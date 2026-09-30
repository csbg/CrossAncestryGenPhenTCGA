## Review =====================================================================

## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
    library(ComplexHeatmap)
    library(data.table)
    library(ggnewscale)
    library(patchwork)
    library(yardstick)
    library(circlize)
    library(ggplot2)
    library(msigdbr)
    library(fgsea)
    library(ggh4x)
  }
)

## Script fun. ================================================================

# Make shorter phenotype
vs_newline <- function(x) {
  gsub("\\s*vs\\s*", " vs\n", x)
}

rename_phenotype <- function(x) {
  x <- gsub("\\bLumA\\b", "Luminal A", x)
  x <- gsub("\\bLumB\\b", "Luminal B", x)
  x
}

# Save cor
cor_safe <- function(x, y, n) {

  if (
    n >= 3 &&
    sd(x) > 0 &&
    sd(y) > 0
  ) {
    cor(x, y, method = "pearson")
  } else {
    NA_real_
  }
}

# Save heatmap
save_heatmapsDK <- function(
  plot,
  file,
  height,
  width,
  bg = "transparent",
  dpi = NA,
  unit = "cm",
  ncol = 2
) {
  
  # helper: convert to inches if needed
  to_in <- function(x, unit) {
    if (unit == "cm") return(x / 2.54)
    if (unit == "mm") return(x / 25.4)
    if (unit == "in") return(x)
    stop("Unsupported unit: ", unit)
  }
  
  width_in  <- to_in(width, unit)
  height_in <- to_in(height, unit)
  
  # infer format
  ext <- tools::file_ext(file)
  
  # open device
  if (ext == "svg") {
    svg(file, width = width_in, height = height_in)
    
  } else if (ext == "png") {
    png(
      file,
      width = width_in,
      height = height_in,
      units = "in",
      res = ifelse(is.na(dpi), 300, dpi),
      bg = bg
    )
    
  } else {
    stop(paste("Unsupported format:", ext))
  }
  
  heatmaps <- plot
  n <- length(heatmaps)
  nrow <- ceiling(n / ncol)
  
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  
  for (i in seq_along(heatmaps)) {
    row <- ceiling(i / ncol)
    col <- i %% ncol
    if (col == 0) col <- ncol
    
    pushViewport(
      viewport(
        layout.pos.row = row,
        layout.pos.col = col
      )
    )

    draw(
      heatmaps[[i]], 
      merge_legends = FALSE,              
      legend_grouping = "original",
        
      heatmap_legend_side = "right",       
      annotation_legend_side = "bottom", 
        
      newpage = FALSE
    )

    grid::grid.text(
      letters[i],
      x = unit(2, "mm"),
      y = unit(1, "npc") - unit(2, "mm"),
      just = c("left", "top"),
      gp = gpar(fontsize = 10, fontface = "bold")
    )
    upViewport()
  }
  
  dev.off()
}

## Colors =====================================================================

ancestry_cols <- c(
  "EUR"   = "#0072B2",
  "AFR"   = "#D55E00",
  "EAS"   = "#56B4E9",
  "AMR"   = "#E69F00",
  "SAS"   = "#009E73",
  "ADMIX" = "#999999"
)

phenotype_cols <- c(
  "Basal vs non-Basal"      = "#c67aa3", 
  "Luminal A vs Luminal B"  = "darkred", 
  "Normal vs Primary"       = "#A65628", 
  "Serous vs Endometrioid"  = "#414141", 
  "Classical vs Follicular" = "#1B9E77", 
  "M0 vs MX"                = "#d9c002"
)

## Directories ================================================================

# Result directories (per cance study)
res_dirs <- list(
  BRCA = file.path("results", "tcga", "analysis", "TCGA_BRCA"),
  UCEC = file.path("results", "tcga", "analysis", "TCGA_UCEC"),
  THCA = file.path("results", "tcga", "analysis", "TCGA_THCA")
)

# Figures
fig_dir <- file.path("results_major_review", "tcga", "figures", "subset_prediction_effect")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Tables
tab_dir <- file.path("results_major_review", "tcga", "tables", "subset_prediction_effect")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

## Results ====================================================================

# AUC
auc_stats <- rbindlist(lapply(names(res_dirs), function(study) {

    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    # List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$summary_stats)
      dt[, tech := t]
    }, files, techs)

  # Add phenotype and study
  out <- rbindlist(dts)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

# Loss
pred_loss <- rbindlist(lapply(names(res_dirs), function(study) {

    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    # List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$subsets_stats)
      dt[, tech := t]
    }, files, techs)

    # Compute log losses → fractional summaries → bind
    ll <- lapply(dts, function(dt)
      dt[, .(logloss = mn_log_loss_vec(
        truth = droplevels(true), estimate = prob, event_level = "second")),
        by = .(iteration, coef_id, g_1, g_2, a_1, a_2, tech)]
    )

    # Summaries
    fr <- lapply(ll, function(dt) {
      meta <- unique(dt[, .(g_1, g_2, a_1, a_2, tech)])
      w <- dcast(dt, iteration ~ coef_id, value.var = "logloss")
      w[, frac := relationship_Y / relationship_X]

      # Add p-values
      B <- nrow(w)                                    # nr. of EUR-subsets
      p_left  <- (1 + sum(w$frac <= 1)) / (B + 1)     # tests frac > 1
      p_right <- (1 + sum(w$frac >= 1)) / (B + 1)     # tests frac < 1
      p_two_sided <- min(1, 2 * min(p_left, p_right)) # symmetric test

      sm <- w[, .(
          frac_mean = mean(frac),
          frac_q025 = quantile(frac, 0.025),
          frac_q975 = quantile(frac, 0.975),
          p_value   = p_two_sided
      )]

      cbind(meta, sm)
    })

  # Add phenotype and study
  out <- rbindlist(fr)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

# Probabilities
pred_prob <- rbindlist(lapply(names(res_dirs), function(study) {
    ## Files 
    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    # List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$subsets_stats)
      dt[, tech := t]
    }, files, techs)

    # Output
    out <- rbindlist(dts)
    out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
    out[, study := study] 
    out
}), use.names = TRUE, fill = TRUE)

# Feature importance
pred_feat <- rbindlist(lapply(names(res_dirs), function(study) {

  fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
  files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  techs <- sub("_.*$", "", basename(dirname(files)))

  # List of subset stats + tech
  dts <- Map(function(f,t) {
    dt <- as.data.table(readRDS(f)$feature_stats)
    dt[, tech := t]
    dt[, l1_norm := estimate / sum(abs(estimate)), by = .(iteration)]
    dt
  }, files, techs)

  # Add phenotype and study
  out <- rbindlist(dts)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

# Interaction effect
subset_res <- rbindlist(lapply(names(res_dirs), function(study) {
  file_path <- file.path(res_dirs[[study]], "subset_interaction_effect")
  dge_files <- list.files(file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  # Load all files
  all_files <- lapply(dge_files, readRDS)
  dt <- rbindlist(lapply(all_files, \(x) x$summary_stats), use.names = TRUE, fill = TRUE)
  # Set phenotype
  dt[, phenotype := fcase(
    coef_id == "baseline_1", paste0(g_1),
    coef_id == "baseline_2", paste0(g_2),
    coef_id == "relationship_1", paste(g_1, "vs", g_2),
    coef_id == "relationship_2", paste(g_1, "vs", g_2),
    coef_id == "interaction", paste(g_1, "vs", g_2)
  )]
  dt[, phenotype := gsub("_", "-", phenotype)]
  dt
}), fill = TRUE)


alpha <- 0.1
formats <- c("svg", "png")

## Reviewer 1 =================================================================
### R1.07 ---------------------------------------------------------------------
res_dirs_major_review <- list(
  BRCA = file.path("results_major_review", "tcga", "analysis", "TCGA_BRCA"),
  UCEC = file.path("results_major_review", "tcga", "analysis", "TCGA_UCEC"),
  THCA = file.path("results_major_review", "tcga", "analysis", "TCGA_THCA")
)



### R1.11 ---------------------------------------------------------------------
#### Panel A (Correlation of raw weights vs. L1 norm. weights) ----------------
settings <- list(
  title_gp    = gpar(fontsize = 6),
  labels_gp   = gpar(fontsize = 6),
  grid_width  = unit(0.2, "cm"), 
  grid_height = unit(0.2, "cm")
)

techs <- c("mrna", "meth")

r1.11_panel_A <- list()
for (t in techs) {
  
  # Long
  r1.11_df <- pred_feat[
    study == "THCA" & 
    tech == t
  ]
  
  r1.11_df[, model_id := paste(study, tech, phenotype, a_2, paste0("Iter", iteration), sep = ":")]
  
  # Wide
  r1.11_df_wide <- dcast(
    r1.11_df, 
    feature ~ model_id, 
    value.var = c("estimate", "l1_norm"),
    sep = ":"
  )
  r1.11_df_wide <- na.omit(r1.11_df_wide)
  
  # Matrix
  r1.11_mat <- as.matrix(r1.11_df_wide[, -1])
  r1.11_cor_mat <- cor(r1.11_mat, method = "spearman")
  diag(r1.11_cor_mat) <- NA

  # Meta
  r1.11_meta <- data.frame(model_id = colnames(r1.11_cor_mat), stringsAsFactors = FALSE)
  r1.11_meta$metric    <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[1])
  r1.11_meta$metric    <- ifelse(r1.11_meta$metric == "estimate", "Raw", "L1 norm.")
  r1.11_meta$study     <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[2])
  r1.11_meta$tech      <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[3])
  r1.11_meta$phenotype <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[4])
  r1.11_meta$ancestry  <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[5])
  r1.11_meta$iteration <- sapply( strsplit(r1.11_meta$model_id, ":"), function(x) x[6])

  weight_cols <- c("Raw" = "#4DAF4A", "L1 norm." = "#984EA3")

  # Top annotation 
  anno_col_top <- HeatmapAnnotation(
    `Model weight` = r1.11_meta$metric,
    Comparison = r1.11_meta$phenotype,
    Ancestry = r1.11_meta$ancestry,
    col = list(`Model weight` = weight_cols, Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(0.2, "cm"), 
    annotation_name_gp = gpar(fontsize = 6), 
    annotation_legend_param = list(`Model weight` = settings, Comparison = settings, Ancestry = settings)
  )

  # Left annotation
  anno_row_left <- rowAnnotation(
    `Model weight` = r1.11_meta$metric,
    Comparison = r1.11_meta$phenotype,
    Ancestry = r1.11_meta$ancestry,
    col = list(`Model weight` = weight_cols, Ancestry = ancestry_cols, Comparison = phenotype_cols),
    simple_anno_size = unit(0.2, "cm"), 
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  # Heatmap 
  r1.11_panel_A[[t]] <- Heatmap(
    r1.11_cor_mat,
    name = "Spearman",

    col    = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    na_col = "gray90",

    column_title = if (t == "mrna") {"THCA Expression"} else if (t == "meth") {"THCA Methylation"},

    column_title_side = "top",
    row_title_side    = "right",
    
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    
    top_annotation  = anno_col_top,
    left_annotation = anno_row_left,
    
    show_row_names    = FALSE,
    show_column_names = FALSE,
    
    heatmap_legend_param = list(
      at = c(-1, 0 ,1),
      labels = c("-1", "0", "1"),
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      legend_height = unit(1, "cm"),
      grid_width = unit(0.2, "cm")
    ),
    
    column_names_gp = gpar(fontsize = 6),
    row_names_gp    = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_title_gp    = gpar(fontsize = 6),

    # Ht dimensions
    width  = unit(4.0, "cm"),
    height = unit(4.0, "cm")
  )
}

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = r1.11_panel_A,
    file = file.path(
      fig_dir, 
      paste0("r1.11_panel_A.", ext)
    ),
    height = 7.5,
    width = 16,
    ncol = 2,
    bg = "white",
    dpi = 300
  )
}

#### Panel B (Paper Figure 4E with L1 norm. weights ) -------------------------
feat <- unique(
  c(
    subset_res[
      study == "BRCA" & 
      tech == "mrna" &
      phenotype == "Basal vs non-Basal" &
      coef_id == "interaction" &
      p_adj < alpha
    ][
      , .(feature = head(feature, 15)), 
      by = a_2
    ][["feature"]],
    pred_feat[
      study == "BRCA" & 
      tech == "mrna" &
      phenotype == "Basal vs non-Basal",
      .(
        mean_estimate = mean(l1_norm, na.rm = TRUE),
        abs_estimate  = abs(mean(l1_norm, na.rm = TRUE))
      ),
      by = .(
        coef_id, coef_type, contrast, g_1, g_2, 
        a_1, a_2, study, tech, phenotype, feature
      )
    ][
      order(-abs_estimate),
      .(feature = head(feature, 20)),
      by = a_2
    ][["feature"]]
  )
)

r1.11 <- wrap_plots(
  ggplot(
    data = pred_feat[
      study == "BRCA" & 
      tech == "mrna" &
      phenotype == "Basal vs non-Basal" &
      feature %in% feat,
      .(mean_estimate = mean(l1_norm, na.rm = TRUE)),
      by = .(
        coef_id, coef_type, contrast, g_1, g_2, 
        a_1, a_2, study, tech, phenotype, feature
      )
    ][
      , feature := factor(
        feature,
        levels = feat
      )
    ][
      , M := max(abs(mean_estimate), na.rm = TRUE)
    ][],
    mapping = aes(
      x = a_2, 
      y = feature, 
      fill = mean_estimate
    )
  ) +
  geom_tile() +
  scale_fill_gradient2(
    name = "L1 normalized\nmodel weights",
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    midpoint = 0,
    limits = {
      dt <- pred_feat[
        study == "BRCA" & 
        tech == "mrna" &
        phenotype == "Basal vs non-Basal" &
        feature %in% feat,
        .(mean_estimate = mean(l1_norm, na.rm = TRUE)),
        by = .(
          coef_id, coef_type, contrast, g_1, g_2, 
          a_1, a_2, study, tech, phenotype, feature
        )
      ]
      M <- max(abs(dt$mean_estimate), na.rm = TRUE)
      c(-M, M)
    }
  ) +
  labs(
    x = "Ancestry",
    y = "Feature"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1,
    rotate = 45, 
    show_borders = TRUE
  ) + 
  theme(
    panel.spacing.y = unit(0.15, "lines"),
    plot.margin = margin(r = 0.15, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(size = 4),
  ),
  ggplot(
    data = subset_res[
      study == "BRCA" & 
      tech == "mrna" &
      phenotype == "Basal vs non-Basal" &
      coef_id %in% c(
          "interaction", 
          "relationship_1" , 
          "relationship_2"
        ) &
      feature %in% feat
    ][
      , feature := factor(
        feature,
        levels = feat
      )
    ][],
    mapping = aes(
      x = a_2, 
      y = feature,
      size = -log10(p_adj),
      color = T_obs
    )
  ) +
  geom_point() +
  scale_color_gradient2(
    name = expression(log[2] *  ( "fold-change" )),
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = expression(-log[10] * ( "adj. p-value" )),
    range = c(0.5, 1.8)
  ) +
  facet_grid(
    cols   = vars(coef_id),
    scales = "free_y",
    space  = "free_y",
    labeller = labeller(
      coef_id = c(
        interaction = "Interaction\neffect",
        relationship_1 = "EUR\ncancer effect",
        relationship_2 = "non-EUR\ncancer effect"
      )
    )
  ) + 
  labs(
    x = "Ancestry",
    y = "Feature"
  ) +
  theme_CrossAncestryGenPhen(
    rotate = 45, 
    legend_key = 1, 
    show_borders = TRUE
  ) + 
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ),
  ncol = 2
) + 
plot_layout(
  guides = "collect",
  widths = c(0.22, 1)
) 

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r1.11,
    file = file.path(
      fig_dir, 
      paste0(
        "r1.11.", 
        ext
      )
    ),
    height = 13,
    width = 10,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}

### R1.14 ---------------------------------------------------------------------
#### Panel A ------------------------------------------------------------------
# Overall robustness of model weights across trained models (Spearman)

settings <- list(
  title_gp    = gpar(fontsize = 6),
  labels_gp   = gpar(fontsize = 6),
  grid_width  = unit(0.2, "cm"), 
  grid_height = unit(0.2, "cm")
)

studies <- c("BRCA", "THCA", "UCEC")
techs   <- c("mrna", "meth") 

r1.14_panel_A <- list()
for (st in studies) {

  tech_heatmaps <- list()
  for (tc in techs) {
    sub_df <- pred_feat[study == st & tech == tc]
  
    if (nrow(sub_df) == 0) next
    
    sub_df$phenotype <- rename_phenotype(sub_df$phenotype)
    sub_df[, model_id := paste(study, tech, phenotype, a_2, paste0("Iter", iteration), sep = "_")]
    
    # Wide format
    df_wide <- dcast(sub_df, feature ~ model_id, value.var = "l1_norm")
    df_wide <- na.omit(df_wide)

    # Matrix
    mat <- as.matrix(df_wide[, -1, with = FALSE])
    rownames(mat) <- df_wide$feature
    
    # Correlation
    cor_mat <- cor(mat, method = "spearman")
    diag(cor_mat) <- NA
    
    # Meta
    meta <- data.frame(model_id = colnames(cor_mat))
    meta$phenotype <- sapply(strsplit(meta$model_id, "_"), function(x) x[3])
    meta$ancestry  <- sapply(strsplit(meta$model_id, "_"), function(x) x[4])
    
    # Top annotation 
    if (tc == techs[1]) {
      anno_col_top <- HeatmapAnnotation(
        `Comparison` = meta$phenotype,
        Ancestry = meta$ancestry,
        col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE,
        show_legend = TRUE,
        simple_anno_size = unit(0.2, "cm"), 
        annotation_name_gp = gpar(fontsize = 6), 
        annotation_legend_param = list(`Comparison` = settings, Ancestry = settings)
      )
    } else {
      anno_col_top <- NULL
    }
    
    # Left annotation
    anno_row_left <- rowAnnotation(
      `Comparison` = meta$phenotype,
      Ancestry = meta$ancestry,
      col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
      simple_anno_size = unit(0.2, "cm"), 
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    # Heatmap
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name = "Spearman",

      col    = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
      na_col = "gray90",

      column_title = st,
      row_title    = if (tc == "mrna") {"Expression"} else if (tc == "meth") {"Methylation"},

      column_title_side = "top",
      row_title_side    = "right",
      
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      
      top_annotation  = anno_col_top,
      left_annotation = anno_row_left,
      
      show_row_names    = FALSE,
      show_column_names = FALSE,
      
      heatmap_legend_param = list(
        at = c(-1, 0 ,1),
        labels = c("-1", "0", "1"),
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        legend_height = unit(1, "cm"),
        grid_width = unit(0.2, "cm")
      ),
      
      column_names_gp = gpar(fontsize = 6),
      row_names_gp    = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 6),
      row_title_gp    = gpar(fontsize = 6),

      # Ht dimensions
      width  = unit(2.0, "cm"),
      height = unit(2.0, "cm")
    )
    
    tech_heatmaps[[tc]] <- ht
  }
  
  r1.14_panel_A[[st]] <- Reduce(`%v%`, tech_heatmaps)
}

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = r1.14_panel_A,
    file = file.path(
      fig_dir, 
      paste0("r1.14_panel_A.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 3,
    bg = "white",
    dpi = 300
  )
}

#### Panel B ------------------------------------------------------------------
# Overall robustness of model weights across trained models (Pearson)

settings <- list(
  title_gp    = gpar(fontsize = 6),
  labels_gp   = gpar(fontsize = 6),
  grid_width  = unit(0.2, "cm"), 
  grid_height = unit(0.2, "cm")
)

studies <- c("BRCA", "THCA", "UCEC")
techs   <- c("mrna", "meth") 

r1.14_panel_B <- list()
for (st in studies) {

  tech_heatmaps <- list()
  for (tc in techs) {
    sub_df <- pred_feat[study == st & tech == tc]
  
    if (nrow(sub_df) == 0) next
    
    sub_df$phenotype <- rename_phenotype(sub_df$phenotype)
    sub_df[, model_id := paste(study, tech, phenotype, a_2, paste0("Iter", iteration), sep = "_")]
    
    # Wide format
    df_wide <- dcast(sub_df, feature ~ model_id, value.var = "l1_norm")
    df_wide <- na.omit(df_wide)

    mat <- as.matrix(df_wide[, -1, with = FALSE])
    rownames(mat) <- df_wide$feature
    
    # Correlation
    cor_mat <- cor(mat, method = "pearson")
    diag(cor_mat) <- NA
    
    # Meta
    meta <- data.frame(model_id = colnames(cor_mat))
    meta$phenotype <- sapply(strsplit(meta$model_id, "_"), function(x) x[3])
    meta$ancestry  <- sapply(strsplit(meta$model_id, "_"), function(x) x[4])
    
    # Top annotation 
    if (tc == techs[1]) {
      anno_col_top <- HeatmapAnnotation(
        `Comparison` = meta$phenotype,
        Ancestry = meta$ancestry,
        col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE,
        show_legend = TRUE,
        simple_anno_size = unit(0.2, "cm"), 
        annotation_name_gp = gpar(fontsize = 6), 
        annotation_legend_param = list(`Comparison` = settings, Ancestry = settings)
      )
    } else {
      anno_col_top <- NULL
    }
    
    # Left annotation
    anno_row_left <- rowAnnotation(
      `Comparison` = meta$phenotype,
      Ancestry = meta$ancestry,
      col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
      simple_anno_size = unit(0.2, "cm"), 
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    # Heatmap
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name = "Pearson",

      col    = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
      na_col = "gray90",

      column_title = st,
      row_title    = if (tc == "mrna") {"Expression"} else if (tc == "meth") {"Methylation"},

      column_title_side = "top",
      row_title_side    = "right",
      
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      
      top_annotation  = anno_col_top,
      left_annotation = anno_row_left,
      
      show_row_names    = FALSE,
      show_column_names = FALSE,
      
      heatmap_legend_param = list(
        at = c(-1, 0 ,1),
        labels = c("-1", "0", "1"),
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        legend_height = unit(1, "cm"),
        grid_width = unit(0.2, "cm")
      ),
      
      column_names_gp = gpar(fontsize = 6),
      row_names_gp    = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 6),
      row_title_gp    = gpar(fontsize = 6),

      # Ht dimensions
      width  = unit(2.0, "cm"),
      height = unit(2.0, "cm")
    )
    
    tech_heatmaps[[tc]] <- ht
  }
  
  r1.14_panel_B[[st]] <- Reduce(`%v%`, tech_heatmaps)
}

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = r1.14_panel_B,
    file = file.path(
      fig_dir, 
      paste0("r1.14_panel_B.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 3,
    bg = "white",
    dpi = 300
  )
}

#### Panel C ------------------------------------------------------------------
# Checking consistency of the top 100 features across models and iterations

settings <- list(
  title_gp    = gpar(fontsize = 6),
  labels_gp   = gpar(fontsize = 6),
  grid_width  = unit(0.2, "cm"), 
  grid_height = unit(0.2, "cm")
)

studies <- c("BRCA", "THCA", "UCEC")
techs   <- c("mrna", "meth") 

r1.14_panel_C <- list()
for (st in studies) {

  tech_heatmaps <- list()
  for (tc in techs) {
    sub_df <- pred_feat[study == st & tech == tc]
  
    if (nrow(sub_df) == 0) next
    
    sub_df$phenotype <- rename_phenotype(sub_df$phenotype)
    sub_df[, model_id := paste(study, tech, phenotype, a_2, paste0("Iter", iteration), sep = "_")]
    
    # Top 100 features per model
    model_names  <- unique(sub_df$model_id)
    top_100_list <- list()
    
    for (mod in model_names) {
      mod_data            <- sub_df[model_id == mod][order(-abs(l1_norm))]
      top_features        <- mod_data[l1_norm != 0, feature][1:100]
      top_100_list[[mod]] <- na.omit(top_features)
    }
    
    # Jaccard matrix
    n_models <- length(model_names)
    jaccard_mat <- matrix(NA, nrow = n_models, ncol = n_models, dimnames = list(model_names, model_names))
    
    for (i in 1:n_models) {
      for (j in 1:n_models) {
        if (i == j) {
          jaccard_mat[i, j] <- NA 
        } else {
          set_A <- top_100_list[[model_names[i]]]
          set_B <- top_100_list[[model_names[j]]]
          
          intersection <- length(intersect(set_A, set_B))
          union_set    <- length(union(set_A, set_B))
          
          jaccard_mat[i, j] <- if (union_set > 0) intersection / union_set else 0
        }
      }
    }
    
    # Meta
    meta <- data.frame(model_id = colnames(jaccard_mat))
    meta$phenotype <- sapply(strsplit(meta$model_id, "_"), function(x) x[3])
    meta$ancestry  <- sapply(strsplit(meta$model_id, "_"), function(x) x[4])
    
    # Top annotation 
    if (tc == techs[1]) {
      anno_col_top <- HeatmapAnnotation(
        `Comparison` = meta$phenotype,
        Ancestry     = meta$ancestry,
        col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE,
        show_legend          = TRUE,
        simple_anno_size     = unit(0.2, "cm"), 
        annotation_name_gp   = gpar(fontsize = 6), 
        annotation_legend_param = list(`Comparison` = settings, Ancestry = settings)
      )
    } else {
      anno_col_top <- NULL
    }
    
    # Left annotation
    anno_row_left <- rowAnnotation(
      `Comparison` = meta$phenotype,
      Ancestry = meta$ancestry,
      col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
      simple_anno_size = unit(0.2, "cm"), 
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    # Heatmap
    ht <- ComplexHeatmap::Heatmap(
      jaccard_mat,
      name = "Jaccard\ntop 100",
      
      col    = colorRamp2(c(0, 0.5, 1), c("white", "orange", "red")),
      na_col = "gray90",

      column_title = st,
      row_title    = if (tc == "mrna") {"Expression"} else if (tc == "meth") {"Methylation"},

      column_title_side = "top",
      row_title_side    = "right",
      
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      
      top_annotation  = anno_col_top,
      left_annotation = anno_row_left,
      
      show_row_names    = FALSE,
      show_column_names = FALSE,
      
      heatmap_legend_param = list(
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        legend_height = unit(1, "cm"),
        grid_width = unit(0.2, "cm")
      ),
      
      column_names_gp = gpar(fontsize = 6),
      row_names_gp    = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 6),
      row_title_gp    = gpar(fontsize = 6),

      width  = unit(2.0, "cm"),
      height = unit(2.0, "cm")
    )
    
    tech_heatmaps[[tc]] <- ht
  }
  
  r1.14_panel_C[[st]] <- Reduce(`%v%`, tech_heatmaps)
}

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = r1.14_panel_C,
    file = file.path(
      fig_dir, 
      paste0("r1.14_panel_C.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 3,
    bg = "white",
    dpi = 300
  )
}

#### Panel D ------------------------------------------------------------------
# Checking consistency of the top 50 features across models

settings <- list(
  title_gp    = gpar(fontsize = 6),
  labels_gp   = gpar(fontsize = 6),
  grid_width  = unit(0.2, "cm"), 
  grid_height = unit(0.2, "cm")
)

studies <- c("BRCA", "THCA", "UCEC")
techs   <- c("mrna", "meth") 

r1.14_panel_D <- list()
for (st in studies) {

  tech_heatmaps <- list()
  for (tc in techs) {
    sub_df <- pred_feat[study == st & tech == tc]
  
    if (nrow(sub_df) == 0) next
    
    sub_df$phenotype <- rename_phenotype(sub_df$phenotype)
    sub_df[, model_id := paste(study, tech, phenotype, a_2, paste0("Iter", iteration), sep = "_")]
    
    # Top 50 features per model
    model_names  <- unique(sub_df$model_id)
    top_100_list <- list()
    
    for (mod in model_names) {
      mod_data            <- sub_df[model_id == mod][order(-abs(l1_norm))]
      top_features        <- mod_data[l1_norm != 0, feature][1:50]
      top_100_list[[mod]] <- na.omit(top_features)
    }
    
    # Jaccard matrix
    n_models <- length(model_names)
    jaccard_mat <- matrix(NA, nrow = n_models, ncol = n_models, dimnames = list(model_names, model_names))
    
    for (i in 1:n_models) {
      for (j in 1:n_models) {
        if (i == j) {
          jaccard_mat[i, j] <- NA 
        } else {
          set_A <- top_100_list[[model_names[i]]]
          set_B <- top_100_list[[model_names[j]]]
          
          intersection <- length(intersect(set_A, set_B))
          union_set    <- length(union(set_A, set_B))
          
          jaccard_mat[i, j] <- if (union_set > 0) intersection / union_set else 0
        }
      }
    }
    
    # Meta
    meta <- data.frame(model_id = colnames(jaccard_mat))
    meta$phenotype <- sapply(strsplit(meta$model_id, "_"), function(x) x[3])
    meta$ancestry  <- sapply(strsplit(meta$model_id, "_"), function(x) x[4])
    
    # Top annotation
    if (tc == techs[1]) {
      anno_col_top <- HeatmapAnnotation(
        `Comparison` = meta$phenotype,
        Ancestry = meta$ancestry,
        col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE,
        show_legend = TRUE,
        simple_anno_size = unit(0.2, "cm"), 
        annotation_name_gp = gpar(fontsize = 6), 
        annotation_legend_param = list(`Comparison` = settings, Ancestry = settings)
      )
    } else {
      anno_col_top <- NULL
    }
    
    # Left annotation
    anno_row_left <- rowAnnotation(
      `Comparison` = meta$phenotype,
      Ancestry = meta$ancestry,
      col = list(Ancestry = ancestry_cols, `Comparison` = phenotype_cols),
      simple_anno_size = unit(0.2, "cm"), 
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    # Heatmap
    ht <- ComplexHeatmap::Heatmap(
      jaccard_mat,
      name = "Jaccard\ntop 50",
      
      col    = colorRamp2(c(0, 0.5, 1), c("white", "orange", "red")),
      na_col = "gray90",

      column_title = st,
      row_title    = if (tc == "mrna") {"Expression"} else if (tc == "meth") {"Methylation"},

      column_title_side = "top",
      row_title_side    = "right",
      
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      
      top_annotation  = anno_col_top,
      left_annotation = anno_row_left,
      
      show_row_names    = FALSE,
      show_column_names = FALSE,
      
      heatmap_legend_param = list(
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        legend_height = unit(1, "cm"),
        grid_width = unit(0.2, "cm")
      ),
      
      column_names_gp = gpar(fontsize = 6),
      row_names_gp    = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 6),
      row_title_gp    = gpar(fontsize = 6),

      width  = unit(2.0, "cm"),
      height = unit(2.0, "cm")
    )
    
    tech_heatmaps[[tc]] <- ht
  }
  
  r1.14_panel_D[[st]] <- Reduce(`%v%`, tech_heatmaps)
}

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = r1.14_panel_D,
    file = file.path(
      fig_dir, 
      paste0("r1.14_panel_D.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 3,
    bg = "white",
    dpi = 300
  )
}

#### Panel E (Histogram of model weights for all ancestry-specific genes) ------
# Robustness of ancestry-specific model weights (across significant genes)

sig_features_dt <- subset_res[
  coef_id == "interaction" & p_adj < alpha, 
  .(study, tech, phenotype, a_2, feature)
]

density_sig_df <- merge(
  pred_feat, 
  sig_features_dt, 
  by = c("study", "tech", "phenotype", "a_2", "feature")
)

r1.14_panel_E <- ggplot(
  data = density_sig_df, 
  mapping = aes(
    x = l1_norm,
    y = after_stat(count / tapply(count, PANEL, sum)[PANEL]),
    fill = phenotype
  )
) +
geom_histogram(
) +
scale_fill_manual(
  values = phenotype_cols
) +
scale_x_continuous(
  # limits = c(-1, 1),
  # breaks = c(-1, 0, 1),
  breaks = scales::pretty_breaks(n = 3)
) +
scale_y_continuous(
  limits = c(0, 1),
  breaks = scales::pretty_breaks(n = 3),
  labels = scales::label_percent()
) +
facet_grid(
  cols = vars(tech),
  rows = vars(a_2),
  labeller = labeller(
    tech = c(
      meth = "Methylation",
      mrna = "Expression"
    )
  )
) +
labs(
  x = "L1 normalized weights", 
  y = "Pct. ancestry-specific features",
  fill = "Comparison"
) + 
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_border = TRUE,
  rotate = 45
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines")
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r1.14_panel_E,
    file = file.path(
      fig_dir, 
      paste0(
        "r1.14_panel_E.", 
        ext
      )
    ),
    height = 9,
    width = 16,
    trimmed = FALSE,
    bg = "white",
    dpi = 300
  )
}

# studies <- c("BRCA", "THCA", "UCEC")
# techs   <- c("mrna", "meth")

# plot_list <- list()
# for (t in techs) {
#   for (s in studies) {
    
#     sub_df   <- density_sig_df[study == s & tech == t]
#     plot_key <- paste(t, s, sep = "_")

#     # Title
#     clean_tech <- t
#     if (t == "mrna") {
#       clean_tech <- "Expression"
#     } else if (t == "meth") {
#       clean_tech <- "Methylation"
#     }
#     plot_title <- paste(s, clean_tech, sep = " ")
    
#     # Plot
#     p <- ggplot(
#       data = sub_df, 
#       mapping = aes(
#         x = l1_norm
#       )
#     ) +
#     geom_histogram(
#       mapping = aes(y = after_stat(count / tapply(count, PANEL, sum)[PANEL])),
#       fill = "red"
#     ) +
#     facet_grid(
#       rows = vars(a_2), 
#       cols = vars(vs_newline(phenotype))
#     ) +
#     scale_x_continuous(
#       # limits = c(-1, 1),
#       # breaks = c(-1, 0, 1),
#       breaks = scales::pretty_breaks(n = 3)
#     ) +
#     scale_y_continuous(
#       limits = c(0, 1),
#       breaks = scales::pretty_breaks(n = 3),
#       labels = scales::label_percent()
#     ) +
#     labs(
#       title = plot_title,
#       x = "L1 normalized weights", 
#       y = "Pct. ancestry-specific features"
#     ) + 
#     theme_CrossAncestryGenPhen(
#       legend_key = 1,
#       show_border = TRUE,
#       rotate = 45
#     ) +
#     theme(
#       panel.spacing.x = unit(0.15, "lines"),
#       panel.spacing.y = unit(0.15, "lines")
#     ) +
#     theme(
#       plot.margin = margin(0, 0, 0, 0, unit = "lines")
#     )
      
#     plot_list[[plot_key]] <- p
#   }
# }

# # Patchwork
# ordered_plots <- list()
# for (t in techs) {
#   for (s in studies) {
#     ordered_plots <- c(
#       ordered_plots, 
#       list(plot_list[[paste(t, s, sep = "_")]])
#     )
#   }
# }

# r1.14_panel_E <- wrap_plots(ordered_plots, ncol = length(studies)) +
#  plot_annotation(
#     tag_levels = "a"
#   ) &
#   theme(
#     plot.title = element_text(margin = margin(b = -3)),
#     plot.tag = element_text(size = 10, face = "bold")
#   )

# # Save
# for (ext in formats) {
#   ggsaveDK(
#     plot = r1.14_panel_E,
#     file = file.path(
#       fig_dir, 
#       paste0(
#         "r1.14_panel_E.", 
#         ext
#       )
#     ),
#     height = 12,
#     width = 16,
#     trimmed = FALSE,
#     bg = "white",
#     dpi = 300
#   )
# }

#### Panel F (Boxplots of model weights for top 5 ancestry-specific genes) -----
# Robustness of ancestry-specific model weights (individual genes)
techs <- c("mrna", "meth")

tech_plots <- list()
for (tc in techs) {
  
  top_genes_dt <- subset_res[
    study == "BRCA" & 
    tech == tc &
    phenotype == "Basal vs non-Basal" &
    coef_id == "interaction" &
    p_adj < alpha
  ][order(p_adj), .(feature = head(feature, 10)), by = .(a_2)]
  
  unique_feat <- unique(top_genes_dt$feature)
  
  bg_significance <- subset_res[
    study == "BRCA" & 
    tech == tc &
    phenotype == "Basal vs non-Basal" &
    coef_id == "interaction" &
    feature %in% unique_feat
  , .(
      is_significant = any(p_adj < alpha)
    ), by = .(a_2, feature)]
  
  bg_significance <- bg_significance[is_significant == TRUE]
  
  density_df <- pred_feat[
    study == "BRCA" & 
    tech == tc &
    phenotype == "Basal vs non-Basal" &
    feature %in% unique_feat
  ]

  density_df[, feature := factor(feature, levels = unique_feat)]
  bg_significance[, feature := factor(feature, levels = unique_feat)]
  
  clean_label <- ifelse(tc == "mrna", "Expression", "Methylation")
  
  p <- ggplot() +
    geom_rect(
      data = bg_significance,
      mapping = aes(
        xmin = as.numeric(feature) - 0.5, 
        xmax = as.numeric(feature) + 0.5, 
        ymin = -Inf, 
        ymax = Inf,
        fill = "Ancestry-specific (FDR < 0.1)"
      ),
      alpha = 0.8,
      inherit.aes = FALSE  
    ) +
    geom_boxplot(
      data = density_df,
      mapping = aes(
        x = feature,
        y = l1_norm
      ),
      color = "black",
      linewidth = 0.1,
      outlier.shape = NA
    ) +
    facet_grid(
      rows = vars(a_2)
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("Ancestry-specific (FDR < 0.1)" = "#FFECEC")
    ) +
    labs(
      title = paste("BRCA", clean_label, "Basal vs non-Basal"),
      x = if (tc == "mrna") {"Genes"} else if (tc == "meth") {"Methylation sites"},
      y = "L1 normalized weights"
    ) +
    theme_CrossAncestryGenPhen(
      legend_key = 1,
      show_border = TRUE,
      rotate = 45
    )
  
  tech_plots[[tc]] <- p
}

r1.14_panel_F <- wrap_plots(tech_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 10, face = "bold")
  )

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r1.14_panel_F,
    file = file.path(
      fig_dir, 
      paste0(
        "r1.14_panel_F.", 
        ext
      )
    ),
    height = 9,
    width = 16,
    trimmed = FALSE,
    bg = "white",
    dpi = 300
  )
}


## Reviewer 2 =================================================================
### R2.11 ---------------------------------------------------------------------
#### Panel A (Inconclusive prediction tasks & sample nr.) ---------------------
# Calc. AUC
pred_AUC <- pred_prob[
  ,
  {
    # Observed AUC
    AUC_obs <- roc_auc_vec(
      truth = droplevels(true),
      estimate = prob,
      event_level = "second"
    )

    # Permuted AUCs
    AUC_perm <- replicate(
      100,
      roc_auc_vec(
        truth = sample(droplevels(true)),
        estimate = prob,
        event_level = "second"
      )
    )

    rbindlist(list(
      data.table(type = "Observed", AUC = AUC_obs),
      data.table(type = "Permuted", AUC = AUC_perm)
    ))
  },
  by = .(iteration, coef_id, coef_type, a_1, a_2, g_1, g_2, study, tech, phenotype)
][
  ,
  .(
    AUC_mean = mean(AUC),
    CI_lower = quantile(AUC, 0.025, na.rm = TRUE),
    CI_upper = quantile(AUC, 0.975, na.rm = TRUE),
    n_iter   = uniqueN(iteration)
  ),
  by = .(type, coef_id, coef_type, a_1, a_2, g_1, g_2, study, tech, phenotype)
][
  ,
  `:=`(
    ancestry = fifelse(
      coef_id == "relationship_X", a_1,
      fifelse(coef_id == "relationship_Y", a_2, NA_character_)
    )
  )
][
  ancestry == "EUR", ancestry := "subset-EUR"
][
  ,
  `:=`(
    ancestry = factor(ancestry, levels = c("subset-EUR", setdiff(levels(ancestry), "subset-EUR"))),
    set = fifelse(ancestry == "subset-EUR", "Validation (EUR)", "Inference (non-EUR)")
  )
][
  ,
  set := factor(set, levels = c("Inference (non-EUR)", "Validation (EUR)"))
][
  ,
  phenotype_plot := vs_newline(
    rename_phenotype(phenotype)
  )
][
  ,
  phenotype_plot := factor(
    phenotype_plot,
    levels = c(
      "Normal vs\nPrimary",
      "Basal vs\nnon-Basal",
      "Luminal A vs\nLuminal B"
    )
  )
]

# Sample sizes
sample_n_table <- fread("results/tcga/analysis/summary_subset_prediction_effect/summary_sample_n.csv")
sample_n_table[, phenotype := sub("^[a-zA-Z0-9]+_(.*)_EUR_vs_.*$", "\\1", comp)]
sample_n_table[, phenotype := gsub("non_", "non-", phenotype)]
sample_n_table[, phenotype := gsub("_", " ", phenotype)]

sample_n_table[, g_col := gsub("_", "-", g_col)]
sample_n_table[g_col == "LumA", g_col := "Luminal A"]
sample_n_table[g_col == "LumB", g_col := "Luminal B"]

sample_n_table[, phenotype := vs_newline(rename_phenotype(phenotype))]
sample_n_table[, phenotype := factor(phenotype, levels = c("Normal vs\nPrimary", "Basal vs\nnon-Basal", "Luminal A vs\nLuminal B"))]

# Plots
techs <- c("meth", "mrna")

tech_combined_plots <- list()
for (tc in techs) {

  # Common sample-size axis across technologies
  sample_n_max <- max(
    sample_n_table[
      study == "BRCA" &
      a_col %in% unique(pred_AUC$a_2),
      n
    ],
    na.rm = TRUE
  )
  sample_n_limits <- c(0, sample_n_max * 1.15)
  
  # Calculate performance overlaps (obs. versus perm.)
  perf_overlap <- merge(
    pred_AUC[study == "BRCA" & tech == tc & type == "Observed", .(a_2, set, phenotype_plot, tech, AUC_mean)],
    pred_AUC[study == "BRCA" & tech == tc & type == "Permuted", .(a_2, set, phenotype_plot, tech, CI_upper)],
    by = c("a_2", "set", "phenotype_plot", "tech")
  )
  
  # Combine conditions: Inconclusive if error bar cross OR any group sample size < 5
  inconclusive_c1 <- perf_overlap[CI_upper >= AUC_mean]
  inconclusive_c1[, a_2 := factor(a_2, levels = levels(factor(pred_AUC$a_2)))]
  
  # Format background table for sample sizes (C.1)
  inconclusive_c2 <- copy(inconclusive_c1)
  setnames(inconclusive_c2, old = c("a_2", "phenotype_plot"), new = c("a_col", "phenotype"))
  inconclusive_c2[, a_col := factor(a_col, levels = levels(factor(sample_n_table$a_col)))]

  # Identify observed bars that are inconclusive
  inconclusive_bars <- merge(
    pred_AUC[
      study == "BRCA" &
        tech == tc &
        type == "Observed",
      .(a_2, set, phenotype_plot, AUC_mean)
    ],
    inconclusive_c1[
      , .(a_2, phenotype_plot)
    ],
    by = c("a_2", "phenotype_plot"),
    allow.cartesian = TRUE
  )

  inconclusive_bars[, inconclusive_fill := fifelse(
    set == "Validation (EUR)",
    "Inconclusive Validation",
    "Inconclusive Inference"
  )]

  # Plot model performance
  p_c1 <- ggplot(
    mapping = aes(x = a_2)
  ) +
  geom_col(
    data = pred_AUC[
      study == "BRCA" &
        tech == tc &
        type == "Observed"
    ],
    mapping = aes(
      y = AUC_mean,
      fill = set
    ),
    position = position_dodge(width = 0.9),
    width = 0.9,
    linewidth = 0
  ) +
  geom_col(
    data = inconclusive_bars,
    mapping = aes(
      y = AUC_mean,
      group = set,
      fill = inconclusive_fill,
      alpha = "Inconclusive"
    ),
    position = position_dodge(width = 0.9),
    width = 0.9,
    color = NA
  ) +
  geom_errorbar(
    data = pred_AUC[
      study == "BRCA" &
      tech == tc &
      type == "Permuted"
    ],
    mapping = aes(
      ymin = CI_lower,
      ymax = CI_upper,
      group = set,
      alpha = "Permuted AUC"
    ),
    position = position_dodge(width = 0.9),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_point(
    data = pred_AUC[
      study == "BRCA" &
      tech == tc &
      type == "Permuted"
    ],
    mapping = aes(
      y = AUC_mean,
      group = set,
      alpha = "Permuted AUC"
    ),
    position = position_dodge(width = 0.9),
    shape = 18,
    size = 0.8
  ) +
  scale_fill_manual(
    name = "Performance",
    values = c(
      "Inference (non-EUR)" = "orange",
      "Validation (EUR)" = "#0072B2",
      "Inconclusive Validation" = "grey65",
      "Inconclusive Inference" = "grey85"
    ),
    breaks = c(
      "Validation (EUR)",
      "Inference (non-EUR)"
    ),
    labels = c(
      "Validation (subset-EUR)",
      "Inference (non-EUR)"
    ),
    guide = guide_legend(order = 1)
  ) +
  scale_alpha_manual(
    name = NULL,
    values = c(
      "Permuted AUC" = 1,
      "Inconclusive" = 1
    ),
    breaks = c(
      "Permuted AUC",
      "Inconclusive"
    ),
    guide = guide_legend(
      order = 2,
      override.aes = list(
        shape = c(18, 22),
        fill = c("black", "grey85"),
        color = c("black", "grey85"),
        alpha = 1
      )
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  facet_grid(
    rows = vars(phenotype_plot),
    scales = "free",
    space = "free",
    labeller = labeller(
      phenotype = vs_newline,
      tech = c(
        meth = "Methylation",
        mrna = "Expression"
      )
    )
  ) +
  coord_flip() +
  labs(
    x = "Ancestry",
    y = "Mean ROC AUC"
  ) +
  theme_CrossAncestryGenPhen(
    # rotate = 45,
    legend_key = 1,
    show_borders = TRUE,
    show_facets = TRUE
  ) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    strip.text.y = element_blank()
  )
  
  # Plot sample sizes
  p_c2 <- ggplot(
    data = sample_n_table[
      study == "BRCA" &
      tech == tc &
      a_col %in% unique(pred_AUC$a_2)
    ],
    mapping = aes(
      x = a_col, 
      y = n, 
      fill = g_col
    )
  ) +
  # geom_rect(
  #   data = inconclusive_c2,
  #   mapping = aes(
  #     xmin = as.numeric(a_col) - 0.47, 
  #     xmax = as.numeric(a_col) + 0.47,
  #     ymin = -Inf, ymax = Inf
  #   ),
  #   fill = "grey90",
  #   linewidth = 0.1,
  #   inherit.aes = FALSE
  # ) +
  geom_col(
    linewidth = 0, 
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(label = n), 
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 1.5, 
    inherit.aes = TRUE
  ) +
  scale_fill_brewer(
    name = "Cancer type", 
    palette = "Set2",
    breaks = c(
      "Normal", 
      "Primary", 
      "Basal", 
      "non-Basal", 
      "Luminal A", 
      "Luminal B"
    )
  ) + 
  facet_grid(
    rows = vars(phenotype), 
    scales = "free", 
    space = "free"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = sample_n_limits,
    breaks = c(0, 50, 125),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    x = "Ancestry", 
    y = "Nr. samples"
  ) +
  theme_CrossAncestryGenPhen(
    # rotate = 45,
    legend_key = 1, 
    show_borders = FALSE
  ) +
  theme(
    panel.spacing.x = unit(0.15, "lines"), 
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin = margin(0, 0, 0, 0), 
    plot.margin = margin(0, 0, 0, 0),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank()
  )

  # Tech title
  p_c3 <- ggplot(
    data = sample_n_table[
      study == "BRCA" &
      tech == tc &
      a_col %in% unique(pred_AUC$a_2)
    ],
    mapping = aes(
      x = a_col, 
      y = n, 
      fill = g_col
    )
  ) +
  facet_grid(
    rows = vars(tech), 
    scales = "free", 
    space = "free",
    labeller = labeller(
      tech  = c(
        meth = "Methylation", 
        mrna = "Expression"
      )
    )
  ) +
  coord_flip(clip = "off") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_CrossAncestryGenPhen(
    show_borders = FALSE,
    show_axis = FALSE
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    strip.text.y = element_text(margin = margin(0, 0, 0, 0)),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
  
  # Append
  tech_combined_plots[[tc]] <- wrap_plots(p_c1, p_c2, p_c3, ncol = 3, widths = c(0.9, 0.75, 0.0001))
}

r2.11 <- wrap_plots(
  tech_combined_plots[[1]] &
  theme(
    plot.margin = margin(l = 0.1, t = 0, b = 0.1, r = 0.1, unit = "lines"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ),
  tech_combined_plots[[2]] & 
  theme(
    plot.margin = margin(l = 0.1, t = 0.1, b = 0, r = 0.1, unit = "lines")
  ),
  ncol = 1
) +
plot_layout(
  guides = "collect"
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r2.11,
    file = file.path(
      fig_dir, 
      paste0("r2.11.", ext)
    ),
    height = 9.5,
    width = 8.2,
    trimmed = FALSE,
    bg = "transparent"
  )
}

### R2.12 ---------------------------------------------------------------------
#### Panel A (Ambiguous legends 4B) -------------------------------------------
r2.12 <- ggplot(
  data = pred_loss[, `:=`(
        T_obs = log2(frac_mean),
        p_adj = p.adjust(p_value, method = "BH")
      )
  ],
  mapping = aes(
    x = a_2,
    y = rename_phenotype(phenotype),
    color = T_obs,
    size = -log10(p_adj)
  )
) +
geom_point(
  #shape = 21,
  # color = "black",
  # stroke = 0.1
) +
scale_color_gradient2(
  name = expression(log[2] ~ (frac("non-EUR"["log loss"], "subset-EUR"["log loss"]))),
  low  = "#4575b4",
  mid  = "white",
  high = "#d73027",
  midpoint = 0
) +
scale_size_continuous(
  name = expression(-log[10] * ( "adj. p-value" )),
  range = c(0.5, 2),
  breaks = scales::pretty_breaks(n = 3)
) +
facet_grid(
  cols = vars(tech),
  rows = vars(study),
  scales = "free_y",
  space  = "free_y",
  labeller = labeller(
    tech  = c(
      meth = "Methylation", 
      mrna = "Expression"
    )
  )
) +
labs(
  x = "Ancestry",
  y = "Cancer comparison"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE,
  show_grid = FALSE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r2.12,
    file = file.path(
      fig_dir, 
      paste0(
        "r2.12.", 
        ext
      )
    ),
    height = 5,
    width = 9,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}

### R2.13 ---------------------------------------------------------------------
#### Panel A (Histogram of log2FC loss ratio) ---------------------------------
r2.13.1 <- ggplot(
  data = pred_loss,
  mapping = aes(
    x = frac_mean,
    fill = a_2
  )
) +
geom_histogram(
  binwidth = 1,
  color = "black",
  linewidth = 0.1
) +
geom_vline(
  xintercept = 5,
  linetype = "dashed",
  linewidth = 0.3
) +
scale_fill_manual(
  values = ancestry_cols
) +
facet_grid(
  rows = vars(tech),
  labeller = ggplot2::labeller(
    tech = c(
      meth = "Methylation",
      mrna = "Expression"
    )
  )
) +
labs(
  x = "Fold-change log loss (non-EUR/EUR)",
  y = "Count (cancer comparisons)",
  fill = "Ancestry"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1
)

# Actual loss values
logloss <- rbindlist(lapply(names(res_dirs), function(study) {

    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    # List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$subsets_stats)
      dt[, tech := t]
    }, files, techs)

    # Compute log losses → fractional summaries → bind
    ll <- lapply(dts, function(dt)
      dt[, .(logloss = mn_log_loss_vec(
        truth = droplevels(true), estimate = prob, event_level = "second")),
        by = .(iteration, coef_id, g_1, g_2, a_1, a_2, tech)]
    )

  # Add phenotype and study
  out <- rbindlist(ll)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

r2.13.2.1 <- ggplot(
  data = logloss[
    study %in% c("BRCA") &
    tech %in% c("mrna") &
    phenotype %in% c("Basal vs non-Basal") &
    a_2 %in% c("AMR")
  ],
  mapping = aes(
    x = a_2,
    y = logloss,
    fill = coef_id
  )
) + 
stat_summary(
  fun = "mean",
  geom = "bar",
  position = position_dodge(width = 0.9)
) +
facet_grid(
  cols = vars(tech),
  rows = vars(phenotype),
  labeller = ggplot2::labeller(
    tech = c(
      meth = "BRCA Methylation",
      mrna = "BRCA Expression"
    ),
    phenotype = vs_newline
  )
) +
scale_fill_manual(
    values = c(
      "relationship_Y" = "orange",
      "relationship_X" = "#0072B2"
    ),
    breaks = c(
      "relationship_X",
      "relationship_Y"
    ),
    labels = c(
      "Validation",
      "Inference"
    ),
) +
labs(
  x = "Ancestry",
  y = "Mean log loss",
  fill = "Performance"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_borders = TRUE,
  plot.margin = margin(0, 0, 0, 0)
)

r2.13.2.2 <- ggplot(
  data = logloss[
    study %in% c("UCEC") &
    tech %in% c("meth") &
    phenotype %in% c("Normal vs Primary") &
    a_2 %in% c("AFR")
  ],
  mapping = aes(
    x = a_2,
    y = logloss,
    fill = coef_id
  )
) + 
stat_summary(
  fun = "mean",
  geom = "bar",
  position = position_dodge(width = 0.9)
) +
scale_fill_manual(
    values = c(
      "relationship_Y" = "orange",
      "relationship_X" = "#0072B2"
    ),
    breaks = c(
      "relationship_X",
      "relationship_Y"
    ),
    labels = c(
      "Validation",
      "Inference"
    ),
) +
facet_grid(
  cols = vars(tech),
  rows = vars(phenotype),
  labeller = ggplot2::labeller(
    tech = c(
      meth = "UCEC Methylation",
      mrna = "UCEC Expression"
    ),
    phenotype = vs_newline
  )
) +
labs(
  x = "Ancestry",
  y = "Mean log loss",
  fill = "Performance"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_borders = TRUE,
  plot.margin = margin(0, 0, 0, 0)
)

r2.13.2.3 <- ggplot(
  data = logloss[
    study %in% c("THCA") &
    tech %in% c("meth") &
    phenotype %in% c("Normal vs Primary") &
    a_2 %in% c("AFR")
  ],
  mapping = aes(
    x = a_2,
    y = logloss,
    fill = coef_id
  )
) + 
stat_summary(
  fun = "mean",
  geom = "bar",
  position = position_dodge(width = 0.9)
) +
scale_fill_manual(
    values = c(
      "relationship_Y" = "orange",
      "relationship_X" = "#0072B2"
    ),
    breaks = c(
      "relationship_X",
      "relationship_Y"
    ),
    labels = c(
      "Validation",
      "Inference"
    ),
) +
facet_grid(
  cols = vars(tech),
  rows = vars(phenotype),
  labeller = ggplot2::labeller(
    tech = c(
      meth = "THCA Methylation",
      mrna = "THCA Expression"
    ),
    phenotype = vs_newline
  )
) +
labs(
  x = "Ancestry",
  y = "Mean log loss",
  fill = "Performance"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_borders = TRUE,
  plot.margin = margin(0, 0, 0, 0)
)

##### Reviewer Figure ---------------------------------------------------------
r2.13_panel_A <- wrap_plots(
  r2.13.1,
  wrap_plots(
    r2.13.2.1,
    r2.13.2.2,
    r2.13.2.3,
    ncol = 1
  ) + 
  plot_layout(guides = "collect"),
  ncol = 2
) + 
plot_layout(
  widths = c(1, 0.25)
) + 
plot_annotation(
  tag_levels = list(c("a", "b", "", ""))
) &
theme(
  plot.tag = element_text(size = 10, face = "bold")
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r2.13_panel_A,
    file = file.path(
      fig_dir, 
      paste0("r2.13_panel_A.", ext)
    ),
    height = 9,
    width = 16,
    trimmed = FALSE,
    bg = "white"
  )
}

#### Panel B (Correlation log2FC loss versus delta AUC) -----------------------
frac_loss <- pred_loss[
  ,
  .(
    study,
    phenotype,
    tech,
    a_2,
    frac_mean
  )
][
  auc_stats[
    ,
    .(
      study,
      phenotype,
      tech,
      a_2,
      delta_mean
    )
  ],
  on = .(study, phenotype, tech, a_2)
]

# Correlation
frac_loss_cor <- frac_loss[
  ,
  `:=`(
    ancestry_n = .N,
    ancestry_cor = cor_safe(log2(frac_mean), delta_mean, .N)
  ),
  by = .(tech, a_2)
][
  ,
  `:=`(
    across_n = .N,
    across_cor = cor_safe(log2(frac_mean), delta_mean, .N)
  ),
  by = tech
][
  ,
  rbindlist(list(
    unique(.SD[, .(
      tech,
      cor_type  = "Across",
      cor_value = across_cor,
      n         = across_n
    )]),
    unique(.SD[, .(
      tech,
      cor_type  = a_2,
      cor_value = ancestry_cor,
      n         = ancestry_n
    )])
  ))
][
  ,
  cor_type := factor(
    cor_type,
    levels = c("Across", setdiff(unique(cor_type), "Across"))
  )
][
  ,
  cor_group := fifelse(
    cor_type == "Across",
    "Across",
    "Ancestry"
  )
][
  ,
  cor_group := factor(
    cor_group,
    levels = c("Across", "Ancestry")
  )
]

# Plot
r2.13_panel_B <- wrap_plots(
  
  # Scatter plot
  ggplot(
    data = frac_loss,
    mapping = aes(
      x = log2(frac_mean),
      y = delta_mean,
      color = a_2
    )
  ) +
  geom_polygon(
    data = {
      x <- log2(frac_loss$frac_mean)
      y <- frac_loss$delta_mean
      
      x_lims <- range(x, na.rm = TRUE)
      y_lims <- range(y, na.rm = TRUE)
      
      data.frame(
        x = c(x_lims[1], x_lims[1], x_lims[2]),
        y = c(y_lims[2], y_lims[1], y_lims[2])
      )
    },
    mapping = aes(
      x = x,
      y = y
    ),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.2
  ) +
  geom_point(
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_smooth(
    mapping = aes(
      x = log2(frac_mean),
      y = delta_mean
    ),
    linewidth = 0.7,
    formula = y ~ x,
    method = "lm",
    se = FALSE,
    color = "blue",
    inherit.aes = FALSE
  ) +
  facet_grid(
    rows = vars(tech),
    labeller = ggplot2::labeller(
      tech = c(
        meth = "Methylation",
        mrna = "Expression"
      )
    )
  ) +
  scale_color_manual(
    values = ancestry_cols
  ) +
  labs(
    color = "Ancestry",
    x = "Log2FC of log loss (non-EUR/ EUR)",
    y = "Delta ROC AUC (EUR - non-EUR)"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1,
    show_grid = FALSE,
    show_borders = TRUE
  ) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ),

  # Correlation plot
  ggplot(
    data = frac_loss_cor,
    mapping = aes(
      x = cor_type,
      y = cor_value,
      fill = cor_type
    )
  ) +
  geom_col(
    na.rm = TRUE
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  facet_grid(
    cols = vars(tech),
    labeller = ggplot2::labeller(
      tech = c(
        meth = "Methylation", 
        mrna = "Expression"
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      ancestry_cols,
      Across = "blue"
    )
  ) +
  labs(
    fill = "Ancestry",
    x = "Ancestry",
    y = "Correlation (Pearson)"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1,
    rotate = 45
  ) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ),
  ncol = 2
) +
plot_annotation(
  tag_levels = "a"
) &
theme(
  plot.tag = element_text(face = "bold", size = 10)
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r2.13_panel_B,
    file = file.path(
      fig_dir, 
      paste0(
        "r2.13_panel_B.", 
        ext
      )
    ),
    height = 7,
    width = 16,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}
