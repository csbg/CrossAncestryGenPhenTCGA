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

## Colors =====================================================================

ancestry_cols <- c(
  "EUR" = "#0072B2",
  "AFR" = "#D55E00",
  "EAS" = "#56B4E9",
  "AMR" = "#E69F00",
  "SAS" = "#009E73",
  "ADMIX" = "#999999"
)

## Directories ================================================================

# Result directories (per cance study)
res_dirs <- list(
  BRCA = file.path("results", "tcga", "analysis", "TCGA_BRCA"),
  UCEC = file.path("results", "tcga", "analysis", "TCGA_UCEC"),
  THCA = file.path("results", "tcga", "analysis", "TCGA_THCA")
)

# Figures
fig_dir <- file.path("results", "tcga", "figures", "subset_prediction_effect")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Tables
tab_dir <- file.path("results", "tcga", "tables", "subset_prediction_effect")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

## Results ====================================================================

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
      B <- nrow(w)
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
  p <- unique(out$phenotype)
  out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
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
    p <- unique(out$phenotype)
    out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
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
  }, files, techs)
  
  # Add phenotype and study
  out <- rbindlist(dts)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(out$phenotype)
  out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
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

## Main 4 =====================================================================

alpha <- 0.1
formats <- c("svg", "png")

### Panel B -------------------------------------------------------------------

# Plot
main4_panel_B <- ggplot(
  data = pred_loss[, `:=`(
        T_obs = log2(frac_mean),
        p_adj = p.adjust(p_value, method = "BH")
      )
  ],
  mapping = aes(
    x = a_2,
    y = phenotype,
    fill = T_obs,
    size = -log10(p_adj)
  )
) +
geom_point(
  shape = 21,
  color = "black",
  stroke = 0.1
) +
scale_fill_gradient2(
  name = "Log2FC of loss\n(non-EUR / EUR)  ",
  low  = "#4575b4",
  mid  = "white",
  high = "#d73027",
  midpoint = 0
) +
scale_size_continuous(
  name = expression(P.adj ~ "(-log"[10] * ")  "),
  range = c(2, 3.5),
  breaks = scales::pretty_breaks(n = 3)
) +
facet_grid(
  cols = vars(tech),
  rows = vars(study),
  scales = "free_y",
  space  = "free_y",
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
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
ggsaveDK(
  plot = main4_panel_B,
  file = file.path(
    fig_dir, 
    "main4_panel_B.svg"
  ),
  height = 5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel C -------------------------------------------------------------------

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
    set = fifelse(ancestry == "subset-EUR", "Validation", "Inference")
  )
][
  ,
  set := factor(set, levels = c("Validation", "Inference"))
]

# Plot
main4_panel_C <- ggplot(
  mapping = aes(
    x = a_2
  )
) +
geom_col(
  data = pred_AUC[
    study == "BRCA" &
    type == "Observed"
  ],
  mapping = aes(
    y = AUC_mean,
    fill = set
  ),
  position = position_dodge(
    width = 0.9
  ),
  width = 0.9
) +
geom_errorbar(
  data = pred_AUC[
    study == "BRCA" &
    type == "Permuted"
  ],
  mapping = aes(
    ymin = CI_lower, 
    ymax = CI_upper,
    group = set
  ),
  position = position_dodge(
    width = 0.9
  ),
  width = 0.2,
  linewidth = 0.3
) +
geom_point(
  data = pred_AUC[
    study == "BRCA" &
    type == "Permuted"
  ],
  mapping = aes(
    y = AUC_mean,
    group = set,
    shape = "Permuted AUC",
    color = "Permuted AUC"
  ),
  position = position_dodge(
    width = 0.9
  ),
  size  = 1
) +
scale_shape_manual(
  values = c("Permuted AUC" = 18)
) +
scale_color_manual(
  values = c("Permuted AUC" = "black")
) +
scale_fill_manual(
  values = c(
    "Inference" = "orange",
    "Validation" = "#0072B2"
  )
) +
scale_y_continuous(
  breaks = scales::pretty_breaks(n = 3)
) +
facet_grid(
  cols = vars(phenotype),
  rows = vars(tech),
  scales = "free",
  space = "free",
  labeller = labeller(
    phenotype = vs_newline,
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  x = "Ancestry",
  y = "Mean ROC AUC\nacross subsets",
  color = NULL,
  shape = NULL,
  fill  = NULL
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save
ggsaveDK(
  plot = main4_panel_C,
  file = file.path(
    fig_dir, 
    "main4_panel_C.svg"
  ),
  height = 6,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel D -------------------------------------------------------------------

# Compute overlap
overlap <- pred_feat[
  , .(mean_estimate = mean(estimate, na.rm = TRUE)),
  by = .(
    coef_id, coef_type, contrast, g_1, g_2, 
    a_1, a_2, study, tech, phenotype, feature
  )
][
  subset_res[coef_id == "interaction"],
  on = .(study, tech, phenotype, a_2, feature),
  nomatch = 0
][
  , .(
    n_sig     = sum(p_adj < alpha),
    n_nonzero = sum(mean_estimate != 0),
    n_both    = sum((p_adj < alpha) & (mean_estimate != 0)),
    recall    = sum((p_adj < alpha) & (mean_estimate != 0)) / 
                sum(mean_estimate != 0)
  ),
  by = .(study, tech, phenotype, a_2)
]

# Plot
main4_panel_D <- ggplot(
  mapping = aes(
    x = a_2,
    y = phenotype
  )
) +
geom_tile(
  data = overlap[
    recall == 0
  ],
  mapping = aes(
    fill = "No overlap"
  ),
  color = "black",
  linewidth = 0.1,
  show.legend = TRUE
) +
scale_fill_manual(
  name = NULL,
  values = c("No overlap" = "grey90"),
  guide = guide_legend(
    override.aes = list(color = "black")
  )
) +
ggnewscale::new_scale_fill() +
geom_tile(
  data = overlap[
    recall != 0
  ],
  mapping = aes(
    fill = recall
  ),
  color = "black",
  linewidth = 0.1
) +
scale_fill_gradient(
  name  = "Recall",
  low   = "white",
  high  = "gold",
  guide = guide_colorbar(order = 1)
) +
facet_grid(
  rows = vars(study),
  cols = vars(tech),
  scales = "free",
  space  = "free",
     labeller = labeller(
    tech = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  x = "Ancestry",
  y = "Cancer relationship"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.spacing.y = unit(10, "pt"),
  legend.margin = margin(0, 0, 0, 0)
) 

# Save
ggsaveDK(
  plot = main4_panel_D,
  file = file.path(
    fig_dir, 
    "main4_panel_D.svg"
  ),
  height = 5.5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel E -------------------------------------------------------------------

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
        mean_estimate = mean(estimate, na.rm = TRUE),
        abs_estimate  = abs(mean(estimate, na.rm = TRUE))
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

# Plot
main4_panel_E <- wrap_plots(
  ggplot(
    data = pred_feat[
      study == "BRCA" & 
      tech == "mrna" &
      phenotype == "Basal vs non-Basal" &
      feature %in% feat,
      .(mean_estimate = mean(estimate, na.rm = TRUE)),
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
    name = "Importance\nscore",
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
        .(mean_estimate = mean(estimate, na.rm = TRUE)),
        by = .(
          coef_id, coef_type, contrast, g_1, g_2, 
          a_1, a_2, study, tech, phenotype, feature
        )
      ]
      M <- max(abs(dt$mean_estimate), na.rm = TRUE)
      c(-M, M)
    }
    # labels = function(x) {
    #   zero_idx <- x == 0
    #   exp <- log10(x)
    #   out <- scales::math_format(10^.x)(exp)
    #   out[zero_idx] <- "0"
    #   out
    # }
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
    name = "log2FC  ",
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = expression(P.adj ~ "(-log"[10] * ")  "),
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
) + plot_layout(
  guides = "collect",
  widths = c(0.2, 1)
)

# Save
ggsaveDK(
  plot = main4_panel_E,
  file = file.path(
    fig_dir, 
    "main4_panel_E.svg"
  ),
  height = 12,
  width = 10,
  trimmed = FALSE,
  bg = "transparent"
)

### Supp. 1 -------------------------------------------------------------------

nr_degs <- subset_res[
  coef_id == "interaction",
  .(nr_degs = uniqueN(feature[p_adj <= alpha])),
  by = .(study, phenotype, tech, a_2)
][
  , degs_val := log10(nr_degs + 1)
][
  pred_loss[
    , .(study, phenotype, tech, a_2, loss_val = T_obs)
  ],
  on = .(study, phenotype, tech, a_2)
][
  , alpha := alpha
][]

# Correlation
nr_degs_cor <- nr_degs[
  , `:=`(
    ancestry_n = .N,
    ancestry_cor = cor_safe(degs_val, loss_val, .N)
  ),
  by = .(tech, a_2)
][
  , `:=`(
    across_n = .N,
    across_cor = cor_safe(degs_val, loss_val, .N)
  ),
  by = tech
][
  , rbindlist(list(
      unique(.SD[, .(
        tech,
        alpha,
        cor_type  = "Across",
        cor_value = across_cor,
        n         = across_n
      )]),
      unique(.SD[, .(
        tech,
        alpha,
        cor_type  = a_2,
        cor_value = ancestry_cor,
        n         = ancestry_n
      )])
    ))
][
  , cor_type := factor(
      cor_type,
      levels = c("Across", setdiff(unique(cor_type), "Across"))
    )
][
  , cor_group := fifelse(cor_type == "Across", "Across", "Ancestry")
][
  , cor_group := factor(
    cor_group, 
    levels = c("Across", "Ancestry")
  )
][
  , alpha := alpha
][]

# Plot
main4_supp1 <- wrap_plots(
  ggplot(
    data = nr_degs,
    mapping = aes(
      x = degs_val,
      y = loss_val,
      color = a_2
    )
  ) +
  geom_polygon(
    data = {

      x <- nr_degs$degs_val
      y <- nr_degs$loss_val
      
      lims <- range(c(x, y))
      
      data.frame(
        x = c(lims[1], lims[1], lims[2]),
        y = c(lims[2], lims[1], lims[2])
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
      x = degs_val,
      y = loss_val
    ),
    linewidth = 0.7,
    formula = y ~ x,
    method = "lm",
    se = FALSE,
    color = "blue",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::math_format(10^.x),
    expand = expansion(mult = 0.1)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::math_format(10^.x),
    expand = expansion(mult = 0.1)
  ) +
  facet_grid(
    rows = vars(tech),
    labeller = labeller(
      tech = c(meth = "Methylation", mrna = "Expression"),
    )
  ) +
  scale_color_manual(
    values = ancestry_cols
  ) +
  labs(
    color = "Ancestry",
    x = "# interaction DEGs",
    y = "log2FC of loss",
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
  ggplot(
    data = nr_degs_cor,
    mapping = aes(
      x = cor_type,
      y = cor_value,
      fill = cor_type
    )
  ) +
  geom_col(
    na.rm = TRUE
  ) +
  facet_grid(
    cols = vars(tech),
      labeller = labeller(
        tech = c(meth = "Methylation", mrna = "Expression"),
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
) + plot_annotation(
  tag_levels = "a"
) &
theme(
  plot.tag = element_text(
    face = "bold",
    size = 10
  ),
  plot.tag.position = c(0, 1)
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = main4_supp1,
    file = file.path(
      fig_dir, 
      paste0(
        "main4_supp1.", 
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
