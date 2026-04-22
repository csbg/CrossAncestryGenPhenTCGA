## max. code width ============================================================
## Libraries ==================================================================

library(CrossAncestryGenPhen)
library(ComplexHeatmap)
library(ggVennDiagram)
library(ggnewscale)
library(data.table)
library(yardstick)
library(patchwork)
library(ggplot2)
library(scales)

ancestry_colors <- c(
  "EUR" = "#0072B2",
  "AFR" = "#D55E00",
  "EAS" = "#56B4E9",
  "AMR" = "#E69F00",
  "SAS" = "#009E73",
  "ADMIX" = "#999999"
)

## Directories ================================================================

# Figures
fig_dir <- file.path("figures", "subset_prediction_effect")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Patchwork
patch_dir <- file.path(fig_dir, "patch_fig")
if (!dir.exists(patch_dir)) dir.create(patch_dir, recursive = TRUE)

# Excel sheets
tab_dir <- file.path("tables", "subset_prediction_effect")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

# Result directories (per cance study)
res_dirs <- list(
  BRCA = "results/TCGA_BRCA",
  UCEC = "results/TCGA_UCEC",
  THCA = "results/TCGA_THCA"
)

# Message
cat("Created directories to save data -------------------------------------\n")

## Load results ===============================================================

# Loss
logistic_prediction <- rbindlist(lapply(names(res_dirs), function(study) {

    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    ## List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$subsets_stats)
      dt[, tech := t]
    }, files, techs)

    ## Compute log losses → fractional summaries → bind
    ll <- lapply(dts, function(dt)
      dt[, .(logloss = mn_log_loss_vec(
        truth = droplevels(true), estimate = prob, event_level = "second")),
        by = .(iteration, coef_id, g_1, g_2, a_1, a_2, tech)]
    )

    ## Summaries
    fr <- lapply(ll, function(dt) {
      meta <- unique(dt[, .(g_1, g_2, a_1, a_2, tech)])
      w <- dcast(dt, iteration ~ coef_id, value.var = "logloss")
      w[, frac := relationship_Y / relationship_X]

      ## Add p-values
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

  ## Add phenotype and study
  out <- rbindlist(fr)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(out$phenotype)
  out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

# Probabilities
logistic_probabilities <- rbindlist(lapply(names(res_dirs), function(study) {
    ## Files 
    fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    ## List of subset stats + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f)$subsets_stats)
      dt[, tech := t]
    }, files, techs)

    ## Output
    out <- rbindlist(dts)
    out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
    p <- unique(out$phenotype)
    out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
    out[, study := study] 
    out
}), use.names = TRUE, fill = TRUE)

# Load feature importance
logistic_features <- rbindlist(lapply(names(res_dirs), function(study) {

  fp    <- file.path(res_dirs[[study]], "subset_prediction_effect")
  files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  techs <- sub("_.*$", "", basename(dirname(files)))

  ## List of subset stats + tech
  dts <- Map(function(f,t) {
    dt <- as.data.table(readRDS(f)$feature_stats)
    dt[, tech := t]
  }, files, techs)
  
  ## Add phenotype and study
  out <- rbindlist(dts)
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(out$phenotype)
  out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  out[, study := study] 
  out
}), use.names = TRUE, fill = TRUE)

# Load interactions
subset_interaction <- rbindlist(lapply(names(res_dirs), function(study) {
  ## Files
  fp <- file.path(res_dirs[[study]], "subset_interaction_effect")
  files <- list.files(fp, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  ## Load all files
  all_files <- lapply(files, readRDS)
  out <- rbindlist(lapply(all_files, \(x) x$summary_stats), use.names = TRUE, fill = TRUE)

  ## Add phenotyep study
  out[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(out$phenotype)
  out[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  out
}), fill = TRUE)

## Loss difference ============================================================
### LogFC ---------------------------------------------------------------------

# Define difference (log2FC of fraction)
sig_logistic_prediction <- copy(logistic_prediction)[, `:=`(
      T_obs = log2(frac_mean),
      p_adj = p.adjust(p_value, method = "BH")
    )
]

# Heatmap
gg_heatmap_loss_difference <- ggplot(
  data = sig_logistic_prediction,
  mapping = aes(
    x = a_2,
    y = phenotype,
    fill = T_obs
  )
) +
geom_tile(
  color = "black",
  linewidth = 0.1
) +
geom_point(
  mapping = aes(
    size = -log10(p_adj)
  ),
  shape = 8,
  stroke = 0.1,
  color  = "black"
) +
scale_fill_gradient2(
  name = "Log2FC of loss (non-EUR / EUR)  ",
  low  = "#4575b4",
  mid  = "white",
  high = "#d73027",
  midpoint = 0
) +
scale_size_continuous(
  name = expression(P.adj ~ "(-log"[10] * ")  "),
  range = c(0.5, 1),
  breaks = scales::pretty_breaks(n = 3)
) +
facet_grid(
  cols  = vars(tech),
  rows  = vars(study),
  scales = "free_y",
  space  = "free_y",
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression"),
    study = label_value  
  )
) +
labs(
  x = "Ancestry",
  y = "Cancer relationship"
) +
theme_CrossAncestryGenPhen(
  legend_key   = 1,
  rotate       = 45,
  show_borders = TRUE,
  show_grid    = FALSE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.position = "bottom",
  legend.title.position = "left",
  legend.box = "vertical",
  legend.box.spacing = unit(0, "pt"),
  legend.spacing = unit(0, "pt"),
  legend.spacing.y = unit(0, "pt"),
  legend.spacing.x = unit(0, "pt"),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
) +
guides(
  fill  = guide_colorbar(order = 1),
  shape = guide_legend(order = 2)
)

# Save
ggsaveDK(
  plot = gg_heatmap_loss_difference,
  file = file.path(fig_dir, "gg_heatmap_loss_difference.svg"),
  height = 5,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

# Dotplot
gg_dotplot_loss_difference <- ggplot(
  data = sig_logistic_prediction,
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
    tech  = c(meth = "Methylation", mrna = "Expression"),
    study = label_value  
  )
) +
labs(
  x = "Ancestry",
  y = "Cancer relationship"
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
  legend.position = "right",
  legend.margin = margin(0, 0, 0, 0)
)

# Save
ggsaveDK(
  plot = gg_dotplot_loss_difference,
  file = file.path(fig_dir, "gg_dotplot_loss_difference.svg"),
  height = 5,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

### Corr (loss logFC vs nr. DEGs) ---------------------------------------------

# Define DEGs.
dt_degs_nr <- rbindlist(
  lapply(
    c(0.05, 0.1),
    function(a) {
      subset_interaction[
        ,
        .(
          alpha = a,
          nr_sig_genes = sum(p_adj < a)
        ),
        by = .(study, phenotype, tech, coef_id, a_2)
      ]
    }
  )
)[coef_id == "interaction" & alpha == 0.1 ]

# Correlation
dt_cor <- merge(
  dt_degs_nr[
    ,
    .(study, phenotype, tech, a_2, cor_value = log10(nr_sig_genes + 1))
  ],
  sig_logistic_prediction[
    ,
    .(study, phenotype, tech, a_2, cor_value = T_obs)
  ],
  by = c("study", "phenotype", "tech", "a_2")
)

# Global
dt_cor[, across_cor :=
  if (.N >= 3 && sd(cor_value.x) > 0 && sd(cor_value.y) > 0)
    cor(cor_value.x, cor_value.y, method = "spearman")
  else NA_real_,
  by = .(tech)
]

across_cor <- unique(
  dt_cor[, .(tech, across_cor)]
)[
  , .(
    tech,
    cor_type = "Across",
    cor_value = across_cor
  )
]

# Within ancestry
dt_cor[, ancestry_cor :=
  if (.N >= 3 && sd(cor_value.x) > 0 && sd(cor_value.y) > 0)
    cor(cor_value.x, cor_value.y, method = "spearman")
  else NA_real_,
  by = .(tech, a_2)
]

ancestry_cor <- unique(
  dt_cor[, .(tech, a_2, ancestry_cor)]
)[
  , .(
    tech,
    cor_type = a_2,
    cor_value = ancestry_cor
  )
]

# Combine
cor_long <- rbind(across_cor, ancestry_cor)
cor_long[, cor_type := factor(
  cor_type,
  levels = c("Across", setdiff(sort(unique(cor_type)), "Across"))
)]
cor_long[, cor_group := ifelse(cor_type == "Across", "Across", "Ancestry")]
cor_long[, cor_group := factor(cor_group, levels = c("Across", "Ancestry"))]

# ggplot pointplot (relationship)
gg_pointplot_cor_degs_pred_alpha0.1 <- ggplot(
  data = dt_cor,
  mapping = aes(
    x = cor_value.x,
    y = cor_value.y,
    color = a_2
  )
) +
geom_smooth(
  mapping = aes(
    x = cor_value.x,
    y = cor_value.y
  ),
  method = "lm",
  se = FALSE,
  color = "blue",
  linewidth = 0.7,
  inherit.aes = FALSE
) +
geom_point(
  size = 0.5
) +
scale_y_continuous(
  breaks = scales::pretty_breaks(n = 3),
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
  values = ancestry_colors
) +
labs(
  color = "Ancestry",
  x = "# interaction DEGs (log10p)",
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
)

# ggplot barplot (spearman)
gg_barplot_cor_degs_pred_alpha0.1 <- ggplot(
  data = cor_long,
  mapping = aes(
    x = cor_value,
    y = cor_type,
    fill = cor_type
  )
) +
geom_col() +
facet_grid(
  cols = vars(tech),
    labeller = labeller(
      tech = c(meth = "Methylation", mrna = "Expression"),
  )
) +
scale_fill_manual(
  values = c(
    ancestry_colors,
    Across = "blue"
  )
) +
coord_flip() +
scale_x_continuous(
  limits = c(-1, 1),
  breaks = c(-1, 0, 1)
) +
labs(
  fill = NULL,
  x = "Correlation (Spearman)",
  y = "Ancestry"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Patchwork
p_cor_degs_pred_alpha0.1 <- (
  gg_barplot_cor_degs_pred_alpha0.1 +
    labs(x = "Spearman") +
    theme(plot.margin = margin(b = 15, 0, 0, 0)) +
    theme(legend.margin = margin(0, 0, 0, 0)) +
    theme(legend.position = "right")
) /
(
  gg_pointplot_cor_degs_pred_alpha0.1 +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0, 0, 0))
) +
plot_layout(guides = "collect", heights = c(0.5, 1))

# Save
ggsaveDK(
  plot = p_cor_degs_pred_alpha0.1,
  file = file.path(patch_dir, "p_cor_degs_pred_alpha0.1.svg"),
  height = 8,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

## Classification (AUC) =======================================================

# Calc. AUC
logistic_metrics <- logistic_probabilities[
  ,
  {
    # Observed AUC
    AUC_obs <- roc_auc_vec(
      truth = droplevels(true),
      estimate = prob,
      event_level = "second"
    )

    # Permuted AUCs (shuffle labels, keep probabilities fixed)
    AUC_perm <- replicate(
      100,
      roc_auc_vec(
        truth = sample(droplevels(true)),
        estimate = prob,
        event_level = "second"
      )
    )

    # Combibne
    rbindlist(list(
      data.table(type = "Observed",  AUC = AUC_obs),
      data.table(type = "Permuted",  AUC = AUC_perm)
    ))
  },
  by = .(iteration, coef_id, coef_type, a_1, a_2, g_1, g_2, study, tech, phenotype)
][
  ,
  # Summarize across iterations
  .(
    AUC_mean = mean(AUC),
    CI_lower = quantile(AUC, 0.025, na.rm = TRUE),
    CI_upper = quantile(AUC, 0.975, na.rm = TRUE),
    n_iter   = uniqueN(iteration)
  ),
  by = .(type, coef_id, coef_type, a_1, a_2, g_1, g_2, study, tech, phenotype)
]

# Grouping
logistic_metrics[, 
  ancestry := fifelse(
    coef_id == "relationship_X", a_1,
    fifelse(coef_id == "relationship_Y", a_2, 
    NA_character_
  )
)]

logistic_metrics[ancestry == "EUR", ancestry := "subset-EUR"]
logistic_metrics[, ancestry := factor(ancestry)]
logistic_metrics[, ancestry := factor(ancestry, levels = c("subset-EUR", setdiff(levels(ancestry), "subset-EUR")))]

logistic_metrics[, set := fifelse(ancestry == "subset-EUR", "Validation", "Inference")]
logistic_metrics[, set := factor(set, levels = c("Validation", "Inference"))]

# Per study AUC
AUC_list <- list()
for (s in c("BRCA", "THCA", "UCEC")){

  ## Subset data
  study_data <- logistic_metrics[
    study == s, 
  ]

  # Labeller
  vs_newline <- function(x) {
    gsub("\\s*vs\\s*", " vs\n", x)
  }

  ## ggplot
  p <- ggplot(
    data = study_data,
    mapping = aes(
      x = a_2
    )
  ) +
  geom_col(
    data = study_data[
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
    data = study_data[
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
    data = study_data[
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
    legend.position = "right",
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.spacing.x = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )


  # Save indiv. plot
  ggsaveDK(
    plot = p,
    file = file.path(fig_dir, paste0("gg_barplot_auc_", s, ".svg")),
    height = 6,
    width = 8,
    trimmed = FALSE,
    bg = "transparent"
  )

  ## Append to list
  AUC_list[[s]] <- p
}

## Patchwork
p_barplot_auc <-
  AUC_list[[1]] + AUC_list[[2]] +
  AUC_list[[3]] + guide_area() +
  plot_layout(
    ncol   = 2,
    guides = "collect"
  ) &
  theme(
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Save
ggsaveDK(
  plot = p_barplot_auc,
  file = file.path(patch_dir, "p_barplot_auc.svg"),
  height = 18,
  width = 18,
  trimmed = TRUE,
  bg = "transparent"
)



## Feature importance =========================================================
### fun. ----------------------------------------------------------------------

sets <- c(
  "EUR cancer effect",
  "non-EUR cancer effect",
  "Interaction effect",
  "Importance score"
)

matrix <- CJ(
  `EUR cancer effect` = 0:1,
  `non-EUR cancer effect` = 0:1,
  `Interaction effect` = 0:1,
  `Importance score` = 0:1
)[
  ## Drop empty intersection
  !(`EUR cancer effect` == 0 &
    `non-EUR cancer effect` == 0 &
    `Interaction effect` == 0 &
    `Importance score` == 0)
][]

count_sets <- function(sets, membership_matrix) {

  set_names <- names(sets)

  ## Build membership (long)
  membership <- rbindlist(
    lapply(set_names, function(nm) {
      data.table(feature = sets[[nm]], set = nm)
    }),
    fill = TRUE
  )

  ## Build wide matrix
  wide <- dcast(
    membership,
    feature ~ set,
    fun.aggregate = length,
    fill = 0
  )

  ## Ensure ALL set columns exist BEFORE using .SDcols
  for (nm in set_names) {
    if (!nm %in% names(wide)) {
      wide[, (nm) := 0L]
    }
  }

  ## Now it is SAFE to reference .SDcols
  wide[, (set_names) := lapply(.SD, as.integer), .SDcols = set_names]

  ## Count intersections
  counts <- wide[
    membership_matrix,
    on = set_names,
    .N,
    by = .EACHI
  ][
    membership_matrix,
    on = set_names
  ][
    is.na(N), N := 0L
  ]

  counts[]
}

make_sets <- function(a2_value, ph_value, t_value, s_value) {

  ## Interaction effects
  int_sub <- subset_interaction[
    coef_id == "interaction" &
    a_2 == a2_value &
    phenotype == ph_value &
    tech == t_value &
    study == s_value &
    p_adj < 0.1
  ]

  ## EUR cancer effects
  EUR_sub <- subset_interaction[
    coef_id == "relationship_1" &
    a_2 == a2_value &
    phenotype == ph_value &
    tech == t_value &
    study == s_value &
    p_adj < 0.1
  ]

  ## non-EUR cancer effects
  nonEUR_sub <- subset_interaction[
    coef_id == "relationship_2" &
    a_2 == a2_value &
    phenotype == ph_value &
    tech == t_value &
    study == s_value &
    p_adj < 0.1
  ]

  ## ML feature importance
  imp_sub <- logistic_features_agg[
    coef_type == "relationship" &
    a_2 == a2_value &
    phenotype == ph_value &
    tech == t_value &
    study == s_value &
    abs_mean_estimate > 0
  ]

  ## Extract meta
  meta <- data.table(
    a_2  = a2_value,
    phenotype = ph_value,
    tech = t_value,
    study  = s_value
  )

  list(
    sets = list(
      `EUR cancer effect` = unique(EUR_sub$feature),
      `non-EUR cancer effect` = unique(nonEUR_sub$feature),
      `Interaction effect` = unique(int_sub$feature),
      `Importance score` = unique(imp_sub$feature)
    ),
    meta = meta
  )
}

### Feature importance overlap with interaction effects -----------------------

# Aggreate feature importance
logistic_features_agg <- logistic_features[
  ,
  .(
    mean_estimate = mean(estimate, na.rm = TRUE),
    abs_mean_estimate = abs(mean(estimate, na.rm = TRUE))
  ),
  by = .(coef_id, coef_type, contrast, g_1, g_2, 
    a_1, a_2, study, tech, phenotype, feature
  )
]

# Compute all counts
cancer_relationships <- unique(
  subset_interaction[
    ,
    .(a_2, phenotype, tech, study)
  ]
)

all_counts <- rbindlist(
  lapply(seq_len(nrow(cancer_relationships)), function(i) {

    res <- make_sets(
      a2_value = cancer_relationships$a_2[i],
      ph_value = cancer_relationships$phenotype[i],
      t_value = cancer_relationships$tech[i],
      s_value = cancer_relationships$study[i]
    )

    counts <- count_sets(res$sets, matrix)

    counts[, `:=`(
      a_2 = res$meta$a_2,
      phenotype = res$meta$phenotype,
      tech = res$meta$tech,
      study = res$meta$study,
      intersection = paste0("I", .I)
    )]

    counts
  })
)

# Feature importance overlap with interaction effects
int_imp_counts <- all_counts[
  `Interaction effect` == 1 &
  `Importance score`   == 1
][, 
  `:=`(
      fill_group = fifelse(N == 0, " ", "effect"),
      fill_value = log10(N + 1)
    )
]

gg_heatmap_feature_overlap <- ggplot(
  mapping = aes(
    x = a_2,
    y = phenotype
  )
) +
geom_tile(
  data = int_imp_counts[
    fill_group == " "
  ],
  mapping = aes(
    fill = fill_group
  ),
  color = "black",
  linewidth = 0.1,
  show.legend = TRUE
) +
scale_fill_manual(
  name = " ",
  labels = "No overlap",
  values = c(" " = "grey90"),
  breaks = " ",
  guide = guide_legend(
    order = 2,
    override.aes = list(color = "black")
  )
) +
ggnewscale::new_scale_fill() +
geom_tile(
  data = int_imp_counts[fill_group == "effect"],
  mapping = aes(
    fill = fill_value
  ),
  color = "black",
  linewidth = 0.1
) +
scale_fill_gradient(
  name  = "# overlaping\ngenes",
  low   = "white",
  high  = "gold",
  limits = c(0, NA),
  breaks = c(0, 1, 2),
  labels = math_format(10^.x),
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
  legend.position = "right",
  legend.title.position = "top",
  legend.spacing.y = unit(10, "pt"),
  legend.margin = margin(0, 0, 0, 0)
) 

# Save
ggsaveDK(
  plot = gg_heatmap_feature_overlap,
  file = file.path(fig_dir, "gg_heatmap_feature_overlap.svg"),
  height = 5.5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Individual overlaps -------------------------------------------------------

# Filter for phenotype
filter <- expression(study == "BRCA" & tech == "mrna" & phenotype == "Basal vs non-Basal")
feature_imp <- logistic_features_agg[eval(filter), ]
interaction <- subset_interaction[eval(filter), ]

sig_interaction <- interaction[coef_id == "interaction" & p_adj < 0.1, unique(feature)[1:15], by = a_2]$V1
sig_importance  <- feature_imp[order(-abs_mean_estimate), unique(feature)[1:20], by = a_2]$V1

## Selcet features
ordered_features <- unique(c(sig_interaction, sig_importance))
ordered_features <- ordered_features[!is.na(ordered_features)]

## ggplot (feature importance)
plot_data <- feature_imp[feature %in% ordered_features]
plot_data[, feature := factor(
  feature,
  levels = ordered_features
)]

M <- max(abs(plot_data$mean_estimate), na.rm = TRUE)

plot_data <- dcast(
  plot_data,
  feature ~ a_2,
  value.var = "mean_estimate",
  fill = 0
)

p_feature_imp <- ggplot(
  data = melt(
    plot_data, 
    id.vars = "feature", 
    variable.name = "ancestry", 
    value.name = "importance"
  ),
  mapping = aes(
    x = ancestry, 
    y = feature, 
    fill = importance
  )
) +
geom_tile() +
scale_fill_gradient2(
  name = "Importance\nscore",
  low  = "#4575b4",
  mid  = "white",
  high = "#d73027",
  midpoint = 0,
  limits   = c(-M, M),
  labels = function(x) {
    zero_idx <- x == 0
    exp <- log10(abs(x))
    out <- scales::math_format(10^.x)(exp)
    out[zero_idx] <- "0"
    out
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
  axis.text.y = element_text(size = 4),
  panel.spacing.y = unit(0.15, "lines"),
  plot.margin = margin(r = 0.15, 0, 0, 0, "lines")
) 

## ggplot for limma
dotplot_data <- copy(interaction)[, 
  feature := factor(feature, levels = ordered_features)
]

p_dotplot <- ggplot(
  data = dotplot_data[
    eval(filter) & 
    coef_id %in% c(
      "interaction", 
      "relationship_1" , 
      "relationship_2"
    ) &
    feature %in% ordered_features
  ],
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
  legend.position = "bottom",
  legend.title.position = "top",
  legend.direction = "vertical",
  panel.spacing.x = unit(0.15, "lines"),
  plot.margin = margin(0, 0, 0, 0),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

# Patchwork
p <- p_feature_imp + p_dotplot +
  plot_layout(
    guides = "collect",
    widths = c(0.2, 1)
  ) &
  theme(
    legend.margin = margin(0, 0, 0, 0)
  )

# Save
ggsaveDK(
  plot = p,
  file = file.path(patch_dir, "p_feature_importance_BRCA.svg"),
  height = 12,
  width = 10,
  trimmed = FALSE,
  bg = "transparent"
)


## --- Cummulative interactions ---
effects <- subset_interaction[coef_type == "interaction",][
  logistic_features,
  on = .(g_1, g_2, a_1, a_2, study, tech, phenotype, feature),
  nomatch = 0
][mean_estimate != 0]

## Weighted interaction effect by feature importance
effects[, weighted_interaction := mean_estimate * T_obs]
effects[, feature_rank := frank(
    -abs(mean_estimate),
    ties.method = "first"
  ),
  by = .(study, phenotype, tech, a_2)
]

## Cummulative iteractions show shift in linear prediction
setorder(effects, study, phenotype, tech, a_2, feature_rank)
effects[, cum_interaction := cumsum(weighted_interaction),
  by = .(study, phenotype, tech, a_2)
][, norm_cum_interaction := cum_interaction / sum(abs(weighted_interaction)),
  by = .(study, phenotype, tech, a_2)
]

## ggplot
p_cum <- ggplot(
  data = effects[study == "BRCA" & phenotype == "LumA vs LumB"],
  mapping = aes(
    x = feature_rank, 
    y = norm_cum_interaction
  )
) +
geom_hline(
  yintercept = 0, 
  linetype   = "dashed",
  linewidth  = 0.1
) +
geom_line(
  linewidth = 0.3
) +
facet_grid(
  cols = vars(phenotype, a_2),
  rows = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  ),
  scales = "free"
) +
labs(
  x = "Ranked features",
  y = "Cumulative sum of interaction effects\nweighted by elastic net coefficients"
) +
theme_CrossAncestryGenPhen(
  show_borders = TRUE
)

# Save
ggsaveDK(
  plot = p_cum,
  file = file.path(out_dir, "Cumulative_interactions.svg"),
  height = 5.5,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)



## --- Individual overlaps ---
## Filter for phenotype
filter <- expression(study == "BRCA" & tech == "mrna" & phenotype == "Basal vs non-Basal")
feature_imp <- logistic_features[eval(filter), ]
interaction <- subset_interaction[eval(filter), ]

sig_interaction <- interaction[coef_id == "interaction" & p_adj < 0.1, unique(feature)[1:5], by = a_2]$V1
sig_importance  <- feature_imp[order(-abs_mean_estimate), unique(feature)[1:25], by = a_2]$V1

## Selcet features
ordered_features <- unique(c(sig_interaction, sig_importance))
ordered_features <- ordered_features[!is.na(ordered_features)]

## ggplot (feature importance)
plot_data <- feature_imp[feature %in% ordered_features, ]
plot_data[, feature := factor(feature, levels = ordered_features)]
M <- max(abs(plot_data$mean_estimate), na.rm = TRUE)

plot_data <- dcast(
  plot_data,
  feature ~ a_2,
  value.var = "mean_estimate",
  fill = 0
)

p_feature_imp <- ggplot(
  data = melt(
    plot_data, 
    id.vars = "feature", 
    variable.name = "ancestry", 
    value.name = "importance"
  ),
  mapping = aes(
    x = ancestry, 
    y = feature, 
    fill = importance
  )
) +
geom_tile() +
scale_fill_gradient2(
  name = "Importance  ",
  low  = "#4575b4",
  mid  = "white",
  high = "#d73027",
  midpoint = 0,
  limits   = c(-M, M),
  labels = function(x) {
    zero_idx <- x == 0
    exp <- log10(abs(x))
    out <- scales::math_format(10^.x)(exp)
    out[zero_idx] <- "0"
    out
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
  legend.position = "bottom",
  legend.title.position = "top",
  legend.direction = "vertical",
  axis.text.y = element_text(size = 5),
  panel.spacing.y = unit(0.15, "lines"),
  plot.margin = margin(0, 0, 0, 0)
) 

## ggplot for limma
dotplot <- copy(interaction)[, feature := factor(feature, levels = ordered_features)]

p_dotplot <- ggplot(
  mapping = aes(
    x = a_2, 
    y = feature,
    size = -log10(p_adj),
    fill = T_obs
  )
) +
geom_point(
  data = dotplot[
    eval(filter) & 
    coef_id == "interaction" & 
    feature %in% ordered_features
  ],
  shape = 21,
  stroke = 0
) +
geom_point(
  data = dotplot[
    eval(filter) & 
    coef_id == "interaction" & 
    feature %in% ordered_features & 
    p_adj < 0.1
  ],
  mapping = aes(
    color = p_adj < 0.1
  ),
  shape = 21,
  stroke = 0.2
) +
geom_point(
  data = dotplot[
    eval(filter) & 
    coef_id %in% c("relationship_1", "relationship_2") & 
    feature %in% ordered_features
  ],
  shape  = 21,
  stroke = 0
) +
geom_point(
  data = dotplot[
    eval(filter) & 
    coef_id %in% c("relationship_1", "relationship_2") & 
    feature %in% ordered_features & 
    p_adj < 0.1
  ],
  shape  = 21,
  stroke = 0.2,
  color  = "black"
) +
scale_color_manual(
  name   = "Significant effect\n(alpha < 0.1)  ",
  values = c(`TRUE` = "black", `FALSE` = NA),
  labels = c(`TRUE` = ""),
  na.translate = FALSE
) +
scale_fill_gradient2(
  name = "log2FC  ",
  low      = "#4575b4",
  mid      = "white",
  high     = "#d73027",
  midpoint = 0
) +
scale_size_continuous(
  name   = expression(P.adj ~ "(-log"[10] * ")  "),
  range  = c(0.8, 2),
  breaks = c(0, 2, 5, 10, 12)
) +
guides(
  size = guide_legend(
    override.aes = list(
      fill = "black"
    )
  ),
  color = guide_legend(
    override.aes = list(
      size   = 2,
      stroke = 0.2
    )
  )
) +
facet_grid(
  cols   = vars(coef_id),
  scales = "free_y",
  space  = "free_y",
  labeller = labeller(
    coef_id = c(
      interaction    = "Interaction",
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
  legend.position = "bottom",
  legend.title.position = "top",
  legend.direction = "vertical",
  panel.spacing.x = unit(0.15, "lines"),
  plot.margin = margin(t = 0, b = 0, l = 0.15, r = 0, unit = "lines"),
  # Remove x.axis
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

# Patchwork
p <- p_feature_imp + p_dotplot + plot_layout(widths = c(0.2, 1)) & 
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
  )

# Save
ggsaveDK(
  plot = p,
  file = file.path(out_dir, "Feature_importance.svg"),
  height = 11,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

## --- Upset heatmaps ---


## --- Set matrix ---
sets <- c(
  "EUR cancer effect",
  "non-EUR cancer effect",
  "Interaction effect",
  "Importance score"
)

matrix <- CJ(
  `EUR cancer effect` = 0:1,
  `non-EUR cancer effect` = 0:1,
  `Interaction effect` = 0:1,
  `Importance score` = 0:1
)[
  ## Drop empty intersection
  !(`EUR cancer effect` == 0 &
    `non-EUR cancer effect` == 0 &
    `Interaction effect` == 0 &
    `Importance score` == 0)
][]

## --- Counts per set ---
count_sets <- function(sets, membership_matrix) {

  set_names <- names(sets)

  ## Build membership (long)
  membership <- rbindlist(
    lapply(set_names, function(nm) {
      data.table(feature = sets[[nm]], set = nm)
    }),
    fill = TRUE
  )

  ## Build wide matrix
  wide <- dcast(
    membership,
    feature ~ set,
    fun.aggregate = length,
    fill = 0
  )

  ## Ensure ALL set columns exist BEFORE using .SDcols
  for (nm in set_names) {
    if (!nm %in% names(wide)) {
      wide[, (nm) := 0L]
    }
  }

  ## Now it is SAFE to reference .SDcols
  wide[, (set_names) := lapply(.SD, as.integer), .SDcols = set_names]

  ## Count intersections
  counts <- wide[
    membership_matrix,
    on = set_names,
    .N,
    by = .EACHI
  ][
    membership_matrix,
    on = set_names
  ][
    is.na(N), N := 0L
  ]

  counts[]
}

## --- Compute all counts ----
conditions <- unique(
  subset_interaction[
    ,
    .(a_2, phenotype, tech, study)
  ]
)

all_counts <- rbindlist(
  lapply(seq_len(nrow(conditions)), function(i) {

    res <- make_sets(
      a2_value = conditions$a_2[i],
      ph_value = conditions$phenotype[i],
      t_value  = conditions$tech[i],
      s_value  = conditions$study[i]
    )

    counts <- count_sets(res$sets, matrix)

    counts[, `:=`(
      a_2         = res$meta$a_2,
      phenotype   = res$meta$phenotype,
      tech        = res$meta$tech,
      study       = res$meta$study,
      intersection = paste0("I", .I)
    )]

    counts
  })
)


## ---- UpSet membership matrix (structure only) ----
set_cols <- c(
  "EUR cancer effect",
  "non-EUR cancer effect",
  "Interaction effect",
  "Importance score"
)

# One copy of membership matrix
upset_def <- unique(
  all_counts[, c("intersection", set_cols), with = FALSE]
)

# Enforce  order
upset_def[, intersection := factor(
  intersection,
  levels = paste0("I", 1:15)
)]

# Long format for ggplot
upset_long <- melt(
  upset_def,
  id.vars = "intersection",
  measure.vars = set_cols,
  variable.name = "set",
  value.name = "member"
)[member == 1]

# Fixed vertical order
upset_long[, set := factor(set, levels = rev(set_cols))]

# Data for vertical connecting lines
lines_dt <- upset_long[
  ,
  .(
    y_min = min(as.numeric(set)),
    y_max = max(as.numeric(set))
  ),
  by = intersection
]

# Row background bands
row_bg <- data.table(
  set = factor(rev(set_cols), levels = rev(set_cols)),
  ymin = seq_along(set_cols) - 0.5,
  ymax = seq_along(set_cols) + 0.5
)
# Alternating row colors
row_bg[, fill := rep(c("#F7F7F7", "#FFFFFF"), length.out = .N)]

# Shadow for rows
shadow_dt <- CJ(
  intersection = levels(upset_long$intersection),
  set = factor(rev(set_cols), levels = rev(set_cols))
)

## --- UpSet membership plot ---
p_upset <- ggplot() +
geom_rect(
  data = row_bg,
  aes(
    xmin = -Inf,
    xmax = Inf,
    ymin = ymin,
    ymax = ymax,
    fill = fill
  ),
  inherit.aes = FALSE
) +
scale_fill_identity() +
geom_point(
  data = shadow_dt,
  aes(
    x = intersection,
    y = as.numeric(set)
  ),
  color = "grey90",
  size  = 0.5
) +
geom_segment(
  data = lines_dt,
  aes(
    x    = intersection,
    xend = intersection,
    y    = y_min,
    yend = y_max
  ),
  linewidth = 0.3,
  color = "black"
) +
geom_point(
  data = upset_long,
  aes(
    x = intersection,
    y = as.numeric(set)
  ),
  size  = 0.5,
  color = "black"
) +
scale_y_continuous(
  breaks = seq_along(set_cols),
  labels = rev(
    c(
    "EUR cancer effect\n(alpha < 0.1)",
    "non-EUR cancer effect\n(alpha < 0.1)",
    "Interaction effect\n(alpha < 0.1)",
    "Importance score\n(|β| > 0)"
    )
  ),
  expand = expansion(add = 0.4)
) +
labs(
  x = NULL,
  y = "Gene sets", 
) +
theme_CrossAncestryGenPhen(
  show_axis = FALSE
) +
theme(
  axis.ticks  = element_blank(),
  axis.text.x = element_blank(),
  plot.margin = margin(t = 0, b = 0, l = 0, r = 0)
)

## --- UpSet heatmap ---
p_heat <- ggplot(
  all_counts[phenotype == "Basal vs non-Basal"],
  aes(
    x = intersection,
    y = a_2,
    fill = log10(N + 1)
  )
) +
geom_tile(
  color = "white", 
  linewidth = 0.0
) +
facet_grid(
  rows = vars(tech),
  scales = "free",
  space = "free",
   labeller = labeller(
    tech = c(meth = "Methylation", mrna = "Expression")
  )
) +
scale_fill_gradient(
  low  = "white",
  high = "purple",
  limits = c(0, NA),
  labels = math_format(10^.x)
) +
labs(
  x = NULL,
  y = "Ancestry",
  fill = "Overlap\nsize"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_axis = FALSE,
  show_borders = TRUE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  axis.ticks   = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  plot.margin  = margin(t = 0, b = 0, l = 0, r = 0)
)

# Patchwork
p <- p_heat / p_upset + plot_layout(heights = c(1, 1)) & 
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
  )

# Save
ggsaveDK(
  plot = p,
  file = file.path(out_dir, "Feature_overlap.svg"),
  height = 5.5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

## --- Interaction effect vs Importance scores ---
int_imp_counts <- all_counts[
  `Interaction effect` == 1 &
  `Importance score`   == 1
][, 
  `:=`(
      fill_group = fifelse(N == 0, " ", "effect"),
      fill_value = log10(N + 1)
    )
]


p <- ggplot(
  mapping = aes(
    x    = a_2,
    y    = phenotype
  )
) +
geom_tile(
  data = int_imp_counts[fill_group == " "],
  mapping = aes(
    fill = fill_group
  ),
  color = "black",
  linewidth = 0.1,
  show.legend = TRUE
) +
scale_fill_manual(
  name   = "No overlaps  ",
  values = c(" " = "grey90"),
  guide  = guide_legend(
    override.aes = list(color = "black")
  )
) +
ggnewscale::new_scale_fill() +
geom_tile(
  data = int_imp_counts[fill_group == "effect"],
  mapping = aes(
    fill = fill_value
  ),
  color = "black",
  linewidth = 0.1
) +
scale_fill_gradient(
  name  = "Overlap size  ",
  low   = "white",
  high  = "purple",
  limits = c(0, NA),
  breaks = c(0, 1, 2),
  labels = math_format(10^.x)
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
  y = "Phenotpye"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE
) +
theme(
  panel.spacing.x       = unit(0.15, "lines"),
  panel.spacing.y       = unit(0.15, "lines"),
  legend.position       = "bottom",
  legend.title.position = "left",
  legend.box            = "vertical",
  legend.box.spacing    = unit(0, "pt"),
  legend.spacing        = unit(0, "pt"),
  legend.spacing.y      = unit(0, "pt"),
  legend.spacing.x      = unit(0, "pt"),
  legend.margin         = margin(t = 0, r = 0, b = 0, l = 0),
)

# Save
ggsaveDK(
  plot = p,
  file = file.path(out_dir, "Interaction_overlap.svg"),
  height = 5.5,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

## --- Correlation feature importance scores vs interaction effect ----
imp_dt <- feature_imp[
  coef_type == "relationship" &
  a_1 == "EUR" &
  abs_mean_estimate > 0,
  .(
    feature,
    importance = abs_mean_estimate,
    study,
    tech,
    phenotype
  )
]

int_dt <- interaction[
  coef_type == "interaction",
  .(
    feature,
    a_2,
    beta_int = T_obs,
    p_int    = p_adj,
    study,
    tech,
    phenotype
  )
]

dt <- merge(
  imp_dt,
  int_dt,
  by = c("feature", "study", "tech", "phenotype"),
  allow.cartesian = TRUE
)

dt[, weighted_int := importance * abs(beta_int)]
dt[, signed_weighted_int := importance * beta_int]

impact_by_ancestry <- dt[
  ,
  .(
    net_shift      = sum(signed_weighted_int),
    abs_shift      = sum(abs(signed_weighted_int)),
    n_features     = .N
  ),
  by = .(a_2, study, tech, phenotype)
]


ggplot(
  dt,
  aes(abs(beta_int))
) +
  geom_histogram(bins = 60, fill = "grey70", color = "white") +
  labs(
    x = "|Interaction effect (β_int)|",
    y = "Number of features",
    title = "Ancestry interaction effects are modest at the single-feature level"
  ) +
  theme_minimal()

  ggplot(
  dt,
  aes(
    x = importance,
    y = abs(beta_int)
  )
) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  labs(
    x = "ML feature importance (EUR-trained)",
    y = "|Interaction effect|",
    title = "Model relies on features with non-zero interaction effects"
  ) +
  theme_minimal()

cum_dt <- dt[
  order(-weighted_int)
][
  ,
  cum_weighted_int := cumsum(weighted_int)
][
  ,
  rank := .I
]

ggplot(
  cum_dt,
  aes(rank, cum_weighted_int)
) +
  geom_line(linewidth = 0.7) +
  labs(
    x = "Features ranked by importance-weighted interaction",
    y = "Cumulative interaction-induced predictor shift",
    title = "Many small interaction effects accumulate into large prediction shifts"
  ) +
  theme_minimal()


ggplot(
  dt,
  aes(abs(beta_int), weighted_int, color = p_int < 0.1)
) +
  geom_point(alpha = 0.4) +
  labs(
    x = "|Interaction effect|",
    y = "Importance-weighted interaction impact",
    color = "p_adj < 0.1",
    title = "Statistical significance underestimates predictive impact"
  ) +
  theme_minimal()

# The model fails in non-EUR because it relies on many features whose ancestry-specific interaction effects are 
# individually modest but collectively large once weighted by model importance, leading to a systematic shift in 
# the model’s linear predictor.