## max. code width ============================================================
## Libraries ==================================================================

library(CrossAncestryGenPhen)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(openxlsx)

## Directories ================================================================

# Results
in_dir  <- file.path("results", "permutation_sim")

# Figures
fig_dir <- file.path("figures", "permutation_sim")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Patchwork
patch_dir <- file.path(fig_dir, "patch_fig")
if (!dir.exists(patch_dir)) dir.create(patch_dir, recursive = TRUE)

# Excel sheets
tab_dir <- file.path("tables", "permutation_sim")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

## Load results ===============================================================

res_path <- file.path(in_dir, "all_methods_results.rds")
res      <- readRDS(res_path)
res_mem  <- object.size(res)
message(sprintf("Loaded file: %s uses %d of RAM.", res_path, res_mem))

res[
  kind == "real" &
  coef_id == "interaction" &
  p_adj < 0.05 &
  method %in% c("limma-voom", "subset-cct", "subset_bonfferoni"),
  .N,
  by = .(method, study)
]

res[
  kind == "perm" &
  coef_id == "interaction" &
  method %in% c("limma-voom", "subset-cct", "subset_bonfferoni"),
][
  , .(n_sig = sum(p_adj < 0.05, na.rm = TRUE)),
  by = .(study, coef_id, method, perm_id)
][
  , .(
      mean_n_sig = mean(n_sig),
      sd_n_sig   = sd(n_sig)
    ),
  by = .(study, coef_id, method)
]

## Figure 1A ==================================================================
### fun. ----------------------------------------------------------------------

plot_null_vs_real_degs <- function(
  data,
  sig_thr,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){

  # Split data 
  perm_res <- data[kind == "perm"]
  real_res <- data[kind == "real"]

  # Permutation: count significant genes
  perm_counts <- perm_res[, .(
    n_sig = sum(p_adj < sig_thr, na.rm = TRUE)
  ), by = .(study, coef_id, method, perm_id)]

  # Summarize across permutations
  perm_summary <- perm_counts[, .(
    mean_DEG = mean(log10(n_sig + 1), na.rm = TRUE),
    sd_DEG   = sd(log10(n_sig + 1), na.rm = TRUE)
  ), by = .(study, coef_id, method)]

  # Original: count significant genes 
  real_summary <- real_res[, .(
    real_DEG = mean(log10(sum(p_adj < sig_thr, na.rm = TRUE) + 1))
  ), by = .(study, coef_id, method)]

  # Merge
  plot_dt <- merge(
    perm_summary, real_summary,
    by = c("study", "method", "coef_id"),
    all = TRUE
  )

  # Construct stacked segments
  stacked_dt <- rbindlist(list(
    # lower bar = DEGs in permuted data
    plot_dt[, .(
      study, method, coef_id,
      y = mean_DEG,
      segment = "DEGs in permuted data (FP)"
    )],
    # upper bar = potential TP (real - mean)
    plot_dt[real_DEG > mean_DEG, .(
      study, method, coef_id,
      y = real_DEG - mean_DEG,
      segment = "Potential true positives (TP)"
    )]
  ), use.names = TRUE, fill = TRUE)

  # Ensure stacking order
  stacked_dt[, segment := factor(segment,
    levels = c("Potential true positives (TP)", "DEGs in permuted data (FP)")
  )]

  # Plot
  p <- ggplot(
    data = stacked_dt,
    aes(
      x = method,
      y = y,
      fill = segment
    ) 
  ) +
  geom_col(
    color = "black", 
    linewidth = 0.1, 
    position = "stack"
  ) +
  geom_errorbar(
    data = plot_dt,
    aes(
      x = method,
      ymin = mean_DEG - sd_DEG,
      ymax = mean_DEG + sd_DEG
    ),
    width = 0.2,
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = plot_dt,
    aes(
      x = method,
      y = real_DEG,
      shape = "DEGs in original data (TP + FP)",
      color = "DEGs in original data (TP + FP)"
    ),
    size = point_size,
    inherit.aes = FALSE
  ) +
  facet_grid(
    cols = vars(study)
  ) +
  labs(
    title = title,
    x = ifelse(is.null(x_label), "Method", x_label),
    y = ifelse(is.null(y_label), paste0("Log10 # of identified DEGs (alpha = ", sig_thr, ")" ), y_label),
    fill = "",
    shape = "",
    color = ""
  ) +
  scale_fill_manual(
    values = c(
      "DEGs in permuted data (FP)" = "grey80",
      "Potential true positives (TP)" = "#fc8d62"
    )
  ) +
  scale_shape_manual(
    values = c("DEGs in original data (TP + FP)" = 18)
  ) +
  scale_color_manual(
    values = c("DEGs in original data (TP + FP)" = "#1f78b4")
  ) +
  guides(
    color = guide_legend(
      order = 1, 
      override.aes = list(size = 2),
      ncol = 1
    ),
    shape = guide_legend(
      order = 1, 
      override.aes = list(size = 2),
      ncol = 1
    ),
    fill = guide_legend(
      order = 2, 
      override.aes = list(color = NA), 
      reverse = TRUE,
      ncol = 1
    )
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1, 
    rotate = 45, 
    show_facets = FALSE,
    show_borders = FALSE
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.x = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )

  # Return
  return(p)
}

### TCGA study imbalance ------------------------------------------------------

study_configs <- list(
  BRCA = list(
    file  = "data/TCGA_BRCA_GDC_study.rds",
    tech = "mrna",
    g_col = "pooled_subtype",
    g1    = "non_Basal",
    g2    = "Basal",
    a_col = "pooled_ancestry",
    a1    = "EUR",
    a2    = "AFR",
    title = "Breast Invasive Carcinoma\n(BRCA)",
    fill_label = ""
  ),
  UCEC = list(
    file  = "data/TCGA_UCEC_GDC_study.rds",
    tech = "mrna",
    g_col = "subtype",
    g1    = "Endometrioid",
    g2    = "Serous",
    a_col = "pooled_ancestry",
    a1    = "EUR",
    a2    = "AFR",
    title = "Uterine Endometrial Carcinoma\n(UCEC)",
    fill_label = ""
  )
)

plots <- list()
for (study_name in names(study_configs)) {
  cfg <- study_configs[[study_name]]
  study <- readRDS(cfg$file)[[cfg$tech]]

  study <- filter_phenotype_ancestry(
    X = study$matr,
    M = study$meta,
    g_col = cfg$g_col,
    a_col = cfg$a_col,
    g_levels = c(cfg$g1, cfg$g2),
    a_levels = c(cfg$a1, cfg$a2),
    plot = FALSE,
    verbose = TRUE
  )

  p <- plot_imbalanced_groups(
    MX = study$X$meta,
    MY = study$Y$meta,
    x_var = cfg$a_col,
    fill_var = cfg$g_col,
    title = cfg$title,
    x_label = "Ancestry",
    y_label = "Number patients",
    fill_label = cfg$fill_label
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_CrossAncestryGenPhen(legend_key = 1, rotate = 45) +
  theme(legend.position = "top", legend.title.position = "left")

  plots[[study_name]] <- p
}

### BRCA perm. DEGs -----------------------------------------------------------

# Alpha = 0.05
perm_degs0.05_BRCA_p <- plot_null_vs_real_degs(
  data = res[coef_id == "interaction" & study == "BRCA", ],
  sig_thr = 0.05,
  point_size = 2
)

# Save
ggsaveDK(
  plot = perm_degs0.05_BRCA_p,
  file = file.path(fig_dir, "perm_degs0.05_BRCA.svg"),
  width = 8,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha = 0.1
perm_degs0.1_BRCA_p <- plot_null_vs_real_degs(
  data = res[coef_type == "interaction" & study == "BRCA", ],
  sig_thr = 0.1,
  point_size = 2
)

# Save
ggsaveDK(
  plot = perm_degs0.05_BRCA_p,
  file = file.path(fig_dir, "perm_degs0.1_BRCA.svg"),
  width = 8,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### UCEC perm. DEGs -----------------------------------------------------------

# Alpha = 0.05
perm_degs0.05_UCEC_p <- plot_null_vs_real_degs(
  data = res[coef_type == "interaction" & study == "UCEC", ],
  sig_thr = 0.05,
  point_size = 2
)

# Save
ggsaveDK(
  plot = perm_degs0.05_UCEC_p,
  file = file.path(fig_dir, "perm_degs0.05_UCEC.svg.svg"),
  width = 8,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha = 0.1
perm_degs0.1_UCEC_p <- plot_null_vs_real_degs(
  data = res[coef_type == "interaction" & study == "UCEC", ],
  sig_thr = 0.1,
  point_size = 2
)

# Save
ggsaveDK(
  plot = perm_degs0.1_UCEC_p,
  file = file.path(fig_dir, "perm_degs0.1_UCEC.svg.svg"),
  width = 8,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Tables --------------------------------------------------------------------
# Split data 
perm_res <- res[coef_type == "interaction" & kind == "perm"]
real_res <- res[coef_type == "interaction" & kind == "real"]

# Permutation: count significant genes
perm_counts <- perm_res[, .(
  n_sig = sum(p_adj < 0.05, na.rm = TRUE)
), by = .(study, coef_id, method, perm_id)]

# Summarize across permutations
perm_summary <- perm_counts[, .(
  log10_perm_DEG = mean(log10(n_sig + 1), na.rm = TRUE),
  log10_sd_DEG = sd(log10(n_sig + 1), na.rm = TRUE),
  perm_DEG = mean(n_sig, na.rm = TRUE),
  sd_DEG = sd(n_sig, na.rm = TRUE)
), by = .(study, coef_id, method)]

# Original: count significant genes 
real_summary <- real_res[, .(
  log10_real_DEG = mean(log10(sum(p_adj < 0.05, na.rm = TRUE) + 1)),
  real_DEG = mean(sum(p_adj < 0.05, na.rm = TRUE))
), by = .(study, coef_id, method)]

# Merge
plot_dt <- merge(
  perm_summary, real_summary,
  by = c("study", "method", "coef_id"),
  all = TRUE
)

### Patchwork -----------------------------------------------------------------

# Alpha = 0.05
p_perm_degs_alpha0.05 <- (
    plots$BRCA + plots$UCEC
) / 
(
  (
    perm_degs0.05_BRCA_p + 
     labs(y = "Log10 # of identified DEGs") +
    perm_degs0.05_UCEC_p + 
      labs(y = "Log10 # of identified DEGs")
    ) + 
    plot_layout(guides = "collect") +
    theme(legend.position = "bottom") &
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.x = unit(0, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 1, 0, 0)
    )
) + 
plot_layout(heights = c(0.5, 1))

# Save
ggsaveDK(
  plot = p_perm_degs_alpha0.05,
  file = file.path(patch_dir, "p_perm_degs_alpha0.05.svg"),
  height = 11,
  width = 8,
  bg = "transparent"
)

# Alpha = 0.1
p_perm_degs_alpha0.1 <- (
  plots$BRCA + plots$UCEC
) / 
(
  (
    perm_degs0.1_BRCA_p + 
      labs(y = "Log10 # of identified DEGs") +
    perm_degs0.1_UCEC_p + 
      labs(y = "Log10 # of identified DEGs")
    ) + 
    plot_layout(guides = "collect") +
    theme(legend.position = "bottom") &
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.spacing.x = unit(0, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 1, 0, 0)
    )
) + 
plot_layout(heights = c(0.5, 1))

# Save
ggsaveDK(
  plot = p_perm_degs_alpha0.1,
  file = file.path(patch_dir, "p_perm_degs_alpha0.1.svg"),
  height = 11,
  width = 8,
  bg = "transparent"
)

## Supplementary figure 1A ==================================================== 
### fun. ----------------------------------------------------------------------

plot_pval_histograms <- function(
  data,
  density = FALSE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  bins = 50
) {

  # Filter for perm
  data <- data[data$kind == "perm", ]
  

  # Plot
  if (density) {
    # Density curves 
    if (!"perm_id" %in% names(data)) {
      stop("Need column 'perm_id' for density plots.")
    }

    # Calculate bin width
    binwidth <- 1 / bins

    p <- ggplot(
      data = data, 
      mapping = aes(
        x = p_value
        )
      ) +
      geom_histogram(
        mapping = aes(
            y = log10(after_stat(count) + 1)
        ),
        binwidth = binwidth,
        fill = "gray80",
        color = "black",
        linewidth = 0.1
      ) +
      geom_density(
        mapping = aes(
          y = log10(after_stat(density) * nrow(data) * binwidth + 1),
          group = perm_id, 
          color = factor(perm_id)
        ),
        linewidth = 0.3
      ) +
      facet_grid(
        study ~ method,
        scales = "free_y"
      ) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(
        title = title,
        x = ifelse(is.null(x_label), "P-value", x_label),
        y = ifelse(is.null(y_label), "Log10 Count" , y_label),
        color = "Permutation"
      ) +
      theme_CrossAncestryGenPhen(
        legend.position = "none"
      )

  } else {

    # Pooled histograms
    p <- ggplot(
      data = data, 
      mapping = aes(
          x = p_value
        )
      ) +
      geom_histogram(
        bins = bins,
        fill = "gray80",
        color = "black",
        linewidth = 0.1
      ) +
      facet_grid(
        study ~ method,
        scales = "free_y"
      ) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(
        title = title,
        x = ifelse(is.null(x_label), "P-value", x_label),
        y = ifelse(is.null(y_label), "Count" , y_label),
      ) +
      theme_CrossAncestryGenPhen(
        legend.position = "none"
      )
  }

  # Return
  return(p)
}

### Pvalue dist. --------------------------------------------------------------

null_dist_p <- plot_pval_histograms(
  data = res,
  density = TRUE,
  bins = 30
)

# Save
ggsaveDK(
  plot = null_dist_p,
  file = file.path(fig_dir, "null_dist.svg"),
  width = 16,
  height = 4
)

### Patchwork -----------------------------------------------------------------

sub_fig_1A <- null_dist_p

# Save
ggsaveDK(
  plot = sub_fig_1A,
  file = file.path(fig_dir, "Supplementary_fig_1A.svg"),
  width = 16,
  height = 4
)
