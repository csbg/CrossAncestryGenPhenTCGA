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
in_dir  <- file.path("results", "synthetic_sim")

# Figures
fig_dir <- file.path("figures", "synthetic_sim")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Patchwork
patch_dir <- file.path(fig_dir, "patch_fig")
if (!dir.exists(patch_dir)) dir.create(patch_dir, recursive = TRUE)

# Excel sheets
tab_dir <- file.path("tables", "synthetic_sim")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

## Load results ===============================================================

# Raw results
res_raw_path <- file.path(in_dir, "all_methods_results.rds")
res_raw      <- readRDS(res_raw_path)
res_raw_mem  <- object.size(res_raw)
message(sprintf("Loaded file: %s uses %d RAM.", res_raw_path, res_raw_mem))

# Create the imbalance ratio lables
res_raw[, ratio_label := paste0("br = ", between_ratio, ", wr = ", within_ratio)]
res_raw[, ratio_label := factor(
  ratio_label, 
  levels = unique(ratio_label[order(between_ratio)])
)]

# Filter interaction 
res_degs  <- res_raw[n_degs > 0 & coef_type %in% c("interaction")]
res_0degs <- res_raw[n_degs == 0 & coef_type %in% c("interaction")]

# Compute metric
alpha_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)

g_vars <- c(
  "study", "coef_id", "contrast", "coef_type", "g_1", "g_2", 
  "a_1", "a_2", "method", "n_samples", "n_degs", "log2fc", 
  "total_samples", "between_ratio", "within_ratio", "ratio_label"
)

confusion_matrix <- res_degs[
  ,
  {
    res <- lapply(alpha_grid, function(alpha) {
      rej   <- p_adj <= alpha
      truth <- T_truth != 0
      
      TP <- sum(rej & truth)
      FP <- sum(rej & !truth)
      FN <- sum(!rej & truth)
      TN <- sum(!rej & !truth)
      
      n_rej <- TP + FP
      
      obs_fdr <- if (n_rej == 0) NA_real_ else FP / n_rej
      obs_power <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
      
      data.table(     
        TP     = TP,
        FP     = FP,
        FN     = FN,
        TN     = TN,
        alpha  = alpha,
        FDR    = obs_fdr,
        Power  = obs_power
      )
    })
    rbindlist(res)
  },
  by = g_vars
]

# Save
saveRDS(confusion_matrix, file = file.path(in_dir, "confusion_matrix.rds"))

# Already computed metric
res_path <- file.path(in_dir, "confusion_matrix.rds")
res      <- readRDS(res_path)
res_mem  <- object.size(res)
message(sprintf("Loaded file: %s uses %d RAM.", res_path, res_mem))

## Figure 1A ==================================================================
### fun. ----------------------------------------------------------------------

plot_imbalance_strip <- function(
  between_ratio = c(1, 5, 10),
  within_ratio  = c(1, 5, 10),
  ratios = NULL,
  N_fixed = 100,
  ncols = 1
) {

  # Grid of ratio combinations
  combos <- CJ(br = between_ratio, wr = within_ratio)

  # Ordered ratio labels
  combos[, ratio_label := paste0("br = ", br, ", wr = ", wr)]
  combos[, ratio_label := factor(
    ratio_label, levels = unique(ratio_label[order(br, wr)])
  )]

  # Optional filter for user-specified ratio labels
  if (!is.null(ratios)) {
    combos <- combos[ratio_label %in% ratios]
    if (nrow(combos) == 0) {
      stop("No matching ratio_label found in the provided grid.")
    }
  }

  # Helper: convert ratio → proportion
  r_to_prop <- function(r) r / (r + 1)

  # Compute synthetic counts 
  mini_base <- combos[, {
    p_between <- r_to_prop(br)   # majority share
    p_within  <- r_to_prop(wr)   # case share within ancestry

    n_major <- round(N_fixed * p_between)
    n_minor <- N_fixed - n_major

    n_major_case <- round(n_major * p_within)
    n_major_ctrl <- n_major - n_major_case

    n_minor_case <- round(n_minor * p_within)
    n_minor_ctrl <- n_minor - n_minor_case

    data.table(
      ancestry = factor(
        rep(c("maj", "min"), each = 2L), 
        levels = c("maj", "min")
      ),
      condition = factor(
        rep(c("case", "control"), 2L), 
        levels = c("case", "control")
      ),
      count = c(n_major_case, n_major_ctrl, n_minor_case, n_minor_ctrl)
    )
  }, by = .(ratio_label)]

  # Replicate into dummy facets 
  rep_tbl <- data.table(dummy = factor(seq_len(ncols)))
  rep_tbl[, key := 1L]; mini_base[, key := 1L]
  mini <- merge(
    mini_base, 
    rep_tbl, 
    by = "key", 
    allow.cartesian = TRUE
  )[, key := NULL]

  # Convert counts to percentages
  mini[, percent := 100 * count / sum(count), by = .(ratio_label, dummy)]

  # Fixed order of sub-bars
  mini[, bar := factor(
    paste(ancestry, condition, sep = "_"), 
    levels = c("maj_case", "maj_control", "min_case", "min_control")
  )]
  mini[, xnum := as.numeric(
    factor(ratio_label, levels = unique(ratio_label))
  )]

  # Appearance parameters
  width_stack <- 0.35  
  nudge_amt   <- 0.25

  ggplot() +
    geom_col(
      data = mini[ancestry == "maj"],
      aes(x = xnum - nudge_amt, y = percent, fill = condition),
      width = width_stack, position = "dodge"
    ) +
    geom_col(
      data = mini[ancestry == "min"],
      aes(x = xnum + nudge_amt, y = percent, fill = condition),
      width = width_stack, position = "dodge"
    ) +
    scale_x_continuous(
      breaks = sort(unique(mini$xnum)),
      labels = unique(mini$ratio_label),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0.05, 0.02)),
      breaks = seq(0, 100, by = 25),
      labels = function(x) paste0(x, "%")
    ) +
    scale_fill_manual(
      values = c(
        "case" = "#999999ff", 
        "control" = "#d40000ff"
      )
    ) +
    facet_grid(
      cols = vars(dummy, ratio_label), 
      scales = "free_x"
    ) +
    labs(
      x = "Imbalance ratio", 
      y = "Fraction"
    ) +
    theme_CrossAncestryGenPhen(
      show_facets = FALSE
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
}

### Imbalance structure -------------------------------------------------------

imbalance_structure_p_I <- plot_imbalance_strip()

# Save
ggsaveDK(
  plot = imbalance_structure_p_I,
  file = file.path(fig_dir, "imbalance_structure_I.svg"),
  width = 8,
  height = 4,
  bg = "transparent"
)

### Heatmap FDR ---------------------------------------------------------------

# Mean results 
fdr_power_alpha <- res[
  , 
  .(
    FDR = mean(FDR, na.rm = TRUE),
    Power = mean(Power, na.rm = TRUE)
  ), 
  by = .(between_ratio, within_ratio, ratio_label, method, alpha)
]

# Deviation from nominal FDR
fdr_power_alpha[, FDR_dev := FDR - alpha]

# Alpha = 0.05 (raw)
heatmap_fdr_alpha0.05_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.05],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = FDR
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradientn(
  name = "Observed FDR",
  colours = c("white", "purple", "#2E004E"),
  limits = c(0, 1)
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "FDR\n(alpha 0.05)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_fdr_alpha0.05_p,
  file = file.path(fig_dir, "heatmap_fdr_alpha0.05.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.05 (deviation)
heatmap_fdr_dev_alpha0.05_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.05],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = FDR_dev
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradient2(
  name = "Deviation from\nclaimed FDR",
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  limits = c(
    -max(abs(fdr_power_alpha[alpha == 0.05, FDR - alpha])),
     max(abs(fdr_power_alpha[alpha == 0.05, FDR - alpha]))
  )
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "FDR\n(alpha 0.05)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_fdr_dev_alpha0.05_p,
  file = file.path(fig_dir, "heatmap_fdr_dev_alpha0.05.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.1 (raw)
heatmap_fdr_alpha0.1_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.1],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = FDR
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradientn(
  name = "Observed FDR",
  colours = c("white", "purple", "#2E004E"),
  limits = c(0, 1)
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "FDR\n(alpha 0.1)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_fdr_alpha0.1_p,
  file = file.path(fig_dir, "heatmap_fdr_alpha0.1.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.1 (deviation)
heatmap_fdr_dev_alpha0.1_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.1],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = FDR_dev
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradient2(
  name = "Deviation from\nclaimed FDR",
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  limits = c(
    -max(abs(fdr_power_alpha[alpha == 0.1, FDR - alpha])),
     max(abs(fdr_power_alpha[alpha == 0.1, FDR - alpha]))
  )
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "FDR\n(alpha 0.05)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_fdr_dev_alpha0.1_p,
  file = file.path(fig_dir, "heatmap_fdr_dev_alpha0.1.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

### Heatmap power -------------------------------------------------------------

# Alpha = 0.05
heatmap_power_alpha0.05_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.05],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = Power
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradientn(
    colours = c("white", "purple", "#2E004E"),
    limits = c(0, 1)
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "Power"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_power_alpha0.05_p,
  file = file.path(fig_dir, "heatmap_power_alpha0.05.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.1
heatmap_power_alpha0.1_p <- ggplot(
  data = fdr_power_alpha[alpha == 0.1],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = Power
  )
) +
geom_tile(
  color = "white",
  linewidth = 0.3
) +
scale_fill_gradientn(
    colours = c("white", "purple", "#2E004E"),
    limits = c(0, 1)
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "Power"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_axis = FALSE, 
  show_border = FALSE
)

# Save
ggsaveDK(
  plot = heatmap_power_alpha0.1_p,
  file = file.path(fig_dir, "heatmap_power_alpha0.1.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

### Patchwork -----------------------------------------------------------------

fig_1A_1 <- (
  imbalance_structure_p_I + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  heatmap_fdr_alpha0.05_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank())
) /
(
  heatmap_power_alpha0.05_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(heights = c(0.3, 1, 1))

# Save
ggsaveDK(
  plot = fig_1A_1,
  file = file.path(patch_dir, "p_heatmap_alpha0.05.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

fig_1A_2 <- (
  imbalance_structure_p_I + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  heatmap_fdr_alpha0.1_p + 
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank())
) /
(
  heatmap_power_alpha0.1_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(heights = c(0.3, 1, 1))

# Save
ggsaveDK(
  plot = fig_1A_2,
  file = file.path(patch_dir, "p_heatmap_alpha0.1.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

fig_1A_3 <- (
  imbalance_structure_p_I + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  heatmap_fdr_dev_alpha0.05_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank())
) /
(
  heatmap_power_alpha0.05_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(heights = c(0.3, 1, 1))

# Save
ggsaveDK(
  plot = fig_1A_3,
  file = file.path(patch_dir, "p_heatmap_dev_alpha0.05.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)


fig_1A_4 <- (
  imbalance_structure_p_I + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  heatmap_fdr_dev_alpha0.1_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank())
) /
(
  heatmap_power_alpha0.1_p +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(heights = c(0.3, 1, 1))

# Save
ggsaveDK(
  plot = fig_1A_4,
  file = file.path(patch_dir, "p_heatmap_dev_alpha0.1.svg"),
  width = 8,
  height = 8,
  bg = "transparent"
)

### Datatables ----------------------------------------------------------------

write.xlsx(
  fdr_power_alpha, 
  file = file.path(tab_dir, "Figure_1A_data.xlsx")
)

write.csv(
  fdr_power_alpha, 
  file = file.path(tab_dir, "Figure_1A_data.csv")
)

## Figure 1B ==================================================================
### Imabalance structure ------------------------------------------------------

imbalance_structure_p_II <- plot_imbalance_strip(
  ratios = c("br = 1, wr = 1", "br = 5, wr = 5", "br = 10, wr = 10")
)

# Save
ggsaveDK(
  plot = imbalance_structure_p_II,
  file = file.path(fig_dir, "imbalance_structure_II.svg"),
  width = 8,
  height = 4,
  bg = "transparent"
)

### FDR curve -----------------------------------------------------------------

# FDR per alpha
fdr_summary <- res[
  ,
  .(
    mean_FDR   = 100 * mean(FDR, na.rm = TRUE),       
    sd_FDR     = 100 * sd(FDR, na.rm = TRUE) / sqrt(.N),
    mean_Power = 100 * mean(Power, na.rm = TRUE),     
    sd_Power   = 100 * sd(Power, na.rm = TRUE) / sqrt(.N),
    n_sim       = .N
  ),
  by = .(method, ratio_label, alpha)
][, alpha := 100 * alpha]

# Ideal scenario
id_df <- data.frame(
  alpha = sort(unique(fdr_summary$alpha))
)
id_df$mean_FDR <- id_df$alpha
id_df$line <- "Ideal scenario"

# FDR curve
fdr_curve_p <- ggplot(
  data = fdr_summary[
    ratio_label %in% c("br = 1, wr = 1", "br = 5, wr = 5", "br = 10, wr = 10")
  ],
  mapping = aes(
    x = alpha, 
    y = mean_FDR
  )
) +
geom_ribbon(
  aes(
    ymin = alpha, 
    ymax = ceiling(
      max(fdr_summary$mean_FDR, na.rm = TRUE) / 10
    ) * 10
  ),
  fill = "grey50", 
  alpha = 0.1
) +
ggnewscale::new_scale_color() +
geom_line(
  data = id_df,
  aes(x = alpha, y = mean_FDR, color = line),
  linewidth = 0.4,
  inherit.aes = FALSE
) +
scale_color_manual(
  name = NULL,
  values = c("Ideal scenario" = "gray30")
) +
ggnewscale::new_scale_color() +
geom_line(
  aes(color = method),
  linewidth = 0.3
) +
geom_point(
  aes(color = method),
  size = 0.3
) +
geom_errorbar(
  aes(
    ymin = mean_FDR - sd_FDR,
    ymax = mean_FDR + sd_FDR,
    color = method
  ),
  width = 0.05,
  linewidth = 0.3
) +
scale_color_brewer(
  name = "Method",
  palette = "Set2"
) +
scale_x_continuous(
  trans = pseudo_log_trans(base = 10),
  breaks = c(0, 1, 5, 10, 20),
  labels = c("0", "1", "5", "10", "20") 
) +
scale_y_continuous(
  trans = pseudo_log_trans(base = 10),
  breaks = c(0, 5, 20, 50, 100),
  labels = c("0", "5", "20", "50", "100") 
) +
facet_grid(
  cols = vars(ratio_label)
) +
labs(
  x = "Claimed FDR (%)",
  y = "Observed FDR (%)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  show_facets = FALSE, 
  show_borders = TRUE
)

# Save
ggsaveDK(
  plot = fdr_curve_p,
  file = file.path(fig_dir, "fdr_curve.svg"),
  width = 8,
  height = 4,
  bg = "transparent"
)

### Patchwork -----------------------------------------------------------------

fig_1B <- (
  imbalance_structure_p_II + 
    theme(axis.line = ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks.x = element_blank())
) /
(
 fdr_curve_p +
  theme(legend.spacing.y = unit(0.3, "line"))
) +
plot_layout(heights = c(0.5, 1))

# Save
ggsaveDK(
  plot = fig_1B,
  file = file.path(patch_dir, "Figure_1B.svg"),
  width = 8,
  height = 4,
  bg = "transparent"
)

write.xlsx(
  fdr_summary, 
  file = file.path(tab_dir, "Figure_1B_data.xlsx")
)

write.csv(
  fdr_summary, 
  file = file.path(tab_dir, "Figure_1B_data.csv")
)

## Supplementary Figure 1A ====================================================
### Confusion matrix ----------------------------------------------------------

# Prepare long format
grouping_vars <- c(
  "study", "coef_id", "contrast", "coef_type", "g_1", "g_2", 
  "a_1", "a_2", "method", "n_samples", "n_degs", "log2fc", 
  "total_samples", "between_ratio", "within_ratio", "ratio_label"
)

res_long <- melt(
  res,
  id.vars = c(grouping_vars, "alpha"),
  variable.name = "metric",
  value.name = "value"
)

# Alpha = 0.05
conf_matr_alpha0.05_p <- ggplot(
  data = res_long[metric %in% c("TP", "FP", "TN", "FN") & alpha == 0.05, ],
  mapping = aes(
    x = ratio_label,
    y = value,
    color = method
  )
) +
geom_boxplot(
  position = "dodge", 
  outlier.size = 0.1,
  linewidth = 0.3
) +
scale_color_brewer(palette = "Set2") +
facet_grid(
  rows = vars(metric),
  cols = vars(ratio_label),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Value",
  color = "Method"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE,
  strip.text.x = element_blank(),
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines")
)

# Save
ggsaveDK(
  plot = conf_matr_alpha0.05_p,
  file = file.path(fig_dir, "conf_matr_alpha0.05.svg"),
  width = 16,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.1
conf_matr_alpha0.1_p <- ggplot(
  data = res_long[metric %in% c("TP", "FP", "TN", "FN") & alpha == 0.1, ],
  mapping = aes(
    x = ratio_label,
    y = value,
    color = method
  )
) +
geom_boxplot(
  position = "dodge", 
  outlier.size = 0.1,
  linewidth = 0.3
) +
scale_color_brewer(palette = "Set2") +
facet_grid(
  rows = vars(metric),
  cols = vars(ratio_label),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Value",
  color = "Method"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_borders = TRUE,
  strip.text.x = element_blank(),
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines")
)

# Save
ggsaveDK(
  plot = conf_matr_alpha0.1_p,
  file = file.path(fig_dir, "conf_matr_alpha0.1.svg"),
  width = 16,
  height = 8,
  bg = "transparent"
)

### FDR & Power ---------------------------------------------------------------

# Alpha = 0.05
fdr_power_alpha0.05_p <- ggplot(
  data = res_long[metric %in% c("FDR", "Power") & alpha == 0.05], 
  mapping = aes(
    x = ratio_label, 
    y = value, 
    color = method
  )
) +
scale_color_brewer(palette = "Set2") +
geom_hline(
  data = data.frame(
    metric = "FDR", 
    yint = 0.05
  ),
  mapping = aes(
    yintercept = yint, 
    linetype = " 0.05"
  ),
  linewidth = 0.3,
  color = "blue"
) +
scale_linetype_manual(
  values = c(
    " 0.05" = "dashed"
  )
) +
geom_boxplot(
  position = "dodge", 
  outlier.size = 0.1,
  linewidth = 0.3
) +
facet_grid(
  rows = vars(metric),
  cols = vars(ratio_label),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Value",
  color = "Method",
  linetype = "Target FDR"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  strip.text.x = element_blank()
) + guides(
  linetype = guide_legend(order = 1),
  color    = guide_legend(order = 2)
)

# Save
ggsaveDK(
  plot = fdr_power_alpha0.05_p,
  file = file.path(fig_dir, "fdr_power_alpha0.05.svg"),
  width = 16,
  height = 8,
  bg = "transparent"
)

# Alpha = 0.1
fdr_power_alpha0.1_p <- ggplot(
  data = res_long[metric %in% c("FDR", "Power") & alpha == 0.1], 
  mapping = aes(
    x = ratio_label, 
    y = value, 
    color = method
  )
) +
scale_color_brewer(palette = "Set2") +
geom_hline(
  data = data.frame(
    metric = "FDR", 
    yint = 0.1
  ),
  mapping = aes(
    yintercept = yint, 
    linetype = " 0.1"
  ),
  linewidth = 0.3,
  color = "blue"
) +
scale_linetype_manual(
  values = c(
    " 0.1" = "dashed"
  )
) +
geom_boxplot(
  position = "dodge", 
  outlier.size = 0.1,
  linewidth = 0.3
) +
facet_grid(
  rows = vars(metric),
  cols = vars(ratio_label),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Value",
  color = "Method",
  linetype = "Target FDR"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  strip.text.x = element_blank()
) + guides(
  linetype = guide_legend(order = 1),
  color    = guide_legend(order = 2)
)

# Save
ggsaveDK(
  plot = fdr_power_alpha0.1_p,
  file = file.path(fig_dir, "fdr_power_alpha0.1.svg"),
  width = 16,
  height = 8,
  bg = "transparent"
)

### Patchwork -----------------------------------------------------------------

# Alpha = 0.05
sub_fig_1A_1 <- (
  imbalance_structure_p_I + 
     theme(axis.ticks.x = element_blank()) +
    theme(plot.margin = margin(t = 0, r = 0, b = 2, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  conf_matr_alpha0.05_p + 
    theme(plot.margin = margin(t = 0, r = 0, b = 4, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.ticks.x = element_blank())
) /
(
  fdr_power_alpha0.05_p + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(
  heights = c(0.3, 1, 0.8),
  guides = "collect"
)

# Save
ggsaveDK(
  plot = sub_fig_1A_1,
  file = file.path(fig_dir, "Supplementary_figure_1A_1.svg"),
  width = 16,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha = 0.1
sub_fig_1A_2 <- (
  imbalance_structure_p_I + 
     theme(axis.ticks.x = element_blank()) +
    theme(plot.margin = margin(t = 0, r = 0, b = 2, l = 0)) +
    theme(panel.spacing.x = unit(0.5, "mm")) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  conf_matr_alpha0.1_p + 
    theme(plot.margin = margin(t = 0, r = 0, b = 4, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.ticks.x = element_blank())
) /
(
  fdr_power_alpha0.1_p + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank())
) +
plot_layout(
  heights = c(0.3, 1, 0.8),
  guides = "collect"
)

# Save
ggsaveDK(
  plot = sub_fig_1A_2,
  file = file.path(fig_dir, "Supplementary_figure_1A_2.svg"),
  width = 16,
  height = 8,
  trimmed = FALSE,
  bg = "transparent"
)

write.xlsx(
  res_long, 
  file = file.path(tab_dir, "Supplementary_figure_1A_1|2_data.xlsx")
)

write.csv(
  res_long, 
  file = file.path(tab_dir, "Supplementary_figure_1A_1|2_data.csv")
)


## Supplementary Figure 1B ====================================================
### FDR vs. Power curve -------------------------------------------------------

fdr_power_curve_p <- ggplot(
  data = fdr_summary, 
  mapping = aes(
    x = mean_FDR, 
    y = mean_Power, 
    color = method
  )
) +
geom_line(
  linewidth = 0.5
) +
# geom_point(
#   size = 0.5
# ) +
scale_color_brewer(palette = "Set2") +
# scale_x_continuous(
#   trans = pseudo_log_trans(base = 10),
#   breaks = c(0, 5, 20, 80)
# ) +
# scale_y_continuous(
#   trans = pseudo_log_trans(base = 10),
#   breaks = c(0, 5, 10, 30, 50)
# ) +
facet_grid(
  cols = vars(ratio_label)
) +
labs(
  x = "Observed FDR (%)",
  y = "Power (%)",
  color = "Method"
) +
theme_CrossAncestryGenPhen(
  legend = 1,
  show_facets = FALSE,
  show_grid = TRUE
)

# Save
ggsaveDK(
  plot = fdr_power_curve_p,
  file = file.path(fig_dir, "fdr_power_curve.svg"),
  width = 16,
  height = 4
)

### 0 DEG pval dist. ----------------------------------------------------------

null_dist_p <- ggplot(
  data = res_0degs, 
  mapping = aes(
    x = p_value
  )
) +
geom_histogram(
  mapping = aes(
    y = log10(after_stat(count) + 1)
  ),
  binwidth = 1/20,
  fill = "gray80",
  color = "black",
  linewidth = 0.1
) +
geom_density(
  aes(
    y = log10(after_stat(density) * nrow(res0) * 1/20 + 1),
    color = ratio_label
  ),
  linewidth = 0.3
) +
facet_grid(
  cols = vars(method), 
  scales = "free_y"
) +
scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
labs(
  x = "P-value",
  y = "Log10 Count",
  color = "Imabalance ratios"
) +
theme_CrossAncestryGenPhen()

# Save
ggsaveDK(
  plot = null_dist_p,
  file = file.path(fig_dir, "null_dist.pdf"),
  width = 16,
  height = 4
)