## Review =====================================================================

## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(patchwork)
    library(openxlsx)
  }
)

## Directories ================================================================

# Figures
fig_dir <- file.path("results_major_review", "sims", "figures", "sim_synthetic")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Tables
tab_dir <- file.path("results_major_review", "sims", "tables", "sim_synthetic")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

alpha <- 0.5
formats <- c("svg", "png")

## Results ====================================================================
in_dir  <- file.path("results", "sims", "synthetic_sim")

# Already computed metric
res_path <- file.path(in_dir, "confusion_matrix.rds")
res      <- readRDS(res_path)
res_mem  <- object.size(res)
message(sprintf("Loaded file: %s uses %d RAM.", res_path, res_mem))

## R1.12 ----------------------------------------------------------------------
## Panel A --------------------------------------------------------------------
### Original paper figure -----------------------------------------------------
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

imbalance_structure_p_I <- plot_imbalance_strip()


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

# Heatmap deviation
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

# Heatmap power
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

main1_panel_C <- (
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

### Reviewed paper figure -----------------------------------------------------
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
  }, by = .(ratio_label, br)]

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
  
  # xnum wird innerhalb jeder Facette vergeben (1, 2, 3...)
  mini[, xnum := as.numeric(factor(ratio_label, levels = unique(ratio_label))), by = .(br)]

  # Erzeugt einen Lookup-Vector für die X-Achsen-Beschriftung
  label_lookup <- setNames(as.character(unique(mini$ratio_label)), unique(mini$xnum))

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
      # Nutzt eine Funktion basierend auf dem Lookup-Vector, um Längenkonflikte pro Facette zu umgehen
      labels = function(x) label_lookup[as.character(x)],
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
      cols = vars(br), 
      scales = "free_x",
      space = "free_x",
      labeller = labeller(br = function(x) paste0("Ancestry 1/", x))
    ) +
    labs(
      x = "Imbalance ratio", 
      y = "Fraction"
    ) +
    theme_CrossAncestryGenPhen(
      show_facets = TRUE
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
}

r1.12.1 <- plot_imbalance_strip()

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

# Alpha = 0.05 (deviation)
r1.12.2 <- ggplot(
  data = fdr_power_alpha[alpha == 0.05],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = FDR_dev
  )
) +
geom_tile(
  # color = "white",
  # linewidth = 0.1
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
facet_grid(
  cols = vars(between_ratio),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "FDR\n(alpha 0.05)"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45,
  show_facets = FALSE,
  show_axis = FALSE, 
  show_border = FALSE
)

r1.12.3 <- ggplot(
  data = fdr_power_alpha[alpha == 0.05],
  mapping = aes(
    x = ratio_label,
    y = method,
    fill = Power
  )
) +
geom_tile(
  # color = "white",
  # linewidth = 0.1
) +
scale_fill_gradientn(
    colours = c("white", "purple", "#2E004E"),
    limits = c(0, 1)
) +
scale_x_discrete(
  # Mappt das originale Faktor-Label auf deinen gewünschten Cancer-String
  labels = function(x) {
    # Extrahiert den numerischen Wert nach "wr = " und baut den String
    wr_vals <- sub(".*wr = ([0-9]+).*", "\\1", x)
    paste0("Cancer 1/", wr_vals)
  }
) +
facet_grid(
  cols = vars(between_ratio),
  scales = "free"
) +
labs(
  x = "Imbalance ratio",
  y = "Method",
  fill = "Power"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 45, 
  show_facets = FALSE,
  show_axis = FALSE, 
  show_border = FALSE
)

main1_panel_C_r1.12 <- (
  r1.12.1 + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.ticks.length = ggplot2::unit(0.5, "mm")) +
    theme(panel.spacing.x = unit(0.3, "mm")) +
    theme(axis.line= ggplot2::element_line(color = "black", linewidth = 0.3)) +
    theme(axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3))
) /
(
  r1.12.2 +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(panel.spacing.x = unit(-0.05, "cm")) +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank())
) /
(
  r1.12.3 +
    theme(legend.position = "right") +
    theme(legend.box.just = "left") +
    theme(legend.justification = "left") +
    theme(panel.spacing.x = unit(-0.05, "cm")) +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    theme(plot.margin = margin(t = 4, r = 0, b = 0, l = 0))
    # theme(axis.text.x = element_blank())
) +
plot_layout(heights = c(0.3, 1, 1))

# Save for paper
for (ext in formats) {
  ggsaveDK(
    plot = main1_panel_C_r1.12,
    file = file.path(
      fig_dir, 
      paste0("main2_panel_C_reviewed.", ext)
    ),
  width = 8,
  height = 8,
    bg = "transparent"
  )
}

# Patchwork
r1.12_panel_A <- wrap_plots(
  main1_panel_C, 
  main1_panel_C_r1.12, 
  ncol = 2
) +
plot_annotation(
  tag_levels = list(c("a", "", "", "b", "", ""))
) & 
theme(
  plot.tag = element_text(size = 10, face = "bold")
)

# Save
for (ext in formats) {
  ggsaveDK(
    plot = r1.12_panel_A,
    file = file.path(
      fig_dir, 
      paste0("r1.12_panel_A.", ext)
    ),
    height = 10,
    width = 16,
    bg = "transparent"
  )
}
