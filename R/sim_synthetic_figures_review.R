## max. code width ============================================================


## Libraries ==================================================================
suppressPackageStartupMessages({
  library(CrossAncestryGenPhen)
  library(data.table)
  library(ggnewscale)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(openxlsx)
})

## utils.R ====================================================================
source("R/utils.R")

## Parameters =================================================================
run_id <- "sims_n_10"

target_alpha <- 0.05
formats <- c("svg", "png")

# Data.table threads
setDTthreads(threads = 0)



## Paper ======================================================================
message("\nCreating paper figures...")
### paper directories ---------------------------------------------------------
fig_dir <- file.path("results_major_review", run_id, "figures", "sim_synthetic")
tab_dir <- file.path("results_major_review", run_id, "tables", "sim_synthetic")

# Figures
fig_dir_main <- file.path(fig_dir, "main")
fig_dir_supp <- file.path(fig_dir, "supp")
if (!dir.exists(fig_dir_main)) dir.create(fig_dir_main, recursive = TRUE)
if (!dir.exists(fig_dir_supp)) dir.create(fig_dir_supp, recursive = TRUE)

# Tables
tab_dir_supp <- file.path(tab_dir, "supp")
tab_dir_fig  <- file.path(tab_dir, "fig")
if (!dir.exists(tab_dir_fig)) dir.create(tab_dir_fig, recursive = TRUE)
if (!dir.exists(tab_dir_supp)) dir.create(tab_dir_supp, recursive = TRUE)



## Results ====================================================================
### res directories  ----------------------------------------------------------
res_dirs <- list(
  BRCA = file.path("results_major_review", run_id, "analysis", "TCGA_BRCA"),
  UCEC = file.path("results_major_review", run_id, "analysis", "TCGA_UCEC")
)
sum_dir <- file.path("results_major_review", run_id, "analysis", "summary_sim_synthetic")
stopifnot(length(unique(dirname(c(unlist(res_dirs), sum_dir)))) == 1)
# Create dir
if (!dir.exists(sum_dir)) dir.create(sum_dir, recursive = TRUE)

cfg_dirs <- list(
  BRCA = file.path("configs", "sims", "TCGA_BRCA"),
  UCEC = file.path("configs", "sims", "TCGA_UCEC")
)

### confusion matrix ----------------------------------------------------------
if (!file.exists(file = file.path(sum_dir, "summary_conf_mat.rds"))) {

  message("Preparing 'conf_mat' from individual results.")
  conf_mat <- rbindlist(lapply(names(res_dirs), function(study) {

    fp    <- file.path(res_dirs[[study]], "sim_synthetic")
    files <- list.files(fp, "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
    techs <- sub("_.*$", "", basename(dirname(files)))

    cfg   <- list.files(cfg_dirs[[study]])
    message(sprintf("%s found %d/%d files at: %s", study, length(files), length(cfg), fp))

    # List of files + tech
    dts <- Map(function(f,t) {
      dt <- as.data.table(readRDS(f))
      dt[, tech := t]
    }, files, techs)

    # Create confusion matrix
    conf_mats <- lapply(dts, function(dt) {

      # Imbalance 
      dt[, ratio_label := paste0("br = ", between_ratio, ", wr = ", within_ratio)]
      dt[, method_label := paste0(method, " (", strata, ")")]

      dt[, ratio_label := factor(ratio_label, levels = unique(ratio_label[order(between_ratio)]))]

      # Interaction only
      dt <- dt[n_degs > 0 & coef_type %in% c("interaction")]
      # dt <- dt[coef_type %in% c("interaction")]

      # Confusion matrix
      mat <- dt[, {
          res <- lapply(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2), function(alpha) {

            rej   <- p_adj <= alpha & !is.na(p_adj)
            truth <- T_truth != 0
            
            TP <- sum(rej & truth)
            FP <- sum(rej & !truth)
            FN <- sum(!rej & truth)
            TN <- sum(!rej & !truth)
            
            n_rej <- TP + FP
            
            fdr   <- if (n_rej == 0) NA_real_ else FP / n_rej
            power <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
            
            data.table(     
              TP     = TP,
              FP     = FP,
              FN     = FN,
              TN     = TN,
              alpha  = alpha,
              FDR    = fdr,
              Power  = power
            )
          })
          rbindlist(res)
        }, 
        by = c(
          "method_label", "method", "n_samples", "n_degs", "log2fc", "total_samples",
          "between_ratio", "within_ratio", "coef_id", "coef_type", "contrast",
          "g_1", "g_2", "a_1", "a_2", "strata", "ratio_label"
        )
      ]
    })

    out <- rbindlist(conf_mats)
    out[, study := study] 
    out
  }), use.names = TRUE, fill = TRUE)

  # Save
  saveRDS(
    object = conf_mat,
    file = file.path(
      sum_dir,
      "summary_conf_mat.rds"
    )
  )

} else {
  message("Loading 'conf_mat' from: ", sum_dir)
  conf_mat <- readRDS(file = file.path(sum_dir, "summary_conf_mat.rds"))
}


## Main Figure 2 ==============================================================
### panel C -------------------------------------------------------------------
mean_conf <- conf_mat[, .(
    FDR   = mean(FDR, na.rm = TRUE),
    Power = mean(Power, na.rm = TRUE)
  ), 
  by = .(between_ratio, within_ratio, ratio_label, method, strata, alpha)
][, `:=`(
  FDR_dev      = FDR - alpha,           # Absolute deviation
  FDR_dev_scal = (FDR - alpha) / alpha  # Scaled (relative)
)][, strata := factor(strata, levels = c("Imbalanced", "Ancestry balanced", "Ancestry + cancer balanced"))]

# Plot 
main2_panel_C <- wrap_plots(
  plot_imbalance_strip() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.length = ggplot2::unit(0.5, "mm"),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    plot.margin = margin(0, 0, 0, 0)
  ),
  ggplot() +
  geom_tile(
    data = mean_conf[
      alpha == target_alpha &
      !is.na(FDR_dev_scal)
    ],
    mapping = aes(
      x = ratio_label,
      y = method,
      fill = FDR_dev_scal
    ),
        color = "grey90",
    linewidth = 0.1
  ) +
  scale_fill_gradient2(
    name = "FDR deviation\n(x target FDR)",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    # transform = "pseudo_log",
    breaks = scales::pretty_breaks(n = 3),
    labels = function(x) {ifelse(x == 0, "0", sprintf("%+.fx", x))}
  ) +
  new_scale_fill() +
  geom_tile(
    data = mean_conf[
      alpha == target_alpha &
      is.na(FDR_dev_scal)
    ],
    mapping = aes(
      x = ratio_label, 
      y = method, 
      fill = "No positives"
    ),
    color = "grey90",
    linewidth = 0.1
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("No positives" = "grey80")
  ) +
  facet_grid(
    rows = vars(strata),
    cols = vars(between_ratio),
    labeller = labeller(
      strata = function(label) {
        gsub(" balanced", "\nbalanced", label)
      }
    ),
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Imbalance ratio",
    y = "Method",
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1, 
    rotate = 45,
    show_axis = FALSE, 
    show_border = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.00, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  ),
  ggplot(
    data = mean_conf[alpha == target_alpha],
    mapping = aes(
      x = ratio_label,
      y = method,
      fill = Power
    )
  ) +
  geom_tile(
    color = "grey90",
    linewidth = 0.1
  ) +
  geom_text(
    mapping = aes(
      label = sprintf("%.1f", Power)
    ),
    color = "black",
    size = 1 
  ) +
  scale_fill_gradientn(
    colours = c("white", "purple", "#2E004E"),
    limits = c(0, 1)
  ) +
  scale_x_discrete(
    labels = function(x) {
      wr_vals <- sub(".*wr = ([0-9]+).*", "\\1", x)
      paste0("Cancer 1/", wr_vals)
    }
  ) +
  facet_grid(
    rows = vars(strata),
    cols = vars(between_ratio),
    labeller = labeller(
      strata = function(label) {
        gsub(" balanced", "\nbalanced", label)
      }
    ),
    scales = "free",
    space = "free"
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
  ) +
  theme(
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.00, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(t = 4, 0, 0, 0)
  ),
  ncol = 1
) + plot_layout(heights = c(0.3, 1, 1))

# Save
for (ext in formats) {
  ggsaveDK(
    plot = main2_panel_C,
    file = file.path(
      fig_dir_main, 
      paste0("main2_panel_C.", ext)
    ),
    width = 10,
    height = 12,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main2_panel_C,
  file = file.path(
    fig_dir_main, 
    "main2_panel_C.rds"
  )
)

### panel D -------------------------------------------------------------------
main2_panel_D <- wrap_plots(
  plot_imbalance_strip(
    ratios = c("br = 1, wr = 1", "br = 5, wr = 5", "br = 10, wr = 10")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.length = ggplot2::unit(0.5, "mm"),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    plot.margin = margin(0, 0, 0, 0)
  ),
  ggplot(
    data = conf_mat[
      ratio_label %in% c("br = 1, wr = 1", "br = 5, wr = 5", "br = 10, wr = 10")
    ][, strata := factor(strata, levels = c("Imbalanced", "Ancestry balanced", "Ancestry + cancer balanced"))]
  ) +
  geom_ribbon(
    data = ~ {
      unique_alpha <- sort(unique(.x$alpha))
      max_fdr <- max(.x$FDR, na.rm = TRUE) * 1.05
      
      data.frame(
        alpha = unique_alpha,
        ymin  = unique_alpha,
        ymax  = max_fdr
      )
    },
    mapping = aes(
      x = alpha, 
      ymin = ymin, 
      ymax = ymax
    ),
    fill = "grey50", 
    alpha = 0.1,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = ~ data.frame(alpha = sort(unique(.x$alpha))),
    mapping = aes(
      x = alpha, 
      y = alpha, 
      color = "Ideal scenario"
    ),
    linetype = "dashed",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Ideal scenario" = "gray30")
  ) +
  ggnewscale::new_scale_color() +
  stat_summary(
    mapping = aes(
      x = alpha,
      y = FDR,
      color = method,
      group = method
    ),
    fun = mean,
    geom = "line",
    linewidth = 0.3
  ) +
  stat_summary(
    mapping = aes(
      x = alpha,
      y = FDR,
      color = method
    ),
    fun = mean,
    geom = "point",
    size = 0.3
  ) +
  scale_color_brewer(
    name = "Method",
    palette = "Set2"
  ) +
  scale_x_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 0.01, 0.05, 0.1, 0.2),
    labels = scales::label_percent(accuracy = 1, suffix = "")
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 0.05, 0.2, 0.5, 1),
    labels = scales::label_percent(accuracy = 1, suffix = "")
  ) +
  facet_grid(
    rows = vars(strata),
    cols = vars(between_ratio),
    labeller = labeller(
      strata = function(label) {
        gsub(" balanced", "\nbalanced", label)
      }
    ),
    space = "free"
  ) +
  labs(
    x = "Claimed FDR (%)",
    y = "Realized FDR (%)"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1, 
    show_facets = TRUE, 
    show_borders = TRUE
  ) +
  theme(
    strip.text.x    = element_blank(),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin   = margin(0, 0, 0, 0),
    plot.margin     = margin(t = 4, 0, 0, 0)
  ),
  ncol = 1
) + plot_layout(heights = c(0.3, 1))

# Save
for (ext in formats) {
  ggsaveDK(
    plot = main2_panel_D,
    file = file.path(
      fig_dir_main, 
      paste0("main2_panel_D.", ext)
    ),
    width = 10,
    height = 8,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main2_panel_D,
  file = file.path(
    fig_dir_main, 
    "main2_panel_D.rds"
  )
)




