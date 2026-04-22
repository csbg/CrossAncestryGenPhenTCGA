## max. code width ============================================================
## Libraries ==================================================================

library(CrossAncestryGenPhen)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(yardstick)
library(patchwork)
library(ggnewscale)
library(openxlsx)
library(msigdbr)
library(ggplot2)
library(ggbreak)
library(scales)
library(edgeR)
library(Rtsne)
library(fgsea)
library(ggh4x)
library(TRADEtools)

## Reusables ==================================================================
### fun. ----------------------------------------------------------------------

vs_newline <- function(x) {
  gsub("\\s*vs\\s*", " vs\n", x)
}

### Colors --------------------------------------------------------------------

ancestry_colors <- c(
  "EUR" = "#0072B2",
  "AFR" = "#D55E00",
  "EAS" = "#56B4E9",
  "AMR" = "#E69F00",
  "SAS" = "#009E73",
  "ADMIX" = "#999999"
)

subtype_colors <- c(
  "LumA" = "#8DB6CD",
  "LumB" = "#e31a1c",
  "Her2" = "#33a02c",
  "Basal" = "#ff7f00",
  "Normal" = "#CD00CD"
)

effect_group_colors <- c(
  "Divergent effect" = "#D62728",  
  "Minor interaction effect" = "black", 
  "non-EUR enhances\ncancer effect" = "#1F77B4",  
  "non-EUR reverts\ncancer effect"  = "#AEC7E8",  
  "Cancer enhances\nancestry effect"= "#98DF8A",  
  "Cancer reverts\nancestry effect" = "#2CA02C"
)

### Databases ----------------------------------------------------------------

# Prepare probe mapping
probe2gene <- fread("download/GDCportal/TCGA_BRCA/probe2gene.csv")
probe2gene <- setnames(probe2gene, "probe", "feature")

# Prepare database
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene <- as.data.table(hallmark)[, .(gs_name, gene_symbol)]
term2gene[, gs_name := {
    x <- sub("^HALLMARK_", "", gs_name)   
    x <- gsub("_", " ", x)         
    x <- tools::toTitleCase(tolower(x))
    x
}]

# FGSEA database
fgsea_pathways <- split(term2gene$gene_symbol, term2gene$gs_name)

# Message
cat("Created databases for enrichment -----------------------------------\n\n")

## Directories ================================================================

# Figures
fig_dir <- file.path("figures", "subset_interaction_effect")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Patchwork
patch_dir <- file.path(fig_dir, "patch_fig")
if (!dir.exists(patch_dir)) dir.create(patch_dir, recursive = TRUE)

# Excel sheets
tab_dir <- file.path("tables", "subset_interaction_effect")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

# Result directories (per cance study)
res_dirs <- list(
  BRCA = "results/TCGA_BRCA",
  UCEC = "results/TCGA_UCEC",
  THCA = "results/TCGA_THCA"
)

# Message
cat("Created directories to save data -----------------------------------\n\n")

## Load results ===============================================================
### Subset interaction results ------------------------------------------------

subset_interaction <- rbindlist(lapply(names(res_dirs), function(study) {
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

# Message
cat("Loaded subset results ----------------------------------------------\n\n")

### Iteration interaction results ---------------------------------------------

iteration_interaction <- rbindlist(lapply(names(res_dirs), function(study) {
  
  file_path <- file.path(res_dirs[[study]], "subset_interaction_effect")
  dge_files <- list.files(file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  # Load files + extract tech from folder name
  all_stats <- lapply(dge_files, function(f) {
    x <- readRDS(f)
    folder <- basename(dirname(f))             
    tech   <- strsplit(folder, "_")[[1]][1]    
    stats  <- as.data.table(x$subsets_stats)    
    stats[, tech := tech]                     
    stats                                      
  })
  # Bind results
  dt <- rbindlist(all_stats, use.names = TRUE, fill = TRUE)[coef_id == "interaction"]
  # Prepare phenotype
  dt[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(dt$phenotype)
  dt[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  dt[, study := study]
  dt
}), fill = TRUE)

## Figure 1A ==================================================================
### fun. ----------------------------------------------------------------------

plot_nr_degs_heatmap <- function(
  data,
  coef,
  legend1,
  legend2,
  fill_c = "purple"
) {

  p <- ggplot(
      data = data,
      mapping = aes(
        x = a_2, 
        y = phenotype
      )
    ) +
    geom_tile(
      data = data[
        coef_id %in% coef & 
        fill_group == " "
      ],
      aes(fill = fill_group),
      color = "black",
      linewidth = 0.1,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(" " = "grey90"),
      breaks = " ",
      labels = legend1,
      guide = guide_legend(
        title = NULL,
        override.aes = list(color = "black")
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = data[
        coef_id %in% coef & 
        fill_group == "effect"
      ],
      aes(fill = fill_value),
      color = "black",
      linewidth = 0.1
    ) +
    scale_fill_gradientn(
      name = legend2,
      colours = c("white", fill_c),
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 2),
      labels = scales::math_format(10^.x),
      guide = guide_colorbar(
        order = 1,
        title.hjust = 0
      )
    ) +
    facet_grid(
      cols = vars(tech),
      rows = vars(study),
      scales = "free_y",
      space = "free_y",
      labeller = labeller(
        tech = c(meth = "Methylation", mrna = "Expression"),
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
      legend.position = "bottom",
      legend.title.position = "left",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.box.spacing = unit(0, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    )
}

plot_nr_degs_barplot <- function(
  data
) {

  ggplot(
    data = data,
    mapping = aes(
      x = vs_newline(phenotype),
      y = log10(nr_sig_genes + 1),
      fill = a_2
    )
  ) +
  geom_col(
    position = position_dodge2(
      preserve = "single",
      width = 1
    ),
    color = "black",
    linewidth = 0.1,
  ) +
  geom_text(
    aes(
      label = nr_sig_genes,
      y = log10(nr_sig_genes + 1),
      color = a_2
    ),
    position = position_dodge2(
      preserve = "single",
      width = 1
    ),
    show.legend = FALSE,
    # angle = 90,
    # hjust = -0.4,
    vjust = -0.5,
    size  = 2,
    na.rm = TRUE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.5)),
    labels = scales::math_format(10^.x),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  facet_grid(
    cols = vars(tech),
    scales = "free_y",
    space  = "free_y",
    labeller = labeller(
      tech = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  theme_CrossAncestryGenPhen(
    base_size = 5,
    legend_key = 1,
    rotate = 45
  ) +
  theme(
    strip.text = element_text(margin = margin(1, 1, 1, 1)),
    plot.title = element_text(margin = margin(b = 1)),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.spacing = unit(0, "pt"),
    legend.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    axis.text.x = element_text(margin = margin(t = 0)),
    axis.title.x = element_text(margin = margin(t = 0)),
    axis.text.y = element_text(margin = margin(t = 0)),
    axis.title.y = element_text(margin = margin(t = 0))
  )
}

plot_volcano <- function(
  data
){

  ggplot(
    data = data,
    mapping = aes(
      x = T_obs,
      y = -log10(p_adj)
    )
  ) + 
  geom_point(
    size = 0.1,
    alpha = 0.5
  ) +
  facet_grid(
    cols = vars(tech),
    rows = vars(phenotype, a_2),
    # scales = "free_y",
    # space  = "free_y",
    labeller = labeller(
      tech = c(meth = "Methylation", mrna = "Expression"),
      phenotype = vs_newline
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(
    x = "Log2FC",
    y = expression(-log[10]("adj. p-value"))
  ) +
  theme_CrossAncestryGenPhen(
    base_size = 5,
    show_border = TRUE,
    show_grid = TRUE
  ) +
  theme(
    strip.text = element_text(margin = margin(1, 1, 1, 1)),
    plot.title = element_text(margin = margin(b = 1)),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    axis.text.x = element_text(margin = margin(t = 0)),
    axis.title.x = element_text(margin = margin(t = 0))
  )
}

FGSEA_one <- function(
  ranks,
  metadata,
  pathways
) {

  res <- fgsea(
    pathways = pathways,
    stats = ranks
  )

  # If returned nothing
  if (is.null(res) || nrow(res) == 0) return(NULL)

  # Convert to data.table
  dt <- as.data.table(res)

  # Add metadata columns
  dt[, study := metadata$study]
  dt[, phenotype := metadata$phenotype]
  dt[, tech := metadata$tech]
  dt[, a_2 := metadata$a_2]

  return(dt)
}

plot_fgsea_dotplot <- function(
  data
){

  ggplot(
    data = data,
    mapping = aes(
      x = a_2,
      y = reorder(pathway, NES),
      size  = -log10(padj),
      color = NES
    )
  ) +
  geom_point() +
  scale_color_gradient2(
    name = "Normalized\nenrichment",
    low  = "blue",
    mid  = "white",
    high = "red"
  ) +
  scale_size_continuous(
    name = expression(p.adj ~ "-log"[10] * ""),
    range = c(0.8, 2),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  facet_grid(
    cols = vars(phenotype),
    rows = vars(tech),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      phenotype = vs_newline,
      tech = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    x = "Ancestry",
    y = "Hallmark pathway"
  ) +
  theme_CrossAncestryGenPhen(
    base_size = 6,
    legend_key = 1,
    rotate = 45,
    show_border = TRUE,
    show_grid = TRUE
  ) +
  theme(
    plot.margin = margin(r = 1, 0, 0, 0),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.position  = "right",
    legend.direction = "vertical"
  )
}

TRADE_one <- function(
  stats,
  metadata
){

  res <- TRADE(
    mode = "univariate",
    results1 = stats
  )

  # If returned nothing
  if (is.null(res))  return(NULL)

  # data.table
  ti <- res$distribution_summary$transcriptome_wide_impact
  ti_scaled <- ti / dim(stats)[1]
  dt <- data.table(
    ti = ti,
    ti_scaled = ti_scaled
  )

  # Add metadata columns
  dt[, study := metadata$study]
  dt[, phenotype := metadata$phenotype]
  dt[, tech := metadata$tech]
  dt[, a_2 := metadata$a_2]

  return(dt)
}

plot_trade_heatmap <- function(
  data,
  legend
) {

  p <- ggplot(
      data = data,
      mapping = aes(
        x = a_2, 
        y = phenotype
      )
    ) +
    geom_tile(
      aes(fill = ti),
      color = "black",
      linewidth = 0.1,
      show.legend = TRUE
    ) +
    scale_fill_gradientn(
      name = legend,
      colours = c("white", "gold"),
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 2),
      guide = guide_colorbar(
        order = 1,
        title.hjust = 0
      )
    ) +
    facet_grid(
      cols = vars(tech),
      rows = vars(study),
      scales = "free_y",
      space = "free_y",
      labeller = labeller(
        tech = c(meth = "Methylation", mrna = "Expression"),
        study = label_value
      )
    ) +
    labs(
      x = "Ancestry",
      y = "Phenotype"
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
}

### [Interaction effect] Heatmap nr. sig. genes -------------------------------

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
      ][
        ,
        `:=`(
          fill_group = fifelse(nr_sig_genes == 0, " ", "effect"),
          fill_value = log10(nr_sig_genes + 1)
        )
      ]
    }
  )
)

write.xlsx(
  dt_degs_nr, 
  file = file.path(tab_dir, "dt_degs_nr.xlsx")
)

write.csv(
  dt_degs_nr, 
  file = file.path(tab_dir, "dt_degs_nr.csv")
)

# Alpha 0.05
gg_heatmap_degs_nr_inter_alpha0.05 <- plot_nr_degs_heatmap(
  data = dt_degs_nr[alpha == 0.05],
  coef = "interaction",
  legend1 = "No effect ",
  legend2 = "# DEGs with\ninteraction effect",
  fill_c = "forestgreen"
) +
theme(legend.position = "right") +
theme(legend.title.position = "top") +
theme(legend.spacing.y = unit(10, "pt")) 

# Save
ggsaveDK(
  plot = gg_heatmap_degs_nr_inter_alpha0.05,
  file = file.path(fig_dir, "heatmap_nr_degs_inter_alpha0.05.svg"),
  height = 5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha 0.1
gg_heatmap_degs_nr_inter_alpha0.1 <- plot_nr_degs_heatmap(
  data = dt_degs_nr[alpha == 0.1],
  coef = "interaction",
  legend1 = "No effect",
  legend2 = "# DEGs with\ninteraction effect",
  fill_c = "forestgreen"
) +
theme(legend.position = "right") +
theme(legend.title.position = "top") +
theme(legend.spacing.y = unit(10, "pt")) 

# Save
ggsaveDK(
  plot = gg_heatmap_degs_nr_inter_alpha0.1,
  file = file.path(fig_dir, "heatmap_nr_degs_inter_alpha0.1.svg"),
  height = 5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### [Interaction effect] Bar plot nr. sig. genes ------------------------------

# Alpha = 0.05
barplots_nr_degs_inter_alpha0.05 <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_nr_degs_barplot(
      data = dt_degs_nr[
        coef_id == "interaction" &
        study == s &
        alpha == 0.05
      ]
    )
    
    # Color + labs
    p <- p +
    scale_fill_manual(
      values = ancestry_colors
    ) +
    scale_color_manual(
      values = ancestry_colors
    ) +
    labs(
      title = s,
      x = NULL,
      y = "# genes\n(alpha 0.05)",
      fill = NULL
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("barplot_nr_degs_inter_alpha0.05_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 5,
      width = 6.5,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

# Alpha = 0.1
barplots_nr_degs_inter_alpha0.1 <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_nr_degs_barplot(
      data = dt_degs_nr[
        coef_id == "interaction" &
        study == s &
        alpha == 0.1
      ]
    ) 
    
    # Color + labs
    p <- p +
    scale_fill_manual(
      values = ancestry_colors
    ) +
    scale_color_manual(
      values = ancestry_colors
    ) +
    labs(
      title = s,
      x = NULL,
      y = "# genes\n(alpha 0.1)",
      fill = NULL
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("barplot_nr_degs_inter_alpha0.1_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 5,
      width = 6.5,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Interaction effect] Volcano plots ----------------------------------------

volcanoplots_inter <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_volcano(
      data = subset_interaction[
        coef_id == "interaction" & 
        study == s
      ]
    )

    # Title
    p <- p + ggtitle(s)

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("volcanoplot_inter_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 8,
      width = 8,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Interaction effect] FGSEA analysis (TOP 5 pahtways) ----------------------

fgsea_res_inter <- copy(subset_interaction)[
  coef_id == "interaction"
][
  tech == "meth",
  feature := probe2gene[.SD, on = "feature", gene]
][
  feature %in% term2gene$gene_symbol
][,
  {

    # Ranking features by logFC
    ranks <- setNames(T_obs, feature)
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- ranks[!is.na(names(ranks))]

    # Sort decreasing (required)
    ranks <- sort(ranks, decreasing = TRUE)

    meta <- .SD[1, .(study, phenotype, tech, a_2)]
    res  <- FGSEA_one(
      ranks = ranks,
      metadata = meta,
      pathways = fgsea_pathways
    )

    .(res = list(res))

  },
  by = .(study, phenotype, tech, a_2)
]
fgsea_res_inter <- rbindlist(fgsea_res_inter$res)

# Save
write.xlsx(
  fgsea_res_inter[,
    `:=`(
      leadingEdge = sapply(leadingEdge, paste, collapse = ";"),
      leadingEdge_size = lengths(leadingEdge)
    )
  ], 
  file = file.path(tab_dir, "fgsea_interaction_effects.xlsx")
)

write.csv(
  fgsea_res_inter[,
    `:=`(
      leadingEdge = sapply(leadingEdge, paste, collapse = ";"),
      leadingEdge_size = lengths(leadingEdge)
    )
  ], 
  file = file.path(tab_dir, "fgsea_interaction_effects.csv")
)

gg_dotplot_fgsea_interaction <- lapply(
  names(res_dirs),
  function(s) {

    # Plot
    p <- plot_fgsea_dotplot(
      data = fgsea_res_inter[
        padj < 0.05 &
        study == s
      ][
        order(-abs(NES))
      ][
        ,
        head(.SD, 5),
        by = .(phenotype, tech, a_2)
      ]
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("dotplot_fgsea_inter_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 9,
      width = 9,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

# Barplot
gg_barplot_fgsea_interaction <- lapply(
  names(res_dirs),
  function(s) {

    # Plot
    p <- ggplot(
      data = fgsea_res_inter[
        study == s & 
        padj < 0.05
      ][
        order(-abs(NES))
      ][,
        head(.SD, 5),
        by = .(phenotype, tech, a_2)
      ],
      mapping = aes(
        x = NES,
        y = reorder(pathway, NES),
        fill = a_2
      )
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      color = "grey80"
    ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.1
    ) +
    scale_fill_manual(
      values = ancestry_colors
    ) +
    facet_grid(
      cols = vars(phenotype),
      rows = vars(tech),
      scales = "free_y",
      space = "free",
      labeller = labeller(
        phenotype = vs_newline,
        tech = c(meth = "Methylation", mrna = "Expression")
      )
    ) +
    labs(
    fill = "Ancestry",
    x = "Normalized\nenrichment score",
    y = "Hallmark pathway"
    ) +
    theme_CrossAncestryGenPhen(
      legend_key = 1,
      show_border = TRUE,
      show_grid = TRUE
    ) +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      panel.spacing.y = unit(0.15, "lines"),
      legend.margin = margin(0, 0, 0, 0)
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("gg_barplot_fgsea_interaction_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 10,
      width = 10,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Interaction effect] TRADE analysis ---------------------------------------

# trade_res_inter <- copy(subset_interaction)[
#   coef_id == "interaction"
# ][,
#   `:=`(
#     log2FoldChange = T_obs,
#     lfcSE = SE,
#     pvalue = p_value
#   )
# ][,
#   {
#     stats <- .SD[, .(
#       log2FoldChange,
#       lfcSE,
#       pvalue
#     )]
#     rownames(stats) <- .SD$feature

#     meta <- .SD[1, .(study, phenotype, tech, a_2)]
#     res  <- TRADE_one(
#       stats = stats,
#       metadata = meta
#     )

#     .(res = list(res))

#   },
#   by = .(study, phenotype, tech, a_2)
# ]
# trade_res_inter <- rbindlist(trade_res_inter$res)

# # Plot (heatmap)
# heatmap_trade_inter <- plot_trade_heatmap(
#   data = trade_res_inter,
#   legend = "TI"
# )

# # Save
# ggsaveDK(
#   plot = heatmap_trade_inter,
#   file = file.path(fig_dir, "heatmap_trade_inter.svg"),
#   height = 5,
#   width = 6.5,
#   trimmed = TRUE,
#   bg = "transparent"
# )

### [Baseline effect] Heatmap nr. sig. genes ----------------------------------

# Alpha = 0.05
gg_heatmap_degs_dr_base_alpha0.05 <- plot_nr_degs_heatmap(
  data = dt_degs_nr[alpha == 0.05],
  coef = c("baseline_1", "baseline_2"),
  legend1 = "No effect",
  legend2 = "# DEGs with\nbaseline effect"
)

# Save
ggsaveDK(
  plot = gg_heatmap_degs_dr_base_alpha0.05,
  file = file.path(fig_dir, "heatmap_nr_degs_base_alpha0.05.svg"),
  height = 6.5,
  width = 6.5,
  trimmed = TRUE,
  bg = "transparent"
)

# Alpha = 0.1
gg_heatmap_degs_nr_base_alpha0.1 <- plot_nr_degs_heatmap(
  data = dt_degs_nr[alpha == 0.1],
  coef = c("baseline_1", "baseline_2"),
  legend1 = "No effect",
  legend2 = "# DEGs with\nbaseline effect"
)

# Save
ggsaveDK(
  plot = gg_heatmap_degs_nr_base_alpha0.1,
  file = file.path(fig_dir, "heatmap_nr_degs_base_alpha0.1.svg"),
  height = 6.5,
  width = 6.5,
  trimmed = TRUE,
  bg = "transparent"
)

### [Baseline effect] Bar plot nr. sig. genes =================================

# Alpha = 0.05
barplots_nr_degs_base_alpha0.05 <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_nr_degs_barplot(
      data = dt_degs_nr[
        coef_id %in% c("baseline_1", "baseline_2") &
        study == s &
        alpha == 0.05
      ]
    )
    
    # Color + labs
    p <- p +
    scale_fill_manual(
      values = ancestry_colors
    ) +
    scale_color_manual(
      values = ancestry_colors
    ) +
    labs(
      title = s,
      x = NULL,
      y = "# genes\n(alpha 0.05)",
      fill = NULL
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("barplot_nr_degs_base_alpha0.05_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 5,
      width = 6.5,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

# Alpha = 0.1
barplots_nr_degs_base_alpha0.1 <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_nr_degs_barplot(
      data = dt_degs_nr[
        coef_id %in% c("baseline_1", "baseline_2") &
        study == s &
        alpha == 0.1
      ]
    )
    
    # Color + labs
    p <- p +
    scale_fill_manual(
      values = ancestry_colors
    ) +
    scale_color_manual(
      values = ancestry_colors
    ) +
    labs(
      title = s,
      x = NULL,
      y = "# genes\n(alpha 0.1)",
      fill = NULL
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("barplot_nr_degs_base_alpha0.1_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 5,
      width = 6.5,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Baseline effect] Volcano plots ===========================================

volcanoplots_base <- lapply(
  names(res_dirs),
  function(s) {
    
    # Plot
    p <- plot_volcano(
      data = subset_interaction[
        coef_id %in% c("baseline_1", "baseline_2") &
        study == s
      ]
    )

    # Title
    p <- p + ggtitle(s)

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("volcanoplot_base_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 16,
      width = 8,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Baseline effect] FGSEA analysis (TOP 5) -------------------------------

fgsea_res_base <- copy(subset_interaction)[
  coef_id %in% c("baseline_1", "baseline_2")
][
  tech == "meth",
  feature := probe2gene[.SD, on = "feature", gene]
][
  feature %in% term2gene$gene_symbol
][,
  {

    # Ranking features by logFC
    ranks <- setNames(T_obs, feature)
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- ranks[!is.na(names(ranks))]

    # Sort decreasing (required)
    ranks <- sort(ranks, decreasing = TRUE)

    meta <- .SD[1, .(study, phenotype, tech, a_2)]
    res  <- FGSEA_one(
      ranks = ranks,
      metadata = meta,
      pathways = fgsea_pathways
    )

    .(res = list(res))

  },
  by = .(study, phenotype, tech, a_2)
]
fgsea_res_base <- rbindlist(fgsea_res_base$res)

# Plot
dotplot_fgsea_base <- lapply(
  names(res_dirs),
  function(s) {

    # Plot
    p <- plot_fgsea_dotplot(
      data = fgsea_res_base[
        padj < 0.05 &
        study == s
      ][
        order(-abs(NES))
      ][
        ,
        head(.SD, 5),
        by = .(phenotype, tech)
      ]
    )

    # Outdir
    file_p <- file.path(
      fig_dir,
      paste0("dotplot_fgsea_base_", s, ".svg")
    )

    # Save
    ggsaveDK(
      plot = p,
      file = file_p,
      height = 6,
      width = 12,
      trimmed = TRUE,
      bg = "transparent"
    )
  }
)

### [Baseline effect] TRADE analysis ---------------------------------------

trade_res_base <- copy(subset_interaction)[
  coef_id %in% c("baseline_1", "baseline_2")
][,
  `:=`(
    log2FoldChange = T_obs,
    lfcSE = SE,
    pvalue = p_value
  )
][,
  {
    stats <- .SD[, .(
      log2FoldChange,
      lfcSE,
      pvalue
    )]
    rownames(stats) <- .SD$feature

    meta <- .SD[1, .(study, phenotype, tech, a_2)]
    res  <- TRADE_one(
      stats = stats,
      metadata = meta
    )

    .(res = list(res))

  },
  by = .(study, phenotype, tech, a_2)
]
trade_res_base <- rbindlist(trade_res_base$res)

# Plot (heatmap)
heatmap_trade_base <- plot_trade_heatmap(
  data = trade_res_base,
  legend = "TI"
)

# Save
ggsaveDK(
  plot = heatmap_trade_base,
  file = file.path(fig_dir, "heatmap_trade_base.svg"),
  height = 5,
  width = 6.5,
  trimmed = TRUE,
  bg = "transparent"
)

### [Baseline & interaction effect] Nr. sig genes correlation -----------------

dt_alpha <- rbindlist(
  lapply(c(0.05, 0.1), function(a) {
    tmp <- copy(subset_interaction)
    tmp[, alpha := a]
    tmp
  })
)

baseline_counts <- dt_alpha[
  coef_type == "baseline",
  .(baseline_genes = uniqueN(feature[p_adj <= alpha])),
  by = .(study, tech, a_2, g_1, g_2, alpha)
]

interaction_counts <- dt_alpha[
  coef_type == "interaction",
  .(interaction_genes = uniqueN(feature[p_adj <= alpha])),
  by = .(study, tech, a_2, g_1, g_2, phenotype, alpha)
]

counts <- merge(
  baseline_counts,
  interaction_counts,
  by = c("study", "tech", "a_2", "g_1", "g_2", "alpha"),
  all = TRUE
)

# Add correlation of DEG nr.
counts[, `:=`(
  x = log10(baseline_genes + 1),
  y = log10(interaction_genes + 1)
)]

counts[, across_cor :=
  if (.N >= 3 && sd(x) > 0 && sd(y) > 0)
    cor(x, y, method = "pearson")
  else NA_real_,
  by = .(tech, alpha)
]

counts[, ancestry_cor :=
   if (.N >= 3 && sd(x) > 0 && sd(y) > 0)
    cor(x, y, method = "pearson")
  else NA_real_,
  by = .(tech, alpha, a_2)
]

across_cor <- unique(
  counts[, .(tech, alpha, across_cor)]
)[
  , .(
    tech,
    alpha,
    cor_type = "Across",
    cor_value = across_cor
  )
]

ancestry_cor <- unique(
  counts[, .(tech, alpha, a_2, ancestry_cor)]
)[
  , .(
    tech,
    alpha,
    cor_type = a_2,
    cor_value = ancestry_cor
  )
]

# Combine across & per ancestry correlation
cor_long <- rbind(across_cor, ancestry_cor)
cor_long[, cor_type := factor(
  cor_type,
  levels = c("Across", setdiff(sort(unique(cor_type)), "Across"))
)]
cor_long[, cor_group := ifelse(cor_type == "Across", "Across", "Ancestry")]
cor_long[, cor_group := factor(cor_group, levels = c("Across", "Ancestry"))]

# Alpha = 0.05
# Point plot
pointplot_nr_degs_corr_alpha0.05 <- ggplot(
  data = counts[alpha == 0.05],
  mapping = aes(
    x = log10(baseline_genes + 1),
    y = log10(interaction_genes + 1),
    color = a_2
  )
) +
geom_polygon(
  data = {
    d <- counts[alpha == 0.05]
    
    x <- log10(d$baseline_genes + 1)
    y <- log10(d$interaction_genes + 1)
    
    lims <- range(c(x, y))
    
    data.frame(
      x = c(lims[1], lims[1], lims[2]),
      y = c(lims[2], lims[1], lims[2])
    )
  },
  aes(x = x, y = y),
  inherit.aes = FALSE,
  fill = "grey80",
  alpha = 0.2
) +
geom_point(
  size = 0.5
) +
geom_smooth(
  mapping = aes(
    x = log10(baseline_genes + 1),
    y = log10(interaction_genes + 1)
  ),
  linewidth = 0.7,
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
  values = ancestry_colors
) +
labs(
  color = "Ancestry",
  x = "# baseline DEGs\n(alpha 0.05)",
  y = "# interaction DEGs\n(alpha 0.05)",
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

# Save (svg)
ggsaveDK(
  plot = pointplot_nr_degs_corr_alpha0.05,
  file = file.path(fig_dir, "pointplot_nr_degs_corr_alpha0.05.svg"),
  height = 8,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

# Barplot
barplot_nr_degs_corr_alpha0.05 <- ggplot(
  data = cor_long[alpha == 0.05],
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
scale_x_continuous(
  limits = c(-1, 1),
  breaks = c(-1, 0, 1)
) +
labs(
  fill = NULL,
  x = "Correlation # DEGs\n(Pearson)",
  y = "Ancestry"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save (svg)
ggsaveDK(
  plot = barplot_nr_degs_corr_alpha0.05,
  file = file.path(fig_dir, "barplot_nr_degs_corr_alpha0.05.svg"),
  height = 8,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

# Alpha = 0.1
# Point plot
pointplot_nr_degs_corr_alpha0.1 <- ggplot(
  data = counts[alpha == 0.1],
  mapping = aes(
    x = log10(baseline_genes + 1),
    y = log10(interaction_genes + 1),
    color = a_2
  )
) +
geom_polygon(
  data = {
    d <- counts[alpha == 0.1]
    
    x <- log10(d$baseline_genes + 1)
    y <- log10(d$interaction_genes + 1)
    
    lims <- range(c(x, y))
    
    data.frame(
      x = c(lims[1], lims[1], lims[2]),
      y = c(lims[2], lims[1], lims[2])
    )
  },
  aes(x = x, y = y),
  inherit.aes = FALSE,
  fill = "grey80",
  alpha = 0.2
) +
geom_point(
  size = 0.5
) +
geom_smooth(
  mapping = aes(
    x = log10(baseline_genes + 1),
    y = log10(interaction_genes + 1)
  ),
  method = "lm",
  se = FALSE,
  color = "blue",
  linewidth = 0.7,
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
  values = ancestry_colors
) +
labs(
  color = "Ancestry",
  x = "# baseline DEGs\n(alpha 0.1)",
  y = "# interaction DEGs\n(alpha 0.1)",
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

# Save (svg)
ggsaveDK(
  plot = pointplot_nr_degs_corr_alpha0.1,
  file = file.path(fig_dir, "pointplot_nr_degs_corr_alpha0.1.svg"),
  height = 8,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

# Barplot
barplot_nr_degs_corr_alpha0.1 <- ggplot(
  data = cor_long[alpha == 0.1],
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
scale_x_continuous(
  limits = c(-1, 1),
  breaks = c(-1, 0, 1)
) +
labs(
  fill = NULL,
  x = "Correlation # DEGs\n(Pearson)",
  y = "Ancestry"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save (svg)
ggsaveDK(
  plot = barplot_nr_degs_corr_alpha0.1,
  file = file.path(fig_dir, "barplot_nr_degs_corr_alpha0.1.svg"),
  height = 8,
  width = 8,
  trimmed = TRUE,
  bg = "transparent"
)

### Patchwork -----------------------------------------------------------------
#### Nr. degs -----------------------------------------------------------------

# Alpha = 0.05
p_degs_nr_alpha0.05 <- (
  gg_heatmap_degs_dr_base_alpha0.05 +
    theme(legend.position = "right") +
    theme(legend.title.position = "top")
    # theme(axis.title.x = element_blank()) +
    # theme(axis.text.x = element_blank()) 
) /
(
  gg_heatmap_degs_nr_inter_alpha0.05 +
    theme(legend.position = "right") +
    theme(legend.title.position = "top")
    # theme(strip.text.x = element_blank())
) +
plot_layout(
  heights = c(1, 0.65)
) &
theme(
  legend.spacing.y = unit(10, "pt"),
) 

# Save (svg)
ggsaveDK(
  plot = p_degs_nr_alpha0.05,
  file = file.path(patch_dir, "p_degs_nr_alpha0.05.svg"),
  height = 12,
  width = 9,
  trimmed = TRUE,
  bg = "transparent"
)

# Alpha = 0.1
p_nr_degs_alpha0.1 <- (
  gg_heatmap_degs_nr_base_alpha0.1 +
    theme(legend.position = "right") +
    theme(legend.title.position = "top")
    # theme(axis.title.x = element_blank()) +
    # theme(axis.text.x = element_blank()) 
) /
(
  gg_heatmap_degs_nr_inter_alpha0.1 +
    theme(legend.position = "right") +
    theme(legend.title.position = "top")
    # theme(strip.text.x = element_blank())
    # theme(axis.text.x = element_text(size = 5))
) +
plot_layout(
  heights = c(1, 0.65)
) &
theme(
  legend.spacing.y = unit(10, "pt"),
) 

# Save (svg)
ggsaveDK(
  plot = p_nr_degs_alpha0.1,
  file = file.path(patch_dir, "p_degs_nr_alpha0.1.svg"),
  height = 12,
  width = 9,
  trimmed = TRUE,
  bg = "transparent"
)

#### Nr. degs correlation -----------------------------------------------------

# Alpha = 0.05
p_cor_degs_alpha0.05 <- (
  barplot_nr_degs_corr_alpha0.05 +
    labs(x = "Correlation # DEGs") +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(plot.margin = margin(b = 15, 0, 0, 0)) +
    theme(legend.margin = margin(0, 0, 0, 0)) +
    theme(legend.position = "top")
) /
(
  pointplot_nr_degs_corr_alpha0.05 +
    labs(
      x = "# baseline DEGs",
      y = "# interaction DEGs",
    ) +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0, 0, 0))
) +
plot_layout(heights = c(0.5, 1))

# Save (svg)
ggsaveDK(
  plot = p_cor_degs_alpha0.05,
  file = file.path(patch_dir, "p_degs_cor_alpha0.05.svg"),
  height = 10,
  width = 4,
  trimmed = TRUE,
  bg = "transparent"
)

# Alpha = 0.1
p_cor_degs_alpha0.1 <- (
  barplot_nr_degs_corr_alpha0.1 +
    labs(x = "Correlation # DEGs") +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    theme(legend.margin = margin(0, 0, 0, 0)) +
    theme(legend.position = "top")
) /
(
  pointplot_nr_degs_corr_alpha0.1 +
    labs(
      x = "# baseline DEGs",
      y = "# interaction DEGs",
    ) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0, 0, 0))
) +
plot_layout(heights = c(0.3, 1))

# Save (svg)
ggsaveDK(
  plot = p_cor_degs_alpha0.1,
  file = file.path(patch_dir, "p_degs_cor_alpha0.1.svg"),
  height = 10,
  width = 4,
  trimmed = TRUE,
  bg = "transparent"
)

#### [Interaction effect] Volcano ---------------------------------------------

# Volcano plots
p_volcano_inter.1 <- (
  (
    (
    volcanoplots_inter[[1]]
    ) + 
    (
      volcanoplots_inter[[2]] + 
      theme(axis.title.y = element_blank())
    ) +
    (
      volcanoplots_inter[[3]] +
      theme(axis.title.y = element_blank())
    )
  ) 
)

# Nr. sig. degs. barplots
p_volcano_inter.2 <- (
  (
    (
    barplots_nr_degs_inter_alpha0.05[[1]] +
      theme(plot.title = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    ) + 
    (
      barplots_nr_degs_inter_alpha0.05[[2]] + 
      theme(plot.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    ) +
    (
      barplots_nr_degs_inter_alpha0.05[[3]] +
      theme(plot.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    )
  ) /
  (
    (
      barplots_nr_degs_inter_alpha0.1[[1]] +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    ) +
     (
      barplots_nr_degs_inter_alpha0.1[[2]] + 
        theme(axis.title.y = element_blank()) +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    ) +
    (
      barplots_nr_degs_inter_alpha0.1[[3]] +
        theme(axis.title.y = element_blank()) +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    )
  ) 
)

# Combination
p_volcano_inter <- wrap_plots(
  c(p_volcano_inter.1, p_volcano_inter.2), 
  ncol = 1
  ) + 
  plot_layout(heights = c(1, 0.2))

# Save (svg)
ggsaveDK(
  plot = p_volcano_inter,
  file = file.path(patch_dir, "p_degs_volcano_inter.svg"),
  height = 18,
  width = 18,
  trimmed = TRUE,
  bg = "transparent"
)

#### [Interaction effect] FGSEA (BRCA | Basal vs non-Basal) -------------------

# Plot
gg_barplot_fgsea_interaction_BRCA <- ggplot(
  data = fgsea_res_inter[
    study == "BRCA" & 
    phenotype == "Basal vs non-Basal" &
    padj < 0.05
  ][
    order(-abs(NES))
  ][,
    head(.SD, 7),
    by = .(phenotype, tech, a_2)
  ],
  mapping = aes(
    x = NES,
    y = reorder(pathway, NES),
    fill = a_2
  )
) +
geom_col(
  # mapping = aes(
  #   fill = NES > 0
  # ),
  position = position_dodge(width = 0.8),
  width = 0.7,
  color = "black",
  linewidth = 0.1
) +
# geom_vline(
#   xintercept = 0,
#   linewidth = 0.1,
#   color = "grey80"
# ) +
scale_fill_manual(
  values = ancestry_colors
) +
facet_grid(
  rows = vars(tech),
  scales = "free_y",
  space = "free",
  labeller = labeller(
    tech = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
 fill = "Ancestry",
 x = "Normalized\nenrichment score",
 y = "Hallmark pathway"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  show_border = TRUE,
  show_grid = FALSE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save
ggsaveDK(
  plot = gg_barplot_fgsea_interaction_BRCA,
  file = file.path(patch_dir, "gg_barplot_fgsea_interaction_BRCA.svg"),
  height = 6,
  width = 7,
  trimmed = TRUE,
  bg = "transparent"
)

fgsea_res_inter[
  study == "BRCA" & 
  phenotype == "Basal vs non-Basal" &
  padj < 0.05
][
  order(-abs(NES))
][,
  head(.SD, 5),
  by = .(phenotype, tech, a_2)
][, .(pathway, tech, a_2, NES)]

#### [Baseline effect] Volcano ------------------------------------------------

# Volcano plots
p_volcano_base.1 <- (
  (
    (
    volcanoplots_base[[1]]
    ) + 
    (
      volcanoplots_base[[2]] + 
      theme(axis.title.y = element_blank())
    ) +
    (
      volcanoplots_base[[3]] +
      theme(axis.title.y = element_blank())
    )
  ) 
)

# Nr. sig. degs. barplots
p_volcano_base.2 <- (
  (
    (
    barplots_nr_degs_base_alpha0.05[[1]] +
      theme(plot.title = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    ) + 
    (
      barplots_nr_degs_base_alpha0.05[[2]] + 
      theme(plot.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    ) +
    (
      barplots_nr_degs_base_alpha0.05[[3]] +
      theme(plot.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(legend.position = "none")
    )
  ) /
  (
    (
      barplots_nr_degs_base_alpha0.1[[1]] +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    ) +
     (
      barplots_nr_degs_base_alpha0.1[[2]] + 
        theme(axis.title.y = element_blank()) +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    ) +
    (
      barplots_nr_degs_base_alpha0.1[[3]] +
        theme(axis.title.y = element_blank()) +
        theme(strip.text.x = element_blank()) +
        theme(plot.title = element_blank())
    )
  ) 
)

# Combination
p_volcano_base <- wrap_plots(
    c(p_volcano_base.1, p_volcano_base.2), 
    ncol = 1
  ) + 
  plot_layout(heights = c(1, 0.2))

# Save
ggsaveDK(
  plot = p_volcano_base,
  file = file.path(patch_dir, "p_degs_volcano_base.png"),
  height = 18,
  width = 18,
  trimmed = TRUE,
  bg = "transparent"
)

# Save (png)
ggsaveDK(
  plot = p_volcano_base,
  file = file.path(patch_dir, "Supplementary_figure_1A_2.png"),
  height = 18,
  width = 18,
  trimmed = TRUE,
  bg = "transparent"
)


##  Interaction effect groups =================================================
### fun. ----------------------------------------------------------------------

group_effects <- function(x) {
  
  sign2 <- function(x) fifelse(x > 0, 1, fifelse(x < 0, -1, 0))

  x <- copy(x)

  x[, `:=`(
    mut         = abs(T_obs_baseline_1),       # ancestry effect (analog of mutation)
    stim_EUR    = abs(T_obs_relationship_1),   # EUR cancer effect
    stim_nonEUR = abs(T_obs_relationship_2),   # non-EUR cancer effect
    int         = abs(T_obs_interaction),

    # Signs for effect direction
    s_mut       = sign2(T_obs_baseline_1),
    s_EUR       = sign2(T_obs_relationship_1),
    s_nonEUR    = sign2(T_obs_relationship_2),
    s_int       = sign2(T_obs_interaction)
  )]

  x[, effect_group := fcase(

    # (I) De novo effect
    int > mut & int > stim_EUR & int > stim_nonEUR,
      "Divergent effect",

    ## (II) Minor interaction effect
    int <= mut & int <= stim_EUR & int <= stim_nonEUR,
      "Minor interaction effect",

    # (III) Ancestry modifies cancer effect (choose EUR or non-EUR stim)
    # SIGN MATCH → enhances
    int > mut &
      (
        (int <= stim_EUR    & s_int == s_EUR) |
        (int <= stim_nonEUR & s_int == s_nonEUR)
      ),
      "non-EUR enhances\ncancer effect",

    # SIGN DIFFERENT → reverts
    int > mut &
      (
        (int <= stim_EUR    & s_int != s_EUR) |
        (int <= stim_nonEUR & s_int != s_nonEUR)
      ),
      "non-EUR reverts\ncancer effect",

    # (IV) Cancer modifies ancestry effect
    # SIGN MATCH → enhances
    (int > stim_EUR | int > stim_nonEUR) &
      int <= mut &
      s_int == s_mut,
      "Cancer enhances\nancestry effect",

    # SIGN DIFFERENT → reverts
    (int > stim_EUR | int > stim_nonEUR) &
      int <= mut &
      s_int != s_mut,
      "Cancer reverts\nancestry effect",

    default = NA_character_
  )]

  # Order categories
  order <- c(
    "Divergent effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR reverts\ncancer effect",
    "Cancer enhances\nancestry effect",
    "Cancer reverts\nancestry effect",
    "Minor interaction effect"
  )

  x[, effect_group := factor(effect_group, levels = order)]

  ## Cancer dominance
  x[, cancer_dom :=
    fifelse(
      abs(T_obs_relationship_1) > abs(T_obs_relationship_2),
      "Cancer effect\nEUR > non-EUR",
      "Cancer effect\nEUR < non-EUR"
    )
  ]

  return(x)
}

summarize_effects <- function(
  x, 
  alpha
) {

  # True tested design (preserve structure)
  full_design <- unique(
    x[
      coef_id == "interaction", 
      .(study, tech, phenotype, a_2)
    ]
  )

  # Filter significant interaction effects
  x_sig_keys <- x[
    coef_id == "interaction" & p_adj < alpha,
    .(study, tech, phenotype, a_1, a_2, g_1, g_2, feature)
  ]

  # Only use interaction phenotypes
  cols_x <- setdiff(names(x), "phenotype")

  x_sig <- x[, ..cols_x][
    x_sig_keys,
    on = .(study, tech, a_1, a_2, g_1, g_2, feature),
    nomatch = 0
  ]

  # If nothing significant → return empty but structured result
  if (nrow(x_sig) == 0) {
    return(list(
      x_sig = data.table(),
      x_sum = data.table()
    ))
  }

  # Reshape + classify
  exp_formula <- expression(
    study + tech + phenotype + a_1 + a_2 + g_1 + g_2  + feature  ~ coef_id
  )

  x_sig <- dcast(
    x_sig,
    eval(exp_formula),
    value.var = c("T_obs", "p_adj")
  )

  x_sig <- group_effects(x_sig)
  x_sig[, alpha := alpha]

  # Summarize significant ones
  x_sum <- x_sig[
    , .(N = .N),
    by = .(study, phenotype, tech, a_2, effect_group, cancer_dom)
  ]

  # Keep only effect groups that actually occur
  effect_levels <- sort(unique(na.omit(x_sig$effect_group)))

  dominance_levels <- c(
    "Cancer effect\nEUR > non-EUR",
    "Cancer effect\nEUR < non-EUR"
  )

  # Expand ONLY within tested combinations
  full_grid <- full_design[
    ,
    CJ(
      effect_group = effect_levels,
      cancer_dom = dominance_levels
    ),
    by = .(study, phenotype, tech, a_2)
  ]

  # Merge summary onto valid tested grid
  x_sum <- merge(
    full_grid,
    x_sum,
    by = c(
      "study", "phenotype", "tech", "a_2",
      "effect_group", "cancer_dom"
    ),
    all.x = TRUE
  )

  # Replace missing counts with 0
  x_sum[is.na(N), N := 0]

  # Compute percentages
  x_sum[
    , pct := if (sum(N) == 0) 0 else N / sum(N) * 100,
    by = .(study, phenotype, tech, a_2, cancer_dom)
  ]

  # Restore factor ordering
  x_sum[, `:=`(
    effect_group = factor(effect_group, levels = effect_levels),
    cancer_dom = factor(
      cancer_dom,
      levels = dominance_levels
    ),
    a_2 = factor(a_2, levels = unique(full_design$a_2)),
    alpha = alpha
  )]

  return(
    list(
      x_sig = x_sig,
      x_sum = x_sum
    )
  )
}

plot_effect_groups_pct_barplot <- function(
  data
){

  ggplot(
    data = data,
    mapping = aes(
      x = a_2, 
      y = pct, 
      fill = effect_group
    )
  ) +
  geom_col() +
  # geom_col(
  #   aes(color = pct > 0),
  #   linewidth = 0.1
  # ) +
  # scale_color_manual(
  #   values = c("FALSE" = NA, "TRUE" = "black"),
  #   guide = "none"
  # ) +
  facet_grid(
    rows = vars(cancer_dom),
    cols = vars(phenotype, tech),
    labeller = labeller(
      phenotype = vs_newline,
      tech  = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  theme_CrossAncestryGenPhen(
    rotate = 45, 
    show_borders = TRUE
  ) +
  theme(
    strip.text = element_text(margin = margin(1, 1, 1, 1)),
    plot.title = element_text(margin = margin(b = 1)),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.box.spacing = unit(0, "pt"),
    legend.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )
}

### Interaction effect groups (educational) -----------------------------------

# Cancer effect bigger in EUR
dt_effect_groups_edu_cancer_bigger_in_EUR <- data.table(
  ancestry = c(
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR"
  ),
  phenotype = c(
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer"
  ),
  effect_group = c(
    "Divergent effect", 
    "Divergent effect", 
    "Divergent effect", 
    "Divergent effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR reverts\ncancer effect", 
    "non-EUR reverts\ncancer effect",
    "non-EUR reverts\ncancer effect",
    "non-EUR reverts\ncancer effect",
    "Cancer enhances\nancestry effect", 
    "Cancer enhances\nancestry effect",
    "Cancer enhances\nancestry effect",
    "Cancer enhances\nancestry effect",
    "Cancer reverts\nancestry effect", 
    "Cancer reverts\nancestry effect", 
    "Cancer reverts\nancestry effect", 
    "Cancer reverts\nancestry effect"
  ),
  value = c(
    0, 0.1, 0, -1, 
    0, -0.7, 0, -1.4, 
    0, 1.9, 0, 0.1, 
    0, 0.1, -1, -2, 
    0, 0, 2, 0.1
  )
)
dt_effect_groups_edu_cancer_bigger_in_EUR[, cancer_dom := "Cancer effect\nEUR > non-EUR"]

# Cancer smaller in EUR
dt_effect_groups_edu_cancer_smaller_in_EUR <- data.table(
  ancestry = c(
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR",
    "EUR", "EUR", "non-EUR", "non-EUR"
  ),
  phenotype = c(
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer",
    "Control", "Cancer", "Control", "Cancer"
  ),
  effect_group = c(
    "Divergent effect", 
    "Divergent effect", 
    "Divergent effect", 
    "Divergent effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR enhances\ncancer effect",
    "non-EUR reverts\ncancer effect", 
    "non-EUR reverts\ncancer effect",
    "non-EUR reverts\ncancer effect",
    "non-EUR reverts\ncancer effect",
    "Cancer enhances\nancestry effect", 
    "Cancer enhances\nancestry effect",
    "Cancer enhances\nancestry effect",
    "Cancer enhances\nancestry effect",
    "Cancer reverts\nancestry effect",
    "Cancer reverts\nancestry effect",
    "Cancer reverts\nancestry effect",
    "Cancer reverts\nancestry effect"
  ),
  value = c(
    0, -0.1, 0, 1, 
    0, 0.8, 0, 1.6, 
    0, -1, 0, 0.1, 
    0, 0, 0.5, 2, 
    0, 0, -2, -0.2
  )
)
dt_effect_groups_edu_cancer_smaller_in_EUR[, cancer_dom := "Cancer effect\nEUR < non-EUR"]

# Combine
dt_effect_groups_edu <- rbind(
  dt_effect_groups_edu_cancer_bigger_in_EUR, 
  dt_effect_groups_edu_cancer_smaller_in_EUR
)

dt_effect_groups_edu[, `:=`(
  phenotype = factor(
    phenotype,
    levels = c("Control", "Cancer")
  ),
  cancer_dom = factor(
    cancer_dom,
    levels = c(
      "Cancer effect\nEUR > non-EUR",
      "Cancer effect\nEUR < non-EUR"
    )
  ),
  effect_group = factor(
    effect_group,
    levels = c(
      "Divergent effect",
      "non-EUR enhances\ncancer effect",
      "non-EUR reverts\ncancer effect",
      "Cancer enhances\nancestry effect",
      "Cancer reverts\nancestry effect"
    )
  ),
  color_group = ifelse(
    ancestry == "EUR",
    "EUR_black",
    as.character(effect_group)
  )
)]

# Data for end-of-line labels (Cancer only)
dt_effect_groups_edu_labels <- dt_effect_groups_edu[
  phenotype == "Cancer" &
  !effect_group %in% c(
    "Minor interaction effect",
    "Cancer enhances\nancestry effect"
  )
]

dt_effect_groups_edu_labels[, max_value := max(value), by = .(effect_group, cancer_dom)]
dt_effect_groups_edu_labels[, direction := ifelse(value == max_value, 1, -1)]
dt_effect_groups_edu_labels[, nudge_y := 0.5 * direction]

# Plot
gg_effect_groups_edu <- ggplot(
  data = dt_effect_groups_edu[
    !effect_group %in% c(
      "Minor interaction effect", 
      "Cancer enhances\nancestry effect"
      )
    ],
  mapping = aes(
    x = phenotype,
    y = value,
    color = color_group,
    group = ancestry
  )
) +
geom_hline(
  yintercept = 0,
  color = "gray80",
  linewidth = 0.5
) +
geom_line(linewidth = 1.5) +
scale_color_manual(
  values = setNames(
    effect_group_colors, 
    sub(
      "Minor interaction effect", 
      "EUR_black", 
      names(effect_group_colors)
    )
  )
) +
geom_text(
  data = dt_effect_groups_edu_labels,
  aes(
    label = ancestry,
    color = color_group,
    y = value + nudge_y
  ),
  nudge_x = -0.05,
  hjust = 0,
  size = 2,
  show.legend = FALSE
) +
scale_x_discrete(
  expand = expansion(add = c(0.2, 0.8))
) +
scale_y_continuous(
  expand = expansion(mult = c(0.05, 0.15))
) +
facet_grid(
  cols = vars(effect_group),
  rows = vars(cancer_dom)
) + 
labs(
  y = "Differntial expression (log2FC)",
  x = "Phenotype"
) +
theme_CrossAncestryGenPhen(
  show_borders = TRUE
) +
theme(
  legend.position = "none", 
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  # plot.margin = margin(0, 0, 0, 0),
  # strip.text = element_text(margin = margin(0, 0, 0, 0))
)

# Save
ggsaveDK(
  plot = gg_effect_groups_edu,
  file = file.path(fig_dir, "education_inter_groups.svg"),
  height = 5,
  width = 9.5,
  trimmed = TRUE,
  bg = "transparent"
)

### Interaction effect groups (actual nr.) ------------------------------------

# Group interaction effects
dt_effect_groups_pct_alpha0.05 <- summarize_effects(subset_interaction, alpha = 0.05)
dt_effect_groups_pct_alpha0.1  <- summarize_effects(subset_interaction, alpha = 0.1)

# Save table
write.xlsx(
  rbind(
    dt_effect_groups_pct_alpha0.05$x_sum,
    dt_effect_groups_pct_alpha0.1$x_sum
  ), 
  file = file.path(tab_dir, "dt_effect_groups_pct.xlsx")
)

write.csv(
  rbind(
    dt_effect_groups_pct_alpha0.05$x_sum,
    dt_effect_groups_pct_alpha0.1$x_sum
  ),
  file = file.path(tab_dir, "dt_effect_groups_pct.csv")
)

# Alpha = 0.05
gg_effect_groups_pct_barplot_alpha0.05 <- lapply(
  names(res_dirs),
  function(s) {

  # Plot
  p <- plot_effect_groups_pct_barplot(
    data = dt_effect_groups_pct_alpha0.05$x_sum[
      study == s
    ]
  )

  # Color + title
  p <- p +
    scale_fill_manual(
      values = effect_group_colors
    ) +
    labs(
      x = "Ancestry",
      y = "Share of genes per effect group\n (% of genes with interaction effect)",
      fill = NULL
    ) +
    ggtitle(s)

  # Outdir
  file_p <- file.path(
    fig_dir,
     paste0("barplot_pct_effect_groups_alpha0.05_", s, ".svg")
  )

  # Save
  ggsaveDK(
    plot = p,
    file = file_p,
    height = 8,
    width = 8,
    bg = "transparent"
  )
})

# Alpha = 0.1
gg_effect_groups_pct_barplot_alpha0.1 <- lapply(
  names(res_dirs),
  function(s) {

  # Plot
  p <- plot_effect_groups_pct_barplot(
    data = dt_effect_groups_pct_alpha0.1$x_sum[
      study == s
    ]
  )

  # Color + title
  p <- p +
    scale_fill_manual(
      values = effect_group_colors
    ) +
    labs(
      x = "Ancestry",
      y = "Share of genes per effect group\n(% of genes with interaction effect)",
      fill = NULL
    ) +
    ggtitle(s)

  # Outdir
  file_p <- file.path(
    fig_dir,
     paste0("barplot_pct_effect_groups_alpha0.1_", s, ".svg")
  )

  # Save
  ggsaveDK(
    plot = p,
    file = file_p,
    height = 8,
    width = 8,
    bg = "transparent"
  )
})

### Patchwork -----------------------------------------------------------------
#### Percentage effect groups (across studies) --------------------------------

# Alpha = 0.05
p_effect_groups_pct_barplot_across_alpha0.05 <- ggplot(
  data = dt_effect_groups_pct_alpha0.05$x_sum[
    , .(N = sum(N)), 
    by = .(tech, cancer_dom, a_2, effect_group)
  ][
    , pct := 100 * N / sum(N), 
    by = .(tech, cancer_dom, a_2)
  ],
  mapping = aes(
    x = a_2, 
    y = pct, 
    fill = effect_group
  )
) +
geom_col() +
scale_fill_manual(
  values = effect_group_colors
) +
facet_grid(
  rows = vars(cancer_dom),
  cols = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  fill = NULL,
  x = "Ancestry",
  y = "Share of genes per effect group\n (% of genes with interaction effect)"
) +
theme_CrossAncestryGenPhen(
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
  plot = p_effect_groups_pct_barplot_across_alpha0.05,
  file = file.path(patch_dir, "p_inter_groups_across_alpha0.05.svg"),
  height = 6,
  width = 8,
  bg = "transparent"
)

# Alpha = 0.1
p_effect_groups_pct_barplot_across_alpha0.1 <- ggplot(
  data = dt_effect_groups_pct_alpha0.1$x_sum[
    , .(N = sum(N)), 
    by = .(tech, cancer_dom, a_2, effect_group)
  ][
    , pct := 100 * N / sum(N), 
    by = .(tech, cancer_dom, a_2)
  ],
  mapping = aes(
    x = a_2, 
    y = pct, 
    fill = effect_group
  )
) +
geom_col() +
scale_fill_manual(
  values = effect_group_colors
) +
facet_grid(
  rows = vars(cancer_dom),
  cols = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  fill = NULL,
  x = "Ancestry",
  y = "Share of genes per effect group\n (% of genes with interaction effect)"
) +
theme_CrossAncestryGenPhen(
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
  plot = p_effect_groups_pct_barplot_across_alpha0.1,
  file = file.path(patch_dir, "p_inter_groups_across_alpha0.1.svg"),
  height = 6,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

#### Percentage effect groups (per study) -------------------------------------

# Alpha = 0.05
p_effect_groups_pct_barplot_per_alpha0.05 <- (
  gg_effect_groups_edu +
    theme(axis.title.x = element_blank())
  ) / 
  (
  (
    gg_effect_groups_pct_barplot_alpha0.05[[1]] +
      labs(y = "BRCA\n% of genes per effect group") +
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(strip.text.y = element_blank()) +
      theme(plot.margin = margin(r = 3)) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))

      # theme(strip.text.y = element_blank()) +
      # theme(axis.title.x = element_blank()) +
      # theme(axis.text.x = element_blank()) +
      # theme(plot.margin = margin(b = 10))
  ) /
  (
    gg_effect_groups_pct_barplot_alpha0.05[[2]] +
      labs(y = "UCEC\n% of genes per effect group") +
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(strip.text.y = element_blank()) +
      theme(plot.margin = margin(r = 3)) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))

      # theme(axis.title.x = element_blank()) +
      # theme(axis.text.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      # theme(strip.text.y = element_blank()) +
      # theme(plot.margin = margin(l = 3, b = 10))
  ) /
  (
    gg_effect_groups_pct_barplot_alpha0.05[[3]] +
      labs(y = "THCA\n% of genes per effect group") +
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(strip.text.y = element_blank()) +
      theme(plot.margin = margin(r = 3)) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))

      # theme(axis.title.x = element_blank()) +
      # theme(axis.text.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      # theme(plot.margin = margin(l = 3, b = 10))
  )
) + plot_layout(heights = c(0.3, 1)) 

# Save
ggsaveDK(
  plot = p_effect_groups_pct_barplot_per_alpha0.05,
  file = file.path(patch_dir, "p_inter_groups_per_study_alpha0.05.svg"),
  height = 18,
  width = 18,
  bg = "transparent"
)

# Alpha = 0.1
p_effect_groups_pct_barplot_per_alpha0.1 <- (
  gg_effect_groups_edu +
    theme(axis.title.x = element_blank()) 
  ) / 
  (
  (
    gg_effect_groups_pct_barplot_alpha0.1[[1]] +
      labs(y = "BRCA\n% of genes per effect group") + 
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))
      
      # theme(strip.text.y = element_blank()) +
      # theme(plot.title = element_blank()) +
      # theme(strip.text.x = element_blank())
  ) /
  (
    gg_effect_groups_pct_barplot_alpha0.1[[2]] +
      labs(y = "UCEC\n% of genes per effect group") +
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))

      # theme(strip.text.y = element_blank()) +
      # theme(plot.title = element_blank()) +
      # theme(strip.text.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      # theme(plot.margin = margin(l = 3))
  ) /
  (
    gg_effect_groups_pct_barplot_alpha0.1[[3]] +
      labs(y = "THCA\n% of genes per effect group") +
      labs(title = NULL) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      theme(strip.text = element_text(margin = margin(0, 0, 0, 0)))

      # theme(plot.title = element_blank()) +
      # theme(strip.text.x = element_blank()) +
      # theme(axis.title.y = element_blank()) +
      # theme(axis.text.y = element_blank()) +
      # theme(plot.margin = margin(l = 3))
  )
) + plot_layout(heights = c(0.3, 1)) 

# Save
ggsaveDK(
  plot = p_effect_groups_pct_barplot_per_alpha0.1,
  file = file.path(patch_dir, "p_inter_groups_per_study_alpha0.1.svg"),
  height = 18,
  width = 18,
  bg = "transparent"
)

#### Percentage effect groups (BRCA | Basal vs non-Basal) ---------------------

# Alpha 0.05
p_effect_groups_pct_barplot_BRCA_alpha0.05 <- ggplot(
  data = dt_effect_groups_pct_alpha0.05$x_sum[
    study == "BRCA" &
    phenotype == "Basal vs non-Basal"
  ],
  mapping = aes(
    x = a_2, 
    y = pct, 
    fill = effect_group
  )
) +
geom_col() +
scale_fill_manual(
  values = effect_group_colors
) +
facet_grid(
  rows = vars(cancer_dom),
  cols = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  fill = NULL,
  x = "Ancestry",
  y = "Share of genes per effect group\n (% of genes with interaction effect)"
) +
theme_CrossAncestryGenPhen(
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
  plot = p_effect_groups_pct_barplot_BRCA_alpha0.05,
  file = file.path(patch_dir, "p_inter_groups_BRCA_alpha0.05.svg"),
  height = 6,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha = 0.1
p_effect_groups_pct_barplot_BRCA_alpha0.1 <- ggplot(
  data = dt_effect_groups_pct_alpha0.1$x_sum[
    study == "BRCA" &
    phenotype == "Basal vs non-Basal"
  ],
  mapping = aes(
    x = a_2, 
    y = pct, 
    fill = effect_group
  )
) +
geom_col() +
scale_fill_manual(
  values = effect_group_colors
) +
facet_grid(
  rows = vars(cancer_dom),
  cols = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  fill = NULL,
  x = "Ancestry",
  y = "Share of genes per effect group\n (% of genes with interaction effect)"
) +
theme_CrossAncestryGenPhen(
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
  plot = p_effect_groups_pct_barplot_BRCA_alpha0.1,
  file = file.path(patch_dir, "p_inter_groups_BRCA_alpha0.1.svg"),
  height = 6,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

#### Statistics effect groups (BRCA | Basal vs non-Basal) ---------------------

# Alpha 0.1
vec_effect_groups_sig_features_alpha0.1 <- subset_interaction[
  study == "BRCA" & 
  tech == "mrna" & 
  phenotype == "Basal vs non-Basal" & 
  a_2 %in% c("AFR", "ADMIX", "EAS", "AMR") &
  coef_id == "interaction" & 
  p_adj < 0.1, 
  feature
]

# order of sig. features
vec_effect_groups_sig_features_order_alpha0.1 <- dt_effect_groups_pct_alpha0.1$x_sig[
  study == "BRCA" & 
  tech == "mrna" & 
  phenotype == "Basal vs non-Basal" & 
  a_2 %in% c("AFR", "ADMIX", "EAS", "AMR")
][
  order(effect_group, feature), 
  unique(feature)
] 

# Only the interaction effect groups
dt_effect_groups_tile_alpha0.1 <- copy(dt_effect_groups_pct_alpha0.1$x_sig)[,
  feature := factor(
    feature, 
    levels = vec_effect_groups_sig_features_order_alpha0.1
  )
][,
  facet := "Expression"
]

# Level by interaction effect group
dt_effect_groups_stats_alpha0.1 <- copy(subset_interaction)[,
  feature := factor(
    feature, 
    levels = vec_effect_groups_sig_features_order_alpha0.1
  )
]

# Tileplot
gg_effect_groups_tileplot_alpha0.1 <- ggplot(
  data = dt_effect_groups_tile_alpha0.1[
    study == "BRCA" & 
    tech == "mrna" & 
    phenotype == "Basal vs non-Basal" & 
    a_2 %in% c("AFR","ADMIX","EAS","AMR")
  ],
  aes(
    x = feature,
    y = "Effect group",              
    fill = effect_group
  )
) +
geom_tile() +
scale_fill_manual(
  values = effect_group_colors
) +
# facet_grid(
#   cols = vars(facet)
# ) +
labs(
  y = "Effect group"
) +
# scale_y_discrete(breaks = NULL) +
theme_CrossAncestryGenPhen(
  show_axis = FALSE, 
  rotate = 90
) +
theme(
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  plot.margin = margin(0, 0, 0, 0),
  legend.position = "none"
)

# Dotplot
gg_effect_groups_dotplot_alpha0.1 <- ggplot(
  data = dt_effect_groups_stats_alpha0.1[
    study == "BRCA" & 
    tech == "mrna" & 
    phenotype == "Basal vs non-Basal" & 
    a_2 %in% c("AFR","ADMIX","EAS","AMR") & 
    coef_id %in% c("relationship_1", "relationship_2", "interaction") & 
    feature %in% vec_effect_groups_sig_features_order_alpha0.1
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
  midpoint = 0,
  breaks = scales::pretty_breaks(n = 3)
) +
scale_size_continuous(
  name = expression(P.adj ~ "-log"[10] * ""),
  # breaks = scales::pretty_breaks(n = 4),
  range = c(0.5, 1.4)
) +
facet_grid(
  rows  = vars(coef_id),
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
coord_flip() +
labs(
  x = "Ancestry",
  y = "Feature"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1, 
  rotate = 90,
  show_borders = TRUE
) + 
theme(
  legend.position = "right",
  legend.title.position = "top",
  legend.margin = margin(0, 0, 0, 0),
  panel.spacing.y = unit(0.15, "lines"),
  panel.spacing.x = unit(0.15, "lines"),
  axis.text.x = element_text(size = 4),
  plot.margin = margin(0, 0, 0, 0)
)

# Patchwork
p_effect_groups_stats_dotplot_BRCA_alpha0.1 <- (
  gg_effect_groups_tileplot_alpha0.1 / 
  gg_effect_groups_dotplot_alpha0.1
) + 
plot_layout(heights = c(0.05, 1))

# Save
ggsaveDK(
  plot = p_effect_groups_stats_dotplot_BRCA_alpha0.1,
  file = file.path(patch_dir, "p_effect_groups_stats_dotplot_BRCA_alpha0.1.svg"),
  height = 6,
  width = 11,
  trimmed = FALSE,
  bg = "transparent"
)

#### Correlation effect size (BRCA | Basal vs non-Basal) ----------------------
dt_effect_sizes_alpha0.1 <- dcast(
  subset_interaction[
    study == "BRCA" & 
    tech == "mrna" & 
    phenotype == "Basal vs non-Basal" & 
    a_2 %in% c("AFR", "ADMIX", "EAS", "AMR") &
    coef_id == "interaction" &
    feature %in% vec_effect_groups_sig_features_alpha0.1
  ],
  feature ~ a_2, 
  value.var = "T_obs"
)

# Add effect group
dt_effect_sizes_alpha0.1 <- merge(
  dt_effect_groups_pct_alpha0.1$x_sig[
    study == "BRCA" & 
    tech == "mrna" & 
    phenotype == "Basal vs non-Basal" & 
    a_2 %in% c("AFR", "ADMIX", "EAS", "AMR"),
    .(feature, effect_group)
  ],
  dt_effect_sizes_alpha0.1,
  by = "feature"
)

# Correlation
split_list <- split(
  dt_effect_sizes_alpha0.1, 
  dt_effect_sizes_alpha0.1$effect_group
)
split_list <- split_list[sapply(split_list, nrow) > 0]

heatmaps <- lapply(names(split_list), function(grp) {
  
  sub <- split_list[[grp]]
  
  mat <- as.matrix(sub[, c("ADMIX", "AFR", "AMR", "EAS"), with = FALSE])
  rownames(mat) <- sub$feature
  
  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  diag(cor_mat) <- NA
  
  ht <- Heatmap(
    cor_mat,
    name = "Pearson",
    column_title = grp,
    col = colorRamp2(
      c(-1, 0, 1), 
      c("blue", "white", "red")
    ),
    rect_gp = gpar(col = "white"),
    na_col = "gray90",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    # Dims
    width = unit(1, "cm"),
    height = unit(1, "cm"),
    # Legend params
    heatmap_legend_param = list(
      at = c(-1, 0 ,1),
      labels = c("-1", "0", "1"),
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      legend_height = unit(1, "cm"),
      grid_width = unit(0.15, "cm")
    ),
    # Dendogram
    row_dend_width = unit(0.3, "cm"),
    column_dend_height = unit(0.3, "cm"),
    # Font
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6)
  )

  # Clean group
  grp_clean <- gsub("[ \n]", "_", grp)

  # File name
  file <- file.path(
    fig_dir,
    paste0(
      "cmplx_heatmap_degs_corr_",
      grp_clean,
      "_BRCA_alpha0.1.pdf"
    )
  )

  # Save
  pdf(
    file,
    width = 1.5, 
    height = 1.5
  )
  draw(ht)
  dev.off()
})

#### Single features (BRCA | Basal vs non-Basal) ------------------------------
data <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")
data <- data[["mrna"]]

# Filter sample type
data <- filter_sample_type(
  X = data$matr, 
  M = data$meta, 
  s_col ="SAMPLE_TYPE", 
  s_levels = "Primary"
)

# Extract components
matr <- as.data.table(data$matr, keep.rownames = "SAMPLE_ID")
meta <- as.data.table(data$meta)

# Features
vec_features_alpha0.1 <- c(
  "TPD52L1",   # De novo | EAS
  "ZFP69",     # De novo | EAS
  "FBXW5",     # De novo | EAS
  "PRR36",     # De novo | EAS
  "PTPN9",     # De novo | EAS
  "MEST",      # De novo | AFR
  "ACSL6",     # non-EUR enhances | EAS
  "TNFAIP1",   # non-EUR enhances | EAS
  "LRRC34",    # non-EUR enhances | EAS
  "TNF",       # non-EUR enhances | EAS
  "BBS2",      # non-EUR enhances | AFR
  "CDKN2B",    # non-EUR enhances | AFR
  "ACSS1",     # non-EUR enhances | AFR
  "RNASEH2B",  # Cancer reverts | EAS
  "CDKAL1",    # Cancer reverts | EAS
  "MCF2L-AS1", # Cancer reverts | EAS
  "RHCE",      # Cancer reverts | EAS
  "PDXK",      # non-EUR reverts | EAS
  "TMEM30B",   # non-EUR reverts | EAS
  "EIF4E3"     # non-EUR reverts | AFR
)

# Subset matrix to features
matr_feature_alpha0.1 <- matr[, c("SAMPLE_ID", vec_features_alpha0.1), with = FALSE]
setkey(meta, SAMPLE_ID); setkey(matr_feature_alpha0.1, SAMPLE_ID)
matr_feature_alpha0.1 <- meta[matr_feature_alpha0.1]

# Merge with stats results
dt_expr_alpha0.1 <- melt(
  matr_feature_alpha0.1,
  id.vars = c(
    "SAMPLE_ID", 
    "POOLED_SUBTYPE",
    "POOLED_GENETIC_ANCESTRY"
  ),  
  measure.vars = vec_features_alpha0.1,
  variable.name = "feature",
  value.name = "expression"
)[
  , `:=`(
    a_2 = POOLED_GENETIC_ANCESTRY
  )
][
  , `:=` (
    study = "BRCA",
    tech = "mrna",
    phenotype = "Basal vs non-Basal"
  )
][
  , c("POOLED_GENETIC_ANCESTRY") := NULL
][
  , zscore := scale(expression), 
  by = feature
]

# Stats for features
dt_expr_stats_alpha0.1 <- merge(
  dt_expr_alpha0.1,
  subset_interaction[
    coef_id == "interaction" &
    study == "BRCA" &
    tech == "mrna" &
    phenotype == "Basal vs non-Basal" &
    feature %in% vec_features_alpha0.1,
  .(study, tech, phenotype, a_2, feature, p_adj)
  ],
  by = c("study", "tech", "phenotype", "a_2", "feature"),
  all.x = TRUE
)

# Plot
dt_plot_alpha0.1 <- dt_expr_stats_alpha0.1[
  POOLED_SUBTYPE %in% c("Basal", "non_Basal") &
  !is.na(a_2)
][, 
  if (uniqueN(POOLED_SUBTYPE) == 2) .SD, by = a_2
][, 
  a_2 := factor(
    a_2,
    levels = c("EUR", setdiff(unique(a_2), "EUR"))
  )
][,
  feature := factor(
    feature,
    levels = vec_features_alpha0.1
  )
]

# Add interaction effect groups
dt_sig_alpha0.1 <- dt_effect_groups_pct_alpha0.1$x_sig[
  dt_expr_stats_alpha0.1[
    POOLED_SUBTYPE %in% c("Basal", "non_Basal") & !is.na(a_2)
  ][
    , if (uniqueN(POOLED_SUBTYPE) == 2) .SD, by = a_2
  ][
    p_adj < 0.1,
    .(study, tech, phenotype, feature, a_2)
  ],
  on = .(study, tech, phenotype, feature, a_2),
  nomatch = NA 
][, 
  a_2 := factor(
    a_2, 
    levels = c("EUR", setdiff(unique(dt_expr_stats_alpha0.1$a_2), "EUR"))
  )
][, 
  x := as.numeric(a_2)
][, 
  .(study, tech, phenotype, feature, a_2, effect_group, x)
]

# Zoomed facets
coef_val <- 0.5

# Compute whisker-filtered min/max per group
dt_whiskers <- dt_plot_alpha0.1[
  , {
    qs <- quantile(zscore, c(0.25, 0.75), na.rm = TRUE)
    q1 <- qs[1]
    q3 <- qs[2]
    iqr <- q3 - q1

    lower_bound <- q1 - coef_val * iqr
    upper_bound <- q3 + coef_val * iqr

    # keep only values within whiskers
    z_in <- zscore[zscore >= lower_bound & zscore <= upper_bound]

    .(
      wmin = min(z_in, na.rm = TRUE),
      wmax = max(z_in, na.rm = TRUE)
    )
  },
  by = .(feature, a_2, POOLED_SUBTYPE)
]

ylims <- dt_whiskers[
  , .(
    ymin = min(wmin, na.rm = TRUE),
    ymax = max(wmax, na.rm = TRUE)
  ),
  by = feature
]

ylims[, {
  range <- ymax - ymin
  
  pad_bottom <- 0.0 * range   
  pad_top    <- 0.0 * range  
  
  .(
    lims = list(c(
      ymin - pad_bottom,
      ymax + pad_top
    ))
  )
}, by = feature] -> ylims
ylims[, feature := factor(feature, levels = vec_features_alpha0.1)]
setorder(ylims, feature)

yscales <- lapply(ylims$lims, function(lim) {
  scale_y_continuous(
    limits = lim,
    breaks = pretty(lim, n = 3)
  )
})

dt_sig_alpha0.1 <- merge(
  dt_sig_alpha0.1,
  ylims,
  by = "feature",
  all.x = TRUE
)[
  , `:=`(
    ymin = lims[[1]][1],
    ymax = lims[[1]][2]
  ),
  by = feature
]

# Plot
gg_effect_groups_zscore_boxplot_BRCA_alpha0.1 <- ggplot(
  data = dt_plot_alpha0.1[, 
    feature := factor(
      feature, 
      levels = vec_features_alpha0.1
    )
  ],
  mapping = aes(
    x = a_2,
    y = zscore,
  )
) +
geom_rect(
  data = unique(
    dt_sig_alpha0.1[, .(feature, x, effect_group)]
  )[, feature := factor(feature, levels = vec_features_alpha0.1)],
  mapping = aes(
    xmin = x - 0.5,
    xmax = x + 0.5,
    fill = effect_group
  ),
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.3,
  color = NA,
  show.legend = FALSE,
  inherit.aes = FALSE
) +
scale_fill_manual(
  name = "Effect group",
  values = effect_group_colors,
) +
ggnewscale::new_scale_fill() +
geom_boxplot(
  mapping = aes(
    fill = POOLED_SUBTYPE
  ),
  outlier.shape = NA,
  coef = coef_val,
  linewidth = 0.1
) +
facet_wrap(
  ~ feature,
  scales = "free_y",
  nrow = 4
  #space = "free"
) +
facetted_pos_scales(
  y = yscales
) +
labs(
  x = "Ancestry",
  y = "Z-score",
  fill = "Phenotype"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_border = TRUE
) +
theme(
  panel.spacing.x = unit(0.05, "lines"),
  panel.spacing.y = unit(0.05, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  # plot.margin = margin(0, 0, 0, 0),
  strip.text = element_text(margin = margin(0, 0, 0, 0))
)

# Save
ggsaveDK(
  plot = gg_effect_groups_zscore_boxplot_BRCA_alpha0.1,
  file = file.path(fig_dir, "gg_effect_groups_zscore_boxplot_BRCA_alpha0.1.svg"),
  height = 6,
  width = 11,
  trimmed = TRUE,
  bg = "transparent"
)

##  Overrepresentation analysis ===============================================
### fun. ----------------------------------------------------------------------

ORA_one <- function(
  features, 
  metadata,
  term2gene,
  universe
  ) {

  res <- enricher(
    gene = features,
    TERM2GENE = term2gene,
    universe = universe,
    pAdjustMethod = "BH"
  )

  # If enrichment returned nothing
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  # Convert S4 enrichResult → data.table
  dt <- as.data.table(res@result)

  # Add metadata columns
  dt[, study        := metadata$study]
  dt[, phenotype    := metadata$phenotype]
  dt[, tech         := metadata$tech]
  dt[, a_2          := metadata$a_2]
  dt[, effect_group := metadata$effect_group]
  dt[, cancer_dom   := metadata$cancer_dom]

  return(dt)
}

plot_ora_dotplot <- function(
  data
){

  ggplot(
    data = data,
    mapping = aes(
      x = a_2, 
      y = reorder(Description, FoldEnrichment),,
      size  = -log10(p.adjust), 
      color = log10(FoldEnrichment)
    )
  ) +
  geom_point(
  ) +
  scale_color_continuous(
    name   = "Fold\nenrichment",
    low    = "blue",
    high   = "red",
    label  =  math_format(10^.x),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_size_continuous(
    name = expression(P.adj ~ "-log"[10] * ""),
    range = c(0.8, 2)
  ) +
  facet_grid(
    cols = vars(phenotype, tech),
    rows = vars(cancer_dom, effect_group),
    scales = "free_y",
    labeller = labeller(
      phenotype = vs_newline,
      tech = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    x = "Ancestry",
    y = "Hallmark pathway"
  ) +
  theme_CrossAncestryGenPhen(
    base_size = 6,
    legend_key = 1, 
    rotate = 45, 
    show_border = TRUE, 
    show_grid = TRUE
  ) +
  theme(
    plot.margin = margin(r = 1, 0, 0, 0),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0)),
    strip.text.y = element_blank(),
    legend.position  = "right",
    legend.direction = "vertical"
  )
}

plot_effect_groups_education <- function(
  data,
  labels
){
  ggplot(
    data = data,
    mapping = aes(
      x = phenotype,
      y = value,
      color = color_group,
      group = ancestry
    )
  ) +
  geom_hline(
    yintercept = 0,
    color = "gray80",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(
    values = setNames(
      effect_group_colors, 
      sub(
        "Minor interaction effect", 
        "EUR_black", 
        names(effect_group_colors)
      )
    )
  ) +
  geom_text(
    data = labels,
    aes(
      label = ancestry,
      color = color_group,
      y = value + nudge_y
    ),
    nudge_x = -0.05,
    hjust = 0,
    size = 2,
    show.legend = FALSE
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.2, 1.8))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  facet_grid(
    rows = vars(cancer_dom, effect_group),
    drop = TRUE
  ) + 
  labs(
    y = "Differntial expression (log2FC)",
    x = "Phenotype"
  ) +
  theme_CrossAncestryGenPhen(
    show_borders = FALSE,
    show_axis = FALSE
  ) +
  theme(
    legend.position = "none", 
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
}

### ORA analysis --------------------------------------------------------------

# Alpha = 0.05
dt_ora_alpha0.05 <- copy(dt_effect_groups_pct_alpha0.05$x_sig)[
  tech == "meth", 
  feature := probe2gene[.SD, on = "feature", gene]
][
  effect_group != "Minor interaction effect",
  feature,
  by = .(study, phenotype, tech, a_2, effect_group, cancer_dom)
][
  feature %in% term2gene$gene_symbol
][,
  {
    # ORA run
    meta <- .SD[1, .(study, phenotype, tech, a_2, effect_group, cancer_dom)]
    res  <- ORA_one(
      features = feature, 
      metadata = meta,
      term2gene = term2gene,
      universe = unique(term2gene$gene_symbol)
      )
    .(res = list(res))
  },
  by = .(study, phenotype, tech, a_2, effect_group)
]

# Add alpha
dt_ora_alpha0.05 <- rbindlist(
  dt_ora_alpha0.05$res, 
  fill = TRUE
)[, alpha := 0.05]

# Plot
gg_ora_dotplot_alpha0.05 <- lapply(
  names(res_dirs),
  function(s) {

  # Subset ORA results
  study_data <- dt_ora_alpha0.05[
    study == s &
    p.adjust < 0.1
  ]

  # Present effect groups 
  present_effect_groups <- unique(
    study_data[, .(effect_group, cancer_dom)]
  )

  # Subset educational groups
  subset_effect_groups_edu <- dt_effect_groups_edu[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Subset ancestry labels
  subset_ancestry_labels <- dt_effect_groups_edu_labels[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Effect plot
  p_edu <- plot_effect_groups_education(
    data = subset_effect_groups_edu,
    labels = subset_ancestry_labels
  )

  # ORA plot
  p_ora <- plot_ora_dotplot(
    data = study_data
  )

  # Patchwork
  p <- p_ora + p_edu + 
    plot_layout(
      guides = "collect",
      widths = c(1, 0.3)
    ) &
    theme(
      #axis.title.y = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(5, "pt")
    )

  # Outdir
  file_p <- file.path(
    fig_dir,
    paste0("dotplot_ora_", s,"_alpha0.05.svg")
  )

  # Save
  ggsaveDK(
    plot = p,
    file = file_p,
    height = 8,
    width = 12,
    bg = "transparent"
  )
})

# Alpha = 0.1
dt_ora_alpha0.1 <- copy(dt_effect_groups_pct_alpha0.1$x_sig)[
  tech == "meth", 
  feature := probe2gene[.SD, on = "feature", gene]
][
  effect_group != "Minor interaction effect",
  feature,
  by = .(study, phenotype, tech, a_2, effect_group, cancer_dom)
][
  feature %in% term2gene$gene_symbol
][,
  {
    # ORA run
    meta <- .SD[1, .(study, phenotype, tech, a_2, effect_group, cancer_dom)]
    res  <- ORA_one(
      features = feature, 
      metadata = meta,
      term2gene = term2gene,
      universe = unique(term2gene$gene_symbol)
      )
    .(res = list(res))
  },
  by = .(study, phenotype, tech, a_2, effect_group)
]

# Add alpha
dt_ora_alpha0.1 <- rbindlist(
  dt_ora_alpha0.1$res, 
  fill = TRUE
)[, alpha := 0.1]

dt_ora_alpha0.1[
  study == "BRCA" &
  tech == "mrna" &
  phenotype == "Basal vs non-Basal" &
  p.adjust < 0.1,
  .(Description, a_2, GeneRatio, geneID, p.adjust, effect_group)
]

# Plots
gg_ora_dotplot_alpha0.1 <- lapply(
  names(res_dirs),
  function(s) {

  # Subset ORA results
  study_data <- dt_ora_alpha0.1[
    study == s &
    p.adjust < 0.1
  ]

  # Present effect groups 
  present_effect_groups <- unique(
    study_data[, .(effect_group, cancer_dom)]
  )

  # Subset educational groups
  subset_effect_groups_edu <- dt_effect_groups_edu[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Subset ancestry labels
  subset_ancestry_labels <- dt_effect_groups_edu_labels[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Effect plot
  p_edu <- plot_effect_groups_education(
    data = subset_effect_groups_edu,
    labels = subset_ancestry_labels
  )

  # ORA plot
  p_ora <- plot_ora_dotplot(
    data = study_data
  )

  # Patchwork
  p <- p_ora + p_edu + 
    plot_layout(
      guides = "collect",
      widths = c(1, 0.3)
    ) &
    theme(
      #axis.title.y = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(5, "pt")
    )


  # Outdir
  file_p <- file.path(
    fig_dir,
    paste0("dotplot_ora_", s,"_alpha0.1.svg")
  )

  # Save
  ggsaveDK(
    plot = p,
    file = file_p,
    height = 8,
    width = 12,
    bg = "transparent"
  )
})

# Tables
write.xlsx(
  rbind(
    dt_ora_alpha0.05, 
    dt_ora_alpha0.1
  ), 
  file = file.path(tab_dir, "ORA_res.xlsx")
)

write.csv(
  rbind(
    dt_ora_alpha0.05, 
    dt_ora_alpha0.1
  ),
  file = file.path(tab_dir, "ORA_res.csv")
)

### Patchwork -----------------------------------------------------------------
#### ORA analysis (BRCA | Basal vs non-Basal) ---------------------------------

# Alpha 0.05
p_effect_groups_ora_barplot_BRCA_alpha0.05 <- lapply(
  c("BRCA"),
  function(s) {

  # Subset ORA results
  study_data <- dt_ora_alpha0.05[
    study == s &
    phenotype == "Basal vs non-Basal" &
    tech == "mrna" &
    p.adjust < 0.1
  ]


  # Present effect groups 
  present_effect_groups <- unique(
    study_data[, .(effect_group, cancer_dom)]
  )

  # Subset educational groups
  subset_effect_groups_edu <- dt_effect_groups_edu[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Subset ancestry labels
  subset_ancestry_labels <- dt_effect_groups_edu_labels[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Effect plot
  p_edu <- ggplot(
    data = subset_effect_groups_edu,
    mapping = aes(
      x = phenotype,
      y = value,
      color = color_group,
      group = ancestry
    )
  ) +
  geom_hline(
    yintercept = 0,
    color = "gray80",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(
      effect_group_colors, 
      sub(
        "Minor interaction effect", 
        "EUR_black", 
        names(effect_group_colors)
      )
    )
  ) +
  scale_y_discrete(
    expand = expansion(add = c(2, 2))
  ) +
  facet_grid(
    rows = vars(cancer_dom, effect_group),
    drop = TRUE
  ) + 
  labs(
    y = "Differntial expression (log2FC)",
    x = "Phenotype"
  ) +
  theme_CrossAncestryGenPhen(
    show_borders = FALSE,
    show_axis = FALSE
  ) +
  theme(
    legend.position = "none", 
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

  # ORA plot
  p_ora <- ggplot(
    data = study_data,
    mapping = aes(
      x = FoldEnrichment, 
      y = reorder(Description, FoldEnrichment),
      fill = a_2,
    )
  ) +
  geom_col(
     position = position_dodge(width = 0.8),
     width = 0.7,
     color = "black",
     linewidth = 0.1
  ) +
  scale_fill_manual(
    values = ancestry_colors
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 3)
  ) +
  facet_grid(
    # cols = vars(tech),
    rows = vars(cancer_dom, effect_group),
    scales = "free_y",
    space = "free",
    labeller = labeller(
      cancer_dom = function(x) "",
      effet_group = function(x) "",
      tech = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    fill = "Ancestry",
    x = "Fold enrichment",
    y = "Hallmark pathway"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1, 
    show_border = TRUE, 
    show_grid = TRUE
  ) +
  theme(
    plot.margin = margin(r = 1, 0, 0, 0),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.text.y = element_blank(),
    strip.placement = "inside",
    legend.position  = "right",
    axis.text.y = element_text(size = 6)
  )

  # Patchwork
  p <- p_ora + p_edu + 
    plot_layout(
      guides = "collect",
      widths = c(1, 0.5)
    ) &
    theme(
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(5, "pt")
    )
  
  # Return
  return(p)
})

# Save
ggsaveDK(
  plot = p_effect_groups_ora_barplot_BRCA_alpha0.05,
  file = file.path(patch_dir, "p_effect_groups_ora_barplot_BRCA_alpha0.05.svg"),
  height = 6,
  width = 7,
  trimmed = FALSE,
  bg = "transparent"
)

# Alpha 0.1
p_effect_groups_ora_barplot_BRCA_alpha0.1 <- lapply(
  c("BRCA"),
  function(s) {

  # Subset ORA results
  study_data <- dt_ora_alpha0.1[
    study == s &
    phenotype == "Basal vs non-Basal" &
    tech == "mrna" &
    p.adjust < 0.1
  ]


  # Present effect groups 
  present_effect_groups <- unique(
    study_data[, .(effect_group, cancer_dom)]
  )

  # Subset educational groups
  subset_effect_groups_edu <- dt_effect_groups_edu[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Subset ancestry labels
  subset_ancestry_labels <- dt_effect_groups_edu_labels[
    present_effect_groups,
    on = .(effect_group, cancer_dom),
    nomatch = 0
  ]

  # Effect plot
  p_edu <- ggplot(
    data = subset_effect_groups_edu,
    mapping = aes(
      x = phenotype,
      y = value,
      color = color_group,
      group = ancestry
    )
  ) +
  geom_hline(
    yintercept = 0,
    color = "gray80",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(
      effect_group_colors, 
      sub(
        "Minor interaction effect", 
        "EUR_black", 
        names(effect_group_colors)
      )
    )
  ) +
  scale_y_discrete(
    expand = expansion(add = c(2, 2))
  ) +
  facet_grid(
    rows = vars(cancer_dom, effect_group),
    drop = TRUE
  ) + 
  labs(
    y = "Differntial expression (log2FC)",
    x = "Phenotype"
  ) +
  theme_CrossAncestryGenPhen(
    show_borders = FALSE,
    show_axis = FALSE
  ) +
  theme(
    legend.position = "none", 
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

  # ORA plot
  p_ora <- ggplot(
    data = study_data,
    mapping = aes(
      x = FoldEnrichment, 
      y = reorder(Description, FoldEnrichment),
      fill = a_2,
    )
  ) +
  geom_col(
     position = position_dodge(width = 0.8),
     width = 0.7,
     color = "black",
     linewidth = 0.1
  ) +
  scale_fill_manual(
    values = ancestry_colors
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 3)
  ) +
  facet_grid(
    # cols = vars(tech),
    rows = vars(cancer_dom, effect_group),
    scales = "free_y",
    space = "free",
    labeller = labeller(
      cancer_dom = function(x) "",
      effet_group = function(x) "",
      tech = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    fill = "Ancestry",
    x = "Fold enrichment",
    y = "Hallmark pathway"
  ) +
  theme_CrossAncestryGenPhen(
    legend_key = 1, 
    show_border = TRUE, 
    show_grid = FALSE
  ) +
  theme(
    plot.margin = margin(r = 1, 0, 0, 0),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.text.y = element_blank(),
    legend.position  = "right",
    legend.title.position = "top"
  )

  # Patchwork
  p <- p_ora + p_edu + 
    plot_layout(
      guides = "collect",
      widths = c(1, 1)
    ) &
    theme(
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(5, "pt")
    )
  
  # Return
  return(p)
})

# Save
ggsaveDK(
  plot = p_effect_groups_ora_barplot_BRCA_alpha0.1,
  file = file.path(patch_dir, "p_effect_groups_ora_barplot_BRCA_alpha0.1.svg"),
  height = 7,
  width = 7,
  trimmed = FALSE,
  bg = "transparent"
)

## Clusters ===================================================================
### TSNE ----------------------------------------------------------------------
# Seed
set.seed(1)

# Load data 
study_file <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")
techs <- c("meth", "mrna")
attrs <- c(
  "POOLED_GENETIC_ANCESTRY", 
  "SUBTYPE"
  # "SAMPLE_TYPE", 
  # "SEX"
)

attr_labels <- c(
  SAMPLE_TYPE = "Sample type",
  SUBTYPE = "Molecular subtype",
  POOLED_GENETIC_ANCESTRY = "Ancestry",
  SEX = "Sex"
)

# Custom palettes for each attribute
attr_colors <- list(
  SAMPLE_TYPE = sample_type_colors,
  SUBTYPE = subtype_colors,
  POOLED_GENETIC_ANCESTRY = ancestry_colors,
  SEX = sex_colors
)

# Prepare data
tsne_data <- lapply(techs, function(tech) {

  mat  <- study_file[[tech]]$matr
  meta <- study_file[[tech]]$meta
  meta <- meta[match(rownames(mat), meta$SAMPLE_ID), ]

  # Filter
  if (tech == "mrna") {
    dge  <- edgeR::DGEList(counts = t(mat))
    keep <- edgeR::filterByExpr(dge, group = factor(rep("All", ncol(dge))))
    mat  <- t(dge$counts[keep,,drop=FALSE])
    mat  <- t(edgeR::cpm(t(mat), log=TRUE))
  }

  if (tech == "meth") {
    v    <- apply(mat,2,var,na.rm=TRUE)
    v[is.na(v)] <- 0
    keep <- v > quantile(v,0.20)
    mat  <- mat[,keep,drop=FALSE]
    mat  <- beta_to_mval(mat)
  }
  
  # Cluster
  tsne <- Rtsne(
    mat, 
    perplexity = 50, 
    dims = 2
  )

  # Return
  list(
    coords = data.frame(
      TSNE1 = tsne$Y[,1], 
      TSNE2 = tsne$Y[,2]
    ),
    meta = meta,
    mat = mat
  )
})
names(tsne_data) <- techs

# Plot
tsne_plots <- lapply(techs, function(tech) {
  coords <- tsne_data[[tech]]$coords
  meta   <- tsne_data[[tech]]$meta

  panels  <- list()
  legends <- list()
  for (a in attrs) {

    keep <- !is.na(meta[[a]])
    keep <- keep & meta$SAMPLE_TYPE == "Primary"
    df   <- cbind(coords[keep, , drop = FALSE], color = factor(meta[[a]][keep]))

    # Plot
    p_full <- ggplot(
      data = df, 
      mapping = aes(
        x = TSNE1, 
        y = TSNE2, 
        color = color
        )
      ) +
      geom_point(
        size = 0.8
      ) +
      scale_color_manual(
        values = attr_colors[[a]], 
        drop = FALSE
      ) +
      labs(
        title = attr_labels[[a]], 
        color = attr_labels[[a]], 
        x = "t-SNE 1", 
        y = "t-SNE 2"
      ) +
      coord_fixed(ratio = 1) +
      theme_CrossAncestryGenPhen(
        legend_key = 1
      ) +
      theme(
        # Mimick void theme
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        legend.title = element_blank()
      )

    panels[[a]]  <- p_full
  }

  # Save
  ggsaveDK(
    file = file.path(fig_dir, paste0("pointplot_cluster_", tech, "_BRCA.svg")),
    plot = wrap_plots(panels, ncol = 2),
    width = 8, 
    height = 8,
    bg = "transparent",
    trimmed = FALSE
  )

  # Return
  return(NULL)
})


## --- Single feature expression ---
data <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")
tech <- "mrna"
data <- data[[tech]]
# Filter sample type
data <- filter_sample_type(X = data$matr, M = data$meta, s_col ="SAMPLE_TYPE", s_levels = "Primary")
# Filter by ancestry and phenotype
g_col <- "POOLED_SUBTYPE"; g1 <- "Basal"; g2 <- "non_Basal"
a_col <- "POOLED_GENETIC_ANCESTRY"; a1 <- "EUR"; a2 <- "EAS"
study <- filter_phenotype_ancestry(X = data$matr, M = data$meta, g_col = g_col, a_col = a_col, g_levels = c(g1, g2), a_levels = c(a1, a2), plot = FALSE)

# Select featurs
DT <- sig_interaction_grp[
  tech == "mrna" &
  phenotype == "Basal vs non-Basal" &
  a_2 == "EAS" &
  effect_group != "Minor interaction effect", feature
]

n_groups <- nrow(top1)
remaining_slots <- max(0, 9 - n_groups)
DT_remaining <- DT[!feature %in% top1$feature]
top_extra <- DT_remaining[1:remaining_slots]

feat_sel <- rbind(top1, top_extra)[1:9]
feat_p <- feat_sel$feature

feat_p <- c(
  "ZFP69",    # De novo
  "PTPN9",    # De novo
  "TNFAIP1",  # non-EUR enhance
  "LRRC34",   # non-EUR enhances
  "TNF",      # non-EUR enhances
  "PDXK",     # non-EUR reverts
  "TMEM30B",  # non-EUR reverts
  "RNASEH2B", # Cancer reverts
  "MCF2L-AS1" # Cancer reverts
)
# Boxplots
p <- plot_feature(
  X = study$X$matr,
  Y = study$Y$matr, 
  MX = study$X$meta, 
  MY = study$Y$meta, 
  x_var = a_col,
  fill_var = g_col,
  show_outliers = FALSE,
  features = feat_p,
  x_label = "Ancestry",
  y_label = "Z-score"
)


df <- p$data
# Compute per-feature zoom limits
ylims <- lapply(split(df$value, df$feature), function(v) {
  quantile(v, c(0.005, 0.995), na.rm = TRUE)
})
# Build facet-specific y scales with pretty breaks
yscales <- lapply(ylims, function(lim) {
  scale_y_continuous(
    limits = lim,
    breaks = pretty(lim, n = 3)   
  )
})
# Apply to plot
p2 <- p + facetted_pos_scales(y = yscales)


ggsaveDK(
  plot = p2,
  file = file.path(out_dir, "Interaction_effect_zscore.svg"),
  height = 6,
  width = 6,
  trimmed = FALSE,
  guides = FALSE,
  bg = "transparent"
)


## --- Enrichment of interaction effect groups ---
# Map probes to gene
probe2gene <- as.data.table(fread("download/GDCportal/TCGA_BRCA/probe2gene.csv"))
setnames(probe2gene, "probe", "feature")
sig_interaction_grp[tech == "meth", feature := probe2gene[.SD, on = "feature", gene]]

# Enrichment for one gene vector returning a CLEAN data.table
enrich_one <- function(genes, metadata) {

  res <- enricher(
    gene          = genes,
    TERM2GENE     = term2gene,
    universe      = unique(term2gene$gene_symbol),
    pAdjustMethod = "BH"
  )

  # if enrichment returned nothing
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  # convert S4 enrichResult → data.table
  dt <- as.data.table(res@result)

  # add metadata columns
  dt[, study        := metadata$study]
  dt[, phenotype    := metadata$phenotype]
  dt[, tech         := metadata$tech]
  dt[, a_2          := metadata$a_2]
  dt[, effect_group := metadata$effect_group]

  return(dt)
}


# Prepare database
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
TERM2GENE <- as.data.table(hallmark)[, .(gs_name, gene_symbol)]
TERM2GENE[, gs_name := {
    x <- sub("^HALLMARK_", "", gs_name)   
    x <- gsub("_", " ", x)         
    x <- tools::toTitleCase(tolower(x))
    x
}]


# Prepare enrich data
gene_sets <- sig_interaction_grp[effect_group != "Minor interaction effect", feature, by = .(study, phenotype, tech, a_2, effect_group)]
gene_sets <- gene_sets[feature %in% unique(TERM2GENE$gene_symbol)]
# Enrich
tmp <- gene_sets[,
  {
    # extract metadata
    meta <- .SD[1, .(study, phenotype, tech, a_2, effect_group)]
    res <- enrich_one(feature, meta)
    .(res = list(res))
  },
  by = .(study, phenotype, tech, a_2, effect_group)
]
clean_enriched <- rbindlist(tmp$res, fill = TRUE)

# Save enrichment analysis
enrich_dir <- file.path(out_dir, "enrichment_analysis")
if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)
saveRDS(clean_enriched, file.path(enrich_dir, "enrichment_results.rds"))

# Plotting universe
p_data <- clean_enriched[p.adjust < 0.05, ]
color_limits <- range(log10(p_data$FoldEnrichment), na.rm = TRUE)
size_limits  <- range(-log10(p_data$p.adjust), na.rm = TRUE)

# Per study enrichment
enrich_list <- list()
for (s in c("BRCA", "THCA", "UCEC")){

  ## Labeller
  vs_newline <- function(x) {
    gsub("\\s*vs\\s*", " vs\n", x)
  }

  ## ggplot
  p <- ggplot(
    data = p_data[study == s, ],
    mapping = aes(
        x = a_2, 
        y = Description, 
        size  = -log10(p.adjust), 
        color = log10(FoldEnrichment)
      )
    ) +
    geom_point(
    ) +
    scale_color_continuous(
      name   = "Fold enrichment  ",
      low    = "blue",
      high   = "red",
      limits = color_limits,
      label  =  math_format(10^.x),
      breaks = scales::pretty_breaks(n = 3)
    ) +
    scale_size_continuous(
      name   = expression(P.adj ~ "(-log"[10] * ")  "),
      limits = size_limits,
      range  = c(0.8, 2)
    ) +
    facet_grid(
      cols = vars(phenotype, tech),
      rows = vars(effect_group),
      scales = "free_y",
      labeller = labeller(
        phenotype = vs_newline,
        tech = c(meth = "Methylation", mrna = "Expression")
      )
    ) +
    labs(
      x = "Ancestry",
      y = "Hallmark Pathway"
    ) +
    theme_CrossAncestryGenPhen(
      base_size = 5,
      legend_key = 1, 
      rotate = 45, 
      show_border = TRUE, 
      show_grid = TRUE
    ) +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      panel.spacing.y = unit(0.15, "lines"),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "bottom", 
      legend.title.position = "left",
      legend.box            = "vertical",
      legend.box.spacing    = unit(0, "pt"),
      legend.spacing        = unit(0, "pt"),
      legend.spacing.y      = unit(0, "pt"),
      legend.spacing.x      = unit(0, "pt"),
      legend.margin         = margin(t = 0, r = 0, b = 0, l = 0),
      strip.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
    )

    ## Indiv. study dir
    enrich_dir <- file.path(out_dir, "enrichment_analysis")
    if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)

    ## Save indiv. plot
    ggsaveDK(
      plot = p,
      file = file.path(enrich_dir, paste0("enrichment_", s, ".svg")),
      height = 8,
      width = 8,
      trimmed = FALSE,
      guides = TRUE,
      bg = "transparent"
    )

    ## Append to list
    enrich_list[[s]] <- p
}

## Patchwork
enrichment_plot <-
  enrich_list[[1]] + enrich_list[[2]] +
  enrich_list[[3]] + guide_area() +
  plot_layout(
    ncol   = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom", 
    legend.title.position = "left",
    legend.box            = "vertical",
    legend.box.spacing    = unit(0, "pt"),
    legend.spacing        = unit(0, "pt"),
    legend.spacing.y      = unit(0, "pt"),
    legend.spacing.x      = unit(0, "pt"),
    legend.margin         = margin(t = 0, r = 0, b = 0, l = 0),
  )


ggsaveDK(
  plot = enrichment_plot,
  file = file.path(out_dir, "Interaction_effect_enrichment.svg"),
  height = 16,
  width = 16,
  trimmed = FALSE,
  guides = TRUE,
  bg = "transparent"
)


### =================================================================================================================================
## --- Heatmap of interaction groups ---
data <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")
toi <- "mrna"
aoi <- "EAS"
data <- data[[toi]]
# Filter sample type
data <- filter_sample_type(X = data$matr, M = data$meta, s_col ="SAMPLE_TYPE", s_levels = "Primary")
data <- filter_phenotype_ancestry(X = data$matr, M = data$meta, g_col = "POOLED_SUBTYPE", a_col = "POOLED_GENETIC_ANCESTRY", g_levels = c("Basal", "non_Basal"), a_levels = c("EUR", aoi), plot = FALSE)
# Raw data
meta <- as.data.table(rbind(data$X$meta, data$Y$meta))
matr <- rbind(data$X$matr, data$Y$matr)
anno <- sig_interaction_grp_
# Filter only for EAS
anno <- anno[tech == toi & a_2 == aoi, ]
# Filter to grouped genes
matr <- matr[, colnames(matr) %in% unique(anno$feature), drop = FALSE]
# Reorder annotation to match matrix gene order
anno <- anno[match(colnames(matr), feature)]
stopifnot(all(anno$feature == colnames(matr)))
# Scale expression
matr <- t(scale(t(matr)))
matr[is.na(matr)] <- 0
# Row annotation
gene_ha <- rowAnnotation(
  "Effect group" = anno$effect_group,
  col = list("Effect group" = effect_group_colors),
  annotation_name_gp = text_gp,
  annotation_name_side = "top",
  show_annotation_name = FALSE,
  show_legend = FALSE
)
# Col annotation
sample_ha <- HeatmapAnnotation(
  Ancestry = meta$POOLED_GENETIC_ANCESTRY,
  Subtype  = meta$POOLED_SUBTYPE,
  col = list(
    Ancestry = ancestry_colors[unique(meta$POOLED_GENETIC_ANCESTRY)],
    Subtype = c(Basal = "#f8766dff", non_Basal = "#00bfc4ff")
  ),
  annotation_name_side = "right",
  annotation_name_gp = text_gp,
  show_annotation_name = FALSE,
  show_legend = FALSE
)
# Column ordering: ancestry → subtype
# Row ordering: effect-group
col_order <- order(factor(meta$POOLED_GENETIC_ANCESTRY, levels = c("EAS", "EUR")), factor(meta$POOLED_SUBTYPE, levels = c("Basal", "non_Basal")))
row_order <- order(factor(anno$effect_group, levels = effect_group_order))
col_fun <- colorRamp2(c(-3, 0, 3), c("#4575b4", "white", "#d73027"))
# Heatmap
ht <- Heatmap(
  t(matr),
  name = "Z-score",
  row_title = "All genes with interaction effect",
  column_title = "Samples",
  column_title_side = "bottom",
  col = col_fun,
  top_annotation = sample_ha,
  left_annotation = gene_ha,
  row_order = row_order,
  column_order = col_order, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_title_gp = text_gp,
  row_names_gp = gpar(fontsize = unit(7, "pt"), fontfamily = "Arial", fontface = "plain"),
  column_title_gp = text_gp
)

svg(file.path(out_dir, "expression_heatmap.svg"), width  = 10, height = 4.5)
draw(
  ht, 
  show_annotation_legend = FALSE, 
  show_heatmap_legend = TRUE,
  padding = unit(c(0, 0, 0, 0), "mm")
)
dev.off()



## --- Individual genes ----
filter_BRCA <- expression(study == "BRCA" & toi == "mrna" & phenotype == "Basal vs non-Basal" & a_2 == "EAS")
feat_to_plot <- subset_interaction[eval(filter_BRCA) & p_adj < 0.1, ]
feat_to_plot <- feat_to_plot[order(p_adj)][1:8, feature]
pval_dist <- iteration_interaction[eval(filter_BRCA) & feature %in% feat_to_plot]


feat_info <- subset_interaction[eval(filter_BRCA) & feature %in% feat_to_plot, .(feature, p_adj)]
facet_labels <- feat_info[, setNames(paste0(feature, "\n(p.adj = ", signif(p_adj, 1), ")"), feature)]
pval_dist[, facet_label := facet_labels[feature]]


## ggplot
p <- plot_pvalue_distribution(
  data = pval_dist,
  x_var = "p_value",
  fill_var = "T_obs",
  facet_col = "feature",
  facet_levels = feat_to_plot,
  x_label = "P-value",
  y_label = "Count",
  fill_label = "log2FC",
  fill_bins = 9,
  bins = 15
) +
scale_x_continuous(
  breaks = c(0, 0.5, 1),
  labels = c("0", "0.5", "1")
) +
scale_y_continuous(
  breaks = c(0, 5, 10),
  labels = c("0", "5", "10")
) +
facet_wrap(~ feature, ncol = 2) +
theme_CrossAncestryGenPhen(legend_key = 1, show_borders = TRUE) +
theme(panel.spacing = unit(0.5, "mm") )

ggsaveDK(
  plot = p,
  file = file.path(out_dir, paste0(paste(feat_to_plot, collapse = c("_")), "_pval_dist.svg")),
  height = 6,
  width = 6,
  trimmed = FALSE,
  guides = TRUE,
  bg = "transparent"
)



## --- Clusters ----
study_file <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")
techs <- c("meth", "mrna")
attrs <- c("POOLED_GENETIC_ANCESTRY", "SUBTYPE")
attr_labels <- c(
  SAMPLE_TYPE             = "Sample type",
  SUBTYPE                 = "Molecular subtype",
  POOLED_GENETIC_ANCESTRY = "Ancestry",
  SEX                     = "Sex"
)

# custom palettes for each attribute
attr_colors <- list(
  SAMPLE_TYPE = sample_type_colors,
  SUBTYPE = subtype_colors,
  POOLED_GENETIC_ANCESTRY = ancestry_colors,
  SEX = sex_colors
)

tsne_data <- lapply(techs, function(tech) {

  mat  <- study_file[[tech]]$matr
  meta <- study_file[[tech]]$meta
  meta <- meta[match(rownames(mat), meta$SAMPLE_ID), ]

  if (tech == "mrna") {
    dge  <- edgeR::DGEList(counts = t(mat))
    keep <- edgeR::filterByExpr(dge, group = factor(rep("All", ncol(dge))))
    mat  <- t(dge$counts[keep,,drop=FALSE])
    mat  <- t(edgeR::cpm(t(mat), log=TRUE))
  }

  if (tech == "meth") {
    v    <- apply(mat,2,var,na.rm=TRUE)
    v[is.na(v)] <- 0
    keep <- v > quantile(v,0.20)
    mat  <- mat[,keep,drop=FALSE]
    mat  <- beta_to_mval(mat)
  }

  tsne <- Rtsne(mat, perplexity=50, dims=2)

  list(
    coords = data.frame(TSNE1=tsne$Y[,1], TSNE2=tsne$Y[,2]),
    meta   = meta,
    mat    = mat
  )
})
names(tsne_data) <- techs

tsne_results <- lapply(techs, function(tech) {
  coords <- tsne_data[[tech]]$coords
  meta   <- tsne_data[[tech]]$meta

  panels  <- list()
  legends <- list()
  for (a in attrs) {

    keep <- !is.na(meta[[a]])
    if (a == "SUBTYPE") keep <- keep & meta$SAMPLE_TYPE != "Normal"
    df   <- cbind(coords[keep, , drop = FALSE], color = factor(meta[[a]][keep]))

    # --- full plot ---
    p_full <- ggplot(df, aes(TSNE1, TSNE2, color = color)) +
      geom_point(size = 0.08) +
      scale_color_manual(values = attr_colors[[a]], drop = FALSE) +
      labs(title = attr_labels[[a]], color = attr_labels[[a]], x = "t-SNE 1", y = "t-SNE 2") +
      theme_CrossAncestryGenPhen()

    # number of categories for this attribute
    n_levels <- length(unique(na.omit(meta[[a]])))
    if (n_levels > 3) p_leg <- p_full + guides(color = guide_legend(ncol = 2)) else p_leg <- p_full

    legends[[a]] <- cowplot::get_legend(p_leg)
    panels[[a]]  <- p_full + theme_void() + theme(legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5))

  }

  combined_plot <- wrap_plots(panels, ncol = 2) & theme(plot.margin = margin(3, 3, 3, 3))
  combined_legend <- cowplot::plot_grid(plotlist = legends, col = 3, align = "h")

  ggsaveDK(
    file   = file.path(out_dir, paste0("BRCA_", tech, "_cluster.svg")),
    plot   = combined_plot,
    width  = 4, 
    height = 2,
    bg = "transparent",
    trimmed = FALSE
  )

  ggsaveDK(
    file   = file.path(out_dir, paste0("BRCA_", tech, "_cluster_legend.svg")),
    plot   = combined_legend,
    width  = 10, 
    height = 10,
    bg = "white",
    trimmed = FALSE
  )


  list(
    plot   = combined_plot,
    legend = combined_legend
  )
})


## =========================== Subset correlation effect ===========================
subset_correlation <- rbindlist(lapply(names(results), function(study) {

  file_path <- file.path(results[[study]], "subset_correlation_effect")
  dge_files <- list.files(path = file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  ## Load all files
  all_files <- lapply(dge_files, readRDS)
  dt <- lapply(all_files, \(x) x$summary_stats)
  dt <- rbindlist(dt, use.names = TRUE, fill = TRUE)
  ## Set phenotype
  dt[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(dt$phenotype)
  dt[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  dt[, study := study] 
  ## Significance
  dt[, sig := ifelse(delta_q025 > 0 | delta_q975 < 0, "significant", "ns")]

  dt
}), fill = TRUE)

## --- Heatmap ---
p_correlation_effect <- ggplot(
  subset_correlation,
  aes(
    x    = a_2,
    y    = phenotype,
    fill = delta_mean
  )
) +
geom_tile(
  color = "black", 
  linewidth = 0.1
) +
geom_text(
  aes(label = ifelse(sig == "significant", "*", "")),
  color = "black",
  size  = 2
) +
scale_fill_gradient2(
  name     = expression(Delta ~ "Pearson (EUR - non-EUR)  "),
  low      = "#4575b4",
  mid      = "white",
  high     = "#d73027",
  midpoint = 0
) +
facet_grid(
  rows  = vars(tech),
  cols  = vars(study),
  scales = "free_y",
  space  = "free_y"
) +
labs(
  x = "Ancestry",
  y = "Phenotype"
) +
theme_CrossAncestryGenPhen(
  rotate = 45,
  show_borders = TRUE,
  show_grid = FALSE
) +
theme(
  legend.position = "bottom",
  legend.title.position = "left",
  plot.caption = element_text(size = 6, hjust = 0)
)

ggsaveDK(
  plot = p_correlation_effect,
  file = file.path(out_dir, "Fig2B_correlation_effect.svg"),
  height = 7,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

## --- Correlation of one phenotype ----
p_indv_correlation_effect <- ggplot(
data = subset_correlation[phenotype == "Basal vs non-Basal"], 
mapping = aes(
    x = a_2, 
    y = delta_mean,
    fill = sig      
  )
) +
geom_col(
  color = "black",
  linewidth = 0.1,
) +
geom_errorbar(
  mapping = aes(
    ymin = delta_q025, 
    ymax = delta_q975
  ),
  width = 0.2,
  linewidth = 0.3
) +
geom_hline(
  yintercept = 0,
  color = "black",
  linewidth = 0.1
) +
geom_text(
  data = subset_correlation[
    phenotype == "Basal vs non-Basal" & sig == "significant"
  ],
  aes(
    x = a_2,
    y = delta_q975 + (0.02 * max(abs(subset_correlation$delta_mean))),  # offset
    label = "*"
  ),
  size = 2
) +
scale_fill_manual(
  name = "",
  values = c(
    "ns" = "gray80",
    "significant" = "red"
  ),
  breaks = c("significant"),
  labels = c("significant" = "95% quantile CI excludes 0")
) +
facet_grid(
  tech ~ phenotype,
  space = "free",
  scales = "free"
) +
labs(
  x = "Ancestry",
  y = expression(Delta ~ "Pearson (EUR - non-EUR)")
) +
theme_CrossAncestryGenPhen(
  rotate = 45,
  legend_key = 1,
  show_borders = TRUE,
  show_grid = FALSE
) +
theme(legend.position = "none")

ggsaveDK(
  plot = p_indv_correlation_effect,
  file = file.path(out_dir, "Fig2B_indiv_correlation_effect.svg"),
  height = 5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)


















plot_list <- list()
for (study in names(results)){

  ## --- Results ---
  res <- results[[study]]

  ## --- Subset-interaction effect ---
  analysis  <- "subset_interaction_effect"
  file_path <- file.path(res, analysis)


  ## --- Extract ---
  dge_files <- list.files(path = file_path, pattern = "dge_res\\.rds$", recursive = TRUE,full.names = TRUE)
  all_files <- lapply(dge_files, readRDS)

  summary_list <- lapply(all_files, function(x) x$summary_stats)


  ## --- Prepare dt ---
  dt <- rbindlist(summary_list, use.names = TRUE)
  dt <- dt[coef_id == "interaction"]

  # Direction of effect
  dt[, direction := ifelse(T_obs > 0, "up_in_a2", ifelse(T_obs < 0, "up_in_a1", "no_change"))]


  # Count
  counts <- dt[, .(
      total_genes = .N,
      total_sig   = sum(p_adj < 0.1),                                
      sig_up_a1   = sum(direction == "up_in_a1" & p_adj < 0.1),
      sig_up_a2   = sum(direction == "up_in_a2" & p_adj < 0.1)
  ), by = .(contrast, g_1, g_2, a_1, a_2, method, tech, study)]
  
  counts[, pct_up_a1 := fifelse(total_sig == 0, 0, sig_up_a1 / total_sig)]
  counts[, pct_up_a2 := fifelse(total_sig == 0, 0, sig_up_a2 / total_sig)]

  counts[, pct_signed :=
    ifelse(sig_up_a1 > 0,  pct_up_a1,
    ifelse(sig_up_a2 > 0, -pct_up_a2, 0))]

  counts_long <- melt(
    counts,
    id.vars = c("contrast","g_1","g_2","a_1","a_2", "method","tech","study", "total_genes","total_sig", "sig_up_a1","sig_up_a2"),
    measure.vars = c("pct_up_a1","pct_up_a2"),
    variable.name = "direction",
    value.name = "pct"
  )

  counts_long[, pct_signed := ifelse(direction == "pct_up_a1",  pct, ifelse(direction == "pct_up_a2", -pct, 0))]


  ## --- Add phenotype ---
  counts_long[, phenotype := paste(g_1, "vs", g_2)]
  phens  <- unique(counts_long$phenotype)
  normal <- phens[grepl("^Normal", phens)]
  others <- setdiff(phens, normal)
  counts_long[, phenotype := factor(phenotype, levels = c(normal, others))]
  
  # techs and phenotypes ordered
  techs       <- sort(unique(counts_long$tech))
  phenotypes  <- levels(counts_long$phenotype)

  column_widths <- nchar(phenotypes)
  column_widths <- column_widths / sum(column_widths)


  # Loop over phenotypes
  plot_grid <- list()
  for (ti in techs) {
    for (ph in phenotypes) {


      ## --- Determine row/column position ---
      row_index <- which(techs == ti)
      col_index <- which(phenotypes == ph)
      last_col <- length(phenotypes)
      last_row <- length(techs)


      ## --- Subset for this tile ---
      dt_pct_sub   <- counts_long[tech == ti & phenotype == ph]

      dt_total_sub <- counts_long[tech == ti & paste(g_1,"vs",g_2)==ph]
      dt_total_sub <- dt_total_sub[, .SD[1], by = .(contrast, g_1, g_2, a_1, a_2, method, tech, study)]


      ## --- Top panel ---
      p_pct <- ggplot(
        data = dt_pct_sub,
        aes(
          x = a_2, 
          y = pct_signed, 
          fill = direction
        )
      ) +
      geom_col(
        color = "black", 
        linewidth = 0.1
      ) +
      geom_hline(
        yintercept = 0, 
        color = "black", 
        linewidth = 0.1
      ) +
      facet_grid(
        tech ~ phenotype,
        scales = "free",
        space = "free"
      ) +
      scale_y_continuous(
        limits = c(-1,1),
        labels = \(x) scales::percent(abs(x), accuracy=1)
      ) +
      scale_fill_manual(
        values = c(
          "pct_up_a1"="red",
          "pct_up_a2"="blue"
        ),
        labels = c(
          "pct_up_a1" = "Upregulated in EUR",
          "pct_up_a2" = "Downregulated in EUR"
        )
      ) +
      labs(
        title = ph,
        x = NULL, 
        y = "% sig. genes"
      ) +
      theme_CrossAncestryGenPhen(
        rotate = 45, 
        show_borders = TRUE, 
        show_grid = FALSE
      )

    
      ## --- Bottom panel ---
      p_total <- ggplot(
        data = dt_total_sub,
        mapping = aes(
          x = a_2, 
          y = total_sig
        )
      ) +
      geom_col(
        fill = "grey80", 
        color = "black", 
        linewidth = 0.1
      ) +
      facet_grid(
        tech ~ phenotype,
        scales = "free",
        space = "free"
      ) +
      scale_y_continuous(
        limits = c(0, NA),
        breaks = function(x) {
          pretty_vals <- scales::breaks_pretty(n = 4)(x)
          pretty_ints <- pretty_vals[!is.na(pretty_vals)]
          pretty_ints <- pretty_ints[pretty_ints >= 0]  # keep >=0 only
          pretty_ints <- unique(floor(pretty_ints))     # force integers safely
          pretty_ints
        }
      ) +
      labs(
        x = "Ancestry", 
        y = "Sig. genes\n(alpha < 0.1)"
      ) +
      theme_CrossAncestryGenPhen(
        rotate = 45, 
        show_borders = FALSE
      )


      ## --- MODIFY p_pct (TOP PANEL) ---
      p_pct <- p_pct + theme(strip.text.x = element_blank())

      # GGTITLE: shown only in first row
      if (row_index != 1) {
        p_pct <- p_pct + theme(
          plot.title = element_blank()
        )
      }

      # TECH STRIP: shown only in last column
      if (col_index != last_col) {
        p_pct <- p_pct + theme(
          strip.text.y = element_blank()
        )
      }

      # Y-AXIS only in FIRST COLUMN
      if (col_index != 1) {
        p_pct <- p_pct + theme(
          axis.title.y = element_blank()
        )
      }

      # TOP PANEL NEVER shows x-axis
      p_pct <- p_pct + theme(
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )


      ## --- MODIFY p_total (BOTTOM PANEL) ---
      # PHENOTYPE STRIP: always removed (bottom panel)
      p_total <- p_total + theme(strip.text.x = element_blank())

      # TECH STRIP: shown only in last column
      if (col_index != last_col) {
        p_total <- p_total + theme(
          strip.text.y = element_blank()
        )
      }

      # Y-AXIS only in FIRST COLUMN
      if (col_index != 1) {
        p_total <- p_total + theme(
          axis.title.y = element_blank()
        )
      }

      # X-AXIS only in LAST ROW
      if (row_index != last_row) {
        p_total <- p_total + theme(
          axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      }

      ## --- Combine vertically ---
      p_pct   <- p_pct + theme(plot.margin = margin(r = 10, b = 5))
      p_total <- p_total + theme(plot.margin = margin(r = 10, b = 5))
      tile <- wrap_plots(p_pct / p_total + plot_layout(heights = c(1.5, 1)))

      ## Store
      plot_grid[[paste(ti, ph, sep="__")]] <- tile
    }
  }

  ## --- Stack the rows ---
  rows <- lapply(techs, function(ti) {
    wrap_plots(
      plot_grid[paste(ti, phenotypes, sep="__")],
      ncol = length(phenotypes),
      widths = column_widths
    )
  })

  ## --- Stack the columns ---
  study_plot <- wrap_plots(rows, ncol = 1) +
    plot_annotation(
      theme = theme(
        plot.margin = margin(0, 0, 0, 0),
        plot.title  = element_text(
        hjust = 0.5,
        size  = 6,
        )
      )
    ) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.height  = ggplot2::unit(6 * 0.352778, "mm"),
      legend.key.width   = ggplot2::unit(6 * 0.352778, "mm"),
      legend.key.spacing = ggplot2::unit(0.5, "mm"),
      legend.margin      = ggplot2::margin(t = 0, b = 0, l = 2, r = 2, unit = "mm"),
      legend.box.spacing = ggplot2::unit(1, "mm"),
      legend.background  = ggplot2::element_rect(fill = NA, color = NA)
    )
  
  ## --- Attach to final list ---
  plot_list[[study]] <- study_plot

  space_ph <- 3
  width <- space_ph * length(phenotypes)

  ## --- Save ---
  ggsaveDK(
    plot = study_plot,
    file = file.path(out_dir, paste0(study, "_up_down.svg")),
    height = 8,
    width = 8,
    trimmed = FALSE,
    bg = "transparent"
  )
}
