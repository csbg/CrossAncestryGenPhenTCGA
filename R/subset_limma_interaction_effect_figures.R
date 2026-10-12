## max. code width ============================================================


## Libraries ==================================================================
suppressPackageStartupMessages({
  library(CrossAncestryGenPhen)
  library(ComplexHeatmap)
  library(data.table)
  library(ggnewscale)
  library(openxlsx)
  library(patchwork)
  library(yardstick)
  library(circlize)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
  library(ggh4x)
})

## utils.R ====================================================================
source("R/utils.R")

## Settings ===================================================================
alpha   <- 0.1
formats <- c("svg", "png")

## Databases ==================================================================

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



# Paper =======================================================================
message("\nCreating paper figures...")

## Results ====================================================================
### directories               -------------------------------------------------
res_dirs <- list(
  BRCA = file.path("results", "tcga", "analysis", "TCGA_BRCA"),
  UCEC = file.path("results", "tcga", "analysis", "TCGA_UCEC"),
  THCA = file.path("results", "tcga", "analysis", "TCGA_THCA")
)
sum_dir <- file.path("results", "tcga", "analysis", "summary_subset_limma_interaction_effect")
stopifnot(length(unique(dirname(c(unlist(res_dirs), sum_dir)))) == 1)

cfg_dirs <- list(
  BRCA = file.path("configs", "tcga", "TCGA_BRCA"),
  UCEC = file.path("configs", "tcga", "TCGA_UCEC"),
  THCA = file.path("configs", "tcga", "TCGA_THCA")
)

# Dirs
fig_dir <- file.path("results", "tcga", "figures", "subset_limma_interaction_effect")
tab_dir <- file.path("results", "tcga", "tables", "subset_limma_interaction_effect")

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

### ancestry-specific effects -------------------------------------------------
if (!file.exists(file = file.path(sum_dir, "summary_limma_res.rds"))) {

  message("Loading 'limma_res' from individual results.")
  limma_res <- rbindlist(lapply(names(res_dirs), function(study) {
    file_path <- file.path(res_dirs[[study]], "subset_limma_interaction_effect")
    dge_files <- list.files(file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)

    cfg <- list.files(cfg_dirs[[study]])
    message(sprintf("%s found %d/%d files.", study, length(dge_files), length(cfg)))

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

  # Save
  saveRDS(
    object = limma_res,
    file = file.path(
      sum_dir,
      "summary_limma_res.rds"
    )
  )

} else {
  message("Loading 'limma_res' from: ", sum_dir)
  limma_res <- readRDS(file = file.path(sum_dir, "summary_limma_res.rds"))
}


## Main Figure 3 ==============================================================
### panel A -------------------------------------------------------------------
nr_degs <- limma_res[
  ,
  .(
    alpha = alpha,
    nr_degs = sum(p_adj < alpha)
  ),
  by = .(
    study, phenotype, 
    tech, coef_id, a_2
  )
]

# Plot
main3_panel_A <- ggplot(
  mapping = aes(
    x = a_2, 
    y = rename_phenotype(phenotype)
  )
) +
geom_tile(
  data = nr_degs[
    coef_id %in% "interaction" & 
    alpha == alpha &
    nr_degs == 0
  ],
  mapping = aes(
    fill = "No effect"
  ),
  color = "black",
  linewidth = 0.1,
  show.legend = TRUE
) +
scale_fill_manual(
  name = NULL,
  values = c("No effect" = "grey90"),
  guide = guide_legend(
    override.aes = list(color = "black")
  )
) +
ggnewscale::new_scale_fill() +
geom_tile(
  data = nr_degs[
    coef_id %in% "interaction" &
    alpha == alpha &
    nr_degs != 0
  ],
  mapping = aes(
    fill = log10(nr_degs + 1)
  ),
  color = "black",
  linewidth = 0.1
) +
scale_fill_gradientn(
  name = "# DEGs with\ninteraction effect",
  colours = c("white", "forestgreen"),
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
    tech = c(
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
  legend.spacing.y = unit(10, "pt"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_A,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_A.", ext)
    ),
    height = 5,
    width = 8,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_A,
  file = file.path(
    fig_dir_main, 
    "main3_panel_A.rds"
  )
)

### panel B -------------------------------------------------------------------
effect_grp_EUR <- data.table(
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

effect_grp_EUR[, 
  cancer_dom := "Cancer effect\nEUR > non-EUR"
]

# Cancer effect smaller in EUR
effect_grp_nonEUR <- data.table(
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

effect_grp_nonEUR[,
  cancer_dom := "Cancer effect\nEUR < non-EUR"
]

# Combine
effect_grp_dt <- rbind(
  effect_grp_EUR, 
  effect_grp_nonEUR
)

effect_grp_dt[, `:=`(
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
effect_grp_labels <- effect_grp_dt[
  phenotype == "Cancer" &
  !effect_group %in% c(
    "Minor interaction effect",
    "Cancer enhances\nancestry effect"
  )
]

effect_grp_labels[, max_value := max(value), by = .(effect_group, cancer_dom)]
effect_grp_labels[, direction := ifelse(value == max_value, 1, -1)]
effect_grp_labels[, nudge_y := 0.5 * direction]

# Plot
main3_panel_B <- ggplot(
  data = effect_grp_dt[
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
    effect_grp_cols, 
    sub(
      "Minor interaction effect", 
      "EUR_black", 
      names(effect_grp_cols)
    )
  )
) +
geom_text(
  data = effect_grp_labels,
  aes(
    label = ancestry,
    color = color_group,
    y = value + nudge_y
  ),
  nudge_x = -0.05,
  hjust = 0,
  size = 1.8,
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
  panel.spacing.y = unit(0.15, "lines")
)

# Save
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_B,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_B.", ext)
    ),
    height = 5,
    width = 9.5,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_B,
  file = file.path(
    fig_dir_main, 
    "main3_panel_B.rds"
  )
)

### panel C -------------------------------------------------------------------
effect_grps_pct <- summarize_effects(limma_res, alpha = alpha)

# Plot
main3_panel_C <- ggplot(
  data = effect_grps_pct$x_sum[
    study == "BRCA" &
    phenotype == "Basal vs non-Basal"
  ],
  mapping = aes(
    x = a_2, 
    y = pct, 
    fill = effect_group
  )
) +
geom_col(
  na.rm = TRUE
) +
geom_text(
  data = effect_grps_pct$x_sum[
    study == "BRCA" & 
    phenotype == "Basal vs non-Basal",
    .(all_zero = all(pct == 0)),
    by = .(cancer_dom, tech, a_2)
  ][all_zero == TRUE],
  mapping = aes(
    x = a_2, 
    label = "No differences"
  ),
  y = 50,
  angle = 90,
  hjust = 0.5,
  vjust = 0.5,
  size = 1.8,
  inherit.aes = FALSE
) +
scale_fill_manual(
  values = effect_grp_cols
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
  y = "Share of genes per effect group\n(% of genes with interaction effect)"
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
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_C,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_C.", ext)
    ),
    height = 6,
    width = 8,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_C,
  file = file.path(
    fig_dir_main, 
    "main3_panel_C.rds"
  )
)

### panel D -------------------------------------------------------------------
main3_panel_D <- wrap_plots(
  ggplot(
    data = copy(effect_grps_pct$x_sig)[
      study == "BRCA" & 
      tech == "mrna" & 
      phenotype == "Basal vs non-Basal"
    ][order(effect_group, feature)][, 
      feature := factor(feature, levels = feature)
    ],
    mapping = aes(
      x = feature,
      y = "Effect group",              
      fill = effect_group
    )
  ) +
  geom_tile() +
  scale_fill_manual(
    values = effect_grp_cols
  ) +
  labs(
    y = "Effect group"
  ) +
  theme_CrossAncestryGenPhen(
    show_axis = FALSE, 
    rotate = 90
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    plot.margin = margin(0, 0, 0, 0)
  ),
  ggplot(
    data = limma_res[
      study == "BRCA" & 
      tech == "mrna" & 
      phenotype == "Basal vs non-Basal" & 
      coef_id %in% c("relationship_1", "relationship_2", "interaction") &
      feature %in% limma_res[
        study == "BRCA" & 
        tech == "mrna" & 
        phenotype == "Basal vs non-Basal" & 
        coef_id == "interaction" & 
        p_adj < alpha,
        feature
      ]
    ][, 
      feature := factor(
      feature,
      levels = effect_grps_pct$x_sig[
        study == "BRCA" & 
        tech == "mrna" & 
        phenotype == "Basal vs non-Basal"
      ][order(effect_group, feature), unique(feature)]
    )
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
    axis.text.x = element_text(size = 4),
    panel.spacing.y = unit(0.15, "lines"),
    panel.spacing.x = unit(0.15, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  ),
  ncol = 1
) + plot_layout(heights = c(0.05, 1))

# Save
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_D,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_D.", ext)
    ),
    height = 6,
    width = 11,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_D,
  file = file.path(
    fig_dir_main, 
    "main3_panel_D.rds"
  )
)

### panel E -------------------------------------------------------------------
mrna <- readRDS("download/GDCportal/TCGA_BRCA/TCGA_BRCA_omics_layers.rds")[["mrna"]]
mrna <- filter_sample_type(
  X = mrna$matr, 
  M = mrna$meta, 
  s_col = "SAMPLE_TYPE", 
  s_levels = "Primary",
  verbose = FALSE
)

# Keys & Features
keys <- c(
  "SAMPLE_ID", 
  "POOLED_SUBTYPE", 
  "POOLED_GENETIC_ANCESTRY"
)

feats <- c(
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

# Merge statistics and effect groups
data_expr <- effect_grps_pct$x_sig[
  study == "BRCA" &
  tech == "mrna" &
  phenotype == "Basal vs non-Basal" &
  feature %in% feats,
  .(study, tech, phenotype, a_2, feature, effect_group)
][
  limma_res[
    coef_id == "interaction" &
    study == "BRCA" &
    tech == "mrna" &
    phenotype == "Basal vs non-Basal" &
    feature %in% feats,
    .(study, tech, phenotype, a_2, feature, p_adj)
  ][
    melt(
      as.data.table(mrna$meta)[
        as.data.table(mrna$matr, keep.rownames = keys[[1]])[
          , c(keys[[1]], feats), with = FALSE
        ],
        on = keys[[1]]
      ],
      id.vars = keys,
      measure.vars = feats,
      variable.name = "feature",
      value.name = "expression"
    )[
      , `:=`(
        POOLED_SUBTYPE = fifelse(
          POOLED_SUBTYPE == "non_Basal",
          "non-Basal",
          POOLED_SUBTYPE
        ),
        a_2 = POOLED_GENETIC_ANCESTRY,
        study = "BRCA",
        tech = "mrna",
        phenotype = "Basal vs non-Basal",
        POOLED_GENETIC_ANCESTRY = NULL
      )
    ][
      !is.na(POOLED_SUBTYPE) &
      !is.na(a_2)
    ][, 
      if (uniqueN(POOLED_SUBTYPE) == 2) .SD, by = a_2
    ][
      , zscore := scale(expression), by = feature
    ],
    on = .(study, tech, phenotype, a_2, feature)
  ],
  on = .(study, tech, phenotype, a_2, feature)
][
  , `:=`(
    a_2 = factor(a_2, levels = c("EUR", setdiff(unique(a_2), "EUR"))),
    feature = factor(feature, levels = feats),
    sig = fifelse(p_adj < alpha, "*", "n.s.")
  )
][, 
  x := as.numeric(a_2)
][]

# Compute zoomed ylims
coef_val <- 0.5
ylims <- data_expr[
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
][
  , .(
    ymin = min(wmin, na.rm = TRUE),
    ymax = max(wmax, na.rm = TRUE)
  ), by = feature
][
  , {
    range <- ymax - ymin
    
    pad_bottom <- 0.0 * range   
    pad_top    <- 0.0 * range  
    
    .(
      lims = list(c(
        ymin - pad_bottom,
        ymax + pad_top
      ))
    )
  }, by = feature
][
  , feature := factor(
    feature, 
    levels = feats
  )
][order(feature)]

yscales <- lapply(ylims$lims, function(lim) {
  scale_y_continuous(
    limits = lim,
    breaks = pretty(lim, n = 3)
  )
})

# Plot
main3_panel_E <- ggplot(
  data = data_expr,
  mapping = aes(
    x = a_2,
    y = zscore,
  )
) +
geom_rect(
  data = unique(
    data_expr[
      p_adj < alpha, 
      .(feature, x, effect_group)
    ]
  ),
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
  values = effect_grp_cols,
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
) +
facetted_pos_scales(
  y = yscales
) + labs(
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
  strip.text = element_text(margin = margin(0, 0, 0, 0))
)

# Save
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_E,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_E.", ext)
    ),
    height = 6,
    width = 11,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_E,
  file = file.path(
    fig_dir_main, 
    "main3_panel_E.rds"
  )
)

### panel F -------------------------------------------------------------------
data_fgsea <- rbindlist(
  limma_res[
    coef_id == "interaction"
  ][
    tech == "meth",
    feature := tstrsplit(
      probe2gene[.SD, on = "feature", gene],
      ";",
      fixed = TRUE,
      keep = 1L
    )
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
      res  <- FGSEA(
        ranks = ranks,
        metadata = meta,
        pathways = fgsea_pathways
      )

      .(res = list(res))

    },
    by = .(study, phenotype, tech, a_2)
  ][["res"]]
)

# Plot
main3_panel_F <- ggplot(
  data = data_fgsea[
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
    x = a_2,
    y = reorder(pathway, NES),
    color = NES,
    size = -log10(padj)
  )
) +
geom_point() +
scale_size_continuous(
  name = expression(P.adj ~ "-log"[10] * ""),
  breaks = scales::pretty_breaks(n = 4),
  range = c(0.5, 1.8)
) +
scale_color_gradient2(
  high = "purple",
  mid = "white",
  low = "forestgreen"
) +
facet_grid(
  rows = vars(tech),
  scales = "free_y",
  space = "free",
  labeller = labeller(
    tech = c(
      meth = "Methylation", 
      mrna = "Expression"
    )
  )
) +
labs(
 color = "NES",
 x = "Ancestry",
 y = "Hallmark pathway"
) +
theme_CrossAncestryGenPhen(
  legend_key = 1,
  rotate = 45,
  show_border = TRUE,
  show_grid = FALSE
) +
theme(
  panel.spacing.x = unit(0.15, "lines"),
  panel.spacing.y = unit(0.15, "lines"),
  legend.margin = margin(0, 0, 0, 0)
)

# Save
for (ext in formats){
  ggsaveDK(
    plot = main3_panel_F,
    file = file.path(
      fig_dir_main, 
      paste0("main3_panel_F.", ext)
    ),
    height = 6,
    width = 7,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_panel_F,
  file = file.path(
    fig_dir_main, 
    "main3_panel_F.rds"
  )
)


## Supplementary Figure 3 =====================================================
### supp. 1 -------------------------------------------------------------------
nr_degs <- limma_res[
  coef_type == "interaction",
  .(i_count = uniqueN(feature[p_adj <= alpha])),
  by = .(study, tech, a_2, g_1, g_2, phenotype)
][
  limma_res[
    coef_type == "baseline",
    .(b_count = uniqueN(feature[p_adj <= alpha])),
    by = .(study, tech, a_2, g_1, g_2)
  ],
  on = .(study, tech, a_2, g_1, g_2)
][
  , `:=`(
    lb = log10(b_count + 1),
    li = log10(i_count + 1)
  )
][
  , alpha := alpha
][]

# Correlation
nr_degs_cor <- nr_degs[
  , `:=`(
    ancestry_n = .N,
    ancestry_cor = cor_safe(lb, li, .N)
  ),
  by = .(tech, a_2)
][
  , `:=`(
    across_n = .N,
    across_cor = cor_safe(lb, li, .N)
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
main3_supp1 <- wrap_plots(
  ggplot(
    data = nr_degs,
    mapping = aes(
      x = lb,
      y = li,
      color = a_2
    )
  ) +
  geom_polygon(
    data = {

      x <- nr_degs$lb
      y <- nr_degs$li
      
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
      x = lb,
      y = li
    ),
    linewidth = 0.7,
    method = "lm",
    formula = y ~ x,
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
    x = "Nr. of significant genes within cancer subtype\n(the ancestry baseline effect)",
    y = "Nr. of significant ancestry-specific genes\n(the interaction effect)",
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
    y = "Pearson correlation of\nnr. baseline genes vs. nr. interaction genes"
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
    plot = main3_supp1,
    file = file.path(
      fig_dir_supp, 
      paste0("main3_supp1.", ext)
    ),
    height = 8,
    width = 16,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_supp1,
  file = file.path(
    fig_dir_supp, 
    "main3_supp1.rds"
  )
)

### supp. 2 -------------------------------------------------------------------
main3_supp2 <- wrap_plots(
  ggplot(
    data = effect_grps_pct$x_sum[
      , .(x_sum = sum(N, na.rm = TRUE)),
      by = .(tech, a_2, effect_group)
    ],
    mapping = aes(
      x = effect_group,
      y = x_sum,
      fill = a_2
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = ancestry_cols
  ) +
  facet_grid(
    rows = vars(tech),
    scales = "free",
    space = "free",
    labeller = labeller(
      tech  = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    fill = "Ancestry",
    x = "Effect group",
    y = "Nr. ancestry-specific genes\nper effect group"
  ) +
  theme_CrossAncestryGenPhen(
    rotate = 45, 
    show_borders = TRUE
  ) +
  theme(
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ),
  ggplot(
    data = effect_grps_pct$x_sum[
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
  geom_col(
    na.rm = TRUE
  ) +
  geom_text(
    data = effect_grps_pct$x_sum[,
      .(all_zero = all(pct == 0)),
      by = .(cancer_dom, tech, a_2)
    ][all_zero == TRUE],
    mapping = aes(
      x = a_2, 
      label = "No differences"
    ),
    y = 50,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    size = 1.8,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = effect_grp_cols
  ) +
  facet_grid(
    rows = vars(cancer_dom),
    cols = vars(tech),
    labeller = labeller(
      tech  = c(meth = "Methylation", mrna = "Expression")
    )
  ) +
  labs(
    fill = "Effect group",
    x = "Ancestry",
    y = "Share of genes per effect group\n(% of genes with interaction effect)"
  ) +
  theme_CrossAncestryGenPhen(
    rotate = 45, 
    show_borders = TRUE
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
    plot = main3_supp2,
    file = file.path(
      fig_dir_supp, 
      paste0("main3_supp2.", ext)
    ),
    height = 8,
    width = 16,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_supp2,
  file = file.path(
    fig_dir_supp, 
    "main3_supp2.rds"
  )
)

### supp. 3 -------------------------------------------------------------------
logFC_cor <- effect_grps_pct$x_sig[
  study == "BRCA" & 
  tech == "mrna" & 
  phenotype == "Basal vs non-Basal",
  .(feature, effect_group)
][
  dcast(
    limma_res[
      study == "BRCA" & 
      tech == "mrna" & 
      phenotype == "Basal vs non-Basal" & 
      coef_id == "interaction" &
      feature %in% limma_res[
        study == "BRCA" & 
        tech == "mrna" & 
        phenotype == "Basal vs non-Basal" & 
        coef_id == "interaction" & 
        p_adj < alpha, 
        feature
      ]
    ],
    feature ~ a_2,
    value.var = "T_obs"
  ),
  on = "feature"
][effect_group != "Minor interaction effect", ]

# Plot
split_list <- split(logFC_cor, logFC_cor$effect_group)
split_list <- split_list[sapply(split_list, nrow) > 0]

main3_supp3 <- lapply(names(split_list), function(grp) {
  
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
    width = unit(2.5, "cm"),
    height = unit(2.5, "cm"),
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

  ht
})

# Save
for (ext in formats) {
  save_heatmapsDK(
    plot = main3_supp3,
    file = file.path(
      fig_dir_supp, 
      paste0("main3_supp3.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 2,
    bg = "transparent",
    dpi = 300
  )
}

saveRDS(
  object = main3_supp3,
  file = file.path(
    fig_dir_supp, 
    "main3_supp3.rds"
  )
)

### supp. 4 -------------------------------------------------------------------
main3_supp4 <- ggplot(
  data = data_fgsea[
    padj < 0.05,
    .(n_pathways = uniqueN(pathway)),
    by = .(tech, a_2)
  ],
  mapping = aes(
    x = a_2,
    y = n_pathways,
    fill = a_2
  )
) +
geom_col(
  na.rm = TRUE
) +
scale_fill_manual(
  values = ancestry_cols
) +
facet_grid(
  cols = vars(tech),
  labeller = labeller(
    tech  = c(meth = "Methylation", mrna = "Expression")
  )
) +
labs(
  fill = "Ancestry",
  x = "Ancestry",
  y = "Nr. enriched pathways"
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
for (ext in formats) {
  ggsaveDK(
    plot = main3_supp4,
    file = file.path(
      fig_dir_supp, 
      paste0("main3_supp4.", ext)
    ),
    height = 7,
    width = 8,
    trimmed = FALSE,
    bg = "transparent"
  )
}

saveRDS(
  object = main3_supp4,
  file = file.path(
    fig_dir_supp, 
    "main3_supp4.rds"
  )
)


## Supplementary Table 3 ======================================================
### ancestry-specific genes ---------------------------------------------------
for (t in c("mrna", "meth")) {

  wb <- createWorkbook()
  dt <- limma_res[tech == t]
  dt[, interaction_id := paste(g_1, g_2, a_1, a_2, sep = "__")]

  groups <- split(dt, dt$interaction_id)
  for (nm in names(groups)) {

    dt <- groups[[nm]]
    dt <- dt[p_adj < alpha]

    sheet <- paste(unique(dt$g_1), "vs",unique(dt$g_2), unique(dt$a_2), sep = "_")
    sheet <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", sheet)
    sheet <- substr(sheet, 1, 31)
    addWorksheet(wb, sheet)

    writeData(wb, sheet = sheet, x = dt)
  }

  # Save
  saveWorkbook(
    wb,
    file.path(
      tab_dir_supp,
      paste0("significant_genes_", t, ".xlsx")
    ),
    overwrite = TRUE
  )
}


