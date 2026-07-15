## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
    library(ComplexHeatmap)
    library(data.table)
    library(ggnewscale)
    library(patchwork)
    library(openxlsx)
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

# Interaction effect groups (Fortelny et al. 2024)
group_effects <- function(x) {
  
  sign2 <- function(x) fifelse(x > 0, 1, fifelse(x < 0, -1, 0))

  x <- copy(x)

  x[, `:=`(
    mut         = abs(T_obs_baseline_1),       # ancestry effect (analog of mutation)
    stim_EUR    = abs(T_obs_relationship_1),   # EUR cancer effect
    stim_nonEUR = abs(T_obs_relationship_2),   # non-EUR cancer effect
    int         = abs(T_obs_interaction),

    # Signs for effect direction
    s_mut    = sign2(T_obs_baseline_1),
    s_EUR    = sign2(T_obs_relationship_1),
    s_nonEUR = sign2(T_obs_relationship_2),
    s_int    = sign2(T_obs_interaction)
  )]

  x[, effect_group := fcase(

    # (I) Divergent effect
    int > mut & int > stim_EUR & int > stim_nonEUR,
      "Divergent effect",

    # (II) Minor interaction effect
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

  # Return
  return(
    list(
      x_sig = x_sig,
      x_sum = x_sum
    )
  )
}

# Gene set enrichment
FGSEA <- function(
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

  # Return
  return(dt)
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

    draw(heatmaps[[i]], newpage = FALSE)

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
  "EUR" = "#0072B2",
  "AFR" = "#D55E00",
  "EAS" = "#56B4E9",
  "AMR" = "#E69F00",
  "SAS" = "#009E73",
  "ADMIX" = "#999999"
)

effect_grp_cols <- c(
  "Divergent effect" = "#D62728",  
  "Minor interaction effect" = "black", 
  "non-EUR enhances\ncancer effect" = "#1F77B4",  
  "non-EUR reverts\ncancer effect"  = "#AEC7E8",  
  "Cancer enhances\nancestry effect"= "#98DF8A",  
  "Cancer reverts\nancestry effect" = "#2CA02C"
)

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

## Directories ================================================================

# Result directories (per cance study)
res_dirs <- list(
  BRCA = file.path("results", "tcga", "analysis", "TCGA_BRCA"),
  UCEC = file.path("results", "tcga", "analysis", "TCGA_UCEC"),
  THCA = file.path("results", "tcga", "analysis", "TCGA_THCA")
)

# Figures
fig_dir <- file.path("results", "tcga", "figures", "subset_interaction_effect")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Tables
tab_dir <- file.path("results", "tcga", "tables", "subset_interaction_effect")
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

## Results ====================================================================

# Subset interaction results 
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

## Main 3 =====================================================================

alpha <- 0.1
formats <- c("svg", "png")

### Excel sheet ---------------------------------------------------------------

for (tech_name in c("mrna", "meth")) {

  wb <- createWorkbook()

  dt_tech <- subset_res[
    tech == tech_name
  ]

  dt_tech[, interaction_id := paste(
    g_1,
    g_2,
    a_1,
    a_2,
    sep = "__"
  )]

  groups <- split(
    dt_tech,
    dt_tech$interaction_id
  )

  for (nm in names(groups)) {

    dt <- groups[[nm]]
    dt <- dt[p_adj < 0.1]

    sheet_name <- paste(
      unique(dt$g_1),
      "vs",
      unique(dt$g_2),
      unique(dt$a_2),
      sep = "_"
    )

    sheet_name <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", sheet_name)
    sheet_name <- substr(sheet_name, 1, 31)

    addWorksheet(wb, sheet_name)

    writeData(
      wb,
      sheet = sheet_name,
      x = dt
    )
  }

  saveWorkbook(
    wb,
    file.path(
      tab_dir,
      paste0("significant_genes_", tech_name, ".xlsx")
    ),
    overwrite = TRUE
  )
}


### Panel A -------------------------------------------------------------------

# Calc. number of DEGs
nr_degs <- subset_res[
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

# Save
write.csv(
  nr_degs,
  file = file.path(
    tab_dir, 
    "main3_panel_A.csv"
  ),
  row.names = FALSE
)

# Save Excel
openxlsx::write.xlsx(
  nr_degs,
  file = file.path(
    tab_dir,
    "main3_panel_A.xlsx"
  )
)

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
ggsaveDK(
  plot = main3_panel_A,
  file = file.path(
    fig_dir, 
    "main3_panel_A.svg"
  ),
  height = 5,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel B -------------------------------------------------------------------

# Cancer effect bigger in EUR
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
ggsaveDK(
  plot = main3_panel_B,
  file = file.path(
    fig_dir, 
    "main3_panel_B.svg"
  ),
  height = 5,
  width = 9.5,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel C -------------------------------------------------------------------

# Group effects 
effect_grps_pct <- summarize_effects(subset_res, alpha = alpha)

# Save
write.csv(
  effect_grps_pct$x_sum,
  file = file.path(
    tab_dir, 
    "main3_panel_C.csv"
  ),
  row.names = FALSE
)

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
ggsaveDK(
  plot = main3_panel_C,
  file = file.path(
    fig_dir, 
    "main3_panel_C.svg"),
  height = 6,
  width = 8,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel D -------------------------------------------------------------------

# Plot
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
    data = subset_res[
      study == "BRCA" & 
      tech == "mrna" & 
      phenotype == "Basal vs non-Basal" & 
      coef_id %in% c("relationship_1", "relationship_2", "interaction") &
      feature %in% subset_res[
        study == "BRCA" & 
        tech == "mrna" & 
        phenotype == "Basal vs non-Basal" & 
        coef_id == "interaction" & 
        p_adj < 0.1,
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
ggsaveDK(
  plot = main3_panel_D,
  file = file.path(
    fig_dir, 
    "main3_panel_D.svg"
  ),
  height = 6,
  width = 11,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel E -------------------------------------------------------------------

# BRCA study
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
  subset_res[
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
    sig = fifelse(p_adj < 0.1, "*", "n.s.")
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
      p_adj < 0.1, 
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
ggsaveDK(
  plot = main3_panel_E,
  file = file.path(
    fig_dir, 
    "main3_panel_E.svg"
  ),
  height = 6,
  width = 11,
  trimmed = FALSE,
  bg = "transparent"
)

### Panel F -------------------------------------------------------------------

# Fast gene set enrichment analysis
data_fgsea <- rbindlist(
  subset_res[
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
    x = NES,
    y = reorder(pathway, NES),
    fill = a_2
  )
) +
geom_col(
  position = position_dodge(
    width = 0.8
  ),
  width = 0.7,
  color = "black",
  linewidth = 0.1,
  na.rm = TRUE
) +
scale_fill_manual(
  values = ancestry_cols
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
  plot = main3_panel_F,
  file = file.path(
    fig_dir, 
    "main3_panel_F.svg"
  ),
  height = 6,
  width = 7,
  trimmed = FALSE,
  bg = "transparent"
)

# Plot
main3_panel_F.1 <- ggplot(
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
ggsaveDK(
  plot = main3_panel_F.1,
  file = file.path(
    fig_dir, 
    "main3_panel_F.1.svg"
  ),
  height = 6,
  width = 7,
  trimmed = FALSE,
  bg = "transparent"
)

### Supp. 1 -------------------------------------------------------------------

# DEGs per coef
nr_degs <- subset_res[
  coef_type == "interaction",
  .(i_count = uniqueN(feature[p_adj <= alpha])),
  by = .(study, tech, a_2, g_1, g_2, phenotype)
][
  subset_res[
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
      fig_dir, 
      paste0(
        "main3_supp1.", 
        ext
      )
    ),
    height = 8,
    width = 16,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}

### Supp. 2 -------------------------------------------------------------------

# Plot
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
      fig_dir, 
      paste0(
        "main3_supp2.", 
        ext
      )
    ),
    height = 8,
    width = 16,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}

### Supp. 3 -------------------------------------------------------------------

# Correlation of effect sizes
logFC_cor <- effect_grps_pct$x_sig[
  study == "BRCA" & 
  tech == "mrna" & 
  phenotype == "Basal vs non-Basal",
  .(feature, effect_group)
][
  dcast(
    subset_res[
      study == "BRCA" & 
      tech == "mrna" & 
      phenotype == "Basal vs non-Basal" & 
      coef_id == "interaction" &
      feature %in% subset_res[
        study == "BRCA" & 
        tech == "mrna" & 
        phenotype == "Basal vs non-Basal" & 
        coef_id == "interaction" & 
        p_adj < 0.1, 
        feature
      ]
    ],
    feature ~ a_2,
    value.var = "T_obs"
  ),
  on = "feature"
]

# Plot
split_list <- split(logFC_cor, logFC_cor$effect_group)
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
    plot = heatmaps,
    file = file.path(
      fig_dir, 
      paste0("main3_supp3.", ext)
    ),
    height = 9,
    width = 16,
    ncol = 2,
    bg = "transparent",
    dpi = 300
  )
}

### Supp. 4 -------------------------------------------------------------------

# Nr. enriched pathways
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
      fig_dir, 
      paste0(
        "main3_supp4.", 
        ext
      )
    ),
    height = 7,
    width = 8,
    trimmed = FALSE,
    bg = "transparent",
    dpi = 300
  )
}