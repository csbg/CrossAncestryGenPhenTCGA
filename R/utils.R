## max. code width ============================================================

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

effect_grp_cols <- c(
  "Divergent effect"                = "#D62728",  
  "Minor interaction effect"        = "black", 
  "non-EUR enhances\ncancer effect" = "#1F77B4",  
  "non-EUR reverts\ncancer effect"  = "#AEC7E8",  
  "Cancer enhances\nancestry effect"= "#98DF8A",  
  "Cancer reverts\nancestry effect" = "#2CA02C"
)

weight_cols <- c(
  "Raw"    = "#4DAF4A", 
  "L1 norm." = "#984EA3"
)

## Script fun. ================================================================
vs_newline <- function(x) {
  gsub("\\s*vs\\s*", " vs\n", x)
}

rename_phenotype <- function(x) {
  x <- gsub("\\bLumA\\b", "Luminal A", x)
  x <- gsub("\\bLumB\\b", "Luminal B", x)
  x
}

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

run_one <- function(
  X, 
  Y, 
  MX, 
  MY,
  g_col,
  a_col,
  covariates,
  n_iter,
  seed,
  perm_id, 
  kind
) {
  # Collect all DGE methods
  dge_methods <- list()
  
  # Ancestry balanced
  message("\t\tlimma  (ancestry)")
  subset_res_ <- subset_limma_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    covariates = covariates,
    match = FALSE,
    use_voom = TRUE, 
    n_iter = n_iter,
    method = "mean", 
    seed = seed,
    verbose = FALSE
  )

  # All p-value agg. methods
  subset_res <- subset_res_$methods_stats
  for (m in names(subset_res)) {
    m_name <- paste0("subset-", m, " (ancestry)")

    method_res <- setDT(subset_res[[m]])
    method_res$method <- paste0("subset-", m)
    method_res$strata <- "Ancestry balanced"

    # Check and add
    stopifnot(!anyNA(method_res$T_obs))
    dge_methods[[m_name]] <- method_res
  }

  # Ancestry + cancer balance
  message("\t\tlimma  (ancestry + cancer)")
  subset_res_ <- subset_limma_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    covariates = covariates,
    match = TRUE,
    use_voom = TRUE,
    n_iter = n_iter,
    method = "mean",
    seed = seed,
    verbose = FALSE
  )

  # All p-value agg. methods
  subset_res <- subset_res_$methods_stats
  for (m in names(subset_res)) {
    m_name <- paste0("subset-", m, " (ancestry + cancer)")

    method_res <- setDT(subset_res[[m]])
    method_res$method <- paste0("subset-", m)
    method_res$strata <- "Ancestry + cancer balanced"

    # Check and add
    stopifnot(!anyNA(method_res$T_obs))
    dge_methods[[m_name]] <- method_res
  }
  
  # Limma full
  message("\t\tlimma  (imbalanced)")
  limma_full_res_ <- limma_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    covariates = covariates,
    use_voom = TRUE,
    verbose = FALSE
  )

  limma_full_res_$method <- "limma-voom"
  limma_full_res_$strata <- "Imbalanced"

  # Check and add
  stopifnot(!anyNA(limma_full_res_$T_obs))
  dge_methods[["limma-voom"]] <- setDT(limma_full_res_)
  
  # EdgeR full
  message("\t\tedger  (imbalanced)")
  edgeR_full_res_ <- edgeR_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    covariates = covariates,
    verbose = FALSE
  )

  edgeR_full_res_$method <- "edgeR-QLFTest"
  edgeR_full_res_$strata <- "Imbalanced"

  # Check and add
  stopifnot(!anyNA(edgeR_full_res_$T_obs))
  dge_methods[["edgeR-QLFTest"]] <- setDT(edgeR_full_res_)
  
  # DESeq full
  message("\t\tdeseq2 (imbalanced)")
  deseq_full_res_ <- DESeq_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    covariates = covariates,
    verbose = FALSE
  )

  deseq_full_res_$method <- "DESeq2-Wald"
  deseq_full_res_$strata <- "Imbalanced"

  # Check and add
  stopifnot(!anyNA(deseq_full_res_$T_obs))
  dge_methods[["DESeq2-Wald"]] <- setDT(deseq_full_res_)


  # Combine results for all methods
  dge_res <- rbindlist(lapply(names(dge_methods), function(m) {

    dge_method <- dge_methods[[m]]

    dge_method[, `:=`(
      method_label = m,
      perm_id = perm_id,
      kind = kind
    )]
    dge_method[]
  }), fill = TRUE)
  
  # Return
  return(dge_res[])
}

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
    )
}