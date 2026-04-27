## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
  }
)

## Configs ====================================================================
configs <- list.files(
  fil.path("configs", "tcga"), 
  pattern = "\\.json$", 
  recursive  = TRUE, 
  full.names = TRUE
)

## Parameters =================================================================
seed <- 42

# Parallel env.
future::plan(multisession, workers = 2)

## Loop: Ancestry-cancer combination ==========================================
for (config_file in configs) {

  # Read config
  cfg <- jsonlite::read_json(config_file)

  # Output directory
  out_dir <- file.path(
    cfg$out_dir, 
    "subset_correlation_effect", 
    paste0(
      cfg$tech, "_", 
      cfg$g1, "_vs_", 
      cfg$g2, "_", 
      cfg$a1, "_vs_", 
      cfg$a2
    )
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Output file
  result_file <- file.path(out_dir, "dge_res.rds")
  if (file.exists(result_file)) {
    message("\nSkipping: '", config_file, "' — results already exist!")
    message("Results file: ", result_file)
    next  
  }

  # Read data
  cat("Processing study:", cfg$study, cfg$tech, "\n")
  data <- readRDS(cfg$file)
  data <- data[[cfg$tech]]

  # Filter sample type
  data <- filter_sample_type(
    X = data$matr,
    M = data$meta,
    s_col = cfg$s_col,
    s_levels = c(cfg$s1, cfg$s2),
    verbose = FALSE
  )

  # Filter phenotype & ancestry
  data <- filter_phenotype_ancestry(
    X = data$matr, 
    M = data$meta,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    g_levels = c(cfg$g1, cfg$g2),
    a_levels = c(cfg$a1, cfg$a2),
    covariates = unlist(cfg$covariates),
    omit_na = TRUE,
    auto_cast = TRUE,
    plot = FALSE,
    verbose = FALSE
  )

  # QC plot
  ggsaveDK(
    plot = data$p, 
    file = file.path(out_dir, "imbalance.svg"),
    trimmed = FALSE,
    height = 4, 
    width = 6
  )
  
  # Filter features
  if (cfg$tech == "mrna"){
    data <- filter_by_expression(
      X = data$X$matr,
      Y = data$Y$matr,
      MX = data$X$meta,
      MY = data$Y$meta,
      g_col = cfg$g_col,
      a_col = cfg$a_col,
      any_group = FALSE,
      verbose = FALSE,
      plot = FALSE
    )
  } else if (cfg$tech == "meth") {
    data <- filter_by_methylation(
      X = data$X$matr,
      Y = data$Y$matr,
      MX = data$X$meta,
      MY = data$Y$meta,
      g_col = cfg$g_col,
      a_col = cfg$a_col,
      any_group = FALSE,
      verbose = FALSE,
      plot = FALSE
    )
  }

  # QC plot
  ggsaveDK(
    plot = data$p, 
    file = file.path(out_dir, "mean_variance.svg"), 
    trimmed = FALSE,
    height = 4, 
    width = 10
  )

  # Cluster
  p <- plot_tsne_cluster(
    X = data$X$matr,
    Y = data$Y$matr,
    MX = data$X$meta,
    MY = data$Y$meta,
    color_var = cfg$a_col,
    shape_var = cfg$g_col,
    cpm = if (cfg$tech  == "mrna") TRUE else FALSE,
    mval = if (cfg$tech == "meth") TRUE else FALSE,
    seed = seed
  )

  # QC plot
  ggsaveDK(
    plot = p, 
    file = file.path(out_dir, "tsne_cluster.svg"), 
    trimmed = FALSE,
    height = 4, 
    width = 10
  )

  # Sample density
  p <- plot_sample_density(
    X = data$X$matr,
    Y = data$Y$matr,
    cpm = if (cfg$tech == "mrna") TRUE else FALSE,
    mval = if (cfg$tech == "meth") TRUE else FALSE
  )

  # QC plot
  ggsaveDK(
    plot = p, 
    file = file.path(out_dir, "sample_density.svg"), 
    trimmed = FALSE,
    height = 6, 
    width = 10
  )

  # Cor.
  res <- subset_limma_correlation_effect(
    X = if (cfg$tech == "mrna") data$X$matr else beta_to_mval(data$X$matr),
    Y = if (cfg$tech == "mrna") data$Y$matr else beta_to_mval(data$Y$matr),
    MX = data$X$meta,
    MY = data$Y$meta,
    g_col = cfg$g_col,
    a_col = cfg$a_col,
    covariates = unlist(cfg$covariates),
    use_voom = if (cfg$tech == "mrna") TRUE else FALSE,
    method = "pearson",
    n_iter = 10,
    seed = seed,
    verbose = TRUE
  )

  res$summary_stats$method <- paste0("subset-pearson")
  res$summary_stats$tech   <- cfg$tech
  res$summary_stats$study  <- cfg$study

  # Save
  saveRDS(res, file = result_file)
}

