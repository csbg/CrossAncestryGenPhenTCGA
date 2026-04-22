library(CrossAncestryGenPhen)
library(future)
library(data.table)
 
# ======= Parallel environment =======
future::plan(multisession, workers = 2)


# ======= DGE configurations =======
configs <- list.files("configs", pattern = "\\.json$", recursive  = TRUE, full.names = TRUE)

# ======= DGE loop =======
for (config_file in configs) {
  
  ## --- Load config ---
  cfg <- jsonlite::read_json(config_file)
  message("\n================================================================================================================================")


  ## --- Make dir ---
  out_dir <- file.path(cfg$out_dir, "subset_prediction_effect", paste0(cfg$tech, "_", cfg$g1, "_vs_", cfg$g2, "_", cfg$a1, "_vs_", cfg$a2))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


  ## --- Define results file ---
  result_file <- file.path(out_dir, "dge_res.rds")


  ## --- Check if results already exist ---
  if (file.exists(result_file)) {
    message("\nSkipping: '", config_file, "' — results already exist!")
    message("Results file: ", result_file)
    next  
  }

  ## --- Message ---
  message("\nStudy: ", cfg$study, "  tech: ", cfg$tech, "  save: ", out_dir)


  ## --- Load file ---
  data <- readRDS(cfg$file)


  ## --- Filter tech ---
  data <- data[[cfg$tech]]


  ## --- Filter sample type ---
  s_col <- cfg$s_col; s1 <- cfg$s1; s2 <- cfg$s2

  data <- filter_sample_type(
    X = data$matr,
    M = data$meta,
    s_col = cfg$s_col,
    s_levels = c(s1, s2)
  )


  ## --- Filter by ancestry and phenotype ---
  g_col <- cfg$g_col; g1 <- cfg$g1; g2 <- cfg$g2
  a_col <- cfg$a_col; a1 <- cfg$a1; a2 <- cfg$a2
  # Covariates
  covariates <- unlist(cfg$covariates)
  # Reproducibility
  seed <- cfg$seed

  study <- filter_phenotype_ancestry(
    X = data$matr,
    M = data$meta,
    g_col = g_col,
    a_col = a_col,
    g_levels = c(g1, g2),
    a_levels = c(a1, a2),
    covariates = covariates,
    omit_na = TRUE,
    auto_cast = TRUE,
    plot = FALSE,
    verbose = TRUE
  )

  # Save
  ggsaveDK(
    plot = study$p, 
    file = file.path(out_dir, "imbalance.svg"),
    trimmed = FALSE,
    height = 4, 
    width = 6
  )

  ## --- Filter features ---
  if (cfg$tech == "mrna"){
    study <- filter_by_expression(
      X = study$X$matr,
      Y = study$Y$matr,
      MX = study$X$meta,
      MY = study$Y$meta,
      g_col = g_col,
      a_col = a_col,
      any_group = FALSE,
      verbose = TRUE,
      plot = FALSE
    )
  } else if (cfg$tech == "meth") {
    study <- filter_by_methylation(
      X = study$X$matr,
      Y = study$Y$matr,
      MX = study$X$meta,
      MY = study$Y$meta,
      g_col = g_col,
      a_col = a_col,
      any_group = FALSE,
      verbose = TRUE,
      plot = FALSE
    )
  }

  # Save
  ggsaveDK(
    plot = study$p, 
    file = file.path(out_dir, "mean_variance.svg"), 
    trimmed = FALSE,
    height = 4, 
    width = 10
  )


  ## --- Cluster ---
  p <- plot_tsne_cluster(
    X = study$X$matr,
    Y = study$Y$matr,
    MX = study$X$meta,
    MY = study$Y$meta,
    color_var = a_col,
    shape_var = g_col,
    cpm = if (cfg$tech  == "mrna") TRUE else FALSE,
    mval = if (cfg$tech == "meth") TRUE else FALSE,
    seed = seed
  )

  # Save
  ggsaveDK(
    plot = p, 
    file = file.path(out_dir, "tsne_cluster.svg"), 
    trimmed = FALSE,
    height = 4, 
    width = 10
  )


  ## --- Sample density ---
  p <- plot_sample_density(
    X = study$X$matr,
    Y = study$Y$matr,
    cpm = if (cfg$tech == "mrna") TRUE else FALSE,
    mval = if (cfg$tech == "meth") TRUE else FALSE
  )

  # Save
  ggsaveDK(
    plot = p, 
    file = file.path(out_dir, "sample_density.svg"), 
    trimmed = FALSE,
    height = 6, 
    width = 10
  )


  ## --- Subset ---
  method = "auc"
  res <- subset_logistic_prediction_effect(
    X = study$X$matr,
    Y = study$Y$matr,
    MX = study$X$meta,
    MY = study$Y$meta,
    g_col = g_col,
    a_col = a_col,
    n_folds = 5,
    n_models = 5,
    method = method,
    n_iter = 10,
    seed = seed,
    verbose = TRUE
  )
  res$summary_stats$method  <- paste0("subset-", method)
  res$summary_stats$tech    <- cfg$tech
  res$summary_stats$study   <- cfg$study


  ## --- Save ---
  saveRDS(res, file = result_file)
}

