## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
  }
)

## Configs ====================================================================
configs <- list.files(
  file.path("configs", "tcga"), 
  pattern = "\\.json$", 
  recursive  = TRUE, 
  full.names = TRUE
)

## Parameters =================================================================
seed <- 42

# Parallel env.
future::plan(future::multisession, workers = 2)

## Loop: Ancestry-cancer combination ==========================================
all_n_list <- list()
for (config_file in configs) {
  
  # Read config
  cfg <- jsonlite::read_json(config_file)
  comp <-  paste0(
    cfg$tech, "_", 
    cfg$g1, "_vs_",  
    cfg$g2, "_", 
    cfg$a1, "_vs_", 
    cfg$a2
  )

  # Output directory
  out_dir <- file.path(
    cfg$out_dir, 
    "subset_prediction_effect", 
    comp
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

  # Sample sizes (n)
  .get_n <- function(meta, g_col, a_col) {
    
    g_groups <- unique(na.omit(meta[[g_col]]))
    g_counts <- sapply(g_groups, function(g) sum(meta[[g_col]] == g, na.rm = TRUE))
    ancestry_label <- unique(na.omit(meta[[a_col]]))[1]
    
    data.frame(
      a_col = ancestry_label,
      g_col = as.character(g_groups),
      n     = as.numeric(g_counts)
    )
  }

  n_X <- .get_n(meta = data$X$meta, g_col = cfg$g_col, a_col = cfg$a_col)
  n_Y <- .get_n(meta = data$Y$meta, g_col = cfg$g_col, a_col = cfg$a_col)
  final_n <- rbind(n_X, n_Y)
  
  # Add meta
  final_n$study <- cfg$study
  final_n$tech  <- cfg$tech
  final_n$comp  <- comp


  # Save
  write.csv(
    final_n, 
    file = file.path(out_dir, "sample_n.csv"), 
    row.names = FALSE
  )

  all_n_list[[config_file]] <- final_n

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

  # Pred.
  res <- subset_logistic_prediction_effect(
    X = data$X$matr,
    Y = data$Y$matr,
    MX = data$X$meta,
    MY = data$Y$meta,
    g_col = cfg$g_col,
    a_col = cfg$a_col,
    n_folds = 5,
    n_models = 5,
    method = "auc",
    n_iter = 10,
    seed = seed,
    verbose = TRUE
  )

  res$summary_stats$method <- paste0("subset-auc")
  res$summary_stats$tech   <- cfg$tech
  res$summary_stats$study  <- cfg$study

  # Save
  saveRDS(res, file = result_file)
}

# Summary
if (length(all_n_list) > 0) {
  
  all_n   <- do.call(rbind, all_n_list)

  out_dir <- file.path("results/tcga/analysis/summary_subset_prediction_effect")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  write.csv(
    all_n,
    file = file.path(
      out_dir,
      "summary_sample_n.csv"
    ), 
    row.names = FALSE
  )
}

