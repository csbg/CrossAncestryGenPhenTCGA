## max. code width ============================================================


## Libraries ==================================================================
suppressPackageStartupMessages({
  library(CrossAncestryGenPhen)
  library(data.table)
  library(patchwork)
  library(openxlsx)
  library(ggplot2)
  library(future)
})

## utils.R ====================================================================
source("R/utils.R")

## Configs ====================================================================
configs <- list.files(
  file.path("configs", "sims"), 
  pattern = "\\.json$", 
  recursive  = TRUE, 
  full.names = TRUE
)

## Parameters =================================================================
seed   <- 42
n_iter <- 20

sim_grid <- setDT(
    sim_make_grid(
    # Simulations
    n_samples = c(1000),
    n_degs    = c(500, 1000, 2000),
    log2fc    = c(1, 2, 3),
    # Imbalancing
    total_samples = c(150, 500, 1000),
    between_ratio = c(1, 5, 10),
    within_ratio  = c(1, 5, 10),
    # Other
    zero_rules = TRUE,
    verbose    = FALSE
  )
)

# sim_grid <- setDT(
#     sim_make_grid(
#     # Simulations
#     n_samples = c(20),
#     n_degs = c(500),
#     log2fc = c(1),
#     # Imbalancing
#     total_samples = c(10), # 150 is the smallest possible sample size with 1/10 imbalance
#     between_ratio = c(1),
#     within_ratio = c(1),
#     zero_rules = TRUE,
#     verbose = FALSE   
#   )
# )

## Loop (NB sim. -> imbalance -> DGE) =========================================
for (config_file in configs) {

  # Read config
  cfg  <- jsonlite::read_json(config_file)
  comp <- paste0(
    cfg$tech, "_", 
    cfg$g1, "_vs_",  
    cfg$g2, "_", 
    cfg$a1, "_vs_", 
    cfg$a2
  )

  # Output directory
  out_dir <- sub("^results", "results_major_review", cfg$out_dir)
  out_dir <- file.path(
    out_dir,
    paste0("sim_synthetic_nsubsets_", n_iter),
    comp
  )

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Output file
  result_file <- file.path(out_dir, "dge_res.rds")
  if (file.exists(result_file)) {
    message("\nSkipping: ", config_file, " — results already exist!")
    message("Results file: ", result_file)
    next  
  } else {
    message("\nProcessing: '", config_file, "'")
    message("Results file: ", result_file)
  }

  # Read data
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
    plot = FALSE, 
    verbose = FALSE
  )

  # Filter by expression
  data <- filter_by_expression(
    X = data$X$matr, 
    Y = data$Y$matr,
    MX = data$X$meta, 
    MY = data$Y$meta,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    any_group = FALSE,
    plot = FALSE, 
    verbose = FALSE
  )

  # Estimate NB params.
  NB_params_X <- estimate_nbinom_params(data$X$matr, plot = FALSE, verbose = FALSE)
  NB_params_Y <- estimate_nbinom_params(data$Y$matr, plot = FALSE, verbose = FALSE)

  # Split 'sim_grid' into sim params vs imbalance params
  sim_params_grid <- unique(sim_grid[, .(n_samples, n_degs, log2fc)])
  imb_params_grid <- unique(sim_grid[, .(total_samples, between_ratio, within_ratio)])

  # Outer loop: Simulate once per sim config (world population)
  res_sim_list <- list()
  for (j in seq_len(nrow(sim_params_grid))) {

    # Define 'world population' params
    n_samples <- sim_params_grid$n_samples[j]
    n_degs    <- sim_params_grid$n_degs[j]
    log2fc    <- sim_params_grid$log2fc[j]
    message(sprintf("\nsim. id: %d [n: %s degs: %s logfc: %s]", j, n_samples, n_degs, log2fc))

    # Simulate (NB dist.)
    sim <- sim_4group_expression(
      estimates_X = NB_params_X,
      estimates_Y = NB_params_Y,
      g_col = "condition",
      g_levels = c("Control", "Case"),
      a_col = "ancestry",
      a_levels = c("EUR_sim", "AFR_sim"),
      n_samples = n_samples,
      n_degs = n_degs,
      log2fc = log2fc,
      mean_method = "libnorm_mle",
      disp_method = "mle",
      drop_zeros = TRUE,
      plot = FALSE,
      seed = seed,
      verbose = FALSE
    )

    # Extract truth (TRUE interaction effects!!!!)
    cols_check <- c("coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "feature")
    sim_truth  <- setDT(sim$degs)
    truth_keep <- sim_truth[, c(.SD, .(T_truth = T_obs)), .SDcols = cols_check]
    if (anyNA(truth_keep$T_truth)) stop("Simulation logFC have NA values.")

    # Inner loop: Run all imbalance configs on same simulation
    for (k in seq_len(nrow(imb_params_grid))) {
      
      # Define imbalance params
      total_samples <- imb_params_grid$total_samples[k]
      between_ratio <- imb_params_grid$between_ratio[k]
      within_ratio  <- imb_params_grid$within_ratio[k]
      message(sprintf("\timb. id: %d [n: %s a_col: %s g_col: %s]", k, total_samples, between_ratio, within_ratio))

      # Imbalance 'world population' (biased sampling)
      imb <- sim_imbalanced_ancestry(
        X = sim$X$matr,
        Y = sim$Y$matr,
        MX = sim$X$meta,
        MY = sim$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        total_samples = total_samples,
        between_ratio = between_ratio,
        within_major_ratio = within_ratio,
        within_minor_ratio = within_ratio,
        seed = seed,
        replace = FALSE,
        verbose = FALSE
      )

      # Benchmark DGE methods
      dge_methods <- list()

      # Ancestry balanced
      message("\t\tlimma  (ancestry)")
      subset_res_ <- subset_limma_interaction_effect(
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        match = FALSE,
        covariates = NULL,
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
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        match = TRUE,
        covariates = NULL,
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
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        covariates = NULL,
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
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        covariates = NULL,
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
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        covariates = NULL,
        verbose = FALSE
      )

      deseq_full_res_$method <- "DESeq2-Wald"
      deseq_full_res_$strata <- "Imbalanced"

      # Check and add
      stopifnot(!anyNA(deseq_full_res_$T_obs))
      dge_methods[["DESeq2-Wald"]] <- setDT(deseq_full_res_)

      # Merge DEG res with truth
      for (m in names(dge_methods)) {

        dge_method <- dge_methods[[m]]

        # Merge
        merged <- merge(
          dge_method, 
          truth_keep, 
          by = cols_check, 
          all.x = TRUE, 
          sort = FALSE
        )

        # P-value NAs are expected as DESeq2 does internal filtering
        if (anyNA(merged$T_truth)) stop(sprintf("Method '%s' has NA values in 'T_truth'.", m))
        if (anyNA(merged$T_obs)) stop(sprintf("Method '%s' has NA values in 'T_obs'.", m))

        # Save per sim config
        res_sim_list[[length(res_sim_list) + 1]] <- data.table(
          n_samples = n_samples,
          n_degs = n_degs,
          log2fc = log2fc,
          total_samples = total_samples,
          between_ratio = between_ratio,
          within_ratio = within_ratio,
          method_label = m,
          merged
        )
      }
    }
  }

  # Simulation DGE res (combine)
  dge_res <- rbindlist(res_sim_list, use.names = TRUE, fill = TRUE)
  stopifnot(!anyNA(dge_res$T_truth))

  # Save
  saveRDS(dge_res, file = result_file)
}
