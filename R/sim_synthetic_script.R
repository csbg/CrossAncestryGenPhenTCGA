## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
    library(data.table)
    library(patchwork)
    library(ggplot2)
    library(future)
  }
)

## script fun. ================================================================
check_alignment <- function(dt_method, dt_truth, cols_check) {

  # Check same number of rows
  if (nrow(dt_method) != nrow(dt_truth)) {
    stop("Row count mismatch: method=", nrow(dt_method), " vs truth=", nrow(dt_truth))
  }

  # Check columns exist
  missing_cols <- setdiff(cols_check, intersect(names(dt_method), names(dt_truth)))
  if (length(missing_cols) > 0) {
    stop("Missing key columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check identical ordering + values
  equal_keys <- dt_method[, ..cols_check] == dt_truth[, ..cols_check]
  if (!all(equal_keys)) {
    bad_rows <- which(!apply(equal_keys, 1, all))
    example  <- head(bad_rows)
    stop("Data misaligned in ", length(bad_rows), " rows. Example indices: ", paste(example, collapse = ", "))
  }

  # Optional hash check for speed and safety
  hash_m <- digest::digest(dt_method[, ..cols_check])
  hash_t <- digest::digest(dt_truth[, ..cols_check])
  if (hash_m != hash_t) stop("Hash mismatch: possible ordering issue")

  # Row difference
  anti_miss1 <- data.table::fsetdiff(dt_method[, ..cols_check], dt_truth[, ..cols_check])
  anti_miss2 <- data.table::fsetdiff(dt_truth[, ..cols_check], dt_method[, ..cols_check])
  if (nrow(anti_miss1) > 0) {
    message("→ Rows in method but not in truth: ", nrow(anti_miss1))
    stop("Sim truth and DGE results misaligned — see printed diagnostics for details.")
    }
  if (nrow(anti_miss2) > 0) {
    message("→ Rows in truth but not in method: ", nrow(anti_miss2))
    stop("Sim truth and DGE results misaligned — see printed diagnostics for details.")
  }
  
  # Return
  return(TRUE)
}

## Configs ====================================================================
configs <- list.files(
  file.path("configs", "sims"), 
  pattern = "\\.json$", 
  recursive  = TRUE, 
  full.names = TRUE
)

## Parameters =================================================================
seed   <- 42
n_jobs <- length(configs)

# Parallel env.
future::plan(multisession, workers = n_jobs)

# Simulation configs 
# sim_grid <- setDT(
#     sim_make_grid(
#     # Simulations
#     n_samples = c(1000),
#     n_degs = c(0, 500, 1000),
#     log2fc = c(0, 1, 2),
#     # Imbalancing
#     total_samples = c(100, 500, 1000),
#     between_ratio = c(1, 5, 10),
#     within_ratio = c(1, 5, 10),
#     zero_rules = TRUE,  # Only run one negative control simulation (don't expand grid where n_degs = 0 to more then one row)
#     verbose = TRUE   
#   )
# )

sim_grid <- setDT(
    sim_make_grid(
    # Simulations
    n_samples = c(50),
    n_degs = c(500),
    log2fc = c(1),
    # Imbalancing
    total_samples = c(40),
    between_ratio = c(1),
    within_ratio = c(1),
    zero_rules = TRUE,
    verbose = FALSE   
  )
)

## Loop (NB sim. -> imbalance -> DGE) =========================================
for (config_file in configs) {

  # Read config
  cfg <- jsonlite::read_json(config_file)

  # Output directory
  out_dir <- file.path(
    cfg$out_dir, 
    "synthetic_sim", 
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
    cat("Skipping: '", config_file, "' — results already exist!\n")
    cat("Results file: ", result_file, "\n")
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
  NB_params_X <- estimate_nbinom_params(
    data$X$matr, 
    plot = FALSE,
    verbose = FALSE
  )

  NB_params_Y <- estimate_nbinom_params(
    data$Y$matr, 
    plot = FALSE,
    verbose = FALSE
  )

  # QC plot
  ggsaveDK(
    plot = NB_params_X$plot / NB_params_Y$plot,
    file = file.path(out_dir, "nbinom_params.svg"),
    height = 5,
    width = 10,
    units = "cm"
  )

  # Split 'sim_grid' into sim params vs imbalance params
  sim_params_grid <- unique(sim_grid[, .(n_samples, n_degs, log2fc)])
  imb_params_grid <- unique(sim_grid[, .(total_samples, between_ratio, within_ratio)])

  # Outer loop: Simulate once per sim config (world population)
  res_sim_list <- list()
  qc_sim_list  <- list()
  for (j in seq_len(nrow(sim_params_grid))) {

    # Define 'world population' params
    cat(sprintf(" - sim. id: %d\n", j))
    n_samples <- sim_params_grid$n_samples[j]
    n_degs    <- sim_params_grid$n_degs[j]
    log2fc    <- sim_params_grid$log2fc[j]

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

    # QC plot 
    qc_data  <- rbind(sim$X$matr, sim$Y$matr)
    qc_title <- sprintf( "N: %d nDEG: %d log2fc: %.1f", n_samples, n_degs, log2fc)
    qc_theme <- theme(plot.title = element_text(size = 6, hjust = 0.5))

    qc_sim_list[[j]] <- wrap_plots(
      sim$plot /
        (
          plot_mean_variance_trend(qc_data, point_size = 0.5) +
          plot_mean_variance_density(qc_data, point_size = 0.5)
        )
      ) +
      plot_annotation(
        title = qc_title,
        theme = qc_theme
      )

    # Extract truth (TRUE interaction effects!!!!)
    cols_check <- c("coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "feature")
    sim_truth  <- setDT(sim$degs)
    truth_keep <- sim_truth[, c(.SD, .(T_truth = T_obs)), .SDcols = cols_check]

    # Inner loop: Run all imbalance configs on same simulation
    for (k in seq_len(nrow(imb_params_grid))) {
      
      # Define imbalance params
      cat(sprintf("    imb. id: %d\n", k))
      total_samples <- imb_params_grid$total_samples[k]
      between_ratio <- imb_params_grid$between_ratio[k]
      within_ratio  <- imb_params_grid$within_ratio[k]

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

      # Subset-limma
      subset_res_ <- subset_limma_interaction_effect(
        X = imb$X$matr,
        Y = imb$Y$matr,
        MX = imb$X$meta,
        MY = imb$Y$meta,
        g_col = "condition",
        a_col = "ancestry",
        covariates = NULL,
        use_voom = TRUE,
        n_iter = 10,
        method = "mean",
        seed = seed,
        verbose = FALSE
      )

      # All p-value agg. methods
      subset_res <- subset_res_$methods_stats
      for (m in names(subset_res)) {
        m_name <- paste0("subset-", m)
        method_res  <- setDT(subset_res[[m]])
        check_alignment(method_res, sim_truth, cols_check)
        dge_methods[[m_name]] <- method_res
      }

      # Limma full
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
      check_alignment(setDT(limma_full_res_), sim_truth, cols_check)
      dge_methods[["limma-voom"]] <- setDT(limma_full_res_)

      # EdgeR full
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
      check_alignment(setDT(edgeR_full_res_), sim_truth, cols_check)
      dge_methods[["edgeR-QLFTest"]] <- setDT(edgeR_full_res_)

      # DESeq full
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
      check_alignment(setDT(deseq_full_res_), sim_truth, cols_check)
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

        # Save per sim config
        res_sim_list[[length(res_sim_list) + 1]] <- data.table(
          study = cfg$study,
          method = m,
          n_samples = n_samples,
          n_degs = n_degs,
          log2fc = log2fc,
          total_samples = total_samples,
          between_ratio = between_ratio,
          within_ratio = within_ratio,
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

  # Simulation QC plots (combine)
  ggsaveDK(
    plot = wrap_plots(
      lapply(qc_sim_effects_list, wrap_elements),
      ncol = nrow(sim_params_grid)
    ),
    file = file.path(out_dir, "sim_effects.svg"),
    height = 8 * nrow(sim_params_grid),
    width = 8 * nrow(sim_params_grid),
    units = "cm"
  )
}
