## max. code width ============================================================


## Libraries ==================================================================
suppressPackageStartupMessages({
  library(CrossAncestryGenPhen)
  library(data.table)
  library(jsonlite)
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
n_perm <- 20

# Max perm. retires
n_retries <- 10

# Parallel env.
future::plan(multisession, workers = 2)



## Loop (orig. data -> perm. data) ============================================
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
    paste0("sim_permutation_nsubsets_", n_iter),
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
    covariates = NULL,
    omit_na = TRUE,
    auto_cast = TRUE,
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
  
  # Run orig. data (unshuffled)
  message("run orig. data")
  message(sprintf("\torig. id: 1 (seed: %d)", seed))
  orig_res <- run_one(
    X = data$X$matr, 
    Y = data$Y$matr, 
    MX = data$X$meta, 
    MY = data$Y$meta,
    g_col = cfg$g_col,
    a_col = cfg$a_col,
    covariates = unlist(cfg$covariates),
    n_iter = n_iter,
    seed = seed,
    perm_id = NA, 
    kind = "real"
  )

  # Run perm. data (shuffled)
  message("run perm. data")
  perm_res_list <- vector("list", n_perm)

  success_count  <- 0
  attempt_count  <- 0
  failed_retries <- 0
  skip_config    <- FALSE

  while (success_count < n_perm) {
    attempt_count <- attempt_count + 1

    # Set new seed
    seed_new <- seed + attempt_count
    set.seed(seed_new)

    MX_perm <- data$X$meta
    MY_perm <- data$Y$meta
    
    # Permute
    MX_perm[[cfg$g_col]] <- sample(MX_perm[[cfg$g_col]])
    MY_perm[[cfg$g_col]] <- sample(MY_perm[[cfg$g_col]])
  
    possible_error <- tryCatch({
      res_one <- run_one(
        X = data$X$matr, 
        Y = data$Y$matr, 
        MX = MX_perm,
        MY = MY_perm,
        g_col = cfg$g_col,
        a_col = cfg$a_col, 
        covariates = unlist(cfg$covariates),
        n_iter = n_iter,
        seed = seed_new,
        perm_id = success_count + 1,
        kind = "perm"
      )

      # Update counter
      success_count  <- success_count + 1
      failed_retries <- 0

      # Add 
      message(sprintf("\tperm. id: %d (success with seed: %d after %d attempts)", success_count, seed_new, attempt_count))
      perm_res_list[[success_count]] <- res_one

    }, error = function(e) {
      failed_retries <<- failed_retries + 1
      message(sprintf("\tAttempt %d (seed: %d) failed due to collinearity (failed retries row: %d/%d). Trying next...", attempt_count, seed_new, failed_retries, n_retries))
      
      # Skript break after n retries
      if (failed_retries >= n_retries) {
        message(sprintf("\tWarning: Skipping '%s' because %d consecutive permutations failed.", config_file, n_retries))
        skip_config <<- TRUE
      }
      return(NULL) 
    })

    if (skip_config) break
  }

  if (skip_config) next

  # Combine orig. res & perm. res
  perm_res <- do.call(rbind, perm_res_list)
  res <- rbind(orig_res, perm_res)

  # Save
  saveRDS(res, file = result_file)
}
