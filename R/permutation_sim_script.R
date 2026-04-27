## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
    library(data.table)
    library(jsonlite)
    library(future)
  }
)

## script fun. ================================================================
run_one <- function(
  X, 
  Y, 
  MX, 
  MY,
  seed,
  cfg, 
  study, 
  perm_id, 
  kind
) {
  methods <- list()
  
  # Limma subset
  limma_subset_res_ <- subset_limma_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    covariates = cfg$covariates,
    use_voom = TRUE, 
    n_iter = 10,
    method = "mean", 
    seed = seed,
    verbose = FALSE
  )
  limma_subset_res <- limma_subset_res_$methods_stats
  for (m in names(limma_subset_res)) {
    method_name <- paste0("subset-", m)
    methods[[method_name]] <- setDT(limma_subset_res[[m]])
  }
  
  # limma full
  limma_full_res_ <- limma_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    covariates = cfg$covariates,
    use_voom = TRUE, 
    verbose = FALSE
  )
  methods[["limma-voom"]] <- setDT(limma_full_res_)
  

  # edgeR full
  edgeR_full_res_ <- edgeR_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    covariates = cfg$covariates,
    verbose = FALSE
  )
  methods[["edgeR-QLFTest"]] <- setDT(edgeR_full_res_)
  

  # Deseq full
  deseq_full_res_ <- DESeq_interaction_effect(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY,
    g_col = cfg$g_col, 
    a_col = cfg$a_col,
    covariates = cfg$covariates,
    verbose = FALSE
  )
  methods[["DESeq-Wald"]] <- setDT(deseq_full_res_)


  # Combine results for all methods
  combined_res <- rbindlist(lapply(names(methods), function(method_name) {
    res <- methods[[method_name]]
    res[, `:=`(
      study   = study,
      perm_id = perm_id,
      kind    = kind,
      method  = method_name
    )]
    res[]
  }), fill = TRUE)
  
  # Return
  return(combined_res[])
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
n_perm <- 10
n_jobs <- length(configs)

# Parallel env.
future::plan(multisession, workers = n_jobs)

## Loop (orig. data -> perm. data) ============================================
for (config_file in configs) {

  # Read config
  cfg <- jsonlite::read_json(config_file)

  # Output directory
  out_dir <- file.path(
    cfg$out_dir, 
    "permutation_sim", 
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
  cat(" - run orig. data\n")
  orig_res <- run_one(
    X = data$X$matr, 
    Y = data$Y$matr, 
    MX = data$X$meta, 
    MY = data$Y$meta,
    cfg = cfg,
    seed = seed,
    study = cfg$study, 
    perm_id = NA, 
    kind = "real"
  )

  # Run perm. data (shuffled)
  cat(" - run perm. data\n")
  perm_res_list <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {

    # Set new seed
    seed_new <- seed + i
    set.seed(seed_new)
    cat(sprintf("   perm. id: %d (seed: %d)\n", i, seed_new))

    MX_perm <- data$X$meta
    MY_perm <- data$Y$meta
    
    # Permute 'g_col'
    MX_perm[[cfg$g_col]] <- sample(MX_perm[[cfg$g_col]])
    MY_perm[[cfg$g_col]] <- sample(MY_perm[[cfg$g_col]])
  
    perm_res_list[[i]] <- run_one(
      X = data$X$matr, 
      Y = data$Y$matr, 
      MX = MX_perm,
      MY = MY_perm,
      seed = seed_new,
      cfg = cfg, 
      study = cfg$study,
      perm_id = i,
      kind = "perm"
    )
  }

  # Combine orig. res & perm. res
  perm_res <- do.call(rbind, perm_res_list)
  res <- rbind(real_res, perm_res)

  # Save
  saveRDS(res, file = result_file)
}
