## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(library(CrossAncestryGenPhen))

## Toy data ===================================================================
set.seed(42)

# Example data 
n_features <- 100
n_EUR <- 600 
n_AFR <- 40  

id_EUR <- paste0("EUR_", seq_len(n_EUR))
id_AFR <- paste0("AFR_", seq_len(n_AFR))

lambda_features <- sample(
  c(
    rpois(20, 3),    # ~20 low-count features
    rpois(60, 20),   # ~60 medium-count features
    rpois(20, 120)   # ~20 high-count features
  ),
  size = n_features,
  replace = FALSE
)

# Expression matrices for EUR and AFR ancestries
X <- sapply(lambda_features, function(lam) rpois(n_EUR, lambda = lam))
Y <- sapply(lambda_features, function(lam) rpois(n_AFR, lambda = lam))
colnames(X) <- colnames(Y) <- paste0("Feature_", seq_len(n_features))


# Metadata for EUR and AFR ancestries
# EUR: overrepresented compared to AFR
MX <- data.frame(
  condition = factor(c(rep("Control", 400), rep("Case", 200)), levels = c("Control", "Case")),
  ancestry = "EUR",
  sex      = sample(c("Male", "Female"), n_EUR, replace = TRUE),
  age      = round(rnorm(n_EUR, mean = 50, sd = 12))
)

# AFR: underrepresented compared to EUR
MY <- data.frame(
  condition = factor(c(rep("Control", 10), rep("Case", 30)), levels = c("Control", "Case")),
  ancestry = "AFR",
  sex      = sample(c("Male", "Female"), n_AFR, replace = TRUE),
  age      = round(rnorm(n_AFR, mean = 48, sd = 11))
)

# Rownames of matrix must be smaple ids
rownames(X) <- rownames(MX) <- id_EUR
rownames(Y) <- rownames(MY) <- id_AFR


# Spike in effect in AFR Case (strong effect, but only for half of AFR cases)
# AFR Case: strong effect, half of cases
afr_case_idx  <- which(MY$condition == "Case")
afr_spike_idx <- sample(afr_case_idx, length(afr_case_idx) / 2)
Y[afr_spike_idx, 1:4] <- Y[afr_spike_idx, 1:4] + rpois(length(afr_spike_idx) * 4, lambda = 80)

# EUR Case: weaker effect, 20% of cases
eur_case_idx  <- which(MX$condition == "Case")
eur_spike_idx <- sample(eur_case_idx, length(eur_case_idx) * 0.2)
X[eur_spike_idx, 1:4] <- X[eur_spike_idx, 1:4] + rpois(length(eur_spike_idx) * 4, lambda = 20)

## Splitting ==================================================================
split_stratified_ancestry_sets <- function(
  X, 
  Y, 
  MX, 
  MY,
  g_col, 
  a_col,
  match_mutual = FALSE,
  seed = NULL,
  verbose = TRUE
) {
    
  ## --- Seed ---
  if (!is.null(seed)) set.seed(seed)

  ## --- Input checks ---
  assert_input(
    X = X, 
    Y = Y,
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col
  )

  ## --- Factor setup ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  if (length(g_levels) != 2 || length(unique(c(a_1, a_2))) != 2) {
    stop("[split_stratified_ancestry_sets] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }

  # Vektoren direkt aus den Dataframes extrahieren (vermeidet wiederholten Spaltenzugriff)
  vec_g_X <- MX[[g_col]]
  vec_g_Y <- MY[[g_col]]

  ## --- Target Count Calculation ---
  count_X <- table(vec_g_X)
  count_Y <- table(vec_g_Y)

  if (match_mutual) {
    min_overall   <- min(c(count_X, count_Y))
    target_counts <- setNames(rep(min_overall, length(count_Y)), names(count_Y))
    message("Enforcing 'a_col x g_col' balance.")

  } else {
    target_counts <- count_Y
    message("Enforcing 'a_col' balance.")
  }

  ## --- Feasibility check ---
  strata_names <- names(target_counts)
  insufficient <- strata_names[target_counts[strata_names] > count_X[strata_names]]
  missing      <- setdiff(strata_names, names(count_X))

  if (length(missing) > 0 || length(insufficient) > 0) {
    stop("[split_stratified_ancestry_sets] X cannot fulfill the requested distribution layout.\n",
         "Missing strata: ", paste(missing, collapse = ", "), "\n",
         "Insufficient strata: ", paste(insufficient, collapse = ", "))
  }

  ## --- Process Y (Inference & Remaining RY) ---
  ids_Y <- rownames(Y)
  
  if (match_mutual) {
    sampled_ids_Y <- vector("list", length(strata_names))
    for (i in seq_along(strata_names)) {
      stratum <- strata_names[i]
      idx <- which(vec_g_Y == stratum)
      sampled_ids_Y[[i]] <- ids_Y[sample(idx, size = target_counts[stratum], replace = FALSE)]
    }
    sampled_ids_Y <- unlist(sampled_ids_Y, use.names = FALSE)
    
    mask_Y_subset <- ids_Y %in% sampled_ids_Y
    
    Y_matr  <- Y[mask_Y_subset, , drop = FALSE]
    Y_meta  <- MY[mask_Y_subset, , drop = FALSE]
    RY_matr <- Y[!mask_Y_subset, , drop = FALSE]
    RY_meta <- MY[!mask_Y_subset, , drop = FALSE]
  } else {
    Y_matr  <- Y
    Y_meta  <- MY
    RY_matr <- NULL
    RY_meta <- NULL
  }

  ## --- Process X (Subset X & Remaining RX) ---
  ids_X <- rownames(X)
  sampled_ids_X <- vector("list", length(strata_names))
  for (i in seq_along(strata_names)) {
    stratum <- strata_names[i]
    idx <- which(vec_g_X == stratum)
    sampled_ids_X[[i]] <- ids_X[sample(idx, size = target_counts[stratum], replace = FALSE)]
  }
  sampled_ids_X <- unlist(sampled_ids_X, use.names = FALSE)

  mask_X_subset <- ids_X %in% sampled_ids_X

  X_matr  <- X[mask_X_subset, , drop = FALSE]
  X_meta  <- MX[mask_X_subset, , drop = FALSE]
  RX_matr <- X[!mask_X_subset, , drop = FALSE]
  RX_meta <- MX[!mask_X_subset, , drop = FALSE]

  ## --- Verbose summary ---
  if (verbose) {
    fmt_counts <- function(M_sub, g_col) {
      if (is.null(M_sub) || nrow(M_sub) == 0) return("N/A")
      tab <- table(M_sub[[g_col]])
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    message("\nStratified split:")
    message(sprintf("%-20s  N: %-4d %s features: %-4d", paste0("Remaining RX (", a_1, "):"), nrow(RX_matr), fmt_counts(RX_meta, g_col), ncol(RX_matr)))
    if (match_mutual) {message(sprintf("%-20s  N: %-4d %s features: %-4d", paste0("Remaining RY (", a_2, "):"), nrow(RY_matr), fmt_counts(RY_meta, g_col), ncol(RY_matr)))}
    message(sprintf("%-20s  N: %-4d %s features: %-4d", paste0("Subset    X  (", a_1, "):"), nrow(X_matr), fmt_counts(X_meta, g_col), ncol(X_matr)))
    if (match_mutual) {message(sprintf("%-20s  N: %-4d %s features: %-4d", paste0("Subset    Y  (", a_2, "):"), nrow(Y_matr), fmt_counts(Y_meta, g_col), ncol(Y_matr)))} else {
      message(sprintf("%-20s  N: %-4d %s features: %-4d", paste0("Inference Y  (", a_2, "):"), nrow(Y_matr), fmt_counts(Y_meta, g_col), ncol(Y_matr)))
    }
  }

  ## --- Return ---
  return(
    list(
      RX = list(matr = RX_matr, meta = RX_meta, ids = rownames(RX_matr)),
      RY = list(matr = RY_matr, meta = RY_meta, ids = if(!is.null(RY_matr)) rownames(RY_matr) else NULL),
      X  = list(matr = X_matr,  meta = X_meta,  ids = rownames(X_matr)),
      Y  = list(matr = Y_matr,  meta = Y_meta,  ids = rownames(Y_matr)),
      strata_info = list(usable = strata_names, missing = missing, insufficient = insufficient)
    )
  )
}


split <- split_stratified_ancestry_sets(
  X = X,
  Y = Y,
  MX = MX,
  MY = MY,
  g_col = "condition",
  a_col = "ancestry",
  seed = 42,
  match_mutual = FALSE,
  verbose = TRUE
)

plot_stratified_sets(
  MX = split$X$meta,
  MY = split$Y$meta,
  MR = split$RX$meta,
  x_var = "ancestry",
  fill_var = "condition",
  title = "Stratified ancestry sets",
  x_label = "Ancestry",
  y_label = "Nr. of patients"
)
