# ==============================================================================
# SPECTRAL NORMALIZATION MODULE
# ==============================================================================
# Provides normalization methods for NMR spectral data:
# - CSN (Constant Sum Normalization): normalizes to total area = 1
# - PQN (Probabilistic Quotient Normalization): corrects dilution effects
#
# NOTE: Pareto scaling is NOT included because it's useless for STOCSY/RANSY.
# Pearson correlation (used in STOCSY) already centers and normalizes internally,
# so r(x_raw, y_raw) = r(x_pareto, y_pareto). Pareto is only useful for PCA/OPLS-DA
# where covariance (not correlation) is computed.
#
# Reference: Dieterle et al. (2006) Anal Chem 78:4281-4290 (PQN)
# ==============================================================================

# ------------------------------------------------------------------------------
# CSN - Constant Sum Normalization
# ------------------------------------------------------------------------------
# Divides each sample (column) by its total sum, so each sample has total area = 1
# This corrects for differences in total concentration between samples.
#
# @param mat Numeric matrix (rows = ppm points, columns = samples)
# @return Normalized matrix with same dimensions
# ------------------------------------------------------------------------------
normalize_csn <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("normalize_csn: input must be a numeric matrix")
  }

  # Calculate sum for each column (sample)
  col_sums <- colSums(mat, na.rm = TRUE)

  # Avoid division by zero
  col_sums[col_sums == 0 | is.na(col_sums)] <- 1

  # Divide each column by its sum
  mat_norm <- sweep(mat, 2, col_sums, FUN = "/")

  return(mat_norm)
}

# ------------------------------------------------------------------------------
# PQN - Probabilistic Quotient Normalization
# ------------------------------------------------------------------------------
# Standard procedure:
# 1. Apply CSN as pre-normalization
# 2. Calculate reference spectrum (median across samples at each ppm)
# 3. For each sample, calculate quotient = sample/reference at each ppm
# 4. Median of quotients = dilution factor for that sample
# 5. Divide ORIGINAL matrix by these dilution factors
#
# @param mat Numeric matrix (rows = ppm points, columns = samples)
# @param min_valid_points Minimum number of valid quotients required (default 10)
# @return Normalized matrix with same dimensions
# ------------------------------------------------------------------------------
normalize_pqn <- function(mat, min_valid_points = 10) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("normalize_pqn: input must be a numeric matrix")
  }

  n_samples <- ncol(mat)
  n_points <- nrow(mat)

  # Step 1: Apply CSN as pre-normalization (standard PQN procedure)
  mat_csn <- normalize_csn(mat)

  # Step 2: Calculate reference spectrum (median at each ppm point)
  reference <- apply(mat_csn, 1, median, na.rm = TRUE)

  # Step 3-4: Calculate dilution factors for each sample
  dilution_factors <- numeric(n_samples)

  for (j in seq_len(n_samples)) {
    sample_vals <- mat_csn[, j]

    # Only calculate quotient where both sample and reference > 0
    valid_idx <- which(sample_vals > 0 & reference > 0 &
                        !is.na(sample_vals) & !is.na(reference))

    if (length(valid_idx) >= min_valid_points) {
      quotients <- sample_vals[valid_idx] / reference[valid_idx]
      dilution_factors[j] <- median(quotients, na.rm = TRUE)
    } else {
      # Not enough valid points - use factor of 1 (no correction)
      dilution_factors[j] <- 1
    }
  }

  # Handle zero or NA factors
  dilution_factors[dilution_factors == 0 | is.na(dilution_factors)] <- 1

  # Step 5: Divide ORIGINAL matrix by dilution factors
  mat_pqn <- sweep(mat, 2, dilution_factors, FUN = "/")

  return(mat_pqn)
}

# ------------------------------------------------------------------------------
# Main orchestration function
# ------------------------------------------------------------------------------
# Applies normalization to correct dilution/concentration effects between samples.
# PQN and CSN are mutually exclusive - if both TRUE, PQN is used.
#
# @param mat Numeric matrix (rows = ppm points, columns = samples)
# @param use_pqn Boolean - apply PQN normalization
# @param use_csn Boolean - apply CSN normalization
# @param log_func Optional logging function (takes single string argument)
# @return List with:
#   - matrix: the transformed matrix
#   - methods_applied: character vector of methods applied
# ------------------------------------------------------------------------------
apply_normalization <- function(mat,
                                 use_pqn = FALSE,
                                 use_csn = FALSE,
                                 log_func = NULL) {

  # Helper for logging
  log_msg <- function(msg) {
    if (!is.null(log_func) && is.function(log_func)) {
      log_func(msg)
    }
  }

  methods_applied <- character(0)
  result_mat <- mat

  # Validate input
  if (!is.matrix(mat) || !is.numeric(mat)) {
    log_msg("[NORMALIZATION] Warning: input is not a numeric matrix, skipping normalization")
    return(list(
      matrix = mat,
      methods_applied = "none (invalid input)"
    ))
  }

  # Handle mutually exclusive PQN/CSN
  if (use_pqn && use_csn) {
    log_msg("[NORMALIZATION] Warning: both PQN and CSN requested - using PQN (they are mutually exclusive)")
    use_csn <- FALSE
  }

  # Apply normalization (PQN or CSN)
  if (use_pqn) {
    log_msg("[NORMALIZATION] Applying PQN (Probabilistic Quotient Normalization)...")
    result_mat <- normalize_pqn(result_mat)
    methods_applied <- c(methods_applied, "PQN")
    log_msg("[NORMALIZATION] PQN applied successfully")
  } else if (use_csn) {
    log_msg("[NORMALIZATION] Applying CSN (Constant Sum Normalization)...")
    result_mat <- normalize_csn(result_mat)
    methods_applied <- c(methods_applied, "CSN")
    log_msg("[NORMALIZATION] CSN applied successfully")
  }

  # Summary
  if (length(methods_applied) == 0) {
    log_msg("[NORMALIZATION] No normalization applied (raw intensities)")
    methods_applied <- "none"
  } else {
    log_msg(paste0("[NORMALIZATION] Methods applied: ", paste(methods_applied, collapse = " + ")))
  }

  return(list(
    matrix = result_mat,
    methods_applied = paste(methods_applied, collapse = " + ")
  ))
}
