# Spectral normalization utilities

# Constant Sum Normalization
normalize_csn <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("normalize_csn: input must be a numeric matrix")
  }

  # Calculate column sums
  col_sums <- colSums(mat, na.rm = TRUE)

  # Avoid division by zero
  col_sums[col_sums == 0 | is.na(col_sums)] <- 1

  # Normalize each column
  mat_norm <- sweep(mat, 2, col_sums, FUN = "/")

  return(mat_norm)
}

# Probabilistic Quotient Normalization
normalize_pqn <- function(mat, min_valid_points = 10) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("normalize_pqn: input must be a numeric matrix")
  }

  n_samples <- ncol(mat)
  n_points <- nrow(mat)

  # Apply CSN pre-normalization
  mat_csn <- normalize_csn(mat)

  # Calculate the reference spectrum
  reference <- apply(mat_csn, 1, median, na.rm = TRUE)

  # Calculate dilution factors
  dilution_factors <- numeric(n_samples)

  for (j in seq_len(n_samples)) {
    sample_vals <- mat_csn[, j]

    # Keep valid sample/reference pairs
    valid_idx <- which(sample_vals > 0 & reference > 0 &
                        !is.na(sample_vals) & !is.na(reference))

    if (length(valid_idx) >= min_valid_points) {
      quotients <- sample_vals[valid_idx] / reference[valid_idx]
      dilution_factors[j] <- median(quotients, na.rm = TRUE)
    } else {
      dilution_factors[j] <- 1
    }
  }

  # Replace invalid dilution factors
  dilution_factors[dilution_factors == 0 | is.na(dilution_factors)] <- 1

  # Normalize the original matrix
  mat_pqn <- sweep(mat, 2, dilution_factors, FUN = "/")

  return(mat_pqn)
}

# Normalization orchestration
apply_normalization <- function(mat,
                                 use_pqn = FALSE,
                                 use_csn = FALSE,
                                 log_func = NULL) {

  # Logging helper
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

  # Enforce exclusive normalization mode
  if (use_pqn && use_csn) {
    log_msg("[NORMALIZATION] Warning: both PQN and CSN requested - using PQN (they are mutually exclusive)")
    use_csn <- FALSE
  }

  # Apply normalization
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

  # Log summary
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
