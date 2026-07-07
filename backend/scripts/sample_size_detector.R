# =====================================================================
# sample_size_detector.R — Adaptive annotation mode based on sample count
# =====================================================================
# Detects the number of samples in spectrum matrix and determines:
# - n >= 5: STOCSY + RANSY + Intensity ratios (full mode)
# - 2 <= n < 5: Intensity ratios only + warning (degraded mode)
# - n = 1: Intensity ratios only (single sample mode)
#
# Usage:
#   source("sample_size_detector.R")
#   mode_info <- detect_annotation_mode(spectrum_matrix)
#   if (!mode_info$use_stocsy) {
#     log_message(mode_info$warning)
#     # Skip STOCSY/RANSY, use intensity_ratio_annotation()
#   }
# =====================================================================

#' Detect annotation mode based on sample count
#'
#' @param spectrum_matrix data.table with 'ppm' column + intensity columns
#' @return list(n_samples, use_stocsy, use_ransy, mode, warning)
detect_annotation_mode <- function(spectrum_matrix) {
  intensity_cols <- setdiff(names(spectrum_matrix), c("ppm", "ppm_matrix", "ppm_F1", "ppm_F2"))
  n_samples <- length(intensity_cols)

  mode_info <- list(
    n_samples = n_samples,
    use_stocsy = FALSE,
    use_ransy = FALSE,
    mode = "unknown",
    warning = "",
    min_samples_stocsy = 5,
    min_samples_ransy = 5
  )

  if (n_samples >= 5) {
    mode_info$use_stocsy <- TRUE
    mode_info$use_ransy <- TRUE
    mode_info$mode <- "full"
    mode_info$warning <- sprintf("Full annotation mode (n=%d samples): STOCSY + RANSY + Intensity ratios", n_samples)
  } else if (n_samples >= 2) {
    mode_info$use_stocsy <- FALSE
    mode_info$use_ransy <- FALSE
    mode_info$mode <- "ratio_only"
    mode_info$warning <- sprintf(
      "WARNING: Only %d sample(s) detected (minimum 5 required for STOCSY/RANSY). Using INTENSITY RATIOS ONLY. Results may be less reliable.",
      n_samples
    )
  } else if (n_samples == 1) {
    mode_info$use_stocsy <- FALSE
    mode_info$use_ransy <- FALSE
    mode_info$mode <- "single_sample"
    mode_info$warning <- "SINGLE SAMPLE MODE: Using intensity ratios only. No statistical validation possible."
  } else {
    mode_info$mode <- "error"
    mode_info$warning <- "ERROR: No intensity columns found in spectrum matrix."
  }

  return(mode_info)
}


#' Intensity ratio-based annotation (fallback for n < 5)
#' This function mimics the EXACT SAME logic as mode normal, but skips STOCSY/RANSY calculations
#'
#' @param matched_ppms character/numeric vector of matched ppm values from spectrum
#' @param molecule_ppm character/numeric vector OR string with "|" separator
#' @param tol_match tolerance for matching
#' @param ratio_min minimum intensity ratio threshold (default 0.01)
#' @param spectrum_matrix data.table with ppm + intensity columns
#' @param spectrum_mat matrix of intensities (for optimization)
#' @return list of annotation results (SAME STRUCTURE as mode normal)
intensity_ratio_annotation <- function(matched_ppms, molecule_ppm, tol_match = 0.1,
                                       ratio_min = 0.01, spectrum_matrix, spectrum_mat) {
  if (length(matched_ppms) == 0 || length(molecule_ppm) == 0) {
    return(NULL)
  }

  # Detect 1D vs 2D spectrum
  is_2D <- "ppm_matrix" %in% names(spectrum_matrix)
  ppm_col <- if (is_2D) spectrum_matrix$ppm_matrix else spectrum_matrix$ppm

  # Step 1: Calculate average intensities and sort by descending intensity
  idx_match <- match(matched_ppms, ppm_col)
  idx_match <- idx_match[!is.na(idx_match)]
  if (!length(idx_match)) {
    return(NULL)
  }

  avg_int <- rowMeans(spectrum_mat[idx_match, , drop = FALSE], na.rm = TRUE)
  matched_ppms_clean <- matched_ppms[!is.na(match(matched_ppms, ppm_col))]

  # Sort by descending intensity
  ord_desc <- order(avg_int, decreasing = TRUE)
  ppms_desc <- matched_ppms_clean[ord_desc]
  int_desc <- avg_int[ord_desc]

  # Step 2: Select pivot (most intense signal)
  if (length(ppms_desc) == 0) {
    return(NULL)
  }
  pivot <- ppms_desc[1]
  pivot_int <- int_desc[1]

  if (is.na(pivot_int) || pivot_int == 0) {
    return(NULL)
  }

  # Step 3: Filter by intensity ratio (SAME AS MODE NORMAL)
  ratios <- int_desc / pivot_int
  idx_valide <- which(int_desc <= pivot_int & ratios >= ratio_min & !is.na(int_desc))

  if (!length(idx_valide)) {
    return(NULL)
  }

  sub_chain_vec <- sort(unique(ppms_desc[idx_valide]))

  # Step 4: Match sub_chain to molecule (EXACT SAME LOGIC AS MODE NORMAL)
  # Detect if 2D by checking column types in molecule_ppm input
  is_2D_mol <- any(grepl("\\|", as.character(molecule_ppm)))

  pairs_dt <- NULL

  if (is_2D_mol) {
    # === 2D ANNOTATION (HSQC/COSY) ===
    # Parse molecule_ppm if it's in "F1|F2" format
    mol_parts <- strsplit(as.character(molecule_ppm), "\\|")
    unique_mppm <- unique(data.table(
      shift_value_F1 = sapply(mol_parts, function(x) if (length(x) >= 1) as.numeric(x[1]) else NA_real_),
      shift_value_F2 = sapply(mol_parts, function(x) if (length(x) >= 2) as.numeric(x[2]) else NA_real_)
    ))
    unique_mppm <- unique_mppm[!is.na(shift_value_F1) & !is.na(shift_value_F2)]

    if (!nrow(unique_mppm)) {
      return(NULL)
    }

    # Merge with spectrum_matrix to get detailed info for sub_chain
    avail <- spectrum_matrix[ppm_matrix %in% sub_chain_vec, .(ppm_matrix, ppm_F1, ppm_F2)]
    if (!nrow(avail)) {
      return(NULL)
    }

    # RANSY pairings (1-1) - SAME CODE AS MODE NORMAL but without CV cost
    pairs_dt <- data.table(mol_F1 = numeric(0), mol_F2 = numeric(0), ppm_ransy = character(0))

    # Tolerance handling (for 2D, tol_match might be pre-calculated Euclidean)
    tol_F1 <- tol_match
    tol_F2 <- tol_match

    for (i in seq_len(nrow(unique_mppm))) {
      m1 <- unique_mppm$shift_value_F1[i]
      m2 <- unique_mppm$shift_value_F2[i]
      cand <- avail[abs(ppm_F1 - m1) <= tol_F1 & abs(ppm_F2 - m2) <= tol_F2]

      if (nrow(cand) > 0) {
        # Distance-based matching only (no CV cost in degraded mode)
        cand[, dist_norm := sqrt(((ppm_F1 - m1) / tol_F1)^2 + ((ppm_F2 - m2) / tol_F2)^2)]
        best_r <- cand[which.min(dist_norm), ppm_matrix]
        pairs_dt <- rbind(pairs_dt, data.table(mol_F1 = m1, mol_F2 = m2, ppm_ransy = best_r))
        avail <- avail[ppm_matrix != best_r]
      }
    }

    matched_unique <- unique(pairs_dt$ppm_ransy)
    coverage_ratio <- if (nrow(unique_mppm) > 0) length(matched_unique) / nrow(unique_mppm) * 100 else 0
  } else {
    # === 1D ANNOTATION (1H/13C) ===
    unique_mppm <- sort(unique(as.numeric(molecule_ppm)))
    avail <- spectrum_matrix[ppm %in% sub_chain_vec, .(ppm)]
    if (!nrow(avail)) {
      return(NULL)
    }

    pairs_dt <- data.table(ppm_molecule = numeric(0), ppm_ransy = numeric(0))

    for (u in unique_mppm) {
      cand <- avail[abs(ppm - u) <= tol_match]
      if (nrow(cand) > 0) {
        best_r <- cand[which.min(abs(ppm - u)), ppm]
        pairs_dt <- rbind(pairs_dt, data.table(ppm_molecule = u, ppm_ransy = best_r))
        avail <- avail[ppm != best_r]
      }
    }

    coverage_ratio <- if (length(unique_mppm) > 0) nrow(pairs_dt) / length(unique_mppm) * 100 else 0
  }

  if (!is.finite(coverage_ratio)) coverage_ratio <- 0

  # Step 5: Calculate intensity ratios string (SAME FORMAT AS MODE NORMAL)
  idx_sub <- match(sub_chain_vec, ppm_col)
  idx_sub_valid <- which(!is.na(idx_sub))
  if (!length(idx_sub_valid)) {
    return(NULL)
  }

  sub_chain_clean <- sub_chain_vec[idx_sub_valid]
  avg_mat_all <- rowMeans(spectrum_mat[idx_sub[idx_sub_valid], , drop = FALSE], na.rm = TRUE)
  max_int_all <- max(avg_mat_all, na.rm = TRUE)
  if (!is.finite(max_int_all) || max_int_all == 0) max_int_all <- 1

  # Separator depends on 1D vs 2D
  sep <- if (is_2D) "|" else ";"
  intensity_ratios_str <- paste(round(avg_mat_all / max_int_all, 3), collapse = sep)

  # Return result in EXACT SAME FORMAT AS MODE NORMAL
  return(list(list(
    base_ppm          = pivot,
    sub_chain         = sub_chain_clean,
    coverage_sub      = coverage_ratio,
    intensity_ratios  = intensity_ratios_str,
    pivot_ppm         = pivot,
    RANSY_signals_all = sub_chain_clean, # All signals passing ratio filter = RANSY without CV
    matched_for_cov   = if (is_2D_mol) unique(pairs_dt$ppm_ransy) else pairs_dt$ppm_ransy,
    matched_pairs     = pairs_dt, # SAME FORMAT AS MODE NORMAL
    coverage_ransy    = coverage_ratio,
    RANSY_variabilite = "N/A (ratio mode)",
    avg_variabilite   = NA_real_,
    avg_error_ppm     = NA_real_, # Calculated later in CSV generation code
    groupe_ppms       = matched_ppms_clean,
    annotation_mode   = "intensity_ratio_only"
  )))
}


#' Intensity ratio-based annotation for 2D (HSQC/COSY) - EXACT SAME FORMAT as mode normal
#'
#' @param matched_ppms character vector of matched ppm_matrix values
#' @param molecule_shifts data.table with shift_value_F1 and shift_value_F2 columns
#' @param tol_match_F1 tolerance for F1 dimension
#' @param tol_match_F2 tolerance for F2 dimension
#' @param ratio_min minimum intensity ratio threshold
#' @param spectrum_matrix data.table with ppm_matrix, ppm_F1, ppm_F2 columns
#' @param spectrum_mat matrix of intensities
#' @return list of annotation results (SAME STRUCTURE as mode normal)
intensity_ratio_annotation_2D <- function(matched_ppms, molecule_shifts, tol_match_F1, tol_match_F2,
                                          ratio_min = 0.01, spectrum_matrix, spectrum_mat) {
  if (length(matched_ppms) == 0 || !nrow(molecule_shifts)) {
    return(NULL)
  }

  # Step 1: Calculate average intensities and sort by descending intensity
  idx_match <- match(matched_ppms, spectrum_matrix$ppm_matrix)
  idx_match <- idx_match[!is.na(idx_match)]
  if (!length(idx_match)) {
    return(NULL)
  }

  avg_int <- rowMeans(spectrum_mat[idx_match, , drop = FALSE], na.rm = TRUE)
  matched_ppms_clean <- matched_ppms[!is.na(match(matched_ppms, spectrum_matrix$ppm_matrix))]

  # Sort by descending intensity
  ord_desc <- order(avg_int, decreasing = TRUE)
  ppms_desc <- matched_ppms_clean[ord_desc]
  int_desc <- avg_int[ord_desc]

  # Step 2: Select pivot (most intense signal)
  if (length(ppms_desc) == 0) {
    return(NULL)
  }
  pivot <- ppms_desc[1]
  pivot_int <- int_desc[1]

  if (is.na(pivot_int) || pivot_int == 0) {
    return(NULL)
  }

  # Step 3: Filter by intensity ratio (SAME AS MODE NORMAL)
  ratios <- int_desc / pivot_int
  idx_valide <- which(int_desc <= pivot_int & ratios >= ratio_min & !is.na(int_desc))

  if (!length(idx_valide)) {
    return(NULL)
  }

  sub_chain_vec <- sort(unique(ppms_desc[idx_valide]))

  # Step 4: Build unique_mppm from molecule_shifts (EXACT SAME AS MODE NORMAL)
  unique_mppm <- unique(molecule_shifts[, .(shift_value_F1, shift_value_F2)])

  if (!nrow(unique_mppm)) {
    return(NULL)
  }

  # Merge with spectrum_matrix to get ppm_F1 and ppm_F2 for sub_chain
  avail <- spectrum_matrix[ppm_matrix %in% sub_chain_vec, .(ppm_matrix, ppm_F1, ppm_F2)]
  if (!nrow(avail)) {
    return(NULL)
  }

  # Step 5: RANSY pairings (1-1) - EXACT SAME CODE AS MODE NORMAL (lines 582-594)
  pairs_dt <- data.table(mol_F1 = numeric(0), mol_F2 = numeric(0), ppm_ransy = character(0))

  for (i in seq_len(nrow(unique_mppm))) {
    m1 <- unique_mppm$shift_value_F1[i]
    m2 <- unique_mppm$shift_value_F2[i]
    cand <- avail[abs(ppm_F1 - m1) <= tol_match_F1 & abs(ppm_F2 - m2) <= tol_match_F2]

    if (nrow(cand) > 0) {
      # Distance-based matching only (no CV cost in degraded mode)
      cand[, dist_norm := sqrt(((ppm_F1 - m1) / tol_match_F1)^2 + ((ppm_F2 - m2) / tol_match_F2)^2)]
      best_r <- cand[which.min(dist_norm), ppm_matrix]
      pairs_dt <- rbind(pairs_dt, data.table(mol_F1 = m1, mol_F2 = m2, ppm_ransy = best_r))
      avail <- avail[ppm_matrix != best_r]
    }
  }

  # Step 6: Calculate coverage (SAME AS MODE NORMAL)
  matched_unique <- unique(pairs_dt$ppm_ransy)
  coverage_ratio <- if (nrow(unique_mppm) > 0) length(matched_unique) / nrow(unique_mppm) * 100 else 0
  if (!is.finite(coverage_ratio)) coverage_ratio <- 0

  # Step 7: Calculate intensity ratios string (SAME FORMAT AS MODE NORMAL)
  idx_sub <- match(sub_chain_vec, spectrum_matrix$ppm_matrix)
  idx_sub_valid <- which(!is.na(idx_sub))
  if (!length(idx_sub_valid)) {
    return(NULL)
  }

  sub_chain_clean <- sub_chain_vec[idx_sub_valid]
  avg_mat_all <- rowMeans(spectrum_mat[idx_sub[idx_sub_valid], , drop = FALSE], na.rm = TRUE)
  max_int_all <- max(avg_mat_all, na.rm = TRUE)
  if (!is.finite(max_int_all) || max_int_all == 0) max_int_all <- 1

  intensity_ratios_str <- paste(round(avg_mat_all / max_int_all, 3), collapse = "|")

  # Return result in EXACT SAME FORMAT AS MODE NORMAL
  return(list(list(
    base_ppm          = pivot,
    sub_chain         = sub_chain_clean,
    coverage_sub      = coverage_ratio,
    intensity_ratios  = intensity_ratios_str,
    pivot_ppm         = pivot,
    RANSY_signals_all = sub_chain_clean,
    matched_for_cov   = matched_unique,
    matched_pairs     = pairs_dt, # SAME FORMAT: mol_F1, mol_F2, ppm_ransy
    coverage_ransy    = coverage_ratio,
    RANSY_variabilite = "N/A (ratio mode)",
    avg_variabilite   = NA_real_,
    avg_error_ppm_2D  = NA_real_,
    groupe_ppms       = matched_ppms_clean,
    annotation_mode   = "intensity_ratio_only"
  )))
}


#' NULL coalescing operator (for default parameters)
`%||%` <- function(a, b) if (is.null(a)) b else a
