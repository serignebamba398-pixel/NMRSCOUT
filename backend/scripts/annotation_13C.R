#!/usr/bin/env Rscript

# =====================================================================
# 13C NMR Annotation — Multi-Pivot STOCSY Strategy (V2)
# Pipeline: Match → 3 Anchor Pivots → STOCSY → Ratio → RANSY → Combinatorial Crossing
# Compatible with existing frontend/backend (same parameter names and CSV + new columns)
# =====================================================================

suppressWarnings(suppressMessages({}))

# -------------------------------
# 0) Arguments & basic helpers
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript annotation_13C_stocsy.R <work_dir> <prediction_13C_csv_file> <spectrum_matrix_csv_file> [params_json]\n")
  quit(status = 1)
}
work_dir <- args[1]
prediction_13C_file <- args[2]
spectrum_matrix_file <- args[3]

# Default parameters (backend compatibility + new multi-pivot params)
params <- list(
  tol_correspondance = 2,
  tol_match = 2,
  correlation_threshold = 0.90,
  correlation_method = "pearson",
  coverage_threshold_group = 0,
  coverage_threshold_subchain = 0,
  ratio_min = 0.1,
  ransy_cv_threshold = 0.10,
  coverage_threshold_ransy = 0,
  max_missing_ransy = 50,
  nb_cores = 5,
  progress_style = 3,
  batch_size = 1,
  memory_cleanup_interval = 50,
  enable_checkpoints = TRUE,
  max_memory_mb = 2048,
  max_pivots = 50,
  # === NEW multi-pivot parameters ===
  enable_multi_pivot = TRUE,
  max_pivots_per_anchor = 20,
  max_matches_per_shift = 20,
  intersection_coverage_threshold = 50,
  ransy_score_good_threshold = 60,
  anchor_min_distance_ppm = 0.5,
  triple_weight = 1.3
)

# -------------------------------
# 1) Thread/ENV sanitization
# -------------------------------
sanitize_threads <- function(x) {
  xi <- suppressWarnings(as.integer(x))
  if (is.na(xi) || xi < 1L) NA_integer_ else xi
}
unset_if_invalid_int <- function(var) {
  val <- Sys.getenv(var, unset = "")
  thr <- sanitize_threads(val)
  if (is.na(thr)) Sys.unsetenv(var)
}
invisible(lapply(c("OMP_NUM_THREADS", "R_NTHREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"), unset_if_invalid_int))

# -------------------------------
# 2) SSE logs + memory
# -------------------------------
robust_mem_mb <- function() {
  m <- NA_real_
  if (!is.finite(m) && file.exists("/proc/self/status")) {
    m <- tryCatch(as.numeric(system("grep VmRSS /proc/self/status | awk '{print $2}'", intern = TRUE)[1]) / 1024, error = function(e) NA_real_)
  }
  if (!is.finite(m) && Sys.info()[["sysname"]] == "Darwin") {
    m <- tryCatch(as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024, error = function(e) NA_real_)
  }
  if (!is.finite(m)) m <- NA_real_
  m
}
log_json <- function(obj) {
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n")
  flush.console()
}
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  mem_mb <- robust_mem_mb()
  log_json(list(type = "log", message = sprintf("[%s] %s (Mem: %sMB)", timestamp, message, ifelse(is.na(mem_mb), "NA", round(mem_mb, 1)))))
}
log_progress <- function(current, total, message = "") {
  log_json(list(
    type = "progress", current = current, total = total,
    percent = if (total > 0) round((current / total) * 100, 1) else 0,
    message = message, memory_mb = round(robust_mem_mb(), 1)
  ))
}
log_error <- function(error) log_json(list(type = "error", message = as.character(error)))
log_result <- function(result) log_json(list(type = "result", data = result))
log_total_sources <- function(n) log_json(list(type = "total_sources", count = n))
log_memory <- function(context = "") {
  log_message(paste("Memory used", context, ":", ifelse(is.na(robust_mem_mb()), "NA", round(robust_mem_mb(), 1)), "MB"))
}

# -------------------------------
# 3) Parameter merge
# -------------------------------
if (length(args) >= 4) {
  params_file <- args[4]
  custom_params <- list()
  if (file.exists(file.path(work_dir, params_file))) {
    custom_params <- tryCatch(jsonlite::fromJSON(file.path(work_dir, params_file)), error = function(e) {
      log_error(paste("JSON parsing error:", e$message))
      list()
    })
  } else {
    custom_params <- tryCatch(jsonlite::fromJSON(params_file), error = function(e) {
      log_error(paste("JSON parsing error:", e$message))
      list()
    })
  }
  for (name in names(custom_params)) if (name %in% names(params)) params[[name]] <- custom_params[[name]]
  if (!is.null(custom_params$tol_correspondence)) {
    params$tol_correspondance <- custom_params$tol_correspondence
  }
}

nb_cores <- sanitize_threads(params$nb_cores)
if (is.na(nb_cores)) nb_cores <- 1L
available <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (is.na(available) || available < 1L) available <- 1L
nb_cores <- max(1L, min(as.integer(nb_cores), 3L, as.integer(available)))

# -------------------------------
# 4) Libs + threads data.table
# -------------------------------
suppressWarnings(suppressMessages({
  ok <- TRUE
  tryCatch(
    {
      library(data.table)
      library(doSNOW)
      library(foreach)
      library(parallel)
      library(jsonlite)
    },
    error = function(e) {
      ok <<- FALSE
      log_error(paste("ERROR: Unable to load libraries:", e$message))
    }
  )
  if (!ok) quit(status = 1)
  data.table::setDTthreads(nb_cores)
}))
if (.Platform$OS.type == "windows") {
  lim <- suppressWarnings(as.numeric(params$max_memory_mb))
  if (is.finite(lim) && lim > 256) try(suppressWarnings(memory.limit(size = lim)), silent = TRUE)
}
options(datatable.optimize = 2, datatable.verbose = FALSE, datatable.alloccol = 1024)
gc()

# -------------------------------
# 5) Parameter extraction
# -------------------------------
tol_correspondance <- as.numeric(params$tol_correspondance)
tol_match <- as.numeric(params$tol_match)
correlation_threshold <- as.numeric(params$correlation_threshold)
correlation_method <- if (is.null(params$correlation_method)) "pearson" else tolower(as.character(params$correlation_method))
if (!correlation_method %in% c("pearson", "spearman")) correlation_method <- "pearson"
coverage_threshold_group <- as.numeric(params$coverage_threshold_group)
coverage_threshold_subchain <- as.numeric(params$coverage_threshold_subchain)
ratio_min <- as.numeric(params$ratio_min)
ransy_cv_threshold <- as.numeric(params$ransy_cv_threshold)
coverage_threshold_ransy <- as.numeric(params$coverage_threshold_ransy)
max_missing_ransy <- as.integer(params$max_missing_ransy)
batch_size <- as.integer(params$batch_size)
memory_cleanup_interval <- as.integer(params$memory_cleanup_interval)
enable_checkpoints <- isTRUE(params$enable_checkpoints)
max_pivots <- as.integer(params$max_pivots)
if (is.na(max_pivots) || max_pivots < 1) max_pivots <- 50L

# New multi-pivot parameters
enable_multi_pivot <- isTRUE(params$enable_multi_pivot)
max_pivots_per_anchor <- as.integer(params$max_pivots_per_anchor)
if (is.na(max_pivots_per_anchor) || max_pivots_per_anchor < 1) max_pivots_per_anchor <- 15L
max_matches_per_shift <- as.integer(params$max_matches_per_shift)
if (is.na(max_matches_per_shift) || max_matches_per_shift < 1) max_matches_per_shift <- 5L
intersection_coverage_threshold <- as.numeric(params$intersection_coverage_threshold)
if (is.na(intersection_coverage_threshold)) intersection_coverage_threshold <- 50
ransy_score_good_threshold <- as.numeric(params$ransy_score_good_threshold)
if (is.na(ransy_score_good_threshold)) ransy_score_good_threshold <- 60
anchor_min_distance_ppm <- as.numeric(params$anchor_min_distance_ppm)
if (is.na(anchor_min_distance_ppm)) anchor_min_distance_ppm <- 5
triple_weight <- as.numeric(params$triple_weight)
if (is.na(triple_weight) || triple_weight < 1) triple_weight <- 1.3

# Normalization parameters (default: no normalization for backward compatibility)
use_pqn <- isTRUE(params$use_pqn)
use_csn <- isTRUE(params$use_csn)
normalization_applied <- "none"  # Will be updated after applying normalization

# Global reservoir
intensity_columns <- character(0)
spectrum_mat <- NULL
ppm_idx_lookup <- NULL

# -------------------------------
# 6) Data reading
# -------------------------------
tryCatch(
  {
    log_message("=== 13C ANNOTATION START - MULTI-PIVOT STOCSY V2 ===")
    log_message(paste("Working directory:", work_dir))
    log_message(paste("Threads data.table:", data.table::getDTthreads()))
    log_memory("initial")

    setwd(work_dir)
    log_message("Optimized file reading...")

    prediction_13C <- data.table::fread(
      prediction_13C_file,
      dec = ".", sep = ",", header = TRUE, encoding = "UTF-8",
      showProgress = FALSE, verbose = FALSE
    )
    log_memory("after reading 13C predictions")

    spectrum_matrix <- data.table::fread(
      spectrum_matrix_file,
      dec = ".", sep = ",", header = TRUE, encoding = "UTF-8",
      showProgress = FALSE, verbose = FALSE
    )
    log_memory("after reading matrix")

    log_message(paste("13C Predictions:", nrow(prediction_13C), "rows"))
    log_message(paste("Spectrum matrix:", nrow(spectrum_matrix), "rows"))

    setDT(prediction_13C)
    setDT(spectrum_matrix)

    if (!"Source" %in% names(prediction_13C) && "Filename" %in% names(prediction_13C)) prediction_13C[, Source := Filename]
    if (!"shift_value" %in% names(prediction_13C) && "Delta_ppm" %in% names(prediction_13C)) prediction_13C[, shift_value := Delta_ppm]
    if (!("ECFP4" %in% names(prediction_13C))) prediction_13C[, ECFP4 := ""]

    prediction_13C[, shift_value := as.numeric(shift_value)]
    prediction_13C <- prediction_13C[!is.na(shift_value)]
    if (!"ppm" %in% names(spectrum_matrix)) stop("Missing 'ppm' column in matrix.")
    spectrum_matrix[, ppm := as.numeric(ppm)]
    spectrum_matrix <- spectrum_matrix[!is.na(ppm)]
    setorder(spectrum_matrix, ppm)

    intensity_columns <<- setdiff(names(spectrum_matrix), "ppm")
    if (length(intensity_columns) == 0) stop("No intensity columns (only 'ppm').")
    log_message(paste("Intensity columns:", length(intensity_columns)))

    # 6.5) Sample size detection
    script_dir <- tryCatch(
      {
        dirname(sys.frame(1)$ofile)
      },
      error = function(e) {
        args_full <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args_full, value = TRUE)
        if (length(file_arg) > 0) {
          dirname(sub("^--file=", "", file_arg[1]))
        } else {
          getwd()
        }
      }
    )

    detector_path <- file.path(script_dir, "sample_size_detector.R")
    if (!file.exists(detector_path)) {
      log_error(paste("Cannot find sample_size_detector.R at:", detector_path))
      quit(status = 1)
    }
    source(detector_path)

    mode_info <- detect_annotation_mode(spectrum_matrix)
    log_message(mode_info$warning)

    USE_STOCSY <<- mode_info$use_stocsy
    USE_RANSY <<- mode_info$use_ransy
    N_SAMPLES <<- mode_info$n_samples
    ANNOTATION_MODE <<- mode_info$mode

    if (!USE_STOCSY) {
      log_message(sprintf(
        "DEGRADED MODE: Only %d sample(s) - STOCSY/RANSY disabled, using intensity ratios only",
        N_SAMPLES
      ))
    }

    spectrum_matrix <<- spectrum_matrix
    spectrum_mat <<- as.matrix(spectrum_matrix[, ..intensity_columns])

    # Apply spectral normalization if requested
    if (use_pqn || use_csn) {
      norm_script <- file.path(dirname(script_path), "spectral_normalization.R")
      if (file.exists(norm_script)) {
        source(norm_script)
        norm_result <- apply_normalization(
          mat = spectrum_mat,
          use_pqn = use_pqn,
          use_csn = use_csn,
          log_func = log_message
        )
        spectrum_mat <<- norm_result$matrix
        normalization_applied <<- norm_result$methods_applied
      } else {
        log_message(paste("[NORMALIZATION] Warning: spectral_normalization.R not found at", norm_script))
        normalization_applied <<- "none (file not found)"
      }
    } else {
      log_message("[NORMALIZATION] No normalization applied (raw intensities)")
    }

    ppm_idx_lookup <<- setNames(seq_len(nrow(spectrum_matrix)), as.character(spectrum_matrix$ppm))
    log_message(paste("PPM index lookup table created:", length(ppm_idx_lookup), "entries"))
    log_message(paste("Multi-pivot mode:", ifelse(enable_multi_pivot, "ENABLED", "DISABLED (legacy)")))

    tableau_final <- prediction_13C
    gc()
    log_memory("after initial cleanup")
  },
  error = function(e) {
    log_error(paste("Data reading error:", e$message))
    quit(status = 1)
  }
)

# -------------------------------
# 7) CSV sanitization helper
# -------------------------------
sanitize_for_csv <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\r\n\t]", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

# -------------------------------
# 7.5) Fast ppm index lookup
# -------------------------------
fast_ppm_idx <- function(ppms) {
  if (is.null(ppm_idx_lookup) || length(ppms) == 0L) return(integer(0))
  keys <- as.character(ppms)
  idx <- ppm_idx_lookup[keys]
  idx[!is.na(idx)]
}

# -------------------------------
# 8) Matching with max_matches_per_shift control
# -------------------------------
trouver_correspondances <- function(shifts, tol) {
  if (length(shifts) == 0L || nrow(spectrum_matrix) == 0L) {
    return(data.table(ppm_matrix = numeric(), ppm_molecule = numeric(), distance = numeric()))
  }
  temp <- data.table(shift = as.numeric(shifts))
  m <- spectrum_matrix[, .(ppm, ppm_min = ppm - tol, ppm_max = ppm + tol)]
  res <- m[temp,
    on = .(ppm_min <= shift, ppm_max >= shift),
    allow.cartesian = TRUE, nomatch = 0,
    .(ppm_matrix = ppm, ppm_molecule = shift, distance = abs(shift - ppm))
  ]
  if (nrow(res)) setorder(res, distance)
  res
}

# Build filtered matched_ppms: max N most INTENSE matches per predicted shift, deduplicated
# With minimum spacing to avoid picking multiple points from the same peak
build_matched_ppms <- function(shifts, tol, max_per_shift, min_spacing = 0.01) {
  res <- trouver_correspondances(shifts, tol)
  if (nrow(res) == 0) return(list(matched_ppms = numeric(0), matches_dt = res))

  # Lookup mean intensity for each matched experimental signal
  unique_matrix_ppms <- unique(res$ppm_matrix)
  int_lookup <- data.table(
    ppm_matrix = unique_matrix_ppms,
    mean_intensity = {
      idx <- fast_ppm_idx(unique_matrix_ppms)
      if (length(idx) > 0) rowMeans(spectrum_mat[idx, , drop = FALSE], na.rm = TRUE)
      else rep(NA_real_, length(unique_matrix_ppms))
    }
  )
  res <- merge(res, int_lookup, by = "ppm_matrix", all.x = TRUE)

  # For each predicted shift: greedy selection by intensity with min spacing
  res_filtered <- res[, {
    dt_grp <- .SD[order(-mean_intensity)]
    if (nrow(dt_grp) == 0) {
      dt_grp[0]
    } else {
      selected_idx <- 1L
      for (ri in seq_len(nrow(dt_grp))[-1]) {
        if (all(abs(dt_grp$ppm_matrix[ri] - dt_grp$ppm_matrix[selected_idx]) >= min_spacing)) {
          selected_idx <- c(selected_idx, ri)
        }
        if (length(selected_idx) >= max_per_shift) break
      }
      dt_grp[selected_idx]
    }
  }, by = ppm_molecule]

  res_filtered[, mean_intensity := NULL]

  matched_ppms <- sort(unique(res_filtered$ppm_matrix))
  list(matched_ppms = matched_ppms, matches_dt = res_filtered)
}

# -------------------------------
# 9) STOCSY — VECTORIZED version
# -------------------------------
calculer_STOCSY <- function(ppms, intensity_cols, driver_ppm, r_min, method = "pearson") {
  if (length(ppms) == 0L) return(data.table(ppm = numeric(), r = numeric()))

  idx <- fast_ppm_idx(ppms)
  if (!length(idx)) return(data.table(ppm = numeric(), r = numeric()))

  mat <- spectrum_mat[idx, , drop = FALSE]
  row_ppm <- spectrum_matrix$ppm[idx]

  driver_idx <- which.min(abs(row_ppm - driver_ppm))
  if (!length(driver_idx) || is.na(driver_idx) || driver_idx < 1L) {
    return(data.table(ppm = numeric(), r = numeric()))
  }

  driver_v <- mat[driver_idx, ]
  valid_cols <- which(is.finite(driver_v) & !is.na(driver_v))
  if (length(valid_cols) < 3L) return(data.table(ppm = numeric(), r = numeric()))

  mat2 <- mat[, valid_cols, drop = FALSE]
  drv <- as.numeric(driver_v[valid_cols])

  if (tolower(method) == "spearman") {
    mat2 <- t(apply(mat2, 1, rank, na.last = "keep"))
    drv <- rank(drv, na.last = "keep")
  }

  mat_centered <- mat2 - rowMeans(mat2, na.rm = TRUE)
  drv_centered <- drv - mean(drv, na.rm = TRUE)
  numerator <- as.numeric(mat_centered %*% drv_centered)
  mat_sd <- sqrt(rowSums(mat_centered^2, na.rm = TRUE))
  drv_sd <- sqrt(sum(drv_centered^2, na.rm = TRUE))
  denominator <- mat_sd * drv_sd
  r_vec <- numerator / denominator
  r_vec[!is.finite(r_vec)] <- NA_real_

  dt <- data.table(ppm = row_ppm, r = r_vec)
  dt <- dt[is.finite(r) & abs(r) >= r_min]
  if (!nrow(dt)) return(data.table(ppm = numeric(), r = numeric()))
  dt[order(ppm)]
}

# -------------------------------
# 10) Group filtering by coverage
# -------------------------------
gfilter_groups_par_coverage <- function(groups, molecule_ppm, tol, cov_thresh) {
  out <- list()
  up <- sort(unique(molecule_ppm))
  if (length(up) == 0) return(out)
  for (grp in groups) {
    cand <- sort(as.numeric(grp))
    matched <- integer()
    avail <- cand
    for (x in up) {
      in_r <- avail[abs(avail - x) <= tol]
      if (length(in_r) > 0) {
        best <- in_r[which.min(abs(in_r - x))]
        matched <- c(matched, best)
        avail <- setdiff(avail, best)
      }
    }
    cov <- length(matched) / length(up) * 100
    if (cov >= cov_thresh) out[[length(out) + 1]] <- grp
  }
  out
}

# -------------------------------
# 11) RANSY — VECTORIZED version
# -------------------------------
calculer_RANSY <- function(ppms, intensity_cols, cv_th, driver_ppm) {
  if (length(ppms) == 0L) return(NULL)

  idx <- fast_ppm_idx(ppms)
  if (!length(idx)) return(NULL)

  mat <- spectrum_mat[idx, , drop = FALSE]
  row_ppm <- spectrum_matrix$ppm[idx]

  driver_idx <- which.min(abs(row_ppm - driver_ppm))
  if (!length(driver_idx) || is.na(driver_idx) || driver_idx < 1L) return(NULL)

  driver_v <- mat[driver_idx, ]
  valid_cols <- which(driver_v > 0 & !is.na(driver_v))
  if (length(valid_cols) < 2L) return(NULL)

  mat2 <- mat[, valid_cols, drop = FALSE]
  ratio_mat <- sweep(mat2, 2, driver_v[valid_cols], "/")

  n_cols <- ncol(ratio_mat)
  row_means <- rowMeans(ratio_mat, na.rm = TRUE)
  row_sq_means <- rowMeans(ratio_mat^2, na.rm = TRUE)
  row_var <- (row_sq_means - row_means^2) * n_cols / (n_cols - 1)
  row_sd <- sqrt(pmax(row_var, 0))
  cvs <- ifelse(row_means > 0 & is.finite(row_means), row_sd / row_means, Inf)
  cvs[!is.finite(cvs)] <- Inf

  data.table(ppm = row_ppm, cv = cvs)[cv <= cv_th]
}

# -------------------------------
# 12) Substring coverage helpers
# -------------------------------
get_overlap_groups_1D <- function(ppms, overlap_tol) {
  n <- length(ppms)
  groups <- seq_len(n)
  if (n <= 1) return(groups)
  ppms_s <- sort(ppms)
  for (i in seq_len(n - 1)) {
    if (abs(ppms_s[i + 1] - ppms_s[i]) <= overlap_tol) {
      groups[groups == groups[i + 1]] <- groups[i]
    }
  }
  groups
}

calculer_coverage_sub_ratio <- function(sub_chain, molecule_ppm, tol_match, overlap_tol = 0.5) {
  # NOTE: overlap_tol parameter kept for API compatibility but NO LONGER USED
  # Score is now calculated per individual signal, not per group
  if (length(sub_chain) == 0 || length(molecule_ppm) == 0) return(0)
  unique_mppm <- sort(unique(molecule_ppm))
  avail_sub <- sort(sub_chain)
  nb_signals <- length(unique_mppm)
  signals_matches <- 0L

  for (u in unique_mppm) {
    candidats <- avail_sub[abs(avail_sub - u) <= tol_match]
    if (length(candidats) > 0) {
      meilleur <- candidats[which.min(abs(candidats - u))]
      avail_sub <- setdiff(avail_sub, meilleur)
      signals_matches <- signals_matches + 1L
    }
  }
  signals_matches / nb_signals * 100
}

# -------------------------------
# 12.5) EXACT intersection (for crossing sub-groups)
# All ransy_ppms come from the SAME spectrum_matrix → exact match
# -------------------------------
exact_intersect <- function(ppms_a, ppms_b) {
  if (length(ppms_a) == 0 || length(ppms_b) == 0) return(numeric(0))
  sort(intersect(ppms_a, ppms_b))
}

# -------------------------------
# 12.6) Anchor selection
# -------------------------------
select_anchors <- function(molecule_ppm, min_distance = 5) {
  sorted_ppm <- sort(unique(molecule_ppm))
  n <- length(sorted_ppm)
  if (n == 0) return(list(anchors = NULL, anchors_close = FALSE))
  if (n == 1) return(list(anchors = sorted_ppm[1], anchors_close = FALSE))
  if (n == 2) return(list(anchors = sorted_ppm, anchors_close = abs(diff(sorted_ppm)) < min_distance))

  anchor_debut <- sorted_ppm[1]
  anchor_milieu <- sorted_ppm[ceiling(n / 2)]
  anchor_fin <- sorted_ppm[n]

  anchors_close <- (abs(anchor_debut - anchor_milieu) < min_distance) ||
    (abs(anchor_milieu - anchor_fin) < min_distance)

  if (anchors_close) {
    anchor_milieu <- sorted_ppm[round(n / 2)]
    anchors_close <- (abs(anchor_debut - anchor_milieu) < min_distance) ||
      (abs(anchor_milieu - anchor_fin) < min_distance)
  }

  list(anchors = c(anchor_debut, anchor_milieu, anchor_fin), anchors_close = anchors_close)
}

# -------------------------------
# 13) Single-pivot STOCSY+Ratio+RANSY pipeline (extracted as reusable function)
# Returns a sub-group result or NULL if pivot fails any filter
# -------------------------------
run_single_pivot_pipeline <- function(pivot_ppm, anchor_ppm, group_index,
                                      matched_ppms, molecule_ppm, tol_match,
                                      correlation_threshold, correlation_method,
                                      coverage_threshold_group, coverage_threshold_subchain,
                                      ratio_min, ransy_cv_threshold, coverage_threshold_ransy,
                                      max_missing, overlap_tol = 0.5) {
  unique_mppm <- sort(unique(molecule_ppm))

  # STOCSY
  dt_stocsy <- calculer_STOCSY(
    ppms = matched_ppms,
    intensity_cols = intensity_columns,
    driver_ppm = pivot_ppm,
    r_min = correlation_threshold,
    method = correlation_method
  )
  if (is.null(dt_stocsy) || !nrow(dt_stocsy)) return(NULL)
  stocsy_ppms <- sort(unique(dt_stocsy$ppm))

  # Coverage group filter
  if (coverage_threshold_group > 0) {
    grps_filtered <- gfilter_groups_par_coverage(list(stocsy_ppms), unique_mppm, tol_match, coverage_threshold_group)
    if (length(grps_filtered) == 0) return(NULL)
  }

  # Intensity ratio filter
  # In multi-pivot mode: symmetric filter — keep signals with reasonable ratio
  # to the pivot in EITHER direction (signal/pivot OR pivot/signal >= ratio_min).
  # This avoids excluding signals more intense than the pivot, which would
  # create incompatible sub-groups across different pivot intensities.
  # In legacy mode: keep original behavior (signal <= pivot intensity).
  pivot_idx <- fast_ppm_idx(pivot_ppm)
  if (!length(pivot_idx)) return(NULL)
  pivot_int <- mean(spectrum_mat[pivot_idx, , drop = FALSE], na.rm = TRUE)
  if (!is.finite(pivot_int) || pivot_int == 0) return(NULL)

  idx_grp <- fast_ppm_idx(stocsy_ppms)
  if (!length(idx_grp)) return(NULL)
  grp_clean <- stocsy_ppms[stocsy_ppms %in% spectrum_matrix$ppm[idx_grp]]
  grp_int <- rowMeans(spectrum_mat[idx_grp, , drop = FALSE], na.rm = TRUE)
  ratios <- grp_int / pivot_int

  if (enable_multi_pivot) {
    # Symmetric: keep if ratio OR inverse ratio >= ratio_min
    inverse_ratios <- pivot_int / grp_int
    idx_valide <- which(ratios >= ratio_min | inverse_ratios >= ratio_min)
  } else {
    # Legacy: only keep signals weaker than pivot
    idx_valide <- which(grp_int <= pivot_int & ratios >= ratio_min)
  }
  if (!length(idx_valide)) return(NULL)
  sub_chain_vec <- sort(unique(grp_clean[idx_valide]))

  # Coverage sub (per signal, no grouping)
  nb_signals <- length(unique_mppm)
  avail_final <- sort(sub_chain_vec)
  nb_signals_sub_matches <- 0L
  for (u in unique_mppm) {
    candidats <- avail_final[abs(avail_final - u) <= tol_match]
    if (length(candidats) > 0) {
      meilleur <- candidats[which.min(abs(candidats - u))]
      avail_final <- setdiff(avail_final, meilleur)
      nb_signals_sub_matches <- nb_signals_sub_matches + 1L
    }
  }
  coverage_sub_final <- nb_signals_sub_matches / nb_signals * 100
  if (!is.finite(coverage_sub_final)) coverage_sub_final <- 0
  if (coverage_sub_final < coverage_threshold_subchain) return(NULL)

  # RANSY
  dt_ransy <- calculer_RANSY(sub_chain_vec, intensity_columns, ransy_cv_threshold, pivot_ppm)
  if (is.null(dt_ransy) || !nrow(dt_ransy)) return(NULL)
  ransy_pass <- sort(unique(dt_ransy$ppm))

  # RANSY pairings (per signal, no grouping)
  pairs_dt <- data.table(ppm_molecule = numeric(0), ppm_ransy = numeric(0))
  avail_ransy <- sort(ransy_pass)
  nb_signals_ransy_matches <- 0L
  for (u in unique_mppm) {
    cand_r <- avail_ransy[abs(avail_ransy - u) <= tol_match]
    if (length(cand_r) > 0) {
      best_r <- cand_r[which.min(abs(cand_r - u))]
      avail_ransy <- setdiff(avail_ransy, best_r)
      nb_signals_ransy_matches <- nb_signals_ransy_matches + 1L
      pairs_dt <- rbind(pairs_dt, data.table(ppm_molecule = u, ppm_ransy = best_r))
    }
  }
  matched_unique <- unique(pairs_dt$ppm_ransy)
  coverage_ransy <- nb_signals_ransy_matches / nb_signals * 100
  if (coverage_ransy > coverage_sub_final) coverage_ransy <- coverage_sub_final

  missing_cnt <- nb_signals - nb_signals_ransy_matches
  if (missing_cnt > max_missing || coverage_ransy < coverage_threshold_ransy) return(NULL)

  # Score
  avg_cv <- round(mean(dt_ransy[ppm %in% matched_unique, cv], na.rm = TRUE), 4)
  ransy_score <- coverage_ransy

  # Intensity ratios
  idx_sub <- fast_ppm_idx(sub_chain_vec)
  if (!length(idx_sub)) return(NULL)
  sub_chain_clean <- sub_chain_vec[sub_chain_vec %in% spectrum_matrix$ppm[idx_sub]]
  avg_mat_all <- rowMeans(spectrum_mat[idx_sub, , drop = FALSE], na.rm = TRUE)
  max_int_all <- max(avg_mat_all, na.rm = TRUE)
  if (!is.finite(max_int_all) || max_int_all == 0) max_int_all <- 1
  intensity_ratios_str <- paste(round(avg_mat_all / max_int_all, 3), collapse = ";")
  ransy_var_str <- paste(paste(dt_ransy$ppm, round(dt_ransy$cv, 4), sep = ":"), collapse = "|")

  avg_error_ppm <- if (nrow(pairs_dt) > 0) {
    round(mean(abs(pairs_dt$ppm_molecule - pairs_dt$ppm_ransy), na.rm = TRUE), 4)
  } else {
    NA_real_
  }

  list(
    pivot_ppm = pivot_ppm,
    anchor_ppm = anchor_ppm,
    group_index = group_index,
    stocsy_ppms = stocsy_ppms,
    sub_chain = sub_chain_clean,
    ransy_ppms = ransy_pass,
    matched_pairs = pairs_dt,
    matched_for_cov = matched_unique,
    coverage_sub = coverage_sub_final,
    coverage_ransy = coverage_ransy,
    ransy_score = ransy_score,
    avg_cv = avg_cv,
    avg_error_ppm = avg_error_ppm,
    intensity_ratios = intensity_ratios_str,
    RANSY_variabilite = ransy_var_str,
    groupe_ppms = stocsy_ppms
  )
}

# -------------------------------
# 13.5) Multi-pivot crossing logic — Score-based (no cascade)
# Strategy: compute ALL intersections (triple + pair), score each on filtered
# signals, weight triples higher, best weighted score wins.
# -------------------------------
cross_sub_groups <- function(sub_groups, molecule_ppm, tol_match,
                             inter_cov_threshold, ransy_good_threshold,
                             triple_weight = 1.3, overlap_tol = 0.5) {
  valid_group_indices <- which(sapply(sub_groups, function(g) length(g) > 0))
  n_valid <- length(valid_group_indices)

  if (n_valid == 0) return(NULL)

  # Helper: choose sub-group with best ransy_score
  pick_best_sg <- function(sg_list) {
    scores <- sapply(sg_list, function(x) x$ransy_score)
    sg_list[[which.max(scores)]]
  }

  # =====================================================================
  # score_intersection: Given intersection signals + representative pivot,
  # recalculate score on filtered signals only.
  # NOTE: No grouping - score calculated per individual signal
  # =====================================================================
  score_intersection <- function(inter_ppms, representative_sg, molecule_ppm,
                                 tol_match, overlap_tol) {
    # NOTE: overlap_tol parameter kept for API compatibility but NO LONGER USED
    if (length(inter_ppms) == 0) {
      return(list(score = 0, coverage = 0, inter_ppms = numeric(0),
                  pairs = data.table(ppm_molecule = numeric(0), ppm_ransy = numeric(0)),
                  avg_error = NA_real_))
    }

    filtered <- sort(inter_ppms)
    unique_mppm <- sort(unique(molecule_ppm))
    nb_signals <- length(unique_mppm)

    # Recalculate coverage on intersection signals (per signal, no grouping)
    avail_r <- sort(filtered)
    nb_matches <- 0L
    recalc_pairs <- data.table(ppm_molecule = numeric(0), ppm_ransy = numeric(0))

    for (u in unique_mppm) {
      cand_r <- avail_r[abs(avail_r - u) <= tol_match]
      if (length(cand_r) > 0) {
        best_r <- cand_r[which.min(abs(cand_r - u))]
        avail_r <- setdiff(avail_r, best_r)
        nb_matches <- nb_matches + 1L
        recalc_pairs <- rbind(recalc_pairs, data.table(ppm_molecule = u, ppm_ransy = best_r))
      }
    }

    coverage <- nb_matches / nb_signals * 100
    avg_error <- if (nrow(recalc_pairs) > 0) {
      round(mean(abs(recalc_pairs$ppm_molecule - recalc_pairs$ppm_ransy), na.rm = TRUE), 4)
    } else NA_real_

    # Score = coverage of predicted shifts by intersection signals
    inter_score <- coverage

    list(score = inter_score, coverage = coverage, inter_ppms = filtered,
         pairs = recalc_pairs, avg_error = avg_error)
  }

  # =====================================================================
  # build_final: Apply intersection filtering to representative sub-group
  # CRITICAL: Only keep signals that match a PREDICTED shift of the molecule.
  # Extra signals that correlate but don't match any prediction are removed.
  # =====================================================================
  build_final <- function(representative_sg, scored, conf_level, n_conv,
                          w_pivots, excl_grp = NA, excl_reason = NA) {
    sg <- representative_sg

    if (length(scored$inter_ppms) > 0 && n_conv >= 2) {
      filtered_ransy <- scored$inter_ppms

      # FINAL FILTER: Only keep intersection signals that match a predicted shift
      # This removes fortuitous signals from other molecules that passed STOCSY/RANSY
      # because they belong to a dominant compound whose signals correlate well.
      unique_mppm <- sort(unique(molecule_ppm))
      kept_signals <- numeric(0)
      avail_inter <- sort(filtered_ransy)
      for (mp in unique_mppm) {
        cands <- avail_inter[abs(avail_inter - mp) <= tol_match]
        if (length(cands) > 0) {
          best <- cands[which.min(abs(cands - mp))]
          kept_signals <- c(kept_signals, best)
          avail_inter <- setdiff(avail_inter, best)
        }
      }
      filtered_ransy <- sort(unique(kept_signals))

      # Update signals with prediction-filtered data
      sg$ransy_ppms_unfiltered <- sg$ransy_ppms
      sg$ransy_ppms <- filtered_ransy
      sg$matched_pairs <- scored$pairs[ppm_ransy %in% filtered_ransy]
      sg$matched_for_cov <- unique(sg$matched_pairs$ppm_ransy)
      sg$avg_error_ppm <- if (nrow(sg$matched_pairs) > 0) {
        round(mean(abs(sg$matched_pairs$ppm_molecule - sg$matched_pairs$ppm_ransy), na.rm = TRUE), 4)
      } else NA_real_

      # Recalculate coverage on filtered signals
      sg$coverage_ransy <- scored$coverage

      # Filter sub_chain
      sg$sub_chain <- sort(intersect(sg$sub_chain, filtered_ransy))
      sg$coverage_sub <- calculer_coverage_sub_ratio(sg$sub_chain, molecule_ppm, tol_match)

      # Recalculate intensity ratios on filtered signals only
      idx_filt <- fast_ppm_idx(filtered_ransy)
      if (length(idx_filt) > 0) {
        avg_filt <- rowMeans(spectrum_mat[idx_filt, , drop = FALSE], na.rm = TRUE)
        max_filt <- max(avg_filt, na.rm = TRUE)
        if (!is.finite(max_filt) || max_filt == 0) max_filt <- 1
        sg$intensity_ratios <- paste(round(avg_filt / max_filt, 3), collapse = ";")
      }

      sg$ransy_score <- scored$score
    }

    sg$confidence_level <- conf_level
    sg$n_groups_converged <- n_conv
    sg$intersection_coverage <- scored$coverage
    sg$intersection_ppms <- scored$inter_ppms
    sg$winning_pivots <- w_pivots
    sg$excluded_group <- excl_grp
    sg$exclusion_reason <- excl_reason
    sg$note <- NA_character_
    sg
  }

  # =====================================================================
  # With one group, use the best solo result without cross-group comparison.
  # =====================================================================
  if (n_valid == 1) {
    g <- valid_group_indices[1]
    best <- pick_best_sg(sub_groups[[g]])
    best$confidence_level <- "LOW"
    best$n_groups_converged <- 1L
    best$intersection_coverage <- 0
    best$intersection_ppms <- numeric(0)
    best$winning_pivots <- best$pivot_ppm
    best$excluded_group <- NA_integer_
    best$exclusion_reason <- NA_character_
    best$note <- "single_group_only"
    return(best)
  }

  # =====================================================================
  # Collect ALL candidate intersections: TRIPLE + PAIR
  # Score each, weight triples higher. Best weighted score wins.
  # =====================================================================
  all_candidates <- list()
  cand_count <- 0

  # --- TRIPLE intersections (only if 3+ groups) ---
  if (n_valid >= 3) {
    g1 <- valid_group_indices[1]
    g2 <- valid_group_indices[2]
    g3 <- valid_group_indices[3]
    n1 <- length(sub_groups[[g1]])
    n2 <- length(sub_groups[[g2]])
    n3 <- length(sub_groups[[g3]])

    for (i in seq_len(n1)) {
      for (j in seq_len(n2)) {
        for (k in seq_len(n3)) {
          sg1 <- sub_groups[[g1]][[i]]
          sg2 <- sub_groups[[g2]][[j]]
          sg3 <- sub_groups[[g3]][[k]]

          inter_12 <- exact_intersect(sg1$ransy_ppms, sg2$ransy_ppms)
          inter_123 <- exact_intersect(inter_12, sg3$ransy_ppms)
          if (length(inter_123) == 0) next

          three_sg <- list(sg1, sg2, sg3)
          rep_sg <- pick_best_sg(three_sg)

          scored <- score_intersection(inter_123, rep_sg, molecule_ppm,
                                       tol_match, overlap_tol)
          # Triple weight bonus
          weighted_score <- scored$score * triple_weight

          cand_count <- cand_count + 1
          all_candidates[[cand_count]] <- list(
            type = "triple", n_conv = 3L,
            weighted_score = weighted_score, raw_score = scored$score,
            scored = scored, representative = rep_sg,
            w_pivots = c(sg1$pivot_ppm, sg2$pivot_ppm, sg3$pivot_ppm),
            excl_grp = NA_integer_, excl_reason = NA_character_
          )
        }
      }
    }
  }

  # --- PAIR intersections (all pairs of valid groups) ---
  if (n_valid >= 2) {
    group_pair_combos <- combn(valid_group_indices, 2, simplify = FALSE)

    for (gp in group_pair_combos) {
      ga <- gp[1]
      gb <- gp[2]
      excluded_groups <- setdiff(valid_group_indices, gp)

      for (i in seq_len(length(sub_groups[[ga]]))) {
        for (j in seq_len(length(sub_groups[[gb]]))) {
          sg_a <- sub_groups[[ga]][[i]]
          sg_b <- sub_groups[[gb]][[j]]

          inter_ab <- exact_intersect(sg_a$ransy_ppms, sg_b$ransy_ppms)
          if (length(inter_ab) == 0) next

          rep_sg <- if (sg_a$ransy_score >= sg_b$ransy_score) sg_a else sg_b

          scored <- score_intersection(inter_ab, rep_sg, molecule_ppm,
                                       tol_match, overlap_tol)
          # Pair: no weight bonus (× 1.0)
          weighted_score <- scored$score

          cand_count <- cand_count + 1
          all_candidates[[cand_count]] <- list(
            type = "pair", n_conv = 2L,
            weighted_score = weighted_score, raw_score = scored$score,
            scored = scored, representative = rep_sg,
            w_pivots = c(sg_a$pivot_ppm, sg_b$pivot_ppm),
            excl_grp = if (length(excluded_groups) > 0) excluded_groups[1] else NA_integer_,
            excl_reason = NA_character_
          )
        }
      }
    }
  }

  # =====================================================================
  # No valid intersection found → fallback to best solo sub-group
  # =====================================================================
  if (length(all_candidates) == 0) {
    all_sg <- unlist(sub_groups[valid_group_indices], recursive = FALSE)
    best <- pick_best_sg(all_sg)
    best$confidence_level <- "LOW"
    best$n_groups_converged <- 1L
    best$intersection_coverage <- 0
    best$intersection_ppms <- numeric(0)
    best$winning_pivots <- best$pivot_ppm
    best$excluded_group <- NA_integer_
    best$exclusion_reason <- NA_character_
    best$note <- "no_intersection_found_fallback_solo"
    return(best)
  }

  # =====================================================================
  # Pick best weighted_score → this is the winner
  # =====================================================================
  w_scores <- sapply(all_candidates, function(x) x$weighted_score)
  best_cand <- all_candidates[[which.max(w_scores)]]

  # Confidence level
  if (best_cand$type == "triple") {
    conf <- "HIGH"
  } else {
    conf <- if (best_cand$raw_score >= ransy_good_threshold) "MEDIUM_HIGH" else "MEDIUM"
  }

  final <- build_final(
    representative_sg = best_cand$representative,
    scored = best_cand$scored,
    conf_level = conf,
    n_conv = best_cand$n_conv,
    w_pivots = best_cand$w_pivots,
    excl_grp = best_cand$excl_grp,
    excl_reason = best_cand$excl_reason
  )

  return(final)
}


# -------------------------------
# 14) Main pipeline function (multi-pivot or legacy)
# -------------------------------
former_groupes_stocsy <- function(matched_ppms, molecule_ppm, tol_match,
                                  correlation_threshold, correlation_method = "pearson",
                                  coverage_threshold_group, coverage_threshold_subchain,
                                  ratio_min, ransy_cv_threshold, coverage_threshold_ransy,
                                  max_missing, max_pivots = 200) {
  # Adaptive mode: n < 5 → intensity ratios only
  if (exists("USE_STOCSY") && !USE_STOCSY) {
    return(intensity_ratio_annotation(
      matched_ppms = matched_ppms, molecule_ppm = molecule_ppm,
      tol_match = tol_match, ratio_min = ratio_min,
      spectrum_matrix = spectrum_matrix, spectrum_mat = spectrum_mat
    ))
  }

  if (length(matched_ppms) == 0 || length(molecule_ppm) == 0) return(NULL)

  # ================================================================
  # MULTI-PIVOT MODE (new)
  # ================================================================
  if (enable_multi_pivot) {
    # Step 1: Select anchors
    anchor_info <- select_anchors(molecule_ppm, anchor_min_distance_ppm)
    if (is.null(anchor_info$anchors)) return(NULL)
    anchors <- anchor_info$anchors

    # Step 2: Build pivot groups — with DEDUPLICATION
    # Each experimental signal can only be a pivot in ONE group (the one
    # whose anchor is closest). This prevents overlapping groups from
    # producing identical sub-groups when anchors are close.
    sub_groups <- list()
    used_pivots <- numeric(0)  # track pivots already assigned to a group

    for (a_idx in seq_along(anchors)) {
      anchor <- anchors[a_idx]
      candidates <- matched_ppms[abs(matched_ppms - anchor) <= tol_correspondance]
      # Exclude pivots already assigned to a previous group
      candidates <- setdiff(candidates, used_pivots)

      # Sort by INTENSITY (descending) — intense pivots give better STOCSY and RANSY
      # Then GREEDY selection with minimum spacing: each new pivot must be at least
      # min_pivot_spacing ppm away from ALL already selected pivots in this group.
      # This avoids selecting multiple points from the same peak.
      min_pivot_spacing <- 0.01  # ppm, avoids redundant points on same peak (13C)
      if (length(candidates) > 0) {
        cand_idx <- fast_ppm_idx(candidates)
        if (length(cand_idx) > 0) {
          cand_int <- rowMeans(spectrum_mat[cand_idx, , drop = FALSE], na.rm = TRUE)
          candidates_sorted <- candidates[order(cand_int, decreasing = TRUE)]

          # Greedy selection: take most intense, then next most intense that is far enough
          selected <- candidates_sorted[1]
          for (ci in seq_along(candidates_sorted)[-1]) {
            cp <- candidates_sorted[ci]
            if (all(abs(cp - selected) >= min_pivot_spacing)) {
              selected <- c(selected, cp)
            }
            if (length(selected) >= max_pivots_per_anchor) break
          }
          candidates <- selected
        } else {
          candidates <- head(candidates, max_pivots_per_anchor)
        }
      }

      if (length(candidates) == 0) {
        log_message(sprintf("  Anchor %d (%.3f): 0 candidates (after dedup)", a_idx, anchor))
        sub_groups[[a_idx]] <- list()
        next
      }

      # Mark these pivots as used
      used_pivots <- c(used_pivots, candidates)

      log_message(sprintf("  Anchor %d (%.3f): %d candidates, testing pivots...",
                          a_idx, anchor, length(candidates)))

      # Step 3: Run pipeline for each pivot in this group
      sg_results <- list()
      for (pivot in candidates) {
        result <- run_single_pivot_pipeline(
          pivot_ppm = pivot, anchor_ppm = anchor, group_index = a_idx,
          matched_ppms = matched_ppms, molecule_ppm = molecule_ppm,
          tol_match = tol_match,
          correlation_threshold = correlation_threshold,
          correlation_method = correlation_method,
          coverage_threshold_group = coverage_threshold_group,
          coverage_threshold_subchain = coverage_threshold_subchain,
          ratio_min = ratio_min, ransy_cv_threshold = ransy_cv_threshold,
          coverage_threshold_ransy = coverage_threshold_ransy,
          max_missing = max_missing, overlap_tol = 0.5
        )
        if (!is.null(result)) sg_results[[length(sg_results) + 1]] <- result
      }
      sub_groups[[a_idx]] <- sg_results
      log_message(sprintf("  Anchor %d (%.3f): %d/%d pivots survived all filters",
                          a_idx, anchor, length(sg_results), length(candidates)))
      if (length(sg_results) > 0) {
        best_sg <- sg_results[[which.max(sapply(sg_results, function(x) x$ransy_score))]]
        log_message(sprintf("    Best sub-group: %d ransy_ppms, coverage=%.1f%%",
                            length(best_sg$ransy_ppms), best_sg$coverage_ransy))
      }
    }

    # Log summary before crossing
    n_valid_grps <- sum(sapply(sub_groups, function(g) length(g) > 0))
    log_message(sprintf("  Multi-pivot: %d/%d groups have surviving pivots", n_valid_grps, length(anchors)))

    # Step 4: Combinatorial crossing
    final <- cross_sub_groups(sub_groups, molecule_ppm, tol_match,
                              intersection_coverage_threshold, ransy_score_good_threshold,
                              triple_weight = triple_weight, overlap_tol = 0.5)
    if (is.null(final)) return(NULL)

    # Wrap in list for compatibility (pipeline expects list of results)
    return(list(final))
  }

  # ================================================================
  # LEGACY MODE (original algorithm, pivots by descending intensity)
  # ================================================================
  idx_match <- match(matched_ppms, spectrum_matrix$ppm)
  idx_match <- idx_match[!is.na(idx_match)]
  if (!length(idx_match)) return(NULL)

  matched_ppms_clean <- matched_ppms[!is.na(match(matched_ppms, spectrum_matrix$ppm))]
  avg_int <- rowMeans(spectrum_mat[idx_match, , drop = FALSE], na.rm = TRUE)
  ord_desc <- order(avg_int, decreasing = TRUE)
  ppms_desc <- matched_ppms_clean[ord_desc]
  int_desc <- avg_int[ord_desc]
  baseline_candidates <- head(ppms_desc, min(max_pivots, length(ppms_desc)))
  unique_mppm <- sort(unique(molecule_ppm))

  all_results <- list()
  for (idx in seq_along(baseline_candidates)) {
    pivot <- baseline_candidates[idx]
    result <- run_single_pivot_pipeline(
      pivot_ppm = pivot, anchor_ppm = pivot, group_index = 1,
      matched_ppms = matched_ppms_clean, molecule_ppm = molecule_ppm,
      tol_match = tol_match,
      correlation_threshold = correlation_threshold,
      correlation_method = correlation_method,
      coverage_threshold_group = coverage_threshold_group,
      coverage_threshold_subchain = coverage_threshold_subchain,
      ratio_min = ratio_min, ransy_cv_threshold = ransy_cv_threshold,
      coverage_threshold_ransy = coverage_threshold_ransy,
      max_missing = max_missing, overlap_tol = 0.5
    )
    if (!is.null(result)) {
      result$confidence_level <- NA_character_
      result$n_groups_converged <- NA_integer_
      result$intersection_coverage <- NA_real_
      result$intersection_ppms <- numeric(0)
      result$winning_pivots <- pivot
      result$excluded_group <- NA_integer_
      result$exclusion_reason <- NA_character_
      result$note <- NA_character_
      all_results[[length(all_results) + 1]] <- result
    }
  }

  if (!length(all_results)) return(NULL)
  all_results
}

# -------------------------------
# 15) Intensity average
# -------------------------------
get_mean_intensity_for_ppm <- function(p) {
  if (!is.finite(p)) return(NA_real_)
  idx <- which.min(abs(spectrum_matrix$ppm - p))
  if (!length(idx) || is.na(idx) || idx < 1L) return(NA_real_)
  mean(spectrum_mat[idx, , drop = FALSE], na.rm = TRUE)
}

# -------------------------------
# 16) Default lines
# -------------------------------
creer_resultat_defaut <- function(src, ecfp_str, mppm, match_ppms) {
  data.table(
    Source = src, ECFP4 = ecfp_str,
    ppm_molecule = paste(unique(mppm), collapse = "|"),
    ppm_matrix = paste(match_ppms, collapse = "|"),
    groupe_correlation = "", baseline_ppm = NA_real_, sous_chaine = "",
    coverage_sub = 0, coverage_sub_ratio = 0, pivot_ppm = NA_real_,
    RANSY_signals = "", ppm_molecule_matched = "",
    matched_atom_indices = "",
    RANSY_mean_intensities = "",
    coverage_ransy = 0, RANSY_variabilite = "",
    avg_variabilite = 0, avg_error_ppm = NA_real_, intensity_ratios = "",
    confidence_level = NA_character_, n_groups_converged = NA_integer_,
    intersection_coverage = NA_real_, intersection_ppms = "",
    anchor_shifts = "", winning_pivots = "",
    pivot_group_index = NA_integer_, excluded_group = NA_integer_,
    exclusion_reason = NA_character_,
    error = "no_valid_results"
  )
}

# -------------------------------
# 17) Checkpoints
# -------------------------------
checkpoint_file <- file.path(work_dir, "checkpoint.rds")
resume_from_checkpoint <- function() {
  if (file.exists(checkpoint_file)) {
    log_message("Checkpoint detected - attempting resume...")
    out <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (is.null(out)) {
      log_message("Corrupted checkpoint - full restart")
      file.remove(checkpoint_file)
    }
    return(out)
  }
  NULL
}
save_checkpoint <- function(results_so_far, processed_sources, batch_idx) {
  if (!enable_checkpoints) return(invisible())
  tryCatch(
    {
      saveRDS(list(results = results_so_far, processed_sources = processed_sources,
                   batch_idx = batch_idx, timestamp = Sys.time()), checkpoint_file)
      log_message(paste("Checkpoint saved - batch", batch_idx))
    },
    error = function(e) log_message(paste("Checkpoint save error:", e$message))
  )
}

# -------------------------------
# 18) Batch pipeline
# -------------------------------
tryCatch(
  {
    sources <- unique(tableau_final$Source)
    total_sources <- length(sources)
    log_total_sources(total_sources)

    batch_size <- if (total_sources < 100) 10
                  else if (total_sources < 500) 25
                  else if (total_sources < 2000) 50
                  else 100

    log_message("=== PARAMETERS ===")
    log_message(paste("Total sources:", total_sources))
    log_message(paste("Batch size:", batch_size))
    log_message(paste("Max workers (nb_cores):", nb_cores))
    log_message(paste("Multi-pivot:", enable_multi_pivot))
    if (enable_multi_pivot) {
      log_message(paste("Max pivots/anchor:", max_pivots_per_anchor))
      log_message(paste("Max matches/shift:", max_matches_per_shift))
      log_message(paste("Intersection threshold:", intersection_coverage_threshold, "%"))
      log_message(paste("RANSY good threshold:", ransy_score_good_threshold))
    } else {
      log_message(paste("Legacy max pivots:", max_pivots))
    }

    checkpoint_data <- resume_from_checkpoint()
    if (!is.null(checkpoint_data)) {
      log_message("=== RESUME FROM CHECKPOINT ===")
      results_list <- checkpoint_data$results
      start_batch <- checkpoint_data$batch_idx + 1
    } else {
      log_message("=== FULL START ===")
      results_list <- list()
      start_batch <- 1
    }

    cores <- min(nb_cores, 3L, ceiling(max(1, total_sources) / 50))
    cores <- max(1L, cores)
    log_message(paste("Workers used (cluster):", cores))
    use_parallel <- cores > 1 && total_sources > 10
    cl <- NULL
    if (use_parallel) {
      cl <- tryCatch(
        {
          cluster <- parallel::makeCluster(cores, type = "PSOCK")
          doSNOW::registerDoSNOW(cluster)
          parallel::clusterEvalQ(cluster, {
            suppressMessages(library(data.table))
            options(datatable.optimize = 2, datatable.verbose = FALSE)
            data.table::setDTthreads(1L)
            gc()
            TRUE
          })
          parallel::clusterExport(cluster, varlist = c(
            "tableau_final", "spectrum_matrix", "spectrum_mat", "intensity_columns",
            "tol_match", "tol_correspondance", "ratio_min", "ransy_cv_threshold",
            "coverage_threshold_group", "coverage_threshold_subchain",
            "coverage_threshold_ransy", "max_missing_ransy", "correlation_threshold",
            "correlation_method", "max_pivots",
            "enable_multi_pivot", "max_pivots_per_anchor", "max_matches_per_shift",
            "intersection_coverage_threshold", "ransy_score_good_threshold",
            "anchor_min_distance_ppm",
            "triple_weight",
            "use_pqn", "use_csn", "normalization_applied",
            "USE_STOCSY", "USE_RANSY", "N_SAMPLES", "ANNOTATION_MODE",
            "intensity_ratio_annotation",
            "trouver_correspondances", "build_matched_ppms",
            "calculer_STOCSY", "gfilter_groups_par_coverage",
            "calculer_RANSY", "calculer_coverage_sub_ratio", "former_groupes_stocsy",
            "run_single_pivot_pipeline", "cross_sub_groups",
            "select_anchors", "exact_intersect", "get_overlap_groups_1D",
            "creer_resultat_defaut", "get_mean_intensity_for_ppm",
            "sanitize_for_csv", "fast_ppm_idx", "ppm_idx_lookup"
          ), envir = environment())
          log_message(paste("Parallel cluster initialized with", cores, "workers"))
          cluster
        },
        error = function(e) {
          log_message(paste("Cluster error:", e$message, "- switching to single-thread"))
          NULL
        }
      )
      if (is.null(cl)) use_parallel <- FALSE
    }
    `%op%` <- if (use_parallel) `%dopar%` else `%do%`
    if (!is.null(cl)) {
      on.exit({ tryCatch({ parallel::stopCluster(cl); log_message("Cluster closed cleanly") },
                          error = function(e) {}) }, add = TRUE)
    }

    total_batches <- ceiling(total_sources / batch_size)
    if (!is.finite(total_batches) || total_batches < 1) total_batches <- 1L
    log_progress(0, 100, "Initializing 13C")

    for (batch_idx in seq(from = start_batch, to = total_batches)) {
      batch_start_time <- Sys.time()
      start_idx <- (batch_idx - 1) * batch_size + 1
      end_idx <- min(batch_idx * batch_size, total_sources)
      batch_sources <- if (total_sources > 0) sources[start_idx:end_idx] else character(0)
      pct_start <- round(((start_idx - 1) / max(1, total_sources)) * 100, 1)

      log_message(paste("=== BATCH", batch_idx, "/", total_batches, "==="))
      log_message(paste("Sources:", length(batch_sources), "(", start_idx, "-", end_idx, ")"))
      log_progress(pct_start, 100, sprintf("Batch %d/%d — sources %d-%d", batch_idx, total_batches, start_idx, end_idx))
      log_memory(paste("batch start", batch_idx))

      batch_results <- tryCatch(
        {
          foreach(i = seq_along(batch_sources), .packages = "data.table",
                  .combine = "c", .init = list(), .errorhandling = "pass") %op% {
            src <- batch_sources[i]
            tryCatch(
              {
                dt_raw <- tableau_final[Source == src]
                mppm <- dt_raw$shift_value
                ecfp_vals <- unique(dt_raw$ECFP4)
                ecfp_str <- paste(ecfp_vals, collapse = "|")
                if (length(mppm) == 0) return(list(creer_resultat_defaut(src, ecfp_str, numeric(), numeric())))

                # Use filtered matching for multi-pivot, full matching for legacy
                if (enable_multi_pivot) {
                  match_info <- build_matched_ppms(mppm, tol_correspondance, max_matches_per_shift, min_spacing = 0.01)
                  match_ppms <- match_info$matched_ppms
                  matches_dt <- match_info$matches_dt
                } else {
                  matches_dt <- trouver_correspondances(mppm, tol_correspondance)
                  match_ppms <- sort(unique(matches_dt$ppm_matrix))
                }

                if (length(match_ppms) == 0) return(list(creer_resultat_defaut(src, ecfp_str, mppm, numeric())))

                all_results <- former_groupes_stocsy(
                  matched_ppms = match_ppms, molecule_ppm = mppm, tol_match = tol_match,
                  correlation_threshold = correlation_threshold,
                  correlation_method = correlation_method,
                  coverage_threshold_group = coverage_threshold_group,
                  coverage_threshold_subchain = coverage_threshold_subchain,
                  ratio_min = ratio_min, ransy_cv_threshold = ransy_cv_threshold,
                  coverage_threshold_ransy = coverage_threshold_ransy,
                  max_missing = max_missing_ransy, max_pivots = max_pivots
                )

                if (is.null(all_results) || length(all_results) == 0) {
                  return(list(data.table(
                    Source = src, ECFP4 = ecfp_str,
                    ppm_molecule = paste(unique(mppm), collapse = "|"),
                    ppm_matrix = paste(match_ppms, collapse = "|"),
                    groupe_correlation = "", baseline_ppm = NA_real_,
                    sous_chaine = "", coverage_sub = 0, coverage_sub_ratio = 0,
                    pivot_ppm = NA_real_, RANSY_signals = "", ppm_molecule_matched = "",
                    matched_atom_indices = "", RANSY_mean_intensities = "",
                    coverage_ransy = 0, RANSY_variabilite = "", avg_variabilite = 0,
                    avg_error_ppm = NA_real_, intensity_ratios = "",
                    confidence_level = NA_character_, n_groups_converged = NA_integer_,
                    intersection_coverage = NA_real_, intersection_ppms = "",
                    anchor_shifts = "", winning_pivots = "",
                    pivot_group_index = NA_integer_, excluded_group = NA_integer_,
                    exclusion_reason = NA_character_,
                    error = "no_ransy_valid"
                  )))
                }

                # Compute anchor_shifts string once
                anchor_str <- ""
                if (enable_multi_pivot) {
                  anchor_info <- select_anchors(mppm, anchor_min_distance_ppm)
                  if (!is.null(anchor_info$anchors)) anchor_str <- paste(anchor_info$anchors, collapse = "|")
                }

                rows <- rbindlist(lapply(all_results, function(best) {
                  coverage_sub_ratio <- calculer_coverage_sub_ratio(best$sub_chain, mppm, tol_match)
                  ransy_all <- as.numeric(best$ransy_ppms)
                  ransy_means <- vapply(ransy_all, get_mean_intensity_for_ppm, numeric(1))

                  pairs <- best$matched_pairs
                  if (is.null(pairs) || !nrow(pairs)) pairs <- data.table(ppm_molecule = numeric(0), ppm_ransy = numeric(0))
                  pairs <- pairs[!is.na(ppm_ransy)]
                  mppm_matched <- pairs$ppm_molecule

                  matched_atom_indices <- ""
                  if (length(mppm_matched) && "C_index" %in% names(dt_raw)) {
                    idx_parts <- vapply(mppm_matched, function(sv) {
                      cand <- dt_raw[abs(shift_value - sv) <= tol_match]
                      if (nrow(cand)) paste(unique(cand$C_index), collapse = "~") else ""
                    }, character(1))
                    matched_atom_indices <- paste(idx_parts[idx_parts != ""], collapse = "|")
                  }

                  # Extract new fields safely
                  conf_level <- if (!is.null(best$confidence_level)) best$confidence_level else NA_character_
                  n_conv <- if (!is.null(best$n_groups_converged)) as.integer(best$n_groups_converged) else NA_integer_
                  inter_cov <- if (!is.null(best$intersection_coverage)) best$intersection_coverage else NA_real_
                  inter_ppms_str <- if (!is.null(best$intersection_ppms) && length(best$intersection_ppms) > 0) {
                    paste(best$intersection_ppms, collapse = "|")
                  } else ""
                  w_pivots_str <- if (!is.null(best$winning_pivots)) {
                    paste(best$winning_pivots, collapse = "|")
                  } else ""
                  grp_idx <- if (!is.null(best$group_index)) as.integer(best$group_index) else NA_integer_
                  excl_g <- if (!is.null(best$excluded_group) && !all(is.na(best$excluded_group))) {
                    as.integer(best$excluded_group[1])
                  } else NA_integer_
                  excl_r <- if (!is.null(best$exclusion_reason) && !all(is.na(best$exclusion_reason))) {
                    as.character(best$exclusion_reason)
                  } else NA_character_

                  data.table(
                    Source = src, ECFP4 = ecfp_str,
                    ppm_molecule = paste(unique(mppm), collapse = "|"),
                    ppm_matrix = paste(sort(unique(matches_dt$ppm_matrix)), collapse = "|"),
                    groupe_correlation = paste(best$groupe_ppms, collapse = "-"),
                    baseline_ppm = best$pivot_ppm,
                    sous_chaine = paste(best$sub_chain, collapse = "|"),
                    coverage_sub = best$coverage_sub,
                    coverage_sub_ratio = coverage_sub_ratio,
                    pivot_ppm = best$pivot_ppm,
                    RANSY_signals = paste(ransy_all, collapse = "|"),
                    ppm_molecule_matched = paste(mppm_matched, collapse = "|"),
                    matched_atom_indices = matched_atom_indices,
                    RANSY_mean_intensities = sanitize_for_csv(paste(round(ransy_means, 6), collapse = ";")),
                    coverage_ransy = best$coverage_ransy,
                    RANSY_variabilite = sanitize_for_csv(best$RANSY_variabilite),
                    avg_variabilite = best$avg_cv,
                    avg_error_ppm = ifelse(is.null(best$avg_error_ppm) || is.na(best$avg_error_ppm), NA_real_, best$avg_error_ppm),
                    intensity_ratios = sanitize_for_csv(best$intensity_ratios),
                    confidence_level = conf_level,
                    n_groups_converged = n_conv,
                    intersection_coverage = inter_cov,
                    intersection_ppms = inter_ppms_str,
                    anchor_shifts = anchor_str,
                    winning_pivots = w_pivots_str,
                    pivot_group_index = grp_idx,
                    excluded_group = excl_g,
                    exclusion_reason = excl_r,
                    error = NA_character_
                  )
                }), fill = TRUE)

                list(rows)
              },
              error = function(e) {
                log_message(sprintf("Error for source %s: %s", src, e$message))
                list(creer_resultat_defaut(
                  src = src,
                  ecfp_str = tryCatch(paste(unique(tableau_final[Source == src]$ECFP4), collapse = "|"), error = function(e2) ""),
                  mppm = tryCatch(tableau_final[Source == src]$shift_value, error = function(e2) numeric()),
                  match_ppms = numeric()
                ))
              }
            )
          }
        },
        error = function(e) {
          log_error(paste("Batch error", batch_idx, ":", e$message))
          list()
        }
      )

      if (length(batch_results) > 0) results_list <- c(results_list, batch_results)

      if (batch_idx %% 2 == 0) { gc(); log_memory(paste("after batch", batch_idx)) }
      if (enable_checkpoints && batch_idx %% 5 == 0) {
        processed_sources_so_far <- if (total_sources > 0) sources[1:end_idx] else character(0)
        save_checkpoint(results_list, processed_sources_so_far, batch_idx)
      }

      batch_duration <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
      log_message(paste("Batch", batch_idx, "completed in", round(batch_duration, 1), "seconds"))
      pct_end <- round((end_idx / max(1, total_sources)) * 100, 1)
      if (batch_idx < total_batches) {
        remaining_batches <- total_batches - batch_idx
        eta_minutes <- round((remaining_batches * batch_duration) / 60, 1)
        log_message(paste("ETA:", eta_minutes, "minutes remaining"))
      }
      log_progress(pct_end, 100, sprintf("Batch %d/%d completed (%d/%d sources)", batch_idx, total_batches, end_idx, total_sources))
    }

    # FINAL ASSEMBLY
    log_progress(95, 100, "Final assembly...")
    log_message("=== FINAL ASSEMBLY ===")
    log_memory("before assembly")
    if (length(results_list) > 0) {
      resultats <- data.table::rbindlist(results_list, fill = TRUE)

      resultats[, ransy_score := ifelse(is.na(coverage_ransy), 0, coverage_ransy)]

      # Penalize molecules with only 1 matched signal — trivial match, unreliable
      resultats[, n_matched := sapply(ppm_molecule_matched, function(x) {
        if (is.na(x) || x == "") return(0L)
        length(strsplit(as.character(x), "\\|")[[1]])
      })]
      resultats[n_matched <= 1, ransy_score := 0]

      # Multi-criteria sort: confidence > ransy_score (= coverage_ransy)
      confidence_order <- c("HIGH" = 4, "MEDIUM_HIGH" = 3, "MEDIUM" = 2, "LOW" = 1)
      resultats[, conf_num := fifelse(is.na(confidence_level), 0L, as.integer(confidence_order[confidence_level]))]
      data.table::setorder(resultats, Source, -conf_num, -ransy_score)
      log_message("Ranking by confidence_level then ransy_score (= coverage_ransy)")

      best_per_source13C <- resultats[, .SD[1], by = Source]
      best_per_source13C[, Source := sub("^NMR_", "", Source)]
      best_per_source13C[, Source := sub("\\.sdf$", "", Source)]
      # Remove temp columns
      best_per_source13C[, c("conf_num", "n_matched") := NULL]
    } else {
      log_error("No results to assemble")
      quit(status = 1)
    }

    output_file <- file.path(work_dir, "13C_With_Ransy_and_ECFP4.csv")
    data.table::fwrite(best_per_source13C, output_file, na = "NA", sep = ";", quote = TRUE, bom = TRUE)

    if (enable_checkpoints && file.exists(checkpoint_file)) {
      file.remove(checkpoint_file)
      log_message("Final checkpoint removed")
    }

    log_message("=== FINAL 13C RESULTS ===")
    log_message(paste("File saved:", output_file))
    log_message(paste("Sources processed:", nrow(best_per_source13C)))

    resultats_valides <- best_per_source13C[is.na(error)]
    resultats_erreur <- best_per_source13C[!is.na(error)]
    log_message(paste("Success:", nrow(resultats_valides)))
    log_message(paste("Errors:", nrow(resultats_erreur)))
    if (nrow(resultats_erreur) > 0) {
      error_types <- table(resultats_erreur$error)
      for (err_type in names(error_types)) log_message(paste("  -", err_type, ":", error_types[[err_type]]))
    }

    # Confidence level distribution
    if (nrow(resultats_valides) > 0 && "confidence_level" %in% names(resultats_valides)) {
      conf_dist <- resultats_valides[!is.na(confidence_level), .N, by = confidence_level]
      if (nrow(conf_dist) > 0) {
        log_message("=== CONFIDENCE DISTRIBUTION ===")
        for (r in seq_len(nrow(conf_dist))) {
          log_message(paste("  -", conf_dist$confidence_level[r], ":", conf_dist$N[r]))
        }
      }
    }

    log_memory("final")

    log_result(list(
      total_sources = length(sources),
      successful_sources = nrow(resultats_valides),
      error_sources = nrow(resultats_erreur),
      output_file = output_file,
      annotation_type = "13C",
      version = "multi_pivot_stocsy_v2",
      batches_processed = ceiling(length(sources) / batch_size),
      memory_optimized = TRUE,
      multi_pivot_enabled = enable_multi_pivot,
      annotation_mode = if (exists("ANNOTATION_MODE")) ANNOTATION_MODE else "unknown",
      n_samples = if (exists("N_SAMPLES")) N_SAMPLES else NA,
      stocsy_enabled = if (exists("USE_STOCSY")) USE_STOCSY else TRUE,
      ransy_enabled = if (exists("USE_RANSY")) USE_RANSY else TRUE,
      normalization_applied = if (exists("normalization_applied")) normalization_applied else "none",
      summary = list(
        avg_coverage_ransy = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$coverage_ransy, na.rm = TRUE), 2) else 0,
        avg_coverage_sub   = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$coverage_sub, na.rm = TRUE), 2) else 0,
        avg_variabilite    = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$avg_variabilite, na.rm = TRUE), 4) else 0,
        avg_error_ppm      = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$avg_error_ppm, na.rm = TRUE), 4) else 0
      )
    ))

    log_progress(100, 100, "13C completed")
  },
  error = function(e) {
    log_error(paste("Critical error in pipeline:", e$message))
    quit(status = 1)
  }
)

log_message("=== 13C ANNOTATION WITH MULTI-PIVOT STOCSY V2 (OK) ===")