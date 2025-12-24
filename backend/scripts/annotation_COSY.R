#!/usr/bin/env Rscript

# =====================================================================
# Annotation COSY 2D NMR — With STOCSY instead of correlations (OPTIMIZED VERSION)
# Pipeline: Match → Intensity sort → pivots → STOCSY → Group coverage → Substring → RANSY
# Compatible with existing frontend/backend (same parameter names)
# =====================================================================

suppressWarnings(suppressMessages({}))

# -------------------------------
# 0) Arguments & basic helpers
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: Rscript annotation_COSY_stocsy.R <work_dir> <prediction_COSY_csv_file> <spectrum_matrix_csv_file> [params_json]\n")
    quit(status = 1)
}
work_dir <- args[1]
prediction_COSY_file <- args[2]
spectrum_matrix_file <- args[3]

# Default parameters (backend compatibility)
params <- list(
    tol_correspondance_F1 = 0.2,
    tol_correspondance_F2 = 0.2,
    tol_match_F1 = 0.2,
    tol_match_F2 = 0.2,
    correlation_threshold = 0.80, # |r| threshold for STOCSY
    coverage_threshold_group = 0,
    coverage_threshold_subchain = 0,
    ratio_min = 0.01,
    ransy_cv_threshold = 0.15,
    coverage_threshold_ransy = 0,
    max_missing_ransy = 50,
    nb_cores = 6,
    progress_style = 3,
    batch_size = 1,
    memory_cleanup_interval = 50,
    enable_checkpoints = TRUE,
    max_memory_mb = 2048,
    max_pivots = 100,
    alpha_cost = 0.3
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
    # memory.size() no longer supported in R >= 4.2, use NA directly
    # if (.Platform$OS.type == "windows") m <- tryCatch(utils::memory.size(), error = function(e) NA_real_)
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
# 3) Merge parameters
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
}

nb_cores <- sanitize_threads(params$nb_cores)
if (is.na(nb_cores)) nb_cores <- 1L
available <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (is.na(available) || available < 1L) available <- 1L
nb_cores <- max(1L, min(as.integer(nb_cores), 3L, as.integer(available)))

# -------------------------------
# 4) Libs + threads data.table
# -------------------------------
suppressMessages({
    ok <- TRUE
    tryCatch(
        {
            suppressWarnings(suppressMessages({
                library(data.table)
                library(doSNOW)
                library(foreach)
                library(parallel)
                library(jsonlite)
            }))
        },
        error = function(e) {
            ok <<- FALSE
            log_error(paste("ERROR: Unable to load libraries:", e$message))
        }
    )
    if (!ok) quit(status = 1)
    data.table::setDTthreads(nb_cores)
})
if (.Platform$OS.type == "windows") {
    lim <- suppressWarnings(as.numeric(params$max_memory_mb))
    if (is.finite(lim) && lim > 256) try(suppressWarnings(memory.limit(size = lim)), silent = TRUE)
}
options(datatable.optimize = 2, datatable.verbose = FALSE, datatable.alloccol = 1024)
gc()

# -------------------------------
# 5) Parameter extraction
# -------------------------------
tol_correspondance_F1 <- as.numeric(params$tol_correspondance_F1)
tol_correspondance_F2 <- as.numeric(params$tol_correspondance_F2)
tol_match_F1 <- as.numeric(params$tol_match_F1)
tol_match_F2 <- as.numeric(params$tol_match_F2)
correlation_threshold <- as.numeric(params$correlation_threshold)
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
if (is.na(max_pivots) || max_pivots < 1) max_pivots <- 200L
alpha_cost <- as.numeric(params$alpha_cost)

# Global reservoir
intensity_columns <- character(0)
spectrum_mat <- NULL # global 2D intensity matrix

# -------------------------------
# 6) Data reading
# -------------------------------
tryCatch(
    {
        log_message("=== COSY ANNOTATION START — STOCSY ===")
        log_message(paste("Working directory:", work_dir))
        log_message(paste("data.table threads:", data.table::getDTthreads()))
        log_memory("initial")

        setwd(work_dir)
        log_message("Optimized file reading...")

        prediction_COSY <- data.table::fread(prediction_COSY_file, dec = ".", sep = ",", header = TRUE, encoding = "UTF-8", showProgress = FALSE, verbose = FALSE)
        log_memory("after reading COSY predictions")

        spectrum_matrix <- data.table::fread(spectrum_matrix_file, dec = ".", sep = ",", header = TRUE, encoding = "UTF-8", showProgress = FALSE, verbose = FALSE)
        log_memory("after reading matrix")

        log_message(paste("COSY predictions:", nrow(prediction_COSY), "rows"))
        log_message(paste("Spectrum matrix:", nrow(spectrum_matrix), "rows"))

        setDT(prediction_COSY)
        setDT(spectrum_matrix)

        # Mapping compat
        if (!"Source" %in% names(prediction_COSY) && "Filename" %in% names(prediction_COSY)) prediction_COSY[, Source := Filename]
        if (!"shift_value_F1" %in% names(prediction_COSY) && "Delta1" %in% names(prediction_COSY)) prediction_COSY[, shift_value_F1 := Delta1]
        if (!"shift_value_F2" %in% names(prediction_COSY) && "Delta2" %in% names(prediction_COSY)) prediction_COSY[, shift_value_F2 := Delta2]
        if (!("ECFP4" %in% names(prediction_COSY))) prediction_COSY[, ECFP4 := ""]

        prediction_COSY[, shift_value_F1 := as.numeric(shift_value_F1)]
        prediction_COSY[, shift_value_F2 := as.numeric(shift_value_F2)]
        prediction_COSY <- prediction_COSY[!is.na(shift_value_F1) & !is.na(shift_value_F2)]

        if (!all(c("ppm_F1", "ppm_F2") %in% names(spectrum_matrix))) stop("Missing 'ppm_F1' and/or 'ppm_F2' columns.")
        spectrum_matrix[, ppm_F1 := as.numeric(ppm_F1)]
        spectrum_matrix[, ppm_F2 := as.numeric(ppm_F2)]
        spectrum_matrix <- spectrum_matrix[!is.na(ppm_F1) & !is.na(ppm_F2)]
        setorder(spectrum_matrix, ppm_F1, ppm_F2)

        spectrum_matrix[, ppm_matrix := paste(ppm_F1, ppm_F2, sep = "|")]
        setkey(spectrum_matrix, ppm_matrix)

        intensity_columns <<- setdiff(names(spectrum_matrix), c("ppm_F1", "ppm_F2", "ppm_matrix"))
        if (length(intensity_columns) == 0) stop("No intensity columns.")
        log_message(paste("Intensity columns:", length(intensity_columns)))

        # Global intensity matrix (all COSY rows x samples)
        spectrum_mat <<- as.matrix(spectrum_matrix[, ..intensity_columns])

        tableau_final <- prediction_COSY
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
# 8) Matching (ALL signals within tolerance)
# -------------------------------
trouver_correspondances_2D <- function(shifts, tol_F1, tol_F2) {
    if (!nrow(shifts) || !nrow(spectrum_matrix)) {
        return(data.table(ppm_matrix = character(), distance = numeric()))
    }
    temp <- copy(shifts)
    spec <- spectrum_matrix[, .(ppm_F1, ppm_F2, ppm_matrix,
        F1_min = ppm_F1 - tol_F1, F1_max = ppm_F1 + tol_F1,
        F2_min = ppm_F2 - tol_F2, F2_max = ppm_F2 + tol_F2
    )]
    res <- spec[temp,
        on = .(
            F1_min <= shift_value_F1, F1_max >= shift_value_F1,
            F2_min <= shift_value_F2, F2_max >= shift_value_F2
        ),
        allow.cartesian = TRUE, nomatch = 0,
        .(ppm_matrix, shift_value_F1, shift_value_F2,
            distance = sqrt((shift_value_F1 - ppm_F1)^2 + (shift_value_F2 - ppm_F2)^2)
        )
    ]
    if (!nrow(res)) {
        return(res)
    }
    setorder(res, distance)
    res
}

# -------------------------------
# 9) 2D STOCSY (OPTIMIZED)
# -------------------------------
calculer_STOCSY_2D <- function(ppms, intensity_cols, driver_ppm, r_min) {
    if (length(ppms) == 0L) {
        return(data.table(ppm_matrix = character(), r = numeric()))
    }

    idx <- match(ppms, spectrum_matrix$ppm_matrix)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
        return(data.table(ppm_matrix = character(), r = numeric()))
    }

    mat <- spectrum_mat[idx, , drop = FALSE]
    row_ppm <- spectrum_matrix$ppm_matrix[idx]

    driver_idx <- which(row_ppm == driver_ppm)[1]
    if (!length(driver_idx) || is.na(driver_idx) || driver_idx < 1L) {
        return(data.table(ppm_matrix = character(), r = numeric()))
    }

    driver_v <- mat[driver_idx, ]
    valid_cols <- which(is.finite(driver_v) & !is.na(driver_v))

    if (length(valid_cols) < 3L) {
        return(data.table(ppm_matrix = character(), r = numeric()))
    }

    mat2 <- mat[, valid_cols, drop = FALSE]
    drv <- driver_v[valid_cols]

    r_vec <- apply(mat2, 1, function(x) {
        suppressWarnings(cor(as.numeric(x), drv, use = "pairwise.complete.obs", method = "pearson"))
    })
    r_vec[!is.finite(r_vec)] <- NA_real_

    dt <- data.table(ppm_matrix = row_ppm, r = r_vec)
    dt <- dt[is.finite(r) & abs(r) >= r_min]

    if (!nrow(dt)) {
        return(data.table(ppm_matrix = character(), r = numeric()))
    }
    dt[order(ppm_matrix)]
}

# -------------------------------
# 10) Group filtering by coverage (FIXED)
# -------------------------------
gfilter_groups_par_coverage_2D <- function(groups, molecule_shifts, tol_F1, tol_F2, cov_thresh) {
    out <- list()
    uniq_ms <- unique(molecule_shifts[, .(shift_value_F1, shift_value_F2)])
    if (!nrow(uniq_ms)) {
        return(out)
    }
    for (grp in groups) {
        coords <- spectrum_matrix[ppm_matrix %in% grp, .(ppm_F1, ppm_F2)]
        if (!nrow(coords)) next
        covered <- vapply(seq_len(nrow(uniq_ms)), function(i) {
            m1 <- uniq_ms$shift_value_F1[i]
            m2 <- uniq_ms$shift_value_F2[i]
            any(abs(coords$ppm_F1 - m1) <= tol_F1 & abs(coords$ppm_F2 - m2) <= tol_F2)
        }, logical(1))
        cov <- if (length(covered)) mean(covered) * 100 else 0
        if (cov >= cov_thresh) out[[length(out) + 1]] <- grp
    }
    out
}

# -------------------------------
# 11) 2D RANSY (OPTIMIZED)
# -------------------------------
calculer_RANSY_2D <- function(ppms, intensity_cols, cv_th, driver_ppm) {
    if (length(ppms) == 0L) {
        return(NULL)
    }

    idx <- match(ppms, spectrum_matrix$ppm_matrix)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
        return(NULL)
    }

    mat <- spectrum_mat[idx, , drop = FALSE]
    row_ppm <- spectrum_matrix$ppm_matrix[idx]

    driver_idx <- which(row_ppm == driver_ppm)[1]
    if (!length(driver_idx) || is.na(driver_idx) || driver_idx < 1L) {
        return(NULL)
    }
    driver_v <- mat[driver_idx, ]
    valid_cols <- which(driver_v > 0 & !is.na(driver_v))
    if (length(valid_cols) < 2L) {
        return(NULL)
    }

    mat2 <- mat[, valid_cols, drop = FALSE]
    ratio_mat <- sweep(mat2, 2, driver_v[valid_cols], "/")

    cvs <- apply(ratio_mat, 1, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) {
            return(Inf)
        }
        m <- mean(x)
        s <- sd(x)
        if (!is.finite(m) || m == 0) Inf else s / m
    })
    data.table(ppm_matrix = row_ppm, cv = cvs)[cv <= cv_th]
}

# -------------------------------
# 12) Substring coverage (FIXED)
# -------------------------------
calculer_coverage_sub_ratio_2D <- function(sub_chain, molecule_shifts, tol_F1, tol_F2) {
    if (length(sub_chain) == 0 || !nrow(molecule_shifts)) {
        return(0)
    }
    uniq_ms <- unique(molecule_shifts[, .(shift_value_F1, shift_value_F2)])
    matched_u <- rep(FALSE, nrow(uniq_ms))
    avail_sub <- spectrum_matrix[ppm_matrix %in% sub_chain, .(ppm_F1, ppm_F2, ppm_matrix)]
    for (i in seq_len(nrow(uniq_ms))) {
        m1 <- uniq_ms$shift_value_F1[i]
        m2 <- uniq_ms$shift_value_F2[i]
        cand <- avail_sub[abs(ppm_F1 - m1) <= tol_F1 & abs(ppm_F2 - m2) <= tol_F2]
        if (nrow(cand) > 0) {
            cand[, dist := sqrt((ppm_F1 - m1)^2 + (ppm_F2 - m2)^2)]
            best <- cand[which.min(dist), ppm_matrix]
            avail_sub <- avail_sub[ppm_matrix != best]
            matched_u[i] <- TRUE
        }
    }
    sum(matched_u) / nrow(uniq_ms) * 100
}

# -------------------------------
# 13) STOCSY pipeline for each molecule (OPTIMIZED)
# -------------------------------
former_groupes_stocsy_2D <- function(matched_ppms, molecule_shifts, tol_match_F1, tol_match_F2,
                                     correlation_threshold, coverage_threshold_group,
                                     coverage_threshold_subchain, ratio_min,
                                     ransy_cv_threshold, coverage_threshold_ransy,
                                     max_missing, max_pivots = 200, alpha = 0.3) {
    if (length(matched_ppms) == 0 || !nrow(molecule_shifts)) {
        return(NULL)
    }

    # Step 1: average intensities via spectrum_mat
    idx_match <- match(matched_ppms, spectrum_matrix$ppm_matrix)
    idx_valid <- which(!is.na(idx_match))
    if (!length(idx_valid)) {
        return(NULL)
    }

    matched_ppms_clean <- matched_ppms[idx_valid]
    avg_int <- rowMeans(spectrum_mat[idx_match[idx_valid], , drop = FALSE], na.rm = TRUE)

    ord_desc <- order(avg_int, decreasing = TRUE, na.last = TRUE)
    ppms_desc <- matched_ppms_clean[ord_desc]
    int_desc <- avg_int[ord_desc]

    keep <- which(!is.na(int_desc) & int_desc > 0)
    if (!length(keep)) {
        return(NULL)
    }
    ppms_desc <- ppms_desc[keep]
    int_desc <- int_desc[keep]

    # Step 2: pivots
    baseline_candidates <- head(ppms_desc, min(max_pivots, length(ppms_desc)))

    # Target deduplication F1/F2
    unique_mppm <- unique(molecule_shifts[, .(shift_value_F1, shift_value_F2)])

    # Step 3: STOCSY for each pivot
    groups_stocsy <- list()
    for (pivot in baseline_candidates) {
        dt_stocsy <- calculer_STOCSY_2D(
            ppms = matched_ppms_clean,
            intensity_cols = intensity_columns,
            driver_ppm = pivot,
            r_min = correlation_threshold
        )
        if (is.null(dt_stocsy) || !nrow(dt_stocsy)) next
        stocsy_ppms <- sort(unique(dt_stocsy$ppm_matrix))
        groups_stocsy[[length(groups_stocsy) + 1]] <- stocsy_ppms
    }

    if (length(groups_stocsy) == 0) {
        return(NULL)
    }

    # Step 4: group coverage filtering
    if (coverage_threshold_group > 0) {
        groups_stocsy <- gfilter_groups_par_coverage_2D(groups_stocsy, unique_mppm, tol_match_F1, tol_match_F2, coverage_threshold_group)
    }
    if (length(groups_stocsy) == 0) {
        return(NULL)
    }

    # Step 5: substrings + RANSY
    all_results <- list()

    for (idx in seq_along(groups_stocsy)) {
        grp <- groups_stocsy[[idx]]
        pivot <- baseline_candidates[idx]
        pivot_int <- int_desc[ppms_desc == pivot][1]

        if (is.na(pivot_int) || pivot_int == 0) next

        # Group intensities via spectrum_mat
        idx_grp <- match(grp, spectrum_matrix$ppm_matrix)
        idx_grp_valid <- which(!is.na(idx_grp))
        if (!length(idx_grp_valid)) next
        grp_clean <- grp[idx_grp_valid]
        grp_int <- rowMeans(spectrum_mat[idx_grp[idx_grp_valid], , drop = FALSE], na.rm = TRUE)

        ratios <- grp_int / pivot_int
        idx_valide <- which(grp_int <= pivot_int & ratios >= ratio_min & !is.na(grp_int))

        if (!length(idx_valide)) next
        sub_chain_vec <- sort(unique(grp_clean[idx_valide]))

        # Coverage sub
        coverage_sub_final <- calculer_coverage_sub_ratio_2D(sub_chain_vec, unique_mppm, tol_match_F1, tol_match_F2)
        if (!is.finite(coverage_sub_final)) coverage_sub_final <- 0
        if (coverage_sub_final < coverage_threshold_subchain) next

        # RANSY
        dt_ransy <- calculer_RANSY_2D(sub_chain_vec, intensity_columns, ransy_cv_threshold, pivot)
        if (is.null(dt_ransy) || !nrow(dt_ransy)) next

        avail <- merge(dt_ransy, spectrum_matrix[, .(ppm_matrix, ppm_F1, ppm_F2)], by = "ppm_matrix")
        if (!nrow(avail)) next

        max_cv_all <- max(avail$cv, na.rm = TRUE)
        if (!is.finite(max_cv_all) || max_cv_all == 0) max_cv_all <- 1

        # RANSY pairings (1-1) with distance+CV cost
        pairs_dt <- data.table(mol_F1 = numeric(0), mol_F2 = numeric(0), ppm_ransy = character(0))

        for (i in seq_len(nrow(unique_mppm))) {
            m1 <- unique_mppm$shift_value_F1[i]
            m2 <- unique_mppm$shift_value_F2[i]
            cand <- avail[abs(ppm_F1 - m1) <= tol_match_F1 & abs(ppm_F2 - m2) <= tol_match_F2]

            if (nrow(cand) > 0) {
                cand[, dist_norm := sqrt(((ppm_F1 - m1) / tol_match_F1)^2 + ((ppm_F2 - m2) / tol_match_F2)^2)]
                cand[, cost := alpha * (cv / max_cv_all) + (1 - alpha) * dist_norm]
                best_r <- cand[which.min(cost), ppm_matrix]
                pairs_dt <- rbind(pairs_dt, data.table(mol_F1 = m1, mol_F2 = m2, ppm_ransy = best_r))
                avail <- avail[ppm_matrix != best_r]
            }
        }

        matched_unique <- unique(pairs_dt$ppm_ransy)
        coverage_ransy <- if (nrow(unique_mppm) > 0) length(matched_unique) / nrow(unique_mppm) * 100 else 0
        if (coverage_ransy > coverage_sub_final) coverage_ransy <- coverage_sub_final

        missing_cnt <- nrow(unique_mppm) - length(matched_unique)
        if (missing_cnt > max_missing || coverage_ransy < coverage_threshold_ransy) next

        # Intensities in substring (via spectrum_mat)
        idx_sub <- match(sub_chain_vec, spectrum_matrix$ppm_matrix)
        idx_sub_valid <- which(!is.na(idx_sub))
        if (!length(idx_sub_valid)) next
        sub_chain_clean <- sub_chain_vec[idx_sub_valid]
        avg_mat_all <- rowMeans(spectrum_mat[idx_sub[idx_sub_valid], , drop = FALSE], na.rm = TRUE)
        max_int_all <- max(avg_mat_all, na.rm = TRUE)
        if (!is.finite(max_int_all) || max_int_all == 0) max_int_all <- 1
        intensity_ratios_str <- paste(round(avg_mat_all / max_int_all, 3), collapse = "|")
        ransy_var_str <- paste(paste(dt_ransy$ppm_matrix, round(dt_ransy$cv, 4), sep = ":"), collapse = "|")

        all_results[[length(all_results) + 1]] <- list(
            base_ppm          = pivot,
            sub_chain         = sub_chain_clean,
            coverage_sub      = coverage_sub_final,
            intensity_ratios  = intensity_ratios_str,
            pivot_ppm         = pivot,
            RANSY_signals_all = dt_ransy$ppm_matrix,
            matched_for_cov   = matched_unique,
            matched_pairs     = pairs_dt,
            coverage_ransy    = coverage_ransy,
            RANSY_variabilite = ransy_var_str,
            avg_variabilite   = round(mean(dt_ransy[ppm_matrix %in% matched_unique, cv], na.rm = TRUE), 4),
            groupe_ppms       = grp
        )
    }

    if (!length(all_results)) {
        return(NULL)
    }
    all_results
}

# -------------------------------
# 14) 2D intensity average (OPTIMIZED)
# -------------------------------
get_mean_intensity_for_ppm2d <- function(pm) {
    if (!is.character(pm) || length(pm) != 1L) {
        return(NA_real_)
    }
    idx <- match(pm, spectrum_matrix$ppm_matrix)
    if (is.na(idx) || idx < 1L) {
        return(NA_real_)
    }
    mean(spectrum_mat[idx, , drop = FALSE], na.rm = TRUE)
}

# -------------------------------
# 15) Default lines
# -------------------------------
creer_resultat_defaut_2D <- function(src, ecfp_str, dt_mol, match_ppms) {
    data.table(
        Source = src, ECFP4 = ecfp_str,
        ppm_molecule = if (nrow(dt_mol)) paste(dt_mol$shift_value_F1, dt_mol$shift_value_F2, sep = "|", collapse = "||") else "",
        ppm_matrix = paste(match_ppms, collapse = "||"),
        groupe_correlation = "",
        baseline_ppm = NA_character_, sous_chaine = "",
        coverage_sub = 0, coverage_sub_ratio = 0,
        pivot_ppm = NA_character_,
        RANSY_signals = "", RANSY_mean_intensities = "",
        coverage_ransy = 0, RANSY_variabilite = "",
        avg_variabilite = 0, intensity_ratios = "",
        ppm_molecule_matched = "",
        error = "no_valid_results"
    )
}

# -------------------------------
# 16) Checkpoints
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
    if (!enable_checkpoints) {
        return(invisible())
    }
    tryCatch(
        {
            saveRDS(list(results = results_so_far, processed_sources = processed_sources, batch_idx = batch_idx, timestamp = Sys.time()), checkpoint_file)
            log_message(paste("Checkpoint saved - batch", batch_idx))
        },
        error = function(e) {
            log_message(paste("Checkpoint save error:", e$message))
        }
    )
}

# -------------------------------
# 17) Batch pipeline
# -------------------------------
tryCatch(
    {
        sources <- unique(tableau_final$Source)
        total_sources <- length(sources)
        log_total_sources(total_sources)
        
        # Adaptive batch size according to molecule count
        batch_size <- if (total_sources < 100) {
          10
        } else if (total_sources < 500) {
          25
        } else if (total_sources < 2000) {
          50
        } else {
          100
        }

        log_message("=== PARAMETERS ===")
        log_message(paste("Total sources:", total_sources))
        log_message(paste("Batch size:", batch_size))
        log_message(paste("Max workers (nb_cores):", nb_cores))
        log_message(paste("Max pivots per molecule:", max_pivots))

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
                        "tol_match_F1", "tol_match_F2", "tol_correspondance_F1", "tol_correspondance_F2",
                        "ratio_min", "ransy_cv_threshold", "correlation_threshold",
                        "coverage_threshold_group", "coverage_threshold_subchain",
                        "coverage_threshold_ransy", "max_missing_ransy", "max_pivots", "alpha_cost",
                        "trouver_correspondances_2D", "calculer_STOCSY_2D", "gfilter_groups_par_coverage_2D",
                        "calculer_RANSY_2D", "calculer_coverage_sub_ratio_2D", "former_groupes_stocsy_2D",
                        "creer_resultat_defaut_2D", "get_mean_intensity_for_ppm2d",
                        "sanitize_for_csv"
                    ), envir = environment())
                    log_message(paste("Parallel cluster initialized with", cores, "workers"))
                    cluster
                },
                error = function(e) {
                    log_message(paste("Cluster error:", e$message, "- switch to single-thread"))
                    NULL
                }
            )
            if (is.null(cl)) use_parallel <- FALSE
        }
        `%op%` <- if (use_parallel) `%dopar%` else `%do%`
        if (!is.null(cl)) {
            on.exit(
                {
                    tryCatch(
                        {
                            parallel::stopCluster(cl)
                            log_message("Cluster closed cleanly")
                        },
                        error = function(e) {}
                    )
                },
                add = TRUE
            )
        }

        total_batches <- ceiling(total_sources / batch_size)
        if (!is.finite(total_batches) || total_batches < 1) total_batches <- 1L
        log_progress(0, 100, "Initializing COSY")

        for (batch_idx in seq(from = start_batch, to = total_batches)) {
            batch_start_time <- Sys.time()
            start_idx <- (batch_idx - 1) * batch_size + 1
            end_idx <- min(batch_idx * batch_size, total_sources)
            batch_sources <- if (total_sources > 0) sources[start_idx:end_idx] else character(0)
            
            # Percentage calculation based on processed sources
            pct_start <- round(((start_idx - 1) / max(1, total_sources)) * 100, 1)

            log_message(paste("=== BATCH", batch_idx, "/", total_batches, "==="))
            log_message(paste("Sources:", length(batch_sources), "(", start_idx, "-", end_idx, ")"))
            log_progress(pct_start, 100, sprintf("Batch %d/%d — sources %d-%d", batch_idx, total_batches, start_idx, end_idx))
            log_memory(paste("batch start", batch_idx))

            batch_results <- tryCatch(
                {
                    foreach(i = seq_along(batch_sources), .packages = "data.table", .combine = "c", .init = list(), .errorhandling = "pass") %op% {
                        src <- batch_sources[i]
                        tryCatch(
                            {
                                dt_mol <- tableau_final[Source == src, .(shift_value_F1, shift_value_F2, ECFP4, H1, H2, Bonds, Group1, Group2)]
                                ecfp_str <- paste(unique(dt_mol$ECFP4), collapse = "|")

                                if (!nrow(dt_mol)) {
                                    return(list(creer_resultat_defaut_2D(src, ecfp_str, dt_mol, character())))
                                }

                                matches <- trouver_correspondances_2D(dt_mol, tol_correspondance_F1, tol_correspondance_F2)
                                match_ppms <- if (nrow(matches)) unique(matches$ppm_matrix) else character()
                                if (!nrow(matches)) {
                                    return(list(creer_resultat_defaut_2D(src, ecfp_str, dt_mol, character())))
                                }

                                all_results <- former_groupes_stocsy_2D(
                                    matched_ppms = match_ppms,
                                    molecule_shifts = dt_mol,
                                    tol_match_F1 = tol_match_F1,
                                    tol_match_F2 = tol_match_F2,
                                    correlation_threshold = correlation_threshold,
                                    coverage_threshold_group = coverage_threshold_group,
                                    coverage_threshold_subchain = coverage_threshold_subchain,
                                    ratio_min = ratio_min,
                                    ransy_cv_threshold = ransy_cv_threshold,
                                    coverage_threshold_ransy = coverage_threshold_ransy,
                                    max_missing = max_missing_ransy,
                                    max_pivots = max_pivots,
                                    alpha = alpha_cost
                                )

                                if (is.null(all_results) || length(all_results) == 0) {
                                    return(list(data.table(
                                        Source = src, ECFP4 = ecfp_str,
                                        ppm_molecule = paste(dt_mol$shift_value_F1, dt_mol$shift_value_F2, sep = "|", collapse = "||"),
                                        ppm_matrix = paste(match_ppms, collapse = "||"),
                                        groupe_correlation = "",
                                        baseline_ppm = NA_character_, sous_chaine = "", coverage_sub = 0, coverage_sub_ratio = 0,
                                        pivot_ppm = NA_character_, RANSY_signals = "", ppm_molecule_matched = "",
                                        RANSY_mean_intensities = "", coverage_ransy = 0, RANSY_variabilite = "",
                                        avg_variabilite = 0, intensity_ratios = "", error = "no_ransy_valid"
                                    )))
                                }

                                rows <- rbindlist(lapply(all_results, function(best) {
                                    coverage_sub_ratio <- calculer_coverage_sub_ratio_2D(best$sub_chain, dt_mol, tol_match_F1, tol_match_F2)
                                    ransy_all <- as.character(best$RANSY_signals_all)
                                    ransy_means <- vapply(ransy_all, get_mean_intensity_for_ppm2d, numeric(1))

                                    pairs <- best$matched_pairs
                                    if (is.null(pairs) || !nrow(pairs)) pairs <- data.table(mol_F1 = numeric(0), mol_F2 = numeric(0), ppm_ransy = character(0))
                                    ppm_molecule_matched <- if (nrow(pairs)) paste(paste(pairs$mol_F1, pairs$mol_F2, sep = "|"), collapse = "||") else ""

                                    data.table(
                                        Source = src, ECFP4 = ecfp_str,
                                        ppm_molecule = paste(dt_mol$shift_value_F1, dt_mol$shift_value_F2, sep = "|", collapse = "||"),
                                        ppm_matrix = paste(match_ppms, collapse = "||"),
                                        groupe_correlation = paste(best$groupe_ppms, collapse = "-"),
                                        baseline_ppm = best$base_ppm,
                                        sous_chaine = paste(best$sub_chain, collapse = "|"),
                                        coverage_sub = best$coverage_sub,
                                        coverage_sub_ratio = coverage_sub_ratio,
                                        pivot_ppm = best$pivot_ppm,
                                        RANSY_signals = paste(ransy_all, collapse = "||"),
                                        ppm_molecule_matched = ppm_molecule_matched,
                                        RANSY_mean_intensities = sanitize_for_csv(paste(round(ransy_means, 6), collapse = ";")),
                                        coverage_ransy = best$coverage_ransy,
                                        RANSY_variabilite = sanitize_for_csv(best$RANSY_variabilite),
                                        avg_variabilite = best$avg_variabilite,
                                        intensity_ratios = sanitize_for_csv(best$intensity_ratios),
                                        error = NA_character_
                                    )
                                }), fill = TRUE)

                                list(rows)
                            },
                            error = function(e) {
                                log_message(sprintf("Error for source %s: %s", src, e$message))
                                dt_mol_safe <- tryCatch(tableau_final[Source == src, .(shift_value_F1, shift_value_F2, ECFP4, H1, H2, Bonds, Group1, Group2)], error = function(e2) data.table())
                                ecfp_safe <- tryCatch(paste(unique(dt_mol_safe$ECFP4), collapse = "|"), error = function(x) "")
                                list(creer_resultat_defaut_2D(src, ecfp_safe, dt_mol_safe, character()))
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

            if (batch_idx %% 2 == 0) {
                gc()
                log_memory(paste("after batch", batch_idx))
            }
            if (enable_checkpoints && batch_idx %% 5 == 0) {
                processed_sources_so_far <- if (total_sources > 0) sources[1:end_idx] else character(0)
                save_checkpoint(results_list, processed_sources_so_far, batch_idx)
            }

            batch_duration <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
            log_message(paste("Batch", batch_idx, "completed in", round(batch_duration, 1), "seconds"))
            
            # Percentage calculation based on processed sources
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
            data.table::setorder(resultats, Source, -coverage_ransy)
            best_per_sourceCOSY <- resultats[, .SD[1], by = Source]
            best_per_sourceCOSY[, Source := sub("^NMR_", "", Source)]
            best_per_sourceCOSY[, Source := sub("\\.sdf$", "", Source)]
        } else {
            log_error("No results to assemble")
            quit(status = 1)
        }

        output_file <- file.path(work_dir, "COSY_With_Ransy_and_ECFP4.csv")
        data.table::fwrite(best_per_sourceCOSY, output_file, na = "NA", sep = ";", quote = TRUE, bom = TRUE)

        if (enable_checkpoints && file.exists(checkpoint_file)) {
            file.remove(checkpoint_file)
            log_message("Final checkpoint removed")
        }

        log_message("=== FINAL COSY RESULTS ===")
        log_message(paste("File saved:", output_file))
        log_message(paste("Sources processed:", nrow(best_per_sourceCOSY)))

        resultats_valides <- best_per_sourceCOSY[is.na(error)]
        resultats_erreur <- best_per_sourceCOSY[!is.na(error)]
        log_message(paste("Success:", nrow(resultats_valides)))
        log_message(paste("Errors:", nrow(resultats_erreur)))
        if (nrow(resultats_erreur) > 0) {
            error_types <- table(resultats_erreur$error)
            for (err_type in names(error_types)) log_message(paste("  -", err_type, ":", error_types[[err_type]]))
        }
        log_memory("final")

        log_result(list(
            total_sources = length(sources),
            successful_sources = nrow(resultats_valides),
            error_sources = nrow(resultats_erreur),
            output_file = output_file,
            annotation_type = "COSY",
            version = "stocsy_pipeline_optimise_dt",
            batches_processed = ceiling(length(sources) / batch_size),
            memory_optimized = TRUE,
            summary = list(
                avg_coverage_ransy = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$coverage_ransy, na.rm = TRUE), 2) else 0,
                avg_coverage_sub   = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$coverage_sub, na.rm = TRUE), 2) else 0,
                avg_variabilite    = if (nrow(resultats_valides) > 0) round(mean(resultats_valides$avg_variabilite, na.rm = TRUE), 4) else 0
            )
        ))

        log_progress(100, 100, "COSY completed")
    },
    error = function(e) {
        log_error(paste("Critical error in pipeline:", e$message))
        quit(status = 1)
    }
)

log_message("=== COSY ANNOTATION WITH STOCSY (OK) ===")
