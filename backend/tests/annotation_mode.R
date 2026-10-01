# Run from the repository root: Rscript backend/tests/annotation_mode.R
# Evaluate the real argument/parameter parsing, without launching annotation workers.
expressions <- parse("backend/scripts/annotation_COSY.R")
check_mode <- function(parameters, expected) {
    scope <- new.env(parent = globalenv())
    scope$commandArgs <- function(trailingOnly = FALSE) {
        c(tempdir(), "prediction.csv", "matrix.csv", parameters)
    }
    for (expr in expressions) {
        if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
            identical(expr[[2]], as.name("nb_cores"))) break
        eval(expr, envir = scope)
    }
    stopifnot(identical(scope$annotation_type, expected))
    stopifnot(identical(scope$output_basename, paste0(expected, "_With_Ransy_and_ECFP4.csv")))
}
check_mode('{}', 'COSY')
check_mode('{"annotation_type":"COSY"}', 'COSY')
check_mode('{"annotation_type":"TOCSY"}', 'TOCSY')
check_mode('{"annotation_type":"tocsy"}', 'TOCSY')
cat("OK: COSY default and explicit COSY/TOCSY output filenames\n")
