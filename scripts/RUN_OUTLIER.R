# Fixed SVM sensitivity analysis: all training spectra versus within-farm MBD quality control.
sourced_files <- unlist(lapply(sys.frames(), function(frame) {
  if (exists("ofile", envir = frame, inherits = FALSE))
    get("ofile", envir = frame, inherits = FALSE)
  else character(0)
}), use.names = FALSE)
command_files <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
active_file <- character(0)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  active_file <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                          error = function(e) character(0))
}
script_files <- c(rev(sourced_files), command_files, active_file)
script_files <- script_files[!is.na(script_files) & nzchar(script_files)]
candidates <- unique(c(Sys.getenv("FTIR_PROJECT_ROOT", unset = ""),
                       dirname(script_files), file.path(dirname(script_files), ".."), getwd(), file.path(getwd(), "FTIR_unified_code")))
candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
valid <- vapply(candidates, function(dir) {
  file.exists(file.path(dir, "R", "FTIR_FUNCTIONS.R")) &&
    all(file.exists(file.path(dir, "data_synthetic", c("FTIR_simulated_spectra_dataset_1.xlsx", "FTIR_simulated_spectra_dataset_2.xlsx"))))
}, logical(1))
if (!any(valid)) {
  stop("FTIR project not found. Open the repository root or set FTIR_PROJECT_ROOT.")
}
BASE_DIR <- normalizePath(candidates[which(valid)[1L]], winslash = "/")
DATA_DIR <- file.path(BASE_DIR, "data_synthetic")
setwd(BASE_DIR)
source(file.path(BASE_DIR, "R", "FTIR_FUNCTIONS.R"), local = FALSE)

required <- c("readxl", "dplyr", "stringr", "signal", "pracma",
              "roahd", "caret", "e1071")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing R packages: ", paste(missing, collapse = ", "))
}

SEED <- 69L
OUTER_FOLDS <- 10L
INNER_FOLDS <- 5L
PVE <- 0.9999
SG_DEGREE <- 4L
FENCE_MULTIPLIER <- 1
ALLOWED_EXCEED_FRACTION <- 0

inputs <- list(
  Sarda = readxl::read_excel(file.path(DATA_DIR, "FTIR_simulated_spectra_dataset_1.xlsx")),
  `Valle del Belice` = readxl::read_excel(file.path(DATA_DIR, "FTIR_simulated_spectra_dataset_2.xlsx"))
)
normalized <- normalizza_dataset(inputs$Sarda, inputs$`Valle del Belice`)
inputs <- list(Sarda = normalized$d1, `Valle del Belice` = normalized$d2)

# Each row is an independent sample; retain the original animal ID as metadata.
for (breed in names(inputs)) {
  dat <- inputs[[breed]]
  dat$Original_ID <- dat$Matricola
  dat$Matricola <- paste0(if (breed == "Sarda") "SW" else "VDB",
                         "_sample_", seq_len(nrow(dat)))
  dat$Altezza <- factor(as.character(dat$Altezza), levels = ALTITUDE_LEVELS)
  if (anyNA(dat$Altezza) || anyNA(dat$Azienda)) {
    stop("Missing altitude class or farm for ", breed)
  }
  if (length(get_spectral_columns(dat, "240", "1299")) != 1060L) {
    stop("Expected 1,060 spectral points for ", breed)
  }
  spec <- as.matrix(dat[, get_spectral_columns(dat, "240", "1299"), drop = FALSE])
  if (any(!is.finite(spec))) stop("Non-finite spectral values in ", breed)
  inputs[[breed]] <- dat
}

message("Fixed chemical-domain linear SVM (C = 0.1, 1, 10); smoothed spectra + first derivatives; PVE = 99.99%")
message("Ten matched outer folds, five inner folds, training-only preprocessing and MBD QC.")
results <- lapply(names(inputs), function(breed) {
  run_outlier_sensitivity_nested(
    data = inputs[[breed]], breed = breed,
    seed = SEED, outer_folds = OUTER_FOLDS, inner_folds = INNER_FOLDS,
    pve = PVE, polynomial_degree = SG_DEGREE,
    fence_multiplier = FENCE_MULTIPLIER,
    allowed_exceed_fraction = ALLOWED_EXCEED_FRACTION
  )
})
names(results) <- names(inputs)

outlier_comparison <- do.call(rbind, lapply(results, `[[`, "summary"))
outlier_fold_metrics <- do.call(rbind, lapply(results, `[[`, "fold_metrics"))
outlier_predictions <- do.call(rbind, lapply(results, `[[`, "predictions"))
outlier_farm_qc <- do.call(rbind, lapply(results, `[[`, "farm_qc"))

print(outlier_comparison, row.names = FALSE)
output_dir <- file.path(BASE_DIR, "OUTPUT_OUTLIER")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(outlier_comparison, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(outlier_fold_metrics, file.path(output_dir, "outer_fold_metrics.csv"), row.names = FALSE)
utils::write.csv(outlier_predictions, file.path(output_dir, "outer_test_predictions.csv"), row.names = FALSE)
utils::write.csv(outlier_farm_qc, file.path(output_dir, "farm_qc_counts.csv"), row.names = FALSE)
message("Outlier sensitivity outputs saved in: ", output_dir)

