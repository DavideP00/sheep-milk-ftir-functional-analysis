# Altitude-associated production-system classification (SVM and ENET).
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
if (!any(valid)) stop("FTIR project not found. Open the repository root in RStudio or set FTIR_PROJECT_ROOT.")
BASE_DIR <- normalizePath(candidates[which(valid)[1L]], winslash = "/")
DATA_DIR <- file.path(BASE_DIR, "data_synthetic")
setwd(BASE_DIR)
source(file.path(BASE_DIR, "R", "FTIR_FUNCTIONS.R"), local = FALSE)
library(pracma, exclude = "pdist2")
library(MASS, exclude = "select")

library(readxl)
library(dplyr)

# Primary configuration: sample-level nested CV only.
MODELS <- normalise_model_names(c("svm", "enet"))

RUN_SARDA <- TRUE
RUN_VDB <- TRUE
CONF_LEVEL <- 0.95

RUN_PAPER_SAMPLE_CV <- TRUE


PARALLEL_ENABLED <- TRUE
N_WORKERS <- 7L
MAX_WORKERS <- 7L
FUTURE_GLOBALS_MAX_GB <- 8

PARALLEL_CFG <- q1_configure_parallel(
  enabled = PARALLEL_ENABLED,
  workers = N_WORKERS,
  max_workers = MAX_WORKERS,
  globals_max_gb = FUTURE_GLOBALS_MAX_GB
)

q1_clear_sflda_cache()

# All spectra are retained; preprocessing and tuning are training-only.
FIXED_SG_WINDOW <- NULL
POLYNOMIAL_DEGREE <- 4

DOMAIN_GRID <- c("full", "chemical", "discriminative")

DERIVATIVE_GRID <- c("m0", "m0_d1", "m0_d1_d2")

PVE_GRID <- c(0.95, 0.99, 0.9999)

DISCR_THRESHOLD_GRID <- c(
  0.001,
  0.10, 0.15, 0.20,
  0.30, 0.40
)

SFLDA_BASE_INTERVALS <- ftir_domain_intervals("chemical")
SFLDA_PARAM_GRID <- q1_make_sflda_param_grid()
SFLDA_CV_FOLDS <- 3L

MODEL_GRIDS_USED <- setNames(lapply(MODELS, default_model_grid), MODELS)

PAPER_SEED <- 69L
PAPER_OUTER_K_FOLDS <- 10L
PAPER_INNER_K_FOLDS <- 5L


raw_sw <- read_excel(file.path(DATA_DIR, "FTIR_simulated_spectra_dataset_1.xlsx"))
raw_vdb <- read_excel(file.path(DATA_DIR, "FTIR_simulated_spectra_dataset_2.xlsx"))

norm <- normalizza_dataset(raw_sw, raw_vdb)
dati_sw <- norm$d1
dati_vdb <- norm$d2

dati_sw$Altezza <- factor(
  dati_sw$Altezza,
  levels = ALTITUDE_LEVELS,
  ordered = TRUE
)

dati_vdb$Altezza <- factor(
  dati_vdb$Altezza,
  levels = ALTITUDE_LEVELS,
  ordered = TRUE
)

cat("\n============================================================\n")
cat("ALTITUDE CLASSIFICATION | SVM + ENET\n")
cat("============================================================\n")
cat("Models: ", paste(MODELS, collapse = ", "), "\n", sep = "")
cat("RUN_PAPER_SAMPLE_CV: ", RUN_PAPER_SAMPLE_CV, "\n", sep = "")
cat("Run Sarda: ", RUN_SARDA, "\n", sep = "")
cat("Run Valle del Belice: ", RUN_VDB, "\n", sep = "")
cat("Domains: ", paste(DOMAIN_GRID, collapse = ", "), "\n", sep = "")
cat("Derivative configs: ", paste(DERIVATIVE_GRID, collapse = ", "), "\n", sep = "")
cat("PVE grid: ", paste(PVE_GRID, collapse = ", "), "\n", sep = "")
cat("Discriminative threshold grid: ",
    paste(DISCR_THRESHOLD_GRID, collapse = ", "), "\n", sep = "")
cat("SFMLDA grid: ", nrow(SFLDA_PARAM_GRID),
    " tau/lambda combinations per chemical interval\n", sep = "")
cat("Feature configurations: ",
    length(DERIVATIVE_GRID) * length(PVE_GRID) * 2L +
      length(DERIVATIVE_GRID) * length(PVE_GRID) *
      length(DISCR_THRESHOLD_GRID),
    " (expected 72)\n", sep = "")
cat("SVM model configs: ", nrow(MODEL_GRIDS_USED$svm),
    " (expected 7)\n", sep = "")
cat("ENET model configs: ", nrow(MODEL_GRIDS_USED$enet),
    " (expected 9)\n", sep = "")
cat("Outliers removed in primary analyses: 0\n")
cat("Scaling: training-only min-max\n")
cat("Windows parallel workers: ", PARALLEL_CFG$workers, "\n", sep = "")
cat("Sarda: ", nrow(dati_sw), " samples from ",
    length(unique(dati_sw$Azienda)), " farms\n", sep = "")
cat("Valle del Belice: ", nrow(dati_vdb), " samples from ",
    length(unique(dati_vdb$Azienda)), " farms\n", sep = "")

if (RUN_PAPER_SAMPLE_CV) {
  cat("Sample-level design: outer ", PAPER_OUTER_K_FOLDS,
      "-fold + inner ", PAPER_INNER_K_FOLDS,
      "-fold stratified sample CV\n", sep = "")
}



farm_description <- print_farm_description(dati_sw, dati_vdb)

paper_sw <- list()
paper_vdb <- list()

if (RUN_PAPER_SAMPLE_CV) {

  cat("\n\n############################################################\n")
  cat("STARTING UPDATED SAMPLE-LEVEL NESTED CV\n")
  cat("############################################################\n")

# Reuse model-independent SFMLDA fits across SVM and ENET.
  q1_clear_sflda_cache()

  for (MODEL in MODELS) {

    MODEL_GRID <- MODEL_GRIDS_USED[[MODEL]]

    cat("\n\n============================================================\n")
    cat("SAMPLE-LEVEL MODEL: ", toupper(MODEL), "\n", sep = "")
    cat("============================================================\n")

    if (RUN_SARDA) {
      paper_sw[[MODEL]] <- run_paper_sample_nested_model_discriminative(
        data = dati_sw,
        model_type = MODEL,
        model_grid = MODEL_GRID,
        domain_grid = DOMAIN_GRID,
        derivative_configs = DERIVATIVE_GRID,
        pve_grid = PVE_GRID,
        threshold_grid = DISCR_THRESHOLD_GRID,
        sflda_intervals = SFLDA_BASE_INTERVALS,
        sflda_param_grid = SFLDA_PARAM_GRID,
        sflda_cv_folds = SFLDA_CV_FOLDS,
        outer_k_folds = PAPER_OUTER_K_FOLDS,
        inner_k_folds = PAPER_INNER_K_FOLDS,
        polynomial_degree = POLYNOMIAL_DEGREE,
        fixed_sg_window = FIXED_SG_WINDOW,
        conf_level = CONF_LEVEL,
        seed = PAPER_SEED,
        parallel_cfg = PARALLEL_CFG,
        verbose = TRUE
      )

      print_paper_model_discriminative_results(
        paper_sw[[MODEL]],
        "Sarda",
        MODEL
      )
    }

    if (RUN_VDB) {
      paper_vdb[[MODEL]] <- run_paper_sample_nested_model_discriminative(
        data = dati_vdb,
        model_type = MODEL,
        model_grid = MODEL_GRID,
        domain_grid = DOMAIN_GRID,
        derivative_configs = DERIVATIVE_GRID,
        pve_grid = PVE_GRID,
        threshold_grid = DISCR_THRESHOLD_GRID,
        sflda_intervals = SFLDA_BASE_INTERVALS,
        sflda_param_grid = SFLDA_PARAM_GRID,
        sflda_cv_folds = SFLDA_CV_FOLDS,
        outer_k_folds = PAPER_OUTER_K_FOLDS,
        inner_k_folds = PAPER_INNER_K_FOLDS,
        polynomial_degree = POLYNOMIAL_DEGREE,
        fixed_sg_window = FIXED_SG_WINDOW,
        conf_level = CONF_LEVEL,
        seed = PAPER_SEED,
        parallel_cfg = PARALLEL_CFG,
        verbose = TRUE
      )

      print_paper_model_discriminative_results(
        paper_vdb[[MODEL]],
        "Valle del Belice",
        MODEL
      )
    }
  }

  if (length(MODELS) >= 2L) {

    if (RUN_SARDA) {
      paper_comparison_sw <- compare_joint_paper_models(
        paper_sw,
        conf_level = CONF_LEVEL,
        seed = PAPER_SEED
      )

      print_paired_model_comparison(
        paper_comparison_sw,
        paste0(
          "PAIRED MODEL COMPARISON | UPDATED PAPER SAMPLE-CV | ",
          "JOINT DOMAIN SELECTION | Sarda"
        )
      )
    }

    if (RUN_VDB) {
      paper_comparison_vdb <- compare_joint_paper_models(
        paper_vdb,
        conf_level = CONF_LEVEL,
        seed = PAPER_SEED
      )

      print_paired_model_comparison(
        paper_comparison_vdb,
        paste0(
          "PAIRED MODEL COMPARISON | UPDATED PAPER SAMPLE-CV | ",
          "JOINT DOMAIN SELECTION | Valle del Belice"
        )
      )
    }
  }

  print_analysis_reproducibility(
    settings = list(
      MODELS = MODELS,
      ANALYSIS = "updated paper sample-level nested CV",
      PRIMARY_DESIGN = paste0(
        "outer ", PAPER_OUTER_K_FOLDS,
        "-fold stratified sample CV + inner ", PAPER_INNER_K_FOLDS,
        "-fold stratified sample CV"
      ),
      RUN_PAPER_SAMPLE_CV = RUN_PAPER_SAMPLE_CV,
      RUN_SARDA = RUN_SARDA,
      RUN_VDB = RUN_VDB,
      CONF_LEVEL = CONF_LEVEL,
      PAPER_SEED = PAPER_SEED,
      PAPER_OUTER_K_FOLDS = PAPER_OUTER_K_FOLDS,
      PAPER_INNER_K_FOLDS = PAPER_INNER_K_FOLDS,
      DOMAIN_GRID = DOMAIN_GRID,
      DERIVATIVE_GRID = DERIVATIVE_GRID,
      PVE_GRID = PVE_GRID,
      DISCR_THRESHOLD_GRID = DISCR_THRESHOLD_GRID,
      SFLDA_CV_FOLDS = SFLDA_CV_FOLDS,
      SFLDA_PARAM_GRID = SFLDA_PARAM_GRID,
      POLYNOMIAL_DEGREE = POLYNOMIAL_DEGREE,
      FIXED_SG_WINDOW = FIXED_SG_WINDOW,
      OUTLIER_REMOVAL = "none",
      SCALING = "training-only min-max",
      INNER_SELECTION_METRIC = "mean inner sample accuracy",
      PARALLEL_BACKEND = if (PARALLEL_CFG$enabled) {
        "future::multisession (Windows PSOCK)"
      } else {
        "sequential"
      },
      PARALLEL_WORKERS = PARALLEL_CFG$workers,
      SFMLDA_COMPROMISE = paste0(
        "tau/lambda tuned once on each complete outer sample-training partition; ",
        "fixed during inner sample CV; beta profiles/intervals refitted on each ",
        "inner training partition"
      ),
      CROSS_MODEL_SFLDA_CACHE = "shared between SVM and ENET within sample-level design"
    ),
    model_grids = MODEL_GRIDS_USED
  )
}

q1_shutdown_parallel()

cat("\n\n############################################################\n")
cat("ALTITUDE CLASSIFICATION COMPLETED\n")
cat("Models: ", paste(MODELS, collapse = ", "), "\n", sep = "")
cat("Sample-level run: ", RUN_PAPER_SAMPLE_CV, "\n", sep = "")
cat("No outliers were removed.\n")
cat("No result files were written automatically.\n")
cat("Objects retained in memory:\n")

if (RUN_PAPER_SAMPLE_CV) {
  cat("  paper_sw, paper_vdb\n")
  if (length(MODELS) >= 2L) {
    cat("  paper_comparison_sw, paper_comparison_vdb\n")
  }
}



cat("############################################################\n")
