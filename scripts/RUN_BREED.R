# Supporting exploratory breed classification: chemical intervals, SVM and nested sample-level CV.

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
FILE_SARDA <- "FTIR_simulated_spectra_dataset_1.xlsx"
FILE_VALLE_DEL_BELICE <- "FTIR_simulated_spectra_dataset_2.xlsx"
OUTPUT_DIR <- file.path(getwd(), "OUTPUT_BREED_NCV")

# Reported analysis: fixed chemical domain, SVM and 99.99% PVE.
ANALYSIS_MODES <- c("chemical")
MODELS <- c("svm")
DERIVATIVE_GRID <- c("m0", "m0_d1")
PVE_GRID <- c(0.9999)
OUTER_FOLDS <- 10L
INNER_FOLDS <- 5L
SEED <- 69L
POLYNOMIAL_DEGREE <- 4L
FIXED_SG_WINDOW <- NULL

PARALLEL_INNER <- TRUE
N_WORKERS <- 5L
VERBOSE <- TRUE

BREED_LEVELS <- c("Sarda", "Valle del Belice")
ALTITUDE_LEVELS <- BREED_LEVELS
SUPPORTED_MODELS <- c("svm", "enet")

required <- c("readxl", "dplyr", "stringr", "signal", "pracma", "caret",
              "e1071", "glmnet")
if (PARALLEL_INNER) required <- c(required, "future", "future.apply")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install the missing packages: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
stopifnot(all(ANALYSIS_MODES %in% c("chemical", "full", "joint")),
          all(MODELS %in% c("svm", "enet")),
          all(DERIVATIVE_GRID %in% c("m0", "m0_d1", "m0_d1_d2")),
          all(PVE_GRID > 0 & PVE_GRID < 1))
if (!length(ANALYSIS_MODES) || !length(MODELS)) stop("Select at least one analysis and one model.")

read_breed_data <- function(file_path, breed_label, prefix) {
  if (!file.exists(file_path)) stop("Data file not found: ", file_path)
  df <- as.data.frame(readxl::read_excel(file_path, .name_repair = "minimal"),
                      check.names = FALSE)
  if (ncol(df) < 4L) stop("The file must contain three metadata columns and the FTIR spectral columns: ", file_path)
  idx <- suppressWarnings(as.integer(sub("^VAR_", "", names(df)[-(1:3)])))
  if (anyNA(idx) || anyDuplicated(idx)) {
    stop("Unrecognized or duplicated spectral columns in: ", file_path,
         ". Expected VAR_240...VAR_1299 (or 240...1299).")
  }
  needed <- 240:1299
  if (!setequal(idx, needed)) {
    stop("Expected FTIR indices 240:1299; missing: ",
         paste(setdiff(needed, idx), collapse = ", "),
         "; extra: ", paste(setdiff(idx, needed), collapse = ", "))
  }
  raw_spec <- df[, -(1:3), drop = FALSE]
  raw_spec <- raw_spec[, match(needed, idx), drop = FALSE]
  X <- vapply(raw_spec, function(z) suppressWarnings(as.numeric(as.character(z))),
              numeric(nrow(df)))
  if (!is.matrix(X)) X <- matrix(X, nrow = nrow(df))
  if (any(!is.finite(X))) stop("Missing, infinite or nonnumeric spectral values in: ", file_path)
  colnames(X) <- as.character(needed)
  farm <- as.character(df[[1]])
  if (anyNA(farm) || any(!nzchar(farm))) stop("Missing farm identifiers for some samples in: ", file_path)
  meta <- data.frame(
    Azienda = paste0(prefix, "_", farm),
    Altezza = as.character(df[[2]]),
    Matricola = paste0(prefix, "_", seq_len(nrow(df))),
    Original_ID = as.character(df[[3]]),
    Razza = factor(rep(breed_label, nrow(df)), levels = BREED_LEVELS),
    stringsAsFactors = FALSE
  )
  spec <- as.data.frame(X, check.names = FALSE)
  names(spec) <- as.character(needed)
  out <- cbind(meta, spec)
  names(out)[(ncol(meta)+1L):ncol(out)] <- as.character(needed)
  out
}

all_data <- rbind(
  read_breed_data(file.path(DATA_DIR, FILE_SARDA), "Sarda", "SW"),
  read_breed_data(file.path(DATA_DIR, FILE_VALLE_DEL_BELICE), "Valle del Belice", "VDB")
)
all_data$Razza <- factor(as.character(all_data$Razza), levels = BREED_LEVELS)
if (anyDuplicated(all_data$Matricola)) stop("Identificativi duplicati.")
if (any(table(all_data$Razza) < OUTER_FOLDS)) stop("Insufficient samples in at least one class.")
if (nrow(all_data) != 910L || any(table(all_data$Razza) != c(460L, 450L))) {
  warning("The input does not contain exactly 460 Sarda and 450 Valle del Belice samples; check the data.")
}
cat("\nSAMPLES BY BREED:\n"); print(table(all_data$Razza))
cat("FARMS BY BREED:\n")
print(tapply(all_data$Azienda, all_data$Razza, function(z) length(unique(z))))

macro_f1_binary <- function(y, p) f1_macro(y, p, levels = BREED_LEVELS)
class_statistics <- function(y, p) {
  cm <- table(True = factor(as.character(y), levels = BREED_LEVELS),
              Predicted = factor(as.character(p), levels = BREED_LEVELS))
  do.call(rbind, lapply(seq_along(BREED_LEVELS), function(j) {
    tp <- cm[j,j]; fp <- sum(cm[,j]) - tp; fn <- sum(cm[j,]) - tp
    precision <- if (tp+fp) tp/(tp+fp) else 0
    recall <- if (tp+fn) tp/(tp+fn) else 0
    f1 <- if (precision+recall) 2*precision*recall/(precision+recall) else 0
    data.frame(class = BREED_LEVELS[j], n = sum(cm[j,]), precision = precision,
               recall = recall, f1 = f1)
  }))
}

fit_binary_model <- function(x, y, model, param, seed) {
  set.seed(seed)
  xx <- as.matrix(x)
  yy <- factor(as.character(y), levels = BREED_LEVELS)
  if (model == "svm") {
    if (param$kernel == "linear") {
      fit <- e1071::svm(xx, yy, kernel = "linear", cost = param$cost,
                        scale = FALSE, probability = FALSE)
    } else {
      gamma <- param$gamma_multiplier / max(1L, ncol(xx))
      fit <- e1071::svm(xx, yy, kernel = "radial", cost = param$cost,
                        gamma = gamma, scale = FALSE, probability = FALSE)
    }
  } else if (model == "enet") {
# Use a binary response for the optional ENET classifier.
    fit <- glmnet::glmnet(xx, yy, family = "binomial",
                          alpha = param$alpha, lambda = param$lambda,
                          standardize = FALSE)
  } else stop("Unsupported model: ", model)
  list(fit = fit, model = model, param = param)
}

predict_binary_model <- function(object, x) {
  xx <- as.matrix(x)
  if (object$model == "svm") {
    z <- predict(object$fit, xx)
  } else {
    z <- predict(object$fit, newx = xx, s = object$param$lambda, type = "class")
  }
  factor(as.character(drop(z)), levels = BREED_LEVELS)
}

DOMAIN_GRID <- unique(c(
  if ("full" %in% ANALYSIS_MODES || "joint" %in% ANALYSIS_MODES) "full",
  if ("chemical" %in% ANALYSIS_MODES || "joint" %in% ANALYSIS_MODES) "chemical"
))
VIEW_GRID <- unique(unlist(lapply(DERIVATIVE_GRID, derivative_views_from_config)))
MODEL_GRIDS <- setNames(lapply(MODELS, default_model_grid), MODELS)

build_inner_cache <- function(outer_train, inner_train_indices) {
  worker <- function(j) {
    tr <- inner_train_indices[[j]]
    te <- setdiff(seq_len(nrow(outer_train)), tr)
    pp <- preprocess_grouped_split(
      outer_train[tr, , drop = FALSE], outer_train[te, , drop = FALSE],
      col_start = "240", col_end = "1299",
      polynomial_degree = POLYNOMIAL_DEGREE,
      fixed_window = FIXED_SG_WINDOW
    )
    build_preprocessed_fpca_cache(
      train_views = pp$train, test_views = pp$test,
      domains = DOMAIN_GRID, views_needed = VIEW_GRID,
      max_pve = max(PVE_GRID),
      id_col = "Matricola", group_col = "Azienda", target_col = "Razza"
    )
  }
  if (PARALLEL_INNER) {
    future.apply::future_lapply(seq_along(inner_train_indices), worker,
                                future.seed = TRUE)
  } else {
    lapply(seq_along(inner_train_indices), worker)
  }
}

# Fit all transformations using each inner-training partition only.
tune_inner <- function(inner_cache, fold_seed) {
  rows <- list()
  z <- 0L
  for (dom in DOMAIN_GRID) for (der in DERIVATIVE_GRID) for (pve in PVE_GRID) {
    features <- lapply(inner_cache, assemble_fpca_features,
                       domain = dom, derivative_config = der, pve = pve,
                       target_col = "Razza", scale_minmax = TRUE)
    if (any(vapply(features, is.null, logical(1)))) next
    nf <- mean(vapply(features, function(q) q$n_features, numeric(1)))
    for (model in MODELS) {
      grid <- MODEL_GRIDS[[model]]
      for (k in seq_len(nrow(grid))) {
        par <- grid[k, , drop = FALSE]
        acc <- numeric(length(features)); f1 <- numeric(length(features))
        ok <- TRUE
        for (j in seq_along(features)) {
          ft <- features[[j]]
          pr <- tryCatch({
            fitted <- fit_binary_model(ft$x_train, ft$y_train, model, par,
                                       seed = fold_seed + j + 100L*k)
            predict_binary_model(fitted, ft$x_test)
          }, error = function(e) {
            if (VERBOSE) message("  Invalid configuration: ", conditionMessage(e))
            NULL
          })
          if (is.null(pr) || anyNA(pr)) {ok <- FALSE; break}
          acc[j] <- mean(as.character(pr) == as.character(ft$y_test))
          f1[j] <- macro_f1_binary(ft$y_test, pr)
        }
        if (!ok) next
        z <- z + 1L
        rows[[z]] <- data.frame(domain = dom, derivative = der, pve = pve,
            model = model, param_id = par$param_id,
            parameters = model_param_label(model, par),
            inner_accuracy = mean(acc), inner_accuracy_sd = stats::sd(acc),
            inner_macro_f1 = mean(f1), n_features = nf)
      }
    }
  }
  if (!length(rows)) stop("All inner-CV configurations failed.")
  do.call(rbind, rows)
}

choose_best <- function(tuning, model, mode) {
  x <- tuning[tuning$model == model, , drop = FALSE]
  if (mode != "joint") x <- x[x$domain == mode, , drop = FALSE]
  if (!nrow(x)) stop("No valid candidate for ", model, " / ", mode)
  ix <- order(-x$inner_accuracy, x$inner_accuracy_sd, x$n_features,
              match(x$domain, c("chemical", "full")), x$pve, x$param_id)
  x[ix[1L], , drop = FALSE]
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (PARALLEL_INNER) {
  future::plan(future::multisession, workers = min(N_WORKERS, INNER_FOLDS))
  options(future.globals.maxSize = 8 * 1024^3)
}
outer_train_folds <- make_stratified_sample_folds(all_data$Razza,
                                                  k = OUTER_FOLDS, seed = SEED)
fold_results <- list(); prediction_results <- list(); tuning_results <- list()
for (fold in seq_along(outer_train_folds)) {
  tr <- outer_train_folds[[fold]]
  te <- setdiff(seq_len(nrow(all_data)), tr)
  train <- all_data[tr, , drop = FALSE]
  test <- all_data[te, , drop = FALSE]
  fold_seed <- SEED + fold * 10000L
  message("\n=========== OUTER FOLD ", fold, "/", OUTER_FOLDS,
          " (training=", nrow(train), ", test=", nrow(test), ") ===========")

  inner_train <- make_stratified_sample_folds(train$Razza,
                                             k = INNER_FOLDS, seed = fold_seed + 1L)
  inner_cache <- build_inner_cache(train, inner_train)
  tuning <- tune_inner(inner_cache, fold_seed)
  tuning$outer_fold <- fold
  tuning_results[[fold]] <- tuning
  rm(inner_cache); invisible(gc(FALSE))

# Refit the selected configuration without using the outer test set.
  outer_pp <- preprocess_grouped_split(
    train, test, col_start = "240", col_end = "1299",
    polynomial_degree = POLYNOMIAL_DEGREE, fixed_window = FIXED_SG_WINDOW
  )
  outer_cache <- build_preprocessed_fpca_cache(
    outer_pp$train, outer_pp$test, domains = DOMAIN_GRID,
    views_needed = VIEW_GRID, max_pve = max(PVE_GRID),
    id_col = "Matricola", group_col = "Azienda", target_col = "Razza"
  )
  for (model in MODELS) for (mode in ANALYSIS_MODES) {
    best <- choose_best(tuning, model, mode)
    ft <- assemble_fpca_features(outer_cache, domain = best$domain,
                derivative_config = best$derivative, pve = best$pve,
                target_col = "Razza", scale_minmax = TRUE)
    par <- MODEL_GRIDS[[model]]
    par <- par[par$param_id == best$param_id, , drop = FALSE]
    fitted <- fit_binary_model(ft$x_train, ft$y_train, model, par,
                               seed = fold_seed + 999L)
    pred <- predict_binary_model(fitted, ft$x_test)
    truth <- factor(as.character(ft$y_test), levels = BREED_LEVELS)
    if (anyNA(pred)) stop("Invalid predictions in outer fold ", fold)
    key <- paste(mode, model, sep = "__")
    fold_results[[length(fold_results) + 1L]] <- data.frame(
      mode = mode, model = model, fold = fold, n_train = nrow(train),
      n_test = length(truth), sg_window = outer_pp$window_length,
      selected_domain = best$domain, selected_derivatives = best$derivative,
      selected_pve = best$pve, selected_params = best$parameters,
      selected_features = ft$n_features, inner_accuracy = best$inner_accuracy,
      accuracy = mean(as.character(truth) == as.character(pred)),
      macro_f1 = macro_f1_binary(truth, pred))
    prediction_results[[length(prediction_results) + 1L]] <- data.frame(
      mode = mode, model = model, fold = fold,
      sample_id = as.character(ft$test_meta$Matricola),
      farm = as.character(ft$test_meta$Azienda),
      truth = as.character(truth), prediction = as.character(pred))
    if (VERBOSE) message("  ", key, " | chosen=", best$domain,
          " | ", best$derivative, " | PVE=", best$pve,
          " | FPCA scores=", ft$n_features,
          " | acc=", round(tail(fold_results, 1L)[[1L]]$accuracy, 4),
          " | macroF1=", round(tail(fold_results, 1L)[[1L]]$macro_f1, 4))
  }
# Fit scaling and models on training data only; preserve completed outer-fold results.
  saveRDS(list(fold = fold, folds = fold_results,
               predictions = prediction_results, tuning = tuning_results),
          file.path(OUTPUT_DIR, sprintf("checkpoint_fold_%02d.rds", fold)))
  rm(outer_cache, outer_pp, tuning); invisible(gc(FALSE))
}
if (PARALLEL_INNER) future::plan(future::sequential)

folds <- do.call(rbind, fold_results)
predictions <- do.call(rbind, prediction_results)
tuning_table <- do.call(rbind, tuning_results)
summary_rows <- list(); class_rows <- list(); confusion_rows <- list()
class_fold_rows <- list()
for (mode in ANALYSIS_MODES) for (model in MODELS) {
  f <- folds[folds$mode == mode & folds$model == model, , drop = FALSE]
  p <- predictions[predictions$mode == mode & predictions$model == model, , drop = FALSE]
  if (nrow(p) != nrow(all_data) || anyDuplicated(p$sample_id)) {
    stop("Audit fallito: ciascun campione deve essere predetto esattamente una volta: ",
         mode, " / ", model)
  }
  ai <- mean_t_ci(f$accuracy)
  fi <- mean_t_ci(f$macro_f1)
  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    mode = mode, model = model, n = nrow(p),
    mean_accuracy = ai["mean"], accuracy_95lo = ai["lower"],
    accuracy_95hi = ai["upper"], sd_accuracy = sd(f$accuracy),
    mean_macro_f1 = fi["mean"], macro_f1_95lo = fi["lower"],
    macro_f1_95hi = fi["upper"], sd_macro_f1 = sd(f$macro_f1),
    pooled_accuracy = mean(p$truth == p$prediction),
    chosen_full_folds = sum(f$selected_domain == "full"),
    chosen_chemical_folds = sum(f$selected_domain == "chemical"))
  cl <- class_statistics(p$truth, p$prediction)
  cl$mode <- mode; cl$model <- model
  for (ff in sort(unique(p$fold))) {
    pf <- p[p$fold == ff, , drop = FALSE]
    cf <- class_statistics(pf$truth, pf$prediction)
    cf$mode <- mode; cf$model <- model; cf$fold <- ff
    class_fold_rows[[length(class_fold_rows)+1L]] <- cf
  }
  class_rows[[length(class_rows)+1L]] <- cl
  cm <- as.data.frame(table(True = factor(p$truth, levels = BREED_LEVELS),
                            Predicted = factor(p$prediction, levels = BREED_LEVELS)))
  names(cm)[3L] <- "count"
  cm$mode <- mode; cm$model <- model
  confusion_rows[[length(confusion_rows)+1L]] <- cm
}
summary_table <- do.call(rbind, summary_rows)
class_table <- do.call(rbind, class_rows)
confusion_table <- do.call(rbind, confusion_rows)
class_fold_table <- do.call(rbind, class_fold_rows)
for (i in seq_len(nrow(class_table))) {
  r <- class_table[i, , drop = FALSE]
  z <- class_fold_table[class_fold_table$mode == r$mode &
                        class_fold_table$model == r$model &
                        class_fold_table$class == r$class, , drop = FALSE]
  for (metric in c("precision", "recall", "f1")) {
    ci <- mean_t_ci(z[[metric]])
    class_table[i, paste0(metric, "_fold95lo")] <- ci["lower"]
    class_table[i, paste0(metric, "_fold95hi")] <- ci["upper"]
  }
}

paired_rows <- list()
for (model in MODELS) if (all(c("full", "chemical") %in% ANALYSIS_MODES)) {
  a <- folds[folds$model == model & folds$mode == "full", , drop = FALSE]
  b <- folds[folds$model == model & folds$mode == "chemical", , drop = FALSE]
  a <- a[order(a$fold), , drop = FALSE]
  b <- b[order(b$fold), , drop = FALSE]
  if (!identical(a$fold, b$fold)) stop("Outer folds are not identical.")
  for (metric in c("accuracy", "macro_f1")) {
    d <- a[[metric]] - b[[metric]]
    ci <- mean_t_ci(d, clip_01 = FALSE)
    paired_rows[[length(paired_rows)+1L]] <- data.frame(
      model = model, metric = metric,
      contrast = "full minus chemical", mean_difference = ci["mean"],
      descriptive_95lo = ci["lower"], descriptive_95hi = ci["upper"])
  }
}
paired_table <- if (length(paired_rows)) do.call(rbind, paired_rows) else data.frame()

write.csv(folds, file.path(OUTPUT_DIR, "outer_fold_metrics.csv"), row.names = FALSE)
write.csv(predictions, file.path(OUTPUT_DIR, "outer_test_predictions.csv"), row.names = FALSE)
write.csv(tuning_table, file.path(OUTPUT_DIR, "inner_cv_all_configurations.csv"), row.names = FALSE)
write.csv(summary_table, file.path(OUTPUT_DIR, "performance_summary.csv"), row.names = FALSE)
write.csv(class_table, file.path(OUTPUT_DIR, "class_specific_metrics.csv"), row.names = FALSE)
write.csv(class_fold_table, file.path(OUTPUT_DIR, "class_specific_by_fold.csv"), row.names = FALSE)
write.csv(paired_table, file.path(OUTPUT_DIR, "paired_full_vs_chemical.csv"), row.names = FALSE)
write.csv(confusion_table, file.path(OUTPUT_DIR, "confusion_counts.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
saveRDS(list(config = list(modes = ANALYSIS_MODES, models = MODELS,
           derivative_grid = DERIVATIVE_GRID, pve_grid = PVE_GRID,
           outer_folds = OUTER_FOLDS, inner_folds = INNER_FOLDS, seed = SEED),
             folds = folds, predictions = predictions, tuning = tuning_table,
             summary = summary_table, confusion = confusion_table,
             classes = class_table, paired_domains = paired_table),
        file.path(OUTPUT_DIR, "complete_breed_ncv_results.rds"))
cat("\n================ SUMMARY (values from 0 to 1) ================\n")
print(summary_table, row.names = FALSE)
cat("\nResults saved in: ", normalizePath(OUTPUT_DIR), "\n", sep = "")
cat("Note: 95% intervals describe variation across sample-level folds, not independent farms.\n")
