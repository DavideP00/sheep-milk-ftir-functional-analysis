# Shared functions for altitude classification, breed classification and outlier sensitivity.

`%||%` <- function(x, y) if (is.null(x)) y else x

ALTITUDE_LEVELS <- c("Plain", "Hill", "Mountain")
SUPPORTED_MODELS <- c("multinom", "ord_cumul", "ord_adj", "enet", "svm")

# Input standardization and dataset metadata.
normalizza_dataset <- function(data1, data2) {
  requireNamespace("dplyr")
  requireNamespace("stringr")

  prepara_base <- function(df) {
    df |>
      dplyr::rename(Azienda = 1, Altezza = 2, Matricola = 3) |>
      dplyr::rename_with(~ stringr::str_remove(.x, "VAR_"), dplyr::starts_with("VAR_")) |>
      dplyr::mutate(dplyr::across(4:dplyr::last_col(), as.numeric)) |>
      dplyr::mutate(
        Altezza = dplyr::case_when(
          stringr::str_detect(Altezza, stringr::regex("Pianura", ignore_case = TRUE)) ~ "Plain",
          stringr::str_detect(Altezza, stringr::regex("Montagna", ignore_case = TRUE)) ~ "Mountain",
          stringr::str_detect(Altezza, stringr::regex("Collina", ignore_case = TRUE)) ~ "Hill",
          TRUE ~ as.character(Altezza)
        ),
        Altezza = factor(Altezza, levels = ALTITUDE_LEVELS, ordered = TRUE)
      )
  }

  d1 <- prepara_base(data1)
  d2 <- prepara_base(data2)

  map <- dplyr::bind_rows(
    dplyr::select(d1, Azienda, Altezza),
    dplyr::select(d2, Azienda, Altezza)
  ) |>
    dplyr::distinct(Azienda, Altezza) |>
    dplyr::arrange(Altezza, Azienda) |>
    dplyr::group_by(Altezza) |>
    dplyr::mutate(Codice_Anonimo = paste0(substr(as.character(Altezza), 1, 1), dplyr::row_number())) |>
    dplyr::ungroup() |>
    dplyr::select(Azienda, Codice_Anonimo)

  anonymise <- function(df) {
    df |>
      dplyr::left_join(map, by = "Azienda") |>
      dplyr::select(-Azienda) |>
      dplyr::rename(Azienda = Codice_Anonimo) |>
      dplyr::relocate(Azienda, .before = Altezza)
  }

  list(d1 = anonymise(d1), d2 = anonymise(d2))
}

f1_macro <- function(y_true, y_pred, levels = ALTITUDE_LEVELS) {
  yt <- as.character(y_true)
  yp <- as.character(y_pred)
  f1 <- vapply(levels, function(cl) {
    tp <- sum(yt == cl & yp == cl, na.rm = TRUE)
    fp <- sum(yt != cl & yp == cl, na.rm = TRUE)
    fn <- sum(yt == cl & yp != cl, na.rm = TRUE)
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  }, numeric(1))
  mean(f1)
}

get_spectral_columns <- function(data, col_start = NULL, col_end = NULL) {
  cols <- grep("^[0-9]+(\\.[0-9]+)?$", names(data), value = TRUE)
  if (length(cols) == 0L) stop("No numeric spectral columns found.")
  nums <- as.numeric(cols)
  ord <- order(nums)
  cols <- cols[ord]
  nums <- nums[ord]

  if (!is.null(col_start)) {
    keep <- nums >= as.numeric(col_start)
    cols <- cols[keep]
    nums <- nums[keep]
  }
  if (!is.null(col_end)) {
    keep <- nums <= as.numeric(col_end)
    cols <- cols[keep]
  }
  cols
}

ftir_domain_intervals <- function(domain = c("chemical", "full")) {
  domain <- match.arg(domain)
  if (domain == "full") return(list(c(240, 1299)))
  list(
    c(240, 297), c(298, 361), c(362, 387), c(388, 439),
    c(440, 453), c(725, 776), c(777, 959)
  )
}

derivative_views_from_config <- function(
    config = c("m0", "m0_d1", "m0_d2", "m0_d1_d2")) {
  config <- match.arg(config)
  switch(
    config,
    m0 = "m0",
    m0_d1 = c("m0", "m1"),
    m0_d2 = c("m0", "m2"),
    m0_d1_d2 = c("m0", "m1", "m2")
  )
}

assert_group_target_nesting <- function(data,
                                        group_col = "Azienda",
                                        target_col = "Altezza") {
  requireNamespace("dplyr")
  chk <- data |>
    dplyr::distinct(.data[[group_col]], .data[[target_col]]) |>
    dplyr::count(.data[[group_col]], name = "n_targets")
  bad <- chk[chk$n_targets != 1L, , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop("Each farm must belong to exactly one altitude class: ",
         paste(bad[[group_col]], collapse = ", "))
  }
  invisible(TRUE)
}

mean_t_ci <- function(x, conf_level = 0.95, clip_01 = TRUE) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  m <- mean(x)
  if (n < 2L) return(c(mean = m, lower = NA_real_, upper = NA_real_))
  se <- stats::sd(x) / sqrt(n)
  q <- stats::qt((1 + conf_level) / 2, df = n - 1L)
  lo <- m - q * se
  hi <- m + q * se
  if (clip_01) {
    lo <- max(0, lo)
    hi <- min(1, hi)
  }
  c(mean = m, lower = lo, upper = hi)
}

# Class-specific metrics and paired model comparisons.

normalise_model_names <- function(models) {
  x <- tolower(trimws(as.character(models)))
  x <- gsub("[- ]+", "_", x)
  aliases <- c(
    "elastic_net" = "enet",
    "elasticnet" = "enet",
    "glmnet" = "enet",
    "multinomial" = "multinom",
    "ordinal_cumulative" = "ord_cumul",
    "ordinal_adjacent" = "ord_adj"
  )
  hit <- x %in% names(aliases)
  x[hit] <- unname(aliases[x[hit]])
  bad <- setdiff(unique(x), SUPPORTED_MODELS)
  if (length(bad) > 0L) {
    stop("Unsupported model(s): ", paste(bad, collapse = ", "),
         ". Supported: ", paste(SUPPORTED_MODELS, collapse = ", "))
  }
  unique(x)
}

class_metrics_from_confusion <- function(cm, levels = ALTITUDE_LEVELS) {
  cm2 <- matrix(
    0,
    nrow = length(levels),
    ncol = length(levels),
    dimnames = list(True = levels, Predicted = levels)
  )
  if (!is.null(cm) && length(cm) > 0L) {
    rr <- intersect(rownames(cm), levels)
    cc <- intersect(colnames(cm), levels)
    cm2[rr, cc] <- cm[rr, cc, drop = FALSE]
  }

  total <- sum(cm2)
  out <- lapply(levels, function(cl) {
    tp <- unname(cm2[cl, cl])
    fn <- sum(cm2[cl, ]) - tp
    fp <- sum(cm2[, cl]) - tp
    tn <- total - tp - fn - fp
    den_p <- tp + fp
    den_r <- tp + fn
    precision <- if (den_p > 0) tp / den_p else NA_real_
    recall <- if (den_r > 0) tp / den_r else NA_real_
    f1 <- if (is.finite(precision) && is.finite(recall) &&
              (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else if (is.finite(precision) && is.finite(recall) &&
               precision == 0 && recall == 0) {
      0
    } else {
      NA_real_
    }
    data.frame(
      class = cl,
      support = den_r,
      predicted_n = den_p,
      tp = tp, fp = fp, fn = fn, tn = tn,
      precision = precision,
      recall = recall,
      f1 = f1,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

class_metrics_from_vectors <- function(truth, pred, levels = ALTITUDE_LEVELS) {
  cm <- table(
    True = factor(as.character(truth), levels = levels),
    Predicted = factor(as.character(pred), levels = levels)
  )
  class_metrics_from_confusion(cm, levels)
}

sample_fold_class_uncertainty <- function(predictions,
                                          conf_level = 0.95,
                                          levels = ALTITUDE_LEVELS) {
  folds <- sort(unique(predictions$fold))
  fold_rows <- lapply(folds, function(f) {
    d <- predictions[predictions$fold == f, , drop = FALSE]
    m <- class_metrics_from_vectors(d$truth, d$pred, levels)
    m$fold <- f
    m
  })
  by_fold <- dplyr::bind_rows(fold_rows)

  pooled <- class_metrics_from_vectors(predictions$truth, predictions$pred, levels)

  summary_rows <- lapply(levels, function(cl) {
    z <- by_fold[by_fold$class == cl, , drop = FALSE]
    out <- pooled[pooled$class == cl, c(
      "class", "support", "predicted_n", "precision", "recall", "f1"
    ), drop = FALSE]
    names(out)[names(out) == "precision"] <- "pooled_precision"
    names(out)[names(out) == "recall"] <- "pooled_recall"
    names(out)[names(out) == "f1"] <- "pooled_f1"

    for (metric in c("precision", "recall", "f1")) {
      vals <- z[[metric]]
      ci <- mean_t_ci(vals, conf_level = conf_level, clip_01 = TRUE)
      out[[paste0(metric, "_fold_mean")]] <- unname(ci["mean"])
      out[[paste0(metric, "_ci_lower")]] <- unname(ci["lower"])
      out[[paste0(metric, "_ci_upper")]] <- unname(ci["upper"])
      out[[paste0(metric, "_sd")]] <- if (sum(is.finite(vals)) > 1L) {
        stats::sd(vals, na.rm = TRUE)
      } else NA_real_
      out[[paste0(metric, "_n_valid_folds")]] <- sum(is.finite(vals))
    }
    out
  })

  list(
    by_fold = by_fold,
    summary = dplyr::bind_rows(summary_rows),
    ci_note = paste0(
      "Class-specific point estimates are pooled out-of-fold metrics. ",
      "95% CIs are descriptive t intervals across outer sample folds; ",
      "sample-level folds are not independent farms."
    )
  )
}

paired_signflip_pvalue <- function(differences,
                                   exact_max_n = 18L,
                                   mc_reps = 100000L,
                                   seed = 123) {
  d <- differences[is.finite(differences)]
  n <- length(d)
  if (n == 0L) return(NA_real_)
  if (all(abs(d) < 1e-15)) return(1)
  obs <- abs(mean(d))

  if (n <= exact_max_n) {
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
    perm_stats <- abs(rowMeans(sweep(signs, 2, d, FUN = "*")))
    return(mean(perm_stats >= (obs - 1e-15)))
  }

  set.seed(seed)
  ge <- 0L
  done <- 0L
  chunk <- 5000L
  while (done < mc_reps) {
    m <- min(chunk, mc_reps - done)
    signs <- matrix(sample(c(-1, 1), size = m * n, replace = TRUE),
                    nrow = m, ncol = n)
    st <- abs(rowMeans(sweep(signs, 2, d, FUN = "*")))
    ge <- ge + sum(st >= (obs - 1e-15))
    done <- done + m
  }
  (ge + 1) / (mc_reps + 1)
}

paired_difference_summary <- function(a,
                                      b,
                                      conf_level = 0.95,
                                      seed = 123) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]
  b <- b[ok]
  d <- a - b
  ci <- mean_t_ci(d, conf_level = conf_level, clip_01 = FALSE)
  data.frame(
    n_pairs = length(d),
    mean_model_a = if (length(a)) mean(a) else NA_real_,
    mean_model_b = if (length(b)) mean(b) else NA_real_,
    mean_difference_a_minus_b = unname(ci["mean"]),
    difference_ci_lower = unname(ci["lower"]),
    difference_ci_upper = unname(ci["upper"]),
    median_difference_a_minus_b = if (length(d)) stats::median(d) else NA_real_,
    signflip_p_value = paired_signflip_pvalue(d, seed = seed),
    stringsAsFactors = FALSE
  )
}

compare_paper_models <- function(model_results,
                                 conf_level = 0.95,
                                 seed = 123) {
  models <- names(model_results)
  if (length(models) < 2L) return(data.frame())
  pairs <- combn(models, 2, simplify = FALSE)
  domains <- Reduce(intersect, lapply(model_results, names))
  out <- list()

  for (dom in domains) {
    for (pp in pairs) {
      a_name <- pp[1]
      b_name <- pp[2]
      a_pred <- model_results[[a_name]][[dom]]$predictions
      b_pred <- model_results[[b_name]][[dom]]$predictions

      key_a <- paste(a_pred$fold, a_pred$sample_id, sep = "::")
      key_b <- paste(b_pred$fold, b_pred$sample_id, sep = "::")
      if (!setequal(key_a, key_b)) {
        stop("Sample outer partitions differ between ", a_name, " and ", b_name,
             " for domain ", dom, ".")
      }

      a <- model_results[[a_name]][[dom]]$fold_metrics
      b <- model_results[[b_name]][[dom]]$fold_metrics
      m <- merge(
        a[, c("fold", "accuracy", "macro_f1")],
        b[, c("fold", "accuracy", "macro_f1")],
        by = "fold",
        suffixes = c("_a", "_b"),
        all = FALSE
      )

      for (metric in c("accuracy", "macro_f1")) {
        rr <- paired_difference_summary(
          m[[paste0(metric, "_a")]],
          m[[paste0(metric, "_b")]],
          conf_level = conf_level,
          seed = seed + length(out) + 1L
        )
        rr$design <- "paper_sample_cv"
        rr$domain <- dom
        rr$metric <- metric
        rr$model_a <- a_name
        rr$model_b <- b_name
        rr$interpretation <- "positive difference = model_a higher"
        rr$inference_note <- paste0(
          "Paired by identical outer sample fold; sign-flip p-value and CI are descriptive ",
          "because sample folds are not independent farms."
        )
        out[[length(out) + 1L]] <- rr
      }
    }
  }
  dplyr::bind_rows(out)
}

print_paired_model_comparison <- function(x, title) {
  cat("\n\n================ ", title, " ================\n", sep = "")
  if (is.null(x) || nrow(x) == 0L) {
    cat("At least two models are required for a paired comparison.\n")
    return(invisible(x))
  }
  print(x, row.names = FALSE)
  invisible(x)
}

# Training-only, dense-grid FPCA.
fit_dense_fpca_pair_fast <- function(train_df,
                                     test_df,
                                     spectral_cols,
                                     max_pve = 0.9999) {
  Xtr <- as.matrix(train_df[, spectral_cols, drop = FALSE])
  Xte <- as.matrix(test_df[, spectral_cols, drop = FALSE])

  storage.mode(Xtr) <- "double"
  storage.mode(Xte) <- "double"

  if (nrow(Xtr) < 3L || ncol(Xtr) < 3L) return(NULL)
  if (anyNA(Xtr) || anyNA(Xte)) stop("NA in spectral matrix before FPCA.")

  mu <- colMeans(Xtr)
  Xtrc <- sweep(Xtr, 2, mu, "-")
  Xtec <- sweep(Xte, 2, mu, "-")

  n <- nrow(Xtrc)
  p <- ncol(Xtrc)
  total_var <- sum(Xtrc^2) / (n - 1)
  if (!is.finite(total_var) || total_var <= 0) return(NULL)

  if (p <= n) {
    C <- crossprod(Xtrc) / (n - 1)
    eg <- eigen(C, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    keep <- vals > max(vals, na.rm = TRUE) * .Machine$double.eps * max(n, p)
    vals <- vals[keep]
    V_all <- eg$vectors[, keep, drop = FALSE]

    cum <- cumsum(vals) / total_var
    kmax <- which(cum >= max_pve)[1]
    if (is.na(kmax)) kmax <- length(vals)
    vals <- vals[seq_len(kmax)]
    V <- V_all[, seq_len(kmax), drop = FALSE]
  } else {
    G <- tcrossprod(Xtrc) / (n - 1)
    eg <- eigen(G, symmetric = TRUE)
    vals_all <- pmax(eg$values, 0)
    keep <- vals_all > max(vals_all, na.rm = TRUE) * .Machine$double.eps * max(n, p)
    vals_all <- vals_all[keep]
    U_all <- eg$vectors[, keep, drop = FALSE]

    cum_all <- cumsum(vals_all) / total_var
    kmax <- which(cum_all >= max_pve)[1]
    if (is.na(kmax)) kmax <- length(vals_all)

    vals <- vals_all[seq_len(kmax)]
    U <- U_all[, seq_len(kmax), drop = FALSE]

    denom <- sqrt((n - 1) * vals)
    V <- crossprod(Xtrc, U)
    V <- sweep(V, 2, denom, "/")

    vn <- sqrt(colSums(V^2))
    good <- is.finite(vn) & vn > 0
    V <- V[, good, drop = FALSE]
    vals <- vals[good]
  }

  if (ncol(V) == 0L) return(NULL)
  train_scores <- Xtrc %*% V
  test_scores <- Xtec %*% V
  cum_pve <- cumsum(vals) / total_var

  list(
    train_scores_all = train_scores,
    test_scores_all = test_scores,
    cum_pve = pmin(cum_pve, 1),
    mean_train = mu,
    loadings_train = V
  )
}

fit_fpca_view_pair_cache <- function(train_df,
                                     test_df,
                                     intervals,
                                     max_pve = 0.9999) {
  spec <- get_spectral_columns(train_df)
  nums <- as.numeric(spec)
  out <- vector("list", length(intervals))

  for (i in seq_along(intervals)) {
    rng <- intervals[[i]]
    cols <- spec[nums >= rng[1] & nums <= rng[2]]
    if (length(cols) < 3L) next
    fit <- fit_dense_fpca_pair_fast(train_df, test_df, cols, max_pve = max_pve)
    if (is.null(fit)) next
    fit$interval <- rng
    out[[i]] <- fit
  }

  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) return(NULL)
  out
}

k_from_pve <- function(cum_pve, pve) {
  k <- which(cum_pve >= pve)[1]
  if (is.na(k)) k <- length(cum_pve)
  as.integer(k)
}

build_preprocessed_fpca_cache <- function(train_views,
                                          test_views,
                                          domains,
                                          views_needed,
                                          max_pve,
                                          id_col = "Matricola",
                                          group_col = "Azienda",
                                          target_col = "Altezza") {
  meta_cols <- intersect(c(id_col, group_col, target_col), names(train_views[[views_needed[1]]]))
  train_meta <- train_views[[views_needed[1]]][, meta_cols, drop = FALSE]
  test_meta <- test_views[[views_needed[1]]][, meta_cols, drop = FALSE]

  for (v in views_needed) {
    if (!identical(as.character(train_views[[v]][[id_col]]), as.character(train_meta[[id_col]]))) {
      stop("Training row order differs across spectral views.")
    }
    if (!identical(as.character(test_views[[v]][[id_col]]), as.character(test_meta[[id_col]]))) {
      stop("Test row order differs across spectral views.")
    }
  }

  cache <- list(train_meta = train_meta, test_meta = test_meta, domains = list())
  for (domain in domains) {
    ints <- ftir_domain_intervals(domain)
    cache$domains[[domain]] <- list()
    for (v in views_needed) {
      cache$domains[[domain]][[v]] <- fit_fpca_view_pair_cache(
        train_views[[v]], test_views[[v]], ints, max_pve = max_pve
      )
      if (is.null(cache$domains[[domain]][[v]])) {
        stop("FPCA failed for domain=", domain, ", view=", v)
      }
    }
  }
  cache
}

assemble_fpca_features <- function(cache,
                                   domain,
                                   derivative_config,
                                   pve,
                                   target_col = "Altezza",
                                   scale_minmax = TRUE) {
  views <- derivative_views_from_config(derivative_config)
  tr_blocks <- list()
  te_blocks <- list()
  nm <- character()

  for (v in views) {
    bundles <- cache$domains[[domain]][[v]]
    if (is.null(bundles)) return(NULL)
    for (j in seq_along(bundles)) {
      b <- bundles[[j]]
      k <- k_from_pve(b$cum_pve, pve)
      idx <- seq_len(k)
      tr_blocks[[length(tr_blocks) + 1L]] <- b$train_scores_all[, idx, drop = FALSE]
      te_blocks[[length(te_blocks) + 1L]] <- b$test_scores_all[, idx, drop = FALSE]
      nm <- c(nm, paste0(v, "_I", j, "_FPC", idx))
    }
  }

  if (length(tr_blocks) == 0L) return(NULL)
  Xtr <- do.call(cbind, tr_blocks)
  Xte <- do.call(cbind, te_blocks)
  colnames(Xtr) <- nm
  colnames(Xte) <- nm

  Xtr <- as.data.frame(Xtr, check.names = FALSE)
  Xte <- as.data.frame(Xte, check.names = FALSE)

  if (scale_minmax) {
    sc <- fit_minmax_grouped(Xtr)
    Xtr <- apply_minmax_grouped(Xtr, sc)
    Xte <- apply_minmax_grouped(Xte, sc)
  } else {
    sc <- NULL
  }

  list(
    x_train = Xtr,
    y_train = cache$train_meta[[target_col]],
    x_test = Xte,
    y_test = cache$test_meta[[target_col]],
    train_meta = cache$train_meta,
    test_meta = cache$test_meta,
    n_features = ncol(Xtr),
    scaler = sc
  )
}

fit_minmax_grouped <- function(x_train) {
  x <- as.data.frame(x_train, check.names = FALSE)
  mins <- vapply(x, min, numeric(1), na.rm = TRUE)
  maxs <- vapply(x, max, numeric(1), na.rm = TRUE)
  ranges <- maxs - mins
  keep <- is.finite(mins) & is.finite(maxs) & is.finite(ranges) & ranges > 0
  if (!any(keep)) stop("No non-constant FPCA features remain.")
  list(min = mins[keep], range = ranges[keep], keep = names(ranges)[keep])
}

apply_minmax_grouped <- function(x, scaler) {
  x <- as.data.frame(x, check.names = FALSE)
  x <- x[, scaler$keep, drop = FALSE]
  z <- sweep(as.matrix(x), 2, scaler$min, "-")
  z <- sweep(z, 2, scaler$range, "/")
  as.data.frame(z, check.names = FALSE)
}

# Peak-width estimation for training-derived Savitzky-Golay windows.
calc_fwhm_from_matrix <- function(X, width_range = NULL) {
  requireNamespace("pracma")
  ref <- apply(X, 2, stats::median, na.rm = TRUE)
  threshold <- (max(ref, na.rm = TRUE) - min(ref, na.rm = TRUE)) * 0.05
  peaks <- pracma::findpeaks(
    ref, nups = 5,
    threshold = min(ref, na.rm = TRUE) + threshold
  )
  if (is.null(peaks)) return(NULL)

  widths <- numeric(nrow(peaks))
  for (i in seq_len(nrow(peaks))) {
    idx_max <- peaks[i, 2]
    val_max <- peaks[i, 1]
    left_peak <- peaks[i, 3]
    right_peak <- peaks[i, 4]
    baseline <- min(ref[left_peak:right_peak], na.rm = TRUE)
    half <- baseline + (val_max - baseline) / 2
    l <- idx_max
    r <- idx_max
    while (l > 1L && ref[l] > half) l <- l - 1L
    while (r < length(ref) && ref[r] > half) r <- r + 1L
    widths[i] <- r - l
  }
  widths <- widths[is.finite(widths)]
  if (!is.null(width_range)) {
    widths <- widths[widths > width_range[1] & widths < width_range[2]]
  }
  if (length(widths) == 0L) NULL else widths
}

estimate_sgolay_window_train <- function(train_data,
                                         col_start = "240",
                                         col_end = "1299",
                                         polynomial_degree = 4,
                                         fallback_window = 7,
                                         width_range = NULL) {
  spec <- get_spectral_columns(train_data, col_start, col_end)
  X <- as.matrix(train_data[, spec, drop = FALSE])
  widths <- tryCatch(calc_fwhm_from_matrix(X, width_range = width_range),
                     error = function(e) NULL)
  n_window <- if (is.null(widths)) fallback_window else floor(stats::median(widths, na.rm = TRUE))
  if (!is.finite(n_window)) n_window <- fallback_window
  if (n_window %% 2L == 0L) n_window <- n_window + 1L
  if (n_window <= polynomial_degree) {
    n_window <- polynomial_degree + 2L
    if (n_window %% 2L == 0L) n_window <- n_window + 1L
  }
  as.integer(n_window)
}

apply_sgolay_m0_m1_m2 <- function(data,
                                     col_start = "240",
                                     col_end = "1299",
                                     window_length,
                                     polynomial_degree = 4) {
  requireNamespace("signal")
  if (polynomial_degree < 2L) {
    stop("Savitzky-Golay polynomial_degree must be >= 2 to compute the second derivative.")
  }

  spec <- get_spectral_columns(data, col_start, col_end)
  X <- as.matrix(data[, spec, drop = FALSE])
  meta <- data[, setdiff(names(data), spec), drop = FALSE]

  filt <- function(m) {
    Xm <- t(vapply(seq_len(nrow(X)), function(i) {
      signal::sgolayfilt(X[i, ], p = polynomial_degree, n = window_length, m = m)
    }, FUN.VALUE = numeric(ncol(X))))
    colnames(Xm) <- spec
    cbind(meta, as.data.frame(Xm, check.names = FALSE))
  }

  list(m0 = filt(0L), m1 = filt(1L), m2 = filt(2L))
}

apply_sgolay_m0_m1 <- function(...) {
  apply_sgolay_m0_m1_m2(...)
}

preprocess_grouped_split <- function(raw_train,
                                     raw_test,
                                     col_start = "240",
                                     col_end = "1299",
                                     polynomial_degree = 4,
                                     fixed_window = NULL,
                                     fallback_window = 7,
                                     width_range = NULL) {
  w <- fixed_window %||% estimate_sgolay_window_train(
    raw_train, col_start, col_end, polynomial_degree, fallback_window,
    width_range = width_range
  )
  list(
    train = apply_sgolay_m0_m1_m2(raw_train, col_start, col_end, w, polynomial_degree),
    test = apply_sgolay_m0_m1_m2(raw_test, col_start, col_end, w, polynomial_degree),
    window_length = w
  )
}

empty_model_grid_row <- function(n) {
  data.frame(
    param_id = seq_len(n),
    kernel = rep(NA_character_, n),
    cost = rep(NA_real_, n),
    gamma_multiplier = rep(NA_real_, n),
    alpha = rep(NA_real_, n),
    lambda = rep(NA_real_, n),
    decay = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
}

default_model_grid <- function(model_type) {
  if (!model_type %in% SUPPORTED_MODELS) stop("Unsupported model: ", model_type)

  if (model_type == "svm") {
    lin <- expand.grid(
      kernel = "linear",
      cost = c(0.1, 1, 10),
      gamma_multiplier = NA_real_,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    rad <- expand.grid(
      kernel = "radial",
      cost = c(1, 10),
      gamma_multiplier = c(0.25, 1),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    z <- dplyr::bind_rows(lin, rad)
    out <- empty_model_grid_row(nrow(z))
    out$kernel <- z$kernel
    out$cost <- z$cost
    out$gamma_multiplier <- z$gamma_multiplier
    return(out)
  }

  if (model_type == "enet") {
    z <- expand.grid(
      alpha = c(0, 0.5, 1),
      lambda = c(0.001, 0.01, 0.1),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    out <- empty_model_grid_row(nrow(z))
    out$alpha <- z$alpha
    out$lambda <- z$lambda
    return(out)
  }

  if (model_type == "multinom") {
    z <- data.frame(decay = c(0, 0.001, 0.01, 0.1))
    out <- empty_model_grid_row(nrow(z))
    out$decay <- z$decay
    return(out)
  }

  empty_model_grid_row(1L)
}

model_param_label <- function(model_type, row) {
  if (model_type == "svm") {
    if (row$kernel == "linear") return(paste0("kernel=linear,C=", row$cost))
    return(paste0("kernel=radial,C=", row$cost, ",gamma_mult=", row$gamma_multiplier))
  }
  if (model_type == "enet") return(paste0("alpha=", row$alpha, ",lambda=", row$lambda))
  if (model_type == "multinom") return(paste0("decay=", row$decay))
  if (model_type == "ord_cumul") return("ridge cumulative logit; lambda selected by training AIC")
  if (model_type == "ord_adj") return("bias-reduced adjacent-category logit")
  ""
}

fit_model_candidate <- function(x, y, model_type, param_row, seed = 123) {
  if (!model_type %in% SUPPORTED_MODELS) stop("Unsupported model: ", model_type)
  set.seed(seed)
  yu <- factor(as.character(y), levels = ALTITUDE_LEVELS)
  yo <- factor(as.character(y), levels = ALTITUDE_LEVELS, ordered = TRUE)
  X <- as.matrix(x)

  if (model_type == "svm") {
    requireNamespace("e1071")
    if (param_row$kernel == "linear") {
      mod <- e1071::svm(x = X, y = yu, kernel = "linear", cost = param_row$cost,
                        scale = FALSE, probability = FALSE)
      return(list(type = model_type, model = mod, param = param_row))
    }
    gamma <- as.numeric(param_row$gamma_multiplier) / max(1, ncol(X))
    mod <- e1071::svm(x = X, y = yu, kernel = "radial", cost = param_row$cost,
                      gamma = gamma, scale = FALSE, probability = FALSE)
    return(list(type = model_type, model = mod, param = param_row, gamma = gamma))
  }

  if (model_type == "enet") {
    requireNamespace("glmnet")
    mod <- glmnet::glmnet(
      x = X, y = yu, family = "multinomial",
      alpha = param_row$alpha, lambda = param_row$lambda,
      standardize = FALSE
    )
    return(list(type = model_type, model = mod, param = param_row))
  }

  if (model_type == "multinom") {
    requireNamespace("nnet")
    df <- data.frame(X, y = yu, check.names = TRUE)
    mod <- nnet::multinom(
      y ~ ., data = df, decay = param_row$decay,
      trace = FALSE, MaxNWts = 10000, maxit = 1000
    )
    return(list(type = model_type, model = mod, param = param_row,
                columns = colnames(df)[colnames(df) != "y"]))
  }

  if (model_type == "ord_cumul") {
    requireNamespace("ordinalNet")
    mod <- ordinalNet::ordinalNet(
      x = X, y = yo, family = "cumulative", link = "logit",
      parallelTerms = FALSE, nonparallelTerms = TRUE,
      alpha = 0, warn = FALSE
    )
    idx <- which.min(mod$aic)
    return(list(type = model_type, model = mod, param = param_row, lambda_idx = idx))
  }

  if (model_type == "ord_adj") {
    requireNamespace("brglm2")
    df <- data.frame(X, y = yo, check.names = TRUE)
    mod <- brglm2::bracl(
      y ~ ., data = df, parallel = FALSE,
      control = brglm2::brglm_control(slowit = 0.5, maxit = 1000, epsilon = 1e-4)
    )
    return(list(type = model_type, model = mod, param = param_row,
                columns = colnames(df)[colnames(df) != "y"]))
  }
}

predict_model_candidate <- function(fit, newx) {
  X <- as.matrix(newx)

  if (fit$type == "svm") {
    p <- predict(fit$model, X)
    return(factor(as.character(p), levels = ALTITUDE_LEVELS))
  }

  if (fit$type == "enet") {
    p <- predict(fit$model, newx = X, type = "class")
    p <- drop(p)
    return(factor(as.character(p), levels = ALTITUDE_LEVELS))
  }

  if (fit$type == "multinom") {
    p <- predict(fit$model, newdata = data.frame(X, check.names = TRUE))
    return(factor(as.character(p), levels = ALTITUDE_LEVELS))
  }

  if (fit$type == "ord_cumul") {
    allp <- predict(fit$model, newx = X, type = "class")
    raw <- if (is.matrix(allp)) allp[, fit$lambda_idx] else allp
    if (is.numeric(raw)) raw <- ALTITUDE_LEVELS[raw]
    return(factor(as.character(raw), levels = ALTITUDE_LEVELS, ordered = TRUE))
  }

  if (fit$type == "ord_adj") {
    probs <- predict(fit$model, newdata = data.frame(X, check.names = TRUE), type = "probs")
    raw <- colnames(probs)[apply(probs, 1, which.max)]
    return(factor(raw, levels = ALTITUDE_LEVELS, ordered = TRUE))
  }

  stop("Unknown model type.")
}

# Nested sample-level cross-validation and training-only feature construction.

make_stratified_sample_folds <- function(y, k = 10, seed = 69) {
  requireNamespace("caret")
  y <- factor(as.character(y), levels = ALTITUDE_LEVELS)
  tab <- table(y)
  if (any(tab == 0L)) stop("All altitude classes must be present.")
  k_eff <- min(as.integer(k), as.integer(min(tab)))
  if (k_eff < 2L) stop("Not enough observations per class for stratified CV.")
  set.seed(seed)
  caret::createFolds(y, k = k_eff, list = TRUE, returnTrain = TRUE)
}

build_inner_sample_preprocess_cache <- function(raw_outer_train,
                                                inner_train_folds,
                                                col_start = "240",
                                                col_end = "1299",
                                                polynomial_degree = 4,
                                                fixed_sg_window = NULL) {
  n <- nrow(raw_outer_train)
  out <- vector("list", length(inner_train_folds))
  for (i in seq_along(inner_train_folds)) {
    tr <- inner_train_folds[[i]]
    va <- setdiff(seq_len(n), tr)
    out[[i]] <- list(
      processed = preprocess_grouped_split(
        raw_outer_train[tr, , drop = FALSE],
        raw_outer_train[va, , drop = FALSE],
        col_start = col_start,
        col_end = col_end,
        polynomial_degree = polynomial_degree,
        fixed_window = fixed_sg_window
      )
    )
  }
  out
}

build_inner_sample_fpca_cache <- function(inner_preprocess_cache,
                                          domain,
                                          pve_grid,
                                          target_col = "Altezza",
                                          id_col = "Matricola",
                                          group_col = "Azienda") {
  max_pve <- max(pve_grid)
  lapply(inner_preprocess_cache, function(sp) {
    build_preprocessed_fpca_cache(
      train_views = sp$processed$train,
      test_views = sp$processed$test,
      domains = domain,
      views_needed = c("m0", "m1"),
      max_pve = max_pve,
      id_col = id_col,
      group_col = group_col,
      target_col = target_col
    )
  })
}

evaluate_sample_candidate <- function(inner_fpca_cache,
                                      domain,
                                      pve,
                                      model_type,
                                      param_row,
                                      target_col = "Altezza",
                                      seed = 123) {
  acc <- numeric(length(inner_fpca_cache))
  f1 <- numeric(length(inner_fpca_cache))
  nfeat <- numeric(length(inner_fpca_cache))

  for (j in seq_along(inner_fpca_cache)) {
    feat <- assemble_fpca_features(
      inner_fpca_cache[[j]],
      domain = domain,
      derivative_config = "m0_d1",
      pve = pve,
      target_col = target_col,
      scale_minmax = TRUE
    )
    if (is.null(feat)) return(NULL)

    fit <- tryCatch(
      fit_model_candidate(feat$x_train, feat$y_train, model_type, param_row,
                          seed = seed + j),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    pred <- tryCatch(predict_model_candidate(fit, feat$x_test), error = function(e) NULL)
    if (is.null(pred)) return(NULL)

    truth <- factor(as.character(feat$y_test), levels = ALTITUDE_LEVELS)
    acc[j] <- mean(as.character(pred) == as.character(truth))
    f1[j] <- f1_macro(truth, pred)
    nfeat[j] <- feat$n_features
  }

  if (!all(is.finite(acc))) return(NULL)
  data.frame(
    domain = domain,
    derivative_config = "m0_d1",
    pve = pve,
    model = model_type,
    param_id = param_row$param_id,
    param_label = model_param_label(model_type, param_row),
    mean_inner_sample_accuracy = mean(acc),
    sd_inner_sample_accuracy = stats::sd(acc),
    mean_inner_macro_f1 = mean(f1),
    mean_n_features = mean(nfeat),
    stringsAsFactors = FALSE
  )
}

tune_sample_inner_pipeline <- function(inner_fpca_cache,
                                       domain,
                                       pve_grid,
                                       model_type,
                                       model_grid,
                                       target_col = "Altezza",
                                       seed = 123,
                                       verbose = TRUE) {
  rows <- list()
  z <- 0L
  total <- length(pve_grid) * nrow(model_grid)

  for (pve in pve_grid) {
    for (m in seq_len(nrow(model_grid))) {
      z <- z + 1L
      pr <- model_grid[m, , drop = FALSE]
      if (verbose) {
        message("    candidate ", z, "/", total,
                " | PVE=", pve,
                " | ", model_param_label(model_type, pr))
      }
      rr <- evaluate_sample_candidate(
        inner_fpca_cache = inner_fpca_cache,
        domain = domain,
        pve = pve,
        model_type = model_type,
        param_row = pr,
        target_col = target_col,
        seed = seed + z * 100L
      )
      if (!is.null(rr)) rows[[length(rows) + 1L]] <- rr
    }
  }

  tuning <- dplyr::bind_rows(rows)
  if (nrow(tuning) == 0L) stop("All inner sample-level configurations failed.")

  tuning <- tuning |>
    dplyr::arrange(
      dplyr::desc(mean_inner_sample_accuracy),
      sd_inner_sample_accuracy,
      mean_n_features,
      pve,
      param_id
    )

  list(best = tuning[1, , drop = FALSE], tuning = tuning)
}

summarise_paper_sample_domain <- function(fold_metrics,
                                          predictions,
                                          conf_level = 0.95) {
  ci_acc <- mean_t_ci(fold_metrics$accuracy, conf_level)
  ci_f1 <- mean_t_ci(fold_metrics$macro_f1, conf_level)
  class_unc <- sample_fold_class_uncertainty(
    predictions = predictions,
    conf_level = conf_level,
    levels = ALTITUDE_LEVELS
  )

  list(
    summary = data.frame(
      model = unique(fold_metrics$model),
      domain = unique(fold_metrics$domain),
      mean_accuracy = unname(ci_acc["mean"]),
      accuracy_ci_lower = unname(ci_acc["lower"]),
      accuracy_ci_upper = unname(ci_acc["upper"]),
      sd_accuracy = stats::sd(fold_metrics$accuracy),
      mean_macro_f1 = unname(ci_f1["mean"]),
      macro_f1_ci_lower = unname(ci_f1["lower"]),
      macro_f1_ci_upper = unname(ci_f1["upper"]),
      sd_macro_f1 = stats::sd(fold_metrics$macro_f1),
      pooled_sample_accuracy = mean(predictions$truth == predictions$pred),
      ci_level = conf_level,
      ci_note = paste0(
        "Descriptive t CI across sample-level outer folds; ",
        "the folds are not independent farms."
      ),
      stringsAsFactors = FALSE
    ),
    pooled_confusion_counts = table(
      True = factor(predictions$truth, levels = ALTITUDE_LEVELS),
      Predicted = factor(predictions$pred, levels = ALTITUDE_LEVELS)
    ),
    class_specific_metrics = class_unc$summary,
    class_specific_fold_metrics = class_unc$by_fold,
    class_specific_ci_note = class_unc$ci_note,
    pve_selection_counts = table(fold_metrics$selected_pve),
    model_parameter_selection_counts = table(fold_metrics$selected_param_label),
    sg_window_summary = summary(fold_metrics$sg_window)
  )
}

run_paper_sample_nested_model <- function(data,
                                          model_type,
                                          model_grid = default_model_grid(model_type),
                                          domains = c("full", "chemical"),
                                          pve_grid = c(0.95, 0.98, 0.99, 0.9999),
                                          outer_k_folds = 10,
                                          inner_k_folds = 5,
                                          col_start = "240",
                                          col_end = "1299",
                                          polynomial_degree = 4,
                                          fixed_sg_window = NULL,
                                          conf_level = 0.95,
                                          seed = 69,
                                          target_col = "Altezza",
                                          id_col = "Matricola",
                                          group_col = "Azienda",
                                          verbose = TRUE) {
  if (!model_type %in% SUPPORTED_MODELS) stop("Unsupported model: ", model_type)
  if (anyDuplicated(data[[id_col]]) > 0L) stop("Sample IDs must be unique.")

  outer_folds <- make_stratified_sample_folds(data[[target_col]], outer_k_folds, seed)
  fold_rows <- setNames(vector("list", length(domains)), domains)
  pred_rows <- setNames(vector("list", length(domains)), domains)
  tuning_logs <- setNames(vector("list", length(domains)), domains)
  for (d in domains) {
    fold_rows[[d]] <- vector("list", length(outer_folds))
    pred_rows[[d]] <- vector("list", length(outer_folds))
    tuning_logs[[d]] <- vector("list", length(outer_folds))
  }

  n <- nrow(data)
  for (f in seq_along(outer_folds)) {
    tr <- outer_folds[[f]]
    te <- setdiff(seq_len(n), tr)
    raw_train <- data[tr, , drop = FALSE]
    raw_test <- data[te, , drop = FALSE]
    fold_seed <- seed + f * 10000L

    if (verbose) {
      message("\n============================================================")
      message("PAPER SAMPLE OUTER FOLD ", f, "/", length(outer_folds),
              " | model=", model_type,
              " | train n=", nrow(raw_train), " | test n=", nrow(raw_test))
    }

    inner_folds <- make_stratified_sample_folds(
      raw_train[[target_col]], inner_k_folds, fold_seed + 1L
    )
    inner_pp <- build_inner_sample_preprocess_cache(
      raw_outer_train = raw_train,
      inner_train_folds = inner_folds,
      col_start = col_start,
      col_end = col_end,
      polynomial_degree = polynomial_degree,
      fixed_sg_window = fixed_sg_window
    )

    outer_pp <- preprocess_grouped_split(
      raw_train, raw_test,
      col_start = col_start,
      col_end = col_end,
      polynomial_degree = polynomial_degree,
      fixed_window = fixed_sg_window
    )

    for (domain in domains) {
      if (verbose) message("  DOMAIN: ", domain)

      inner_fpca <- build_inner_sample_fpca_cache(
        inner_preprocess_cache = inner_pp,
        domain = domain,
        pve_grid = pve_grid,
        target_col = target_col,
        id_col = id_col,
        group_col = group_col
      )

      tuned <- tune_sample_inner_pipeline(
        inner_fpca_cache = inner_fpca,
        domain = domain,
        pve_grid = pve_grid,
        model_type = model_type,
        model_grid = model_grid,
        target_col = target_col,
        seed = fold_seed,
        verbose = verbose
      )
      best <- tuned$best
      tuning_logs[[domain]][[f]] <- tuned$tuning

      if (verbose) {
        message("  SELECTED | PVE=", best$pve,
                " | ", best$param_label,
                " | mean inner sample acc=",
                round(best$mean_inner_sample_accuracy, 4))
      }

      outer_cache <- build_preprocessed_fpca_cache(
        train_views = outer_pp$train,
        test_views = outer_pp$test,
        domains = domain,
        views_needed = c("m0", "m1"),
        max_pve = as.numeric(best$pve),
        id_col = id_col,
        group_col = group_col,
        target_col = target_col
      )
      feat <- assemble_fpca_features(
        outer_cache,
        domain = domain,
        derivative_config = "m0_d1",
        pve = as.numeric(best$pve),
        target_col = target_col,
        scale_minmax = TRUE
      )

      pr <- model_grid[model_grid$param_id == best$param_id, , drop = FALSE]
      fit <- fit_model_candidate(
        feat$x_train, feat$y_train, model_type, pr,
        seed = fold_seed + 999L
      )
      pred <- predict_model_candidate(fit, feat$x_test)
      truth <- factor(as.character(feat$y_test), levels = ALTITUDE_LEVELS)
      acc <- mean(as.character(pred) == as.character(truth))
      mf1 <- f1_macro(truth, pred)

      fold_rows[[domain]][[f]] <- data.frame(
        fold = f,
        domain = domain,
        model = model_type,
        accuracy = acc,
        macro_f1 = mf1,
        n_test = length(truth),
        selected_pve = as.numeric(best$pve),
        selected_param_id = best$param_id,
        selected_param_label = best$param_label,
        inner_mean_sample_accuracy = best$mean_inner_sample_accuracy,
        n_features = feat$n_features,
        sg_window = outer_pp$window_length,
        stringsAsFactors = FALSE
      )

      pred_rows[[domain]][[f]] <- data.frame(
        fold = f,
        sample_id = as.character(feat$test_meta[[id_col]]),
        farm = if (group_col %in% names(feat$test_meta)) {
          as.character(feat$test_meta[[group_col]])
        } else NA_character_,
        truth = as.character(truth),
        pred = as.character(pred),
        stringsAsFactors = FALSE
      )

      if (verbose) {
        message("  OUTER SAMPLE-FOLD | acc=", round(acc, 4),
                " | macro-F1=", round(mf1, 4))
      }
    }
  }

  out <- list()
  for (domain in domains) {
    fm <- dplyr::bind_rows(fold_rows[[domain]])
    pp <- dplyr::bind_rows(pred_rows[[domain]])
    sm <- summarise_paper_sample_domain(fm, pp, conf_level)
    out[[domain]] <- c(
      list(fold_metrics = fm, predictions = pp, tuning_logs = tuning_logs[[domain]]),
      sm
    )
  }
  out
}

# Optional leave-one-farm-out validation utilities.

split_indices_for_farm <- function(data, test_farm, group_col = "Azienda") {
  ist <- as.character(data[[group_col]]) == as.character(test_farm)
  list(train = which(!ist), test = which(ist))
}

build_inner_preprocess_cache <- function(raw_outer_train,
                                         inner_plan,
                                         group_col = "Azienda",
                                         col_start = "240",
                                         col_end = "1299",
                                         polynomial_degree = 4,
                                         fixed_sg_window = NULL) {
  out <- vector("list", nrow(inner_plan))
  for (i in seq_len(nrow(inner_plan))) {
    val_farm <- inner_plan$group[i]
    idx <- split_indices_for_farm(raw_outer_train, val_farm, group_col)
    out[[i]] <- list(
      validation_farm = val_farm,
      validation_target = inner_plan$target[i],
      processed = preprocess_grouped_split(
        raw_outer_train[idx$train, , drop = FALSE],
        raw_outer_train[idx$test, , drop = FALSE],
        col_start, col_end, polynomial_degree, fixed_sg_window
      )
    )
  }
  out
}

build_inner_fpca_cache <- function(inner_preprocess_cache,
                                   domain_grid,
                                   derivative_configs,
                                   pve_grid,
                                   target_col = "Altezza",
                                   id_col = "Matricola",
                                   group_col = "Azienda",
                                   verbose = TRUE) {
  views_needed <- unique(unlist(lapply(derivative_configs, derivative_views_from_config)))
  max_pve <- max(pve_grid)
  out <- vector("list", length(inner_preprocess_cache))

  for (i in seq_along(inner_preprocess_cache)) {
    if (verbose) message("  building FPCA cache for inner farm ",
                         inner_preprocess_cache[[i]]$validation_farm, " ...")
    pp <- inner_preprocess_cache[[i]]$processed
    out[[i]] <- build_preprocessed_fpca_cache(
      train_views = pp$train,
      test_views = pp$test,
      domains = domain_grid,
      views_needed = views_needed,
      max_pve = max_pve,
      id_col = id_col, group_col = group_col, target_col = target_col
    )
    out[[i]]$validation_farm <- inner_preprocess_cache[[i]]$validation_farm
    out[[i]]$validation_target <- inner_preprocess_cache[[i]]$validation_target
  }
  out
}

make_farm_description_table <- function(data,
                                        breed_name,
                                        farm_col = "Azienda",
                                        altitude_col = "Altezza",
                                        id_col = "Matricola") {
  required <- c(farm_col, altitude_col, id_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in ", breed_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  data |>
    dplyr::mutate(
      .farm = as.character(.data[[farm_col]]),
      .altitude = as.character(.data[[altitude_col]]),
      .id = as.character(.data[[id_col]])
    ) |>
    dplyr::group_by(.farm, .altitude) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_unique_IDs = dplyr::n_distinct(.id),
      duplicate_ID_rows = n_samples - n_unique_IDs,
      .groups = "drop"
    ) |>
    dplyr::transmute(
      Breed = breed_name,
      Farm = .farm,
      Altitude_class = .altitude,
      n_samples = n_samples,
      n_unique_IDs = n_unique_IDs,
      duplicate_ID_rows = duplicate_ID_rows
    )
}

print_farm_description <- function(dati_sw,
                                   dati_vdb,
                                   farm_col = "Azienda",
                                   altitude_col = "Altezza",
                                   id_col = "Matricola") {
  farm_table <- dplyr::bind_rows(
    make_farm_description_table(
      dati_sw, "Sarda", farm_col, altitude_col, id_col
    ),
    make_farm_description_table(
      dati_vdb, "Valle del Belice", farm_col, altitude_col, id_col
    )
  ) |>
    dplyr::mutate(
      Altitude_class = factor(
        Altitude_class,
        levels = ALTITUDE_LEVELS,
        ordered = TRUE
      )
    ) |>
    dplyr::arrange(Breed, Altitude_class, Farm)

  cat("\n\n============================================================\n")
  cat("FARM DESCRIPTION TABLE\n")
  cat("============================================================\n")
  print(as.data.frame(farm_table), row.names = FALSE)

  cat("\nTotal farms: ", nrow(farm_table), "\n", sep = "")

  cat("\nFarms by breed and altitude:\n")
  tmp_farms <- farm_table |>
    dplyr::count(Breed, Altitude_class, name = "n_farms")
  print(as.data.frame(tmp_farms), row.names = FALSE)

  cat("\nSamples by breed and altitude:\n")
  tmp_samples <- farm_table |>
    dplyr::group_by(Breed, Altitude_class) |>
    dplyr::summarise(
      n_samples = sum(n_samples),
      .groups = "drop"
    )
  print(as.data.frame(tmp_samples), row.names = FALSE)

  if (any(farm_table$duplicate_ID_rows > 0L)) {
    cat("\nNOTE: duplicate_ID_rows > 0 indicates repeated values of the ID column '",
        id_col, "'. Interpret this only according to what that ID represents in the source data.\n",
        sep = "")
  } else {
    cat("\nNo duplicated values of '", id_col, "' were detected within farms.\n", sep = "")
  }

  invisible(farm_table)
}

selection_numeric_summary <- function(x, label = "n_features") {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(data.frame(
      variable = label,
      n = 0L,
      mean = NA_real_, sd = NA_real_, median = NA_real_,
      min = NA_real_, max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    variable = label,
    n = length(x),
    mean = mean(x),
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    median = stats::median(x),
    min = min(x),
    max = max(x),
    stringsAsFactors = FALSE
  )
}

print_selection_frequency <- function(x, title) {
  cat("\n", title, ":\n", sep = "")
  print(table(x, useNA = "ifany"))
}

print_paper_results <- function(res, breed_name, model_type) {
  cat("\n\n================ PAPER SAMPLE-CV | ", breed_name, " | ", model_type,
      " ================\n", sep = "")
  cat("No outlier removal. SG/FPCA/min-max are training-only; PVE/model parameters are inner-tuned.\n")

  for (dom in names(res)) {
    rr <- res[[dom]]
    fm <- rr$fold_metrics

    cat("\nDOMAIN: ", dom, "\n", sep = "")
    cat("Outer sample-fold metrics:\n")
    print(fm)

    cat("\nSummary + descriptive 95% CI across outer sample folds:\n")
    print(rr$summary)

    cat("\nPooled confusion counts:\n")
    print(rr$pooled_confusion_counts)

    cat("\nClass-specific Precision / Recall / F1 + uncertainty:\n")
    print(rr$class_specific_metrics, row.names = FALSE)
    cat("NOTE: ", rr$class_specific_ci_note, "\n", sep = "")

    cat("\nFPC feature-count summary across outer folds:\n")
    print(selection_numeric_summary(fm$n_features, "n_FPC_features"), row.names = FALSE)

    print_selection_frequency(fm$selected_pve, "Selected PVE frequency")
    print_selection_frequency(fm$sg_window, "Training-derived SG-window frequency")
    print_selection_frequency(fm$selected_param_label, "Selected model-parameter frequency")
  }
}

print_analysis_reproducibility <- function(settings,
                                           model_grids = NULL) {
  cat("\n\n============================================================\n")
  cat("ANALYSIS SETTINGS / REPRODUCIBILITY\n")
  cat("============================================================\n")

  for (nm in names(settings)) {
    cat("\n", nm, ":\n", sep = "")
    print(settings[[nm]])
  }

  if (!is.null(model_grids)) {
    cat("\nMODEL TUNING GRIDS USED:\n")
    for (nm in names(model_grids)) {
      cat("\nModel: ", nm, "\n", sep = "")
      print(model_grids[[nm]], row.names = FALSE)
    }
  }

  cat("\n============================================================\n")
  cat("R / PACKAGE SESSION INFORMATION\n")
  cat("============================================================\n\n")
  print(sessionInfo())

  invisible(NULL)
}

# Numerical helpers and sparse functional discriminant analysis.
vector_norm = function(x) {
  return(sqrt(sum(x^2)))
}

pdist2 = function(x,y) {
  dist.vec = matrix(nrow = nrow(x), ncol = ncol(x))
  for (i in c(1:nrow(x))) {
    dist.vec[i] = sqrt((x[i] - y)^2)
  }
  return(dist.vec)
}

wthresh = function(x,sorh,t){

  if (sorh == 's') {
    tmp = (abs(x) - t)
    tmp = (tmp+abs(tmp))/2
    y = sign(x) * tmp
  } else if (sorh == 'h') {
    y = x * (abs(x) > t)
  } else {
    print('Invalid argument value.')
  }
  return(y)
}

coordlasso = function(A, delta, lambda, tol, maxiter) {

  p = nrow(A)

  b_ini = solve(A)%*%delta
  b_ini = b_ini/vector_norm(b_ini)

  bnew = b_ini

  converged = FALSE
  step = 0
  dA = as.matrix(diag(A))

  while ((step < maxiter) & !converged) {
    bold = bnew
    step = step + 1

    for (j in c(1:p)) {
      x = delta[j] - A[j,]%*%bold + dA[j]*bold[j]

      x = wthresh(x,'s',lambda)
      x = x/dA[j]

      bold[j] = x
    }
    conv_criterion = vector_norm(bnew - bold)/vector_norm(bold)
    converged = conv_criterion < tol

    bnew = bold

    if (vector_norm(bnew) == 0) break

  }

  b = bnew

  if (vector_norm(b) != 0) {
    b = b / vector_norm(b)
  }

  return(b)
}

SDelta_binary = function(xtr0cv,xtr1cv) {

  ntr0cv = nrow(xtr0cv)
  ntr1cv = nrow(xtr1cv)
  p = ncol(xtr1cv)

  S = (ntr0cv - 1)*cov(xtr0cv) + (ntr1cv - 1)*cov(xtr1cv)
  S = S/(ntr0cv + ntr1cv - 2)
  mu0 = as.matrix(colMeans(xtr0cv))
  mu1 = as.matrix(colMeans(xtr1cv))

  return(list('S' = S, 'mu0' = mu0, 'mu1' = mu1))
}

SDelta_multi = function(xtr) {

  nclass = length(xtr)

  p = ncol(xtr[[1]])

  ntr = c()

  S = matrix(0,nrow = p, ncol= p)

  for (n in seq(nclass)){

    S = S + (nrow(xtr[[n]]) - 1)*cov(xtr[[n]])

    ntr = c(ntr, nrow(xtr[[n]]))

  }

  S = S/(sum(ntr) - nclass)

  mu_list = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    mu_list[[n]] = as.matrix(colMeans(xtr[[n]]))
  }

  A = do.call(cbind,mu_list)
  A_cov = cov(t(A))

  eigen_A_cov = eigen(A_cov)

  delta_list = vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    delta_list[[n]] = eigen_A_cov$vectors[,n]
  }

  return(list('S' = S, 'delta_list' = delta_list))
}

is_constant = function(xb){

  const = c()

  for (n in seq(ncol(xb))){
    const = c(const, var(xb[,n]))
  }

  if (sum(const) == 0){
    return(NULL)
  } else{
    return(xb[,which(const != 0), drop = FALSE])
  }
}

# Sparse functional discriminant analysis.
SFLDA = function(data, y, tau, lambda){

  if (!is.list(data)){
    stop("data should be given as list of matrices")
  }

  if (length(data) == 0){
    stop('list is empty')
  }

  if (!is.vector(y) && !is.factor(y)){
    stop("y should be given as a vector or factor")
  }

  y = as.character(y)

  nclass = length(unique(y))

  if (nclass == 1){
    stop('number of class equals 1')
  } else if(nclass == 2){
    if (length(data) == 1){

      result = binary_univariate_SFLDA(data,y, tau, lambda)
    }
    else{

      result = binary_multivariate_SFLDA(data,y, tau, lambda)
    }
  } else{
    if (length(data) == 1){

      result = multiclass_univariate_SFLDA(data,y, tau, lambda)
    } else{

      result = multiclass_multivariate_SFLDA(data,y, tau, lambda)
    }
  }

  return(result)

}

binary_univariate_SFLDA = function(data, y, tau, lambda){
  data = data[[1]]

  y_class = unique(y)

  xtr0 = data[y == y_class[1],]

  xtr1 = data[y == y_class[2],]

  SD = SDelta_binary(xtr0, xtr1)

  S = SD$S
  mu0 = SD$mu0
  mu1 = SD$mu1

  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')

  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = 'xb1'

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  return(list('type' = 'bu','beta' = b, 'y_class' = y_class, 'md'=md, 'constant' = constant))
}

binary_multivariate_SFLDA = function(data, y, tau, lambda){
  n_mat = length(data)

  p_list = c()

  for (i in seq(n_mat)){
    p_list = c(p_list, ncol(data[[i]]))
  }

  data = do.call(cbind, data)

  y_class = unique(y)

  xtr0 = data[y == y_class[1],]

  xtr1 = data[y == y_class[2],]

  SD = SDelta_binary(xtr0, xtr1)

  S = SD$S
  mu0 = SD$mu0
  mu1 = SD$mu1

  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  s = 0

  for (j in p_list[-n_mat]){
    D[s +j,] = rep(0,p)

    s = s + j
  }

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')

  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)

  beta = vector('list', length = n_mat)

  s = 1

  for (j in seq(n_mat)){

    beta[[j]] = b[seq(s,s + p_list[j]-1),,drop = FALSE]

    s = s + p_list[j]

  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = 'xb1'

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  return(list('type' = 'bm', 'beta' = beta,'y_class' = y_class,  'md'=md, 'constant' = constant))

}

multiclass_univariate_SFLDA = function(data, y, tau, lambda){
  data = data[[1]]

  y_class = unique(y)

  nclass = length(y_class)

  xtr = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    xtr[[n]] = data[y == y_class[n],]
  }

  SD = SDelta_multi(xtr)

  S = SD$S

  delta_list = SD$delta_list

  p = length(delta_list[[1]])
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  S = S / norm(S, type = 'F')

  beta_list =  vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    beta_list[[n]] = coordlasso((1-tau)*S + tau*DD, delta_list[[n]], lambda, 1e-8, 500)
  }

  b = do.call(cbind, beta_list)

  if (pracma::Rank(b) == nclass - 1){
    b = gramSchmidt(b)$Q
  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = paste('xb',as.character(seq(nclass - 1)),sep = '')

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb = is_constant(xb)
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  return(list('type' = 'mu','beta' = b, 'y_class' = y_class, 'md' = md, 'constant' =  constant))

}

multiclass_multivariate_SFLDA = function(data, y, tau, lambda){
  n_mat = length(data)

  p_list = c()

  for (i in seq(n_mat)){
    p_list = c(p_list, ncol(data[[i]]))
  }

  data = do.call(cbind, data)

  y_class = unique(y)

  nclass = length(y_class)

  xtr = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    xtr[[n]] = data[y == y_class[n],]
  }

  SD = SDelta_multi(xtr)

  S = SD$S

  delta_list = SD$delta_list

  p = length(delta_list[[1]])
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  s = 0

  for (j in p_list[-n_mat]){
    D[s +j,] = rep(0,p)

    s = s + j
  }

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  S = S / norm(S, type = 'F')

  beta_list =  vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    beta_list[[n]] = coordlasso((1-tau)*S + tau*DD, delta_list[[n]], lambda, 1e-8, 500)
  }

  b = do.call(cbind, beta_list)

  if (pracma::Rank(b) == nclass - 1){
    b = gramSchmidt(b)$Q
  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = paste('xb',as.character(seq(nclass - 1)),sep = '')

  if (is.null(is_constant(xb))){

    md = table(y)
    constant = TRUE

  } else {
    xb = is_constant(xb)
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  beta = vector('list', length = n_mat)

  s = 1

  for (j in seq(n_mat)){

    beta[[j]] = b[seq(s,s + p_list[j]-1),]

    s = s + p_list[j]

  }

  return(list('type' = 'mm','beta' = beta,'y_class' = y_class,'md'= md, 'constant' = constant))

}

predict_class = function(model,x_test){

  type = model$type

  if (type == 'bu'){

    x_test = x_test[[1]]

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)

    } else{

      beta = model$beta

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = 'xb1'

      pred_class = predict(model$md, pX)$class
    }

  }

  if (type == 'bm'){

    x_test = do.call(cbind, x_test)

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)
    }else {

      beta = model$beta

      beta = do.call(rbind, beta)

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = 'xb1'

      pred_class = predict(model$md, pX)$class

    }

  }

  if (type == 'mu'){

    x_test = x_test[[1]]

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)

    } else{

      beta = model$beta

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = paste('xb',as.character(seq(ncol(pX))),sep = '')

      pred_class = predict(model$md, pX)$class
    }

  }

  if (type == 'mm'){

    x_test = do.call(cbind,x_test)

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)
    }else {

      beta = model$beta

      beta = do.call(rbind, beta)

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = paste('xb',as.character(seq(ncol(pX))),sep = '')

      pred_class = predict(model$md, pX)$class

    }

  }

  return(pred_class)

}

# Cache SFMLDA fits and reuse discriminant profiles across selection thresholds.

q1_format_seconds <- function(x) {
  if (!is.finite(x) || x < 0) return("NA")
  x <- as.integer(round(x))
  h <- x %/% 3600L
  m <- (x %% 3600L) %/% 60L
  s <- x %% 60L
  sprintf("%02d:%02d:%02d", h, m, s)
}

q1_new_progress <- function(total, label = "PROGRESS") {
  e <- new.env(parent = emptyenv())
  e$total <- max(1L, as.integer(total))
  e$label <- as.character(label)
  e$start <- proc.time()[[3]]
  e$last_done <- 0L
  e
}

q1_progress_update <- function(progress, done, extra = NULL, force = FALSE) {
  done <- min(progress$total, max(0L, as.integer(done)))
  if (!force && done <= progress$last_done) return(invisible(NULL))
  elapsed <- proc.time()[[3]] - progress$start
  pct <- 100 * done / progress$total
  eta <- if (done > 0L) elapsed * (progress$total - done) / done else NA_real_
  suffix <- if (!is.null(extra) && nzchar(extra)) paste0(" | ", extra) else ""
  message(sprintf(
    "[%s] %d/%d (%.1f%%) | elapsed %s | ETA %s%s",
    progress$label, done, progress$total, pct,
    q1_format_seconds(elapsed), q1_format_seconds(eta), suffix
  ))
  progress$last_done <- done
  invisible(NULL)
}

q1_default_workers <- function(max_workers = 6L) {
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(cores) || is.na(cores) || cores < 2L) {
    cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  }
  if (!is.finite(cores) || is.na(cores)) cores <- 2L
  max(1L, min(as.integer(max_workers), as.integer(cores) - 1L))
}

q1_configure_parallel <- function(enabled = TRUE,
                                  workers = NULL,
                                  max_workers = 6L,
                                  globals_max_gb = 8) {
  requireNamespace("future")
  requireNamespace("future.apply")

  if (is.null(workers)) workers <- q1_default_workers(max_workers)
  workers <- max(1L, as.integer(workers))
  enabled <- isTRUE(enabled) && workers > 1L

  options(future.globals.maxSize = globals_max_gb * 1024^3)

  if (enabled) {
    future::plan(future::multisession, workers = workers)
    message("[PARALLEL] Windows-safe multisession backend enabled | workers=", workers)
  } else {
    future::plan(future::sequential)
    workers <- 1L
    message("[PARALLEL] Sequential backend enabled | workers=1")
  }

  list(enabled = enabled, workers = workers)
}

q1_shutdown_parallel <- function() {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  invisible(NULL)
}

q1_parallel_lapply <- function(X, FUN, parallel_cfg, ...) {
  if (length(X) == 0L) return(list())
  if (isTRUE(parallel_cfg$enabled) && parallel_cfg$workers > 1L) {
    future.apply::future_lapply(
      X, FUN, ...,
      future.seed = TRUE,
      future.packages = c("MASS", "pracma")
    )
  } else {
    lapply(X, FUN, ...)
  }
}

q1_make_sflda_param_grid <- function() {
  expand.grid(
    tau = c(0.001, 0.01),
    lambda = 10^seq(log10(0.00625), log10(0.3), length.out = 10),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

tune_sflda_intervals_fast <- function(data,
                                      intervals_list,
                                      target_col = "Altezza",
                                      k_folds = 3,
                                      seed = 123,
                                      param_grid = q1_make_sflda_param_grid(),
                                      parallel_cfg = list(enabled = FALSE, workers = 1L),
                                      batch_multiplier = 2L,
                                      label = "SFMLDA",
                                      verbose = TRUE,
                                      return_all = FALSE) {
  requireNamespace("caret")
  requireNamespace("dplyr")

  if (!target_col %in% names(data)) stop("SFMLDA: target column not found: ", target_col)
  if (!is.list(intervals_list) || length(intervals_list) == 0L) {
    stop("SFMLDA: intervals_list must be a non-empty list.")
  }
  if (!all(c("tau", "lambda") %in% names(param_grid))) {
    stop("SFMLDA param_grid must contain columns tau and lambda.")
  }

  y <- factor(as.character(data[[target_col]]), levels = ALTITUDE_LEVELS)
  if (any(table(y) == 0L)) stop("SFMLDA requires all altitude classes in the outer training set.")

  set.seed(seed)
  k_eff <- min(as.integer(k_folds), as.integer(min(table(y))))
  if (k_eff < 2L) stop("SFMLDA: insufficient observations per class for CV.")
  folds <- caret::createFolds(y, k = k_eff, list = TRUE, returnTrain = TRUE)

  spec_cols <- get_spectral_columns(data)
  spec_num <- as.numeric(spec_cols)

  interval_data <- lapply(seq_along(intervals_list), function(i) {
    rng <- intervals_list[[i]]
    cols <- spec_cols[spec_num >= rng[1] & spec_num <= rng[2]]
    if (length(cols) < 3L) return(NULL)
    X <- as.matrix(data[, cols, drop = FALSE])
    storage.mode(X) <- "double"
    list(interval_id = i, range = rng, X = X)
  })
  valid_ids <- which(!vapply(interval_data, is.null, logical(1)))
  if (length(valid_ids) == 0L) stop("SFMLDA: no valid chemical intervals.")

  grid <- as.data.frame(param_grid, stringsAsFactors = FALSE)
  grid$param_id <- seq_len(nrow(grid))
  tasks <- expand.grid(
    interval_id = valid_ids,
    param_id = grid$param_id,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  worker_fun <- function(task_row_index) {
    task <- tasks[task_row_index, , drop = FALSE]
    iid <- as.integer(task$interval_id)
    pid <- as.integer(task$param_id)
    X <- interval_data[[iid]]$X
    tau <- as.numeric(grid$tau[pid])
    lambda <- as.numeric(grid$lambda[pid])

    fold_scores <- vapply(seq_along(folds), function(fi) {
      tr <- folds[[fi]]
      va <- setdiff(seq_len(nrow(X)), tr)
      mod <- tryCatch(
        SFLDA(list(X[tr, , drop = FALSE]), y[tr], tau = tau, lambda = lambda),
        error = function(e) NULL
      )
      if (is.null(mod)) return(NA_real_)
      pred <- tryCatch(
        predict_class(mod, list(X[va, , drop = FALSE])),
        error = function(e) NULL
      )
      if (is.null(pred)) return(NA_real_)
      pred <- factor(as.character(pred), levels = ALTITUDE_LEVELS)
      f1_macro(y[va], pred, levels = ALTITUDE_LEVELS)
    }, numeric(1))

    data.frame(
      interval_id = iid,
      param_id = pid,
      tau = tau,
      lambda = lambda,
      f1_score = if (all(!is.finite(fold_scores))) NA_real_ else mean(fold_scores, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  total <- nrow(tasks)
  progress <- q1_new_progress(total, label)
  batch_size <- if (isTRUE(parallel_cfg$enabled)) {
    max(1L, as.integer(parallel_cfg$workers) * max(1L, as.integer(batch_multiplier)))
  } else {
    1L
  }
  batches <- split(seq_len(total), ceiling(seq_len(total) / batch_size))
  results <- vector("list", length(batches))
  done <- 0L

  if (verbose) {
    message("[", label, "] tuning tau/lambda | intervals=", length(valid_ids),
            " | grid=", nrow(grid), " | CV folds=", length(folds),
            " | total jobs=", total,
            " | workers=", parallel_cfg$workers)
  }

  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    results[[b]] <- q1_parallel_lapply(ids, worker_fun, parallel_cfg)
    done <- done + length(ids)
    if (verbose) q1_progress_update(progress, done)
  }

  all_res <- dplyr::bind_rows(unlist(results, recursive = FALSE))
  all_res$interval_start <- vapply(all_res$interval_id, function(i) interval_data[[i]]$range[1], numeric(1))
  all_res$interval_end <- vapply(all_res$interval_id, function(i) interval_data[[i]]$range[2], numeric(1))

  ok <- all_res[is.finite(all_res$f1_score), , drop = FALSE]
  if (nrow(ok) == 0L) stop("All SFMLDA tau/lambda tuning jobs failed.")

  best <- ok |>
    dplyr::arrange(interval_id, dplyr::desc(f1_score), param_id) |>
    dplyr::group_by(interval_id) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(interval_start, interval_end, tau, lambda, f1_score)

  if (verbose) {
    message("[", label, "] completed. Best tau/lambda by chemical interval:")
    print(best, row.names = FALSE)
  }

  if (return_all) list(best = best, all = all_res) else best
}

fit_sflda_profiles_fixed_params_fast <- function(data,
                                                 param_sflda,
                                                 target_col = "Altezza",
                                                 parallel_cfg = list(enabled = FALSE, workers = 1L),
                                                 label = "SFMLDA-PROFILE",
                                                 verbose = FALSE) {
  if (is.null(param_sflda) || nrow(param_sflda) == 0L) return(list())
  req <- c("interval_start", "interval_end", "tau", "lambda")
  if (!all(req %in% names(param_sflda))) {
    stop("SFMLDA profile parameters missing: ", paste(setdiff(req, names(param_sflda)), collapse = ", "))
  }

  y <- factor(as.character(data[[target_col]]), levels = ALTITUDE_LEVELS)
  spec_cols <- get_spectral_columns(data)
  spec_num <- as.numeric(spec_cols)

  worker_fun <- function(i) {
    row <- param_sflda[i, , drop = FALSE]
    cols <- spec_cols[spec_num >= row$interval_start & spec_num <= row$interval_end]
    if (length(cols) < 3L) return(NULL)
    X <- as.matrix(data[, cols, drop = FALSE])
    storage.mode(X) <- "double"
    mod <- tryCatch(
      SFLDA(list(X), y, tau = as.numeric(row$tau), lambda = as.numeric(row$lambda)),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NULL)
    beta <- if (!is.null(dim(mod$beta))) as.matrix(mod$beta) else matrix(mod$beta, ncol = 1L)
    if (nrow(beta) != length(cols)) return(NULL)
    mag <- sqrt(rowSums(beta^2))
    list(
      interval_start = as.numeric(row$interval_start),
      interval_end = as.numeric(row$interval_end),
      tau = as.numeric(row$tau),
      lambda = as.numeric(row$lambda),
      f1_score = if ("f1_score" %in% names(row)) as.numeric(row$f1_score) else NA_real_,
      freqs = as.numeric(cols),
      beta = beta,
      magnitude = mag
    )
  }

  ids <- seq_len(nrow(param_sflda))
  if (verbose) {
    message("[", label, "] fitting ", length(ids), " cached beta profiles | workers=", parallel_cfg$workers)
  }
  out <- q1_parallel_lapply(ids, worker_fun, parallel_cfg)
  out <- Filter(Negate(is.null), out)
  if (verbose) message("[", label, "] beta profiles ready: ", length(out), "/", length(ids))
  out
}

q1_pad_merge_active_indices <- function(active_idx,
                                        n_rows,
                                        min_width = 5L,
                                        max_gap = 1L) {
  if (length(active_idx) == 0L) return(list())
  breaks <- which(diff(active_idx) > (1L + max_gap))
  starts <- c(active_idx[1L], active_idx[breaks + 1L])
  ends <- c(active_idx[breaks], active_idx[length(active_idx)])

  ranges <- lapply(seq_along(starts), function(i) {
    s <- starts[i]; e <- ends[i]
    width <- e - s + 1L
    if (width < min_width) {
      add <- min_width - width
      s <- s - floor(add / 2)
      e <- e + ceiling(add / 2)
      if (s < 1L) {
        e <- e + (1L - s)
        s <- 1L
      }
      if (e > n_rows) {
        s <- s - (e - n_rows)
        e <- n_rows
      }
      s <- max(1L, s)
      e <- min(n_rows, e)
    }
    c(as.integer(s), as.integer(e))
  })

  ranges <- ranges[order(vapply(ranges, `[`, integer(1), 1L))]
  merged <- list()
  cur <- ranges[[1L]]
  if (length(ranges) > 1L) {
    for (i in 2:length(ranges)) {
      z <- ranges[[i]]
      if (z[1L] <= cur[2L] + 1L + max_gap) {
        cur[2L] <- max(cur[2L], z[2L])
      } else {
        merged[[length(merged) + 1L]] <- cur
        cur <- z
      }
    }
  }
  merged[[length(merged) + 1L]] <- cur
  merged
}

intervals_from_sflda_profiles <- function(profiles,
                                          threshold,
                                          min_width = 5L,
                                          max_gap = 1L) {
  if (is.null(profiles) || length(profiles) == 0L) return(list())
  out <- list()
  for (p in profiles) {
    mag <- p$magnitude
    if (length(mag) == 0L || !any(is.finite(mag))) next
    cutoff <- max(mag, na.rm = TRUE) * as.numeric(threshold)
    active <- which(mag >= cutoff)
    ranges <- q1_pad_merge_active_indices(active, length(mag), min_width, max_gap)
    if (length(ranges) == 0L) next
    for (r in ranges) {
      out[[length(out) + 1L]] <- c(p$freqs[r[1L]], p$freqs[r[2L]])
    }
  }
  out
}

q1_interval_signature <- function(intervals) {
  if (is.null(intervals) || length(intervals) == 0L) return("EMPTY")
  paste(vapply(intervals, function(z) paste0(z[1L], "-", z[2L]), character(1)), collapse = "|")
}

q1_intervals_to_string <- function(intervals) {
  if (is.null(intervals) || length(intervals) == 0L) return(NA_character_)
  paste(vapply(intervals, function(z) paste0("[", z[1L], ",", z[2L], "]"), character(1)), collapse = "; ")
}

# Nested selection of chemically informed and discriminative intervals.

Q1_STATIC_DOMAINS <- c("full", "chemical")
Q1_ALL_DOMAINS <- c("full", "chemical", "discriminative")

Q1_SFLDA_OUTER_PARAM_CACHE <- new.env(parent = emptyenv())

q1_clear_sflda_cache <- function() {
  rm(list = ls(envir = Q1_SFLDA_OUTER_PARAM_CACHE, all.names = TRUE),
     envir = Q1_SFLDA_OUTER_PARAM_CACHE)
  invisible(NULL)
}

q1_sflda_cache_key <- function(data, view, seed, sflda_cv_folds, sflda_intervals, sflda_param_grid) {
  ids <- if ("Matricola" %in% names(data)) as.character(data$Matricola) else seq_len(nrow(data))
  y <- if ("Altezza" %in% names(data)) as.character(data$Altezza) else rep("", nrow(data))
  farms <- if ("Azienda" %in% names(data)) as.character(data$Azienda) else rep("", nrow(data))
  payload <- list(
    view = view, seed = seed, cv = sflda_cv_folds,
    ids = ids, y = y, farms = farms, intervals = sflda_intervals, grid = sflda_param_grid
  )
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(payload, algo = "xxhash64")
  } else {
    paste(view, seed, sflda_cv_folds, nrow(data), paste(ids, collapse = "|"), sep = "::")
  }
}

q1_validate_model <- function(model_type) {
  model_type <- normalise_model_names(model_type)
  if (length(model_type) != 1L) stop("Provide exactly one model_type.")
  if (!model_type %in% SUPPORTED_MODELS) stop("Unsupported model: ", model_type)
  invisible(model_type)
}

q1_make_feature_grid <- function(domain_grid,
                                 derivative_configs,
                                 pve_grid,
                                 threshold_grid) {
  domain_grid <- unique(as.character(domain_grid))
  bad <- setdiff(domain_grid, Q1_ALL_DOMAINS)
  if (length(bad) > 0L) stop("Unknown domain(s): ", paste(bad, collapse = ", "))

  derivative_configs <- unique(as.character(derivative_configs))
  invisible(lapply(derivative_configs, derivative_views_from_config))

  static_domains <- intersect(domain_grid, Q1_STATIC_DOMAINS)
  rows <- list()

  if (length(static_domains) > 0L) {
    z <- expand.grid(
      domain = static_domains,
      derivative_config = derivative_configs,
      pve = pve_grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    z$threshold <- NA_real_
    rows[[length(rows) + 1L]] <- z
  }

  if ("discriminative" %in% domain_grid) {
    z <- expand.grid(
      domain = "discriminative",
      derivative_config = derivative_configs,
      pve = pve_grid,
      threshold = threshold_grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- z
  }

  out <- dplyr::bind_rows(rows)
  out$feature_id <- seq_len(nrow(out))
  out
}

q1_tune_outer_sflda_params <- function(outer_processed_train,
                                       views_needed,
                                       sflda_intervals,
                                       sflda_param_grid,
                                       sflda_cv_folds,
                                       parallel_cfg,
                                       seed,
                                       label_prefix,
                                       target_col = "Altezza",
                                       verbose = TRUE,
                                       use_cache = TRUE) {
  out <- list()
  for (vi in seq_along(views_needed)) {
    v <- views_needed[vi]
    key <- q1_sflda_cache_key(
      data = outer_processed_train[[v]],
      view = v,
      seed = seed + vi * 1000L,
      sflda_cv_folds = sflda_cv_folds,
      sflda_intervals = sflda_intervals,
      sflda_param_grid = sflda_param_grid
    )

    if (isTRUE(use_cache) && exists(key, envir = Q1_SFLDA_OUTER_PARAM_CACHE, inherits = FALSE)) {
      if (verbose) message("[", label_prefix, "] SFMLDA tau/lambda cache HIT | view=", v,
                           " | reused across models")
      out[[v]] <- get(key, envir = Q1_SFLDA_OUTER_PARAM_CACHE, inherits = FALSE)
      next
    }

    if (verbose) {
      message("\n[", label_prefix, "] SFMLDA tau/lambda tuning | view=", v,
              " | OUTER TRAIN ONLY")
    }
    ans <- tune_sflda_intervals_fast(
      data = outer_processed_train[[v]],
      intervals_list = sflda_intervals,
      target_col = target_col,
      k_folds = sflda_cv_folds,
      seed = seed + vi * 1000L,
      param_grid = sflda_param_grid,
      parallel_cfg = parallel_cfg,
      label = paste0(label_prefix, ":", v),
      verbose = verbose,
      return_all = FALSE
    )
    out[[v]] <- ans
    if (isTRUE(use_cache)) assign(key, ans, envir = Q1_SFLDA_OUTER_PARAM_CACHE)
  }
  out
}

q1_build_static_inner_cache <- function(inner_preprocess_cache,
                                        static_domains,
                                        views_needed,
                                        max_pve,
                                        target_col = "Altezza",
                                        id_col = "Matricola",
                                        group_col = "Azienda",
                                        verbose = TRUE,
                                        label = "STATIC-FPCA") {
  if (length(static_domains) == 0L) return(vector("list", length(inner_preprocess_cache)))
  out <- vector("list", length(inner_preprocess_cache))
  progress <- q1_new_progress(length(inner_preprocess_cache), label)

  for (i in seq_along(inner_preprocess_cache)) {
    pp <- inner_preprocess_cache[[i]]$processed
    out[[i]] <- build_preprocessed_fpca_cache(
      train_views = pp$train,
      test_views = pp$test,
      domains = static_domains,
      views_needed = views_needed,
      max_pve = max_pve,
      id_col = id_col,
      group_col = group_col,
      target_col = target_col
    )
    if (!is.null(inner_preprocess_cache[[i]]$validation_farm)) {
      out[[i]]$validation_farm <- inner_preprocess_cache[[i]]$validation_farm
      out[[i]]$validation_target <- inner_preprocess_cache[[i]]$validation_target
    }
    if (verbose) q1_progress_update(progress, i)
  }
  out
}

q1_build_discriminative_split_cache <- function(processed_split,
                                                sflda_params_by_view,
                                                thresholds,
                                                views_needed,
                                                max_pve,
                                                parallel_cfg,
                                                target_col = "Altezza",
                                                id_col = "Matricola",
                                                group_col = "Azienda",
                                                profile_label = "DISC-PROFILES",
                                                verbose = FALSE) {
  train_views <- processed_split$train
  test_views <- processed_split$test

  meta_cols <- intersect(c(id_col, group_col, target_col), names(train_views[[views_needed[1L]]]))
  train_meta <- train_views[[views_needed[1L]]][, meta_cols, drop = FALSE]
  test_meta <- test_views[[views_needed[1L]]][, meta_cols, drop = FALSE]

  profiles_by_view <- list()
  for (v in views_needed) {
    params <- sflda_params_by_view[[v]]
    if (is.null(params) || nrow(params) == 0L) {
      profiles_by_view[[v]] <- list()
      next
    }
    profiles_by_view[[v]] <- fit_sflda_profiles_fixed_params_fast(
      data = train_views[[v]],
      param_sflda = params,
      target_col = target_col,
      parallel_cfg = parallel_cfg,
      label = paste0(profile_label, ":", v),
      verbose = verbose
    )
  }

  view_threshold <- list()
  for (v in views_needed) {
    profiles <- profiles_by_view[[v]]
    by_threshold <- vector("list", length(thresholds))
    names(by_threshold) <- as.character(thresholds)

    interval_sets <- lapply(thresholds, function(th) {
      intervals_from_sflda_profiles(profiles, th, min_width = 5L, max_gap = 1L)
    })
    signatures <- vapply(interval_sets, q1_interval_signature, character(1))
    unique_sigs <- unique(signatures)
    fpca_by_sig <- list()

    for (sig in unique_sigs) {
      ids <- which(signatures == sig)
      ints <- interval_sets[[ids[1L]]]
      if (identical(sig, "EMPTY") || length(ints) == 0L) {
        fpca_by_sig[[sig]] <- NULL
      } else {
        fpca_by_sig[[sig]] <- fit_fpca_view_pair_cache(
          train_views[[v]],
          test_views[[v]],
          intervals = ints,
          max_pve = max_pve
        )
      }
    }

    for (ti in seq_along(thresholds)) {
      by_threshold[[ti]] <- list(
        intervals = interval_sets[[ti]],
        bundles = fpca_by_sig[[signatures[ti]]],
        signature = signatures[ti]
      )
    }
    view_threshold[[v]] <- by_threshold
  }

  out <- vector("list", length(thresholds))
  names(out) <- as.character(thresholds)

  for (ti in seq_along(thresholds)) {
    cache <- list(
      train_meta = train_meta,
      test_meta = test_meta,
      domains = list(discriminative = list()),
      discriminative_intervals = list(),
      threshold = thresholds[ti]
    )
    for (v in views_needed) {
      cache$domains$discriminative[[v]] <- view_threshold[[v]][[ti]]$bundles
      cache$discriminative_intervals[[v]] <- view_threshold[[v]][[ti]]$intervals
    }
    out[[ti]] <- cache
  }

  out
}

q1_build_discriminative_inner_cache <- function(inner_preprocess_cache,
                                                sflda_params_by_view,
                                                thresholds,
                                                views_needed,
                                                max_pve,
                                                parallel_cfg,
                                                target_col = "Altezza",
                                                id_col = "Matricola",
                                                group_col = "Azienda",
                                                verbose = TRUE,
                                                label = "DISC-CACHE") {
  out <- vector("list", length(inner_preprocess_cache))
  progress <- q1_new_progress(length(inner_preprocess_cache), label)

  for (i in seq_along(inner_preprocess_cache)) {
    if (verbose) {
      extra <- if (!is.null(inner_preprocess_cache[[i]]$validation_farm)) {
        paste0("validation farm=", inner_preprocess_cache[[i]]$validation_farm)
      } else {
        paste0("inner split=", i)
      }
      message("[", label, "] fitting cached SFMLDA profiles + threshold FPCA | ", extra)
    }
    out[[i]] <- q1_build_discriminative_split_cache(
      processed_split = inner_preprocess_cache[[i]]$processed,
      sflda_params_by_view = sflda_params_by_view,
      thresholds = thresholds,
      views_needed = views_needed,
      max_pve = max_pve,
      parallel_cfg = parallel_cfg,
      target_col = target_col,
      id_col = id_col,
      group_col = group_col,
      profile_label = paste0(label, ":split", i),
      verbose = FALSE
    )
    if (!is.null(inner_preprocess_cache[[i]]$validation_farm)) {
      attr(out[[i]], "validation_farm") <- inner_preprocess_cache[[i]]$validation_farm
      attr(out[[i]], "validation_target") <- inner_preprocess_cache[[i]]$validation_target
    }
    if (verbose) q1_progress_update(progress, i)
  }
  out
}

q1_get_cache_for_feature <- function(static_cache,
                                     discriminative_cache,
                                     split_index,
                                     domain,
                                     threshold) {
  if (domain %in% Q1_STATIC_DOMAINS) return(static_cache[[split_index]])
  if (!identical(domain, "discriminative")) stop("Unknown domain: ", domain)
  if (!is.finite(threshold)) return(NULL)
  discriminative_cache[[split_index]][[as.character(threshold)]]
}

q1_build_feature_bundles <- function(static_cache,
                                     discriminative_cache,
                                     feature_row,
                                     target_col = "Altezza") {
  n_splits <- max(length(static_cache), length(discriminative_cache))
  out <- vector("list", n_splits)

  for (j in seq_len(n_splits)) {
    cache <- q1_get_cache_for_feature(
      static_cache = static_cache,
      discriminative_cache = discriminative_cache,
      split_index = j,
      domain = as.character(feature_row$domain),
      threshold = as.numeric(feature_row$threshold)
    )
    if (is.null(cache)) return(NULL)

    feat <- tryCatch(
      assemble_fpca_features(
        cache = cache,
        domain = as.character(feature_row$domain),
        derivative_config = as.character(feature_row$derivative_config),
        pve = as.numeric(feature_row$pve),
        target_col = target_col,
        scale_minmax = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(feat)) return(NULL)
    out[[j]] <- feat
  }
  out
}

q1_evaluate_model_sample_bundles <- function(feature_bundles,
                                             model_type,
                                             param_row,
                                             seed) {
  acc <- numeric(length(feature_bundles))
  f1 <- numeric(length(feature_bundles))
  nf <- numeric(length(feature_bundles))

  for (j in seq_along(feature_bundles)) {
    feat <- feature_bundles[[j]]
    fit <- tryCatch(
      fit_model_candidate(feat$x_train, feat$y_train, model_type, param_row, seed = seed + j),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    pred <- tryCatch(predict_model_candidate(fit, feat$x_test), error = function(e) NULL)
    if (is.null(pred)) return(NULL)
    truth <- factor(as.character(feat$y_test), levels = ALTITUDE_LEVELS)
    acc[j] <- mean(as.character(pred) == as.character(truth))
    f1[j] <- f1_macro(truth, pred)
    nf[j] <- feat$n_features
  }

  data.frame(
    mean_inner_sample_accuracy = mean(acc),
    sd_inner_sample_accuracy = stats::sd(acc),
    mean_inner_macro_f1 = mean(f1),
    mean_n_features = mean(nf),
    stringsAsFactors = FALSE
  )
}

q1_tune_joint_inner_model <- function(static_cache,
                                    discriminative_cache,
                                    feature_grid,
                                    model_type,
                                    model_grid,
                                    design = "sample",
                                    validation_targets = NULL,
                                    target_col = "Altezza",
                                    seed = 123,
                                    parallel_cfg = list(enabled = FALSE, workers = 1L),
                                    verbose = TRUE,
                                    label = "INNER-TUNING") {
  if (!identical(design, "sample")) stop("Only sample-level CV is supported in this release.")
  model_type <- q1_validate_model(model_type)

  total_jobs <- nrow(feature_grid) * nrow(model_grid)
  progress <- q1_new_progress(total_jobs, label)
  rows <- list()
  done <- 0L

  for (fi in seq_len(nrow(feature_grid))) {
    fg <- feature_grid[fi, , drop = FALSE]

    if (verbose) {
      th_text <- if (is.finite(fg$threshold)) paste0(" | threshold=", fg$threshold) else ""
      message("\n[", label, "] feature configuration ", fi, "/", nrow(feature_grid),
              " | domain=", fg$domain,
              " | deriv=", fg$derivative_config,
              " | PVE=", fg$pve, th_text)
    }

    bundles <- q1_build_feature_bundles(
      static_cache = static_cache,
      discriminative_cache = discriminative_cache,
      feature_row = fg,
      target_col = target_col
    )

    if (is.null(bundles)) {
      done <- done + nrow(model_grid)
      if (verbose) q1_progress_update(progress, done, extra = "feature config unavailable")
      next
    }

    param_worker <- function(mi) {
      pr <- model_grid[mi, , drop = FALSE]
      met <- q1_evaluate_model_sample_bundles(
        bundles, model_type, pr, seed = seed + fi * 100000L + mi * 100L
      )
      if (is.null(met)) return(NULL)
      cbind(
        data.frame(
          feature_id = fg$feature_id,
          domain = as.character(fg$domain),
          derivative_config = as.character(fg$derivative_config),
          pve = as.numeric(fg$pve),
          threshold = as.numeric(fg$threshold),
          model = model_type,
          param_id = pr$param_id,
          param_label = model_param_label(model_type, pr),
          stringsAsFactors = FALSE
        ),
        met
      )
    }

    res <- q1_parallel_lapply(seq_len(nrow(model_grid)), param_worker, parallel_cfg)
    res <- Filter(Negate(is.null), res)
    if (length(res) > 0L) rows[[length(rows) + 1L]] <- dplyr::bind_rows(res)

    done <- done + nrow(model_grid)
    if (verbose) q1_progress_update(progress, done)
  }

  tuning <- dplyr::bind_rows(rows)
  if (nrow(tuning) == 0L) stop("All joint inner-CV configurations failed.")

  if (design == "sample") {
    tuning <- tuning |>
      dplyr::arrange(
        dplyr::desc(mean_inner_sample_accuracy),
        sd_inner_sample_accuracy,
        mean_n_features,
        pve,
        domain,
        derivative_config,
        threshold,
        param_id
      )
  } else {
    tuning <- tuning |>
      dplyr::arrange(
        dplyr::desc(class_balanced_inner_farm_accuracy),
        sd_inner_farm_accuracy,
        mean_n_features,
        pve,
        domain,
        derivative_config,
        threshold,
        param_id
      )
  }

  list(best = tuning[1L, , drop = FALSE], tuning = tuning)
}

q1_build_final_outer_feature <- function(outer_processed,
                                         best,
                                         sflda_params_by_view,
                                         views_needed_all,
                                         max_pve,
                                         parallel_cfg,
                                         target_col = "Altezza",
                                         id_col = "Matricola",
                                         group_col = "Azienda") {
  domain <- as.character(best$domain)
  deriv <- as.character(best$derivative_config)
  pve <- as.numeric(best$pve)
  selected_views <- derivative_views_from_config(deriv)

  if (domain %in% Q1_STATIC_DOMAINS) {
    cache <- build_preprocessed_fpca_cache(
      train_views = outer_processed$train,
      test_views = outer_processed$test,
      domains = domain,
      views_needed = selected_views,
      max_pve = pve,
      id_col = id_col,
      group_col = group_col,
      target_col = target_col
    )
    feat <- assemble_fpca_features(
      cache, domain, deriv, pve,
      target_col = target_col,
      scale_minmax = TRUE
    )
    return(list(feat = feat, intervals = list()))
  }

  th <- as.numeric(best$threshold)
  dcache <- q1_build_discriminative_split_cache(
    processed_split = outer_processed,
    sflda_params_by_view = sflda_params_by_view,
    thresholds = th,
    views_needed = selected_views,
    max_pve = pve,
    parallel_cfg = parallel_cfg,
    target_col = target_col,
    id_col = id_col,
    group_col = group_col,
    profile_label = "FINAL-DISC-PROFILE",
    verbose = FALSE
  )[[1L]]

  feat <- assemble_fpca_features(
    dcache, "discriminative", deriv, pve,
    target_col = target_col,
    scale_minmax = TRUE
  )
  list(feat = feat, intervals = dcache$discriminative_intervals)
}

q1_joint_sample_summary <- function(fold_metrics,
                                    predictions,
                                    model_type,
                                    conf_level = 0.95) {
  ci_acc <- mean_t_ci(fold_metrics$accuracy, conf_level)
  ci_f1 <- mean_t_ci(fold_metrics$macro_f1, conf_level)
  class_unc <- sample_fold_class_uncertainty(predictions, conf_level, ALTITUDE_LEVELS)

  list(
    overall = data.frame(
      model = model_type,
      mean_accuracy = unname(ci_acc["mean"]),
      accuracy_ci_lower = unname(ci_acc["lower"]),
      accuracy_ci_upper = unname(ci_acc["upper"]),
      sd_accuracy = stats::sd(fold_metrics$accuracy),
      mean_macro_f1 = unname(ci_f1["mean"]),
      macro_f1_ci_lower = unname(ci_f1["lower"]),
      macro_f1_ci_upper = unname(ci_f1["upper"]),
      sd_macro_f1 = stats::sd(fold_metrics$macro_f1),
      pooled_sample_accuracy = mean(predictions$truth == predictions$pred),
      n_outer_folds = nrow(fold_metrics),
      ci_level = conf_level,
      ci_note = "Descriptive t CI across sample-level outer folds; folds are not independent farms.",
      stringsAsFactors = FALSE
    ),
    pooled_confusion_counts = table(
      True = factor(predictions$truth, levels = ALTITUDE_LEVELS),
      Predicted = factor(predictions$pred, levels = ALTITUDE_LEVELS)
    ),
    class_specific_metrics = class_unc$summary,
    class_specific_ci_note = class_unc$ci_note,
    domain_selection_counts = table(fold_metrics$selected_domain),
    derivative_selection_counts = table(fold_metrics$selected_derivative),
    pve_selection_counts = table(fold_metrics$selected_pve),
    threshold_selection_counts = table(fold_metrics$selected_threshold, useNA = "ifany"),
    model_parameter_selection_counts = table(fold_metrics$selected_param_label)
  )
}

q1_consensus_interval_ranges <- function(interval_lists,
                                         eligible,
                                         min_frequency = 0.5,
                                         grid = 240:1299) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) return(data.frame())
  counts <- integer(length(grid))

  for (i in eligible_ids) {
    ints <- interval_lists[[i]]
    if (is.null(ints) || length(ints) == 0L) next
    mark <- logical(length(grid))
    for (z in ints) mark[grid >= z[1L] & grid <= z[2L]] <- TRUE
    counts <- counts + as.integer(mark)
  }

  freq <- counts / length(eligible_ids)
  active <- which(freq >= min_frequency)
  if (length(active) == 0L) return(data.frame())
  ranges <- q1_pad_merge_active_indices(active, length(grid), min_width = 1L, max_gap = 0L)
  dplyr::bind_rows(lapply(ranges, function(r) {
    ids <- r[1L]:r[2L]
    data.frame(
      start_index = grid[r[1L]],
      end_index = grid[r[2L]],
      min_selection_frequency = min(freq[ids]),
      max_selection_frequency = max(freq[ids]),
      n_discriminative_outer_folds = length(eligible_ids),
      stringsAsFactors = FALSE
    )
  }))
}

q1_print_interval_report <- function(metrics,
                                     intervals_m0,
                                     intervals_m1,
                                     intervals_m2,
                                     title) {
  cat("\n", title, "\n", sep = "")
  disc <- metrics$selected_domain == "discriminative"
  cat("Discriminative domain selected in ", sum(disc), "/", nrow(metrics), " outer folds.\n", sep = "")
  if (!any(disc)) return(invisible(NULL))

  for (i in which(disc)) {
    cat("  fold ", metrics$fold[i],
        " | derivative=", metrics$selected_derivative[i],
        " | threshold=", metrics$selected_threshold[i],
        " | m0: ", q1_intervals_to_string(intervals_m0[[i]]),
        " | m1: ", q1_intervals_to_string(intervals_m1[[i]]),
        " | m2: ", q1_intervals_to_string(intervals_m2[[i]]), "\n", sep = "")
  }

  cat("\nConsensus m0 ranges selected in >=50% of discriminative outer folds:\n")
  c0 <- q1_consensus_interval_ranges(intervals_m0, disc, 0.5)
  if (nrow(c0) == 0L) cat("  none\n") else print(c0, row.names = FALSE)

  m1_eligible <- disc & metrics$selected_derivative %in% c("m0_d1", "m0_d1_d2")
  cat("\nConsensus m1 ranges selected in >=50% of discriminative outer folds that used d1:\n")
  c1 <- q1_consensus_interval_ranges(intervals_m1, m1_eligible, 0.5)
  if (nrow(c1) == 0L) cat("  none / m1 not selected\n") else print(c1, row.names = FALSE)

  m2_eligible <- disc & metrics$selected_derivative %in% c("m0_d2", "m0_d1_d2")
  cat("\nConsensus m2 ranges selected in >=50% of discriminative outer folds that used d2:\n")
  c2 <- q1_consensus_interval_ranges(intervals_m2, m2_eligible, 0.5)
  if (nrow(c2) == 0L) cat("  none / m2 not selected\n") else print(c2, row.names = FALSE)

  invisible(NULL)
}

# Nested altitude classification with joint domain selection.
run_paper_sample_nested_model_discriminative <- function(data,
                                                       model_type,
                                                       model_grid = default_model_grid(model_type),
                                                       domain_grid = c("full", "chemical", "discriminative"),
                                                       derivative_configs = c("m0", "m0_d1", "m0_d2", "m0_d1_d2"),
                                                       pve_grid = c(0.95, 0.98, 0.99, 0.9999),
                                                       threshold_grid = c(0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
                                                       sflda_intervals = ftir_domain_intervals("chemical"),
                                                       sflda_param_grid = q1_make_sflda_param_grid(),
                                                       sflda_cv_folds = 3,
                                                       outer_k_folds = 10,
                                                       inner_k_folds = 5,
                                                       col_start = "240",
                                                       col_end = "1299",
                                                       polynomial_degree = 4,
                                                       fixed_sg_window = NULL,
                                                       conf_level = 0.95,
                                                       seed = 69,
                                                       target_col = "Altezza",
                                                       id_col = "Matricola",
                                                       group_col = "Azienda",
                                                       parallel_cfg = list(enabled = FALSE, workers = 1L),
                                                       verbose = TRUE) {
  model_type <- q1_validate_model(model_type)
  if (anyDuplicated(data[[id_col]]) > 0L) stop("Sample IDs must be unique.")

  outer_folds <- make_stratified_sample_folds(data[[target_col]], outer_k_folds, seed)
  static_domains <- intersect(domain_grid, Q1_STATIC_DOMAINS)
  views_needed <- unique(unlist(lapply(derivative_configs, derivative_views_from_config)))
  max_pve <- max(pve_grid)
  feature_grid <- q1_make_feature_grid(domain_grid, derivative_configs, pve_grid, threshold_grid)

  fold_rows <- vector("list", length(outer_folds))
  pred_rows <- vector("list", length(outer_folds))
  tuning_logs <- vector("list", length(outer_folds))
  sflda_logs <- vector("list", length(outer_folds))
  intervals_m0 <- vector("list", length(outer_folds))
  intervals_m1 <- vector("list", length(outer_folds))
  intervals_m2 <- vector("list", length(outer_folds))
  outer_progress <- q1_new_progress(length(outer_folds), "PAPER-OUTER")
  n <- nrow(data)

  for (f in seq_along(outer_folds)) {
    tr <- outer_folds[[f]]
    te <- setdiff(seq_len(n), tr)
    raw_train <- data[tr, , drop = FALSE]
    raw_test <- data[te, , drop = FALSE]
    fold_seed <- seed + f * 10000L

    if (verbose) {
      message("\n============================================================")
      message("PAPER SAMPLE OUTER FOLD ", f, "/", length(outer_folds),
              " | ", toupper(model_type), " | train n=", nrow(raw_train), " | test n=", nrow(raw_test))
      q1_progress_update(outer_progress, f - 1L, force = TRUE)
    }

    outer_pp <- preprocess_grouped_split(
      raw_train, raw_test,
      col_start = col_start,
      col_end = col_end,
      polynomial_degree = polynomial_degree,
      fixed_window = fixed_sg_window
    )

    sflda_params <- list()
    if ("discriminative" %in% domain_grid) {
      sflda_params <- q1_tune_outer_sflda_params(
        outer_processed_train = outer_pp$train,
        views_needed = views_needed,
        sflda_intervals = sflda_intervals,
        sflda_param_grid = sflda_param_grid,
        sflda_cv_folds = sflda_cv_folds,
        parallel_cfg = parallel_cfg,
        seed = fold_seed + 100L,
        label_prefix = paste0("PAPER-F", f, "-SFMLDA"),
        target_col = target_col,
        verbose = verbose
      )
    }
    sflda_logs[[f]] <- sflda_params

    inner_folds <- make_stratified_sample_folds(raw_train[[target_col]], inner_k_folds, fold_seed + 1L)
    inner_pp <- build_inner_sample_preprocess_cache(
      raw_outer_train = raw_train,
      inner_train_folds = inner_folds,
      col_start = col_start,
      col_end = col_end,
      polynomial_degree = polynomial_degree,
      fixed_sg_window = fixed_sg_window
    )

    static_cache <- q1_build_static_inner_cache(
      inner_preprocess_cache = inner_pp,
      static_domains = static_domains,
      views_needed = views_needed,
      max_pve = max_pve,
      target_col = target_col,
      id_col = id_col,
      group_col = group_col,
      verbose = verbose,
      label = paste0("PAPER-F", f, "-STATIC-FPCA")
    )

    discr_cache <- vector("list", length(inner_pp))
    if ("discriminative" %in% domain_grid) {
      discr_cache <- q1_build_discriminative_inner_cache(
        inner_preprocess_cache = inner_pp,
        sflda_params_by_view = sflda_params,
        thresholds = threshold_grid,
        views_needed = views_needed,
        max_pve = max_pve,
        parallel_cfg = parallel_cfg,
        target_col = target_col,
        id_col = id_col,
        group_col = group_col,
        verbose = verbose,
        label = paste0("PAPER-F", f, "-DISC-CACHE")
      )
    }

    tuned <- q1_tune_joint_inner_model(
      static_cache = static_cache,
      discriminative_cache = discr_cache,
      feature_grid = feature_grid,
      model_type = model_type,
      model_grid = model_grid,
      design = "sample",
      target_col = target_col,
      seed = fold_seed + 500L,
      parallel_cfg = parallel_cfg,
      verbose = verbose,
      label = paste0("PAPER-F", f, "-INNER")
    )
    best <- tuned$best
    tuning_logs[[f]] <- tuned$tuning

    if (verbose) {
      message("\n[PAPER-F", f, "] SELECTED | domain=", best$domain,
              " | derivative=", best$derivative_config,
              " | PVE=", best$pve,
              if (is.finite(best$threshold)) paste0(" | threshold=", best$threshold) else "",
              " | ", best$param_label,
              " | inner acc=", round(best$mean_inner_sample_accuracy, 4))
    }

    final <- q1_build_final_outer_feature(
      outer_processed = outer_pp,
      best = best,
      sflda_params_by_view = sflda_params,
      views_needed_all = views_needed,
      max_pve = max_pve,
      parallel_cfg = parallel_cfg,
      target_col = target_col,
      id_col = id_col,
      group_col = group_col
    )
    feat <- final$feat
    if (is.null(feat)) stop("Final outer feature construction failed in sample fold ", f)

    pr <- model_grid[model_grid$param_id == best$param_id, , drop = FALSE]
    fit <- fit_model_candidate(feat$x_train, feat$y_train, model_type, pr, seed = fold_seed + 999L)
    pred <- predict_model_candidate(fit, feat$x_test)
    truth <- factor(as.character(feat$y_test), levels = ALTITUDE_LEVELS)
    acc <- mean(as.character(pred) == as.character(truth))
    mf1 <- f1_macro(truth, pred)

    intervals_m0[[f]] <- final$intervals$m0 %||% list()
    intervals_m1[[f]] <- final$intervals$m1 %||% list()
    intervals_m2[[f]] <- final$intervals$m2 %||% list()

    fold_rows[[f]] <- data.frame(
      fold = f,
      model = model_type,
      accuracy = acc,
      macro_f1 = mf1,
      n_test = length(truth),
      selected_domain = as.character(best$domain),
      selected_derivative = as.character(best$derivative_config),
      selected_pve = as.numeric(best$pve),
      selected_threshold = as.numeric(best$threshold),
      selected_param_id = best$param_id,
      selected_param_label = best$param_label,
      inner_mean_sample_accuracy = best$mean_inner_sample_accuracy,
      inner_mean_macro_f1 = best$mean_inner_macro_f1,
      n_features = feat$n_features,
      sg_window = outer_pp$window_length,
      stringsAsFactors = FALSE
    )

    pred_rows[[f]] <- data.frame(
      fold = f,
      sample_id = as.character(feat$test_meta[[id_col]]),
      farm = as.character(feat$test_meta[[group_col]]),
      truth = as.character(truth),
      pred = as.character(pred),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message("[PAPER-F", f, "] OUTER ACC=", round(acc, 4), " | MACRO-F1=", round(mf1, 4))
      q1_progress_update(outer_progress, f, force = TRUE)
    }
  }

  fm <- dplyr::bind_rows(fold_rows)
  pp <- dplyr::bind_rows(pred_rows)
  list(
    model = model_type,
    fold_metrics = fm,
    predictions = pp,
    summary = q1_joint_sample_summary(fm, pp, model_type, conf_level),
    tuning_logs = tuning_logs,
    sflda_params_by_outer_fold = sflda_logs,
    selected_intervals_m0 = intervals_m0,
    selected_intervals_m1 = intervals_m1,
    selected_intervals_m2 = intervals_m2,
    feature_grid = feature_grid
  )
}

q1_collect_sflda_outer_params <- function(sflda_logs) {
  rows <- list()
  for (f in seq_along(sflda_logs)) {
    z <- sflda_logs[[f]]
    if (is.null(z) || length(z) == 0L) next
    for (v in names(z)) {
      d <- z[[v]]
      if (is.null(d) || nrow(d) == 0L) next
      d$fold <- f
      d$view <- v
      rows[[length(rows) + 1L]] <- d
    }
  }
  dplyr::bind_rows(rows)
}

q1_print_sflda_outer_summary <- function(sflda_logs, title = "SFMLDA OUTER-TRAIN PARAMETER SUMMARY") {
  d <- q1_collect_sflda_outer_params(sflda_logs)
  cat("\n", title, ":\n", sep = "")
  if (nrow(d) == 0L) {
    cat("  no SFMLDA tuning results\n")
    return(invisible(NULL))
  }

  sm <- d |>
    dplyr::group_by(view, interval_start, interval_end) |>
    dplyr::summarise(
      n_outer_folds = dplyr::n(),
      mean_tau = mean(tau, na.rm = TRUE),
      min_tau = min(tau, na.rm = TRUE),
      max_tau = max(tau, na.rm = TRUE),
      mean_lambda = mean(lambda, na.rm = TRUE),
      min_lambda = min(lambda, na.rm = TRUE),
      max_lambda = max(lambda, na.rm = TRUE),
      mean_internal_f1 = mean(f1_score, na.rm = TRUE),
      .groups = "drop"
    )
  print(sm, n = Inf)
  invisible(sm)
}

compare_joint_paper_models <- function(model_results,
                                       conf_level = 0.95,
                                       seed = 123) {
  models <- names(model_results)
  if (length(models) < 2L) return(data.frame())
  pairs <- combn(models, 2, simplify = FALSE)
  out <- list()

  for (pp in pairs) {
    a_name <- pp[1L]
    b_name <- pp[2L]
    a_res <- model_results[[a_name]]
    b_res <- model_results[[b_name]]

    key_a <- paste(a_res$predictions$fold, a_res$predictions$sample_id, sep = "::")
    key_b <- paste(b_res$predictions$fold, b_res$predictions$sample_id, sep = "::")
    if (!setequal(key_a, key_b)) {
      stop("Sample outer partitions differ between ", a_name, " and ", b_name, ".")
    }

    m <- merge(
      a_res$fold_metrics[, c("fold", "accuracy", "macro_f1")],
      b_res$fold_metrics[, c("fold", "accuracy", "macro_f1")],
      by = "fold", suffixes = c("_a", "_b"), all = FALSE
    )

    for (metric in c("accuracy", "macro_f1")) {
      rr <- paired_difference_summary(
        m[[paste0(metric, "_a")]],
        m[[paste0(metric, "_b")]],
        conf_level = conf_level,
        seed = seed + length(out) + 1L
      )
      rr$design <- "paper_sample_cv_joint_domain_selection"
      rr$metric <- metric
      rr$model_a <- a_name
      rr$model_b <- b_name
      rr$interpretation <- "positive difference = model_a higher"
      rr$inference_note <- paste0(
        "Paired by identical outer sample fold; each model independently selects domain/derivative/PVE/",
        "discriminative threshold and model hyperparameters in inner CV. CI and sign-flip p-value are descriptive."
      )
      out[[length(out) + 1L]] <- rr
    }
  }
  dplyr::bind_rows(out)
}

print_paper_model_discriminative_results <- function(res, breed_name, model_type = res$model) {
  fm <- res$fold_metrics
  sm <- res$summary

  cat("\n\n================ PAPER SAMPLE-CV | ", breed_name,
      " | ", toupper(model_type), " | JOINT DOMAIN SELECTION ================\n", sep = "")
  cat("No outlier removal. SG/FPCA/min-max are training-only.\n")
  cat("SFMLDA tau/lambda are tuned once on each OUTER training set; beta profiles/intervals are refitted on INNER training sets.\n")

  cat("\nOuter sample-fold metrics:\n")
  print(fm, row.names = FALSE)

  cat("\nSummary + descriptive 95% CI across outer sample folds:\n")
  print(sm$overall, row.names = FALSE)

  cat("\nPooled confusion counts:\n")
  print(sm$pooled_confusion_counts)

  cat("\nClass-specific Precision / Recall / F1 + uncertainty:\n")
  print(sm$class_specific_metrics, row.names = FALSE)
  cat("NOTE: ", sm$class_specific_ci_note, "\n", sep = "")

  cat("\nFPC feature-count summary across outer folds:\n")
  print(selection_numeric_summary(fm$n_features, "n_FPC_features"), row.names = FALSE)

  print_selection_frequency(fm$selected_domain, "Selected domain frequency")
  print_selection_frequency(fm$selected_derivative, "Selected derivative frequency")
  print_selection_frequency(fm$selected_pve, "Selected PVE frequency")
  print_selection_frequency(fm$selected_threshold, "Selected discriminative threshold frequency (NA = non-discriminative domain)")
  print_selection_frequency(fm$sg_window, "Training-derived SG-window frequency")
  print_selection_frequency(fm$selected_param_label, paste0("Selected ", toupper(model_type), "-parameter frequency"))

  q1_print_sflda_outer_summary(res$sflda_params_by_outer_fold)

  q1_print_interval_report(
    metrics = fm,
    intervals_m0 = res$selected_intervals_m0,
    intervals_m1 = res$selected_intervals_m1,
    intervals_m2 = res$selected_intervals_m2,
    title = "DISCRIMINATIVE INTERVALS BY OUTER SAMPLE FOLD"
  )
  invisible(res)
}

q1_build_inner_grouped_preprocess_cache <- function(raw_outer_train,
                                                     grouped_plan,
                                                     group_col = "Azienda",
                                                     col_start = "240",
                                                     col_end = "1299",
                                                     polynomial_degree = 4,
                                                     fixed_sg_window = NULL) {
  folds <- sort(unique(grouped_plan$inner_fold))
  out <- vector("list", length(folds))

  for (ii in seq_along(folds)) {
    f <- folds[ii]
    val_farms <- grouped_plan$group[grouped_plan$inner_fold == f]
    is_val <- as.character(raw_outer_train[[group_col]]) %in% val_farms
    tr <- which(!is_val)
    va <- which(is_val)
    if (length(tr) == 0L || length(va) == 0L) stop("Empty train/validation split in grouped inner fold ", f)

    out[[ii]] <- list(
      validation_farm = paste(val_farms, collapse = "+"),
      validation_target = "grouped",
      validation_farms = val_farms,
      inner_fold = f,
      processed = preprocess_grouped_split(
        raw_outer_train[tr, , drop = FALSE],
        raw_outer_train[va, , drop = FALSE],
        col_start, col_end, polynomial_degree, fixed_sg_window
      )
    )
  }
  out
}

get_qc_spectral_columns <- function(data, col_start = 240, col_end = 1299) {
  nm <- names(data)
  is_spec <- grepl("^[0-9]+(\\.[0-9]+)?$", nm)
  vals <- suppressWarnings(as.numeric(nm))
  nm[is_spec & !is.na(vals) & vals >= col_start & vals <= col_end]
}

fit_train_global_functional_qc <- function(train_smoothed,
                                           f_value = 1,
                                           central_prop = 0.50,
                                           col_start = 240,
                                           col_end = 1299) {
  if (!requireNamespace("roahd", quietly = TRUE)) stop("Package 'roahd' is required.")
  spec_cols <- get_qc_spectral_columns(train_smoothed, col_start, col_end)
  if (length(spec_cols) < 5) stop("Too few spectral columns for functional QC.")

  X <- as.matrix(train_smoothed[, spec_cols, drop = FALSE])
  storage.mode(X) <- "double"
  grid <- as.numeric(spec_cols)

  fd <- roahd::fData(grid, X)
  depths <- roahd::MBD(fd)
  ord <- order(depths, decreasing = TRUE)
  n_central <- max(2L, ceiling(nrow(X) * central_prop))
  central <- X[ord[seq_len(n_central)], , drop = FALSE]

  bag_min <- apply(central, 2, min, na.rm = TRUE)
  bag_max <- apply(central, 2, max, na.rm = TRUE)
  iqr_functional <- bag_max - bag_min

  list(
    spec_cols = spec_cols,
    grid = grid,
    lower = bag_min - f_value * iqr_functional,
    upper = bag_max + f_value * iqr_functional,
    f_value = f_value,
    central_prop = central_prop,
    n_train = nrow(X)
  )
}

apply_train_global_functional_qc <- function(data,
                                             qc_model,
                                             allowed_exceed_fraction = 0) {
  miss <- setdiff(qc_model$spec_cols, names(data))
  if (length(miss) > 0) stop("QC columns missing from the dataset.")

  X <- as.matrix(data[, qc_model$spec_cols, drop = FALSE])
  storage.mode(X) <- "double"
  below <- sweep(X, 2, qc_model$lower, FUN = "<")
  above <- sweep(X, 2, qc_model$upper, FUN = ">")
  exceed_fraction <- rowMeans(below | above, na.rm = TRUE)
  is_outlier <- exceed_fraction > allowed_exceed_fraction

  list(is_outlier = is_outlier, exceed_fraction = exceed_fraction)
}

apply_train_only_farm_qc_to_views <- function(split_orig,
                                               split_d1 = NULL,
                                               split_d2 = NULL,
                                               f_value = 1,
                                               allowed_exceed_fraction = 0,
                                               remove_test = FALSE,
                                               farm_col = "Azienda",
                                               id_col = "Matricola",
                                               min_train_per_farm = 6L,
                                               col_start = 240,
                                               col_end = 1299) {
  if (isTRUE(remove_test)) {
    stop("This sensitivity design keeps ALL test spectra: remove_test must be FALSE.")
  }
  if (is.null(split_orig) || length(split_orig) != 2L) {
    stop("Provide the two-part (train/test) smoothed-spectra split.")
  }
  tr0 <- split_orig[[1]]
  te0 <- split_orig[[2]]
  needed <- c(farm_col, id_col)
  if (!all(needed %in% names(tr0)) || !all(needed %in% names(te0))) {
    stop("Farm or sample identifier missing from the smoothed spectra.")
  }
  if (anyNA(tr0[[farm_col]]) || anyNA(te0[[farm_col]]) ||
      anyNA(tr0[[id_col]]) || anyNA(te0[[id_col]])) {
    stop("Missing farm or sample identifiers in the analysis data.")
  }

  assert_aligned <- function(sp, view_label) {
    if (is.null(sp)) return(invisible(TRUE))
    if (length(sp) != 2L) stop("Invalid split for view: ", view_label)
    for (j in seq_len(2L)) {
      base <- split_orig[[j]]
      v <- sp[[j]]
      if (nrow(v) != nrow(base) ||
          !identical(as.character(v[[farm_col]]), as.character(base[[farm_col]])) ||
          !identical(as.character(v[[id_col]]), as.character(base[[id_col]]))) {
        stop("View misaligned with m0 (farm/sample/order): ", view_label)
      }
    }
    invisible(TRUE)
  }
  assert_aligned(split_d1, "d1")
  assert_aligned(split_d2, "d2")

  train_farms <- as.character(tr0[[farm_col]])
  test_farms <- as.character(te0[[farm_col]])
  train_flags <- rep(FALSE, nrow(tr0))
  test_flags <- rep(FALSE, nrow(te0))
  farm_rows <- vector("list", length(unique(train_farms)))
  farm_names <- sort(unique(train_farms))

  for (fi in seq_along(farm_names)) {
    farm <- farm_names[[fi]]
    idx_tr <- which(train_farms == farm)
    idx_te <- which(test_farms == farm)
    if (length(idx_tr) < min_train_per_farm) {
      stop("Farm ", farm, " has only ", length(idx_tr),
           " training spectra; minimum required is ", min_train_per_farm, ".")
    }
    farm_train <- tr0[idx_tr, , drop = FALSE]
    model <- fit_train_global_functional_qc(
      train_smoothed = farm_train,
      f_value = f_value,
      col_start = col_start,
      col_end = col_end
    )
    flagged_tr <- apply_train_global_functional_qc(
      farm_train, model, allowed_exceed_fraction
    )$is_outlier
    if (anyNA(flagged_tr)) stop("Undefined QC decision in farm ", farm)
    n_after <- length(idx_tr) - sum(flagged_tr)
    if (n_after < min_train_per_farm) {
      stop("Farm ", farm, " would retain only ", n_after,
           " training spectra after QC. Review f_value and the data.")
    }
    train_flags[idx_tr] <- flagged_tr

    if (length(idx_te) > 0L) {
      flagged_te <- apply_train_global_functional_qc(
        te0[idx_te, , drop = FALSE], model, allowed_exceed_fraction
      )$is_outlier
      if (anyNA(flagged_te)) stop("Undefined test QC decision in farm ", farm)
      test_flags[idx_te] <- flagged_te
    }

    farm_rows[[fi]] <- data.frame(
      Azienda = farm,
      train_before = length(idx_tr),
      train_removed = sum(flagged_tr),
      train_after = n_after,
      test_n = length(idx_te),
      test_flagged = if (length(idx_te)) sum(test_flags[idx_te]) else 0L,
      test_removed = 0L,
      stringsAsFactors = FALSE
    )
  }
  missing_test_farms <- setdiff(unique(test_farms), farm_names)
  if (length(missing_test_farms)) {
    stop("No corresponding training data for test farms: ",
         paste(missing_test_farms, collapse = ", "))
  }

  keep_train <- !train_flags
  apply_keep <- function(sp) {
    if (is.null(sp)) return(NULL)
    list(sp[[1]][keep_train, , drop = FALSE], sp[[2]])
  }
  list(
    orig = apply_keep(split_orig),
    d1   = apply_keep(split_d1),
    d2   = apply_keep(split_d2),
    info = list(
      mode = "train_within_farm_qc",
      train_n_before = nrow(tr0),
      train_removed = sum(train_flags),
      train_n_after = sum(keep_train),
      test_n_before = nrow(te0),
      test_flagged = sum(test_flags),
      test_removed = 0L,
      test_n_used = nrow(te0),
      remove_test = FALSE,
      f_value = f_value,
      allowed_exceed_fraction = allowed_exceed_fraction,
      by_farm = do.call(rbind, farm_rows)
    )
  )
}

# Fixed-configuration, sample-level sensitivity comparison using nested SVM tuning.
run_outlier_sensitivity_nested <- function(data, breed, seed = 69L,
                                           outer_folds = 10L,
                                           inner_folds = 5L,
                                           pve = 0.9999,
                                           polynomial_degree = 4L,
                                           fence_multiplier = 1,
                                           allowed_exceed_fraction = 0) {
  intervals <- ftir_domain_intervals("chemical")
  model_grid <- default_model_grid("svm")
  model_grid <- model_grid[model_grid$kernel == "linear", , drop = FALSE]
  scenarios <- c("no_removal", "within_farm_qc")
  outer_plan <- make_stratified_sample_folds(data$Altezza, outer_folds, seed)
  n <- nrow(data)
  fold_results <- vector("list", length(outer_plan) * length(scenarios))
  predictions <- vector("list", length(fold_results))
  qc_details <- vector("list", length(fold_results))
  count <- 0L

  prepare_views <- function(preprocessed, scenario) {
    if (scenario == "no_removal") {
      farms <- as.character(preprocessed$train$m0$Azienda)
      farm_ids <- sort(unique(farms))
      farm_counts <- data.frame(
        Azienda = farm_ids,
        train_before = vapply(farm_ids, function(farm) sum(farms == farm), integer(1)),
        train_removed = 0L,
        train_after = vapply(farm_ids, function(farm) sum(farms == farm), integer(1)),
        test_n = vapply(farm_ids, function(farm) {
          sum(as.character(preprocessed$test$m0$Azienda) == farm)
        }, integer(1)),
        test_flagged = 0L, test_removed = 0L,
        stringsAsFactors = FALSE
      )
      info <- list(train_removed = 0L, test_removed = 0L,
                   test_flagged = 0L, by_farm = farm_counts)
      return(list(train = preprocessed$train[c("m0", "m1")],
                  test = preprocessed$test[c("m0", "m1")], qc = info))
    }
    qc <- apply_train_only_farm_qc_to_views(
      split_orig = list(preprocessed$train$m0, preprocessed$test$m0),
      split_d1 = list(preprocessed$train$m1, preprocessed$test$m1),
      f_value = fence_multiplier,
      allowed_exceed_fraction = allowed_exceed_fraction,
      remove_test = FALSE,
      farm_col = "Azienda", id_col = "Matricola",
      min_train_per_farm = 6L,
      col_start = 240, col_end = 1299
    )
    list(train = list(m0 = qc$orig[[1L]], m1 = qc$d1[[1L]]),
         test = list(m0 = qc$orig[[2L]], m1 = qc$d1[[2L]]),
         qc = qc$info)
  }

  extract_features <- function(views) {
    cache <- build_preprocessed_fpca_cache(
      train_views = views$train, test_views = views$test,
      domains = "chemical", views_needed = c("m0", "m1"),
      max_pve = pve, id_col = "Matricola",
      group_col = "Azienda", target_col = "Altezza"
    )
    features <- assemble_fpca_features(
      cache, domain = "chemical", derivative_config = "m0_d1",
      pve = pve, target_col = "Altezza", scale_minmax = TRUE
    )
    if (is.null(features) || any(table(features$y_train) == 0L)) {
      stop("Missing training features or an altitude class after preprocessing.")
    }
    features
  }

  for (outer_id in seq_along(outer_plan)) {
    train_idx <- outer_plan[[outer_id]]
    test_idx <- setdiff(seq_len(n), train_idx)
    outer_train <- data[train_idx, , drop = FALSE]
    outer_test <- data[test_idx, , drop = FALSE]
    inner_plan <- make_stratified_sample_folds(
      outer_train$Altezza, inner_folds, seed + 1000L * outer_id
    )

    # Derive the Savitzky-Golay window from unfiltered training data in each split.
    inner_preprocessed <- lapply(inner_plan, function(inner_train_idx) {
      inner_test_idx <- setdiff(seq_len(nrow(outer_train)), inner_train_idx)
      preprocess_grouped_split(
        outer_train[inner_train_idx, , drop = FALSE],
        outer_train[inner_test_idx, , drop = FALSE],
        polynomial_degree = polynomial_degree,
        fallback_window = 7L, width_range = NULL
      )
    })
    outer_preprocessed <- preprocess_grouped_split(
      outer_train, outer_test, polynomial_degree = polynomial_degree,
      fallback_window = 7L, width_range = NULL
    )

    for (scenario in scenarios) {
      message(breed, " | outer fold ", outer_id, "/", length(outer_plan),
              " | ", scenario)
      inner_features <- lapply(inner_preprocessed, function(pp) {
        extract_features(prepare_views(pp, scenario))
      })
      tuning <- lapply(seq_len(nrow(model_grid)), function(param_id) {
        param <- model_grid[param_id, , drop = FALSE]
        fold_accuracy <- vapply(seq_along(inner_features), function(inner_id) {
          feat <- inner_features[[inner_id]]
          fit <- fit_model_candidate(
            feat$x_train, feat$y_train, "svm", param,
            seed = seed + 1000L * outer_id + 100L * param_id + inner_id
          )
          pred <- predict_model_candidate(fit, feat$x_test)
          mean(as.character(pred) == as.character(feat$y_test))
        }, numeric(1))
        data.frame(param_id = param$param_id,
                   mean_accuracy = mean(fold_accuracy),
                   sd_accuracy = stats::sd(fold_accuracy),
                   mean_features = mean(vapply(inner_features, `[[`, numeric(1),
                                               "n_features")))
      })
      tuning <- do.call(rbind, tuning)
      tuning <- tuning[order(-tuning$mean_accuracy,
                             tuning$sd_accuracy,
                             tuning$mean_features,
                             tuning$param_id), , drop = FALSE]
      best <- model_grid[model_grid$param_id == tuning$param_id[1L], , drop = FALSE]

      views <- prepare_views(outer_preprocessed, scenario)
      feat <- extract_features(views)
      fit <- fit_model_candidate(
        feat$x_train, feat$y_train, "svm", best,
        seed = seed + 10000L * outer_id + best$param_id
      )
      pred <- predict_model_candidate(fit, feat$x_test)
      truth <- factor(as.character(feat$y_test), levels = ALTITUDE_LEVELS)
      if (length(pred) != length(test_idx) ||
          !identical(as.character(feat$test_meta$Matricola),
                     as.character(outer_test$Matricola)) ||
          views$qc$test_removed != 0L) {
        stop("An outer-test sample was removed or its order changed.")
      }

      count <- count + 1L
      fold_results[[count]] <- data.frame(
        breed = breed, scenario = scenario, outer_fold = outer_id,
        accuracy = mean(as.character(pred) == as.character(truth)),
        macro_f1 = f1_macro(truth, pred),
        n_train_before = nrow(outer_train),
        n_train_after = nrow(feat$x_train),
        n_train_removed = views$qc$train_removed,
        n_test = length(test_idx), n_test_removed = views$qc$test_removed,
        n_test_flagged = views$qc$test_flagged,
        n_fpca_scores = feat$n_features,
        sg_window = outer_preprocessed$window_length,
        svm_kernel = best$kernel, svm_cost = best$cost,
        svm_gamma_multiplier = best$gamma_multiplier,
        inner_accuracy = tuning$mean_accuracy[1L],
        stringsAsFactors = FALSE
      )
      predictions[[count]] <- data.frame(
        breed = breed, scenario = scenario, outer_fold = outer_id,
        sample_id = as.character(outer_test$Matricola),
        observed = as.character(truth), predicted = as.character(pred),
        stringsAsFactors = FALSE
      )
      farms <- views$qc$by_farm
      farms$breed <- breed
      farms$scenario <- scenario
      farms$outer_fold <- outer_id
      qc_details[[count]] <- farms
    }
  }
  metrics <- do.call(rbind, fold_results)
  preds <- do.call(rbind, predictions)
  farms <- do.call(rbind, qc_details)

  for (scenario in scenarios) {
    ids <- preds$sample_id[preds$scenario == scenario]
    if (length(ids) != n || anyDuplicated(ids) ||
        !setequal(ids, as.character(data$Matricola))) {
      stop("Each sample must appear exactly once in each scenario's outer test sets.")
    }
  }
  a <- metrics[metrics$scenario == scenarios[1L], , drop = FALSE]
  b <- metrics[metrics$scenario == scenarios[2L], , drop = FALSE]
  if (!identical(a$sg_window, b$sg_window) ||
      !identical(a$n_test, b$n_test) ||
      any(metrics$n_test_removed != 0L)) {
    stop("Sensitivity scenarios must have matched folds/windows and no test exclusion.")
  }

  summary <- do.call(rbind, lapply(scenarios, function(scenario) {
    z <- metrics[metrics$scenario == scenario, , drop = FALSE]
    data.frame(breed = breed, scenario = scenario,
               mean_accuracy = mean(z$accuracy), sd_accuracy = stats::sd(z$accuracy),
               mean_macro_f1 = mean(z$macro_f1), sd_macro_f1 = stats::sd(z$macro_f1),
               mean_train_removed_per_fold = mean(z$n_train_removed),
               total_train_removal_occurrences = sum(z$n_train_removed),
               total_test_removed = sum(z$n_test_removed),
               stringsAsFactors = FALSE)
  }))
  list(summary = summary, fold_metrics = metrics,
       predictions = preds, farm_qc = farms)
}
