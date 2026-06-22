# Shared constants and strict validators for the Monte Carlo pipeline.

MC_SAMPLER_VERSION <- "mc_sampler_v2"
MC_T_HORIZON <- 480L
MC_HAZARD_LO <- 1e-3
MC_HAZARD_HI <- 1e-1

MC_FAIR_COLS <- c("ecs", "tcr", "rcp", "pulse", "cool", "term", "start", "term_delta")
MC_POST_COLS <- c("theta", "alpha", "delta", "prob", "mortality_srm", "forctoTg",
                  "TgtoUSD", "mortality_ozone", "vsl", "vsl_eta", "dg")
MC_METADATA_COLS <- c("sampler_version", "sampling_method", "seed", "draw_index",
                      "term_delay", "censored_termination")
MC_REQUIRED_ID_COLS <- c("ID", MC_FAIR_COLS, MC_POST_COLS, MC_METADATA_COLS)

mc_strip_csv_index <- function(x) {
  if ("X" %in% names(x)) x <- x[setdiff(names(x), "X")]
  x
}

mc_is_power_of_two <- function(n) {
  n <- as.numeric(n)
  is.finite(n) && n >= 1 && n == floor(n) && bitwAnd(as.integer(n), as.integer(n) - 1L) == 0L
}

mc_geometric_delay <- function(u, hazard) {
  if (length(u) != length(hazard)) stop("u and hazard must have the same length.")
  if (any(!is.finite(u)) || any(!is.finite(hazard))) stop("u and hazard must be finite.")
  if (any(hazard <= 0 | hazard >= 1)) stop("hazard must be in (0, 1).")
  u <- pmin(pmax(u, 0), 1 - .Machine$double.eps)
  as.integer(pmax(ceiling(log1p(-u) / log1p(-hazard)), 1L))
}

mc_term_delta <- function(pulse, term_delay, t_horizon = MC_T_HORIZON) {
  as.integer(pmin(as.integer(pulse) + as.integer(term_delay), as.integer(t_horizon)))
}

mc_censored_termination <- function(pulse, term_delay, t_horizon = MC_T_HORIZON) {
  as.integer(pulse) + as.integer(term_delay) >= as.integer(t_horizon)
}

mc_stop <- function(context, errors) {
  stop(paste0(context, " failed strict Monte Carlo validation:\n- ",
              paste(errors, collapse = "\n- ")), call. = FALSE)
}

mc_assert_no_missing <- function(x, cols, context, max_rows = 5L) {
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0) {
    stop(context, " missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  bad_cols <- cols[vapply(cols, function(cname) any(is.na(x[[cname]])), logical(1))]
  if (length(bad_cols) == 0) return(invisible(TRUE))

  bad <- which(Reduce(`|`, lapply(bad_cols, function(cname) is.na(x[[cname]]))))
  details <- paste(utils::capture.output(print(utils::head(x[bad, intersect(names(x), cols), drop = FALSE], max_rows))),
                   collapse = "\n")
  stop(context, " contains missing values in: ", paste(bad_cols, collapse = ", "),
       "\nExamples:\n", details, call. = FALSE)
}

mc_assert_unique_key <- function(x, key_cols, context) {
  missing_cols <- setdiff(key_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(context, " missing key columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  key_df <- x[key_cols]
  dup <- duplicated(key_df) | duplicated(key_df, fromLast = TRUE)
  if (!any(dup)) return(invisible(TRUE))
  details <- paste(utils::capture.output(print(utils::head(key_df[dup, , drop = FALSE], 10))),
                   collapse = "\n")
  stop(context, " contains duplicate key rows for: ", paste(key_cols, collapse = ", "),
       "\nExamples:\n", details, call. = FALSE)
}

mc_assert_key_set_equal <- function(expected, observed, key_cols, context) {
  missing_cols <- setdiff(key_cols, intersect(names(expected), names(observed)))
  if (length(missing_cols) > 0) {
    stop(context, " cannot compare key sets; missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  expected_key <- unique(expected[key_cols])
  observed_key <- unique(observed[key_cols])
  expected_s <- do.call(paste, c(expected_key, sep = "\r"))
  observed_s <- do.call(paste, c(observed_key, sep = "\r"))
  missing_s <- setdiff(expected_s, observed_s)
  extra_s <- setdiff(observed_s, expected_s)
  if (length(missing_s) == 0 && length(extra_s) == 0) return(invisible(TRUE))
  msg <- character(0)
  if (length(missing_s) > 0) {
    msg <- c(msg, paste0("missing expected keys: ", paste(utils::head(missing_s, 10), collapse = " | ")))
  }
  if (length(extra_s) > 0) {
    msg <- c(msg, paste0("unexpected extra keys: ", paste(utils::head(extra_s, 10), collapse = " | ")))
  }
  stop(context, " key-set mismatch: ", paste(msg, collapse = "; "), call. = FALSE)
}

validate_id_montecarlo <- function(x, context = "id_montecarlo.csv",
                                   require_current_version = TRUE,
                                   allow_empty = FALSE) {
  x <- mc_strip_csv_index(as.data.frame(x, stringsAsFactors = FALSE))
  errors <- character(0)

  if (nrow(x) == 0 && !allow_empty) errors <- c(errors, "file has no realizations")

  missing_cols <- setdiff(MC_REQUIRED_ID_COLS, names(x))
  if (length(missing_cols) > 0) {
    mc_stop(context, c(errors, paste0("missing required columns: ", paste(missing_cols, collapse = ", "))))
  }

  mc_assert_no_missing(x, MC_REQUIRED_ID_COLS, context)

  if (require_current_version && any(as.character(x$sampler_version) != MC_SAMPLER_VERSION)) {
    errors <- c(errors, paste0("sampler_version must be ", MC_SAMPLER_VERSION))
  }
  if (any(!as.character(x$sampling_method) %in% c("sobol", "montecarlo"))) {
    errors <- c(errors, "sampling_method must be 'sobol' or 'montecarlo'")
  }

  num_cols <- setdiff(MC_REQUIRED_ID_COLS, c("rcp", "sampler_version", "sampling_method", "censored_termination"))
  nums <- lapply(x[num_cols], function(v) suppressWarnings(as.numeric(v)))
  bad_num <- names(nums)[vapply(nums, function(v) any(is.na(v) | !is.finite(v)), logical(1))]
  if (length(bad_num) > 0) errors <- c(errors, paste0("non-finite numeric values in: ", paste(bad_num, collapse = ", ")))

  if (length(bad_num) == 0) {
    pulse <- as.integer(round(as.numeric(x$pulse)))
    term_delta <- as.integer(round(as.numeric(x$term_delta)))
    term_delay <- as.integer(round(as.numeric(x$term_delay)))
    expected_term <- mc_term_delta(pulse, term_delay)
    expected_censored <- mc_censored_termination(pulse, term_delay)
    censored <- as.logical(x$censored_termination)

    if (any(term_delta > MC_T_HORIZON)) errors <- c(errors, paste0("term_delta exceeds horizon ", MC_T_HORIZON))
    if (any(term_delta <= pulse)) errors <- c(errors, "term_delta must be strictly greater than pulse")
    if (any(term_delay < 1)) errors <- c(errors, "term_delay must be >= 1")
    if (any(term_delta != expected_term)) errors <- c(errors, "term_delta must equal pmin(pulse + term_delay, horizon)")
    if (any(is.na(censored) | censored != expected_censored)) {
      errors <- c(errors, "censored_termination must equal pulse + term_delay >= horizon")
    }
    if (any(as.numeric(x$prob) <= 0 | as.numeric(x$prob) >= 1)) {
      errors <- c(errors, "prob must be a per-year hazard in (0, 1)")
    }
  }

  if (any(duplicated(as.integer(round(as.numeric(x$ID)))))) errors <- c(errors, "duplicate ID values")
  if (any(duplicated(as.integer(round(as.numeric(x$draw_index)))))) errors <- c(errors, "duplicate draw_index values")

  if (length(errors) > 0) mc_stop(context, errors)
  invisible(x)
}
