# Shared constants and strict validators for the Monte Carlo pipeline.

# Canonical numeric parameters now live in all_parameters.R (the single control
# file). Source it once here (unless a caller already did) and expose the
# historical MC_* names as aliases. Located robustly so this works under Rscript,
# RStudio "Source", or an interactive console at/under/above the project.
if (!exists("T_HORIZON")) {
  .mc_find_root <- function() {
    starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                unlist(lapply(sys.frames(), function(f) f$ofile)), getwd())
    for (s in starts[nzchar(starts)]) {
      d <- if (dir.exists(s)) s else dirname(s)
      repeat {
        if (file.exists(file.path(d, "all_parameters.R"))) return(d)
        if (file.exists(file.path(d, "Paper_SAI", "all_parameters.R"))) return(file.path(d, "Paper_SAI"))
        if (identical(dirname(d), d)) break
        d <- dirname(d)
      }
    }
    stop("montecarlo_utils.R: cannot locate all_parameters.R (the control file).", call. = FALSE)
  }
  source(file.path(.mc_find_root(), "all_parameters.R"))
}

MC_SAMPLER_VERSION <- SAMPLER_VERSION
MC_T_HORIZON <- as.integer(T_HORIZON)
MC_HAZARD_LO <- HAZARD_LO
MC_HAZARD_HI <- HAZARD_HI

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
  # Column subset via data.frame semantics so this works for both data.frame and
  # data.table inputs (DT[<character>] would be a keyed-join, not a column select).
  key_df <- as.data.frame(x, stringsAsFactors = FALSE)[key_cols]
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
  # Column subset via data.frame semantics so this works for both data.frame and
  # data.table inputs (DT[<character>] would be a keyed-join, not a column select).
  expected_key <- unique(as.data.frame(expected, stringsAsFactors = FALSE)[key_cols])
  observed_key <- unique(as.data.frame(observed, stringsAsFactors = FALSE)[key_cols])
  expected_s <- do.call(paste, c(expected_key, sep = "\r"))
  observed_s <- do.call(paste, c(observed_key, sep = "\r"))
  missing_s <- setdiff(expected_s, observed_s)
  extra_s <- setdiff(observed_s, expected_s)
  if (length(missing_s) == 0 && length(extra_s) == 0) return(invisible(TRUE))
  # Keys are built by paste()-ing the key columns with a "\r" separator (chosen so
  # it can't collide with real values); turn that into a visible "/" for the
  # message, otherwise the carriage returns overwrite the text in the console.
  fmt_keys <- function(s) gsub("\r", "/", utils::head(s, 10), fixed = TRUE)
  msg <- character(0)
  if (length(missing_s) > 0) {
    msg <- c(msg, paste0("missing expected keys: ", paste(fmt_keys(missing_s), collapse = " | ")))
  }
  if (length(extra_s) > 0) {
    msg <- c(msg, paste0("unexpected extra keys: ", paste(fmt_keys(extra_s), collapse = " | ")))
  }
  stop(context, " key-set mismatch: ", paste(msg, collapse = "; "), call. = FALSE)
}

mc_assert_complete_time_window <- function(x, key_cols, context,
                                           t_col = "t",
                                           start_col = "pulse_time",
                                           end_t = MC_T_HORIZON,
                                           max_rows = 10L) {
  missing_cols <- setdiff(c(key_cols, t_col, start_col), names(x))
  if (length(missing_cols) > 0) {
    stop(context, " missing columns for time-window validation: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(context, " requires package data.table for time-window validation", call. = FALSE)
  }

  DT <- data.table::as.data.table(x)
  summary <- DT[, {
    t_raw <- suppressWarnings(as.numeric(get(t_col)))
    start_raw <- suppressWarnings(as.numeric(get(start_col)))
    bad_t <- any(is.na(t_raw) | !is.finite(t_raw) | t_raw != floor(t_raw))
    bad_start <- any(is.na(start_raw) | !is.finite(start_raw) | start_raw != floor(start_raw))
    t_vals <- sort(unique(as.integer(t_raw[!is.na(t_raw) & is.finite(t_raw) & t_raw == floor(t_raw)])))
    start_vals <- unique(as.integer(start_raw[!is.na(start_raw) & is.finite(start_raw) & start_raw == floor(start_raw)]))
    start_ok <- !bad_start && length(start_vals) == 1L && start_vals <= as.integer(end_t)
    expected <- if (start_ok) seq.int(start_vals, as.integer(end_t)) else integer(0)
    missing_t <- setdiff(expected, t_vals)
    extra_t <- setdiff(t_vals, expected)
    duplicate_count <- length(t_raw) - length(t_vals)
    list(
      start_values = if (length(start_vals) == 0L) NA_character_ else paste(start_vals, collapse = ","),
      observed_years = length(t_vals),
      expected_years = length(expected),
      min_t = if (length(t_vals) == 0L) NA_integer_ else min(t_vals),
      max_t = if (length(t_vals) == 0L) NA_integer_ else max(t_vals),
      missing_count = length(missing_t),
      extra_count = length(extra_t),
      duplicate_count = duplicate_count,
      invalid_t = bad_t,
      invalid_start = bad_start || length(start_vals) != 1L || any(start_vals > as.integer(end_t)),
      first_missing_t = if (length(missing_t) == 0L) NA_character_ else paste(utils::head(missing_t, 10), collapse = ","),
      first_extra_t = if (length(extra_t) == 0L) NA_character_ else paste(utils::head(extra_t, 10), collapse = ",")
    )
  }, by = key_cols]

  bad <- summary[
    invalid_t | invalid_start | missing_count > 0L |
      extra_count > 0L | duplicate_count > 0L
  ]
  if (nrow(bad) == 0L) return(invisible(TRUE))

  details <- paste(utils::capture.output(print(utils::head(bad, max_rows))),
                   collapse = "\n")
  stop(context, " has incomplete time windows; expected every ",
       paste(key_cols, collapse = "/"), " group to cover ",
       start_col, ":", end_t, " exactly once.\nExamples:\n",
       details, call. = FALSE)
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
