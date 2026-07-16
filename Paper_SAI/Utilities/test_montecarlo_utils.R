# Locate the Paper_SAI folder robustly (Rscript / RStudio "Source" / interactive
# console at/under/above the project) and load the shared validators.
.find_paper_root <- function() {
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
  stop("Cannot locate the Paper_SAI folder (all_parameters.R).", call. = FALSE)
}
source(file.path(.find_paper_root(), "Utilities", "montecarlo_utils.R"))

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

expect_error <- function(expr, pattern) {
  ok <- FALSE
  msg <- NULL
  tryCatch(
    force(expr),
    error = function(e) {
      msg <<- conditionMessage(e)
      ok <<- grepl(pattern, msg)
    }
  )
  if (!ok) stop("Expected error matching '", pattern, "', got: ", msg, call. = FALSE)
}

set.seed(11)
n <- 200000
hazard <- exp(runif(n, log(MC_HAZARD_LO), log(MC_HAZARD_HI)))
delay <- mc_geometric_delay(runif(n), hazard)
assert_true(abs(mean(delay == 1) - mean(hazard)) < 0.0015,
            "Geometric inverse-CDF does not match first-year hazard.")

pulse <- c(5L, 30L, 80L)
term_delay <- c(10L, 450L, 600L)
term_delta <- mc_term_delta(pulse, term_delay)
assert_true(identical(term_delta, c(15L, 480L, 480L)), "Censoring term_delta failed.")
assert_true(identical(mc_censored_termination(pulse, term_delay), c(FALSE, TRUE, TRUE)),
            "Censoring flag failed.")

valid <- data.frame(
  ID = 1:2,
  ecs = c(30, 31),
  tcr = c(18, 19),
  rcp = c("RCP45", "RCP6"),
  pulse = c(5, 30),
  cool = c(10, 20),
  term = c(2400, 2500),
  start = c(2025, 2050),
  term_delta = c(15, 480),
  theta = c(10, 11),
  alpha = c(0.005, 0.006),
  delta = c(0.02, 0.03),
  prob = c(0.01, 0.02),
  mortality_srm = c(7400, 7600),
  forctoTg = c(2, 3),
  TgtoUSD = c(1, 2),
  mortality_ozone = c(112, 113),
  vsl = c(1e7, 1.1e7),
  vsl_eta = c(1, 0.9),
  dg = c(0.015, 0.016),
  sampler_version = MC_SAMPLER_VERSION,
  sampling_method = "montecarlo",
  seed = 123,
  draw_index = 1:2,
  term_delay = c(10, 500),
  censored_termination = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
validate_id_montecarlo(valid, "valid fixture")

legacy <- valid[setdiff(names(valid), "sampler_version")]
expect_error(validate_id_montecarlo(legacy, "legacy fixture"), "missing required columns")

dupe <- valid
dupe$draw_index[2] <- dupe$draw_index[1]
expect_error(validate_id_montecarlo(dupe, "duplicate fixture"), "duplicate draw_index")

bad_term <- valid
bad_term$term_delta[1] <- bad_term$pulse[1]
expect_error(validate_id_montecarlo(bad_term, "bad term fixture"), "term_delta must be strictly greater")

complete_window <- data.frame(
  ID = c(1, 1, 1, 2, 2),
  gas = c("co2", "co2", "co2", "ch4", "ch4"),
  pulse_time = c(1, 1, 1, 2, 2),
  t = c(1, 2, 3, 2, 3)
)
mc_assert_complete_time_window(complete_window, c("ID", "gas"),
                               "complete time-window fixture", end_t = 3)

missing_window <- complete_window[!(complete_window$ID == 1 & complete_window$t == 2), ]
expect_error(
  mc_assert_complete_time_window(missing_window, c("ID", "gas"),
                                 "missing time-window fixture", end_t = 3),
  "incomplete time windows"
)

duplicate_window <- rbind(complete_window, complete_window[1, ])
expect_error(
  mc_assert_complete_time_window(duplicate_window, c("ID", "gas"),
                                 "duplicate time-window fixture", end_t = 3),
  "incomplete time windows"
)

# mc_assert_unique_key must select columns, not trigger a data.table keyed join.
# Every caller in Analyze_montecarlo.R passes a data.table, so cover that input
# type for unkeyed, keyed and plain data.frame cases.
unique_keys_dt <- data.table::data.table(ID = c("1", "2"), gas = c("co2", "ch4"), v = c(1, 2))
assert_true(mc_assert_unique_key(unique_keys_dt, c("ID", "gas"), "unkeyed data.table fixture"),
            "unique keys in an unkeyed data.table must pass")
assert_true(mc_assert_unique_key(unique_keys_dt[, .(ID, v)], "ID", "single-column key fixture"),
            "a single-column key on a data.table must pass")
assert_true(mc_assert_unique_key(as.data.frame(unique_keys_dt), c("ID", "gas"), "data.frame fixture"),
            "unique keys in a data.frame must pass")

duplicate_keys_dt <- data.table::data.table(ID = c("1", "1"), gas = c("co2", "co2"), v = c(1, 2))
expect_error(
  mc_assert_unique_key(duplicate_keys_dt, c("ID", "gas"), "duplicate key fixture"),
  "contains duplicate key rows"
)

keyed_duplicates_dt <- data.table::copy(duplicate_keys_dt)
data.table::setkeyv(keyed_duplicates_dt, "ID")
expect_error(
  mc_assert_unique_key(keyed_duplicates_dt, c("ID", "gas"), "keyed duplicate fixture"),
  "contains duplicate key rows"
)

expect_error(
  mc_assert_unique_key(unique_keys_dt, c("ID", "nope"), "missing key column fixture"),
  "missing key columns"
)

cat("Monte Carlo utility tests passed.\n")
