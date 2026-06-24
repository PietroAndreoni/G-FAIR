source("montecarlo_utils.R")

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

cat("Monte Carlo utility tests passed.\n")
