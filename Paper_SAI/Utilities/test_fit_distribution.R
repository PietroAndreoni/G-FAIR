# =============================================================================
# test_fit_distribution.R
#
# Correctness harness for fit_distribution() defined in Generate_montecarlo.R.
# It loads the LIVE function (extracted from the script, single source of truth)
# and exercises every input combination, classifying each outcome as:
#   PASS       - behaves as intended
#   IMPOSSIBLE - inputs cannot be satisfied by the requested family (the
#                function returns a documented least-discrepancy compromise)
#   BUG        - a self-consistent input is NOT recovered, or an expected error
#                does not fire  (=> exits non-zero)
#
# Run:  Rscript test_fit_distribution.R
# =============================================================================

# ---- load the live functions (single source of truth) -----------------------
# Locate the Paper_SAI folder robustly (Rscript / RStudio "Source" / interactive
# console at/under/above the project) and load the live functions.
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
source(file.path(.find_paper_root(), "Utilities", "distribution_utils.R"))

# ---- ground-truth theoretical moments/quantiles -----------------------------
ln_inputs <- function(mu, s) list(
  median = exp(mu),
  mean   = exp(mu + s^2 / 2),
  sd     = sqrt((exp(s^2) - 1)) * exp(mu + s^2 / 2),
  q33 = exp(mu + s * qnorm(.33)), q66 = exp(mu + s * qnorm(.66)),
  q5  = exp(mu + s * qnorm(.05)), q95 = exp(mu + s * qnorm(.95))
)
nm_inputs <- function(mu, s) list(
  median = mu, mean = mu, sd = s,
  q33 = mu + s * qnorm(.33), q66 = mu + s * qnorm(.66),
  q5  = mu + s * qnorm(.05), q95 = mu + s * qnorm(.95)
)

# branch -> which inputs it consumes
branches <- list(
  lognormal = list(
    "mean&sd"        = c("mean","sd"),
    "mean&median"    = c("mean","median"),
    "median&sd"      = c("median","sd"),
    "median&q33&q66" = c("median","q33","q66"),
    "median&q5&q95"  = c("median","q5","q95"),
    "q33&q66"        = c("q33","q66"),
    "q5&q95"         = c("q5","q95")
  ),
  normal = list(
    "mean&sd"        = c("mean","sd"),
    "median&sd"      = c("median","sd"),
    "median&q33&q66" = c("median","q33","q66"),
    "median&q5&q95"  = c("median","q5","q95"),
    "q33&q66"        = c("q33","q66"),
    "q5&q95"         = c("q5","q95")
  )
)

results <- list()
add <- function(family, branch, scenario, detail, verdict)
  results[[length(results) + 1]] <<- data.frame(
    family = family, branch = branch, scenario = scenario,
    detail = detail, verdict = verdict, stringsAsFactors = FALSE)

quiet <- function(expr) suppressWarnings(expr)

# ---- TEST 1: branch recovery from self-consistent inputs (BUG detector) -----
truth <- list(lognormal = c(mu = 1.0, s = 0.5), normal = c(mu = 5.0, s = 2.0))
tol_par <- 1e-6
for (fam in names(branches)) {
  g  <- truth[[fam]]; mu0 <- unname(g["mu"]); s0 <- unname(g["s"])
  inp <- if (fam == "lognormal") ln_inputs(mu0, s0) else nm_inputs(mu0, s0)
  for (br in names(branches[[fam]])) {
    args <- c(list(distribution = fam, return_params = TRUE), inp[branches[[fam]][[br]]])
    got  <- tryCatch(quiet(do.call(fit_distribution, args)), error = function(e) e)
    if (inherits(got, "error")) {
      add(fam, br, "consistent", paste("unexpected error:", conditionMessage(got)), "BUG")
    } else {
      ok <- abs(unname(got$mu) - mu0) < tol_par && abs(unname(got$sigma) - s0) < tol_par
      add(fam, br, "consistent",
          sprintf("mu=%.5f (exp %.5f), sigma=%.5f (exp %.5f)", got$mu, mu0, got$sigma, s0),
          if (isTRUE(ok)) "PASS" else "BUG")
    }
  }
}

# ---- TEST 2: sampling sanity - rng output matches the fitted theoretical ----
set.seed(20240101)
for (fam in names(branches)) {
  g  <- truth[[fam]]; mu0 <- unname(g["mu"]); s0 <- unname(g["s"])
  inp <- if (fam == "lognormal") ln_inputs(mu0, s0) else nm_inputs(mu0, s0)
  br  <- "median&q5&q95"
  args <- c(list(distribution = fam, n = 2e5), inp[branches[[fam]][[br]]])
  draws <- quiet(do.call(fit_distribution, args))
  qf <- if (fam == "lognormal") function(p) qlnorm(p, mu0, s0) else function(p) qnorm(p, mu0, s0)
  emp <- quantile(draws, c(.05, .5, .95), names = FALSE); theo <- qf(c(.05, .5, .95))
  rel <- max(abs(emp - theo) / pmax(abs(theo), 1e-9))
  add(fam, br, "sampling(n=2e5)",
      sprintf("emp med=%.4g q5=%.4g q95=%.4g | max rel.err=%.1f%%", emp[2], emp[1], emp[3], 100 * rel),
      if (rel < 0.05) "PASS" else "BUG")
}

# ---- TEST 3: over-determined / inconsistent inputs (IMPOSSIBLE detector) ----
# lognormal feasibility: q5*q95 == median^2 ; normal: q5+q95 == 2*median
report_over <- function(fam, label, median, q5, q95) {
  got <- tryCatch(quiet(fit_distribution(fam, median = median, q5 = q5, q95 = q95, return_params = TRUE)),
                  error = function(e) e)
  if (inherits(got, "error")) { add(fam, "median&q5&q95", label, paste("error:", conditionMessage(got)), "BUG"); return() }
  if (fam == "lognormal") { rq5 <- qlnorm(.05, got$mu, got$sigma); rq95 <- qlnorm(.95, got$mu, got$sigma); ratio <- q5 * q95 / median^2 }
  else                    { rq5 <- qnorm(.05, got$mu, got$sigma);  rq95 <- qnorm(.95, got$mu, got$sigma);  ratio <- (q5 + q95) / (2 * median) }
  consistent <- abs(ratio - 1) < 1e-6
  add(fam, "median&q5&q95", label,
      sprintf("req q5/q95=%.4g/%.4g -> got %.4g/%.4g (consistency=%.3f)", q5, q95, rq5, rq95, ratio),
      if (consistent) "PASS" else "IMPOSSIBLE")
}
# the real Generate calls
report_over("lognormal", "ecs",            3,    2,    5)
report_over("lognormal", "tcr",            1.8,  1.2,  2.4)
report_over("lognormal", "theta",          10,   3,    30)
report_over("lognormal", "mortality_srm",  7400, 2300, 16000)
report_over("normal",    "mortality_ozone",11250,5000, 17500)
# a deliberately consistent control for each family
report_over("lognormal", "consistent-ctrl", exp(1), exp(1 + .5*qnorm(.05)), exp(1 + .5*qnorm(.95)))
report_over("normal",    "consistent-ctrl", 5, 5 + 2*qnorm(.05), 5 + 2*qnorm(.95))

# ---- TEST 4: error / edge paths (expect graceful stop) ----------------------
expect_stop <- function(fam, label, expr) {
  r <- tryCatch({ quiet(expr); "RETURNED" }, error = function(e) "STOPPED")
  add(fam, "error-path", label, r, if (r == "STOPPED") "PASS" else "BUG")
}
expect_stop("lognormal", "mean<median",   fit_distribution("lognormal", mean = 2, median = 5, return_params = TRUE))
expect_stop("lognormal", "negative q5",   fit_distribution("lognormal", median = 3, q5 = -1, q95 = 5, return_params = TRUE))
expect_stop("lognormal", "only mean",     fit_distribution("lognormal", mean = 3, return_params = TRUE))
expect_stop("normal",    "only median",   fit_distribution("normal", median = 3, return_params = TRUE))
expect_stop("either",    "bad family",    fit_distribution("poisson", mean = 1, sd = 1, return_params = TRUE))

# ---- TEST 5: ignored-input warning fires (and value is unaffected) ----------
ws <- character(0)
val <- withCallingHandlers(
  fit_distribution("lognormal", mean = 99, q5 = exp(1 + .5*qnorm(.05)), q95 = exp(1 + .5*qnorm(.95)), return_params = TRUE),
  warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
fired <- any(grepl("ignored input", ws))
add("lognormal", "ignored-input", "mean+q5+q95",
    sprintf("mu=%.4f (mean ignored), warning fired=%s", val$mu, fired),
    if (fired && isTRUE(all.equal(val$mu, 1, tolerance = 1e-6))) "PASS" else "BUG")

# ---- TEST 6: fit_quantile_dist honors over-determined quantiles -------------
# The 3-parameter (split-normal) families should match median+q5+q95 that the
# 2-parameter normal/lognormal cannot.
for (cs in list(list("ecs", 3, 2, 5), list("mortality_srm", 7400, 2300, 16000))) {
  lbl <- cs[[1]]; md <- cs[[2]]; q5 <- cs[[3]]; q95 <- cs[[4]]
  fqd <- quiet(fit_quantile_dist(median = md, q5 = q5, q95 = q95))
  hits <- c(fqd$q(0.05), fqd$q(0.5), fqd$q(0.95))
  relerr <- max(abs(hits / c(q5, md, q95) - 1))
  add("multi-family", "fit_quantile_dist", lbl,
      sprintf("best=%s residual=%.1e | q5/med/q95 -> %.4g/%.4g/%.4g (max rel.err=%.1e)",
              fqd$best, fqd$residual, hits[1], hits[2], hits[3], relerr),
      if (relerr < 1e-3) "PASS" else "BUG")
}

# ---- TEST 7: truncated inverse-CDF stays inside [lo,hi] ---------------------
set.seed(7); u <- runif(5e4)
th <- qlnorm_fit_trunc(u, lo = 0, hi = 90, median = 10, q5 = 3, q95 = 30)
add("truncation", "qlnorm_fit_trunc", "theta (0,90]",
    sprintf("range=[%.3f, %.3f]  frac_outside=%.4f", min(th), max(th), mean(th <= 0 | th > 90)),
    if (all(th > 0 & th <= 90)) "PASS" else "BUG")
oz <- qnorm_fit_trunc(u, lo = 0, median = 11250, q5 = 5000, q95 = 17500)
add("truncation", "qnorm_fit_trunc", "ozone >=0",
    sprintf("min=%.2f  frac<0=%.4f", min(oz), mean(oz < 0)),
    if (all(oz >= 0)) "PASS" else "BUG")

# ---- TEST 8: split-normal internal consistency ------------------------------
si <- abs(integrate(function(x) dsnorm(x, 2, 1, 3), -40, 60)$value - 1)
pq <- abs(psnorm(qsnorm(0.37, 2, 1, 3), 2, 1, 3) - 0.37)
add("split-normal", "dsnorm/psnorm/qsnorm", "integrate & inverse",
    sprintf("|integral-1|=%.1e  |psnorm(qsnorm)-p|=%.1e", si, pq),
    if (si < 1e-4 && pq < 1e-6) "PASS" else "BUG")

# ---- report -----------------------------------------------------------------
res <- do.call(rbind, results)
cat("\n================ fit_distribution test report ================\n")
print(res, row.names = FALSE, right = FALSE)
n_bug  <- sum(res$verdict == "BUG")
n_imp  <- sum(res$verdict == "IMPOSSIBLE")
n_pass <- sum(res$verdict == "PASS")
cat(sprintf("\nPASS=%d  IMPOSSIBLE=%d  BUG=%d\n", n_pass, n_imp, n_bug))
if (n_imp > 0) {
  cat("\nIMPOSSIBLE rows are inputs that no single distribution of that family can\n",
      "satisfy (lognormal needs q5*q95=median^2; normal needs q5+q95=2*median).\n",
      "fit_distribution returns a least-discrepancy compromise and warns.\n", sep = "")
}
if (n_bug > 0) { cat("\nFAILURE: code bug(s) detected.\n"); quit(status = 1) }
cat("\nAll branches correct (no code bugs).\n")
