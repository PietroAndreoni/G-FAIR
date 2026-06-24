# this scripts generates the realizations for the montecarlo parameters
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Generate_montecarlo.R [-o <res>] [-n <n_scenarios>] [-w <overwrite_data>] [-p <run_parallel>] [-s <start_job>] [-e <end_job>] [--hpc <run_hpc>] [--base <main_scenario>] [--angle <angle>] [--method <sampling_method>] [--seed <seed>] [--diagnostics <diagnostics>]

Options:
-o <res>                     Path of the input/outputs
-n <n_scenarios>             Number of new scenarios to generate
-w <overwrite_data>          T/F to overwrite old realizations
--seed <seed>                seed number (for reproducibility)
--base <main_scenario>       T/F if to run the main scenario only (no montecarlo over policy parameters)
--angle <angle>              theta value for --base=T (default: 10; use na to keep theta uncertain)
--method <sampling_method>   "sobol" (quasi-random, default) or "montecarlo" (pseudo-random)
--diagnostics <diagnostics>  T/F write a PDF comparing the drawn vs theoretical distributions (default F)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Generate_montecarlo')

# Shared distribution-fitting / sampling helpers: fit_distribution,
# qlnorm_fit/qnorm_fit, truncated inverse-CDF samplers, split-normal and
# fit_quantile_dist. Kept in one file so Generate and the test harness agree.
.utils <- "distribution_utils.R"
.sp <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(.sp) == 1 && file.exists(file.path(dirname(.sp), .utils)))
  .utils <- file.path(dirname(.sp), .utils)
source(.utils)

.mc_utils <- "montecarlo_utils.R"
if (length(.sp) == 1 && file.exists(file.path(dirname(.sp), .mc_utils)))
  .mc_utils <- file.path(dirname(.sp), .mc_utils)
source(.mc_utils)

# Single control file for every hard-coded input / distribution (also pulled in
# transitively by montecarlo_utils.R; sourced explicitly here for clarity).
.all_params <- "all_parameters.R"
if (length(.sp) == 1 && file.exists(file.path(dirname(.sp), .all_params)))
  .all_params <- file.path(dirname(.sp), .all_params)
source(.all_params)

# logical
overwrite_data = ifelse(is.null(opts[["w"]]), T, as.logical(opts["w"]) )
main_scenario = ifelse(is.null(opts[["base"]]), F, as.logical(opts["base"]) )
diagnostics = ifelse(is.null(opts[["diagnostics"]]), F, as.logical(opts["diagnostics"]) )

# numeric
n_scenarios = ifelse(is.null(opts[["n"]]), N_SCENARIOS_DEFAULT, as.numeric(opts["n"]) )
seed = ifelse(is.null(opts[["seed"]]), DEFAULT_SEED, as.integer(opts["seed"]) )
if (length(n_scenarios) != 1 || is.na(n_scenarios) || n_scenarios < 1 || n_scenarios != floor(n_scenarios)) {
  stop("-n must be a positive integer")
}
n_scenarios <- as.integer(n_scenarios)
if (length(seed) != 1 || is.na(seed)) stop("--seed must be an integer")

# strings
res = ifelse(is.null(opts[["o"]]), "Montecarlo", as.character(opts["o"]) )
sampling_method = ifelse(is.null(opts[["method"]]), SAMPLING_METHOD, as.character(opts["method"]) )
if (!sampling_method %in% c("sobol","montecarlo")) stop("--method must be 'sobol' or 'montecarlo'")

# theta for the deterministic main scenario (--base=T): a fixed angle, or "na"
# to keep theta uncertain (sampled). Ignored when --base=F.
angle_opt = ifelse(is.null(opts[["angle"]]), as.character(MAIN_THETA), as.character(opts["angle"]))
keep_theta_uncertain <- str_to_lower(str_trim(angle_opt)) == "na"
if (!keep_theta_uncertain) {
  base_theta <- suppressWarnings(as.numeric(angle_opt))
  if (length(base_theta) != 1 || is.na(base_theta)) stop("--angle must be a number or 'na'")
  if (base_theta < THETA_TRUNC["lo"] || base_theta > THETA_TRUNC["hi"])
    stop("--angle must be between ", THETA_TRUNC["lo"], " and ", THETA_TRUNC["hi"], ", or 'na'")
} else {
  base_theta <- NA_real_
}

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res)
}

cat("Generating data... \n")

set.seed(seed)

if (overwrite_data==T | !file.exists(paste0(res,"/id_montecarlo.csv"))) {
    data <- data.frame()
  } else {
    data <- mc_strip_csv_index(as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv"), stringsAsFactors = FALSE)))
    validate_id_montecarlo(data, paste0(res, "/id_montecarlo.csv"))
    write.csv(data, file=paste0(res,"/id_montecarlo_copy.csv"), row.names = FALSE) # store a copy of previous version database
  }

max_id <- if (nrow(data) == 0) 0L else max(as.integer(data$ID), na.rm = TRUE)
start_draw_index <- if (nrow(data) == 0) 1L else max(as.integer(data$draw_index), na.rm = TRUE) + 1L
end_draw_index <- start_draw_index + as.integer(n_scenarios) - 1L

# Joint lognormal sampling for climate parameters
# as in FAIR v1.3
# Target: Corr(ECS, TCR) = 0.81 on the original (non-log) scale
ecs_par <- fit_distribution(median=ECS_DIST$median, q5=ECS_DIST$q5, q95=ECS_DIST$q95, return_params=T)
tcr_par <- fit_distribution(median=TCR_DIST$median, q5=TCR_DIST$q5, q95=TCR_DIST$q95, return_params=T)
rho_ecs_tcr <- ECS_TCR_CORR
rho_log <- log(1 + rho_ecs_tcr * sqrt((exp(ecs_par$sigma^2) - 1) * (exp(tcr_par$sigma^2) - 1))) /
  (ecs_par$sigma * tcr_par$sigma)

if (!is.finite(rho_log) || abs(rho_log) > 1) {
  stop("Inconsistent ECS/TCR marginals and requested correlation (", rho_ecs_tcr, ").")
}

cov_matrix <- matrix(c(1, rho_log, rho_log, 1), nrow=2, byrow=TRUE)
chol_cov <- chol(cov_matrix)
n_draws <- n_scenarios

# Discrete grids (equally-weighted) sampled jointly with the climate parameters.
rcp_choices        <- RCP_CHOICES
pulse_choices      <- PULSE_CHOICES
cool_choices       <- COOL_CHOICES
term_choices       <- TERM_CHOICES
start_choices      <- START_CHOICES

# Termination is modelled as a constant per-year hazard h = prob:
#   P(termination at year t | not yet terminated) = h.
# h is drawn log-uniform on [hazard_lo, hazard_hi]; this range is chosen so the
# implied median termination delay, ln(0.5)/ln(1-h), is long enough to span the
# analysed run (~6930 yr at h=1e-4 down to ~7 yr at h=1e-1).
hazard_lo <- HAZARD_LO
hazard_hi <- HAZARD_HI

# FAIR analysis horizon, in model-year units (t=1 is calendar 2020). FAIR solves
# t=1..1000, but post-processing keeps only t<=480 (the projection window
# t_proj=2020..2500; see sanitize_dt in Analyze_montecarlo.R). The termination
# time term_delta is the model-year SRM is switched off (t.val gt term_delta -> 0
# in experiments/srm.gms); censoring it at t_horizon means a termination in the
# last analysed year leaves no post-termination period, i.e. "no termination
# within the run". Keep this in sync with the t<=T_HORIZON cut in Analyze.
t_horizon <- T_HORIZON

# Map uniform (0,1) draws onto an equally-weighted discrete grid (inverse CDF).
map_discrete <- function(u, choices) {
  k <- length(choices)
  choices[pmin(floor(u * k) + 1L, k)]
}


# Each realization is a point in a 19-dimensional unit hypercube:
#  FAIR parameters (drive a FAIR run, encoded in the gdx file name):
#   1: ecs, 2: tcr (-> correlated joint lognormal)
#   3: rcp, 4: pulse, 5: cool, 6: term, 7: start, 8: term_delta (geometric; hazard from col 12)
#  Post-processing parameters (applied to the gdx output, not part of any run):
#   9: theta, 10: alpha, 11: delta, 12: prob (per-year termination hazard), 13: mortality_srm, 14: forctoTg,
#   15: TgtoUSD, 16: mortality_ozone, 17: vsl, 18: vsl_eta, 19: dg
# (the column->parameter mapping is SOBOL_COLUMN_MAP in all_parameters.R)
n_dim <- N_SOBOL_DIM

if (sampling_method == "sobol") {
  if (!mc_is_power_of_two(n_scenarios)) {
    warning("Sobol sample size is not a power of two; balance guarantees are strongest at powers of two.", call. = FALSE)
  }
  # Owen-scrambled Sobol sequence (randomized QMC). qrng delegates Owen
  # scrambling to spacefillr, so both packages are required. Scrambling keeps
  # the low-discrepancy structure while randomizing the net: `seed` therefore
  # selects an unbiased design (repeated seeds give Monte-Carlo error bars on
  # the QMC estimate) and every coordinate lands strictly inside (0,1).
  for (pkg in c("spacefillr","qrng")) {
    if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="http://cloud.r-project.org")
  }
  unit_draws_all <- qrng::sobol(n=end_draw_index, d=n_dim, randomize="Owen", seed=seed)
  unit_draws <- unit_draws_all[seq.int(start_draw_index, end_draw_index), , drop = FALSE]
} else if (sampling_method == "montecarlo") {
  # Plain Monte Carlo: independent pseudo-random uniforms. Generate through
  # end_draw_index and keep only the append tail so repeated appends with the
  # same seed never reuse earlier uniforms.
  unit_draws_all <- matrix(runif(end_draw_index * n_dim), ncol=n_dim)
  unit_draws <- unit_draws_all[seq.int(start_draw_index, end_draw_index), , drop = FALSE]
} else {
  stop("--method must be 'sobol' or 'montecarlo'")
}
unit_draw_indices <- seq.int(start_draw_index, end_draw_index)

# Transform uniforms into the correlated joint lognormal (ecs, tcr) marginals.
draw_joint_ecs_tcr <- function(u_mat) {
  u_mat <- pmin(pmax(u_mat, 1e-10), 1 - 1e-10)  # guard qnorm() against 0/1
  z_joint <- qnorm(u_mat) %*% chol_cov   # independent normals -> correlated
  list(
    ecs = exp(ecs_par$mu + ecs_par$sigma * z_joint[,1]),
    tcr = exp(tcr_par$mu + tcr_par$sigma * z_joint[,2])
  )
}

draw_replacement_ecs_tcr <- function(draw_indices, iter) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  redraw_u <- matrix(NA_real_, nrow = length(draw_indices), ncol = 2)
  for (j in seq_along(draw_indices)) {
    redraw_seed <- as.integer((as.numeric(seed) +
      as.numeric(draw_indices[j]) * 10007 + iter * 104729) %% .Machine$integer.max)
    set.seed(redraw_seed)
    redraw_u[j, ] <- runif(2)
  }
  draw_joint_ecs_tcr(redraw_u)
}

joint_draw <- draw_joint_ecs_tcr(unit_draws[, 1:2, drop=FALSE])
ecs_draw <- joint_draw$ecs
tcr_draw <- joint_draw$tcr

# Enforce ecs > tcr for every stored realization (after rounding to tenths).
# The few rejected draws are refreshed with pseudo-random joint normals: this is
# negligible relative to the quasi-random bulk and keeps the constraint exact.
invalid <- which(round(ecs_draw * 10, 0) <= round(tcr_draw * 10, 0))
iter <- 0
while (length(invalid) > 0) {
  iter <- iter + 1
  if (iter > 10000) {
    stop("Unable to obtain n_draws with ecs > tcr after repeated redraws.")
  }
  redraw <- draw_replacement_ecs_tcr(unit_draw_indices[invalid], iter)
  ecs_draw[invalid] <- redraw$ecs
  tcr_draw[invalid] <- redraw$tcr
  invalid <- which(round(ecs_draw * 10, 0) <= round(tcr_draw * 10, 0))
}

# Per-year termination hazard (col 12) and the termination delay it induces
# (col 8). Given hazard h, the delay D since deployment is geometric on yearly
# support {1, 2, ...}; use the discrete inverse CDF, not a rounded continuous
# approximation. prob stores the same raw hazard used to draw term_delay.
hazard <- exp(qunif(unit_draws[,12], log(hazard_lo), log(hazard_hi)))
term_delay <- mc_geometric_delay(unit_draws[,8], hazard)

# Absolute termination time (model-year index) = deployment year + delay,
# censored at the FAIR horizon. term_delta == t_horizon is the "no termination
# within the run" outcome. pulse is drawn here too so term_delta can fold it in.
pulse_draw <- map_discrete(unit_draws[,4], pulse_choices)
term_delta_draw <- mc_term_delta(pulse_draw, term_delay, t_horizon)

new_data <- data.frame(
  # --- FAIR-run parameters ---
  ecs = round(ecs_draw * 10, 0),
  tcr = round(tcr_draw * 10, 0),
  rcp = map_discrete(unit_draws[,3], rcp_choices),
  pulse = pulse_draw,
  cool = map_discrete(unit_draws[,5], cool_choices),
  term = map_discrete(unit_draws[,6], term_choices),
  start = map_discrete(unit_draws[,7], start_choices),
  term_delta = term_delta_draw,   # absolute termination year (deployment + geometric delay, censored at horizon)
  term_delay = term_delay,
  # --- post-processing parameters (all distributions live in all_parameters.R) ---
  theta           = round(qlnorm_fit_trunc(unit_draws[,9], lo=THETA_TRUNC[["lo"]], hi=THETA_TRUNC[["hi"]],
                          median=THETA_DIST$median, q5=THETA_DIST$q5, q95=THETA_DIST$q95), THETA_ROUND),
  alpha           = round(qlnorm_fit(unit_draws[,10], median=ALPHA_DIST$median, sd=ALPHA_DIST$sd), ALPHA_ROUND),
  delta           = round(qunif(unit_draws[,11], DELTA_UNIF[["min"]], DELTA_UNIF[["max"]]), DELTA_ROUND),  # widened +/- half a step so grid edges get a full bin
  prob            = hazard,   # per-year termination hazard (log-uniform), drives term_delta
  mortality_srm   = round(qlnorm_fit(unit_draws[,13], median=MORTALITY_SRM_DIST$median,
                          q5=MORTALITY_SRM_DIST$q5, q95=MORTALITY_SRM_DIST$q95), MORTALITY_SRM_ROUND),  # lognormal>0: no clamp needed
  forctoTg        = round(1 / qunif(unit_draws[,14], FORCTOTG_INV_UNIF[["min"]], FORCTOTG_INV_UNIF[["max"]]), FORCTOTG_ROUND),
  TgtoUSD         = round(qunif(unit_draws[,15], TGTOUSD_UNIF[["min"]], TGTOUSD_UNIF[["max"]]), TGTOUSD_ROUND),
  mortality_ozone = round(qnorm_fit_trunc(unit_draws[,16], lo=MORTALITY_OZONE_TRUNC[["lo"]], median=MORTALITY_OZONE_DIST$median,
                          q5=MORTALITY_OZONE_DIST$q5, q95=MORTALITY_OZONE_DIST$q95) * MORTALITY_OZONE_SCALE, MORTALITY_OZONE_ROUND),
  vsl             = round(qunif(unit_draws[,17], VSL_UNIF_MILLIONS[["min"]], VSL_UNIF_MILLIONS[["max"]]), VSL_ROUND) * VSL_SCALE,
  vsl_eta         = round(qunif(unit_draws[,18], VSL_ETA_UNIF[["min"]], VSL_ETA_UNIF[["max"]]), VSL_ETA_ROUND),  # widened +/- half a step so grid edges get a full bin
  dg              = round(qnorm_fit(unit_draws[,19], median=DG_DIST$median, sd=DG_DIST$sd), DG_ROUND),
  sampler_version = MC_SAMPLER_VERSION,
  sampling_method = sampling_method,
  seed            = seed,
  censored_termination = mc_censored_termination(pulse_draw, term_delay, t_horizon),
  stringsAsFactors = FALSE
)

new_data <- new_data %>%
    mutate(term=ifelse(cool==0,COOL0_TERM_SENTINEL,term),
           start=ifelse(cool==0,COOL0_START_SENTINEL,start) )

if (main_scenario==T) {
    new_data <- rbind(new_data %>%
                    mutate(rcp=MAIN_RCP, cool=MAIN_COOL, term=MAIN_TERM, start=MAIN_START, pulse=MAIN_PULSE_SET[1] ),
                  new_data %>%
                    mutate(rcp=MAIN_RCP, cool=MAIN_COOL, term=MAIN_TERM, start=MAIN_START, pulse=MAIN_PULSE_SET[2] ) )
    # Fix the post-processing parameters for the deterministic main scenario
    # (previously applied in Analyze_montecarlo.R).
    new_data <- new_data %>% mutate(vsl = MAIN_VSL, delta = MAIN_DELTA, vsl_eta = MAIN_VSL_ETA)
    if (!keep_theta_uncertain) new_data <- new_data %>% mutate(theta = base_theta) }

new_data <- new_data %>%
  mutate(
    term_delta = mc_term_delta(pulse, term_delay, t_horizon),
    censored_termination = mc_censored_termination(pulse, term_delay, t_horizon),
    ID = max_id + seq_len(n()),
    draw_index = start_draw_index + seq_len(n()) - 1L
  ) %>%
  select(ID, all_of(MC_FAIR_COLS), all_of(MC_POST_COLS), all_of(MC_METADATA_COLS))

validate_id_montecarlo(new_data, "new Monte Carlo realizations")
data <- data %>% bind_rows(new_data)
validate_id_montecarlo(data, paste0(res, "/id_montecarlo.csv"))
  
write.csv(data,file=paste0(res,"/id_montecarlo.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Optional diagnostics: compare the freshly drawn marginals (new_data) against
# the theoretical distribution each parameter targets. The theoretical curve
# folds in the transforms the sampler actually applies -- truncation (theta,
# ozone, via the truncated inverse-CDF) and grid binning (coarse rounded params)
# -- so the only residual, annotated deviations are integer rounding and the
# joint ecs>tcr rejection (which no single marginal can capture).
# ---------------------------------------------------------------------------
if (diagnostics == T) {
  cat("Writing distribution diagnostics... \n")
  pdf_path <- paste0(res, "/diagnostics_distributions.pdf")
  pdf(pdf_path, width = 9, height = 11)

  # fitted (mu,sigma); warnings already surfaced during sampling, so muffle here
  ln <- function(...) suppressWarnings(fit_distribution("lognormal", ..., return_params = TRUE))
  nm <- function(...) suppressWarnings(fit_distribution("normal",    ..., return_params = TRUE))

  # one parameter -> a density overlay + a theoretical Q-Q plot
  diag_panel <- function(x, dfun, qfun, title, markers = NULL, note = NULL) {
    x <- x[is.finite(x)]
    xs <- seq(min(x), max(x), length.out = 512)
    yt <- dfun(xs)
    h  <- hist(x, breaks = 40, plot = FALSE)
    hist(x, breaks = 40, freq = FALSE, col = "grey85", border = "grey60",
         main = title, xlab = "value",
         ylim = c(0, max(c(h$density, yt), na.rm = TRUE)))
    lines(xs, yt, col = "red", lwd = 2)
    if (!is.null(markers)) abline(v = markers, col = "blue", lty = 2)
    legend("topright", c("drawn","theoretical","target q"),
           col = c("grey60","red","blue"), lwd = c(6,2,1), lty = c(1,1,2),
           bty = "n", cex = 0.7)
    summ <- sprintf("emp med/q5/q95 = %.4g / %.4g / %.4g | theo = %.4g / %.4g / %.4g",
                    median(x), quantile(x, .05), quantile(x, .95),
                    qfun(.5), qfun(.05), qfun(.95))
    mtext(summ, side = 1, line = 2.6, cex = 0.5)
    if (!is.null(note)) mtext(note, side = 3, line = -1, cex = 0.55, col = "darkgreen")
    # theoretical Q-Q
    pp  <- ppoints(min(length(x), 2000))
    plot(qfun(pp), quantile(x, pp, names = FALSE), pch = ".", cex = 2, col = "grey30",
         main = paste("Q-Q:", title), xlab = "theoretical quantile", ylab = "empirical quantile")
    abline(0, 1, col = "red", lwd = 1.5)
  }

  # truncated density / quantile on [lo,hi] from a base d/p/q (matches the
  # truncated inverse-CDF sampler: redistributes tail mass, no boundary spike)
  dtrunc <- function(dfun, pfun, lo, hi) function(x) ifelse(x >= lo & x <= hi, dfun(x) / (pfun(hi) - pfun(lo)), 0)
  qtrunc <- function(pfun, qfun, lo, hi) function(pp) qfun(pfun(lo) + pp * (pfun(hi) - pfun(lo)))

  # discrete-aware panel for coarsely-rounded (grid) parameters: observed vs the
  # theoretical mass binned onto the grid, P(round(X)=k) = F(k+step/2)-F(k-step/2)
  diag_discrete <- function(x, cdf, step, title, note = NULL) {
    vals <- sort(unique(x))
    obs  <- as.numeric(table(factor(x, levels = vals))) / length(x)
    exq  <- pmax(0, cdf(vals + step / 2) - cdf(vals - step / 2)); exq <- exq / sum(exq)
    ymax <- max(obs, exq)
    plot(vals, obs, type = "h", lwd = 6, col = "grey70", ylim = c(0, ymax),
         main = title, xlab = "grid value", ylab = "probability")
    points(vals, exq, col = "red", pch = 19, cex = 0.9)
    legend("topright", c("observed", "expected (binned theory)"),
           col = c("grey70", "red"), pch = c(NA, 19), lwd = c(6, NA), bty = "n", cex = 0.7)
    if (!is.null(note)) mtext(note, side = 3, line = -1, cex = 0.55, col = "darkgreen")
    plot(exq, obs, pch = 19, col = "grey30", xlim = c(0, ymax), ylim = c(0, ymax),
         main = paste("obs vs exp:", title), xlab = "expected prob", ylab = "observed prob")
    abline(0, 1, col = "red", lwd = 1.5)
  }

  # --- joint ecs/tcr panel (correlation + constraint) ---
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  ex <- new_data$ecs / 10; tx <- new_data$tcr / 10
  smoothScatter(ex, tx, xlab = "ECS", ylab = "TCR",
                main = sprintf("Joint ECS-TCR  (empirical corr = %.3f, target %.2f)", cor(ex, tx), ECS_TCR_CORR))
  abline(0, 1, col = "red", lty = 2)
  legend("bottomright", "ECS = TCR (constraint)", col = "red", lty = 2, bty = "n", cex = 0.8)

  # --- per-parameter marginals (3 params = 3 rows of [density | Q-Q] per page) ---
  par(mfrow = c(3, 2), mar = c(4.5, 4, 3, 1))

  p <- do.call(ln, ECS_DIST)
  diag_panel(ex, function(x) dlnorm(x, p$mu, p$sigma), function(q) qlnorm(q, p$mu, p$sigma),
             "ecs (lognormal)", markers = c(ECS_DIST$q5, ECS_DIST$median, ECS_DIST$q95),
             note = "drawn rounded to 0.1 & constrained ECS>TCR: lower-tail deviation expected")
  p <- do.call(ln, TCR_DIST)
  diag_panel(tx, function(x) dlnorm(x, p$mu, p$sigma), function(q) qlnorm(q, p$mu, p$sigma),
             "tcr (lognormal)", markers = c(TCR_DIST$q5, TCR_DIST$median, TCR_DIST$q95),
             note = "drawn rounded to 0.1 & constrained ECS>TCR: upper-tail deviation expected")
  p <- do.call(ln, THETA_DIST)
  diag_panel(new_data$theta,
             dtrunc(function(x) dlnorm(x, p$mu, p$sigma), function(x) plnorm(x, p$mu, p$sigma), THETA_TRUNC[["lo"]], THETA_TRUNC[["hi"]]),
             qtrunc(function(x) plnorm(x, p$mu, p$sigma), function(q) qlnorm(q, p$mu, p$sigma), THETA_TRUNC[["lo"]], THETA_TRUNC[["hi"]]),
             "theta (lognormal truncated 0-90)", markers = c(THETA_DIST$q5, THETA_DIST$median, THETA_DIST$q95),
             note = "truncated inverse-CDF on (0,90]; drawn rounded to 1")
  p <- do.call(ln, ALPHA_DIST)
  diag_panel(new_data$alpha, function(x) dlnorm(x, p$mu, p$sigma), function(q) qlnorm(q, p$mu, p$sigma),
             "alpha (lognormal, median+sd)")
  p <- do.call(ln, MORTALITY_SRM_DIST)
  diag_panel(new_data$mortality_srm, function(x) dlnorm(x, p$mu, p$sigma), function(q) qlnorm(q, p$mu, p$sigma),
             "mortality_srm (lognormal)",
             markers = c(MORTALITY_SRM_DIST$q5, MORTALITY_SRM_DIST$median, MORTALITY_SRM_DIST$q95),
             note = "inputs inconsistent (q5*q95 != median^2): compromise fit")
  p <- do.call(nm, MORTALITY_OZONE_DIST)
  mo <- p$mu * MORTALITY_OZONE_SCALE; so <- p$sigma * MORTALITY_OZONE_SCALE
  diag_panel(new_data$mortality_ozone,
             dtrunc(function(x) dnorm(x, mo, so), function(x) pnorm(x, mo, so), MORTALITY_OZONE_TRUNC[["lo"]], MORTALITY_OZONE_TRUNC[["hi"]]),
             qtrunc(function(x) pnorm(x, mo, so), function(q) qnorm(q, mo, so), MORTALITY_OZONE_TRUNC[["lo"]], MORTALITY_OZONE_TRUNC[["hi"]]),
             "mortality_ozone (normal /100, truncated >=0)",
             markers = c(MORTALITY_OZONE_DIST$q5, MORTALITY_OZONE_DIST$median, MORTALITY_OZONE_DIST$q95) * MORTALITY_OZONE_SCALE,
             note = "truncated inverse-CDF at 0, scaled /100, rounded to 1")
  p <- do.call(nm, DG_DIST)
  diag_panel(new_data$dg, function(x) dnorm(x, p$mu, p$sigma), function(q) qnorm(q, p$mu, p$sigma),
             "dg (normal, median+sd)")

  # fine-grid uniform parameters (rounding far below the spread -> continuous overlay)
  diag_panel(new_data$TgtoUSD,
             function(x) dunif(x, TGTOUSD_UNIF[["min"]], TGTOUSD_UNIF[["max"]]),
             function(q) qunif(q, TGTOUSD_UNIF[["min"]], TGTOUSD_UNIF[["max"]]),
             "TgtoUSD ~ U(0.75,3)", note = "drawn rounded to 0.01 (effectively continuous)")
  .ftg_lo <- FORCTOTG_INV_UNIF[["min"]]; .ftg_hi <- FORCTOTG_INV_UNIF[["max"]]
  diag_panel(new_data$forctoTg,
             function(x) 1 / ((.ftg_hi - .ftg_lo) * x^2),
             function(q) 1 / (.ftg_hi - q * (.ftg_hi - .ftg_lo)),
             "forctoTg = 1/U(0.2,1.5)", note = "reciprocal-uniform; drawn rounded to 0.01")

  # prob = per-year termination hazard ~ LogUniform(hazard_lo,hazard_hi): a
  # log-uniform shows up as a flat density on the log10 scale, so check it there.
  diag_panel(log10(new_data$prob),
             function(x) dunif(x, log10(hazard_lo), log10(hazard_hi)),
             function(q) qunif(q, log10(hazard_lo), log10(hazard_hi)),
             sprintf("log10(prob): hazard ~ LogUniform(%g,%g)", hazard_lo, hazard_hi),
             note = "per-year termination hazard; drives term_delta, drawn to 4 sig figs")

  # Termination delay (years of SRM before termination), the quantity the hazard
  # acts on directly: term_delta itself is the absolute year (pulse + delay,
  # censored at the horizon), so we validate the underlying delay here. Each
  # realization's delay is geometric with that realization's own hazard h, and h
  # is log-uniform, so the marginal is the log-uniform mixture of geometrics.
  # The delay is shown censored at t_horizon to mirror the run.
  .hgrid   <- exp(seq(log(hazard_lo), log(hazard_hi), length.out = 512))
  d_delay  <- function(d) {
    vapply(d, function(z) {
      k <- pmax(as.integer(round(z)), 1L)
      mean(.hgrid * (1 - .hgrid)^(k - 1L))
    }, numeric(1))
  }
  p_delay  <- function(d) vapply(d, function(z) mean(1 - (1 - .hgrid)^floor(z)), numeric(1))
  .dg <- seq_len(t_horizon)
  .pg <- p_delay(.dg)
  q_delay  <- function(pp) approx(.pg, .dg, xout = pp, rule = 2)$y
  delay_plot <- pmin(new_data$term_delay, t_horizon)
  diag_panel(delay_plot, d_delay, q_delay,
             "termination delay (LogUnif-hazard mixture of geometrics)",
             note = sprintf("%.1f%% censored at t_horizon=%d (no termination in run): hist spike & Q-Q flatten there",
                            100 * mean(delay_plot == t_horizon), t_horizon))

  # coarse-grid uniform parameters (few distinct values -> discrete-aware panel).
  # Sampling ranges widened by +/- half a grid step so every grid value (incl.
  # the edges) gets an equal-width bin; the expected-mass cdf uses the same
  # widened bounds so expected stays flat across all grid points.
  diag_discrete(new_data$delta,    function(x) punif(x, DELTA_UNIF[["min"]], DELTA_UNIF[["max"]]), DELTA_GRID_STEP, "delta ~ U(0.01,0.07) [grid 0.01, edges equalized]")
  diag_discrete(new_data$vsl_eta,  function(x) punif(x, VSL_ETA_UNIF[["min"]], VSL_ETA_UNIF[["max"]]), VSL_ETA_GRID_STEP,  "vsl_eta ~ U(0.4,1) [grid 0.1, edges equalized]")
  diag_discrete(new_data$vsl / VSL_SCALE, function(x) punif(x, VSL_UNIF_MILLIONS[["min"]], VSL_UNIF_MILLIONS[["max"]]), 1, "vsl/1e6 ~ round(U(7.5,13.6)) [grid 1]",
                note = "banker's rounding at .5 may shift edge bins slightly")

  dev.off()
  cat("Diagnostics written to", pdf_path, "\n")
}
