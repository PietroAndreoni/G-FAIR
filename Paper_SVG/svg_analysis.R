# =============================================================================
# svg_analysis.R
#
# Social value of SAI (SVG) from the `svg` GAMS experiment (experiments/svg.gms).
#
# That experiment now writes ONE gdx per scenario configuration, holding a single
# report(vrn,tc,rep,t) parameter plus the exogenous paths and calibration scalars.
#   vrn = experiment variant     tc = commitment length / pulse delay [years]
#   rep = reported quantity      t  = model period (t = 1 is 2020)
#
#   base     (tc 0)   background SAI schedule, no emission pulse. Its temperature
#                     path is the `target_temp` the other runs are held to.
#   ghgpulse (tc 0)   same SAI schedule + a one-year CO2 pulse at `pulse_time`,
#                     SAI NOT adjusted.
#   srmpulse (tc 0)   the *extra* SAI that exactly offsets that pulse, with
#                     emissions back at their baseline level: the baseline world
#                     plus a marginal amount of SAI and nothing else. Never
#                     terminated.
#   term     (tc > 0) that same marginal SAI, switched off tc years after the
#                     pulse -- the actual value of a finite commitment, rebound
#                     included.
#   scc      (tc > 0) background SAI only, CO2 pulse emitted tc years later:
#                     the SCC schedule along the deployment.
#
# With damages D(t) [trillion USD] and discount factor 1/(1+r)^(t-t0):
#
#   NPV_pulse = NPV[ D_ghgpulse - D_base ]    damage caused by the pulse
#   NPV_sai   = NPV[ D_srmpulse - D_base ]    damage caused by the SAI
#   MB_sai    = -NPV_sai                      marginal *benefit* of the SAI
#
#   SCC       = NPV_pulse / pulse size in tCO2          [USD/tCO2]
#   SVG       = MB_sai / NPV_pulse                      [-]  marginal benefit of
#                                                       the SAI normalised by the SCC
#   SVG_tCO2  = MB_sai / SCC = SVG * pulse_tCO2         [tCO2]
#
# Numerator and denominator refer to the SAME pulse, so the ratio is unit-free and
# independent of the pulse size and of the discounting reference year.
#
# Damage accounting (Model/srm_impacts.gms):
#   DAMFRAC_TEMP = a0*TATM^b0                     realised, SAI-cooled temperature
#   DAMFRAC_SRM  = a0*TATM_GHG^b0*(1-EFF_SRM)     imperfect-offsetting penalty,
#                = 2*a0*(Tecs/forc2x)^2*TOT_FORC*FORC_SRM*(1-cos(srm_angle)) [b0=2]
#   DAMAGES      = (DAMFRAC_TEMP + DAMFRAC_SRM)*y + VLL + COST_SRM
# The penalty is zero without SAI and at srm_angle = 0, and non-negative always.
# svg_damage_components() re-derives that identity from the report and checks it.
#
# ---------------------------------------------------------------------------
# Usage
#   Rscript Paper_SVG/svg_analysis.R                 # full run, writes ./output
#   source("Paper_SVG/svg_analysis.R")               # same, functions stay available
#   svg_report(svg_social_value(rate = 0.05))        # another discount rate
#   plot_gdx_diff("tatm")                            # generic difference plot
#   plot_gdx_diff("damages", scenarios = c("term50", "srmpulse"))
#   options(svg.no.main = TRUE); source(...)         # load the functions only
#
# Model time convention: t = 1 is calendar 2020 (calendar year = 2019 + t), as in
# Model/initialization.gms.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gdxtools)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
SVG_DISCOUNT_RATE <- 0.03    # default discount rate
SVG_T_MAX         <- 1000L   # last model period kept: the full FAIR horizon, which
                             # the longest pulse delays (tc = 500) need. Beyond
                             # ~2500 emissions and GDP are flat and, at r >= 1%,
                             # the tail is worth < 1e-6 of the total.
SVG_BASE_YEAR     <- 2019L   # calendar year = SVG_BASE_YEAR + t

# Scenario tags parsed out of the gdx file name, e.g.
#   RCP45_EXPsvg_GASco2_ECS30_TCR18_PT5_RC10_EC2400_BC2025.gdx
SVG_EXPERIMENT   <- "svg"
SVG_TAG_PATTERNS <- c(
  rcp        = "(?<=RCP)[A-Za-z0-9]+?(?=_)",
  gas        = "(?<=GAS)[a-z0-9]+",
  ecs        = "(?<=ECS)[0-9]+",
  tcr        = "(?<=TCR)[0-9]+",
  pulse_time = "(?<=PT)[0-9]+",
  cool_rate  = "(?<=RC)[0-9]+",
  geo_end    = "(?<=EC)[0-9]+",
  geo_start  = "(?<=BC)[0-9]+"
)
SVG_SCENARIO_KEYS <- c("rcp", "gas", "ecs", "tcr", "pulse_time", "cool_rate",
                       "geo_end", "geo_start")

# Variant labels, as written by experiments/svg.gms.
SVG_BASE     <- "base"
SVG_GHGPULSE <- "ghgpulse"
SVG_SRMPULSE <- "srmpulse"
SVG_TERM     <- "term"
SVG_SCC      <- "scc"


# -----------------------------------------------------------------------------
# Locations / GAMS
# -----------------------------------------------------------------------------

# The project root is the directory holding FAIR.gms. Resolved from the script
# path (Rscript --file / RStudio "Source") or from the working directory, so the
# script runs from anywhere.
svg_project_root <- function() {
  starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
              unlist(lapply(sys.frames(), function(f) f$ofile)), getwd())
  for (s in starts[nzchar(starts)]) {
    d <- if (dir.exists(s)) s else dirname(s)
    repeat {
      if (file.exists(file.path(d, "FAIR.gms"))) return(normalizePath(d, "/", FALSE))
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Cannot locate the project root (the folder holding FAIR.gms).", call. = FALSE)
}

# Directory this script lives in (Rscript --file, RStudio "Source", or the working
# directory as a fallback), so its outputs travel with it if the folder is moved.
svg_script_dir <- function() {
  p <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
         unlist(lapply(sys.frames(), function(f) f$ofile)))
  p <- p[nzchar(p)]
  if (length(p)) normalizePath(dirname(p[1]), "/", FALSE) else normalizePath(getwd(), "/", FALSE)
}

SVG_ROOT       <- svg_project_root()
SVG_RESULTS    <- file.path(SVG_ROOT, "Results")   # where svg.gms unloads its gdx
SVG_OUTPUT_DIR <- file.path(svg_script_dir(), "output")

# Point the GDX library at the GAMS installation. Once per session; `gams` on the
# PATH by default, override with the GAMS_DIR environment variable.
svg_init_gams <- function(gams_dir = Sys.getenv("GAMS_DIR", unset = "")) {
  if (isTRUE(getOption("svg.gams.ready"))) return(invisible(TRUE))
  if (!nzchar(gams_dir)) {
    exe <- Sys.which("gams")
    if (!nzchar(exe)) stop("GAMS not found on the PATH; set GAMS_DIR to the GAMS folder.", call. = FALSE)
    gams_dir <- dirname(exe)
  }
  # igdx() comes from gdxtools (which vendors the gdxrrw loader).
  invisible(utils::capture.output(gdxtools::igdx(gams_dir)))
  if (!gdxtools::igdx(silent = TRUE, returnStr = FALSE))
    stop("Could not load the GDX library from '", gams_dir, "'.", call. = FALSE)
  options(svg.gams.ready = TRUE)
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# Locating the run
# -----------------------------------------------------------------------------

# Index every experiment gdx of a results folder: experiment name, scenario tags,
# modification time.
svg_gdx_index <- function(folder = SVG_RESULTS) {
  files <- list.files(folder, pattern = "\\.gdx$", full.names = TRUE)
  files <- files[grepl("_EXP", basename(files), fixed = TRUE)]
  if (!length(files)) stop("No experiment gdx found in ", folder, call. = FALSE)
  idx <- data.table(path = normalizePath(files, "/", FALSE), file = basename(files))
  grab <- function(x, pattern) {
    m <- regexpr(pattern, x, perl = TRUE)
    out <- rep(NA_character_, length(x))
    out[m > 0] <- regmatches(x, m)
    out
  }
  idx[, experiment := sub("\\.gdx$", "", grab(file, "(?<=EXP)[^_]+"))]
  for (tag in names(SVG_TAG_PATTERNS)) idx[, (tag) := grab(file, SVG_TAG_PATTERNS[[tag]])]
  # An unresolved GAMS macro in the name (ECS%ecs%) means the run was launched
  # without that command-line argument: not a usable scenario.
  idx <- idx[!grepl("%", file, fixed = TRUE)]
  idx[, mtime := file.mtime(path)]
  idx[]
}

# Resolve the svg gdx to analyse.
#   filters: scenario tags to pin, e.g. list(geo_start = "2025", ecs = "30").
# The newest matching run wins; the alternatives are listed so a stale one never
# takes over silently. The returned path carries the scenario tags as attributes.
svg_run <- function(folder = SVG_RESULTS, filters = list(), quiet = FALSE) {
  idx <- svg_gdx_index(folder)[experiment == SVG_EXPERIMENT]
  if (!nrow(idx))
    stop("No '", SVG_EXPERIMENT, "' gdx in ", folder,
         ". Run: gams FAIR.gms --experiment=svg --ecs=.. --tcr=..", call. = FALSE)
  for (nm in names(filters)) {
    if (!nm %in% names(idx)) stop("Unknown scenario tag '", nm, "'.", call. = FALSE)
    idx <- idx[get(nm) == as.character(filters[[nm]])]
  }
  if (!nrow(idx)) stop("No svg gdx matches those filters.", call. = FALSE)
  setorder(idx, -mtime)
  if (nrow(idx) > 1 && !quiet) {
    fmt <- function(r) paste(sprintf("%s=%s", SVG_SCENARIO_KEYS, r), collapse = " ")
    message("svg_run(): ", nrow(idx), " svg runs found; using the newest\n  ",
            fmt(unlist(idx[1, ..SVG_SCENARIO_KEYS])), "\n  others (select with e.g. ",
            "filters = list(geo_start = \"2020\")):\n  ",
            paste(apply(idx[-1, ..SVG_SCENARIO_KEYS], 1, fmt), collapse = "\n  "))
  }
  structure(idx[1, path], names = SVG_EXPERIMENT,
            scenario = as.list(idx[1, ..SVG_SCENARIO_KEYS]),
            class = c("svg_run", "character"))
}

print.svg_run <- function(x, ...) {
  cat("svg run:", basename(unclass(x)), "\n  ",
      paste(sprintf("%s=%s", names(attr(x, "scenario")), unlist(attr(x, "scenario"))),
            collapse = " "), "\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# gdx reading
# -----------------------------------------------------------------------------

# Read one symbol from one or more gdx.
#   files : named character vector of gdx paths; the names become `source`
#   field : gdx field -- "l" (level, default), "m" (marginal), "lo", "up"
# Returns a data.table: source, <set dimensions...>, t (integer), value.
svg_read <- function(files, symbol, field = "l", t_max = SVG_T_MAX) {
  svg_init_gams()
  files <- stats::setNames(as.character(files), names(files) %||% basename(files))
  got <- gdxtools::batch_extract(symbol, files = unname(files), field = field)[[symbol]]
  if (is.null(got) || !nrow(got))
    stop("Symbol '", symbol, "' not found (or empty) in the gdx.", call. = FALSE)
  dt <- as.data.table(got)
  dims <- setdiff(names(dt), c("value", "gdx"))
  if (!length(dims)) stop("Symbol '", symbol, "' is a scalar: use svg_scalar().", call. = FALSE)
  # Every reported symbol carries time as its last index.
  setnames(dt, if ("t" %in% dims) "t" else dims[length(dims)], "t")
  dt[, t := as.integer(as.character(t))]
  dt[, source := names(files)[match(gdx, unname(files))]]
  dt[, gdx := NULL]
  if (is.finite(t_max)) dt <- dt[t <= t_max]
  setcolorder(dt, c("source", setdiff(names(dt), c("source", "t", "value")), "t", "value"))
  dt[]
}

# Read a 0-dimensional parameter (a0, b0, Tecs, CO2toC, ...).
svg_scalar <- function(file, symbol) {
  svg_init_gams()
  v <- gdxtools::batch_extract(symbol, files = unname(file)[1])[[symbol]]
  if (is.null(v) || !nrow(v)) stop("Scalar '", symbol, "' not found in ", file, call. = FALSE)
  as.numeric(v$value[1])
}

# The report parameter, one row per (variant, commitment, period) and one column
# per reported quantity, with the GDP path joined on.
#   scenario : "base", "ghgpulse", "srmpulse", "term10", "scc250", ...
# Records GAMS dropped for being exactly zero are restored as zeros.
svg_report_table <- function(run = NULL, t_max = SVG_T_MAX, folder = SVG_RESULTS,
                             filters = list(), quiet = FALSE) {
  if (is.null(run)) run <- svg_run(folder, filters, quiet = quiet)
  file <- stats::setNames(as.character(run), "svg")
  long <- tryCatch(svg_read(file, "report", t_max = t_max), error = function(e)
    stop("No `report` parameter in ", basename(file),
         ": that gdx predates the single-unload svg.gms. Re-run svg.gms.", call. = FALSE))
  long[, `:=`(tc = as.integer(as.character(tc)), source = NULL)]
  w <- dcast(long, vrn + tc + t ~ rep, value.var = "value", fill = 0)
  y <- svg_read(file, "y", t_max = t_max)[, .(t, y = value)]
  w <- merge(w, y, by = "t")
  w[, scenario := fifelse(tc == 0L, vrn, paste0(vrn, tc))]
  setcolorder(w, c("scenario", "vrn", "tc", "t"))
  setorder(w, vrn, tc, t)
  setattr(w, "run", run)
  w[]
}

# Quantities the report carries, i.e. what gdx_diff() / plot_gdx_diff() can draw.
svg_quantities <- function(rep) setdiff(names(rep), c("scenario", "vrn", "tc", "t", "y"))


# -----------------------------------------------------------------------------
# Damage accounting
# -----------------------------------------------------------------------------

# Damage path of every scenario, by component, in trillion USD per year.
#   temp    : DAMFRAC_TEMP * y -- damage of the realised, SAI-cooled temperature
#   penalty : DAMFRAC_SRM * y  -- imperfect-offsetting damage
#   vll     : VLL              -- value of lives lost to SAI aerosol
#   cost    : COST_SRM         -- SAI deployment cost
#   total   : their sum, which must reproduce the reported DAMAGES (eq_damtot)
svg_damage_components <- function(rep) {
  d <- copy(rep)
  d[, `:=`(temp = damfrac_temp * y, penalty = damfrac_srm * y)]
  d[, total := temp + penalty + vll + cost]
  err <- d[, max(abs(damages - total))]
  if (!is.finite(err) || err > 1e-6 * max(1, d[, max(abs(damages))]))
    stop("Damage components do not reproduce DAMAGES from the gdx (max abs error ",
         signif(err, 3), "): check eq_damtot in Model/srm_impacts.gms.", call. = FALSE)
  d[]
}

# Size and timing of the emission pulse carried by each pulsed variant, read off
# the reported CO2 emissions rather than assumed from the file name.
# W_EMI('co2',.) is in GtCO2 (initialization.gms divides the RCP GtC by CO2toC).
svg_pulses <- function(rep) {
  base <- rep[vrn == SVG_BASE, .(t, base_emi = wemi_co2)]
  d <- merge(rep[vrn != SVG_BASE], base, by = "t")
  d[, demi := wemi_co2 - base_emi]
  out <- d[abs(demi) > 1e-9, .(pulse_t = min(t), pulse_GtCO2 = sum(demi)), by = .(vrn, tc, scenario)]
  out[, `:=`(pulse_year = SVG_BASE_YEAR + pulse_t, pulse_tons = pulse_GtCO2 * 1e9)]
  out[]
}


# -----------------------------------------------------------------------------
# Generic: difference of any reported quantity against a baseline scenario
# -----------------------------------------------------------------------------

# Difference of `variable` between each scenario and the baseline scenario.
#   variable  : a reported quantity, see svg_quantities() -- "tatm", "damages",
#               "forc_srm", ... (names are case-insensitive, so "TATM" works)
#   baseline  : scenario label to difference against
#   relative  : TRUE -> (x - base)/|base| instead of (x - base)
gdx_diff <- function(variable, rep = NULL, baseline = SVG_BASE, scenarios = NULL,
                     relative = FALSE, ...) {
  if (is.null(rep)) rep <- svg_report_table(...)
  v <- tolower(variable)
  if (!v %in% svg_quantities(rep))
    stop("'", variable, "' is not reported. Available: ",
         paste(svg_quantities(rep), collapse = ", "),
         ". Add it to `rep` in experiments/svg.gms and to experiments/svg_store.gms.",
         call. = FALSE)
  if (!baseline %in% rep$scenario)
    stop("Baseline '", baseline, "' is not one of the scenarios.", call. = FALSE)
  b <- rep[scenario == baseline, .(t, base = get(v))]
  d <- rep[scenario != baseline, .(scenario, vrn, tc, t, value = get(v))]
  if (!is.null(scenarios)) {
    unknown <- setdiff(scenarios, unique(rep$scenario))
    if (length(unknown)) stop("No such scenario(s): ", paste(unknown, collapse = ", "),
                              ". Available: ", paste(unique(rep$scenario), collapse = ", "),
                              call. = FALSE)
    d <- d[scenario %in% scenarios]
  }
  out <- merge(d, b, by = "t")
  out[, diff := value - base]
  if (relative) out[, diff := fifelse(base == 0, NA_real_, diff / abs(base))]
  out[, year := SVG_BASE_YEAR + t]
  setattr(out, "variable", v)
  setattr(out, "baseline", baseline)
  setorder(out, scenario, t)[]
}

# Plot the difference of `variable` against the baseline scenario; returns a ggplot.
# With no `scenarios`, the reference runs (ghgpulse, srmpulse) are drawn -- the
# commitment families would otherwise put 36 lines on one panel.
plot_gdx_diff <- function(variable, rep = NULL, baseline = SVG_BASE, scenarios = NULL,
                          relative = FALSE, years = TRUE, t_max = NULL,
                          title = NULL, ylab = NULL, ...) {
  if (is.null(rep)) rep <- svg_report_table(...)
  if (is.null(scenarios)) scenarios <- setdiff(rep[tc == 0L, unique(scenario)], baseline)
  d <- gdx_diff(variable, rep = rep, baseline = baseline, scenarios = scenarios,
                relative = relative)
  if (!is.null(t_max)) d <- d[t <= t_max]
  ggplot(d, aes(x = if (years) year else t, y = diff, colour = scenario)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line(linewidth = 0.8) +
    labs(x = if (years) "Year" else "Model period t",
         y = ylab %||% paste0(if (relative) "Relative difference, " else "Difference, ",
                              attr(d, "variable")),
         colour = NULL,
         title = title %||% paste0(attr(d, "variable"), ": deviation from ", baseline)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}


# -----------------------------------------------------------------------------
# Social value of SAI
# -----------------------------------------------------------------------------

# NPV of a scenario's damage difference from the baseline, discounted to `t0` and
# summed from `from` to the end of the horizon. Returns one row per scenario.
svg_npv <- function(dam, scenarios, rate, t0, from = t0, by_component = FALSE) {
  cols <- c("temp", "penalty", "vll", "cost", "total")
  base <- dam[vrn == SVG_BASE, c("t", cols), with = FALSE]
  setnames(base, cols, paste0("b_", cols))
  d <- merge(dam[scenario %in% scenarios], base, by = "t")
  for (cc in cols) d[, (cc) := get(cc) - get(paste0("b_", cc))]
  d <- d[t >= from]
  d[, discount := (1 + rate)^-(t - t0)]
  out <- d[, lapply(.SD, function(v) sum(v * discount)), by = .(scenario, vrn, tc),
           .SDcols = cols]
  if (by_component) out else out[, .(scenario, vrn, tc, total)]
}

# Social value of SAI.
#   rate : discount rate (default 3%)
#   t0   : discounting reference and first period summed; defaults to the pulse
#          period read off the reported emissions
# Returns a list: $summary (one row), $components (NPV per damage component),
# $damages (annual paths), $pulses, $rep, $run.
svg_social_value <- function(rep = NULL, rate = SVG_DISCOUNT_RATE, t0 = NULL,
                             t_max = SVG_T_MAX, folder = SVG_RESULTS,
                             filters = list(), quiet = FALSE) {
  if (is.null(rep)) rep <- svg_report_table(t_max = t_max, folder = folder,
                                            filters = filters, quiet = quiet)
  dam <- svg_damage_components(rep)
  pulses <- svg_pulses(rep)
  ref <- pulses[scenario == SVG_GHGPULSE]
  if (!nrow(ref)) stop("No '", SVG_GHGPULSE, "' pulse found in the report.", call. = FALSE)
  if (is.null(t0)) t0 <- ref$pulse_t

  npv <- svg_npv(dam, c(SVG_GHGPULSE, SVG_SRMPULSE), rate, t0, by_component = TRUE)
  npv_ghg <- npv[scenario == SVG_GHGPULSE, total]
  npv_srm <- npv[scenario == SVG_SRMPULSE, total]
  mb_sai  <- -npv_srm
  scc     <- npv_ghg * 1e12 / ref$pulse_tons

  q <- dcast(rep[scenario %in% c(SVG_BASE, SVG_SRMPULSE)], t ~ scenario,
             value.var = c("qsrm", "forc_srm"))
  summary <- data.table(
    discount_rate    = rate,
    pulse_t          = t0,
    pulse_year       = SVG_BASE_YEAR + t0,
    t_max            = max(rep$t),
    pulse_GtCO2      = ref$pulse_GtCO2,
    npv_dam_pulse_tn = npv_ghg,                          # trillion USD
    npv_dam_sai_tn   = npv_srm,                          # trillion USD (negative = SAI helps)
    mb_sai_tn        = mb_sai,
    scc_usd_tco2     = scc,                              # both marginal quantities are
    mb_sai_usd_tco2  = mb_sai * 1e12 / ref$pulse_tons,   # per ton of the same pulse
    svg              = mb_sai / npv_ghg,
    svg_tco2         = mb_sai / npv_ghg * ref$pulse_tons,
    dq_srm_max       = q[, max(get(paste0("qsrm_", SVG_SRMPULSE)) - get(paste0("qsrm_", SVG_BASE)))],
    dforc_srm_max    = q[, max(get(paste0("forc_srm_", SVG_SRMPULSE)) - get(paste0("forc_srm_", SVG_BASE)))]
  )
  list(summary = summary, components = npv, damages = dam, pulses = pulses,
       rep = rep, run = attr(rep, "run"))
}


# -----------------------------------------------------------------------------
# Commitment time
# -----------------------------------------------------------------------------

# Social value of SAI as a function of the commitment time T: the marginal benefit
# accumulated from the pulse year up to (pulse year + T), discounted to the pulse
# year, as a share of the SCC of the pulse over the WHOLE horizon.
# NOTE this truncates the accounting window, not the deployment: the SAI in the
# srmpulse run runs to the end of the horizon either way, so this is "value
# accrued within T years". The `term` runs below are the deployment counterpart.
svg_commitment_profile <- function(x = NULL, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  cols <- c("temp", "penalty", "vll", "cost", "total")
  base <- x$damages[vrn == SVG_BASE, c("t", cols), with = FALSE]
  setnames(base, cols, paste0("b_", cols))
  p <- merge(x$damages[scenario == SVG_SRMPULSE], base, by = "t")[order(t)][t >= s$pulse_t]
  q <- merge(x$damages[scenario == SVG_GHGPULSE], base, by = "t")[order(t)][t >= s$pulse_t]
  disc <- function(z) (1 + s$discount_rate)^-(z - s$pulse_t)
  out <- data.table(
    commitment   = p$t - s$pulse_t,
    year         = SVG_BASE_YEAR + p$t,
    mb_tn        = -cumsum((p$total - p$b_total) * disc(p$t)),
    pulse_dam_tn =  cumsum((q$total - q$b_total) * disc(q$t)))
  out[, `:=`(mb_usd_tco2 = mb_tn * 1e12 / (s$pulse_GtCO2 * 1e9),
             frac_scc    = mb_tn / s$npv_dam_pulse_tn,
             pulse_frac  = pulse_dam_tn / s$npv_dam_pulse_tn)]
  setattr(out, "summary", s)
  out[]
}

# Value of a commitment of length T from the runs that actually terminate the
# marginal SAI after T years: the FULL-horizon NPV of the damage difference, so
# the post-termination rebound is priced in. Same discounting reference and pulse
# size as svg_social_value(), so it sits directly on the theoretical curve.
svg_terminated_value <- function(x = NULL, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  scen <- x$rep[vrn == SVG_TERM, unique(scenario)]
  if (!length(scen)) return(NULL)
  npv <- svg_npv(x$damages, scen, s$discount_rate, s$pulse_t)
  npv[, .(commitment = tc, mb_tn = -total,
          mb_usd_tco2 = -total * 1e12 / (s$pulse_GtCO2 * 1e9),
          frac_scc = -total / s$npv_dam_pulse_tn)][order(commitment)]
}

# SCC(t) along the background SAI deployment, from the runs that emit the CO2
# pulse tc years later. Each pulse is discounted to its OWN emission year, the
# standard SCC convention; `scc_pv_usd_tco2` brings it back to the reference pulse
# year, and `delay_value_usd_tco2` is the value of postponing the pulse by tc.
svg_scc_schedule <- function(x = NULL, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  scen <- x$rep[vrn == SVG_SCC, unique(scenario)]
  if (!length(scen)) return(NULL)
  out <- rbindlist(lapply(scen, function(sc) {
    p <- x$pulses[scenario == sc]
    n <- svg_npv(x$damages, sc, s$discount_rate, t0 = p$pulse_t, from = p$pulse_t)
    data.table(delay = p$tc, pulse_t = p$pulse_t, pulse_year = p$pulse_year,
               pulse_GtCO2 = p$pulse_GtCO2, npv_tn = n$total,
               scc_usd_tco2 = n$total * 1e12 / p$pulse_tons)
  }))
  out[, scc_pv_usd_tco2 := scc_usd_tco2 / (1 + s$discount_rate)^delay]
  out[, delay_value_usd_tco2 := s$scc_usd_tco2 - scc_pv_usd_tco2]
  setorder(out, delay)[]
}


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

# Marginal benefit of SAI against commitment time, in USD per tCO2 of the pulse,
# with the share of the full-horizon SCC on the right axis (the same curve over a
# constant, so one line carries both).
#   actual : svg_terminated_value() points, drawn on top of the theoretical line
#   delay  : svg_scc_schedule(), drawn as the value of merely postponing the pulse
#   panels : TRUE stacks the two readings instead of using a secondary axis
plot_svg_commitment <- function(prof = NULL, max_T = NULL, panels = FALSE,
                                actual = NULL, delay = NULL, ...) {
  if (is.null(prof)) prof <- svg_commitment_profile(...)
  s <- attr(prof, "summary")
  keep <- function(z) if (is.null(z) || is.null(max_T)) z else z[commitment <= max_T]
  d <- keep(prof)
  lab_abs <- "Marginal benefit of SAI [USD/tCO2 of the pulse]"
  lab_frc <- "Share of the full-horizon SCC [-]"
  sub <- sprintf("r = %s%%, denominator = SCC over t <= %d (%d)",
                 format(100 * s$discount_rate, digits = 3), s$t_max,
                 SVG_BASE_YEAR + s$t_max)
  if (panels) {
    long <- rbind(d[, .(commitment, panel = lab_abs, value = mb_usd_tco2)],
                  d[, .(commitment, panel = lab_frc, value = frac_scc)])
    long[, panel := factor(panel, levels = c(lab_abs, lab_frc))]
    return(ggplot(long, aes(commitment, value)) +
             geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
             geom_line(linewidth = 0.8, colour = "#1f78b4") +
             facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "left") +
             labs(x = paste0("Commitment time T [years from ", s$pulse_year, "]"),
                  y = NULL, title = "Social value of SAI by commitment time", subtitle = sub) +
             theme_bw(base_size = 11) +
             theme(strip.background = element_blank(), strip.placement = "outside",
                   panel.grid.minor = element_blank()))
  }
  cols <- c("theoretical (accounting window truncated at T)" = "#1f78b4",
            "actual (marginal SAI terminated at T)"          = "#e31a1c",
            "postponing the pulse by T (from SCC(t))"        = "#33a02c")
  p <- ggplot(d, aes(commitment, mb_usd_tco2)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line(aes(colour = names(cols)[1]), linewidth = 0.8) +
    labs(x = paste0("Commitment time T [years from ", s$pulse_year, "]"),
         colour = NULL, title = "Social value of SAI by commitment time", subtitle = sub) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  if (!is.null(delay))
    p <- p + geom_line(data = keep(data.table(commitment = delay$delay,
                                              mb_usd_tco2 = delay$delay_value_usd_tco2)),
                       aes(colour = names(cols)[3]), linewidth = 0.8, linetype = "dashed")
  if (!is.null(actual))
    p <- p + geom_point(data = keep(actual), aes(colour = names(cols)[2]), size = 2.4)
  p + scale_colour_manual(values = cols, breaks = names(cols)) +
    scale_y_continuous(name = lab_abs,
                       sec.axis = sec_axis(~ . / s$scc_usd_tco2, name = lab_frc))
}

# SCC(t) along the background SAI deployment: each point is a CO2 pulse emitted in
# that year, discounted to its own emission year.
plot_scc_schedule <- function(sched, s = NULL) {
  ggplot(sched, aes(pulse_year, scc_usd_tco2)) +
    geom_line(linewidth = 0.8, colour = "#6a3d9a") +
    geom_point(size = 2, colour = "#6a3d9a") +
    labs(x = "Year of the emission pulse", y = "SCC [USD/tCO2, discounted to the pulse year]",
         title = "SCC schedule under the background SAI deployment",
         subtitle = if (is.null(s)) NULL else
           sprintf("r = %s%%", format(100 * s$discount_rate, digits = 3))) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
}


# -----------------------------------------------------------------------------
# Console report
# -----------------------------------------------------------------------------
svg_report <- function(x) {
  s <- x$summary
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat("Social value of SAI  |  r = ", format(100 * s$discount_rate, digits = 3),
      "%  |  t <= ", s$t_max, " (", SVG_BASE_YEAR + s$t_max, ")\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")
  # Both marginal quantities are NPVs per ton of the same pulse, so the ratio of
  # the last two columns is the SVG whichever unit it is read in.
  cat(sprintf("  CO2 pulse                      %12.2f GtCO2 in %d\n", s$pulse_GtCO2, s$pulse_year))
  cat(sprintf("  NPV damage of the pulse        %12.4f tn USD   %10.2f USD/tCO2  (SCC)\n",
              s$npv_dam_pulse_tn, s$scc_usd_tco2))
  cat(sprintf("  NPV damage of the extra SAI    %12.4f tn USD   %10.2f USD/tCO2\n",
              s$npv_dam_sai_tn, -s$mb_sai_usd_tco2))
  cat(sprintf("  Marginal benefit of the SAI    %12.4f tn USD   %10.2f USD/tCO2\n",
              s$mb_sai_tn, s$mb_sai_usd_tco2))
  cat(sprintf("  SVG = MB(SAI) / SCC(pulse)     %12.4f  [-]\n", s$svg))
  cat(sprintf("      = %.4g tCO2 of avoided emissions\n", s$svg_tco2))
  cat(sprintf("  Extra SAI: peak dQ_SRM %.4g, peak dFORC_SRM %.4g W/m2\n",
              s$dq_srm_max, s$dforc_srm_max))
  cat("\n  NPV by damage component [trillion USD, difference from baseline]:\n")
  print(x$components)
  prof <- svg_commitment_profile(x)
  pay <- prof[mb_tn > 0, min(commitment)]
  cat("\n  Value by commitment time T [years from ", s$pulse_year, "]:\n", sep = "")
  cat(sprintf("    payback (first T with MB > 0): %s\n",
              if (is.finite(pay)) paste(pay, "years") else "never"))
  for (sh in c(0.5, 0.9, 0.99)) {
    hit <- prof[frac_scc >= sh * s$svg, min(commitment)]
    cat(sprintf("    %2.0f%% of the full-horizon value reached at T = %s\n",
                100 * sh, if (is.finite(hit)) paste(hit, "years") else "never"))
  }
  print(prof[commitment %in% c(10, 25, 50, 100, 200, 400),
             .(T = commitment, year, mb_tn = round(mb_tn, 3),
               mb_usd_tco2 = round(mb_usd_tco2, 1), frac_scc = round(frac_scc, 4))])
  invisible(x)
}


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
svg_main <- function(rate = SVG_DISCOUNT_RATE, t_max = SVG_T_MAX, filters = list(),
                     outdir = SVG_OUTPUT_DIR, plots = TRUE) {
  rep <- svg_report_table(t_max = t_max, filters = filters)
  print(attr(rep, "run"))
  cat("  variants:", paste(rep[, uniqueN(scenario)], "scenarios;"),
      "quantities:", paste(svg_quantities(rep), collapse = ", "), "\n")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  x <- svg_social_value(rep, rate = rate)
  svg_report(x)
  fwrite(x$summary, file.path(outdir, "svg_summary.csv"))
  fwrite(x$components, file.path(outdir, "svg_damage_components.csv"))

  prof   <- svg_commitment_profile(x)
  actual <- svg_terminated_value(x)
  sched  <- svg_scc_schedule(x)
  fwrite(prof, file.path(outdir, "svg_commitment_profile.csv"))

  if (!is.null(actual)) {
    fwrite(actual, file.path(outdir, "svg_terminated_value.csv"))
    cat("\nActual commitment experiments (marginal SAI terminated after T years):\n")
    print(merge(actual, prof[, .(commitment, theory_usd_tco2 = mb_usd_tco2)], by = "commitment"
                )[, .(T = commitment, actual_usd_tco2 = round(mb_usd_tco2, 1),
                      theory_usd_tco2 = round(theory_usd_tco2, 1),
                      actual_frac_scc = round(frac_scc, 4),
                      actual_over_theory = round(mb_usd_tco2 / theory_usd_tco2, 3))])
  }
  if (!is.null(sched)) {
    fwrite(sched, file.path(outdir, "svg_scc_schedule.csv"))
    cat("\nSCC schedule under the background SAI deployment:\n")
    print(sched[, .(delay, pulse_year, scc_usd_tco2 = round(scc_usd_tco2, 1),
                    scc_pv_usd_tco2 = round(scc_pv_usd_tco2, 1),
                    delay_value_usd_tco2 = round(delay_value_usd_tco2, 1))])
  }

  # Sensitivity of the ratio to the discount rate.
  rates <- sort(unique(c(0.01, 0.02, 0.03, 0.05, rate)))
  sens <- rbindlist(lapply(rates, function(rr) svg_social_value(rep, rate = rr)$summary))
  fwrite(sens, file.path(outdir, "svg_discount_sensitivity.csv"))
  cat("\nDiscount-rate sensitivity:\n")
  print(sens[, .(r = discount_rate, scc_usd_tco2, mb_sai_usd_tco2, svg)])

  if (plots) {
    ggsave(file.path(outdir, "svg_commitment.png"),
           plot_svg_commitment(prof, actual = actual, delay = sched),
           width = 9, height = 5.5, dpi = 200)
    ggsave(file.path(outdir, "svg_commitment_zoom.png"),
           plot_svg_commitment(prof, actual = actual, delay = sched, max_T = 200),
           width = 9, height = 5.5, dpi = 200)
    ggsave(file.path(outdir, "svg_commitment_panels.png"),
           plot_svg_commitment(prof, panels = TRUE), width = 8, height = 7, dpi = 200)
    if (!is.null(sched))
      ggsave(file.path(outdir, "scc_schedule.png"), plot_scc_schedule(sched, x$summary),
             width = 8, height = 5, dpi = 200)
    for (v in c("tatm", "damages", "forc_srm", "qsrm"))
      ggsave(file.path(outdir, paste0("diff_", v, ".png")),
             plot_gdx_diff(v, rep = rep, t_max = 480), width = 8, height = 5, dpi = 200)
  }
  cat("\nWrote CSVs", if (plots) " and figures" else "", " to ", outdir, "\n", sep = "")
  invisible(list(rep = rep, value = x, profile = prof, actual = actual, schedule = sched))
}

# Runs on `Rscript Paper_SVG/svg_analysis.R` and on source(); the functions above
# stay available either way. options(svg.no.main = TRUE) loads them without running.
if (!isTRUE(getOption("svg.no.main"))) {
  svg_results <- svg_main()
}
