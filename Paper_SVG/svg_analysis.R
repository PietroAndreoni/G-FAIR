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
#   ghgstore (tc > 0) the CO2 mirror of `term`: background SAI only, the pulse
#                     REMOVED at `pulse_time` and re-emitted tc years later, so
#                     the cooling is a finite commitment bought with carbon
#                     instead of aerosol.
#
# With damages D(t) [trillion USD] and the discount factor DF(t|t0) of the curve
# in svg_discounting.R -- by default the Ramsey rate r(t) = rho + eta*g(t) built
# from the model's own per-capita income path, which declines from about 3%/yr
# today to rho once growth has run out; a constant rate is the same object with
# a flat r, for which DF(t|t0) = (1+r)^-(t-t0):
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
#   svg_report(svg_social_value(rate = 0.05))        # a constant rate instead
#   svg_report(svg_social_value(rate = svg_ramsey_curve(rep$t, rep$ypc, rho = .01)))
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
# Structural settings (calibration lives in svg_config.R, sourced below)
# -----------------------------------------------------------------------------
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
SVG_GHGSTORE <- "ghgstore"
SVG_SRMDELAY <- "srmdelay"


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

# Directory this script lives in, so its companions (svg_config.R, ...) and its
# outputs travel with it if the folder is moved. The file being sourced wins over
# --file=: when another script source()s this one under Rscript, --file points at
# that outer script, not at this file.
svg_script_dir <- function() {
  p <- c(unlist(lapply(sys.frames(), function(f) f$ofile)),
         sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  p <- p[nzchar(p)]
  if (length(p)) normalizePath(dirname(p[1]), "/", FALSE) else normalizePath(getwd(), "/", FALSE)
}

SVG_ROOT       <- svg_project_root()
SVG_RESULTS    <- file.path(SVG_ROOT, "Results")   # where svg.gms unloads its gdx
SVG_DIR        <- svg_script_dir()
SVG_OUTPUT_DIR <- file.path(SVG_DIR, "output")

# Calibration, grids and tolerances (single control file), then the hazard and
# termination modules that depend on them.
for (.f in c("svg_config.R", "svg_discounting.R", "svg_termination.R",
             "svg_hazard.R")) {
  .p <- file.path(SVG_DIR, .f)
  if (!file.exists(.p)) stop("Missing ", .f, " next to svg_analysis.R.", call. = FALSE)
  source(.p)
}
rm(.f, .p)

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
  # Per-capita income: exogenous, one path shared by every scenario, and the
  # input the Ramsey discount curve is built from (svg_discounting.R).
  ypc <- tryCatch(svg_read(file, "ypc", t_max = t_max)[, .(t, ypc = value)],
                  error = function(e) NULL)
  if (!is.null(ypc)) w <- merge(w, ypc, by = "t")
  w[, scenario := fifelse(tc == 0L, vrn, paste0(vrn, tc))]
  setcolorder(w, c("scenario", "vrn", "tc", "t"))
  setorder(w, vrn, tc, t)

  # Drop any scenario GAMS did not solve cleanly. svg.gms records the status
  # instead of aborting, so one bad state cannot discard a whole run -- but it
  # must not reach the analysis either.
  st <- svg_solve_status(run)
  if (!is.null(st)) {
    bad <- st[ok == FALSE]
    if (nrow(bad)) {
      if (!quiet)
        message("svg_report_table(): dropping ", nrow(bad), " unsolved scenario(s): ",
                paste(bad[, paste0(vrn, tc)], collapse = ", "))
      w <- w[!bad[, .(vrn, tc)], on = c("vrn", "tc")]
    }
    setattr(w, "solve_status", st)
  }
  setattr(w, "run", run)
  w[]
}

# Does this gdx hold every scenario family the analysis needs? A file left over
# from an earlier version of experiments/svg.gms exists and reads fine but is
# missing whole variants, so "the file is there" is not a usable skip test.
SVG_REQUIRED_VARIANTS <- c(SVG_BASE, SVG_GHGPULSE, SVG_SRMPULSE, SVG_TERM,
                           SVG_SCC, SVG_GHGSTORE, SVG_SRMDELAY)

svg_gdx_complete <- function(path, need = SVG_REQUIRED_VARIANTS) {
  if (!file.exists(path)) return(FALSE)
  st <- tryCatch(svg_solve_status(path), error = function(e) NULL)
  !is.null(st) && all(need %in% st$vrn)
}

# Solver outcome of every stored scenario. solve_status(vrn,tc,*) does not carry
# time as its last index, so it is read directly rather than through svg_read().
svg_solve_status <- function(run) {
  svg_init_gams()
  got <- gdxtools::batch_extract("solve_status", files = unname(as.character(run))[1])$solve_status
  if (is.null(got) || !nrow(got)) return(NULL)   # gdx predates solve_status
  st <- as.data.table(got)
  setnames(st, names(st)[3], "stat")
  w <- dcast(st, vrn + tc ~ stat, value.var = "value")
  w[, tc := as.integer(as.character(tc))]
  w[, ok := solvestat == SVG_OK_SOLVESTAT & modelstat %in% SVG_OK_MODELSTAT]
  w[]
}

# Quantities the report carries, i.e. what gdx_diff() / plot_gdx_diff() can draw.
svg_quantities <- function(rep)
  setdiff(names(rep), c("scenario", "vrn", "tc", "t", "y", "ypc"))


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
# `ghgstore` is the one variant with two deviations -- a removal, then the
# re-emission tc years later -- so the pulse size is the FIRST of them (negative,
# it is a removal) and `net_GtCO2` is what is left in the atmosphere at the end.
# Summing would net those two to zero and any per-tonne normalisation would then
# divide by it. Every other variant has a single deviation, for which first == sum.
svg_pulses <- function(rep) {
  base <- rep[vrn == SVG_BASE, .(t, base_emi = wemi_co2)]
  d <- merge(rep[vrn != SVG_BASE], base, by = "t")
  d[, demi := wemi_co2 - base_emi]
  out <- d[abs(demi) > 1e-9][order(t),
           .(pulse_t = first(t), pulse_GtCO2 = first(demi), net_GtCO2 = sum(demi)),
           by = .(vrn, tc, scenario)]
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
#   rate : a discount curve (svg_discounting.R) or a constant rate
svg_npv <- function(dam, scenarios, rate, t0, from = t0, by_component = FALSE) {
  cols <- c("temp", "penalty", "vll", "cost", "total")
  base <- dam[vrn == SVG_BASE, c("t", cols), with = FALSE]
  setnames(base, cols, paste0("b_", cols))
  d <- merge(dam[scenario %in% scenarios], base, by = "t")
  for (cc in cols) d[, (cc) := get(cc) - get(paste0("b_", cc))]
  d <- d[t >= from]
  d[, discount := svg_df(rate, t, t0)]
  out <- d[, lapply(.SD, function(v) sum(v * discount)), by = .(scenario, vrn, tc),
           .SDcols = cols]
  if (by_component) out else out[, .(scenario, vrn, tc, total)]
}

# Social value of SAI.
#   rate : NULL for the configured curve (Ramsey off the run's own ypc path), a
#          number for a constant rate, or a curve from svg_discounting.R
#   t0   : discounting reference and first period summed; defaults to the pulse
#          period read off the reported emissions
# Returns a list: $summary (one row), $components (NPV per damage component),
# $damages (annual paths), $pulses, $disc (the curve everything downstream uses),
# $rep, $run.
svg_social_value <- function(rep = NULL, rate = NULL, t0 = NULL,
                             t_max = SVG_T_MAX, folder = SVG_RESULTS,
                             filters = list(), quiet = FALSE) {
  if (is.null(rep)) rep <- svg_report_table(t_max = t_max, folder = folder,
                                            filters = filters, quiet = quiet)
  disc <- svg_resolve_discount(rate, rep)
  rate <- disc
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
    # `discount_rate` stays a single number -- the rate in force at the pulse
    # period -- so anything with one slot to print still has one; the curve
    # itself travels as x$disc, and `discount_label` says what it is.
    discount_mode    = disc$mode,
    discount_rate    = svg_rate_scalar(disc, t0),
    discount_rho     = disc$rho %||% NA_real_,
    discount_eta     = disc$eta %||% NA_real_,
    discount_label   = svg_discount_label(disc, t0),
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
       disc = disc, rep = rep, run = attr(rep, "run"))
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
  disc <- function(z) svg_df(x$disc, z, s$pulse_t)
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
  npv <- svg_npv(x$damages, scen, x$disc, s$pulse_t)
  npv[, .(commitment = tc, mb_tn = -total,
          mb_usd_tco2 = -total * 1e12 / (s$pulse_GtCO2 * 1e9),
          frac_scc = -total / s$npv_dam_pulse_tn)][order(commitment)]
}

# The carbon counterpart of svg_terminated_value(), from the `ghgstore` runs: the
# pulse is REMOVED at the reference pulse period and re-emitted T years later, so
# this is the value of holding a tonne of CO2 out of the atmosphere for T years,
# with the re-emission rebound priced in. Background SAI only -- no marginal
# aerosol -- so it is the same commitment question asked of the other lever.
# Same discounting reference, same pulse size and the same tc grid as the SAI
# runs, so the two sit on one axis without rescaling. Once T is long enough that
# the re-emission is discounted away the value converges to the SCC from below --
# from below even in the limit, because damages are convex, so removing a tonne
# saves slightly less than adding one costs.
#
# In a LINEAR model this equals the "postponing the pulse by T" line that
# svg_scc_schedule() constructs from two separate runs; here it is simulated in
# one, so the gap between the two is the carbon-cycle and damage nonlinearity.
svg_offset_value <- function(x = NULL, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  scen <- x$rep[vrn == SVG_GHGSTORE, unique(scenario)]
  if (!length(scen)) return(NULL)
  npv <- svg_npv(x$damages, scen, x$disc, s$pulse_t)
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
  # ghgpulse is the same experiment at delay 0, so it anchors the schedule at the
  # reference pulse year rather than starting it at the first delayed pulse.
  scen <- x$rep[vrn %in% c(SVG_GHGPULSE, SVG_SCC), unique(scenario)]
  if (!length(scen)) return(NULL)
  out <- rbindlist(lapply(scen, function(sc) {
    p <- x$pulses[scenario == sc]
    n <- svg_npv(x$damages, sc, x$disc, t0 = p$pulse_t, from = p$pulse_t)
    data.table(delay = p$tc, pulse_t = p$pulse_t, pulse_year = p$pulse_year,
               pulse_GtCO2 = p$pulse_GtCO2, npv_tn = n$total,
               scc_usd_tco2 = n$total * 1e12 / p$pulse_tons)
  }))
  # Back to the reference pulse year along the same curve the NPVs used, from the
  # actual periods rather than from `delay`: under a declining rate the factor is
  # not a power of a single number.
  out[, scc_pv_usd_tco2 := scc_usd_tco2 * svg_df(x$disc, pulse_t, s$pulse_t)]
  out[, delay_value_usd_tco2 := s$scc_usd_tco2 - scc_pv_usd_tco2]
  setorder(out, delay)[]
}


# -----------------------------------------------------------------------------
# State-dependent social value  SV*(G | s,h)   (spec sections 2.10 - 2.13)
# -----------------------------------------------------------------------------

# SCC and marginal damages on the full annual grid, interpolated from the state
# grid the delayed-pulse runs provide. Linear interpolation with flat
# extrapolation: the schedule is smooth and monotone over 10-50 year steps, and a
# spline would be free to oscillate between them. MD comes from the schedule's
# own identity dSCC/dt = r*SCC - MD, which is the spec's section 8.3 check used
# forwards.
svg_scc_annual <- function(sched, t_grid, rate = SVG_DISCOUNT_RATE) {
  scc <- stats::approx(sched$pulse_t, sched$scc_usd_tco2, xout = t_grid, rule = 2)$y
  dscc <- c(diff(scc), 0)                       # forward difference, dt = 1
  # r is r(t) under a declining curve, so the identity is applied period by period.
  data.table(t = t_grid, scc = scc, r = svg_rate_at(rate, t_grid),
             md = svg_rate_at(rate, t_grid) * scc - dscc)
}

# Annual marginal damage flow MD(t) [USD/tCO2/yr]: the damage the reference pulse
# adds in year t, straight from the model. This is the quantity the SCC integrates,
# so it is the honest MD schedule; the difference quotient of the interpolated SCC
# in svg_scc_annual() exists only for the dSCC/dt = r*SCC - MD identity check.
svg_md_flow <- function(x) {
  p <- x$pulses[scenario == SVG_GHGPULSE]
  b <- x$damages[vrn == SVG_BASE, .(t, b_total = total)]
  d <- merge(x$damages[scenario == SVG_GHGPULSE, .(t, total)], b, by = "t")
  d[, .(t, md = (total - b_total) * SVG_TN_USD / (p$pulse_GtCO2 * SVG_GTCO2_TO_TCO2))]
}

# One row per starting state t0 along the background pathway, with everything the
# social-value integral and its diagnostics need. States are the tc grid: tc = 0
# is the reference pulse period (srmpulse / ghgpulse), tc > 0 the delayed ones
# (srmdelay / scc).
svg_states <- function(x = NULL, beta = SVG_BETA, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  sched <- svg_scc_schedule(x)
  if (is.null(sched)) return(NULL)
  have <- x$rep[vrn %in% c(SVG_SRMPULSE, SVG_SRMDELAY), unique(tc)]
  st <- sched[delay %in% have]
  if (!nrow(st)) return(NULL)
  base <- x$damages[vrn == SVG_BASE]
  bg <- base[, .(t, G = qsrm, M = tatm_ghg - tatm, T = tatm, T_ghg = tatm_ghg,
                 forc_srm, tot_forc)]
  out <- merge(st[, .(tc = delay, t0 = pulse_t, year = pulse_year, pulse_GtCO2,
                      scc0 = scc_usd_tco2)],
               bg, by.x = "t0", by.y = "t")
  # Model-implied marginal efficacy of SAI at this state (see the header note):
  # eta = 1 - (T_GHG/TATM)*(1-cos(theta)); reported as a diagnostic only.
  theta <- svg_scalar(x$run, "srm_angle") * pi / 180
  out[, eta := 1 - (T_ghg / T) * (1 - cos(theta))]
  # Annual marginal damage flow, taken from the reference pulse's own damage
  # stream rather than by differentiating the interpolated SCC: the state grid is
  # 10-50 years wide, so a difference quotient across it is dominated by the
  # interpolation and can even come out negative.
  mdf <- svg_md_flow(x)
  out[, md0 := mdf$md[match(t0, mdf$t)]]
  out[, deltaL0 := svg_marginal_termination_loss(M, scc0, beta = beta)]
  setorder(out, tc)[]
}

# Social value of the marginal intervention at every state, for every hazard.
# The operating flow is the model's own damage difference for the marginal SAI
# bought at that state (srmpulse at tc = 0, srmdelay after), expressed per tonne
# of the pulse it offsets so it is directly comparable with the SCC.
#   keep_paths : also return the full SV(H) curve per state (spec section 7.3)

# Annual operating benefit and termination loss of the marginal intervention
# bought at one state, both in USD/tCO2 of the pulse it offsets, from t0 to the
# end of the horizon. This is the input the commitment integral consumes.
svg_state_flows <- function(x, tc = 0L, states = NULL, ann = NULL) {
  tc_val <- tc          # the argument would otherwise shadow the tc column
  if (is.null(states)) states <- svg_states(x)
  st <- states[tc == tc_val]
  if (!nrow(st)) return(NULL)
  if (is.null(ann)) ann <- svg_deltaL_annual(x)
  base <- x$damages[vrn == SVG_BASE, .(t, b_total = total)]
  vrn_i <- if (tc_val == 0L) SVG_SRMPULSE else SVG_SRMDELAY
  d <- x$damages[vrn == vrn_i & tc == tc_val, .(t, total)]
  if (!nrow(d)) return(NULL)
  d <- merge(d, base, by = "t")[t >= st$t0][order(t)]
  list(t0 = st$t0, state = st,
       # Benefit is avoided damage, so the sign flips; per tonne of the pulse.
       benefit = -(d$total - d$b_total) * SVG_TN_USD /
                  (st$pulse_GtCO2 * SVG_GTCO2_TO_TCO2),
       deltaL = ann$deltaL[match(d$t, ann$t)],
       t = d$t)
}

# DeltaL(t) on the annual grid: the interpolated SCC schedule and the background
# masked warming, through the reduced form in svg_termination.R.
svg_deltaL_annual <- function(x, beta = SVG_BETA) {
  sched <- svg_scc_schedule(x)
  base <- x$damages[vrn == SVG_BASE, .(t, M = tatm_ghg - tatm)]
  ann <- svg_scc_annual(sched, base$t, rate = x$disc)
  ann <- merge(ann, base, by = "t")
  ann[, deltaL := svg_marginal_termination_loss(M, scc, beta = beta)][]
}

# SV(H) at one state under one hazard, with the optimum attached. This is the
# curve the commitment figure draws when termination is included.
svg_commitment_hazard <- function(x = NULL, h = SVG_HAZARD_MAIN,
                                  beta = SVG_BETA, tc = 0L, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  fl <- svg_state_flows(x, tc, states = svg_states(x, beta = beta),
                        ann = svg_deltaL_annual(x, beta = beta))
  if (is.null(fl)) return(NULL)
  sv <- svg_commitment_value(fl$benefit, fl$deltaL, h = h, rate = x$disc, t0 = fl$t0)
  setattr(sv, "opt", svg_optimize_commitment(sv))
  setattr(sv, "h", h)
  setattr(sv, "beta", beta)
  setattr(sv, "state", fl$state)
  sv[]
}

svg_value_by_state <- function(x = NULL, hazards = SVG_HAZARDS,
                               beta = SVG_BETA, keep_paths = FALSE, ...) {
  if (is.null(x)) x <- svg_social_value(...)
  s <- x$summary
  states <- svg_states(x, beta = beta)
  if (is.null(states)) return(NULL)
  ann <- svg_deltaL_annual(x, beta = beta)

  rows <- list(); paths <- list()
  for (i in seq_len(nrow(states))) {
    tc_i <- states$tc[i]; t0 <- states$t0[i]
    fl <- svg_state_flows(x, tc_i, states = states, ann = ann)
    if (is.null(fl)) next
    ben <- fl$benefit; dl <- fl$deltaL
    for (h in hazards) {
      sv <- svg_commitment_value(ben, dl, h = h, rate = x$disc, t0 = t0)
      opt <- svg_optimize_commitment(sv)
      rows[[length(rows) + 1L]] <- data.table(
        tc = tc_i, t0 = t0, year = states$year[i], h = h,
        G = states$G[i], M = states$M[i], T = states$T[i], eta = states$eta[i],
        scc0 = states$scc0[i], md0 = states$md0[i], deltaL0 = states$deltaL0[i],
        H_star = opt$H_star, H_class = opt$class, SV_star = opt$SV_star,
        SV_star_over_scc = opt$SV_star / states$scc0[i],
        SV_at_last_H = opt$sv_last, last_H = opt$last_H,
        tail_gain = opt$tail_gain, tail_converged = opt$tail_converged)
      if (keep_paths)
        paths[[length(paths) + 1L]] <- sv[, .(tc = tc_i, t0, h, H, sv,
                                              operating, termination)]
    }
  }
  res <- rbindlist(rows)
  if (keep_paths) setattr(res, "paths", rbindlist(paths))
  setattr(res, "summary", s)
  setorder(res, h, tc)[]
}


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

# Marginal benefit of SAI against commitment time, in USD per tCO2 of the pulse,
# with the share of the full-horizon SCC on the right axis (the same curve over a
# constant, so one line carries both).
#   actual : svg_terminated_value() points, drawn on top of the theoretical line
#   delay  : svg_scc_schedule(), drawn as the value of merely postponing the pulse
#   offset : svg_offset_value(), the same commitment question asked of carbon --
#            a tonne removed now and re-emitted at T -- so the two levers are read
#            off one axis. It is the simulated counterpart of the `delay` line.
#   hazard : svg_commitment_hazard(), the curve once the termination hazard is
#            priced in; its optimum (H*, SV*) is marked and labelled
#   panels : TRUE stacks the two readings instead of using a secondary axis
plot_svg_commitment <- function(prof = NULL, max_T = NULL, panels = FALSE,
                                actual = NULL, delay = NULL, hazard = NULL,
                                offset = NULL, ...) {
  if (is.null(prof)) prof <- svg_commitment_profile(...)
  s <- attr(prof, "summary")
  keep <- function(z) if (is.null(z) || is.null(max_T)) z else z[commitment <= max_T]
  d <- keep(prof)
  lab_abs <- "Marginal benefit of SAI [USD/tCO2 of the pulse]"
  lab_frc <- "Share of the full-horizon SCC [-]"
  sub <- sprintf("%s; denominator = SCC over t <= %d (%d)",
                 s$discount_label, s$t_max, SVG_BASE_YEAR + s$t_max)
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
  haz_lab <- if (is.null(hazard)) "with termination hazard" else
    sprintf("with termination hazard (h = %s/yr)", format(attr(hazard, "h")))
  if (!is.null(hazard)) {
    o <- attr(hazard, "opt")
    sub <- paste0(sub, sprintf("\nh = %s/yr, beta = %s: H* = %d yr (%s), SV* = %.1f USD/tCO2 = %.0f%% of SCC",
                               format(attr(hazard, "h")),
                               format(attr(hazard, "beta") %||% SVG_BETA), o$H_star,
                               o$class, o$SV_star, 100 * o$SV_star / s$scc_usd_tco2))
  }
  # Categorical hues in fixed order, never cycled: a series keeps its colour
  # whether or not the others are drawn. The fifth is near-black rather than the
  # next ColorBrewer-Paired step, which was checked and rejected -- orange sits
  # at OKLab dE 0.7 from the green under protanopia, i.e. the same colour. Black
  # is >= 28 from all four in every vision mode, so it costs the palette nothing.
  # Every series also carries a mark or linetype of its own, so identity never
  # rests on hue alone (the red/green pair is dE 3.0 under deuteranopia and is
  # separated by geometry: bare points vs a dashed line).
  cols <- stats::setNames(c("#1f78b4", "#e31a1c", "#33a02c", "#6a3d9a", "#111111"),
                          c("theoretical (accounting window truncated at T)",
                            "actual (marginal SAI terminated at T)",
                            "postponing the pulse by T (from SCC(t))",
                            haz_lab,
                            "carbon offset (CO2 removed at 0, re-emitted at T)"))
  ttl <- if (is.null(offset)) "Social value of SAI by commitment time" else
    "Social value of a finite commitment: marginal SAI vs temporary carbon storage"
  p <- ggplot(d, aes(commitment, mb_usd_tco2)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line(aes(colour = names(cols)[1]), linewidth = 0.8) +
    labs(x = paste0("Commitment time T [years from ", s$pulse_year, "]"),
         colour = NULL, title = ttl, subtitle = sub) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  if (!is.null(delay))
    p <- p + geom_line(data = keep(data.table(commitment = delay$delay,
                                              mb_usd_tco2 = delay$delay_value_usd_tco2)),
                       aes(colour = names(cols)[3]), linewidth = 0.8, linetype = "dashed")
  # Points, not a line, and drawn straight after the `delay` line: this is the
  # simulated version of that constructed curve and lands on top of it, so a line
  # would simply hide it. Same idiom as `actual` over the theoretical line --
  # line = constructed, marks = the run that was actually solved. Triangles keep
  # it apart from the `actual` circles without relying on hue.
  if (!is.null(offset))
    p <- p + geom_point(data = keep(offset), aes(colour = names(cols)[5]),
                        size = 2.2, shape = 17)
  if (!is.null(actual))
    p <- p + geom_point(data = keep(actual), aes(colour = names(cols)[2]), size = 2.4)
  if (!is.null(hazard)) {
    hz <- keep(data.table(commitment = hazard$H, mb_usd_tco2 = hazard$sv))
    opt <- attr(hazard, "opt")
    p <- p + geom_line(data = hz, aes(colour = names(cols)[4]), linewidth = 0.9)
    # The optimum is the point of this exercise: mark it, drop guides to both
    # axes, and say what it is.
    if (is.null(max_T) || opt$H_star <= max_T) {
      om <- data.table(commitment = opt$H_star, mb_usd_tco2 = opt$SV_star)
      # Label to whichever side keeps it inside the panel.
      x_hi <- max(hz$commitment)
      right <- opt$H_star > 0.6 * x_hi
      p <- p +
        geom_segment(data = om, aes(x = commitment, xend = commitment,
                                    y = -Inf, yend = mb_usd_tco2),
                     linetype = "dotted", linewidth = 0.4, colour = "#6a3d9a") +
        geom_segment(data = om, aes(x = -Inf, xend = commitment,
                                    y = mb_usd_tco2, yend = mb_usd_tco2),
                     linetype = "dotted", linewidth = 0.4, colour = "#6a3d9a") +
        geom_point(data = om, aes(colour = names(cols)[4]), size = 3.2, shape = 18) +
        annotate("text", x = opt$H_star, y = opt$SV_star,
                 label = sprintf(if (right) "H* = %d yr, SV* = %.1f $/tCO2  " else
                                            "  H* = %d yr, SV* = %.1f $/tCO2",
                                 opt$H_star, opt$SV_star),
                 hjust = if (right) 1 else 0, vjust = -0.9, size = 3.1, colour = "#6a3d9a")
    }
  }
  p + scale_colour_manual(values = cols, breaks = names(cols)) +
    guides(colour = guide_legend(nrow = if (is.null(offset)) 2 else 3, byrow = TRUE)) +
    scale_y_continuous(name = lab_abs,
                       sec.axis = sec_axis(~ . / s$scc_usd_tco2, name = lab_frc))
}

# Line-only SV(H) sensitivity figure. Each panel moves one parameter around the
# configured baseline while holding the other three fixed. The rho, eta and beta
# panels all retain the headline termination hazard; the h panel adds h = 0 as
# the deterministic reference to the configured hazard grid.
svg_commitment_sensitivity <- function(rep = NULL, max_T = SVG_FIG_MAX_T,
                                       rho = SVG_RHO_SENS, eta = SVG_ETA_SENS,
                                       beta = SVG_BETA_SENS,
                                       h = sort(unique(c(0, SVG_HAZARDS))), ...) {
  if (is.null(rep)) rep <- svg_report_table(...)

  one_curve <- function(parameter, value, rho_i = SVG_RHO, eta_i = SVG_ETA,
                        beta_i = SVG_BETA, h_i = SVG_HAZARD_MAIN) {
    disc <- svg_resolve_discount(rep = rep, rho = rho_i, eta = eta_i, mode = "ramsey")
    x <- svg_social_value(rep, rate = disc)
    z <- svg_commitment_hazard(x, h = h_i, beta = beta_i, tc = 0L)
    if (is.null(z)) stop("Cannot construct the commitment-hazard sensitivity curve.",
                         call. = FALSE)
    z[H <= max_T, .(parameter, value, H, sv)]
  }

  curves <- rbindlist(c(
    lapply(rho, function(z) one_curve("rho", z, rho_i = z)),
    lapply(eta, function(z) one_curve("eta", z, eta_i = z)),
    lapply(beta, function(z) one_curve("beta", z, beta_i = z)),
    lapply(h, function(z) one_curve("h", z, h_i = z))))

  baseline <- c(rho = SVG_RHO, eta = SVG_ETA, beta = SVG_BETA,
                h = SVG_HAZARD_MAIN)
  curves[, case := fifelse(value < baseline[parameter], "Lower",
                           fifelse(value > baseline[parameter], "Upper", "Baseline"))]
  curves[, case := factor(case, levels = c("Lower", "Baseline", "Upper"))]

  fmt <- function(z) paste(format(z, digits = 3, trim = TRUE), collapse = ", ")
  fmt_pct <- function(z) paste0(format(100 * z, digits = 3, trim = TRUE), "%",
                                collapse = ", ")
  panel_labels <- c(
    rho  = paste0("rho: ", fmt_pct(sort(unique(rho)))),
    eta  = paste0("eta: ", fmt(sort(unique(eta)))),
    beta = paste0("beta: ", fmt(sort(unique(beta)))),
    h    = paste0("h: ", fmt(sort(unique(h))), " /yr"))
  curves[, panel := factor(panel_labels[parameter], levels = unname(panel_labels))]
  setorder(curves, parameter, value, H)
  setattr(curves, "max_T", max_T)
  curves[]
}

plot_svg_commitment_sensitivity <- function(sens) {
  max_T <- attr(sens, "max_T") %||% max(sens$H)
  ggplot(sens, aes(H, sv, colour = case, linetype = case, linewidth = case,
                   group = interaction(parameter, value))) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line() +
    facet_wrap(~panel, ncol = 2) +
    scale_colour_manual(values = c(Lower = "#1f78b4", Baseline = "#111111",
                                   Upper = "#e31a1c"), drop = FALSE) +
    scale_linetype_manual(values = c(Lower = "dashed", Baseline = "solid",
                                     Upper = "dotdash"), drop = FALSE) +
    scale_linewidth_manual(values = c(Lower = 0.75, Baseline = 1.0, Upper = 0.75),
                           drop = FALSE) +
    scale_x_continuous(limits = c(0, max_T), expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(ylim = c(-50, 50)) +
    labs(x = "Commitment H [years]", y = "Social value [USD/tCO2]",
         colour = NULL, linetype = NULL, linewidth = NULL,
         title = "Social value versus commitment: parameter sensitivity",
         subtitle = paste0("One parameter at a time; rho, eta and beta panels use ",
                           "h = ", format(SVG_HAZARD_MAIN), "/yr; ",
                           "the h panel includes the no-hazard reference")) +
    guides(colour = guide_legend(override.aes = list(linewidth = c(0.75, 1, 0.75)))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

# The discount curve itself: the Ramsey rate and the per-capita growth rate it is
# built from, with the constant rate it replaces drawn for reference. Every other
# figure in this folder is a weighted sum over this line, so it is worth seeing.
plot_discount_curve <- function(disc, max_year = 2300, ref = SVG_DISCOUNT_RATE) {
  d <- svg_discount_table(disc)[year <= max_year]
  long <- rbind(d[, .(year, series = "discount rate r(t)", value = r)],
                d[, .(year, series = "per-capita growth g(t)", value = g)])
  long <- long[!is.na(value)]
  ggplot(long, aes(year, 100 * value, colour = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_hline(yintercept = 100 * ref, linetype = "dashed", linewidth = 0.4,
               colour = "grey40") +
    annotate("text", x = max_year, y = 100 * ref, vjust = -0.6, hjust = 1,
             size = 3.1, colour = "grey40",
             label = sprintf("constant r = %s%%", format(100 * ref, digits = 3))) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("discount rate r(t)" = "#6a3d9a",
                                   "per-capita growth g(t)" = "#33a02c")) +
    labs(x = "Year", y = "[%/yr]", colour = NULL,
         title = "Discount curve applied to every social-value integral",
         subtitle = svg_discount_label(disc)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

# SCC(t) along the background SAI deployment: each point is a CO2 pulse emitted in
# that year, discounted to its own emission year.
plot_scc_schedule <- function(sched, s = NULL) {
  ggplot(sched, aes(pulse_year, scc_usd_tco2)) +
    geom_line(linewidth = 0.8, colour = "#6a3d9a") +
    geom_point(size = 2, colour = "#6a3d9a") +
    labs(x = "Year of the emission pulse", y = "SCC [USD/tCO2, discounted to the pulse year]",
         title = "SCC schedule under the background SAI deployment",
         subtitle = if (is.null(s)) NULL else s$discount_label) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
}


# -----------------------------------------------------------------------------
# Console report
# -----------------------------------------------------------------------------
svg_report <- function(x) {
  s <- x$summary
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat("Social value of SAI  |  ", s$discount_label,
      "  |  t <= ", s$t_max, " (", SVG_BASE_YEAR + s$t_max, ")\n", sep = "")
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
svg_main <- function(rate = NULL, t_max = SVG_T_MAX,
                     filters = list(cool_rate = SVG_REFERENCE_RC),
                     outdir = SVG_OUTPUT_DIR, plots = TRUE) {
  rep <- svg_report_table(t_max = t_max, filters = filters)
  print(attr(rep, "run"))
  cat("  variants:", paste(rep[, uniqueN(scenario)], "scenarios;"),
      "quantities:", paste(svg_quantities(rep), collapse = ", "), "\n")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  x <- svg_social_value(rep, rate = rate)
  print(x$disc)
  svg_report(x)
  fwrite(x$summary, file.path(outdir, "svg_summary.csv"))
  fwrite(x$components, file.path(outdir, "svg_damage_components.csv"))
  fwrite(svg_discount_table(x$disc), file.path(outdir, "svg_discount_curve.csv"))

  prof   <- svg_commitment_profile(x)
  actual <- svg_terminated_value(x)
  offset <- svg_offset_value(x)
  sched  <- svg_scc_schedule(x)
  hz     <- svg_commitment_hazard(x, h = SVG_HAZARD_MAIN, tc = 0L)
  fwrite(prof, file.path(outdir, "svg_commitment_profile.csv"))
  if (!is.null(hz)) {
    fwrite(hz, file.path(outdir, "svg_commitment_hazard.csv"))
    o <- attr(hz, "opt")
    cat(sprintf("\nWith the termination hazard (h = %s/yr, beta = %s, n = %s):\n",
                format(SVG_HAZARD_MAIN), format(SVG_BETA), format(SVG_P * SVG_Q)))
    cat(sprintf("  optimal commitment H* = %d years (%s), SV* = %.2f USD/tCO2 = %.1f%% of SCC\n",
                o$H_star, o$class, o$SV_star, 100 * o$SV_star / x$summary$scc_usd_tco2))
    cat(sprintf("  vs %.2f USD/tCO2 with no hazard and no planned end\n",
                x$summary$mb_sai_usd_tco2))
  }
  states <- svg_value_by_state(x)
  if (!is.null(states)) {
    fwrite(states, file.path(outdir, "social_value_by_state.csv"))
    cat("\nSV*(G | s,h) by starting state:\n")
    print(states[h == SVG_HAZARD_MAIN,
                 .(year, G = round(G, 3), M = round(M, 3), scc0 = round(scc0, 1),
                   deltaL0 = round(deltaL0, 1), H_star, H_class,
                   SV_star = round(SV_star, 2),
                   SV_over_scc = round(SV_star_over_scc, 4))])
  }

  if (!is.null(actual)) {
    fwrite(actual, file.path(outdir, "svg_terminated_value.csv"))
    cat("\nActual commitment experiments (marginal SAI terminated after T years):\n")
    print(merge(actual, prof[, .(commitment, theory_usd_tco2 = mb_usd_tco2)], by = "commitment"
                )[, .(T = commitment, actual_usd_tco2 = round(mb_usd_tco2, 1),
                      theory_usd_tco2 = round(theory_usd_tco2, 1),
                      actual_frac_scc = round(frac_scc, 4),
                      actual_over_theory = round(mb_usd_tco2 / theory_usd_tco2, 3))])
  }
  if (!is.null(offset)) {
    fwrite(offset, file.path(outdir, "svg_offset_value.csv"))
    cat("\nCarbon-offset experiments",
        "(CO2 removed at the pulse period, re-emitted after T years):\n")
    # Beside the SAI commitment of the same length when that family is present:
    # same pulse, same units, so the ratio is the like-for-like comparison.
    o <- if (is.null(actual)) copy(offset)[, sai_usd_tco2 := NA_real_] else
      merge(offset, actual[, .(commitment, sai_usd_tco2 = mb_usd_tco2)],
            by = "commitment", all.x = TRUE)
    print(o[, .(T = commitment, offset_usd_tco2 = round(mb_usd_tco2, 1),
                sai_usd_tco2 = round(sai_usd_tco2, 1),
                offset_frac_scc = round(frac_scc, 4),
                offset_over_sai = round(mb_usd_tco2 / sai_usd_tco2, 3))])
  }
  if (!is.null(sched)) {
    fwrite(sched, file.path(outdir, "svg_scc_schedule.csv"))
    cat("\nSCC schedule under the background SAI deployment:\n")
    print(sched[, .(delay, pulse_year, scc_usd_tco2 = round(scc_usd_tco2, 1),
                    scc_pv_usd_tco2 = round(scc_pv_usd_tco2, 1),
                    delay_value_usd_tco2 = round(delay_value_usd_tco2, 1))])
  }

  # Sensitivity of the ratio to the discounting assumption: the Ramsey
  # parameters one at a time around the default, and the constant rates the
  # earlier convention used, so the two conventions are read off one table.
  sens <- rbindlist(c(
    lapply(SVG_RHO_SENS, function(z)
      svg_social_value(rep, rate = svg_resolve_discount(rep = rep, rho = z,
                                                        mode = "ramsey"))$summary),
    lapply(SVG_ETA_SENS, function(z)
      svg_social_value(rep, rate = svg_resolve_discount(rep = rep, eta = z,
                                                        mode = "ramsey"))$summary),
    lapply(SVG_DISCOUNT_RATES_SENS, function(rr)
      svg_social_value(rep, rate = rr)$summary)))
  sens <- unique(sens, by = c("discount_mode", "discount_rho", "discount_eta",
                              "discount_rate"))
  fwrite(sens, file.path(outdir, "svg_discount_sensitivity.csv"))
  cat("\nDiscounting sensitivity:\n")
  print(sens[, .(mode = discount_mode, rho = discount_rho, eta = discount_eta,
                 r_at_pulse = round(discount_rate, 4),
                 scc_usd_tco2 = round(scc_usd_tco2, 1),
                 mb_sai_usd_tco2 = round(mb_sai_usd_tco2, 1),
                 svg = round(svg, 4))])

  commitment_sens <- svg_commitment_sensitivity(rep, max_T = SVG_FIG_MAX_T)
  fwrite(commitment_sens[, .(parameter, value, case, H, sv)],
         file.path(outdir, "svg_commitment_sensitivity.csv"))

  if (plots) {
    ggsave(file.path(outdir, "svg_commitment.png"),
           plot_svg_commitment(prof, actual = actual, delay = sched, hazard = hz,
                               offset = offset),
           width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
    ggsave(file.path(outdir, "svg_commitment_zoom.png"),
           plot_svg_commitment(prof, actual = actual, delay = sched, hazard = hz,
                               offset = offset, max_T = SVG_FIG_MAX_T),
           width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
    ggsave(file.path(outdir, "svg_commitment_panels.png"),
           plot_svg_commitment(prof, panels = TRUE), width = 8, height = 7, dpi = 200)
    ggsave(file.path(outdir, "svg_commitment_sensitivity.png"),
           plot_svg_commitment_sensitivity(commitment_sens),
           width = 9, height = 7, dpi = SVG_FIG_DPI)
    ggsave(file.path(outdir, "discount_curve.png"), plot_discount_curve(x$disc),
           width = 8, height = 5, dpi = SVG_FIG_DPI)
    if (!is.null(sched))
      ggsave(file.path(outdir, "scc_schedule.png"), plot_scc_schedule(sched, x$summary),
             width = 8, height = 5, dpi = 200)
    for (v in c("tatm", "damages", "forc_srm", "qsrm"))
      ggsave(file.path(outdir, paste0("diff_", v, ".png")),
             plot_gdx_diff(v, rep = rep, t_max = 480), width = 8, height = 5, dpi = 200)
  }
  cat("\nWrote CSVs", if (plots) " and figures" else "", " to ", outdir, "\n", sep = "")
  invisible(list(rep = rep, value = x, profile = prof, actual = actual,
                 offset = offset, schedule = sched,
                 commitment_sensitivity = commitment_sens))
}

# Runs on `Rscript Paper_SVG/svg_analysis.R` and on source(); the functions above
# stay available either way. options(svg.no.main = TRUE) loads them without running.
if (!isTRUE(getOption("svg.no.main"))) {
  svg_results <- svg_main()
}
