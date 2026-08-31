# =============================================================================
# svg_scenarios.R
#
# Runs the background cooling-ramp grid of the social-value experiment: one GAMS
# call per ramp, each producing the single svg gdx that Paper_SVG/svg_analysis.R
# reads. Follows the runner pattern of Paper_SAI/Sampling/Run_montecarlo.R --
# derive the expected result file name, skip what already exists, capture the
# solve status, and never fail silently.
#
#   Rscript Paper_SVG/svg_scenarios.R            # full grid (SVG_COOLING_RC)
#   Rscript Paper_SVG/svg_scenarios.R mvp        # SVG_COOLING_RC_MVP only
#   Rscript Paper_SVG/svg_scenarios.R 0 20       # just these ramps, in deg/millennium
#   Rscript Paper_SVG/svg_scenarios.R plan       # print the plan, run nothing
#   Rscript Paper_SVG/svg_scenarios.R force      # re-run even if the gdx exists
#
# rc = 0 is the no-background-SAI reference: the only aerosol in the run is the
# marginal pulse itself, which is the G = 0 end of the SV*(G) curve.
#
# A command-line --rate_of_cooling survives the $setglobal default at the top of
# experiments/svg.gms (verified: GAMS gives double-dash parameters precedence),
# which is what makes one experiment file serve the whole grid.
# =============================================================================

options(svg.no.main = TRUE)
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "svg_analysis.R"))

# Expected gdx name for a ramp, mirroring the execute_unload in experiments/svg.gms.
svg_expected_gdx <- function(rc, cfg = SVG_RUN_DEFAULTS, ramp = SVG_RAMP) {
  file.path(SVG_RESULTS,
            paste0(cfg$rcp, "_EXPsvg_GAS", cfg$gas, "_ECS", cfg$ecs, "_TCR", cfg$tcr,
                   "_PT", cfg$pulse_time, "_RC", rc,
                   "_EC", ramp$end_rampdown, "_BC", ramp$start_rampup, ".gdx"))
}

# GAMS command for one ramp. Every scenario axis is passed explicitly so the grid
# is reproducible from svg_config.R alone.
svg_gams_command <- function(rc, cfg = SVG_RUN_DEFAULTS, ramp = SVG_RAMP) {
  gams <- Sys.which("gams")
  if (!nzchar(gams)) stop("GAMS not found on the PATH.", call. = FALSE)
  paste0("\"", gams, "\" FAIR.gms --experiment=svg",
         " --gas=", cfg$gas, " --rcp=", cfg$rcp,
         " --ecs=", cfg$ecs, " --tcr=", cfg$tcr,
         " --pulse_time=", cfg$pulse_time,
         " --rate_of_cooling=", rc,
         " --start_rampup=", ramp$start_rampup, " --end_rampup=", ramp$end_rampup,
         " --start_rampdown=", ramp$start_rampdown, " --end_rampdown=", ramp$end_rampdown,
         " lo=2")
}

# Run the grid. Returns the run manifest (spec section 13.2) and writes it next
# to the other outputs.
svg_run_grid <- function(rc_grid = SVG_COOLING_RC, skip_existing = TRUE,
                         dry_run = FALSE, outdir = SVG_OUTPUT_DIR) {
  man <- data.table(rc = rc_grid, s_degC_per_decade = svg_rc_to_s(rc_grid))
  man[, gdx := vapply(rc, svg_expected_gdx, character(1))]
  # A gdx counts as done only if it holds every scenario family: files written by
  # an earlier version of svg.gms exist but lack whole variants.
  man[, exists_before := vapply(gdx, svg_gdx_complete, logical(1))]
  man[, command := vapply(rc, svg_gams_command, character(1))]
  todo <- if (skip_existing) man[exists_before == FALSE] else man

  cat("svg scenario grid:", nrow(man), "ramp(s), s =",
      paste(man$s_degC_per_decade, collapse = ", "), "degC/decade\n")
  cat("  to run:", nrow(todo), "| already present:", nrow(man) - nrow(todo), "\n")
  cat("  each run is ~153 GAMS solves (~2.5 min); estimated",
      round(nrow(todo) * 2.5), "min total\n")
  if (dry_run) {
    cat("\nplan only, nothing executed. Commands:\n")
    cat(paste0("  ", todo$command, collapse = "\n"), "\n")
    return(invisible(man[]))
  }

  # FAIR.gms resolves experiments/ and Results/ relative to the working directory.
  old <- setwd(SVG_ROOT); on.exit(setwd(old), add = TRUE)
  man[, `:=`(status = NA_integer_, elapsed_s = NA_real_, produced = FALSE)]
  for (i in seq_len(nrow(man))) {
    if (skip_existing && man$exists_before[i]) {
      cat("[skip]", basename(man$gdx[i]), "\n"); man[i, produced := TRUE]; next
    }
    cat("[run ] rc =", man$rc[i], "(s =", man$s_degC_per_decade[i], ")... ")
    t0 <- Sys.time()
    st <- suppressWarnings(system(man$command[i], intern = FALSE))
    man[i, `:=`(status = st, elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                # `gdx` is already the i-th row here: indexing it again read
                # past the end and reported every ramp after the first as failed.
                produced = file.exists(gdx))]
    cat(if (man$produced[i]) "ok" else "FAILED",
        sprintf(" (%.0f s, exit %d)\n", man$elapsed_s[i], st))
  }
  bad <- man[produced == FALSE]
  if (nrow(bad))
    warning(nrow(bad), " ramp(s) produced no gdx: rc = ", paste(bad$rc, collapse = ", "),
            ". Check FAIR.lst.", call. = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fwrite(man[, .(rc, s_degC_per_decade, gdx = basename(gdx), exists_before,
                 status, elapsed_s, produced)],
         file.path(outdir, "run_manifest.csv"))
  cat("\nManifest written to", file.path(outdir, "run_manifest.csv"), "\n")
  invisible(man[])
}

if (!isTRUE(getOption("svg.no.grid"))) {
  .args <- commandArgs(TRUE)
  # Bare numbers select the ramps to run, so a single ramp can be added or
  # repaired without re-running the whole grid.
  .pick <- suppressWarnings(as.numeric(.args))
  .pick <- .pick[!is.na(.pick)]
  .rc   <- if (length(.pick)) .pick
           else if ("mvp" %in% .args) SVG_COOLING_RC_MVP else SVG_COOLING_RC
  svg_grid_manifest <- svg_run_grid(.rc,
                                    skip_existing = !("force" %in% .args),
                                    dry_run = "plan" %in% .args)
}
