# =============================================================================
# svg_figures.R
#
# Cross-pathway products: the social value of the marginal intervention as a
# function of the background deployment it is added to, SV*(G | s,h), and the
# optimal commitment H*(G | s,h) that goes with it (spec sections 2.13, 5.13).
#
#   Rscript Paper_SVG/svg_figures.R              # every complete svg run in Results/
#
# One GAMS run per cooling ramp s; each already carries a whole grid of starting
# states, so the curve below needs no extra solves. Values are per tonne of the
# pulse the intervention offsets, and normalised by the SCC of the SAME pathway
# at the SAME state -- never by a universal reference SCC (spec section 11.1).
# =============================================================================

options(svg.no.main = TRUE)
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "svg_analysis.R"))

# Every complete svg gdx in the results folder, newest first per cooling ramp.
svg_all_runs <- function(folder = SVG_RESULTS, filters = list()) {
  idx <- svg_gdx_index(folder)[experiment == SVG_EXPERIMENT]
  for (nm in names(filters)) idx <- idx[get(nm) == as.character(filters[[nm]])]
  idx <- idx[vapply(path, svg_gdx_complete, logical(1))]
  if (!nrow(idx)) stop("No complete svg run found in ", folder,
                       ". Run Paper_SVG/svg_scenarios.R first.", call. = FALSE)
  setorder(idx, -mtime)
  idx <- idx[, .SD[1], by = cool_rate]          # newest run per ramp
  setorder(idx, cool_rate)
  lapply(seq_len(nrow(idx)), function(i)
    structure(idx$path[i], names = SVG_EXPERIMENT,
              scenario = as.list(idx[i, ..SVG_SCENARIO_KEYS]),
              class = c("svg_run", "character")))
}

# State table across every ramp, tagged with the ramp rate. `rate` is NULL for
# the configured discount curve (Ramsey, built from each run's own ypc path --
# the same path in every run, since it is exogenous).
svg_states_all <- function(runs = NULL, hazards = SVG_HAZARDS, rate = NULL) {
  if (is.null(runs)) runs <- svg_all_runs()
  rbindlist(lapply(runs, function(run) {
    rc <- as.numeric(attr(run, "scenario")$cool_rate)
    x  <- svg_social_value(svg_report_table(run, quiet = TRUE), rate = rate)
    v  <- svg_value_by_state(x, hazards = hazards)
    if (is.null(v)) return(NULL)
    v[, `:=`(rc = rc, s = svg_rc_to_s(rc), scc_ref = x$summary$scc_usd_tco2,
             svg_nohazard = x$summary$svg,
             disc_label = x$summary$discount_label)][]
  }), fill = TRUE)
}

# SV*(G | s,h): value of the marginal intervention against the deployment it is
# added to. One line per ramp; the same G under a different future pathway is a
# different state, so the ramps are NOT pooled into one curve (spec section 6.9).
#   x : "G" (deployment, as the spec asks) or "M" (masked warming). The ramp
#       plateaus after end_rampup, so many later states share one G and the
#       G-curve turns vertical there; M keeps rising and separates them.
plot_sv_curve <- function(states, h = SVG_HAZARD_MAIN, normalised = FALSE, x = "G") {
  h_val <- h            # the argument would otherwise shadow the h column
  d <- states[h == h_val]
  yv <- if (normalised) "SV_star_over_scc" else "SV_star"
  xlab <- if (x == "G") "Background SAI deployment G at the starting state [TgS/yr]"
          else "Masked warming M at the starting state [degC]"
  ggplot(d, aes(get(x), get(yv), colour = factor(s), group = s)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    labs(x = xlab,
         y = if (normalised) "SV* / SCC(t0)  [-]" else "SV* [USD/tCO2]",
         colour = "cooling ramp s\n[degC/decade]",
         title = "Social value of marginal SAI along the background deployment",
         subtitle = sprintf("h = %s/yr, beta = %s, %s\neach point is one starting state, normalised by its own pathway's SCC",
                            format(h), format(SVG_BETA),
                            d$disc_label[1] %||% "")) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
}

# H*(G | s,h). Horizons pinned at the simulation end are drawn as open symbols so
# an "infinite" commitment is never confused with a genuinely interior optimum.
plot_h_star <- function(states, h = SVG_HAZARD_MAIN) {
  h_val <- h
  d <- copy(states[h == h_val])
  d[, kind := fifelse(H_class == "finite", "finite optimum",
              fifelse(H_class == "zero", "H* = 0 (no commitment)",
              fifelse(H_class == "infinite", "unbounded (tail converged)", "boundary (unresolved)")))]
  ggplot(d, aes(G, H_star, colour = factor(s), group = s)) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = kind), size = 2.4, fill = "white") +
    scale_shape_manual(values = c("finite optimum" = 16, "H* = 0 (no commitment)" = 4,
                                  "unbounded (tail converged)" = 21,
                                  "boundary (unresolved)" = 8)) +
    labs(x = "Background SAI deployment G at the starting state [TgS/yr]",
         y = "H* [years]", colour = "cooling ramp s\n[degC/decade]", shape = NULL,
         title = "Optimal commitment horizon along the background deployment",
         subtitle = sprintf("h = %s/yr, beta = %s", format(h), format(SVG_BETA))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
}

# Background pathway diagnostics: deployment, masked warming, temperature, SCC.
plot_backgrounds <- function(runs = NULL) {
  if (is.null(runs)) runs <- svg_all_runs()
  d <- rbindlist(lapply(runs, function(run) {
    rc <- as.numeric(attr(run, "scenario")$cool_rate)
    rp <- svg_report_table(run, quiet = TRUE)[vrn == SVG_BASE]
    data.table(s = svg_rc_to_s(rc), year = SVG_BASE_YEAR + rp$t,
               G = rp$qsrm, M = rp$tatm_ghg - rp$tatm, T = rp$tatm)
  }))
  long <- melt(d[year <= 2300], id.vars = c("s", "year"),
               measure.vars = c("G", "M", "T"), variable.name = "q")
  levels(long$q) <- c("G [TgS/yr]", "masked warming M [degC]", "temperature T [degC]")
  ggplot(long, aes(year, value, colour = factor(s))) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~q, ncol = 1, scales = "free_y", strip.position = "left") +
    labs(x = "Year", y = NULL, colour = "cooling ramp s\n[degC/decade]",
         title = "Background SAI pathways") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          panel.grid.minor = element_blank())
}

svg_figures_main <- function(outdir = SVG_OUTPUT_DIR, h = SVG_HAZARD_MAIN) {
  runs <- svg_all_runs()
  cat("svg runs found:\n"); invisible(lapply(runs, print))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  states <- svg_states_all(runs)
  fwrite(states, file.path(outdir, "social_value_by_state.csv"))
  h_val <- h
  cat("\nSV*(G | s, h =", format(h), ") by ramp and starting state:\n")
  print(states[h == h_val, .(s, year, G = round(G, 3), M = round(M, 3),
                           scc0 = round(scc0, 1), H_star, H_class,
                           SV_star = round(SV_star, 2),
                           SV_over_scc = round(SV_star_over_scc, 3))])

  ggsave(file.path(outdir, "sv_star_by_deployment.png"), plot_sv_curve(states, h),
         width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  ggsave(file.path(outdir, "sv_star_by_deployment_normalised.png"),
         plot_sv_curve(states, h, normalised = TRUE),
         width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  ggsave(file.path(outdir, "sv_star_by_maskedwarming.png"),
         plot_sv_curve(states, h, x = "M"),
         width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  ggsave(file.path(outdir, "sv_star_by_maskedwarming_normalised.png"),
         plot_sv_curve(states, h, normalised = TRUE, x = "M"),
         width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  ggsave(file.path(outdir, "h_star_by_deployment.png"), plot_h_star(states, h),
         width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  ggsave(file.path(outdir, "background_pathways.png"), plot_backgrounds(runs),
         width = 8, height = 7, dpi = SVG_FIG_DPI)
  # The curve every value above is integrated against; identical across ramps,
  # since ypc is exogenous.
  ggsave(file.path(outdir, "discount_curve.png"),
         plot_discount_curve(svg_resolve_discount(
           rep = svg_report_table(runs[[1]], quiet = TRUE))),
         width = 8, height = 5, dpi = SVG_FIG_DPI)

  # One commitment figure per ramp, each with its own hazard optimum marked.
  for (run in runs) {
    rc <- attr(run, "scenario")$cool_rate
    x  <- svg_social_value(svg_report_table(run, quiet = TRUE))
    ggsave(file.path(outdir, paste0("svg_commitment_RC", rc, ".png")),
           plot_svg_commitment(svg_commitment_profile(x),
                               actual = svg_terminated_value(x),
                               offset = svg_offset_value(x),
                               delay  = svg_scc_schedule(x),
                               hazard = svg_commitment_hazard(x, h = h, tc = 0L)),
           width = SVG_FIG_WIDTH, height = SVG_FIG_HEIGHT, dpi = SVG_FIG_DPI)
  }
  cat("\nFigures written to", outdir, "\n")
  invisible(states)
}

if (!isTRUE(getOption("svg.no.figures"))) svg_figures_result <- svg_figures_main()
