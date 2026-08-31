# =============================================================================
# svg_validate.R
#
# Validation suite for the social-value pipeline (implementation spec section 8).
# Runs as one script and prints a pass/fail table:
#
#   Rscript Paper_SVG/svg_validate.R
#
# Three kinds of check:
#   * unit tests on the closed-form pieces (units, termination loss, survival);
#   * an analytic regression -- the integrator is fed synthetic constant-growth
#     schedules whose SV has a closed form, which is what ties this code to the
#     simplified calculations the project started from;
#   * accounting identities on the real model output.
#
# Exits non-zero if anything fails, so it can gate a production run.
# =============================================================================

options(svg.no.main = TRUE)
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "svg_analysis.R"))

.tests <- list()
svg_check <- function(name, got, want = TRUE, tol = 1e-8, note = "") {
  ok <- if (is.logical(want) && length(want) == 1L && is.logical(got))
          isTRUE(got) == isTRUE(want)
        else all(abs(got - want) <= tol * pmax(1, abs(want)))
  .tests[[length(.tests) + 1L]] <<- data.table(
    test = name, ok = ok,
    got = paste(signif(as.numeric(got), 8), collapse = ", "),
    want = paste(signif(as.numeric(want), 8), collapse = ", "), note = note)
  invisible(ok)
}
svg_skip <- function(name, note) {
  .tests[[length(.tests) + 1L]] <<- data.table(test = name, ok = NA, got = "-",
                                               want = "-", note = note)
}

r <- SVG_DISCOUNT_RATE

# -----------------------------------------------------------------------------
# 8.1 Unit tests
# -----------------------------------------------------------------------------
svg_check("units: TgSO2 -> TgS", svg_tgso2_to_tgs(64), 32, note = "molar 64/32")
svg_check("units: round trip", svg_tgs_to_tgso2(svg_tgso2_to_tgs(3.7)), 3.7)

svg_check("termination: DeltaL = 0 at M = 0",
          svg_marginal_termination_loss(0, 100), 0)
svg_check("termination: DeltaL/SCC = 8 at M = 1",
          svg_marginal_termination_loss(1, 100) / 100, 8,
          note = "beta = 4, p = 2, q = 1")
svg_check("termination: linear in M at n = 2",
          svg_marginal_termination_loss(2.5, 100) / svg_marginal_termination_loss(1, 100),
          2.5)
svg_check("termination: beta scales DeltaL",
          svg_marginal_termination_loss(1, 100, beta = 8) /
            svg_marginal_termination_loss(1, 100, beta = 4), 2)

svg_check("hazard: survival exp(-1) at h*T = 1 (h = 0.01)", svg_survival(100, 0.01), exp(-1))
svg_check("hazard: survival exp(-1) at h*T = 1 (h = 0.001)", svg_survival(1000, 0.001), exp(-1))
svg_check("hazard: survival is 1 at h = 0", svg_survival(500, 0), 1)

# h = 0 must remove every termination charge and every survival weight, leaving
# the plain discounted sum: this is what makes the hazard curve comparable with
# the deterministic commitment profile.
.f <- rep(2, 50); .l <- rep(1000, 50)
.sv0 <- svg_commitment_value(.f, .l, h = 0, rate = r)
svg_check("hazard: h = 0 gives the plain discounted sum",
          .sv0$sv[nrow(.sv0)], sum(.f * (1 + r)^-(0:49)))
svg_check("hazard: h = 0 charges no termination", max(abs(.sv0$termination)), 0)
svg_check("commitment: SV(0) = 0", .sv0$sv[1], 0)

# Deposition separability: the sulfate-mortality term is a component of the same
# damage total, so removing it must move the NPV by exactly its own NPV.
svg_check("deposition: SV is linear in the VLL component",
          svg_commitment_value(.f - 0.5, .l, h = 0, rate = r)$sv[51],
          .sv0$sv[51] - sum(rep(0.5, 50) * (1 + r)^-(0:49)))

# -----------------------------------------------------------------------------
# 8.1b Discount curve (svg_discounting.R)
# -----------------------------------------------------------------------------
# A flat curve must be indistinguishable from the scalar it replaces: this is
# what guarantees the constant-rate results did not move when the curve was
# introduced.
.flat <- svg_rate_curve(1:200, rep(r, 200))
svg_check("discount: flat curve equals the closed form",
          max(abs(svg_df(.flat, 1:150, 20) - (1 + r)^-((1:150) - 20))), 0, tol = 1e-12)
svg_check("discount: constant curve equals the closed form",
          svg_df(svg_constant_curve(r), 77, 5), (1 + r)^-72)
svg_check("discount: DF is 1 at the reference period", svg_df(.flat, 30, 30), 1)

# Ramsey construction on a synthetic economy growing at a constant gc: the rate
# is rho + eta*gc everywhere and the factor is the corresponding geometric one.
.gc <- 0.015; .rho <- SVG_RHO; .eta <- SVG_ETA
.ram <- svg_ramsey_curve(1:300, 100 * (1 + .gc)^(1:300), rho = .rho, eta = .eta)
svg_check("ramsey: r = rho + eta*g on constant growth",
          max(abs(.ram$r - (.rho + .eta * .gc))), 0, tol = 1e-12)
svg_check("ramsey: factor matches the implied constant rate",
          svg_df(.ram, 101, 1), (1 + .rho + .eta * .gc)^-100, tol = 1e-10)
svg_check("ramsey: zero growth leaves r = rho",
          svg_rate_scalar(svg_ramsey_curve(1:50, rep(100, 50)), 10), .rho)

# Chaining: a declining curve must compose, DF(t2|t0) = DF(t2|t1)*DF(t1|t0).
.dec <- svg_rate_curve(1:400, seq(0.04, 0.005, length.out = 400))
svg_check("discount: factors chain across periods",
          svg_df(.dec, 300, 5), svg_df(.dec, 300, 120) * svg_df(.dec, 120, 5),
          tol = 1e-12)
svg_check("discount: a declining curve discounts the tail less than its start rate",
          svg_df(.dec, 400, 1) > (1 + .dec$r[1])^-399)
# Rates are read where they are asked for, not averaged.
svg_check("discount: r(t) recovered from the curve",
          svg_rate_at(.dec, c(1, 200, 400)), .dec$r[c(1, 200, 400)])

# The hazard integrator must see the curve at the right place on the time axis:
# the same commitment started later is discounted by the later, lower rates.
.f1 <- rep(1, 100)
svg_check("hazard: t0 shifts the commitment along a declining curve",
          svg_commitment_value(.f1, rep(0, 100), h = 0, rate = .dec, t0 = 200)$sv[101] >
            svg_commitment_value(.f1, rep(0, 100), h = 0, rate = .dec, t0 = 1)$sv[101])
svg_check("hazard: t0 is irrelevant under a constant rate",
          svg_commitment_value(.f1, rep(0, 100), h = 0, rate = .flat, t0 = 60)$sv[101],
          svg_commitment_value(.f1, rep(0, 100), h = 0, rate = r)$sv[101], tol = 1e-12)

# -----------------------------------------------------------------------------
# 8.2 Analytic regression: constant-growth schedules
# -----------------------------------------------------------------------------
# With F(k) = F0*(1+g)^k, DeltaL(k) = L0*(1+g)^k, discounting (1+r)^-k, survival
# exp(-hk) and per-interval event probability exp(-hk)*(1-exp(-h)):
#   SV(H) = [F0 - L0*(1-exp(-h))] * (1-q^H)/(1-q),   q = (1+g)*exp(-h)/(1+r)
.g <- 0.01; .h <- 0.001; .F0 <- 3; .L0 <- 500; .H <- 200
.k <- 0:(.H - 1)
.sv <- svg_commitment_value(.F0 * (1 + .g)^.k, .L0 * (1 + .g)^.k, h = .h, rate = r)
.q <- (1 + .g) * exp(-.h) / (1 + r)
.closed <- (.F0 - .L0 * (1 - exp(-.h))) * (1 - .q^.H) / (1 - .q)
svg_check("analytic: SV(H) matches the closed form", .sv$sv[.H + 1], .closed, tol = 1e-10)

# Same, with no hazard: the geometric series with q = (1+g)/(1+r).
.sv2 <- svg_commitment_value(.F0 * (1 + .g)^.k, .L0 * (1 + .g)^.k, h = 0, rate = r)
.q2 <- (1 + .g) / (1 + r)
svg_check("analytic: no-hazard SV(H) matches the closed form",
          .sv2$sv[.H + 1], .F0 * (1 - .q2^.H) / (1 - .q2), tol = 1e-10)

# Optimiser: a flow that turns negative at a known period must give that H*.
.flow <- c(rep(1, 30), rep(-1, 70))
.opt <- svg_optimize_commitment(svg_commitment_value(.flow, rep(0, 100), h = 0, rate = r))
svg_check("optimiser: interior optimum found", .opt$H_star, 30)
svg_check("optimiser: classified finite", .opt$class == "finite")
.opt0 <- svg_optimize_commitment(svg_commitment_value(rep(-1, 50), rep(0, 50), h = 0, rate = r))
svg_check("optimiser: H* = 0 when the flow starts negative", .opt0$H_star, 0)
svg_check("optimiser: classified zero", .opt0$class == "zero")
.optinf <- svg_optimize_commitment(svg_commitment_value(rep(1, 400), rep(0, 400), h = 0, rate = r))
svg_check("optimiser: converged tail classified infinite", .optinf$class == "infinite")

# -----------------------------------------------------------------------------
# 8.3 SCC identity on synthetic input
# -----------------------------------------------------------------------------
# SCC(t) = A*(1+g)^t implies MD = r*SCC - dSCC/dt = A*(r-g)*(1+g)^t.
.t <- 1:300
.sched <- data.table(pulse_t = .t, scc_usd_tco2 = 100 * (1 + .g)^.t)
.ann <- svg_scc_annual(.sched, .t, rate = r)
svg_check("SCC identity: MD = r*SCC - dSCC/dt on synthetic input",
          max(abs(.ann$md[1:299] - (r - .g) * 100 * (1 + .g)^.t[1:299])), 0, tol = 1e-6,
          note = "exact for a geometric SCC")

# -----------------------------------------------------------------------------
# Model-based checks (skipped when no complete run is present)
# -----------------------------------------------------------------------------
.runs <- tryCatch(svg_gdx_index()[experiment == SVG_EXPERIMENT], error = function(e) NULL)
.ok_run <- if (is.null(.runs)) NULL else .runs[vapply(path, svg_gdx_complete, logical(1))]

if (is.null(.ok_run) || !nrow(.ok_run)) {
  svg_skip("model checks", "no complete svg gdx found; run Paper_SVG/svg_scenarios.R")
} else {
  .run <- structure(.ok_run[order(-mtime)][1, path], names = SVG_EXPERIMENT,
                    scenario = as.list(.ok_run[order(-mtime)][1, ..SVG_SCENARIO_KEYS]),
                    class = c("svg_run", "character"))
  .rep <- svg_report_table(.run, quiet = TRUE)
  .x   <- svg_social_value(.rep)
  .s   <- .x$summary

  # eq_damtot: svg_damage_components() already stops if this fails, so reaching
  # here is the pass.
  svg_check("model: damage components reproduce DAMAGES", TRUE,
            note = basename(unclass(.run)))
  svg_check("model: pulse is one period", nrow(.x$pulses[scenario == SVG_GHGPULSE]), 1)
  svg_check("model: G >= 0 throughout", .rep[vrn == SVG_BASE, min(qsrm)] >= -1e-9)
  svg_check("model: masked warming M >= 0 throughout",
            .rep[vrn == SVG_BASE, min(tatm_ghg - tatm)] >= -1e-9)

  # The curve the run is analysed under is the configured one, built from the
  # model's own ypc path, and it really does decline to rho.
  svg_check("discount: run uses the configured mode", .x$disc$mode == SVG_DISCOUNT_MODE,
            note = .s$discount_label)
  if (identical(SVG_DISCOUNT_MODE, "ramsey")) {
    svg_check("ramsey: r(t) settles at rho once growth stops",
              .x$disc$r[length(.x$disc$r)], SVG_RHO, tol = 1e-12)
    # r(t) edges up for the first decade -- population growth falls from 2020
    # while GDP growth holds flat to 2030, so per-capita growth rises before it
    # rolls over -- and is non-increasing from that peak on.
    .rr <- .x$disc$r
    .pk <- which.max(.rr)
    svg_check("ramsey: r(t) is non-increasing after its peak",
              max(diff(.rr[.pk:length(.rr)])) <= 1e-12,
              note = sprintf("peak %.3f%% in %d, then %.3f%%", 100 * .rr[.pk],
                             SVG_BASE_YEAR + .x$disc$t[.pk], 100 * .rr[length(.rr)]))
    # A declining curve weights the tail more, so both marginal quantities are
    # larger than under the flat rate the curve replaced.
    .const <- svg_social_value(.rep, rate = SVG_DISCOUNT_RATE)$summary
    svg_check("ramsey: SCC exceeds the constant-rate SCC",
              .s$scc_usd_tco2 > .const$scc_usd_tco2,
              note = sprintf("%.1f vs %.1f USD/tCO2", .s$scc_usd_tco2,
                             .const$scc_usd_tco2))
  }

  # Absolute vs normalised value must be the same statement.
  svg_check("accounting: SV = normalised SV * SCC(t0)",
            .s$svg * .s$scc_usd_tco2, .s$mb_sai_usd_tco2, tol = 1e-10)

  # At the M = 0 state the marginal intervention adds nothing to a programme that
  # is not yet masking anything, so its termination externality is zero even at h > 0.
  .st <- svg_states(.x)
  if (!is.null(.st)) {
    svg_check("accounting: DeltaL = 0 at the M = 0 state",
              .st[tc == 0, deltaL0], 0, tol = 1e-9)
    svg_check("model: eta = cos(theta) at M = 0", .st[tc == 0, eta],
              cos(svg_scalar(.run, "srm_angle") * pi / 180), tol = 1e-6)

    # h -> 0 must reproduce the deterministic marginal benefit of the never-
    # terminated intervention: the two calculations meet where the hazard vanishes.
    .hz0 <- svg_commitment_hazard(.x, h = 0, tc = 0L)
    svg_check("accounting: h = 0 reproduces the deterministic MB(SAI)",
              .hz0$sv[nrow(.hz0)], .s$mb_sai_usd_tco2, tol = 1e-8)

    # Reported optima must really be optima on the grid.
    .v <- svg_value_by_state(.x, hazards = SVG_HAZARD_MAIN)
    .bad <- 0L
    for (i in seq_len(nrow(.v))) {
      sv <- svg_commitment_hazard(.x, h = .v$h[i], tc = .v$tc[i])
      j <- which(sv$H == .v$H_star[i])
      nb <- sv$sv[intersect(c(j - 1, j + 1), seq_len(nrow(sv)))]
      if (any(nb > .v$SV_star[i] + 1e-9)) .bad <- .bad + 1L
    }
    svg_check("optimiser: no neighbouring horizon beats H*", .bad, 0)
  }

  # Carbon-offset family. Two invariants that hold whatever the calibration:
  # holding the tonne out of the atmosphere for longer can only be worth more,
  # and even permanent storage cannot beat the SCC of the pulse it removes --
  # damages are convex, so removing a tonne saves slightly less than adding one
  # costs. Both would break if the removal and the re-emission were mis-timed.
  .off <- svg_offset_value(.x)
  if (is.null(.off)) {
    svg_skip("model: carbon-offset checks", "no ghgstore family in this gdx")
  } else {
    svg_check("model: carbon-offset value is non-decreasing in T",
              min(diff(.off$mb_usd_tco2)) >= -1e-6,
              note = sprintf("%d commitment lengths", nrow(.off)))
    svg_check("model: carbon-offset value stays within the SCC",
              max(.off$mb_usd_tco2) <= .s$scc_usd_tco2 * (1 + 1e-6),
              note = sprintf("max %.1f vs SCC %.1f USD/tCO2",
                             max(.off$mb_usd_tco2), .s$scc_usd_tco2))
  }

  # Tail convergence: the reported values must not depend on where the horizon
  # is cut (spec 3.3), tested against a cut just inside SVG_T_MAX. Under the
  # Ramsey curve the rate settles at rho = 0.8% from 2200, so the tail is worth
  # far more than it was at a flat 3% and the test horizon has to sit near the
  # end of the simulation: a cut at t = 480 (2499) still moves the SCC by ~2%,
  # which is why SVG_T_MAX is the full FAIR horizon and not something shorter.
  .short <- svg_social_value(svg_report_table(.run, t_max = SVG_T_TAIL_TEST,
                                              quiet = TRUE))
  .lab <- sprintf(" (t_max %d vs %d)", SVG_T_TAIL_TEST, SVG_T_MAX)
  svg_check(paste0("tail: SCC(t0) converged", .lab),
            abs(.short$summary$scc_usd_tco2 / .s$scc_usd_tco2 - 1) <= SVG_SCC_TOL,
            note = sprintf("%.4f%% change", 100 * (.short$summary$scc_usd_tco2 / .s$scc_usd_tco2 - 1)))
  svg_check(paste0("tail: SVG converged", .lab),
            abs(.short$summary$svg / .s$svg - 1) <= SVG_SV_TOL,
            note = sprintf("%.4f%% change", 100 * (.short$summary$svg / .s$svg - 1)))
}

# -----------------------------------------------------------------------------
res <- rbindlist(.tests)
cat("\n", strrep("-", 78), "\nsvg validation suite\n", strrep("-", 78), "\n", sep = "")
print(res[, .(test, result = fifelse(is.na(ok), "SKIP", fifelse(ok, "pass", "FAIL")),
              got, want, note)], nrows = 100)
n_fail <- res[ok == FALSE, .N]; n_skip <- res[is.na(ok), .N]
cat(sprintf("\n%d passed, %d failed, %d skipped\n", res[ok == TRUE, .N], n_fail, n_skip))
if (n_fail > 0) quit(save = "no", status = 1)
