# =============================================================================
# svg_config.R
#
# SINGLE CONTROL FILE for the social-value-of-SAI pipeline, in the same spirit as
# Paper_SAI/all_parameters.R: every calibration choice, grid and tolerance lives
# here with a one-line explanation, and nothing scientific is buried as a literal
# in the analysis code.
#
# Sourced by svg_analysis.R (which defines the file locations), so the objects
# below are available to every module. Base R only; sourcing has no side effects.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. DISCOUNTING AND HORIZON
# -----------------------------------------------------------------------------
# The GAMS model carries no discounting of its own: the whole curve is applied in
# post-processing, by svg_discounting.R.
#
# "ramsey"   r(t) = SVG_RHO + SVG_ETA * g(t), with g(t) the per-capita growth
#            rate of the economy taken from the model's own ypc path. Declining:
#            about 3%/yr today, SVG_RHO once growth has run out (2200 in this
#            calibration).
# "constant" a flat SVG_DISCOUNT_RATE, the earlier convention.
SVG_DISCOUNT_MODE <- "ramsey"

SVG_RHO <- 0.008   # pure rate of time preference [1/yr]
SVG_ETA <- 1.45    # elasticity of marginal utility of consumption [-]

SVG_DISCOUNT_RATE <- 0.03    # the flat rate of the "constant" mode, and the
                             # comparison point the sensitivity table reports
SVG_DISCOUNT_RATES_SENS <- c(0.01, 0.02, 0.03, 0.05)   # constant-rate sensitivity grid

# Ramsey sensitivity: one parameter moved at a time around (SVG_RHO, SVG_ETA).
SVG_RHO_SENS <- c(0.001, 0.008, 0.02)
SVG_ETA_SENS <- c(1.0, 1.45, 2.0)

SVG_T_MAX     <- 1000L   # last model period kept: the full FAIR horizon, which the
                         # longest pulse delays (tc = 500 -> emission in 2524) need.
                         # Under the Ramsey curve the rate settles at SVG_RHO from
                         # 2200, so the tail is worth much more than it was at a
                         # flat 3%: a cut at t = 480 still moves the SCC by ~2%,
                         # while the last 100 periods of the full horizon are worth
                         # 0.03% of it. Nothing shorter than this is safe.
SVG_T_TAIL_TEST <- 900L  # horizon the convergence check in svg_validate.R cuts at
SVG_BASE_YEAR <- 2019L   # calendar year = SVG_BASE_YEAR + t  (t = 1 is 2020)

# Discounting convention within a period. "start" applies (1+r)^-(t-t0), which is
# what the SCC and the deterministic commitment profile already use, so the h = 0
# case reproduces them exactly. "mid" is the midpoint form of spec section 2.11.
SVG_DISCOUNT_AT <- "start"


# -----------------------------------------------------------------------------
# 2. TERMINATION HAZARD AND LOSS  (spec sections 1.1, 2.8)
# -----------------------------------------------------------------------------
# Constant Poisson hazard of an abrupt, unplanned collapse of the SAI programme
# while the marginal intervention is active. One-in-a-millennium and
# one-in-a-century mean recurrence times. There is no separate failure hazard,
# and no termination cost is charged at the PLANNED end of the commitment.
SVG_HAZARDS <- c(0.001, 0.01)

# Headline hazard: the one used for the SV-vs-commitment figure and the reported
# optimum. The other value in SVG_HAZARDS is still carried through the state
# table as a sensitivity.
SVG_HAZARD_MAIN <- 0.001

# Reduced-form termination damage. With p the convexity of damages in the rebound
# velocity and q the scaling of that velocity with masked warming, n = p*q and
#   L_total(M,t) = beta * SCC(t) * (M_ref/alpha) * (M/M_ref)^n
#   DeltaL(t)    = alpha * dL_total/dM = beta * n * SCC(t) * (M/M_ref)^(n-1)
# so DeltaL is free of alpha, is zero at M = 0, and equals 8*SCC at M = 1 degC
# under the central calibration.
SVG_BETA  <- 4.0    # termination damage at M_ref, in multiples of the SCC benchmark
SVG_P     <- 2.0    # damage convexity in rebound velocity
SVG_Q     <- 1.0    # rebound velocity scaling in masked warming
SVG_M_REF <- 1.0    # reference masked warming [degC]

# Sensitivity grids (post-processing only, no extra model runs).
SVG_BETA_SENS <- c(2, 4, 8)
SVG_P_SENS    <- c(1.5, 2, 3)


# -----------------------------------------------------------------------------
# 3. BACKGROUND COOLING-RAMP GRID  (spec section 3.1)
# -----------------------------------------------------------------------------
# GAMS takes --rate_of_cooling in deg/millennium, so s [degC/decade] = rc/100:
# rc = 10 is the 0.1 degC/decade pathway used so far. rc = 0 is the no-SAI
# reference, which defines the G = 0 state of the SV*(G) curve.
SVG_COOLING_RC <- c(0, 5, 10, 15, 20, 25, 30)

# The ramp the unlabelled headline outputs of svg_analysis.R refer to. Without
# it svg_main() takes whichever gdx solved last, so adding a ramp to the grid
# would silently move svg_commitment.png onto a different pathway.
SVG_REFERENCE_RC <- 10
svg_rc_to_s <- function(rc) rc / 100          # deg/millennium -> degC/decade
svg_s_to_rc <- function(s)  round(s * 100)

# Ramp geometry passed to GAMS with every run, so the grid is reproducible from
# configuration alone. The ramp-down is placed beyond the analysis horizon: these
# are sustained-deployment pathways.
SVG_RAMP <- list(start_rampup = 2025, end_rampup = 2100,
                 start_rampdown = 3000, end_rampdown = 3100)

# Climate and scenario settings shared by every run of the grid.
SVG_RUN_DEFAULTS <- list(rcp = "RCP45", gas = "co2", ecs = 30, tcr = 18, pulse_time = 5)

# Test subset: run and validate these two before committing to the full grid.
SVG_COOLING_RC_MVP <- c(10, 20)


# -----------------------------------------------------------------------------
# 4. CONVERGENCE AND CLASSIFICATION TOLERANCES  (spec sections 2.12, 3.3)
# -----------------------------------------------------------------------------
# A commitment whose optimum sits at the last available horizon is only called
# "infinite" if the tail has converged: the value still to come over the final
# SVG_TAIL_FRAC of the horizon must be below SVG_TAIL_TOL of |SV*|.
SVG_TAIL_FRAC <- 0.10
SVG_TAIL_TOL  <- 0.01
SVG_SCC_TOL   <- 0.005   # SCC tail-convergence tolerance (0.5%)
SVG_SV_TOL    <- 0.01    # SV* tail-convergence tolerance (1%)

# A stored scenario is usable only if GAMS reported solvestat 1 and modelstat 1 or 2.
SVG_OK_SOLVESTAT <- 1
SVG_OK_MODELSTAT <- c(1, 2)


# -----------------------------------------------------------------------------
# 5. UNITS  (spec section 2.7, 5.5)
# -----------------------------------------------------------------------------
# Q_SRM is the model's SAI deployment, in TgS/yr: eq_mortsrm multiplies it by
# mort_srm = 7e3 deaths per unit, matching the ~7400 deaths per TgS/yr the
# project calibration uses. Nothing in the model is expressed in TgSO2, but the
# conversion is defined here so any external coefficient can be brought in safely.
SVG_TGSO2_PER_TGS <- 64 / 32    # molar mass ratio SO2:S
svg_tgso2_to_tgs <- function(x) x / SVG_TGSO2_PER_TGS
svg_tgs_to_tgso2 <- function(x) x * SVG_TGSO2_PER_TGS

SVG_GTCO2_TO_TCO2 <- 1e9        # W_EMI('co2') is GtCO2/yr
SVG_TN_USD        <- 1e12       # damages, y, VLL and COST_SRM are trillion USD


# -----------------------------------------------------------------------------
# 6. FIGURES
# -----------------------------------------------------------------------------
SVG_FIG_MAX_T   <- 200     # commitment-time axis cap on the zoomed figures [years]
SVG_FIG_DPI     <- 200
SVG_FIG_WIDTH   <- 9
SVG_FIG_HEIGHT  <- 5.5
