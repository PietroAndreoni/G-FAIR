# =============================================================================
# all_parameters.R
#
# SINGLE, USER-FRIENDLY CONTROL FILE for the SRM-substitution Monte Carlo
# pipeline. Every hard-coded input, probability distribution, quantile, sd,
# discrete grid, economic constant and plotting choice that currently lives
# scattered across the generation / analysis / figure scripts is collected here
# with a recognizable name and a one-line explanation of what it does.
#
# Intended use: `source("all_parameters.R")` at the top of each script and refer
# to the names below instead of repeating literals. (The scripts are NOT yet
# wired to this file -- see the "WIRING" note at the bottom and the inconsistency
# report that accompanied this file. Until they are wired in, this file is the
# documented source of truth and the place to change a value once.)
#
# Conventions used throughout the pipeline:
#   * Time is in FAIR model-year units: t = 1 is calendar year 2020, so
#     calendar_year = 2019 + t  (t = 480  ->  2499; "2020..2500 window").
#   * ECS and TCR are stored as integer tenths of a kelvin (3.0 K -> 30).
#   * Distributions are sampled by inverse-CDF from a shared Sobol sequence in
#     Generate_montecarlo.R via the helpers in distribution_utils.R:
#       fit_distribution(), qlnorm_fit(), qnorm_fit(),
#       qlnorm_fit_trunc(), qnorm_fit_trunc().
#   * Each distribution list below maps directly onto a fit_distribution() call:
#       list(median=, q5=, q95=)        -> median + 5th/95th percentile
#       list(median=, sd=)              -> median + standard deviation
#     plus, where relevant, truncation bounds (lo/hi) and a `round` digit count.
#
# Base R only; sourcing this file has no side effects beyond defining objects.
# Lines tagged `# FLAG:` mark a value that disagrees with a comment/label/sibling
# script elsewhere -- see the inconsistency report. Values here are set to what
# the LIVE code currently does, so wiring this in is behaviour-preserving.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SAMPLER / RUN-CONTROL META  (Generate_montecarlo.R, montecarlo_utils.R)
# -----------------------------------------------------------------------------
SAMPLER_VERSION   <- "mc_sampler_v2"  # stamped on every realization; bump on any
                                      # change to the sampling logic below.
SAMPLING_METHOD   <- "sobol"          # "sobol" (Owen-scrambled QMC) or "montecarlo"
DEFAULT_SEED      <- 123L             # RNG seed (reproducibility / QMC net selection)

N_SCENARIOS_DEFAULT <- 8192L          # default number of draws per Generate call.
# NOTE: changed from 8096 -> 8192 (2^13). 8096 is not a power of two and tripped
#       the Sobol balance warning on every default run.

# Dimension of the unit hypercube each realization occupies, and the FIXED column
# -> parameter assignment of the Sobol/uniform draw matrix. Changing the order
# here means re-aligning Generate_montecarlo.R; documented so the mapping is
# explicit and never silently reshuffled.
N_SOBOL_DIM <- 19L
SOBOL_COLUMN_MAP <- c(
  ecs = 1, tcr = 2,            # 1-2 correlated joint lognormal (climate)
  rcp = 3, pulse = 4, cool = 5, term = 6, start = 7,
  term_delta = 8,             # 8 geometric termination delay (hazard from col 12)
  theta = 9, alpha = 10, delta = 11,
  prob = 12,                  # 12 per-year termination hazard (drives col 8)
  mortality_srm = 13, forctoTg = 14, TgtoUSD = 15,
  mortality_ozone = 16, vsl = 17, vsl_eta = 18, dg = 19
)


# -----------------------------------------------------------------------------
# 1. CLIMATE SENSITIVITY  (joint lognormal ECS-TCR, FAIR v1.3 style)
#    Generate_montecarlo.R lines ~94-104, 378-385
# -----------------------------------------------------------------------------
# Equilibrium climate sensitivity [K]: lognormal from median + 5/95 percentiles.
ECS_DIST <- list(median = 3.0, q5 = 2.0, q95 = 5.0)
# Transient climate response [K]: lognormal from median + 5/95 percentiles.
TCR_DIST <- list(median = 1.8, q5 = 1.2, q95 = 2.4)
# Target Pearson correlation between ECS and TCR on the ORIGINAL (non-log) scale.
ECS_TCR_CORR <- 0.81
# Hard constraint enforced after rounding to tenths: every stored draw has ECS > TCR.
ECS_TCR_ENFORCE_ECS_GT_TCR <- TRUE
ECS_TCR_ROUND_TENTHS <- TRUE          # store as integer tenths of K (3.0 -> 30)


# -----------------------------------------------------------------------------
# 2. SRM / SAI POLICY DISCRETE GRIDS  (equally-weighted, sampled jointly)
#    Generate_montecarlo.R lines ~109-113, 272-274
# -----------------------------------------------------------------------------
RCP_CHOICES   <- c("RCP3PD", "RCP45", "RCP6", "RCP85")  # baseline emission scenario
PULSE_CHOICES <- seq(5, 80, by = 1)        # pulse_time: deployment year offset (yr from 2020)
COOL_CHOICES  <- seq(0, 40, by = 1)        # rate_of_cooling (cool==0 => "no SRM", see below)
TERM_CHOICES  <- seq(2200, 2600, by = 100) # geo_end:  calendar year SAI ramp-down ends
START_CHOICES <- seq(2025, 2100, by = 25)  # geo_start: calendar year SAI ramp-up starts

# Sentinel values written when cool==0 (degenerate "no deployment") so term/start
# do not alias a real policy. NB Figure_2.R uses DIFFERENT sentinels (2500 / max+100)
# when it re-derives the no-SRM rows -- see inconsistency report.
COOL0_TERM_SENTINEL  <- 2700
COOL0_START_SENTINEL <- 2700


# -----------------------------------------------------------------------------
# 3. TERMINATION HAZARD & ANALYSIS HORIZON
#    montecarlo_utils.R lines 4-6; Generate_montecarlo.R lines ~115-130, 233-240
# -----------------------------------------------------------------------------
# Termination is a constant per-year hazard h (= the stored `prob`). The delay
# from deployment to termination is geometric on {1,2,...}; h is drawn
# LOG-UNIFORM on [HAZARD_LO, HAZARD_HI]. term_delta = pulse + delay, censored at
# T_HORIZON (a draw landing at the horizon == "no termination within the run").
HAZARD_LO <- 1e-4   # lower per-year hazard bound (median delay ~6930 yr at 1e-4)
HAZARD_HI <- 1e-1   # upper per-year hazard bound (median delay ~7 yr at 1e-1)
# RESOLVED: 1e-4 is the live sampled bound; the stale "1e-3" text in
#       Generate_montecarlo.R (comment + diagnostics title) was corrected to 1e-4.
#       Value left at 1e-4 (behaviour-preserving). Set to 1e-3 here if you instead
#       want the looser termination distribution implied by the old comments.

# FAIR projection horizon in model-year units (t = 1 is 2020; 480 -> 2020..2500).
# FAIR integrates t = 1..1000 but post-processing keeps only t <= T_HORIZON.
T_HORIZON <- 480L
# RESOLVED: the Analyze scripts now reference T_HORIZON instead of literal 480.
BASE_YEAR <- 2020L  # calendar year of t = 1 (calendar_year = BASE_YEAR - 1 + t)


# -----------------------------------------------------------------------------
# 4. POST-PROCESSING UNCERTAIN PARAMETERS
#    Sampled once in Generate_montecarlo.R (lines ~254-264), read by Analyze.
#    Each entry: the distribution, any truncation, the rounding applied to the
#    stored draw, and what the parameter physically controls.
# -----------------------------------------------------------------------------

# theta -- SAI injection / mixing angle [degrees]. Lognormal truncated to (0,90].
THETA_DIST  <- list(median = 10, q5 = 3, q95 = 30)
THETA_TRUNC <- c(lo = 0, hi = 90)
THETA_ROUND <- 0L            # rounded to whole degrees

# alpha -- climate-damage coefficient [fraction of GDP per K^2]. Lognormal from
# median + sd, with sd expressed as median * (382/179) (keeps the original ratio).
ALPHA_DIST  <- list(median = 0.00575, sd = 0.00575 * 382 / 179)
ALPHA_ROUND <- 4L

# delta -- pure rate of time preference / discount rate [1/yr]. Uniform on a
# 0.01-step grid {0.01,...,0.07}. The SAMPLING range is widened by +/- half a
# step (0.005, 0.075) so the edge bins 0.01 and 0.07 receive full weight.
DELTA_UNIF        <- c(min = 0.005, max = 0.075)  # widened sampling bounds
DELTA_GRID_STEP   <- 0.01                          # nominal grid resolution
DELTA_GRID_RANGE  <- c(0.01, 0.07)                 # nominal (pre-widening) range
DELTA_ROUND       <- 2L

# prob -- per-year termination hazard. Same object as Section 3 (log-uniform on
# [HAZARD_LO, HAZARD_HI]); stored raw (no rounding beyond its 4 sig figs).

# mortality_srm -- excess mortality attributable to SAI sulfate [deaths per Tg-S/yr].
# Lognormal from median + 5/95 percentiles, rounded to integer.
MORTALITY_SRM_DIST  <- list(median = 7400, q5 = 2300, q95 = 16000)
MORTALITY_SRM_ROUND <- 0L
# FLAG: these quantiles are mutually inconsistent for a lognormal
#       (q5*q95 = 3.68e7 != median^2 = 5.48e7) -> fit_distribution() returns a
#       least-discrepancy COMPROMISE and warns. Intended, but verify the targets.

# forctoTg -- radiative forcing -> SAI sulfur loading conversion. Drawn as the
# RECIPROCAL of a uniform: forctoTg = 1 / U(0.2, 1.5)  [Tg-S per (W/m^2)], 2 dp.
FORCTOTG_INV_UNIF <- c(min = 0.2, max = 1.5)
FORCTOTG_ROUND    <- 2L

# TgtoUSD -- direct deployment cost per unit sulfur [USD per Tg-S/yr]. Uniform.
TGTOUSD_UNIF  <- c(min = 0.75, max = 3)
TGTOUSD_ROUND <- 2L

# mortality_ozone -- ozone-related mortality per unit CH4 concentration change.
# Normal truncated at >= 0, then SCALED by 1/100 and rounded to integer.
MORTALITY_OZONE_DIST  <- list(median = 11250, q5 = 5000, q95 = 17500)
MORTALITY_OZONE_TRUNC <- c(lo = 0, hi = Inf)
MORTALITY_OZONE_SCALE <- 1 / 100
MORTALITY_OZONE_ROUND <- 0L

# vsl -- value of a statistical life [USD]. Uniform in millions, rounded to whole
# millions then multiplied by 1e6  => {8,...,14} x 1e6 for most draws.
VSL_UNIF_MILLIONS <- c(min = 7.5, max = 13.6)
VSL_SCALE         <- 1e6
VSL_ROUND         <- 0L   # rounding applied (in millions) before scaling

# vsl_eta -- income elasticity of VSL [-]. Uniform on a 0.1-step grid {0.4,...,1.0};
# sampling range widened by +/- half a step (0.35, 1.05) to equalize edge bins.
VSL_ETA_UNIF       <- c(min = 0.35, max = 1.05)  # widened sampling bounds
VSL_ETA_GRID_STEP  <- 0.1
VSL_ETA_GRID_RANGE <- c(0.4, 1.0)                # nominal (pre-widening) range
VSL_ETA_ROUND      <- 1L

# dg -- annual GDP growth rate used to scale damages forward [1/yr]. Normal.
DG_DIST  <- list(median = 0.015, sd = 0.005)
DG_ROUND <- 3L


# -----------------------------------------------------------------------------
# 5. DETERMINISTIC MAIN-SCENARIO OVERRIDES  (--base=T)
#    Generate_montecarlo.R lines ~276-284
# -----------------------------------------------------------------------------
# When the main (non-Monte-Carlo) scenario is requested, the policy axes are
# pinned and a fixed pair of pulse times is run.
MAIN_RCP        <- "RCP45"
MAIN_COOL       <- 10
MAIN_TERM       <- 2400
MAIN_START      <- 2025
MAIN_PULSE_SET  <- c(5, 30)   # one row each
# Post-processing parameters fixed for the deterministic scenario:
MAIN_VSL        <- 10 * 1e6
MAIN_DELTA      <- 0.02
MAIN_VSL_ETA    <- 1
MAIN_THETA      <- 10         # default --angle; "na" keeps theta uncertain (0-90 valid)


# -----------------------------------------------------------------------------
# 6. ECONOMIC / DAMAGE-FUNCTION CONSTANTS  (Analyze_montecarlo.R)
# -----------------------------------------------------------------------------
# Initial gross world product [USD]: 2019 World Bank GDP (constant 2015 USD)
# rescaled to 2020 USD by the 1.085 deflator.
GWP_INITIAL <- 85.28 * 1e12 * 1.085

# Ratio of global to US per-capita GDP in 2020 (1.82 x 10937/63515), adjusted to
# match McDuffie 2023; converts US-denominated VSL/ozone terms to a global basis.
GLOBAL_TO_US_PC <- 1.82 * 10937 / 63515

# Cap on cumulative GDP growth multiple in compute_gwpt(): gwpt <= GWP_INITIAL * GWP_MAX.
GWP_MAX <- 15
# RESOLVED: the stale (1+dg)^(280-1) comment above compute_gwpt() was corrected
#       to describe the actual cap, GWP_INITIAL * GWP_MAX (=15x). Value unchanged.

# GDP growth-rate decay schedule in compute_gwpt(): dg applies until DG_DECAY_T1,
# halves to dg/2 until DG_DECAY_T2, then quarters to dg/4 thereafter
# (model-year breakpoints; 80 -> 2100, 180 -> 2200).
DG_DECAY_T1 <- 80L
DG_DECAY_T2 <- 180L


# -----------------------------------------------------------------------------
# 7. FAIR RUN CONFIGURATION  (Run_montecarlo.R)
# -----------------------------------------------------------------------------
GASES_FULL <- c("ch4", "co2")  # gases simulated in the full Monte Carlo
GASES_MAIN <- c("ch4")         # gases simulated for --base=T

# SAI ramp window offsets [years] used when composing the GAMS call:
#   start_rampdown = geo_end   - RAMP_DOWN_YEARS
#   end_rampdown   = geo_end
#   start_rampup   = geo_start
#   end_rampup     = geo_start + RAMP_UP_YEARS
RAMP_DOWN_YEARS <- 100
RAMP_UP_YEARS   <- 100

# HPC submission (bsub) defaults -- environment-specific, edit per cluster.
HPC_QUEUE          <- "p_short"
HPC_PROJECT        <- "0638"
HPC_MEMORY         <- "2G"
HPC_GAMS_PATH      <- "/work/cmcc/pa12520/gams40.4_linux_x64_64_sfx"  # igdx() on HPC
SOLVE_OK_STATUS    <- 112   # GAMS exit status treated as "solved" (not problematic)


# -----------------------------------------------------------------------------
# 8. PULSE-EFFECTS ANALYSIS  (Analyze_montecarlo_pulse_effects.R)
# -----------------------------------------------------------------------------
TRAJ_PROBS    <- c(0.05, 0.25, 0.5, 0.75, 0.95)  # percentiles for coherent-run selection
RMSE_WINDOW   <- 200    # years from pulse used for RMSE / monotonic checks
PULSE_CHUNK_N <- 400L   # matched scenario-groups per batch
QUANTILE_TYPE <- 8L     # stats::quantile type used everywhere in the pulse analysis


# -----------------------------------------------------------------------------
# 9. FIGURE / PLOTTING SETTINGS  (Figure_1.R, Figure_2.R, Figure_3*.R)
# -----------------------------------------------------------------------------
# Result folders the figures read from (currently hard-coded per figure).
RESULTS_FOLDER_MAIN     <- "Results_1903"               # Figure_1, Figure_2
RESULTS_FOLDER_FIG3     <- "Results_base_1903_angle30"  # Figure_3
RESULTS_FOLDER_FIG3_PAT <- "^Results_base_1903"         # Figure_3_SI folder glob

# Ribbon/line percentiles and time window shared by the trajectory figures.
FIG_PERCENTILES  <- c(0.05, 0.25, 0.5, 0.75, 0.95)
FIG_INNER_RIBBON <- c(0.25, 0.75)
FIG_OUTER_RIBBON <- c(0.05, 0.95)
FIG_TREL_WINDOW  <- 200   # x-axis cap [years from pulse]

# Outlier trimming applied before density/sensitivity plots (Figure_2/3): keep
# rows whose listed columns sit within the [lo, hi] quantiles.
FIG_OUTLIER_COLS <- c("ecs", "tcr", "alpha", "theta", "mortality_srm",
                      "mortality_ozone", "dg", "TgtoUSD", "forctoTg")
FIG_OUTLIER_QLO  <- 0.001
FIG_OUTLIER_QHI  <- 0.999

# KDE bandwidth for the normalized-cost density (Figure_2): adjust = log(C * sd * n^p).
FIG_KDE_BW_CONST <- 1.06
FIG_KDE_BW_POW   <- -1/5

# gsaot global sensitivity analysis settings (Figure_2).
GSA_M     <- 15   # number of partitions
GSA_BOOT  <- TRUE
GSA_R     <- 100  # bootstrap replicates

# Figure_3 MACC density scaling and display window.
DENSITY_COST_UNIT <- 100               # density shown as % per $100/ton
FIG3_XLIM         <- c(0, 10000)       # abatement-cost axis [USD/ton CH4]
FIG3_YLIM         <- c(0, 50)          # Figure_3 only
FIG3_MACC_YEARS   <- c(2025, 2050)     # pulse years shown as MACC curves

# CH4 cost-axis unit conversions on Figure_3 ($/tonCH4 -> $/tonCO2eq -> $/tonCeq):
CH4_GWP100        <- 25     # AR4 GWP-100 of CH4 (also /25 tick on the x-axis)
C_PER_CO2         <- 12/44  # carbon mass fraction of CO2 (tonCeq = tonCO2 * 12/44)
USD_DEFLATOR_2010_2020 <- 1.18  # 2010->2020 USD deflator on Harmsen MACC costs
                                # (Harmsen data is in 2010 USD; both Fig_3 & Fig_3_SI)

# AR4 GWP-100 values for non-CO2 gases (Figure_3_SI.R; Figure_3.R relies on this
# vector being present in the environment -- see FLAG below).
AR4_GWP100 <- c(c2f6 = 12200, c6f14 = 9300, cf4 = 7390,
                hfc125 = 124, hfc134a = 1430, hfc143a = 4470, hfc152a = 124,
                hfc227ea = 3220, hfc23 = 14800, hfc236fa = 675, hfc245ca = 693,
                hfc32 = 675, hfc4310 = 1640, ch4 = 25, n2o = 298)
# RESOLVED: Figure_3.R previously used an undefined `ar4gwp` (defined only in
#       Figure_3_SI.R / the Harmsen input script) and errored standalone. Both
#       figures now source this file and use AR4_GWP100.

# Harmsen non-CO2 MACC input files. Figure_3 and Figure_3_SI now both read the
# same .parquet inputs (Figure_3_SI was switched off the .lz4 variant by hand).
HARMSEN_BASELINE_FIG3    <- "input/data/harmsen_nonco2_baseline.parquet"
HARMSEN_MACC_FIG3        <- "input/data/harmsen_nonco2_macc.parquet"
HARMSEN_BASELINE_FIG3_SI <- "input/data/harmsen_nonco2_baseline.parquet"
HARMSEN_MACC_FIG3_SI     <- "input/data/harmsen_nonco2_macc.parquet"
MACC_CO2_FILE            <- "input/data/macc_ed_full_2022.lz4.parquet"
# NOTE: MACC_CO2_FILE is referenced by both figures' macc_co2 block but is not
#       present in the repo tree -- confirm it exists before running.

# Figure_2 sensitivity-panel axis transforms / limits (documented, not all wired):
FIG2_TERM_YEAR_BASE <- 2020      # x = 2020 + term  (termination calendar year)
FIG2_TERM_XLIM      <- c(2020, 2400)
FIG2_DELTA_PCT      <- 100       # delta plotted as delta*100 [%]
FIG2_ALPHA_PCT      <- 100       # alpha plotted as alpha*100 [%GDP/K^2]
FIG2_ECS_TENTHS     <- 10        # ecs plotted as ecs/10 [K]
FIG2_PANEL_QLO      <- 0.05      # per-panel x-axis trimming quantiles
FIG2_PANEL_QHI      <- 0.95
FIG2_TAIL_QUANTILE  <- 0.95      # "bad tail" threshold for the contrast matrix


# -----------------------------------------------------------------------------
# WIRING NOTE
# -----------------------------------------------------------------------------
# To make this the ONLY place these values are decided, each script should:
#   1. source("all_parameters.R") near the top (after the self-locating logic
#      already used for distribution_utils.R / montecarlo_utils.R), and
#   2. replace the literal with the name above, e.g. in Generate_montecarlo.R:
#         ecs_par <- do.call(function(...) fit_distribution(..., return_params=TRUE), ECS_DIST)
#         pulse_choices <- PULSE_CHOICES
#         theta <- round(qlnorm_fit_trunc(unit_draws[,SOBOL_COLUMN_MAP["theta"]],
#                        lo=THETA_TRUNC["lo"], hi=THETA_TRUNC["hi"],
#                        median=THETA_DIST$median, q5=THETA_DIST$q5, q95=THETA_DIST$q95),
#                        THETA_ROUND)
#   montecarlo_utils.R currently DEFINES MC_T_HORIZON / MC_HAZARD_* and is sourced
#   by everything; if this file is wired in there, keep the MC_* aliases pointing
#   at T_HORIZON / HAZARD_LO / HAZARD_HI so the validators and tests stay intact.
# =============================================================================
