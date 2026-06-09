#### Global Warming Potential of methane from the FaIR pulse experiments.
####
#### GWP_CH4(H) = AGWP_CH4(H) / AGWP_CO2(H), where each AGWP is the time integral
#### over [pulse_time, pulse_time + H] of the radiative-forcing difference that a
#### 1000-tonne pulse of the gas produces relative to the base scenario. Because
#### both the CH4 and the CO2 pulse are 1000 t of their respective gas (see
#### experiments/pulse.gms), the ratio of integrals is already the per-mass GWP.
####
#### Two forcing definitions are reported side by side:
####   GWP_direct  uses ONLY the gas's own band forcing -- FORCING('ch4') in the
####               numerator and FORCING('co2') in the denominator. It excludes the
####               indirect methane channels (CO2 from CH4 oxidation, which lands in
####               FORCING('co2'); tropospheric O3, FORCING('o3trop'); and
####               stratospheric water vapour, FORCING('h2o') = 0.12*FORCING('ch4')).
####   GWP_total   uses the model's full radiative forcing summed across all gases,
####               so every indirect channel is included. This reproduces the GWP
####               column of compute_emission_metrics() in gdx_difference_profiles.R.
####
#### The arithmetic, GDX extraction and forcing-difference helpers are reused from
#### gdx_difference_profiles.R (sourced below): extract_expression_delta() evaluates
#### "experiment minus base" on a FORCING expression (whole symbol or a
#### FORCING[ghg=...] slice), and resolve_gdx() turns a bare scenario name into a
#### path inside --folder.
####
#### CLI usage from the repository root (defaults target Results_coremoval):
####   Rscript Plots/methane_gwp_forcing.R
####   Rscript Plots/methane_gwp_forcing.R --folder=Results_coremoval \
####     --ch4=RCP45_EXPpulse_GASch4_ECS3.24_TCR1.79_PT5_IChistorical_run \
####     --co2=RCP45_EXPpulse_GASco2_ECS3.24_TCR1.79_PT5_IChistorical_run \
####     --base=RCP45_EXPbase_ECS3.24_TCR1.79_IChistorical_run \
####     --horizons=20,50,100 --pulse_time=5
####
#### Scenario names are given WITHOUT the .gdx extension; the extension and the
#### --folder location are added for you. An optional --output=<file.csv> writes the
#### result table.

suppressPackageStartupMessages({
  library(tidyverse)
})

# Reuse the GDX/forcing helpers from the sibling script. Resolve its path from this
# file's own location so the source works regardless of the caller's directory.
this_file <- local({
  cmd <- commandArgs(FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit) > 0) sub("^--file=", "", hit[[1]]) else "Plots/methane_gwp_forcing.R"
})
source(file.path(dirname(this_file), "gdx_difference_profiles.R"))

DEFAULT_FOLDER <- "Results"
DEFAULT_BASE <- "RCP45_EXPbase_ECS32_TCR18_IChistorical_run"
DEFAULT_CH4 <- "RCP45_EXPpulse_GASch4_ECS32_TCR18_PT5_IChistorical_run"
DEFAULT_CO2 <- "RCP45_EXPpulse_GASco2_ECS32_TCR18_PT5_IChistorical_run"
DEFAULT_HORIZONS <- c(20, 50, 100)

#### Compute methane GWP at each horizon under both forcing definitions.
#### Returns one row per horizon with the AGWPs (W/m2*yr per 1000-t pulse) and the
#### dimensionless GWP for the direct-only and the total-forcing definitions.
compute_methane_gwp <- function(ch4_gdx=paste0(DEFAULT_FOLDER,"/",DEFAULT_CH4),
                                co2_gdx=paste0(DEFAULT_FOLDER,"/",DEFAULT_CO2),
                                base_gdx=paste0(DEFAULT_FOLDER,"/",DEFAULT_BASE),
                                folder = ".",
                                horizons = DEFAULT_HORIZONS,
                                pulse_time = 5,
                                tstep = 1,
                                year_offset = 2019) {
  ch4_path <- resolve_gdx(ch4_gdx, folder = folder)
  co2_path <- resolve_gdx(co2_gdx, folder = folder)
  base_path <- resolve_gdx(base_gdx, folder = folder)

  # Integral of a delta profile over the closed window [t_start, t_end] using
  # tstep as the rectangle-rule weight; mirrors compute_emission_metrics().
  integrate_window <- function(df, t_start, t_end) {
    d <- df %>% filter(t >= t_start, t <= t_end)
    sum(d$value, na.rm = TRUE) * tstep
  }

  # Forcing differences (scenario minus base), arranged on t.
  forcing_direct_ch4 <- extract_expression_delta("FORCING[ghg=ch4]", ch4_path, base_path) %>% arrange(t)
  forcing_direct_co2 <- extract_expression_delta("FORCING[ghg=co2]", co2_path, base_path) %>% arrange(t)
  forcing_total_ch4  <- extract_expression_delta("FORCING", ch4_path, base_path) %>% arrange(t)
  forcing_total_co2  <- extract_expression_delta("FORCING", co2_path, base_path) %>% arrange(t)

  map_dfr(horizons, function(H) {
    t_end <- pulse_time + H

    agwp_direct_ch4 <- integrate_window(forcing_direct_ch4, pulse_time, t_end)
    agwp_direct_co2 <- integrate_window(forcing_direct_co2, pulse_time, t_end)
    agwp_total_ch4  <- integrate_window(forcing_total_ch4,  pulse_time, t_end)
    agwp_total_co2  <- integrate_window(forcing_total_co2,  pulse_time, t_end)

    tibble(
      horizon = H,
      pulse_time = pulse_time,
      end_t = t_end,
      end_year = year_offset + t_end,
      agwp_direct_ch4 = agwp_direct_ch4,
      agwp_direct_co2 = agwp_direct_co2,
      GWP_direct = agwp_direct_ch4 / agwp_direct_co2,
      agwp_total_ch4 = agwp_total_ch4,
      agwp_total_co2 = agwp_total_co2,
      GWP_total = agwp_total_ch4 / agwp_total_co2
    )
  })
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  folder <- args$folder %||% DEFAULT_FOLDER
  horizons <- if (is.null(args$horizons)) {
    DEFAULT_HORIZONS
  } else {
    as.numeric(str_split(args$horizons, ",", simplify = TRUE))
  }

  result <- compute_methane_gwp(
    ch4_gdx = args$ch4 %||% DEFAULT_CH4,
    co2_gdx = args$co2 %||% DEFAULT_CO2,
    base_gdx = args$base %||% DEFAULT_BASE,
    folder = folder,
    horizons = horizons,
    pulse_time = as.numeric(args$pulse_time %||% "5")
  )

  message("Methane GWP (CH4 pulse vs CO2 pulse, both relative to base)")
  message("Folder: ", folder)
  print(
    result %>%
      select(horizon, GWP_direct, GWP_total, agwp_direct_ch4, agwp_direct_co2,
             agwp_total_ch4, agwp_total_co2),
    width = Inf
  )

  if (!is.null(args$output)) {
    readr::write_csv(result, args$output)
    message("Wrote ", args$output)
  }
}
