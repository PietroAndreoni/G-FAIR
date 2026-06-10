#### Scenario-matrix emission metrics for the coremoval experiment.
####
#### For every coremoval scenario in --folder (all capture_co2 x methane_source
#### combinations) this computes GWP, GTP and the relative social value at each
#### horizon, ALWAYS using the CO2 emission pulse as the reference (denominator)
#### for GWP and GTP, and as the reference for the social-value ratio. Every
#### quantity is the marginal difference from the common base scenario.
####
#### GWP(H)  = AGWP_scenario(H) / AGWP_co2pulse(H)     (total forcing integral ratio)
#### GTP(H)  = dT_scenario(pulse+H) / dT_co2pulse(pulse+H)
#### social_value_fraction(H) = NPV climate-damage difference of the scenario,
####            accumulated to pulse+H, divided by the same for the CO2 pulse.
####
#### All machinery is reused from gdx_difference_profiles.R::compute_emission_metrics().

suppressPackageStartupMessages({ library(tidyverse) })

this_file <- local({
  cmd <- commandArgs(FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit) > 0) sub("^--file=", "", hit[[1]]) else "Plots/matrix_metrics.R"
})
source(file.path(dirname(this_file), "gdx_difference_profiles.R"))

args        <- parse_args(commandArgs(trailingOnly = TRUE))
folder      <- args$folder    %||% "Results_matrix"
base_gdx    <- args$base      %||% "RCP45_EXPbase_ECS32.4_TCR17.9_IChistorical_run"
co2_gdx     <- args$co2       %||% "RCP45_EXPpulse_GASco2_ECS32.4_TCR17.9_PT5_IChistorical_run"
pulse_time  <- as.numeric(args$pulse_time %||% "5")
horizons    <- if (is.null(args$horizons)) c(20, 50, 100) else as.numeric(str_split(args$horizons, ",", simplify = TRUE))

# All coremoval scenarios present in the folder (capture_co2 x methane_source).
removal_files <- list.files(folder, pattern = "^RCP45_EXPremoval.*\\.gdx$", full.names = FALSE)
scenarios     <- str_remove(removal_files, "\\.gdx$")

metrics <- compute_emission_metrics(
  experiment_gdx = scenarios,
  co2_gdx        = co2_gdx,
  base_gdx       = base_gdx,
  folder         = folder,
  horizons       = horizons,
  pulse_time     = pulse_time
)

# Decode capture_co2 and methane_source from the scenario name, and keep both the
# headline ratios and the raw components that build them, so the CSV is fully
# self-documenting and editable.
full <- metrics %>%
  mutate(
    capture_co2    = str_match(scenario, "EXPremoval([a-z0-9]+)_PT")[, 2],
    methane_source = str_match(scenario, "SRC([a-z]+)$")[, 2]
  ) %>%
  transmute(
    capture_co2,
    methane_source,
    horizon_years          = horizon,
    end_year,
    # --- headline metrics (all referenced to the CO2 pulse) ---
    GWP,
    GTP,
    social_value_fraction          = social_cost_fraction,
    social_value_fraction_infinite = social_cost_fraction_infinite,
    # --- raw components (experiment vs base, and the CO2-pulse reference) ---
    agwp_experiment_Wm2yr = agwp_experiment,   # integral of forcing diff over [pulse, pulse+H]
    agwp_co2_Wm2yr        = agwp_co2,
    dT_experiment_K       = agtp_experiment,   # temperature diff at pulse+H
    dT_co2_K              = agtp_co2,
    npv_damage_experiment_USD = social_cost_experiment,  # NPV damage diff accumulated to pulse+H
    npv_damage_co2_USD        = social_cost_co2,
    scenario
  ) %>%
  arrange(capture_co2, methane_source, horizon_years)

options(pillar.sigfig = 4, width = 220)
cat("\n=== Scenario-matrix metrics (reference = CO2 pulse; differences vs base) ===\n\n")
print(as.data.frame(full %>% select(capture_co2, methane_source, horizon_years,
                                    GWP, GTP, social_value_fraction,
                                    social_value_fraction_infinite)), digits = 4)

out_csv <- file.path(folder, "matrix_metrics.csv")
readr::write_csv(full, out_csv)
cat("\nWrote", out_csv, "\n")
