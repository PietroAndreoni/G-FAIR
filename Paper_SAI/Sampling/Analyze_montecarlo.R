install_witchtools <- function(){
  if (!"remotes" %in% rownames(installed.packages()))
    install.packages("remotes", repos="http://cloud.r-project.org")
  remotes::install_github("witch-team/witchtools")
}
if (!"witchtools" %in% rownames(installed.packages())) {
  install_witchtools()
  if (!requireNamespace("witchtools")) stop("Package witchtools not found")
}
if (packageVersion("witchtools") < "0.4.0") {
  cat("Need a more recent version of witchtools. Try to update.\n")
  install_witchtools()
  cat(paste("Installed version:",packageVersion("witchtools"),"\n"))
  if (packageVersion("witchtools") < "0.4.0") {
    stop("Please install witchtools version >= 0.4.0.\n")
  }
}

library(witchtools)

# Packages
pkgs <- c('data.table','stringr','docopt')
res <- lapply(pkgs, require_package)
require_gdxtools()

start.time <- Sys.time()

# Locate the Paper_SAI root (holds all_parameters.R), load the control file and
# the shared validators / MC_* constants in Utilities/. Works from any wd.
# Locate the Paper_SAI folder (holds all_parameters.R) robustly so the script
# works under Rscript (--file), RStudio "Source" (sys.frame $ofile), and an
# interactive console whose working dir is at/under/above the project.
.find_paper_root <- function() {
  starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
              unlist(lapply(sys.frames(), function(f) f$ofile)), getwd())
  for (s in starts[nzchar(starts)]) {
    d <- if (dir.exists(s)) s else dirname(s)
    repeat {
      if (file.exists(file.path(d, "all_parameters.R")))
        return(normalizePath(d, "/", FALSE))
      if (file.exists(file.path(d, "Paper_SAI", "all_parameters.R")))
        return(normalizePath(file.path(d, "Paper_SAI"), "/", FALSE))
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Cannot locate all_parameters.R (Paper_SAI control file); set the working ",
       "directory to the project root or the Paper_SAI folder.", call. = FALSE)
}
source(file.path(.find_paper_root(), "all_parameters.R"))
source(file.path(PAPER_ROOT, "Utilities", "montecarlo_utils.R"))

# extract names function rewritten using data.table/vector ops
extract_names_dt <- function(files_vec) {
  DT <- data.table(gdx = files_vec)
  DT[, file := str_remove_all(gdx, ".*(?=RCP)|.gdx")]
  DT[, gas := str_extract(file, "(?<=GAS).+?(?=_)")]
  DT[, rcp := str_extract(file, "(?<=RCP).+?(?=_)")]
  DT[, cool_rate := str_extract(str_remove(file, "RCP"), "(?<=RC).+?(?=_)")]
  DT[, pulse_time := str_extract(gdx, "(?<=PT).+?(?=_)")]
  DT[, geo_end := str_extract(str_remove(file, "ECS"), "(?<=EC).+?(?=_)")]
  DT[, geo_start := str_extract(file, "(?<=BC).+?(?=_)")]
  DT[, ecs := str_extract(file, "(?<=ECS).+?(?=_)")]
  DT[, tcr := str_extract(file, "(?<=TCR).+?(?=_)")]
  DT[, term := str_extract(file, "(?<=TER).+?(?=_)")]
  DT[, experiment := str_extract(file, "(?<=EXP).+?(?=_)")]
  return(DT)
}


sanitize_dt <- function(DT) {
  setDT(DT)
  # ensure t numeric and <= T_HORIZON, drop gdx
  if("t" %in% names(DT)) DT[, t := as.numeric(t)]
  if("gdx" %in% names(DT)) DT <- DT[, !"gdx", with = FALSE]
  DT <- DT[ t <= T_HORIZON ]
  DT <- unique(DT)
  return(DT)
}


# fast replace numeric NA with 0 in-place
replace_num_na0 <- function(DT) {
  numcols <- names(DT)[vapply(DT, is.numeric, logical(1))]
  for (cname in numcols) set(DT, which(is.na(DT[[cname]])), cname, 0)
  invisible(NULL)
}

# vectorized gwpt computation
# Cumulative GDP path, capped at GWP_INITIAL * GWP_MAX (default 15x). Growth uses
# dg until model-year DG_DECAY_T1, dg/2 until DG_DECAY_T2, then dg/4 thereafter.
compute_gwpt <- function(DT, dg_col = "dg", t_col = "t", gwp_max = GWP_MAX) {
  t1 <- DG_DECAY_T1; t2 <- DG_DECAY_T2
  DT[t <= t1, gwpt := pmin(gwp * gwp_max,
                    gwp * (1 + dg)^(t-1))]
  DT[t > t1 & t <= t2, gwpt := pmin(gwp * gwp_max,
                  gwp * (1 + dg)^(t1-1) * (1 + dg/2)^(t-t1))]
  DT[t > t2, gwpt := pmin(gwp * gwp_max,
                           gwp * (1 + dg)^(t1-1) * (1 + dg/2)^(t2-t1) * (1 + dg/4)^(t-t2))]
}

npv_aggregator <- function(DT, keep_names = all_names) {
  DT = copy(DT)
  compute_gwpt(DT)
  DT[, dam := gwpt * (alpha * (pmax(0,temp)^2 - pmax(0,temp_srm)^2))]
  DT[, direct_cost := srm_masking * forctoTg * TgtoUSD]
  DT[, srm_pollution := vsl * globaltouspc * (gwpt/gwp)  ^ (vsl_eta) * srm_masking * forctoTg * mortality_srm]
  DT[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp)  ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
  DT[, imperfect_masking := 2 * gwpt * alpha * pmax(0,temp_base) * (1 - cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1 + srm/forc)]
  mc_assert_no_missing(DT, c("dam", "direct_cost", "srm_pollution",
                             "tropoz_pollution", "imperfect_masking"),
                       "npv_aggregator computed costs")
  # aggregator
  res <- DT[, .(
    dirnpv = sum(direct_cost / (1 + delta)^(t - as.numeric(pulse_time))),
    srmpnpv = sum(srm_pollution / (1 + delta)^(t - as.numeric(pulse_time))),
    ozpnpv = sum(tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time))),
    masknpv = sum(imperfect_masking / (1 + delta)^(t - as.numeric(pulse_time))),
    damnpv = sum(dam / (1 + delta)^(t - as.numeric(pulse_time)))
  ), by = keep_names]
  # convert to long format
  res_long <- melt(
    res,
    id.vars = keep_names,
    variable.name = "source",
    value.name = "cost")
  return(res_long)
}

# Total (NOT marginal) cost of the srm-only world (the "srmnopulse" run: SRM
# deployed at its baseline schedule, no GHG pulse), summed over every damage/cost
# component and evaluated across the FULL FAIR horizon. Returns, per id, the FIRST
# model-year in which that total reaches or exceeds gross world product (gwpt) --
# the point beyond which the world is no longer economically meaningful -- or NA
# for a realization that never breaches.
#
# Realizations are NOT rejected on this basis. The caller keeps the draw and
# truncates its time series at first_breach_t, so the NPVs are taken over the
# economically meaningful window only (see the truncation block below). The breach
# is permanent in practice -- once total damage reaches GDP it stays there for the
# rest of the horizon -- so a single cut point loses nothing that a
# breach-by-breach mask would keep.
#
# Every field here is read from the srm-only world, as intended:
#   temp_srm, concch4_base  <- experiment "srm"  (baseline SRM active, no pulse)
#   temp_base               <- experiment "base" (no SRM; counterfactual for masking)
#   srm                     <- forcing_srm, the baseline SAI deployment schedule.
# forcing_srm is a GAMS *parameter* set once per scenario and never re-assigned
# across experiments (termination acts on the SRM *variable*, which is identically
# 0 in the srm run), so taking it from the selected experiment's gdx yields the
# srm-run values. That makes the total gas- and termination-independent, so we key
# on the srm climate scenario only and pool over gases (the earliest breach across
# gases governs the draw).
srmnopulse_first_breach_t <- function(DT, id_col, climate_keys) {
  D <- copy(DT)
  compute_gwpt(D)
  D[, srm_dam           := gwpt * alpha * pmax(0, temp_srm)^2]
  D[, srm_dir_cost      := srm * forctoTg * TgtoUSD]
  D[, srm_pollution_tot := vsl * globaltouspc * (gwpt / gwp)^(vsl_eta) * srm * forctoTg * mortality_srm]
  D[, srm_tropoz        := vsl * globaltouspc * (gwpt / gwp)^(vsl_eta) * mortality_ozone * concch4_base]
  D[, srm_masking_tot   := gwpt * alpha * pmax(0, temp_base)^2 * (1 - cos(theta * pi / 180))]
  D[, srm_total_damage  := srm_dam + srm_dir_cost + srm_pollution_tot + srm_tropoz + srm_masking_tot]
  mc_assert_no_missing(D, c("gwpt", "srm_total_damage"), "srmnopulse_first_breach_t")
  D[, .(first_breach_t = {
          breached <- as.numeric(t)[srm_total_damage >= gwpt]
          if (length(breached) == 0L) NA_real_ else min(breached)
        }), by = c(id_col, climate_keys)]
}

# First year in which the NO-SRM world's own total damage reaches GDP. That world
# carries no SAI at all, so its total is climate damage plus ozone mortality from the
# ambient CH4 concentration -- none of the deployment-cost, sulfate-pollution or
# imperfect-masking terms that enter the srm world's total. Mirrors
# srmnopulse_first_breach_t in every other respect (full horizon, NA = never breaches).
# Both temp_base and concch4_base are read from experiment "base", so this is the
# no-SRM world judged by its own standard rather than the srm world's.
nosrm_first_breach_t <- function(DT, id_col) {
  D <- copy(DT)
  compute_gwpt(D)
  D[, nosrm_dam    := gwpt * alpha * pmax(0, temp_base)^2]
  D[, nosrm_tropoz := vsl * globaltouspc * (gwpt / gwp)^(vsl_eta) * mortality_ozone * concch4_base]
  D[, nosrm_total  := nosrm_dam + nosrm_tropoz]
  mc_assert_no_missing(D, c("gwpt", "nosrm_total"), "nosrm_first_breach_t")
  D[, .(first_breach_t = {
          breached <- as.numeric(t)[nosrm_total >= gwpt]
          if (length(breached) == 0L) NA_real_ else min(breached)
        }), by = id_col]
}

# Drop every timestep from a draw's first breach onward. `breach` carries one
# first_breach_t per ID (NA = never breaches, keep the full horizon). The breach
# year itself goes: total damage already meets GDP there, so it is the first
# non-economic step rather than the last economic one.
# ID is character in some tables (expected_chunk_keys) and integer in others, so
# match on as.character() rather than a join. An ID absent from `breach` gets NA and
# is kept whole, which is also the right answer for "never breaches".
truncate_at_breach <- function(DT, breach) {
  cut_t <- stats::setNames(as.numeric(breach$first_breach_t), as.character(breach$ID))
  ct <- cut_t[as.character(DT$ID)]
  DT[is.na(ct) | as.numeric(DT$t) < ct]
}

# Build the shared dropped-scenarios file consumed by Analyze_montecarlo_pulse_effects.R
# from the comprehensive per-draw discard list. A srm-climate scenario is fully
# dropped only when EVERY draw behind it (across term and both gases) was discarded --
# whether by the srmnopulse damage filter or by FAIR-run infeasibility. Totals come
# from id_montecarlo, so the decision is correct across chunk boundaries and --skip
# resumes. Also collapses any duplicate discard rows accumulated across resumes.
write_dropped_scenarios <- function() {
  if (!file.exists(discarded_runs_path)) return(invisible(NULL))
  disc <- unique(fread(discarded_runs_path))
  fwrite(disc, file = discarded_runs_path)                 # de-duplicate in place
  for (k in climate_keys) set(disc, j = k, value = as.character(disc[[k]]))
  totals <- id_montecarlo[, .(n_total = uniqueN(ID)), by = climate_keys]
  for (k in climate_keys) set(totals, j = k, value = as.character(totals[[k]]))
  # A climate scenario counts as fully dropped only when its draws were genuinely
  # discarded (damage filter / fair_infeasible). not_run draws were never submitted,
  # and the damage_truncated* reasons mark draws that are KEPT (just shortened), so
  # none of these may mark a scenario as dropped for the pulse-effects plots.
  disc_real <- if ("reason" %in% names(disc)) {
    disc[!reason %in% c("not_run", "damage_truncated", "damage_truncated_nosrm")]
  } else disc
  discounts <- disc_real[, .(n_disc = uniqueN(ID)), by = climate_keys]
  merged <- merge(discounts, totals, by = climate_keys, all.x = TRUE)
  dropped_scen <- merged[!is.na(n_total) & n_disc >= n_total, ..climate_keys]
  dropped_scen_path <- file.path(output_folder, dropped_scen_output)
  fwrite(dropped_scen, file = dropped_scen_path)
  cat("Discard filter: wrote", nrow(dropped_scen), "fully-dropped scenario(s) to",
      dropped_scen_path, "and", nrow(disc), "discarded run row(s) to",
      discarded_runs_path, "\n")
}

make_scc_diagnostics <- function(DT, key_cols, value_cols, context) {
  mc_assert_no_missing(DT, c(key_cols, "t", value_cols),
                       paste0(context, " diagnostics inputs"))
  diag <- DT[, list(
    observed_years = .N,
    min_t = min(as.numeric(t)),
    max_t = max(as.numeric(t))
  ), by = key_cols]

  for (cname in value_cols) {
    stats <- DT[, list(
      stat_min = min(get(cname), na.rm = TRUE),
      stat_max = max(get(cname), na.rm = TRUE),
      stat_mean = mean(get(cname), na.rm = TRUE),
      stat_negative_years = sum(get(cname) < 0, na.rm = TRUE),
      stat_zero_years = sum(get(cname) == 0, na.rm = TRUE),
      stat_nonfinite_years = sum(!is.finite(get(cname)))
    ), by = key_cols]
    setnames(
      stats,
      c("stat_min", "stat_max", "stat_mean",
        "stat_negative_years", "stat_zero_years", "stat_nonfinite_years"),
      paste0(cname, c("_min", "_max", "_mean",
                      "_negative_years", "_zero_years", "_nonfinite_years"))
    )
    diag <- merge(diag, stats, by = key_cols, all = TRUE)
  }
  diag
}

# NOTE: the post-processing uncertainties (theta, alpha, delta, prob,
# mortality_srm, forctoTg, TgtoUSD, mortality_ozone, vsl, vsl_eta, dg) are no
# longer sampled here. They are drawn once (via Sobol) in Generate_montecarlo.R
# and read from id_montecarlo.csv below. The fit_distribution() helper therefore
# lives only in Generate_montecarlo.R now.

prepare_join_table <- function(filter_experiment) {
  
  temp_exp <- copy(TATM)[experiment == filter_experiment]  # make copy to avoid modifying original
  temp_exp[, c("file", "experiment") := NULL]               # remove columns
  setnames(temp_exp, "value", "temp", skip_absent = TRUE)
  
  backsrm_exp <- copy(background_srm)[experiment == filter_experiment]  # make copy to avoid modifying original
  backsrm_exp[, c("file", "experiment") := NULL]               # remove columns
  setnames(backsrm_exp, "value", "srm", skip_absent = TRUE)
  
  dsrm_exp <- copy(SRM)[experiment == filter_experiment]  # make copy to avoid modifying original
  dsrm_exp[, c("file", "experiment") := NULL]               # remove columns
  setnames(dsrm_exp, "value", "srm_masking", skip_absent = TRUE)
  
  tot_forcing_exp <- copy(tot_forcing)[experiment == filter_experiment]  # make copy to avoid modifying original
  tot_forcing_exp[, c("file", "experiment") := NULL]               # remove columns
  setnames(tot_forcing_exp, "value", "forc", skip_absent = TRUE)
  
  concch4_exp <- copy(CONC)[ghg == "ch4" & experiment == filter_experiment]
  concch4_exp[, c("file", "experiment","ghg") := NULL]
  setnames(concch4_exp, "value", "concch4_pulse", skip_absent = TRUE)
  
  concch4_srm <- copy(CONC)[ghg == "ch4" & experiment == "srm"]
  concch4_srm[, c("file", "experiment","term","ghg") := NULL]
  setnames(concch4_srm, "value", "concch4_base", skip_absent = TRUE)
  
  cols <- c("t",base_scenarios,"value")
  temp_base <- copy(TATM)[experiment == "base", ..cols]
  setnames(temp_base, "value", "temp_base", skip_absent = TRUE)
  
  temp_srm <- copy(TATM)[experiment == "srm"]
  temp_srm[, c("file", "experiment","term") := NULL]
  setnames(temp_srm, "value", "temp_srm", skip_absent = TRUE)
  
  temp_srmpulse <- copy(TATM)[experiment == "srmpulse"]
  temp_srmpulse[, c("file", "experiment","term") := NULL]
  setnames(temp_srmpulse, "value", "temp_srmpulse", skip_absent = TRUE)
  
  cols <- c("t",pulse_scenarios,"value")
  temp_ghgpulse <- copy(TATM)[experiment == "pulse", ..cols]
  setnames(temp_ghgpulse, "value", "temp_ghgpulse", skip_absent = TRUE)
  
  base <- copy(full_grid)
  base <- base[tot_forcing_exp, on = c("t", scenario_names)]
  base <- base[temp_exp, on = c("t", scenario_names)]
  base <- backsrm_exp[base, on = c("t",scenario_names)]
  # forcing_srm is a GAMS parameter, so zero records are omitted from the gdx:
  # the SRM background forcing is genuinely 0 outside the deployment window.
  base[is.na(srm), srm := 0]
  base <- dsrm_exp[base, on = c("t",scenario_names)]
  base <- concch4_exp[base, on = c("t", scenario_names)]
  base <- concch4_srm[base, on = setdiff(c("t", scenario_names),"term") ]
  base <- temp_srm[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_srmpulse[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_base[base, on = c("t", base_scenarios)]
  base <- temp_ghgpulse[base, on = c("t",pulse_scenarios), nomatch = NULL]
  
  # The terminated run (srmpulsemaskedterm) carries a real, term-specific
  # temperature path, so its `term` (= termination time) MUST be kept: each path
  # stays matched to its own termination time. The other experiments are
  # termination-independent (their gdx has term = NA), so we drop the placeholder
  # term and broadcast that single path across every termination time present in
  # the chunk. Broadcasting the terminated run instead would mis-pair the
  # overshoot temperature with the wrong term (and double-count) whenever two
  # realizations share a FAIR scenario but differ in termination time.
  if (filter_experiment != "srmpulsemaskedterm") {
    base[, term := NULL] # term-independent run: drop placeholder, re-add below
    base <- merge(base, all_scenarios[!is.na(term)], by = setdiff(scenario_names, c("term")), allow.cartesian = TRUE)
  }
  # allow.cartesian: several realizations may share the same FAIR scenario (and
  # thus the same gdx) while differing in the post-processing parameters.
  base <- merge(base, id_montecarlo, by = setdiff(scenario_names,c("gas")), allow.cartesian = TRUE)
  mc_assert_no_missing(
    base,
    c("ID", "gas", "t", "temp", "srm", "srm_masking", "forc",
      "concch4_pulse", "concch4_base", "temp_srm", "temp_srmpulse",
      "temp_base", "temp_ghgpulse", post_cols),
    paste0("prepare_join_table(", filter_experiment, ")")
  )
  mc_assert_unique_key(base, c("ID", "gas", "t"), paste0("prepare_join_table(", filter_experiment, ")"))
  return(base)
}

# Usage
'Launch script to analyze montecarlo scenarios (produces a csv file in the same folder)

Usage:
  Analyze_montecarlo.R [-i <input>] [-o <results>] [--hpc <run_hpc>] [-p <plot_results>] [--chunk <chunk>] [--skip <skip>] [--res <output_folder>] [--termination <termination>]

Options:
-i <input>             Path where the montecarlo id are
-o <output>           Where to save output
--res <results_folder>  name of the results folder (default: Results_montecarlo). For multiple folders separate with -
--hpc <run_hpc>        T/F if running from Juno (T) or local (F)
--chunk <chunk>        how many scenarios to run together
--skip <skip>          skip or rerun existing scenarios
--termination <termination>  T: NPC of the terminated run (SRM off at term_delta); F: NPC of the run where SRM is never terminated (default T)
' -> doc

opts <- docopt(doc, version = 'Montecarlo')

# Folder names are user-specifiable; their locations are fixed. The GAMS .gdx are
# read from Results/ (under Paper_SAI); id_montecarlo.csv is read from, and the
# npc/scc CSV output is written to, the working folder under Sampling/.
res <- ifelse(is.null(opts[["res"]]), RUN_RESULTS_DEFAULT, as.character(opts["res"]) )
res <- file.path(RUN_RESULTS_PARENT, str_split(res,"-")[[1]])

input_folder <- file.path(MC_WORK_PARENT, ifelse(is.null(opts[["i"]]), MC_WORK_DEFAULT, as.character(opts["i"])) )
output_folder <- file.path(MC_WORK_PARENT, ifelse(is.null(opts[["o"]]), MC_WORK_DEFAULT, as.character(opts["o"])) )

# Which world to price: the terminated run (SRM switched off at the realization's
# term_delta) or the run where SRM is never terminated. Termination probability
# is NOT applied here -- it is embedded in the frequency of termination times
# drawn in Generate_montecarlo.R (term_delta ~ geometric hazard; a draw at the
# FAIR horizon = no termination within the run). Distinct output names so the two
# modes don't append into the same file.
run_termination = ifelse(is.null(opts[["termination"]]), T, as.logical(opts["termination"]) )
name_output <- if (run_termination) "npc_output.csv" else "npc_noterm_output.csv"
scc_diag_output <- paste0("scc_diag_", str_remove(name_output, "npc_"))
sccnosrm_diag_output <- paste0("sccnosrm_diag_", str_remove(name_output, "npc_"))

# Discard bookkeeping. discarded_runs.csv is the comprehensive per-draw list (both
# gases, reason-tagged: damage_filter / fair_infeasible / fair_infeasible_sibling)
# kept for later statistical analysis; dropped_srmnopulse_scenarios.csv is the
# derived climate-scenario drop list consumed by Analyze_montecarlo_pulse_effects.R.
# Both are termination-independent, so they use fixed names across termination modes.
discarded_runs_output <- "discarded_runs.csv"
discarded_runs_path   <- file.path(output_folder, discarded_runs_output)
dropped_scen_output   <- "dropped_srmnopulse_scenarios.csv"


# Make sure the output folder exists (create it if not)
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

if (any(!dir.exists(res)) ) stop("some of the folder specified do not exsist")
if (!dir.exists(input_folder) ) stop("no id folder")

run_hpc = ifelse(is.null(opts[["hpc"]]), F, as.logical(opts["hpc"]) )
N = ifelse(is.null(opts[["chunk"]]), 100, as.integer(opts["chunk"]) )
skip_scenarios = ifelse(is.null(opts[["skip"]]), T, as.logical(opts["skip"]) )

if(run_hpc==F) {igdx(dirname(Sys.which("gams")))} else {
  igdx(HPC_GAMS_PATH)
  logfile <- file(file.path(output_folder,"r_console.log"), open = "wt")
  sink(logfile)
  sink(logfile, type = "message")  # capture warnings/errors too
}

## climate damage function parameters (defined in all_parameters.R)
gwp <- GWP_INITIAL      # initial world gdp: 2019 World Bank (constant 2015 USD -> 2020)
globaltouspc <- GLOBAL_TO_US_PC  # ratio of global to US per-capita GDP (2020, McDuffie 2023)
#forctoUSD <- forctoTg * TgtoUSD # US$/(W/m^2)
#ozone_rftoconc <- 50 / 0.263  # 50 ppb as https://www.sciencedirect.com/science/article/pii/S2542519622002601?ref=pdf_download&fr=RR-2&rr=9a07c2002fcb708b

# ----------------------------
# Files discovery
# ----------------------------
filelist <- unlist(lapply(res, function(folder) list.files(folder, pattern = ".gdx$", full.names = TRUE)))
filelist <- filelist[str_detect(filelist, "EXP")]


cat("Sanitizing the data...\n")
sanitized_names <- extract_names_dt(filelist)

# Load id_montecarlo from folders
id_list <- lapply(input_folder, function(folder) {
  mc_strip_csv_index(read.csv(file.path(folder, "id_montecarlo.csv"), stringsAsFactors = FALSE))
})
id_montecarlo <- rbindlist(id_list, fill = TRUE)
setDT(id_montecarlo)
validate_id_montecarlo(id_montecarlo, "combined id_montecarlo.csv inputs")
# The FAIR (scenario) parameters drive a FAIR run and are encoded in the gdx
# file names; the post-processing parameters are now sampled once in
# Generate_montecarlo.R and read here straight from id_montecarlo.csv.
fair_cols <- MC_FAIR_COLS
post_cols <- MC_POST_COLS
metadata_cols <- MC_METADATA_COLS
id_cols <- c("ID", fair_cols, post_cols, metadata_cols)
missing_cols <- setdiff(id_cols, names(id_montecarlo))
if (length(missing_cols) > 0) {
  stop("id_montecarlo csv missing columns: ", paste(missing_cols, collapse=", "),
       ". Regenerate it with the current Generate_montecarlo.R.")
}
id_montecarlo <- id_montecarlo[, ..id_cols]

# FAIR/scenario columns are matched against the gdx file names as strings.
for (cname in setdiff(fair_cols, "rcp")) {
  id_montecarlo[, (cname) := as.character(as.integer(round(as.numeric(get(cname)))))]
}
id_montecarlo[, ID := as.character(as.integer(round(as.numeric(ID))))]
id_montecarlo[, rcp := str_remove(as.character(rcp), "RCP")]
# Post-processing columns stay numeric for the cost/damage computations.
for (cname in post_cols) id_montecarlo[, (cname) := as.numeric(get(cname))]
for (cname in c("seed", "draw_index", "term_delay")) id_montecarlo[, (cname) := as.numeric(get(cname))]
id_montecarlo[, censored_termination := as.logical(censored_termination)]
setnames(id_montecarlo, c("pulse","cool","term_delta","start","term"), c("pulse_time","cool_rate","term","geo_start","geo_end"))

# id_montecarlo carries three kinds of column: FAIR/scenario inputs, post-processing
# inputs (post_cols), and provenance metadata (metadata_cols: sampler_version,
# sampling_method, seed, draw_index, term_delay, censored_termination). The metadata
# is not a model input, so it must stay in id_montecarlo.csv and never leak into the
# analysis outputs. Split the table here:
#   * id_montecarlo_full keeps everything and feeds ONLY the discard record, so
#     discarded_runs.csv is a self-contained per-draw dump (inputs + provenance).
#   * id_montecarlo is the "minimum" version used for EVERY downstream merge, so
#     npc/scc outputs carry only inputs + computed results.
# Keyed by ID, so the metadata can always be re-attached from id_montecarlo.csv.
id_montecarlo_full <- copy(id_montecarlo)
id_montecarlo <- id_montecarlo[, setdiff(names(id_montecarlo), metadata_cols), with = FALSE]

cat("Loaded", nrow(id_montecarlo), "montecarlo realizations from id_montecarlo.csv \n")

# Post-processing parameters (incl. the deterministic main-scenario overrides
# for vsl/delta/vsl_eta/theta) are sampled and fixed in Generate_montecarlo.R.

all_names <- c("gas", names(id_montecarlo))

# extract scenario names
scenario_names <- setdiff(names(sanitized_names), c("gdx","file","experiment"))
pulse_scenarios <- setdiff(scenario_names, c("cool_rate","geo_end","geo_start","term"))
base_scenarios <- setdiff(pulse_scenarios, c("gas","pulse_time"))

# srm-only world climate scenario (termination- and gas-independent): the key the
# srmnopulse damage filter aggregates on, and which the pulse-effects script joins.
climate_keys <- c("ecs", "tcr", "rcp", "cool_rate", "pulse_time", "geo_start", "geo_end")
# The comprehensive per-draw discard list is a self-contained dump: every column of
# id_montecarlo (FAIR inputs + post-processing inputs + provenance metadata) plus the
# gas and the reason tag. Discarded draws never reach npc/scc, so this file is the only
# record of their parameters -- which is why, unlike the outputs, it DOES keep metadata.
# `breach_t` is the model timestep at which the relevant world's total damage first
# reaches GDP -- the cut point applied to that draw -- in the SAME units as the
# `pulse` column carried over from id_montecarlo, so the two are directly comparable:
#   reason = damage_truncated       -> srm-world cut, applied to npc_output/scc_output
#   reason = damage_truncated_nosrm -> no-SRM-world cut, applied to sccnosrm_output
#   reason = damage_filter          -> srm-world breach at/before pulse: nothing to price
#   reason = fair_infeasible/not_run-> NA (no breach involved)
# Every record carries the column (NA where not applicable) because fwrite(append=TRUE)
# writes rows positionally and would silently misalign a ragged schema.
append_discards <- function(x) {
  if (!"breach_t" %in% names(x)) x <- copy(x)[, breach_t := NA_real_]
  rec <- unique(x[, .(ID = as.character(ID), gas = as.character(gas),
                      reason = as.character(reason), breach_t = as.numeric(breach_t))])
  rec <- merge(rec, id_montecarlo_full, by = "ID", all.x = TRUE, allow.cartesian = TRUE)
  setcolorder(rec, c(names(id_montecarlo_full), "gas", "reason", "breach_t"))
  # fwrite(append=TRUE) writes rows positionally without checking the header, so a file
  # left over from an older schema (e.g. one predating breach_t, kept across a --skip=T
  # resume) would be silently corrupted. Fail loudly instead.
  if (file.exists(discarded_runs_path)) {
    existing_cols <- names(fread(discarded_runs_path, nrows = 0L))
    if (!identical(existing_cols, names(rec))) {
      stop("Existing ", discarded_runs_output, " has columns [",
           paste(existing_cols, collapse = ", "), "] but this run appends [",
           paste(names(rec), collapse = ", "), "]. Delete the file or rerun with --skip=F.",
           call. = FALSE)
    }
  }
  fwrite(rec, file = discarded_runs_path, append = file.exists(discarded_runs_path))
}


## pre-filter by experiment so we don't repeat this in the loop
dt_srmpt   <- sanitized_names[experiment == "srmpulsemaskedterm"]
dt_srmptm  <- sanitized_names[experiment == "srmpulsemasked"]
dt_srmp    <- sanitized_names[experiment == "srmpulse"]
dt_srm     <- sanitized_names[experiment == "srm"]
dt_pulse   <- sanitized_names[experiment == "pulse"]
dt_base    <- sanitized_names[experiment == "base"]

## define column sets used as join keys
k_all     <- scenario_names
k_no_term <- setdiff(scenario_names, "term")
k_pulse   <- setdiff(scenario_names, c("term", "cool_rate", "geo_start", "geo_end"))
k_base    <- setdiff(scenario_names,
                     c("term", "cool_rate", "geo_start", "geo_end", "gas", "pulse_time"))

## set keys once for fast joins inside the loop
setkeyv(dt_srmpt,   k_all)
setkeyv(dt_srmptm,  k_no_term)
setkeyv(dt_srmp,    k_no_term)
setkeyv(dt_srm,     k_no_term)
setkeyv(dt_pulse,   k_pulse)
setkeyv(dt_base,    k_base)

completed_runs <- data.table()
if(file.exists(file.path(output_folder, name_output))) { 
  existing_output <- fread(file = file.path(output_folder, name_output))
  if (skip_scenarios == F) {
    file.remove(file.path(output_folder, name_output))
    file.remove(file.path(output_folder, paste0("scc_",str_remove(name_output,"npc_"))) )
    file.remove(file.path(output_folder, paste0("sccnosrm_",str_remove(name_output,"npc_"))) )
    file.remove(file.path(output_folder, scc_diag_output))
    file.remove(file.path(output_folder, sccnosrm_diag_output))
    # The discard list is a per-run output too: leaving it in place would append this
    # run's records to a previous run's (stale draws, and a stale column schema).
    if (file.exists(discarded_runs_path)) file.remove(discarded_runs_path)
  } else {
    if (!all(c("ID", "gas") %in% names(existing_output))) {
      stop("Existing ", name_output, " is legacy output without ID/gas keys. Rerun with --skip=F to replace it.")
    }
    completed_runs <- unique(existing_output[, .(ID = as.character(ID), gas = as.character(gas))])
  }
} else {skip_scenarios <- F}

# Keep the discard list in step with the main output: a fresh/replace run (--skip=F)
# rebuilds it from scratch; a resumed run (--skip=T) appends to it and
# write_dropped_scenarios() de-duplicates at the end.
if (!skip_scenarios && file.exists(discarded_runs_path)) file.remove(discarded_runs_path)

expected_gases <- sort(unique(na.omit(sanitized_names$gas)))
if (length(expected_gases) == 0L) stop("No gas-specific GDX files found in results folders.")

id_scenarios <- unique(id_montecarlo[, c("ID", setdiff(scenario_names, "gas")), with = FALSE])
full_expected <- id_scenarios[, .(gas = expected_gases), by = names(id_scenarios)]
setcolorder(full_expected, c("ID", scenario_names))

# Map every (draw, gas) to its GDX experiment files over the FULL grid (independent
# of --skip) so FAIR-run infeasibility is judged across BOTH gases: a draw is usable
# only if every gas produced a complete experiment set. If any gas is missing /
# infeasible, the whole draw is discarded (both gases) -- if one gas completed but
# the other did not, both are dropped -- and recorded to the shared discard list.
runs_all <- copy(full_expected)
runs_all[, f_srmpulsemaskedterm := dt_srmpt[runs_all, on = k_all, mult = "first"]$gdx]
runs_all[, f_srmpulsemasked     := dt_srmptm[runs_all, on = k_no_term, mult = "first"]$gdx]
runs_all[, f_srmpulse           := dt_srmp[runs_all, on = k_no_term, mult = "first"]$gdx]
runs_all[, f_srm                := dt_srm[runs_all, on = k_no_term, mult = "first"]$gdx]
runs_all[, f_pulse              := dt_pulse[runs_all, on = k_pulse, mult = "first"]$gdx]
runs_all[, f_base               := dt_base[runs_all, on = k_base, mult = "first"]$gdx]

file_cols <- grep("^f_", names(runs_all), value = TRUE)
runs_all[, complete_run := stats::complete.cases(runs_all[, ..file_cols])]
runs_all[, draw_ok := all(complete_run), by = "ID"]

# Monte Carlo jobs are launched in id_montecarlo order (ID == draw_index == Sobol
# order) as a contiguous prefix, so the largest ID that produced ANY gdx marks how
# far the launcher actually got. Draws past this frontier were never submitted and
# are tagged "not_run" (rather than "fair_infeasible", i.e. submitted-but-failed),
# so an incomplete/partial launch doesn't masquerade as thousands of infeasible runs.
# Only the SRM-family experiments are draw-specific evidence of a launch: f_base is
# keyed on (rcp,ecs,tcr) alone and f_pulse adds only (gas,pulse_time), so both are
# shared by many later, never-submitted draws and would inflate the frontier.
srm_file_cols <- intersect(c("f_srm", "f_srmpulse", "f_srmpulsemasked", "f_srmpulsemaskedterm"), file_cols)
runs_all[, any_file := Reduce(`|`, lapply(.SD, Negate(is.na))), .SDcols = srm_file_cols]
attempted_frontier <- suppressWarnings(max(as.integer(runs_all[any_file == TRUE, ID])))
if (!is.finite(attempted_frontier)) attempted_frontier <- 0L
runs_all[, any_file := NULL]

infeasible_runs <- runs_all[draw_ok == FALSE]
if (nrow(infeasible_runs) > 0L) {
  infeasible_runs[, missing_experiments := apply(.SD, 1, function(z) paste(names(z)[is.na(z)], collapse = ", ")), .SDcols = file_cols]
  inf_rec <- infeasible_runs[, c("ID", scenario_names), with = FALSE]
  # not_run: draw beyond the launch frontier, i.e. never submitted. Otherwise the
  # FAIR run was submitted but did not produce a complete set for every gas: the
  # failing gas -> fair_infeasible; a gas that completed but is dropped only because
  # its sibling failed -> fair_infeasible_sibling.
  inf_rec[, reason := fifelse(as.integer(ID) > attempted_frontier, "not_run",
                       fifelse(infeasible_runs$complete_run, "fair_infeasible_sibling", "fair_infeasible"))]
  n_not_run <- uniqueN(inf_rec[reason == "not_run", ID])
  n_failed  <- uniqueN(inf_rec[reason != "not_run", ID])
  # Genuine solve failures (submitted but incomplete) warrant a warning; never-run
  # draws are just an incomplete launch, so they only get an informational note.
  if (n_failed > 0L) {
    details <- paste(capture.output(print(head(
      unique(infeasible_runs[as.integer(ID) <= attempted_frontier,
             c("ID", scenario_names, "missing_experiments"), with = FALSE]), 20))), collapse = "\n")
    warning(n_failed, " Monte Carlo draw(s) submitted but missing a complete GDX ",
            "experiment set for at least one gas; BOTH gases dropped (fair_infeasible):\n", details)
  }
  if (n_not_run > 0L) {
    cat("Note:", n_not_run, "Monte Carlo draw(s) beyond the launch frontier (ID >",
        attempted_frontier, ") were never submitted; recorded as not_run.\n")
  }
  # DROP the whole draw regardless (leaving it in breaks the downstream key-set
  # assertions) and record it for later analysis / the pulse-effects script.
  append_discards(inf_rec[, .(ID, gas, reason)])
}

feasible_runs <- runs_all[draw_ok == TRUE]
feasible_runs[, c("complete_run", "draw_ok") := NULL]
if (nrow(feasible_runs) == 0L) {
  stop("No Monte Carlo draws have a complete set of GDX experiment files for every gas; nothing to analyze.")
}

# Apply the --skip (resume) filter to the feasible draws only.
required_runs <- copy(feasible_runs)
if (skip_scenarios == TRUE && nrow(completed_runs) > 0L) {
  required_runs <- required_runs[!completed_runs, on = c("ID", "gas")]
}
if (nrow(required_runs) == 0L) {
  cat("No new Monte Carlo runs to analyze after applying existing-output skip.\n")
  write_dropped_scenarios()
  quit(save = "no", status = 0)
}

all_experiments <- unique(required_runs[, ..scenario_names])

cat("End data preprocessing...\n")

for (n_chunk in seq(1,nrow(all_experiments), by = N+1 )) { 
  
exp_chunk <- all_experiments[n_chunk:min(n_chunk+N,nrow(all_experiments))]
runs_chunk <- required_runs[exp_chunk, on = scenario_names, nomatch = 0L]
if (nrow(runs_chunk) == 0L) next
expected_chunk_keys <- unique(runs_chunk[, .(ID = as.character(ID), gas = as.character(gas))])

files_loop <- unique(na.omit(unlist(runs_chunk[, ..file_cols], use.names = FALSE)))
if (length(files_loop) == 0) next     # skip loop iteration

cat("Loading the data for experiments", n_chunk, "to", n_chunk+N, "from a total of", nrow(all_experiments), "experiments...\n")

# create full grid to handle missing data
all_t <- seq(1, T_HORIZON)  # or tot_forcing$t
all_scenarios <- unique(sanitized_names[gdx %in% files_loop, ..scenario_names])  # all scenario combinations

# Create the full grid
full_grid <- CJ(t = all_t, 
                scenario = 1:nrow(all_scenarios), unique = TRUE)

# Merge scenario columns back
full_grid <- cbind(full_grid[, -"scenario", with = FALSE], all_scenarios[full_grid$scenario])

# Batch extract using gdxtools, then join sanitized_names and sanitize
TATM <- setDT(gdxtools::batch_extract("TATM", files = files_loop)$TATM)
TATM <- merge(TATM, sanitized_names, by.x = "gdx", by.y = "gdx", all = FALSE)
TATM <- sanitize_dt(TATM)

W_EMI <- setDT(gdxtools::batch_extract("W_EMI", files = files_loop)$W_EMI)
W_EMI <- merge(W_EMI, sanitized_names, by = "gdx", all = FALSE)
W_EMI <- sanitize_dt(W_EMI)

FORC <- setDT(gdxtools::batch_extract("FORCING", files = files_loop)$FORCING)
FORC <- merge(FORC, sanitized_names, by = "gdx", all = FALSE)
FORC <- sanitize_dt(FORC)

CONC <- setDT(gdxtools::batch_extract("CONC", files = files_loop)$CONC)
CONC <- merge(CONC, sanitized_names, by = "gdx", all = FALSE)
CONC <- sanitize_dt(CONC)

SRM <- setDT(gdxtools::batch_extract("SRM", files = files_loop)$SRM)
SRM <- merge(SRM, sanitized_names, by = "gdx", all = FALSE)
SRM <- sanitize_dt(SRM)

background_srm <- setDT(gdxtools::batch_extract("forcing_srm", files = files_loop)$forcing_srm)
background_srm <- merge(background_srm, sanitized_names, by = "gdx", all = FALSE)
background_srm <- sanitize_dt(background_srm)

# total forcing aggregation (equivalent to dplyr group_by + summarise)
tot_forcing <- FORC[, .(value = sum(value)), by = c("t", "file",scenario_names, "experiment")]

cat("Analyzing the data... \n")

# NPC of a single world over the whole post-deployment horizon [pulse_time, 480]:
# either the terminated run (srmpulsemaskedterm: SRM off for t > term_delta) or
# the run where SRM is never terminated (srmpulsemasked). No probability blend --
# the termination distribution is already carried by the term_delta draws.
selected_experiment <- if (run_termination) "srmpulsemaskedterm" else "srmpulsemasked"
dt_sel <- prepare_join_table(selected_experiment)

# --- srmnopulse damage truncation (full horizon) ------------------------------
# Locate, per draw, the first year in which the srm-only world's total damage
# reaches GDP. Draws are NOT rejected on this basis: they are flagged in
# discarded_runs.csv (reason "damage_truncated") and every timestep from the breach
# onward is dropped, so the NPVs below run over each draw's economically meaningful
# window only. `srm_breach` is carried to the scc blocks so all three outputs
# (npc, scc, scc_nosrm) are truncated at the SAME year for a given draw: Figure_2
# forms per-draw ratios of these NPVs (npc_srm / scc_nosrm, npc_srm / scc_srm), which
# are only meaningful if every term spans the same horizon.
#
# The one case that still has to be rejected: a draw breaching at or before its own
# pulse_time has no valid window left at all (the NPV would be an empty sum), so it
# keeps reason "damage_filter" and is dropped as before.
srm_breach <- srmnopulse_first_breach_t(dt_sel, "ID", climate_keys)
# Collapse to exactly one row per ID (earliest breach wins) so the lookup in
# truncate_at_breach cannot resolve an ID to an arbitrary one of several rows.
srm_breach <- srm_breach[, .(first_breach_t = if (all(is.na(first_breach_t))) NA_real_
                             else min(first_breach_t, na.rm = TRUE)), by = ID]
mc_assert_unique_key(srm_breach, "ID", paste0("srmnopulse breach table chunk ", n_chunk))
pulse_by_id <- unique(dt_sel[, .(ID, pulse_time = as.numeric(pulse_time))])
srm_breach <- merge(srm_breach, pulse_by_id, by = "ID", all.x = TRUE)

empty_ids <- srm_breach[!is.na(first_breach_t) & first_breach_t <= pulse_time, unique(ID)]
trunc_ids <- srm_breach[!is.na(first_breach_t) & first_breach_t >  pulse_time, unique(ID)]

if (length(trunc_ids) > 0L) {
  trunc_rec <- unique(dt_sel[ID %in% trunc_ids, c("ID", scenario_names), with = FALSE])
  trunc_rec[, reason := "damage_truncated"]
  trunc_rec <- merge(trunc_rec, srm_breach[, .(ID, breach_t = first_breach_t)],
                     by = "ID", all.x = TRUE)
  append_discards(trunc_rec[, .(ID, gas, reason, breach_t)])
  cat("srmnopulse truncation (chunk", n_chunk, "):", length(trunc_ids),
      "draw(s) flagged and truncated at first year where total srm-run damage >= GDP",
      "(kept in all outputs).\n")
}
if (length(empty_ids) > 0L) {
  dmg_rec <- unique(dt_sel[ID %in% empty_ids, c("ID", scenario_names), with = FALSE])
  dmg_rec[, reason := "damage_filter"]
  dmg_rec <- merge(dmg_rec, srm_breach[, .(ID, breach_t = first_breach_t)],
                   by = "ID", all.x = TRUE)
  append_discards(dmg_rec[, .(ID, gas, reason, breach_t)])
  cat("srmnopulse filter (chunk", n_chunk, "): dropping",
      expected_chunk_keys[ID %in% empty_ids, .N], "run(s) across",
      length(empty_ids), "draw(s) breaching at or before their own pulse_time.\n")
  dt_sel <- dt_sel[!ID %in% empty_ids]
  expected_chunk_keys <- expected_chunk_keys[!ID %in% empty_ids]
  srm_breach <- srm_breach[!ID %in% empty_ids]
}
if (nrow(expected_chunk_keys) == 0L) {
  cat("All runs in chunk", n_chunk, "dropped by the srmnopulse filter; skipping chunk.\n")
  next
}

dt_sel <- dt_sel[ t >= as.numeric(pulse_time) ]
mc_assert_key_set_equal(expected_chunk_keys, unique(dt_sel[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("NPC selected input chunk ", n_chunk))
mc_assert_complete_time_window(dt_sel, c("ID", "gas"),
                               paste0("NPC selected input chunk ", n_chunk))
# Truncate AFTER the completeness assertion: the assertion validates that the gdx
# delivered a full [pulse_time, T_HORIZON] window, which is still what we require
# of the raw data. The truncation below is a deliberate, per-draw shortening of
# that window and must not trip it.
dt_sel <- truncate_at_breach(dt_sel, srm_breach)
res_sel <- npv_aggregator(dt_sel, all_names)

# pulse_size & scc (reuse earlier approach but modularized)
pulse_size <- W_EMI[ ghg == gas & experiment %in% c("srm","srmpulse"), .(t, value, experiment, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[, .(pulse_size = value[ experiment == "srmpulse" ] - value[ experiment == "srm" ]), by = .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[ pulse_size != 0 ]
pulse_size[ , pulse_size := ifelse(gas == "co2", pulse_size * 1e9, pulse_size * 1e6)]
pulse_size[, t:=NULL]
mc_assert_no_missing(pulse_size, c("pulse_size", "gas", "rcp", "ecs", "tcr", "cool_rate", "pulse_time", "geo_start", "geo_end"),
                     paste0("pulse_size chunk ", n_chunk))

# scc calculation (with SRM)
scc <- merge(
  TATM[experiment == "srm", .(t, temp_srm = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  TATM[experiment == "srmpulse", .(t, temp_srmpulse = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE
)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "srmpulse", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, concch4_pulse = value)], by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "srm", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, concch4_base = value)],  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all.x = TRUE)
scc <- merge(scc, tot_forcing[experiment == "srm", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, forc_srm = value)],
             by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all.x = TRUE)
scc <- merge(scc, tot_forcing[experiment == "srmpulse", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, forc_srmpulse = value)],
             by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all.x = TRUE)
mc_assert_no_missing(scc, c("temp_srm", "temp_srmpulse", "concch4_pulse", "concch4_base",
                            "forc_srm", "forc_srmpulse"),
                     paste0("SCC-with-SRM climate inputs chunk ", n_chunk))
scc <- merge(scc, id_montecarlo[, !c("theta","term","prob","mortality_srm","forctoTg"), with = FALSE], by = intersect(names(scc), names(id_montecarlo)), allow.cartesian = TRUE)
mc_assert_no_missing(scc, c("ID", "alpha", "delta", "mortality_ozone", "vsl", "vsl_eta", "dg"),
                     paste0("SCC-with-SRM Monte Carlo inputs chunk ", n_chunk))
scc <- scc[ t >= as.numeric(pulse_time) ]
# The id_montecarlo join keys on the SRM scenario without `term`, so it can pull
# in realizations whose termination time places their full scenario in a later
# chunk. Restrict to this chunk's realizations (see the no-SRM note below) so
# nothing is double-counted across chunks.
scc <- scc[expected_chunk_keys, on = c("ID", "gas"), nomatch = 0L]
scc[, `:=`(
  dtemp_srmpulse_srm = temp_srmpulse - temp_srm,
  dforc_srmpulse_srm = forc_srmpulse - forc_srm
)]
mc_assert_key_set_equal(expected_chunk_keys, unique(scc[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("SCC-with-SRM input chunk ", n_chunk))
mc_assert_complete_time_window(scc, c("ID", "gas"),
                               paste0("SCC-with-SRM input chunk ", n_chunk))
scc <- truncate_at_breach(scc, srm_breach)   # same cut point as the npc block
scc_by <- c("gas", intersect(names(scc), names(id_montecarlo)))
scc_diag <- make_scc_diagnostics(
  scc,
  scc_by,
  c("forc_srm", "forc_srmpulse", "dforc_srmpulse_srm",
    "temp_srm", "temp_srmpulse", "dtemp_srmpulse_srm"),
  paste0("SCC-with-SRM chunk ", n_chunk)
)
compute_gwpt(scc)
scc[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp) ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
scc[, dam := gwpt * ( alpha * (pmax(0,temp_srmpulse)^2 - pmax(0,temp_srm)^2) ) ]
mc_assert_no_missing(scc, c("tropoz_pollution", "dam"), paste0("SCC-with-SRM computed costs chunk ", n_chunk))
scc_agg <- scc[, .(damnpv = sum( dam / (1 + delta)^(t - as.numeric(pulse_time))),
                    ozpnpv = sum( tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)))),
                by = scc_by ]
scc_agg <- merge(scc_agg, pulse_size, by = intersect(names(scc_agg), names(pulse_size)), all.x = TRUE)
scc_agg[, scc := (damnpv + ozpnpv) / pulse_size]
mc_assert_no_missing(scc_agg, c("ID", "gas", "pulse_size", "scc"), paste0("SCC-with-SRM output chunk ", n_chunk))
mc_assert_unique_key(scc_agg, c("ID", "gas"), paste0("SCC-with-SRM output chunk ", n_chunk))
mc_assert_key_set_equal(expected_chunk_keys, unique(scc_agg[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("SCC-with-SRM output chunk ", n_chunk))

# scc calculation (w/o SRM)
scc <- merge(
  TATM[experiment == "base", .(t, temp_base = value, rcp, ecs, tcr)],
  TATM[experiment == "pulse", .(t, temp_pulse = value, gas, rcp, ecs, tcr, pulse_time)],
  by = c("t","rcp","ecs","tcr"), all = FALSE
)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "pulse", .(t, gas, rcp, ecs, tcr, pulse_time, concch4_pulse = value)], by = c("t","gas","rcp","ecs","tcr","pulse_time"), all = FALSE)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "base", .(t, rcp, ecs, tcr, concch4_base = value)],  by = c("t","rcp","ecs","tcr"), all.x = TRUE)
scc <- merge(scc, tot_forcing[experiment == "base", .(t, rcp, ecs, tcr, forc_base = value)],
             by = c("t","rcp","ecs","tcr"), all.x = TRUE)
scc <- merge(scc, tot_forcing[experiment == "pulse", .(t, gas, rcp, ecs, tcr, pulse_time, forc_pulse = value)],
             by = c("t","gas","rcp","ecs","tcr","pulse_time"), all.x = TRUE)
mc_assert_no_missing(scc, c("temp_base", "temp_pulse", "concch4_pulse", "concch4_base",
                            "forc_base", "forc_pulse"),
                     paste0("SCC-no-SRM climate inputs chunk ", n_chunk))
scc <- merge(scc, id_montecarlo[, c("ID","ecs","tcr","rcp","pulse_time","alpha","delta",
                                    "mortality_ozone","vsl","vsl_eta","dg"), with = FALSE],
             by = intersect(names(scc), names(id_montecarlo)), allow.cartesian = TRUE)
mc_assert_no_missing(scc, c("ID", "alpha", "delta", "mortality_ozone", "vsl", "vsl_eta", "dg"),
                     paste0("SCC-no-SRM Monte Carlo inputs chunk ", n_chunk))
# Breach year of the no-SRM world judged by its OWN total damage, computed here while
# the table still spans the full horizon (the srm-world equivalent is likewise taken
# before any pulse_time restriction). A draw whose no-SRM world breaches at or before
# its pulse_time has no economically meaningful window and no definable scc_nosrm, so
# it is excluded from this block's expected key set rather than emitted as a zero.
nosrm_breach <- nosrm_first_breach_t(scc, "ID")
nosrm_breach <- merge(nosrm_breach, pulse_by_id, by = "ID", all.x = TRUE)
mc_assert_unique_key(nosrm_breach, "ID", paste0("no-SRM breach table chunk ", n_chunk))
nosrm_empty_ids <- nosrm_breach[!is.na(first_breach_t) & first_breach_t <= pulse_time, unique(ID)]
expected_nosrm_keys <- expected_chunk_keys[!ID %in% nosrm_empty_ids]

nosrm_trunc_ids <- nosrm_breach[!is.na(first_breach_t) & first_breach_t > pulse_time, unique(ID)]
# Build the record from expected_chunk_keys, not from `scc`: the id_montecarlo join
# above pulls in every realization sharing this base scenario, including IDs belonging
# to other chunks, and only this chunk's keys should be flagged here.
nosrm_rec <- expected_chunk_keys[ID %in% nosrm_trunc_ids, .(ID, gas)]
if (nrow(nosrm_rec) > 0L) {
  nosrm_rec[, reason := "damage_truncated_nosrm"]
  # expected_chunk_keys carries ID as character, nosrm_breach as numeric -- match on
  # character so the breach year actually lands instead of merging to all-NA.
  nosrm_rec <- merge(nosrm_rec,
                     nosrm_breach[, .(ID = as.character(ID), breach_t = first_breach_t)],
                     by = "ID", all.x = TRUE)
  mc_assert_no_missing(nosrm_rec, "breach_t", paste0("no-SRM discard record chunk ", n_chunk))
  append_discards(nosrm_rec[, .(ID, gas, reason, breach_t)])
  cat("no-SRM truncation (chunk", n_chunk, "):", uniqueN(nosrm_rec$ID),
      "draw(s) truncated at their own no-SRM breach year in sccnosrm_output (kept).\n")
}
if (length(nosrm_empty_ids) > 0L) {
  cat("no-SRM filter (chunk", n_chunk, "):", length(nosrm_empty_ids),
      "draw(s) have no definable scc_nosrm (no-SRM world breaches at/before pulse_time);",
      "omitted from sccnosrm_output only.\n")
}

scc <- scc[ t >= as.numeric(pulse_time) ]
# The no-SRM SCC is keyed only on the base FAIR scenario (ecs/tcr/rcp/pulse_time),
# so the id_montecarlo join above pulls in every realization sharing that base
# scenario -- including IDs whose SRM dimensions (cool_rate/geo/term) place their
# full scenario in a different chunk. Keep only this chunk's realizations so the
# same ID isn't aggregated and appended again from its own chunk (duplicate rows)
# and the key-set assertion below holds.
scc <- scc[expected_chunk_keys, on = c("ID", "gas"), nomatch = 0L]
scc[, `:=`(
  dtemp_pulse_base = temp_pulse - temp_base,
  dforc_pulse_base = forc_pulse - forc_base
)]
mc_assert_key_set_equal(expected_chunk_keys, unique(scc[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("SCC-no-SRM input chunk ", n_chunk))
mc_assert_complete_time_window(scc, c("ID", "gas"),
                               paste0("SCC-no-SRM input chunk ", n_chunk))
# Cut at the no-SRM world's OWN breach year, not the srm world's: this block prices a
# world with no SAI in it, so the srm world's viability is not the relevant standard.
# The two cut points genuinely differ -- the no-SRM world is hotter and usually breaches
# first, but for draws with costly SAI (high mortality_srm / theta) the srm world's
# deployment and pollution terms tip it over earlier instead.
# CONSEQUENCE: npc/scc and scc_nosrm are now integrated over different horizons for the
# same draw, so the per-draw ratios in Figure_2 (npc_srm / scc_nosrm, npc_srm / scc_srm)
# divide NPVs spanning different windows. That is a deliberate choice -- see the SI note
# on draw-dependent horizons -- but it is the reason those ratios are not like-for-like.
scc <- truncate_at_breach(scc, nosrm_breach)
scc_nosrm_by <- c("gas", intersect(names(scc), names(id_montecarlo)))
scc_nosrm_diag <- make_scc_diagnostics(
  scc,
  scc_nosrm_by,
  c("forc_base", "forc_pulse", "dforc_pulse_base",
    "temp_base", "temp_pulse", "dtemp_pulse_base"),
  paste0("SCC-no-SRM chunk ", n_chunk)
)
compute_gwpt(scc)
scc[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp) ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
scc[, dam := gwpt * ( alpha * (pmax(0,temp_pulse)^2 - pmax(0,temp_base)^2) )]
mc_assert_no_missing(scc, c("tropoz_pollution", "dam"), paste0("SCC-no-SRM computed costs chunk ", n_chunk))
scc_agg2 <- scc[, .(damnpv = sum( dam / (1 + delta)^(t - as.numeric(pulse_time))),
                    ozpnpv = sum( tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)))),
                by = scc_nosrm_by ]
scc_agg2 <- merge(scc_agg2, pulse_size[, c("gas","ecs","tcr","rcp","pulse_time","pulse_size"), with = FALSE], by = intersect(names(scc_agg2), names(pulse_size)), all.x = TRUE)
scc_agg2[, scc_nosrm := (damnpv + ozpnpv) / pulse_size]
mc_assert_no_missing(scc_agg2, c("ID", "gas", "pulse_size", "scc_nosrm"), paste0("SCC-no-SRM output chunk ", n_chunk))
mc_assert_unique_key(scc_agg2, c("ID", "gas"), paste0("SCC-no-SRM output chunk ", n_chunk))
# Assert against the narrowed key set: draws whose no-SRM world breaches at or before
# their pulse_time have no scc_nosrm and are legitimately absent here, while every other
# draw must still be present. Narrowing the expectation keeps this a strict equality
# check rather than weakening it to a subset test.
mc_assert_key_set_equal(expected_nosrm_keys, unique(scc_agg2[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("SCC-no-SRM output chunk ", n_chunk))

# combine results: a single experiment's NPC, no probability weighting
combined <- res_sel[, c(all_names, "source", "cost"), with = FALSE]
mc_assert_no_missing(combined, c("ID", "gas", "source", "cost"), paste0("NPC long output chunk ", n_chunk))
if (nrow(combined) == 0L) {
  stop(
    paste0(
      "No aggregated NPV rows were produced for chunk starting at ", n_chunk,
      ". This usually means the matched GDX files for this chunk are incomplete or prepare_join_table() returned no joined rows."
    )
  )
}
combined_wide <- dcast(
  combined,
  formula = as.formula(paste(paste(all_names, collapse = " + "), "~ source")),
  value.var = "cost"
)
required_sources <- c("dirnpv", "srmpnpv", "ozpnpv", "masknpv", "damnpv")
missing_sources <- setdiff(required_sources, names(combined_wide))
if (length(missing_sources) > 0L) {
  stop(
    paste0(
      "Missing aggregated cost columns after dcast for chunk starting at ", n_chunk,
      ": ", paste(missing_sources, collapse = ", "),
      ". This happens when one or more source categories are absent in the chunk input."
    )
  )
}
combined_wide <- merge(combined_wide, unique(pulse_size), by = intersect(names(combined_wide), names(pulse_size)), all.x = TRUE)
combined_wide[, npc_srm := (dirnpv + srmpnpv + ozpnpv + masknpv + damnpv) / pulse_size]
mc_assert_no_missing(combined_wide, c("ID", "gas", "pulse_size", "npc_srm"), paste0("NPC output chunk ", n_chunk))
mc_assert_unique_key(combined_wide, c("ID", "gas"), paste0("NPC output chunk ", n_chunk))
mc_assert_key_set_equal(expected_chunk_keys, unique(combined_wide[, .(ID = as.character(ID), gas = as.character(gas))]),
                        c("ID", "gas"), paste0("NPC output chunk ", n_chunk))

# ----------------------------
# Save final output (combined table)
# ----------------------------
fwrite(combined_wide, file = file.path(output_folder, name_output), append = TRUE)
fwrite(scc_agg, file = file.path(output_folder, paste0("scc_",str_remove(name_output,"npc_") )), append = TRUE)
fwrite(scc_agg2, file = file.path(output_folder, paste0("sccnosrm_",str_remove(name_output,"npc_"))), append = TRUE)
fwrite(scc_diag, file = file.path(output_folder, scc_diag_output), append = TRUE)
fwrite(scc_nosrm_diag, file = file.path(output_folder, sccnosrm_diag_output), append = TRUE)

}

# Now that every chunk's discards (damage filter) and the pre-loop discards
# (FAIR infeasibility) have been recorded, build the shared climate-scenario drop
# list for Analyze_montecarlo_pulse_effects.R and de-duplicate the discard list.
write_dropped_scenarios()

cat("Saved output to:", file.path(output_folder, name_output), "\n")

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("This code run for exactly",time.taken, "seconds")

if(run_hpc==T) {
  sink(type = "message")
  sink()
  close(logfile) }
