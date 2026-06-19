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
  # ensure t numeric and <=480, drop gdx
  if("t" %in% names(DT)) DT[, t := as.numeric(t)]
  if("gdx" %in% names(DT)) DT <- DT[, !"gdx", with = FALSE]
  DT <- DT[ t <= 480 ]
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
compute_gwpt <- function(DT, dg_col = "dg", t_col = "t", gwp_max = 15) {
  # gwpt = pmin(gwp * (1+dg)^(280-1), gwp * (1+dg)^(t-1))
  DT[t <= 80, gwpt := pmin(gwp * gwp_max,
                    gwp * (1 + dg)^(t-1))]
  DT[t > 80 & t <= 180, gwpt := pmin(gwp * gwp_max,
                  gwp * (1 + dg)^(80-1) * (1 + dg/2)^(t-80))]
  DT[t > 180, gwpt := pmin(gwp * gwp_max,
                           gwp * (1 + dg)^(80-1) * (1 + dg/2)^(180-80) * (1 + dg/4)^(t-180))]
}

npv_aggregator <- function(DT, keep_names = all_names) {
  DT = copy(DT)
  compute_gwpt(DT)
  DT[, dam := gwpt * (alpha * (pmax(0,temp)^2 - pmax(0,temp_srm)^2))]
  DT[, direct_cost := srm_masking * forctoTg * TgtoUSD]
  DT[, srm_pollution := vsl * globaltouspc * (gwpt/gwp)  ^ (vsl_eta) * srm_masking * forctoTg * mortality_srm]
  DT[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp)  ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
  DT[, imperfect_masking := 2 * gwpt * alpha * pmax(0,temp_base) * (1 - cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1 + srm/forc)]
  # aggregator
  res <- DT[, .(
    dirnpv = sum(direct_cost / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    srmpnpv = sum(srm_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    ozpnpv = sum(tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    masknpv = sum(imperfect_masking / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    damnpv = sum(dam / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE)
  ), by = keep_names]
  # convert to long format
  res_long <- melt(
    res,
    id.vars = keep_names,
    variable.name = "source",
    value.name = "cost")
  return(res_long)
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
  base <- dsrm_exp[base, on = c("t",scenario_names)]
  base <- concch4_exp[base, on = c("t", scenario_names)]
  base <- concch4_srm[base, on = setdiff(c("t", scenario_names),"term") ]
  base <- temp_srm[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_srmpulse[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_base[base, on = c("t", base_scenarios)]
  base <- temp_ghgpulse[base, on = c("t",pulse_scenarios), nomatch = NULL]
  
  replace_num_na0(base)
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
  return(base)
}

# Usage
'Launch script to analyze montecarlo scenarios (produces a csv file in the same folder)

Usage:
  Analyze_montecarlo.R [-i <input>] [-o <results>] [--hpc <run_hpc>] [-p <plot_results>] [--chunk <chunk>] [--skip <skip>] [--res <output_folder>] [--termination <termination>]

Options:
-i <input>             Path where the montecarlo id are
-o <results>           Where to save results
--res <output_folder>  name of the output folder (default: Results_montecarlo). For multiple folders separate with -
--hpc <run_hpc>        T/F if running from Juno (T) or local (F)
--chunk <chunk>        how many scenarios to run together
--skip <skip>          skip or rerun existing scenarios
--termination <termination>  T: NPC of the terminated run (SRM off at term_delta); F: NPC of the run where SRM is never terminated (default T)
' -> doc

opts <- docopt(doc, version = 'Montecarlo')

res <- ifelse(is.null(opts[["res"]]), "Results_montecarlo", as.character(opts["res"]) )
res <- str_split(res,"-")[[1]]

input_folder <- ifelse(is.null(opts[["i"]]), "Montecarlo", as.character(opts["i"]) )
output_folder <- ifelse(is.null(opts[["o"]]), "Montecarlo", as.character(opts["o"]) )

# Which world to price: the terminated run (SRM switched off at the realization's
# term_delta) or the run where SRM is never terminated. Termination probability
# is NOT applied here -- it is embedded in the frequency of termination times
# drawn in Generate_montecarlo.R (term_delta ~ geometric hazard; a draw at the
# FAIR horizon = no termination within the run). Distinct output names so the two
# modes don't append into the same file.
run_termination = ifelse(is.null(opts[["termination"]]), T, as.logical(opts["termination"]) )
name_output <- if (run_termination) "npc_output.csv" else "npc_noterm_output.csv"


# Make sure the output folder exists (create it if not)
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

if (any(!dir.exists(res)) ) stop("some of the folder specified do not exsist")
if (!dir.exists(input_folder) ) stop("no id folder")

run_hpc = ifelse(is.null(opts[["hpc"]]), F, as.logical(opts["hpc"]) )
N = ifelse(is.null(opts[["chunk"]]), 100, as.integer(opts["chunk"]) )
skip_scenarios = ifelse(is.null(opts[["skip"]]), T, as.logical(opts["skip"]) )

if(run_hpc==F) {igdx()} else {
  igdx("/work/cmcc/pa12520/gams40.4_linux_x64_64_sfx")
  logfile <- file(file.path(output_folder,"r_console.log"), open = "wt")
  sink(logfile)
  sink(logfile, type = "message")  # capture warnings/errors too
}

## climate damage function parameters
gwp <- 105*1e12 # initial world gdp
#forctoTg <- 1/0.2 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
globaltouspc <- 1.82 * 10937 / 63515 # ratio between global GDPc and USgdpc in 2020, adjusted to match McDuffie 2023 
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
  read.csv(file.path(folder, "id_montecarlo.csv"), stringsAsFactors = FALSE)
})
id_montecarlo <- rbindlist(id_list, fill = TRUE)
setDT(id_montecarlo)
# The FAIR (scenario) parameters drive a FAIR run and are encoded in the gdx
# file names; the post-processing parameters are now sampled once in
# Generate_montecarlo.R and read here straight from id_montecarlo.csv.
fair_cols <- c("ecs","tcr","rcp","pulse","cool","term","start","term_delta")
post_cols <- c("theta","alpha","delta","prob","mortality_srm","forctoTg",
               "TgtoUSD","mortality_ozone","vsl","vsl_eta","dg")
missing_cols <- setdiff(c(fair_cols, post_cols), names(id_montecarlo))
if (length(missing_cols) > 0) {
  stop("id_montecarlo csv missing columns: ", paste(missing_cols, collapse=", "),
       ". Regenerate it with the current Generate_montecarlo.R.")
}
id_montecarlo <- id_montecarlo[, c(fair_cols, post_cols), with = FALSE]

# FAIR/scenario columns are matched against the gdx file names as strings.
for (cname in setdiff(fair_cols, "rcp")) {
  id_montecarlo[, (cname) := as.character(as.integer(round(as.numeric(get(cname)))))]
}
id_montecarlo[, rcp := str_remove(as.character(rcp), "RCP")]
# Post-processing columns stay numeric for the cost/damage computations.
for (cname in post_cols) id_montecarlo[, (cname) := as.numeric(get(cname))]
setnames(id_montecarlo, c("pulse","cool","term_delta","start","term"), c("pulse_time","cool_rate","term","geo_start","geo_end"))
id_montecarlo <- unique(id_montecarlo)

cat("Loaded", nrow(id_montecarlo), "montecarlo realizations from id_montecarlo.csv \n")

# Post-processing parameters (incl. the deterministic main-scenario overrides
# for vsl/delta/vsl_eta/theta) are sampled and fixed in Generate_montecarlo.R.

all_names <- c("gas", names(id_montecarlo))

# extract scenario names
scenario_names <- setdiff(names(sanitized_names), c("gdx","file","experiment"))
pulse_scenarios <- setdiff(scenario_names, c("cool_rate","geo_end","geo_start","term"))
base_scenarios <- setdiff(pulse_scenarios, c("gas","pulse_time"))


## all_experiments: srmpulsemaskedterm rows, scenario columns only
all_experiments <- unique(
  sanitized_names[experiment == "srmpulsemaskedterm", ..scenario_names]
)

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

if(file.exists(file.path(output_folder, name_output))) { 
  existing_output <- fread(file = file.path(output_folder, name_output))
  existing_output <- unique(existing_output[,..scenario_names])
  existing_output[, names(existing_output) := lapply(.SD, as.character)]
  if (skip_scenarios == F) {
    file.remove(file.path(output_folder, name_output))
    file.remove(file.path(output_folder, paste0("scc_",name_output)) ) 
    file.remove(file.path(output_folder, paste0("sccnosrm_",name_output)) ) }
  } else {skip_scenarios <- F}

cat("End data preprocessing...\n")

for (n_chunk in seq(1,nrow(all_experiments), by = N+1 )) { 
  
  files_loop <- c()
  
  ## current experiment row
  for (n_loop in seq(n_chunk,min(n_chunk+N,nrow(all_experiments))) ) {
  
  exp_row <- all_experiments[n_loop]
  
  if ( skip_scenarios==T) if (nrow(existing_output[exp_row, on = k_all, nomatch = 0L]) > 0 ) next 
  
  ## 1) srmpulsemaskedterm: match on all scenario columns
  f1 <- dt_srmpt[exp_row, on = k_all, mult = "first"]$gdx
  
  ## 2-4) srmpulsemasked / srmpulse / srm: drop term
  exp_no_term <- exp_row[, ..k_no_term]
  
  f2 <- dt_srmptm[exp_no_term, on = k_no_term, mult = "first"]$gdx
  f3 <- dt_srmp[exp_no_term,   on = k_no_term, mult = "first"]$gdx
  f4 <- dt_srm[exp_no_term,    on = k_no_term, mult = "first"]$gdx
  
  ## 5) pulse: drop term, cool_rate, geo_start, geo_end
  exp_pulse <- exp_row[, ..k_pulse]
  f5 <- dt_pulse[exp_pulse, on = k_pulse, mult = "first"]$gdx
  
  ## 6) base: drop term, cool_rate, geo_start, geo_end, gas, pulse_time
  exp_base <- exp_row[, ..k_base]
  f6 <- dt_base[exp_base, on = k_base, mult = "first"]$gdx
  
  files_loop_mini <- c(f1, f2, f3, f4, f5, f6)
  
  if (length(files_loop_mini) != 6 || any(is.na(files_loop_mini))) next else files_loop <- append(files_loop,files_loop_mini) 
  
  }

if (length(files_loop) == 0) next     # skip loop iteration
  
files_loop <- unique(files_loop)

cat("Loading the data for experiments", n_chunk, "to", n_chunk+N, "from a total of", nrow(all_experiments), "experiments...\n")

# create full grid to handle missing data
all_t <- seq(1,480)  # or tot_forcing$t
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
dt_sel <- dt_sel[ t >= as.numeric(pulse_time) ]
res_sel <- npv_aggregator(dt_sel, all_names)

# pulse_size & scc (reuse earlier approach but modularized)
pulse_size <- W_EMI[ ghg == gas & experiment %in% c("srm","srmpulse"), .(t, value, experiment, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[, .(pulse_size = value[ experiment == "srmpulse" ] - value[ experiment == "srm" ]), by = .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[ pulse_size != 0 ]
pulse_size[ , pulse_size := ifelse(gas == "co2", pulse_size * 1e9, pulse_size * 1e6)]
pulse_size[, t:=NULL]

# scc calculation (with SRM)
scc <- merge(
  TATM[experiment == "srm", .(t, temp_srm = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  TATM[experiment == "srmpulse", .(t, temp_srmpulse = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE
)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "srmpulse", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, concch4_pulse = value)], by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "srm", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, concch4_base = value)],  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all.x = TRUE)
replace_num_na0(scc)
scc <- merge(scc, id_montecarlo[, !c("theta","term","prob","mortality_srm","forctoTg"), with = FALSE], by = intersect(names(scc), names(id_montecarlo)), allow.cartesian = TRUE)
scc <- scc[ t >= as.numeric(pulse_time) ]
compute_gwpt(scc)
scc[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp) ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
scc[, dam := gwpt * ( alpha * (pmax(0,temp_srmpulse)^2 - pmax(0,temp_srm)^2) )]
scc_agg <- scc[, .(damnpv = sum( dam / (1 + delta)^(t - as.numeric(pulse_time)),na.rm = TRUE),
                    ozpnpv = sum( tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE)), 
                by = c("gas",intersect(names(scc), names(id_montecarlo))) ]
scc_agg <- merge(scc_agg, pulse_size, by = intersect(names(scc_agg), names(pulse_size)), all.x = TRUE)
scc_agg[, scc := (damnpv + ozpnpv) / pulse_size]
scc_agg <- unique(scc_agg)

# scc calculation (w/o SRM)
scc <- merge(
  TATM[experiment == "base", .(t, temp_base = value, rcp, ecs, tcr)],
  TATM[experiment == "pulse", .(t, temp_pulse = value, gas, rcp, ecs, tcr, pulse_time)],
  by = c("t","rcp","ecs","tcr"), all = FALSE
)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "pulse", .(t, gas, rcp, ecs, tcr, pulse_time, concch4_pulse = value)], by = c("t","gas","rcp","ecs","tcr","pulse_time"), all = FALSE)
scc <- merge(scc, CONC[ ghg == "ch4" & experiment == "base", .(t, rcp, ecs, tcr, concch4_base = value)],  by = c("t","rcp","ecs","tcr"), all.x = TRUE)
replace_num_na0(scc)
scc <- merge(scc, id_montecarlo[, c("ecs","tcr","rcp","pulse_time","alpha","delta","mortality_ozone","vsl","vsl_eta","dg"), with = FALSE], by = intersect(names(scc), names(id_montecarlo)), allow.cartesian = TRUE)
scc <- scc[ t >= as.numeric(pulse_time) ]
compute_gwpt(scc)
scc[, tropoz_pollution := vsl * globaltouspc * (gwpt/gwp) ^ (vsl_eta) * mortality_ozone * (concch4_pulse - concch4_base)]
scc[, dam := gwpt * ( alpha * (pmax(0,temp_pulse)^2 - pmax(0,temp_base)^2) )]
scc_agg2 <- scc[, .(damnpv = sum( dam / (1 + delta)^(t - as.numeric(pulse_time)),na.rm = TRUE),
                    ozpnpv = sum( tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE)), 
                by = c("gas",intersect(names(scc), names(id_montecarlo))) ]
scc_agg2 <- merge(scc_agg2, pulse_size[, c("gas","ecs","tcr","rcp","pulse_time","pulse_size"), with = FALSE], by = intersect(names(scc_agg2), names(pulse_size)), all.x = TRUE)
scc_agg2[, scc_nosrm := (damnpv + ozpnpv) / pulse_size]
scc_agg2 <- unique(scc_agg2)

# combine results: a single experiment's NPC, no probability weighting
combined <- res_sel[, c(all_names, "source", "cost"), with = FALSE]
replace_num_na0(combined)
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

# ----------------------------
# Save final output (combined table)
# ----------------------------
fwrite(combined_wide, file = file.path(output_folder, name_output), append = TRUE)
fwrite(scc_agg, file = file.path(output_folder, paste0("scc_",str_remove(name_output,"npc_") )), append = TRUE)
fwrite(scc_agg2, file = file.path(output_folder, paste0("sccnosrm_",str_remove(name_output,"npc_"))), append = TRUE)

}

cat("Saved output to:", file.path(output_folder, name_output), "\n")

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("This code run for exactly",time.taken, "seconds")

if(run_hpc==T) {
  sink(type = "message")
  sink()
  close(logfile) }
