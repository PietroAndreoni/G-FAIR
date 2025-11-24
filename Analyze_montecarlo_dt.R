library(data.table)
library(stringr)
library(gdxtools)
library(EnvStats) # for rtri

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
  return(DT)
}


# fast replace numeric NA with 0 in-place
replace_num_na0 <- function(DT) {
  numcols <- names(DT)[vapply(DT, is.numeric, logical(1))]
  for (cname in numcols) set(DT, which(is.na(DT[[cname]])), cname, 0)
  invisible(NULL)
}

# vectorized gwpt computation
compute_gwpt <- function(DT, gwp_col = "gwp", dg_col = "dg", t_col = "t") {
  # gwpt = pmin(gwp * (1+dg)^(280-1), gwp * (1+dg)^(t-1))
  DT[, gwpt := pmin(get(gwp_col) * (1 + get(dg_col))^(280-1),
                    get(gwp_col) * (1 + get(dg_col))^(get(t_col)-1))]
}

npv_aggregator <- function(DT, keep_names = all_names) {
  compute_gwpt(DT)
  DT[, dam := gwp * (alpha * (temp^2 - temp_srm^2))]
  DT[, direct_cost := srm_masking * forctoUSD]
  DT[, srm_pollution := vsl * srm_masking * forctoTg * mortality_srm]
  DT[, tropoz_pollution := vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc]
  DT[, imperfect_masking := 2 * gwp * alpha * temp_base * (1 - cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1 + srm/forc)]
  # aggregator
  res <- DT[, .(
    dirnpv = sum(direct_cost / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    srmpnpv = sum(srm_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    ozpnpv = sum(tropoz_pollution / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    masknpv = sum(imperfect_masking / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE),
    damnpv = sum(dam / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE)
  ), by = keep_names]
  res[, costnpv := dirnpv + srmpnpv + ozpnpv + masknpv + damnpv]
  return(res)
}

drawln <- function(median,std,n=1,plot=F) {
  location <- log(median)#log(m^2 / sqrt(s^2 + m^2))
  y <- (1 + sqrt(1 + 4 * (std^2 / median^2))) / 2
  shape <- sqrt(log(y))
  if (plot==T) {
    draws <- rlnorm(n = 10000, location, shape)
    print(paste("mean: ",round(mean(draws),5) ) )
    print(paste("median: ",round(median(draws), 5) ) )
    print(paste("std: ",round(sd(draws), 5) ) )
    plot(density(draws[draws > 0 & draws < 10*std])) }
  
  return(rlnorm(n, location, shape))
}


prepare_join_table <- function(filter_experiment, include_term = TRUE) {
  
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
  
  ozone_exp <- copy(FORC)[ghg == "o3trop" & experiment == filter_experiment]
  ozone_exp[, c("file", "experiment","ghg") := NULL]
  setnames(ozone_exp, "value", "ozone_pulse", skip_absent = TRUE)
  
  ozone_srm <- copy(FORC)[ghg == "o3trop" & experiment == "srm"]
  ozone_srm[, c("file", "experiment","term","ghg") := NULL]
  setnames(ozone_srm, "value", "ozone_base", skip_absent = TRUE)
  
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
  base <- ozone_exp[base, on = c("t", scenario_names)]
  base <- ozone_srm[base, on = setdiff(c("t", scenario_names),"term") ]
  base <- temp_srm[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_srmpulse[base, on = setdiff(c("t", scenario_names),"term")]
  base <- temp_base[base, on = c("t", base_scenarios)]
  base <- temp_ghgpulse[base, on = c("t",pulse_scenarios), nomatch = NULL]
  
  replace_num_na0(base)
  base[, term := NULL] # remove term (readded later)
  base <- merge(base, id_montecarlo, by = setdiff(scenario_names,c("term","gas")))
  return(base)
}

# Usage
'Launch script to analyze montecarlo scenarios (produces a csv file in the same folder)

Usage:
  Analyze_montecarlo.R [-i <res>] [-o <output_folder>] [--hpc <run_hpc>] [-p <plot_results>] [--seed <seed>] 

Options:
-i <res>              Path where the results are (default: Results_montecarlo). For multiple folders separate with -
-o <output_folder>   Where to save results
--hpc <run_hpc>          T/F if running from Juno (T) or local (F) 
-p <plot_results>     T/F to plot data and save data analysis
--seed <seed>         seed number (for reproducibility)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Montecarlo')

res <- ifelse(is.null(opts[["i"]]), "Results_montecarlo-Results_montecarlo_2-Results_montecarlo_3", as.character(opts["i"]) )
res <- str_split(res,"-")[[1]]

output_folder <- ifelse(is.null(opts[["o"]]), "Results_output", as.character(opts["o"]) )

# Make sure the output folder exists (create it if not)
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

if (any(!dir.exists(res)) ) stop("some of the folder specified do not exsist")

run_hpc = ifelse(is.null(opts[["hpc"]]), T, as.logical(opts["h"]) )
plot_results = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )
seed = ifelse(is.null(opts[["seed"]]), 123, as.integer(opts["seed"]) )

if(run_hpc==F) {igdx()} else {igdx("/work/cmcc/pa12520/gams40.4_linux_x64_64_sfx")}

## climate damage function parameters
gwp <- 105*1e12 # initial world gdp
forctoTg <- 1/0.2 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
TgtoUSD <- 2250*10^6 # from https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7/pdf  
forctoUSD <- forctoTg * TgtoUSD # US$/(W/m^2)
ozone_rftoconc <- 50 / 0.263  # 50 ppb as https://www.sciencedirect.com/science/article/pii/S2542519622002601?ref=pdf_download&fr=RR-2&rr=9a07c2002fcb708b

# ----------------------------
# Files discovery
# ----------------------------
filelist <- unlist(lapply(res, function(folder) list.files(folder, pattern = ".gdx$", full.names = TRUE)))
filelist <- filelist[str_detect(filelist, "EXP")]


cat("Sanitizing the data...\n")
all_scenarios <- extract_names_dt(filelist)
# deduplicate by file
all_scenarios <- unique(all_scenarios, by = "file")
cat("Removing", length(filelist) - nrow(all_scenarios), "duplicate scenarios \n")

# keep only experiments that have required companion experiments
check1 <- copy(all_scenarios)
# filter out base and pulse
check1 <- check1[ !experiment %in% c("base","pulse") ]
setkeyv(check1, c("gas","rcp","cool_rate","pulse_time","geo_start","geo_end","ecs","tcr"))
check1 <- check1[, .(tot = .N, diff = any(duplicated(experiment))), by = key(check1)]
check1 <- check1[ tot == 4 & diff == FALSE ]
check1 <- check1[, c("tot","diff") := NULL ]


# Now build check2: keep rows from all_scenarios that match check1 and have pulse and base
if (nrow(check1) == 0) stop("No complete scenario groups found")

# Join pulses and base
pulses <- all_scenarios[ experiment == "pulse", .(gas, rcp, ecs, tcr, pulse_time)]
pulses <- unique(pulses, by =  c("gas","rcp","ecs","tcr","pulse_time"))
base <- all_scenarios[ experiment == "base", .(rcp, ecs, tcr)]
base <- unique(base, by =  c("rcp","ecs","tcr"))
check2 <- merge(check1, pulses, by = c("gas","rcp","ecs","tcr","pulse_time"), all = FALSE)
check2 <- merge(check2, base, by = c("rcp","ecs","tcr"), all = FALSE)

# sanitized_names: include matching rows from all_scenarios
sanitized_names <- rbindlist( list(merge(all_scenarios, check2, by = c("gas","rcp","cool_rate","pulse_time","geo_start","geo_end","ecs","tcr"), all = FALSE),
                    merge(all_scenarios[experiment == "pulse"], pulses, by = c("gas","rcp","ecs","tcr","pulse_time"), all = FALSE),
                    merge(all_scenarios[experiment == "base"], base, by = c("rcp","ecs","tcr"), all = FALSE)),
                    use.names=TRUE)

cat(sprintf("Careful! %d scenarios removed.\n", nrow(all_scenarios) - nrow(sanitized_names)))

# add pulse and base unique rows
files <- sanitized_names$gdx

# extract scenario names
scenario_names <- setdiff(names(sanitized_names), c("gdx","file","experiment"))
base_scenarios <- names(base)
pulse_scenarios <- names(pulses)

# create full grid to handle missing data
all_t <- seq(1,480)  # or tot_forcing$t
all_scenarios <- unique(sanitized_names[, ..scenario_names])  # all scenario combinations

# Create the full grid
full_grid <- CJ(t = all_t, 
                scenario = 1:nrow(all_scenarios), unique = TRUE)

# Merge scenario columns back
full_grid <- cbind(full_grid[, -"scenario", with = FALSE], all_scenarios[full_grid$scenario])




cat("Loading the data...\n")

# Batch extract using gdxtools, then join sanitized_names and sanitize
TATM <- setDT(gdxtools::batch_extract("TATM", files = files)$TATM)
TATM <- merge(TATM, sanitized_names, by.x = "gdx", by.y = "gdx", all = FALSE)
TATM <- sanitize_dt(TATM)


W_EMI <- setDT(gdxtools::batch_extract("W_EMI", files = files)$W_EMI)
W_EMI <- merge(W_EMI, sanitized_names, by = "gdx", all = FALSE)
W_EMI <- sanitize_dt(W_EMI)


FORC <- setDT(gdxtools::batch_extract("FORCING", files = files)$FORCING)
FORC <- merge(FORC, sanitized_names, by = "gdx", all = FALSE)
FORC <- sanitize_dt(FORC)


SRM <- setDT(gdxtools::batch_extract("SRM", files = files)$SRM)
SRM <- merge(SRM, sanitized_names, by = "gdx", all = FALSE)
SRM <- sanitize_dt(SRM)


background_srm <- setDT(gdxtools::batch_extract("forcing_srm", files = files)$forcing_srm)
background_srm <- merge(background_srm, sanitized_names, by = "gdx", all = FALSE)
background_srm <- sanitize_dt(background_srm)


# total forcing aggregation (equivalent to dplyr group_by + summarise)
tot_forcing <- FORC[, .(value = sum(value)), by = c("t", "file",scenario_names, "experiment")]


# Load id_montecarlo from folders
id_list <- lapply(res, function(folder) {
  read.csv(file.path(folder, "id_montecarlo.csv"), stringsAsFactors = FALSE)
})
id_montecarlo <- rbindlist(id_list, fill = TRUE)
setDT(id_montecarlo)
# select, mutate types & rename as original script
keep_cols <- c("ecs","tcr","rcp","pulse","cool","term","start","term_delta")
if (!all(keep_cols %in% names(id_montecarlo))) stop("id_montecarlo csv missing columns")
id_montecarlo <- id_montecarlo[, .(ecs, tcr, rcp, pulse, cool, term, start, term_delta)]
# mutate/rename
id_montecarlo[, term := as.integer(term)]
# convert integer columns to character as in original
int_cols <- names(id_montecarlo)[vapply(id_montecarlo, is.integer, logical(1))]
for (cname in int_cols) id_montecarlo[, (cname) := as.character(get(cname))]
id_montecarlo[, rcp := str_remove(rcp, "RCP")]
setnames(id_montecarlo, c("pulse","cool","term_delta","start","term"), c("pulse_time","cool_rate","term","geo_start","geo_end"))
id_montecarlo <- unique(id_montecarlo)

cat("Proucing montecarlo realizations...")
n_scenarios <- nrow(id_montecarlo)
set.seed(seed)
# add uncertainties
id_montecarlo[, theta := round(drawln(15,7,n_scenarios),0)]
id_montecarlo[, alpha := round(drawln(0.00575,0.00575*150/(230-100),n_scenarios),5)]
id_montecarlo[, delta := round(runif(n_scenarios,0.001,0.07),3)]
id_montecarlo[, prob := round(runif(n_scenarios,0,1),2)]
id_montecarlo[, mortality_srm := round(drawln(7400,(16000-2300)/1.96,n_scenarios),0)]
id_montecarlo[, mortality_ozone := round(EnvStats::rtri(n_scenarios, min = 100, max = 107, mode = 104),1)]
id_montecarlo[, vsl := round(runif(n_scenarios,1,10),0) * 1e6]
id_montecarlo[, dg := round(rnorm(n_scenarios,mean=0.015,sd=0.005),4)]

all_names <- c("gas", names(id_montecarlo))


cat("Analyzing the data...")
# build and aggregate for pre, post_noterm, post_term using above functions
dt_noterm <- prepare_join_table("srmpulsemasked")
dt_pre <- dt_noterm[ t >= as.numeric(pulse_time) & t <= as.numeric(term) ]
res_pre <- npv_aggregator(dt_pre, all_names)

dt_post_noterm <- dt_noterm[ t >= as.numeric(pulse_time) & t > as.numeric(term) ]
res_post_noterm <- npv_aggregator(dt_post_noterm, all_names)

dt_term <- prepare_join_table("srmpulsemaskedterm")
dt_post_term <- dt_term[ t >= as.numeric(pulse_time) & t > as.numeric(term) ]
res_post_term <- npv_aggregator(dt_post_term, all_names)

# pulse_size & scc (reuse earlier approach but modularized)
pulse_size <- W_EMI[ ghg == gas & experiment %in% c("srm","srmpulse"), .(t, value, experiment, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[, .(pulse_size = value[ experiment == "srmpulse" ] - value[ experiment == "srm" ]), by = .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)]
pulse_size <- pulse_size[ pulse_size != 0 ]
pulse_size[ , pulse_size := ifelse(gas == "co2", pulse_size * 1e9, pulse_size * 1e6)]

# scc calculation (modular)
scc <- merge(
  TATM[experiment == "srm", .(t, temp_srm = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  TATM[experiment == "srmpulse", .(t, temp_srmpulse = value, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end)],
  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE
)
scc <- merge(scc, FORC[ ghg == "o3trop" & experiment == "srmpulse", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, ozone_pulse = value)], by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all = FALSE)
scc <- merge(scc, FORC[ ghg == "o3trop" & experiment == "srm", .(t, gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, ozone_base = value)],  by = c("t","gas","rcp","ecs","tcr","cool_rate","pulse_time","geo_start","geo_end"), all.x = TRUE)
replace_num_na0(scc)
scc <- merge(scc, id_montecarlo[, !c("theta","term","prob","mortality_srm"), with = FALSE], by = intersect(names(scc), names(id_montecarlo)), allow.cartesian = TRUE)
scc <- scc[ t >= as.numeric(pulse_time) ]
compute_gwpt(scc)
scc[, tropoz_pollution := vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc]
scc[, dam := gwp * ( alpha * ((temp_srmpulse)^2 - (temp_srm)^2) )]
scc_agg <- scc[, .(damnpv = sum( (dam + tropoz_pollution) / (1 + delta)^(t - as.numeric(pulse_time)), na.rm = TRUE)), by = setdiff(all_names, c("gas","theta","term","prob","mortality_srm"))]
scc_agg <- merge(scc_agg, pulse_size, by = intersect(names(scc_agg), names(pulse_size)), all.x = TRUE)
scc_agg[, scc := damnpv / pulse_size]

# combine results
combined <- merge(res_pre[, c(all_names, "costnpv"), with = FALSE], res_post_noterm[, c(all_names, "costnpv"), with = FALSE], by = all_names, all = TRUE, suffixes = c("_pre","_postnoterm"))
combined <- merge(combined, res_post_term[, c(all_names, "costnpv"), with = FALSE], by = all_names, all = TRUE)
setnames(combined, c("costnpv_pre","costnpv_postnoterm","costnpv"), c("costpre","costpostnoterm","costpostterm"))
replace_num_na0(combined)
combined[, costnpv := costpre + as.numeric(prob) * costpostterm + (1 - as.numeric(prob)) * costpostnoterm]
combined[, c("costpre","costpostterm","costpostnoterm") := NULL]
combined <- merge(combined, pulse_size, by = intersect(names(combined), names(pulse_size)), all.x = TRUE)
combined <- merge(combined, scc_agg, by = intersect(names(combined), names(scc_agg)), all.x = TRUE)
combined[, npc_srm := costnpv / pulse_size]

# ----------------------------
# Save final output (combined table)
# ----------------------------
fwrite(combined, file = file.path(output_folder, "output_analysis_dt.csv"))
cat("Saved output to:", file.path(output_folder, "output_analysis_dt.csv"), "\n")

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("This code run for exactly",time.taken, "seconds")