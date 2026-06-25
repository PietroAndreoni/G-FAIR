# this scripts runs the montecarlo scenarios in batches 
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Run_montecarlo.R [--res <res>] [-i <load_data>] [-q <which_queue>] [-m <max_scenarios>] [-x <skip_problematic>] [-w <overwrite_data>] [-p <run_parallel>] [-s <start_job>] [-e <end_job>] [--hpc <run_hpc>] [--base <main_scenario>]

Options:
--res <res>                  Path of the results (default: Results)
-i <load_data>               Path(s) and name of .csv with existing realizations 
-p <run_parallel>            T/F if to run in parallel or in series (for Juno)
-s <start_job>               Number of line to start calling the scenarios 
-e <end_job>                 End of line to start calling the scenarios
-q <which_queue>             Select queue to send the jobs to (for parallel solving)
-m <max_scenarios>           Maximum scenarios to send to solve
-x <skip_problematic>        T to skip known problematic scenarios (default)
--hpc <run_hpc>              T/F if to run on Juno or local (windows)
--base <main_scenario>       T/F to run the base scenario      
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Run_montecarlo')

# Locate the Paper_SAI root (holds all_parameters.R), load the control file and
# the shared validators / MC_* constants in Utilities/. Works from any wd.
.sp <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- if (length(.sp) == 1) dirname(.sp) else getwd()
while (!file.exists(file.path(.root, "all_parameters.R")) && dirname(.root) != .root)
  .root <- dirname(.root)
if (!file.exists(file.path(.root, "all_parameters.R")))
  stop("Cannot locate all_parameters.R (Paper_SAI root).")
source(file.path(.root, "all_parameters.R"))
source(file.path(PAPER_ROOT, "Utilities", "montecarlo_utils.R"))

# logical
run_parallel = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )
run_hpc = ifelse(is.null(opts[["hpc"]]), F, as.logical(opts["hpc"]) )
skip_problematic = ifelse(is.null(opts[["x"]]),T, as.logical(opts["x"]) )
main_scenario = ifelse(is.null(opts[["base"]]), F, as.logical(opts["base"]) )

# numeric
max_scenarios = ifelse(is.null(opts[["m"]]), 100, as.numeric(opts["m"]) )
start_job = ifelse(is.null(opts[["s"]]), 1, as.numeric(opts["s"]) )
end_job = ifelse(is.null(opts[["e"]]), 100000, as.numeric(opts["e"]) )

# strings
# Folder names are user-specifiable; their locations are fixed. The GAMS .gdx
# results go to Results/ (under Paper_SAI); id_montecarlo.csv is read from the
# working folder under Sampling/.
res_name = ifelse(is.null(opts[["res"]]), RUN_RESULTS_DEFAULT, as.character(opts["res"]) )
res = file.path(RUN_RESULTS_PARENT, res_name)
input_name = ifelse(is.null(opts[["i"]]), MC_WORK_DEFAULT, as.character(opts["i"]) )
input = file.path(MC_WORK_PARENT, input_name)
which_queue = ifelse(is.null(opts[["q"]]), HPC_QUEUE, as.character(opts["q"]) )

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res, recursive = TRUE)
}

problematic_existing <- data.frame()
if (file.exists(paste0(res,"/problematic_scenarios.csv"))) {
  problematic_existing <- mc_strip_csv_index(data.frame(read.csv(paste0(res,"/problematic_scenarios.csv"), stringsAsFactors = FALSE)))
}

if (!file.exists(paste0(input,"/id_montecarlo.csv"))) stop("Please generate new data if no pre-existing are available")

data <- mc_strip_csv_index(as.data.frame(read.csv(paste0(input,"/id_montecarlo.csv"), stringsAsFactors = FALSE)))
validate_id_montecarlo(data, paste0(input, "/id_montecarlo.csv"))

if (skip_problematic==T && nrow(problematic_existing) > 0) {
  problematic_keys <- problematic_existing[setdiff(names(problematic_existing), c("error","gas"))] %>% unique()
  join_cols <- intersect(names(data), names(problematic_keys))
  if (length(join_cols) > 0) data <- data %>% anti_join(problematic_keys, by = join_cols)
}
cat("Launching jobs... \n")

filelist <- list.files(path=paste0(res,"/"),pattern="*.gdx")
known_results <- tools::file_path_sans_ext(filelist)

scenarios_launched <- 0
problematic_data <- data.frame()
if(main_scenario==T) gases <- GASES_MAIN else gases <- GASES_FULL

last_job <- min(end_job,nrow(data))
if (start_job > last_job) stop("No scenarios to run in the requested start/end range after filtering.")

for (i in seq.int(start_job,last_job) ) {
  
for (gas in gases) {
    
bsub <- paste("bsub", "-q",which_queue, "-n 1",
              paste0("-P ", HPC_PROJECT), paste0("-J scenariosrmpulse", i,"_gas",gas), paste0("-K -M ", HPC_MEMORY))

gams <- paste0("gams FAIR.gms --experiment=srm",
               " --gas=",gas,
              " --ecs=",data[i,]$ecs,
              " --tcr=",data[i,]$tcr,
              " --rcp=",data[i,]$rcp,
              " --pulse_time=",data[i,]$pulse,
              " --rate_of_cooling=",data[i,]$cool,
              " --start_rampdown=",data[i,]$term-RAMP_DOWN_YEARS,
              " --end_rampdown=",data[i,]$term,
              " --start_rampup=",data[i,]$start,
              " --end_rampup=",data[i,]$start+RAMP_UP_YEARS,
              " --termination_time=",data[i,]$term_delta,
              ' --results_folder="',res,'"')  # quoted: the resolved path may contain spaces

results_name <-  paste0(data[i,]$rcp,
                        "_EXPsrmpulsemaskedterm_TER",data[i,]$term_delta,
                        "_GAS",gas,
                        "_ECS",data[i,]$ecs,
                        "_TCR",data[i,]$tcr,
                        "_PT",data[i,]$pulse,
                        "_RC",data[i,]$cool,
                        "_EC",data[i,]$term,
                        "_BC",data[i,]$start,
                        "_IChistorical_run")

cat("Checking scenario",i,"with gas",gas,"...\n")

if (!results_name %in% known_results ) {
  
  if (run_parallel==T) { 
    bsub <- str_remove(bsub, "-K ") 
    scenarios_launched <- scenarios_launched + 1 }
  
  command <- paste(bsub, gams)
  
  if (run_hpc==F) {command <- gams}
  
  if (run_parallel==T) {cat("Attention! This scenarios was not solved before. \n Launching scenario",i,"with gas",gas,"...\n")} else {cat("Attention! This scenarios was not solved before. \n Solving scenario",i,"with gas",gas,"...\n")}
    
  ret <- system(command = command, intern = TRUE)
  known_results <- c(known_results, results_name)
  
    if(!is.null(attr(ret,"status")) )  {
    
        if(attr(ret,"status")!=SOLVE_OK_STATUS) {
        
          cat("Careful! scenario",i,"with gas",gas,"couldn't solve...\n")
          problematic_data <- problematic_data %>% 
          bind_rows(data.frame(data[i,], error=attr(ret,"status"), gas=gas, stringsAsFactors = FALSE)) }
      
        }
  }

if(run_parallel==T & scenarios_launched >= max_scenarios) {
  if(nrow(problematic_data)!=0) {
    cat("Careful! Some scenarios couldn't solve. Information is stored in problematic_scenarios.csv")
    write.csv(problematic_data,file=paste0(res,"/problematic_scenarios.csv"), row.names = FALSE) }
  stop("The maximum number of scenarios (",max_scenarios,") has been launched for this batch") }
}

}

if(nrow(problematic_data)!=0) {
  cat("Careful! Some scenarios couldn't solve. Information is stored in problematic_scenarios.csv")
  write.csv(problematic_data,file=paste0(res,"/problematic_scenarios.csv"), row.names = FALSE) }


