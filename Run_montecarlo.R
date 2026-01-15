# this scripts runs the montecarlo scenarios in batches 
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Run_montecarlo.R [-o <res>] [-l <load_data>] [-q <which_queue>] [-m <max_scenarios>] [-x <rerun_problem>] [-w <overwrite_data>] [-p <run_parallel>] [-s <start_job>] [-e <end_job>] [--hpc <run_hpc>] 

Options:
-o <res>                     Path of the results (default: Results)
-l <load_data>               Path(s) and name of .csv with existing realizations 
-p <run_parallel>            T/F if to run in parallel or in series (for Juno)
-s <start_job>               Number of line to start calling the scenarios 
-e <end_job>                 End of line to start calling the scenarios
-q <which_queue>             Select queue to send the jobs to (for parallel solving)
-m <max_scenarios>           Maximum scenarios to send to solve
-x <rerun_problem>           T to avoid rerunning problematic scenarios (default)
--hpc <run_hpc>              T/F if to run on Juno or local (windows)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Run_montecarlo')

# logical
run_parallel = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )
run_hpc = ifelse(is.null(opts[["hpc"]]), F, as.logical(opts["hpc"]) )
rerun_problem = ifelse(is.null(opts[["x"]]),T, as.logical(opts["x"]) )

# numeric
max_scenarios = ifelse(is.null(opts[["m"]]), 100, as.numeric(opts["m"]) )
start_job = ifelse(is.null(opts[["s"]]), 1, as.numeric(opts["s"]) )
end_job = ifelse(is.null(opts[["e"]]), 100000, as.numeric(opts["e"]) )

# strings
res = ifelse(is.null(opts[["o"]]), "Results_montecarlo", as.character(opts["o"]) )
input = ifelse(is.null(opts[["l"]]), "Montecarlo", as.character(opts["l"]) )
which_queue = ifelse(is.null(opts[["q"]]), "p_short", as.character(opts["q"]) )

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res)
}

if (file.exists(paste0(res,"/problematic_scenarios.csv"))) {
  problematic_data <- data.frame(read.csv(paste0(res,"/problematic_scenarios.csv"))) %>% select(-X)
}

if (!file.exists(paste0(input,"/id_montecarlo.csv"))) stop("Please generate new data if no pre-existing are available")

data <- as.data.frame(read.csv(paste0(input,"/id_montecarlo.csv"))) %>% select(-X)

if (rerun_problem==F) data <- data %>% anti_join(problematic_data %>% select(-error,-gas) %>% unique())
cat("Launching jobs... \n")

filelist <- list.files(path=paste0(res,"/"),pattern="*.gdx")

scenarios_launched <- 0
gases <- c("ch4","co2")

for (gas in gases) {
  
for (i in seq(start_job,min(end_job,nrow(data))) ) {
bsub <- paste("bsub", "-q",which_queue, "-n 1",
              "-P 0638", paste0("-J scenariosrmpulse", i,"_gas",gas), "-K -M 64G")

gams <- paste0("gams FAIR.gms --experiment=srm",
               " --gas=",gas,
              " --ecs=",data[i,]$ecs,
              " --tcr=",data[i,]$tcr,
              " --rcp=",data[i,]$rcp,
              " --pulse_time=",data[i,]$pulse,
              " --rate_of_cooling=",data[i,]$cool,
              " --start_rampdown=",data[i,]$term-100,
              " --end_rampdown=",data[i,]$term,
              " --start_rampup=",data[i,]$start,
              " --end_rampup=",data[i,]$start+100,
              " --termination_time=",data[i,]$term_delta,
              " --results_folder=",res)

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

if (!any(str_detect(str_remove(filelist,".gdx"),results_name)) ) {
  
  if (run_parallel==T) { 
    bsub <- str_remove(bsub, "-K ") 
    scenarios_launched <- scenarios_launched + 1 }
  
  command <- paste(bsub, gams)
  
  if (run_hpc==F) {command <- gams}
  
  if (run_parallel==T) {cat("Attention! This scenarios was not solved before. \n Launching scenario",i,"with gas",gas,"...\n")} else {cat("Attention! This scenarios was not solved before. \n Solving scenario",i,"with gas",gas,"...\n")}
    
  ret <- system(command = command, intern = TRUE)
  
    if(!is.null(attr(ret,"status")) )  {
    
        if(attr(ret,"status")!=112) {
        
          cat("Careful! scenario",i,"with gas",gas,"couldn't solve...\n")
          problematic_data <- problematic_data %>% 
          bind_rows(c(data[i,],error=attr(ret,"status"),gas=gas)) } 
      
        }
  }

if(run_parallel==T & scenarios_launched >= max_scenarios) {
  if(nrow(problematic_data)!=0) {
    cat("Careful! Some scenarios couldn't solve. Information is stored in problematic_scenarios.csv")
    write.csv(problematic_data,file=paste0(res,"/problematic_scenarios.csv")) }
  stop("The maximum number of scenarios (",max_scenarios,") has been launched for this batch") }
}

}

if(nrow(problematic_data)!=0) {
  cat("Careful! Some scenarios couldn't solve. Information is stored in problematic_scenarios.csv")
  write.csv(problematic_data,file=paste0(res,"/problematic_scenarios.csv")) }


