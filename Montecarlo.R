# this scripts generates the realizations for the montecarlo parameters
# generates a script with all the scenarios to run (excluding those with results)
# and runs the scenarios in batches 
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Montecarlo.R [-g <generate_data>] [-o <res>] [-q <which_queue>] [-n <n_scenarios>] [-m <max_scenarios>] [-x <rerun_problem>] [-w <overwrite_data>] [-p <run_parallel>] [-r <run_scenarios>] [-h <run_hpc>] [-s <start_job>] [-e <end_job>] [--seed <seed>] 

Options:
-o <res>              Path where the results are (default: Results)
-n <n_scenarios>      Number of new scenarios to generate
-g <generate_data>    T/F to generate new realizations
-w <overwrite_data>   T/F to overwrite old realizations
-r <run_scenarios>    T/F to run or just generate scenarios (in .sh script)
-p <run_parallel>     T/F if to run in parallel or in series (for Juno)
-h <run_hpc>          T/F if to run on Juno or local (windows)
-s <start_job>        Number of line to start calling the scenarios 
-e <end_job>          End of line to start calling the scenarios
-q <which_queue>      Select queue to send the jobs to (for parallel solving)
-m <max_scenarios>    Maximum scenarios to send to solve
-x <rerun_problem>    T to avoid rerunning problematic scenarios (default)
--seed <seed>   seed number (for reproducibility)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Montecarlo')

drawln <- function(median,std,plot=F) {
  location <- log(median)#log(m^2 / sqrt(s^2 + m^2))
  y <- (1 + sqrt(1 + 4 * (std^2 / median^2))) / 2
  shape <- sqrt(log(y))
  if (plot==T) {
    draws <- rlnorm(n = 10000, location, shape)
    print(paste("mean: ",round(mean(draws),5) ) )
    print(paste("median: ",round(median(draws), 5) ) )
    print(paste("std: ",round(sd(draws), 5) ) )
    plot(density(draws[draws > 0 & draws < 10*std])) }
  
  return(rlnorm(n = 1, location, shape))
}

# logical
generate_data = ifelse(is.null(opts[["g"]]), T, as.logical(opts["g"]) )
overwrite_data = ifelse(is.null(opts[["w"]]), F, as.logical(opts["w"]) )
run_scenarios = ifelse(is.null(opts[["r"]]), T, as.logical(opts["r"]) )
run_parallel = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )
run_hpc = ifelse(is.null(opts[["h"]]), T, as.logical(opts["h"]) )
rerun_problem = ifelse(is.null(opts[["x"]]),T, as.logical(opts["x"]) )

# numeric
n_scenarios = ifelse(is.null(opts[["n"]]), 10, as.numeric(opts["n"]) )
seed = ifelse(is.null(opts[["seed"]]), 123, as.integer(opts["seed"]) )
max_scenarios = ifelse(is.null(opts[["m"]]), 100, as.numeric(opts["m"]) )
start_job = ifelse(is.null(opts[["s"]]), 1, as.numeric(opts["s"]) )
end_job = ifelse(is.null(opts[["e"]]), 100000, as.numeric(opts["e"]) )

# strings
res = ifelse(is.null(opts[["o"]]), "Results_montecarlo", as.character(opts["o"]) )
which_queue = ifelse(is.null(opts[["q"]]), "p_short", as.character(opts["q"]) )

# Define the path to the .ssh file
sh_file <- paste0(res,"/montecarlo.sh" ) 
scenarios_launched <- 0
problematic_data <- data.frame()

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res)
}

# Make sure the file exists (create it if not)
if (!file.exists(sh_file)) {
  file.create(sh_file)
} else {
  file.remove(sh_file)
  file.create(sh_file)}

if (file.exists(paste0(res,"/problematic_scenarios.csv"))) {
  problematic_data <- data.frame(read.csv(paste0(res,"/problematic_scenarios.csv"))) %>% select(-X)
}
cat("Generating data... \n")

if(generate_data==F & overwrite_data==T) stop("Invalid options! You cannot overwrite old data without generating new ones")

if (generate_data==F & !file.exists(paste0(res,"/id_montecarlo.csv"))) stop("Please generate new data if no pre-existing are available")

set.seed(seed)

if (generate_data==T) {

if (overwrite_data==T | !file.exists(paste0(res,"/id_montecarlo.csv"))) {
  data <- data.frame()
} else {
  data <-  as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv")) ) %>% select(-X)
  write.csv(data,file=paste0(res,"/id_montecarlo_copy.csv")) # story a copy of previous version database
  }

max_id <- nrow(data)

for (i in seq(1,n_scenarios,by=1)) {
# climate equilibrium sensitivity
  ecs <- round(drawln(3.14,2)*10,0)

# tcr
  tcr <- round(drawln(1.85,0.8)*10,0)

# rcp
  rcp <- sample(c("RCP3PD","RCP45","RCP6","RCP85"),1)

# time of the pulse
  pt <- sample(c(2,10,20,30,40,50,60,70,80,90,100),1)

# cooling rate 
  cool <- round(runif(1,min=0,max=40),0)

# termination year 
  term <- sample(seq(2300,2600,by=100),1)

# start year 
  start <- sample(seq(2025,2100,by=25),1)

#### parameters post-solve

# theta
  theta <- round(drawln(15,7),0)

# alpha
  alpha <- round(drawln(0.00575,0.00575*150/(230-100)),5)

# delta 
  delta <- sample(seq(0,0.1,by=0.01),1)

# prob 
  prob <- sample(seq(0,1,by=0.05),1)

# time of termination event 
  time_term <- sample(seq(10,500,by=10),1)
  
  
  data <- data %>% 
    bind_rows(data.frame(ID=i+max_id,
                       ecs=ecs,
                       tcr=tcr,
                       rcp=rcp,
                       pulse=pt,
                       cool=cool,
                       term=term,
                       start=start,
                       theta=theta,
                       alpha=alpha,
                       delta=delta,
                       term_delta=time_term,
                       prob=prob) ) }

write.csv(data,file=paste0(res,"/id_montecarlo.csv")) } else {
  
data <- as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv"))) %>% select(-X)

}
 
# select unique scenarios 
data_srmpulse <- data %>% 
  select(ecs, tcr, rcp, pulse, cool, term, start, term_delta) %>% unique() 
if (rerun_problem==F) data_srmpulse <- data_srmpulse %>% anti_join(problematic_data %>% select(-error,-gas) %>% unique())
cat("Launching jobs... \n")

filelist <- list.files(path=paste0(res,"/"),pattern="*.gdx")

for (gas in c("ch4","co2")) {
  
for (i in seq(start_job,min(end_job,nrow(data_srmpulse))) ) {
bsub <- paste("bsub", "-q",which_queue, "-n 1",
              "-P 0638", paste0("-J scenariosrmpulse", i,"_gas",gas), "-K -M 64G")

gams <- paste0("gams FAIR.gms --experiment=srm",
               " --gas=",gas,
              " --ecs=",data_srmpulse[i,]$ecs,
              " --tcr=",data_srmpulse[i,]$tcr,
              " --rcp=",data_srmpulse[i,]$rcp,
              " --pulse_time=",data_srmpulse[i,]$pulse,
              " --rate_of_cooling=",data_srmpulse[i,]$cool,
              " --start_rampdown=",data_srmpulse[i,]$term-100,
              " --end_rampdown=",data_srmpulse[i,]$term,
              " --start_rampup=",data_srmpulse[i,]$start,
              " --end_rampup=",data_srmpulse[i,]$start+100,
              " --termination_time=",data_srmpulse[i,]$term_delta,
              " --results_folder=",res)

results_name <-  paste0(data_srmpulse[i,]$rcp,
                        "_EXPsrmpulsemaskedterm_TER",data_srmpulse[i,]$term_delta,
                        "_GAS",gas,
                        "_ECS",data_srmpulse[i,]$ecs,
                        "_TCR",data_srmpulse[i,]$tcr,
                        "_PT",data_srmpulse[i,]$pulse,
                        "_RC",data_srmpulse[i,]$cool,
                        "_EC",data_srmpulse[i,]$term,
                        "_BC",data_srmpulse[i,]$start,
                        "_IChistorical_run")

cat("Checking scenario",i,"with gas",gas,"...\n")

if (!any(str_detect(str_remove(filelist,".gdx"),results_name)) ) {
  if (run_parallel==T) { 
    bsub <- str_remove(bsub, "-K ") 
    scenarios_launched <- scenarios_launched + 1 }
  command <- paste(bsub, gams)
  write(str_remove(command, "-K "), file = sh_file, append = TRUE)
  if (run_hpc==F) {command <- gams}
  if (run_parallel==T) {cat("Attention! This scenarios was not solved before. \n Launching scenario",i,"with gas",gas,"...\n")} else {cat("Attention! This scenarios was not solved before. \n Solving scenario",i,"with gas",gas,"...\n")}
  if (run_scenarios==T) {
    ret <- system(command = command, intern = TRUE)
    if(!is.null(attr(ret,"status")) )  {
      cat("Careful! scenario",i,"with gas",gas,"couldn't solve...\n")
      problematic_data <- problematic_data %>% 
        bind_rows(c(data_srmpulse[i,],error=attr(ret,"status"),gas=gas)) }
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


