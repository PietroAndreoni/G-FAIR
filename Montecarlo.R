# this scripts generates the realizations for the montecarlo parameters
# generates a script with all the scenarios to run (excluding those with results)
# and runs the scenarios in batches 
require(dplyr)
require(stringr)

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

generate_data <- F # generate new data
overwrite_data <- F #overwrite old data
n_scenarios <- 1000 # number of (new) scenarios to generate
run_batch <- T # run scenarios 
run_parallel <- F # launches all scenarios in parallel (careful not to flood Juno..)
run_hpc <- T # run from juno (F for local machine)
start_job <- 250 # beginning of 
end_job <- 300
res <- "Results_montecarlo"

# Define the path to the .ssh file
sh_file <- "Montecarlo.sh"   # or any other file you want to modify

# Make sure the file exists (create it if not)
if (!file.exists(sh_file)) {
  file.create(sh_file)
} else {
  file.remove(sh_file)
  file.create(sh_file)}

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res)
}

if (generate_data==T) {

if (overwrite_data==T) {data <- data.frame()} else {data <-  as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv")) ) %>% select(-X)}

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
  select(ecs, tcr, rcp, pulse, cool, term, start, term_delta) %>% unique() %>% 
  mutate()

data_pulse <- data %>% 
  select(ecs, tcr, rcp, pulse) %>% 
  unique() 

filelist <- list.files(path=paste0(res,"/"),pattern="*.gdx")

for (gas in c("ch4","co2")) {
  
for (i in seq(start_job,min(end_job,nrow(data_srmpulse))) ) {
bsub <- paste("bsub", "-q p_short", "-n 1",
              "-P 0638", paste0("-J scenariosrmpulse", i,"_gas",gas), "-K -M 64G")

gams <- paste0("gams FAIR.gms --experiment=srm",
               " --gas=",gas,
              " --ecs=",data_srmpulse[i,]$ecs,
              " --tcr=",data_srmpulse[i,]$tcr,
              " --rcp=",data_srmpulse[i,]$rcp,
              " --pulse_time=",data_srmpulse[i,]$pulse,
              " --rate_of_cooling=",data_srmpulse[i,]$cool,
              " --start_rampdown=",data_srmpulse[i,]$term,
              " --end_rampdown=",data_srmpulse[i,]$term+100,
              " --start_rampup=",data_srmpulse[i,]$start,
              " --end_rampup=",data_srmpulse[i,]$start+100,
              " --termination_time=",data_srmpulse[i,]$term_delta,
              " --results_folder=",res)

results_name <-  paste0(data_srmpulse[i,]$rcp,
                        "_EXPsrmpulsemasked_",
                        "_GAS",gas,
                        "ECS_",data_srmpulse[i,]$ecs,
                        "_TCR",data_srmpulse[i,]$tcr,
                        "_PT",data_srmpulse[i,]$pulse,
                        "_RC",data_srmpulse[i,]$cool,
                        "_EC",data_srmpulse[i,]$term,
                        "_BC",data_srmpulse[i,]$start,
                        "_IChistorical_run")

if (!any(str_detect(filelist,results_name)) ) {
  if (run_parallel==T) {bsub <- str_remove(bsub, "-K ")}
  command <- paste(bsub, gams)
  write(str_remove(command, "-K "), file = sh_file, append = TRUE)
  if (run_hpc==F) {command <- gams}
  if (run_batch==T) {ret <- system(command = command, intern = TRUE)}
  }

}

for (i in seq(start_job,min(end_job,nrow(data_pulse))) ) {
  bsub <- paste("bsub", "-q p_short", "-n 1",
                "-P 0638", paste0("-J scenariosrmpulse", i,"_gas",gas), "-M -K 64G")
  
  gams <- paste0("gams FAIR.gms --experiment=pulse",
                 " --gas=",gas,
                 " --ecs=",data_pulse[i,]$ecs,
                 " --tcr=",data_pulse[i,]$tcr,
                 " --rcp=",data_pulse[i,]$rcp,
                 " --pulse_time=",data_pulse[i,]$pulse,
                 " --results_folder=",res)
  
  results_name <-  paste0(data_srmpulse[i,]$rcp,
                          "_GAS",gas,
                          "_EXPpulse",
                          "_ECS_",data_srmpulse[i,]$ecs,
                          "_TCR",data_srmpulse[i,]$tcr,
                          "_PT",data_srmpulse[i,]$pulse,
                          "_IChistorical_run")
  
  if (!any(str_detect(filelist,results_name)) ) {
    if (run_parallel==T) {bsub <- str_remove(bsub, "-K ")}
    command <- paste(bsub, gams)
    write(str_remove(command, "-K "), file = sh_file, append = TRUE)
    if (run_hpc==F) {command <- gams}
    if (run_batch==T) {ret <- system(command = command, intern = TRUE)}
    }
  
}
  
}




