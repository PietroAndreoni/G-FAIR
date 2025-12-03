# this scripts generates the realizations for the montecarlo parameters
# generates a script with all the scenarios to run (excluding those with results)
# and runs the scenarios in batches 
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Montecarlo.R [-g <generate_data>] [-o <res>] [-q <which_queue>] [-n <n_scenarios>] [-m <max_scenarios>] [-x <rerun_problem>] [-w <overwrite_data>] [-p <run_parallel>] [-r <run_scenarios>] [-s <start_job>] [-e <end_job>] [--hpc <run_hpc>] [--seed <seed>] 

Options:
-o <res>              Path where the results are (default: Results)
-n <n_scenarios>      Number of new scenarios to generate
-g <generate_data>    T/F to generate new realizations
-w <overwrite_data>   T/F to overwrite old realizations
-r <run_scenarios>    T/F to run or just generate scenarios (in .sh script)
-p <run_parallel>     T/F if to run in parallel or in series (for Juno)
-s <start_job>        Number of line to start calling the scenarios 
-e <end_job>          End of line to start calling the scenarios
-q <which_queue>      Select queue to send the jobs to (for parallel solving)
-m <max_scenarios>    Maximum scenarios to send to solve
-x <rerun_problem>    T to avoid rerunning problematic scenarios (default)
--seed <seed>         seed number (for reproducibility)
--hpc <run_hpc>       T/F if to run on Juno or local (windows)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Montecarlo')

fit_distribution <- function(distribution = "lognormal",
                             n =1, #number of required realizations
                             median = NULL,
                             mean   = NULL,
                             sd     = NULL,
                             q33    = NULL,
                             q66    = NULL,
                             q5     = NULL,
                             q95    = NULL,
                             probs_33_66 = c(0.33, 0.66),
                             probs_5_95  = c(0.05, 0.95),
                             plot=F) {
  
  mu    <- NA_real_
  sigma <- NA_real_
  
  if (distribution == "lognormal") {
    # Helper to check positivity (lognormal only makes sense for > 0)
    check_pos <- function(x, name) {
      if (!is.null(x) && any(x <= 0)) {
        stop(name, " must be > 0 for a lognormal distribution.")
      }
    }
    
    check_pos(median, "median")
    check_pos(mean,   "mean")
    check_pos(sd,     "sd")
    check_pos(q33,    "q33")
    check_pos(q66,    "q66")
    check_pos(q5,     "q5")
    check_pos(q95,    "q95")
    
    ## 1) Case: mean and sd of X
    if (!is.null(mean) && !is.null(sd)) {
      if (sd <= 0) stop("sd must be > 0.")
      v      <- sd^2
      sigma2 <- log(1 + v / mean^2)
      if (sigma2 <= 0) stop("Inconsistent mean and sd for a lognormal.")
      sigma  <- sqrt(sigma2)
      mu     <- log(mean) - 0.5 * sigma2
      
      ## 2) Case: mean and median of X
    } else if (!is.null(mean) && !is.null(median)) {
      mu     <- log(median)
      sigma2 <- 2 * (log(mean) - mu)
      if (sigma2 <= 0) stop("Inconsistent mean and median for a lognormal.")
      sigma  <- sqrt(sigma2)
      
      ## 3) NEW: median and sd of X
    } else if (!is.null(median) && !is.null(sd)) {
      if (sd <= 0) stop("sd must be > 0.")
      mu <- log(median)
      v  <- sd^2
      
      # For a lognormal:
      # var = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
      # Let t = exp(sigma^2); then var = (t - 1) * t * exp(2*mu)
      # => (t^2 - t) = var / exp(2*mu)
      # Solve t^2 - t - A = 0, where A = var / median^2
      A    <- v / median^2
      disc <- 1 + 4 * A
      t    <- (1 + sqrt(disc)) / 2   # positive root
      sigma2 <- log(t)
      if (sigma2 <= 0) stop("Inconsistent median and sd for a lognormal.")
      sigma  <- sqrt(sigma2)
      
      ## 4) Case: median and 33–66% quantiles
    } else if (!is.null(median) && !is.null(q33) && !is.null(q66)) {
      mu <- log(median)
      p1 <- probs_33_66[1]
      p2 <- probs_33_66[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      sigma_33 <- (log(q33) - mu) / z1
      sigma_66 <- (log(q66) - mu) / z2
      sigma    <- mean(c(sigma_33, sigma_66))
      
      ## 5) Case: median and 5–95% quantiles
    } else if (!is.null(median) && !is.null(q5) && !is.null(q95)) {
      mu <- log(median)
      p1 <- probs_5_95[1]
      p2 <- probs_5_95[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      sigma_5  <- (log(q5)  - mu) / z1
      sigma_95 <- (log(q95) - mu) / z2
      sigma    <- mean(c(sigma_5, sigma_95))
      
      ## 6) Case: 33–66% quantiles only
    } else if (!is.null(q33) && !is.null(q66)) {
      p1 <- probs_33_66[1]
      p2 <- probs_33_66[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      L1 <- log(q33)
      L2 <- log(q66)
      
      sigma <- (L2 - L1) / (z2 - z1)
      mu    <- L1 - sigma * z1
      
      ## 7) Case: 5–95% quantiles only
    } else if (!is.null(q5) && !is.null(q95)) {
      p1 <- probs_5_95[1]
      p2 <- probs_5_95[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      L1 <- log(q5)
      L2 <- log(q95)
      
      sigma <- (L2 - L1) / (z2 - z1)
      mu    <- L1 - sigma * z1
      
    } else {
      stop("Not enough information for lognormal. Provide one of:
      (mean & sd),
      (mean & median),
      (median & sd),
      (median & q33 & q66),
      (median & q5 & q95),
      (q33 & q66), or
      (q5 & q95).")
    }
    
    # Implied stats for lognormal X
    out_mean   <- exp(mu + 0.5 * sigma^2)
    out_median <- exp(mu)
    out_sd     <- sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
    
    
    if (plot==T) {
      draws <- rlnorm(n = 10000, mu, sigma)
      print(paste("mean: ",round(mean(draws),5) ) )
      print(paste("median: ",round(median(draws), 5) ) )
      print(paste("std: ",round(sd(draws), 5) ) )
      print(paste("66%: ",round(quantile(draws,0.66), 5) ) )
      print(paste("33%: ",round(quantile(draws,0.33), 5) ) )
      print(paste("95%: ",round(quantile(draws,0.95), 5) ) )
      print(paste("5%: ",round(quantile(draws,0.05), 5) ) )
      plot(density(draws[draws > 0 & draws < 10*sd(draws)])) }
    
    return(rlnorm(n, mu, sigma)) 
    
    
  } else if (distribution == "normal") {
    
    if (!is.null(sd) && sd <= 0) stop("sd must be > 0.")
    
    ## 1) mean & sd
    if (!is.null(mean) && !is.null(sd)) {
      mu    <- mean
      sigma <- sd
      
      ## 2) median & sd
    } else if (!is.null(median) && !is.null(sd)) {
      mu    <- median
      sigma <- sd
      
      ## 3) median and 33–66% quantiles
    } else if (!is.null(median) && !is.null(q33) && !is.null(q66)) {
      mu <- median
      p1 <- probs_33_66[1]
      p2 <- probs_33_66[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      sigma_33 <- (q33 - mu) / z1
      sigma_66 <- (q66 - mu) / z2
      sigma    <- mean(c(sigma_33, sigma_66))
      
      ## 4) median and 5–95% quantiles
    } else if (!is.null(median) && !is.null(q5) && !is.null(q95)) {
      mu <- median
      p1 <- probs_5_95[1]
      p2 <- probs_5_95[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      sigma_5  <- (q5  - mu) / z1
      sigma_95 <- (q95 - mu) / z2
      sigma    <- mean(c(sigma_5, sigma_95))
      
      ## 5) 33–66% quantiles only
    } else if (!is.null(q33) && !is.null(q66)) {
      p1 <- probs_33_66[1]
      p2 <- probs_33_66[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      sigma <- (q66 - q33) / (z2 - z1)
      mu    <- q33 - sigma * z1
      
      ## 6) 5–95% quantiles only
    } else if (!is.null(q5) && !is.null(q95)) {
      p1 <- probs_5_95[1]
      p2 <- probs_5_95[2]
      z1 <- qnorm(p1)
      z2 <- qnorm(p2)
      
      sigma <- (q95 - q5) / (z2 - z1)
      mu    <- q5 - sigma * z1
      
    }  else {
      stop("Not enough information for lognormal. Provide one of:
      (mean & sd),
      (mean & median),
      (median & sd),
      (median & q33 & q66),
      (median & q5 & q95),
      (q33 & q66), or
      (q5 & q95).")
    }
    
    # Return also some implied summary stats for convenience
    out_mean <- exp(mu + 0.5 * sigma^2)
    out_median <- exp(mu)
    out_sd <- sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
    
    
    if (plot==T) {
      draws <- rnorm(n = 10000, mu, sigma)
      print(paste("mean: ",round(mean(draws),5) ) )
      print(paste("median: ",round(median(draws), 5) ) )
      print(paste("std: ",round(sd(draws), 5) ) )
      print(paste("66%: ",round(quantile(draws,0.66), 5) ) )
      print(paste("33%: ",round(quantile(draws,0.33), 5) ) )
      print(paste("95%: ",round(quantile(draws,0.95), 5) ) )
      print(paste("5%: ",round(quantile(draws,0.05), 5) ) )
      plot(density(draws[draws > 0 & draws < 10*sd(draws)])) }
    
    return(rnorm(n, mu, sigma)) 
    
    
  }  else {stop("choose normal or lognormal distribution")}
  
}

# logical
generate_data = ifelse(is.null(opts[["g"]]), T, as.logical(opts["g"]) )
overwrite_data = ifelse(is.null(opts[["w"]]), T, as.logical(opts["w"]) )
run_scenarios = ifelse(is.null(opts[["r"]]), T, as.logical(opts["r"]) )
run_parallel = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )
run_hpc = ifelse(is.null(opts[["hpc"]]), F, as.logical(opts["hpc"]) )
rerun_problem = ifelse(is.null(opts[["x"]]),T, as.logical(opts["x"]) )

# numeric
n_scenarios = ifelse(is.null(opts[["n"]]), 5, as.numeric(opts["n"]) )
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
  ecs <- round(fit_distribution(median=3,q5=2,q95=5)*10,0)

# tcr
  tcr <- round(fit_distribution(median=1.8,q5=1,q95=2.5)*10,0)

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
                       term_delta=time_term) ) }

data <- data %>% 
  mutate(term_delta=pulse+term_delta,
         term=ifelse(cool==0,2700,term),
         start=ifelse(cool==0,2700,start) ) %>% 
  unique()

write.csv(data,file=paste0(res,"/id_montecarlo.csv")) } else {
  
data <- as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv"))) %>% select(-X)

}

if (rerun_problem==F) data <- data %>% anti_join(problematic_data %>% select(-error,-gas) %>% unique())
cat("Launching jobs... \n")

filelist <- list.files(path=paste0(res,"/"),pattern="*.gdx")

for (gas in c("ch4","co2")) {
  
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
  write(str_remove(command, "-K "), file = sh_file, append = TRUE)
  if (run_hpc==F) {command <- gams}
  if (run_parallel==T) {cat("Attention! This scenarios was not solved before. \n Launching scenario",i,"with gas",gas,"...\n")} else {cat("Attention! This scenarios was not solved before. \n Solving scenario",i,"with gas",gas,"...\n")}
  if (run_scenarios==T) {
    ret <- system(command = command, intern = TRUE)
    if(!is.null(attr(ret,"status")) )  {
      if(attr(ret,"status")!=112) {
        cat("Careful! scenario",i,"with gas",gas,"couldn't solve...\n")
        problematic_data <- problematic_data %>% 
          bind_rows(c(data[i,],error=attr(ret,"status"),gas=gas)) } 
      }
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


