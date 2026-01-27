# this scripts generates the realizations for the montecarlo parameters
require(dplyr)
require(stringr)

# Usage
'Launch montecarlo script for SRM substitution pulse analysis

Usage:
  Generate_montecarlo.R [-o <res>] [-n <n_scenarios>] [-w <overwrite_data>] [-p <run_parallel>] [-s <start_job>] [-e <end_job>] [--hpc <run_hpc>] [--base <main_scenario>] 

Options:
-o <res>                     Path of the input/outputs
-n <n_scenarios>             Number of new scenarios to generate
-w <overwrite_data>          T/F to overwrite old realizations
--seed <seed>                seed number (for reproducibility)
--base <main_scenario>       T/F if to run the main scenario only (no montecarlo over policy parameters)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Generate_montecarlo')

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
overwrite_data = ifelse(is.null(opts[["w"]]), T, as.logical(opts["w"]) )
main_scenario = ifelse(is.null(opts[["base"]]), F, as.logical(opts["base"]) )

# numeric
n_scenarios = ifelse(is.null(opts[["n"]]), 10000, as.numeric(opts["n"]) )
seed = ifelse(is.null(opts[["seed"]]), 123, as.integer(opts["seed"]) )

# strings
res = ifelse(is.null(opts[["o"]]), "Montecarlo", as.character(opts["o"]) )

# Make sure the file exists (create it if not)
if (!dir.exists(res)) {
  dir.create(res)
}

cat("Generating data... \n")

set.seed(seed)

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
    pt <- sample(c(5,10,20,30,40,50,60,70,80),1)
    
    # cooling rate 
    cool <- sample(seq(0,40,by=5),1)
    
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
  
if (main_scenario==T) {
    data <- data %>% 
      mutate(rcp="RCP45", cool=10, term=2400, start=2025, pulse=5 ) }
  
write.csv(data,file=paste0(res,"/id_montecarlo.csv")) 

