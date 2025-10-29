# this scripts generates the realization for the montecarlo parameters
require(dplyr)
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

generate_data <- T
overwrite_data <- T
res <- "Results_montecarlo"

if (generate_data==T) {

if (overwrite_data==T) {data <- data.frame()} else {data <-  as.data.frame(read.csv("id_montecarlo.csv")) %>% select(-X)}

max_id <- nrow(data)

for (i in seq(1,10000,by=1)) {
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
  delta <- sample(seq(0,0.07,by=0.01),1)

# prob 
  prob <- sample(seq(0,0.1,by=0.01),1)

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
                       prob=prob) ) }

write.csv(data,file="id_montecarlo.csv") } else {
  
data <- as.data.frame(read.csv("id_montecarlo.csv")) %>% select(-X)

}
 
# select unique scenarios 
data_srmpulse <- data %>% 
  select(ecs, tcr, rcp, pulse, cool, term, start) %>% unique() %>% 
  mutate()

for (i in seq(1,nrow(data_srmpulse)) ) {
bsub <- paste("bsub", "-q p_short", "-n 1",
              "-P  0638", paste0("-J scenariosrmpulse_", i), "-K -M 64G")

gams <- paste0("gams FAIR.gms --experiment=srm",
              " --ecs=",data_srmpulse[i,]$ecs,
              " --tcr=",data_srmpulse[i,]$tcr,
              " --rcp=",data_srmpulse[i,]$rcp,
              " --pulse_time=",data_srmpulse[i,]$pulse,
              " --rate_of_cooling=",data_srmpulse[i,]$cool,
              " --start_rampdown=",data_srmpulse[i,]$term,
              " --end_rampdown=",data_srmpulse[i,]$term+100,
              " --start_rampup=",data_srmpulse[i,]$start,
              " --end_rampup=",data_srmpulse[i,]$start+100)
command <- paste(bsub, gams)
# command <- gams
ret <- system(command = command, intern = TRUE)

for (t in seq(10,100,by=10)) {
  bsub <- paste("bsub", "-q p_short", "-n 1",
                "-P 0638", paste0("-J scenariosrmpulse_", i,"TERM",t), "-K -M 64G")
  
  gams <- paste0("gams FAIR.gms --experiment=srm",
                 " --ecs=",data_srmpulse[i,]$ecs,
                 " --tcr=",data_srmpulse[i,]$tcr,
                 " --rcp=",data_srmpulse[i,]$rcp,
                 " --pulse_time=",data_srmpulse[i,]$pulse,
                 " --rate_of_cooling=",data_srmpulse[i,]$cool,
                 " --start_rampdown=",data_srmpulse[i,]$term,
                 " --end_rampdown=",data_srmpulse[i,]$term+100,
                 " --start_rampup=",data_srmpulse[i,]$start,
                 " --end_rampup=",data_srmpulse[i,]$start+100,
                 " --termination_time=",t,
                 " --results_folder=",res)
  command <- paste(bsub, gams)
  # command <- gams
  ret <- system(command = command, intern = TRUE)
}

}

data_pulse <- data %>% 
  select(ecs, tcr, rcp, pulse) %>% 
  unique() 

for (i in seq(1,nrow(data_pulse)) ) {
  bsub <- paste("bsub", "-q p_short", "-n 1",
                "-P 0588", paste0("-J scenariopulse_", i), "-K -M 64G")
  
  gams <- paste0("gams FAIR.gms --experiment=pulse",
                 " --ecs=",data_pulse[i,]$ecs,
                 " --tcr=",data_pulse[i,]$tcr,
                 " --rcp=",data_pulse[i,]$rcp,
                 " --pulse_time=",data_pulse[i,]$pulse)
  command <- paste(bsub, gams)
  # command <- gams
  ret <- system(command = command, intern = TRUE)
}

