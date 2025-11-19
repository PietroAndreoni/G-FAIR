library(gdxtools)
library(dplyr)
library(stringr)

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

extract_names <- function(.x) {
  .x %>% 
    as_tibble() %>% 
    mutate(file=str_remove_all(gdx,".*(?=RCP)|.gdx"),
           gas=str_extract(file,"(?<=GAS).+?(?=_)"),
           rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
           cool_rate=str_extract(str_remove(file,"RCP"),"(?<=RC).+?(?=_)"),
           pulse_time=str_extract(gdx,"(?<=PT).+?(?=_)"),
           geo_end=str_extract(str_remove(file,"ECS"),"(?<=EC).+?(?=_)"),
           geo_start=str_extract(file,"(?<=BC).+?(?=_)"),
           ecs=str_extract(file,"(?<=ECS).+?(?=_)"),
           tcr=str_extract(file,"(?<=TCR).+?(?=_)"),
           term=str_extract(file,"(?<=TER).+?(?=_)"),
           experiment=str_extract(file,"(?<=EXP).+?(?=_)")) %>%
    ungroup()}

sanitize <- function(.x) {
.x %>% as_tibble() %>% 
  mutate(t=as.numeric(t) ) %>% 
    filter(t<=480)}


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

res <- ifelse(is.null(opts[["i"]]), "Results_montecarlo", as.character(opts["i"]) )
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

### extract scenarios 
filelist <- unlist(lapply(res, function(folder) {
  files <- list.files(paste0(folder, "/"), pattern = ".gdx", full.names = TRUE)
  return(files)}))

filelist <- filelist[stringr::str_detect(filelist,"EXP")]

cat("Sanitizing the data...\n")
### make sure to include only files with the full scenario matrix
all_scenarios <- data.frame(filelist) %>% 
  rename(gdx=filelist) %>% 
  extract_names(.) 

# check that: (1) srm, srmpulse, srmpulsemasked, srmpulsemaskedterm are equal in number
check1 <- all_scenarios %>% 
  filter(!experiment %in% c("base","pulse") ) %>% 
  group_by(gas,rcp,cool_rate,pulse_time,geo_start,geo_end,ecs,tcr) %>% 
  summarise(tot=n(), diff=any(duplicated(experiment))) %>%
  ungroup() %>% 
  filter(tot==4 & diff==F)

# check that each of the srm experiments has the pulse and the base experiments associated
check2 <-  inner_join(all_scenarios %>% select(-term),check1 %>% select(-tot,-diff)) %>% 
  inner_join( all_scenarios %>% filter(experiment=="pulse") %>% select(gas,rcp,ecs,tcr,pulse_time)) %>% 
  inner_join( all_scenarios %>% filter(experiment=="base") %>% select(rcp,ecs,tcr))
  
# include only relevant scenarios
sanitized_names <- inner_join(all_scenarios,check2) %>% 
  bind_rows(inner_join(all_scenarios %>% filter(experiment=="pulse"),check2 %>% select(gas,rcp,ecs,tcr,pulse_time) ) %>%  unique()) %>% 
  bind_rows(inner_join(all_scenarios %>% filter(experiment=="base"),check2 %>% select(rcp,ecs,tcr) ) %>%  unique())

# discard duplicates in case any (possible in multiple folders)
sanitized_names <- sanitized_names %>% 
  group_by(file) %>% 
  filter(n()==1) 

files <- sanitized_names$gdx 

cat("Careful!", nrow(anti_join(all_scenarios,sanitized_names)), "scenarios out of",nrow(all_scenarios) ,"have been removed.\n")


cat("Loading the data...\n")

TATM <- gdxtools::batch_extract("TATM", files=files)$TATM %>%
  inner_join(sanitized_names) %>% sanitize()
  
W_EMI <- gdxtools::batch_extract("W_EMI", files=files)$W_EMI %>%
  inner_join(sanitized_names) %>% sanitize() 

FORC <- gdxtools::batch_extract("FORCING", files=files)$FORCING  %>% 
  inner_join(sanitized_names) %>% sanitize() 

SRM <- gdxtools::batch_extract("SRM", files=files)$SRM  %>%
  inner_join(sanitized_names) %>% sanitize()

background_srm <- gdxtools::batch_extract("forcing_srm", files=files)$forcing_srm  %>%
  inner_join(sanitized_names) %>% sanitize()

tot_forcing <- FORC %>% 
  group_by(t,gas,rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,term,experiment) %>% 
  summarise(value=sum(value)) %>% ungroup() 


cat("Loading the IDs from Montecarlo id file...")

id_montecarlo <- do.call(rbind, lapply(res, function(folder) {
    data <- read.csv(paste0(folder, "/id_montecarlo.csv"))
#    data$folder <- folder  # add folder info as a new column (optional but useful)
    return(data)
  }))

id_montecarlo <- id_montecarlo %>% 
  select(ecs,tcr,rcp,pulse,cool,term,start,term_delta) %>% 
  as_tibble() %>% 
  mutate(term=as.integer(term)) %>% 
  mutate_if(is.integer, as.character) %>% 
  mutate(rcp=str_remove(rcp,"RCP")) %>% 
  rename(pulse_time=pulse,cool_rate=cool,geo_end=term,geo_start=start,term=term_delta) %>% 
  unique()
  
n_scenarios <- nrow(id_montecarlo)

set.seed(seed)

#### add uncertaint post-solve parameters
# theta
id_montecarlo$theta <- round(drawln(15,7,n_scenarios),0)

# alpha
id_montecarlo$alpha <- round(drawln(0.00575,0.00575*150/(230-100),n_scenarios),5)

# delta 
id_montecarlo$delta <- round(runif(n_scenarios,0.001,0.07),3)

# prob 
id_montecarlo$prob <- round(runif(n_scenarios,0,1),2)

# additional SAI mortality (thousands death x tGS)
id_montecarlo$mortality_srm <- round(drawln(7400,(16000-2300)/1.96,n_scenarios),0)

# pollution costs (ozone mortality)
id_montecarlo$mortality_ozone <- round(EnvStats::rtri(n_scenarios, min = 100, max = 107, mode = 104),1)

# value of statistical life
id_montecarlo$vsl <- round(runif(n_scenarios,1,10),0)*10^6

# gdp growth
id_montecarlo$dg <- round(rnorm(n_scenarios,mean=0.015,sd=0.005),4)

# this exclude
all_names <- c("gas",names(id_montecarlo))
scenario_names <- names(sanitized_names)
scenario_names <- scenario_names[scenario_names != c("gdx","file","experiment")] 

cat("Calculating net present cost: part1...\n")

damnpv_pre <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemasked") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemasked")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemasked")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemasked" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_pulse = value) %>% 
               filter(experiment=="srmpulsemasked") %>% 
               select_at(c("t",scenario_names,"ozone_pulse"))) %>%
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_base = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(c("t",scenario_names,"ozone_base"))) %>%
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(c("t",scenario_names,"temp_srm")) ) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select_at(c("t",scenario_names,"temp_srmpulse"))) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  select(-term) %>% 
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t<=as.numeric(term)) %>% 
  rowwise() %>% mutate(gwpt = min(gwp* (1 +dg)^(280-1),gwp* (1 + dg)^(t-1)) ) %>% 
  mutate(direct_cost = srm_masking * forctoUSD,
         srm_pollution = vsl * srm_masking * forctoTg * mortality_srm,
         tropoz_pollution = vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc,
         imperfect_masking = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1+srm/forc) ) %>%
  group_by_at( all_names ) %>%
  summarise(dirnpv = sum( direct_cost / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            srmpnpv = sum( srm_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            ozpnpv = sum( tropoz_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            masknpv = sum( imperfect_masking / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>% 
  ungroup() %>% mutate(costnpv=dirnpv + srmpnpv + ozpnpv + masknpv)

cat("Calculating net present cost: part2...\n")

damnpv_post_noterm <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemasked") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemasked")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemasked")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemasked" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_pulse = value) %>% 
               filter(experiment=="srmpulsemasked") %>% 
               select_at(c("t",scenario_names,"ozone_pulse"))) %>%
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_base = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(c("t",scenario_names,"ozone_base"))) %>%
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(c("t",scenario_names,"temp_srm")) ) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select_at(c("t",scenario_names,"temp_srmpulse"))) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  select(-term) %>% 
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t>as.numeric(term)) %>% 
  rowwise() %>% mutate(gwpt = min(gwp* (1 +dg)^(280-1),gwp* (1 + dg)^(t-1)) ) %>% 
  mutate(direct_cost = srm_masking * forctoUSD,
         srm_pollution = vsl * srm_masking * forctoTg * mortality_srm,
         tropoz_pollution = vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc,
         imperfect_masking = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1+srm/forc) ) %>%
  group_by_at( all_names ) %>%
  summarise(dirnpv = sum( direct_cost / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            srmpnpv = sum( srm_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            ozpnpv = sum( tropoz_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            masknpv = sum( imperfect_masking / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>% 
  ungroup() %>% mutate(costnpv=dirnpv + srmpnpv + ozpnpv + masknpv)

cat("Calculating net present cost: part3...\n")

damnpv_post_term <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemaskedterm") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemaskedterm")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemaskedterm")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemaskedterm" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_pulse = value) %>% 
               filter(experiment=="srmpulsemasked") %>% 
               select_at(setdiff(c("t", scenario_names, "ozone_pulse"), "term")) ) %>%
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_base = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(setdiff(c("t", scenario_names, "ozone_base"), "term")) ) %>%
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(setdiff(c("t", scenario_names, "temp_srm"), "term")) ) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select_at(setdiff(c("t", scenario_names, "temp_srmpulse"), "term")) ) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t>as.numeric(term)) %>% 
  rowwise() %>% mutate(gwpt = min(gwp* (1 +dg)^(280-1),gwp* (1 + dg)^(t-1)) ) %>% 
  mutate(direct_cost = srm_masking * forctoUSD,
         srm_pollution = vsl * srm_masking * forctoTg * mortality_srm,
         tropoz_pollution = vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc,
         dam = gwp*(alpha * (temp**2-temp_srm**2)),
         imperfect_masking = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1+srm/forc) ) %>%
  group_by_at( all_names ) %>%
  summarise(dirnpv = sum( direct_cost / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            srmpnpv = sum( srm_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            ozpnpv = sum( tropoz_pollution / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            masknpv = sum( imperfect_masking / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>% 
  ungroup() %>% mutate(costnpv=dirnpv + srmpnpv + ozpnpv + masknpv + damnpv)

cat("Calculating net present cost: putting things together...\n")

pulse_size <- W_EMI %>% 
  filter(ghg==gas) %>% 
  filter(experiment %in% c("srm","srmpulse") ) %>% 
  group_by_at(c("t",scenario_names)) %>% 
  mutate(pulse_size=value-value[experiment=="srm"]) %>% 
  filter(experiment=="srmpulse" & pulse_size!=0) %>% 
  ungroup() %>% 
  select_at(setdiff(c(scenario_names,"pulse_size"),"term")) %>% 
  mutate(pulse_size=ifelse(gas=="co2",pulse_size*1e9,pulse_size*1e6))

scc <- TATM %>%  rename(temp_srm = value) %>% 
  filter(experiment=="srm") %>% 
  select_at(c("t",scenario_names,"temp_srm")) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select_at(c("t",scenario_names,"temp_srmpulse")) ) %>%
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_pulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select_at(setdiff(c("t", scenario_names, "ozone_pulse"), "term")) ) %>%
  inner_join(FORC %>% filter(ghg=="o3trop") %>% rename(ozone_base = value) %>% 
               filter(experiment=="srm") %>% 
               select_at(setdiff(c("t", scenario_names, "ozone_base"), "term")) ) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  select(-term) %>% 
  inner_join(id_montecarlo %>% select_at(setdiff(all_names,c("gas","theta","term","prob","mortality_srm" ) ) ) ) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time)) %>% 
  rowwise() %>% mutate(gwpt = min(gwp* (1 +dg)^(280-1),gwp* (1 + dg)^(t-1)) ) %>% 
  mutate(tropoz_pollution = vsl * mortality_ozone * (ozone_pulse - ozone_base) * ozone_rftoconc,
         dam = gwp*( alpha * ((temp_srmpulse)**2-temp_srm**2)) ) %>% 
  group_by_at(setdiff(all_names,c("gas","theta","term","prob","mortality_srm" ) ) ) %>%
  summarise(damnpv = sum( (dam + tropoz_pollution) / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>%
  full_join(pulse_size) %>% 
  mutate(scc=damnpv/pulse_size ) %>% 
  select(setdiff(c(all_names,"scc"),c("theta","term","prob","mortality_srm" ) ) )

damnpv <- damnpv_pre %>% select_at(c(all_names,"costnpv")) %>%  rename(costpre=costnpv) %>% 
  full_join(damnpv_post_noterm %>% select_at(c(all_names,"costnpv")) %>%  rename(costpostnoterm=costnpv)) %>% 
  full_join(damnpv_post_term %>%  select_at(c(all_names,"costnpv")) %>%  rename(costpostterm=costnpv)) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(costnpv=costpre+as.numeric(prob)*costpostterm+(1-as.numeric(prob))*costpostnoterm) %>% 
  select(-costpre,-costpostterm,-costpostnoterm) %>% 
  full_join(pulse_size) %>% 
  full_join(scc) %>% 
  mutate(npc_srm=costnpv/pulse_size)

cat("Saving output...\n")

write.csv(damnpv,file=paste0(output_folder,"/output_analysis.csv"))

### produce plots (optional)
if(plot_results==T) {

damnpv <- bind_rows(lapply(list.files(output_folder, pattern = "\\.csv$", full.names = TRUE), read.csv)) %>% select(-X)
  
damnorm <- damnpv %>% 
  filter(scc>0 & npc_srm>0 & delta>0) %>% 
  group_by(gas) %>% 
  mutate(npc_norm=npc_srm/median(npc_srm,na.rm=TRUE))

damnorm_diff <- damnorm %>% 
  filter(scc>0 & npc_srm>0 & delta>0) %>%
  select(-scc,-npc_srm,-pulse_size,-costnpv) %>% 
  pivot_wider(values_from="npc_norm",names_from="gas") %>% 
  unnest(ch4,co2) %>% 
  filter(!is.null(ch4) & !is.null(co2)) %>% 
  mutate(norm_diff=co2-ch4)

damnorm %>% 
  group_by(gas) %>% 
  summarise(p66=quantile(npc_norm,0.66,na.rm=TRUE), 
            p75=quantile(npc_norm,0.75,na.rm=TRUE),
            p90=quantile(npc_norm,0.9,na.rm=TRUE), 
            p95=quantile(npc_norm,0.95,na.rm=TRUE), 
            p99=quantile(npc_norm,0.99,na.rm=TRUE), 
            p999=quantile(npc_norm,0.999,na.rm=TRUE)) 

stats::ks.test(damnorm %>% filter(gas=="co2") %>% pull(npc_norm),
               damnorm %>% filter(gas=="ch4") %>% pull(npc_norm))
  
ggplot(damnorm) +
  geom_density(aes(x=npc_norm,
                   color=gas),
               adjust=10,
               linewidth=1) +
  geom_point(aes(x=npc_norm,
                 color=gas,y=0),
             shape=108) +
  coord_cartesian(xlim = c(0, quantile(damnorm$npc_norm, 0.95, na.rm=TRUE)) ) +
  xlab("Normalized present cost") + ylab("density") 

ggplot(damnorm %>% filter(gas=="ch4")) +
  geom_density(aes(x=npc_srm),adjust=5,color="red") +
  coord_cartesian(xlim = c(0, quantile( (damnorm %>% filter(gas=="ch4"))$npc_srm, 0.95, na.rm=TRUE)) ) +
  xlab("Net present cost [$/tonCH4]") + ylab("")

ggplot(damnorm) +
  geom_smooth(aes(x=delta,
                  y=npc_norm,
                  color=gas), 
              method="loess") +
  geom_density(aes(x=delta,
                   y=after_stat(scaled)) ) +
  xlab("Normalized present cost") + 
  ylab("density") + 
  coord_cartesian(ylim=c(0,10),
                  xlim=c(0,45))


ggplot(damnorm) +
  geom_smooth(aes(x=alpha,
                  y=npc_norm,
                  color=gas), method="loess" ) +
  geom_density(aes(x=alpha,
                   y=after_stat(scaled)) ) +
  xlab("Normalized present cost") + 
  ylab("density") + 
  coord_cartesian(ylim=c(0,10),
                  xlim=c(0,0.03))


ggplot(damnorm) +
  geom_smooth(aes(x=term-pulse_time,y=npc_norm,color=gas), method="loess" ) +
  geom_density(aes(x=term-pulse_time,y = after_stat(scaled)) ) +
  xlab("Normalized present cost") + 
  ylab("density") + 
  coord_cartesian(ylim=c(0,6))
           
                  

ggplot(damnorm) +
  geom_smooth(aes(x=term-pulse_time,y=npc_norm,color=gas), method="loess" ) +
  geom_density(aes(x=term-pulse_time,y = after_stat(scaled)) ) +
  xlab("Normalized present cost") + 
  ylab("density") + 
  coord_cartesian(ylim=c(0,6),
                  xlim=c(0,300))

ggplot(damnorm) +
  geom_smooth(aes(x=prob,y=npc_norm,color=gas), method="loess" ) +
  geom_density(aes(x=prob,y=after_stat(scaled)) ) +
  xlab("Normalized present cost") + 
  ylab("density") + 
  coord_cartesian(ylim=c(0,6),
                  xlim=c(0,1))

ggplot(damnorm) +
  geom_smooth(aes(x=tcr,y=npc_norm,color=gas)) +
  xlab("Normalized present cost") + ylab("density")

library(gsaot)
damnpv <- damnpv %>%  
  ungroup() %>% 
  select(-X) %>% 
  mutate(across(c(tcr, ecs, term, theta, pulse_time, geo_start, geo_end, term, delta, prob), as.numeric))

gsoat_data <- damnpv %>% filter(gas=="ch4" & scc>0 & npc_srm>0 & delta>0) 

stat_analysis_ch4 <- ot_indices_1d(gsoat_data %>% 
                                     select(rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                   gsoat_data %>% pull(npc_srm ), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_ch4 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

gsoat_data <- damnpv %>% filter(gas=="co2" & scc>0 & npc_srm>0 & delta>0) %>% ungroup()

stat_analysis_co2 <- ot_indices_1d(gsoat_data %>% 
                                     select(rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_co2 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")


stat_analysis_norm <- ot_indices_1d(damnorm %>% 
                                     select(gas, rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                    damnorm %>% pull(npc_norm ), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_norm <- irrelevance_threshold(damnorm %>% pull(npc_norm), M= 15, solver="1d")


ggplot(as_tibble(stat_analysis_norm$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_norm$indices) +
  geom_errorbar(aes(x=input,
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [ch4]") + xlab("") + 
  theme(legend.position = "none")

ggplot(as_tibble(stat_analysis_ch4$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_ch4$indices) +
  geom_errorbar(aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [ch4]") + xlab("") + 
  theme(legend.position = "none")

ggplot(as_tibble(stat_analysis_co2$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_co2$indices) +
  geom_errorbar(aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [co2]") + xlab("") + 
  theme(legend.position = "none")

}



stat_analysis_diff <- ot_indices_1d(damnorm_diff %>% 
                                     select(rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                    damnorm_diff %>% pull(norm_diff), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_diff <- irrelevance_threshold(damnorm_diff %>% pull(norm_diff), M= 15, solver="1d")

ggplot(as_tibble(stat_analysis_diff$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_diff$indices) +
  geom_errorbar(aes(x=input,
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [difference]") + xlab("") + 
  theme(legend.position = "none")
