library(gdxtools)
library(dplyr)
library(stringr)

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
  Analyze_montecarlo.R [-i <res>] [-o <output_folder>] [-h <run_hpc>] [-p <plot_results>] 

Options:
-i <res>              Path where the results are (default: Results_montecarlo). For multiple folders separate with -
-o <output_folder>   Where to save results
-h <run_hpc>          T/F if running from Juno (T) or local (F) 
-p <plot_results>     T/F to plot data and save data analysis
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Montecarlo')

res = ifelse(is.null(opts[["i"]]), "Results_montecarlo", as.character(opts["i"]) )
res <- str_split(res,"-")[[1]]

output_folder = ifelse(is.null(opts[["o"]]), "Results_output", as.character(opts["o"]) )

# Make sure the output folder exists (create it if not)
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

if (any(!dir.exists(res)) ) stop("some of the folder specified do not exsist")
  
run_hpc = ifelse(is.null(opts[["h"]]), T, as.logical(opts["h"]) )
plot_results = ifelse(is.null(opts[["p"]]), F, as.logical(opts["p"]) )

if(run_hpc==F) {igdx()} else {igdx("/work/cmcc/pa12520/gams40.4_linux_x64_64_sfx")}

## climate damage function parameters
gwp <- 105*1e12 # initial world gdp
cost <- 0.001 # 0.1% GDP per W/m2 as Belaiia et al
forctoTg <- 1/0.2 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
TgtoUSD <- 2250*10^6 # from https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7/pdf  
forctoUSD <- forctoTg * TgtoUSD # US$/(W/m^2)
gwpmax <- gwp* (1 + 0.022)^(180-1)


### extract scenarios 
filelist <- unlist(lapply(res, function(folder) {
  files <- list.files(paste0(folder, "/"), pattern = ".gdx", full.names = TRUE)
  return(files)}))

filelist <- filelist[stringr::str_detect(filelist,"EXP")]
files <- c(paste0(res,"/",filelist)) 

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
  select(-X,-ID) %>% 
  as_tibble() %>% 
  mutate(term=as.integer(term) ) %>% 
  mutate_if(is.integer, as.character) %>% 
  mutate(theta=as.numeric(theta),rcp=str_remove(rcp,"RCP")) %>% 
  rename(pulse_time=pulse,cool_rate=cool,geo_end=term,geo_start=start,term=term_delta)
  
cat("Calculating net present cost: part1...\n")

damnpv_pre <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemasked") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemasked")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemasked")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemasked" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srmpulse)) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  select(-term) %>% 
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t<=as.numeric(term)) %>% 
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)), deltatemp=temp_srmpulse-temp ) %>% 
  mutate(impl = srm_masking * forctoUSD,
         dir = gwp * srm_masking * cost,
         dam = gwp*( alpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1+srm/forc) ) %>%
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,gas,delta,alpha,theta,term,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

cat("Calculating net present cost: part2...\n")

damnpv_post_noterm <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemasked") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemasked")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemasked")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemasked" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srmpulse)) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  select(-term) %>% 
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t>as.numeric(term)) %>% 
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)), deltatemp=temp_srmpulse-temp ) %>% 
  mutate(impl = srm_masking * forctoUSD,
         dir = gwp * srm_masking * cost,
         dam = gwp*( alpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10 / 3.71 * (1+srm/forc) ) %>%
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,gas,delta,alpha,theta,term,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

cat("Calculating net present cost: part3...\n")

damnpv_post_term <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemaskedterm") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemaskedterm")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemaskedterm")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemaskedterm" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,ecs,tcr,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srmpulse)) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  inner_join(id_montecarlo) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time) & t>as.numeric(term)) %>% 
  rowwise() %>% 
  mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha * (temp**2-temp_srm**2)),
         mask = 2 * gwp * alpha * temp_ghg *  (1-cos(theta * pi/180)) * srm_masking * as.numeric(ecs)/10/3.71 * (1+srm/forc)) %>% 
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,gas,term,delta,alpha,theta,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv + damnpv)

cat("Calculating net present cost: putting things together...\n")

pulse_size <- W_EMI %>% 
  filter(ghg==gas) %>% 
  filter(experiment %in% c("srm","srmpulse") ) %>% 
  group_by(ghg,t,rcp,ecs,tcr,gas,cool_rate,pulse_time,geo_start,geo_end) %>% 
  mutate(pulse_size=value-value[experiment=="srm"]) %>% 
  filter(experiment=="srmpulse" & pulse_size!=0) %>% 
  ungroup() %>% 
  select(-experiment,-ghg,-t,-value,-gdx,-file,-term) %>% 
  mutate(pulse_size=ifelse(gas=="co2",pulse_size*1e9,pulse_size*1e6))

scc <- TATM %>%  rename(temp_srm = value) %>% 
  filter(experiment=="srm") %>% 
  select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srm) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_start,geo_end,temp_srmpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  inner_join(id_montecarlo %>% select(-term,-prob,-theta)) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time)) %>% 
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(dam = gwp*( alpha * ((temp_srmpulse)**2-temp_srm**2)) ) %>% 
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,gas,delta,alpha) %>%
  summarise(damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>%
  full_join(pulse_size) %>% 
  mutate(scc=damnpv/pulse_size ) %>% 
  select(rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,gas,delta,alpha,scc)

damnpv <- damnpv_pre %>% select(gas,rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpre=costnpv) %>% 
  full_join(damnpv_post_noterm %>% select(gas,rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpostnoterm=costnpv)) %>% 
  full_join(damnpv_post_term %>% select(gas,rcp,ecs,tcr,cool_rate,pulse_time,geo_start,geo_end,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpostterm=costnpv)) %>% 
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
damnpv %>% 
  group_by(gas) %>% 
  summarise(med=median(npc_srm,na.rm=TRUE),sd=sd(npc_srm,na.rm=TRUE) )

damnorm <- damnpv %>% group_by(gas) %>% mutate(npc_norm=npc_srm/median(npc_srm,na.rm=TRUE))

ggplot(damnorm %>% filter(gas=="co2")) +
  geom_boxplot(aes(x=npc_norm,color=gas),outlier.shape = NA) +
  facet_wrap(gas~.,scales="free") +
  coord_cartesian(xlim = quantile(damnorm$npc_norm, c(0.03, 0.97), na.rm=TRUE))

library(gsaot)
gsoat_data <- damnpv %>% filter(gas=="ch4") %>% ungroup()
stat_analysis_ch4 <- ot_indices_1d(gsoat_data %>% 
                                     select(rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                   gsoat_data %>% pull(npc_srm ), 
                                   M= 15, 
                                   boot = T, 
                                   R = 100)
lowerbound_ch4 <- lower_bound(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

gsoat_data <- damnpv %>% filter(gas=="co2") %>% ungroup()
stat_analysis_co2 <- ot_indices_1d(gsoat_data %>% 
                                     select(rcp, ecs, tcr, cool_rate, pulse_time, geo_start, geo_end, delta, prob, alpha, theta, term),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15, 
                                   boot = T, 
                                   R = 100)
lowerbound_co2 <- lower_bound(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

ggplot(as_tibble(stat_analysis_ch4$indices_ci)) + 
  geom_bar(aes(x=input,
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_ch4$indices) +
  ylab("importance [ch4]") + xlab("")

ggplot(as_tibble(stat_analysis_co2$indices_ci)) + 
  geom_bar(aes(x=input,
               y=original,
               fill=input),stat="identity",color="black") +
  geom_hline(yintercept=lowerbound$indices) +
  ylab("importance [co2]") + xlab("") }
