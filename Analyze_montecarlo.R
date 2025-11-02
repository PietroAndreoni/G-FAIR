library(gdxtools)
library(dplyr)
library(stringr)

extract_names <- function(.x,res) {
  .x %>% 
    as_tibble() %>% 
    mutate(file=str_remove_all(gdx,paste0(res,"/|.gdx")),
           path=res,
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
  Analyze_montecarlo.R [-o <res>] 

Options:
-o <res>              Path where the results are (default: Results_montecarlo)
' -> doc

library(docopt)
opts <- docopt(doc, version = 'Montecarlo')

res = ifelse(is.null(opts[["o"]]), "Results_montecarlo_mini", as.character(opts["o"]) )

### extract scenarios with no probability of termination
filelist <- list.files(paste0(res,"/"),pattern=".gdx")
filelist <- filelist[stringr::str_detect(filelist,"EXP")]
files <- c(paste0(res,"/",filelist)) 

cat("Loading the data...")

TATM <- gdxtools::batch_extract("TATM", files=files)$TATM    

sanitized_names <- TATM %>% filter(t==2) %>% 
  extract_names(.,res) %>% select(-t,-value) %>% unique()

TATM <- TATM  %>%
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
  group_by(t,gas,rcp,ecs,tcr,cool_rate,pulse_time,geo_end,term,experiment) %>% 
  summarise(value=sum(value)) %>% ungroup() 


cat("Loading the IDs from Montecarlo id file...")

id_montecarlo <-  as.data.frame(read.csv(paste0(res,"/id_montecarlo.csv")) ) %>% 
  select(-X,-ID) %>% 
  as_tibble() %>% 
  mutate(term=as.integer(term+100) ) %>% 
  mutate_if(is.integer, as.character) %>% 
  mutate(theta=as.numeric(theta),rcp=str_remove(rcp,"RCP")) %>% 
  rename(pulse_time=pulse,cool_rate=cool,geo_end=term,geo_start=start,term=term_delta)
  
cat("Calculating net present cost: part1...")

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
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
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
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,delta,alpha,theta,term,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

cat("Calculating net present cost: part2...")

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
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
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
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,delta,alpha,theta,term,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

cat("Calculating net present cost: part3...")

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
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
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
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,term,delta,alpha,theta,prob) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv + damnpv)

cat("Calculating net present cost: putting things together...")

scc <- TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srm) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,ecs,tcr,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  inner_join(id_montecarlo %>% select(-term,-prob,-theta)) %>% 
  ungroup() %>% 
  filter(t>=as.numeric(pulse_time)) %>% 
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(dam = gwp*( alpha * ((temp_srmpulse)**2-temp_srm**2)) ) %>% 
  group_by(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,delta,alpha) %>%
  summarise(damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE) ) %>%
  full_join(pulse_size) %>% 
  mutate(scc=damnpv/pulse_size ) %>% 
  select(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,delta,alpha,scc)

pulse_size <- W_EMI %>% 
  filter(ghg==gas) %>% 
  filter(experiment %in% c("srm","srmpulse") ) %>% 
  group_by(ghg,t,rcp,ecs,tcr,gas,cool_rate,pulse_time,geo_end) %>% 
  mutate(pulse_size=value-value[experiment=="srm"]) %>% 
  filter(experiment=="srmpulse" & pulse_size!=0) %>% 
  ungroup() %>% 
  select(-experiment,-ghg,-t,-value,-gdx,-file,-path,-term) %>% 
  mutate(pulse_size=ifelse(gas=="co2",pulse_size*1e9,pulse_size*1e6))

damnpv <- damnpv_pre %>% select(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpre=costnpv) %>% 
  full_join(damnpv_post_noterm %>% select(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpostnoterm=costnpv)) %>% 
  full_join(damnpv_post_term %>% select(rcp,ecs,tcr,cool_rate,pulse_time,geo_end,gas,term,delta,alpha,theta,prob,costnpv) %>%  rename(costpostterm=costnpv)) %>% 
  mutate(costnpv=costpre+as.numeric(prob)*costpostterm+(1-as.numeric(prob))*costpostnoterm) %>% 
  select(-costpre,-costpostterm,-costpostnoterm) %>% 
  full_join(pulse_size) %>% 
  full_join(scc) %>% 
  mutate(npc_srm=costnpv/pulse_size)

cat("Saving output...")

write.csv(damnpv,file=paste0(res,"/output_analysis.csv"))
