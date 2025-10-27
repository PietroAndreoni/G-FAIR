eff <- function(g,c,theta) { pmin(1,(1 + (g/c)**2 - g/c * 2 * cos(theta* pi/180)) ) }

sanitize <- function(.x) {
  .x %>% 
    as_tibble() %>% 
    mutate(t=as.numeric(t),
           file=str_remove_all(gdx,"../Results/|.gdx"),
           gas=str_extract(file,"(?<=GAS).+?(?=_)"),
           rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
           cool_rate=str_extract(str_remove(file,"RCP"),"(?<=RC).+?(?=_)"),
           pulse_time=str_extract(file,"(?<=PT).+?(?=_)"),
           geo_end=str_extract(file,"(?<=EC).*"),
           term=str_extract(file,"(?<=TT).+?(?=_)"),
           experiment=str_extract(file,"(?<=EXP).+?(?=_)")) %>%
    ungroup() %>%
    filter(t<=480) %>% 
    complete(t) }


### extract scenarios with no probability of termination
filelist <- list.files("../Results/",pattern=".gdx")
filelist <- filelist[stringr::str_detect(filelist,"EXP")]
filelist <- filelist[!stringr::str_detect(filelist,"EXPsrmpulsemaskedterm")]

TATM <- gdxtools::batch_extract("TATM", files=c(paste0("../Results/",filelist)) )$TATM %>% 
  sanitize()       

W_EMI <- gdxtools::batch_extract("W_EMI", files=c(paste0("../Results/",filelist)) )$W_EMI %>% 
  sanitize()    

FORC <- gdxtools::batch_extract("FORCING", files=c(paste0("../Results/",filelist)) )$FORCING %>% 
  sanitize()    

SRM <- gdxtools::batch_extract("SRM", files=c(paste0("../Results/",filelist)) )$SRM  %>% 
  sanitize() 

background_srm <- gdxtools::batch_extract("forcing_srm", files=c(paste0("../Results/",filelist)) )$forcing_srm %>% 
  sanitize()

tot_forcing <- FORC %>% 
  group_by(t,gas,rcp,cool_rate,pulse_time,geo_end,experiment) %>% 
  summarise(value=sum(value)) %>% ungroup()

tot_forcing %>% 
  filter(experiment=="srm") %>% 
  select(t,cool_rate,geo_end,rcp,value) %>% 
  unique() %>% 
  full_join(background_srm %>%
              filter(experiment=="srm") %>%  
              select(t,rcp,cool_rate,geo_end,value) %>%  
              unique() %>% 
              rename(srm=value)) %>% 
  filter(t<=480) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  unique() %>% 
  ggplot() +
  geom_line(aes(x=t+2019,y=value-srm,color=rcp,group=interaction(cool_rate,geo_end,rcp)))

background_srm %>% select(t,cool_rate,geo_end,value) %>% 
  unique() %>% 
  filter(t<=480) %>% 
  ggplot() +
  geom_line(aes(x=t+2019,y=value,color=cool_rate,linetype=geo_end))

cross_join(data.frame(g=seq(0,1,by=0.1)),data.frame(theta=seq(0,90,by=10))) %>% mutate(eff=eff(g,1,theta)) %>% 
  ggplot() +
  geom_line(aes(x=g,y=eff,color=as.factor(theta))) + scale_color_viridis_d()

totcosttime <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemasked") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemasked")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemasked")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemasked" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,pulse_time,cool_rate,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)), deltatemp=temp_srmpulse-temp ) 



#### extract scenarios with termination
filelist <- list.files("../Results/",pattern=".gdx")
filelist <- filelist[stringr::str_detect(filelist,"EXP")]
filelist <- filelist[stringr::str_detect(filelist,"RCP45")]
filelist <- filelist[stringr::str_detect(filelist,"ECna|EXPbase|EXPpulse")]
filelist <- filelist[stringr::str_detect(filelist,"PT2_|EXPbase")]
filelist <- filelist[!stringr::str_detect(filelist,"GASn2o")]

TATM <- gdxtools::batch_extract("TATM", files=c(paste0("../Results/",filelist)) )$TATM %>% 
  sanitize() 

SRM <- gdxtools::batch_extract("SRM", files=c(paste0("../Results/",filelist)) )$SRM %>% 
  sanitize() 

dtnoterm <- TATM %>% filter(experiment %in% c("srmpulse","srm")) %>% 
  group_by(t,gas) %>% 
  mutate(dt=value[experiment=="srmpulse"]-value[experiment=="srm"]) %>% 
  ungroup() %>% 
  filter(experiment=="srmpulse") %>% 
  select(-experiment,-file,-gdx,-term) 

FORC <- gdxtools::batch_extract("FORCING", files=c(paste0("../Results/",filelist)) )$FORCING %>% 
  sanitize()    

background_srm <- gdxtools::batch_extract("forcing_srm", files=c(paste0("../Results/",filelist)) )$forcing_srm %>% 
  sanitize()

tot_forcing <- FORC %>% 
  group_by(t,gas,rcp,cool_rate,pulse_time,geo_end,term,experiment) %>% 
  summarise(value=sum(value)) %>% ungroup()

TATM %>% 
  filter(experiment!="base") %>% 
  group_by(t,gas) %>% 
  mutate(dt=value-value[experiment=="srm"]) %>% 
  ungroup() %>% 
  filter(experiment=="srmpulsemaskedterm") %>% 
  select(-experiment,-file,-gdx) %>% 
  bind_rows(dtnoterm) %>% 
  ggplot() +
  geom_line(aes(x=t,y=dt,color=term,group=term)) +
  facet_wrap(gas~.,scales="free") + xlim(c(0,200))

cross_join(data.frame(t=seq(1,500,by=1)),
           data.frame(delta=c(seq(0.01,0.06,by=0.01) ) ) ) %>% 
  mutate(rr=1/(1+delta)^t) %>% 
  ggplot(aes(x=t,y=rr,color=delta,group=delta)) + geom_line() + scale_color_viridis_c()

cross_join(data.frame(t=seq(1,100,by=1)),
           data.frame(tt=seq(10,100,by=10)) ) %>% 
  filter(t>tt) %>% 
  cross_join( data.frame(p=seq(0,0.1,by=0.01) ) ) %>% 
  group_by(t,p) %>% 
  mutate(rr=sum(p * (1-p) ^ (tt/10-1)) ) %>% 
  ggplot(aes(x=t,y=rr,color=p,group=p)) + geom_line() + scale_color_viridis_c()

totcosttime_term <- tot_forcing %>% rename(forc=value) %>% filter(experiment=="srmpulsemaskedterm") %>% 
  full_join(TATM %>%  rename(temp = value) %>% filter(experiment=="srmpulsemaskedterm")) %>%
  full_join(background_srm %>% rename(srm=value) %>% filter(experiment=="srmpulsemaskedterm")) %>% 
  full_join(SRM %>% rename(srm_masking=value) %>% filter(experiment=="srmpulsemaskedterm" & srm_masking>=0)) %>% 
  select(-file,-gdx,-experiment) %>% 
  inner_join(TATM %>%  rename(temp_ghg = value) %>% 
               filter(experiment=="base") %>% 
               select(t,rcp,temp_ghg)) %>%
  inner_join(TATM %>%  rename(temp_srm = value) %>% 
               filter(experiment=="srm") %>% 
               select(t,gas,rcp,pulse_time,cool_rate,geo_end,temp_srm)) %>%
  inner_join(TATM %>%  rename(temp_srmpulse = value) %>% 
               filter(experiment=="srmpulse") %>% 
               select(t,gas,rcp,pulse_time,cool_rate,geo_end,temp_srmpulse)) %>%
  inner_join(TATM %>%  rename(temp_ghgpulse = value) %>% 
               filter(experiment=="pulse") %>% 
               select(t,gas,rcp,pulse_time,temp_ghgpulse)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) 
