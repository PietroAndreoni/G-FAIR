#### load data for figure 1
require(tidyverse)
require(gdxtools)

experiment <- "pulse"
ghg <- c("ch4","co2")
respath <- "../Paper and figures/"
tstart <- "historical_run"
rcp <- "RCP45"

## climate damage function parameters
alpha <- 0.595*1.25/100
powerdam <- 2
gwp <- 105*1e12
cost <- 0.001
powercost <- 1
eff <- 0.15

##
tsec <- gdxtools::batch_extract("save_delta",
                                files=c(paste0("../Results/",rcp,"_EXP",experiment,"_GAS",ghg,"_IC",tstart,".gdx"),
                                        paste0("../Results/",rcp,"_EXP",experiment,"masked","_GAS",ghg,"_IC",tstart,".gdx")))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
         file=str_remove_all(file,"../Results/|.gdx"),
         gas=str_extract(file,"(?<=GAS).+?(?=_)"),
         rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
         experiment=str_extract(file,"(?<=EXP).+?(?=_)"),
         initial_conditions=str_extract(file,"(?<=IC).*"),
         masking=ifelse(str_detect(file,"masked"),"yes","no")) %>%
  mutate(experiment=str_remove(experiment,"masking")) %>%
  ungroup() %>%
  complete(t) 

temp <- gdxtools::batch_extract("TATM",
                                files=unique(paste0("../Results/",rcp,"_EXPsimulation_IC",tstart,".gdx")))$TATM %>%
  as.data.frame() %>% 
  rename(file=gdx) %>%
  mutate(t=as.numeric(t),
         file=str_remove_all(file,"../Results/|.gdx"),
         rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
         experiment=str_extract(file,"(?<=EXP).+?(?=_)"),
         initial_conditions=str_extract(file,"(?<=IC).*")) %>%
  rename(temp = value)

temp0 <- (temp %>% 
            filter(t==2))$temp

# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
forctoTg <- 1/0.2

# from https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7/pdf  
TgtoUSD <- 2250*10^6

forctoUSD <- forctoTg * TgtoUSD # US$/(W/m^2)

gwpmax <- gwp* (1 + 0.022)^(80-1)

pulse_size <- tsec %>% filter(Variable=="emi" & masking=="yes") %>% 
  select(gas,value) %>% 
  rename(pulse_size=value) %>% 
  mutate(pulse_size=ifelse(gas=="co2",pulse_size*1e9,pulse_size*1e6))

totcosttime <- tsec %>% 
  filter(Variable=="forc" & ghg=="sai" & masking=="yes") %>%
  mutate(value=-value) %>%
  inner_join(tsec %>% 
               filter(Variable=="T" & ghg=="co2" & masking=="no") %>%
               rename(deltatemp = value) %>% 
               select(t,gas,deltatemp)) %>%
  inner_join(temp %>% 
               select(t,rcp,temp,initial_conditions)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=c(seq(0.002,0.06,by=0.002)))) %>%
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(impl = value*forctoUSD,
         dir = gwp*value*cost,
         dam = (gwp*alpha*((temp+deltatemp)**powerdam-temp**powerdam) ) ) 

damt <- totcosttime %>%
  group_by(file,gas,delta) %>%
  mutate(forcnorm=value/max(value)) %>%
  ungroup() %>% filter(forcnorm>0.001) %>%
  cross_join(data.frame(tt=seq(10,1000,by=10))) %>%
  filter(t>tt) %>%
  group_by(file,gas,delta,tt) %>%
  summarise(damnpv = sum( ((1-eff)*dam) / (1+delta)^t) ) %>%
  ungroup() %>% cross_join(data.frame(prob=seq(0,0.05,by=0.002))) %>%
  group_by(file,gas,delta,prob) %>%
  summarise(termnpv = sum(prob * (1-prob) ^ (tt/10-1) * damnpv )  ) 

damres <- totcosttime %>%
  group_by(file,gas,delta) %>%
  mutate(forcnorm=value/max(value)) %>%
  ungroup() %>% filter(forcnorm>0.001) %>%
  group_by(file,gas,delta) %>%
  summarise(damnpv = sum( (eff*dam) / (1+delta)^t) ) 

totcost <- totcosttime %>% 
  group_by(file,gas,delta) %>%
  summarise(implnpv = sum( (impl+dir)/(1+delta)^(t-1) ) ) %>% 
  ungroup() %>% 
  full_join(damt) %>%
  full_join(damres) %>%
  full_join(pulse_size) %>%
  mutate(costnpv=(termnpv+implnpv+damnpv)/pulse_size ) %>%
  group_by(file,gas,prob) %>%
  mutate(costnormd = costnpv/costnpv[delta==0.02]) %>%
  group_by(file,gas,delta) %>%
  mutate(costnormp = costnpv/costnpv[prob==0]) %>%
  group_by(file,gas) %>%
  mutate(costnorm = (costnpv)/costnpv[prob==0 & delta==0.02])



scc <- tsec %>% 
  filter(Variable=="T" & ghg=="co2" & masking=="no") %>%
  rename(deltatemp = value) %>% 
  select(t,gas,deltatemp) %>%
  inner_join(temp %>% 
               select(t,rcp,temp,initial_conditions)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=c(seq(0,0.009,by=0.001),seq(0.01,0.05,by=0.01)))) %>%
  cross_join(data.frame(sensdam= c(seq(0.5,2,by=0.5),5))) %>%
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  group_by(gas,delta,sensdam) %>%
  summarise(scc = sum( gwpt * ( (1- 1 / (1+alpha*sensdam*(temp+deltatemp)^powerdam)) - (1- 1 / (1+alpha*sensdam*(temp)^powerdam )) ) /(1+delta)^(t-1) ) ) 

GWP <- tsec %>% 
  filter(Variable=="forc" & masking=="no") %>%
  group_by(file,gas) %>%
  summarise(gwp=sum(value)) %>%
  ungroup() %>% mutate(gwp=gwp/gwp[gas=="co2"])

####### figure 1
temperature1 <- ggplot(tsec %>% 
                        filter(Variable=="T" & ghg=="co2" & masking=="no" & t<76 )) +
  geom_line(aes(x=2019+t,
                y=value,
                color=gas),linewidth=1) + 
  xlab('') +
  ylab('Temperature (K)') + ylim(c(0,0.025)) + 
  theme_classic() +
  theme(legend.position="none") 

temperature2 <- ggplot(tsec %>% 
                         filter(Variable=="T" & ghg=="co2" & masking=="no" & t>200 )) +
  geom_line(aes(x=2019+t,
                y=value,
                color=gas),linewidth=1) +  
  xlab('') + ylab('')+ ylim(c(0,0.025)) + 
  theme_classic() +
  theme(axis.line.y=element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y =  element_blank(), 
        legend.position = "none")

temperature <- ggpubr::ggarrange(temperature1,temperature2)

forcing1 <- ggplot(tsec %>% 
                    filter(Variable=="forc" & t<76 & masking=="yes" & ghg=="sai") %>%
                    group_by(t,gas,masking) %>%
                    summarise(value=sum(-value))) +
  geom_line(aes(x=2019+t,
                y=value,
                color=gas),
            linewidth=1)  + 
  xlab('') +
  ylab('Radiative forcing variation (W*m^-2)') + 
  ylim(c(-1e-8,0.065)) + 
  theme_classic() +
  theme(legend.position="none")

forcing2 <- ggplot(tsec %>% 
                     filter(Variable=="forc" & t>200 & masking=="yes" & ghg=="sai") %>%
                     group_by(t,gas,masking) %>%
                     summarise(value=sum(-value))) +
  geom_line(aes(x=2019+t,
                y=value,
                color=gas),
            linewidth=1)  + 
  xlab(' ') +
  ylab('')+   
  ylim(c(-1e-8,0.065)) + 
  theme_classic() +
  theme(axis.line.y=element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y =  element_blank())

forcing <- ggpubr::ggarrange(forcing1+theme(legend.position = "none"),forcing2+theme(legend.position = "none"))

fig1 <- ggpubr::ggarrange(temperature,forcing)
ggsave("figure_1.svg",width=18,height=9,path=respath,plot=fig1)

fig2a <- ggplot(totcost %>% 
                  filter(delta > 0.01) ) +
  geom_line(data=.%>% filter(prob==0),aes(x=delta*100,
                y=costnormd,
                color=gas ),
            linewidth=1) +
  geom_ribbon(data=.%>% 
              group_by(delta,gas) %>% 
              summarise(max=max(costnormd),min=min(costnormd)),
            aes(x=delta*100,
                ymin=min,
                ymax=max,
                fill=gas,
                color=gas ),
            linewidth=0.1,alpha=0.2) +
  ylab('Relative cost') +
  xlab('Discount rate (%)') + 
  theme_classic() +
  theme(legend.position = "bottom") +
  ylim(c(0,3))


fig2b <- ggplot(totcost %>% filter(delta>=0.01)) +
  geom_line(data=.%>% filter(delta==0.02),
            aes(x=prob*100,
                y=costnormp,
                color=gas), 
            linewidth=1) +
  geom_ribbon(data=.%>% 
                group_by(gas,prob) %>% 
                summarise(max=max(costnormp),min=min(costnormp)),
              aes(x=prob*100,
                  ymin=min,
                  ymax=max,
                  fill=gas,
                  color=gas ),
              linewidth=0.1,alpha=0.2) +
  ylab('Relative cost') +
  xlab('Probability of termination (%/decade)') + 
  theme_classic() +
  theme(legend.position = "bottom") +
  ylim(c(0,3))

fig2 <- ggpubr::ggarrange(fig2a,fig2b,nrow=1,common.legend = TRUE)
ggsave("figure_2.svg",width=11,height=5.5,path=respath,plot=fig2)

totcost %>% filter(delta == 0.02 & prob==0)

totcost %>% filter(delta == 0.01 & prob==0) 

totcost %>% filter(delta == 0.05 & prob==0) 

totcost %>% 
  filter(delta >= 0.01) %>% 
  ggplot() + 
  #geom_raster(aes(x=prob,y=delta,fill=abs(log(costnorm) ) ) ) + 
  geom_contour_filled(aes(x=prob*100,y=delta*100,z=abs(log(costnorm) ) ) )  +
  geom_contour(aes(x=prob*100,y=delta*100,z=abs(log(costnorm) ) ) )  +
  facet_wrap(gas~.,scales="free") +
  scale_fill_viridis_d() + scale_color_viridis_c() + 
  ylab('Discount rate (%)') +
  xlab('Probability of termination (%/decade)') + theme_classic()

totcost %>% 
  filter(delta >= 0.01 & gas=="co2") %>% 
  ggplot() + 
  #geom_raster(aes(x=prob,y=delta,fill=abs(log(costnorm) ) ) ) + 
  geom_contour_filled(aes(x=prob*100,y=delta*100,z=costnpv) )  +
  geom_contour(aes(x=prob*100,y=delta*100,z=costnpv) )  +
  facet_wrap(gas~.,scales="free") +
  scale_fill_viridis_d() + scale_color_viridis_c() + 
  ylab('Discount rate (%)') +
  xlab('Probability of termination (%/decade)') + theme_classic()

totcost %>% 
  filter(delta >= 0.01 & gas=="ch4") %>% 
  ggplot() + 
  #geom_raster(aes(x=prob,y=delta,fill=abs(log(costnorm) ) ) ) + 
  geom_contour_filled(aes(x=prob*100,y=delta*100,z=costnpv) )  +
  geom_contour(aes(x=prob*100,y=delta*100,z=costnpv) )  +
  facet_wrap(gas~.,scales="free") +
  scale_fill_viridis_d() + scale_color_viridis_c() + 
  ylab('Discount rate (%)') +
  xlab('Probability of termination (%/decade)') + theme_classic()
