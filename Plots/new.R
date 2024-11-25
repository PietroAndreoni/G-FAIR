#### load data for figure 1
require(tidyverse)
require(gdxtools)

experiment <- "pulse"
ghg <- c("ch4","co2")
respath <- "../Paper and figures/"
tstart <- "historical_run"
rcp <- "RCP45"

## climate damage function parameters
powerdam <- 2
gwp <- 105*1e12
cost <- 0.001
powercost <- 1
pulse_size <- 1e6
climsens <- 0.84 # K/(W/M^2)
eff <- 0.85

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
  cross_join(data.frame(delta=c(seq(0,0.009,by=0.001),seq(0.01,0.07,by=0.01)))) %>%
  cross_join(data.frame(alpha= c(0.0074,0.0023))) %>%
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(impl = value*forctoUSD,
         dir = gwp*value*cost,
         dam =  gwp*(1-eff)*alpha*((climsens*value)^2 + 2*deltatemp*climsens*value) ) 

totcost <- totcosttime %>%
  group_by(file,gas,delta,alpha) %>%
  summarise(implnpv = sum( impl/(1+delta)^(t-1) ),
            dirnpv = sum( dir/(1+delta)^(t-1) ),
            damnpv = sum( dam/(1+delta)^(t-1) ) ) %>% 
  ungroup() %>% 
  mutate(cost=(damnpv+implnpv+dirnpv)/1e6 ) %>%
  group_by(file,gas,alpha) %>%
  mutate(costnorm = cost/cost[delta==0.02]) 

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
                        filter(Variable=="T" & ghg=="co2" & masking=="no" & t<76 ) %>% 
                        mutate(value=case_when(gas=="ch4"~value*800,
                                               gas=="n2o"~value,
                                               gas=="co2"~value*37200) )) +
  geom_line(aes(x=2024+t,
                y=value,
                color=gas),linewidth=1) + 
  xlab('time [years]') +
  ylab('Temperature variation [K]') + ylim(c(0,0.03)) + theme(legend.position="none")

temperature2 <- ggplot(tsec %>% 
                         filter(Variable=="T" & ghg=="co2" & masking=="no" & t>200 ) %>% 
                         mutate(value=case_when(gas=="ch4"~value*800,
                                                gas=="n2o"~value,
                                                gas=="co2"~value*37200) )) +
  geom_line(aes(x=2024+t,
                y=value,
                color=gas),linewidth=1) + 
  xlab('') + ylab('')+ ylim(c(0,0.03)) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y =  element_blank(), legend.position = "none")

temperature <- ggpubr::ggarrange(temperature1,temperature2)

forcing1 <- ggplot(tsec %>% 
                    filter(Variable=="forc" & t<76 & masking=="yes" & ghg=="sai") %>%
                    group_by(t,gas,masking) %>%
                    summarise(value=sum(-value)) %>% 
                    mutate(value=case_when(gas=="ch4"~value*800,
                                           gas=="n2o"~value,
                                           gas=="co2"~value*37200) ) ) +
  geom_line(aes(x=2024+t,
                y=value,
                color=gas),
            linewidth=1)  + 
  xlab('time [years]') +
  ylab('Radiative forcing variation [W*m^-2]') + ylim(c(0,0.13)) + 
  theme(legend.position="none",
        plot.margin = unit(c(0,0,0,0), 'lines'))

forcing2 <- ggplot(tsec %>% 
                     filter(Variable=="forc" & t>200 & masking=="yes" & ghg=="sai") %>%
                     group_by(t,gas,masking) %>%
                     summarise(value=sum(-value)) %>% 
                     mutate(value=case_when(gas=="ch4"~value*800,
                                            gas=="n2o"~value,
                                            gas=="co2"~value*37200) ) ) +
  geom_line(aes(x=2024+t,
                y=value,
                color=gas),
            linewidth=1)  + 
  xlab(' ') +
  ylab('')+ ylim(c(0,0.13))  + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y =  element_blank(),
        plot.margin = unit(c(0,0,0,0), 'lines'))

forcing <- ggpubr::ggarrange(forcing1,forcing2+theme(legend.position = "none"))

fig1 <- ggpubr::ggarrange(temperature,forcing)
ggsave("figure_1.png",width=12,height=6,path=respath,plot=fig1)

fig2 <- ggplot() +
  geom_line(data=totcost %>% 
              filter(delta > 0.002 & alpha==0.0074),
            aes(x=delta*100,
                y=cost,
                color=gas),
            linewidth=1) +
  ylab('Relative cost [frac]') +
  xlab('Discount rate [%]') + theme(legend.position = "bottom") + facet_wrap(gas~.,scales="free")
ggsave("figure_2.png",width=10,height=10,path=respath,plot=fig2)


##### social cost of carbon 
ggplot() +
  geom_line(data=scc %>% 
              filter(sensdam==1 & delta>=0.01),
            aes(x=delta*100,
                y=scc/pulse_size,
                color=gas),
            linewidth=1) +  facet_wrap(.~gas,scales="free_y")
