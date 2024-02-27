#### load data for figure 1
require(tidyverse)
require(gdxtools)

experiment <- "pulse"
ghg <- c("n2o","ch4","co2")
respath <- "../Papers and figures/"
tstart <- "historical_run"
rcp <- "RCP45"
alpha <- 0.595*1.25/100
power <- 2
eff <- 1
gwp <- 105*1e12
tsec <- gdxtools::batch_extract("save_delta",
                        files=paste0("../Results/",rcp,"_EXP",experiment,"_GAS",ghg,"_IC",tstart,".gdx"))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
         file=str_remove_all(file,"../Results/|.gdx"),
         gas=str_extract(file,"(?<=GAS).+?(?=_)"),
         rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
         experiment=str_extract(file,"(?<=EXP).+?(?=_)"),
         initial_conditions=str_extract(file,"(?<=IC).*") ) 

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

totcosttime <- tsec %>% 
  filter(Variable=="forc") %>%
  group_by(t,file,gas) %>%
  summarise(value=sum(value)) %>%
  full_join(tsec %>% 
              filter(Variable=="T" & ghg=="co2") %>%
              rename(deltatemp = value) %>% 
              select(t,gas,deltatemp)) %>%
  inner_join(temp %>% 
              select(t,rcp,temp,initial_conditions)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=seq(0,0.1,by=0.01))) %>%
  cross_join(data.frame(senscost= c(0.1,0.5,1,2,10,100))) %>%
  mutate(dir = value*forctoUSD*senscost,
         dam = ( gwp*(1- 1 / (1+(1-eff)*alpha*(temp+deltatemp)^power)) -  gwp*(1- 1 / (1+(1-eff)*alpha*(temp)^power )) )  ) 
  
totcost <- totcosttime %>%
  group_by(file,delta,gas,senscost) %>%
  summarise(dirnpv = sum( dir/(1+delta)^(t-1) ),
            damnpv = sum(  dam / (1+delta)^(t-1) ) ) %>% 
  ungroup() %>% mutate(cost=dirnpv+damnpv) %>%
  group_by(file,gas,senscost) %>%
  mutate(costnorm = cost/cost[delta==0.01]) 


GWP <- tsec %>% 
  filter(Variable=="forc") %>%
  group_by(file,gas) %>%
  summarise(gwp=sum(value)) %>%
  ungroup() %>% mutate(gwp=gwp/gwp[gas=="co2"])

####### figure 1
temperature <- ggplot(tsec %>% 
       filter(Variable=="T" & ghg=="co2") %>% 
         mutate(value=case_when(gas=="ch4"~value/10,
                                gas=="n2o"~value/100,
                                gas=="co2"~value) )) +
  geom_line(aes(x=t,
                y=value/1000,
                color=gas),linewidth=1) + 
  geom_text(data=data.frame(gas=c("ch4","n2o","co2"),
                            text=c("x10","x100"," "),
                            value=c(2e-10,3e-9,0)),
            aes(x=250,
                y=value/1000,
                label=text,
                color=gas)) +
  xlab('time [years]') +
  ylab('Temperature variation [K]')

forcing <- ggplot(tsec %>% 
         filter(Variable=="forc") %>%
         group_by(t,file,gas) %>%
         summarise(value=sum(value)) %>%
         mutate(value=case_when(gas=="ch4"~value/10,
                                gas=="n2o"~value/100,
                                gas=="co2"~value) ) ) +
  geom_line(aes(x=t,
                y=value/1000,
                color=gas),linewidth=1)  + 
  geom_text(data=data.frame(gas=c("ch4","n2o","co2"),
                            text=c("x10","x100"," "),
                            value=c(2e-10,3e-9,0)),
            aes(x=250,
                y=value/1000,
                label=text,
                color=gas)) +
  xlab('time [years]') +
  ylab('Radiative forcing variation [W*m^-2]')

fig1 <- ggpubr::ggarrange(temperature,forcing,common.legend=TRUE)
ggsave("figure_1.png",width=10,height=6,path=respath)


#### figure 2
absolute <- ggplot(totcost %>% 
                     filter( delta<=0.05 & delta>=0.01 & senscost == 1) %>%
                     mutate(cost=case_when(gas=="n2o"~cost/126,
                                           gas=="ch4"~cost/3.28,
                                           .default = cost) ) ) +
  geom_line(aes(x=delta*100,
                y=cost/1000,
                color=gas),
            linewidth=1) + 
  geom_text(data=data.frame(gas=c("ch4","n2o","co2"),
                            text=c("x3.28","x126"," "),
                            value=c(3.8,1.7,0),
                            intercept=c(2.5,4,0)),
            aes(x=intercept,
                y=value,
                label=text,
                color=gas)) +
  ylab('Net Present Cost, CO2 and CH4 [US$]') +
  xlab('Discount rate [%]') + xlim(c(1,5))

ggplot(totcost %>% 
         filter( delta %in% c(0,0.02,0.04) & senscost %in% c(0.1,1,10) ) %>%
         mutate(cost=case_when(gas=="n2o"~cost,
                               gas=="ch4"~cost,
                               .default = cost) ) ) +
  geom_bar(aes(x=as.factor(delta),
                y=cost/1000,
                fill=gas),
            linewidth=1, 
           stat="identity",position="dodge") + 
  ylab('Net Present Cost [US$/ton]') +
  xlab('Discount rate [%]') +
  facet_wrap(gas~.,scales="free")

relative <- ggplot(totcost %>% 
                     filter(delta<=0.05 & senscost == 1)) +
  geom_line(aes(x=delta*100,
                y=costnorm,
                color=gas),linewidth=1) +
  ylab('Relative cost [frac]') +
  xlab('Discount rate [%]') +
  ylim(c(0,1)) + xlim(c(1,5))

fig2 <- ggpubr::ggarrange(absolute,relative,common.legend=TRUE)
ggsave("figure_2.png",width=10,height=6,path=respath)
