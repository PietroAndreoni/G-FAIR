#### load data for figure 1
require(tidyverse)
require(gdxtools)

experiment <- "pulse"
ghg <- c("n2o","ch4","co2")
respath <- "../Paper and figures/"
tstart <- "historical_run"
rcp <- "RCP45"

## climate damage function parameters
alpha <- 0.595*1.25/100
powerdam <- 2
gwp <- 105*1e12
cost <- 0.001
powercost <- 1
pulse_size <- 1e6

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
  cross_join(data.frame(delta=c(seq(0,0.009,by=0.001),seq(0.01,0.05,by=0.01)))) %>%
  cross_join(data.frame(senscost= c(seq(0,2,by=0.5),5))) %>%
  cross_join(data.frame(eff= seq(0.75,0.95,by=0.1))) %>%
  cross_join(data.frame(sensdam= c(seq(0.5,2,by=0.5),5))) %>%
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(impl = value*forctoUSD,
         dir = gwp*value*senscost*cost^powercost,
         dam = ( gwpt*(1- 1 / (1+(1-eff)*alpha*sensdam*(temp+deltatemp)^powerdam)) -  gwpt*(1- 1 / (1+(1-eff)*alpha*sensdam*(temp)^powerdam )) )  ) 
  
totcost <- totcosttime %>%
  group_by(file,gas,delta,senscost,eff,sensdam) %>%
  summarise(implnpv = sum( impl/(1+delta)^(t-1) ),
            dirnpv = sum( dir/(1+delta)^(t-1) ),
            damnpv = sum(  dam / (1+delta)^(t-1) ) ) %>% 
  ungroup() %>% 
  mutate(cost=dirnpv+damnpv+implnpv) %>%
  group_by(file,gas,senscost,eff,sensdam) %>%
  mutate(costnorm = cost/cost[delta==0.01]) 

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
temperature <- ggplot(tsec %>% 
       filter(Variable=="T" & ghg=="co2" & masking=="no") %>% 
         mutate(value=case_when(gas=="ch4"~value*10,
                                gas=="n2o"~value,
                                gas=="co2"~value*100) )) +
  geom_line(aes(x=t,
                y=value/1000,
                color=gas),linewidth=1) + 
  # geom_text(data=data.frame(gas=c("ch4","n2o","co2"),
  #                           text=c("x10","x100"," "),
  #                           value=c(2e-10,3e-9,0)),
  #           aes(x=250,
  #               y=value,
  #               label=text,
  #               color=gas)) +
  xlab('time [years]') +
  ylab('Temperature variation [K]')

forcing <- ggplot(tsec %>% 
         filter(Variable=="forc" & t<1000 & masking=="yes" & ghg=="sai") %>%
         group_by(t,gas,masking) %>%
         summarise(value=sum(-value)) %>% 
         bind_rows(tsec %>% 
                     filter(Variable=="forc" & t<1000 & masking=="no" & ghg!="sai") %>%
                     group_by(t,gas,masking) %>%
                     summarise(value=sum(value)))  %>%
         mutate(value=case_when(gas=="ch4"~value*10,
                                gas=="n2o"~value,
                                gas=="co2"~value*100) ) ) +
  geom_line(aes(x=t,
                y=value,
                color=gas,
                linetype=masking),linewidth=1)  + 
  # geom_text(data=data.frame(gas=c("ch4","n2o","co2"),
  #                           text=c("x10","x100"," "),
  #                           value=c(2e-10,3e-9,0)),
  #           aes(x=250,
  #               y=value/1000,
  #               label=text,
  #               color=gas)) +
  xlab('time [years]') +
  ylab('Radiative forcing variation [W*m^-2]')

fig1 <- ggpubr::ggarrange(temperature,forcing,common.legend=TRUE)
ggsave("figure_1.png",width=10,height=6,path=respath,plot=fig1)


#### figure 2
absolute1 <- ggplot(totcost %>% 
                     filter( delta >= 0.001 & 
                               delta <= 0.01 &
                              senscost==1 & 
                              eff==0.85 & 
                              sensdam==1 ) %>%
                     pivot_longer(c(implnpv,dirnpv,damnpv), 
                        names_to="Cost component") ) +
  geom_bar(aes(x=as.factor(delta*100),
               y=value/pulse_size,
               fill=`Cost component`),
           linewidth=1,
           color="black",
           stat="identity",
           position="stack") +
  geom_text(aes(x=as.factor(delta*100),
                y=value/pulse_size,
                color=`Cost component`,
                label=paste0(round(value/pulse_size,1),"$/ton")),
            vjust = -1,
            size = 2,
            stat="identity",
            position="stack") +
  ylab('Net Present Cost [US$/ton]') +
  xlab('Discount rate [%]') +
  scale_fill_manual(values=c("blue","grey","lightblue"),labels=c("implnpv"="Deployment",
                                                                 "dirnpv"="Side-effects",
                                                                 "damnpv"="Imperfect masking")) +
  scale_color_manual(values=c("blue","grey","lightblue"),labels=c("implnpv"="Deployment",
                                                                  "dirnpv"="Side-effects",
                                                                  "damnpv"="Imperfect masking")) +
  facet_wrap(.~gas,scales="free_y") + 
  theme(legend.position = "top")
ggsave("figure_2bis.png",width=10,height=5,path=respath)

absolute2 <- ggplot(totcost %>% 
                     filter(delta >= 0.01 &
                              senscost==1 & 
                              eff==0.85 & 
                              sensdam==1 ) %>%
                     pivot_longer(c(implnpv,dirnpv,damnpv), 
                                  names_to="Cost component") ) +
  geom_bar(aes(x=as.factor(delta*100),
               y=value/pulse_size,
               fill=`Cost component`),
           linewidth=1,
           color="black",
           stat="identity",
           position="stack") +
  geom_text(aes(x=as.factor(delta*100),
               y=value/pulse_size,
               color=`Cost component`,
               label=paste0(round(value/pulse_size,1)," $/ton")),
           vjust = -1,
           stat="identity",
           position="stack") +
  ylab('Net Present Cost [US$/ton]') +
  xlab('Discount rate [%]') +
  scale_fill_manual(values=c("blue","grey","lightblue"),labels=c("implnpv"="Deployment",
                                                                 "dirnpv"="Side-effects",
                                                                 "damnpv"="Imperfect masking")) +
  scale_color_manual(values=c("blue","grey","lightblue"),labels=c("implnpv"="Deployment",
                                                                 "dirnpv"="Side-effects",
                                                                 "damnpv"="Imperfect masking")) +
  facet_wrap(.~gas,scales="free_y") + 
  theme(legend.position = "top")

relative <- ggplot() +
  geom_line(data=totcost %>% 
              filter(delta > 0.005 & 
                       eff==0.85 & 
                       senscost==1 & 
                       sensdam==1),
            aes(x=delta*100,
                y=costnorm,
                color=gas),
            linewidth=1) +
  geom_ribbon(data=totcost %>% 
              filter(delta %in% c(0.01,0.02,0.03,0.04,0.05)) %>%
                group_by(delta,gas) %>%
                summarise(max=max(costnorm),min=min(costnorm)),
            aes(x=delta*100,
                ymin=min,
                ymax=max,
                color=gas,
                fill=gas),
            linewidth=0.5,
            alpha=0.2) +
  ylab('Relative cost [frac]') +
  xlab('Discount rate [%]') + theme(legend.position = "bottom")
void <- ggplot() + theme_void()
fig2 <- ggpubr::ggarrange(absolute2,ggpubr::ggarrange(void,relative,void,nrow=1,widths=c(0.2,0.6,0.2)),nrow=2,heights = c(0.6,0.4))
ggsave("figure_2.png",width=10,height=10,path=respath)


##### social cost of carbon 
ggplot() +
  geom_line(data=scc %>% 
              filter(sensdam==1),
            aes(x=delta*100,
                y=scc/pulse_size,
                color=gas),
            linewidth=1) +  facet_wrap(.~gas,scales="free_y")

