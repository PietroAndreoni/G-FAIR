#this script requires pulsemasked.R to be run before
ghg <- c("co2")
tremoval <- seq(20,300,by=20)

tsecrem <- gdxtools::batch_extract("save_delta",
                                files=unique(c(paste0("../Results/",rcp,"_EXP",experiment,"_REM",tremoval,"_GAS",ghg,"_IC",tstart,".gdx"),
                                               paste0("../Results/",rcp,"_EXP",experiment,"masked","_REM",tremoval,"_GAS",ghg,"_IC",tstart,".gdx"))))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
    file=str_remove_all(file,"../Results/|.gdx"),
    gas=str_extract(file,"(?<=GAS).+?(?=_)"),
    rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
    experiment=str_extract(file,"(?<=EXP).+?(?=_)"),
    tremoval=as.numeric(str_extract(file,"(?<=REM).+?(?=_)")),
    initial_conditions=str_extract(file,"(?<=IC).*"),
    masking=ifelse(str_detect(file,"masked"),"yes","no")) %>%
  mutate(experiment=str_remove(experiment,"masking")) %>%
  ungroup() %>%
  complete(t) 



temperature <- ggplot(tsecrem %>% 
                        filter(Variable=="T" & 
                                 masking=="no" & 
                                 ghg=="co2" & 
                                 tremoval %in% c(20,60,120,200) )) +
  geom_line(aes(x=t,
                y=value/pulse_size,
                color=as.factor(tremoval)), 
            linewidth=1) +
  geom_line(data=tsec %>% 
              filter(Variable=="T" & 
                       masking=="no" & 
                       ghg=="co2" & 
                       gas=="co2" ),
            aes(x=t,
                y=value/pulse_size ), 
            linewidth=1, 
            linetype=2)  + 
  scale_color_viridis_d() +
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Temperature variation [K]')

forcing <- ggplot(tsecrem %>%
                    filter(Variable=="forc" & t<1000 & masking=="yes" & ghg=="sai" & 
                             tremoval %in% c(20,60,120,200)) %>%
                    group_by(t,gas,masking,tremoval) %>%
                    summarise(value=sum(-value)) %>% 
                    bind_rows(tsecrem %>% 
                                filter(Variable=="forc" & t<1000 & masking=="no" & ghg!="sai" & 
                                         tremoval %in% c(20,60,120,200)) %>%
                                group_by(t,gas,masking,tremoval) %>%
                                summarise(value=sum(value))) ) +
  geom_line(aes(x=t,
                y=value/pulse_size,
                color=as.factor(tremoval),
                linetype=masking),
            linewidth=1)  +
  geom_line(data=tsec %>%
              filter(Variable=="forc" & gas=="co2" & ghg!="sai" & masking=="no") %>%
              group_by(t,file,gas) %>%
              summarise(value=sum(value) ),
            aes(x=t,
                y=value/pulse_size ),
            linewidth=1,
            linetype=2)  +
  scale_color_viridis_d() +
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Radiative forcing variation [W*m^-2]')
fig3 <- ggpubr::ggarrange(temperature,forcing,common.legend=TRUE)
ggsave("figure_3.png",width=10,height=6,path=respath,plot=fig3)

###### totcost with removal 
totcostremtime <- tsecrem %>% 
  filter(Variable=="forc" & ghg=="sai" & masking=="yes") %>%
  mutate(value=-value) %>%
  inner_join(tsecrem %>% 
               filter(Variable=="T" & ghg=="co2" & masking=="no") %>%
               rename(deltatemp = value) %>% 
               select(t,gas,tremoval,deltatemp)) %>%
  inner_join(temp %>% 
               select(t,rcp,temp,initial_conditions)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=c(seq(0.01,0.05,by=0.01)))) %>%
  cross_join(data.frame(senscost= 1)) %>%
  cross_join(data.frame(eff= 0.85)) %>%
  cross_join(data.frame(sensdam= 1)) %>%
  cross_join(data.frame(removal=c(seq(500,0,by=-100)))) %>%
  rowwise() %>% mutate(gwpt = min(gwpmax,gwp* (1 + 0.022)^(t-1)) ) %>%
  mutate(impl = value*forctoUSD,
         dir = gwp*value*senscost*cost^powercost,
         dam = ( gwpt*(1- 1 / (1+(1-eff)*alpha*sensdam*(temp+deltatemp)^powerdam)) -  gwpt*(1- 1 / (1+(1-eff)*alpha*sensdam*(temp)^powerdam )) )  ) 

totcostrem <- totcostremtime %>%
  group_by(file,delta,tremoval,removal,senscost,sensdam,eff) %>%
  summarise(implnpv = sum( impl/(1+delta)^(t-1) ),
            dirnpv = sum( dir/(1+delta)^(t-1) ),
            damnpv = sum( dam / (1+delta)^(t-1) ) ) %>%
  inner_join(totcostremtime %>%
               group_by(file,delta,tremoval,removal,senscost,sensdam,eff) %>%
               filter(t==tremoval) %>%
               summarise(removalnpv = pulse_size*removal/(1+delta)^(tremoval-1) ) ) %>%
  mutate(cost = dirnpv + damnpv + implnpv + removalnpv) %>% 
  group_by(file,tremoval,removal,senscost,sensdam,eff) %>%
  mutate(costnorm = cost/cost[delta==0.01]) 
  
####### pareto fronts
fronts <- totcostrem %>%
  ungroup() %>% 
  inner_join(totcost %>% ungroup() %>%
              filter(gas=="co2") %>%
              dplyr::select(-file,-costnorm, -damnpv, -dirnpv, -implnpv) %>%
              rename(benchmark=cost) ) %>%
  filter(cost<benchmark*1.1) %>%
  group_by(file,delta,removal,senscost,eff) %>%
  filter(tremoval==min(tremoval)) %>%
  group_by(delta,tremoval,senscost,eff) %>%
  filter(removal==max(removal)) 
  

ggplot(fronts %>% filter(delta %in% c(0.01,0.03,0.05) & 
                           senscost %in% c(1) & 
                           sensdam==1 &
                           eff==0.85 ) ) +
  geom_line(aes(x=tremoval,
                  y=removal,
                  color=as.factor(paste0(delta*100," %")),
                  linetype=as.factor(senscost) ), 
            linewidth=1 ) +
  ylab('Removal cost [US$/tonCO2]') +
  xlab('Time of removal [years]') +
  scale_linetype_manual(values=c(2,1,3)) +
  guides(color=guide_legend(title = "Discount rate"),
         linetype=guide_legend(title = "Direct SAI costs [x central spec]"))
ggsave("figure_4.png",width=7,height=6,path=respath)


ggplot(totcostrem %>% filter(removal %in% c(100,200,500) & 
                               delta %in% c(0.01,0.03) & 
                               senscost==1 &
                               eff==0.85 & sensdam==1 )) +
  geom_line(aes(x=tremoval,
                y=cost/pulse_size,
                color=as.factor(removal),
                linetype=as.factor(delta))) +
  geom_hline(data=totcost %>%
               filter(gas=="co2" & 
                        delta %in% c(0.01,0.03) & 
                        senscost==1 &
                        eff==0.8 ),
             aes(yintercept=cost/pulse_size,
                 linetype=as.factor(delta))) +
  ylab('Net Present Cost [US$]') +
  xlab('Time of removal')  + 
  guides(color=guide_legend(title = "Cost of removal"),
         linetype=guide_legend(title = "Discount rate"))
ggsave("figure_5.png",width=7,height=6,path=respath)


ggplot(totcostrem %>% 
         filter( delta %in% c(0,0.01,0.04) & 
                   senscost==1 & 
                   eff == 0.8 & 
                   removal==100) %>%
         pivot_longer(c(implnpv,dirnpv,damnpv,removalnpv)) %>%
         bind_rows(totcost %>%   
                     filter( delta %in% c(0,0.01,0.04) & 
                                           senscost==1 & 
                                           eff == 0.8) %>%
                      pivot_longer(c(implnpv,dirnpv,damnpv)) %>%
                     mutate(tremoval=NA) ) %>%
         bind_rows(totcost %>%   
                     filter( delta %in% c(0,0.01,0.04) & 
                               senscost==1 & 
                               eff == 0) %>%
                     pivot_longer(c(damnpv)) %>%
                     mutate(tremoval=NaN) ) ) +
  geom_bar(aes(x=as.factor(tremoval),
               y=value/1000,
               fill=name),
           linewidth=1, 
           stat="identity",position="stack") + 
  ylab('Net Present Cost [US$/ton]') +
  xlab('Discount rate [%]') +
  facet_wrap(delta~.,scales="free")


ggplot(totcostrem %>% 
         filter( delta %in% c(0.01,0.04,0.05) & 
                   senscost==1 &
                   sensdam==1 &
                   eff %in% c(0.85) & 
                   removal %in% c(100,200,500) ) ) +
  geom_line(aes(x=tremoval,
               y=cost/pulse_size,
               color=as.factor(delta),
               linetype=as.factor(eff))) + 
  ylab('Net Present Cost [US$/ton]') +
  xlab('Time of removal [years]') +
  facet_wrap(removal~.,)


relative <- ggplot(totcostrem %>% 
                     filter(delta %in% c(0.01,0.02,0.03,0.04,0.05) & 
                              tremoval %in% c(20,60,120) )) +
  geom_line(data=.%>%filter(removal==500 & senscost == 1 & eff==0.8 & sensdam==1),
            aes(x=delta*100,
                y=costnorm,
                color=as.factor(tremoval)),
            linewidth=1) +
  geom_ribbon(data=.%>%
                group_by(delta,tremoval) %>%
                summarise(max=max(costnorm),min=min(costnorm)),
              aes(x=delta*100,
                  ymin=min,
                  ymax=max,
                  fill=as.factor(tremoval),
                  color=as.factor(tremoval)),
              linewidth=0.5,
              alpha=0.2) +
  geom_line(data=totcost %>% 
              filter(delta %in% c(0.01,0.02,0.03,0.04,0.05) & senscost == 1 & eff==0.8 & sensdam==1 & gas !="n2o"),
            aes(x=delta*100,
                y=costnorm,
                color=gas),
            linewidth=1) +
  geom_ribbon(data=totcost %>% 
                filter(delta %in% c(0.01,0.02,0.03,0.04,0.05) & gas !="n2o") %>%
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
  xlab('Discount rate [%]') 

ggplot(totcostremtime %>% 
         filter(tremoval %in% c(20,60,120,200) & 
                  delta==0.02 & 
                  removal==100 & 
                  eff==0.85 & 
                  senscost==1) %>%
         mutate(rem = ifelse(t==tremoval,removal*pulse_size,0)) ) +
  geom_line(aes(x=t,
                y=(dir+dam+rem+impl)/pulse_size,
                color=as.factor(tremoval) ), 
            linewidth=1)  + 
  geom_line(data=totcosttime %>% 
              filter(senscost==1 & 
                       sensdam==1 &
                       eff==0.85 & 
                       delta==0.02 & 
                       gas=="co2"),
            aes(x=t,
                y=(dir+dam+impl)/pulse_size ), 
            linewidth=1,linetype=2)  + 
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Cost of strategy [US$/yr]')

