experiment <- "pulse"
ghg <- c("co2")
tremoval <- seq(20,300,by=20)
tstart <- "historical_run"


tsecrem <- gdxtools::batch_extract("save_delta",
                                files=unique(paste0("../Results/",rcp,"_EXP",experiment,"_REM",tremoval,"_GAS",ghg,"_IC",tstart,".gdx")))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
    file=str_remove_all(file,"../Results/|.gdx"),
    gas=str_extract(file,"(?<=GAS).+?(?=_)"),
    rcp=str_extract(file,"(?<=RCP).+?(?=_)"),
    experiment=str_extract(file,"(?<=EXP).+?(?=_)"),
    tremoval=as.numeric(str_extract(file,"(?<=REM).+?(?=_)")),
    initial_conditions=str_extract(file,"(?<=IC).*")) 


###### totcost with removal 
totcostrem <- tsecrem %>% 
  filter(Variable=="forc" & gas=="co2") %>%
  group_by(t,file,tremoval) %>%
  summarise(value=sum(value)) %>%
  full_join(tsecrem %>% 
              filter(Variable=="T" & ghg=="co2") %>%
              rename(deltatemp = value) %>% 
              select(t,gas,tremoval,deltatemp)) %>%
  inner_join(temp %>% 
               select(t,rcp,temp,initial_conditions)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=seq(0.01,0.1,by=0.02))) %>%
  cross_join(data.frame(removal=c(seq(1000,0,by=-20)))) %>%
  cross_join(data.frame(senscost= c(0.1,0.5,1,2,10,100))) %>%
#  rowwise() %>% mutate(value=max(value,0)) %>%
  group_by(file,delta,tremoval,removal,senscost) %>%
#  filter(t<tremoval) %>%
  summarise(dir = sum( value*forctoUSD*senscost /(1+delta)^(t-1) ),
            dam = sum( ( (1-eff)*alpha*(temp+deltatemp-temp0)^power -  (1-eff)*alpha*(temp-temp0)^power ) / (1+delta)^(t-1) ) ) %>% 
  ungroup() %>% 
  mutate(cost = dir + dam + 1000*removal/(1+delta)^(tremoval-1)) %>%
  unique()
  

temperature <- ggplot(tsecrem %>% 
                        filter(Variable=="T" & ghg=="co2" & tremoval %in% c(20,60,120,200) )) +
  geom_line(aes(x=t,y=value/1000,color=as.factor(tremoval)), linewidth=1) +
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Temperature variation [K]')

forcing <- ggplot(tsecrem %>% 
                    filter(Variable=="forc" & tremoval %in% c(20,60,120,200)) %>%
                    group_by(t,file,tremoval) %>%
                    summarise(value=sum(value) ) ) +
  geom_line(aes(x=t,y=value/1000,color=as.factor(tremoval) ), linewidth=1)  + 
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Radiative forcing variation [W*m^-2]')
fig3 <- ggpubr::ggarrange(temperature,forcing,common.legend=TRUE)
ggsave("figure_3.png",width=10,height=6)

deltaplot <- c(0.03,0.06,0.09)

####### pareto fronts
fronts <- totcostrem %>%
  ungroup() %>% 
  full_join(totcost %>% ungroup() %>%
              filter(str_detect(file,"co2")) %>%
              dplyr::select(-file,-costnorm, -dam, -dir) %>%
              rename(benchmark=cost) ) %>%
  filter(cost<benchmark*1.1) %>%
  group_by(file,delta,removal,senscost) %>%
  filter(tremoval==min(tremoval)) %>%
  group_by(delta,tremoval,senscost) %>%
  filter(removal==max(removal) & removal!=1000) 
  

ggplot(fronts %>% filter(delta %in% c(0.01,0.03,0.05) & senscost %in% c(0.1,1,10)   ) ) +
  # stat_smooth(aes(x=tremoval,
  #                 y=removal,
  #                 color=as.factor(paste0(delta*100," %")),
  #                 linetype=as.factor(senscost) ),
  #             method="lm",
  #             formula=y ~ poly(x, 4, raw=TRUE),
  #             se=FALSE,
  #             fullrange=TRUE ) +
  geom_line(aes(x=tremoval,
                  y=removal,
                  color=as.factor(paste0(delta*100," %")),
                  linetype=as.factor(senscost) ), 
            linewidth=1 ) +
  xlim(c(5,300)) + ylim(c(0,1000))+
  ylab('Removal cost [US$/tonCO2]') +
  xlab('Time of removal [years]') + 
  guides(color=guide_legend(title = "Discount rate"))
ggsave("figure_4.png",width=7,height=6)


ggplot(totcostrem %>% filter(tremoval & removal %in% c(500,200,100,0) & delta==0.01 & senscost==10 )) +
  geom_line(aes(x=tremoval,
                y=cost/1000,
                color=as.factor(removal),
                linetype=as.factor(delta))) +
  geom_hline(data=totcost %>%
               filter(str_detect(file,"co2") & delta==0.01 & senscost==10),
             aes(yintercept=cost/1000,
                 linetype=as.factor(delta))) +
  ylab('Net Present Cost [US$]') +
  xlab('Time of removal')  + 
  guides(color=guide_legend(title = "Cost of removal"),
         linetype=guide_legend(title = "Discount rate"))
ggsave("figure_5.png",width=7,height=6)

