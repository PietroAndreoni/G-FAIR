experiment <- "pulse_removal"
ghg <- c("co2")
tremoval <- c(seq(10,100,by=10),seq(100,200,by=25),250,300)
tstart <- "2020"

tsecrem <- gdxtools::batch_extract("save_delta",
                                files=unique(paste0("../",experiment,tremoval,"_",ghg,"_",tstart,".gdx")))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
         file=str_remove(file,".gdx"),
         tremoval=str_sub(file,17,19)) %>%
  mutate(tremoval=as.numeric(str_remove(tremoval,"_")))

###### totcost with removal 
totcostrem <- tsecrem %>% 
  filter(Variable=="forc" & str_detect(file,"co2")) %>%
  group_by(t,file,tremoval) %>%
  summarise(value=sum(value)) %>%
  ungroup() %>% 
  cross_join(data.frame(delta=seq(0.01,0.1,by=0.01))) %>%
  cross_join(data.frame(removal=c(seq(500,0,by=-10)))) %>%
#  rowwise() %>% mutate(value=max(value,0)) %>%
  group_by(file,delta,tremoval,removal) %>%
  filter(t<tremoval) %>%
  reframe(cost = sum(value*forctoTg*TgtoUSD/(1+delta)^(t-1)) + 1000*removal/(1+delta)^(tremoval-1)) %>% 
  unique()

temperature <- ggplot(tsecrem %>% 
                        filter(Variable=="T" & ghg=="co2" & tremoval %in% c(10,50,100,200) )) +
  geom_line(aes(x=t,y=value/1000,color=as.factor(tremoval)), linewidth=1) +
  guides(color = guide_legend(title = "Time of removal [years]"))  +
  xlab('time [years]') +
  ylab('Temperature variation [K]')

forcing <- ggplot(tsecrem %>% 
                    filter(Variable=="forc" & tremoval %in% c(10,50,100,200)) %>%
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
              dplyr::select(-file,-costnorm) %>%
              rename(benchmark=cost) ) %>%
  filter(cost<benchmark*1.05) %>%
  group_by(file,delta,removal) %>%
  filter(tremoval==min(tremoval)) %>%
  group_by(delta,tremoval) %>%
  filter(removal==max(removal) & removal!=500) 
  

ggplot(fronts %>% filter(delta %in% c(0.01,0.03,0.05,0.09) ) ) +
  stat_smooth(aes(x=tremoval,y=removal,color=as.factor(paste0(delta*100," %") )),
              method="lm",
              formula=y ~ poly(x, 4, raw=TRUE),
              se=FALSE,
              fullrange=TRUE ) +
  xlim(c(5,300)) + ylim(c(0,500))+
  ylab('Removal cost [US$/tonCO2]') +
  xlab('Time of removal [years]') + 
  guides(color=guide_legend(title = "Discount rate"))
ggsave("figure_4.png",width=7,height=6)


ggplot(totcostrem %>% filter(tremoval > 20 & removal %in% c(500,200,100) & delta %in% deltaplot )) +
  geom_line(aes(x=tremoval,
                y=cost/1000,
                color=as.factor(removal),
                linetype=as.factor(delta))) +
  geom_hline(data=totcost %>%
               filter(str_detect(file,"co2") & delta%in% deltaplot),
             aes(yintercept=cost/1000,
                 linetype=as.factor(delta))) +
  ylab('Net Present Cost [US$]') +
  xlab('Time of removal')  + 
  guides(color=guide_legend(title = "Cost of removal"),
         linetype=guide_legend(title = "Discount rate"))
ggsave("figure_5.png",width=7,height=6)

