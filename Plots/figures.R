#### load data for figure 1
require(tidyverse)
require(gdxtools)

experiment <- "pulse"
ghg <- c("n2o","ch4","co2")
tstart <- "2020"

tsec <- gdxtools::batch_extract("save_delta",
                        files=paste0("../",experiment,"_",ghg,"_",tstart,".gdx"))$save_delta %>%
  as.data.frame() %>% 
  rename(Variable=V3,file=gdx) %>%
  mutate(t=as.numeric(t),
         file=str_remove(file,".gdx"),
         file=str_remove(file,"../"),
         gas=str_sub(file,7,9)) 

# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
forctoTg <- 1/0.2

# from https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7/pdf  
TgtoUSD <- 2250*10^6

totcost <- tsec %>% 
  filter(Variable=="forc") %>%
  group_by(t,file,gas) %>%
  summarise(value=sum(value)) %>%
  ungroup() %>% cross_join(data.frame(delta=seq(0.01,0.2,by=0.01))) %>%
  group_by(file,delta,gas) %>%
  summarise(cost = sum(value*forctoTg*TgtoUSD/(1+delta)^(t-1))) %>%
  group_by(file,gas) %>%
  mutate(costnorm = cost/cost[delta==0.01]) 


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
ggsave("figure_1.png",width=10,height=6)


#### figure 2
absolute <- ggplot(totcost %>% 
                     filter(delta<=0.1  ) %>%
                     mutate(cost=ifelse(gas=="n2o",cost/100,cost))) +
  geom_line(aes(x=delta*100,
                y=cost/1000,
                color=gas),linewidth=1) +
  ylab('Net Present Cost, CO2 and CH4 [US$]') +
  xlab('Discount rate [%]') + 
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 100, name = "Net Present Cost, N20 [US$]") )

relative <- ggplot(totcost %>% 
                     filter(delta<=0.1)) +
  geom_line(aes(x=delta*100,
                y=costnorm,
                color=gas),linewidth=1) +
  ylab('Relative cost [frac]') +
  xlab('Discount rate [%]') +
  ylim(c(0,1))

fig2 <- ggpubr::ggarrange(absolute,relative,common.legend=TRUE)
ggsave("figure_2.png",width=10,height=6)
