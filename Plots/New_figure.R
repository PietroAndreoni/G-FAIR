require(tidyverse)

## climate damage function parameters
alpha <- 0.595*1.25/100 # from Howard and Stern
gwp <- 105*1e12 # initial world gdp
cost <- 0.001 # 0.1% GDP per W/m2 as Belaiia et al
ecs <- 3.24/3.71 # consistent with FAIR model
forctoTg <- 1/0.2 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897825/ 
TgtoUSD <- 2250*10^6 # from https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7/pdf  
forctoUSD <- forctoTg * TgtoUSD # US$/(W/m^2)
gwpmax <- gwp* (1 + 0.022)^(180-1)

damnpv_norm <- totcosttime %>% 
  mutate(xalpha=1,
         xtheta=10,
         impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha* xalpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * xalpha * temp_ghg *  (1-cos(xtheta * pi/180)) * srm_masking * ecs * (1+srm/forc),
         mask_belaiia = gwp * alpha * xalpha * ( temp_ghgpulse ** 2 * eff(srm_masking+srm,srm_masking+forc,xtheta * pi/180) - temp_ghg ** 2 * eff(srm,forc,xtheta * pi/180) ) ) %>%
  cross_join(data.frame(delta=seq(0.01,0.06,by=0.01))) %>%
  filter(t>=as.numeric(pulse_time)) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha,xtheta) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

pulse_size <- W_EMI %>% 
  filter(ghg==gas) %>% 
  filter(experiment %in% c("srm","srmpulse") ) %>% 
  group_by(ghg,t,rcp,gas,cool_rate,pulse_time,geo_end) %>% 
  mutate(pulse_size=value-value[experiment=="srm"]) %>% 
  filter(experiment=="srmpulse" & pulse_size!=0) %>% 
  ungroup() %>% 
  select(-experiment,-ghg,-t,-value,-gdx,-file) %>% 
  mutate(pulse_size=ifelse(gas=="co2",pulse_size*1e9,pulse_size*1e6))

f1a <- ggplot(damnpv_norm %>%
         full_join(pulse_size) %>% 
         filter(gas %in% c("co2","ch4") & xalpha %in% c(1) & xtheta %in% c(10) ) %>% 
         group_by(gas) %>% 
         mutate(costnpvnorm=costnpv/costnpv[delta==0.02 & cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10]) ) +
  geom_ribbon(data=. %>% 
                group_by(gas,delta) %>%  
                summarise(min=min(costnpvnorm), max=max(costnpvnorm)),
              aes(x=delta,
                   ymin=min,
                   ymax=max,
                   fill=gas,color=gas),linewidth=0.5,alpha=0.4 ) +
  geom_line(data=. %>% filter(cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" &  xalpha==1 & xtheta==10), 
            aes(x=delta,
                y=costnpvnorm, 
                color=gas),linewidth=1 ) +
  ylim(c(0,5))


damnpv_alpha<- totcosttime %>% 
  cross_join(data.frame(xalpha=c(0.2,0.5,1,2,5)) ) %>% 
  mutate(xtheta=10, 
         delta=0.02,
         impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha* xalpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * xalpha * temp_ghg *  (1-cos(xtheta * pi/180)) * srm_masking * ecs * (1+srm/forc),
         mask_belaiia = gwp * alpha * xalpha * ( temp_ghgpulse ** 2 * eff(srm_masking+srm,srm_masking+forc,xtheta * pi/180) - temp_ghg ** 2 * eff(srm,forc,xtheta * pi/180) ) ) %>%
  filter(t>=as.numeric(pulse_time)) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha,xtheta) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

f1b <- ggplot(damnpv_alpha %>%
         filter(gas %in% c("co2","ch4")) %>% 
         group_by(gas) %>% 
         mutate(costnpvnorm=costnpv/costnpv[cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10] ) )  +
  geom_ribbon(data=. %>% 
                group_by(gas,xalpha,xtheta) %>%  
                summarise(min=min(costnpvnorm), max=max(costnpvnorm)),
              aes(x=xalpha,
                  ymin=min,
                  ymax=max,
                  fill=gas,
                  color=gas),
              linewidth=0.5,
              alpha=0.4 ) +
  geom_line(data=. %>% filter(cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10), 
            aes(x=xalpha,
                y=costnpvnorm, 
                color=gas),
            linewidth=1) +
  ylim(c(0,5))
  


damnpv_theta <- totcosttime %>% 
  cross_join(data.frame(xtheta=c(0,5,10,15,20,25,30)) ) %>% 
  mutate(xalpha=1,
         delta=0.02,
         impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha* xalpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * xalpha * temp_ghg *  (1-cos(xtheta * pi/180)) * srm_masking * ecs * (1+srm/forc),
         mask_belaiia = gwp * alpha * xalpha * ( temp_ghgpulse ** 2 * eff(srm_masking+srm,srm_masking+forc,xtheta * pi/180) - temp_ghg ** 2 * eff(srm,forc,xtheta * pi/180) ) ) %>%
  filter(t>=as.numeric(pulse_time)) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha,xtheta) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)

f1c <- ggplot(damnpv_theta %>%
         filter(gas %in% c("co2","ch4") & delta==0.02) %>% 
         group_by(gas) %>% 
         mutate(costnpvnorm=costnpv/costnpv[delta==0.02 & cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10] ) )  +
  geom_ribbon(data= . %>% 
                group_by(gas,xalpha,xtheta) %>%  
                summarise(min=min(costnpvnorm), max=max(costnpvnorm)),
              aes(x=xtheta,
                  ymin=min,
                  ymax=max,
                  fill=gas,
                  color=gas),
              linewidth=0.5,
              alpha=0.4 ) +
  geom_line(data= . %>% filter(delta==0.02 & cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10), 
            aes(x=xtheta,
                y=costnpvnorm, 
                color=gas),
            linewidth=1)  +
  ylim(c(0,5))

ggpubr::ggarrange(f1a,f1b,f1c,nrow=1,common.legend=TRUE)

scc <- totcosttime %>% 
  cross_join(data.frame(xalpha=c(0.2,0.5,1,2,5)) ) %>% 
  mutate(dam = gwp*(alpha* xalpha * ((temp+deltatemp)**2-temp**2)) ) %>% 
  cross_join(data.frame(delta=seq(0.01,0.06,by=0.01)) ) %>%
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha) %>%
  summarise(damnpv = sum( dam / (1+delta)^t, na.rm = TRUE) ) %>%
  full_join(pulse_size) %>% 
  mutate(scc=damnpv/pulse_size ) %>% 
  select(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha,scc)

ggplot(scc %>%
         full_join(pulse_size) %>% 
         filter(delta==0.02  & gas %in% c("co2","ch4"))) +
  geom_boxplot(aes(x=as.factor(xalpha),
                y=scc, 
                color=rcp) ) +
  facet_wrap(gas~.,scales="free") 

ggplot(damnpv_theta %>% full_join(pulse_size) %>% 
         filter(gas %in% c("co2","ch4")) %>% 
         mutate(costnpvnorm=(costnpv/pulse_size)) )  +
  geom_ribbon(data= . %>% 
                group_by(gas,xalpha,xtheta) %>%  
                summarise(min=min(costnpvnorm), max=max(costnpvnorm)),
              aes(x=xtheta,
                  ymin=min,
                  ymax=max,
                  fill=gas,
                  color=gas),
              linewidth=0.5,
              alpha=0.4 ) +
  geom_line(data= . %>% filter(delta==0.02 & cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10), 
            aes(x=xtheta,
                y=costnpvnorm, 
                color=gas),
            linewidth=1) + facet_wrap(gas~.,scales="free")

damnpv_all <- totcosttime %>% 
  cross_join(data.frame(xalpha=c(0.2,0.5,1,2,5)) ) %>% 
  cross_join(data.frame(xtheta=c(0,5,10,15,20,25,30)) ) %>% 
  mutate(impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha* xalpha * ((temp+deltatemp)**2-temp**2)),
         mask = 2 * gwp * alpha * xalpha * temp_ghg *  (1-cos(xtheta * pi/180)) * srm_masking * ecs * (1+srm/forc),
         mask_belaiia = gwp * alpha * xalpha * ( temp_ghgpulse ** 2 * eff(srm_masking+srm,srm_masking+forc,xtheta * pi/180) - temp_ghg ** 2 * eff(srm,forc,xtheta * pi/180) ) ) %>%
  filter(t>=as.numeric(pulse_time)) %>% 
  cross_join(data.frame(delta=seq(0.01,0.06,by=0.01)) ) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,xalpha,xtheta) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv)
