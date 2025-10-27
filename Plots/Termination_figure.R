damnpv_term <- totcosttime_term %>%
  cross_join(data.frame(xalpha=c(0.2,0.5,1,2,5)) ) %>% 
  cross_join(data.frame(xtheta=c(10,30)) ) %>% 
  mutate(impl = srm_masking*forctoUSD,
         dir = gwp*srm_masking*cost,
         dam = gwp*(alpha* xalpha * (temp**2-temp_srm**2)),
         mask = 2 * gwp * alpha * xalpha * temp_ghg *  (1-cos(xtheta * pi/180)) * srm_masking * ecs * (1+srm/forc),
         mask_belaiia = gwp * alpha * xalpha * ( temp_ghgpulse ** 2 * eff(srm_masking+srm,srm_masking+forc,xtheta * pi/180) - temp_ghg ** 2 * eff(srm,forc,xtheta * pi/180) ) ) %>% 
  cross_join(data.frame(delta=seq(0.01,0.06,by=0.01))) %>%
  filter(t>=as.numeric(pulse_time)) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,term,delta,xalpha,xtheta) %>%
  summarise(masknpv = sum( mask / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            implnpv = sum( impl / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            dirnpv = sum( dir / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE),
            damnpv = sum( dam / (1+delta)^(t-as.numeric(pulse_time)), na.rm = TRUE)) %>% 
  ungroup() %>%  mutate(costnpv = masknpv + implnpv + dirnpv+damnpv,term=as.numeric(term))

costnpv_term <- damnpv_term  %>% filter(term<=100) %>% 
  cross_join(data.frame(prob=seq(0,0.05,by=0.005))) %>%
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,prob,xalpha,xtheta) %>%
  arrange(term) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,prob,xalpha,xtheta) %>%
  mutate(dt = term - lag(term), 
         p = prob/10 * (1-prob/10) ^ (term-1)) %>% 
  mutate(dt = ifelse(is.na(dt),10,dt)) %>% 
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,prob,xalpha,xtheta) %>%
  summarise(termnpv =  sum( dt * p * costnpv )) 
  
costnpv_noterm <- damnpv_term  %>% filter(as.numeric(term)==1000) %>% 
  cross_join(data.frame(prob=seq(0,0.05,by=0.005))) %>%
  group_by(rcp,cool_rate,pulse_time,geo_end,gas,delta,prob) %>%
  mutate(termnpv_base = (1-prob/10) ^ (100) * costnpv) %>% 
  select(rcp,gas,cool_rate,pulse_time,geo_end,prob,,xalpha,xtheta,termnpv_base)

costnpv <- full_join(costnpv_term,costnpv_noterm) %>% 
  mutate(termnpv = termnpv + termnpv_base)

ggplot(costnpv %>% filter(xtheta==10 & xalpha==1 & delta==0.02 & xtheta==10 & xalpha==1) %>% 
         group_by(gas) %>% 
         mutate(costnpvnorm=termnpv/termnpv[prob==0 & delta==0.02 & xtheta==10 & xalpha==1]) )  +
  geom_ribbon(data= . %>% 
                group_by(gas,prob) %>%  
                summarise(min=min(costnpvnorm), max=max(costnpvnorm)),
              aes(x=prob,
                  ymin=min,
                  ymax=max,
                  fill=gas,
                  color=gas),
              linewidth=0.5,
              alpha=0.4 ) +
  geom_line(data= . %>% filter(delta==0.02 & cool_rate=="10" & geo_end=="2300" & rcp=="45" & pulse_time=="2" & xalpha==1 & xtheta==10), 
            aes(x=prob,
                y=costnpvnorm, 
                color=gas),
            linewidth=1) + 
  ylim(c(0,5))
