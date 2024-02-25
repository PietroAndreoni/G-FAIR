alpha <- 0.23 #[%GDP/K^2]
sens <- 0.84 #[0.84 K/(W*m^2)]
direct_costs <- 0.1 #[%GDP/(W*m^2)]
power <- 2 #exponent of direct costs
Cmax <- 4#
cost <- 0.2

study <- cross_join(data.frame(C=seq(0,Cmax,0.5)),
                    data.frame(g=seq(0,1,0.1))) %>%
  cross_join(data.frame(theta=seq(0,60,by=10)*pi/180)) %>%
  mutate(direct = direct_costs*(g*C)^power,
         perfect_masking = alpha * (sens*C - sens*g*C)^2,
         indirect = alpha * (sens*C)^2 * (1 + g^2 - 2*g*cos(theta)),
         abatement = cost * (Cmax - C)^2,
         total = direct + indirect,
         efficacy = (1 + g^2 - 2*g*cos(theta)) ) %>%
  pivot_longer(c(total,perfect_masking,direct,indirect,efficacy)) %>%
  group_by(C,theta) %>%
  mutate(valuenorm=ifelse(name!="efficacy",value/value[name=="total" & g==0],value))

ggplot(study %>% 
         filter(theta==30*pi/180 & C %in% c(1,4) & name!="efficacy")) +
  geom_line(aes(x=g,
                y=valuenorm,
                color=name,
                linetype=as.factor(C) ),
            linewidth=1)

ggplot(study %>% 
         filter(C==1 & name=="efficacy")) +
  geom_line(aes(x=g,
                y=value,
                color=as.factor(theta*180/pi) ),
            linewidth=1)
