
damnpv <- bind_rows(lapply(list.files(output_folder, pattern = "\\.csv$", full.names = TRUE), read.csv)) %>% select(-X)

damnpv <- read.csv("Results_output/output_analysis_dt.csv")

damnorm <- damnpv %>% select(-t) %>% 
  mutate(geo_start=ifelse(cool_rate==0,2500,geo_start),
         geo_end =ifelse(cool_rate==0,max(geo_end)+100,geo_end )) %>% 
  mutate_if(is.integer, as.numeric) %>%
  group_by(gas) %>% 
  mutate(npc_std=(npc_srm-median(npc_srm,na.rm=TRUE))/mad(npc_srm,na.rm=TRUE),
         npc_norm=npc_srm/median(npc_srm,na.rm=TRUE)) %>% 
  group_by(gas) %>% 
  filter(npc_srm>quantile(npc_srm,0.001) & npc_srm<quantile(npc_srm,0.999)) 

damnpv %>% 
  group_by(gas) %>% 
  summarise(npc_mad = mad(npc_srm,na.rm=TRUE), 
            npc_med = median(npc_srm,na.rm=TRUE),
            p99 = median(npc_srm,na.rm=TRUE))
         
damnorm %>% 
  ggplot() + geom_qq(aes(sample=npc_std,color=gas)) +
  geom_abline(slope=1,intercept=0) 

stats::ks.test(damnorm %>% filter(gas=="co2") %>% pull(npc_std),
               damnorm %>% filter(gas=="ch4") %>% pull(npc_std))


th_co2 <- quantile(damnorm %>% filter(gas=="co2") %>% pull(npc_srm), 0.95)
th_ch4 <- quantile(damnorm %>% filter(gas=="ch4") %>% pull(npc_srm), 0.95)

quantile(damnorm %>% filter(gas=="co2") %>% pull(npc_srm), 0.99)/quantile(damnorm %>% filter(gas=="co2") %>% pull(npc_srm), 0.95)
quantile(damnorm %>% filter(gas=="ch4") %>% pull(npc_srm), 0.99)/quantile(damnorm %>% filter(gas=="ch4") %>% pull(npc_srm), 0.95)

fit_co2 <- POT::fitgpd(damnorm %>% filter(gas=="co2") %>% pull(npc_srm), threshold = th_co2)
fit_ch4 <- POT::fitgpd(damnorm %>% filter(gas=="ch4") %>% pull(npc_srm), threshold = th_ch4)

fit_co2$fitted.values["shape"] + fit_co2$std.err["shape"]   # ξ for CO2
fit_ch4$fitted.values["shape"] + fit_ch4$std.err["shape"]   # ξ for CH4


adjust_opt <- 1.06 * sd(damnorm$npc_norm)*(length(damnorm$npc_norm))^(-1/5)

density_plot <- ggplot(damnorm) +
  geom_density(aes(x=npc_norm,
                   color=gas),
               adjust=adjust_opt,
               linewidth=1) +
geom_point(aes(x=npc_norm,
               color=gas,
               y=-(as.numeric(as.factor(gas))-1)/100),
           shape=108) +
  coord_cartesian(xlim = c(quantile(damnorm$npc_norm, 0.01, na.rm=TRUE), 
                           quantile(damnorm$npc_norm, 0.99, na.rm=TRUE)) ) +
  xlab("Normalized present cost") + ylab("density") 

ggplot(damnorm %>% filter(gas=="ch4")) +
  geom_density(aes(x=npc_srm),adjust=5,color="red") +
  coord_cartesian(xlim = c(0, quantile( (damnorm %>% filter(gas=="ch4"))$npc_srm, 0.95, na.rm=TRUE)) ) +
  xlab("Net present cost [$/tonCH4]") + ylab("")

ggplot(damnorm %>% filter(gas=="co2")) +
  geom_density(aes(x=scc,color=gas),adjust=5) +
  coord_cartesian(xlim = c(0, quantile((damnorm %>% filter(gas=="co2"))$scc, 0.95, na.rm=TRUE)) ) +
  xlab("Net present cost [$/tonCH4]") + ylab("") 

ggplot(damnorm) +
  geom_boxplot(aes(x=gas,y=scc,color=gas)) +
  xlab("Net present cost [$/tonCH4]") + ylab("") +
  facet_wrap(gas~.,scales="free")

library(gsaot)
gsoat_data <- damnorm %>% filter(gas=="ch4") %>% ungroup()

stat_analysis_ch4 <- ot_indices_1d(gsoat_data %>% 
                                     select(-gas,-pulse_size,-costnpv,-damnpv,-scc,-npc_srm, -npc_norm),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_ch4 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

gsoat_data <- damnorm %>% filter(gas=="co2") %>% ungroup()

stat_analysis_co2 <- ot_indices_1d(gsoat_data %>% 
                                     select(-gas,-pulse_size,-costnpv,-damnpv,-scc,-npc_srm, -npc_norm),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_co2 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

input_categories <- c("Climate","Climate","Socio-economic","Socio-economic","background SAI",
                      "background SAI","background SAI","Impacts","Normative","Impacts","Normative",
                      "Socio-economic","Normative","Impacts","Normative",
                      "Impacts")
names(input_categories) <- stat_analysis_ch4$indices_ci$input

importance_ch4 <- ggplot(as_tibble(stat_analysis_ch4$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input_categories[input]),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_ch4$indices) +
  geom_errorbar(aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [ch4]") + xlab("") + 
  theme(legend.position = "top") + scale_fill_viridis_d()

importance_co2 <- ggplot(as_tibble(stat_analysis_co2$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input_categories[input]),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_co2$indices) +
  geom_errorbar(aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),
                stat="identity",position="dodge",color="black") +
  ylab("importance [co2]") + xlab("")  + scale_fill_viridis_d()

importances <- ggpubr::ggarrange(importance_co2,importance_ch4,common.legend=TRUE,nrow=1)


dr <- ggplot(damnorm) +
  geom_smooth(aes(x=delta*100,y=npc_std,color=gas),
              method="loess") +
  xlab("Discount rate [%]") + ylab("")+
  coord_cartesian(ylim = c(0,5) ) 

nm <- ggplot(damnorm) +
  geom_smooth(aes(x=term,y=npc_std,color=gas),
              method="loess") +
  xlab("Termination year") + ylab("") +
  coord_cartesian(ylim = c(0,5) ) 

ecs <- ggplot(damnorm) +
  geom_smooth(aes(x=ecs/10,y=npc_std,color=gas),
              method="loess") +
  xlab("Equilibrium climate sensitvity [K]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$ecs/10, 0.01, na.rm=TRUE), 
                           quantile(damnorm$ecs/10, 0.99, na.rm=TRUE)),
                  ylim = c(0,5) ) 


alpha <- ggplot(damnorm) +
  geom_smooth(aes(x=alpha*100,y=npc_std,color=gas),
              method="loess") +
  geom_density(aes(x=alpha*100)) +
  xlab("Climate damages [%gdp/K^2]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$alpha*100, 0.01, na.rm=TRUE), 
                           quantile(damnorm$alpha*100, 0.99, na.rm=TRUE)),
                  ylim = c(0,5) ) 

theta <- ggplot(damnorm) +
  geom_smooth(aes(x=theta,y=npc_std,color=gas),
              method="loess") +
  geom_density(aes(x=theta)) +
  xlab("SRM efficacy angle [°]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$theta, 0.01, na.rm=TRUE), 
                           quantile(damnorm$theta, 0.99, na.rm=TRUE)),
                  ylim = c(0,5) ) 

tcr <- ggplot(damnorm) +
  geom_smooth(aes(x=tcr/10,y=npc_std,color=gas),
              method="loess") +
  xlab("Transient climate response [K]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$tcr/10, 0.01, na.rm=TRUE), 
                           quantile(damnorm$tcr/10, 0.99, na.rm=TRUE)),
                  ylim = c(0,5)) 

pulse_time <- ggplot(damnorm %>% filter(pulse_time < term)) +
  geom_smooth(aes(x=pulse_time+2020,y=npc_std,color=gas),
              method="loess") +
  xlab("Transient climate response [K]") + ylab("")


stat_analysis_diff <- ot_indices_1d(damnorm_diff %>% 
                                      select(-damnpv,-norm_diff,-co2,-ch4),
                                    damnorm_diff %>% pull(norm_diff), 
                                    M= 15,
                                    boot = T,
                                    R = 100)
lowerbound_diff <- irrelevance_threshold(damnorm_diff %>% pull(norm_diff), M= 15, solver="1d")

ggplot(as_tibble(stat_analysis_diff$indices_ci)  ) + 
  geom_bar(aes(x=reorder(input, -original),
               y=original,
               fill=input),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_diff$indices) +
  geom_errorbar(aes(x=input,
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [difference]") + xlab("") + 
  theme(legend.position = "none")
