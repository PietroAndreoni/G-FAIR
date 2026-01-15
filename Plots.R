library(tidyverse)
output_folder <- "C:/Users/Andreoni/OneDrive - The University of Chicago/G-FAIR/Results_2112"
damnpv <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "npc_output")), read.csv)) 
scc <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "sccnosrm_output")), read.csv)) %>% rename(scc=scc_nosrm)
scc_srm <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "scc_output")), read.csv)) %>% rename(scc_srm=scc)
all_cols <- names(damnpv)[1:19]

# filter scc and npc>0 and scenarios with both CH4 and CO2
damnpv <- damnpv %>% 
  inner_join(scc %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  inner_join(scc_srm %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  filter(npc_srm>0 & scc>0 & scc_srm>0) %>% 
  group_by_at(setdiff(all_cols,c("gas")) ) %>% 
  filter(n()==2 & !any(duplicated(gas))) %>% 
  ungroup() %>% unique()

scc %>% 
  select(gas,ecs,tcr,rcp,alpha,delta,vsl,vsl_eta,dg,mortality_ozone,damnpv,pulse_size) %>% 
  mutate(scc= damnpv/pulse_size) %>% select(-damnpv,-pulse_size) %>% 
  unique() %>% pivot_wider(names_from="gas",values_from="scc") %>% 
  mutate(ratio=ch4/co2) %>% 
  ggplot() + geom_density(aes(x=ratio))
  
frac_damages <- damnpv %>% 
  pivot_longer(c(dirnpv,srmpnpv,ozpnpv,masknpv,damnpv),
               names_to="source",values_to="npc_partial") %>% 
  group_by_at(all_cols) %>% 
  mutate(npc_partial=npc_partial/pulse_size ) %>% 
  group_by_at(all_cols) %>%
  mutate(npc_frac=npc_partial/sum(npc_partial,na.rm=TRUE))

fraction <- frac_damages %>% filter(npc_frac!=0  & npc_frac>0) %>% 
  ggplot()+
  geom_density(aes(x=npc_frac,y=after_stat(scaled),color=source),adjust=1,linewidth=1.5) +
  facet_wrap(gas~.,) + ggpubr::theme_pubr() +
  xlab("Fraction of total cost") + ylab("Density")

damnorm <- damnpv %>% unique() %>% 
  mutate(geo_start=ifelse(cool_rate==0,2500,geo_start),
         geo_end =ifelse(cool_rate==0,max(geo_end)+100,geo_end )) %>% 
  mutate_if(is.integer, as.numeric) %>%
  group_by(gas) %>% 
  mutate(npc_std=(npc_srm-median(npc_srm,na.rm=TRUE))/mad(npc_srm,na.rm=TRUE),
         npc_norm=npc_srm/median(npc_srm,na.rm=TRUE))

damnpv %>% 
  group_by(gas) %>% 
  summarise(npc_mad = mad(npc_srm,na.rm=TRUE),
            p25 = quantile(npc_srm,0.25,na.rm=TRUE), 
            p75 = quantile(npc_srm,0.75,na.rm=TRUE),  
            p95 = quantile(npc_srm,0.95,na.rm=TRUE),
            npc_med = median(npc_srm,na.rm=TRUE))
    
scc %>% 
  group_by(gas) %>% 
  summarise(npc_mad = mad(scc,na.rm=TRUE), 
            p25 = quantile(scc,0.25,na.rm=TRUE), 
            p75 = quantile(scc,0.75,na.rm=TRUE),
            npc_med = median(scc,na.rm=TRUE))

stats::ks.test(damnorm %>% filter(gas=="co2") %>% pull(npc_std),
               damnorm %>% filter(gas=="ch4") %>% pull(npc_std))

quantile(damnpv %>% filter(gas=="co2") %>% pull(npc_srm), 0.99)/quantile(damnpv %>% filter(gas=="co2") %>% pull(npc_srm), 0.95)
quantile(damnpv %>% filter(gas=="ch4") %>% pull(npc_srm), 0.99)/quantile(damnpv %>% filter(gas=="ch4") %>% pull(npc_srm), 0.95)

adjust_opt <- log(1.06 * sd(damnorm$npc_std)*(length(damnorm$npc_std))^(-1/5))

density_plot <- ggplot(damnorm) +
  geom_density(aes(x=log10(npc_norm),
                   color=gas),
               adjust=1,
               linewidth=1.5) +
geom_point(aes(x=log10(npc_norm),
               color=gas,
               y=-(as.numeric(as.factor(gas))-1)/100),
           shape=108) +
  coord_cartesian(xlim = c(-2.5,2.5) ) +
  xlab("Log 10 of normalized present cost") + ylab("Density") +
  ggpubr::theme_pubr(legend="none")

library(gsaot)
gsoat_data <- damnorm %>% ungroup() %>% filter(gas=="ch4") %>% #filter(npc_srm>quantile(npc_srm,0.01)  & npc_srm<quantile(npc_srm,0.99)) %>% 
  select(-gas,-pulse_size,-ozpnpv,-srmpnpv,-masknpv,-dirnpv,-damnpv,-npc_std, -npc_norm, -scc, -scc_srm)

stat_analysis_ch4 <- ot_indices_1d(gsoat_data %>% select(-npc_srm),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_ch4 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

gsoat_data <- damnorm %>% ungroup() %>% filter(gas=="co2") %>%   
  select(-gas,-pulse_size,-ozpnpv,-srmpnpv,-masknpv,-dirnpv,-damnpv,-npc_std, -npc_norm, -scc, -scc_srm)

stat_analysis_co2 <- ot_indices_1d(gsoat_data %>% 
                                     select(-npc_srm),
                                   gsoat_data %>% pull(npc_srm), 
                                   M= 15,
                                   boot = T,
                                   R = 100)
lowerbound_co2 <- irrelevance_threshold(gsoat_data %>% pull(npc_srm), M= 15, solver="1d")

input_categories <- c("ecs"="Climate",
                      "tcr"="Climate",
                      "rcp"="Socio-economic",
                      "pulse_time"="Socio-economic",
                      "geo_start"="SAI policy",
                      "geo_end"="SAI policy",
                      "cool_rate"="SAI policy",
                      "alpha"="Impacts",
                      "delta"="Normative",
                      "theta"="SAI physics",
                      "prob"="SAI policy",
                      "dg"="Socio-economic",
                      "term"="SAI policy",
                      "mortality_ozone"="Impacts",
                      "vsl"="Normative",
                      "vsl_eta"="Normative",
                      "mortality_srm"="Impacts",
                      "forctoTg"="SAI physics")

importance_ch4 <- ggplot( ) + 
  geom_bar(data=tibble("input"=names(stat_analysis_ch4$indices),
                     "original"=stat_analysis_ch4$indices), aes(x=reorder(input, -original),
               y=original,
               fill=input_categories[input]),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_ch4$indices) +
  geom_errorbar(data=as_tibble(stat_analysis_ch4$indices_ci),
                aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [ch4]") + xlab("") + 
  theme(legend.position = "top") + scale_fill_viridis_d(name="") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+ ggpubr::theme_pubr()

importance_co2 <- ggplot( ) + 
  geom_bar(data=tibble("input"=names(stat_analysis_co2$indices),
                       "original"=stat_analysis_co2$indices), aes(x=reorder(input, -original),
                y=original,
                fill=input_categories[input]),stat="identity",color="black") + 
  geom_hline(yintercept=lowerbound_co2$indices) +
  geom_errorbar(data=as_tibble(stat_analysis_co2$indices_ci),
                aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),stat="identity",position="dodge",color="black") +
  ylab("importance [co2]") + xlab("") + 
  theme(legend.position = "top") + scale_fill_viridis_d(name="") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+ ggpubr::theme_pubr()

require(patchwork)
importances <- (importance_co2+importance_ch4)+
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

gsoat_data <- damnorm %>% ungroup() %>% 
  group_by_at(setdiff(all_cols,c("gas","npc_norm")) ) %>%  
  summarise(diff=log10(npc_norm[gas=="co2"])-log10(npc_norm[gas=="ch4"]) ) %>% ungroup()
#  summarise(diff=npc_norm[gas=="co2"]-npc_norm[gas=="ch4"] ) %>% ungroup()

stat_analysis_diff <- ot_indices_1d(gsoat_data %>% 
                                      select(-diff),
                                    gsoat_data %>% pull(diff), 
                                    M= 15,
                                    boot = T,
                                    R = 100)
lowerbound_diff <- irrelevance_threshold(gsoat_data %>% pull(diff), M= 15, solver="1d")

data <- data.frame(i=seq(1,1000),irr=NA,sample=NA)
for(i in seq(1,1000)) {
data[i,"sample"] <- pmin(nrow(gsoat_data),sample(seq(1000,10000,by=100),1))
data[i,"irr"] <- irrelevance_threshold(sample(gsoat_data %>% pull(diff),data[i,"sample"]), M= 15, solver="1d")$indices
data[i,"irr_tot"] <- irrelevance_threshold(gsoat_data %>% pull(diff), M= 15, solver="1d")$indices

 }

hist(gsoat_data %>% pull(diff))
ggplot(data) +
  geom_point(aes(x=sample,y=irr))

importance_diff <- ggplot() + 
  geom_bar(data=tibble("input"=names(stat_analysis_diff$indices),
                       "original"=stat_analysis_diff$indices), aes(x=reorder(input, -original),
                                                                   y=original,
                                                                   fill=input_categories[input]),stat="identity",color="black") + 
#  geom_hline(yintercept=lowerbound_diff$indices) +
  geom_errorbar(data=as_tibble(stat_analysis_diff$indices_ci),
                aes(x=reorder(input, -original),
                    ymin=low.ci,
                    ymax=high.ci),
                stat="identity",position="dodge",color="black") +
  ylab("importance [log10(ch4)-log10(co2)]") + xlab("") + 
  theme(legend.position = "top") + scale_fill_viridis_d(name="") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+ ggpubr::theme_pubr()

damnorm %>% 
  ggplot() + geom_qq(aes(sample=npc_std,color=gas)) +
  geom_abline(slope=1,intercept=0)+ ggpubr::theme_pubr()

dr <- ggplot(damnorm) +
  geom_smooth(aes(x=delta*100,y=log10(npc_norm),color=gas),
              method="loess") +
  xlab("Discount rate [%]") + ylab("Log 10 of normalized cost") +
  coord_cartesian(ylim = c(-0.5,1) ) + ggpubr::theme_pubr(legend="none")

alpha <- ggplot(damnorm) +
  geom_smooth(aes(x=alpha*100,y=log(npc_norm),color=gas),
              method="loess") +
  geom_density(aes(x=alpha*100,y=after_stat(scaled))) +
  xlab("Climate damages [%gdp/K^2]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$alpha*100, 0.01, na.rm=TRUE), 
                           quantile(damnorm$alpha*100, 0.95, na.rm=TRUE)),
                  ylim = c(-0.5,1) ) + 
  ggpubr::theme_pubr(legend="none")

theta <- ggplot(damnorm) +
  geom_smooth(aes(x=theta,y=log10(npc_norm),color=gas),
              method="loess") +
  geom_density(aes(x=theta,y=after_stat(scaled))) +
  xlab("SAI angle [°]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$theta, 0.01, na.rm=TRUE), 
                           quantile(damnorm$theta, 0.99, na.rm=TRUE)),
                  ylim = c(-0.5,1) ) + 
  ggpubr::theme_pubr(legend="none")

term <- ggplot(damnorm) +
  geom_smooth(aes(x=2020+term,y=log10(npc_norm),color=gas),
              method="loess") +
  xlab("Year of termination") + ylab("") +
  coord_cartesian(xlim = c(2020,2400),
                  ylim = c(-0.5,1) ) + 
  ggpubr::theme_pubr(legend="none")

ecs <- ggplot(damnorm) +
  geom_smooth(aes(x=ecs/10,y=log10(npc_norm),color=gas),
              method="loess") +
  geom_density(aes(x=ecs/10,y=after_stat(scaled))) +
  xlab("Climate equilibrium sensitivity [K]") + ylab("") +
  coord_cartesian(xlim = c(quantile(damnorm$ecs/10, 0.01, na.rm=TRUE), 
                           quantile(damnorm$ecs/10, 0.99, na.rm=TRUE)),
                  ylim = c(-0.5,1)) + ggpubr::theme_pubr(legend="none")

void <- ggplot() + theme_void()
fig2 <- (void + density_plot + void + plot_layout(widths = c(0.2,1,0.2)) ) / (dr + term) / (theta + alpha) + plot_layout(heights=c(1,0.6,0.6))

ggsave("fig_2.png",fig2,width=12,height=12,dpi=300)
ggsave("fig_3.png",fig3,width=8,height=7,dpi=300)
ggsave("extfig_gsa.png",importances,width=12,height=6,dpi=300)
ggsave("extfig_gsadiff.png",importance_diff,width=7,height=6,dpi=300)
ggsave("extfig_fra.png",fraction,width=12,height=6,dpi=300)


damnorm %>% filter(gas=="co2" & scc<100) %>% 
  mutate(nosrm = npc_srm / scc, srm = npc_srm / scc_srm ) %>% 
  group_by(gas) %>% 
  filter(nosrm<quantile(nosrm,0.99,na.rm=TRUE) & 
         scc<quantile(scc,0.99,na.rm=TRUE)) %>% 
  ggplot() +
  geom_point(aes(x=scc,y=nosrm,color=gas),alpha=0.1) + 
  geom_point(data=. %>% filter(nosrm<=1),
             aes(x=scc,y=nosrm,color=gas)) +
  facet_wrap(gas~.,scales="free_x")

damnorm %>% 
  mutate(nosrm = npc_srm / scc, srm = npc_srm / scc_srm ) %>% 
  group_by(gas) %>% 
  filter(nosrm<quantile(nosrm,0.95,na.rm=TRUE) & 
         srm<quantile(srm,0.95,na.rm=TRUE) & 
         scc<quantile(scc,0.95,na.rm=TRUE)) %>% 
  select_at(c(all_cols,"srm","nosrm")) %>% 
  pivot_longer(c(nosrm,srm),names_to="fracscc") %>% 
  ggplot() +
  geom_density(aes(x=value,
                   color=gas,
                   linetype=fracscc),
               adjust=10,
               linewidth=1) +
  xlab("Normalized present cost") + ylab("density") 
