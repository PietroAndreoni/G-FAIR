# Usage
require(witchtools)

pkgs <- c('tidyverse','stringr','countrycode','arrow','data.table')
res <- lapply(pkgs,require_package)
require_gdxtools()
igdx(dirname(Sys.which('gams'))) # Please have gams in your PATH!

# Single control file for plotting settings / input files / unit conversions.
# Locate the Paper_SAI folder (holds all_parameters.R) robustly so the script
# works under Rscript (--file), RStudio "Source" (sys.frame $ofile), and an
# interactive console whose working dir is at/under/above the project.
.find_paper_root <- function() {
  starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
              unlist(lapply(sys.frames(), function(f) f$ofile)), getwd())
  for (s in starts[nzchar(starts)]) {
    d <- if (dir.exists(s)) s else dirname(s)
    repeat {
      if (file.exists(file.path(d, "all_parameters.R")))
        return(normalizePath(d, "/", FALSE))
      if (file.exists(file.path(d, "Paper_SAI", "all_parameters.R")))
        return(normalizePath(file.path(d, "Paper_SAI"), "/", FALSE))
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Cannot locate all_parameters.R (Paper_SAI control file); set the working ",
       "directory to the project root or the Paper_SAI folder.", call. = FALSE)
}
source(file.path(.find_paper_root(), "all_parameters.R"))

# load data from Harmsen (provided in 2010 $/tonCeq)
baseline <- read_parquet(HARMSEN_BASELINE_FIG3)
macc <- read_parquet(HARMSEN_MACC_FIG3)

output_folder <- RESULTS_FOLDER_FIG3
damnpv <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "npc_output")), read.csv)) 
scc <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "sccnosrm_output")), read.csv)) %>% rename(scc=scc_nosrm)
scc_srm <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "scc_output")), read.csv)) %>% rename(scc_srm=scc)
all_cols <- names(damnpv)[1:21]
remove_outliers <- FIG_OUTLIER_COLS
check_densities <- damnpv %>%
  filter(gas=="ch4") %>%
  filter_at(remove_outliers, ~ (.x <= quantile(.x, FIG_OUTLIER_QHI, na.rm = TRUE) & .x >= quantile(.x, FIG_OUTLIER_QLO, na.rm = TRUE) ) ) %>%
  select_at(setdiff(all_cols,c("gas") ))
check_densities %>% 
  select_at(remove_outliers) %>% 
  pivot_longer(remove_outliers) %>% 
  ggplot() +
  geom_density(aes(x=value) ) +
  facet_wrap(name~.,scales="free")

# filter scc and npc>0 and scenarios with both CH4 and CO2
damnpv <- damnpv %>% 
#  inner_join(check_densities) %>% 
#  inner_join(scc %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
#  inner_join(scc_srm %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  filter(npc_srm>0) %>% 
  ungroup() %>% unique() 

damnpv %>% 
  group_by(pulse_time) %>% 
  summarise(med=median(npc_srm),
            p95=quantile(npc_srm,.95))

# enerdata co2 (provided in 2015 $/tonCO2)
macc_co2 <- read_parquet(MACC_CO2_FILE) %>%
  filter(Variable=="Emissions") %>% 
  group_by(Year,Scenario,Carbon_value) %>%
  summarise(value=sum(Value)) %>% 
  group_by(Year,Scenario) %>%
  mutate(miu=(value[Carbon_value==0]-value)/value[Carbon_value==0] ) %>% 
  rename(cost=Carbon_value,year=Year) %>% 
  select(Scenario,year,cost,miu)
  
macc_by_gas_w <- macc %>% 
  filter(unit=="emissions") %>% 
  full_join(baseline %>% rename(base=value) ) %>%
  group_by(year,e,cost) %>%
  summarize(value=sum(value,na.rm=TRUE),
            base=sum(base,na.rm=TRUE)) %>%
  ungroup() %>% mutate(miu=(base-value)/base, e=tolower(e)) %>%
  mutate(cost=cost*AR4_GWP100[e]*C_PER_CO2*USD_DEFLATOR_2010_2020) %>%
  select(year,e,cost,miu)

data <- damnpv %>%
  filter(gas=="ch4") %>%
  mutate(year=pulse_time+2020) %>%
  filter(year %in% FIG3_MACC_YEARS )

density_cost_unit <- DENSITY_COST_UNIT

fig3 <- ggplot() +
  geom_line(data=macc_by_gas_w %>%
              filter(e=="ch4" & year %in% FIG3_MACC_YEARS ),
            aes(y=miu*100,x=cost),
            color="black",linewidth=1,linetype=2)+
  geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv)/pulse_size,
                             y=after_stat(density * 100 * density_cost_unit) ),
               adjust=2,linetype=1, linewidth=1, color="#6BAED6") +
  # geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv)/pulse_size,y=after_stat(scaled)*10 ),
  #              adjust=2,linetype=1, linewidth=1.5, color="#08306B") +
  geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv+ozpnpv)/pulse_size,
                             y=after_stat(density * 100 * density_cost_unit),
                             color=as.factor(year) ),
               adjust=2,linetype=1, linewidth=1.5, color="#08306B") +
  geom_hline(yintercept=0) +
  geom_point(data=data,
             aes(x=npc_srm,
                 y=0),
             shape=108, 
             color="#08306B") +
  facet_wrap(year~.,) +
  theme_classic() + 
  ylab("Emission reductions (% of baseline)\nDensity (% per $100/ton)") + 
  theme(legend.position = "none") + 
  coord_cartesian(xlim=FIG3_XLIM,ylim=FIG3_YLIM) +
  scale_x_continuous(labels = ~paste(., ./CH4_GWP100, ./CH4_GWP100/C_PER_CO2,sep = "\n"),
                     name = expression(atop("Abatement cost ($/ton" * CH[4] * ")",
                                            "Abatement cost ($/ton" * CO[2] * "eq)",
                                            "Abatement cost ($/tonCeq)"))) #+ facet_wrap(year~.,)
save_figure("fig_3.png",fig3,width=12,height=6,dpi=300)

frac_damages <- damnpv %>% 
  pivot_longer(c(dirnpv,srmpnpv,ozpnpv,masknpv,damnpv),
               names_to="source",values_to="npc_partial") %>% 
  group_by_at(all_cols) %>% 
  mutate(npc_partial=npc_partial/pulse_size ) %>% 
  group_by_at(all_cols) %>%
  mutate(npc_frac=npc_partial/sum(npc_partial,na.rm=TRUE))


fraction <- frac_damages %>% filter(npc_frac!=0  & npc_frac>0) %>% 
  ggplot()+
  geom_density(aes(x=npc_partial,y=after_stat(scaled),color=source),adjust=1,linewidth=1.5) +
  facet_wrap(gas~.,) + ggpubr::theme_pubr() +
  xlab("Fraction of total cost") + ylab("Density") +
coord_cartesian(xlim=c(0,3000)) 

### co2
# macc_co2 %>% filter(year %in% c(2025,2050) & Scenario=="EnerBase" ) %>%
#   ungroup() %>% 
#   select(year,miu,cost) %>% 
#   full_join(damnorm %>% filter(gas=="co2" ) %>% mutate(cost=case_when( npc_srm <= 20 ~ round(npc_srm / 5) * 5,
#                                                                        npc_srm <=110 ~ round(npc_srm / 10) * 10,
#                                                                        npc_srm <=520 ~ round(npc_srm / 20) * 20,
#                                                                        npc_srm <=925 ~ round(npc_srm / 40) * 40,
#                                                                        npc_srm <=1000 ~ round(npc_srm / 50) * 50,
#                                                                        .default= 1000 ) ) ) %>% 
#   ggplot() +
#   geom_density(aes(x=miu,color=as.factor(year)),adjust=5)
# 
# macc_by_gas_w %>% filter(e=="ch4" & year %in% c(2025,2050) ) %>%
#   ungroup() %>% 
#   select(year,miu,cost) %>% mutate(cost=round(cost/25*3.66  / 20) * 20) %>% 
#   full_join(damnorm %>% filter(gas=="ch4" ) %>% mutate(cost=pmax(round(npc_srm / 20) * 20), 4000 )) %>% ggplot() +
#   geom_density(aes(x=miu,color=as.factor(year)),adjust=5)
# 
# 
# data <- damnorm %>% 
#   filter(gas=="co2") %>%
#   mutate(year=pulse_time+2020) %>% 
#   mutate(year=ifelse(year==2022,2025,year) ) %>% 
#   filter(year %in% c(2025,2050))
# 
# fig3_co2 <- ggplot() +
#   # geom_line(data=macc_by_gas_w %>% filter(e=="ch4" & year %in% c(2025,2050) ),
#   #           aes(y=miu*100,x=cost,color=as.factor(year)),linewidth=1.5)+
#   geom_line(data=macc_co2 %>% filter(year %in% c(2025,2050) & Scenario=="EnerBase" ),
#             aes(y=miu*100,x=cost,color=as.factor(year)),linewidth=1.5)+
#   geom_density(data=data,aes(x=npc_srm,y=after_stat(scaled)*100,color=as.factor(year)),adjust=5) +
#   geom_point(data=data,
#              aes(x=npc_srm,
#                  y=0),
#              shape=108) +
#   geom_hline(yintercept=0) +
#   scale_color_manual(values=c("#6BAED6","#08306B")) + 
#   theme_classic() + 
#   ylab("Emission reductions (% of baseline) / density (scaled)") + xlab("Abatement cost ($/tonCH4)") +
#   theme(legend.position = "top") + xlim(c(0,1000))

