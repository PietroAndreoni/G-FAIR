# Usage
require(witchtools)

pkgs <- c('tidyverse','stringr','countrycode','arrow','data.table')
res <- lapply(pkgs,require_package)
require_gdxtools()
igdx(dirname(Sys.which('gams'))) # Please have gams in your PATH!

# load data from Harmsen (provided in 2010 $/tonCeq)
baseline <- read_parquet("input/data/harmsen_nonco2_baseline.lz4.parquet")
macc <- read_parquet("input/data/harmsen_nonco2_macc.lz4.parquet")

output_folder <- "Results_fig3"
damnpv <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "npc_output")), read.csv)) 
scc <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "sccnosrm_output")), read.csv)) %>% rename(scc=scc_nosrm)
scc_srm <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "scc_output")), read.csv)) %>% rename(scc_srm=scc)
all_cols <- names(damnpv)[1:19]

# filter scc and npc>0 and scenarios with both CH4 and CO2
damnpv <- damnpv %>% 
  inner_join(scc %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  inner_join(scc_srm %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  filter(npc_srm>0 & scc>0 & scc_srm>0) %>% 
  ungroup() %>% unique()

# enerdata co2 (provided in 2015 $/tonCO2)
macc_co2 <- read_parquet("input/data/macc_ed_full_2022.lz4.parquet") %>% 
  filter(Variable=="Emissions") %>% 
  group_by(Year,Scenario,Carbon_value) %>%
  summarise(value=sum(Value)) %>% 
  group_by(Year,Scenario) %>%
  mutate(miu=(value[Carbon_value==0]-value)/value[Carbon_value==0] ) %>% 
  rename(cost=Carbon_value,year=Year) %>% 
  select(Scenario,year,cost,miu)
  

map_sectors_to_ipcc <- c("fossil"="coal",
                         "fossil"="oil",
                         "fossil"="gas",
                         "agriculture"="animals / enteric fermentation",
                         "agriculture"="fertilizer use",
                         "agriculture"="wetland rice production",
                         "agriculture"="animal waste",
                         "waste"="landfills",
                         "waste"="domestic sewage",
                         "industry"="adipic acid production",
                         "industry"="nitric acid production",
                         "industry"="pfc",
                         "industry"="hfc",
                         "transport"="transport")

ar4gwp <- c("c2f6"=12200,
            "c6f14"=9300,
            "cf4"=7390,
            "hfc125"=124,   
            "hfc134a"=1430,
            "hfc143a"=4470,
            "hfc152a"=124,
            "hfc227ea"=3220,
            "hfc23"=14800,
            "hfc236fa"=675,
            "hfc245ca"=693,
            "hfc32"=675,   
            "hfc4310"=1640,
            "ch4"=25,
            "n2o"=298)  

correction_factor <- macc %>%
  filter(cost==0) %>%
  select(-cost) %>% rename(miu0=value) %>%
  full_join(baseline %>% rename(base=value)) %>%
  filter(year>=2015 & e %in% c("CH4","N2O") ) %>% 
  group_by(e,sector,year) %>%
  mutate(wbase = sum(base,na.rm=TRUE)) %>%
  group_by(image26,region,e,sector) %>%
  mutate(share = base/wbase[year==2015] / (1-miu0) )

macc_by_gas_w <- macc %>% 
  inner_join(correction_factor) %>%
  group_by(image26,e,sector,cost) %>%
  mutate(value=(1-value)*share*wbase[year==2015] ) %>%
  group_by(year,e,cost) %>%
  summarize(value=sum(value,na.rm=TRUE),
            base=sum(base,na.rm=TRUE)) %>%
  ungroup() %>% mutate(miu=(base-value)/base, e=tolower(e)) %>%
  mutate(cost=cost*ar4gwp[e]*12/44) %>%
  select(year,e,cost,miu)

data <- damnpv %>% 
  filter(gas=="ch4") %>%
  mutate(year=pulse_time+2020) %>% 
  filter(year %in% c(2025) )

fig3 <- ggplot() +
  geom_line(data=macc_by_gas_w %>% 
              filter(e=="ch4" & year %in% c(2025) ),
            aes(y=miu*100,x=cost), 
            color="black",linewidth=1,linetype=2)+
  # geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv)/pulse_size,y=after_stat(scaled)*10 ),
  #              adjust=2,linetype=1, linewidth=1.5, color="#6BAED6") +
  # geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv)/pulse_size,y=after_stat(scaled)*10 ),
  #              adjust=2,linetype=1, linewidth=1.5, color="#08306B") +
  geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv+ozpnpv)/pulse_size,y=after_stat(scaled)*10, color=as.factor(year) ),
               adjust=2,linetype=1, linewidth=1.5, color="black") +
  geom_hline(yintercept=0) +
  geom_point(data=data,
             aes(x=npc_srm,
                 y=0),
             shape=108, 
             color="black") +
  theme_classic() + 
  ylab("Emission reductions (% of baseline)\nDensity (scaled to 100)") + 
  theme(legend.position = "none") + 
  coord_cartesian(xlim=c(0,25000)) +
  scale_x_continuous(labels = ~paste(., ./25, sep = "\n"),
                     name = "Abatement cost ($/tonCH4)\nAbatement cost ($/tonCO2eq)") #+ facet_wrap(year~.,)


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

