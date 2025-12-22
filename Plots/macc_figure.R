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
damnpv <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "npc_")), read.csv)) 
scc <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "sccnosrm_")), read.csv)) %>% rename(scc=scc_nosrm)
scc_srm <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "scc_")), read.csv)) %>% rename(scc_srm=scc)
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
            aes(y=miu*100,x=cost,color=as.factor(year)),
            linewidth=1,linetype=2)+
  geom_density(data=data,aes(x=npc_srm,y=after_stat(scaled)*30, color=as.factor(year) ),
               adjust=2,linetype=1, linewidth=1.5) +
  geom_hline(yintercept=0) +
  scale_color_manual(name="",
                     values=c("#6BAED6","#08306B","black")) + 
  geom_point(data=data,
             aes(x=npc_srm,
                 color=as.factor(year),
                 y=0),
             shape=108) +
  theme_classic() + 
  ylab("Emission reductions (% of baseline)\nDensity (scaled, %)") + 
  theme(legend.position = "none") + 
  coord_cartesian(xlim=c(0,25000)) +
  scale_x_continuous(labels = ~paste(., ./25, sep = "\n"),
                     name = "Abatement cost ($/tonCH4)\nAbatement cost ($/tonCO2eq)") #+ facet_wrap(year~.,)

