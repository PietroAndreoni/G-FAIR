# Usage
require(witchtools)

pkgs <- c('tidyverse','stringr','countrycode','arrow','data.table')
res <- lapply(pkgs,require_package)
require_gdxtools()
igdx(dirname(Sys.which('gams'))) # Please have gams in your PATH!

# load data from Harmsen (provided in 2010 $/tonCeq)
baseline <- read_parquet("input/data/harmsen_nonco2_baseline.lz4.parquet")
macc <- read_parquet("input/data/harmsen_nonco2_macc.lz4.parquet")

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

fig3 <- ggplot() +
  geom_line(data=macc_by_gas_w %>% filter(e=="ch4" & year %in% c(2025,2050) ),
            aes(x=miu*100,y=cost/ar4gwp["ch4"],linetype=as.factor(year)),linewidth=1.5,color="red")+
  geom_line(data=macc_co2 %>% filter(year %in% c(2025,2050) & Scenario=="EnerBase" ),
            aes(x=miu*100,y=cost,linetype=as.factor(year)),linewidth=1.5,color="lightblue")+
  geom_hline(yintercept=0) +
  scale_color_manual(values=c("#6BAED6","#08306B")) + 
  theme_classic() + 
  ylab("Abatement cost ($/tonCH4)") + xlab("Emission reductions (% of baseline)") +
  theme(legend.position = "top") + xlim(c(0,100))
ggsave("figure_3.svg",width=5.5,height=5.5,path=respath,plot=fig3)


fig3 <- ggplot() +
  geom_line(aes(x=miu*100,y=cost,color=as.factor(year)),linewidth=1.5)+
  geom_hline(yintercept=0) +
  scale_color_manual(values=c("#6BAED6","#08306B")) + 
  theme_classic() + 
  ylab("Abatement cost ($/tonCH4)") + xlab("Emission reductions (% of baseline)") +
  theme(legend.position = "none")
