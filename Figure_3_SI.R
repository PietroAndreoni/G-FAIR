# Usage
require(witchtools)

pkgs <- c('tidyverse','stringr','countrycode','arrow','data.table')
res <- lapply(pkgs,require_package)
require_gdxtools()
igdx(dirname(Sys.which('gams'))) # Please have gams in your PATH!

# load data from Harmsen (provided in 2010 $/tonCeq)
baseline <- read_parquet("input/data/harmsen_nonco2_baseline.lz4.parquet")
macc <- read_parquet("input/data/harmsen_nonco2_macc.lz4.parquet")

output_folders <- tibble(
  output_folder = list.files(pattern = "^Results_base_1903", full.names = TRUE)
) %>%
  filter(dir.exists(output_folder)) %>%
  arrange(output_folder) %>%
  mutate(
    folder_name = basename(output_folder),
    angle = case_when(
      str_detect(folder_name, "freeangle") ~ "free",
      str_detect(folder_name, "angle\\d+") ~ str_extract(folder_name, "(?<=angle)\\d+"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(angle))

angle_levels <- output_folders %>%
  pull(angle) %>%
  unique()

read_output_files <- function(folders, pattern) {
  files_to_read <- folders %>%
    mutate(file = map(output_folder, ~ list.files(path = .x, pattern = pattern, full.names = TRUE))) %>%
    unnest(file)

  out <- files_to_read %>%
    mutate(data = map2(file, angle, ~ read.csv(.x) %>% mutate(angle = .y))) %>%
    pull(data) %>%
    bind_rows() %>%
    mutate(angle = factor(angle, levels = angle_levels))

  out
}

damnpv <- read_output_files(output_folders, "npc_output")
scc <- read_output_files(output_folders, "sccnosrm_output") %>% rename(scc = scc_nosrm)
scc_srm <- read_output_files(output_folders, "scc_output") %>% rename(scc_srm = scc)
all_cols <- c(names(damnpv)[1:20], "angle")
remove_outliers <- c("ecs","tcr","alpha","theta","mortality_srm","mortality_ozone","dg","TgtoUSD","forctoTg")
check_densities <- damnpv %>%
  filter(gas=="ch4") %>%
  filter_at(remove_outliers, ~ (.x <= quantile(.x, 0.999, na.rm = TRUE) & .x >= quantile(.x, 0.001, na.rm = TRUE) ) ) %>%
  select_at(setdiff(all_cols,c("gas") ))
check_densities %>%
  select_at(remove_outliers) %>%
  pivot_longer(remove_outliers) %>%
  ggplot() +
  geom_density(aes(x=value) ) +
  facet_wrap(name~.,scales="free")

# filter scc and npc>0 and scenarios with both CH4 and CO2
damnpv <- damnpv %>%
  inner_join(check_densities) %>%
#  inner_join(scc %>% select(-damnpv,-ozpnpv,-pulse_time)) %>%
#  inner_join(scc_srm %>% select(-damnpv,-ozpnpv,-pulse_time)) %>%
  filter(npc_srm>0) %>%
  ungroup() %>% unique()

damnpv %>%
  group_by(pulse_time, angle) %>%
  summarise(med=median(npc_srm),
            p95=quantile(npc_srm,.95),
            .groups = "drop")

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
  filter(year %in% c(2025,2050) )

density_cost_unit <- 100

fig3_si <- ggplot() +
  # geom_line(data=macc_by_gas_w %>%
  #             filter(e=="ch4" & year %in% c(2025,2050) ),
  #           aes(y=miu*100,x=cost),
  #           color="black",linewidth=1,linetype=2)+
  geom_density(data=data,
               aes(x=(masknpv+damnpv+dirnpv+srmpnpv)/pulse_size,
                   y=after_stat(density * 100 * density_cost_unit),
                   color=angle,
                   group=angle),
               adjust=2, linewidth=1, linetype=2) +
  # geom_density(data=data,aes(x=(masknpv+damnpv+dirnpv+srmpnpv)/pulse_size,y=after_stat(scaled)*10 ),
  #              adjust=2,linetype=1, linewidth=1.5, color="#08306B") +
  geom_density(data=data,
               aes(x=(masknpv+damnpv+dirnpv+srmpnpv+ozpnpv)/pulse_size,
                   y=after_stat(density * 100 * density_cost_unit),
                   color=angle,
                   group=angle),
               adjust=2, linewidth=1.5, linetype=1) +
  geom_hline(yintercept=0) +
  geom_point(data=data,
             aes(x=npc_srm,
                 color=angle,
                 y=0),
             shape=108) +
  facet_wrap(year~.,) +
  theme_classic() +
  ylab("Emission reductions (% of baseline)\nDensity (%)") +
  theme(legend.position = "top") +
  coord_cartesian(xlim=c(0,10000)) +
  scale_linetype_discrete(name = "Angle") +
  scale_x_continuous(labels = ~paste(., ./25, sep = "\n"),
                     name = expression(atop("Abatement cost ($/ton" * CH[4] * ")",
                                            "Abatement cost ($/ton" * CO[2] * "eq)"))) #+ facet_wrap(year~.,)

ggsave("fig_3_SI.png",fig3_si,width=12,height=6,dpi=300)

