# Load data
# 100 year GWP (AR4)
# HFCs	HFC-23	HFC-32	HFC-43_10 	HFC-125	HFC-134a	HFC-143a	HFC-152a	HFC-227ea	HFC-236fa	HFC245ca
#       14800	  675	    1640	      3500	  1430	    4470	    124	      3220	    9810	    693
# PFCs + SF6	CF4	  C2F6	SF6	  C6F14
#             7390	12200	22800	9300
#
library(tidyverse)
library(readxl)
library(arrow)

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

# load CH4, N2O MACC (relative to emission factor in 2015)
# costs in USD2010/tCe (using AR4 GWP)
file <- "input/data/1-s2.0-S2352340919306882-mmc1.xlsx"
sheets <- readxl::excel_sheets(file)
macc_sheets <- sheets[str_detect(sheets,"^(SSP2 CH4_|SSP2 N2O_)")]
macc_sheets <- str_remove(macc_sheets,"SSP2 ")
macc_list <- NULL
for(s in macc_sheets) {
  print(s)
  dd <- read_excel(file,sheet = s) |>
    rename(year = t) |>
    rename(cost = `Implementation costs`) |>
    pivot_longer(cols = 3:28, names_to = "region") |>
    mutate(esector = s, .before = year) |>
    separate(esector, into = c("e","sector"), sep = "_")
  
  macc_list <- c(macc_list,list(dd))
  
}
ch4n2o_macc <- bind_rows(macc_list)

# Harmonize names
ch4n2o_macc <- ch4n2o_macc |> 
  mutate(sector = case_when(
    sector == "adip acid" ~ "adipic acid production",
    sector == "ent fermentation" ~ "animals / enteric fermentation",
    sector == "fertilizer" ~ "fertilizer use",
    sector == "manure" ~ "animal waste",
    sector == "ngas" ~ "gas",
    sector == "nitr acid" ~ "nitric acid production",
    sector == "oilp" ~ "oil",
    sector == "rice" ~ "wetland rice production",
    sector == "sewage" ~ "domestic sewage",
    TRUE ~ sector
  ))

# load CH4, N2O sectoral baseline (SSP2)
file <- "input/data/1-s2.0-S2352340919306882-mmc1.xlsx"
sheets <- readxl::excel_sheets(file)
nonco2_baseline <- read_excel(file,sheet = "SSP2 CH4 N2O baseline emissions") |>
    rename(year = Year, region = Region) |>
    pivot_longer(cols = 3:16, names_to = "esector") |>
    mutate(esector = str_replace(esector, " from", "")) |>
    separate(esector, into = c("e","sector"), sep = " ", extra = "merge")

nonco2_baseline <- nonco2_baseline |> 
  mutate(sector = case_when(
    sector == "Domestic sewage" ~ "domestic sewage",
    TRUE ~ sector
  ))

# produce correction factor from % of 2015 activity to % of baseline emissions
correction_factor <- ch4n2o_macc %>%
  filter(cost==0) %>%
  select(-cost) %>% rename(miu0=value) %>%
  full_join(nonco2_baseline %>% rename(base=value)) %>%
  filter(year >= 2015) %>%
  group_by(e,sector,year) %>%
  mutate(wbase = sum(base,na.rm=TRUE)) %>%
  group_by(region,e,sector) %>%
  mutate(share = base/wbase[year==2015] / (1-miu0) )

# convert mius in ton of abated emissions
nonco2_macc <- ch4n2o_macc %>% 
  mutate(unit="activity") %>%
  bind_rows( ch4n2o_macc %>%
               inner_join(correction_factor) %>%
               group_by(region,e,sector,cost) %>%
               mutate(value=(1-value)*share*wbase[year==2015], unit="emissions" ) ) %>% 
  select(e,sector,year,cost,region,unit,value) 

# Load Baseline BASELINE HFCs (ktHFC)
file <- "input/data/1-s2.0-S2352340919306882-mmc2.xlsx"
dd <- read_excel(file,sheet = "Baseline HFC", skip = 3)
hfc_types <- c("hfc23","hfc32","hfc4310","hfc125","hfc134a","hfc143a",
               "hfc152a","hfc227ea","hfc236fa","hfc245ca")
names(dd) <- c("year","region",hfc_types)
baseline_hfc <- dd |>
  pivot_longer(cols = 3:12, names_to = "e") |>
  mutate(sector = "hfc", .before = "value")

# Load Baseline BASELINE PFCs and SF6 (ktPFC)
file <- "input/data/1-s2.0-S2352340919306882-mmc2.xlsx"
dd <- read_excel(file,sheet = "Baseline PFC + SF6", skip = 3)
names(dd) <- c("year","region","cf4","c2f6","sf6","c6f14")
baseline_pfc <- dd |>
  pivot_longer(cols = 3:6, names_to = "e")  |>
  mutate(sector = "pfc", .before = "value")

# Merge baselines, all species in MtCeq (CH4,N2O) / ktCeq
nonco2_baseline <- bind_rows(nonco2_baseline, 
                             baseline_hfc %>% mutate(value=value*ar4gwp[e]/3.66), 
                             baseline_pfc %>% mutate(value=value*ar4gwp[e]/3.66)) 

# Reductions in Kt HFC
# costs in USD2005/tCe (using AR4 GWP)
file <- "input/data/1-s2.0-S2352340919306882-mmc2.xlsx"
cc_hfc <- read_excel(file,sheet = "cost_curve_absolute_HFC_all", skip = 2) |>
  rename(year = 1, e = 2, cost = 3) |>
  mutate(e = hfc_types[e]) |>
  select(-`class_ 27`,-`class_ 28`) |>
  pivot_longer(cols = 4:29, names_to = "region") |>
  mutate(sector = "hfc")

# Reductions in Kt HFC
# costs in USD2005/tCe (using AR4 GWP)
file <- "input/data/1-s2.0-S2352340919306882-mmc2.xlsx"
pfc_types <- c("cf4","c2f6","c6f14")
cc_pfc <- read_excel(file,sheet = "cost_curve_absolute_PFC_all", skip = 3) |>
  rename(year = 1, e = 2, cost = 3) |>
  mutate(e = pfc_types[e]) |>
  select(-`class_ 27`,-`class_ 28`) |>
  pivot_longer(cols = 4:29, names_to = "region") |>
  mutate(sector = "pfc")

# Reductions in % of baseline
# costs in USD2005/tCe (using AR4 GWP)
file <- "input/data/1-s2.0-S2352340919306882-mmc2.xlsx"
cc_sf6 <- read_excel(file,sheet = "cost_curve_relative_SF6", skip = 2) |>
  rename(year = 1, cost = 2) |>
  select(-`Global average`) |>
  pivot_longer(cols = 3:28, names_to = "region") |>
  mutate(sector = "sf6") |>
  mutate(e = "sf6")

nonco2_macc <- bind_rows(nonco2_macc, 
                         cc_hfc %>% mutate(unit="emissions", 
                                           value = value*ar4gwp[e]/3.66), 
                         cc_pfc %>% mutate(unit="emissions", 
                                           value = value*ar4gwp[e]/3.66), 
                         cc_sf6 %>% mutate(unit="emissions") %>% 
                           inner_join(baseline_pfc %>% 
                                        rename(base=value) %>% 
                                        filter(e=="sf6")) %>% 
                           mutate(value=value/100*base) %>% select(-base)) 

# Add region code
image_regions <- tibble(image26 = unique(witchtools::region_mappings[['image26']]$image26)) |>
  mutate(region = case_when(
    image26 == "bra" ~ "Brazil",             
    image26 == "can" ~ "Canada",
    image26 == "rcam" ~ "Central America",
    image26 == "ceu" ~ "Central Europe",
    image26 == "chn" ~ "China+",
    image26 == "eaf" ~ "East Africa",
    image26 == "india" ~ "India",
    image26 == "indo" ~ "Indonesia",
    image26 == "jap" ~ "Japan",
    image26 == "stan" ~ "Kazachstan",      
    image26 == "kor" ~ "Korea",
    image26 == "mex" ~ "Mexico",
    image26 == "me" ~ "Middle East",
    image26 == "naf" ~ "North Africa",
    image26 == "oce" ~ "Oceania",
    image26 == "rsaf" ~ "Rest of South Africa",
    image26 == "rsas" ~ "Rest of South Asia",   
    image26 == "rsam" ~ "Rest of South-America",
    image26 == "rus" ~ "Russia",
    image26 == "saf" ~ "South Africa",      
    image26 == "seas" ~ "South East Asia",
    image26 == "tur" ~ "Turkey",
    image26 == "ukr" ~ "Ukraine",
    image26 == "usa" ~ "USA",
    image26 == "waf" ~ "West Africa",
    image26 == "weu" ~ "West Europe",
    TRUE ~ "ERROR"
  ))

nonco2_macc <- image_regions |>
  full_join(nonco2_macc)

nonco2_baseline <- image_regions |>
  full_join(nonco2_baseline)

write_parquet(nonco2_macc, "input/data/harmsen_nonco2_macc.parquet")
write_parquet(nonco2_baseline, "input/data/harmsen_nonco2_baseline.parquet")
