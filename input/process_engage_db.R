# Packages
require(witchtools)
pkgs <- c('MonoPoly','colf','tidyverse','data.table','gdxtools')
res <- lapply(pkgs, require_package)

exclude_removal <- TRUE
load_raw_data <- FALSE
compute_removal_curves <- TRUE

#-------------------------------------- PART 1 ----------------------------------------

#----------------------- 1.1) Build EMISSION tables AR6 ---------------------
#if required, recreate from this file (in the backup RICE50+ folder): engage_raw_data.parquet
if(load_raw_data==TRUE) {
engage_raw_data <- arrow::read_parquet(file.path('input','data','engage_raw_data.parquet')) %>% 
  dplyr::select(Model,Scenario,Region,Variable,Unit,paste(seq(2010,2100,by = 5)) ) %>% 
  pivot_longer(c(paste(seq(2010,2100,by = 5))), names_to = "year")

data_emissions_ar6 <- engage_raw_data %>% filter(Variable %in% c("Price|Carbon",
                                                              "Emissions|CH4",
                                                              "Emissions|N2O",
                                                              "Emissions|F-Gases",
                                                              "Diagnostics|MAGICC6|Temperature|Global Mean") |
                                                str_detect(Variable, "^Emissions\\|CO2") )

data_removal_ar6 <- engage_raw_data %>% filter(str_detect(Variable,"Carbon Sequestration") &
                                        !str_detect(Variable,"CCS\\|Fossil") &
                                        !str_detect(Variable,"Land Use"))

data_gdp_ar6 <- engage_raw_data %>% filter(Variable %in% c("GDP|PPP","Population"))

arrow::write_parquet(data_emissions_ar6, file.path('Data','Input','data_engage_ar6.parquet'))
arrow::write_parquet(data_removal_ar6, file.path('Data','Input','data_removal_ar6.parquet')) 
arrow::write_parquet(data_gdp_ar6, file.path('Data','Input','data_gdp_ar6.parquet')) 
}

emissions <- arrow::read_parquet(file.path('input','data','data_engage_ar6.parquet'))
removal <- arrow::read_parquet(file.path('input','data','data_removal_ar6.parquet'))

map_removal_to_emissions <- c("Carbon Sequestration|CCS|Biomass"="Emissions|CO2|Energy and Industrial Processes",
                              "Carbon Sequestration|Direct Air Capture"="Emissions|CO2|Energy and Industrial Processes",
                              "Carbon Sequestration|Enhanced Weathering"="Emissions|CO2|Energy and Industrial Processes",
                              "Carbon Sequestration|CCS|Biomass|Energy|Supply|Electricity"="Emissions|CO2|Energy|Supply|Electricity",
                              "Carbon Sequestration|CCS|Biomass|Energy|Demand|Industry"="Emissions|CO2|Energy|Demand|Industry",
                              "Carbon Sequestration|CCS|Biomass|Energy|Supply|Gases"="Emissions|CO2|Energy|Supply|Gases",
                              "Carbon Sequestration|CCS|Biomass|Energy|Supply|Liquids"="Emissions|CO2|Energy|Supply|Liquids",
                              "Carbon Sequestration|CCS|Biomass|Energy|Supply|Hydrogen"="Emissions|CO2|Energy|Supply|Liquids", #arbitrary, but ending up in other energy transformation
                              "Carbon Sequestration|CCS|Biomass|Energy|Supply|Other"="Emissions|CO2|Energy|Supply|Other sector", #arbitrary, but ending up in other energy transformation
                              "Carbon Sequestration|Direct Air Capture"="Emissions|CO2|Industrial Processes") # this is really arbitrary!

total_removal <- removal %>%
  mutate(value=ifelse(str_detect(Model,"REMIND") & str_detect(Variable,"Direct Air Capture"),-value,value)) %>%
  filter(Variable %in% names(map_removal_to_emissions) ) %>%
  mutate(Variable=map_removal_to_emissions[Variable]) %>% 
  filter(!is.na(value)) %>%
  group_by(Model,Scenario,Region,Variable,Unit,year) %>%
  summarise(removal=sum(value))
  

sectors_to_emissions <- c("Total_CO2"="Emissions|CO2|Energy and Industrial Processes",
                          "Power"="Emissions|CO2|Energy|Supply|Electricity",
                          "Industry_fuel"="Emissions|CO2|Energy|Demand|Industry",
                          "Industry_process"="Emissions|CO2|Industrial Processes",
                          "Transport"="Emissions|CO2|Energy|Demand|Transportation",
                          "Buildings"="Emissions|CO2|Energy|Demand|Residential and Commercial",
                          "Agriculture"="Emissions|CO2|Energy|Demand|AFOFI",
                          "Other_energy_transformation"="Emissions|CO2|Energy|Supply|Liquids",
                          "Other_energy_transformation"="Emissions|CO2|Energy|Supply|Gases",
                          "Other_energy_transformation"="Emissions|CO2|Energy|Supply|Heat",
                          "Other_energy_transformation"="Emissions|CO2|Energy|Supply|Solids",
                          "Other_energy_transformation"="Emissions|CO2|Energy|Supply|Other sector",
                          "Methane"="Emissions|CH4")

emissions_to_sectors <- setNames(names(sectors_to_emissions), sectors_to_emissions)

if(exclude_removal==TRUE) {
gross_emissions <- emissions %>% 
  filter(Variable %in% unname(sectors_to_emissions) | Variable=="Emissions|CH4" & !is.na(value)) %>%
  full_join(total_removal) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  mutate(value=value+removal) %>% 
  ungroup() %>% mutate(Variable = emissions_to_sectors[Variable]) %>%
  group_by(Model,Scenario,Region,Variable,Unit,year) %>%
  summarise(value=sum(value,na.rm=TRUE)) %>%
  ungroup() 
} else {
gross_emissions <- emissions %>% 
  filter(Variable %in% unname(sectors_to_emissions) | Variable=="Emissions|CH4" & !is.na(value)) %>%
  ungroup() %>% mutate(Variable = emissions_to_sectors[Variable]) %>%
  group_by(Model,Scenario,Region,Variable,Unit,year) %>%
  summarise(value=sum(value,na.rm=TRUE)) %>%
  ungroup() 
}

emi_cprice_all <- gross_emissions %>%
  rename(Emissions=value,Sector=Variable) %>%
  full_join(emissions %>% 
              filter(Variable=="Price|Carbon") %>%
              select(-Variable,-Unit) %>%
              rename(Cprice=value)) 

baseline_scenario_ar6 = 'EN_NPi2100' # = Current Policy scenario

# Calculate abatement_perc(=miu) ((emi_bau - emi) / emi_bau)
macc_with_bau <- emi_cprice_all %>%
  group_by(Model,year,Region,Sector) %>%
  filter(any(Scenario==baseline_scenario_ar6) ) %>% 
  group_by(Model,year,Region,Sector) %>%
  mutate(abatement_perc=(Emissions[Scenario==baseline_scenario_ar6]-Emissions)/Emissions[Scenario==baseline_scenario_ar6] ) %>%
  filter(Scenario!='EN_NoPolicy' & str_detect(Scenario,"f") )

# Prepare regional mapping for the R5 to R10 disaggregation of POLES emissions and Cprice
r5_r10_region_mapping <- data.frame(r5 = c(toupper(c("r5lam", "r5asia", "r5maf", "r5oecd", "r5ref", "r5maf", "r5asia", "r5oecd", "r5oecd", "r5asia")),"World"),
  r10 = c("R10LATIN_AM", "R10INDIA+", "R10AFRICA", "R10EUROPE", "R10REF_ECON", "R10MIDDLE_EAST", "R10REST_ASIA", "R10PAC_OECD", "R10NORTH_AM", "R10CHINA+","World"))

macc_with_bau <- macc_with_bau %>%
  filter(Region %in% setdiff(unique(r5_r10_region_mapping$r5),"World") ) %>%
  full_join(r5_r10_region_mapping %>% rename(Region=r5), relationship="many-to-many") %>%
  mutate(Region=r10) %>%
  select(-r10) %>%
  rbind(macc_with_bau %>%
          filter(Region %in% unique(r5_r10_region_mapping$r10)))

# function for fitting POLYNOMIAL124 CONSTRAINED:  y ~ ax + bx^2 + cx^4   (forcing Y positive, monotonically increasing and convex)
constrained_poly1246_fit <- function(dd){
  colf_nlxb(y ~ miu + I(miu^2) + I(miu^4) + I(miu^6) - 1, data = dd, lower = c( 0, 0, 0, 0)) 
}

constrained_poly124_fit <- function(dd){
  colf_nlxb(y ~ miu + I(miu^2) + I(miu^4) - 1, data = dd, lower = c( 0, 0, 0)) 
}

constrained_poly14_fit <- function(dd){
  colf_nlxb(y ~ miu + I(miu^4) - 1, data = dd, lower = c( 0, 0)) 
}

fit_macc_by_model <- macc_with_bau %>% 
  as_tibble() %>%
  rename(y=Cprice,miu=abatement_perc) %>%
  filter(!is.na(miu) & !is.na(y) & y<=1000 & miu>=0 & miu<=1 & year >= 2025) %>% # exlcude high carbon prices and miu
  group_by(Sector,Model,Region,year) %>%
  do(value=coef(constrained_poly1246_fit(.))) %>% # extrapolate the best-fitting coefficients of the curve
  unnest_wider(value) %>% 
  rename(c1=param_miu,c2=param_I.miu.2.,c4=param_I.miu.4.,c6=param_I.miu.6.) 

fit_macc_by_model <- fit_macc_by_model %>%
  pivot_longer(c(c1,c2,c4,c6)) %>%
  group_by(Sector,Model,Region,name) %>%
  do(year = seq(2025,2100,by=5), value = Hmisc::approxExtrap(.$year, .$value, xout = seq(2025,2100,by=5), method = "linear")$y ) %>% # linear extrapolation
  ungroup() %>% 
  unnest(c(year,value)) %>%
  pivot_wider() 

macc_by_model <- data.frame(miu=seq(0,1,by=0.02)) %>%
  cross_join(fit_macc_by_model) %>%
  mutate(cprice=c1*miu+c2*miu^2+c4*miu^4+c6*miu^6)

fit_macc_by_percentile <- macc_by_model %>% 
  mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                              year>2035 & year <= 2055 ~ "2040-2055",
                              year>2055 & year <= 2075 ~ "2060-2080",
                              year>2075 ~ "2080-2100")) %>%
  group_by(Region,Sector,yearbundle,miu) %>%
  reframe(Quantile=c("prob10","prob25","prob50","prob75","prob90"),
          y=unname(quantile(cprice,c(0.1,0.25,0.5,0.75,0.9)))) %>%
  group_by(Sector,Region,Quantile,yearbundle) %>%
  do(value=coef(constrained_poly1246_fit(.))) %>%
  unnest_wider(value) %>% 
  rename(c1=param_miu,c2=param_I.miu.2.,c4=param_I.miu.4.,c6=param_I.miu.6.,limits10=Region) %>%
  mutate(limits10=tolower(limits10))
  
macc_by_percentile <- data.frame(miu=seq(0,1,by=0.02)) %>%
  cross_join(fit_macc_by_percentile) %>%
  mutate(cprice=c1*miu+c2*miu^2+c4*miu^4+c6*miu^6) %>%
  select(Sector,limits10,Quantile,yearbundle,miu,cprice)

### Figure 1: calibration of hard-to-abate in 2050
ggplot() + 
  geom_point(data=macc_with_bau %>% 
               mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                                           year>2035 & year <= 2055 ~ "2040-2055",
                                           year>2055 & year <= 2075 ~ "2060-2080",
                                           year>2075 ~ "2080-2100")) %>%
               filter(Region=="World" & Sector=="Total_CO2"),
             aes(x=abatement_perc,y=Cprice) ) +
  geom_line(data=macc_by_percentile %>% 
               filter(limits10=="world" & Sector=="Total_CO2"),
             aes(x=miu,y=cprice,color=Quantile) ) + 
  ylim(c(0,1000)) + xlim(c(0,1)) +
  xlab("Fraction of controlled emissions") + ylab("Carbon price [$/tonCO2]") +
  facet_wrap(yearbundle~.,) + ggpubr::theme_pubr()

ggplot() + 
  geom_point(data=macc_with_bau %>% 
               mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                                           year>2035 & year <= 2055 ~ "2040-2055",
                                           year>2055 & year <= 2075 ~ "2060-2080",
                                           year>2075 ~ "2080-2100")) %>%
               filter(Region=="World" & Sector=="Methane"),
             aes(x=abatement_perc,y=Cprice*25*16/32) ) +
  geom_line(data=macc_by_percentile %>% 
              filter(limits10=="world" & Sector=="Methane"),
            aes(x=miu,y=cprice*25*16/32,color=Quantile) ) + 
  ylim(c(0,1000*25*16/32)) + xlim(c(0,1)) +
  xlab("Fraction of controlled emissions") + ylab("Carbon price [$/tonCh4]") +
  facet_wrap(yearbundle~.,) + ggpubr::theme_pubr()


ggplot() + 
  geom_point(data=macc_with_bau %>% 
               mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                                           year>2035 & year <= 2055 ~ "2040-2055",
                                           year>2055 & year <= 2075 ~ "2060-2080",
                                           year>2075 ~ "2080-2100")) %>%
               filter(Region=="World" & Sector=="Methane"),
             aes(x=abatement_perc,y=Cprice), color="red" ) +
  geom_point(data=macc_with_bau %>% 
               mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                                           year>2035 & year <= 2055 ~ "2040-2055",
                                           year>2055 & year <= 2075 ~ "2060-2080",
                                           year>2075 ~ "2080-2100")) %>%
               filter(Region=="World" & Sector=="Total_CO2"),
             aes(x=abatement_perc,y=Cprice), color="blue" ) +
  ylim(c(0,1000)) + xlim(c(0,1)) +
  xlab("Fraction of controlled emissions") + ylab("Carbon price [$/tonCh4]") +
  facet_grid(yearbundle~Model,) + ggpubr::theme_pubr()
  
ggplot() + 
geom_line(data=macc_by_percentile %>% 
            filter(limits10=="world" & Sector=="Total_CO2" & Quantile=="prob50"),
          aes(x=miu,y=cprice,color=yearbundle,group=yearbundle) ) 

ggplot() + 
  geom_point(data=macc_with_bau %>% 
               filter(year %in% c(2030,2100) & Region=="World" & Sector=="Total_CO2"),
             aes(x=abatement_perc,y=Cprice) ) +
  geom_line(data=macc_by_model %>% 
              filter(year %in% c(2030,2100) & Region=="World" & Sector=="Total_CO2"),
            aes(x=miu,y=cprice,color=Model) ) + 
  ylim(c(0,1000)) + xlim(c(0,1)) +
  xlab("Fraction of controlled emissions") + ylab("Carbon price [$/tonCO2]") +
  facet_wrap(year~.,) + theme_pubr()

if(compute_removal_curves==TRUE) {

constrained_poly0124_fit <- function(dd){
  colf_nlxb(y ~ miu + I(miu^2) + I(miu^4), data = dd, lower = c(0, 0, 0, 0)) }

constrained_poly02_fit <- function(dd){
  colf_nlxb(y ~ I(miu^2), data = dd, lower = c( 0, 0)) 
}

removal_types <- c("Carbon Sequestration|CCS|Biomass"="BECCS",
                   "Carbon Sequestration|Direct Air Capture"="DACCS",
                   "Carbon Sequestration|Enhanced Weathering"="EW")
  
rem_cprice_all <- removal %>%
  filter(Variable %in% names(removal_types) ) %>%
  mutate(value=ifelse(str_detect(Model,"REMIND") & str_detect(Variable,"Direct Air Capture"),-value,value)) %>%
  mutate(Variable=removal_types[Variable]) %>%
  rbind(removal %>%
            filter(Variable %in% names(removal_types) ) %>% 
           group_by(Model,Scenario,Region,Unit,year) %>%
           summarise(value=sum(value),Variable="Total_cdr") ) %>% 
  rename(Removal=value,Sector=Variable) %>%
  full_join(emissions %>% 
              filter(Variable=="Price|Carbon") %>%
              select(-Variable,-Unit) %>%
              rename(Cprice=value)) %>%
  filter( !str_detect(Scenario,"f") & str_detect(Scenario,"NPi2020_") )

rem_cprice_all <- rem_cprice_all %>%
  filter(Region %in% setdiff(unique(r5_r10_region_mapping$r5),"World") ) %>%
  full_join(r5_r10_region_mapping %>% rename(Region=r5), relationship="many-to-many") %>%
  mutate(Region=r10) %>% 
  select(-r10) %>%
  rbind(rem_cprice_all %>%
          filter(Region %in% unique(r5_r10_region_mapping$r10)))

fit_rem_by_model <- rem_cprice_all %>% 
  as_tibble() %>%
  rename(y=Cprice,miu=Removal) %>%
  filter(!is.na(miu) & !is.na(y) & y<=1000 & miu>0 & year >= 2025) %>% # exlcude high carbon prices and miu
  group_by(Sector,Model,Region,year) %>%
  do(value=coef(constrained_poly02_fit(.))) %>% # extrapolate the best-fitting coefficients of the curve
  unnest_wider(value) %>% 
  rename(c0=param_X.Intercept.,c2=param_I.miu.2.) 

fit_rem_by_model <- fit_rem_by_model %>%
  pivot_longer(c(c0,c2)) %>%
  group_by(Sector,Model,Region,name) %>%
  do(year = seq(2025,2100,by=5), value = Hmisc::approxExtrap(.$year, .$value, xout = seq(2025,2100,by=5), method = "linear")$y ) %>% # linear extrapolation
  ungroup() %>% 
  unnest(c(year,value)) %>%
  pivot_wider() 

rem_by_model <- data.frame(miu=seq(0,20000,by=200)) %>%
  cross_join(fit_rem_by_model) %>%
  mutate(cprice=c0+c2*miu^2)

fit_rem_by_percentile <- rem_by_model %>%
  mutate(yearbundle=case_when(year<=2035 ~ "2025-2035",
                              year>2035 & year <= 2055 ~ "2040-2055",
                              year>2055 & year <= 2075 ~ "2060-2080",
                              year>2075 ~ "2080-2100")) %>%
  group_by(Region,Sector,yearbundle,miu) %>%
  reframe(Quantile=c("prob10","prob25","prob33","prob50","prob66","prob75","prob90"),
          y=unname(quantile(cprice,c(0.1,0.25,0.33,0.5,0.66,0.75,0.9)))) %>%
  group_by(Sector,Region,Quantile,yearbundle) %>%
  do(value=coef(constrained_poly02_fit(.))) %>%
  unnest_wider(value) %>% 
  rename(c0=param_X.Intercept.,c2=param_I.miu.2.) %>%
  rename(limits10=Region) %>%
  mutate(limits10=tolower(limits10))

rem_by_percentile <- data.frame(miu=seq(0,20000,by=50)) %>%
  cross_join(fit_rem_by_percentile) %>%
  mutate(cprice=c0+c2*miu^2) %>%
  select(Sector,limits10,Quantile,yearbundle,miu,cprice)

}



gdp <- arrow::read_parquet(file.path('Data','Input','data_gdp_ar6.parquet'))

baseline_socioecon <- gdp %>% select(-Unit) %>%
  filter(Scenario==baseline_scenario_ar6 & year>=2020 & Model %in% c("WITCH 5.0","IMAGE 3.0","MESSAGEix-GLOBIOM_1.1","REMIND-MAgPIE 2.1-4.2") ) %>%
  pivot_wider(names_from="Variable",values_from="value") %>%
  inner_join(emissions %>% 
               filter(Variable=="Emissions|CO2|Energy and Industrial Processes") %>% 
               rename(Emissions=value) %>% 
               select(-Variable,-Unit)) %>%
  pivot_longer(c("Emissions","GDP|PPP","Population")) %>%
  group_by(Model,Region,name) %>%
  do(year = seq(2020,2100,by=5), value = Hmisc::approxExtrap(.$year, .$value, xout = seq(2020,2100,by=5), method = "linear")$y ) %>% # linear extrapolation
  ungroup() %>% 
  unnest(c(year,value)) %>%
  group_by(Region,year,name) %>%
  summarise(value=median(value,na.rm=TRUE)) %>% 
  group_by(Region,name) %>%
  summarise(growth = mean((lead(value)-value)/lead(value)/5, na.rm=TRUE)) %>%
  mutate(limits10=tolower(Region)) %>% select(-Region)
  

#save GDX files
invisible(write.gdx(file.path('Data', 'data_mod_macc.gdx'), 
                    params = list(macc_ed_coef=rbind(fit_macc_by_model %>% rename(Specification=Model,limits10=Region), 
                                  fit_macc_by_percentile %>% rename(Specification=Quantile)),
                                  baseline_socioecon %>% rename(Specification=Quantile)) ) ) 
