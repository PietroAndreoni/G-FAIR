require(dplyr)
require(stringr)
require(janitor)

# extract emissions
files <- paste0("data/", list.files(path="./data",pattern = "\\.csv$")) %>%
  str_subset("EMISSIONS") 
emissions_rcps <- tibble()
for (file in files){
emissions_rcps <- emissions_rcps %>%
  readr::bind_rows(read_csv(file) %>%
  slice(36:n()) %>%
  row_to_names(1) %>%
  rename(year=`v YEARS/GAS >`) %>%
  clean_names() %>%
  mutate_all(~as.numeric(.)) %>%
  mutate(rcp=str_remove_all(file,"data/|_EMISSIONS.csv")) ) 
}

#units for emissions
units_emissions <- readr::read_csv(files[1]) %>%
  slice(36:35) %>%
  row_to_names(1) %>%
  select(-1) %>%
  clean_names() %>%
  select(-fossil_co2, -other_co2 ) %>%
  mutate(co2 = "GtC/yr")

emissions_rcps <- emissions_rcps %>%
  mutate(co2 = fossil_co2   + other_co2 ) %>%
  select(-fossil_co2, -other_co2 ) %>%
  pivot_longer(!c(year,rcp),names_to="ghg") 
  
# extract forcing 
files <- paste0("data/", list.files(path="./data",pattern = "\\.csv$")) %>%
  str_subset("RADFORCING") 
forcing_rcps <- tibble()
for (file in files){
  forcing_rcps <- forcing_rcps %>%
    bind_rows(readr::read_csv(file) %>%
                slice(58:n()) %>%
                row_to_names(1) %>%
                rename(year=`v YEARS/GAS >`) %>%
                clean_names() %>%
                mutate_all(~as.numeric(.)) %>%
                mutate(rcp=str_remove_all(file,"data/|_MIDYEAR_RADFORCING.csv")) ) 
}

forcing_rcps <- forcing_rcps %>%
  pivot_longer(!c(year,rcp),names_to="source") %>%
  mutate(source=str_remove(source,"_rf"))

# extract fraction of methane fossil emissions 
methane_fraction <- paste0("data/", list.files(path="./data",pattern = "\\.csv$")) %>%
  str_subset("fossilCH4_fraction") %>%
  readr::read_csv() %>%
  slice(4:n()) %>%
  row_to_names(1) %>%
  mutate_all(~as.numeric(.)) %>%
  pivot_longer(!c(Year),names_to="rcp") %>%
  mutate(rcp=str_remove_all(rcp,"[^[:alnum:]]")) %>%
  rename(year=Year)
  

natural_emissions <- paste0("data/", list.files(path="./data",pattern = "\\.csv$")) %>%
  str_subset("natural") %>%
  readr::read_csv() %>% 
  slice(4:n()) %>%
  separate(col=1,into=c("year","ch4","n2o"), sep="[ \t]+") %>%
  mutate_all(~as.numeric(.)) %>%
  pivot_longer(!c(year),names_to="ghg") 

### save into gdx file
require(gdxtools)
write.gdx("RCPs_consolidated.gdx",params=list("Forcing"=forcing_rcps,
                                              "Emissions"=emissions_rcps,
                                              "fossilch4_frac"=methane_fraction,
                                              "natural_emissions"=natural_emissions))

          