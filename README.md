# G-FAIR

This folder contains a GAMS translation of FAIR v1.3, containing all relevant equations pertaining to CO2, CH4 and N2O.
Other greenhouses gases are exogenous. This version of the model also contains a representation of solar geoengineering and can solve 
For the original FAIR v1.3 model, see: https://gmd.copernicus.org/articles/11/2273/2018/ 

# VERSION 1.0 

For referees: this version contains the replication package for the paper "Solar geoengineering is a better substitute for methane than carbon" (In review as of March 17 2025).
To replicate the relevant results:
(1) run input/consolidate_data_ar5.R in the input folder. 
(2) run "run_experiments.bat" in the main folder or, from a GAMS command line (e.g. GAMS studio):
--experiment_ghg=pulse --gas=ch4 --sai=1
--experiment_ghg=pulse --gas=co2 --sai=1
--experiment_ghg=pulse --gas=ch4 
--experiment_ghg=pulse --gas=co2 
Results replication requires an active license to run. Results were obtained with GAMS 40 and CONOPT.
(3) To replicate the plots and data used in the paper, run "Plots/main_figures.R" in the Plots folder.
Data analysis replications does not require any active license to run. Plots were obtained with RStudio 2023.03.0+386 "Cherry Blossom" Release and packages "gdxtools" (v0.7.0) and "tidyverse" (v2.0.0).
Data can be replicated directly as relevant result files are included in the V1.0 release.
