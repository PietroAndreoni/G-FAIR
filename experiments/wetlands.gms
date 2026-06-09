$setglobal pulse_time 5 #year of restoration (1 is 2020)
$setglobal proj_duration 40 #projection duration in years
$setglobal wetland_area_km2 1000 #wetland area restored in square kilometers
$setglobal CO2perkm 2000 #peak CO2 uptake in tCO2 per km2 per year
$setglobal CH4perkm 80 #peak CH4 emissions in tCH4 per km2 per year
$setglobal carbon_peak_time 5 #years after restoration when carbon uptake peaks
$setglobal carbon_tau 30 #time constant for carbon uptake decline after peak
$setglobal carbon_baseline_frac 0.33 #long-run carbon uptake as fraction of peak uptake
$setglobal methane_tau 5 #time constant for methane emissions decline after restoration
$setglobal methane_baseline_frac 0.2 #long-run methane emissions as fraction of peak emissions
$setglobal methane_emissions yes #whether to include methane emissions from wetlands

scalar wetland_peak_co2_uptake;
scalar wetland_peak_ch4_emissions;

parameter emi_prepulse(ghg,t);
parameter wetland_age(t);
parameter wetland_carbon_profile(t);
parameter wetland_methane_profile(t);
parameter wetland_co2_uptake(t);
parameter wetland_ch4_emissions(t);

emi_prepulse(ghg,t) = W_EMI.l(ghg,t);
wetland_age(t) = t.val - %pulse_time%;
wetland_peak_co2_uptake = %wetland_area_km2% * %CO2perkm% / 1e9;
wetland_peak_ch4_emissions = %wetland_area_km2% * %CH4perkm% / 1e6;

***** wetland restoration emission profiles
wetland_carbon_profile(t) = 0;
wetland_carbon_profile(t)$(wetland_age(t) ge 0 and wetland_age(t) le %carbon_peak_time% and wetland_age(t) lt %proj_duration%) =
        wetland_age(t) / %carbon_peak_time%;
wetland_carbon_profile(t)$(wetland_age(t) gt %carbon_peak_time% and wetland_age(t) lt %proj_duration%) =
        %carbon_baseline_frac% + (1 - %carbon_baseline_frac%) * exp(-(wetland_age(t) - %carbon_peak_time%) / %carbon_tau%);

wetland_co2_uptake(t) =
        wetland_peak_co2_uptake * wetland_carbon_profile(t);
W_EMI.fx('co2',t) = emi_prepulse('co2',t) - wetland_co2_uptake(t);

$ifthen.m "%methane_emissions%"=="yes"
wetland_methane_profile(t) = 0;
wetland_methane_profile(t)$(wetland_age(t) ge 0 and wetland_age(t) lt %proj_duration%) =
        %methane_baseline_frac% + (1 - %methane_baseline_frac%) * exp(-wetland_age(t) / %methane_tau%);
wetland_ch4_emissions(t) =
        wetland_peak_ch4_emissions * wetland_methane_profile(t);

W_EMI.fx('ch4',t) = emi_prepulse('ch4',t) + wetland_ch4_emissions(t);

* methane emissions from wetlands are biogenic
FF_CH4.fx(t)$((emi_prepulse('ch4',t) + wetland_ch4_emissions(t)) gt 0) =
        1 - ((1-FF_CH4.l(t)) * emi_prepulse('ch4',t) + wetland_ch4_emissions(t))
          / (emi_prepulse('ch4',t) + wetland_ch4_emissions(t));
$endif.m

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXPwetlands_ECS%ecs%_TCR%tcr%_PT%pulse_time%_PROJ%proj_duration%_AREA%wetland_area_km2%_CO2KM%CO2perkm%_CH4KM%CH4perkm%_CP%carbon_peak_time%_CT%carbon_tau%_CB%carbon_baseline_frac%_MT%methane_tau%_MB%methane_baseline_frac%.gdx";
