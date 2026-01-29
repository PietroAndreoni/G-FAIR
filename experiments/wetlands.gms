$setglobal pulse_size 100 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)
$setglobal proj_duration 20 #projection duration in years
$setglobal CO2toCH4_long 0.1 #fraction of ch4 emitted per unit of co2 removed (long term)
$setglobal CO2toCH4_pulse 0.25 #fraction of ch4 emitted per unit of co2 removed (short term)

parameter emi_prepulse(ghg,t);
emi_prepulse(ghg,t) = W_EMI.l(ghg,t);

***** emission pulse
W_EMI.fx('co2',t)$(t.val ge %pulse_time% and t.val lt %pulse_time%  + %proj_duration%) = W_EMI.l('co2',t) - %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC);
W_EMI.fx('ch4',t)$(t.val eq %pulse_time%) = W_EMI.l('ch4',t) + %CO2toCH4_pulse% * %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC) * 16/44 ;
W_EMI.fx('ch4',t)$(t.val gt %pulse_time% and t.val lt %pulse_time%  + %proj_duration%) = W_EMI.l('ch4',t) + %CO2toCH4_long% * %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC) * 16/44;

* methane emissions are biogenic
FF_CH4.fx(t)$(t.val eq %pulse_time%) = 1 - ((1-FF_CH4.l(t)) * emi_prepulse('ch4',t) + %CO2toCH4_pulse% * %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC) * 16/44 ) / W_EMI.l('ch4',t);
FF_CH4.fx(t)$(t.val gt %pulse_time% and t.val lt %pulse_time% + %proj_duration%) = 1 - ((1-FF_CH4.l(t)) * emi_prepulse('ch4',t) + %CO2toCH4_long% * %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC) * 16/44 ) / W_EMI.l('ch4',t);

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXPwetlands_ECS%ecs%_TCR%tcr%_PT%pulse_time%_PROJ%proj_duration%_PR%CO2toCH4_pulse%_LR%CO2toCH4_long%.gdx";
