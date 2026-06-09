$set gas %1
$setglobal pulse_size 1000 #tons of emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

parameter emi_prepulse(ghg,t);
emi_prepulse(ghg,t) = W_EMI.l(ghg,t);

***** emission pulse
W_EMI.fx('%gas%','%pulse_time%') = W_EMI.l('%gas%','%pulse_time%') + %pulse_size%* (1e-9$(sameas('%gas%','co2')) + 1e-6$(sameas('%gas%','ch4')) + 1e-6$(sameas('%gas%','n2o')));

$ifthen.source %methane_source%=="fossil"
FF_CH4.fx('%pulse_time%') = 1 - ( ( 1 - FF_CH4.l('%pulse_time%')) * emi_prepulse('ch4','%pulse_time%') - %pulse_size%/100 * emissions_rcp('2005','%rcp%','%gas%') ) / W_EMI.l('ch4','%pulse_time%');
$elseif.source %methane_source%=="biogenic"
FF_CH4.fx('%pulse_time%') = 1 - ( ( 1 - FF_CH4.l('%pulse_time%')) * emi_prepulse('ch4','%pulse_time%') + %pulse_size%/100 * emissions_rcp('2005','%rcp%','%gas%') ) / W_EMI.l('ch4','%pulse_time%');
$else.source 
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%');
$endif.source

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXPpulse_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_IC%initial_conditions%.gdx";
