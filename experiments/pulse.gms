$set gas %1

$setglobal pulse_size 1 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

parameter dT(t),dW(t,ghg);

dT(t) = TATM.l(t);
dW(t,ghg) = W_EMI.l(ghg,t);

***** emission pulse, co2
W_EMI.fx('%gas%','%pulse_time%') = W_EMI.l('%gas%','%pulse_time%') + %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

dT(t) = TATM.l(t) - dT(t);
dW(t,ghg) = W_EMI.l(ghg,t) - dW(t,ghg);

execute_unload "Results/%rcp%_EXPpulse_GAS%gas%_IC%initial_conditions%_PT%pulse_time%_PS%pulse_size%.gdx";
