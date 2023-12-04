***** steady state
active('co2') = yes;
active('ch4') = yes;
active('n20') = yes;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.fx(ghg,t) = 0;

solve fair using nlp minimizing OBJ;

execute_unload "steady_state_ic%initial_conditions%.gdx";
