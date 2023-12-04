***** emission pulse, co2
active('co2') = yes;
active('ch4') = no;
active('n20') = no;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.fx(ghg,t) = 0;
W_EMI.fx('co2',tfirst) = 37.2; #GtCO2/yr in 2020

solve fair using nlp minimizing OBJ;

execute_unload "emission_pulse_co2_ic%initial_conditions%.gdx";
