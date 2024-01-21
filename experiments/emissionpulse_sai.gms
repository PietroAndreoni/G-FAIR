***** emission pulse, co2
active('co2') = no;
active('ch4') = no;
active('n20') = no;
active('sai') = yes;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.fx(ghg,t) = 0;
W_EMI.fx('sai',tfirst) = 100; #GtCO2/yr in 2020

solve fair using nlp minimizing OBJ;

execute_unload "emission_pulse_sai_ic%initial_conditions%.gdx";
