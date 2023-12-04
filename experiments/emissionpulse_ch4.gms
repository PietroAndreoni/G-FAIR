***** 1bis: emission pulse, ch4
active('co2') = no;
active('ch4') = yes;
active('n2o') = no;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
CONC.lo(ghg,t)$(not tfirst(t)) = 0;
CONC.up(ghg,t)$(not tfirst(t)) = +inf;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.fx(ghg,t) = 0;
W_EMI.fx('ch4',tfirst) = 37200/28; #GtCO2/yr in 2020

solve fair using nlp minimizing OBJ;

execute_unload "emission_pulse_ch4_ic%initial_conditions%.gdx";
