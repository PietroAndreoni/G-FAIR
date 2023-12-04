***** 1: constconc, co2
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

CONC.fx('co2',t) = conc0('co2');
W_EMI.lo('co2',t) = -100; #GtCO2/yr in 2020
W_EMI.up('co2',t) = +100; #GtCO2/yr in 2020

solve fair using nlp minimizing OBJ;

execute_unload "const_conc_co2_ic%initial_conditions%.gdx";
