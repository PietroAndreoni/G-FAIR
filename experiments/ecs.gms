**** 2: equilibrium climate sensitivity
active('co2') = yes;
active('ch4') = no;
active('n2o') = no;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
CONC.lo(ghg,t)$(not tfirst(t)) = 0;
CONC.up(ghg,t)$(not tfirst(t)) = +inf;

$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$endif.ic

W_EMI.fx('co2',tfirst) =  (2*conc0('co2')-CONC.fx('co2','1')) / (1e18 / atmosphere_mass * ghg_mm('co2')  / atmosphere_mm); #GtCO2/yr in 2020

solve fair using nlp minimizing OBJ;

execute_unload "ECS_ic%initial_conditions%.gdx";