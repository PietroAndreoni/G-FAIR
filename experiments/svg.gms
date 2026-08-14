$setglobal rate_of_cooling 10 #in deg/millennium, by default half of current warming (0.1 deg/decade). Weird unit is for reporting 
$setglobal start_rampup 2025
$setglobal end_rampup 2100
$setglobal start_rampdown 2300
$setglobal end_rampdown 2400
$setglobal pulse_size 100 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

* save initial state 
parameter wemi_0(ghg,t), ffch4_0(t);
wemi_0(ghg,t) = W_EMI.l(ghg,t);
ffch4_0(t) = FF_CH4.l(t);

** Scenario 1: with SRM, without emission pulse

* build SRM strategy (linear ramp-up and down, flat in between)
background_srm(t)$(2020 + t.val ge %start_rampup% and 2020 + t.val le %end_rampup%) = %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampup%) / 1e3 ; # 0.1 deg/decade
background_srm(t)$(2020 + t.val gt %end_rampup% and 2020 + t.val le %start_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3;
background_srm(t)$(2020 + t.val gt %start_rampdown% and 2020 + t.val le %end_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3 - %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampdown%) / 1e3;
background_srm(t)$(background_srm(t) le 0) = 0;
solve fair using nlp minimizing OBJ;

target_temp(t) = TATM.l(t);

execute_unload "%results_folder%/%rcp%_EXP%experiment%baseline_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%";

** Scenario 2: with SRM, with emission pulse (but no masking)
W_EMI.fx('co2','%pulse_time%') = W_EMI.l('co2','%pulse_time%') + %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC);

solve fair using nlp minimizing OBJ;

execute_unload "%results_folder%/%rcp%_EXP%experiment%ghgpulse_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%";

** Scenario 3: with additional SRM, no emission pulse
SRM.lo(t) = - background_srm(t); SRM.up(t) = +inf; # reduce the SAI by at most the full amount
SRM.fx(t)$(t.val lt %pulse_time%) = 0; # and cannot change SAI before the pulse

* solve 3 times as this is the only non-simulation problem
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

SRM.fx(t) = SRM.l(t); 
W_EMI.fx(ghg,t) = wemi_0(ghg,t);

solve fair using nlp minimizing OBJ;

execute_unload "%results_folder%/%rcp%_EXP%experiment%srmpulse_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%";
