$set gas %1
$setglobal rate_of_cooling 10 #in deg/millennium, by default half of current warming (0.1 deg/decade). Weird unit is for reporting 
$setglobal start_rampdown "na"
$setglobal end_rampdown "na"
$setglobal rate_of_cooling 0
$setglobal pulse_size 10 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

** build SRM strategy (linear ramp-up and down, flat in between)
forcing_srm(t)$(2020 + t.val ge %start_rampup% and 2020 + t.val le %end_rampup%) = %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampup%) / 1e3 ; # 0.1 deg/decade
forcing_srm(t)$(2020 + t.val gt %end_rampup% and 2020 + t.val le %start_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3;
forcing_srm(t)$(2020 + t.val gt %start_rampdown% and 2020 + t.val le %end_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3 - %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampdown%) / 1e3;
forcing_srm(t)$(forcing_srm(t) le 0) = 0;
solve fair using nlp minimizing OBJ;

** set target for experiment with SRM masking (3rd run)
target_temp(t) = TATM.l(t);

** Scenario 1: with SRM, without emission pulse
execute_unload "%results_folder%/%rcp%_EXP%experiment%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

** Scenario 2: with SRM, with emission pulse (but no masking)
W_EMI.fx('%gas%','%pulse_time%') = W_EMI.l('%gas%','%pulse_time%') + %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

solve fair using nlp minimizing OBJ;

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulse_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

** Scenario 3: with SRM, with emission pulse and masking
SRM.lo(t) = - forcing_srm(t); SRM.up(t) = +inf; # reduce the SAI by at most the full amount
SRM.fx(t)$(t.val lt %pulse_time%) = 0; # and cannot change SAI before the pulse

* solve 3 times as this is the only non-simulation problem
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulsemasked_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

$ifthen.tterm set termination_time 

* fix strategy befre tterm
SRM.fx(t)$(t.val le %termination_time%) = SRM.l(t); 
* stop masking at tterm
SRM.fx(t)$(t.val gt %termination_time%) = 0; 

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulsemaskedterm_TER%termination_time%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

$endif.tterm

$ifthen.trem set removal_time 

* remove the same amount pulsed at t2 at tremoval
W_EMI.fx('%gas%','%removal_time%') = W_EMI.l('%gas%','%removal_time%') - %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Pulse model is not solving";

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulsemaskedrem_REM%removal_time%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

$endif.trem

