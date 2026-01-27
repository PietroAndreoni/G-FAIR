$set gas %1
$setglobal rate_of_cooling 10 #in deg/millennium, by default half of current warming (0.1 deg/decade). Weird unit is for reporting 
$setglobal start_rampdown "na"
$setglobal end_rampdown "na"
$setglobal rate_of_cooling 0
$setglobal pulse_size 100 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)
$setglobal methane_source "mix" # "fossil", "biogenic", "mix"

* save initial state 
parameter wemi_0(ghg,t), ffch4_0(t);
wemi_0(ghg,t) = W_EMI.l(ghg,t);
ffch4_0(t) = FF_CH4.l(t);

** Scenario 0: no SRM, emission pulse

* pulse emissions
W_EMI.fx('%gas%','%pulse_time%') = W_EMI.l('%gas%','%pulse_time%') + %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

$ifthen.source %methane_source%=="fossil"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%') + %pulse_size%/100 * ( 0 + ( emissions_rcp('2005','%rcp%','%gas%') / W_EMI.l('ch4','%pulse_time%') )$(sameas('%gas%','ch4')) );
$elseif.source %methane_source%=="biogenic"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%') - %pulse_size%/100 * ( 0 + ( emissions_rcp('2005','%rcp%','%gas%') / W_EMI.l('ch4','%pulse_time%') )$(sameas('%gas%','ch4')) );
$elseif.source %methane_source%=="mix"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%');
$endif.source

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXPpulse_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_IC%initial_conditions%.gdx";

* remove pulse and adjust methane fraction
W_EMI.fx(ghg,t) = wemi_0(ghg,t);
FF_CH4.fx(t) = ffch4_0(t);

** Scenario 1: with SRM, without emission pulse

* build SRM strategy (linear ramp-up and down, flat in between)
forcing_srm(t)$(2020 + t.val ge %start_rampup% and 2020 + t.val le %end_rampup%) = %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampup%) / 1e3 ; # 0.1 deg/decade
forcing_srm(t)$(2020 + t.val gt %end_rampup% and 2020 + t.val le %start_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3;
forcing_srm(t)$(2020 + t.val gt %start_rampdown% and 2020 + t.val le %end_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3 - %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampdown%) / 1e3;
forcing_srm(t)$(forcing_srm(t) le 0) = 0;
solve fair using nlp minimizing OBJ;

* set target for experiment with SRM masking (3rd run)
target_temp(t) = TATM.l(t);

execute_unload "%results_folder%/%rcp%_EXP%experiment%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

** Scenario 2: with SRM, with emission pulse (but no masking)
W_EMI.fx('%gas%','%pulse_time%') = W_EMI.l('%gas%','%pulse_time%') + %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

$ifthen.source %methane_source%=="fossil"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%') + %pulse_size%/100 * ( 0 + emissions_rcp('2005','%rcp%','%gas%') / W_EMI.l('ch4','%pulse_time%') )$(sameas('%gas%','ch4'));
$elseif.source %methane_source%=="biogenic"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%') - %pulse_size%/100 * ( 0 + emissions_rcp('2005','%rcp%','%gas%') / W_EMI.l('ch4','%pulse_time%') )$(sameas('%gas%','ch4'));
$elseif.source %methane_source%=="mix"
FF_CH4.fx('%pulse_time%') = FF_CH4.l('%pulse_time%');
$endif.source


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

** Scenario 4: with SRM, with emission pulse and masking, terminated at termination_time
$ifthen.tterm set termination_time 

* fix strategy befre tterm
SRM.fx(t)$(t.val le %termination_time%) = SRM.l(t); 
* stop masking at tterm
SRM.fx(t)$(t.val gt %termination_time%) = 0; 

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulsemaskedterm_TER%termination_time%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

$endif.tterm

** Scenario 5 (optional): with SRM, with emission pulse removed at time removal_time
$ifthen.trem set removal_time 

* remove the same amount pulsed at t2 at tremoval
W_EMI.fx('%gas%','%removal_time%') = W_EMI.l('%gas%','%removal_time%') - %pulse_size%/100 *  ( (emissions_rcp('2005','%rcp%','%gas%')/CO2toC)$(sameas('%gas%','co2')) + emissions_rcp('2005','%rcp%','%gas%')$(not sameas('%gas%','co2')));

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Pulse model is not solving";

execute_unload "%results_folder%/%rcp%_EXP%experiment%pulsemaskedrem_REM%removal_time%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%_IC%initial_conditions%";

$endif.trem

