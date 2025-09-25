$set gas %1

$setglobal start_rampup 2030
$setglobal end_rampup 2100
$setglobal start_rampdown 2200
$setglobal end_rampdown 2300
$setglobal rate_of_cooling 10 #in deg/millennium, by default half of current warming (0.1 deg/decade). Weird unit is for reporting 
$setglobal pulse_size 1 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

** build SRM strategy (linear ramp-up and down, flat in between)
forcing_srm(t)$(2020 + t.val ge %start_rampup% and 2020 + t.val le %end_rampup%) = %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - 2030) / 1e3 ; # 0.1 deg/decade
forcing_srm(t)$(2020 + t.val gt %end_rampup% and 2020 + t.val le %start_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - 2030) / 1e3;
forcing_srm(t)$(2020 + t.val gt %start_rampdown% and 2020 + t.val le %end_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - 2030) / 1e3 - %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampdown%) / 1e3;
forcing_srm(t)$(forcing_srm(t) le 0) = 0;

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_base(ghg,t,'emi') = W_EMI.l(ghg,t);
save_base(ghg,t,'conc') = CONC.l(ghg,t);
save_base(ghg,t,'forc') = FORCING.l(ghg,t);
save_base(ghg,t,'T') = TATM.l(t);
save_base(ghg,t,'IRF') = IRF.l(t);
save_base(ghg,t,'srm') = forcing_srm(t)+SRM.l(t);


***** emission pulse, co2
W_EMI.fx('%gas%','%pulse_time%') = (1 + %pulse_size%/100 ) *W_EMI.l('%gas%','%pulse_time%');

$if set srm SRM.lo(t) = -inf; SRM.up(t) = +inf;
target_temp(t) = TATM.l(t);

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'emi') = W_EMI.l(ghg,t)-save_base(ghg,t,'emi');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
$if set srm save_delta(ghg,t,'srm') = forcing_srm(t)+SRM.l(t)-save_base(ghg,t,'srm');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');
save_delta(ghg,t,'IRF') = IRF.l(t)-save_base(ghg,t,'IRF');

$if not set srm execute_unload "Results/%rcp%_EXP%experiment%_GAS%gas%_PT%pulse_time%_PS%pulse_size%_RC%rate_of_cooling%_IC%initial_conditions%";
$if set srm execute_unload "Results/%rcp%_EXP%experiment%masked_GAS%gas%_PT%pulse_time%_PS%pulse_size%_RC%rate_of_cooling%_IC%initial_conditions%";

$ifthen.trem set removal_time 

* remove the same amount pulsed at t2 at tremoval
W_EMI.fx('%gas%','%removal_time%') = W_EMI.l('%gas%','%removal_time%') - %pulse_size%/100 * W_EMI.l('%gas%','%pulse_time%');

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');
save_delta(ghg,t,'IRF') = IRF.l(t)-save_base(ghg,t,'IRF');

$if not set srm execute_unload "Results/%rcp%_EXP%experiment_ghg%_GAS%gas%_REM%removal_time%_PT%pulse_time%_PS%pulse_size%_RC%rate_of_cooling%_IC%initial_conditions%";
$if set srm execute_unload "Results/%rcp%_EXP%experiment_ghg%masked_REM%removal_time%_GAS%gas%_REM%removal_time%_PT%pulse_time%_PS%pulse_size%_RC%rate_of_cooling%_IC%initial_conditions%";

$endif.trem

