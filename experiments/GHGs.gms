$set exp %1
$set gas %2


***** emission pulse, co2
W_EMI.fx('%gas%',tsecond) = 2*W_EMI.l('%gas%',tsecond);

$if set sai SRM.lo(t) = -inf; SRM.up(t) = +inf;
target_temp(t) = TATM.l(t);

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'emi') = W_EMI.l(ghg,t)-save_base(ghg,t,'emi');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
$if set sai save_delta(ghg,t,'srm') = SRM.l(t)-save_base(ghg,t,'srm');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');
save_delta(ghg,t,'IRF') = IRF.l(t)-save_base(ghg,t,'IRF');

$if not set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%_GAS%gas%_IC%initial_conditions%";
$if set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%masked_GAS%gas%_IC%initial_conditions%";

$ifthen.trem set tremoval 

* remove the same amount pulsed at t2 at tremoval
W_EMI.fx('%gas%','%tremoval%') = W_EMI.l('%gas%','%tremoval%') - W_EMI.l('%gas%',t)$(ord(t) eq 2);

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');
save_delta(ghg,t,'IRF') = IRF.l(t)-save_base(ghg,t,'IRF');

$if not set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%_REM%tremoval%_GAS%gas%_IC%initial_conditions%";
$if set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%masked_REM%tremoval%_GAS%gas%_IC%initial_conditions%";

$endif.trem

