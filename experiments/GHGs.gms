$set exp %1
$set gas %2
***** emission pulse, co2

$if set sai active('sai') = yes;

$ifthen.exp %exp% =="pulse"
W_EMI.fx('%gas%',tsecond) = W_EMI.l('%gas%',tsecond) + 1e-3$(sameas('%gas%','co2')) + 1$(not sameas('%gas%','co2'));
$elseif.exp %exp% =="const"
W_EMI.fx('%gas%',t)$(ord(t) le 125) = 10*44/12;
$elseif.exp %exp% =="linear"
W_EMI.fx('%gas%',t)$(ord(t) ge 2) = 37 + (ord(t) - 1) * 0.5;
$endif.exp

$if set sai W_EMI.up('sai',t)$(not tfirst(t)) = +inf;
$if set sai FORCING.lo('sai',t) = -inf;
target_temp(t) = TATM.l(t);

solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'emi') = W_EMI.l(ghg,t)-save_base(ghg,t,'emi');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');
save_delta(ghg,t,'IRF') = IRF.l(t)-save_base(ghg,t,'IRF');

$if not set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%_GAS%gas%_IC%initial_conditions%";
$if set sai execute_unload "Results/%rcp%_EXP%experiment_ghg%masked_GAS%gas%_IC%initial_conditions%";

$ifthen.trem set tremoval 
W_EMI.fx('%gas%','%tremoval%') = W_EMI.l('%gas%','%tremoval%') - (1e-6$(sameas('%gas%','co2')) + 1e-3$(not sameas('%gas%','co2')));
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

