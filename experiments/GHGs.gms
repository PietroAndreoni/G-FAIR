$set exp %1
$set gas %2
***** emission pulse, co2
active('%gas%') = yes;
*active(ghg)$(not sameas(ghg,'%gas%')) = no;

$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc_2020(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = conc_preindustrial(ghg);
$elseif.ic %initial_conditions%=="historical_run"
CONC.fx(ghg,t)$(not active(ghg)) = CONC.l(ghg,'270');
$endif.ic

$ifthen.exp %exp% =="pulse"
W_EMI.fx('%gas%',tsecond) = W_EMI.l('%gas%',tsecond) + 1e-6$(sameas('%gas%','co2')) + 1e-3$(not sameas('%gas%','co2'));
$elseif.exp %exp% =="const"
W_EMI.fx('%gas%',t)$(ord(t) le 125) = 10*44/12;
$elseif.exp %exp% =="linear"
W_EMI.fx('%gas%',t)$(ord(t) ge 2) = 37 + (ord(t) - 1) * 0.5;
$endif.exp


