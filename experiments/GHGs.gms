$set exp %1
$set gas %2

***** emission pulse, co2
active('%gas%') = yes;
active(ghg)$(not sameas(ghg,'%gas%')) = no;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$elseif.ic %initial_conditions%=="historical_run"
CONC.fx(ghg,t)$(not active(ghg)) = CONC.l(ghg,'255');
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.fx(ghg,t) = 0; #GtCO2/yr in 2020

$ifthen.exp %exp% =="pulse"
W_EMI.fx('%gas%',tsecond) = 37;
$elseif.exp %exp% =="const"
W_EMI.fx('%gas%',t)$(ord(t) le 125) = 10*44/12;
$elseif.exp %exp% =="linear"
W_EMI.fx('%gas%',t)$(ord(t) ge 2) = 37 + (ord(t) - 1) * 0.5;
$endif.exp


