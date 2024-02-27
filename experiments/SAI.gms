$set exp %1

***** emission pulse, co2
active('sai') = yes;
active('co2') = yes;
active('ch4') = yes;
active('n2o') = yes;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
$ifthen.ic %initial_conditions%=="2020"
CONC.fx(ghg,t)$(not active(ghg)) = conc0(ghg);
$elseif.ic %initial_conditions%=="preindustrial"
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);
$elseif.ic %initial_conditions%=="historical_run"
CONC.fx(ghg,t)$(not active(ghg)) = CONC.l(ghg,'255');
$endif.ic

FF_CH4.fx(t) = 0;

W_EMI.up(ghg,t)$(active(ghg) and not sameas(ghg,'sai'))  = +inf;
W_EMI.lo(ghg,t)$(active(ghg) and not sameas(ghg,'sai'))  = 0;
W_EMI.lo('co2',t)$(active("co2"))  = -inf; 
W_EMI.fx('sai',t) = 0; 

$ifthen.exp %exp% =="pulse"
W_EMI.fx('sai',tsecond) = 12;
$elseif.exp %exp% =="constant"
W_EMI.fx('sai',t) = 12;
$elseif.exp %exp% =="stepwise"
W_EMI.fx('sai',t)$(ord(t) le 100) = 12;
$elseif.exp %exp% =="linear"
W_EMI.fx('sai',t)$(ord(t) le 100) = 12 + (ord(t) - 1) * 0.12;
$endif.exp


