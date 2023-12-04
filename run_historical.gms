set t_hist /1750*2020/;
set thisttot(t_hist,t);
thisttot(t_hist,t) = yes$(t_hist.val-1749 eq t.val);

parameter q_emi_valid_primap(*,*,*);
parameter q_emi_valid_oscar(*,*,*);

$gdxin "input/data_historical_values.gdx"
$load q_emi_valid_primap,q_emi_valid_oscar
$gdxin

parameter emi_hist(ghg,t_hist);
parameter fossilch4_frac(t_hist);

emi_hist('co2',t_hist) = 44 / 12 *(q_emi_valid_primap('co2ffi',t_hist,'world')) ;
emi_hist('ch4',t_hist) = 1e3 / 25 * 44 / 12 *(q_emi_valid_primap('ch4_ffi',t_hist,'world')+q_emi_valid_primap('ch4_wst',t_hist,'world')+q_emi_valid_primap('ch4_agr',t_hist,'world')) ;
emi_hist('n20',t_hist) = 1e3 / 298 * 44 / 12 *(q_emi_valid_primap('n2o_ffi',t_hist,'world')+q_emi_valid_primap('n2o_wst',t_hist,'world')+q_emi_valid_primap('n2o_agr',t_hist,'world')) ;
fossilch4_frac(t_hist) = 1e3 / 25 * 44 / 12 * q_emi_valid_primap('ch4_ffi',t_hist,'world')/emi_hist('ch4',t_hist);

CONC.FX(ghg,tfirst) = preindustrial_conc(ghg);
CUMEMI.fx(tfirst) = 0;
C_ATM.fx(tfirst) = catmeq; 
RES.fx(box,tfirst) = 0;
TATM.FX(tfirst) = 0;
TSLOW.fx(tfirst) = 0;
TFAST.fx(tfirst) = 0;
IRF.fx(tfirst) = irf0;

W_EMI.fx(ghg,t)= sum(t_hist,emi_hist(ghg,t_hist)$thisttot(t_hist,t)) ;
FF_CH4.fx(t) = sum(t_hist,fossilch4_frac(t_hist)$thisttot(t_hist,t));
FF_CH4.fx(t) = 0;

active(ghg) = yes;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
CONC.fx(ghg,t)$(not active(ghg)) = preindustrial_conc(ghg);

*natural_emissions('ch4',t) = 191;
*natural_emissions('n20',t) = 8.99;

solve fair using nlp minimizing OBJ;
execute_unload "historical.gdx";

