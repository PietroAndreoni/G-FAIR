set t_rcp /1765*2500/;
set thisttot(t_rcp,t);
thisttot(t_rcp,t) = yes$(t_rcp.val-1764 eq t.val);

parameter emissions_rcp(*,*,*);
parameter forcing_rcp(*,*,*);
parameter fossilch4_frac(*,*);
parameter natemi_hist(*,*);

$gdxin "input/RCPs_consolidated.gdx"
$load emissions_rcp=Emissions,forcing_rcp=Forcing,fossilch4_frac,natemi_hist=natural_emissions
$gdxin


CONC.FX(ghg,tfirst) = conc_preindustrial(ghg);
CUMEMI.fx(tfirst) = 0;
C_ATM.fx(tfirst) = catm_preindustrial; 
RES.fx(box,tfirst) = 0;
TATM.FX(tfirst) = 0;
TSLOW.fx(tfirst) = 0;
TFAST.fx(tfirst) = 0;
IRF.fx(tfirst) = irf_preindustrial;

W_EMI.fx(ghg,t)= sum(t_rcp,emissions_rcp(t_rcp,'%rcp%',ghg)$thisttot(t_rcp,t)) / CO2toC;
FF_CH4.fx(t) = sum(t_rcp,fossilch4_frac(t_rcp,'%rcp%')$thisttot(t_rcp,t));
natural_emissions(ghg,t) = sum(t_rcp,natemi_hist(t_rcp,ghg)$thisttot(t_rcp,t));
active(ghg) = yes;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;
CONC.fx(ghg,t)$(not active(ghg)) = conc_preindustrial(ghg);

*natural_emissions('ch4',t) = 191;
*natural_emissions('n20',t) = 8.99;

solve fair using nlp minimizing OBJ;
execute_unload "historical.gdx";

