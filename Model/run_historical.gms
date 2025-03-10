parameter emissions_rcp(*,*,*);
parameter fossilch4_frac(*,*);
parameter natemi_hist(*,*);

$gdxin "input/RCPs_consolidated.gdx"
$load emissions_rcp=Emissions,fossilch4_frac,natemi_hist=natural_emissions
$gdxin

W_EMI.fx(ghg,t)= sum(t_rcp,emissions_rcp(t_rcp,'%rcp%',ghg)$thisttot(t_rcp,t))/tstep;
W_EMI.fx('co2',t) = W_EMI.l('co2',t) / CO2toC;
FF_CH4.fx(t) = sum(t_rcp,fossilch4_frac(t_rcp,'%rcp%')$thisttot(t_rcp,t))/tstep;
natural_emissions(ghg,t) = sum(t_rcp,natemi_hist(t_rcp,ghg)$thisttot(t_rcp,t))/tstep;
wemi_pre(pre,t) = sum(t_rcp,emissions_rcp(t_rcp,'%rcp%',pre)$thisttot(t_rcp,t))/tstep;
active(ghg) = yes;

*** initial conditions 
CONC.FX(ghg,tfirst) = conc_preindustrial(ghg);
CUMEMI.fx(tfirst) = 0;
C_ATM.fx(tfirst) = catm_preindustrial; 
RES.fx(box,tfirst) = 0; #emshare(box) * ( W_EMI.l('co2',tfirst) * CO2toC ) * emitoconc('c') * tstep;
TATM.FX(tfirst) = 0;
TSLOW.fx(tfirst) = 0;
TFAST.fx(tfirst) = 0;
IRF.fx(tfirst) = irf_preindustrial;

W_EMI.fx(ghg,t)$(not active(ghg)) = 0;

** fix forcing instead of emissions for non active species
FORCING.fx(ghg,t)$(not active(ghg)) = sum(t_rcp,forcing_rcp(t_rcp,'%rcp%',ghg)$thisttot(t_rcp,t))/tstep;
FORCING.fx(ghg,t)$(ord(t) ge card(t_rcp) and not active(ghg)) = forcing_rcp('2500','%rcp%',ghg);

solve fair using nlp minimizing OBJ;
execute_unload "Results/historical_%rcp%.gdx";
