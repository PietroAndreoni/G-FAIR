
*Initial conditions
$ifthen.ic %initial_conditions%=="2020"
CONC.FX(ghg,tfirst) = conc_2020(ghg);
CUMEMI.fx(tfirst) = cumemi_2020;
C_ATM.fx(tfirst) = catm_2020; 
RES.fx(box,tfirst) = res_2020(box);
TATM.FX(tfirst) = tatm0;
TSLOW.fx(tfirst) = tslow0;
TFAST.fx(tfirst) = tfast0;
IRF.fx(tfirst) = irf_preindustrial + irC * (cumemi_2020 - (catm_2020-catm_preindustrial) ) * CO2toC + irT * tatm0;
FF_CH4.fx(t) = 0;
FF_CH4.fx(tfirst) = FF_CH4.l('255');
target_temp(t) = tatm0;

W_EMI.fx(ghg,t)= sum(t_proj,emissions_rcp(t_proj,'%rcp%',ghg)$tprojtot(t_proj,t))/tstep;
W_EMI.fx(ghg,t)$(ord(t) ge card(t_proj)) = emissions_rcp('2500','%rcp%',ghg);
W_EMI.fx('co2',t) = W_EMI.l('co2',t) / CO2toC;
forcing_exogenous(t)= sum((t_proj,sources),forcing_rcp(t_proj,'%rcp%',sources)$tprojtot(t_proj,t))/tstep;
forcing_exogenous(t)$(ord(t) ge card(t_proj)) = sum(sources,forcing_rcp('2500','%rcp%',sources));

** fix forcing instead of emissions for non active species
FORCING.fx(ghg,t)$(not active(ghg)) = sum(t_proj,forcing_rcp(t_proj,'%rcp%',ghg)$tprojtot(t_proj,t))/tstep;
FORCING.fx(ghg,t)$(ord(t) ge card(t_proj) and not active(ghg)) = forcing_rcp('2500','%rcp%',ghg);

$elseif.ic %initial_conditions%=="historical_run"
CONC.FX(ghg,tfirst) =  CONC.l(ghg,'255');
CUMEMI.fx(tfirst) = CUMEMI.l('255');
C_ATM.fx(tfirst) = C_ATM.l('255'); 
RES.fx(box,tfirst) = RES.l(box,'255');
TATM.FX(tfirst) = TATM.l('255');
TSLOW.fx(tfirst) = TSLOW.l('255');
TFAST.fx(tfirst) = TFAST.l('255');
IRF.fx(tfirst) = irf_preindustrial + irC * (CUMEMI.l('255') - (C_ATM.l('255')-catm_preindustrial) ) * CO2toC + irT * TATM.l('255');
target_temp(t) = TATM.l('255');

W_EMI.fx(ghg,t)= sum(t_proj,emissions_rcp(t_proj,'%rcp%',ghg)$tprojtot(t_proj,t))/tstep;
W_EMI.fx(ghg,t)$(ord(t) ge card(t_proj)) = emissions_rcp('2500','%rcp%',ghg);
W_EMI.fx('co2',t) = W_EMI.l('co2',t) / CO2toC;
FF_CH4.fx(t) = sum(t_proj,fossilch4_frac(t_proj,'%rcp%')$tprojtot(t_proj,t))/tstep;
FF_CH4.fx(t)$(ord(t) ge card(t_proj)) = fossilch4_frac('2500','%rcp%');

wemi_pre(pre,t) = sum(t_proj,emissions_rcp(t_proj,'%rcp%',pre)$tprojtot(t_proj,t))/tstep;

forcing_exogenous(t)= sum((t_proj,sources),forcing_rcp(t_proj,'%rcp%',sources)$tprojtot(t_proj,t))/tstep;
forcing_exogenous(t)$(ord(t) ge card(t_proj)) = sum(sources,forcing_rcp('2500','%rcp%',sources));

** fix forcing instead of emissions for non active species
FORCING.fx(ghg,t)$(not active(ghg)) = sum(t_proj,forcing_rcp(t_proj,'%rcp%',ghg)$tprojtot(t_proj,t))/tstep;
FORCING.fx(ghg,t)$(ord(t) ge card(t_proj) and not active(ghg)) = forcing_rcp('2500','%rcp%',ghg);

$elseif.ic %initial_conditions%=="preindustrial"
CONC.FX(ghg,tfirst) = conc_preindustrial(ghg);
CUMEMI.fx(tfirst) = 0;
C_ATM.fx(tfirst) = catm_preindustrial; 
RES.fx(box,tfirst) = 0;
TATM.FX(tfirst) = 0;
TSLOW.fx(tfirst) = 0;
TFAST.fx(tfirst) = 0;
IRF.fx(tfirst) = irf_preindustrial;
FF_CH4.fx(t) = 0;
target_temp(t) = 0;
W_EMI.fx(ghg,t) = 0;
forcing_exogenous(t) = 0;
$endif.ic