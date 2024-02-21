** Equals old FAIR with recalibrated parameters for revised F2xco2 and Millar model.
** Deletes nonnegative reservoirs. See explanation below
$eolcom #
$onMulti

$setglobal initial_conditions 'historical_run'
$setglobal gas "co2"
$setglobal rcp "RCP45"


set t /1*1000/;
alias (t,tt);

sets     tfirst(t),tsecond(t),tlast(t);
tfirst(t) = yes$(ord(t) eq 1);
tsecond(t) = yes$(ord(t) eq 2);
tlast(t) = yes$(ord(t) eq card(t));

scalar tstep "time-step of the model" /1/;
scalar delta "" /1e-3/; 

set box "boxes for co2 concentration module"
                                /     "geological processes",
                                      "deep ocean",
                                      "biosphere",
                                      "ocean mixed layer" /;

set ghg "Greenhouse gases" /'co2', 'ch4', 'n20'/;

set cghg(ghg) "Core greenhouse gases";
set oghg(ghg) "Other well-mixed greenhouse gases";
set active(ghg) "active greenhouse gases (if not, assumed constant concentration at initial levels)";
cghg("co2") = yes;
cghg("ch4") = yes;
cghg("n20") = yes;
oghg(ghg)$(not cghg(ghg)) = yes;

SCALARS
        yr_2020     "Calendar year that corresponds to model year zero"         /2020/

        tsloweq    "Thermal equilibration parameter for box 1 (m^2 per KW)"         /0.324/
        tfasteq    "Thermal equilibration parameter for box 2 (m^2 per KW)"        /0.44/
        dslow      "Thermal response timescale for deep ocean (year)"               /236/
        dfast      "Thermal response timescale for upper ocean (year)"              /4.07/
 
        irf_preindustrial     "Pre-industrial IRF100 (year)"                                        /35/
        irC      "Increase in IRF100 with cumulative carbon uptake (years per GtC)"  /0.019/
        irT      "Increase in IRF100 with warming (years per degree K)"                /4.165/
        atmosphere_mass "Mass of atmosphere (kg)"                                     /5.1352e18/ 
        atmosphere_mm   "Molecular mass of atmosphere (kg mol-1)"                     /28.97/
        CO2toC          "CO2 to carbon conversion factor"                             /0.2727/
        Tecs            "equilbrium climate sensitivity (K)" /2.75/
        Ttcr            "Transient climate response (K)" /1.6/
        forc2x          "Forcing for 2xCO2 (Wm-2)" /3.71/
        scaling_forc2x  "Scaling factor for CO2 forcing to ensure consistency with user-specified 2xforcing" /1.0/
        emi_2020            "Yearly emissions in 2020 (GtCO2)"
        cumemi_2020         "Initial cumulative emissions in 2020 (GtCO2)"      
        catm_2020           "Initial concentration in atmosphere in 2020 (GtCO2)"       
        catm_preindustrial          "Equilibrium concentration atmosphere  (GtCO2)"            
        tslow0          "Initial temperature box 1 change in 2020 (K from 1765)"  /0.1477  /
        tfast0          "Initial temperature box 2 change in 2020 (K from 1765)"  /1.099454/
        tatm0           "Initial atmospheric temperature change in 2020"          /1.24715 /;
 
PARAMETERS         emshare(box) "Carbon emissions share into Reservoir i"  
                   taubox(box)    "Decay time constant for reservoir *  (year)"
                   taughg(ghg)    "Decay time constant for ghg *  (year)"
                   forcing_coeff(ghg) "Concentration to forcing for other green-house gases [W/]"
                   natural_emissions(ghg,t) "Emissions from natural sources for non co2 gasses"
                   target_temp(t) "Target temperature";

natural_emissions(ghg,t) = 0;
 
taubox("geological processes") = 1000000;
taubox("deep ocean") = 394.4;
taubox("biosphere") = 36.53;
taubox("ocean mixed layer") = 4.304;

taughg("ch4") = 9.3;
taughg("n20") = 121;

emshare("geological processes") = 0.2173;
emshare("deep ocean") = 0.224;
emshare("biosphere") = 0.2824;
emshare("ocean mixed layer") = 0.2763;

forcing_coeff(ghg) = 0;
target_temp(t) = 0;

** INITIAL CONDITIONS TO BE CALIBRATED TO HISTORY
** CALIBRATION;

PARAMETER ghg_mm(*) "Molecular mass of greenhouse gases (kg mol-1)";
ghg_mm('co2') = 44.01;
ghg_mm('ch4') = 16.04;
ghg_mm('n20') = 44.013;
ghg_mm('c') = 12.01;
ghg_mm('n2') = 28.013;
CO2toC = ghg_mm('c') / ghg_mm('co2');

PARAMETER emitoconc(*) "Conversion factor from emissions to concentration for greenhouse gas i (Gt to ppm/ Mt to ppb)";
emitoconc(ghg) = 1e18 / atmosphere_mass * ghg_mm(ghg)  / atmosphere_mm;
emitoconc('c') = 1e18 / atmosphere_mass * ghg_mm('c')  / atmosphere_mm;
emitoconc('n2o') = emitoconc('n2o') * ghg_mm('n2') / ghg_mm('n20'); #n20 is expressed in n2 equivalent

PARAMETER res_2020(box)  "Initial concentration in Reservoir 0 in 2020 (GtCO2)";
        
PARAMETER conc_2020(ghg)  "Initial concentration of greenhouse gas i in 2020 (ppm/ppb)";
conc_2020('co2') = 410.8;
conc_2020('ch4') = 1866.0;
conc_2020('n20') = 331.1;

PARAMETER conc_preindustrial(ghg) "Pre-industrial concentration of greenhouse gas i (ppm/ppb)";
conc_preindustrial('co2') = 278.05;
conc_preindustrial('ch4') = 722.0;
conc_preindustrial('n20') = 270.0;

catm_preindustrial = conc_preindustrial('co2') / emitoconc('co2');
catm_2020 = conc_2020('co2') / emitoconc('co2'); 
cumemi_2020 = 1717.8; #from 1750, source global carbon budget 2022
emi_2020 = 37.1; #GtCO2 
res_2020(box) = emshare(box) * (catm_2020-catm_preindustrial);
scaling_forc2x = ( -2.4e-7 * sqr( conc_preindustrial('co2') ) +  7.2e-4 * conc_preindustrial('co2') -  1.05e-4 * ( 2*conc_preindustrial('n20') ) + 5.36 ) * log( 2 ) / forc2x;

PARAMETER inertia(ghg);
inertia(ghg) = 0;

VARIABLES
*Note: Stock variables correspond to levels at the END of the period
        W_EMI(ghg,t)   "Global missions of greenhouse gas i (GtCO2/MtCH4/MtN20 per year)"
        CONC(ghg,t)    "Concentration of greenhouse gas i (ppm/ppb from 1765)"
        DCONC(ghg,t)   "Difference in concentration of non co2 greenhouse gases due to emissions at time t (ppm/ppb)"
        FORCING(ghg,t) "Increase in radiative forcing due to ghg i (watts per m2 from 1765)"
        EMO(t)         "CO2 emissions from methane oxidation (GtC per year)"
        FF_CH4(t)      "Fraction of fossil methane emissions"
        RES(box,t)     "Carbon concentration in Reservoir i (GtC from 1765)"
        TATM(t)        "Increase temperature of atmosphere (degrees L from 1765)"     
        TSLOW(t)       "Increase temperature from slow response (degrees K from 1765)"
        TFAST(t)       "Increase temperature from fast response (degrees K from 1765)"
        CUMEMI(t)      "Total co2 emitted (GtC from 1765)"
        C_SINKS(t)     "Accumulated carbon in ocean and other sinks (GtC)"
        C_ATM(t)       "Accumulated carbon in atmoshpere (GtC)"
        IRF(t)         "IRF100 at time t"
        ALPHA(t)       "Carbon decay time scaling factor"
        OBJ;

VARIABLES QSLOW, QFAST;

**** IMPORTANT PROGRAMMING NOTE. Earlier implementations has reservoirs as non-negative.
**** However, these are not physical but mathematical solutions.
**** So, they need to be unconstrained so that can have negative emissions.
POSITIVE VARIABLES   CONC, IRF, alpha, FF_CH4, CUMEMI;

EQUATIONS       
        eq_reslom           "Reservoir i law of motion"
        eq_concco2          "Atmospheric concentration equation"
        eq_catm             "Atmospheric carbon equation"
        eq_cumemi           "Total emitted carbon"
        eq_csinks           "Accumulated carbon in sinks equation"
        eq_deltaconcco2     "Delta in concentration of co2"
        eq_deltaconcghg     "Delta in concentration of other GHGs"
        eq_concghg          "Concentration equation for other GHGs"
        eq_methoxi          "Methane oxidation equation"
        eq_forcco2          "CO2 forcing equation"
        eq_forcch4          "CH4 forcing equation"
        eq_forcn20          "N2O forcing equation"
        eq_forcoghg         "Other GHG forcing equation"
        eq_tatm             "Temperature-climate equation for atmosphere"
        eq_tslow            "Temperature box 1 law of motion"
        eq_tfast            "Temperature box 2 law of motion"
        eq_irflhs           "Left-hand side of IRF100 equation"
        eq_irfrhs           "Right-hand side of IRF100 equation"
        eq_inertiaup
        eq_inertiadown
        eq_obj;

$if set sai $batinclude "SAI.gms"

** Four box model for CO2 emission-to-concentrations (FAIR formulation)
eq_reslom(box,t+1)..   RES(box,t+1) =E= RES(box,t) * exp( - tstep / ( taubox(box) * ALPHA(t) ) ) +
                                        emshare(box) * ( W_EMI('co2',t+1) * CO2toC + EMO(t+1) ) * emitoconc('c') * tstep;

eq_concco2(t)$(active('co2'))..            CONC('co2',t) =E=  conc_preindustrial('co2') + sum(box, RES(box,t) )/CO2toC;

eq_deltaconcco2(t+1)$(active('co2'))..     DCONC('co2',t+1) =E= CONC('co2',t+1) - CONC('co2',t);

eq_catm(t)..                               C_ATM(t)  =E=  catm_preindustrial + CONC('co2',t) / emitoconc('co2');
        
eq_cumemi(t+1)..                           CUMEMI(t+1) =E=  CUMEMI(t) +  ( W_EMI('co2',t) + EMO(t) )*tstep;

eq_csinks(t)..                             C_SINKS(t) =E=  CUMEMI(t) - ( C_ATM(t) -  catm_preindustrial );
    
** Single box model for non-CO2 GHGs  
eq_deltaconcghg(ghg,t)$(not sameas(ghg,'co2') and active(ghg))..   DCONC(ghg,t) =E= ( W_EMI(ghg,t) + natural_emissions(ghg,t) ) * emitoconc(ghg)  * tstep;
    
eq_concghg(ghg,t+1)$(not sameas(ghg,'co2') and active(ghg))..      CONC(ghg,t+1) =E= CONC(ghg,t) * ( exp(-tstep/taughg(ghg)) ) + ( DCONC(ghg,t+1) + DCONC(ghg,t) ) / 2;

** methanize oxidation to CO2
eq_methoxi(t)..           EMO(t) =E= 1e-3 * ghg_mm('c') / ghg_mm('ch4') * 0.61 * FF_CH4(t) * (CONC('ch4',t) - conc_preindustrial('ch4')) * (1 - exp(-1/taughg('ch4')) ) ;

** forcing for the three main greenhouse gases (CO2, CH4, N2O) 
eq_forcco2(t)..         FORCING('co2',t) =E=  ( -2.4e-7 * sqr( CONC('co2',t) - conc_preindustrial('co2') ) +
                                                7.2e-4 * ( sqrt( sqr( CONC('co2',t) - conc_preindustrial('co2') ) + sqr(delta) ) - delta ) -
                                                1.05e-4 * ( CONC('n20',t) +  conc_preindustrial('n20') ) + 5.36 ) *
                                                log( CONC('co2',t) / conc_preindustrial('co2') ) / scaling_forc2x;
 
eq_forcch4(t)..         FORCING('ch4',t) =E=  ( -6.5e-7 * (CONC('ch4',t) +  conc_preindustrial('ch4')) -
                                                4.1e-6 * (CONC('n20',t) +  conc_preindustrial('n20')) + 0.043 ) * 
                                                ( sqrt(CONC('ch4',t)) - sqrt(conc_preindustrial('ch4')) );

eq_forcn20(t)..         FORCING('n20',t) =E=  ( -4.0e-6 * (CONC('co2',t) +  conc_preindustrial('co2')) +
                                                2.1e-6 * (CONC('n20',t) +  conc_preindustrial('n20')) -
                                                2.45e-6 * (CONC('ch4',t) +  conc_preindustrial('ch4')) + 0.117 ) * 
                                                ( sqrt(CONC('n20',t)) - sqrt(conc_preindustrial('n20')) );

** forcing for other well-mixed greenhouse gases (F-gases, SOx, BC, OC, NH3, CO, NMVOC, NOx)  
eq_forcoghg(oghg,t)..     FORCING(oghg,t) =E=  (CONC(oghg,t) - conc_preindustrial(oghg)) * forcing_coeff(oghg);

** forcing to temperature 
eq_tslow(t+1)..  TSLOW(t+1) =E=  TSLOW(t) * exp(-tstep/dslow) + QSLOW * sum(ghg, FORCING(ghg,t) * ( 1 - exp(-tstep/dslow) ) );

eq_tfast(t+1)..  TFAST(t+1) =E=  TFAST(t) * exp(-tstep/dfast) + QFAST * sum(ghg, FORCING(ghg,t) * ( 1 - exp(-tstep/dfast) ) );

eq_tatm(t)..       TATM(t)  =E=  TSLOW(t) + TFAST(t);

** calculate alphas imposing IRF 
eq_irflhs(t)..    IRF(t)    =E=  sum(box, ( ALPHA(t) * emshare(box) * taubox(box) * ( 1 - exp(-100/(ALPHA(t)*taubox(box)) ) ) ) );

eq_irfrhs(t+1)..    IRF(t+1)    =E=  irf_preindustrial + irC * C_SINKS(t) * CO2toC + irT * TATM(t);

eq_inertiaup(t+1,ghg)$(not inertia(ghg) eq 0).. W_EMI(ghg,t+1) =L= (1+inertia(ghg))*W_EMI(ghg,t);

eq_inertiadown(t+1,ghg)$(not inertia(ghg) eq 0).. W_EMI(ghg,t+1) =G= (1-inertia(ghg))*W_EMI(ghg,t);

eq_obj..          OBJ =E= sqr(CONC('co2','255') - conc_2020('co2')); #sum(t,sqr(TATM(t)-target_temp(t)) );

**  Upper and lower bounds for stability
CONC.LO(cghg,t) = 1e-9;
CONC.LO(oghg,t) = 0;
TATM.LO(t)  = -10;
TATM.UP(t)  = 20;
alpha.up(t) = 100;
alpha.lo(t) = 1e-2;
IRF.up(t) = 97;
FF_CH4.up(t) = 1;
W_EMI.scale(ghg,t) = 1e3;
** Starting guess
ALPHA.l(t) = 0.35;

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
 
model fair / all /;


**** PRE MODEL 1
************ find QSLOW and QFAST given climate specifications
EQUATIONS eq_tecs, eq_ttcr;
** calculate Qi imposing ECS and TCR (more efficient as a parameter)
eq_tecs..          Tecs =E= forc2x * (QSLOW + QFAST); 

eq_ttcr..          Ttcr =E= forc2x * (QSLOW * (1 - dslow/69.7 * (1 - exp(-69.7/dslow)) ) +
                            QFAST * (1 - dfast/69.7 * (1 - exp(-69.7/dfast)) ) ) ; 

model solveqs  / eq_tecs, eq_ttcr /;
solve solveqs using cns; 

QSLOW.fx = QSLOW.l; 
QFAST.fx = QFAST.l;
*************** end QSLOW and QFAST


**** PRE MODEL 2
*************** find "natural" non co2 emissions that grant constant concentrations at steady state (pre-industrial)
EQUATIONS eq_constantconc,eq_constantemi;

VARIABLE NATEMI(ghg);

eq_constantconc..      OBJ =E= sum((ghg,t), sqr(CONC(ghg,t) - conc_preindustrial(ghg)));

eq_constantemi(ghg,t)..       NATEMI(ghg) =E= W_EMI(ghg,t);
    
model constant_concentrations_ghg /eq_deltaconcghg,eq_concghg,eq_constantconc,eq_constantemi/;

CONC.FX(ghg,tfirst) = conc_preindustrial(ghg);
active(ghg)$(not sameas(ghg,'co2')) = yes;
solve constant_concentrations_ghg using nlp minimizing obj;
natural_emissions(ghg,t) = NATEMI.l(ghg);
active(ghg) = no;
*************** end natural emissions

*** solve the model from pre-industrial era to 2020
$batinclude "run_historical.gms"

* Initial conditions
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
FF_CH4.fx(tfirst) = FF_CH4.l(tfirst);
target_temp(t) = tatm0;
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
$endif.ic

parameter save_base(ghg,t,*);
parameter save_delta(ghg,t,*);

*** solve the basic model
active(ghg) = yes;
$if set sai active('sai') = no;
W_EMI.fx(ghg,t) = 0;
solve fair using nlp minimizing OBJ;
save_base(ghg,t,'conc') = CONC.l(ghg,t);
save_base(ghg,t,'forc') = FORCING.l(ghg,t);
save_base(ghg,t,'T') = TATM.l(t);
execute_unload "simulation.gdx";

***** run some experiments
$ifthen.exp set experiment_ghg 
$if set sai active('sai') = yes;
$batinclude "experiments/GHGs.gms" "%experiment_ghg%" "%gas%"
$if set sai W_EMI.up('sai',t) = +inf;
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');

execute_unload "%experiment_ghg%_%gas%_%initial_conditions%";

$ifthen.trem set tremoval 
W_EMI.fx('%gas%','%tremoval%') = -(1e-6$(sameas('%gas%','co2')) + 1e-3$(not sameas('%gas%','co2')));
solve fair using nlp minimizing OBJ;

save_delta(ghg,t,'conc') = CONC.l(ghg,t)-save_base(ghg,t,'conc');
save_delta(ghg,t,'forc') = FORCING.l(ghg,t)-save_base(ghg,t,'forc');
save_delta(ghg,t,'T') = TATM.l(t)-save_base(ghg,t,'T');

execute_unload "%experiment_ghg%_removal%tremoval%_%gas%_%initial_conditions%";
$endif.trem

$endif.exp
