** Equals old FAIR with recalibrated parameters for revised F2xco2 and Millar model.
** Deletes nonnegative reservoirs. See explanation below
$eolcom #

$setglobal initial_conditions 'preindustrial'
$setglobal gas "co2"
$setglobal experiment "emissionpulse"

set t /1*400/;
alias (t,tt);

sets     tfirst(t), tlast(t);
tfirst(t) = yes$(ord(t) eq 1);
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
        yr0     "Calendar year that corresponds to model year zero"         /2020/

        tsloweq    "Thermal equilibration parameter for box 1 (m^2 per KW)"         /0.324/
        tfasteq    "Thermal equilibration parameter for box 2 (m^2 per KW)"        /0.44/
        dslow      "Thermal response timescale for deep ocean (year)"               /236/
        dfast      "Thermal response timescale for upper ocean (year)"              /4.07/
 
        irf0     "Pre-industrial IRF100 (year)"                                        /32.4/
        irC      "Increase in IRF100 with cumulative carbon uptake (years per GtC)"  /0.019/
        irT      "Increase in IRF100 with warming (years per degree K)"                /4.165/
        atmosphere_mass "Mass of atmosphere (kg)"                                     /5.1352e18/ 
        atmosphere_mm   "Molecular mass of atmosphere (kg mol-1)"                     /28.97/
        CO2toC          "CO2 to carbon conversion factor"                             /0.2727/
        Tecs "equilbrium climate sensitivity (K)" /2.75/
        Ttcr "Transient climate response (K)" /1.6/
        forc2x "Forcing for 2xCO2 (Wm-2)" /3.71/
        scaling_forc2x "Scaling factor for CO2 forcing to ensure consistency with user-specified 2xforcing" /1.0/
        cumemi0   "Initial cumulative emissions in 2020 (GtCO2)"       /2322.833/
        catm0     "Initial concentration in atmosphere in 2020 (GtCO2)"       /3250.547/
        catmeq    "Equilibrium concentration atmosphere  (GtCO2)"            /2156/
        tslow0    "Initial temperature box 1 change in 2020 (C from 1765)"  /0.1477  /
        tfast0    "Initial temperature box 2 change in 2020 (C from 1765)"  /1.099454/
        tatm0     "Initial atmospheric temperature change in 2020"          /1.24715 /;
 
PARAMETERS         emshare(box) "Carbon emissions share into Reservoir i"  
                   taubox(box)    "Decay time constant for reservoir *  (year)"
                   taughg(ghg)    "Decay time constant for ghg *  (year)"
                   forcing_coeff(ghg) "Concentration to forcing for other green-house gases"
                   natural_emissions(ghg,t) "Emissions from natural sources for non co2 gasses";

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
** INITIAL CONDITIONS TO BE CALIBRATED TO HISTORY
** CALIBRATION;

PARAMETER ghg_mm(ghg) "Molecular mass of greenhouse gases (kg mol-1)";
ghg_mm('co2') = 44.01;
ghg_mm('ch4') = 16.04;
ghg_mm('n20') = 44.013;

PARAMETER res0(box)  "Initial concentration in Reservoir 0 in 2020 (GtCO2)";
*res0("geological processes") = 550.341;
*res0("deep ocean") = 376.559;
*res0("biosphere") = 144.958;
*res0("ocean mixed layer") = 22.6838;

        
PARAMETER conc0(ghg)  "Initial concentration of greenhouse gas i in 2020 (ppm/ppb)";
conc0('co2') = 410.8;
conc0('ch4') = 1866.0;
conc0('n20') = 331.1;

PARAMETER preindustrial_conc(ghg) "Pre-industrial concentration of greenhouse gas i (ppm/ppb)";
preindustrial_conc('co2') = 278.05;
preindustrial_conc('ch4') = 722.0;
preindustrial_conc('n20') = 270.0;

catmeq = preindustrial_conc('co2') * atmosphere_mass / 1e18 / ghg_mm('co2') * atmosphere_mm;
catm0 = conc0('co2') * atmosphere_mass / 1e18 / ghg_mm('co2') * atmosphere_mm;
res0(box) = emshare(box) * (catm0-catmeq);
cumemi0 = 1700; #from 1750, source global carbon budget 2022
scaling_forc2x = ( -2.4e-7 * sqr( preindustrial_conc('co2') ) +  7.2e-4 * preindustrial_conc('co2') -  1.05e-4 * ( 2*preindustrial_conc('n20') ) + 5.36 ) * log( 2 ) / forc2x;

VARIABLES
*Note: Stock variables correspond to levels at the END of the period
        W_EMI(ghg,t)   "Global missions of greenhouse gas i (GtCO2/MtCH4/MtN20 per year)"
        CONC(ghg,t)    "Concentration of greenhouse gas i (ppm/ppb from 1765)"
        DCONC(ghg,t)   "Difference in concentration of non co2 greenhouse gases due to emissions at time t (ppm/ppb)"
        FORCING(ghg,t) "Increase in radiative forcing due to ghg i (watts per m2 from 1765)"
        EMO(t)         "CO2 emissions from methane oxidation (GtCO2 per year)"
        FF_CH4(t)      "Fraction of fossil methane emissions"
        RES(box,t)     "Carbon concentration in Reservoir i (GtCO2 from 1765)"
        TATM(t)        "Increase temperature of atmosphere (degrees C from 1765)"     
        TSLOW(t)       "Increase temperature from slow response (degrees C from 1765)"
        TFAST(t)       "Increase temperature from fast response (degrees C from 1765)"
        CUMEMI(t)      "Total co2 emitted (GtCO2 from 1765)"
        C_SINKS(t)        "Accumulated carbon in ocean and other sinks (GtCO2)"
        C_ATM(t)        "Accumulated carbon in atmoshpere (GtCO2)"
        IRF(t)         "IRF100 at time t"
        ALPHA(t)       "Carbon decay time scaling factor"
        OBJ;

VARIABLES QSLOW, QFAST;

**** IMPORTANT PROGRAMMING NOTE. Earlier implementations has reservoirs as non-negative.
**** However, these are not physical but mathematical solutions.
**** So, they need to be unconstrained so that can have negative emissions.
POSITIVE VARIABLES   CONC, RES, IRF, alpha, FF_CH4;

EQUATIONS       
        eq_reslom           "Reservoir i law of motion"
        eq_concco2          "Atmospheric concentration equation"
        eq_catm             "Atmospheric carbon equation"
        eq_cumemi           "Total emitted carbon"
        eq_csinks           "Accumulated carbon in sinks equation"
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
        eq_obj;

** Four box model for CO2 emission-to-concentrations (nordhaus formulation)
eq_reslom(box,t+1)..   RES(box,t+1) =E=  ( emshare(box) * taubox(box) * ALPHA(t) * ( W_EMI('co2',t) + EMO(t) ) ) * 
                                        ( 1 - exp( - tstep / ( taubox(box) * ALPHA(t) ) ) ) + 
                                        RES(box,t) * exp( - tstep / ( taubox(box) * ALPHA(t) ) );


** Four box model for CO2 emission-to-concentrations (own formulation)
*eq_reslom(box,t+1)..   RES(box,t+1) =E=   exp( - tstep / ( taubox(box) * ALPHA(t+1) ) ) *
*                                        ( RES(box,'1') / exp( - tstep / ( taubox(box) * ALPHA('1') ) ) +
*                                        sum(tt$(tt.val le t.val+1), emshare(box) * ( W_EMI('co2',tt) + EMO(tt) ) *
*                                        exp( - tstep / ( taubox(box) * ALPHA(tt) ) ) ) );

eq_catm(t)..                 C_ATM(t)  =E=  catmeq + sum(box, RES(box,t) * tstep );
    
eq_cumemi(t+1)..           CUMEMI(t+1) =E=  CUMEMI(t) +  ( W_EMI('co2',t) + EMO(t) )*tstep;

eq_csinks(t)..             C_SINKS(t) =E=  CUMEMI(t) - ( C_ATM(t) -  catmeq );
    
eq_concco2(t)$(active('co2'))..            CONC('co2',t) =E=  C_ATM(t) * 1e18 / atmosphere_mass * ghg_mm('co2')  / atmosphere_mm;

** Single box model for non-CO2 GHGs  
eq_deltaconcghg(ghg,t)$(not sameas(ghg,'co2') and active(ghg))..   DCONC(ghg,t) =E= ( W_EMI(ghg,t) + natural_emissions(ghg,t) ) * 1e18 / atmosphere_mass * ghg_mm(ghg)  / atmosphere_mm  * tstep;
    
eq_concghg(ghg,t+1)$(not sameas(ghg,'co2') and active(ghg))..      CONC(ghg,t+1) =E= CONC(ghg,t) * ( exp(-tstep/taughg(ghg)) ) + ( DCONC(ghg,t) + DCONC(ghg,t+1) )/2;

** methanize oxidation to CO2
eq_methoxi(t)..           EMO(t) =E= 0.61 * FF_CH4(t) * (CONC('ch4',t) - preindustrial_conc('ch4')) * (1 - exp(-1/taughg('ch4')) ) ;

** forcing for the three main greenhouse gases (CO2, CH4, N2O) 
eq_forcco2(t)..         FORCING('co2',t) =E=  ( -2.4e-7 * sqr( CONC('co2',t) - preindustrial_conc('co2') ) +
                                                7.2e-4 * ( sqrt( sqr( CONC('co2',t) - preindustrial_conc('co2') ) + sqr(delta) ) - delta ) -
                                                1.05e-4 * ( CONC('n20',t) +  preindustrial_conc('n20') ) + 5.36 ) *
                                                log( CONC('co2',t) / preindustrial_conc('co2') ) / scaling_forc2x;
 
eq_forcch4(t)..         FORCING('ch4',t) =E=  ( -6.5e-7 * (CONC('ch4',t) +  preindustrial_conc('ch4')) -
                                                4.1e-6 * (CONC('n20',t) +  preindustrial_conc('n20')) + 0.043 ) * 
                                                ( sqrt(CONC('ch4',t)) - sqrt(preindustrial_conc('ch4')) );

eq_forcn20(t)..         FORCING('n20',t) =E=  ( -4.0e-6 * (CONC('co2',t) +  preindustrial_conc('co2')) +
                                                2.1e-6 * (CONC('n20',t) +  preindustrial_conc('n20')) -
                                                2.45e-6 * (CONC('ch4',t) +  preindustrial_conc('ch4')) + 0.117 ) * 
                                                ( sqrt(CONC('n20',t)) - sqrt(preindustrial_conc('n20')) );

** forcing for other well-mixed greenhouse gases (F-gases, SOx, BC, OC, NH3, CO, NMVOC, NOx)  
eq_forcoghg(oghg,t)..     FORCING(oghg,t) =E=  (CONC(oghg,t) - preindustrial_conc(oghg)) * forcing_coeff(oghg);

** forcing to temperature 
eq_tslow(t+1)..  TSLOW(t+1) =E=  TSLOW(t) * exp(-tstep/dslow) + QSLOW * sum(ghg, FORCING(ghg,t) * ( 1 - exp(-tstep/dslow) ) );

eq_tfast(t+1)..  TFAST(t+1) =E=  TFAST(t) * exp(-tstep/dfast) + QFAST * sum(ghg, FORCING(ghg,t) * ( 1 - exp(-tstep/dfast) ) );

eq_tatm(t)..       TATM(t)  =E=  TSLOW(t) + TFAST(t);

** calculate alphas imposing IRF 
eq_irflhs(t)..    IRF(t)    =E=  sum(box, ( ALPHA(t) * emshare(box) * taubox(box) * ( 1 - exp(-100/(ALPHA(t)*taubox(box)) ) ) ) );

eq_irfrhs(t+1)..  IRF(t+1)    =E=  irf0 + irC * C_SINKS(t) * CO2toC + irT * TATM(t);

eq_obj..          OBJ =E= sum(t,sqr(TATM(t)-tatm0) );

**  Upper and lower bounds for stability
CONC.LO(ghg,t) = 1e-3;
TATM.UP(t)  = 20;
alpha.up(t) = 100;
alpha.lo(t) = 1e-2;
FF_CH4.up(t) = 1;

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

eq_constantconc..      OBJ =E= sum((ghg,t), sqr(CONC(ghg,t) - preindustrial_conc(ghg)));

eq_constantemi(ghg,t)..       NATEMI(ghg) =E= W_EMI(ghg,t);
    
model constant_concentrations_ghg /eq_deltaconcghg,eq_concghg,eq_constantconc,eq_constantemi/;

CONC.FX(ghg,tfirst) = preindustrial_conc(ghg);
active(ghg)$(not sameas(ghg,'co2')) = yes;
solve constant_concentrations_ghg using nlp minimizing obj;
natural_emissions(ghg,t) = NATEMI.l(ghg);
active(ghg) = no;
*************** end natural emissions

$batinclude "run_historical.gms"

* Initial conditions
$ifthen.ic %initial_conditions%=="2020"
CONC.FX(ghg,tfirst) = conc0(ghg);
CUMEMI.fx(tfirst) = cumemi0;
C_ATM.fx(tfirst) = catm0; 
RES.fx(box,tfirst) = res0(box);
TATM.FX(tfirst) = tatm0;
TSLOW.fx(tfirst) = tslow0;
TFAST.fx(tfirst) = tfast0;
IRF.fx(tfirst) = irf0 + irC * (cumemi0 - (catm0-catmeq) ) * CO2toC + irT * tatm0;
$elseif.ic %initial_conditions%=="preindustrial"
CONC.FX(ghg,tfirst) = preindustrial_conc(ghg);
CUMEMI.fx(tfirst) = 0;
C_ATM.fx(tfirst) = catmeq; 
RES.fx(box,tfirst) = 0;
TATM.FX(tfirst) = 0;
TSLOW.fx(tfirst) = 0;
TFAST.fx(tfirst) = 0;
IRF.fx(tfirst) = irf0;
$endif.ic


***** run some experiments
$if set experiment $batinclude "experiments/%experiment%_%gas%.gms"

