** Equals old FAIR with recalibrated parameters for revised F2xco2 and Millar model.
** Deletes nonnegative reservoirs. See explanation below
$eolcom #
$onMulti

$setglobal initial_conditions 'historical_run'
$setglobal gas "co2"
$setglobal rcp "RCP45"
*$setglobal no_oforc

*** time
set t /1*1000/;
alias (t,tt);

sets     tfirst(t),tsecond(t),tlast(t);
tfirst(t) = yes$(ord(t) eq 1);
tsecond(t) = yes$(ord(t) eq 2);
tlast(t) = yes$(ord(t) eq card(t));

scalar tstep "time-step of the model" /1/;

set t_rcp /1765*2500/;
set thisttot(t_rcp,t);
thisttot(t_rcp,t)$( (1765 - tstep + t.val*tstep) ge (t_rcp.val-tstep/2) and (1765 - tstep + t.val*tstep) lt (t_rcp.val+tstep/2) ) = yes;

set t_proj /2020*2500/;
set tprojtot(t_proj,t);
tprojtot(t_proj,t)$( (2020 - tstep + t.val*tstep) ge (t_proj.val-tstep/2) and (2020  - tstep + t.val*tstep) lt (t_proj.val+tstep/2) ) = yes;


scalar delta "Delta for smooth approximation" /1e-8/; 

set box "boxes for co2 concentration module"
                                /     "geological processes",
                                      "deep ocean",
                                      "biosphere",
                                      "ocean mixed layer" /;

set ghg "Greenhouse gases" /'co2', 'ch4', 'n2o', 'h2o', 'o3trop'/;
set pre "Precursor gases (that are not ghgs)" /'co','no_x','nmvoc'/; 
set cghg(ghg) "Core greenhouse gases";
set oghg(ghg) "Other well-mixed greenhouse gases with emission and concentration represenation";
set preghg(ghg) "Other greenhouse gases without emission and concentration represenation";
set active(ghg) "active greenhouse gases (if not, assumed constant concentration at initial levels)";
cghg("co2") = yes;
cghg("ch4") = yes;
cghg("n2o") = yes;
preghg("o3trop") = yes;
preghg("h2o") = yes;
oghg(ghg)$(not cghg(ghg) and not preghg(ghg)) = yes;

SCALARS
        yr_2020     "Calendar year that corresponds to model year zero"         /2020/

        tsloweq    "Thermal equilibration parameter for box 1 (m^2 per KW)"         /0.324/
        tfasteq    "Thermal equilibration parameter for box 2 (m^2 per KW)"        /0.44/
        dslow      "Thermal response timescale for deep ocean (year)"               /236/
        dfast      "Thermal response timescale for upper ocean (year)"              /4.07/
 
        irf_preindustrial     "Pre-industrial IRF100 (%)"                         /35/
        irf_max               "Maximum IRF100 (%)"                                /97/
        irC      "Increase in IRF100 with cumulative carbon uptake (years per GtC)"  /0.019/
        irT      "Increase in IRF100 with warming (years per degree K)"                /4.165/
        atmosphere_mass "Mass of atmosphere (kg)"                                     /5.1352e18/ 
        atmosphere_mm   "Molecular mass of atmosphere (kg mol-1)"                     /28.97/
        CO2toC          "CO2 to carbon conversion factor"                             /0.2727/
        Tecs            "equilbrium climate sensitivity (K)" /3.24/
        Ttcr            "Transient climate response (K)" /1.79/
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
                   forcing_coeff(ghg) "Concentration to forcing for other green-house gases [W/m2]"
                   natural_emissions(ghg,t) "Emissions from natural sources for non co2 gasses"
                   wemi_pre(pre,t)          "Emissions of precursor gases"
                   forcing_exogenous(t) "Exogenous forcing from natural sources and exogenous oghgs [W/m2]"
                   forcing_srm(t) "Exogenous forcing from solar radiation management [W/m2]"
                   target_temp(t) "Target temperature";

PARAMETER ghg_mm(*) "Molecular mass of greenhouse gases (kg mol-1)";

# Conversion between ppb/ppt concentrations and Mt/kt emissions
# in the RCP databases ppb = Mt and ppt = kt so factor always 1e18
PARAMETER emitoconc(*) "Conversion factor from emissions to concentration for greenhouse gas i (Gt to ppm/ Mt to ppb)";

PARAMETER res_2020(box)  "Initial concentration in Reservoir 0 in 2020 (GtCO2)";
        
PARAMETER conc_2020(ghg)  "Initial concentration of greenhouse gas i in 2020 (ppm/ppb)";

PARAMETER conc_preindustrial(ghg) "Pre-industrial concentration of greenhouse gas i (ppm/ppb)";

VARIABLES
*Note: Stock variables correspond to levels at the END of the period
        W_EMI(ghg,t)   "Global missions of greenhouse gas i (GtCO2/MtCH4/MtN20 per year)"
        CONC(ghg,t)    "Concentration of greenhouse gas i (ppm/ppb from 1765)"
        FORCING(ghg,t) "Increase in radiative forcing due to ghg i (watts per m2 from 1765)"
        OXI_CH4(t)         "CO2 emissions from methane oxidation (GtC per year)"
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
        SRM(t)         "Forcing masked with solar radiation management (W/m2)"
        OBJ;

VARIABLES QSLOW, QFAST;

**** IMPORTANT PROGRAMMING NOTE. Earlier implementations has reservoirs as non-negative.
**** However, these are not physical but mathematical solutions.
**** So, they need to be unconstrained so that can have negative emissions.
POSITIVE VARIABLES   CONC, IRF, alpha, FF_CH4, CUMEMI, SRM;

EQUATIONS       
        eq_reslom           "Reservoir i law of motion"
        eq_concco2          "Atmospheric concentration equation"
        eq_catm             "Atmospheric carbon equation"
        eq_cumemi           "Total emitted carbon"
        eq_csinks           "Accumulated carbon in sinks equation"
        eq_concghg          "Concentration equation for other GHGs"
        eq_methoxi          "Methane oxidation equation"
        eq_forcco2          "CO2 forcing equation"
        eq_forcch4          "CH4 forcing equation"
        eq_forcn20          "N2O forcing equation"
        eq_forch2o
        eq_forco3trop
        eq_forcoghg         "Other GHG forcing equation"
        eq_tatm             "Temperature-climate equation for atmosphere"
        eq_tslow            "Temperature box 1 law of motion"
        eq_tfast            "Temperature box 2 law of motion"
        eq_irflhs           "Left-hand side of IRF100 equation"
        eq_irfrhs           "Right-hand side of IRF100 equation"
        eq_obj;

***** parameters initialization
$batinclude "Model/parameters.gms"

** Four box model for CO2 emission-to-concentrations (FAIR formulation)
eq_reslom(box,t+1)$(active('co2'))..   RES(box,t+1) =E= RES(box,t) * exp( - tstep / ( taubox(box) * ALPHA(t) ) ) +
                                        emshare(box) * ( W_EMI('co2',t+1) + OXI_CH4(t+1) ) * emitoconc('co2') * tstep;

eq_concco2(t)$(active('co2'))..            CONC('co2',t) =E=  conc_preindustrial('co2') + sum(box, RES(box,t) );

eq_catm(t)$(active('co2'))..               C_ATM(t)  =E=  CONC('co2',t) / emitoconc('co2');
        
eq_cumemi(t+1)$(active('co2'))..           CUMEMI(t+1) =E=  CUMEMI(t) +  ( W_EMI('co2',t+1) + OXI_CH4(t+1) )*tstep;

eq_csinks(t)$(active('co2'))..             C_SINKS(t) =E=  CUMEMI(t) - ( C_ATM(t) -  catm_preindustrial );
    
** Single box model for non-CO2 GHGs  
eq_concghg(ghg,t+1)$(not sameas(ghg,'co2') and not preghg(ghg) and active(ghg))..      
                        CONC(ghg,t+1) =E= CONC(ghg,t) * exp(-tstep/taughg(ghg)) + 
                        ( (  W_EMI(ghg,t+1) +   W_EMI(ghg,t) ) / 2 + natural_emissions(ghg,t+1) ) * emitoconc(ghg)  * tstep;

** methanize oxidation to CO2
eq_methoxi(t)..         OXI_CH4(t) =E= 1e-3 * ghg_mm('co2') / ghg_mm('ch4') * 0.61 * FF_CH4(t) * (CONC('ch4',t) - conc_preindustrial('ch4')) * (1 - exp(-tstep/taughg('ch4')) ) ;

** forcing for the three main greenhouse gases (CO2, CH4, N2O) 
eq_forcco2(t)..         FORCING('co2',t) =E=  ( -2.4e-7 * sqr( CONC('co2',t) - conc_preindustrial('co2') ) +
                                                7.2e-4 * ( sqrt( sqr( CONC('co2',t) - conc_preindustrial('co2') ) + sqr(delta) ) - delta ) -
                                                1.05e-4 * ( CONC('n2o',t) +  conc_preindustrial('n2o') ) + 5.36 ) *
                                                log( CONC('co2',t) / conc_preindustrial('co2') ) / scaling_forc2x;
 
eq_forcch4(t)..         FORCING('ch4',t) =E=  ( -6.5e-7 * (CONC('ch4',t) +  conc_preindustrial('ch4')) -
                                                4.1e-6 * (CONC('n2o',t) +  conc_preindustrial('n2o')) + 0.043 ) * 
                                                ( sqrt(CONC('ch4',t)) - sqrt(conc_preindustrial('ch4')) );

eq_forcn20(t)..         FORCING('n2o',t) =E=  ( -4.0e-6 * (CONC('co2',t) +  conc_preindustrial('co2')) +
                                                2.1e-6 * (CONC('n2o',t) +  conc_preindustrial('n2o')) -
                                                2.45e-6 * (CONC('ch4',t) +  conc_preindustrial('ch4')) + 0.117 ) * 
                                                ( sqrt(CONC('n2o',t)) - sqrt(conc_preindustrial('n2o')) );

eq_forch2o(t)..         FORCING('h2o',t) =E= 0.12 * FORCING('ch4',t); 

eq_forco3trop(t)..      FORCING('o3trop',t) =E= 1.74e-4 * (CONC('ch4',t) - conc_preindustrial('ch4')) +
                                            9.08e-4 * (wemi_pre('no_x',t)-2) +
                                            8.51e-5 * (wemi_pre('co',t)-170) +
                                            2.25e-4 * (wemi_pre('nmvoc',t)-5) +
                                            ( 0.032 * (exp(-1.35*(TATM(t)) ) - 1) - sqrt( sqr(0.032 * (exp(-1.35*(TATM(t))) - 1) ) + sqr(1e-8)) ) / 2   ;

** forcing for other well-mixed greenhouse gases (F-gases, SOx, BC, OC, NH3, CO, NMVOC, NOx)  
eq_forcoghg(oghg,t)$(active(oghg))..     FORCING(oghg,t) =E=  (CONC(oghg,t) - conc_preindustrial(oghg)) * forcing_coeff(oghg);

** forcing to temperature 
eq_tslow(t+1)..  TSLOW(t+1) =E=  TSLOW(t) * exp(-tstep/dslow) + QSLOW * ( sum(ghg, FORCING(ghg,t) ) + forcing_exogenous(t) - forcing_srm(t) - SRM(t) ) * ( 1 - exp(-tstep/dslow) );

eq_tfast(t+1)..  TFAST(t+1) =E=  TFAST(t) * exp(-tstep/dfast) + QFAST * ( sum(ghg, FORCING(ghg,t) ) + forcing_exogenous(t) - forcing_srm(t) - SRM(t) ) * ( 1 - exp(-tstep/dfast) );

eq_tatm(t)..       TATM(t)  =E=  TSLOW(t) + TFAST(t);

** calculate alphas imposing IRF 
eq_irflhs(t)$(active('co2'))..    IRF(t)    =E= ALPHA(t) * sum(box, emshare(box) * taubox(box) * ( 1 - exp(-100/(ALPHA(t)*taubox(box)) ) ) );

*** IRF max is 97. Smooth GAMS approximation: [f(x) + g(y) - sqrt(sqr(f(x)-g(y)) + sqr(delta))] /2
eq_irfrhs(t)$(active('co2'))..    IRF(t)    =E= ( ( irf_max + ( irf_preindustrial + irC * C_SINKS(t) * CO2toC + irT * TATM(t) ) ) - 
                                                    sqrt( sqr(irf_max - (irf_preindustrial + irC * C_SINKS(t) * CO2toC + irT * TATM(t) ) ) + sqr(delta) ) ) / 2;

eq_obj..          OBJ =E= sum(t,sqr((TATM(t)-target_temp(t)) ) );

**  Upper and lower bounds for stability
CONC.LO(cghg,t) = 1e-9;
CONC.LO(oghg,t) = 0;
TATM.LO(t)  = -10;
TATM.UP(t)  = 20;
ALPHA.lo(t) = 1e-2;
ALPHA.up(t) = 1e3;
IRF.up(t) = 100;
FF_CH4.up(t) = 1;
** Starting guess
ALPHA.l(t) = 0.35;
** Deactivate SRM by default
SRM.fx(t) = 0;

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
 
model fair / all /;
fair.OptFile = 1;

** find QSLOW and QFAST given TCR, ECS, and forc2x parameters 
$batinclude "Input/pre_find_Qs.gms"

** find natural emissions that grant pre-industrial equilbrium of concentrations for ch4 and n2o
$if set neutral_natemi $batinclude "Input/pre_find_natemi.gms"

*** include forcing from natural sources and exogenous 
$batinclude "Model/exogenous_forcing.gms"

*** solve the model from pre-industrial era 
$batinclude "Model/run_historical.gms"

*** initialize the model and emission scenarios
$batinclude "Model/initialization.gms"

$if set no_oforc forcing_exogenous(t) = 0;

parameter save_base(ghg,t,*);
parameter save_delta(ghg,t,*);

*** solve the basic model
active(ghg) = yes;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
execute_unload "Results/%rcp%_EXPsimulation_IC%initial_conditions%.gdx";

***** run some experiments
$if set experiment $batinclude "experiments/%experiment%.gms" "%gas%"
