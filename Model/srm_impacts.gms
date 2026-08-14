# this modules adds economic impacts to the model


Variable 
    DAMFRAC_TEMP(t)      "fraction of GDP lost to climate change damages"
    DAMFRAC_SRM(t)       "fraction of GDP lost to SRM damages"
    EFF_SRM(t)           "efficiency of SRM in offsetting climate change"
    Q_SRM(t)             "quantity of SRM deployed"
    VLL(t)               "value of lost lives due to SRM"
    COST_SRM(t)          "cost of SRM deployment"
    DAMAGES(t)           "total damages including climate change and SRM";

POSITIVE VARIABLE TOT_FORC, Q_SRM;
Q_SRM.up(t) = 99;

PARAMETERS 
a0 "damage function parameter for climate change"
b0 "damage function exponent for climate change"
srm_angle "angle of SRM forcing relative to GHG forcing"
tgtoforc "SRM forcing per unit of SRM"
vsl0 "value of statistical life at base year"
ypc0 "per capita income at base year"
eta_vsl "elasticity of VSL with respect to income"
mort_srm "mortality rate associated with SRM deployment"
csrm "cost per unit of SRM deployment"
eff_decline_srm "efficiency decline rate of SRM as a function of deployment";

*economic variables in trillion USD
*population in billions
*GDP per capita in USD
PARAMETER y(t) "GDP at time t"
          pop(t) "population at time t"
          ypc(t) "per capita income at time t"
          g(t) "growth rate of GDP at time t"
          gpop(t) "growth rate of population at time t";

a0 = 0.00575;
b0 = 2.0;
srm_angle = 10;
tgtoforc = 0.25 * Tecs;
vsl0 = 1e-6;
eta_vsl = 1;
mort_srm = 7e3;
csrm = 1e-3;
eff_decline_srm = 0.01;

pop('1') = 7.0; 
y('1') = 8e2;

g(t) = 0.02;
gpop(t) = 0.01;

loop(t,y(t+1) = y(t) * (1 + g(t)) ** tstep );
y(t) = min(y(t),10*y('1'));
loop(t,pop(t+1) = pop(t) * (1 + gpop(t)) ** tstep );
pop(t) = min(pop(t),11);
ypc(t) = y(t) / pop(t) * 1e3;

ypc0 = y('1') / pop('1') * 1e3;

TOT_FORC.l(t) = 3;

EQUATIONS 
eq_impactcc
eq_totforcghg
eq_impactsrm
eq_srmeff
eq_qsrm
eq_mortsrm
eq_costsrm
eq_damtot;

eq_impactcc(t)..           DAMFRAC_TEMP(t) =E= a0 * power(TATM(t), b0) - a0 * power(TATM('1'), b0);

eq_impactsrm(t)..          DAMFRAC_SRM(t) =E= a0 * power(TATM_GHG(t), b0) * EFF_SRM(t) - a0 * power(TATM(t), b0);

eq_totforcghg(t)..         TOT_FORC(t) =E= delta + ( sum(ghg, FORCING(ghg,t) ) + forcing_exogenous(t)
                                                    + sqrt( sqr( sum(ghg, FORCING(ghg,t) ) + forcing_exogenous(t) ) + sqr(delta) ) ) / 2;

eq_srmeff(t)..             EFF_SRM(t) =E= 1 + sqr( FORC_SRM(t) / TOT_FORC(t) ) - 2 * ( FORC_SRM(t) / TOT_FORC(t)) * cos( srm_angle * 180/3.14159 );

eq_qsrm(t)..               Q_SRM(t) =E= ( background_srm(t) + SRM(t) ) / ( tgtoforc * (1 - eff_decline_srm *  Q_SRM(t) ) );

eq_mortsrm(t)..            VLL(t) =E= vsl0 * (ypc(t) / ypc0) ** eta_vsl * mort_srm * Q_SRM(t);

eq_costsrm(t)..            COST_SRM(t) =E= csrm * Q_SRM(t);

eq_damtot(t)..             DAMAGES(t) =E= ( DAMFRAC_TEMP(t) + DAMFRAC_SRM(t) ) * y(t) + VLL(t) + COST_SRM(t);
