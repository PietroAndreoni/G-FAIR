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
vsl0 = 1e-5;
eta_vsl = 1;
mort_srm = 7e3;
csrm = 1e-2;
eff_decline_srm = 0.01;

pop('1') = 7.0;
y('1') = 8e1;   # trillion USD, 2020 gross world product (~85 tn$); ypc0 ~ 11.4 k$/person

* Calendar year of each model period, same mapping as Model/initialization.gms
* (t = 1 is 2020, so t_proj = yr_2020 - tstep + t*tstep).
PARAMETER yr(t) "calendar year of model period t";
yr(t) = yr_2020 - tstep + t.val * tstep;

SCALARS  g0             "initial GDP growth rate"                          /0.025/
         gpop0          "initial population growth rate"                   /0.01/
         gdp_flat_until "last year of undiminished GDP growth"             /2030/
         gdp_zero_from  "year GDP growth has fallen to zero"               /2200/
         pop_zero_from  "year population growth has fallen to zero"        /2100/;

* GDP growth holds at g0 through gdp_flat_until, then falls linearly to zero at
* gdp_zero_from and stays there.
g(t) = g0 * min(1, max(0, (gdp_zero_from - yr(t)) / (gdp_zero_from - gdp_flat_until) ));

* Population growth falls linearly from gpop0 today to zero at pop_zero_from.
gpop(t) = gpop0 * max(0, (pop_zero_from - yr(t)) / (pop_zero_from - yr('1')) );

* Both paths now plateau on their own (y peaks well below the old 10*y('1') cap
* and pop below the old 11 bn one), so no min() clipping is applied.
loop(t,y(t+1) = y(t) * (1 + g(t)) ** tstep );
loop(t,pop(t+1) = pop(t) * (1 + gpop(t)) ** tstep );
ypc(t) = y(t) / pop(t) * 1e3;

ypc0 = y('1') / pop('1') * 1e3;

TOT_FORC.l(t) = 3;

EQUATIONS 
eq_impactcc
eq_totforcghg
eq_impactsrm
eq_effsrm
eq_qsrm
eq_mortsrm
eq_costsrm
eq_damtot;

eq_impactcc(t)..           DAMFRAC_TEMP(t) =E= a0 * power( TATM(t), b0);

eq_totforcghg(t)..         TOT_FORC(t) =E= delta + ( sum(cghg, FORCING(cghg,t) ) + forcing_exogenous(t)
                                                    + sqrt( sqr( sum(cghg, FORCING(cghg,t) ) + forcing_exogenous(t) ) + sqr(delta) ) ) / 2;

eq_effsrm(t)..             EFF_SRM(t) =E= 1 - (power( Tecs/forc2x * TOT_FORC(t), b0) * (1 + sqr( FORC_SRM(t) / TOT_FORC(t) ) - 2 * ( FORC_SRM(t) / TOT_FORC(t)) * cos( srm_angle * 3.14159 / 180 ) ) - power( Tecs/forc2x * (TOT_FORC(t) - FORC_SRM(t)), b0) ) / power( Tecs/forc2x * TOT_FORC(t), b0);

eq_impactsrm(t)..          DAMFRAC_SRM(t) =E= a0 * power( TATM_GHG(t), b0) * (1 - EFF_SRM(t));

eq_qsrm(t)..               Q_SRM(t) =E= ( background_srm(t) + SRM(t) ) / ( tgtoforc * (1 - eff_decline_srm *  Q_SRM(t) ) );

eq_mortsrm(t)..            VLL(t) =E= vsl0 * (ypc(t) / ypc0) ** eta_vsl * mort_srm * Q_SRM(t);

eq_costsrm(t)..            COST_SRM(t) =E= csrm * Q_SRM(t);

eq_damtot(t)..             DAMAGES(t) =E= ( DAMFRAC_TEMP(t) + DAMFRAC_SRM(t) ) * y(t) + VLL(t) + COST_SRM(t);
