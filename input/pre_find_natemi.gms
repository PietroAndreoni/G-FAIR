**** PRE MODEL 2
*************** find "natural" non co2 emissions that grant constant concentrations at steady state (pre-industrial)
EQUATIONS eq_constantconc,eq_constantemi;

VARIABLE NATEMI(ghg);

eq_constantconc..      OBJ =E= sum((ghg,t), sqr(CONC(ghg,t) - conc_preindustrial(ghg)));

eq_constantemi(ghg,t)..       NATEMI(ghg) =E= W_EMI(ghg,t);
    
model constant_concentrations_ghg /eq_concghg,eq_constantconc,eq_constantemi/;

CONC.FX(ghg,tfirst) = conc_preindustrial(ghg);
solve constant_concentrations_ghg using nlp minimizing obj;
natural_emissions(ghg,t) = NATEMI.l(ghg);
*************** end natural emissions
