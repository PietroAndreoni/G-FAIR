natural_emissions(ghg,t) = 0;
target_temp(t) = 0;
forcing_coeff(cghg) = 0;

taubox("geological processes") = 1000000;
taubox("deep ocean") = 394.4;
taubox("biosphere") = 36.53;
taubox("ocean mixed layer") = 4.304;

taughg("ch4") = 9.3;
taughg("n2o") = 121;

emshare("geological processes") = 0.2173;
emshare("deep ocean") = 0.224;
emshare("biosphere") = 0.2824;
emshare("ocean mixed layer") = 0.2763;

ghg_mm('co2') = 44.01;
ghg_mm('ch4') = 16.04;
ghg_mm('n2o') = 44.013;
ghg_mm('c') = 12.01;
ghg_mm('n2') = 28.013;
CO2toC = ghg_mm('c') / ghg_mm('co2');

emitoconc(ghg) = 1e18 / atmosphere_mass * atmosphere_mm / ghg_mm(ghg) ;
emitoconc('c') = 1e18 / atmosphere_mass * atmosphere_mm / ghg_mm('c');
emitoconc('n2o') = emitoconc('n2o') * ghg_mm('n2o') / ghg_mm('n2') ; #n2o is expressed in n2 equivalent

conc_2020('co2') = 410.8;
conc_2020('ch4') = 1866.0;
conc_2020('n2o') = 331.1;

conc_preindustrial('co2') = 278.05;
conc_preindustrial('ch4') = 722.0;
conc_preindustrial('n2o') = 255.0;

catm_preindustrial = conc_preindustrial('co2') / emitoconc('co2');
catm_2020 = conc_2020('co2') / emitoconc('co2'); 
cumemi_2020 = 1717.8; #from 1750, source global carbon budget 2022
emi_2020 = 37.1; #GtCO2 
res_2020(box) = emshare(box) * (conc_2020('co2')-conc_preindustrial('co2'));
scaling_forc2x = ( -2.4e-7 * sqr( conc_preindustrial('co2') ) +  7.2e-4 * conc_preindustrial('co2') -  1.05e-4 * ( 2*conc_preindustrial('n2o') ) + 5.36 ) * log( 2 ) / forc2x;
