#add stratospheric SAI as short-lived GHG
set ghg /"sai"/;
oghg("sai") = yes;

taughg("sai") = 0.8;
ghg_mm("sai") = 64.066; 

forcing_coeff("sai") = -0.2; #12 tG/S ->  

conc0("sai") = 0;
preindustrial_conc("sai") = 0;
natural_emissions("sai",t) = 0;

equation eq_saiinertiaup, eq_saiinertiadown;

eq_saiinertiaup(t+1).. W_EMI("sai",t+1) =L= 1.2*W_EMI("sai",t);
eq_saiinertiadown(t+1).. W_EMI("sai",t+1) =G= 0.8*W_EMI("sai",t);