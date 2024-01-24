#add stratospheric SAI as short-lived GHG
set ghg /"sai"/;
oghg("sai") = yes;

taughg("sai") = 2;
ghg_mm("sai") = 64.066; 

forcing_coeff("sai") = -0.2; #12 tG/S ->  

conc0("sai") = 0;
preindustrial_conc("sai") = 0;
natural_emissions("sai",t) = 0;
inertia("sai") = 0.2;

W_EMI.fx('sai',tfirst) = 0;