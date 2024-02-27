#add stratospheric SAI as short-lived GHG
set ghg /"sai"/;
oghg("sai") = yes;

taughg("sai") = 0.8;
ghg_mm("sai") = 64.066; 
forcing_coeff("sai") = -0.2; #12 tG/S ->  

conc_2020("sai") = 0;
conc_preindustrial("sai") = 0;
natural_emissions("sai",t) = 0;

W_EMI.fx('sai',tfirst) = 0;