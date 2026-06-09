$setglobal pulse_size 1000 #annual avoided methane as a percentage of 2005 ch4 emissions
$setglobal pulse_time 5 #start year of the intervention (1 is 2020)
$setglobal proj_duration 1 #number of years the constant-rate capture project runs from pulse_time (1 = single-year pulse)
$setglobal capture_co2 "ch4co2oxi" #whether the technology also captures the equivalent amount of co2: ch4 | ch4co2oxi (carbon content of the avoided methane) | ch4co2 (carbon content plus co-captured ambient co2) | co2 (same co2 capture as coremoval but without methane removal, for comparison)
$setglobal methane_source "biogenic" #source of the removed methane: biogenic (atmospheric carbon, raises FF_CH4) | fossil (geological carbon, lowers FF_CH4) | atmospheric (ambient mix, FF_CH4 unchanged)


parameter avoided_ch4(t) "avoided methane per project year (MtCH4)";
parameter captured_co2(t) "captured co2 per project year (GtCO2)";

parameter emi_prepulse(ghg,t);
emi_prepulse(ghg,t) = W_EMI.l(ghg,t);

* the project removes methane (and captures co2) at a constant annual rate over
* proj_duration years starting at pulse_time; proj_duration = 1 is a single pulse.
* avoided_ch4(t)/captured_co2(t) are time profiles, so other shapes can be set here.
set proj_window(t) "years over which the capture project is active";
proj_window(t) = yes$(t.val ge %pulse_time% and t.val le %pulse_time% + %proj_duration% - 1);

avoided_ch4(t) = 0;
avoided_ch4(proj_window) = %pulse_size% * 1e-6;

captured_co2(t) = 0;

***** (a) remove methane emissions over the project window
W_EMI.fx('ch4',t)$proj_window(t) = emi_prepulse('ch4',t) - avoided_ch4(t);

* the source of the removed methane sets how the fossil fraction of the remaining
* methane responds
$ifthen.s "%methane_source%"=="biogenic"
* biogenic: carbon was drawn from atmospheric CO2; the fossil mass is held constant,
* so the fossil fraction of the remaining methane increases
FF_CH4.fx(t)$proj_window(t) =
        1 - ( (1 - FF_CH4.l(t)) * emi_prepulse('ch4',t) - avoided_ch4(t) )
          / ( emi_prepulse('ch4',t) - avoided_ch4(t) );
$elseif.s "%methane_source%"=="fossil"
* fossil: carbon is of geological origin; the biogenic mass is held constant, so
* the fossil fraction of the remaining methane decreases
FF_CH4.fx(t)$proj_window(t) =
        ( FF_CH4.l(t) * emi_prepulse('ch4',t) - avoided_ch4(t) )
          / ( emi_prepulse('ch4',t) - avoided_ch4(t) );
$elseif.s "%methane_source%"=="atmospheric"
* atmospheric: methane is taken from ambient air with its average composition, so
* the fossil fraction of the remaining methane is left unchanged
FF_CH4.fx(t)$proj_window(t) = FF_CH4.l(t);
$endif.s

***** (b) capture the equivalent amount of co2 (carbon content of the avoided methane)
$ifthen.c "%capture_co2%"=="ch4co2oxi"
captured_co2(t)$proj_window(t) = conv_frac * 1e-3 * ghg_mm('co2') / ghg_mm('ch4') * avoided_ch4(t);
W_EMI.fx('co2',t)$proj_window(t) = emi_prepulse('co2',t) - captured_co2(t);
$elseif.c "%capture_co2%"=="ch4co2"
* carbon content of the avoided methane (as CO2) plus the ambient CO2 co-captured
* from the same air stream: per molecule of CH4 removed, the processed air also
* carries CONC(co2)/CONC(ch4) molecules of CO2, captured at 100% efficiency. The
* molar ratio uses the baseline concentrations of each project year (ppm/ppb -> mole fractions).
captured_co2(t)$proj_window(t) =
        conv_frac * 1e-3 * ghg_mm('co2') / ghg_mm('ch4') * avoided_ch4(t)
      + 1e-3 * ghg_mm('co2') / ghg_mm('ch4') * avoided_ch4(t) *
        ( CONC.l('co2',t) * 1e-6 ) / ( CONC.l('ch4',t) * 1e-9 );
W_EMI.fx('co2',t)$proj_window(t) = emi_prepulse('co2',t) - captured_co2(t);
$elseif.c "%capture_co2%"=="co2"
* for comparison, only the ambient CO2 co-capture (as in "coremoval") with the
* methane left in place, to isolate the direct-air-capture effect
captured_co2(t)$proj_window(t) = 1e-3 * ghg_mm('co2') / ghg_mm('ch4') * avoided_ch4(t) *
        ( CONC.l('co2',t) * 1e-6 ) / ( CONC.l('ch4',t) * 1e-9 );
W_EMI.fx('co2',t)$proj_window(t) = emi_prepulse('co2',t) - captured_co2(t);
W_EMI.fx('ch4',t)$proj_window(t) = emi_prepulse('ch4',t);
FF_CH4.fx(t)$proj_window(t) = 1 - ( (1 - FF_CH4.l(t)) * emi_prepulse('ch4',t) )
          / ( emi_prepulse('ch4',t) );
$endif.c

solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

execute_unload "%results_folder%/%rcp%_EXPremoval%capture_co2%_PT%pulse_time%_PS%pulse_size%_PD%proj_duration%_SRC%methane_source%.gdx";
