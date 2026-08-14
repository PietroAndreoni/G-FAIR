$setglobal rate_of_cooling 10 #in deg/millennium, by default half of current warming (0.1 deg/decade). Weird unit is for reporting
$setglobal start_rampup 2025
$setglobal end_rampup 2100
$setglobal start_rampdown 3000
$setglobal end_rampdown 3100
$setglobal pulse_size 100 #percentage increase in emissions
$setglobal pulse_time 5 #year of pulse (1 is 2020)

* save initial state
parameter wemi_0(ghg,t), ffch4_0(t);
wemi_0(ghg,t) = W_EMI.l(ghg,t);
ffch4_0(t) = FF_CH4.l(t);

** Reporting. Every scenario below writes into one parameter and nothing is
** unloaded until the end: the five scenario families span 39 solves, which as
** full gdx dumps would be ~140 Mb of results that are always read together.
** experiments/svg_store.gms holds the per-solve snapshot, so a new reported
** quantity is added in two places only (the `rep` set here and that file).
set tc  "commitment length / pulse delay in years (0 = the reference experiments)"
        / 0, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500 /;
set tc0(tc) "index used by the reference experiments" / 0 /;
set tcl(tc) "lengths used by the commitment loops"
        / 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500 /;

set vrn "experiment variant" /
    base     "background SAI, no emission pulse"
    ghgpulse "background SAI, CO2 pulse at the pulse period, SAI not adjusted"
    srmpulse "marginal SAI that offsets the pulse, never terminated"
    term     "that marginal SAI, terminated tc years after the pulse"
    scc      "background SAI, CO2 pulse emitted tc years after the pulse period" /;

set rep "reported quantity" / damages, damfrac_temp, damfrac_srm, vll, cost,
                              tatm, tatm_ghg, forc_srm, srm, qsrm, wemi_co2 /;

parameter report(vrn,tc,rep,t) "results by variant, commitment length and quantity";

** Scenario 1: with SRM, without emission pulse

* build SRM strategy (linear ramp-up and down, flat in between)
background_srm(t)$(2020 + t.val ge %start_rampup% and 2020 + t.val le %end_rampup%) = %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampup%) / 1e3 ; # 0.1 deg/decade
background_srm(t)$(2020 + t.val gt %end_rampup% and 2020 + t.val le %start_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3;
background_srm(t)$(2020 + t.val gt %start_rampdown% and 2020 + t.val le %end_rampdown%) = %rate_of_cooling% * forc2x / Tecs * (%end_rampup% - %start_rampup%) / 1e3 - %rate_of_cooling% * forc2x / Tecs * (2020 + t.val - %start_rampdown%) / 1e3;
background_srm(t)$(background_srm(t) le 0) = 0;
solve fair using nlp minimizing OBJ;

target_temp(t) = TATM.l(t);

loop(tc0,
$batinclude "experiments/svg_store.gms" base tc0
);

** Scenario 2: with SRM, with emission pulse (but no masking)
W_EMI.fx('co2','%pulse_time%') = W_EMI.l('co2','%pulse_time%') + %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC);

solve fair using nlp minimizing OBJ;

loop(tc0,
$batinclude "experiments/svg_store.gms" ghgpulse tc0
);

** Scenario 3: with additional SRM, no emission pulse
SRM.lo(t) = - background_srm(t); SRM.up(t) = +inf; # reduce the SAI by at most the full amount
SRM.fx(t)$(t.val lt %pulse_time%) = 0; # and cannot change SAI before the pulse

* solve 3 times as this is the only non-simulation problem
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
solve fair using nlp minimizing OBJ;
abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "Base model is not solving";

SRM.fx(t) = SRM.l(t);
W_EMI.fx(ghg,t) = wemi_0(ghg,t);

solve fair using nlp minimizing OBJ;

loop(tc0,
$batinclude "experiments/svg_store.gms" srmpulse tc0
);

** Scenario 4: the scenario-3 SAI, switched off tc years after the pulse (abrupt
** termination; FORC_SRM still decays with tausrm). Baseline emissions throughout,
** so this prices the marginal SAI alone, rebound included.
parameter srm_opt(t) "marginal SAI forcing solved for in scenario 3";
srm_opt(t) = SRM.l(t);

loop(tcl,
  W_EMI.fx(ghg,t) = wemi_0(ghg,t);
  SRM.fx(t) = srm_opt(t)$( t.val le %pulse_time% + tcl.val );

  solve fair using nlp minimizing OBJ;
  abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "svg termination run is not solving";

$batinclude "experiments/svg_store.gms" term tcl
);

** Scenario 5: no marginal SAI, background schedule only, with the CO2 pulse moved
** tc years later. Differenced against the baseline this gives SCC(t) along the
** background SAI deployment.
SRM.fx(t) = 0;

loop(tcl,
  W_EMI.fx(ghg,t) = wemi_0(ghg,t);
  W_EMI.fx('co2',t)$( t.val eq %pulse_time% + tcl.val ) = wemi_0('co2',t) + %pulse_size%/100 * (emissions_rcp('2005','%rcp%','co2')/CO2toC);

  solve fair using nlp minimizing OBJ;
  abort$(not (fair.solvestat eq 1 and (fair.modelstat eq 1 or fair.modelstat eq 2))) "svg delayed-pulse run is not solving";

$batinclude "experiments/svg_store.gms" scc tcl
);

** Single unload, and only the reporting symbols: the model variables a full dump
** would carry are the last solve's alone, and every scenario's are in `report`.
execute_unload "%results_folder%/%rcp%_EXP%experiment%_GAS%gas%_ECS%ecs%_TCR%tcr%_PT%pulse_time%_RC%rate_of_cooling%_EC%end_rampdown%_BC%start_rampup%",
    report, vrn, tc, tc0, tcl, rep,
    y, pop, ypc, g, gpop, background_srm, target_temp, forcing_exogenous,
    a0, b0, srm_angle, tgtoforc, vsl0, ypc0, eta_vsl, mort_srm, csrm, eff_decline_srm,
    Tecs, Ttcr, forc2x, CO2toC, tstep, yr_2020;
