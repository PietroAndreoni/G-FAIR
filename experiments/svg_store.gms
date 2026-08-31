* Snapshot the current solution into report(vrn,tc,rep,t).
*   %1 = variant label (an element of vrn)
*   %2 = the tc index in scope (the set being looped over: tc0 or tcl)
* Add a quantity here and to the `rep` set in svg.gms to have it reported for
* every scenario at once; the analysis reads whatever `rep` contains.
report('%1',%2,'damages',t)      = DAMAGES.l(t);
report('%1',%2,'damfrac_temp',t) = DAMFRAC_TEMP.l(t);
report('%1',%2,'damfrac_srm',t)  = DAMFRAC_SRM.l(t);
report('%1',%2,'vll',t)          = VLL.l(t);
report('%1',%2,'cost',t)         = COST_SRM.l(t);
report('%1',%2,'tatm',t)         = TATM.l(t);
report('%1',%2,'tatm_ghg',t)     = TATM_GHG.l(t);
report('%1',%2,'forc_srm',t)     = FORC_SRM.l(t);
report('%1',%2,'srm',t)          = SRM.l(t);
report('%1',%2,'qsrm',t)         = Q_SRM.l(t);
report('%1',%2,'tot_forc',t)     = TOT_FORC.l(t);
report('%1',%2,'wemi_co2',t)     = W_EMI.l('co2',t);

* Solver outcome of the solve this snapshot belongs to. Recorded rather than
* aborted on, so one bad state in a loop does not discard the whole run; the
* analysis drops any scenario that is not (solvestat 1, modelstat 1 or 2).
solve_status('%1',%2,'solvestat') = fair.solvestat;
solve_status('%1',%2,'modelstat') = fair.modelstat;
