**** PRE MODEL 1
************ find QSLOW and QFAST given climate specifications
EQUATIONS eq_tecs, eq_ttcr;
** calculate Qi imposing ECS and TCR (more efficient as a parameter)
eq_tecs..          Tecs =E= forc2x * (QSLOW + QFAST); 

eq_ttcr..          Ttcr =E= forc2x * (QSLOW * (1 - dslow/69.7 * (1 - exp(-69.7/dslow)) ) +
                            QFAST * (1 - dfast/69.7 * (1 - exp(69.7/dfast)) ) ) ; 

model solveqs  / eq_tecs, eq_ttcr /;
solve solveqs using cns; 

QSLOW.fx = QSLOW.l; 
QFAST.fx = QFAST.l;
*************** end QSLOW and QFAST
