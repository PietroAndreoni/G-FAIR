# =============================================================================
# svg_termination.R
#
# Marginal termination loss of the one-tonne-equivalent SAI intervention, in the
# reduced form of the implementation spec (section 2.8).
#
# The hazard is an abrupt, unplanned collapse of the whole SAI programme. Total
# termination damage is convex in the masked warming the programme is carrying:
#
#   L_total(M,t) = beta * SCC(t) * (M_ref / alpha) * (M(t) / M_ref)^n,   n = p*q
#
# where p is the convexity of damages in the rebound velocity and q the scaling
# of that velocity with masked warming. The marginal intervention adds alpha to
# masked warming, so its own contribution to the loss is
#
#   DeltaL(t) = alpha * dL_total/dM = beta * n * SCC(t) * (M(t) / M_ref)^(n-1)
#
# which is free of alpha, is zero at M = 0, and equals 8*SCC at M = 1 degC under
# the central calibration (beta = 4, p = 2, q = 1). DeltaL is a marginal loss per
# tonne-equivalent; L_total is an absolute global loss. Do not mix them.
#
# The spec allows this imposed form to be replaced later by the derivative of an
# empirically fitted L(M, state) from explicit programme-termination branches of
# the GAMS model (spec sections 2.9 / 3.7). svg_termination_loss_empirical() is
# the hook for that; it is unused until those runs exist.
# =============================================================================


# Marginal termination loss [USD/tCO2], vectorised over M and scc.
#   M    : masked warming at the state, T_GHG - TATM [degC]
#   scc  : social cost of carbon at the same state [USD/tCO2]
svg_marginal_termination_loss <- function(M, scc,
                                          beta = SVG_BETA, p = SVG_P, q = SVG_Q,
                                          M_ref = SVG_M_REF) {
  # On the rc = 0 pathway there is no background aerosol, so M = T_GHG - TATM is
  # zero by construction and comes back off the solver at +-1e-16. Only a
  # genuinely negative M -- SAI that warms -- is worth a warning.
  if (any(M < -1e-9, na.rm = TRUE))
    warning("Negative masked warming passed to svg_marginal_termination_loss().")
  n <- p * q
  beta * n * scc * (pmax(M, 0) / M_ref)^(n - 1)
}

# Total (absolute) termination loss of the programme, for diagnostics. Needs the
# marginal temperature effect of a tonne, alpha [degC/tCO2], which the marginal
# form does not.
svg_total_termination_loss <- function(M, scc, alpha,
                                       beta = SVG_BETA, p = SVG_P, q = SVG_Q,
                                       M_ref = SVG_M_REF) {
  n <- p * q
  beta * scc * (M_ref / alpha) * (pmax(M, 0) / M_ref)^n
}

# Placeholder for the empirical calibration (spec 2.9): given a fitted
# L_total(M) and its derivative, return DeltaL = alpha * dL/dM. Kept so the
# switch is a one-line change in svg_value_by_state() if the programme
# termination branches are ever run.
svg_termination_loss_empirical <- function(M, fit) {
  stop("Empirical termination calibration is not available: run the programme-",
       "termination branches first (spec sections 2.9 / 3.7).", call. = FALSE)
}

# The masked warming at which the expected annual termination cost h*DeltaL
# overtakes the annual operating benefit, given a marginal-damage flow md
# [USD/tCO2/yr] and the SCC at the same state. Beyond it the marginal
# intervention has a negative flow and the optimal commitment collapses to zero.
# With n = 2 this is M* = md / (beta*n*h*scc).
svg_termination_crossover_M <- function(md, scc, h,
                                        beta = SVG_BETA, p = SVG_P, q = SVG_Q,
                                        M_ref = SVG_M_REF) {
  n <- p * q
  if (n != 2) return(NA_real_)   # closed form only for the linear-in-M case
  M_ref * md / (beta * n * h * scc)
}
