# =============================================================================
# svg_hazard.R
#
# Survival weighting, the social value of a commitment of length H, and the
# optimal commitment horizon (implementation spec sections 2.10 - 2.12).
#
#   SV(H;t0) = sum over the first H periods of
#                DF * S * (annual operating flow)          [benefit while active]
#              - DF * (S_start - S_end) * DeltaL           [expected collapse cost]
#
# The operating flow is the model's own damage difference (marginal SAI vs the
# background run), so it already contains the climate benefit, the imperfect-
# offsetting penalty, the sulfate mortality and the deployment cost. Nothing is
# charged at the planned end H: the hazard applies only while the intervention is
# active (spec section 11.4).
#
# Termination is weighted by the exact interval probability S(tau) - S(tau+dt) of
# a constant Poisson hazard, not by h*dt (spec section 2.11), so the numbers are
# right at h = 0.01 and dt = 1 rather than only in the small-dt limit. At h = 0
# every survival weight is 1 and the whole thing collapses to the deterministic
# commitment profile, which is what makes the two directly comparable.
# =============================================================================


# Survival of the intervention to elapsed time tau under a constant hazard.
svg_survival <- function(tau, h) exp(-h * tau)

# Discount factor for a period whose start is tau years after t0.
#   rate : a discount curve (svg_discounting.R) or a constant rate
#   t0   : the model period tau is measured from. Only a declining curve needs
#          it -- with a constant rate the factor depends on tau alone, which is
#          why the default keeps the old signature working.
svg_discount_factor <- function(tau, rate, at = SVG_DISCOUNT_AT, dt = 1, t0 = 0) {
  off <- switch(at, start = 0, mid = dt / 2,
                stop("SVG_DISCOUNT_AT must be \"start\" or \"mid\".", call. = FALSE))
  svg_df(rate, t0 + tau + off, t0)
}

# Expected discounted value of committing for H = 0, 1, ... periods from t0.
#   flow    : annual operating flow [trillion USD/yr], one entry per period from
#             t0 onward, in period order
#   deltaL  : marginal termination loss [trillion USD] for the same periods
#   h       : Poisson hazard [1/yr]; h = 0 gives the deterministic case
#   rate    : a discount curve (svg_discounting.R) or a constant rate
#   t0      : model period the flows start at, which a declining curve needs to
#             know where on itself the commitment sits
# Returns a data.table with one row per horizon H (including H = 0, value 0):
#   contribution   the period's own expected discounted contribution
#   sv             cumulative SV(H)
#   operating, termination   the two halves of the cumulative sum
svg_commitment_value <- function(flow, deltaL, h, rate = SVG_DISCOUNT_RATE,
                                 dt = 1, at = SVG_DISCOUNT_AT, t0 = 0) {
  n <- length(flow)
  if (length(deltaL) != n) stop("flow and deltaL must have the same length.", call. = FALSE)
  tau <- (seq_len(n) - 1) * dt            # elapsed time at the START of each period
  df  <- svg_discount_factor(tau, rate, at, dt, t0 = t0)
  s0  <- svg_survival(tau, h)             # survival into the period
  s1  <- svg_survival(tau + dt, h)        # survival out of it
  # Operating benefit is earned while active; with "start" discounting the weight
  # is the survival into the period, with "mid" the survival at its midpoint.
  s_op <- if (identical(at, "mid")) svg_survival(tau + dt / 2, h) else s0
  operating   <- df * s_op * flow * dt
  termination <- df * (s0 - s1) * deltaL
  data.table(H = 0:n,
             contribution = c(0, operating - termination),
             sv           = c(0, cumsum(operating - termination)),
             operating    = c(0, cumsum(operating)),
             termination  = c(0, cumsum(termination)))
}

# Global optimum over the commitment horizon, including H = 0.
# Classification (spec section 2.12):
#   zero      SV is maximised by not intervening at all
#   finite    an interior maximum
#   infinite  the maximum is at the last horizon AND the tail has converged, so
#             extending the simulation would not move it
#   boundary  the maximum is at the last horizon and the tail has NOT converged:
#             the result is pinned to the simulation end and must not be read as
#             an infinite commitment
svg_optimize_commitment <- function(sv, tail_frac = SVG_TAIL_FRAC, tail_tol = SVG_TAIL_TOL) {
  if (!nrow(sv)) stop("Empty commitment table.", call. = FALSE)
  i <- which.max(sv$sv)
  H_star  <- sv$H[i]
  SV_star <- sv$sv[i]
  last_H  <- sv$H[nrow(sv)]
  tail_n  <- max(1L, ceiling(tail_frac * nrow(sv)))
  tail_gain <- sv$sv[nrow(sv)] - sv$sv[nrow(sv) - tail_n + 1L]
  converged <- abs(tail_gain) <= tail_tol * max(abs(SV_star), .Machine$double.eps)
  cls <- if (H_star == 0) "zero"
         else if (H_star < last_H) "finite"
         else if (converged) "infinite"
         else "boundary"
  list(H_star = H_star, SV_star = SV_star, class = cls,
       tail_gain = tail_gain, tail_converged = converged,
       sv_last = sv$sv[nrow(sv)], last_H = last_H)
}
