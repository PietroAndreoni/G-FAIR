# =============================================================================
# svg_discounting.R
#
# The discount curve of the social-value pipeline: a rate r(t) on the model's
# annual grid plus the cumulative factor that turns it into a discount factor
# between any two periods.
#
# Ramsey (the default, SVG_DISCOUNT_MODE = "ramsey"):
#
#   r(t) = rho + eta * g(t),   g(t) = ypc(t+1)/ypc(t) - 1
#
# with ypc the per-capita income path the model already carries (srm_impacts.gms
# builds it from the exogenous GDP and population growth paths, so it is the same
# in every scenario of a run and the curve can be built once). GDP growth holds
# at g0 until 2030 and falls to zero at 2200, population growth to zero at 2100,
# so g starts near 1.5%/yr and reaches zero at 2200: r(t) declines from about 3%
# today to rho thereafter. That decline is the whole point -- a constant 3% and
# this curve agree over the first decades and part company over the tail, which
# is exactly where the commitment question lives.
#
# Constant discounting (SVG_DISCOUNT_MODE = "constant", and every rate of the
# sensitivity grid) is the same object with a flat r, so nothing downstream needs
# to know which of the two it is holding.
#
# Discount factor between two periods, for a curve given on the grid t_1..t_n:
#
#   DF(t | t0) = exp( logdf(t) - logdf(t0) ),   logdf(t_i) = -sum_{j<i} log(1+r_j)
#
# r_j is the rate applied over the interval [t_j, t_j + 1), logdf is interpolated
# linearly between grid points (so the "mid" discounting convention keeps working)
# and extended at the terminal rate beyond them. With a flat r this reproduces
# (1+r)^-(t-t0) exactly, which is what keeps every constant-rate result identical
# to what it was before this module existed.
#
# Consumers take `rate` and pass it to svg_df() / svg_rate_at(), both of which
# also accept a bare numeric as a constant rate, so calls that pass a scalar --
# the closed-form checks in svg_validate.R among them -- are unaffected.
# =============================================================================


# -----------------------------------------------------------------------------
# Construction
# -----------------------------------------------------------------------------

# Constant-rate curve. Carries no grid: the closed form is exact everywhere.
svg_constant_curve <- function(rate = SVG_DISCOUNT_RATE) {
  structure(list(mode = "constant", const = rate, rate0 = rate,
                 rho = NA_real_, eta = NA_real_),
            class = "svg_disc")
}

# Curve from an arbitrary rate path. r[i] applies over [t[i], t[i] + 1).
svg_rate_curve <- function(t, r, mode = "path", ...) {
  o <- order(t); t <- t[o]; r <- r[o]
  if (anyNA(r)) stop("Discount rate path contains NA.", call. = FALSE)
  # logdf[1] = 0: the reference period cancels out of every ratio anyway.
  structure(c(list(mode = mode, const = NULL, t = t, r = r,
                   logdf = c(0, -cumsum(log1p(r)[-length(r)])), rate0 = r[1]),
              list(...)),
            class = "svg_disc")
}

# Per-capita growth rate of the economy on the model grid, from the model's own
# ypc path: g(t) = ypc(t+1)/ypc(t) - 1, the terminal period carrying the last
# resolved rate (the path has plateaued long before it, by construction).
svg_pc_growth <- function(t, ypc) {
  o <- order(t); t <- t[o]; ypc <- ypc[o]
  if (length(t) < 2L) stop("Need at least two periods to derive g(t).", call. = FALSE)
  if (any(!is.finite(ypc)) || any(ypc <= 0))
    stop("Non-positive or non-finite ypc in the discount curve input.", call. = FALSE)
  g <- diff(ypc) / head(ypc, -1)
  data.table(t = t, g = c(g, g[length(g)]))
}

# Ramsey curve r(t) = rho + eta*g(t) from the per-capita income path.
svg_ramsey_curve <- function(t, ypc, rho = SVG_RHO, eta = SVG_ETA) {
  gg <- svg_pc_growth(t, ypc)
  r <- rho + eta * gg$g
  if (any(r <= -1))
    stop("Ramsey rate <= -100%: check rho, eta and the ypc path.", call. = FALSE)
  svg_rate_curve(gg$t, r, mode = "ramsey", rho = rho, eta = eta, g = gg$g)
}

# What every entry point calls. `rate` may be
#   NULL     -> the configured default (Ramsey, off the report's ypc path)
#   numeric  -> a constant rate
#   svg_disc -> itself
# `rep` is the report table, needed only for the Ramsey path.
svg_resolve_discount <- function(rate = NULL, rep = NULL, rho = SVG_RHO,
                                 eta = SVG_ETA, mode = SVG_DISCOUNT_MODE) {
  if (inherits(rate, "svg_disc")) return(rate)
  if (is.numeric(rate)) return(svg_constant_curve(rate))
  if (identical(mode, "constant")) return(svg_constant_curve(SVG_DISCOUNT_RATE))
  if (!identical(mode, "ramsey"))
    stop("SVG_DISCOUNT_MODE must be \"ramsey\" or \"constant\".", call. = FALSE)
  if (is.null(rep) || !"ypc" %in% names(rep))
    stop("Ramsey discounting needs the ypc path: pass the table from ",
         "svg_report_table(), or a numeric rate.", call. = FALSE)
  p <- unique(rep[, .(t, ypc)])
  if (nrow(p) != uniqueN(rep$t))
    stop("ypc is not one path across scenarios; it must be exogenous.", call. = FALSE)
  svg_ramsey_curve(p$t, p$ypc, rho = rho, eta = eta)
}


# -----------------------------------------------------------------------------
# Use
# -----------------------------------------------------------------------------

# Cumulative log discount factor at (possibly non-integer) periods `t`, linearly
# interpolated inside the grid and extended at the end rates outside it.
svg_logdf <- function(disc, t) {
  tg <- disc$t; L <- disc$logdf; n <- length(tg)
  out <- stats::approx(tg, L, xout = pmin(pmax(t, tg[1]), tg[n]))$y
  lo <- t < tg[1]; hi <- t > tg[n]
  if (any(lo)) out[lo] <- L[1] - log1p(disc$r[1]) * (t[lo] - tg[1])
  if (any(hi)) out[hi] <- L[n] - log1p(disc$r[n]) * (t[hi] - tg[n])
  out
}

# Discount factor of period `t` seen from period `t0`. Vectorised over both.
svg_df <- function(disc, t, t0) {
  if (is.numeric(disc)) return((1 + disc)^-(t - t0))
  if (!is.null(disc$const)) return((1 + disc$const)^-(t - t0))
  exp(svg_logdf(disc, t) - svg_logdf(disc, t0))
}

# The annual rate applying at each period, for the identities that need r itself
# (dSCC/dt = r*SCC - MD) rather than a factor.
svg_rate_at <- function(disc, t) {
  if (is.numeric(disc)) return(rep(disc, length(t)))
  if (!is.null(disc$const)) return(rep(disc$const, length(t)))
  stats::approx(disc$t, disc$r, xout = t, rule = 2)$y
}

# The rate at one period as a scalar: what the summary table stores, so code with
# a single number to print still has one.
svg_rate_scalar <- function(disc, t = 1) as.numeric(svg_rate_at(disc, t)[1])

# Short label for figure subtitles and the console report.
svg_discount_label <- function(disc, t0 = NULL, base_year = SVG_BASE_YEAR) {
  pct <- function(z) paste0(format(100 * z, digits = 3), "%")
  if (is.numeric(disc)) return(paste0("r = ", pct(disc)))
  if (identical(disc$mode, "constant")) return(paste0("r = ", pct(disc$const)))
  n <- length(disc$r)
  # The rate is flat over the tail, so name the year it settles rather than the
  # last model period, which would suggest the decline runs to 3020.
  moving <- which(abs(disc$r - disc$r[n]) > 1e-12)
  y_flat <- base_year + if (length(moving)) disc$t[max(moving)] + 1L else disc$t[1]
  r0 <- if (is.null(t0)) disc$r[1] else svg_rate_scalar(disc, t0)
  y0 <- base_year + if (is.null(t0)) disc$t[1] else t0
  sprintf("Ramsey r(t) = %s + %s g(t): %s in %d, %s from %d",
          pct(disc$rho), format(disc$eta), pct(r0), y0, pct(disc$r[n]), y_flat)
}

print.svg_disc <- function(x, ...) {
  cat("svg discount curve:", svg_discount_label(x), "\n")
  invisible(x)
}

# The rate path as a table, for plotting and for the CSV that documents a run.
svg_discount_table <- function(disc, t = NULL, base_year = SVG_BASE_YEAR) {
  if (is.null(t)) t <- if (is.null(disc$const)) disc$t else seq_len(SVG_T_MAX)
  data.table(t = t, year = base_year + t,
             g = if (is.null(disc$const)) stats::approx(disc$t, disc$g, t, rule = 2)$y
                 else NA_real_,
             r = svg_rate_at(disc, t))
}
