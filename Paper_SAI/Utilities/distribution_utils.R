# =============================================================================
# distribution_utils.R
#
# Shared distribution-fitting / sampling helpers for the Monte Carlo pipeline.
# Sourced by Generate_montecarlo.R and by test_fit_distribution.R so the two
# never drift. Base R only (no package dependencies).
#
# Contents:
#   fit_distribution()       closed-form normal/lognormal fit from summary stats
#   qlnorm_fit / qnorm_fit   inverse-CDF samplers (fit -> quantile of uniforms)
#   qtrunc_inv + *_trunc     truncated inverse-CDF (no boundary spike; keeps the
#                            Sobol low-discrepancy structure)
#   dsnorm/psnorm/qsnorm     split-normal (two-piece normal) distribution
#   fit_quantile_dist()      infer the best-fitting family from quantile points
# =============================================================================

# -----------------------------------------------------------------------------
# Closed-form normal / lognormal fit from a subset of summary statistics.
# -----------------------------------------------------------------------------
fit_distribution <- function(distribution = "lognormal",
                             n      = 1,      # number of required realizations
                             median = NULL,
                             mean   = NULL,
                             sd     = NULL,
                             q33    = NULL,
                             q66    = NULL,
                             q5     = NULL,
                             q95    = NULL,
                             probs_33_66 = c(0.33, 0.66),
                             probs_5_95  = c(0.05, 0.95),
                             plot = FALSE,
                             return_params = FALSE) {

  mu    <- NA_real_
  sigma <- NA_real_

  # Track which inputs the caller supplied, so we can warn about any that the
  # matched branch ends up ignoring (the if/else chains use a fixed precedence).
  supplied <- names(Filter(Negate(is.null),
                           list(median = median, mean = mean, sd = sd,
                                q33 = q33, q66 = q66, q5 = q5, q95 = q95)))
  used <- character(0)

  # Warn when an over-determined (median + two quantiles) request cannot be met
  # by a single distribution, i.e. the two implied sigma estimates disagree.
  warn_inconsistent <- function(sig_lo, sig_hi, where) {
    if (is.finite(sig_lo) && is.finite(sig_hi)) {
      rel <- abs(sig_hi - sig_lo) / mean(c(abs(sig_lo), abs(sig_hi)))
      if (isTRUE(rel > 0.05))
        warning(sprintf("fit_distribution: %s inputs are mutually inconsistent (sigma estimates differ by %.0f%%); returning a least-discrepancy compromise.",
                        where, 100 * rel), call. = FALSE)
    }
  }

  if (distribution == "lognormal") {
    # Lognormal only makes sense for strictly positive inputs.
    check_pos <- function(x, name) {
      if (!is.null(x) && any(x <= 0)) stop(name, " must be > 0 for a lognormal distribution.")
    }
    check_pos(median, "median"); check_pos(mean, "mean"); check_pos(sd, "sd")
    check_pos(q33, "q33"); check_pos(q66, "q66"); check_pos(q5, "q5"); check_pos(q95, "q95")

    ## 1) mean and sd of X
    if (!is.null(mean) && !is.null(sd)) {
      used <- c("mean","sd")
      if (sd <= 0) stop("sd must be > 0.")
      v      <- sd^2
      sigma2 <- log(1 + v / mean^2)
      if (sigma2 <= 0) stop("Inconsistent mean and sd for a lognormal.")
      sigma  <- sqrt(sigma2)
      mu     <- log(mean) - 0.5 * sigma2

      ## 2) mean and median of X
    } else if (!is.null(mean) && !is.null(median)) {
      used <- c("mean","median")
      mu     <- log(median)
      sigma2 <- 2 * (log(mean) - mu)
      if (sigma2 <= 0) stop("Inconsistent mean and median for a lognormal.")
      sigma  <- sqrt(sigma2)

      ## 3) median and sd of X
    } else if (!is.null(median) && !is.null(sd)) {
      used <- c("median","sd")
      if (sd <= 0) stop("sd must be > 0.")
      mu <- log(median)
      v  <- sd^2
      # var = (t-1)*t*exp(2 mu) with t = exp(sigma^2); solve t^2 - t - A = 0.
      A    <- v / median^2
      t    <- (1 + sqrt(1 + 4 * A)) / 2
      sigma2 <- log(t)
      if (sigma2 <= 0) stop("Inconsistent median and sd for a lognormal.")
      sigma  <- sqrt(sigma2)

      ## 4) median and 33/66 quantiles
    } else if (!is.null(median) && !is.null(q33) && !is.null(q66)) {
      used <- c("median","q33","q66")
      mu <- log(median)
      z1 <- qnorm(probs_33_66[1]); z2 <- qnorm(probs_33_66[2])
      sigma_33 <- (log(q33) - mu) / z1
      sigma_66 <- (log(q66) - mu) / z2
      sigma    <- mean(c(sigma_33, sigma_66))
      warn_inconsistent(sigma_33, sigma_66, "lognormal median+q33+q66")

      ## 5) median and 5/95 quantiles
    } else if (!is.null(median) && !is.null(q5) && !is.null(q95)) {
      used <- c("median","q5","q95")
      mu <- log(median)
      z1 <- qnorm(probs_5_95[1]); z2 <- qnorm(probs_5_95[2])
      sigma_5  <- (log(q5)  - mu) / z1
      sigma_95 <- (log(q95) - mu) / z2
      sigma    <- mean(c(sigma_5, sigma_95))
      warn_inconsistent(sigma_5, sigma_95, "lognormal median+q5+q95")

      ## 6) 33/66 quantiles only
    } else if (!is.null(q33) && !is.null(q66)) {
      used <- c("q33","q66")
      z1 <- qnorm(probs_33_66[1]); z2 <- qnorm(probs_33_66[2])
      L1 <- log(q33); L2 <- log(q66)
      sigma <- (L2 - L1) / (z2 - z1)
      mu    <- L1 - sigma * z1

      ## 7) 5/95 quantiles only
    } else if (!is.null(q5) && !is.null(q95)) {
      used <- c("q5","q95")
      z1 <- qnorm(probs_5_95[1]); z2 <- qnorm(probs_5_95[2])
      L1 <- log(q5); L2 <- log(q95)
      sigma <- (L2 - L1) / (z2 - z1)
      mu    <- L1 - sigma * z1

    } else {
      stop("Not enough information for lognormal. Provide one of:
      (mean & sd), (mean & median), (median & sd),
      (median & q33 & q66), (median & q5 & q95), (q33 & q66), or (q5 & q95).")
    }

    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0)
      stop("Invalid lognormal parameters inferred from inputs.")

    extra <- setdiff(supplied, used)
    if (length(extra) > 0)
      warning(sprintf("fit_distribution: ignored input(s) %s; the matched lognormal case used (%s).",
                      paste(extra, collapse = ", "), paste(used, collapse = " & ")), call. = FALSE)

    out_mean   <- exp(mu + 0.5 * sigma^2)
    out_median <- exp(mu)
    out_sd     <- sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))

    if (return_params) return(list(mu = mu, sigma = sigma, mean = out_mean, median = out_median, sd = out_sd))
    if (plot) plot(density(rlnorm(10000, mu, sigma)))
    return(rlnorm(n, mu, sigma))

  } else if (distribution == "normal") {

    if (!is.null(sd) && sd <= 0) stop("sd must be > 0.")

    ## 1) mean & sd
    if (!is.null(mean) && !is.null(sd)) {
      used <- c("mean","sd"); mu <- mean; sigma <- sd

      ## 2) median & sd
    } else if (!is.null(median) && !is.null(sd)) {
      used <- c("median","sd"); mu <- median; sigma <- sd

      ## 3) median and 33/66 quantiles
    } else if (!is.null(median) && !is.null(q33) && !is.null(q66)) {
      used <- c("median","q33","q66")
      mu <- median
      z1 <- qnorm(probs_33_66[1]); z2 <- qnorm(probs_33_66[2])
      sigma_33 <- (q33 - mu) / z1
      sigma_66 <- (q66 - mu) / z2
      sigma    <- mean(c(sigma_33, sigma_66))
      warn_inconsistent(sigma_33, sigma_66, "normal median+q33+q66")

      ## 4) median and 5/95 quantiles
    } else if (!is.null(median) && !is.null(q5) && !is.null(q95)) {
      used <- c("median","q5","q95")
      mu <- median
      z1 <- qnorm(probs_5_95[1]); z2 <- qnorm(probs_5_95[2])
      sigma_5  <- (q5  - mu) / z1
      sigma_95 <- (q95 - mu) / z2
      sigma    <- mean(c(sigma_5, sigma_95))
      warn_inconsistent(sigma_5, sigma_95, "normal median+q5+q95")

      ## 5) 33/66 quantiles only
    } else if (!is.null(q33) && !is.null(q66)) {
      used <- c("q33","q66")
      z1 <- qnorm(probs_33_66[1]); z2 <- qnorm(probs_33_66[2])
      sigma <- (q66 - q33) / (z2 - z1)
      mu    <- q33 - sigma * z1

      ## 6) 5/95 quantiles only
    } else if (!is.null(q5) && !is.null(q95)) {
      used <- c("q5","q95")
      z1 <- qnorm(probs_5_95[1]); z2 <- qnorm(probs_5_95[2])
      sigma <- (q95 - q5) / (z2 - z1)
      mu    <- q5 - sigma * z1

    } else {
      stop("Not enough information for normal. Provide one of:
      (mean & sd), (median & sd), (median & q33 & q66),
      (median & q5 & q95), (q33 & q66), or (q5 & q95).")
    }

    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0)
      stop("Invalid normal parameters inferred from inputs.")

    extra <- setdiff(supplied, used)
    if (length(extra) > 0)
      warning(sprintf("fit_distribution: ignored input(s) %s; the matched normal case used (%s).",
                      paste(extra, collapse = ", "), paste(used, collapse = " & ")), call. = FALSE)

    if (return_params) return(list(mu = mu, sigma = sigma, mean = mu, median = mu, sd = sigma))
    if (plot) plot(density(rnorm(10000, mu, sigma)))
    return(rnorm(n, mu, sigma))

  } else stop("choose normal or lognormal distribution")
}

# -----------------------------------------------------------------------------
# Inverse-CDF samplers: fit (mu,sigma) then evaluate the quantile of uniforms u.
# Driven by the shared Sobol sequence so the whole realization is low-discrepancy.
# -----------------------------------------------------------------------------
qlnorm_fit <- function(u, ...) {
  p <- fit_distribution("lognormal", ..., return_params = TRUE)
  qlnorm(u, p$mu, p$sigma)
}
qnorm_fit <- function(u, ...) {
  p <- fit_distribution("normal", ..., return_params = TRUE)
  qnorm(u, p$mu, p$sigma)
}

# -----------------------------------------------------------------------------
# Truncated inverse-CDF. Maps u in (0,1) into the [lo,hi] truncation of a base
# distribution given its CDF (pfun) and quantile (qfun). Unlike pmin/pmax
# clamping this redistributes the tail mass smoothly (no boundary spike) and
# stays a monotone map of u, preserving the Sobol low-discrepancy structure.
# -----------------------------------------------------------------------------
qtrunc_inv <- function(u, pfun, qfun, lo = -Inf, hi = Inf) {
  plo <- pfun(lo); phi <- pfun(hi)
  qfun(plo + u * (phi - plo))
}
qlnorm_fit_trunc <- function(u, lo = 0, hi = Inf, ...) {
  p <- fit_distribution("lognormal", ..., return_params = TRUE)
  qtrunc_inv(u, function(x) plnorm(x, p$mu, p$sigma),
                function(pp) qlnorm(pp, p$mu, p$sigma), lo, hi)
}
qnorm_fit_trunc <- function(u, lo = -Inf, hi = Inf, ...) {
  p <- fit_distribution("normal", ..., return_params = TRUE)
  qtrunc_inv(u, function(x) pnorm(x, p$mu, p$sigma),
                function(pp) qnorm(pp, p$mu, p$sigma), lo, hi)
}

# -----------------------------------------------------------------------------
# Split-normal (two-piece normal): mode m, left sd sL, right sd sR. A 3-parameter
# family that can match an asymmetric median + q5 + q95 exactly.
# -----------------------------------------------------------------------------
dsnorm <- function(x, m, sL, sR) {
  A <- sqrt(2 / pi) / (sL + sR)
  ifelse(x < m, A * exp(-(x - m)^2 / (2 * sL^2)),
                A * exp(-(x - m)^2 / (2 * sR^2)))
}
psnorm <- function(x, m, sL, sR) {
  ifelse(x < m, (2 * sL / (sL + sR)) * pnorm((x - m) / sL),
                sL / (sL + sR) + (2 * sR / (sL + sR)) * (pnorm((x - m) / sR) - 0.5))
}
qsnorm <- function(p, m, sL, sR) {
  w   <- sL / (sL + sR)                # P(X < m)
  out <- numeric(length(p))
  lo  <- p <= w
  out[lo]  <- m + sL * qnorm(p[lo] * (sL + sR) / (2 * sL))
  out[!lo] <- m + sR * qnorm(0.5 + (p[!lo] * (sL + sR) - sL) / (2 * sR))
  out
}

# -----------------------------------------------------------------------------
# Infer the best-fitting distribution from quantile information.
#
# Supply any >=2 of median / q5 / q95 / q33 / q66. Each candidate family is fit
# by minimizing the relative squared error to ALL supplied (p, q) points; the
# family with the smallest residual is returned. 3-parameter families
# (split-normal) can honor median+q5+q95 exactly where normal/lognormal cannot.
#
# Returns a list with: best (family name), params (interpretable), q()/d()/r()
# closures (quantile / density / inverse-CDF sampler from uniforms), residual,
# and ranking (data.frame of all families with their residual + params).
# -----------------------------------------------------------------------------
fit_quantile_dist <- function(median = NULL, q5 = NULL, q95 = NULL,
                              q33 = NULL, q66 = NULL,
                              families = c("normal","lognormal","split-normal","log-split-normal")) {
  prob_of <- c(median = 0.5, q5 = 0.05, q95 = 0.95, q33 = 0.33, q66 = 0.66)
  pts <- c(median = median, q5 = q5, q95 = q95, q33 = q33, q66 = q66)
  pts <- pts[!vapply(pts, is.null, logical(1))]
  pts <- unlist(pts)
  if (length(pts) < 2) stop("fit_quantile_dist: supply at least two of median/q5/q95/q33/q66.")
  p  <- unname(prob_of[names(pts)])
  q  <- unname(pts)
  o  <- order(p); p <- p[o]; q <- q[o]
  positive <- all(q > 0)
  rel <- function(fit) sqrt(mean(((fit - q) / q)^2))   # relative RMSE across points

  # crude initial spread on the (log-)scale from the widest available pair
  spread <- function(v) {
    zr <- qnorm(max(p)) - qnorm(min(p))
    max((v[length(v)] - v[1]) / zr, 1e-6)
  }
  m0 <- if (!is.null(median)) median else stats::approx(p, q, 0.5, rule = 2)$y

  specs <- list(
    "normal" = list(
      init = c(mu = m0, ls = log(spread(q))),
      Q = function(th, pr) th[1] + exp(th[2]) * qnorm(pr),
      par = function(th) c(mu = unname(th[1]), sigma = exp(unname(th[2]))),
      dens = function(th) { pa <- c(th[1], exp(th[2])); function(x) dnorm(x, pa[1], pa[2]) },
      ok = TRUE),
    "lognormal" = list(
      init = c(mu = log(m0), ls = log(spread(log(pmax(q, 1e-12))))),
      Q = function(th, pr) exp(th[1] + exp(th[2]) * qnorm(pr)),
      par = function(th) c(meanlog = unname(th[1]), sdlog = exp(unname(th[2]))),
      dens = function(th) { pa <- c(th[1], exp(th[2])); function(x) dlnorm(x, pa[1], pa[2]) },
      ok = positive),
    "split-normal" = list(
      init = c(m = m0, lL = log(spread(q)), lR = log(spread(q))),
      Q = function(th, pr) qsnorm(pr, th[1], exp(th[2]), exp(th[3])),
      par = function(th) c(mode = unname(th[1]), sdL = exp(unname(th[2])), sdR = exp(unname(th[3]))),
      dens = function(th) { pa <- c(th[1], exp(th[2]), exp(th[3])); function(x) dsnorm(x, pa[1], pa[2], pa[3]) },
      ok = TRUE),
    "log-split-normal" = list(
      init = c(m = log(m0), lL = log(spread(log(pmax(q, 1e-12)))), lR = log(spread(log(pmax(q, 1e-12))))),
      Q = function(th, pr) exp(qsnorm(pr, th[1], exp(th[2]), exp(th[3]))),
      par = function(th) c(mode_log = unname(th[1]), sdL_log = exp(unname(th[2])), sdR_log = exp(unname(th[3]))),
      dens = function(th) { pa <- c(th[1], exp(th[2]), exp(th[3]))
                            function(x) ifelse(x > 0, dsnorm(log(x), pa[1], pa[2], pa[3]) / x, 0) },
      ok = positive)
  )

  fit_one <- function(sp) {
    obj <- function(th) {
      val <- tryCatch(sp$Q(th, p), error = function(e) rep(NA_real_, length(p)))
      if (any(!is.finite(val))) return(1e12)
      sum(((val - q) / q)^2)
    }
    opt <- tryCatch(optim(sp$init, obj, method = "Nelder-Mead",
                          control = list(maxit = 2000, reltol = 1e-10)),
                    error = function(e) NULL)
    if (is.null(opt) || !is.finite(opt$value)) return(NULL)
    th <- opt$par
    list(residual = rel(sp$Q(th, p)), theta = th,
         Q = function(pr) sp$Q(th, pr), d = sp$dens(th), par = sp$par(th))
  }

  fitted <- list()
  for (fam in families) {
    sp <- specs[[fam]]
    if (is.null(sp) || !isTRUE(sp$ok)) next
    f <- fit_one(sp)
    if (!is.null(f)) fitted[[fam]] <- f
  }
  if (length(fitted) == 0) stop("fit_quantile_dist: no candidate family could be fit.")

  resid <- vapply(fitted, function(f) f$residual, numeric(1))
  ranking <- data.frame(
    family   = names(resid),
    residual = unname(resid),
    params   = vapply(names(resid), function(k)
                 paste(sprintf("%s=%.4g", names(fitted[[k]]$par), fitted[[k]]$par), collapse = ", "),
                 character(1)),
    stringsAsFactors = FALSE)
  ranking <- ranking[order(ranking$residual), ]
  best <- ranking$family[1]
  bf   <- fitted[[best]]
  list(best = best, params = bf$par, residual = bf$residual,
       q = bf$Q, d = bf$d, r = function(u) bf$Q(u), ranking = ranking)
}
