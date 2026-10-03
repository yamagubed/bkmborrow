#' Tuning constant from a clinically elicited tolerance
#'
#' Returns \eqn{k = \log(2)\, V / \delta^2}, so that the borrowing weight
#' \eqn{w = \exp(-k Z^2)} equals 1/2 when the historical and current
#' arm-level RMSTs differ by \eqn{\delta} days and the variance of their
#' difference equals its design value \eqn{V}.
#'
#' @param delta Tolerance in days (the discrepancy receiving half weight).
#' @param v_design Design-stage variance of the difference between the
#'   current and historical arm-level RMST estimates, for example from
#'   \code{\link{design_variance}}.
#' @return The tuning constant \eqn{k}.
#' @examples
#' k_from_delta(delta = 20, v_design = 280)
#' @export
k_from_delta <- function(delta, v_design) {
  if (any(delta <= 0) || any(v_design <= 0)) stop("'delta' and 'v_design' must be positive")
  log(2) * v_design / delta^2
}

#' RMST-guided borrowing weight for one arm
#'
#' Computes \eqn{Z^2 = (\hat\mu_3 - \hat\mu_2)^2 / (V_3 + V_2)} from the
#' current-only and historical-only Bayesian Kaplan-Meier posteriors of the
#' arm's RMST, and returns \eqn{w = \exp(-k Z^2)}.
#'
#' @param cur,hist \code{bkm_counts} objects for the current and historical arm.
#' @param k Non-negative tuning constant (\code{Inf} gives no borrowing).
#' @param prior Shape parameters \eqn{(a_0, b_0)} of the daily beta prior.
#' @param convention RMST convention passed to \code{\link{rmst_moments}}.
#' @return A list with the weight \code{w}, \code{z2}, and the two sets of moments.
#' @export
bkm_weight <- function(cur, hist, k, prior = c(0.001, 1), convention = "left") {
  mc <- rmst_moments(prior[1] + cur$D, prior[2] + cur$Y, convention)
  mh <- rmst_moments(prior[1] + hist$D, prior[2] + hist$Y, convention)
  z2 <- (mc$mean - mh$mean)^2 / (mc$var + mh$var)
  w <- if (is.infinite(k)) 0 else exp(-k * z2)
  list(w = w, z2 = z2, current = mc, historical = mh)
}

#' Bayesian Kaplan-Meier analysis with dynamic historical borrowing
#'
#' Fits the daily-beta Bayesian Kaplan-Meier model to a two-arm current trial,
#' borrowing historical data for the control arm, the treatment arm, or both,
#' with the RMST-guided weight \eqn{w_g = \exp(-k_g Z_g^2)}.
#' With \code{method = "hazard"} (RG-H) the historical daily counts of arm
#' \eqn{g} are multiplied by \eqn{w_g}, so each daily hazard has posterior
#' \eqn{\mathrm{Beta}(a_0 + w_g D_{0gd} + D_{1gd},\ b_0 + w_g Y_{0gd} + Y_{1gd})}.
#' With \code{method = "rmst"} (RG-R) the current and historical arm-level
#' posterior RMSTs are combined as a normal power prior with weight \eqn{w_g}.
#'
#' @param current List with elements \code{control} and \code{treatment}, each
#'   a data frame with \code{time} and \code{status} or a \code{bkm_counts} object.
#' @param historical List with elements \code{control} and/or \code{treatment};
#'   a missing or \code{NULL} element means no historical data for that arm.
#' @param tau Restriction horizon in days.
#' @param method \code{"hazard"} (RG-H, default), \code{"rmst"} (RG-R),
#'   \code{"fixed"} (fixed power prior with weight \code{w} on the daily counts),
#'   or \code{"none"}.
#' @param k Tuning constant(s) for \code{"hazard"} and \code{"rmst"}: a single
#'   value or a named vector \code{c(control = , treatment = )}.
#' @param w Fixed weight(s) for \code{"fixed"} (scalar or named vector).
#' @param prior Daily beta prior shape parameters \eqn{(a_0, b_0)}.
#' @param delta0 Margin for the posterior probability \eqn{P(\Delta > \delta_0)}.
#' @param prob \code{"normal"} (normal approximation using exact moments,
#'   default) or \code{"draws"} (posterior simulation; hazard-level methods only).
#' @param ndraw Number of draws when \code{prob = "draws"}.
#' @param convention RMST convention (see \code{\link{rmst_moments}}).
#' @return An object of class \code{"bkmborrow"}.
#' @examples
#' set.seed(2)
#' cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.7))
#' his <- list(control = sim_arm(60, 1), treatment = sim_arm(60, 0.7))
#' fit <- bkm_borrow(cur, his, tau = 365, k = 0.5)
#' fit
#' @export
bkm_borrow <- function(current, historical = list(), tau,
                       method = c("hazard", "rmst", "fixed", "none"),
                       k = 1, w = 0.5,
                       prior = c(0.001, 1), delta0 = 0,
                       prob = c("normal", "draws"), ndraw = 4000,
                       convention = c("left", "endpoint")) {
  method <- match.arg(method)
  prob <- match.arg(prob)
  convention <- match.arg(convention)
  arms <- c("control", "treatment")
  if (!all(arms %in% names(current))) stop("'current' needs 'control' and 'treatment'")
  cur <- lapply(current[arms], as_counts, tau = tau)
  his <- lapply(arms, function(a) as_counts(historical[[a]], tau = tau))
  names(his) <- arms
  getpar <- function(x, a) if (length(x) == 1L) unname(x) else unname(x[[a]])

  res <- lapply(arms, function(a) {
    h <- his[[a]]
    c0 <- cur[[a]]
    m0 <- rmst_moments(prior[1] + c0$D, prior[2] + c0$Y, convention)
    if (is.null(h) || method == "none") {
      return(list(w = 0, z2 = NA_real_, mean = m0$mean, var = m0$var,
                  comp = list(list(p = 1, A = prior[1] + c0$D, B = prior[2] + c0$Y))))
    }
    mh <- rmst_moments(prior[1] + h$D, prior[2] + h$Y, convention)
    z2 <- (m0$mean - mh$mean)^2 / (m0$var + mh$var)
    if (method == "rmst") {
      kk <- getpar(k, a)
      ww <- if (is.infinite(kk)) 0 else exp(-kk * z2)
      prec <- 1 / m0$var + ww / mh$var
      return(list(w = ww, z2 = z2, mean = (m0$mean / m0$var + ww * mh$mean / mh$var) / prec,
                  var = 1 / prec, comp = NULL))
    }
    ww <- if (method == "fixed") getpar(w, a) else {
      kk <- getpar(k, a)
      if (is.infinite(kk)) 0 else exp(-kk * z2)
    }
    po <- arm_posterior(c0, h, ww, prior)
    mm <- rmst_moments(po$A, po$B, convention)
    list(w = ww, z2 = z2, mean = mm$mean, var = mm$var,
         comp = list(list(p = 1, A = po$A, B = po$B)))
  })
  names(res) <- arms

  dmean <- res$treatment$mean - res$control$mean
  dvar <- res$treatment$var + res$control$var
  draws <- NULL
  if (prob == "normal") {
    sd <- sqrt(dvar)
    pp <- stats::pnorm((dmean - delta0) / sd)
    ci <- dmean + c(-1, 1) * stats::qnorm(0.975) * sd
  } else {
    if (method == "rmst") stop("prob = 'draws' is available for hazard-level methods only")
    draw_arm <- function(r) {
      if (length(r$comp) == 1L) return(rmst_draws(r$comp[[1]]$A, r$comp[[1]]$B, ndraw, convention))
      n1 <- stats::rbinom(1, ndraw, r$comp[[1]]$p)
      c(if (n1 > 0) rmst_draws(r$comp[[1]]$A, r$comp[[1]]$B, n1, convention),
        if (n1 < ndraw) rmst_draws(r$comp[[2]]$A, r$comp[[2]]$B, ndraw - n1, convention))[sample.int(ndraw)]
    }
    draws <- draw_arm(res$treatment) - draw_arm(res$control)
    pp <- mean(draws > delta0)
    ci <- unname(stats::quantile(draws, c(0.025, 0.975)))
  }
  structure(list(
    method = method, tau = tau, prior = prior, convention = convention,
    weights = c(control = res$control$w, treatment = res$treatment$w),
    z2 = c(control = res$control$z2, treatment = res$treatment$z2),
    rmst = rbind(control = c(mean = res$control$mean, var = res$control$var),
                 treatment = c(mean = res$treatment$mean, var = res$treatment$var)),
    diff = c(mean = dmean, sd = sqrt(dvar), lower = ci[1], upper = ci[2]),
    delta0 = delta0, prob = pp, prob_method = prob, draws = draws,
    historical_arms = arms[!vapply(his, is.null, logical(1))]
  ), class = "bkmborrow")
}

#' @export
print.bkmborrow <- function(x, digits = 3, ...) {
  cat("Bayesian Kaplan-Meier analysis, method:", x$method, "\n")
  cat("Horizon:", x$tau, "days; historical arms:",
      if (length(x$historical_arms)) paste(x$historical_arms, collapse = ", ") else "none", "\n")
  cat("Borrowing weights: control =", format(x$weights["control"], digits = digits),
      ", treatment =", format(x$weights["treatment"], digits = digits), "\n\n")
  print(round(cbind(x$rmst, sd = sqrt(x$rmst[, "var"])), digits))
  cat("\nRMST difference (treatment - control):", format(x$diff["mean"], digits = digits + 1),
      " 95% CrI (", format(x$diff["lower"], digits = digits + 1), ",",
      format(x$diff["upper"], digits = digits + 1), ")\n")
  cat("P(Delta >", x$delta0, "| data) =", format(x$prob, digits = digits + 1),
      paste0("[", x$prob_method, "]"), "\n")
  invisible(x)
}

#' Decision with a calibrated threshold
#'
#' @param fit A \code{bkmborrow} fit.
#' @param gamma Posterior probability threshold, e.g. from \code{\link{calibrate_gamma}}.
#' @return Logical: \code{TRUE} if \eqn{P(\Delta > \delta_0 \mid D) > \gamma}.
#' @export
bkm_decision <- function(fit, gamma) unname(fit$prob > gamma)
