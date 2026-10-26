#' Tuning constant from a tolerance in days
#'
#' Returns \eqn{k_g = \log(2)\, \bar V_g / \delta^2}, so that the borrowing
#' weight \eqn{w_g = \exp(-k_g Z_g^2)} equals 1/2 when the historical and
#' current arm-level RMSTs differ by \eqn{\delta} days and
#' \eqn{V_{1g} + V_{0g} = \bar V_g}.
#'
#' @param delta Tolerance \eqn{\delta} in days (the discrepancy receiving half weight).
#' @param v_design Design stage value \eqn{\bar V_g} of \eqn{V_{1g} + V_{0g}},
#'   for example from \code{\link{design_variance}}.
#' @return The tuning constant \eqn{k_g}.
#' @examples
#' k_from_delta(delta = 20, v_design = 280)
#' @export
k_from_delta <- function(delta, v_design) {
  if (any(delta <= 0) || any(v_design <= 0)) stop("'delta' and 'v_design' must be positive")
  log(2) * v_design / delta^2
}

#' RMST-guided borrowing weight for one arm
#'
#' Computes the standardized discrepancy
#' \eqn{Z_g = (\hat\mu_{1g} - \hat\mu_{0g}) / (V_{1g} + V_{0g})^{1/2}} from the
#' current-only (\eqn{s = 1}) and historical-only (\eqn{s = 0}) Bayesian
#' Kaplan-Meier posteriors of the arm's RMST under the reference prior, and
#' returns \eqn{w_g = \exp(-k_g Z_g^2)}.
#'
#' @param cur,hist \code{bkm_counts} objects for the current and historical arm.
#' @param k Non-negative tuning constant \eqn{k_g} (\code{Inf} gives no borrowing).
#' @param prior Shape parameters \eqn{(a_0, b_0)} of the reference prior.
#' @param convention RMST convention passed to \code{\link{rmst_moments}}.
#' @return A list with the weight \code{w}, \code{z2} (\eqn{Z_g^2}), and the
#'   current and historical posterior moments.
#' @export
bkm_weight <- function(cur, hist, k, prior = c(0.001, 1), convention = "left") {
  mc <- rmst_moments(prior[1] + cur$X, prior[2] + cur$Y, convention)
  mh <- rmst_moments(prior[1] + hist$X, prior[2] + hist$Y, convention)
  z2 <- (mc$mean - mh$mean)^2 / (mc$var + mh$var)
  w <- if (is.infinite(k)) 0 else exp(-k * z2)
  list(w = w, z2 = z2, current = mc, historical = mh)
}

#' Bayesian Kaplan-Meier analysis with dynamic historical borrowing
#'
#' Fits the Bayesian Kaplan-Meier model to a two-arm current trial, borrowing
#' historical data for the control arm (configuration C), the treatment arm
#' (configuration T), or both (configuration B), with the RMST-guided weight
#' \eqn{w_g = \exp(-k_g Z_g^2)}.
#' With \code{method = "hazard"} (RG-H), the historical daily counts of arm
#' \eqn{g} are multiplied by \eqn{w_g}, so each current daily hazard has the
#' power-prior posterior
#' \eqn{\mathrm{Beta}(a_0 + w_g X_{0gt} + X_{1gt},\ b_0 + w_g Y_{0gt} + Y_{1gt})}.
#' With \code{method = "rmst"} (RG-R), the current and historical arm-level
#' posterior RMSTs are combined as a normal power prior with weight \eqn{w_g}.
#'
#' The posterior probability that the treatment is beneficial,
#' \eqn{\eta = \Pr(\Delta > \delta_0 \mid \mathcal{D})} with
#' \eqn{\Delta = \mu_{11} - \mu_{10}}, is by default computed by the normal
#' approximation \eqn{\Phi\{(E(\Delta \mid \mathcal{D}) - \delta_0) /
#' \mathrm{Var}(\Delta \mid \mathcal{D})^{1/2}\}} with the exact posterior
#' moments, and the 95\% credible interval is
#' \eqn{E(\Delta \mid \mathcal{D}) \pm 1.96\, \mathrm{Var}(\Delta \mid \mathcal{D})^{1/2}}.
#'
#' @param current List with elements \code{control} and \code{treatment}, each
#'   a data frame with \code{time} and \code{status} or a \code{bkm_counts} object.
#' @param historical List with elements \code{control} and/or \code{treatment};
#'   a missing or \code{NULL} element means no historical data for that arm.
#' @param tau RMST time point \eqn{\tau} in days, common to both trials (see
#'   \code{\link{bkm_counts}}).
#' @param method \code{"hazard"} (RG-H, default), \code{"rmst"} (RG-R),
#'   \code{"fixed"} (fixed power prior with weight \code{w} on the daily counts),
#'   or \code{"none"} (no borrowing).
#' @param k Tuning constant(s) \eqn{k_g} for \code{"hazard"} and \code{"rmst"}:
#'   a single value or a named vector \code{c(control = , treatment = )}.
#' @param w Fixed weight(s) for \code{"fixed"} (scalar or named vector).
#' @param prior Reference prior shape parameters \eqn{(a_0, b_0)}.
#' @param delta0 Margin \eqn{\delta_0} in \eqn{\Pr(\Delta > \delta_0 \mid \mathcal{D})} (default 0).
#' @param posterior \code{"normal"} (normal approximation with exact moments,
#'   default) or \code{"draws"} (posterior sampling; RG-H, fixed power prior
#'   and no borrowing only).
#' @param ndraw Number of draws when \code{posterior = "draws"}.
#' @param convention RMST convention (see \code{\link{rmst_moments}}).
#' @return An object of class \code{"bkmborrow"} with the weights, arm-level
#'   posterior RMST moments, the posterior summary of \eqn{\Delta}
#'   (\code{diff}) and \code{eta}.
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
                       posterior = c("normal", "draws"), ndraw = 4000,
                       convention = c("left", "endpoint")) {
  method <- match.arg(method)
  posterior <- match.arg(posterior)
  convention <- match.arg(convention)
  arms <- c("control", "treatment")
  if (!all(arms %in% names(current))) stop("'current' needs 'control' and 'treatment'")
  cur <- lapply(current[arms], as_counts, tau = tau)
  his <- lapply(arms, function(a) as_counts(historical[[a]], tau = tau))
  names(his) <- arms
  getpar <- function(x, a) if (length(x) == 1L) unname(x) else unname(x[[a]])

  res <- lapply(arms, function(a) {
    h <- his[[a]]
    c1 <- cur[[a]]
    m1 <- rmst_moments(prior[1] + c1$X, prior[2] + c1$Y, convention)
    if (is.null(h) || method == "none") {
      return(list(w = 0, z2 = NA_real_, mean = m1$mean, var = m1$var,
                  A = prior[1] + c1$X, B = prior[2] + c1$Y))
    }
    m0 <- rmst_moments(prior[1] + h$X, prior[2] + h$Y, convention)
    z2 <- (m1$mean - m0$mean)^2 / (m1$var + m0$var)
    ww <- if (method == "fixed") getpar(w, a) else {
      kk <- getpar(k, a)
      if (is.infinite(kk)) 0 else exp(-kk * z2)
    }
    if (method == "rmst") {
      prec <- 1 / m1$var + ww / m0$var
      return(list(w = ww, z2 = z2, mean = (m1$mean / m1$var + ww * m0$mean / m0$var) / prec,
                  var = 1 / prec, A = NULL, B = NULL))
    }
    po <- arm_posterior(c1, h, ww, prior)
    mm <- rmst_moments(po$A, po$B, convention)
    list(w = ww, z2 = z2, mean = mm$mean, var = mm$var, A = po$A, B = po$B)
  })
  names(res) <- arms

  dmean <- res$treatment$mean - res$control$mean
  dvar <- res$treatment$var + res$control$var
  draws <- NULL
  if (posterior == "normal") {
    sd <- sqrt(dvar)
    eta <- stats::pnorm((dmean - delta0) / sd)
    ci <- dmean + c(-1, 1) * stats::qnorm(0.975) * sd
  } else {
    if (method == "rmst") stop("posterior = 'draws' is not available for RG-R")
    dr <- function(r) rmst_draws(r$A, r$B, ndraw, convention)
    draws <- dr(res$treatment) - dr(res$control)
    eta <- mean(draws > delta0)
    ci <- unname(stats::quantile(draws, c(0.025, 0.975)))
  }
  structure(list(
    method = method, tau = tau, prior = prior, convention = convention,
    weights = c(control = res$control$w, treatment = res$treatment$w),
    z2 = c(control = res$control$z2, treatment = res$treatment$z2),
    rmst = rbind(control = c(mean = res$control$mean, var = res$control$var),
                 treatment = c(mean = res$treatment$mean, var = res$treatment$var)),
    diff = c(mean = dmean, sd = sqrt(dvar), lower = ci[1], upper = ci[2]),
    delta0 = delta0, eta = eta, posterior = posterior, draws = draws,
    historical_arms = arms[!vapply(his, is.null, logical(1))]
  ), class = "bkmborrow")
}

#' @export
print.bkmborrow <- function(x, digits = 3, ...) {
  lab <- c(hazard = "RG-H (hazard level)", rmst = "RG-R (RMST level)",
           fixed = "fixed power prior", none = "no borrowing")
  cat("Bayesian Kaplan-Meier analysis:", lab[[x$method]], "\n")
  cat("RMST time point tau:", x$tau, "days; historical arms:",
      if (length(x$historical_arms)) paste(x$historical_arms, collapse = ", ") else "none", "\n")
  cat("Borrowing weights: w_0 (control) =", format(x$weights["control"], digits = digits),
      ", w_1 (treatment) =", format(x$weights["treatment"], digits = digits), "\n\n")
  print(round(cbind(x$rmst, sd = sqrt(x$rmst[, "var"])), digits))
  cat("\nRMST difference Delta (treatment - control):", format(x$diff["mean"], digits = digits + 1),
      " 95% CrI (", format(x$diff["lower"], digits = digits + 1), ",",
      format(x$diff["upper"], digits = digits + 1), ")\n")
  cat("eta = Pr(Delta >", x$delta0, "| data) =", format(x$eta, digits = digits + 1),
      paste0("[", x$posterior, "]"), "\n")
  invisible(x)
}

#' Success criterion
#'
#' Applies the success criterion \eqn{\eta > \gamma_0}: efficacy is concluded
#' if the posterior probability that the treatment is beneficial exceeds the
#' threshold.
#'
#' @param fit A \code{bkmborrow} fit.
#' @param gamma0 Threshold \eqn{\gamma_0}, e.g. from \code{\link{gamma_nb}}.
#' @return Logical: \code{TRUE} if \eqn{\eta > \gamma_0}.
#' @export
bkm_success <- function(fit, gamma0) unname(fit$eta > gamma0)
