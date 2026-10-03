#' Simulate one arm of a trial on a daily grid
#'
#' Event times follow a geometric distribution with daily hazard
#' \eqn{1 - \exp(-\lambda_0 m)}, where \eqn{m} is a hazard multiplier, either
#' constant or a function of the day. Independent geometric dropout is chosen
#' so that its probability by \code{tau} is \code{dropout} in the absence of
#' events; follow-up ends administratively at \code{tau}.
#'
#' @param n Number of participants.
#' @param hazard_mult Hazard multiplier relative to the reference rate: a
#'   number, or a vectorised function of the day returning multipliers.
#' @param tau Administrative follow-up in days.
#' @param base_rate Reference daily rate \eqn{\lambda_0} (default: median one year).
#' @param dropout Probability of dropout by \code{tau} without events.
#' @return A data frame with integer \code{time} and \code{status}.
#' @examples
#' set.seed(4)
#' head(sim_arm(100, 0.7))
#' @export
sim_arm <- function(n, hazard_mult = 1, tau = 365, base_rate = log(2) / 365,
                    dropout = 0.15) {
  pc <- 1 - (1 - dropout)^(1 / tau)
  if (is.function(hazard_mult)) {
    h <- 1 - exp(-base_rate * hazard_mult(seq_len(tau)))
    u <- matrix(stats::runif(n * tau), nrow = tau)
    ev <- u < h
    first <- apply(ev, 2, function(z) { w <- which(z); if (length(w)) w[1] else tau + 1L })
    T <- first
  } else {
    p <- 1 - exp(-base_rate * hazard_mult)
    T <- stats::rgeom(n, p) + 1L
  }
  C <- if (pc > 0) stats::rgeom(n, pc) + 1L else rep(.Machine$integer.max, n)
  event <- T <= C & T <= tau
  data.frame(time = as.integer(ifelse(event, T, pmin(C, tau))),
             status = as.integer(event))
}

#' True RMST of a simulated arm
#'
#' @inheritParams sim_arm
#' @param convention RMST convention (see \code{\link{rmst_moments}}).
#' @return The RMST implied by the generating daily hazards.
#' @export
true_rmst <- function(hazard_mult = 1, tau = 365, base_rate = log(2) / 365,
                      convention = c("left", "endpoint")) {
  convention <- match.arg(convention)
  m <- if (is.function(hazard_mult)) hazard_mult(seq_len(tau)) else rep(hazard_mult, tau)
  S <- cumprod(exp(-base_rate * m))
  if (convention == "left") 1 + sum(S[-tau]) else sum(S)
}

#' Design-stage variance for the tuning constant
#'
#' Estimates, by simulation under planning assumptions, the expected variance
#' \eqn{V_3 + V_2} of the difference between the current and historical
#' arm-level RMST posterior means, for use in \code{\link{k_from_delta}}.
#'
#' @param n_current,n_historical Arm sample sizes.
#' @param hazard_mult Planning hazard multiplier for the arm.
#' @param tau Horizon in days.
#' @param nsim Number of simulated trial pairs.
#' @param prior Daily beta prior.
#' @param ... Further arguments to \code{\link{sim_arm}}.
#' @return Mean of \eqn{V_3 + V_2} over simulated pairs.
#' @examples
#' set.seed(5)
#' design_variance(150, 60, nsim = 200)
#' @export
design_variance <- function(n_current, n_historical, hazard_mult = 1, tau = 365,
                            nsim = 1000, prior = c(0.001, 1), ...) {
  v <- vapply(seq_len(nsim), function(i) {
    d3 <- sim_arm(n_current, hazard_mult, tau, ...)
    d2 <- sim_arm(n_historical, hazard_mult, tau, ...)
    c3 <- bkm_counts(d3$time, d3$status, tau)
    c2 <- bkm_counts(d2$time, d2$status, tau)
    rmst_moments(prior[1] + c3$D, prior[2] + c3$Y)$var +
      rmst_moments(prior[1] + c2$D, prior[2] + c2$Y)$var
  }, numeric(1))
  mean(v)
}
