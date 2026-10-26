#' Simulate one arm of a trial on a daily grid
#'
#' Event times follow a geometric distribution with daily hazard
#' \eqn{1 - \exp(-\lambda_0 m)}, where \eqn{m} is a hazard multiplier, either
#' constant or a function of the day. Independent geometric dropout is chosen
#' so that its probability by the end of follow-up is \code{dropout} in the
#' absence of events, and follow-up ends administratively at day
#' \code{follow_up}. The follow-up is distinct from the RMST time point
#' \eqn{\tau} used in the analysis.
#'
#' @param n Number of participants.
#' @param hazard_mult Hazard multiplier relative to the reference rate: a
#'   number, or a vectorised function of the day returning multipliers.
#' @param follow_up Administrative follow-up in days.
#' @param base_rate Reference daily rate \eqn{\lambda_0} (default: median one year).
#' @param dropout Probability of dropout by \code{follow_up} without events.
#' @return A data frame with integer \code{time} and \code{status}.
#' @examples
#' set.seed(4)
#' head(sim_arm(100, 0.7))
#' @export
sim_arm <- function(n, hazard_mult = 1, follow_up = 365, base_rate = log(2) / 365,
                    dropout = 0.15) {
  follow_up <- as.integer(follow_up)
  pc <- 1 - (1 - dropout)^(1 / follow_up)
  if (is.function(hazard_mult)) {
    h <- 1 - exp(-base_rate * hazard_mult(seq_len(follow_up)))
    u <- matrix(stats::runif(n * follow_up), nrow = follow_up)
    ev <- u < h
    T <- apply(ev, 2, function(z) { w <- which(z); if (length(w)) w[1] else follow_up + 1L })
  } else {
    p <- 1 - exp(-base_rate * hazard_mult)
    T <- stats::rgeom(n, p) + 1L
  }
  C <- if (pc > 0) stats::rgeom(n, pc) + 1L else rep(.Machine$integer.max, n)
  event <- T <= C & T <= follow_up
  data.frame(time = as.integer(ifelse(event, T, pmin(C, follow_up))),
             status = as.integer(event))
}

#' True RMST of a simulated arm
#'
#' @param hazard_mult Hazard multiplier (number or function of the day), as in
#'   \code{\link{sim_arm}}.
#' @param tau RMST time point in days.
#' @param base_rate Reference daily rate.
#' @param convention RMST convention (see \code{\link{rmst_moments}}).
#' @return The RMST to \code{tau} implied by the generating daily hazards.
#' @export
true_rmst <- function(hazard_mult = 1, tau = 365, base_rate = log(2) / 365,
                      convention = c("left", "endpoint")) {
  convention <- match.arg(convention)
  m <- if (is.function(hazard_mult)) hazard_mult(seq_len(tau)) else rep(hazard_mult, tau)
  S <- cumprod(exp(-base_rate * m))
  if (convention == "left") 1 + sum(S[-tau]) else sum(S)
}

#' Design stage variance for the tuning constant
#'
#' Estimates \eqn{\bar V_g}, the design stage value of \eqn{V_{1g} + V_{0g}},
#' the sum of the current-only and historical-only posterior variances of the
#' arm-level RMST, by simulation under planning assumptions. It is used in
#' \code{\link{k_from_delta}}.
#'
#' @param n_current,n_historical Arm sample sizes in the current and
#'   historical trials.
#' @param hazard_mult Planning hazard multiplier for the arm.
#' @param tau RMST time point in days.
#' @param follow_up Follow-up in days (default \code{tau}).
#' @param nsim Number of simulated trial pairs.
#' @param prior Daily beta prior \eqn{(a_0, b_0)}.
#' @param ... Further arguments to \code{\link{sim_arm}} (e.g. \code{dropout}).
#' @return Mean of \eqn{V_{1g} + V_{0g}} over simulated pairs.
#' @examples
#' set.seed(5)
#' design_variance(150, 60, nsim = 200)
#' @export
design_variance <- function(n_current, n_historical, hazard_mult = 1, tau = 365,
                            follow_up = tau, nsim = 1000, prior = c(0.001, 1), ...) {
  v <- vapply(seq_len(nsim), function(i) {
    d1 <- sim_arm(n_current, hazard_mult, follow_up, ...)
    d0 <- sim_arm(n_historical, hazard_mult, follow_up, ...)
    c1 <- bkm_counts(d1$time, d1$status, tau)
    c0 <- bkm_counts(d0$time, d0$status, tau)
    rmst_moments(prior[1] + c1$X, prior[2] + c1$Y)$var +
      rmst_moments(prior[1] + c0$X, prior[2] + c0$Y)$var
  }, numeric(1))
  mean(v)
}
