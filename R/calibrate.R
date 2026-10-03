#' Calibrate the posterior-probability threshold for type I error control
#'
#' Given posterior probabilities \eqn{P(\Delta > \delta_0 \mid D)} simulated
#' under each scenario of a prespecified null set \eqn{\Theta_0}, returns the
#' smallest threshold \eqn{\gamma} such that the rejection rate
#' \eqn{\Pr(P > \gamma)} is at most \code{alpha} in every scenario.
#'
#' With \code{criterion = "estimate"} the empirical rejection rate must not
#' exceed \code{alpha}. With \code{criterion = "upper"} a one-sided
#' Clopper-Pearson upper confidence bound at level \code{conf}, Bonferroni
#' adjusted over the scenarios, must not exceed \code{alpha}. Because the
#' threshold is selected using the same simulated data, its operating
#' characteristics should be confirmed on an independent simulation stream.
#'
#' @param pp A matrix (replicates by null scenarios) or a list of numeric
#'   vectors of posterior probabilities, one per null scenario.
#' @param alpha Target maximum type I error rate.
#' @param criterion \code{"estimate"} or \code{"upper"}.
#' @param conf Confidence level for \code{criterion = "upper"}.
#' @return A list with the calibrated threshold \code{gamma}, the scenario
#'   attaining the maximum, and the per-scenario thresholds and rejection rates.
#' @examples
#' set.seed(3)
#' pp <- cbind(runif(5000), pbeta(runif(5000), 1, 0.8))
#' calibrate_gamma(pp, alpha = 0.025)$gamma
#' @export
calibrate_gamma <- function(pp, alpha = 0.025, criterion = c("estimate", "upper"),
                            conf = 0.95) {
  criterion <- match.arg(criterion)
  if (is.matrix(pp)) pp <- lapply(seq_len(ncol(pp)), function(j) pp[, j])
  S <- length(pp)
  per <- vapply(pp, function(x) {
    x <- x[!is.na(x)]
    R <- length(x)
    cmax <- if (criterion == "estimate") floor(alpha * R) else {
      lev <- 1 - (1 - conf) / S
      cs <- 0:ceiling(alpha * R)
      ub <- ifelse(cs == R, 1, stats::qbeta(lev, cs + 1, R - cs))
      ok <- cs[ub <= alpha]
      if (length(ok)) max(ok) else -1
    }
    if (cmax < 0) return(1)
    xs <- sort(x, decreasing = TRUE)
    # smallest gamma with sum(x > gamma) <= cmax: the (cmax+1)-th largest value
    g <- xs[cmax + 1]
    # ties at g are excluded by the strict inequality; ensure count is respected
    min(g, 1)
  }, numeric(1))
  gamma <- max(per)
  rates <- vapply(pp, function(x) mean(x[!is.na(x)] > gamma), numeric(1))
  list(gamma = gamma, worst = which.max(per), gamma_by_scenario = per,
       rejection = rates, criterion = criterion, alpha = alpha)
}

#' Rejection rate at a threshold
#'
#' @param pp Posterior probabilities (vector, matrix or list as in
#'   \code{\link{calibrate_gamma}}).
#' @param gamma Threshold.
#' @return Rejection rate(s) \eqn{\Pr(P > \gamma)} with Monte Carlo standard errors.
#' @export
rejection_rate <- function(pp, gamma) {
  if (is.matrix(pp)) pp <- lapply(seq_len(ncol(pp)), function(j) pp[, j])
  if (!is.list(pp)) pp <- list(pp)
  t(vapply(pp, function(x) {
    x <- x[!is.na(x)]
    r <- mean(x > gamma)
    c(rate = r, mcse = sqrt(r * (1 - r) / length(x)))
  }, numeric(2)))
}
