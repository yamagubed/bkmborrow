#' Maximum attainable power under bounded harmful drift
#'
#' Upper bound on the power at no drift of any test that uses a current trial
#' and one historical arm \eqn{g}, subject to a type I error of at most
#' \code{alpha} for every harmful drift whose bias in the arm-level RMST is
#' between 0 and \code{bias_max} days. Under a normal approximation with known
#' variances, \eqn{(\hat\Delta_{\mathrm{NB}}, h_g)} is bivariate normal with mean
#' \eqn{(\Delta, \beta_g)}, where \eqn{h_g = \pm(\hat\mu_{0g} - \hat\mu_{1g})} is
#' oriented so that harmful drift has \eqn{\beta_g > 0}. With
#' \eqn{\sigma^2 = V_{10} + V_{11}}, \eqn{\beta^* = V_{1g}\Delta_1/\sigma^2} and
#' \eqn{\tilde\beta = \min(\beta_{\max}, \beta^*)}, the bound is
#' \eqn{\Phi\{(m^\top\Sigma^{-1}m)^{1/2} - z_{1-\alpha}\}} with
#' \eqn{m = (\Delta_1, -\tilde\beta)^\top}. It is attained by fixed-weight
#' borrowing, \eqn{\hat\Delta_{\mathrm{NB}} + \lambda^* h_g}, and equals the
#' power of no borrowing at level \code{alpha} when
#' \eqn{\beta_{\max} \ge \beta^*}.
#'
#' @param delta1 RMST difference \eqn{\Delta_1} under the alternative (days).
#' @param v_current Named vector \code{c(control = , treatment = )} of the
#'   current-trial posterior variances \eqn{V_{10}, V_{11}} of the arm-level RMSTs.
#' @param v_hist Posterior variance \eqn{V_{0g}} of the historical arm-level RMST.
#' @param arm \code{"control"} or \code{"treatment"}: the arm with historical data.
#' @param bias_max Largest harmful bias \eqn{\beta_{\max} \ge 0} (days) in the
#'   historical arm-level RMST over the calibration range, e.g.
#'   \code{true_rmst(1) - true_rmst(1.1)} for a historical control hazard 10\%
#'   higher than the current one.
#' @param alpha Type I error bound (e.g. \eqn{\alpha_{\mathrm{cal}}}).
#' @return A list with the bound \code{power}, the power of no borrowing at
#'   level \code{alpha} (\code{power_nb}), \code{beta_star}, the optimal
#'   \code{lambda} and the equivalent fixed power-prior weight \code{w}.
#' @examples
#' power_bound(25.84, c(control = 52.2, treatment = 44.5), v_hist = 52.2,
#'             arm = "control", bias_max = true_rmst(1) - true_rmst(1.1))
#' @export
power_bound <- function(delta1, v_current, v_hist, arm = c("control", "treatment"),
                        bias_max, alpha = 0.04) {
  arm <- match.arg(arm)
  if (bias_max < 0) stop("'bias_max' must be non-negative")
  V1 <- unname(v_current[[arm]]); s2 <- sum(v_current[c("control", "treatment")])
  S <- matrix(c(s2, -V1, -V1, V1 + v_hist), 2)
  bstar <- V1 * delta1 / s2; bt <- min(bias_max, bstar)
  m <- c(delta1, -bt); l <- solve(S, m)
  lam <- if (abs(l[2]) < 1e-12) 0 else l[2] / l[1]
  list(power = stats::pnorm(sqrt(sum(m * l)) - stats::qnorm(1 - alpha)),
       power_nb = stats::pnorm(delta1 / sqrt(s2) - stats::qnorm(1 - alpha)),
       beta_star = bstar, lambda = lam, w = lam * v_hist / ((1 - lam) * V1))
}
