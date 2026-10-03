#' Exact posterior moments of the RMST
#'
#' For independent daily hazards \eqn{p_d \sim \mathrm{Beta}(A_d, B_d)},
#' computes the exact mean and variance of the restricted mean survival time.
#' With \code{convention = "left"} the RMST is \eqn{\sum_{t=0}^{K-1} S(t)},
#' the area under the step survival function for integer event times; with
#' \code{"endpoint"} it is \eqn{\sum_{t=1}^{K} S(t)}.
#'
#' The computation uses \eqn{E\{S(t)S(u)\} = E\{S(t)^2\}\,E\{S(u)\}/E\{S(t)\}}
#' for \eqn{t \le u}, so it needs work proportional to \eqn{K}.
#'
#' @param A,B Vectors of beta shape parameters, one per day.
#' @param convention \code{"left"} (default) or \code{"endpoint"}.
#' @return A list with \code{mean} and \code{var}.
#' @examples
#' rmst_moments(A = rep(1, 30), B = rep(100, 30))
#' @export
rmst_moments <- function(A, B, convention = c("left", "endpoint")) {
  convention <- match.arg(convention)
  K <- length(A)
  if (length(B) != K) stop("'A' and 'B' differ in length")
  if (any(A <= 0) || any(B <= 0)) stop("beta parameters must be positive")
  q1 <- B / (A + B)
  q2 <- q1 * (B + 1) / (A + B + 1)
  if (convention == "left") {
    m <- c(1, cumprod(q1)[-K])
    m2 <- c(1, cumprod(q2)[-K])
    r <- c(1, cumprod(q2 / q1)[-K])
  } else {
    m <- cumprod(q1)
    m2 <- cumprod(q2)
    r <- cumprod(q2 / q1)
  }
  mu <- sum(m)
  tail <- rev(cumsum(rev(m))) - m
  es2 <- sum(m2) + 2 * sum(r * tail)
  list(mean = mu, var = max(es2 - mu^2, 0))
}

#' Posterior draws of the RMST
#'
#' @inheritParams rmst_moments
#' @param ndraw Number of posterior draws.
#' @return Numeric vector of RMST draws.
#' @examples
#' set.seed(1)
#' summary(rmst_draws(A = rep(1, 30), B = rep(100, 30), ndraw = 1000))
#' @export
rmst_draws <- function(A, B, ndraw = 4000, convention = c("left", "endpoint")) {
  convention <- match.arg(convention)
  K <- length(A)
  S <- rep(1, ndraw)
  total <- if (convention == "left") rep(1, ndraw) else rep(0, ndraw)
  last <- if (convention == "left") K - 1L else K
  for (d in seq_len(last)) {
    S <- S * (1 - stats::rbeta(ndraw, A[d], B[d]))
    total <- total + S
  }
  total
}

arm_posterior <- function(cur, hist, w, prior) {
  A <- prior[1] + cur$D
  B <- prior[2] + cur$Y
  if (!is.null(hist) && w > 0) {
    A <- A + w * hist$D
    B <- B + w * hist$Y
  }
  list(A = A, B = B)
}
