#' Exact posterior moments of the RMST
#'
#' For independent daily hazards \eqn{p_t \sim \mathrm{Beta}(A_t, B_t)},
#' \eqn{t = 1,\dots,\tau}, computes the exact mean and variance of the
#' restricted mean survival time \eqn{\mu = \sum_{t=0}^{\tau-1} S(t)}, where
#' \eqn{S(t) = \prod_{u \le t}(1 - p_u)}. With \eqn{q_{1t} = B_t/(A_t + B_t)},
#' \eqn{q_{2t} = q_{1t}(B_t + 1)/(A_t + B_t + 1)},
#' \eqn{m_1(t) = \prod_{u \le t} q_{1u}} and \eqn{m_2(t) = \prod_{u \le t} q_{2u}},
#' \deqn{E(\mu) = \sum_{t=0}^{\tau-1} m_1(t),\quad
#' E(\mu^2) = \sum_{t=0}^{\tau-1} m_2(t) + 2\sum_{t=0}^{\tau-1}
#' \frac{m_2(t)}{m_1(t)} \sum_{v=t+1}^{\tau-1} m_1(v).}
#' The computation needs \eqn{O(\tau)} operations. With
#' \code{convention = "endpoint"}, the RMST is instead
#' \eqn{\sum_{t=1}^{\tau} S(t)} (Paydarfar et al., 2026).
#'
#' @param A,B Vectors of beta shape parameters, one per day \eqn{t = 1,\dots,\tau}.
#' @param convention \code{"left"} (default) or \code{"endpoint"}.
#' @return A list with \code{mean} and \code{var}.
#' @examples
#' rmst_moments(A = rep(1, 30), B = rep(100, 30))
#' @export
rmst_moments <- function(A, B, convention = c("left", "endpoint")) {
  convention <- match.arg(convention)
  tau <- length(A)
  if (length(B) != tau) stop("'A' and 'B' differ in length")
  if (any(A <= 0) || any(B <= 0)) stop("beta parameters must be positive")
  q1 <- B / (A + B)
  q2 <- q1 * (B + 1) / (A + B + 1)
  if (convention == "left") {
    m1 <- c(1, cumprod(q1)[-tau])
    m2 <- c(1, cumprod(q2)[-tau])
    r <- c(1, cumprod(q2 / q1)[-tau])            # m2 / m1, computed directly
  } else {
    m1 <- cumprod(q1)
    m2 <- cumprod(q2)
    r <- cumprod(q2 / q1)
  }
  mu <- sum(m1)
  tail <- rev(cumsum(rev(m1))) - m1               # sum_{v > t} m1(v)
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
  tau <- length(A)
  S <- rep(1, ndraw)
  total <- if (convention == "left") rep(1, ndraw) else rep(0, ndraw)
  last <- if (convention == "left") tau - 1L else tau
  for (t in seq_len(last)) {
    S <- S * (1 - stats::rbeta(ndraw, A[t], B[t]))
    total <- total + S
  }
  total
}

## internal: beta parameters of the (power-prior) posterior of the current daily hazards
arm_posterior <- function(cur, hist, w, prior) {
  A <- prior[1] + cur$X
  B <- prior[2] + cur$Y
  if (!is.null(hist) && w > 0) {
    A <- A + w * hist$X
    B <- B + w * hist$Y
  }
  list(A = A, B = B)
}
