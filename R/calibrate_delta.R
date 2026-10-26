#' Threshold giving a target type I error without borrowing
#'
#' Simulates current trials under the null hypothesis (\eqn{\Delta = 0}) and
#' returns the threshold \eqn{\gamma_0} at which the success criterion
#' \eqn{\eta > \gamma_0} without borrowing has type I error \code{alpha}.
#'
#' @param n_current Current sample size per arm.
#' @param alpha Target type I error.
#' @param tau RMST time point in days.
#' @param follow_up Follow-up in days (default \code{tau}).
#' @param R Number of simulated trials.
#' @param prior Reference prior \eqn{(a_0, b_0)}.
#' @param seed Random seed.
#' @param ... Further arguments to \code{\link{sim_arm}} (e.g. \code{dropout}).
#' @return The threshold \eqn{\gamma_0}.
#' @examples
#' gamma_nb(150, R = 2000, seed = 1)
#' @export
gamma_nb <- function(n_current, alpha = 0.025, tau = 365, follow_up = tau, R = 1e5,
                     prior = c(0.001, 1), seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  eta <- vapply(seq_len(R), function(r) {
    a <- sim_arm(n_current, 1, follow_up, ...); b <- sim_arm(n_current, 1, follow_up, ...)
    ca <- bkm_counts(a$time, a$status, tau); cb <- bkm_counts(b$time, b$status, tau)
    ma <- rmst_moments(prior[1] + ca$X, prior[2] + ca$Y)
    mb <- rmst_moments(prior[1] + cb$X, prior[2] + cb$Y)
    stats::pnorm((mb$mean - ma$mean) / sqrt(ma$var + mb$var))
  }, numeric(1))
  sort(eta, decreasing = TRUE)[floor(alpha * R) + 1]
}

#' Design stage calibration of the tolerance
#'
#' Selects the largest tolerance \eqn{\delta} (days) such that the maximum
#' type I error of the success criterion \eqn{\eta > \gamma_0} over a
#' prespecified calibration range \eqn{\Omega_{\mathrm{cal}}} of historical
#' drift does not exceed \code{alpha_cal}, with \eqn{\gamma_0} fixed.
#' The drift \eqn{d_g} is the hazard ratio of historical arm \eqn{g} relative to
#' the corresponding current arm. Borrowing inflates the type I error when the
#' historical control has a higher hazard (\eqn{d_0 > 1}) or the historical
#' treatment arm a lower hazard (\eqn{d_1 < 1}), so the default ranges cover
#' these harmful sides; when both historical arms are available
#' (configuration B), all combinations are used.
#' The tuning constants are \eqn{k_g = \log(2)\, \bar V_g / \delta^2}.
#'
#' @param n_current Current sample size per arm.
#' @param n_hist Named vector of historical sample sizes, e.g.
#'   \code{c(control = 150)} (configuration C), \code{c(treatment = 150)}
#'   (configuration T), or both (configuration B).
#' @param v_design Named vector of design stage values \eqn{\bar V_g}
#'   (see \code{\link{design_variance}}) for the arms with historical data.
#' @param gamma0 Threshold \eqn{\gamma_0} of the success criterion (e.g. from
#'   \code{\link{gamma_nb}}).
#' @param drift_control,drift_treatment Calibration ranges of the drift
#'   \eqn{d_0} (historical control) and \eqn{d_1} (historical treatment).
#' @param method \code{"hazard"} (RG-H) or \code{"rmst"} (RG-R).
#' @param delta_grid Candidate tolerances \eqn{\mathcal{T}} in days.
#' @param alpha_cal Maximum type I error \eqn{\alpha_{\mathrm{cal}}} allowed
#'   over the calibration range.
#' @param R Replicates per drift configuration.
#' @param tau RMST time point in days.
#' @param follow_up Follow-up in days (default \code{tau}).
#' @param prior Reference prior \eqn{(a_0, b_0)}.
#' @param seed Random seed.
#' @param cores Number of cores (parallel on all platforms when > 1).
#' @param ... Further arguments to \code{\link{sim_arm}}.
#' @return A list with the selected \code{delta}, the tuning constants
#'   \code{k}, the maximum type I error over the range for each candidate
#'   (\code{max_type1}), the matrix of estimated type I error rates (drift
#'   configurations by tolerances), and the drift configurations.
#' @examples
#' \donttest{
#' cal <- calibrate_delta(150, c(control = 150), v_design = c(control = 213),
#'                        gamma0 = 0.976, drift_control = c(1, 1.1),
#'                        delta_grid = c(5, 10, 20), R = 200, seed = 1)
#' cal$delta
#' }
#' @export
calibrate_delta <- function(n_current, n_hist, v_design, gamma0,
                            drift_control = c(1, 1.05, 1.10),
                            drift_treatment = c(0.90, 0.95, 1),
                            method = c("hazard", "rmst"),
                            delta_grid = c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60),
                            alpha_cal = 0.04, R = 5000, tau = 365, follow_up = tau,
                            prior = c(0.001, 1), seed = NULL, cores = 1, ...) {
  method <- match.arg(method)
  arms <- names(n_hist)
  if (!length(arms) || !all(arms %in% c("control", "treatment")))
    stop("'n_hist' must be named 'control' and/or 'treatment'")
  if (!all(arms %in% names(v_design))) stop("'v_design' must be given for each historical arm")
  grid <- if (length(arms) == 2) expand.grid(control = drift_control, treatment = drift_treatment) else
    stats::setNames(data.frame(if (arms == "control") drift_control else drift_treatment), arms)
  kfac <- log(2) * v_design[arms]
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, R)
  dots <- list(...)
  one <- function(r) {
    set.seed(seeds[r])
    sim <- function(n, h) { d <- do.call(sim_arm, c(list(n, h, follow_up), dots)); bkm_counts(d$time, d$status, tau) }
    mom <- function(cnt, wcnt = NULL, w = 0) {
      A <- prior[1] + cnt$X; B <- prior[2] + cnt$Y
      if (w > 0) { A <- A + w * wcnt$X; B <- B + w * wcnt$Y }
      rmst_moments(A, B)
    }
    cur <- list(control = sim(n_current, 1), treatment = sim(n_current, 1))
    m1 <- lapply(cur, mom)
    out <- matrix(NA_real_, nrow(grid), length(delta_grid))
    for (j in seq_len(nrow(grid))) {
      his <- lapply(arms, function(a) sim(n_hist[[a]], grid[j, a]))
      names(his) <- arms
      m0 <- lapply(his, mom)
      z2 <- vapply(arms, function(a) (m1[[a]]$mean - m0[[a]]$mean)^2 / (m1[[a]]$var + m0[[a]]$var), numeric(1))
      for (i in seq_along(delta_grid)) {
        mm <- lapply(c("control", "treatment"), function(a) {
          if (!a %in% arms) return(m1[[a]])
          w <- exp(-kfac[[a]] / delta_grid[i]^2 * z2[[a]])
          if (method == "hazard") return(mom(cur[[a]], his[[a]], w))
          prec <- 1 / m1[[a]]$var + w / m0[[a]]$var
          list(mean = (m1[[a]]$mean / m1[[a]]$var + w * m0[[a]]$mean / m0[[a]]$var) / prec, var = 1 / prec)
        })
        out[j, i] <- stats::pnorm((mm[[2]]$mean - mm[[1]]$mean) / sqrt(mm[[1]]$var + mm[[2]]$var))
      }
    }
    out
  }
  res <- par_apply(seq_len(R), one, cores)
  rej <- Reduce(`+`, lapply(res, function(p) p > gamma0)) / R
  dimnames(rej) <- list(NULL, paste0("delta=", delta_grid))
  mx <- apply(rej, 2, max)
  ok <- which(mx <= alpha_cal)
  delta <- if (length(ok)) max(delta_grid[ok]) else NA_real_
  k <- if (is.na(delta)) stats::setNames(rep(Inf, length(arms)), arms) else kfac / delta^2
  list(delta = delta, k = k, max_type1 = stats::setNames(mx, delta_grid), type1 = rej,
       drift = grid, method = method, gamma0 = gamma0, alpha_cal = alpha_cal, R = R)
}
