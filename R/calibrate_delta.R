#' Threshold giving a target type I error without borrowing
#'
#' Simulates current trials under the null hypothesis and returns the
#' posterior-probability threshold \eqn{\gamma_0} at which the analysis
#' without borrowing has type I error \code{alpha}.
#'
#' @param n_current Current sample size per arm.
#' @param alpha Target type I error.
#' @param tau Horizon in days.
#' @param R Number of simulated trials.
#' @param prior Daily beta prior.
#' @param seed Random seed.
#' @param ... Further arguments to \code{\link{sim_arm}} (e.g. \code{dropout}).
#' @return The threshold \eqn{\gamma_0}.
#' @examples
#' gamma_nb(150, R = 2000, seed = 1)
#' @export
gamma_nb <- function(n_current, alpha = 0.025, tau = 365, R = 1e5, prior = c(0.001, 1),
                     seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  pp <- vapply(seq_len(R), function(r) {
    a <- sim_arm(n_current, 1, tau, ...); b <- sim_arm(n_current, 1, tau, ...)
    ca <- bkm_counts(a$time, a$status, tau); cb <- bkm_counts(b$time, b$status, tau)
    ma <- rmst_moments(prior[1] + ca$D, prior[2] + ca$Y)
    mb <- rmst_moments(prior[1] + cb$D, prior[2] + cb$Y)
    stats::pnorm((mb$mean - ma$mean) / sqrt(ma$var + mb$var))
  }, numeric(1))
  sort(pp, decreasing = TRUE)[floor(alpha * R) + 1]
}

#' Design-stage calibration of the tolerance
#'
#' Selects the largest tolerance \eqn{\delta} (days) such that the maximum
#' type I error over a prespecified calibration range of historical drift does
#' not exceed \code{alpha_cal}, with the decision threshold fixed at
#' \code{gamma0}. Drift is the hazard ratio of a historical arm relative to the
#' corresponding current arm. Type I error inflation arises when the historical
#' control is worse (drift > 1) or the historical treatment is better
#' (drift < 1), so the default ranges cover these harmful sides; when both
#' historical arms are available, all combinations are used.
#' The tuning constants are \eqn{k_g = \log(2) V_g / \delta^2}.
#'
#' @param n_current Current sample size per arm.
#' @param n_hist Named vector of historical sample sizes, e.g.
#'   \code{c(control = 150)}, \code{c(treatment = 150)} or both.
#' @param v_design Named vector of design variances \eqn{V_g}
#'   (see \code{\link{design_variance}}) for the arms with historical data.
#' @param gamma0 Decision threshold (e.g. from \code{\link{gamma_nb}}).
#' @param drift_control,drift_treatment Calibration ranges of drift hazard
#'   ratios for the historical control and treatment arms.
#' @param method \code{"hazard"} (RG-H) or \code{"rmst"} (RG-R).
#' @param delta_grid Candidate tolerances in days.
#' @param alpha_cal Calibration tolerance for the type I error.
#' @param R Replicates per drift configuration.
#' @param tau Horizon in days.
#' @param prior Daily beta prior.
#' @param seed Random seed.
#' @param cores Number of cores (uses \code{parallel::mclapply} when > 1).
#' @param ... Further arguments to \code{\link{sim_arm}}.
#' @return A list with the selected \code{delta}, the tuning constants
#'   \code{k}, the matrix of estimated type I error rates (drift
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
                            alpha_cal = 0.04, R = 5000, tau = 365, prior = c(0.001, 1),
                            seed = NULL, cores = 1, ...) {
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
  mom <- function(cnt, wcnt = NULL, w = 0) {
    A <- prior[1] + cnt$D; B <- prior[2] + cnt$Y
    if (w > 0) { A <- A + w * wcnt$D; B <- B + w * wcnt$Y }
    rmst_moments(A, B)
  }
  one <- function(r) {
    set.seed(seeds[r])
    cur <- lapply(c(control = 1, treatment = 1), function(h) {
      d <- sim_arm(n_current, h, tau, ...); bkm_counts(d$time, d$status, tau) })
    m0 <- lapply(cur, mom)
    out <- matrix(NA_real_, nrow(grid), length(delta_grid))
    for (j in seq_len(nrow(grid))) {
      his <- lapply(arms, function(a) { d <- sim_arm(n_hist[[a]], grid[j, a], tau, ...); bkm_counts(d$time, d$status, tau) })
      names(his) <- arms
      mh_all <- lapply(his, mom)
      z2 <- vapply(arms, function(a) (m0[[a]]$mean - mh_all[[a]]$mean)^2 / (m0[[a]]$var + mh_all[[a]]$var), numeric(1))
      for (i in seq_along(delta_grid)) {
        mm <- lapply(c("control", "treatment"), function(a) {
          if (!a %in% arms) return(m0[[a]])
          w <- exp(-kfac[[a]] / delta_grid[i]^2 * z2[[a]])
          if (method == "hazard") return(mom(cur[[a]], his[[a]], w))
          mh <- mh_all[[a]]; prec <- 1 / m0[[a]]$var + w / mh$var
          list(mean = (m0[[a]]$mean / m0[[a]]$var + w * mh$mean / mh$var) / prec, var = 1 / prec) })
        out[j, i] <- stats::pnorm((mm[[2]]$mean - mm[[1]]$mean) / sqrt(mm[[1]]$var + mm[[2]]$var))
      }
    }
    out
  }
  res <- if (cores > 1) parallel::mclapply(seq_len(R), one, mc.cores = cores) else lapply(seq_len(R), one)
  rej <- Reduce(`+`, lapply(res, function(p) p > gamma0)) / R
  dimnames(rej) <- list(NULL, paste0("delta=", delta_grid))
  mx <- apply(rej, 2, max)
  ok <- which(mx <= alpha_cal)
  delta <- if (length(ok)) max(delta_grid[ok]) else NA_real_
  k <- if (is.na(delta)) stats::setNames(rep(Inf, length(arms)), arms) else kfac / delta^2
  list(delta = delta, k = k, max_type1 = stats::setNames(mx, delta_grid), type1 = rej,
       drift = grid, method = method, gamma0 = gamma0, alpha_cal = alpha_cal, R = R)
}
