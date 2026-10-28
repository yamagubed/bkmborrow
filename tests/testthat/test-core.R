test_that("counts follow the event-priority and censoring conventions", {
  cnt <- bkm_counts(time = c(2, 2, 3, 5, 7), status = c(1, 0, 1, 0, 1), tau = 6)
  expect_equal(cnt$N, c(5, 5, 3, 2, 2, 1))
  expect_equal(cnt$X, c(0, 1, 1, 0, 0, 0))
  expect_equal(cnt$Y, cnt$N - cnt$X)
  expect_error(bkm_counts(c(1, -1), c(1, 1), 5))
})

test_that("exact RMST moments match brute-force simulation", {
  set.seed(10)
  K <- 40
  A <- runif(K, 0.5, 3); B <- runif(K, 20, 60)
  for (cv in c("left", "endpoint")) {
    m <- rmst_moments(A, B, cv)
    d <- rmst_draws(A, B, 2e5, cv)
    expect_equal(m$mean, mean(d), tolerance = 2e-3)
    expect_equal(m$var, var(d), tolerance = 2e-2)
  }
})

test_that("left and endpoint conventions differ by 1 - S(K)", {
  A <- rep(1, 20); B <- rep(30, 20)
  ml <- rmst_moments(A, B, "left")$mean
  me <- rmst_moments(A, B, "endpoint")$mean
  expect_equal(ml - me, 1 - prod(B / (A + B)))
})

test_that("weights behave as specified", {
  set.seed(11)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.7))
  his <- list(control = sim_arm(60, 1), treatment = sim_arm(60, 0.7))
  f0 <- bkm_borrow(cur, his, tau = 365, method = "none")
  fk <- bkm_borrow(cur, his, tau = 365, k = Inf)
  expect_equal(f0$diff, fk$diff)
  expect_equal(unname(fk$weights), c(0, 0))
  ff <- bkm_borrow(cur, his, tau = 365, method = "fixed", w = 1)
  cc <- lapply(cur, function(x) bkm_counts(x$time, x$status, 365))
  hh <- lapply(his, function(x) bkm_counts(x$time, x$status, 365))
  pooled <- rmst_moments(0.001 + cc$control$X + hh$control$X, 1 + cc$control$Y + hh$control$Y)
  expect_equal(unname(ff$rmst["control", "mean"]), pooled$mean)
  fk1 <- bkm_borrow(cur, his, tau = 365, k = 0.5)
  z2 <- bkm_weight(cc$treatment, hh$treatment, 0.5)$z2
  expect_equal(unname(fk1$weights["treatment"]), exp(-0.5 * z2))
  expect_true(all(fk1$weights > 0 & fk1$weights < 1))
})

test_that("missing historical arms are not borrowed", {
  set.seed(12)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 1))
  fc <- bkm_borrow(cur, list(control = sim_arm(60, 1)), tau = 365, k = 0.5)
  expect_equal(unname(fc$weights["treatment"]), 0)
  expect_identical(fc$historical_arms, "control")
  ft <- bkm_borrow(cur, list(treatment = sim_arm(60, 1)), tau = 365, k = 0.5)
  expect_equal(unname(ft$weights["control"]), 0)
})

test_that("arm-specific k is respected", {
  set.seed(13)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 1))
  his <- list(control = sim_arm(60, 1), treatment = sim_arm(60, 1))
  f <- bkm_borrow(cur, his, tau = 365, k = c(control = Inf, treatment = 0))
  expect_equal(unname(f$weights), c(0, 1))
})

test_that("normal and draw-based posterior probabilities agree", {
  set.seed(14)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.8))
  his <- list(control = sim_arm(60, 1), treatment = sim_arm(60, 0.8))
  fn <- bkm_borrow(cur, his, tau = 365, k = 0.5, variance = "posterior")
  fd <- bkm_borrow(cur, his, tau = 365, k = 0.5, variance = "posterior", posterior = "draws", ndraw = 20000)
  expect_equal(fn$eta, fd$eta, tolerance = 0.02)
  expect_error(bkm_borrow(cur, his, tau = 365, method = "rmst", variance = "posterior", posterior = "draws"))
  expect_error(bkm_borrow(cur, his, tau = 365, posterior = "draws"))
})

test_that("tolerance and success criterion behave as specified", {
  expect_equal(k_from_delta(10, 100), log(2))
  set.seed(15)
  cur <- list(control = sim_arm(100, 1), treatment = sim_arm(100, 0.7))
  f <- bkm_borrow(cur, tau = 365, method = "none")
  expect_identical(bkm_success(f, 0), TRUE)
  expect_identical(bkm_success(f, 1), FALSE)
  expect_equal(unname(f$diff["upper"] - f$diff["mean"]), qnorm(0.975) * unname(f$diff["sd"]))
})

test_that("follow-up is distinct from the RMST time point", {
  set.seed(17)
  d <- sim_arm(200, 1, follow_up = 730)
  expect_true(max(d$time) > 365)
  cnt <- bkm_counts(d$time, d$status, tau = 365)
  expect_length(cnt$X, 365)
  expect_equal(cnt$N[1], 200)
})

test_that("true RMST matches the closed form for a constant hazard", {
  lam <- log(2) / 365
  x <- 1 - exp(-lam * 0.7)
  expect_equal(true_rmst(0.7), sum((1 - x)^(0:364)))
})

test_that("RG-R combines arm-level RMSTs as a normal power prior", {
  set.seed(16)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.8))
  his <- list(control = sim_arm(100, 1), treatment = sim_arm(100, 0.8))
  f <- bkm_borrow(cur, his, tau = 365, method = "rmst", k = 0.5, variance = "posterior")
  cc <- bkm_counts(cur$control$time, cur$control$status, 365)
  hh <- bkm_counts(his$control$time, his$control$status, 365)
  w <- bkm_weight(cc, hh, 0.5)
  prec <- 1 / w$current$var + w$w / w$historical$var
  expect_equal(unname(f$rmst["control", "mean"]),
               (w$current$mean / w$current$var + w$w * w$historical$mean / w$historical$var) / prec)
  expect_equal(unname(f$rmst["control", "var"]), 1 / prec)
  f0 <- bkm_borrow(cur, his, tau = 365, method = "rmst", k = Inf)
  expect_equal(f0$diff, bkm_borrow(cur, his, tau = 365, method = "none")$diff)
})

test_that("calibrate_delta returns an admissible tolerance", {
  cal <- calibrate_delta(100, c(control = 100), v_design = c(control = 250), gamma0 = 0.976,
                         drift_control = c(1, 1.1), delta_grid = c(2, 10, 40), R = 60, seed = 3)
  expect_true(is.na(cal$delta) || cal$delta %in% c(2, 10, 40))
  expect_equal(dim(cal$type1), c(2, 3))
  calr <- calibrate_delta(100, c(control = 100, treatment = 100), v_design = c(control = 250, treatment = 220),
                          gamma0 = 0.976, delta_grid = c(2, 10), R = 20, seed = 4, method = "rmst")
  expect_equal(nrow(calr$drift), 9)
})

test_that("adjusted variance accounts for the data-dependent weight", {
  m1 <- list(mean = 260, var = 50); m0 <- list(mean = 260, var = 40)
  vr <- function(w) 1 / (1 / m1$var + w / m0$var)
  expect_equal(adjusted_variance(m1, m0, 1, 2), vr(1))                    # D = 0
  m0b <- list(mean = 252, var = 40); w <- exp(-0.5 * 64 / 90)
  expect_equal(adjusted_variance(m1, m0b, w, 0), (1 - w * 50 / (40 + w * 50))^2 * 50 + (w * 50 / (40 + w * 50))^2 * 40)
  expect_gt(adjusted_variance(m1, m0b, w, 0.5), vr(w))
  expect_equal(adjusted_variance(m1, m0b, 0, Inf), 50)
  ## Monte Carlo: adjusted variance tracks the sampling variance of the RG-R estimate
  set.seed(18); V1 <- 100; V0 <- 100; k <- 1
  a <- rnorm(4e4, 0, 10); b <- rnorm(4e4, 0, 10); D <- b - a
  w <- exp(-k * D^2 / (V1 + V0)); est <- a + w * V1 / (V0 + w * V1) * D
  va <- mapply(function(x, y, ww) adjusted_variance(list(mean = x, var = V1), list(mean = y, var = V0), ww, k), a, b, w)
  expect_lt(mean(1 / (1 / V1 + w / V0)), var(est))
  expect_gt(mean(va), var(est))
  ## in bkm_borrow
  set.seed(19)
  cur <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.8))
  his <- list(control = sim_arm(100, 1.1), treatment = sim_arm(100, 0.8))
  fa <- bkm_borrow(cur, his, tau = 365, k = 0.5); fp <- bkm_borrow(cur, his, tau = 365, k = 0.5, variance = "posterior")
  expect_equal(fa$diff[["mean"]], fp$diff[["mean"]])
  expect_gte(fa$diff[["sd"]], fp$diff[["sd"]])
  f0 <- bkm_borrow(cur, his, tau = 365, k = Inf)
  expect_equal(f0$diff, bkm_borrow(cur, his, tau = 365, method = "none")$diff)
})

test_that("power bound matches the linear programme and reduces to no borrowing", {
  v <- c(control = 52.4, treatment = 44.6)
  expect_equal(power_bound(25.84, v, 52.4, "control", 7.9)$power, 0.832, tolerance = 1e-3)
  expect_equal(power_bound(25.84, v, 52.4, "control", 4.0)$power, 0.866, tolerance = 1e-3)
  b <- power_bound(25.84, v, 52.4, "control", 20)
  expect_equal(b$power, b$power_nb); expect_equal(b$lambda, 0)
  expect_gt(power_bound(25.84, v, 52.4, "control", 0)$power, power_bound(25.84, v, 52.4, "control", 5)$power)
  ## the attaining fixed-weight test has size alpha at the edge and power equal to the bound
  r <- power_bound(25.84, v, 52.4, "treatment", 3); s2 <- sum(v); V1 <- v[["treatment"]]
  vt <- s2 + 2 * r$lambda * (-V1) + r$lambda^2 * (V1 + 52.4)
  crit <- r$lambda * 3 + qnorm(0.96) * sqrt(vt)
  expect_equal(1 - pnorm((crit - 25.84) / sqrt(vt)), r$power, tolerance = 1e-8)
})
