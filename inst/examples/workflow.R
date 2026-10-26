## Worked example: design stage calibration of the tolerance and analysis with RG-H
## (hypothetical data generated with sim_arm; configuration B, both arms borrowed)
library(bkmborrow)
tau <- 365; n1 <- 310; n0 <- 200; alpha <- 0.025; alpha_cal <- 0.04
cores <- 2   # increase on a multi-core machine

## 1. Threshold gamma0 at which no borrowing has type I error alpha
gamma0 <- gamma_nb(n1, alpha = alpha, tau = tau, R = 50000, seed = 11)

## 2. Design stage values bar V_g (planning hazard ratio 0.7)
set.seed(101)
vd <- c(control   = design_variance(n1, n0, hazard_mult = 1,   tau = tau, nsim = 2000),
        treatment = design_variance(n1, n0, hazard_mult = 0.7, tau = tau, nsim = 2000))

## 3. Calibrate the tolerance over the local harmful-side drift range (both arms jointly)
cal <- calibrate_delta(n1, c(control = n0, treatment = n0), vd, gamma0,
                       drift_control = c(1.00, 1.05, 1.10), drift_treatment = c(0.90, 0.95, 1.00),
                       alpha_cal = alpha_cal, R = 2000, tau = tau, seed = 202, cores = cores)

## 4. Analysis of one hypothetical trial pair
set.seed(303)
current <- list(control = sim_arm(n1, 1), treatment = sim_arm(n1, 0.7))
historical <- list(control = sim_arm(n0, 1.05), treatment = sim_arm(n0, 0.7))
fit <- bkm_borrow(current, historical, tau = tau, method = "hazard", k = cal$k)
fit_nb <- bkm_borrow(current, historical, tau = tau, method = "none")
success <- bkm_success(fit, gamma0)
