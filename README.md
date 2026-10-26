# bkmborrow

Dynamic historical borrowing for Bayesian Kaplan-Meier inference on the restricted mean survival time (RMST).

The daily hazards of each arm get independent `Beta(a0, b0)` priors (the reference prior; Hjort 1990; Pan et al. 2026; Paydarfar et al. 2026). For each arm g with historical data, the weight

    w_g = exp(-k_g Z_g^2),   Z_g = (mu_hat_1g - mu_hat_0g) / sqrt(V_1g + V_0g),

compares the current (s = 1) and historical (s = 0) posterior RMSTs. RG-H multiplies the historical daily event and event-free counts by w_g, so the posterior stays conjugate; RG-R combines the arm-level RMSTs as a normal power prior with the same weight. Historical data may be available for the control arm only (configuration C), the treatment arm only (T), or both arms (B).

The tuning constant is expressed as a tolerance delta in days, k_g = log(2) bar V_g / delta^2, the discrepancy at which history receives half weight. With the threshold gamma0 of the success criterion eta = Pr(Delta > 0 | D) > gamma0 fixed so that no borrowing has type I error alpha, delta is calibrated at the design stage as the largest tolerance whose maximum type I error over a prespecified range of harmful drift is at most alpha_cal.

## Main functions

| Function | Purpose |
|---|---|
| `bkm_counts()` | daily events `X`, event-free counts `Y` and numbers at risk `N` up to the RMST time point tau |
| `rmst_moments()`, `rmst_draws()` | exact posterior RMST moments (O(tau)) and posterior draws |
| `bkm_borrow()` | analysis with RG-H (`"hazard"`), RG-R (`"rmst"`), a fixed power prior or no borrowing; returns `eta` and the 95% credible interval |
| `bkm_success()` | success criterion `eta > gamma0` |
| `k_from_delta()`, `design_variance()` | tuning constant from a tolerance in days; design stage value bar V_g |
| `gamma_nb()` | threshold gamma0 at which no borrowing has type I error alpha |
| `calibrate_delta()` | design stage calibration of the tolerance over a harmful-side drift range |
| `bkm_weight()` | the RMST-guided weight for one arm |
| `sim_arm()`, `true_rmst()` | daily-grid simulation helpers (follow-up is distinct from tau) |

## Example

```r
library(bkmborrow)
set.seed(1)
current    <- list(control = sim_arm(310, 1), treatment = sim_arm(310, 0.7))
historical <- list(control = sim_arm(200, 1.1), treatment = sim_arm(200, 0.75))
k   <- k_from_delta(20, c(control = 132, treatment = 113))
fit <- bkm_borrow(current, historical, tau = 365, k = k)
fit
bkm_success(fit, gamma0 = 0.976)
```

To borrow only the control arm, pass `historical = list(control = ...)`; to borrow only the treatment arm, pass `historical = list(treatment = ...)`. `inst/examples/workflow.R` contains the complete design-and-analysis workflow used in the paper, including calibration.

## Changes in 0.3.0

* Notation aligned with the paper: daily event counts are `X` (was `D`); the posterior probability is `eta` (was `prob`); the argument `prob` is now `posterior`.
* `bkm_success()` replaces `bkm_decision()`.
* `sim_arm()`, `design_variance()`, `gamma_nb()` and `calibrate_delta()` take a `follow_up` argument separate from the RMST time point `tau`.
* `calibrate_gamma()` and `rejection_rate()` (threshold calibration) were removed; the method fixes gamma0 and calibrates the tolerance.
* `calibrate_delta(cores = )` runs in parallel on Windows as well.
