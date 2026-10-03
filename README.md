# bkmborrow

Dynamic historical borrowing for Bayesian Kaplan-Meier (daily-beta) inference on the restricted mean survival time (RMST).

The current trial's daily hazards get independent `Beta(a0, b0)` priors (Hjort 1990; Pan et al. 2026; Paydarfar et al. 2026). Each historical arm's daily event and event-free counts are discounted by one weight

    w_g = exp(-k_g Z_g^2),   Z_g = (RMST_current - RMST_historical) / sqrt(V_current + V_historical),

so the posterior stays conjugate (RG-H). Alternatively the arm-level RMSTs themselves are combined as a normal power prior with the same weight (RG-R). Historical data may be available for the control arm only, the treatment arm only, or both arms. The tuning constant can be set from a clinically elicited tolerance `delta` (days). The posterior-probability threshold is calibrated by simulation to control the maximum type I error over a prespecified conflict set.

## Main functions

| Function | Purpose |
|---|---|
| `bkm_counts()` | daily sufficient statistics from (time, status) |
| `rmst_moments()`, `rmst_draws()` | exact posterior RMST moments (O(K)) and posterior draws |
| `bkm_borrow()` | analysis with RG-H (`"hazard"`, hazard-level), RG-R (`"rmst"`, RMST-level), fixed power or no borrowing |
| `k_from_delta()`, `design_variance()` | tuning constant from an elicited tolerance |
| `gamma_nb()` | threshold at which no borrowing has type I error alpha |
| `calibrate_delta()` | design-stage calibration of the tolerance delta (and k) to a type I error tolerance over a drift range |
| `calibrate_gamma()`, `rejection_rate()` | threshold calibration (alternative design) |
| `sim_arm()`, `true_rmst()` | daily-grid simulation helpers |

## Example

```r
library(bkmborrow)
set.seed(1)
current    <- list(control = sim_arm(150, 1), treatment = sim_arm(150, 0.7))
historical <- list(control = sim_arm(100, 1.1), treatment = sim_arm(100, 0.75))
k   <- k_from_delta(20, c(control = 265, treatment = 226))
fit <- bkm_borrow(current, historical, tau = 365, k = k)
fit
```

To borrow only the control arm, pass `historical = list(control = ...)`; to borrow only the treatment arm, pass `historical = list(treatment = ...)`. `inst/examples/workflow.R` contains a complete design-and-analysis workflow, including calibration.
