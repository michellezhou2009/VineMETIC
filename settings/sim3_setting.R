## settings/sim3_setting.R
##
## Configuration for Simulation Study 3 — comparison of VineSurvIC versus
## the Li (2020) pseudo-MLE under a trivariate Clayton DGP without covariates.
##
## DGP:
##   (T1, T2, D) ~ trivariate Clayton with single parameter alpha
##   Marginals: Weibull(shape=2) with scale 70 / 60 / 85 for T1 / T2 / D
##   Administrative censoring: Ca ~ Exp(rate = 1/350)
##
## Two true tau levels:
##   tau = 0.3  =>  alpha = 2*0.3/(1-0.3) = 6/7  ~= 0.8571
##   tau = 0.6  =>  alpha = 2*0.6/(1-0.6) = 3.0
##
## Methods compared:
##   (A) Li (2020): iterative pseudo-MLE assuming ONE shared alpha
##       for all pairs.  Targets: tau(T1,D) = tau(T2,D) = tau(T1,T2|D).
##   (B) VineSurvIC: vine copula with all-Clayton intercept-only pairs.
##       Three separate edge parameters; each edge tau estimated
##       independently.
##
## Note on identifiability:
##   For trivariate Clayton with parameter alpha, all three bivariate
##   marginals are bivariate Clayton with the SAME alpha (Archimedean
##   exchangeability).  Hence ALL three vine-edge taus have the same
##   true value alpha / (alpha + 2).

## ------------------------------------------------------------------
## Two true tau values
## ------------------------------------------------------------------
tau_vals  <- c(0.3, 0.6)
alpha_vals <- 2 * tau_vals / (1 - tau_vals)   # Clayton: tau = alpha/(alpha+2)

## ------------------------------------------------------------------
## VineSurvIC fit spec: all-Clayton, intercept-only W_formula
## ------------------------------------------------------------------
suppressPackageStartupMessages(library(VineSurvIC))

vine_spec_sim3 <- VineSurvIC::default_vine_j3(
  families   = c("Clayton", "Clayton", "Clayton"),
  links      = c("log", "log", "log"),
  W_formulas = list(~ 1, ~ 1, ~ 1)
)

## ------------------------------------------------------------------
## Print summary
## ------------------------------------------------------------------
cat("\n=== Sim 3 setting ===\n")
cat("DGP copula : trivariate Clayton (single alpha for ALL pairs)\n")
cat("Marginals  : Weibull(shape=2, scale=70/60/85) for T1/T2/D\n")
cat("Censoring  : Exp(rate = 1/350) administrative\n")
cat("True tau   :", paste(tau_vals,   collapse = ", "), "\n")
cat("True alpha :", paste(round(alpha_vals, 4), collapse = ", "), "\n")
cat("Methods    : Li (2020) pseudo-MLE, VineSurvIC vine (Clayton x3)\n\n")

## ------------------------------------------------------------------
## Save
## ------------------------------------------------------------------
setting3 <- list(
  tau_vals   = tau_vals,
  alpha_vals = alpha_vals,
  vine_spec  = vine_spec_sim3
)

saveRDS(setting3, file = "settings/sim3_setting.rds")
cat("Setting saved to settings/sim3_setting.rds\n")
