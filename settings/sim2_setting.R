## sim2_setting.R
##
## Configuration for Simulation Study II — copula-family
## misspecification (R&R reviewer's Major Comment 2(b)).
##
## DGP: identical to Sim 1 — Gumbel C(1,3) at intercept-only tau = 0.7,
##      Clayton C(2,3) at intercept-only tau = 0.4, Gaussian C(1,2|3)
##      with covariate-dependent tau ranging over (-0.3, 0.3) across the
##      (z1, z2) support. Calibrated marginal intercepts (gammas) so the
##      censoring rates land at (10%, 30%, 25%) for (TD, T1, T2).
##
## Four misspecified-fit sub-scenarios. Each scenario is named by its
## fit-time family combination (C(1,3) _ C(2,3) _ C(1,2|3)), so the
## output folder under raw_results/sim2/ is self-documenting.
##
##   Scenario name              C(1,3)    C(2,3)    C(1,2|3)   Type of misspec
##   ---------------------------------------------------------------------------
##   Clayton_Clayton_Gaussian   Clayton✗  Clayton✓  Gaussian✓  only C(1,3) wrong (cleanest one-edge)
##   Clayton_Gumbel_Gaussian    Clayton✗  Gumbel ✗  Gaussian✓  both first-layer wrong (cross-swap)
##   Gumbel_Clayton_Frank       Gumbel ✓  Clayton✓  Frank   ✗  only C(1,2|3) wrong (second layer)
##   Clayton_Clayton_Frank      Clayton✗  Clayton✓  Frank   ✗  combined first + second layer wrong
##
## ✓ = matches truth, ✗ = misspecified vs truth (Gumbel/Clayton/Gaussian).
##
## Why these four (vs other candidates we considered):
##   - C-C-G isolates the effect of getting just C(1,3) wrong (cleanest
##     "single first-layer edge wrong" demonstration).
##   - C-G-G adds the effect of getting both first-layer pairs wrong;
##     also exposes propagation into C(1,2|3) via plug-in.
##   - G-C-Frank: only C(1,2|3) is wrong. Frank can hit negative tau
##     (so it's a "feasible" misspec, no boundary issues), but lives
##     on a totally different link (identity vs fisher_z), which makes
##     the parameter-scale BIAS huge while coverage on C(1,2|3) drops
##     to ~5-40%.
##   - C-C-Frank: combined first-layer + second-layer misspecification.
##     Adds a scenario the other three don't cover: how do the two
##     misspecification effects compose?
##
## SCENARIOS DROPPED (and why):
##   - Gumbel_Clayton_Clayton: Clayton fitted on C(1,2|3) with a truth
##     covering negative tau forces the optimizer to the parameter
##     boundary on most reps, taking many minutes per rep. Same
##     qualitative message as G-C-Frank, just much slower.
##   - Clayton_Clayton_Clayton: same boundary issue. Plus the message
##     is largely a superposition of C-C-G + G-C-Clayton.
##   - Gumbel_Clayton_Gumbel: same positive-only-family boundary issue
##     on C(1,2|3) as G-C-Clayton; redundant.

suppressPackageStartupMessages({
  library(VineSurvIC)
  library(VineCopula)
})

## ------------------------------------------------------------------
## Load Sim 1's DGP. We reuse: betas, vine_spec, copula_pars, gammas,
## c.lwr, c.upr, tt.all. With the same seed master across Sim 1 and
## Sim 2, every rep r produces the identical generated dataset --
## common random numbers (CRN) gives the cleanest variance reduction
## in scenario-to-scenario comparisons.
## ------------------------------------------------------------------
sim1_setting_path <- "settings/sim1_setting.rds"
if (!file.exists(sim1_setting_path)) {
  ## Build Sim 1's setting first if it isn't there.
  source("settings/sim1_setting.R", chdir = FALSE)
}
sim1 <- readRDS(sim1_setting_path)

## ------------------------------------------------------------------
## The five misspecified-fit vine_specs. All use the same regression
## formula (~ z1 + z2) on every edge to match Sim 1; only the family
## list differs. Links are NULL so the package picks each family's
## default registered link (log-1 for Gumbel, log for Clayton,
## fisher_z for Gaussian, identity for Frank).
## ------------------------------------------------------------------
make_spec <- function(f_13, f_23, f_12_3) {
  list(
    j1J  = list(family = f_13,   link = NULL, W_formula = ~ z1 + z2),
    j2J  = list(family = f_23,   link = NULL, W_formula = ~ z1 + z2),
    j12J = list(family = f_12_3, link = NULL, W_formula = ~ z1 + z2)
  )
}

## Scenario name == "{C(1,3)}_{C(2,3)}_{C(1,2|3)}" — same names that
## appear under raw_results/sim2/ so a folder listing is also a
## scenario summary. Order below = recommended run order (lighter
## misspecification first, combined misspecification last).
sim2_specs <- list(
  "Clayton_Clayton_Gaussian" = make_spec("Clayton", "Clayton", "Gaussian"),
  "Clayton_Gumbel_Gaussian"  = make_spec("Clayton", "Gumbel",  "Gaussian"),
  "Gumbel_Clayton_Frank"     = make_spec("Gumbel",  "Clayton", "Frank"),
  "Clayton_Clayton_Frank"    = make_spec("Clayton", "Clayton", "Frank")
)

## ------------------------------------------------------------------
## Print a summary so a sanity check is easy from the console.
## ------------------------------------------------------------------
true_fams <- c(sim1$vine_spec$j1J$family,
               sim1$vine_spec$j2J$family,
               sim1$vine_spec$j12J$family)
cat("\n=== Sim 2 setting summary ===\n")
cat("DGP families (true):  ", paste(true_fams, collapse = " | "), "\n",
    sep = "")
cat("Fit scenarios:\n")
maxw <- max(nchar(names(sim2_specs)))
for (s in names(sim2_specs)) {
  fit_fams <- c(sim2_specs[[s]]$j1J$family,
                sim2_specs[[s]]$j2J$family,
                sim2_specs[[s]]$j12J$family)
  flags <- ifelse(fit_fams == true_fams, " ", "*")
  cat("  ", formatC(s, width = maxw, flag = "-"), ":  ",
      sprintf("%-9s%s | %-9s%s | %-9s%s",
              fit_fams[1], flags[1], fit_fams[2], flags[2],
              fit_fams[3], flags[3]),
      "\n", sep = "")
}
cat("(* = misspecified vs DGP)\n")

## ------------------------------------------------------------------
## Save the setting so simulation drivers can readRDS() it.
## ------------------------------------------------------------------
saveRDS(
  list(
    dgp = list(
      betas       = sim1$betas,
      vine_spec   = sim1$vine_spec,
      copula_pars = sim1$copula_pars,
      gammas      = sim1$gammas,
      c.lwr       = sim1$c.lwr,
      c.upr       = sim1$c.upr,
      tt.all      = sim1$tt.all
    ),
    fit_specs = sim2_specs
  ),
  file = file.path("settings", "sim2_setting.rds")
)
cat("\nSetting saved to LiDA/simulation/settings/sim2_setting.rds\n")
