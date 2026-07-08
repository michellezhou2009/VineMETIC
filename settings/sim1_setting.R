## sim1_setting.R
##
## Heterogeneous-families configuration for Simulation Study I of the
## R&R. Reproduces simulation_new/setting_(1,1).R *exactly*, except
## C(1,2|3) is changed from Frank to Gaussian (reviewer Major 2(a)).
##
##   Pair          Family    Default link   tau target
##   --------------------------------------------------------------
##   C(1,3)        Gumbel    log-1          0.7   (intercept-only)
##   C(2,3)        Clayton   log            0.4   (intercept-only)
##   C(1,2|3)      Gaussian  fisher_z       range (-0.3, 0.3)
##                                          across (z1, z2) support
##
## All other pieces (regression coefficients, gamma slopes, censoring
## targets, c.lwr/c.upr) match setting_(1,1).R verbatim.
##
## The marginal intercepts gammas = list(t1, t2, tD) are calibrated via
## find.gammaD() / find.gammaT() from simulation_new/DataGenFuns.R so
## the censoring rates land at the original targets (10%, 30%, 25%).
## The gammas computed here depend ONLY on C(1,3) and C(2,3) (Gumbel
## and Clayton), not on C(1,2|3); changing the conditional family from
## Frank to Gaussian leaves them unchanged from your saved
## setting_(1,1).RData. We re-derive them from scratch for clarity.

suppressPackageStartupMessages({
  library(VineSurvIC)
  library(VineCopula)
})

## Source the censoring-rate calibration helpers (find.gammaD,
## find.gammaT, quantileT.fun). These are simulation-script-flavored
## helpers; they live in the harness, not the package.
source("utils/calibration.R")

## Backwards-compatible names used by find.gammaT (which calls
## MylinkFun and MyCopIndex internally). These are exported by
## VineSurvIC under MylinkFun / MyCopIndex names but live in the
## internal namespace; expose them here.
MylinkFun  <- VineSurvIC:::MylinkFun
MyCopIndex <- VineSurvIC:::MyCopIndex
MyCop      <- VineSurvIC:::MyCop

## ------------------------------------------------------------------
## Marginal regression coefficients (verbatim from setting_(1,1).R)
## ------------------------------------------------------------------
betas <- list(t1 = c(2, 2),
              t2 = c(2, 2),
              tD = c(2, 2))

## ------------------------------------------------------------------
## Vine spec: families, links, copula-parameter regression formulas
## ------------------------------------------------------------------
vine_spec <- list(
  j1J  = list(family = "Gumbel",   link = NULL, W_formula = ~ z1 + z2),
  j2J  = list(family = "Clayton",  link = NULL, W_formula = ~ z1 + z2),
  j12J = list(family = "Gaussian", link = NULL, W_formula = ~ z1 + z2)
)
copula.fams  <- matrix(NA_character_, 3, 3)
copula.fams[1, 3] <- "Gumbel"
copula.fams[2, 3] <- "Clayton"
copula.fams[1, 2] <- "Gaussian"
copula.links <- matrix(NA_character_, 3, 3)
copula.links[1, 3] <- "log-1"
copula.links[2, 3] <- "log"
copula.links[1, 2] <- "fisher_z"

`%||%` <- function(a, b) if (is.null(a)) b else a

## ------------------------------------------------------------------
## True gamma vectors (intercept, Z1 slope, Z2 slope) per copula edge
## ------------------------------------------------------------------
gamma1D.true <- c(
  MylinkFun(copula.links[1, 3])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[1, 3]), tau = 0.7)),
  1, 0.1)
gamma2D.true <- c(
  MylinkFun(copula.links[2, 3])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[2, 3]), tau = 0.4)),
  0.1, 1)
## Gaussian C(1,2|3): chosen so the implied Kendall's tau ranges
## exactly over (-0.3, 0.3) across the covariate support
## (z1 in [-0.5, 0.5], z2 in {0, 1}). With link = fisher_z and equal
## slopes gamma_1 = gamma_2 = 0.4894, intercept gamma_0 = -0.2447,
## the Fisher-z linear predictor at (z1=-0.5, z2=0) gives rho such
## that tau = -0.3, and at (z1= 0.5, z2=1) gives tau = +0.3.
gamma12.true <- c(-0.2447, 0.4894, 0.4894)

copula_pars <- list(
  "(1,3)" = gamma1D.true,
  "(2,3)" = gamma2D.true,
  "(1,2)" = gamma12.true
)

## ------------------------------------------------------------------
## Censoring distribution (verbatim from setting_(1,1).R)
## ------------------------------------------------------------------
c.lwr <- 1
c.upr <- 6

## ------------------------------------------------------------------
## Calibrated marginal intercepts so censoring rates land at targets
## ------------------------------------------------------------------
cat("Calibrating gamma.tD for 10% censoring on T_D...\n")
gamma.tD <- find.gammaD(cen.rate = 0.1)

cat("Calibrating gamma.t1 for 30% censoring on T_1 ...\n")
gamma.t1 <- find.gammaT(cen.rate = 0.3, gamma.tD = gamma.tD,
                        type = "t1",
                        copula.fam  = copula.fams[1, 3],
                        copula.link = copula.links[1, 3],
                        copula.par  = copula_pars$"(1,3)")

cat("Calibrating gamma.t2 for 25% censoring on T_2 ...\n")
gamma.t2 <- find.gammaT(cen.rate = 0.25, gamma.tD = gamma.tD,
                        type = "t2",
                        copula.fam  = copula.fams[2, 3],
                        copula.link = copula.links[2, 3],
                        copula.par  = copula_pars$"(2,3)")

gammas <- list(t1 = gamma.t1, t2 = gamma.t2, tD = gamma.tD)

## ------------------------------------------------------------------
## Time grid for evaluating baseline survival S_j(t) at probability
## levels pp (matches simulation_new/setting_(1,1).R:68-78). For each
## probability p in pp, we solve S_{j,0}(t) = p for t under the
## calibrated marginal intercept gammas$t_j, giving the "true" times
## tt_{1,p}, tt_{2,p}, tt_{D,p}. The simulation harness then evaluates
## the *estimated* baseline survival at these times and compares to p
## -- the basis of Web Table S1.
## ------------------------------------------------------------------
pp  <- seq(0.05, 0.95, by = 0.01)
tt1 <- vapply(pp, function(p) quantileT.fun(p, gammaT = gammas$t1),
              numeric(1))
tt2 <- vapply(pp, function(p) quantileT.fun(p, gammaT = gammas$t2),
              numeric(1))
ttD <- vapply(pp, function(p) quantileT.fun(p, gammaT = gammas$tD,
                                            upr = 200),
              numeric(1))
tt.all <- list(pp = pp, t1 = tt1, t2 = tt2, tD = ttD)

## ------------------------------------------------------------------
## Implied (alpha, tau) range across the covariate support
## ------------------------------------------------------------------
tau_range <- function(spec, gamma_vec) {
  rec  <- VineSurvIC:::bicop_get(spec$family)
  link <- VineSurvIC:::link_funs(spec$link %||% rec$default_link)
  lp_lwr <- gamma_vec[1] - 0.5 * gamma_vec[2] + 0 * gamma_vec[3]
  lp_upr <- gamma_vec[1] + 0.5 * gamma_vec[2] + 1 * gamma_vec[3]
  c(tau_lwr = BiCopPar2Tau(rec$family_code, link$h.fun(lp_lwr)),
    tau_upr = BiCopPar2Tau(rec$family_code, link$h.fun(lp_upr)))
}
tau.all <- list(
  cop1D   = tau_range(vine_spec$j1J,  gamma1D.true),
  cop2D   = tau_range(vine_spec$j2J,  gamma2D.true),
  cop12gD = tau_range(vine_spec$j12J, gamma12.true)
)

cat("\n=== Sim 1 setting summary ===\n")
cat("Marginals: PH, betas = (2, 2) for each event\n")
cat("Vine families: ",
    paste(c(vine_spec$j1J$family,
            vine_spec$j2J$family,
            vine_spec$j12J$family), collapse = " | "), "\n", sep = "")
cat("True gamma vectors:\n")
print(rbind("(1,3)"   = gamma1D.true,
            "(2,3)"   = gamma2D.true,
            "(1,2|3)" = gamma12.true))
cat("\nCalibrated marginal intercepts (gammas):\n")
print(unlist(gammas))
cat("\nImplied tau range across (z1, z2) support:\n")
print(do.call(rbind, tau.all))

## ------------------------------------------------------------------
## Save the setting so simulation drivers can `readRDS()` it.
## ------------------------------------------------------------------
saveRDS(
  list(betas = betas, vine_spec = vine_spec,
       copula_pars = copula_pars,
       gammas = gammas,
       c.lwr = c.lwr, c.upr = c.upr,
       tt.all = tt.all,
       tau.all = tau.all),
  file = file.path("settings", "sim1_setting.rds")
)
cat("\nSetting saved to LiDA/simulation/settings/sim1_setting.rds\n")
