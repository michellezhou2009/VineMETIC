## Simulation Study 3 driver
##
## DGP: symmetric trivariate Clayton copula C(T1, T2, D; phi) with no
## covariates. Two tau values (0.3, 0.7), sample sizes n in (500, 1000),
## 500 replications each.
##
## C-vine decomposition (star with D as root):
##   C(T1, D)      ~ Clayton(phi_true)      -> gamma0_true  = phi_true
##   C(T2, D)      ~ Clayton(phi_true)      -> gamma0_true  = phi_true
##   C(T1, T2 | D) ~ Clayton(phi/(1+phi))   -> gamma12_true = phi_true / (1 + phi_true)
##
## Three methods:
##   A. Li (2020)   -- iterative pseudo-MLE, single shared Clayton parameter
##   B. vine_diff   -- fit_metic_j3() with independent gamma1, gamma2, gamma12
##   C. vine_same   -- fit_trivic_same_j3() with constrained gamma1 = gamma2 = gamma0
rm(list = ls())

library(VineSurvIC)
library(VineCopula)
library(dplyr)
library(copula)
library(survival)
library(trust)
library(pracma)

## Source Li (2020) helpers and tau12 quadrature
source("utils/li2020_fit.R")
source("utils/helpers.R")

# ---- Configuration ----------------------------------------------------------
tau_grid    <- c(0.3, 0.7)
nrep_target <- 5L
# nsamp_grid  <- c(500L, 1000L)
nsamp_grid  <- c(500L)
seed_master <- 20240122L

cat("Sim 3 launcher\n")
cat("  tau_grid    :", paste(tau_grid,    collapse = ", "), "\n")
cat("  nsamp_grid  :", paste(nsamp_grid,  collapse = ", "), "\n")
cat("  nrep_target :", nrep_target, "\n\n")

# ---- True parameter grid for evaluation time points ------------------------
pp.true <- seq(0.05, 0.95, by = 0.05)
tt.all  <- list(
  T1 = qweibull(p = pp.true, shape = 2, scale = 70),
  T2 = qweibull(p = pp.true, shape = 2, scale = 60)
)

# ---- Reproducible per-rep seeds ---------------------------------------------
set.seed(seed_master)
seeds <- sample.int(.Machine$integer.max, size = nrep_target)

# ---- Vine specs (fixed across all reps/scenarios) --------------------------
## vine_diff: independent Clayton for all three edges
vine_spec_diff <- list(
  j1J  = list(family = "Clayton", link = "identity", W_formula = ~ 1),
  j2J  = list(family = "Clayton", link = "identity", W_formula = ~ 1),
  j12J = list(family = "Clayton", link = "identity", W_formula = ~ 1)
)

## vine_same: shared Clayton for C(T1,D) = C(T2,D); Clayton for C(T1,T2|D)
vine_spec_same <- default_vine_same_j3(
  family_jJ   = "Clayton", link_jJ   = "identity",
  family_j12J = "Clayton", link_j12J = "identity"
)

# ---- Helper: survival at evaluation time points from a cumhaz table --------
## Lambda_df has columns 'time' (event grid tk) and 'est' (cumulative hazard).
## Returns S(t) = exp(-Lambda(t)) for each t in tt.
surv_at_tt <- function(Lambda_df, tt) {
  dL <- diff(c(0, Lambda_df$est))
  Lambda.tt <- sum_I(tt, ">=", Lambda_df$time, dL)
  exp(-Lambda.tt)
}

# ---- Per-rep worker ---------------------------------------------------------
fit_one_rep_sim3 <- function(r, n, seed_r,
                              phi_true, gamma0_true, gamma12_true,
                              tau_true, out_root) {

  out_dir  <- file.path(out_root, sprintf("n%d", n))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, sprintf("rep_%04d.rds", r))

  if (file.exists(out_file))
    return(list(rep = r, n = n, status = "cached"))

  t0 <- Sys.time()

  ## ---- Generate data -------------------------------------------------------
  set.seed(seed_r)
  clayton_cop <- copula::claytonCopula(param = phi_true, dim = 3)
  sam   <- copula::rCopula(n = n, clayton_cop)
  T1    <- qweibull(sam[, 1], shape = 2, scale = 70, lower.tail = FALSE)
  T2    <- qweibull(sam[, 2], shape = 2, scale = 60, lower.tail = FALSE)
  D     <- qweibull(sam[, 3], shape = 2, scale = 85, lower.tail = FALSE)
  Ca    <- rexp(n, rate = 1 / 350)
  C     <- pmin(D, Ca)
  U1    <- pmin(T1, C)
  U2    <- pmin(T2, C)
  delta_D <- as.integer(D  <= Ca)
  delta_1 <- as.integer(T1 <= C)
  delta_2 <- as.integer(T2 <= C)
  cen_rate <- data.frame(
    T1 = 1 - mean(delta_1),
    T2 = 1 - mean(delta_2),
    TD = 1 - mean(delta_D)
  )
  dat <- data.frame(U1, U2, C, delta_1, delta_2, delta_D)

  ## ---- Method A: Li (2020) ------------------------------------------------
  Li2020 <- tryCatch({
    res.Li     <- est_cop_par(copula = "clayton", data = dat)
    phi_est    <- as.numeric(res.Li[1L])
    tau_est    <- BiCopPar2Tau(family = 3L, par = phi_est)

    mydat1 <- npest.star_12(dat)$X
    S.X1   <- data.frame(
      time = dat$U1,
      surv = g.fn(mydat1$S1.star.hat, mydat1$SD.U1.hat, phi_est)
    ) %>% arrange(time)
    S.X2   <- data.frame(
      time = dat$U2,
      surv = g.fn(mydat1$S2.star.hat, mydat1$SD.U2.hat, phi_est)
    ) %>% arrange(time)
    surv_T1 <- c(1, S.X1$surv)[sum_I(tt.all[[1]], ">=", S.X1$time) + 1L]
    surv_T2 <- c(1, S.X2$surv)[sum_I(tt.all[[2]], ">=", S.X2$time) + 1L]

    list(phi = phi_est, tau = tau_est,
         surv_T1 = surv_T1, surv_T2 = surv_T2)
  }, error = function(e) {
    list(phi = NA_real_, tau = NA_real_,
         surv_T1 = rep(NA_real_, length(pp.true)),
         surv_T2 = rep(NA_real_, length(pp.true)))
  })

  ## ---- Method B: vine_diff  -----------------------------------------------
  vc_diff <- tryCatch({
    fit.diff <- fit_metic_j3(
      data       = dat,
      times      = c("U1", "U2", "C"),
      statuses   = c("delta_1", "delta_2", "delta_D"),
      Z_formulas = list(~ 1, ~ 1, ~ 1),
      vine_spec  = vine_spec_diff,
      Gfuns      = c("PH", "PH", "PH"),
      pc.method  = "none",
      verbose    = FALSE
    )
    gamma1  <- fit.diff$summary[["Cop(1,3)"]]$est[1L]
    gamma2  <- fit.diff$summary[["Cop(2,3)"]]$est[1L]
    gamma12 <- fit.diff$summary[["Cop(1,2)"]]$est[1L]
    list(
      gamma1  = gamma1,
      gamma2  = gamma2,
      gamma12 = gamma12,
      tau1D   = BiCopPar2Tau(family = 3L, par = gamma1),
      tau2D   = BiCopPar2Tau(family = 3L, par = gamma2),
      surv_T1 = surv_at_tt(fit.diff$summary$T1$cumhaz, tt.all[[1]]),
      surv_T2 = surv_at_tt(fit.diff$summary$T2$cumhaz, tt.all[[2]])
    )
  }, error = function(e) {
    list(gamma1 = NA_real_, gamma2 = NA_real_, gamma12 = NA_real_,
         tau1D = NA_real_, tau2D = NA_real_,
         surv_T1 = rep(NA_real_, length(pp.true)),
         surv_T2 = rep(NA_real_, length(pp.true)))
  })

  ## ---- Method C: vine_same  -----------------------------------------------
  vc_same <- tryCatch({
    fit.same <- fit_trivic_same_j3(
      data       = dat,
      times      = c("U1", "U2", "C"),
      statuses   = c("delta_1", "delta_2", "delta_D"),
      Z_formulas = list(~ 1, ~ 1, ~ 1),
      vine_spec  = vine_spec_same,
      Gfuns      = c("PH", "PH", "PH"),
      pc.method  = "none",
      verbose    = FALSE
    )
    gamma0  <- fit.same$gamma.est[1L]
    gamma12 <- fit.same$gamma12.est[1L]
    list(
      gamma0  = gamma0,
      gamma12 = gamma12,
      tau0D   = BiCopPar2Tau(family = 3L, par = gamma0),
      surv_T1 = surv_at_tt(fit.same$summary$Lambda1, tt.all[[1]]),
      surv_T2 = surv_at_tt(fit.same$summary$Lambda2, tt.all[[2]])
    )
  }, error = function(e) {
    list(gamma0 = NA_real_, gamma12 = NA_real_, tau0D = NA_real_,
         surv_T1 = rep(NA_real_, length(pp.true)),
         surv_T2 = rep(NA_real_, length(pp.true)))
  })

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("rep %d: %.0f s\n", r, elapsed))
  res <- list(
    rep          = r,
    n            = n,
    tau_true     = tau_true,
    phi_true     = phi_true,
    gamma0_true  = gamma0_true,
    gamma12_true = gamma12_true,
    status       = "ok",
    elapsed      = elapsed,
    cen_rate     = cen_rate,
    Li2020       = Li2020,
    vc_diff      = vc_diff,
    vc_same      = vc_same
  )

  saveRDS(res, file = out_file)
  res
}

# ---- Main loop: (tau, n, rep), sequential -----------------------------------
for (tau_true in tau_grid) {
  phi_true     <- BiCopTau2Par(family = 3L, tau = tau_true)
  gamma0_true  <- phi_true
  gamma12_true <- phi_true / (1 + phi_true)
  out_root     <- file.path("raw_results", "sim3", sprintf("tau%g", tau_true))

  cat(sprintf("\n=== tau_true = %.1f  phi_true = %.4f ===\n",
              tau_true, phi_true))

  for (n in nsamp_grid) {
    cat(sprintf("  n = %d\n", n))

    results <- vector("list", nrep_target)
    for (r in seq_len(nrep_target)) {
      results[[r]] <- fit_one_rep_sim3(
        r            = r,
        n            = n,
        seed_r       = seeds[r],
        phi_true     = phi_true,
        gamma0_true  = gamma0_true,
        gamma12_true = gamma12_true,
        tau_true     = tau_true,
        out_root     = out_root
      )
    }

    ok    <- sum(vapply(results, function(x) identical(x$status, "ok"),     logical(1L)))
    cache <- sum(vapply(results, function(x) identical(x$status, "cached"), logical(1L)))
    err   <- sum(vapply(results, inherits, logical(1L), "error"))
    cat(sprintf("    ok = %d  cached = %d  errors = %d\n", ok, cache, err))
  }
}

cat("\nSim 3 finished. Results in: raw_results/sim3/\n")
