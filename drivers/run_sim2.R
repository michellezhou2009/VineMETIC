rm(list = ls())

library(VineSurvIC)
library(VineCopula)
library(survival)
library(trust)

## ------------------------------------------------------------------
## Configuration -- shared across all scenarios in this loop
## ------------------------------------------------------------------
scenarios_overnight <- c(
  "Gumbel_Clayton_Frank",
  "Clayton_Clayton_Gaussian",
  "Clayton_Gumbel_Gaussian"
)
nrep_target <- 500L
nsamp_grid  <- c(500L, 1000L)

seed_master  <- 20240122L
setting_path <- "settings/sim2_setting.rds"
if (!file.exists(setting_path)) {
  source("settings/sim2_setting.R", chdir = FALSE)
}
setting <- readRDS(setting_path)
dgp     <- setting$dgp

set.seed(seed_master)
seeds <- sample.int(.Machine$integer.max, size = nrep_target)

## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
surv_at_grid <- function(cumhaz, tt) {
  tk     <- cumhaz$time
  Lambda <- cumhaz$est
  Lse    <- cumhaz$se
  idx    <- findInterval(tt, tk)
  Ltt    <- ifelse(idx == 0L, 0, Lambda[pmax(idx, 1L)])
  Lse.tt <- ifelse(idx == 0L, 0, Lse[pmax(idx, 1L)])
  Shat   <- exp(-Ltt)
  data.frame(tt = tt, surv = Shat, se = Shat * Lse.tt)
}

fit_one_rep <- function(r, n, seed_r, dgp, fit_spec, scenario, out_root) {
  out_dir  <- file.path(out_root, sprintf("n%d", n))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, sprintf("rep_%04d.rds", r))
  if (file.exists(out_file)) return(list(rep = r, n = n, status = "cached"))

  t0 <- Sys.time()
  set.seed(seed_r)
  dat <- VineSurvIC::gen_metic_data_j3(
    n            = n,
    betas        = dgp$betas,
    vine_spec    = dgp$vine_spec,
    copula_pars  = dgp$copula_pars,
    gammas       = dgp$gammas,
    c.lwr        = dgp$c.lwr,
    c.upr        = dgp$c.upr,
    seed         = seed_r
  )
  cen_rate <- data.frame(
    T1 = 1 - mean(dat$delta1), T2 = 1 - mean(dat$delta2),
    TD = 1 - mean(dat$deltaD), rep = r, nsamp = n)

  fit <- tryCatch(
    VineSurvIC::fit_metic_j3(
      data       = dat,
      times      = c("time1", "time2", "timeD"),
      statuses   = c("delta1", "delta2", "deltaD"),
      Z_formulas = list(~ z1 + z2, ~ z1 + z2, ~ z1 + z2),
      vine_spec  = fit_spec,
      Gfuns      = c("PH", "PH", "PH"),
      pc.method  = "none",
      verbose    = FALSE
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("rep ", r, ":", elapsed, "\n", sep = "")
  if (inherits(fit, "error")) {
    res <- list(rep = r, n = n, status = "error",
                scenario = scenario,
                error = conditionMessage(fit),
                elapsed = elapsed, cen_rate = cen_rate)
  } else {
    s  <- summary(fit)
    pp <- dgp$tt.all$pp
    survT1 <- surv_at_grid(s$Lambda_T1, dgp$tt.all$t1)
    survT1$true <- pp; survT1$rep <- r; survT1$nsamp <- n
    survT2 <- surv_at_grid(s$Lambda_T2, dgp$tt.all$t2)
    survT2$true <- pp; survT2$rep <- r; survT2$nsamp <- n
    survTD <- surv_at_grid(s$Lambda_T3, dgp$tt.all$tD)
    survTD$true <- pp; survTD$rep <- r; survTD$nsamp <- n

    res <- list(
      rep = r, n = n, status = "ok",
      scenario = scenario,
      fit_families = c(fit_spec$j1J$family,
                       fit_spec$j2J$family,
                       fit_spec$j12J$family),
      elapsed = elapsed, cen_rate = cen_rate,
      beta_T1 = s$beta_T1, beta_T2 = s$beta_T2, beta_T3 = s$beta_T3,
      gamma_13 = s$gamma_13, gamma_23 = s$gamma_23,
      gamma_12_3 = s$gamma_12_3,
      survT1 = survT1, survT2 = survT2, survTD = survTD,
      implied_13 = s$implied_13,
      implied_23 = s$implied_23,
      implied_12_3 = s$implied_12_3)
  }
  saveRDS(res, file = out_file)
  res
}

is_scenario_complete <- function(s, n, nrep_target) {
  out_dir <- file.path("raw_results", "sim2", s, sprintf("n%d", n))
  if (!dir.exists(out_dir)) return(FALSE)
  files <- list.files(out_dir, pattern = "^rep_\\d{4}\\.rds$")
  length(files) >= nrep_target
}

## ------------------------------------------------------------------
## Main loop: scenario x n x rep (sequential)
## ------------------------------------------------------------------
cat("\nSim 2 launcher\n")
cat("  scenarios: ", paste(scenarios_overnight, collapse = ", "),
    "\n", sep = "")
cat("  nsamp    : ", paste(nsamp_grid, collapse = ", "), "\n", sep = "")
cat("  nrep     : ", nrep_target, "\n", sep = "")
cat("  start    : ", format(Sys.time()), "\n\n", sep = "")

for (scenario in scenarios_overnight) {
  fit_spec <- setting$fit_specs[[scenario]]
  if (is.null(fit_spec)) {
    cat(">>> Scenario '", scenario, "' not in sim2_setting.rds; skipping.\n",
        sep = "")
    next
  }
  out_root <- file.path("raw_results", "sim2", scenario)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  cat("\n>>> Scenario: ", scenario, "  (", format(Sys.time()), ")\n",
      sep = "")
  cat("    fit families: ",
      fit_spec$j1J$family, " | ",
      fit_spec$j2J$family, " | ",
      fit_spec$j12J$family, "\n", sep = "")

  for (n in nsamp_grid) {
    if (is_scenario_complete(scenario, n, nrep_target)) {
      cat("    n=", n, " already complete; skipping.\n", sep = "")
      next
    }
    cat("    === ", scenario, ", n = ", n, " ===\n", sep = "")
    t0 <- Sys.time()
    results <- vector("list", nrep_target)
    for (r in seq_len(nrep_target)) {
      results[[r]] <- fit_one_rep(r = r, n = n, seed_r = seeds[r],
                                   dgp = dgp, fit_spec = fit_spec,
                                   scenario = scenario, out_root = out_root)
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok    <- sum(vapply(results, function(x) identical(x$status, "ok"),
                        logical(1)), na.rm = TRUE)
    err   <- sum(vapply(results, function(x) identical(x$status, "error"),
                        logical(1)), na.rm = TRUE)
    cache <- sum(vapply(results, function(x) identical(x$status, "cached"),
                        logical(1)), na.rm = TRUE)
    cat("      ok = ", ok, "  cached = ", cache, "  errors = ", err,
        "  wall = ", round(elapsed / 60, 1), " min\n", sep = "")
  }
}

cat("\nSim 2 finished at ", format(Sys.time()), "\n", sep = "")
