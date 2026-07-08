## gather/gather_sim3.R
##
## Aggregates per-replication results for Simulation Study III and
## produces the Web Table (Sim III, supp): Bias / ESD / RMSE for the
## nonterminal-terminal dependence (tau_1,3 = tau_2,3), the between-
## nonterminal dependence (tau_1,2), and the nonterminal marginal
## survival probabilities S_1(t_M), S_2(t_M), compared across the
## tri-Clayton (Li, 2020), Vine, and pVine methods.
##
## Expects that drivers/compute_tau12_sim3.R has already added a
## tau12 field to each per-rep .rds produced by drivers/run_sim3.R.
## Missing tau12 fields are silently carried as NA.
##
## Output:
##   gather/table3_sim3.csv
##
## Run from the VineMETIC/ working directory:
##   source("gather/gather_sim3.R")

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

tau_grid   <- c(0.3, 0.7)
nsamp_grid <- c(500L, 1000L)
out_root   <- "gather"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

pp.true <- seq(0.05, 0.95, by = 0.05)

## ------------------------------------------------------------------
## 1. Read every per-rep .rds into a long data frame
## ------------------------------------------------------------------
read_one_rep <- function(file) {
  obj <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(obj) || !identical(obj$status, "ok")) return(NULL)

  tau_true     <- obj$tau_true
  phi_true     <- obj$phi_true
  gamma0_true  <- obj$gamma0_true
  gamma12_true <- obj$gamma12_true

  rows <- list()

  li <- obj$Li2020
  rows$Li2020_params <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "Li2020",
    para = c("phi", "tau"), est = c(li$phi, li$tau),
    true = c(phi_true, tau_true), stringsAsFactors = FALSE)
  rows$Li2020_surv <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "Li2020",
    para = c(paste0("surv_T1_", round(pp.true, 2)),
             paste0("surv_T2_", round(pp.true, 2))),
    est = c(li$surv_T1, li$surv_T2), true = c(1 - pp.true, 1 - pp.true),
    stringsAsFactors = FALSE)

  vd <- obj$vc_diff
  rows$diff_params <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "vine_diff",
    para = c("gamma0_T1", "gamma0_T2", "gamma12", "tau1D", "tau2D", "tau12"),
    est = c(vd$gamma1, vd$gamma2, vd$gamma12, vd$tau1D, vd$tau2D,
            if (is.null(vd$tau12)) NA_real_ else vd$tau12),
    true = c(gamma0_true, gamma0_true, gamma12_true, tau_true, tau_true, tau_true),
    stringsAsFactors = FALSE)
  rows$diff_surv <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "vine_diff",
    para = c(paste0("surv_T1_", round(pp.true, 2)),
             paste0("surv_T2_", round(pp.true, 2))),
    est = c(vd$surv_T1, vd$surv_T2), true = c(1 - pp.true, 1 - pp.true),
    stringsAsFactors = FALSE)

  vs <- obj$vc_same
  rows$same_params <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "vine_same",
    para = c("gamma0", "gamma12", "tau0D", "tau12"),
    est = c(vs$gamma0, vs$gamma12, vs$tau0D,
            if (is.null(vs$tau12)) NA_real_ else vs$tau12),
    true = c(gamma0_true, gamma12_true, tau_true, tau_true),
    stringsAsFactors = FALSE)
  rows$same_surv <- data.frame(
    rep = obj$rep, n = obj$n, tau_true = tau_true, method = "vine_same",
    para = c(paste0("surv_T1_", round(pp.true, 2)),
             paste0("surv_T2_", round(pp.true, 2))),
    est = c(vs$surv_T1, vs$surv_T2), true = c(1 - pp.true, 1 - pp.true),
    stringsAsFactors = FALSE)

  dplyr::bind_rows(rows)
}

all_reps <- map_dfr(tau_grid, function(tau) {
  map_dfr(nsamp_grid, function(n) {
    dir_path <- file.path("raw_results", "sim3", sprintf("tau%g", tau), sprintf("n%d", n))
    files <- list.files(dir_path, pattern = "^rep_\\d+\\.rds$", full.names = TRUE)
    cat(sprintf("tau=%.1f  n=%d: %d files\n", tau, n, length(files)))
    map_dfr(files, read_one_rep)
  })
})

tau12_complete <- mean(!is.na(all_reps$est[all_reps$para == "tau12"]))
cat(sprintf("tau12 computed in %.0f%% of rows.\n", 100 * tau12_complete))

## ------------------------------------------------------------------
## 2. Summary: bias / ESD / RMSE per (method, para, tau_true, n)
## ------------------------------------------------------------------
summary3 <- all_reps %>%
  filter(!grepl("^surv_", para)) %>%
  group_by(method, para, tau_true, n) %>%
  summarise(nrep = sum(!is.na(est)), bias = mean(est - true, na.rm = TRUE),
           esd = sd(est, na.rm = TRUE),
           rmse = sqrt(mean((est - true)^2, na.rm = TRUE)), .groups = "drop")

summ_surv <- all_reps %>%
  filter(para %in% c("surv_T1_0.5", "surv_T2_0.5")) %>%
  group_by(method, para, tau_true, n) %>%
  summarise(nrep = sum(!is.na(est)), bias = mean(est - true, na.rm = TRUE),
           esd = sd(est, na.rm = TRUE),
           rmse = sqrt(mean((est - true)^2, na.rm = TRUE)), .groups = "drop")

## ------------------------------------------------------------------
## 3. Table (Sim III, supp): wide CSV matching the manuscript layout
## ------------------------------------------------------------------
nsamp_use <- sort(intersect(nsamp_grid, unique(all_reps$n)))

pull_cell <- function(df, method_val, para_val, tau_val, n_val) {
  row <- df %>% filter(method == method_val, para == para_val,
                       tau_true == tau_val, n == n_val)
  if (nrow(row) == 0L) return(c(Bias = NA_real_, ESD = NA_real_, RMSE = NA_real_))
  c(Bias = row$bias[1], ESD = row$esd[1], RMSE = row$rmse[1])
}

build_row <- function(tau_val, parameter_label, method_label, df_src, para_val, method_val) {
  row <- data.frame(Tau_true = tau_val, Parameter = parameter_label,
                    Method = method_label, stringsAsFactors = FALSE)
  for (n in nsamp_use) {
    cell <- pull_cell(df_src, method_val, para_val, tau_val, n)
    row[[paste0("n", n, "_Bias")]] <- cell["Bias"]
    row[[paste0("n", n, "_ESD")]]  <- cell["ESD"]
    row[[paste0("n", n, "_RMSE")]] <- cell["RMSE"]
  }
  row
}

rows <- list()
for (tau_val in tau_grid) {
  tau_lab <- formatC(tau_val, format = "f", digits = 1)

  block1 <- sprintf("tau_1,3=tau_2,3=%s", tau_lab)
  rows[[length(rows) + 1L]] <- build_row(tau_val, block1, "tri-Clayton", summary3, "tau",   "Li2020")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block1, "Vine(1,3)",  summary3, "tau1D", "vine_diff")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block1, "Vine(2,3)",  summary3, "tau2D", "vine_diff")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block1, "pVine",     summary3, "tau0D", "vine_same")

  block2 <- sprintf("tau_1,2=%s", tau_lab)
  rows[[length(rows) + 1L]] <- build_row(tau_val, block2, "tri-Clayton", summary3, "tau",   "Li2020")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block2, "Vine",       summary3, "tau12", "vine_diff")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block2, "pVine",      summary3, "tau12", "vine_same")

  block3 <- "S_1(t_M)=0.5"
  rows[[length(rows) + 1L]] <- build_row(tau_val, block3, "tri-Clayton", summ_surv, "surv_T1_0.5", "Li2020")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block3, "Vine",       summ_surv, "surv_T1_0.5", "vine_diff")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block3, "pVine",      summ_surv, "surv_T1_0.5", "vine_same")

  block4 <- "S_2(t_M)=0.5"
  rows[[length(rows) + 1L]] <- build_row(tau_val, block4, "tri-Clayton", summ_surv, "surv_T2_0.5", "Li2020")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block4, "Vine",       summ_surv, "surv_T2_0.5", "vine_diff")
  rows[[length(rows) + 1L]] <- build_row(tau_val, block4, "pVine",      summ_surv, "surv_T2_0.5", "vine_same")
}

tab3 <- do.call(rbind, rows)

write.csv(tab3, file.path(out_root, "table3_sim3.csv"), row.names = FALSE)
cat("Wrote: ", file.path(out_root, "table3_sim3.csv"), "\n", sep = "")

cat("\nDone.\n")
