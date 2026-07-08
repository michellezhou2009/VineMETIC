## gather/gather_sim2.R
##
## Aggregate Simulation Study II (misspecification) together with the
## correct-specification results from Simulation Study I, and produce:
##
##   Table 2 (main text)  -- marginal regression coefficient Bias (ESD)
##                            under correct specification and three
##                            copula-misspecification scenarios
##                            -> gather/table2_sim2.csv
##
##   Figure 2 (main text) and Web Figures S2-S4 (supp) -- average bias
##   in the implied, covariate-dependent Kendall's tau for each vine
##   edge, under correct and incorrect copula specifications
##                            -> figures/figure2_tau_bias_C12_given_D.*
##                            -> figures/supp_tau_bias_C13.*
##                            -> figures/supp_tau_bias_C23.*
##                            -> figures/supp_tau_bias_C12_given_D.*
##
## Run from the VineMETIC/ working directory:
##   source("gather/gather_sim2.R")

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(VineCopula)
  library(VineSurvIC)
})

out_root <- "gather"
fig_dir  <- "figures"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

sim1_setting_path <- "settings/sim1_setting.rds"
sim2_setting_path <- "settings/sim2_setting.rds"
if (!file.exists(sim1_setting_path)) stop("Cannot find ", sim1_setting_path)
if (!file.exists(sim2_setting_path)) stop("Cannot find ", sim2_setting_path)

sim1_setting <- readRDS(sim1_setting_path)
sim2_setting <- readRDS(sim2_setting_path)
dgp <- sim2_setting$dgp   ## true DGP, shared by Sim I and Sim II

z1_grid <- seq(-0.5, 0.5, length.out = 101)
z2_grid <- c(0, 1)

cat("Sim 2 gather\n")

## ------------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

make_para_rows <- function(tab, dist, true_vec, rep, nsamp, scenario) {
  if (is.null(tab) || nrow(tab) == 0L) return(NULL)
  data.frame(est = as.numeric(tab$est), se = as.numeric(tab$se),
             dist = dist, para = rownames(tab), true = as.numeric(true_vec),
             rep = rep, nsamp = nsamp, scenario = scenario,
             stringsAsFactors = FALSE)
}

## eta -> copula parameter -> Kendall's tau
eta_to_tau <- function(eta, family, link = NULL) {
  rec <- VineSurvIC:::bicop_get(family)
  link_name <- link %||% rec$default_link
  alpha <- switch(link_name,
    "log"      = exp(eta),
    "log-1"    = 1 + exp(eta),
    "fisher_z" = tanh(eta),
    "identity" = eta,
    { link_obj <- VineSurvIC:::link_funs(link_name); link_obj$h.fun(eta) }
  )
  VineCopula::BiCopPar2Tau(rec$family_code, alpha)
}

gamma_to_tau_grid <- function(gamma, family, link, edge, scenario, rep, nsamp) {
  expand.grid(z1 = z1_grid, z2 = z2_grid) |>
    mutate(eta = gamma[1] + gamma[2] * z1 + gamma[3] * z2,
           tau = eta_to_tau(eta, family = family, link = link),
           edge = edge, scenario = scenario, rep = rep, nsamp = nsamp)
}

## ------------------------------------------------------------------
## 1. Read Simulation Study I (Correct specification) per-rep files
## ------------------------------------------------------------------
true_betas <- list(T1 = dgp$betas$t1, T2 = dgp$betas$t2, TD = dgp$betas$tD)
true_gammas <- list("C(1,D)" = dgp$copula_pars[["(1,3)"]],
                    "C(2,D)" = dgp$copula_pars[["(2,3)"]],
                    "C(1,2|D)" = dgp$copula_pars[["(1,2)"]])

read_sim1 <- function(in_root = "raw_results/sim1") {
  ndirs <- list.files(in_root, pattern = "^n\\d+$", full.names = TRUE)
  rows <- list()
  for (nd in ndirs) {
    files <- list.files(nd, pattern = "^rep_\\d{4}\\.rds$", full.names = TRUE)
    for (f in files) {
      res <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(res) || !identical(res$status, "ok")) next
      rows[[length(rows) + 1L]] <- rbind(
        make_para_rows(res$beta_T1,    "T1",       true_betas$T1,             res$rep, res$n, "Correct specification"),
        make_para_rows(res$beta_T2,    "T2",       true_betas$T2,             res$rep, res$n, "Correct specification"),
        make_para_rows(res$beta_T3,    "TD",       true_betas$TD,             res$rep, res$n, "Correct specification"),
        make_para_rows(res$gamma_13,   "C(1,D)",   true_gammas[["C(1,D)"]],   res$rep, res$n, "Correct specification"),
        make_para_rows(res$gamma_23,   "C(2,D)",   true_gammas[["C(2,D)"]],   res$rep, res$n, "Correct specification"),
        make_para_rows(res$gamma_12_3, "C(1,2|D)", true_gammas[["C(1,2|D)"]], res$rep, res$n, "Correct specification")
      )
    }
  }
  do.call(rbind, rows)
}

sim1_para <- read_sim1()
cat("Sim 1 (correct): ", nrow(sim1_para), " parameter rows loaded.\n", sep = "")

## ------------------------------------------------------------------
## 2. Read Simulation Study II (misspecified scenarios) per-rep files
## ------------------------------------------------------------------
read_sim2 <- function(in_root = "raw_results/sim2") {
  scenario_dirs <- list.files(in_root, full.names = TRUE)
  scenario_dirs <- scenario_dirs[file.info(scenario_dirs)$isdir]
  rows <- list()
  for (sd in scenario_dirs) {
    scenario <- basename(sd)
    ndirs <- list.files(sd, pattern = "^n\\d+$", full.names = TRUE)
    for (nd in ndirs) {
      files <- list.files(nd, pattern = "^rep_\\d{4}\\.rds$", full.names = TRUE)
      for (f in files) {
        res <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(res) || !identical(res$status, "ok")) next
        rows[[length(rows) + 1L]] <- rbind(
          make_para_rows(res$beta_T1,    "T1",       true_betas$T1,             res$rep, res$n, scenario),
          make_para_rows(res$beta_T2,    "T2",       true_betas$T2,             res$rep, res$n, scenario),
          make_para_rows(res$beta_T3,    "TD",       true_betas$TD,             res$rep, res$n, scenario),
          make_para_rows(res$gamma_13,   "C(1,D)",   true_gammas[["C(1,D)"]],   res$rep, res$n, scenario),
          make_para_rows(res$gamma_23,   "C(2,D)",   true_gammas[["C(2,D)"]],   res$rep, res$n, scenario),
          make_para_rows(res$gamma_12_3, "C(1,2|D)", true_gammas[["C(1,2|D)"]], res$rep, res$n, scenario)
        )
      }
    }
  }
  do.call(rbind, rows)
}

sim2_para <- read_sim2()
cat("Sim 2 (misspecified): ", nrow(sim2_para), " parameter rows loaded.\n", sep = "")

para_all <- rbind(sim1_para, sim2_para)

## ------------------------------------------------------------------
## 3. Summary statistics: BIAS / ESD per (scenario, nsamp, dist, para)
## ------------------------------------------------------------------
summary.para <- para_all |>
  group_by(scenario, nsamp, dist, para) |>
  summarize(true = first(true), avg = mean(est), BIAS = avg - true,
           ESD = sd(est), .groups = "drop")

## ------------------------------------------------------------------
## 4. Table 2 (main text): marginal Bias (ESD) by scenario, wide CSV
## ------------------------------------------------------------------
scenario_order <- c("Correct specification", "Gumbel_Clayton_Frank",
                    "Clayton_Clayton_Gaussian", "Clayton_Gumbel_Gaussian")
scenario_labels <- c("Correct specification" = "Correct",
                     "Gumbel_Clayton_Frank" = "GuClFr",
                     "Clayton_Clayton_Gaussian" = "ClClGa",
                     "Clayton_Gumbel_Gaussian" = "ClGuGa")
available_scenarios <- intersect(scenario_order, unique(summary.para$scenario))
scen_labs <- scenario_labels[available_scenarios]

para_order <- data.frame(
  dist = c("T1","T1","T2","T2","TD","TD"),
  para = c("z1","z2","z1","z2","z1","z2"),
  Parameter = c("beta_11","beta_12","beta_21","beta_22","beta_D1","beta_D2"),
  stringsAsFactors = FALSE
)

## Restrict to sample sizes that exist for the misspecified (Sim II)
## scenarios -- Sim II never ran n=2000, so a "Correct specification"-only
## n=2000 block (with the other scenario columns empty) is not useful here.
nsamp_vals <- sort(unique(sim2_para$nsamp))

tab2_rows <- list()
for (n in nsamp_vals) {
  for (i in seq_len(nrow(para_order))) {
    row <- data.frame(n = n, Parameter = para_order$Parameter[i],
                      True = NA_real_, stringsAsFactors = FALSE)
    for (sc in available_scenarios) {
      lab <- scen_labs[[sc]]
      cell <- summary.para[summary.para$scenario == sc &
                           summary.para$nsamp == n &
                           summary.para$dist == para_order$dist[i] &
                           summary.para$para == para_order$para[i], ]
      if (nrow(cell) == 1L) {
        row$True <- cell$true
        row[[paste0(lab, "_Bias")]] <- cell$BIAS
        row[[paste0(lab, "_ESD")]]  <- cell$ESD
      } else {
        row[[paste0(lab, "_Bias")]] <- NA_real_
        row[[paste0(lab, "_ESD")]]  <- NA_real_
      }
    }
    tab2_rows[[length(tab2_rows) + 1L]] <- row
  }
}
tab2 <- do.call(rbind, tab2_rows)
col_order <- c("n", "Parameter", "True",
              unlist(lapply(scen_labs, function(l) paste0(l, c("_Bias","_ESD")))))
tab2 <- tab2[, col_order]

write.csv(tab2, file.path(out_root, "table2_sim2.csv"), row.names = FALSE)
cat("Wrote: ", file.path(out_root, "table2_sim2.csv"), "\n", sep = "")

## ------------------------------------------------------------------
## 5. Tau curves for Figure 2 and Web Figures S2-S4
## ------------------------------------------------------------------
get_family_link <- function(edge_name, scenario) {
  spec <- if (identical(scenario, "Correct specification")) dgp$vine_spec
          else sim2_setting$fit_specs[[scenario]]
  switch(edge_name,
        "C(1,D)"   = list(family = spec$j1J$family,  link = spec$j1J$link),
        "C(2,D)"   = list(family = spec$j2J$family,  link = spec$j2J$link),
        "C(1,2|D)" = list(family = spec$j12J$family, link = spec$j12J$link))
}

edges <- c("C(1,D)", "C(2,D)", "C(1,2|D)")

## True tau grid (one per edge, shared across nsamp)
true_tau <- bind_rows(lapply(nsamp_vals, function(n) {
  bind_rows(lapply(edges, function(e) {
    fl <- get_family_link(e, "Correct specification")
    gamma_to_tau_grid(true_gammas[[e]], fl$family, fl$link, e, "Truth", 0, n)
  }))
})) |> select(nsamp, edge, z1, z2, tau_true = tau)

## Per-rep estimated tau grid, all scenarios and edges
## (restrict to nsamp_vals so Sim I's n=2000 "Correct specification" reps
## -- which Sim II never ran -- don't leak into the tau-bias figures)
tau_rep <- para_all |>
  filter(dist %in% edges, nsamp %in% nsamp_vals) |>
  group_by(scenario, dist, rep, nsamp) |>
  group_modify(function(d, key) {
    d <- d[match(c("(Intercept)", "z1", "z2"), d$para), ]
    gamma <- d$est
    if (any(is.na(gamma))) return(tibble())
    fl <- get_family_link(key$dist, key$scenario)
    gamma_to_tau_grid(gamma, fl$family, fl$link, key$dist,
                      key$scenario, key$rep, key$nsamp) |>
      select(-edge, -scenario, -rep, -nsamp)
  }) |>
  ungroup() |>
  rename(edge = dist)

tau_bias_summary <- tau_rep |>
  left_join(true_tau, by = c("nsamp", "edge", "z1", "z2")) |>
  mutate(bias = tau - tau_true) |>
  group_by(scenario, nsamp, edge, z1, z2) |>
  summarize(mean_bias = mean(bias, na.rm = TRUE), emp_sd = sd(bias, na.rm = TRUE),
           lower = mean_bias - 1.96 * emp_sd, upper = mean_bias + 1.96 * emp_sd,
           .groups = "drop") |>
  filter(scenario %in% available_scenarios) |>
  mutate(scenario = factor(scenario, levels = available_scenarios),
        scenario_lab = factor(scenario_labels[as.character(scenario)],
                              levels = unname(scen_labs)),
        z2_lab = factor(z2, levels = c(0, 1), labels = c("Z[2] == 0", "Z[2] == 1")))

## ------------------------------------------------------------------
## 6. Figures
## ------------------------------------------------------------------
scenario_colors <- c("Correct" = "black", "GuClFr" = "red",
                     "ClClGa" = "blue", "ClGuGa" = "#009E73")
scenario_linetypes <- c("Correct" = "solid", "GuClFr" = "dashed",
                        "ClClGa" = "dotted", "ClGuGa" = "dotdash")

plot_tau_bias <- function(edge_name, show_bands) {
  dat <- tau_bias_summary |> filter(edge == edge_name)
  multi_n <- length(unique(dat$nsamp)) > 1L
  if (multi_n) {
    n_levels <- paste0("n == ", sort(unique(dat$nsamp)))
    dat <- dat |> mutate(n_lab = factor(paste0("n == ", nsamp), levels = n_levels))
    facet_spec <- facet_grid(z2_lab ~ n_lab, labeller = label_parsed)
  } else {
    facet_spec <- facet_wrap(~ z2_lab, nrow = 1, labeller = label_parsed)
  }
  p <- ggplot(dat, aes(x = z1, y = mean_bias, color = scenario_lab,
                       linetype = scenario_lab, fill = scenario_lab,
                       group = scenario_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
    facet_spec +
    scale_color_manual(values = scenario_colors, name = "Fitted copula families") +
    scale_linetype_manual(values = scenario_linetypes, name = "Fitted copula families") +
    scale_fill_manual(values = scenario_colors, guide = "none") +
    labs(x = expression(Z[1]), y = "Average bias in Kendall's tau") +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
         strip.text = element_text(size = 15, face = "bold"),
         axis.title = element_text(size = 15), axis.text = element_text(size = 15),
         legend.title = element_text(size = 15), legend.text = element_text(size = 15))
  if (show_bands) {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.08,
                         color = NA, show.legend = FALSE)
  }
  p + geom_line(linewidth = 1, alpha = 1)
}

fig_w <- 6.0 * length(nsamp_vals)
fig_h <- 3.8 * length(z2_grid)

save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = fig_w, height = fig_h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = fig_w, height = fig_h, dpi = 300)
  cat("Saved figure: ", name, "\n", sep = "")
}

save_fig(plot_tau_bias("C(1,2|D)", show_bands = FALSE) +
        labs(y = expression("Bias in "~tau["1,2|3"])),
        "figure2_tau_bias_C12_given_D")
save_fig(plot_tau_bias("C(1,2|D)", show_bands = TRUE) +
        labs(y = expression("Bias in "~tau["1,2|3"])),
        "supp_tau_bias_C12_given_D")
save_fig(plot_tau_bias("C(1,D)", show_bands = TRUE) +
        labs(y = expression("Bias in "~tau["1,3"])),
        "supp_tau_bias_C13")
save_fig(plot_tau_bias("C(2,D)", show_bands = TRUE) +
        labs(y = expression("Bias in "~tau["2,3"])),
        "supp_tau_bias_C23")

cat("\nDone.\n")
