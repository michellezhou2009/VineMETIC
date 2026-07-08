## gather/gather_sim1.R
##
## Aggregate the per-rep .rds files produced by drivers/run_sim1.R into
## the two tables reported for Simulation Study I:
##   Table 1 (main text)   -- marginal and copula coefficient estimates
##   Web Table S1 (supp)   -- baseline survival probability estimates
##
## Both tables are written as .csv files with the same rows/columns
## shown in the manuscript and supplementary material (absolute BIAS,
## ESD, ASE, ECP; sample sizes n = 500, 1000, 2000 as separate column
## blocks).
##
## Run from the VineMETIC/ working directory:
##   source("gather/gather_sim1.R")

rm(list = ls())

library(VineSurvIC)

in_root  <- "raw_results/sim1"
out_root <- "gather"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

setting_path <- "settings/sim1_setting.rds"
if (!file.exists(setting_path)) {
  stop("Cannot find ", setting_path, ". Source settings/sim1_setting.R first.")
}
setting <- readRDS(setting_path)

cat("Sim 1 gather\n")
cat("  in_root  : ", in_root, "\n", sep = "")
cat("  out_root : ", out_root, "\n", sep = "")

## ------------------------------------------------------------------
## 1. Read every per-rep .rds file under in_root/n*/rep_*.rds
## ------------------------------------------------------------------
ndirs <- list.files(in_root, pattern = "^n\\d+$", full.names = TRUE)
if (length(ndirs) == 0L) stop("No n* subfolders found under '", in_root, "'.")
ndirs <- ndirs[order(as.integer(sub("^.*?n(\\d+)$", "\\1", ndirs)))]

true_betas <- list(T1 = setting$betas$t1, T2 = setting$betas$t2,
                   TD = setting$betas$tD)
true_gammas <- list(
  "C(1,D)" = setting$copula_pars[["(1,3)"]],
  "C(2,D)" = setting$copula_pars[["(2,3)"]],
  "C(1,2)" = setting$copula_pars[["(1,2)"]]
)

make_para_rows <- function(tab, dist, true_vec, rep, nsamp) {
  if (is.null(tab) || nrow(tab) == 0L) return(NULL)
  data.frame(est = as.numeric(tab$est), se = as.numeric(tab$se),
             dist = dist, para = rownames(tab), true = as.numeric(true_vec),
             rep = rep, nsamp = nsamp, stringsAsFactors = FALSE)
}

para.est.all   <- list()
survT1.est.all <- list()
survT2.est.all <- list()
survTD.est.all <- list()
n_ok <- n_err <- 0L

for (nd in ndirs) {
  files <- list.files(nd, pattern = "^rep_\\d{4}\\.rds$", full.names = TRUE)
  for (f in files) {
    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res)) next
    if (!identical(res$status, "ok")) { n_err <- n_err + 1L; next }
    n_ok <- n_ok + 1L

    para.est.all[[length(para.est.all) + 1L]] <- rbind(
      make_para_rows(res$beta_T1,    "T1",     true_betas$T1, res$rep, res$n),
      make_para_rows(res$beta_T2,    "T2",     true_betas$T2, res$rep, res$n),
      make_para_rows(res$beta_T3,    "TD",     true_betas$TD, res$rep, res$n),
      make_para_rows(res$gamma_13,   "C(1,D)", true_gammas[["C(1,D)"]], res$rep, res$n),
      make_para_rows(res$gamma_23,   "C(2,D)", true_gammas[["C(2,D)"]], res$rep, res$n),
      make_para_rows(res$gamma_12_3, "C(1,2)", true_gammas[["C(1,2)"]], res$rep, res$n)
    )
    if (!is.null(res$survT1)) survT1.est.all[[length(survT1.est.all) + 1L]] <- res$survT1
    if (!is.null(res$survT2)) survT2.est.all[[length(survT2.est.all) + 1L]] <- res$survT2
    if (!is.null(res$survTD)) survTD.est.all[[length(survTD.est.all) + 1L]] <- res$survTD
  }
}

para.est.all   <- do.call(rbind, para.est.all)
survT1.est.all <- do.call(rbind, survT1.est.all)
survT2.est.all <- do.call(rbind, survT2.est.all)
survTD.est.all <- do.call(rbind, survTD.est.all)

cat("Loaded ", n_ok, " ok reps and ", n_err, " errored reps.\n", sep = "")

## ------------------------------------------------------------------
## 2. Summary statistics -- absolute BIAS / ESD / ASE / ECP per cell
## ------------------------------------------------------------------
summarize_para <- function(df) {
  df$lwr <- df$est - 1.96 * df$se
  df$upr <- df$est + 1.96 * df$se
  df$covered <- as.integer(df$lwr <= df$true & df$upr >= df$true)
  grp <- with(df, paste(nsamp, dist, para, sep = "|"))
  rows <- lapply(split(df, grp), function(d) {
    data.frame(nsamp = d$nsamp[1], dist = d$dist[1], para = d$para[1],
               true = d$true[1], BIAS = mean(d$est) - d$true[1],
               ESD = sd(d$est), ASE = mean(d$se), ECP = 100 * mean(d$covered),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summarize_surv <- function(df) {
  df$lwr <- df$surv - 1.96 * df$se
  df$upr <- df$surv + 1.96 * df$se
  df$covered <- as.integer(df$lwr <= df$true & df$upr >= df$true)
  grp <- with(df, paste(nsamp, round(true, 4), sep = "|"))
  rows <- lapply(split(df, grp), function(d) {
    data.frame(nsamp = d$nsamp[1], tt = mean(d$tt), true = d$true[1],
               BIAS = mean(d$surv) - d$true[1], ESD = sd(d$surv),
               ASE = mean(d$se), ECP = 100 * mean(d$covered),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summary.para  <- summarize_para(para.est.all)
summary.surv1 <- summarize_surv(survT1.est.all)
summary.surv2 <- summarize_surv(survT2.est.all)
summary.survD <- summarize_surv(survTD.est.all)

## ------------------------------------------------------------------
## 3. Table 1 (main text): wide CSV, one row per parameter,
##    column blocks n500_*/n1000_*/n2000_*
## ------------------------------------------------------------------
nsamp_vals <- sort(unique(summary.para$nsamp))

para_order <- data.frame(
  dist  = c("T1","T1","T2","T2","TD","TD",
            "C(1,D)","C(1,D)","C(1,D)",
            "C(2,D)","C(2,D)","C(2,D)",
            "C(1,2)","C(1,2)","C(1,2)"),
  para  = c("z1","z2","z1","z2","z1","z2",
            "(Intercept)","z1","z2",
            "(Intercept)","z1","z2",
            "(Intercept)","z1","z2"),
  Parameter = c("beta_11","beta_12","beta_21","beta_22","beta_D1","beta_D2",
                "gamma_1D_0","gamma_1D_1","gamma_1D_2",
                "gamma_2D_0","gamma_2D_1","gamma_2D_2",
                "gamma_12_0","gamma_12_1","gamma_12_2"),
  stringsAsFactors = FALSE
)

build_wide <- function(summary_df, key_df, key_cols) {
  wide <- key_df
  wide$True <- NA_real_
  for (n in nsamp_vals) {
    wide[[paste0("n", n, "_BIAS")]] <- NA_real_
    wide[[paste0("n", n, "_ESD")]]  <- NA_real_
    wide[[paste0("n", n, "_ASE")]]  <- NA_real_
    wide[[paste0("n", n, "_ECP")]]  <- NA_real_
  }
  for (i in seq_len(nrow(key_df))) {
    for (n in nsamp_vals) {
      sel <- TRUE
      for (kc in key_cols) sel <- sel & (summary_df[[kc]] == key_df[[kc]][i])
      sel <- sel & (summary_df$nsamp == n)
      row <- summary_df[sel, ]
      if (nrow(row) == 1L) {
        wide$True[i] <- row$true
        wide[[paste0("n", n, "_BIAS")]][i] <- row$BIAS
        wide[[paste0("n", n, "_ESD")]][i]  <- row$ESD
        wide[[paste0("n", n, "_ASE")]][i]  <- row$ASE
        wide[[paste0("n", n, "_ECP")]][i]  <- row$ECP
      }
    }
  }
  wide
}

tab1_wide <- build_wide(summary.para, para_order, c("dist", "para"))
tab1_wide <- tab1_wide[, c("Parameter", "True",
                           unlist(lapply(nsamp_vals, function(n)
                             paste0("n", n, c("_BIAS","_ESD","_ASE","_ECP")))))]

write.csv(tab1_wide, file.path(out_root, "table1_sim1.csv"),
         row.names = FALSE)
cat("Wrote: ", file.path(out_root, "table1_sim1.csv"), "\n", sep = "")

## ------------------------------------------------------------------
## 4. Web Table S1 (supp): baseline survival probabilities
## ------------------------------------------------------------------
surv_key <- data.frame(
  Event = c("S1","S1","S1","S2","S2","S2","SD","SD","SD"),
  true  = c(0.25, 0.50, 0.75, 0.25, 0.50, 0.75, 0.25, 0.50, 0.75),
  stringsAsFactors = FALSE
)

surv_all <- rbind(
  cbind(Event = "S1", summary.surv1),
  cbind(Event = "S2", summary.surv2),
  cbind(Event = "SD", summary.survD)
)
surv_all <- surv_all[round(surv_all$true, 2) %in% c(0.25, 0.50, 0.75), ]

build_surv_wide <- function() {
  wide <- surv_key
  wide$t <- NA_real_
  for (n in nsamp_vals) {
    wide[[paste0("n", n, "_BIAS")]] <- NA_real_
    wide[[paste0("n", n, "_ESD")]]  <- NA_real_
    wide[[paste0("n", n, "_ASE")]]  <- NA_real_
    wide[[paste0("n", n, "_ECP")]]  <- NA_real_
  }
  for (i in seq_len(nrow(surv_key))) {
    for (n in nsamp_vals) {
      row <- surv_all[surv_all$Event == surv_key$Event[i] &
                      round(surv_all$true, 2) == surv_key$true[i] &
                      surv_all$nsamp == n, ]
      if (nrow(row) == 1L) {
        wide$t[i] <- row$tt
        wide[[paste0("n", n, "_BIAS")]][i] <- row$BIAS
        wide[[paste0("n", n, "_ESD")]][i]  <- row$ESD
        wide[[paste0("n", n, "_ASE")]][i]  <- row$ASE
        wide[[paste0("n", n, "_ECP")]][i]  <- row$ECP
      }
    }
  }
  wide
}

webS1_wide <- build_surv_wide()
webS1_wide <- webS1_wide[, c("Event", "true", "t",
                             unlist(lapply(nsamp_vals, function(n)
                               paste0("n", n, c("_BIAS","_ESD","_ASE","_ECP")))))]
names(webS1_wide)[names(webS1_wide) == "true"] <- "True_S"

write.csv(webS1_wide, file.path(out_root, "webtableS1_sim1_survival.csv"),
         row.names = FALSE)
cat("Wrote: ", file.path(out_root, "webtableS1_sim1_survival.csv"), "\n", sep = "")

cat("\nDone.\n")
