## compute_tau12_sim3.R
##
## Post-processing step for Simulation Study 3.
## Reads each per-rep RDS produced by run_sim3.R, computes tau12 for
## vine_diff and vine_same via tau12.est.fun() (quadrature, ~45 s each),
## and saves the result back to the same file.
##
## Run AFTER run_sim3.R has finished.
##
## Run from the VineMETIC/ working directory:
##   source("drivers/compute_tau12_sim3.R")

rm(list = ls())

library(VineSurvIC)
library(VineCopula)

source("utils/helpers.R")   # tau12.est.fun

# ---- Configuration ----------------------------------------------------------
tau_grid   <- c(0.3, 0.7)
nsamp_grid <- c(500L, 1000L)

cat("compute_tau12_sim3.R\n")
cat("  tau_grid   :", paste(tau_grid,   collapse = ", "), "\n")
cat("  nsamp_grid :", paste(nsamp_grid, collapse = ", "), "\n\n")

# ---- Worker: compute tau12 for one RDS file --------------------------------
process_one <- function(file) {
  res <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(res)) return(list(file = file, status = "read_error"))
  if (!identical(res$status, "ok"))
    return(list(file = file, status = paste0("skipped:", res$status)))

  ## Skip if already processed
  if (!is.null(res$vc_diff$tau12) && !is.null(res$vc_same$tau12))
    return(list(file = file, status = "cached"))

  ## ---- tau12 for vine_diff ------------------------------------------------
  ## tau12 = Kendall's tau of the bivariate margin C(T1,T2) implied by the
  ## C-vine with C(T1,D) = Clayton(gamma1), C(T2,D) = Clayton(gamma2),
  ## C(T1,T2|D) = Clayton(gamma12). Requires numerical quadrature (~45 s).
  tau12_diff <- tryCatch(
    tau12.est.fun(
      par1    = res$vc_diff$gamma1,
      par2    = res$vc_diff$gamma2,
      par12   = res$vc_diff$gamma12,
      index1  = 3L, index2 = 3L, index12 = 3L,
      pc.method = "none"
    ),
    error = function(e) NA_real_
  )

  ## ---- tau12 for vine_same ------------------------------------------------
  ## Same vine structure but gamma1 = gamma2 = gamma0 (constrained).
  tau12_same <- tryCatch(
    tau12.est.fun(
      par1    = res$vc_same$gamma0,
      par2    = res$vc_same$gamma0,
      par12   = res$vc_same$gamma12,
      index1  = 3L, index2 = 3L, index12 = 3L,
      pc.method = "none"
    ),
    error = function(e) NA_real_
  )

  res$vc_diff$tau12 <- tau12_diff
  res$vc_same$tau12 <- tau12_same

  saveRDS(res, file = file)
  list(file = file, status = "ok",
       tau12_diff = tau12_diff, tau12_same = tau12_same)
}

# ---- Process all files, sequentially ---------------------------------------
for (tau_true in rev(tau_grid)) {
  for (n in rev(nsamp_grid)) {
    dir_path <- file.path("raw_results", "sim3",
                          sprintf("tau%g", tau_true),
                          sprintf("n%d", n))
    files <- list.files(dir_path, pattern = "^rep_\\d+\\.rds$",
                        full.names = TRUE)
    if (length(files) == 0L) {
      message(sprintf("tau=%.1f n=%d: no files found, skipping.", tau_true, n))
      next
    }

    ## How many still need tau12?
    need <- sum(vapply(files, function(f) {
      obj <- tryCatch(readRDS(f), error = function(e) NULL)
      is.null(obj$vc_diff$tau12) || is.null(obj$vc_same$tau12)
    }, logical(1L)))

    cat(sprintf("tau=%.1f  n=%4d: %d files, %d need tau12\n",
                tau_true, n, length(files), need))

    if (need == 0L) next

    results <- vector("list", length(files))
    for (i in seq_along(files)) {
      results[[i]] <- process_one(files[i])
    }

    ok    <- sum(vapply(results, function(x) identical(x$status, "ok"),     logical(1L)))
    cache <- sum(vapply(results, function(x) identical(x$status, "cached"), logical(1L)))
    err   <- sum(vapply(results, function(x) grepl("error", x$status),      logical(1L)))
    cat(sprintf("  done: ok=%d  cached=%d  errors=%d\n", ok, cache, err))
  }
}

cat("\ncompute_tau12_sim3.R finished.\n")
cat("Run gather/gather_sim3.R to aggregate results.\n")
