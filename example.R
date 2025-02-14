rm(list=ls())
library(pracma)
library(tidyverse)
library(VineCopula)
library(survival)
library(trust)
library(parallel)
library(foreach)
library(doSNOW)

source("fitSPT.R")
source("helpers.R")
source("helpers_bvic.R")
source("helpers_trivic.R")
source("cvine_trivic.R")

ncores = 2
cl = makeCluster(ncores)
registerDoSNOW(cl)

data = readRDS("example.rds")

times = c("time1", "time2", "timeD")
deltas = c("delta1", "delta2", "deltaD")
fmlas = list("T1" = ~ z1 + z2, "T2" = ~ z1 + z2, "D" = ~ z1 + z2)
Gfuns = c("T1" = "PH", "T2" = "PH", "D" = "PH")
copula.fams = matrix(NA, 3, 3); 
copula.fams[1, 2] = "Frank"; 
copula.fams[1, 3] = "Gumbel"; copula.fams[2, 3] = "Clayton"
copula.links = matrix(NA, 3, 3)
copula.links[1, 2] = "identity"; 
copula.links[1, 3] = "log-1"; copula.links[2, 3] = "log"
copula.fmlas = matrix(NA, 3, 3)
copula.fmlas[1, 3] = "~ z1 + z2"; copula.fmlas[2, 3] = "~ z1 + z2"
copula.fmlas[1, 2] = "~ z1 + z2"
system.time({
  out = cvine.trivic(
      data = data, times, deltas, fmlas, Gfuns, copula.fams, 
      copula.links, copula.fmlas, ncores = ncores, 
      pc.method = "foreach")
})
stopCluster(cl)


out$T1$beta
out$T2$beta
out$TD$beta
out$`Cop(1,3)`
out$`Cop(2,3)`
out$`Cop(1,2)`



