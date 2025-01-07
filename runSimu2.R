rm(list=ls())
library(bvic)
library(numDeriv)
library(tidyverse)
library(VineCopula)
library(copula)
library(survival)
library(trust)
library(parallel)
library(foreach)
library(doSNOW)

source("fitSPT.R")
source("helpers.R")
source("helpers_bvic.R")
source("helpers_trivic.R")
source("helpers_Li.R")

seed1 = 20240122
set.seed(seed1)
nrep = 500; seed.all = round(runif(nrep, 10, 1000000))
file.dir = "simu_raw/simu2/"

pp.true = seq(0.01, 0.99, by = 0.01)
tt.all = list(
  T1 = qweibull(p = pp.true, shape = 2, scale = 70),
  T2 = qweibull(p = pp.true, shape = 2, scale = 60)
)
tt.dat = rbind(
  data.frame(time = tt.all[[1]], true = 1 - pp.true, para = "surv.T1"),
  data.frame(time = tt.all[[2]], true = 1 - pp.true, para = "surv.T2")
)

tau12.true = tau0.true = 0.7
copula.index = 3
theta12 = BiCopTau2Par(family = copula.index, tau = tau12.true)
theta0 = BiCopTau2Par(family = copula.index, tau = tau0.true)
gamma0.true = theta0
gamma12.true = gamma0.true / (1 + gamma0.true)
# tau12.est.fun(par1 = gamma0.true, par2 = gamma0.true, 
#               par12 = gamma12.true, index12 = 3, index1 = 3, 
#               index2 = 3)


# Model specification
times = c("U1", "U2", "C")
deltas = c("delta_1", "delta_2", "delta_D")
fmlas = list("T1" = ~ 1, "T2" = ~ 1, "T3" = ~ 1)
Gfuns = c("T1" = "PH", "T2" = "PH", "T3" = "PH")
families = matrix(NA, 3, 3); 
families[1, 2] = families[1, 3] = families[2, 3] = "Clayton"
links = matrix(NA, 3, 3); 
links[1, 2] = links[1, 3] = links[2, 3] = "identity"
copula.fmlas = matrix(NA, 3, 3)
copula.fmlas[1, 2] = copula.fmlas[1, 3] = copula.fmlas[2, 3] = "~ 1"

#n = 500
#source("executeSimu2.R")

n = 500
source("executeSimu2.R")




