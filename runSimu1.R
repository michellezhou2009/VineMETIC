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
source("DataGenFuns.R")
source("metic_trivic.R")

nrep = 500
seed1 = 20240122
set.seed(seed1); seeds = round(runif(nrep, 10, 10000000))

## Data generation setting ----
betas = list(tD = c(2, 2), t1 = c(2, 2), t2 = c(2, 2))
copula.fams = matrix(NA, 3, 3); 
copula.fams[1, 2] = "Frank"; 
copula.fams[1, 3] = "Gumbel"; copula.fams[2, 3] = "Clayton"
copula.links = matrix(NA, 3, 3)
copula.links[1, 2] = "identity"; 
copula.links[1, 3] = "log-1"; copula.links[2, 3] = "log"
gamma1D.true = c(
  MylinkFun(copula.links[1, 3])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[1, 3]), tau = 0.7)
  ), 1, 0.1)
copula.fam = copula.fams[1, 3]; copula.link = copula.links[1, 3]
gamma.true = gamma1D.true
lwr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] - 0.5 * gamma.true[2]))
upr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] + 0.5 * gamma.true[2] + gamma.true[3]))
tau.all = list(cop1D = c(lwr, upr))

gamma2D.true = c(
  MylinkFun(copula.links[2, 3])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[2, 3]), tau = 0.4) 
  ), 0.1, 1)
copula.fam = copula.fams[2, 3]; copula.link = copula.links[2, 3]
gamma.true = gamma2D.true
lwr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] - 0.5 * gamma.true[2]))
upr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] + 0.5 * gamma.true[2] + gamma.true[3]))
tau.all$cop2D = c(lwr, upr)

gamma12.true = c(
  MylinkFun(copula.links[1, 2])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[1, 2]), tau = 0.2) 
  ), 1, 1)
copula.fam = copula.fams[1, 2]; copula.link = copula.links[1, 2]
gamma.true = gamma12.true
lwr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] - 0.5 * gamma.true[2]))
upr = BiCopPar2Tau(family = MyCopIndex(copula.fam),
                   par = MylinkFun(copula.link)$h.fun(
                     gamma.true[1] + 0.5 * gamma.true[2] + gamma.true[3]))
tau.all$cop12 = c(lwr, upr)
gamma.true.all = list(gamma1D = gamma1D.true, gamma2D = gamma2D.true,
                      gamma12 = gamma12.true)

copula.pars = list(
  `(1,3)` = gamma1D.true, `(2,3)` = gamma2D.true, `(1,2)` = gamma12.true
)
c.lwr = 1; c.upr = 6
gamma.tD = find.gammaD(cen.rate = 0.1)
gamma.t1 = find.gammaT(cen.rate = 0.3, gamma.tD = gamma.tD, type = "t1",
                       copula.fam = copula.fams[1, 3], 
                       copula.link = copula.links[1, 3],
                       copula.par = copula.pars$`(1,3)`)
gamma.t2 = find.gammaT(cen.rate = 0.25, gamma.tD = gamma.tD, type = "t2",
                       copula.fam = copula.fams[2, 3], 
                       copula.link = copula.links[2, 3],
                       copula.par = copula.pars$`(2,3)`)
gammas = list(t1 = gamma.t1, t2 = gamma.t2, tD = gamma.tD)

pp = seq(0.05, 0.95, by = 0.01)
tt1 = sapply(pp, function(p){
  quantileT.fun(p = p, gammaT = gammas$t1)
})
tt2 = sapply(pp, function(p){
  quantileT.fun(p = p, gammaT = gammas$t2)
})
ttD = sapply(pp, function(p){
  quantileT.fun(p = p, gammaT = gammas$tD, upr = 200)
})
tt.all = list(pp = pp, t1 = tt1, t2 = tt2, tD = ttD)
save(tt.all, tau.all, gamma.true.all, gammas, file = "setting_(1,1).RData")

# Model specification -----
times = c("time1", "time2", "timeD")
deltas = c("delta1", "delta2", "deltaD")
fmlas = list("T1" = ~ z1 + z2, "T2" = ~ z1 + z2, "D" = ~ z1 + z2)
Gfuns = c("T1" = "PH", "T2" = "PH", "D" = "PH")
families = copula.fams
links = copula.links
copula.fmlas = matrix(NA, 3, 3)
copula.fmlas[1, 3] = "~ z1 + z2"; copula.fmlas[2, 3] = "~ z1 + z2"
copula.fmlas[1, 2] = "~ z1 + z2"

# Run -----
nsamp = 500
dir.nm = "simu_raw/simu1/"
bad.ls = NULL
source("executeSimu1.R")
