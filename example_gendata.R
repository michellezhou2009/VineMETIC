rm(list=ls())
library(tidyverse)
library(VineCopula)
library(parallel)
library(foreach)
library(doSNOW)

ncores = 2
cl = makeCluster(ncores)
registerDoSNOW(cl)

source("helpers.R")
source("DataGenFuns.R")

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
gamma2D.true = c(
  MylinkFun(copula.links[2, 3])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[2, 3]), tau = 0.4) 
  ), 0.1, 1)
gamma12.true = c(
  MylinkFun(copula.links[1, 2])$hinv.fun(
    BiCopTau2Par(family = MyCopIndex(copula.fams[1, 2]), tau = 0.2) 
  ), 1, 1)
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

nsamp = 500
data = data.gen.fun(nsamp = nsamp, copula.fams = copula.fams, 
                    copula.links = copula.links, 
                    copula.pars = copula.pars, 
                    gammas = gammas, seed0 = 20250214)
saveRDS(data, file = "expample.rds")
stopCluster(cl)
