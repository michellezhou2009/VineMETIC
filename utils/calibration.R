## utils/calibration.R
##
## Censoring-rate calibration helpers used by settings/sim1_setting.R
## to build Simulation Study I's data-generating configuration.
##
## Ported from the original simulation harness (simulation_new/DataGenFuns.R,
## via LiDA/simulation/utils/calibration_orig.R). Trimmed to keep only the
## functions actually called by the current simulation scripts:
## `find.gammaD()`, `find.gammaT()`, `quantileT.fun()`. (The original file
## also defined `data.gen.fun()`, `calc.tau()`, and a standalone `survT.fun()`,
## none of which are invoked anywhere in the current Sim I/II/III code --
## `find.gammaD()`/`find.gammaT()` each carry their own internal copy of
## the `survT.fun` logic below, so nothing is lost by dropping the rest.)

find.gammaD = function(cen.rate,
                       betas = list(tD = c(2, 2), t1 = c(2, 2), t2 = c(2, 2)),
                       c.lwr = 1, c.upr = 6){

  survT.fun <- function(tt, betaT, gammaT, z1, z2){
    exp(- exp(gammaT) * tt * exp(z1 * betaT[1] + z2 * betaT[2]))
  }
  beta.tD = betas$tD; beta.t1 = betas$t1; beta.t2 = betas$t2
  obj.fun = function(bb, cen.rate){
    sapply(bb, function(b){
      rate.fun = function(cc, z1, z2){
        survT.fun(cc, betaT = beta.tD, gammaT = b, z1 = z1, z2 = z2) /
          (c.upr - c.lwr)
      }
      pracma::integral2(fun = rate.fun, xmin = c.lwr, ymin = -0.5,
                        xmax = c.upr, ymax = 0.5, z2 = 1)$Q * 0.3 +
        pracma::integral2(fun = rate.fun, xmin = c.lwr, ymin = -0.5,
                          xmax = c.upr, ymax = 0.5, z2 = 0)$Q * 0.7 - cen.rate
    })
  }
  uniroot(obj.fun, lower = - 20, upper = 20, cen.rate = cen.rate)$root
}

find.gammaT = function(
    cen.rate, gamma.tD, type = "t1", copula.fam, copula.par, copula.link,
    betas = list(tD = c(2, 2), t1 = c(2, 2), t2 = c(2, 2)),
    c.lwr = 1, c.upr = 6){
  survT.fun <- function(tt, betaT, gammaT, z1, z2){
    exp(- exp(gammaT) * tt * exp(z1 * betaT[1] + z2 * betaT[2]))
  }
  copula.index = MyCopIndex(copula.fam)
  beta.tD = betas$tD; beta.t1 = betas$t1; beta.t2 = betas$t2
  h.fun = MylinkFun(copula.link)$h.fun
  if (type == "t1") betaT = beta.t1 else betaT = beta.t2
  obj.fun = function(bb, cen.rate){
    sapply(bb, function(b){
      rate.fun = function(cc, z1, z2){
        survT.fun(cc, betaT = betaT, gammaT = b, z1 = z1, z2 = z2) /
          (c.upr - c.lwr)
      }
      prob.A = 1 -(pracma::integral2(fun = rate.fun, xmin = c.lwr, ymin = -0.5,
                                     xmax = c.upr, ymax = 0.5, z2 = 1)$Q * 0.3 +
                     pracma::integral2(fun = rate.fun, xmin = c.lwr, ymin = -0.5,
                                       xmax = c.upr, ymax = 0.5, z2 = 0)$Q * 0.7)
      rate.fun = function(uu, z1, z2, copula.index, copula.par, h.fun){
        tt = - log(uu) * exp(- gamma.tD - beta.tD[1] * z1 - beta.tD[2] * z2)
        ST.tt = survT.fun(tt, betaT = betaT, gammaT = b, z1 = z1, z2 = z2)
        alpha = h.fun(copula.par[1] + copula.par[2] * z1 + copula.par[3] * z2)
        1 - BiCopHfunc2(u1 = ST.tt, u2 = uu, family = copula.index, par = alpha)
      }
      prob.A * (pracma::integral2(rate.fun, xmin = 0, ymin = -0.5,
                                  xmax = 1, ymax = 0.5, z2 = 1,
                                  copula.index = copula.index,
                                  copula.par = copula.par, h.fun =h.fun)$Q * 0.3 +
                  pracma::integral2(rate.fun, xmin = 0, ymin = -0.5,
                                    xmax = 1, ymax = 0.5, z2 = 0,
                                    copula.index = copula.index,
                                    copula.par = copula.par, h.fun =h.fun)$Q * 0.7) -
        (1 - cen.rate)
    })
  }
  uniroot(obj.fun, lower = - 20, upper = 20, cen.rate = cen.rate)$root
}

quantileT.fun <- function(p, gammaT, upr = 100){

  objfun <- function(tt, p, gammaT){
    exp(- exp(gammaT) * tt)  - p
  }
  uniroot(objfun, lower = 0, upper = upr, p = p, gammaT = gammaT)$root
}
