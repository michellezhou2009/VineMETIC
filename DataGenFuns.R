data.gen.fun = 
  function(nsamp, copula.fams, copula.links, copula.pars, gammas, 
           betas = list(tD = c(2, 2), t1 = c(2, 2), t2 = c(2, 2)),
           c.lwr = 1, c.upr = 6, seed0 = 1234){
    
    beta.tD = betas$tD; beta.t1 = betas$t1; beta.t2 = betas$t2
    a1D = copula.pars$`(1,3)`; a2D = copula.pars$`(2,3)`; 
    a12 = copula.pars$`(1,2)`
    link1D = MylinkFun(copula.links[1, 3])$h.fun
    link2D = MylinkFun(copula.links[2, 3])$h.fun
    link12 = MylinkFun(copula.links[1, 2])$h.fun
    cop1D.index = MyCopIndex(copula.fams[1, 3])
    cop2D.index = MyCopIndex(copula.fams[2, 3])
    cop12.index = MyCopIndex(copula.fams[1, 2])
    gamma.tD = gammas$tD; gamma.t1 = gammas$t1; gamma.t2 = gammas$t2
    set.seed(seed0)
    
    survT.fun <- function(tt, betaT, gammaT, z1, z2){
      exp(- exp(gammaT) * tt * exp(z1 * betaT[1] + z2 * betaT[2])) 
    }
    
    aa = runif(nsamp, c.lwr, c.upr)
    z1 = runif(nsamp, -0.5, 0.5)
    z2 = rbinom(nsamp, size = 1, prob = 1 / 3)
    alpha1D = link1D(a1D[1] + a1D[2] * z1 + a1D[3] * z2)
    alpha2D = link2D(a2D[1] + a2D[2] * z1 + a2D[3] * z2)
    alpha12 = link12(a12[1] + a12[2] * z1 + a12[3] * z2)
    
    ttD = exp( - (gamma.tD + z1 * beta.tD[1] + z2 * beta.tD[2]) + 
                 log(- log(1 - runif(nsamp))))
    uuD = survT.fun(ttD, beta.tD, gamma.tD, z1, z2)
    uu12 = BiCopSim(N = nsamp, family = cop12.index, par = alpha12)
    
    uu = foreach(i = 1 : nsamp, .packages = c("VineCopula")) %dopar% {
      solve.u.fun <- function(u1, uD, alpha, copula.index){
        fun0 = function(uu, u1, uD, copula.index, alpha){
          BiCopHfunc2(u1 = uu, u2 = rep(uD, length(uu)), 
                      family = copula.index, par = alpha) - u1
        }
        uniroot(fun0, lower = 0, upper = 1, u1 = u1, uD = uD, 
                copula.index = copula.index, alpha = alpha)$root
      }
      u1 = solve.u.fun(u1 = uu12[i, 1], uD = uuD[i], alpha = alpha1D[i], 
                  copula.index = cop1D.index)
      u2 = solve.u.fun(u1 = uu12[i, 2], uD = uuD[i], alpha = alpha2D[i], 
                       copula.index = cop2D.index)
      c(u1, u2)
    }
    
    uu = foreach(i = 1 : nsamp, .packages = c("VineCopula"), .combine = rbind) %dopar% {
      solve.u.fun <- function(u1, uD, alpha, copula.index){
        fun0 = function(uu, u1, uD, copula.index, alpha){
          BiCopHfunc2(u1 = uu, u2 = rep(uD, length(uu)), 
                      family = copula.index, par = alpha) - u1
        }
        uniroot(fun0, lower = 0, upper = 1, u1 = u1, uD = uD, 
                copula.index = copula.index, alpha = alpha)$root
      }
      u1 = solve.u.fun(u1 = uu12[i, 1], uD = uuD[i], alpha = alpha1D[i], 
                  copula.index = cop1D.index)
      u2 = solve.u.fun(u1 = uu12[i, 2], uD = uuD[i], alpha = alpha2D[i], 
                       copula.index = cop2D.index)
      data.frame(u1 = u1, u2 = u2)
    }
    
    uu1 = uu[,1]; uu2 = uu[, 2]
    tt1 = exp(- (gamma.t1 + z1 * beta.t1[1] + z2 * beta.t1[2]) + log(- log(uu1)))
    tt2 = exp(- (gamma.t2 + z1 * beta.t2[1] + z2 * beta.t2[2]) + log(- log(uu2)))
    # ttD = round(ttD, 4); tt1 = round(tt1, 4); tt2 = round(tt2, 4)
    dat = data.frame(timeD = pmin(ttD, aa), deltaD = 1 * (ttD <= aa), 
               z1 = z1, z2 = z2) %>%
      mutate(time1 = pmin(tt1, timeD), delta1 = 1 * (tt1 <= timeD), 
             time2 = pmin(tt2, timeD), delta2 = 1 * (tt2 <= timeD)) %>%
      select(time1, delta1, time2, delta2, timeD, deltaD, z1, z2)
    dat = dat %>% mutate(time1 = round(time1, 7),
                         time2 = round(time2, 7),
                         timeD = round(timeD, 7))
    dat
  }


calc.tau = function(copula.par, copula.link, copula.index){
  h.fun = MylinkFun(copula.link)$h.fun
  obj.fun = function(z1, z2, copula.par, h.fun){
    alpha = h.fun(copula.par[1] + copula.par[2] * z1 + copula.par[2]) 
    BiCopPar2Tau(family = copula.index, par = alpha)
  }
  pracma::integral(fun = obj.fun, xmin = 1, xmax = 2, z2 = 1, 
                   copula.par =  copula.par, h.fun = h.fun) * 0.3 +
    pracma::integral(fun = obj.fun, xmin = 1, xmax = 2, z2 = 0, 
                     copula.par =  copula.par, h.fun = h.fun) * 0.7
}


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

survT.fun <- function(tt, betaT, gammaT, z1, z2){
  exp(- exp(gammaT) * tt * exp(z1 * betaT[1] + z2 * betaT[2])) 
}
