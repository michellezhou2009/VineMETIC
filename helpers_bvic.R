density_bvic <- function(copula.index){
  out = c();
  out[[1]] = function(u1, u2, para) {
    BiCopCDF(u1, u2, family = copula.index, par = para) # (0, 0)
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopHfunc1(u1, u2, family = copula.index, par = para) # (1, 0)
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopHfunc2(u1, u2, family = copula.index, par = para) # (0, 1)
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para)# (1, 1)
  } # (1, 1)
  return(out)
}

density.Dpar_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.para(u1, u2, para) 
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "par", family = copula.index, par = para) 
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "par", family = copula.index, par = para)
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "par", family = copula.index, par = para) 
  } # (1, 1)
  return(out)
}

density.Du1_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfunc1(u1, u2, family = copula.index, par = para) 
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "u2", family = copula.index, par = para) 
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para)
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "u1", family = copula.index, par = para)
  } # (1, 1)
  return(out)
}

density.Du2_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfunc2(u1, u2, family = copula.index, par = para)
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para) 
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "u2", family = copula.index, par = para) 
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "u2", family = copula.index, par = para) 
  } # (1, 1)
  return(out)
}

density.Dpar.Dpar_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.para2(u1, u2, para) 
  }
  out[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "par", family = copula.index, par = para) 
  }
  out[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "par", family = copula.index, par = para)
  }
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "par", family = copula.index, par = para) 
  }
  return(out)
}

density.Du1.Du1_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "u2", 
                    family = copula.index, par = para) 
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "u2", family = copula.index, par = para) 
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para) 
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "u1", family = copula.index, par = para) 
  } # (1, 1)
  return(out)
}

density.Du2.Du2_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "u2", 
                    family = copula.index, par = para) 
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para) 
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "u2", family = copula.index, par = para)
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "u2", family = copula.index, par = para)
  } # (1, 1)
  return(out)
}

density.Du1.Du2_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para) 
  } # (0, 0)
  out[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para)
  } # (1, 0)
  out[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para)
  } # (0, 1)
  out[[4]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.u1u1u2u2(u1, u2, para)
  } # (1, 1)
  return(out)
}

density.Dpar.Du1_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "par", family = copula.index, par = para) 
  }
  out[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "par1u2",
                     family = copula.index, par = para) 
  }
  out[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "par", family = copula.index, par = para) 
  }
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2, deriv = "par1u1", family = copula.index, par = para)
  }
  return(out)
}

density.Dpar.Du2_bvic <- function(copula.index){
  out = c()
  out[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "par", family = copula.index, par = para) 
  }
  out[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "par", family = copula.index, par = para) 
  }
  out[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "par1u2",
                     family = copula.index, par = para)
  }
  out[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2, deriv = "par1u2", family = copula.index, par = para)
  }
  return(out)
}

funsData_bvic <- function(Funs, d1, d2, u1, u2, copula.index, para){
  funs = Funs(copula.index)
  grp = d1 + 2 * d2 + 1
  out = lapply(1 : 4, function(k){
    index = which(grp == k) 
    if (sum(index) == 0) return(NULL) else {
      ll = funs[[k]](u1 = u1[index], u2 = u2[index], para = para[index])
      data.frame(index = index, ll = ll)
    }
  }) %>% do.call(rbind, .)
  out = out %>% arrange(index, )
  out$ll
}

prepare_bvic = function(theta, thetaT = FALSE, Xi = NULL, deltaT = NULL, 
                        ZmatT = NULL,  copula.link, Wmat, GfunT = "PH", 
                        control){
  
  theta0 = theta
  h.fun = copula.link$h.fun
  n.para = ncol(Wmat); copula.para = theta0[1 : n.para]
  copula.lp = as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
  alphai = h.fun(copula.lp)
  theta0 = theta0[- c(1 : n.para)]
  yes.constraint = min(alphai) < control$lwr | 
    max(alphai) > control$upr
  
  if (!thetaT){
    if (is.null(Xi) | is.null(deltaT) | is.null(ZmatT)) 
      stop("Xi, ZmatT should be provide.")
    N = length(Xi)
    tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
    Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk)) 
    dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
    n.bT = ncol(ZmatT);
    if (n.bT == 0){
      dLambdaT = theta0; bZ.T = rep(0, N) 
    } else{
      betaT = theta0[c(1 : n.bT)]; dLambdaT = theta0[n.bT + c(1 : n.tk)] 
      bZ.T = as.vector(ZmatT %*% betaT); 
    }
    dLambdaT.Xi = rep(0, N)
    dLambdaT.Xi[deltaT == 1] = dLambdaT[match(Xi[deltaT == 1], tk)]
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
    G.all = G.funs(GfunT)
    uiT = exp(- G.all$g.fun(LambdaT.Xi * exp(bZ.T)))
    yes.constraint = yes.constraint | any(dLambdaT < 0) 
    out = list(
      alphai = alphai, yes.constraint = yes.constraint,
      ui = uiT, G.all = G.all, Lambda.Xi = LambdaT.Xi, 
      dLambda.Xi = dLambdaT.Xi, dLambda.tk = dLambdaT,
      bz = bZ.T, Xi.g.tk = Xi.g.tk, dN.tk = dN.tk, Zmat = ZmatT)
  } else {
    out = list(alphai = alphai, yes.constraint = yes.constraint)
  } 
  return(out)
}




llfuns.bvic <- function(u1, u2, d1, d2, copula.index, alphai, yes.constraint, 
                        only.ll = TRUE, theta1 = FALSE, dat_bvic = NULL, 
                        copula.link, Wmat, yes.dllDu2 = FALSE){
  if (yes.constraint) 
    return(list(lln = -Inf, dll = NA, ddlln =NA)) else {
      N = length(u1)
      if (copula.index == 4) { # gumbel
        u1[u1 == 1] = exp(- 1 / N)
        u2[u2 == 1] = exp(- 1 / N)
      }
      dd = funsData_bvic(density_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
                    copula.index = copula.index, para = alphai)
      ll = log(dd)
      if (!theta1){
        dLambda.X1 = dat_bvic$dLambda.Xi; G.all = dat_bvic$G.all
        Lambda.X1 = dat_bvic$Lambda.Xi; bz1 = dat_bvic$bz
        X1.g.tk = dat_bvic$Xi.g.tk; dN1.tk = dat_bvic$dN.tk
        dLambda1 = dat_bvic$dLambda.tk; n.tk1 = length(dLambda1)
        Zmat1 = dat_bvic$Zmat; n.b1 = ncol(Zmat1)
        
        log.dLambda.X1 = rep(0, N); 
        log.dLambda.X1[d1 == 1] = log(dLambda.X1[d1 == 1])
        ll = ll + 
          d1 * (- G.all$g.fun(Lambda.X1 * exp(bz1)) + 
                  log(G.all$dg.fun(Lambda.X1 * exp(bz1))) + bz1 + log.dLambda.X1)
        
      } 
      lln = mean(ll, na.rm = T)
      if (only.ll){
        return(list(lln = lln))
      } else {
        ## dll.Dpar & dll.Dpar.Dpar -----
        alphai.lp = copula.link$hinv.fun(alphai); n.gamma = ncol(Wmat)
        dd.Dpar = funsData_bvic(
          Funs = density.Dpar_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
          copula.index = copula.index, para = alphai)
        ll.Dpar = (dd.Dpar / dd) 
        ll.Dgamma = ll.Dpar *  copula.link$dot.h.fun(alphai.lp) * Wmat
        dll = ll.Dgamma
        
        dd.Dpar.Dpar = funsData_bvic(
          Funs = density.Dpar.Dpar_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
          copula.index = copula.index, para = alphai)
        ll.Dpar.Dpar = dd.Dpar.Dpar / dd - ll.Dpar ^ 2 
        ll.Dgamma.Dgamma = matrix(
          apply((ll.Dpar.Dpar * (copula.link$dot.h.fun(alphai.lp) ^ 2) + 
                   ll.Dpar * copula.link$ddot.h.fun(alphai.lp)) * 
                  Wmat[, rep(1 : n.gamma, n.gamma), drop = F] * 
                  Wmat[, rep(1 : n.gamma, each = n.gamma), drop = F], 
                2, mean, na.rm = T),
          byrow = T, nrow = n.gamma, ncol = n.gamma)
        ddlln = ll.Dgamma.Dgamma
        
        if (!theta1){
          ## ll.Dtheta1 ----
          ebz1 = exp(bz1); 
          dd.Du1 = funsData_bvic(
            Funs = density.Du1_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
            copula.index = copula.index, para = alphai)
          ll.Du1 = dd.Du1 / dd
          ll.DLambda1 =
            - ll.Du1 * u1  * ebz1 * X1.g.tk * G.all$dg.fun(Lambda.X1 * ebz1) -
            d1 * ebz1 * X1.g.tk * G.all$dg.fun(Lambda.X1 * ebz1) +
            d1 * ebz1 * X1.g.tk * G.all$ddg.fun(Lambda.X1 * ebz1) /
            G.all$dg.fun(Lambda.X1) +
            dN1.tk / matrix(dLambda1, byrow = T, nrow = N, ncol = n.tk1)
          if (n.b1 == 0) {
            ll.Db1 = NULL 
          } else {
            ll.Db1 = 
              (- ll.Du1 * u1 * Lambda.X1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) -
                 d1 * Lambda.X1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) +
                 d1 * Lambda.X1 * ebz1 * G.all$ddg.fun(Lambda.X1 * ebz1) /
                 G.all$dg.fun(Lambda.X1 * ebz1) + d1
              ) * Zmat1
          } # end for n.b1 == 0
          ll.Dtheta1 = cbind(ll.Db1, ll.DLambda1)
          dll = cbind(dll, ll.Dtheta1)
          
          ## ll.Dtheta1.Dtheta1 ----
          dd.Du1.Du1 = funsData_bvic(
            Funs = density.Du1.Du1_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
            copula.index = copula.index, para = alphai) 
          ll.Du1.Du1 = dd.Du1.Du1 / dd - ll.Du1 ^ 2
          ll.DLambda1.DLambda1 = crossprod(
            (ll.Du1.Du1 * u1 ^ 2 * ebz1 ^ 2 * G.all$dg.fun(Lambda.X1 * ebz1) ^ 2 + 
               ll.Du1 * u1 * ebz1 ^ 2 * G.all$dg.fun(Lambda.X1 * ebz1) ^ 2 -
               ll.Du1 * u1 * ebz1 ^ 2 * G.all$ddg.fun(Lambda.X1 * ebz1) -
               d1 * ebz1 ^ 2 * G.all$ddg.fun(Lambda.X1 * ebz1) +
               d1 * ebz1 ^ 2 * 
               (G.all$dddg.fun(Lambda.X1 * ebz1) / G.all$dg.fun(Lambda.X1 * ebz1) - 
                  (G.all$ddg.fun(Lambda.X1 * ebz1) / 
                     G.all$dg.fun(Lambda.X1 * ebz1)) ^ 2)) * 
              X1.g.tk,  X1.g.tk) / N - diag(apply(dN1.tk, 2, mean) / dLambda1 ^ 2)
          
          dd.Dpar.Du1 = funsData_bvic(
            Funs = density.Dpar.Du1_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
            copula.index = copula.index, para = alphai) 
          ll.Dpar.Du1 = dd.Dpar.Du1 / dd - ll.Dpar * ll.Du1
          ll.Dpar.DLambda1 = - ll.Dpar.Du1 * u1 * ebz1 * X1.g.tk * 
            G.all$dg.fun(Lambda.X1 * ebz1)
          ll.Dgamma.DLambda1 = matrix(
            apply(ll.Dpar.DLambda1[, rep(1 : n.tk1, n.gamma)] *
                    copula.link$dot.h.fun(alphai.lp) *
                    Wmat[, rep(1 : n.gamma, each = n.tk1)],
                  2, mean, na.rm = T),
            byrow = T, ncol = n.tk1, nrow = n.gamma) 
          
          if (n.b1 == 0){
            ddlln = rbind(
              cbind(ddlln, ll.Dgamma.DLambda1),
              cbind(t(ll.Dgamma.DLambda1), ll.DLambda1.DLambda1)
            )
          } else {
            ll.Db1.Db1 = crossprod(
              (ll.Du1.Du1 * u1 ^ 2 * Lambda.X1 ^ 2 * ebz1 ^ 2 * 
                 G.all$dg.fun(Lambda.X1 * ebz1) ^ 2 +
                 ll.Du1 * u1 * Lambda.X1 ^ 2 * ebz1 ^ 2 * 
                 G.all$dg.fun(Lambda.X1 * ebz1) ^ 2  -
                 ll.Du1 * u1 * Lambda.X1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) - 
                 ll.Du1 * u1 * Lambda.X1 ^2  * ebz1 ^2 * 
                 G.all$ddg.fun(Lambda.X1 * ebz1) -
                 d1 * Lambda.X1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) - 
                 d1 * Lambda.X1 ^2  * ebz1 ^2  * 
                 G.all$ddg.fun(Lambda.X1 * ebz1) +
                 d1 * Lambda.X1 * ebz1 * G.all$ddg.fun(Lambda.X1 * ebz1) /
                 G.all$dg.fun(Lambda.X1 * ebz1) +
                 d1 * Lambda.X1 ^ 2 * ebz1 ^ 2 * 
                 (G.all$dddg.fun(Lambda.X1 * ebz1) / 
                    G.all$dg.fun(Lambda.X1 * ebz1) - 
                    (G.all$ddg.fun(Lambda.X1 * ebz1) / 
                       G.all$dg.fun(Lambda.X1 * ebz1)) ^ 2
                 )) * Zmat1, Zmat1
            ) / N
            ll.Db1.DLambda1 = crossprod(
              (ll.Du1.Du1 * u1 ^ 2 * ebz1 ^ 2 * Lambda.X1 * 
                 G.all$dg.fun(Lambda.X1 * ebz1) ^ 2 + 
                 ll.Du1 * u1 * ebz1 ^ 2 * Lambda.X1 * 
                 G.all$dg.fun(Lambda.X1 * ebz1) ^ 2 - 
                 ll.Du1 * u1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) -
                 ll.Du1 * u1 * ebz1 ^ 2 * Lambda.X1 * 
                 G.all$ddg.fun(Lambda.X1 * ebz1)  - 
                 d1 * ebz1 * G.all$dg.fun(Lambda.X1 * ebz1) - 
                 d1 * ebz1 ^ 2 * Lambda.X1 * G.all$ddg.fun(Lambda.X1 * ebz1) + 
                 d1 * ebz1 * G.all$ddg.fun(Lambda.X1 * ebz1) /
                 G.all$dg.fun(Lambda.X1 * ebz1) +
                 d1 * ebz1 ^ 2 *  Lambda.X1 * 
                 (G.all$dddg.fun(Lambda.X1 * ebz1) / 
                    G.all$dg.fun(Lambda.X1 * ebz1) - 
                    (G.all$ddg.fun(Lambda.X1 * ebz1) / 
                       G.all$dg.fun(Lambda.X1 * ebz1)) ^ 2
                 )) * X1.g.tk, Zmat1) / N
            ll.Dpar.Db1 = - ll.Dpar.Du1 * u1 * ebz1 * Lambda.X1 * Zmat1 * 
              G.all$dg.fun(Lambda.X1 * ebz1)
            ll.Dgamma.Db1 = matrix(
              apply(ll.Dpar.Db1[, rep(1 : n.b1, n.gamma)] * 
                      copula.link$dot.h.fun(alphai.lp) * 
                      Wmat[, rep(1 : n.gamma, each = n.b1)],
                    2, mean, na.rm = T),
              byrow = T, ncol = n.b1, nrow = n.gamma)
            ddlln = rbind(
              cbind(ddlln, ll.Dgamma.Db1, ll.Dgamma.DLambda1),
              cbind(t(ll.Dgamma.Db1), ll.Db1.Db1, t(ll.Db1.DLambda1)),
              cbind(t(ll.Dgamma.DLambda1), ll.Db1.DLambda1, 
                    ll.DLambda1.DLambda1)
            )
          } # end for n.b1 == 0
        } # end for !theta1
        ## DparDu2 -----
        if (yes.dllDu2){
          dd.Du2 = funsData_bvic(
            Funs = density.Du2_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
            copula.index = copula.index, para = alphai) 
          ll.Du2 = dd.Du2 / dd
          dd.Dpar.Du2 = funsData_bvic(
            Funs = density.Dpar.Du2_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
            copula.index = copula.index, para = alphai) 
          ll.Dpar.Du2 = dd.Dpar.Du2 / dd - ll.Dpar * ll.Du2
          ll.Dgamma.Du2 = ll.Dpar.Du2 *  copula.link$dot.h.fun(alphai.lp) * Wmat
          
          if (!theta1){
            dd.Du1.Du2 = funsData_bvic(
              Funs = density.Du1.Du2_bvic, d1 = d1, d2 = d2, u1 = u1, u2 = u2,
              copula.index = copula.index, para = alphai) 
            ll.Du1.Du2 = dd.Du1.Du2 / dd - ll.Du1 * ll.Du2
            ll.DLambda1.Du2 = - ll.Du1.Du2 * u1  * ebz1 * X1.g.tk * 
              G.all$dg.fun(Lambda.X1 * ebz1) 
            if (n.b1 > 0) {
              ll.Db1.Du2 = - ll.Du1.Du2 * u1  * Lambda.X1 * ebz1 * 
                G.all$dg.fun(Lambda.X1 * ebz1) * Zmat1
            } else ll.Db1.Du2 = NULL
            dll.u2 = cbind(ll.Dgamma.Du2, ll.Db1.Du2, ll.DLambda1.Du2)
          } else dll.u2 = ll.Dgamma.Du2
          return(list(lln = lln, dll = dll, ddlln = ddlln, dll.u2 = dll.u2))
        } else {
          return(list(lln = lln, dll = dll, ddlln = ddlln))
        }
      } # end for only.ll
    }
}
