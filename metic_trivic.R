metic.trivic = function(data, times, deltas, fmlas, Gfuns, copula.fams, 
                            copula.links, copula.fmlas, ncores = 1,
                        pc.method = "foreach"){
  
  N = nrow(data); J = length(times)
  
  ## estimation for TD ----
  timeD = times[J]; statusD = deltas[J]; fmlaD = fmlas[[J]]; GfunD = Gfuns[[J]]
  deltaD = data[, statusD]
  # time = timeD; status = statusD; formula = fmlaD; Gfun = GfunD
  myfit = fitSPT(data, time = timeD, status = statusD, formula = fmlaD, 
                 Gfun = GfunD)
  
  summary.Lambda = myfit$Lambda[, c("time", "est", "se")]
  if (is.null(myfit$beta)) 
    summary.beta = NULL else
      summary.beta = myfit$beta[, c("est", "se")]
  ui = predict.fitSPT(myfit, data)
  uD = ui$surv$fit
  Psi.uD = ui$Psi.surv
  tmp = myfit$Psi.theta
  Psi.theta = cbind(tmp$beta, tmp$dLambda)
  
  resD = list(
    summary.beta = summary.beta, summary.Lambda = summary.Lambda,
    ui = uD, Psi.ui = Psi.uD,
    Psi.theta = Psi.theta
  )

  # Estimate Tj and Copula (j,D) ----
  res = lapply(seq(J-1), function(j){
    timej = times[j]; statusj = deltas[j]; fmlaj = fmlas[[j]]; Gfunj = Gfuns[j]
    Xi = data[, timej]; deltaT = data[, statusj]
    copula.link = MylinkFun(links[j, J])
    Wmat = model.matrix(as.formula(copula.fmlas[j, J]), data); n.gamma = ncol(Wmat)
    copula.index = MyCopIndex(families[j, J])  
    control = MyCop(copula.index)
    
    ZmatT = model.matrix(fmlaj, data)[ , - 1, drop = F]; n.bT = ncol(ZmatT)
    if (n.bT == 0) b.ini = NULL else b.ini = rep(0, n.bT)
    tk = sort(unique(Xi[deltaT == 1])); n.tk = length(tk)
    dLambda.ini = rep(0.01, n.tk)
    gamma.ini = copula.link$hinv.fun(BiCopTau2Par(copula.index, tau = 0.1))
    if (n.gamma > 1) gamma.ini = c(gamma.ini, rep(0, n.gamma - 1))
    
    theta.ini = c(gamma.ini, b.ini, dLambda.ini)
    
    objfun <- function(x){
      dat_bvic = prepare_bvic(
        theta = x, thetaT = FALSE, Xi = Xi, deltaT = deltaT, ZmatT = ZmatT, 
        copula.link = copula.link, Wmat = Wmat, GfunT = Gfunj, control = control)
      out = llfuns.bvic(u1 = dat_bvic$ui, u2 = uD, d1 = deltaT, d2 = deltaD, 
                        copula.index = copula.index, alphai = dat_bvic$alphai, 
                        yes.constraint = dat_bvic$yes.constraint,
                        only.ll = FALSE, theta1 = FALSE, dat_bvic = dat_bvic,
                        copula.link = copula.link, Wmat = Wmat) 
      f = out$lln
      if (any(is.na(out$dll))) g = NA else g = apply(out$dll, 2, mean, na.rm = T)
      B = out$ddlln
      list(value = f, gradient = g, hessian = B)
    }
    
    est.res <- trust(objfun, theta.ini, 5, 100, iterlim = 300,
                     minimize= FALSE, blather = T)
    
    theta.est = est.res$argument
    dat_bvic = prepare_bvic(theta = theta.est, thetaT = FALSE, Xi = Xi, 
                            deltaT = deltaT, ZmatT = ZmatT, copula.link = copula.link, 
                            Wmat = Wmat, GfunT = Gfunj, control = control)
    lls = llfuns.bvic(u1 = dat_bvic$ui, u2 = uD, d1 = deltaT, d2 = deltaD, 
                      copula.index = copula.index, alphai = dat_bvic$alphai, 
                      yes.constraint = dat_bvic$yes.constraint,
                      only.ll = FALSE, theta1 = FALSE, dat_bvic = dat_bvic,
                      copula.link = copula.link, Wmat = Wmat, yes.dllDu2 = TRUE)
    Imat = - lls$ddlln
    dll.k = lls$dll.u2; dll.k[is.na(dll.k)] = 0
    Psi =  lls$dll; Psi[is.na(Psi)] = 0
    Psi.theta = Psi + Psi.uD %*% dll.k / N
    Vmat = crossprod(Psi.theta, Psi.theta) / N
    theta.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
    theta.se = sqrt(diag(theta.cov))
    Psi.theta = Psi.theta %*% solve(Imat) 
    
    gamma.est = theta.est[1 : n.gamma]
    gamma.se = theta.se[1 : n.gamma]
    Psi.gamma = Psi.theta[, c (1: n.gamma), drop = F]
    etai.gamma = as.vector(Wmat %*% matrix(gamma.est, byrow = F, ncol = 1))
    alphai = copula.link$h.fun(etai.gamma)
    summary.gamma = data.frame(est = gamma.est, se = gamma.se)
    rownames(summary.gamma) = colnames(Wmat)
    
    thetaT.est = theta.est[- c(1 : n.gamma)]
    thetaT.se = theta.se[- c(1 : n.gamma)]
    thetaT.cov = theta.cov[- c(1 : n.gamma), - c(1 : n.gamma)]
    Psi.thetaT = Psi.theta[, - c (1: n.gamma), drop = F]
    
    if (n.bT == 0){
      bT.est = bT.se = summary.beta = NULL
      dLambdaT.est = thetaT.est
      dLambdaT.cov = thetaT.cov
      bz = rep(0, N)
    } else{
      bT.est = thetaT.est[1 : n.bT]
      bT.se = thetaT.se[1 : n.bT]
      summary.beta = data.frame(est = bT.est, se = bT.se)
      rownames(summary.beta) = colnames(ZmatT)
      dLambdaT.est = thetaT.est[- c(1 : n.bT)]
      dLambdaT.cov = thetaT.cov[- c(1 : n.bT), - c(1 : n.bT)]
      bz = as.vector(ZmatT %*% bT.est); 
    }
    LambdaT.est = cumsum(dLambdaT.est)
    LambdaT.var = sapply(1 : n.tk, function(k){
      a = matrix(0, nrow = 1, ncol = n.tk) ; a[1, 1 : k] = 1
      as.numeric(
        a %*% dLambdaT.cov %*% t(a))
    })
    summary.LambdaT = data.frame(
      time = tk, est = LambdaT.est, se = sqrt(LambdaT.var)
    )
    
    ebz = exp(bz)
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT.est)
    dLambdaT.Xi = rep(0, N)
    dLambdaT.Xi[deltaT == 1] = dLambdaT.est[match(Xi[deltaT == 1], tk)]
    G.all =  G.funs(Gfunj)
    uiT = exp(- G.all$g.fun(LambdaT.Xi * ebz))
    Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T,nrow = N, ncol = n.tk))
    dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
    A = cbind(ZmatT * LambdaT.Xi * ebz, Xi.g.tk * ebz) *
      uiT * G.all$dg.fun(LambdaT.Xi * ebz)
    Psi.uiT = - Psi.theta[, - c(1 : n.gamma)] %*% t(A)
    
    out = list(
      summary.gamma = summary.gamma,
      summary.beta = summary.beta,
      summary.Lambda = summary.LambdaT,
      alphai = alphai, 
      gamma.est = gamma.est, Psi.gamma = Psi.gamma,
      thetaT.est = thetaT.est, Psi.thetaT = Psi.thetaT,
      ui = uiT, Psi.ui = Psi.uiT, 
      delta = deltaT, copula.index = copula.index, 
      Wmat = Wmat, copula.link = copula.link, 
      G.all = G.all, Lambda.Xi = LambdaT.Xi,
      dLambda.Xi = dLambdaT.Xi, dLambda.tk = dLambdaT.est,
      bz = bz, Xi.g.tk = Xi.g.tk, dN.tk = dN.tk, Zmat = ZmatT
    )
    out
  })
  names(res) = paste0("(", seq(J - 1), ",", J, ")")
  
  summary.all = lapply(res, function(out){ 
    list(beta = out$summary.beta, cumhaz = out$summary.Lambda)
  })
  names(summary.all) = paste0("T", seq(J-1))
  summary.all$TD = list(beta = resD$summary.beta, 
                        cumhaz = resD$summary.Lambda) 
  
  gamma.all = lapply(res, function(out){out$summary.gamma})
  names(gamma.all) = paste0("Cop", names(res))
  summary.all = c(summary.all, gamma.all)
  
  Psi.all = lapply(res, function(out){out$Psi.thetaT})
  names(Psi.all) = paste0("T", seq(J-1))
  Psi.all$TD = resD$Psi.theta
  
  Psi.gamma.all = lapply(res, function(out){out$Psi.gamma})
  names(Psi.gamma.all) = paste0("Cop", names(res))
  Psi.all = c(Psi.all, Psi.gamma.all)
  
  ## Estimate Cop(1,2) ----
  link12 = MylinkFun(links[1, 2])
  Wmat12 = model.matrix(as.formula(copula.fmlas[1, 2]), data); 
  n.gamma12 = ncol(Wmat12)
  index12 = MyCopIndex(families[1, 2]); control12 = MyCop(index12)
  gamma12.ini = link12$hinv.fun(BiCopTau2Par(index12, tau = 0.1))
  if (n.gamma12 > 1) gamma12.ini = c(gamma12.ini, rep(0, n.gamma12 - 1))
  
  if (n.gamma12 == 1){
    objfun = function(x){
      dat.trivic = prepare_trivic(theta = x, index12 = index12, 
                                  link12 = link12, Wmat12 = Wmat12, 
                                  control12 = control12, N = N, 
                                  thetaT = TRUE, datT = NULL)
      - lln.trivic.PMLE(
        uu1 = res[[1]]$ui, uu2 = res[[2]]$ui, uu3 = uD, 
        dd1 = res[[1]]$delta, dd2 = res[[2]]$delta, dd3 = deltaD, 
        par12 = dat.trivic$alphai, par1 = res[[1]]$alphai, par2 = res[[2]]$alphai,
        index12 = index12, index1 = res[[1]]$copula.index, 
        index2 = res[[2]]$copula.index, link12 = link12, 
        Wmat12 = Wmat12, yes.constraint = dat.trivic$yes.constraint,
        only.ll =  TRUE, pc.method = pc.method)$lln
    }
    gamma12.est = nlm(f = objfun, p = gamma12.ini)$estimate
  } else {
    objfun = function(x){
      dat.trivic = prepare_trivic(theta = x, index12 = index12, 
                                  link12 = link12, Wmat12 = Wmat12, 
                                  control12 = control12, N = N, 
                                  thetaT = TRUE, datT = NULL)
      out = lln.trivic.PMLE(
        uu1 = res[[1]]$ui, uu2 = res[[2]]$ui, uu3 = uD, 
        dd1 = res[[1]]$delta, dd2 = res[[2]]$delta, dd3 = deltaD, 
        par12 = dat.trivic$alphai, par1 = res[[1]]$alphai, par2 = res[[2]]$alphai,
        index12 = index12, index1 = res[[1]]$copula.index, 
        index2 = res[[2]]$copula.index, link12 = link12, 
        Wmat12 = Wmat12, yes.constraint = dat.trivic$yes.constraint,
        only.ll =  FALSE, pc.method = pc.method)
      f = out$lln
      if (any(is.na(out$dll))) g = NA else g = apply(out$dll, 2, mean, na.rm = T)
      B = out$ddlln
      list(value = f, gradient = g, hessian = B)
    }
    est.res <- trust(objfun, gamma12.ini, 5, 100, iterlim = 300,
                     minimize= FALSE, blather = T)
    gamma12.est = est.res$argument
  }
  
  dat.trivic = prepare_trivic(theta = gamma12.est, index12 = index12, 
                              link12 = link12, Wmat12 = Wmat12, 
                              control12 = control12, N = N, 
                              thetaT = TRUE, datT = NULL)
  lls = lln.trivic.PMLE(
    uu1 = res[[1]]$ui, uu2 = res[[2]]$ui, uu3 = uD, 
    dd1 = res[[1]]$delta, dd2 = res[[2]]$delta, dd3 = deltaD, 
    par12 = dat.trivic$alphai, par1 = res[[1]]$alphai, par2 = res[[2]]$alphai,
    index12 = index12, index1 = res[[1]]$copula.index, 
    index2 = res[[2]]$copula.index, link12 = link12, 
    Wmat12 = Wmat12, yes.constraint = dat.trivic$yes.constraint,
    only.ll =  FALSE, pc.method = pc.method)
  dlls = dll.trivic.PMLE(
    uu1 = res[[1]]$ui, uu2 = res[[2]]$ui, uu3 = uD, 
    dd1 = res[[1]]$delta, dd2 = res[[2]]$delta, dd3 = deltaD, 
    par12 = dat.trivic$alphai, par1 = res[[1]]$alphai, par2 = res[[2]]$alphai,
    index12 = index12, index1 = res[[1]]$copula.index, 
    index2 = res[[2]]$copula.index, yes.constraint = dat.trivic$yes.constraint,
    pc.method = pc.method)
  alpha12 = dat.trivic$alphai; alpha12.lp = link12$hinv.fun(alpha12)
  
  Imat = - lls$ddlln
  
  Psi.u1 = res[[1]]$Psi.ui
  dll = dlls$ll.Dpar12.Du1
  dll[is.na(dll)] = 0
  dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
  Psi.u1 = Psi.u1 %*% dll
  
  
  Psi.u2 = res[[2]]$Psi.ui
  dll = dlls$ll.Dpar12.Du2
  dll[is.na(dll)] = 0
  dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
  Psi.u2 = Psi.u2 %*% dll
  
  Psi.u3 = Psi.uD
  dll = dlls$ll.Dpar12.Du3
  dll[is.na(dll)] = 0
  dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
  Psi.u3 = Psi.u3 %*% dll
  
  alphai = res[[1]]$alphai; copula.link = res[[1]]$copula.link
  Psi.gamma = res[[1]]$Psi.gamma; Wmat = res[[1]]$Wmat
  alphai.lp = copula.link$hinv.fun(alphai) 
  dll = dlls$ll.Dpar12.Dpar1
  dll[is.na(dll)] = 0
  dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
  Psi.par1 = apply(dll, 2, function(x){
    apply(Psi.gamma * copula.link$dot.h.fun(alphai.lp) * Wmat * x, 
          1, sum)
  })
  
  alphai = res[[2]]$alphai; copula.link = res[[2]]$copula.link
  Psi.gamma = res[[2]]$Psi.gamma; Wmat = res[[2]]$Wmat
  alphai.lp = copula.link$hinv.fun(alphai) 
  dll = dlls$ll.Dpar12.Dpar2
  dll[is.na(dll)] = 0
  dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
  Psi.par2 = apply(dll, 2, function(x){
    apply(Psi.gamma * copula.link$dot.h.fun(alphai.lp) * Wmat * x, 
          1, sum)
  })
  
  Psi =  lls$dll
  Psi[is.na(Psi)] = 0
  Psi.gamma = Psi + 
    (Psi.u1 + Psi.u2 + Psi.u3 + Psi.par1 + Psi.par2) / N
  Vmat = crossprod(Psi.gamma, Psi.gamma) / N
  gamma.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  gamma.se = sqrt(diag(gamma.cov))
  Psi.gamma = Psi.gamma %*% solve(Imat)
  summary.gamma = data.frame(est = gamma12.est, se = gamma.se)
  rownames(summary.gamma) = colnames(Wmat12)
  
  summary.all$`Cop(1,2)` = summary.gamma
  Psi.all$`Cop(1,2)` = Psi.gamma
  summary.all$Psi = Psi.all
  return(summary.all) 
}

Surv.baseline = function(Lambda.est, tt){
  tk = Lambda.est$time; Lambda = Lambda.est$est; 
  Lambda.tt = c(0, Lambda)[sum.I(tt, ">=", tk) + 1]
  Lambda.se = Lambda.est$se
  Lambda.tt.se = c(0, Lambda.se)[sum.I(tt, ">=", tk) + 1]
  Shat.tt = exp(- Lambda.tt)
  Shat.tt.se = abs(- Shat.tt) * Lambda.tt.se
  # DLambda = - Shat.tt * (tt >= matrix(tk, byrow = T, nrow = length(tt), ncol = length(tk)))
  # Psi.Shat.tt = Psi.dLambda %*% t(DLambda)
  # Shat.tt.cov = crossprod(Psi.Shat.tt, Psi.Shat.tt) 
  # tmp = sqrt(diag(Shat.tt.cov)) / nrow(Psi.Lambda)
  # summary(Shat.tt.se - tmp)
  data.frame(tt = tt, surv = Shat.tt, se = Shat.tt.se)
}



SurvEst = function(Lambda.est, beta.est, Psi.theta, newdat, conf.int = FALSE,
                   level = 0.95){
  tk = Lambda.est$time; Lambda = Lambda.est$est
  Lambda.tt = c(0, Lambda)[sum.I(newdat$tt, ">=", tk) + 1]
  bz = as.vector(as.matrix(newdat[,c("z1", "z2")]) %*% 
                   matrix(beta.est, byrow = F, ncol = 1))
  ebz = exp(bz)
  Shat.tt = exp(- Lambda.tt * ebz)
  Dbeta = - Shat.tt * Lambda.tt * ebz * as.matrix(newdat[,c("z1", "z2")]) 
  DLambda = - Shat.tt * ebz * 
    (newdat$tt >= matrix(tk, byrow = T, nrow = nrow(newdat), ncol = length(tk)))
  Psi.Shat.tt = Psi.theta %*% t(cbind(Dbeta, DLambda))
  Shat.tt.cov = crossprod(Psi.Shat.tt, Psi.Shat.tt) 
  Shat.tt.se = sqrt(diag(Shat.tt.cov)) / nrow(Psi.theta)
  out = newdat %>% 
    mutate(surv = Shat.tt, se = Shat.tt.se)
  if (conf.int == TRUE) {
    Dbeta = as.matrix(newdat[,c("z1", "z2")])
    DLambda = (1 / Lambda.tt) * 
      (newdat$tt >= matrix(tk, byrow = T, nrow = nrow(newdat), ncol = length(tk)))
    Psi.Shat.cloglog = Psi.theta %*% t(cbind(Dbeta, DLambda))
    Shat.cloglog = log(Lambda.tt) + bz
    Shat.cloglog.se = sqrt(diag(crossprod(Psi.Shat.cloglog, Psi.Shat.cloglog))) /
      nrow(Psi.theta)
    critval = qnorm((1 - level) / 2, lower = F)
    lower = exp(-exp(Shat.cloglog - critval * Shat.cloglog.se))
    upper = exp(-exp(Shat.cloglog + critval * Shat.cloglog.se))
    out = out %>% mutate(lower = lower, upper = upper)
  } 
  return(out)
}