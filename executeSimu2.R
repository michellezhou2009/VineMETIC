ncores = 2
cl = makeCluster(ncores)
registerDoSNOW(cl)
for (k in 1 : 2){
  cat("Rep", k, "\n")
  file.nm = paste0(file.dir, "nsamp", n, "_tau", tau0.true, "_rep", k, ".RData")
  t1 = Sys.time()
  if (!file.exists(file.nm)){
    set.seed(seed.all[k])
    sam = rnacopula(n, onacopula("Clayton", C(theta0, 3, C(theta12, c(1, 2)))))
    T1 = qweibull(sam[, 1], shape = 2, scale = 70, lower.tail = F)
    T2 = qweibull(sam[, 2], shape = 2, scale = 60, lower.tail = F)
    D  = qweibull(sam[, 3], shape = 2, scale = 85, lower.tail = F)
    Ca = rexp(nrow(sam), rate = 1/350)
    C = pmin(D, Ca)
    U1 = pmin(T1, C)
    U2 = pmin(T2, C)
    delta_D = ifelse(D  <= Ca, 1, 0)
    delta_1 = ifelse(T1 <= C, 1, 0)
    delta_2 = ifelse(T2 <= C, 1, 0)
    cen.rate = data.frame("T1" = 1 - mean(delta_1), 
                          "T2" = 1 - mean(delta_2), 
                          "TD" = 1 - mean(delta_D), 
                          nsamp = n, rep = k)
    data = data.frame(U1, U2, C, delta_1, delta_2,delta_D)
    
    ## Li's method ----
    res.Li = est_cop_par(copula = "clayton", data = data)
    gamma0.est = as.numeric(res.Li[1])
    gamma12.est = as.numeric(res.Li[2])
    tau0.est = BiCopPar2Tau(family = 3, par = gamma0.est)
    tau12.est = BiCopPar2Tau(family = 3, par = gamma12.est)
    out.nested = data.frame(
      para = c("gamma0", "gamma12", "tau0", "tau12"), 
      est = c(gamma0.est, gamma12.est, tau0.est, tau12.est), 
      method = "nested",
      true = c(theta0, theta12, tau0.true, tau0.true))
    mydat1 = npest.star_12(data)$X
    S.X1 = data.frame(
      time = data$U1, 
      surv = g.fn(mydat1$S1.star.hat, mydat1$SD.U1.hat, gamma0.est)
    ) %>% arrange(time)
    surv.T1 = c(1, S.X1$surv)[sum.I(tt.all[[1]], ">=", S.X1$time) + 1]
    S.X2 = data.frame(
      time = data$U2, 
      surv = g.fn(mydat1$S2.star.hat, mydat1$SD.U2.hat, gamma0.est)
    ) %>% arrange(time)
    surv.T2 = c(1, S.X2$surv)[sum.I(tt.all[[2]], ">=", S.X2$time) + 1]
    out.nested = rbind(
      out.nested,
      data.frame(
        para = "surv.T1", est = surv.T1, method = "nested", 
        true = 1 - pp.true
      ),
      data.frame(
        para = "surv.T2", est = surv.T2, method = "nested", 
        true = 1 - pp.true
      ))
    
    
    ## Vine copula ----
    N = nrow(data); J = 3
    timeD = times[J]; statusD = deltas[J]; fmlaD = fmlas[[J]]; GfunD = Gfuns[[J]]
    deltaD = data[, statusD]
    myfit = fitSPT(data, time = timeD, status = statusD, formula = fmlaD, 
                   Gfun = GfunD)
    ui = predict.fitSPT(myfit, data)
    uD = ui$surv$fit
    
    # gamma1 not equal to gamma2 ----
    res = lapply(seq(J-1), function(j){
      timej = times[j]; statusj = deltas[j]; fmlaj = fmlas[[j]]; Gfunj = Gfuns[j]
      Xi = data[, timej]; deltaT = data[, statusj]
      copula.link = MylinkFun(links[j, J])
      Wmat = model.matrix(as.formula(copula.fmlas[j, J]), data); n.gamma = ncol(Wmat)
      copula.index = MyCopIndex(families[j, J])  
      control = MyCop(copula.index)
      
      ZmatT = model.matrix(fmlaj, data)[ , - 1, drop = F]; 
      b.ini = NULL
      tk = sort(unique(Xi[deltaT == 1])); n.tk = length(tk)
      dLambda.ini = rep(0.01, n.tk)
      gamma.ini = copula.link$hinv.fun(BiCopTau2Par(copula.index, tau = 0.1))
      theta.ini = c(gamma.ini, b.ini, dLambda.ini)
      
      objfun <- function(x){
        dat_bvic = prepare_bvic(
          theta = x, thetaT = FALSE, Xi = Xi, deltaT = deltaT, ZmatT = ZmatT, 
          copula.link = copula.link, Wmat = Wmat, GfunT = Gfunj, 
          control = control)
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
      
      est.res = trust(objfun, theta.ini, 5, 100, iterlim = 300,
                      minimize= FALSE, blather = T)
      
      theta.est = est.res$argument
      gamma.est = theta.est[1 : n.gamma]
      etai.gamma = as.vector(Wmat %*% matrix(gamma.est, byrow = F, ncol = 1))
      alphai = copula.link$h.fun(etai.gamma)
      dLambdaT.est = theta.est[- c(1 : n.gamma)]
      bz = rep(0, N); ebz = exp(bz)
      LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT.est)
      dLambdaT.Xi = rep(0, N)
      dLambdaT.Xi[deltaT == 1] = dLambdaT.est[match(Xi[deltaT == 1], tk)]
      G.all =  G.funs(Gfunj)
      uiT = exp(- G.all$g.fun(LambdaT.Xi * ebz))
      Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T,nrow = N, ncol = n.tk))
      dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
      
      tt = tt.all[[j]]
      LambdaT.tt = sum.I(tt, ">=", tk, dLambdaT.est)
      surv.tt = exp(- LambdaT.tt)
      
      out = list(
        gamma.est = gamma.est,
        alphai = alphai, ui = uiT, 
        delta = deltaT, copula.index = copula.index, 
        Wmat = Wmat, copula.link = copula.link, 
        G.all = G.all, Lambda.Xi = LambdaT.Xi,
        dLambda.Xi = dLambdaT.Xi, dLambda.tk = dLambdaT.est,
        bz = bz, Xi.g.tk = Xi.g.tk, dN.tk = dN.tk, Zmat = ZmatT,
        surv.tt = surv.tt
      )
      out
    })
    link12 = MylinkFun(links[1, 2])
    Wmat12 = model.matrix(as.formula(copula.fmlas[1, 2]), data); 
    n.gamma12 = ncol(Wmat12)
    index12 = MyCopIndex(families[1, 2]); control12 = MyCop(index12)
    gamma12.ini = link12$hinv.fun(BiCopTau2Par(index12, tau = 0.1))
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
        only.ll = TRUE)
      return(-out$lln)
    }
    gamma12.est = nlm(f = objfun, p = gamma12.ini)$estimate
    gamma1.est = res[[1]]$gamma.est
    gamma2.est = res[[2]]$gamma.est
    tau12.est = tau12.est.fun(par1 = gamma1.est, par2 = gamma2.est, 
                              par12 = gamma12.est, index12 = 3, index1 = 3, 
                              index2 = 3)
    tau1.est = BiCopPar2Tau(family = 3, par = gamma1.est)
    tau2.est = BiCopPar2Tau(family = 3, par = gamma2.est)
    out.vine.diff = data.frame(
      para = c("gamma0", "gamma0", "gamma12","tau1", "tau2", "tau12"), 
      est = c(gamma1.est, gamma2.est, gamma12.est, tau1.est, tau2.est, tau12.est), 
      method = c("vine_different_T1", "vine_different_T2", "vine_different"),
      true = c(gamma0.true, gamma0.true, gamma12.true, tau0.true, tau0.true, tau0.true))
    out.vine.diff = rbind(
      out.vine.diff,
      data.frame(
        para = "surv.T1", est = res[[1]]$surv.tt, method = "vine_diff", 
        true = 1 - pp.true
      ),
      data.frame(
        para = "surv.T2", est = res[[2]]$surv.tt, method = "vine_diff", 
        true = 1 - pp.true
      ))
    
    ## gamma1 == gamma2 ----
    X1 = data[, times[1]]; delta1 = data[, deltas[1]]
    X2 = data[, times[2]]; delta2 = data[, deltas[2]]
    tk1 = sort(unique(X1[delta1 == 1])); n.tk1 = length(tk1)
    tk2 = sort(unique(X2[delta2 == 1])); n.tk2 = length(tk2)
    
    gamma0.ini = BiCopTau2Par(family = copula.index, tau = 0.1)
    dLambda1.ini = rep(0.01, n.tk1)
    dLambda2.ini = rep(0.01, n.tk2)
    theta.ini = c(gamma0.ini, dLambda1.ini, dLambda2.ini)
    objfun.S1 = function(x){
      gamma0 = x[1]; dLambda1 = x[1 + c(1 : n.tk1)]; 
      dLambda2 = x[-c(1 : (1 + n.tk1))]
      dat1 = prepare_bvic(
        theta = c(gamma0, dLambda1), thetaT = FALSE, Xi = X1, 
        deltaT = delta1, ZmatT = res[[1]]$Zmat, 
        copula.link = res[[1]]$copula.link, 
        Wmat = res[[1]]$Wmat, GfunT = Gfuns[1], 
        control = MyCop(res[[1]]$copula.index))
      dat2 = prepare_bvic(
        theta = c(gamma0, dLambda2), thetaT = FALSE, Xi = X2, 
        deltaT = delta2, ZmatT = res[[2]]$Zmat, 
        copula.link = res[[2]]$copula.link, 
        Wmat = res[[2]]$Wmat, GfunT = Gfuns[2], 
        control = MyCop(res[[2]]$copula.index))
      
      if(dat1$yes.constraint | dat2$yes.constraint){
        f = -Inf; g = NA; B = NA
      } else {
        out1 = llfuns.bvic(
          u1 = dat1$ui, u2 = uD, d1 = delta1, d2 = deltaD, 
          copula.index = copula.index, alphai = dat1$alphai, 
          yes.constraint = dat1$yes.constraint,
          only.ll = FALSE, theta1 = FALSE, dat_bvic = dat1,
          copula.link = res[[1]]$copula.link, Wmat = res[[1]]$Wmat)
        
        out2 = llfuns.bvic(
          u1 = dat2$ui, u2 = uD, d1 = delta2, d2 = deltaD, 
          copula.index = copula.index, alphai = dat2$alphai, 
          yes.constraint = dat2$yes.constraint,
          only.ll = FALSE, theta1 = FALSE, dat_bvic = dat2,
          copula.link = res[[2]]$copula.link, Wmat = res[[2]]$Wmat)
        
        f = out1$lln * N + out2$lln * N
        dll1 = apply(out1$dll, 2, sum, na.rm = T)
        dll2 = apply(out2$dll, 2, sum, na.rm = T)
        g = c(dll1[1] + dll2[1], dll1[- 1], dll2[- 1])
        out.T1 = out1$ddlln * N; out.T2 = out2$ddlln * N
        out.gamma = out.T1[1, 1] + out.T2[1, 1]
        out.gamma.dLambda1 = out.T1[1, -1, drop = F]
        out.dLambda1.dLambda1 = out.T1[-1, -1]
        out.gamma.dLambda2 = out.T2[1, -1, drop = F]
        out.dLambda2.dLambda2 = out.T2[-1, -1]
        out.dLambda1.dLambda2 = matrix(0, nrow = n.tk1, ncol = n.tk2)
        B = rbind(
          cbind(out.gamma, out.gamma.dLambda1, out.gamma.dLambda2),
          cbind(t(out.gamma.dLambda1), out.dLambda1.dLambda1, out.dLambda1.dLambda2),
          cbind(t(out.gamma.dLambda2), t(out.dLambda1.dLambda2), out.dLambda2.dLambda2)
        )
      }
      return(list(value = f, gradient = g, hessian = B))
    }
    est.S1 = trust(objfun.S1, theta.ini, 5, 100, iterlim = 300,
                   minimize= FALSE, blather = T)
    theta.est = est.S1$argument
    gamma0.est = theta.est[1]; 
    dLambda1.est = theta.est[1 + c(1 : n.tk1)]
    dLambda2.est = theta.est[- c(1 : (n.tk1 + 1))]
    Lambda.X1 = sum.I(X1, ">=", tk1, dLambda1.est)
    uu1 = exp(- Lambda.X1)
    Lambda.X2 = sum.I(X2, ">=", tk2, dLambda2.est)
    uu2 = exp(- Lambda.X2)
    LambdaT.tt = sum.I(tt.all[[1]], ">=", tk1, dLambda1.est)
    surv.T1 = exp(- LambdaT.tt)
    LambdaT.tt = sum.I(tt.all[[2]], ">=", tk2, dLambda2.est)
    surv.T2 = exp(- LambdaT.tt)
    
    objfun.S2 = function(x){
      dat.trivic = prepare_trivic(theta = x, index12 = index12, 
                                  link12 = link12, Wmat12 = Wmat12, 
                                  control12 = control12, N = N, 
                                  thetaT = TRUE, datT = NULL)
      - lln.trivic.PMLE(
        uu1 = uu1, uu2 = uu2, uu3 = uD, 
        dd1 = delta1, dd2 = delta2, dd3 = deltaD, 
        par12 = dat.trivic$alphai, par1 = res[[1]]$alphai, par2 = res[[2]]$alphai,
        index12 = index12, index1 = res[[1]]$copula.index, 
        index2 = res[[2]]$copula.index, link12 = link12, 
        Wmat12 = Wmat12, yes.constraint = dat.trivic$yes.constraint, 
        only.ll = TRUE)$lln
    }
    gamma12.est = nlm(f = objfun.S2, p = gamma12.ini)$estimate
    tau12.est = tau12.est.fun(par1 = gamma0.est, par2 = gamma0.est, 
                              par12 = gamma12.est, index12 = 3, index1 = 3, 
                              index2 = 3)
    tau0.est = BiCopPar2Tau(family = 3, par = gamma0.est)
    
    out.vine.same = data.frame(
      para = c("gamma0",  "gamma12", "tau0", "tau12"), 
      est = c(gamma0.est, gamma12.est, tau0.est, tau12.est), 
      method = "vine_same",
      true = c(gamma0.true, gamma12.true, tau0.true, tau0.true))
    out.vine.same = rbind(
      out.vine.same,
      data.frame(
        para = "surv.T1", est = surv.T1, method = "vine_same", 
        true = 1 - pp.true
      ),
      data.frame(
        para = "surv.T2", est = surv.T2, method = "vine_same", 
        true = 1 - pp.true
      ))
    out = rbind(out.nested, out.vine.diff, out.vine.same) %>% mutate(rep = k)
    
    save(tt.dat, cen.rate, out, file = file.nm)
  }
  print(Sys.time() - t1)
}
stopCluster(cl)
