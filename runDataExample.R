# Setup ----
rm(list=ls())
library(tidyverse)
library(trust)
library(survival)
library(VineCopula)
library(doSNOW)
library(parallel)
library(foreach)
library(caret)
source("helpers.R")
source("helpers_bvic.R")
source("helpers_trivic.R")
Cop2Link = c(
  "Clayton" = "log", 
  "Gumbel" = "log-1", 
  "Frank" = "identity",
  "Gaussian" = "tanh",
  "fClayton" = "neglog",
  "fGumbel" = "neglog-1"
)
load("workdata.RData")

events = c("firstupdate", "firstResp", "first10Com", "end"); J = length(events)
times = paste0("time_", events)
deltas = paste0("status_", events)
Gfuns = rep("PH", J)
copula.families = list(
  Stage1 = c(firstupdate = "Gumbel", firstResp = "Gumbel", 
             first10Com = "Gumbel"),
  Stage2 = c("(firstupdate,first10Com)" = "Gumbel", 
             "(firstResp,first10Com)" = "Clayton", 
             "(firstupdate,firstResp)" = "Clayton"),
  Stage3 = c("Frank")
)

model.nm = "fullmodel"
cov.nms = c("email_is_verified", "yesResp_1st",
            "log_amount_raised_usd_1st", "log_narrative_length_1st",
            "log_numUpdate_1st", "log_perk_price_mean_1st", "log_comment_number_1st")
surv.fmlas = list(
  "firstupdate" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "firstResp" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "first10Com" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "end" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+")))
)
copula.fmlas = c()
copula.fmlas$Stage1 = list(
  "(firstupdate,end)" =  as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "(firstResp,end)" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "(first10Com,end)" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+")))
)
copula.fmlas$Stage2 = list(
  "(firstupdate,first10Com)" =  as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "(firstResp,first10Com)" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+"))),
  "(firstupdate,firstResp)" = as.formula(paste0(" ~ ", paste0(cov.nms, collapse = "+")))
)
copula.fmlas$Stage3 = ~ 1
data = workdata; initial = NULL

N = nrow(data); J = length(events)
summary.gamma.all = as.list(1 : 3)
summary.beta.all = as.list(1 : J)
theta.est.all  = uu.all = dd.all = Psi.theta.all = Psi.uu.all = as.list(1 : J)

# Stage I ----
print("=== Stage I ===")
## Marginal for D ----
event.nm = events[J]; timeD = times[J]; statusD = deltas[J]; 
fmlaD = update(surv.fmlas[[event.nm]], 
               as.formula(paste0("Surv(", timeD, ",", statusD, ") ~ .")))
zi = model.matrix(fmlaD, data)[, - 1, drop = F]
xi = data[, timeD]; di = data[, statusD]; dd.all[[J]] = di; deltaD = di
tk = sort(unique(xi[di == 1]))

myfitD = coxph(fmlaD, data = data, ties = "breslow", robust = TRUE)
b.est = myfitD$coef; n.b = length(b.est)
Lambda.fit = basehaz(myfitD, centered = FALSE)
Lambda = Lambda.fit$hazard
dLambda = c(Lambda[1], diff(Lambda, lag = 1))
tt = Lambda.fit$time[dLambda != 0]; n.tt = length(tt)
dLambda = dLambda[dLambda != 0]

Lambda.xi = sum.I(xi, ">=", tt, dLambda)  
bzi = as.vector(zi %*% matrix(b.est, byrow = F, ncol = 1))
ui = exp(- Lambda.xi * exp(bzi))
uu.all[[J]] = ui

xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) * di
out =  -  exp(bzi) * xi.g.tt +
  dN.tt / matrix(dLambda, byrow = T, nrow = length(xi), ncol = n.tt)
dll.b = zi * (di - Lambda.xi * exp(bzi))
dll = cbind(dll.b, out)
out = - diag(apply(dN.tt, 2, sum) / dLambda ^ 2)
dll.bb = - crossprod(
  Lambda.xi * exp(bzi) * zi, zi)
dll.dLambdab = - crossprod(exp(bzi) * xi.g.tt , zi)
ddlln = rbind(
  cbind(dll.bb, t(dll.dLambdab)),
  cbind(dll.dLambdab, out)
)
Imat = - ddlln / N; Imat.inv = solve(Imat)
Psi.theta = dll %*% Imat.inv 
theta.cov = crossprod(Psi.theta, Psi.theta) / (N ^ 2)
theta.se = sqrt(diag(theta.cov))
theta.est.all[[J]] = list(
  beta = b.est, tt = tt, dLambda = dLambda, cov = theta.cov, se = theta.se)
print(paste0("Estimation for mariginal of ", events[J]))
summary.beta = 
  data.frame(est = b.est, se = theta.se[1 : n.b]) %>%
  mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) %>%
  mutate(est = round(est, 4), se = round(se, 4),  z = round(z, 4))
print(summary.beta)
summary.beta.all[[J]] = summary.beta

Psi.theta.all[[J]] = Psi.theta
A = xi.g.tt * exp(bzi)
A = cbind(zi * Lambda.xi * exp(bzi), A) 
A = A * ui
Psi.uu.all[[J]] = - Psi.theta %*% t(A)

## Copula for (Tj, D) ----
gamma.S1.all = alpha.S1.all = cuu.all = Psi.cuu.all = 
  Psi.gamma.S1.all = summary.gamma.S1 = as.list(1 : (J - 1)) 
for (j in 1 : (J - 1)){
  event.nm = events[j]; timej = times[j]; statusj = deltas[j]; 
  fmlaT = update(surv.fmlas[[event.nm]], 
                 as.formula(paste0("Surv(", timej, ",", statusj, ") ~ .")))
  Gfunj = Gfuns[j]
  copula.fam = copula.families[[1]][event.nm]
  copula.index = MyCopIndex(copula.fam)  
  control = MyCop(copula.index)
  copula.link = MylinkFun(Cop2Link[copula.fam])
  copula.fmlaj = copula.fmlas[[1]][[paste0("(",event.nm,",end)")]]
  print(paste0(copula.fam, " copula for ", events[J], " and ", event.nm, "----"))
  
  Xi = data[, timej]; deltaT = data[, statusj]
  dd.all[[j]] = deltaT
  Wmat = model.matrix(copula.fmlaj, data); n.gamma = ncol(Wmat)
  ZmatT = model.matrix(fmlaT, data)[ , - 1, drop = F]; n.bT = ncol(ZmatT)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(tk)
  if (is.null(initial)){
    myfit = coxph(fmlaT, data)
    b.ini = myfit$coefficients
    Lambda.fit = basehaz(myfit, centered = FALSE)
    Lambda = Lambda.fit$hazard
    dLambda = c(Lambda[1], diff(Lambda, lag = 1))
    dLambda.ini = dLambda[dLambda != 0]
    thetaT.ini = c(b.ini, dLambda.ini)
    gamma.ini = copula.link$hinv.fun(BiCopTau2Par(copula.index, tau = 0.5))
    gamma.ini = c(gamma.ini, rep(0, n.gamma - 1))
    theta.ini = c(gamma.ini, thetaT.ini)
  } else {
    theta.ini = initial$theta[[j]]
  }
  objfun <- function(x){
    save(x, file = "x.RData")
    dat_bvic = prepare_bvic(
      theta = x, thetaT = FALSE, Xi = Xi, deltaT = deltaT, ZmatT = ZmatT, 
      copula.link = copula.link, Wmat = Wmat, GfunT = Gfunj, control = control)
    out = llfuns.bvic(
      u1 = dat_bvic$ui, u2 = uu.all[[J]], d1 = deltaT, d2 = dd.all[[J]], 
      copula.index = copula.index, alphai = dat_bvic$alphai, 
      yes.constraint = dat_bvic$yes.constraint,
      only.ll = FALSE, theta1 = FALSE, dat_bvic = dat_bvic,
      copula.link = copula.link, Wmat = Wmat) 
    f = out$lln
    if (all(is.na(out$dll))) {f=-Inf; g = NA} else {
      if (!is.null(dim(out$dll))) g = apply(out$dll, 2, mean, na.rm = T)
    }
    if (any(is.na(out$ddlln))) {f = -Inf; B = NA} else B = out$ddlln
    if (all(is.numeric(out$ddlln))) B = out$ddlln else {f = -Inf; B = NA}
    if (all(is.finite(out$ddlln))) B = out$ddlln else {f = -Inf; B = NA}
    print(f)
    list(value = f, gradient = g, hessian = B)
  }
  
  time.comp = system.time({
    est.res <- trust(objfun, theta.ini, 5, 100, iterlim = 300,
                     minimize= FALSE, blather = T)
  })
  theta.est = est.res$argument
  if (length(theta.est) != sum(n.gamma + n.bT + n.tk)) 
    stop("parameter dimension does not match")
  dat_bvic = prepare_bvic(theta = theta.est, thetaT = FALSE, Xi = Xi, 
                          deltaT = deltaT, ZmatT = ZmatT, 
                          copula.link = copula.link, 
                          Wmat = Wmat, GfunT = Gfunj, control = control)
  uu = dat_bvic$ui
  if (copula.fam == "Gumbel" & sum(uu == 1) > 0) uu[uu == 1] = exp(- 1 / N)
  dat_bvic$ui = uu
  gamma.est = theta.est[1 : n.gamma]
  names(gamma.est) = colnames(Wmat)
  beta.est = theta.est[c(1 : n.bT) + n.gamma]
  dLambda.est = theta.est[c(1 : n.tk) + n.bT + n.gamma]
  alpha.S1.all[[j]] = dat_bvic$alphai
  uu.all[[j]] = dat_bvic$ui
  cuu.all[[j]] = BiCopHfunc2(dat_bvic$ui, uu.all[[J]], 
                             family = copula.index, par = dat_bvic$alphai)
  lls = llfuns.bvic(
    u1 = dat_bvic$ui, u2 = uu.all[[J]], d1 = deltaT, d2 = dd.all[[J]], 
    copula.index = copula.index, alphai = dat_bvic$alphai, 
    yes.constraint = dat_bvic$yes.constraint,
    only.ll = FALSE, theta1 = FALSE, dat_bvic = dat_bvic,
    copula.link = copula.link, Wmat = Wmat, yes.dllDu2 = TRUE)
  Imat = - lls$ddlln; Imat.inv = solve(Imat)
  dll.k = lls$dll.u2; dll.k[is.na(dll.k)] = 0
  Psi =  lls$dll; Psi[is.na(Psi)] = 0
  Psi.theta = (Psi + Psi.uu.all[[J]] %*% dll.k / N) %*% Imat.inv
  theta.cov = crossprod(Psi.theta, Psi.theta) / (N ^ 2)
  theta.se = sqrt(diag(theta.cov))
  
  gamma.se = theta.se[1 : n.gamma]; names(gamma.se) = colnames(Wmat)
  gamma.cov = theta.cov[c(1 : n.gamma), c(1 : n.gamma)]
  beta.se = theta.se[c(1 : n.bT) + n.gamma]
  dLambda.se = theta.se[c(1 : n.tk) + n.bT + n.gamma]
  
  Psi.gamma = Psi.theta[, c (1: n.gamma), drop = F]
  Psi.gamma.S1.all[[j]] = Psi.gamma
  
  thetaT.est = theta.est[- c(1 : n.gamma)]
  thetaT.se = theta.se[- c(1 : n.gamma)]
  thetaT.cov = theta.cov[- c(1 : n.gamma), - c(1 : n.gamma)]
  Psi.theta.all[[j]]  = Psi.theta[, - c (1: n.gamma), drop = F]
  print(paste0("Estimation for mariginal of ", event.nm))
  summary.beta = 
    data.frame(est = beta.est, se = beta.se) %>%
    mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) %>%
    mutate(est = round(est, 4), se = round(se, 4),  z = round(z, 4))
  print(summary.beta)
  summary.beta.all[[j]] = summary.beta
  
  print(paste0("Estimation for parameter  for copula between", event.nm, 
               " and ", events[J]))
  summary.gamma = 
    data.frame(est = gamma.est, se = gamma.se) %>%
    mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) %>%
    mutate(est = round(est, 4), se = round(se, 4),  z = round(z, 4))
  print(summary.gamma)
  summary.gamma.S1[[j]] = summary.gamma
  
  ebz = exp(as.vector(ZmatT %*% beta.est))
  LambdaT.Xi = sum.I(Xi, ">=", tk, dLambda.est)
  Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T,nrow = N, ncol = n.tk))
  A = cbind(ZmatT * LambdaT.Xi * ebz, Xi.g.tk * ebz) * dat_bvic$ui
  Psi.uu.all[[j]] = - Psi.theta.all[[j]] %*% t(A)
  
  cuu.duT = BiCopPDF(dat_bvic$ui, uu.all[[J]], family = copula.index, 
                     par = dat_bvic$alphai)
  cuu.duD = BiCopHfuncDeriv(dat_bvic$ui, uu.all[[J]], family = copula.index, 
                            par = dat_bvic$alphai, deriv = "u2")
  cuu.dgamma = BiCopHfuncDeriv(dat_bvic$ui, uu.all[[J]], family = copula.index, 
                               par = dat_bvic$alphai, deriv = "par") * 
    copula.link$dot.h.fun(dat_bvic$alphai.lp) * Wmat
  
  Psi.cuu.all[[j]] = t(t(Psi.uu.all[[j]]) * cuu.duT) + 
    t(t(Psi.uu.all[[J]]) * cuu.duD) + 
    Psi.gamma %*% t(cuu.dgamma)
  
  theta.est.all[[j]] = list(beta = beta.est, tt = tk, dLambda = dLambda.est, 
                            cov = thetaT.cov, se = thetaT.se)
  gamma.S1.all[[j]] = list(
    gamma = gamma.est, cov = gamma.cov, se = gamma.se, 
    convergence = est.res$converged, iterations = est.res$iterations
  )
  
  print(time.comp)
}
summary.gamma.all$stage1 = summary.gamma.S1

# Stage II ----
print("=== Stage II ===")

## (T1,T3|D) & (T2,T3|D) ----
gamma.S2.all = index.S2.all = alpha.S2.all = Psi.gamma.S2.all = 
  Wmat.all = copula.link.all = alpha.lp.S2.all = summary.gamma.S2 = 
  as.list(1 : (J - 2)) 
for (j in 1 : (J - 2)){
  pair.nm = paste0("(", events[j], ",", events[J-1],")")
  cu1 = cuu.all[[j]]; cu2 = cuu.all[[J - 1]]
  d1 = dd.all[[j]]; d2 = dd.all[[J - 1]]
  copula.fam = copula.families[[2]][pair.nm]
  copula.fmla = copula.fmlas[[2]][[pair.nm]]
  copula.index = MyCopIndex(copula.fam)  
  index.S2.all[[j]] = copula.index
  control = MyCop(copula.index)
  copula.link = MylinkFun(Cop2Link[copula.fam])
  copula.link.all[[j]] = copula.link
  h.fun = copula.link$h.fun
  Wmat = model.matrix(copula.fmla, data); n.gamma = ncol(Wmat)
  Wmat.all[[j]] = Wmat
  print(paste0(copula.fam, " copula for ", events[j], " and ", events[J-1], 
               " given ", events[J], "----"))
  print("with model for copula parameter")
  print(copula.fmla)
  
  if (is.null(initial)){
    gamma.ini = copula.link$hinv.fun(BiCopTau2Par(copula.index, tau = 0.5))
    gamma.ini = c(gamma.ini, rep(0, n.gamma - 1))
  } else {
    gamma.ini = initial$gamma.S2[[j]]
  }
  
  objfun <- function(x){
    save(x, file = "x.RData")
    copula.lp = as.vector(Wmat %*% matrix(x, byrow = F, ncol = 1))
    alphai = h.fun(copula.lp)
    yes.constraint = min(alphai) <= control$lwr | max(alphai) >= control$upr
    if (yes.constraint) {
      return(list(value = -Inf, gradient = NA, hessian = NA))
    } else{
      N = length(cu1)
      if (copula.index == 4) { # gumbel
        cu1[cu1 == 1] = exp(- 1 / N)
        cu2[cu2 == 1] = exp(- 1 / N)
      }
      dd = funsData_bvic(density_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
                         copula.index = copula.index, para = alphai)
      ll = log(dd); lln = mean(ll, na.rm = T)
      alphai.lp = copula.link$hinv.fun(alphai); 
      dd.Dpar = funsData_bvic(
        Funs = density.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
        copula.index = copula.index, para = alphai)
      ll.Dpar = (dd.Dpar / dd) 
      ll.Dgamma = ll.Dpar *  copula.link$dot.h.fun(alphai.lp) * Wmat
      dll = ll.Dgamma
      
      dd.Dpar.Dpar = funsData_bvic(
        Funs = density.Dpar.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
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
      f = lln
      if (all(is.na(dll))) g = NA else {
        if (!is.null(dim(dll))) g = apply(dll, 2, mean, na.rm = T)
      }
      B = ddlln
      print(f)
      return(list(value = f, gradient = g, hessian = B))
    }
    
  }
  
  time.comp = system.time({
    est.res <- trust(objfun, gamma.ini, 5, 100, iterlim = 300,
                     minimize= FALSE, blather = T)
  })
  gamma.est = est.res$argument; names(gamma.est) = colnames(Wmat)
  alphai.lp = as.vector(Wmat %*% matrix(gamma.est, byrow = F, ncol = 1))
  alpha.lp.S2.all[[j]] = alphai.lp
  alphai = copula.link$h.fun(alphai.lp)
  alpha.S2.all[[j]] = alphai
  
  dd = funsData_bvic(density_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
                     copula.index = copula.index, para = alphai)
  dd.Dpar = funsData_bvic(
    Funs = density.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
    copula.index = copula.index, para = alphai)
  ll.Dpar = dd.Dpar / dd
  dd.Dpar.Dpar = funsData_bvic(
    Funs = density.Dpar.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
    copula.index = copula.index, para = alphai)
  ll.Dpar.Dpar = dd.Dpar.Dpar / dd - ll.Dpar ^ 2 
  ll.Dgamma.Dgamma = matrix(
    apply((ll.Dpar.Dpar * (copula.link$dot.h.fun(alphai.lp) ^ 2) + 
             ll.Dpar * copula.link$ddot.h.fun(alphai.lp)) * 
            Wmat[, rep(1 : n.gamma, n.gamma), drop = F] * 
            Wmat[, rep(1 : n.gamma, each = n.gamma), drop = F], 
          2, mean, na.rm = T),
    byrow = T, nrow = n.gamma, ncol = n.gamma)
  Imat = - ll.Dgamma.Dgamma; Imat.inv = solve(Imat)
  
  ll.Dpar.Du1 = funsData_bvic(
    Funs = density.Dpar.Du1_bvic, 
    d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
    copula.index = copula.index, para = alphai) / dd
  ll.Du1 = funsData_bvic(
    Funs = density.Du1_bvic, 
    d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
    copula.index = copula.index, para = alphai) / dd
  dll.k = ll.Dpar.Du1 - ll.Du1 * ll.Dpar
  dll.k[is.na(dll.k)] = 0
  dll.k = dll.k * copula.link$dot.h.fun(alphai.lp) * Wmat
  tmp = Psi.cuu.all[[j]]
  tmp[is.na(tmp)] = 0
  Psi.dll.u1 = tmp %*% dll.k / N
  
  ll.Dpar.Du2 = funsData_bvic(
    Funs = density.Dpar.Du2_bvic, 
    d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
    copula.index = copula.index, para = alphai) / dd
  ll.Du2 = funsData_bvic(
    Funs = density.Du2_bvic, 
    d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
    copula.index = copula.index, para = alphai) / dd
  dll.k = ll.Dpar.Du2 - ll.Du2 * ll.Dpar
  dll.k[is.na(dll.k)] = 0
  dll.k = dll.k * copula.link$dot.h.fun(alphai.lp) * Wmat
  Psi.dll.u2 = Psi.cuu.all[[J - 1]] %*% dll.k / N
  
  Psi.par = ll.Dpar * copula.link$dot.h.fun(alphai.lp) * Wmat
  Psi.gamma.S2.all[[j]] = (Psi.par +  Psi.dll.u1 + Psi.dll.u2) %*% Imat.inv
  
  gamma.cov = crossprod(Psi.gamma.S2.all[[j]], Psi.gamma.S2.all[[j]]) / (N ^ 2)
  gamma.se = sqrt(diag(gamma.cov))
  names(gamma.se) = colnames(Wmat)
  summary.gamma = 
    data.frame(est = gamma.est, se = gamma.se) %>%
    mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) %>%
    mutate(est = round(est, 4), se = round(se, 4),  z = round(z, 4))
  print(summary.gamma)
  summary.gamma.S2[[j]] = summary.gamma
  
  gamma.S2.all[[j]] = 
    list(gamma = gamma.est, cov = gamma.cov, se = gamma.se,
         convergence = est.res$converged,
         iterations = est.res$iterations)
  
  print(time.comp)
}

## (T1,T2|D) ----
j = 3
pair.nm = paste0("(", events[1], ",", events[2],")")
cu1 = cuu.all[[1]]; cu2 = cuu.all[[2]]
d1 = dd.all[[1]]; d2 = dd.all[[2]]
copula.fam = copula.families[[2]][pair.nm]
copula.fmla = copula.fmlas[[2]][[pair.nm]]
copula.index = MyCopIndex(copula.fam)  
index.S2.all[[j]] = copula.index
control = MyCop(copula.index)
copula.link = MylinkFun(Cop2Link[copula.fam])
copula.link.all[[j]] = copula.link
h.fun = copula.link$h.fun
Wmat = model.matrix(copula.fmla, data); n.gamma = ncol(Wmat)
Wmat.all[[j]] = Wmat
print(paste0(copula.fam, " copula for ", events[1], " and ", events[2], 
             " given ", events[4], "----"))
print("with model for copula parameter")
print(copula.fmla)

if (is.null(initial)){
  gamma.ini = copula.link$hinv.fun(BiCopTau2Par(copula.index, tau = 0.5))
  gamma.ini = c(gamma.ini, rep(0, n.gamma - 1))
} else {
  gamma.ini = initial$gamma.S2[[j]]
}

objfun <- function(x){
  save(x, file = "x.RData")
  copula.lp = as.vector(Wmat %*% matrix(x, byrow = F, ncol = 1))
  alphai = h.fun(copula.lp)
  yes.constraint = min(alphai) <= control$lwr | max(alphai) >= control$upr
  if (yes.constraint) {
    return(list(value = -Inf, gradient = NA, hessian = NA))
  } else{
    N = length(cu1)
    if (copula.index == 4) { # gumbel
      cu1[cu1 == 1] = exp(- 1 / N)
      cu2[cu2 == 1] = exp(- 1 / N)
    }
    dd = funsData_bvic(density_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
                       copula.index = copula.index, para = alphai)
    ll = log(dd); lln = mean(ll, na.rm = T)
    alphai.lp = copula.link$hinv.fun(alphai); 
    dd.Dpar = funsData_bvic(
      Funs = density.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
      copula.index = copula.index, para = alphai)
    ll.Dpar = (dd.Dpar / dd) 
    ll.Dgamma = ll.Dpar *  copula.link$dot.h.fun(alphai.lp) * Wmat
    dll = ll.Dgamma
    
    dd.Dpar.Dpar = funsData_bvic(
      Funs = density.Dpar.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
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
    f = lln
    if (all(is.na(dll))) g = NA else {
      if (!is.null(dim(dll))) g = apply(dll, 2, mean, na.rm = T)
    }
    B = ddlln
    print(f)
    return(list(value = f, gradient = g, hessian = B))
  }
  
}

time.comp = system.time({
  est.res <- trust(objfun, gamma.ini, 5, 100, iterlim = 300,
                   minimize= FALSE, blather = T)
})
gamma.est = est.res$argument; names(gamma.est) = colnames(Wmat)
alphai.lp = as.vector(Wmat %*% matrix(gamma.est, byrow = F, ncol = 1))
alpha.lp.S2.all[[j]] = alphai.lp
alphai = copula.link$h.fun(alphai.lp)
alpha.S2.all[[j]] = alphai

dd = funsData_bvic(density_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
                   copula.index = copula.index, para = alphai)
dd.Dpar = funsData_bvic(
  Funs = density.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
  copula.index = copula.index, para = alphai)
ll.Dpar = dd.Dpar / dd
dd.Dpar.Dpar = funsData_bvic(
  Funs = density.Dpar.Dpar_bvic, d1 = d1, d2 = d2, u1 = cu1, u2 = cu2,
  copula.index = copula.index, para = alphai)
ll.Dpar.Dpar = dd.Dpar.Dpar / dd - ll.Dpar ^ 2 
ll.Dgamma.Dgamma = matrix(
  apply((ll.Dpar.Dpar * (copula.link$dot.h.fun(alphai.lp) ^ 2) + 
           ll.Dpar * copula.link$ddot.h.fun(alphai.lp)) * 
          Wmat[, rep(1 : n.gamma, n.gamma), drop = F] * 
          Wmat[, rep(1 : n.gamma, each = n.gamma), drop = F], 
        2, mean, na.rm = T),
  byrow = T, nrow = n.gamma, ncol = n.gamma)
Imat = - ll.Dgamma.Dgamma; Imat.inv = solve(Imat)

ll.Dpar.Du1 = funsData_bvic(
  Funs = density.Dpar.Du1_bvic, 
  d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
  copula.index = copula.index, para = alphai) / dd
ll.Du1 = funsData_bvic(
  Funs = density.Du1_bvic, 
  d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
  copula.index = copula.index, para = alphai) / dd
dll.k = ll.Dpar.Du1 - ll.Du1 * ll.Dpar
dll.k[is.na(dll.k)] = 0
dll.k = dll.k * copula.link$dot.h.fun(alphai.lp) * Wmat
Psi.dll.u1 = Psi.cuu.all[[j]] %*% dll.k / N

ll.Dpar.Du2 = funsData_bvic(
  Funs = density.Dpar.Du2_bvic, 
  d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
  copula.index = copula.index, para = alphai) / dd
ll.Du2 = funsData_bvic(
  Funs = density.Du2_bvic, 
  d1 = d1, d2 = d2, u1 = cu1, u2 = cu2, 
  copula.index = copula.index, para = alphai) / dd
dll.k = ll.Dpar.Du2 - ll.Du2 * ll.Dpar
dll.k[is.na(dll.k)] = 0
dll.k = dll.k * copula.link$dot.h.fun(alphai.lp) * Wmat
Psi.dll.u2 = Psi.cuu.all[[J - 1]] %*% dll.k / N

Psi.par = ll.Dpar * copula.link$dot.h.fun(alphai.lp) * Wmat
Psi.gamma.S2.all[[j]] = (Psi.par +  Psi.dll.u1 + Psi.dll.u2) %*% Imat.inv

gamma.cov = crossprod(Psi.gamma.S2.all[[j]], Psi.gamma.S2.all[[j]]) / (N ^ 2)
gamma.se = sqrt(diag(gamma.cov))
names(gamma.se) = colnames(Wmat)
summary.gamma = 
  data.frame(est = gamma.est, se = gamma.se) %>%
  mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) %>%
  mutate(est = round(est, 4), se = round(se, 4),  z = round(z, 4))
print(summary.gamma)
summary.gamma.S2[[j]] = summary.gamma
gamma.S2.all[[j]] = 
  list(gamma = gamma.est, cov = gamma.cov, se = gamma.se,
       convergence = est.res$converged,
       iterations = est.res$iterations)

print(time.comp)
summary.gamma.all$stage2 = summary.gamma.S2

# Stage III ----
print("=== Stage III ===")
par1 = alpha.S2.all[[1]]; par2 = alpha.S2.all[[2]]
index1 = index.S2.all[[1]]; index2 = index.S2.all[[2]]
copula.fmla = copula.fmlas[[3]]
Wmat12 = model.matrix(copula.fmla, data); n.gamma = ncol(Wmat12)
uu1 = cuu.all[[1]]; uu2 = cuu.all[[2]]; uu3 = cuu.all[[3]]
dd1 = dd.all[[1]]; dd2 = dd.all[[2]]; dd3 = dd.all[[3]]
copula.fam = fam12
index12 = MyCopIndex(copula.fam)
control12 = MyCop(index12)
link12 = MylinkFun("identity")
print(paste0(copula.fam, " copula for ", events[1], " and ", events[2],
             " given ", events[J - 1], " and ", events[J],"----"))
print("with model for copula parameter")
print(copula.fmla)

alpha.ini = BiCopTau2Par(index12, tau = 0.5)
gamma.ini = c(alpha.ini, rep(0, n.gamma - 1))

objfun = function(x){
  save(x, file = "x.RData")
  par12 = link12$h.fun(
    as.vector(Wmat12 %*% matrix(x, byrow = F, ncol = 1))
  )
  yes.constraint = min(par12) < control12$lwr |
    max(par12) > control12$upr
  out = lln.trivic.PMLE(uu1, uu2, uu3, dd1, dd2, dd3, par12, par1, par2,
                        index12, index1, index2, link12, Wmat12,
                        yes.constraint, only.ll = TRUE,
                        pc.method = "foreach")
  f = out$lln
  return(-f)
  print(f)
}
time.comp = system.time({
  est.res <- nlm(objfun, p = gamma.ini)
})
gamma.est = gamma.PMLE.all$estimate
alpha12.lp = as.vector(Wmat12 %*% matrix(gamma.est, byrow = F, ncol = 1))
alpha12 = link12$h.fun(alpha12.lp)
gamma.S3 = list(gamma = gamma.est)

yes.constraint = min(alpha12) < control12$lwr |
  max(alpha12) > control12$upr
t1 = Sys.time()
lls = lln.trivic.PMLE(
  uu1 = uu1, uu2 = uu2, uu3 = uu3,
  dd1 = dd1, dd2 = dd2, dd3 = dd3,
  par12 = alpha12, par1 = par1, par2 = par2,
  index12 = index12, index1 = index1,
  index2 = index2, link12 = link12,
  Wmat12 = Wmat12, yes.constraint = yes.constraint,
  only.ll =  FALSE, pc.method = "foreach")
print(Sys.time() - t1)
print("lls done")
t1 = Sys.time()
dlls = dll.trivic.PMLE(
  uu1 = uu1, uu2 = uu2, uu3 = uu3,
  dd1 = dd1, dd2 = dd2, dd3 = dd3,
  par12 = alpha12, par1 = par1, par2 = par2,
  index12 = index12, index1 = index1,
  index2 = index2, yes.constraint = yes.constraint,
  pc.method = "foreach")
print(Sys.time() - t1)
print("dlls done")

Imat = - lls$ddlln

dll = dlls$ll.Dpar12.Du1
dll[is.na(dll)] = 0
dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
Psi.u1 = Psi.cuu.all[[1]] %*% dll

dll = dlls$ll.Dpar12.Du2
dll[is.na(dll)] = 0
dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
tmp = Psi.cuu.all[[2]]
tmp[is.na(tmp)] = 0
Psi.u2 = tmp %*% dll

dll = dlls$ll.Dpar12.Du3
dll[is.na(dll)] = 0
dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
Psi.u3 = Psi.cuu.all[[3]] %*% dll

dll = dlls$ll.Dpar12.Dpar1
dll[is.na(dll)] = 0
dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
Psi.par1 = apply(dll, 2, function(x){
  apply(Psi.gamma.S2.all[[1]] *
          copula.link.all[[1]]$dot.h.fun(alpha.lp.S2.all[[1]]) *
          Wmat.all[[1]] * x,
        1, sum)
})

dll = dlls$ll.Dpar12.Dpar2
dll[is.na(dll)] = 0
dll = dll * link12$dot.h.fun(alpha12.lp) * Wmat12
Psi.par2 = apply(dll, 2, function(x){
  apply(Psi.gamma.S2.all[[2]] *
          copula.link.all[[2]]$dot.h.fun(alpha.lp.S2.all[[2]]) *
          Wmat.all[[2]] * x,
        1, sum)
})
Psi =  lls$dll
Psi[is.na(Psi)] = 0
Psi.gamma = (Psi +
               (Psi.u1 + Psi.u2 + Psi.u3 + Psi.par1 + Psi.par2) / N) %*% solve(Imat)
gamma.cov = crossprod(Psi.gamma, Psi.gamma) / N ^ 2
gamma.se = sqrt(diag(gamma.cov))
summary.gamma =
  data.frame(est = gamma.est, se = gamma.se) %>%
  mutate(z = est /se, pval = 2 * pnorm(abs(z), lower = F)) 
print(summary.gamma)
summary.gamma.all$stage3 = summary.gamma 
gamma.S3 = list(
  gamma = gamma.est, cov = gamma.cov, se = gamma.se,
  convergence = est.res$converged, iterations = est.res$iterations)

# Save results -----
save(Psi.theta.all,file = "fullmodel_Psi_theta.RData")
save(Psi.gamma.S1.all,file = "fullmodel_Psi_gamma_S1.RData")
save(Psi.gamma.S2.all,file = "fullmodel_Psi_gamma_S2.RData")
Psi.gamma.S3 = Psi.gamma
save(Psi.gamma.S3, file = "fullmodel_Psi_gamma_S3.RData")

res.fullmodel = list(
  theta = theta.est.all,
  gamma.S1 = gamma.S1.all, gamma.S2 = gamma.S2.all, gamma.S3 = gamma.S3
)
para.est = c()
para.est$S1 = as.list(1 : J)
names(para.est$S1) = events
for (j in 1 : (J-1)){
  gamma.est = res.fullmodel$gamma.S1[[j]]$gamma
  beta.est = res.fullmodel$theta[[j]]$beta
  dLambda.est = res.fullmodel$theta[[j]]$dLambda
  tt = res.fullmodel$theta[[j]]$tt
  para.est$S1[[j]] = list(gamma = gamma.est, 
                          beta = beta.est, dLambda = dLambda.est, tt = tt)
}
para.est$S1$end = list(beta = res.fullmodel$theta[["end"]]$beta, 
                       dLambda = res.fullmodel$theta[["end"]]$dLambda)
para.est$S2 = as.list(1 : 3)
names(para.est$S2) = names(copula.families$Stage2)
for (j in 1 : 3){
  para.est$S2[[j]] = res.fullmodel$gamma.S2[[j]]$gamma
}
para.est$S3 = res.fullmodel$gamma.S3$gamma

save(res.fullmodel, summary.beta.all, summary.gamma.all,
     para.est, events, times, deltas, surv.fmlas, copula.families, 
     copula.fmlas, data, file = "res.fullmodel.RData")

stopCluster(cl)