# Source Function for pMLE: clayton ####

#' Function to return estimated SD, S1star, S2star
#'
#' Revised from the function of the same name in the `bvic` package
#' (https://github.com/dli-stats/bvic), provided in:
#' Li, D., Hu, X. J., & Wang, R. (2023). Evaluating Association Between
#' Two Event Times with Observations Subject to Informative Censoring.
#' Journal of the American Statistical Association, 118(542), 1282-1294.
#' https://doi.org/10.1080/01621459.2021.1990766
#'
#' @importFrom survival survfit Surv
#' @importFrom stats approxfun
#' @keywords internal
npest.star_12 = function(X) {
  k = ncol(X)

  X$eta_1 = ifelse(X$delta_1==1 | X$delta_D==1, 1, 0)
  X$eta_2 = ifelse(X$delta_2==1 | X$delta_D==1, 1, 0)

  fit.SD = survfit(Surv(X$C, X$delta_D) ~ 1)
  SD.hat.fn = approxfun(x = fit.SD$time, y = fit.SD$surv, yleft = 1, yright = 0, method = "constant")
  SD.hat = SD.hat.fn(X$C)

  fit.S1.star = survfit(Surv(X$U1, X$eta_1) ~ 1)
  fit.S2.star = survfit(Surv(X$U2, X$eta_2) ~ 1)
  S1.star.fn = approxfun(x = fit.S1.star$time, y = fit.S1.star$surv, yleft = 1, yright = 0, method = "constant")
  S2.star.fn = approxfun(x = fit.S2.star$time, y = fit.S2.star$surv, yleft = 1, yright = 0, method = "constant")

  S1.star.hat = S1.star.fn(X$U1)
  S2.star.hat = S2.star.fn(X$U2)

  # Calculate SD.hat at values of U1/U2
  SD.U1.hat = SD.hat.fn(X$U1)
  SD.U2.hat = SD.hat.fn(X$U2)

  X$S1.star.hat = S1.star.hat
  X$S2.star.hat = S2.star.hat
  X$SD.hat = SD.hat
  X$SD.U1.hat = SD.U1.hat
  X$SD.U2.hat = SD.U2.hat

  return(list(
    X = X, SD.hat.fn = SD.hat.fn,
    S1.star.fn = S1.star.fn, S2.star.fn = S2.star.fn
  ))
}

clayton.cc.stage.T1T2.D.lik = function(X0, theta, S1.hat, S2.hat, SD.hat) {

  X = data.frame(X0, theta)
  X$S1 = S1.hat
  X$S2 = S2.hat
  X$SD = SD.hat

  case1.1 = X[which(X$delta_1 == 1 & X$delta_2 == 1 & X$delta_D == 1), ]
  if (nrow(case1.1) != 0) {
    L.part1.1 = with(case1.1, clayton.C21(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C.1(S1, S2, theta) *
                       clayton.C1.(S1, S2, theta) +
                       clayton.C11(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C11(S1, S2, theta))
  } else {
    L.part1.1 = 0
  }

  case1.2 = X[which(X$delta_1 == 1 & X$delta_2 == 0 & X$delta_D == 1), ]
  if (nrow(case1.2) != 0) {
    L.part1.2 = with(case1.2, clayton.C11(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C1.(S1, S2, theta))
  } else {
    L.part1.2 = 0
  }

  case1.3 = X[which(X$delta_1 == 0 & X$delta_2 == 1 & X$delta_D == 1), ]
  if (nrow(case1.3) != 0) {
    L.part1.3 = with(case1.3, clayton.C11(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C.1(S1, S2, theta))
  } else {
    L.part1.3 = 0
  }

  case1.4 = X[which(X$delta_1 == 0 & X$delta_2 == 0 & X$delta_D == 1), ]
  if (nrow(case1.4) != 0) {
    L.part1.4 = with(case1.4, clayton.C.1(clayton.C..(S1, S2, theta), SD, theta))
  } else {
    L.part1.4 = 0
  }

  case2.1 = X[which(X$delta_1 == 1 & X$delta_2 == 1 & X$delta_D == 0), ]
  if (nrow(case2.1) != 0) {
    L.part2.1 = with(case2.1, clayton.C2.(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C.1(S1, S2, theta) *
                       clayton.C1.(S1, S2, theta) +
                       clayton.C1.(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C11(S1, S2, theta))
  } else {
    L.part2.1 = 0
  }

  case2.2 = X[which(X$delta_1 == 1 & X$delta_2 == 0 & X$delta_D == 0), ]
  if (nrow(case2.2) != 0) {
    L.part2.2 = with(case2.2, clayton.C1.(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C1.(S1, S2, theta))
  } else {
    L.part2.2 = 0
  }

  case2.3 = X[which(X$delta_1 == 0 & X$delta_2 == 1 & X$delta_D == 0), ]
  if (nrow(case2.3) != 0) {
    L.part2.3 = with(case2.3, clayton.C1.(clayton.C..(S1, S2, theta), SD, theta) *
                       clayton.C.1(S1, S2, theta))
  } else {
    L.part2.3 = 0
  }

  case2.4 = X[which(X$delta_1 == 0 & X$delta_2 == 0 & X$delta_D == 0), ]
  if (nrow(case2.4) != 0) {
    L.part2.4 = with(case2.4, clayton.C..(clayton.C..(S1, S2, theta), SD, theta))
  } else {
    L.part2.4 = 0
  }

  L.result = c(
    L.part1.1, L.part1.2, L.part1.3, L.part1.4,
    L.part2.1, L.part2.2, L.part2.3, L.part2.4
  )

  suppressWarnings(-sum(log(L.result[L.result != 0]), na.rm = TRUE))
}


cdfExpr.Clayton = function(n) {
  expr = "u1^(-theta) - 1"
  for (i in 2:n) {
    cur = paste0("u", i, "^(-theta) - 1")
    expr = paste(expr, cur, sep = " + ")
  }
  expr = paste("(1 + (", expr, "))^ (-1/theta)")
  parse(text = expr)
}

pdfExpr = function(cdf, n) {
  val = cdf
  for (i in 1:n) {
    val = D(val, paste0("u", i))
  }
  val
}


C.clayton.expr = cdfExpr.Clayton(2)
clayton.C1. = function(u1, u2, theta) {}
body(clayton.C1.) = D(C.clayton.expr, "u1")
clayton.C2. = function(u1, u2, theta) {}
body(clayton.C2.) = D(D(C.clayton.expr, "u1"), "u1")
clayton.C.1 = function(u1, u2, theta) {}
body(clayton.C.1) = D(C.clayton.expr, "u2")
clayton.C11 = function(u1, u2, theta) {}
body(clayton.C11) = D(D(C.clayton.expr, "u1"), "u2")
clayton.C21 = function(u1, u2, theta) {}
body(clayton.C21) = D(D(D(C.clayton.expr, "u1"), "u1"), "u2")
clayton.C.. = function(u1, u2, theta) {}
body(clayton.C..) = C.clayton.expr


g.fn = function(w1, w2, theta) {
  (w1^(-theta) - w2^(-theta) + 1)^(-1 / theta)
}

g.expression = expression((w1^(-theta) - w2^(-theta) + 1)^(-1 / theta))

est_cop_par = function(copula = "clayton", data, niter = 50) {

  mydat = data
  n = nrow(mydat)
  out = npest.star_12(mydat)
  mydat1 = out$X

  iter = 1
  K = niter
  procC.iter.clayton.c.theta = matrix(NA, nrow=K, ncol=1)
  procC.iter.clayton.c.tau = matrix(NA, nrow=K, ncol=1)
  procC.iter.clayton.c.theta[iter,1] = 1
  procC.iter.clayton.c.tau[iter,1] = procC.iter.clayton.c.theta[iter,1] / (procC.iter.clayton.c.theta[iter,1]+2)


  # print(paste(iter, "run:", "theta0:", procC.iter.clayton.c.theta0[iter,1], ";theta:", procC.iter.clayton.c.theta[iter,1] ))

  repeat{
    iter = iter + 1
    if (iter == K) {
      break
    }

    procC.iter.clayton.c.theta[iter,1] = procC.iter.clayton.c.theta[(iter-1), 1]
    procC.iter.clayton.c.tau[iter,1] = procC.iter.clayton.c.tau[(iter-1), 1]

    c.S1.hat = g.fn(mydat1$S1.star.hat, mydat1$SD.U1.hat, procC.iter.clayton.c.theta[iter,1])
    c.S2.hat = g.fn(mydat1$S2.star.hat, mydat1$SD.U2.hat, procC.iter.clayton.c.theta[iter,1])

    out = optim(procC.iter.clayton.c.theta[iter,1], lower = 0, upper = 28,
                clayton.cc.stage.T1T2.D.lik, method = "Brent",
                X0 = mydat1,
                S1.hat = c.S1.hat, S2.hat = c.S2.hat,
                SD.hat = mydat1$SD.hat
    )

    procC.iter.clayton.c.theta[iter,1] = out$par[1]


    procC.iter.clayton.c.tau[iter,1] = procC.iter.clayton.c.theta[iter,1] / (procC.iter.clayton.c.theta[iter,1] + 2)

    if (abs(procC.iter.clayton.c.theta[iter,1] - procC.iter.clayton.c.theta[(iter-1),1]) < 1e-5) {
      break
    }
  }
  return(data.frame(
    theta = procC.iter.clayton.c.theta[max(which(!is.na(procC.iter.clayton.c.theta[,1]))),1]
  ))
}
