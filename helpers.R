sum.I <- function(yy, FUN, Yi, Vi=NULL){
  if (FUN == "<" | FUN == ">=") {yy = -yy; Yi = -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy, Yi),ties.method = 'f')[1 : length(yy)] - 
    rank(yy, ties.method = 'f')    
  if (substring(FUN, 2, 2) == "=") pos = length(Yi) - pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN, 2, 2) == "=") tmpind = order(-Yi) else  tmpind = order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop=F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else return(pos)
}

G.funs = function(Gfun){
  switch(Gfun, 
         "PH" = {
           g.fun = function(x){x}
           dg.fun = function(x){rep(1, length(x))}
           ddg.fun = function(x){rep(0, length(x))}
           dddg.fun = function(x){rep(0, length(x))}
         },
         "PO" = {
           g.fun = function(x){log(1+x)}
           dg.fun = function(x){1 / (1 + x)}
           ddg.fun = function(x){- 1 / (1 + x) ^ 2}
           dddg.fun = function(x){2 / (1 + x) ^ 3}
         })
  list(g.fun = g.fun, dg.fun = dg.fun, ddg.fun = ddg.fun, 
       dddg.fun = dddg.fun)
}

MyCop = function(copula.index){
  switch(as.character(copula.index),
         "1" = {
           copula.lwr = - 1; copula.upr = 1
           tau.alpha = NULL
           Dtau.alpha = NULL
           link.default = NULL
         },
         "3" = {
           copula.lwr = 0; copula.upr = 28
           tau.alpha = function(alpha) alpha / (alpha + 2)
           Dtau.alpha = function(alpha){}
           body(Dtau.alpha) = D(expression(alpha/(alpha + 2)), "alpha")
           link.default = "log"
         },
         "23" = {
           copula.lwr = -28; copula.upr = 0
           tau.alpha = function(alpha) alpha / (alpha + 2)
           Dtau.alpha = function(alpha){}
           body(Dtau.alpha) = D(expression(alpha/(alpha + 2)), "alpha")
           link.default = "log"
         },
         "4" = {
           copula.lwr = 1; copula.upr = 17
           tau.alpha = function(alpha) 1 - 1 / alpha
           Dtau.alpha = function(alpha) 1 / (alpha ^ 2)
           link.default = "log-1"
         },
         "24" = {
           copula.lwr = -17; copula.upr = - 1
           tau.alpha = function(alpha) 1 - 1 / alpha
           Dtau.alpha = function(alpha) 1 / (alpha ^ 2)
           link.default = "log-1"
         },
         "5" = {
           copula.lwr = -35; copula.upr = 35
           tau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               1 - 4 / a + (4 / a ^ 2) *
                 integrate(fun0, lower = 0, upper = a)$value
             })
           }
           Dtau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               4 / a ^ 2 + 4 / (a * (exp(a) - 1)) -
                 (8 / a ^ 3) * integrate(fun0, lower = 0, upper = a)$value
             })
           }
           link.default = "identity"
         }
  )
  list(lwr = copula.lwr, upr = copula.upr,
       tau.alpha = tau.alpha, Dtau.alpha = Dtau.alpha, 
       link.default = link.default)
}

BiCopLink = function(copula.fam){
  switch(
    copula.fam,
    "Clayton" = {link = function(x){log(x)}},
    "fClayton" = {link = function(x){log(-x)}},
    "Gumbel" = {link = function(x){log(x - 1)}},
    "fGumbel" = {link = function(x){log(- x + 1)}},
    "Frank" = {link = function(x) x},
    "Gaussian" = {link = function(x)  (log(1 + x) - log(1 - x)) / 2}
  )
  return(link)
}

MyCopIndex = function(copula.fam){
  switch(copula.fam,
         "Clayton" = {cop.index = 3},
         "Gumbel" = {cop.index = 4},
         "Frank" = {cop.index = 5},
         "Gaussian" = {cop.index = 1},
         "fClayton" = {cop.index = 23},
         "fGumbel" = {cop.index = 24}
  )
  return(cop.index)
}

MylinkFun = function(link.fun){
  switch(link.fun,
         "identity" = {
           copula.link = list(h.fun = function(x) {x},
                              dot.h.fun = function(x) {rep(1, length(x))},
                              ddot.h.fun = function(x) {rep(0, length(x))},
                              hinv.fun = function(x) {x})
         },
         "log" = {
           copula.link = list(h.fun = function(x) {exp(x)},
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)},
                              hinv.fun = function(x){log(x)})
         },
         "neglog" = {
           copula.link = list(h.fun = function(x) {- exp(x)},
                              dot.h.fun = function(x) {- exp(x)},
                              ddot.h.fun = function(x) {- exp(x)},
                              hinv.fun = function(x){log(- x)})
         },
         "log-1" = {
           copula.link = list(h.fun = function(x) {exp(x) + 1},
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)},
                              hinv.fun = function(x) log(x - 1))
         },
         "neglog-1" = {
           copula.link = list(h.fun = function(x) {- exp(x) - 1},
                              dot.h.fun = function(x) {- exp(x)},
                              ddot.h.fun = function(x) {- exp(x)},
                              hinv.fun = function(x) log(- x - 1))
         },
         "tanh" = {
           copula.link = list(h.fun = function(x) {tanh(x)},
                              dot.h.fun = function(x) {1 - tanh(x) ^ 2},
                              ddot.h.fun = function(x) {
                                - 2 * tanh(x) * (1 - tanh(x) ^ 2)
                                },
                              hinv.fun = function(x) {
                                (log(1 + x) - log(1 - x)) / 2
                              })
         }
  )
  copula.link
}
MyCopula = function(copula.index){
  
  if (copula.index != 1){
    switch(
      as.character(copula.index),
      "3" = { # Clayton
        C.exp = expression((u1^(-para) + u2^(-para) -1)^(-1/para))
      },
      "23" = { # Clayton
        C.exp = expression(u2 - ((1 - u1) ^(para) + u2 ^(para) -1)^(1/para))
      },
      "4" = { # Gumbel
        C.exp = expression(
          exp(-((-log(u1))^(para) + (-log(u2))^(para))^(1/para))
          )
      },
      "24" = {
        C.exp = expression(
          u2 - exp(-((-log(1 - u1))^(-para) + (-log(u2))^(-para))^(-1/para))
        )
      },
      "5" = { # Frank
        C.exp = 
          expression((-1/para) * log(1 + ((exp(-para*u1)-1)*(exp(-para*u2)-1))/(exp(-para)-1)))
      }
    )
    C.fun = function(u1, u2, para){};
    body(C.fun) = C.exp
    dC.para = function(u1, u2, para){};
    body(dC.para) = D(C.exp, "para")
    dC.para2 = function(u1, u2, para){};
    body(dC.para2) = D(D(C.exp, "para"), "para")
    dC.u1u1u2u2 = function(u1, u2, para){};
    body(dC.u1u1u2u2) = D(D(D(D(C.exp, "u1"), "u1"), "u2"), "u2")
    
  } else {
    C.fun = function(u1, u2, para){
      BiCopCDF(u1 = u1, u2 = u2, family = 1, par = para)
    }
    dC.para = function(u1, u2, para){
      C.fun1 = function(u_1, u_2, alpha){
        BiCopCDF(u_1, u_2, family = 1, par = alpha)
      }
      sapply(1 : length(u1), function(k){
        numDeriv::jacobian(C.fun1, x = para[k], u_1 = u1[k], u_2 = u2[k]) %>%
          unlist()
      })
    }
    dC.para2 = function(u1, u2, para){
      C.fun1 = function(u_1, u_2, alpha){
        BiCopCDF(u_1, u_2, family = 1, par = alpha)
      }
      sapply(1 : length(u1), function(k){
        pracma::hessian(C.fun1, x0 = para[k], u_1 = u1[k], u_2 = u2[k]) %>% 
          unlist()
      }) 
    }
    dC.u1u1u2u2 = function(u1, u2, para){
      BiCopDeriv2(u1 = u1, u2 = u2, family = 1, par = para, deriv = "u1u2")
    }
  }
  list(C.fun = C.fun, dC.para = dC.para, dC.para2 = dC.para2,
       dC.u1u1u2u2 = dC.u1u1u2u2)
}


MyfCopula = function(copula.index){
  switch(
    as.character(copula.index),
    "3" = { # Clayton
      C.exp = expression((u1^(-para) + u2^(-para) -1)^(-1/para))
    },
    "4" = { # Gumbel
      C.exp = expression(exp(-((-log(u1))^(para) + (-log(u2))^(para))^(1/para)))
    },
    "5" = { # Frank
      C.exp = 
        expression((-1/para) * log(1 + ((exp(-para*u1)-1)*(exp(-para*u2)-1))/(exp(-para)-1)))
    }
  )
  C.fun = function(u1, u2, para){};
  body(C.fun) = C.exp
  dC.para = function(u1, u2, para){};
  body(dC.para) = D(C.exp, "para")
  dC.para2 = function(u1, u2, para){};
  body(dC.para2) = D(D(C.exp, "para"), "para")
  dC.u1 = function(u1, u2, para){};
  body(dC.u1) = D(C.exp, "u1")
  dC.u1u1 = function(u1, u2, para){};
  body(dC.u1u1) = D(D(C.exp, "u1"),"u1")
  dC.u1para = function(u1, u2, para){};
  body(dC.u1para) = D(D(C.exp, "u1"), "para")
  dC.u1u1para = function(u1, u2, para){};
  body(dC.u1u1para) = D(D(D(C.exp, "u1"), "u1"), "para")
  dC.u1para2 = function(u1, u2, para){};
  body(dC.u1para2) = D(D(D(C.exp, "u1"), "para"), "para")
  dC.u2 = function(u1, u2, para){};
  body(dC.u2) = D(C.exp, "u2")
  dC.u2u2 = function(u1, u2, para){};
  body(dC.u2u2) = D(D(C.exp, "u2"), "u2")
  dC.u2para = function(u1, u2, para){};
  body(dC.u2para) = D(D(C.exp, "u2"), "para")
  dC.u2u2para = function(u1, u2, para){};
  body(dC.u2u2para) = D(D(D(C.exp, "u2"), "u2"), "para")
  dC.u2para2 = function(u1, u2, para){};
  body(dC.u2para2) = D(D(D(C.exp, "u2"), "para"), "para")
  dC.u1u2 = function(u1, u2, para){};
  body(dC.u1u2)= D(D(C.exp, "u1"), "u2")
  dC.u1u2para = function(u1, u2, para){};
  body(dC.u1u2para) = D(D(D(C.exp, "u1"), "u2"), "para")
  dC.u1u1u1 = function(u1, u2, para){};
  body(dC.u1u1u1) = D(D(D(C.exp, "u1"), "u1"), "u1")
  dC.u2u2u2 = function(u1, u2, para){};
  body(dC.u2u2u2) = D(D(D(C.exp, "u2"), "u2"), "u2")
  dC.u1u1u2 = function(u1, u2, para){};
  body(dC.u1u1u2) = D(D(D(C.exp, "u1"), "u1"), "u2")
  dC.u1u2u2 = function(u1, u2, para){};
  body(dC.u1u2u2) = D(D(D(C.exp, "u1"), "u2"), "u2")
  dC.u1u1u1u2 = function(u1, u2, para){};
  body(dC.u1u1u1u2) = D(D(D(D(C.exp, "u1"), "u1"), "u1"), "u2")
  dC.u1u2u2u2 = function(u1, u2, para){};
  body(dC.u1u2u2u2) = D(D(D(D(C.exp, "u1"), "u2"), "u2"), "u2")
  dC.u1u1u2u2 = function(u1, u2, para){};
  body(dC.u1u1u2u2) = D(D(D(D(C.exp, "u1"), "u1"), "u2"), "u2")
  dC.u1u1u2para = function(u1, u2, para){};
  body(dC.u1u1u2para) = D(D(D(D(C.exp, "u1"), "u1"), "u2"), "para")
  dC.u1u2u2para = function(u1, u2, para){};
  body(dC.u1u2u2para) = D(D(D(D(C.exp, "u1"), "u2"), "u2"), "para")
  dC.u1u2para2 = function(u1, u2, para){};
  body(dC.u1u2para2) = D(D(D(D(C.exp, "u1"), "u2"), "para"), "para")
  list(C.fun = C.fun, dC.para = dC.para, dC.para2 =dC.para2,
       dC.u1 = dC.u1, dC.u1u1 = dC.u1u1, dC.u1para = dC.u1para,
       dC.u1u1para = dC.u1u1para, dC.u1para2 = dC.u1para2,
       dC.u2 = dC.u2, dC.u2u2 = dC.u2u2, dC.u2para = dC.u2para,
       dC.u2u2para = dC.u2u2para, dC.u2para2 = dC.u2para2,
       dC.u1u2 = dC.u1u2, dC.u1u2para = dC.u1u2para,
       dC.u1u1u1 = dC.u1u1u1, dC.u2u2u2 = dC.u2u2u2,
       dC.u1u1u2 = dC.u1u1u2, dC.u1u2u2 = dC.u1u2u2,
       dC.u1u1u1u2 = dC.u1u1u1u2, dC.u1u2u2u2 = dC.u1u2u2u2,
       dC.u1u1u2u2 = dC.u1u1u2u2,
       dC.u1u1u2para = dC.u1u1u2para, dC.u1u2u2para = dC.u1u2u2para,
       dC.u1u2para2 = dC.u1u2para2)
  
}

alphaTD.fun = function(gamma, newdat, fmla, copula.link, copula.index, 
                       se.fit = FALSE, Psi.gamma = NULL){
  w = model.matrix(fmla, newdat)
  lp = as.vector(w %*% gamma)
  alpha = copula.link$h.fun(lp)
  mycop = MyCop(copula.index)
  tau = mycop$tau.alpha(alpha)
  alpha.out = data.frame(para = "alpha", est = alpha) %>%
    cbind(., newdat)
  tau.out =  data.frame(para = "tau", est = tau) %>%
    cbind(., newdat)
  out = rbind(alpha.out, tau.out) %>% as.data.frame()
  
  if (se.fit){
    a = copula.link$dot.h.fun(lp) * w 
    gamma.cov = crossprod(Psi.gamma, Psi.gamma) / nrow(Psi.gamma) ^ 2
    alpha.cov = a %*% gamma.cov %*% t(a)
    alpha.se = sqrt(diag(alpha.cov))
    tau.se = alpha.se * abs(mycop$Dtau.alpha(alpha))
    out$se = c(alpha.se, tau.se)
  }
  return(out)
}
tau12.est.fun2 = function(par1, par2, par12, index1, index2, index12, 
                          pc.method = "foreach", maxlist = 100){
  switch(pc.method,
         "mclapply" = {
           fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
             out = parallel::mclapply(
               1 : length(u1), function(i){
                 fun0 = function(uu){
                   uu1 = rep(u1[i], length(uu)); uu2 = rep(u2[i], length(uu))
                   cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
                   cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
                   BiCopCDF(cuu13, cuu23, family = index12, par = par12) 
                 }
                 S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
                 fun0 = function(uu){
                   uu1 = rep(u1[i], length(uu)); uu2 = rep(u2[i], length(uu))
                   cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
                   cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
                   BiCopPDF(cuu13, cuu23, family = index12, par = par12) *
                     BiCopPDF(uu1, uu, family = index1, par = par1) *
                     BiCopPDF(uu2, uu, family = index2, par = par2) 
                 }
                 f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
                 S12 * f12
               }, mc.cores = ncores
             ) %>% unlist()
             return(out)
           }
         },
         "foreach" = {
           fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
             out = foreach(
               i = 1 : length(u1),
               .packages = c("pracma", "VineCopula"), 
               .combine = c) %dopar% {
                 fun0 = function(uu){
                   uu1 = rep(u1[i], length(uu)); uu2 = rep(u2[i], length(uu))
                   cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
                   cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
                   BiCopCDF(cuu13, cuu23, family = index12, par = par12) 
                 }
                 S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
                 fun0 = function(uu){
                   uu1 = rep(u1[i], length(uu)); uu2 = rep(u2[i], length(uu))
                   cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
                   cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
                   BiCopPDF(cuu13, cuu23, family = index12, par = par12) *
                     BiCopPDF(uu1, uu, family = index1, par = par1) *
                     BiCopPDF(uu2, uu, family = index2, par = par2) 
                 }
                 f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
                 S12 * f12
               }
             return(out)
           }
         }
         
  )
  
  
  4 * pracma::integral2(
    fun, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    par12 = par12, par1 = par1, par2 = par2,
    index12 = index12, index1 = index1, index2 = index2, 
    maxlist = maxlist)$Q - 1
}




tau12.est.fun = function(par1, par2, par12, index1, index2, index12, 
                         pc.method = "foreach", maxlist = 100){
  obj.fun1 = function(uu1, uu2, par1, par2, par12, index12, index1, index2, 
                      pc.method){
    fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopCDF(cuu13, cuu23, family = index12, par = par12) 
      }
      S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopPDF(cuu13, cuu23, family = index12, par = par12) *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopPDF(uu2, uu, family = index2, par = par2) 
      }
      f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      S12 * f12
    }
    switch(pc.method,
           "foreach" = {
             out = foreach(
               i = 1 : length(uu1),
               .packages = c("pracma", "VineCopula"), 
               .combine = c) %dopar% {
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }
           },
           "mclapply" = {
             out = parallel::mclapply(
               1 : length(uu1), function(i){
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }, mc.cores = ncores
             ) %>% unlist()
           }
    )
    return(out)
  }
  est = 4 * pracma::integral2(
    obj.fun1, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    par12 = par12, par1 = par1, par2 = par2,
    index12 = index12, index1 = index1, index2 = index2,
    pc.method = pc.method, maxlist = maxlist)$Q - 1
  return(est)
  
}

Dtau12.fun = function(par1, par2, par12, index1, index2, index12, 
                      pc.method = "foeach", maxlist = 10){
  Dpar12.fun = function(uu1, uu2, par1, par2, par12, index12, index1, index2, 
                        pc.method){
    fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        MyCopula(index12)$dC.para(cuu13, cuu23, par = par12) 
      }
      S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopDeriv(cuu13, cuu23, family = index12, par = par12, deriv = "par") *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopPDF(uu2, uu, family = index2, par = par2) 
      }
      f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      S12 * f12
    }
    switch(pc.method,
           "foreach" = {
             out = foreach(
               i = 1 : length(uu1),
               .packages = c("pracma", "VineCopula"), 
               .combine = c) %dopar% {
                 source("helpers.R")
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }
           },
           "mclapply" = {
             out = parallel::mclapply(
               1 : length(uu1), function(i){
                 fun(u1 = uu1[i], u2 = uu2[i],
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }, mc.cores = ncores)  %>% unlist()
           })
    return(out)
  }
  Dpar12.fun = function(uu1, uu2){
    fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        MyCopula(index12)$dC.para(cuu13, cuu23, par = par12) 
      }
      S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopDeriv(cuu13, cuu23, family = index12, par = par12, deriv = "par") *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopPDF(uu2, uu, family = index2, par = par2) 
      }
      f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      S12 * f12
    }
    switch(pc.method,
           "foreach" = {
             out = foreach(
               i = 1 : length(uu1),
               .packages = c("pracma", "VineCopula"), 
               .combine = c) %dopar% {
                 source("helpers.R")
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }
           },
           "mclapply" = {
             out = parallel::mclapply(
               1 : length(uu1), function(i){
                 fun(u1 = uu1[i], u2 = uu2[i],
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }, mc.cores = ncores)  %>% unlist()
           })
    return(out)
  }
  int2(Dpar12.fun, a = c(0, 0), b = c(1, 1))
  
  Dpar12.fun(uu1 = seq(0.1, 0.9, by = 0.1), uu2 = seq(0.1, 0.9, by = 0.1),
             par12 = par12, par1 = par1, par2 = par2,
             index12 = index12, index1 = index1, index2 = index2, 
             pc.method = pc.method)
  Dpar12 = 4 * pracma::integral2(
    Dpar12.fun, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    par12 = par12, par1 = par1, par2 = par2,
    index12 = index12, index1 = index1, index2 = index2, 
    pc.method = pc.method, maxlist = maxlist)$Q
  Dpar1.fun = function(uu1, uu2, par1, par2, par12, index12, index1, index2, 
                       pc.method){
    fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopHfunc1(cuu13, cuu23, family = index12, par = par12) *
          BiCopHfuncDeriv(uu1, uu, family = index1, par = par1, deriv = "par") +
          BiCopHfunc2(cuu13, cuu23, family = index12, par = par12) *
          BiCopHfuncDeriv(uu2, uu, family = index2, par = par2, deriv = "par")
      }
      S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopDeriv(cuu13, cuu23, family = index12, par = par12, deriv = "u1") *
          BiCopHfuncDeriv(uu1, uu, family = index1, par = par1, deriv = "par") *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopPDF(uu2, uu, family = index2, par = par2) +
          BiCopPDF(cuu13, cuu23, family = index12, par = par12) *
          BiCopDeriv(uu1, uu, family = index1, par = par1, deriv = "par") *
          BiCopPDF(uu2, uu, family = index2, par = par2) 
      }
      f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      S12 * f12
    }
    switch(pc.method,
           "foreach" = {
             out = foreach(i = 1 : length(uu1),
                           .packages = c("pracma", "VineCopula"), .combine = c) %dopar% {
                             fun(u1 = uu1[i], u2 = uu2[i], 
                                 par12 = par12, par1 = par1, par2 = par2,
                                 index12 = index12, index1 = index1, index2 = index2)
                           }
           },
           "mclapply" = {
             out = parallel::mclapply(
               1 : length(uu1), function(i){
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }, mc.cores = ncores) %>% unlist()
           }
    )
    return(out)
  }
  Dpar1 = 4 * pracma::integral2(
    Dpar1.fun, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    par12 = par12, par1 = par1, par2 = par2,
    index12 = index12, index1 = index1, index2 = index2,
    pc.method = pc.method, maxlist = maxlist)$Q
  
  Dpar2.fun = function(uu1, uu2, par1, par2, par12, index12, index1, index2, 
                       pc.method){
    fun = function(u1, u2, par1, par2, par12, index12, index1, index2){
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopHfunc2(cuu13, cuu23, family = index12, par = par12) *
          BiCopHfuncDeriv(uu2, uu, family = index2, par = par2, deriv = "par") +
          BiCopHfunc2(cuu13, cuu23, family = index12, par = par12) *
          BiCopHfuncDeriv(uu2, uu, family = index2, par = par2, deriv = "par")
      }
      S12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      fun0 = function(uu){
        uu1 = rep(u1, length(uu)); uu2 = rep(u2, length(uu))
        cuu13 = BiCopHfunc2(uu1, uu, family = index1, par = par1)
        cuu23 = BiCopHfunc2(uu2, uu, family = index2, par = par2)
        BiCopDeriv(cuu13, cuu23, family = index12, par = par12, deriv = "u2") *
          BiCopHfuncDeriv(uu2, uu, family = index2, par = par2, deriv = "par") *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopPDF(uu2, uu, family = index2, par = par2) +
          BiCopPDF(cuu13, cuu23, family = index12, par = par12) *
          BiCopPDF(uu1, uu, family = index1, par = par1) *
          BiCopDeriv(uu2, uu, family = index2, par = par2, deriv = "par") 
      }
      f12 = pracma::integral(fun0, xmin = 0, xmax = 1)
      S12 * f12
    }
    switch(pc.method,
           "foreach" = {
             out = foreach(
               i = 1 : length(uu1),
               .packages = c("pracma", "VineCopula"), 
               .combine = c) %dopar% {
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }
           },
           "mclapply" = {
             out = parallel::mclapply(
               1 : length(uu1), function(i){
                 fun(u1 = uu1[i], u2 = uu2[i], 
                     par12 = par12, par1 = par1, par2 = par2,
                     index12 = index12, index1 = index1, index2 = index2)
               }, mc.cores = ncores) %>% unlist()
           }
    )
    return(out)
  }
  Dpar2 = 4 * pracma::integral2(
    Dpar2.fun, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    par12 = par12, par1 = par1, par2 = par2,
    index12 = index12, index1 = index1, index2 = index2,
    pc.method = pc.method, maxlist = maxlist)$Q
  
  c(Dpar12, Dpar1, Dpar2)
}

ST.uD = function(uuD, zT, zD, betaT, betaD, LambdaT, LambdaD, ttT, ttD){
  ebz = exp(sum(zD * betaD))
  index = sum.I(- log(uuD) / ebz, ">=", LambdaD)
  tt = rep(0, length(uuD))
  tt[which(index !=0)] = ttD[index[index!=0]]
  exp(- exp(sum(zT * betaT)) * c(0, LambdaT)[sum.I(tt, ">=", ttT) + 1])
}

