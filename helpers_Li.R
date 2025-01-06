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
g.fn = function(w1, w2, theta) {
  (w1^(-theta) - w2^(-theta) + 1)^(-1 / theta)
}