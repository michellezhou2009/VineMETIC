#' Nonparametric MLE for semi-parametric transformation regression model (SPT)
#'
#' @description
#'
#' @param data a data frame containing variables names in the \code{formula}.
#' @param time a character string specifying the variable name in the \code{data} for the event time.
#' @param status a character string specifying the variable name in the \code{data} for the censoring indicator: \code{1} indicates the time is observed, and \code{0} indicates censored.
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic description of the model to be fitted.
#' @param Gfun a character string specifying the link function, \code{Gfun = "PH"} (default) or \code{Gfun = "PO"} in the SPT model.
#'
#' @import dplyr trust rlang
#'
#' @export
#'
#' @return a list of the following components:
#' \describe{
#'  \item{\code{beta}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for regresson coefficients.}
#'  \item{\code{dLambda}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for the jump size of the baseline function at observed time points.}
#'  \item{\code{Lambda}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for the baseline function at observed time points.}
#'  \item{\code{varcov}}{a list containing \code{model}, the model-based variance-covariance matrix, and \code{robust}, the robust variance-covariance matrix for all the parameters.}
#'  \item{\code{Psi.theta}}{a data frame containing the estimated asymptotic expansion of the NPMLE estimator.}
#'  \item{\code{call}}{a list containing the specified values of the arguments.}
#'  \item{\code{convergence}}{a logical value indicating whether the maximization of the log-likelihood converges or not.}
#'  \item{\code{niter}}{an interger which is the number of iterations for the maximizaiton of the log-likelihood.}
#' }
#'
#' @examples
#' \code{data(BMT, package = "SemiCompRisks")}
#' \code{data = BMT %>%
#'   mutate(g = factor(g, levels = c(2, 3, 1),
#'                     labels = c("AML-low", "AML-high", "ALL")))}
#' \code{fitSPT(data, time = "T1", status = "delta1",
#'               formula = ~ g, Gfun = "PH")$beta}
fitSPT = function(data, time = "time", status = "status", formula = ~ 1,
                  Gfun = "PH", se.fit = TRUE){
  
  if (!rlang::is_formula(formula))
    stop("Argument \"formula\" needs to be a formula.")
  allvars = list(time, status); allnm = c("time", "status")
  tmp.index = sapply(allvars, function(x) is.character(x))
  if (any(tmp.index == FALSE))
    stop(paste0("Values of arguments: \"",
                paste0(allnm[which(tmp.index == FALSE)], collapse = ", "),
                "\" needs to be characters."))
  var.nm = all.vars(formula)
  allvars = c(time, status, var.nm)
  exists.index = allvars %in% colnames(data)
  if (any(exists.index == FALSE))
    stop(paste0("Variables: ",
                paste0(allvars[exists.index == FALSE], collapse = ", ")),
         " do not exists in the data.")
  
  if (Gfun == "PH"){
    out = fitPH(data = data, time = time, status = status, formula = formula, 
                se.fit = se.fit)
  } else {
    N = nrow(data)
    xi = data[, colnames(data) == time]
    di = data[, colnames(data) == status]
    zi = model.matrix(formula, data)[, - 1, drop = F]
    n.b = ncol(zi)
    tt = sort(unique(xi[di == 1])); n.tt = length(tt)
    G.all = G.funs(Gfun)
    
    
    lln.fun <- function(theta, xi, di, zi, G.all){
      g.fun = G.all$g.fun; dg.fun = G.all$dg.fun;
      n.b = ncol(zi);
      tt = sort(unique(xi[di == 1])); n.tt = length(tt)
      if (n.b == 0) {
        dLambda = theta
      } else {
        b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
      }
      if (any(dLambda <= 0)) return(-Inf) else {
        if (n.b == 0) bzi = rep(0, length(xi)) else bzi = as.vector(zi %*% b)
        log.dLambda = rep(0, length(xi));
        log.dLambda[di == 1] = log(dLambda[match(xi[di == 1], tt)])
        Lambda.xi = sum.I(xi, ">=", tt, dLambda)
        
        sum(di * bzi + log.dLambda + di * log(dg.fun(Lambda.xi * exp(bzi))) -
              g.fun(Lambda.xi * exp(bzi)))
      }
    }
    dll.fun <- function(theta, xi, di, zi, G.all){
      g.fun = G.all$g.fun; dg.fun = G.all$dg.fun; ddg.fun = G.all$ddg.fun
      n.b = ncol(zi);
      tt = sort(unique(xi[di == 1])); n.tt = length(tt)
      if (n.b == 0) {
        dLambda = theta
      } else {
        b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
      }
      if (any(dLambda <=0)) return(NA) else{
        if (n.b == 0) bzi = rep(0, length(xi)) else bzi = as.vector(zi %*% b)
        Lambda.xi = sum.I(xi, ">=", tt, dLambda)
        xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
        dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) * di
        out = (
          di * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
            dg.fun(Lambda.xi * exp(bzi)) -
            dg.fun(Lambda.xi * exp(bzi)) * exp(bzi)
        ) * xi.g.tt +
          dN.tt / matrix(dLambda, byrow = T, nrow = length(xi), ncol = n.tt)
        if (n.b > 0) {
          dll.b = zi * (
            di + di * Lambda.xi * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
              dg.fun(Lambda.xi * exp(bzi)) - dg.fun(Lambda.xi * exp(bzi)) *
              Lambda.xi * exp(bzi)
          )
          out = cbind(dll.b, out)
        }
        return(out)
      }
    }
    ddlln.fun <- function(theta, xi, di, zi, G.all){
      g.fun = G.all$g.fun; dg.fun = G.all$dg.fun; ddg.fun = G.all$ddg.fun
      dddg.fun = G.all$dddg.fun
      n.b = ncol(zi);
      tt = sort(unique(xi[di == 1])); n.tt = length(tt)
      if (n.b == 0) {
        dLambda = theta
      } else {
        b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
      }
      if (any(dLambda <=0))
        return(matrix(-Inf, length(theta), length(theta))) else{
          if (n.b == 0) bzi = rep(0, length(xi)) else bzi = as.vector(zi %*% b)
          Lambda.xi = sum.I(xi, ">=", tt, dLambda)
          xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
          dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) *
            di
          out = crossprod(
            (di * exp(bzi) ^ 2 *
               (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) -
                  (ddg.fun(Lambda.xi * exp(bzi)) /
                     dg.fun(Lambda.xi * exp(bzi))) ^ 2) -
               ddg.fun(Lambda.xi * exp(bzi)) * exp(bzi) ^ 2
            ) * xi.g.tt, xi.g.tt)  - diag(apply(dN.tt, 2, sum) / dLambda ^ 2)
          
          if (n.b > 0){
            dll.bb = crossprod(
              (di * Lambda.xi * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
                 dg.fun(Lambda.xi * exp(bzi)) +
                 di * Lambda.xi ^ 2 * exp(bzi) ^ 2 *
                 (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) -
                    (ddg.fun(Lambda.xi * exp(bzi)) /
                       dg.fun(Lambda.xi * exp(bzi))) ^ 2) -
                 ddg.fun(Lambda.xi * exp(bzi)) * Lambda.xi ^ 2 * exp(bzi) ^ 2 -
                 dg.fun(Lambda.xi * exp(bzi)) * Lambda.xi * exp(bzi)) *
                zi, zi)
            dll.dLambdab = crossprod(
              (di * exp(bzi) *
                 ddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) +
                 di * exp(bzi) ^ 2 * Lambda.xi *
                 (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) -
                    (ddg.fun(Lambda.xi * exp(bzi)) /
                       dg.fun(Lambda.xi * exp(bzi))) ^ 2) -
                 exp(bzi) * dg.fun(Lambda.xi * exp(bzi)) -
                 exp(bzi) ^ 2 * Lambda.xi * ddg.fun(Lambda.xi * exp(bzi))) * xi.g.tt
              , zi)
            out = rbind(
              cbind(dll.bb, t(dll.dLambdab)),
              cbind(dll.dLambdab, out)
            )
          }
          return(out)
        }
    }
    if (n.b > 0) b.ini = rep(0, n.b) else b.ini = NULL
    dLambda.ini = rep(0.1, n.tt)
    G.all = G.funs(Gfun)
    theta.ini  = c(b.ini, dLambda.ini)
    objfun<- function(x){
      f = lln.fun(theta = x, xi, di, zi, G.all)
      g = dll.fun(theta = x,xi, di, zi, G.all)
      if (any(!is.na(g))) g = apply(g, 2, sum)
      B = ddlln.fun(theta = x, xi, di, zi, G.all)
      list(value = f, gradient = g, hessian = B)
    }
    usedtime <- system.time({
      res <- trust::trust(objfun, theta.ini, 5, 100, iterlim = 300,
                          minimize= FALSE, blather = T)})
    if (!res$converged) stop("MLE did not converge.")
    
    theta.est = res$argument
    if (n.b > 0){
      b.est = theta.est[c(1 : n.b)]; dLambda.est = theta.est[ - c(1 : n.b)]
    } else {
      b.est = NULL; dLambda.est = theta.est
    }
    if (se.fit){
      Imat = - ddlln.fun(theta.est, xi, di, zi, G.all) / N
      dll = dll.fun(theta.est, xi, di, zi, G.all)
      Psi.theta = dll %*% solve(Imat)
      theta.cov.model = solve(Imat) / N
      theta.cov.robust = crossprod(Psi.theta, Psi.theta) / (N ^ 2)
      if (n.b > 0){
        beta.summary = data.frame(
          est = b.est, 
          se = sqrt(diag(theta.cov.robust)[c(1 : n.b)]),
          model.se = sqrt(diag(theta.cov.model)[c(1 : n.b)])
        ); rownames(beta.summary) = colnames(zi)
        Psi.beta = Psi.theta[ , 1 : n.b, drop = F]
        dLambda.cov.model = theta.cov.model[- c(1 : n.b), - c(1 : n.b)]
        dLambda.cov.robust = theta.cov.robust[- c(1 : n.b), - c(1 : n.b)]
        Psi.dLambda = Psi.theta[, - c(1 : n.b), drop = F]
      } else {
        beta.summary = NULL
        Psi.beta = NULL
        dLambda.cov.model = theta.cov.model
        dLambda.cov.robust = theta.cov.robust
        Psi.dLambda = Psi.theta
      }
      dLambda.summary = data.frame(
        time = tt, est = dLambda.est,
        se = sqrt(diag(dLambda.cov.robust)),
        model.se = sqrt(diag(dLambda.cov.model))
      );
      Psi.Lambda = t(apply(Psi.dLambda, 1, cumsum))
      tmp = matrix(1, nrow = n.tt, ncol = n.tt)
      tmp[upper.tri(tmp, diag = FALSE)] = 0
      Lambda.summary = data.frame(
        time = tt, est = cumsum(dLambda.est),
        se = sqrt(apply(Psi.Lambda ^ 2, 2, sum)) / N,
        model.se = sqrt(diag(tmp %*% dLambda.cov.model %*% t(tmp)))
      )
      out = list(
        call = list(time = time, status = status, G.all = G.all, formula = formula),
        beta = beta.summary, dLambda = dLambda.summary, Lambda = Lambda.summary,
        varcov = list(model = theta.cov.model, robust = theta.cov.robust),
        Psi.theta = list(beta = Psi.beta, dLambda = Psi.dLambda,
                         Lambda = Psi.Lambda))
    } else {
      if (n.b > 0){ 
        beta.summary = data.frame(est = b.est)
        rownames(beta.summary) = colnames(zi)
      } else beta.summary = NULL
      dLambda.summary = data.frame(time = tt, est = dLambda.est)
      Lambda.summary = data.frame(time = tt, est = cumsum(dLambda.est))
      out = list(
        call = list(time = time, status = status, G.all = G.all, formula = formula),
        beta = beta.summary, dLambda = dLambda.summary, Lambda = Lambda.summary)
    }
    
  }
  return(out)
}

fitPH = function(data, time = "time", status = "status", formula = ~ 1, se.fit = TRUE){
  
  N = nrow(data)
  xi = data[, colnames(data) == time]
  di = data[, colnames(data) == status]
  zi = model.matrix(formula, data)[, - 1, drop = F]
  n.b = ncol(zi)
  
  ph.fmla = update(formula, as.formula(paste0("Surv(", time, ",", status, ") ~ .")))
  
  if (n.b == 0) {
    b.est = NULL; bzi = rep(0, N)
    myfit = summary(survfit(ph.fmla, data, type = "fh"))
    Lambda.fit = list(
      time = myfit$time, hazard = myfit$cumhaz
    )
  } else {
    myfit = coxph(ph.fmla, data, method = "breslow")
    b.est = myfit$coef; 
    bzi = as.vector(zi %*% matrix(b.est, byrow = F, ncol = 1))
    Lambda.fit = basehaz(myfit, centered = FALSE)
  }
  Lambda = Lambda.fit$hazard
  dLambda = c(Lambda[1], diff(Lambda, lag = 1))
  tt = Lambda.fit$time[dLambda != 0]; n.tt = length(tt)
  dLambda = dLambda[dLambda != 0]
  
  if (se.fit){
    theta.est = c(b.est, dLambda)
    ## dll
    Lambda.xi = sum.I(xi, ">=", tt, dLambda)
    xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
    dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) * di
    out =  -  exp(bzi) * xi.g.tt +
      dN.tt / matrix(dLambda, byrow = T, nrow = length(xi), ncol = n.tt)
    if (n.b > 0) {
      dll.b = zi * (di - Lambda.xi * exp(bzi))
      out = cbind(dll.b, out)
    }
    dll = out
    
    out = - diag(apply(dN.tt, 2, sum) / dLambda ^ 2)
    if (n.b > 0){
      dll.bb = - crossprod(
        Lambda.xi * exp(bzi) * zi, zi)
      dll.dLambdab = - crossprod(exp(bzi) * xi.g.tt , zi)
      out = rbind(
        cbind(dll.bb, t(dll.dLambdab)),
        cbind(dll.dLambdab, out)
      )
    }
    ddlln = out
    Imat = - ddlln / N
    Psi.theta = dll %*% solve(Imat)
    theta.cov.model = solve(Imat) / N
    theta.cov.robust = crossprod(Psi.theta, Psi.theta) / (N ^ 2)
    if (n.b > 0){
      beta.summary = data.frame(
        est = b.est, 
        se = sqrt(diag(theta.cov.robust)[c(1 : n.b)]),
        model.se = sqrt(diag(theta.cov.model)[c(1 : n.b)])
      ); rownames(beta.summary) = colnames(zi)
      Psi.beta = Psi.theta[ , 1 : n.b, drop = F]
      dLambda.cov.model = theta.cov.model[- c(1 : n.b), - c(1 : n.b)]
      dLambda.cov.robust = theta.cov.robust[- c(1 : n.b), - c(1 : n.b)]
      Psi.dLambda = Psi.theta[, - c(1 : n.b), drop = F]
    } else {
      beta.summary = NULL
      Psi.beta = NULL
      dLambda.cov.model = theta.cov.model
      dLambda.cov.robust = theta.cov.robust
      Psi.dLambda = Psi.theta
    }
    dLambda.summary = data.frame(
      time = tt, est = dLambda,
      se = sqrt(diag(dLambda.cov.robust)),
      model.se = sqrt(diag(dLambda.cov.model))
    );
    Psi.Lambda = t(apply(Psi.dLambda, 1, cumsum))
    tmp = matrix(1, nrow = n.tt, ncol = n.tt)
    tmp[upper.tri(tmp, diag = FALSE)] = 0
    Lambda.summary = data.frame(
      time = tt, est = cumsum(dLambda),
      se = sqrt(apply(Psi.Lambda ^ 2, 2, sum)) / N,
      model.se = sqrt(diag(tmp %*% dLambda.cov.model %*% t(tmp)))
    )
    G.all = G.funs("PH")
    out = list(
      call = list(time = time, status = status, G.all = G.all, formula = formula),
      beta = beta.summary, dLambda = dLambda.summary, Lambda = Lambda.summary,
      varcov = list(model = theta.cov.model, robust = theta.cov.robust),
      Psi.theta = list(beta = Psi.beta, dLambda = Psi.dLambda,
                       Lambda = Psi.Lambda))
  } else {
    if (n.b > 0){ 
      beta.summary = data.frame(est = b.est)
      rownames(beta.summary) = colnames(zi)
    } else beta.summary = NULL
    dLambda.summary = data.frame(time = tt, est = dLambda)
    Lambda.summary = data.frame(time = tt, est = cumsum(dLambda))
    G.all = G.funs("PH")
    out = list(
      call = list(time = time, status = status, G.all = G.all, formula = formula),
      beta = beta.summary, dLambda = dLambda.summary, Lambda = Lambda.summary)
  }
  

  return(out)
}
#' Survival estimation under semi-parametric transformation regression model (SPT) for a given data
#'
#' @description
#'
#' @param obj an object of \code{\link[PMLE4SCR]{fitSPT}}.
#' @param newdata a data frame containing the same variables in the \code{formula} of \code{\link[PMLE4SCR]{fitSPT}}.
#'
#' @import dplyr
#'
#' @export
#'
#' @return a list of the following components:
#' \describe{
#'  \item{\code{surv}}{a data frame containing the estimate, model-based standard eror (SE), and robust SE for the survival probabilities for a given data.}
#'  \item{\code{Psi.surv}}{a data frame containing the estimated asymptotic expansion of the estimator of the survival probabilities.}
#' }
#'
#'
predict.fitSPT= function(obj, newdata, se.fit = TRUE){
  time = obj$call$time; status = obj$call$status
  formula = obj$call$formula
  g.fun = obj$call$G.all$g.fun; dg.fun = obj$call$G.all$dg.fun
  tt = as.vector(obj$dLambda$time)
  dLambda = as.vector(obj$dLambda$est)
  xk = newdata[ , colnames(newdata) == time]
  dk = newdata[ , colnames(newdata) == status]
  
  if (is.null(obj$beta)){
    bzk = rep(0, nrow(newdata))
  } else {
    b = as.vector(obj$beta$est)
    zk = model.matrix(formula, newdata)[ , -1, drop = F]
    bzk = as.vector(zk %*% b)
  }
  Lambda.xk = sum.I(xk, ">=", tt, dLambda)
  uk = exp( - g.fun(Lambda.xk * exp(bzk)))
  if (se.fit){
    xk.g.tt = 1 * (xk >= matrix(tt, byrow = T,
                                nrow = length(xk), ncol = length(tt)))
    Psi.theta = cbind(obj$Psi.theta$beta, obj$Psi.theta$dLambda)
    A = xk.g.tt * exp(bzk)
    if (!is.null(obj$beta)) A = cbind(zk * Lambda.xk * exp(bzk), A) 
    A = A * uk * dg.fun(Lambda.xk * exp(bzk))
    se.uk = sqrt(diag(A %*% obj$varcov$model %*% t(A)))
    Psi.uk = - Psi.theta %*% t(A)
    robust.se.uk = sqrt(apply(Psi.uk ^ 2, 2, sum)) / nrow(Psi.uk)
    
    list(surv = data.frame(fit = uk, model.se = se.uk, se = robust.se.uk),
         Psi.surv = Psi.uk)
  } else {
    list(surv = data.frame(fit = uk))
  }
}

predict.fitPH = function(obj, newdata, se.fit = TRUE){
  time = obj$call$time; status = obj$call$status
  formula = obj$call$formula
  
  tt = as.vector(obj$dLambda$time)
  dLambda = as.vector(obj$dLambda$est)
  xk = newdata[ , colnames(newdata) == time]
  dk = newdata[ , colnames(newdata) == status]
  
  if (is.null(obj$beta)){
    bzk = rep(0, nrow(newdata))
  } else {
    b = as.vector(obj$beta$est)
    zk = model.matrix(formula, newdata)[ , -1, drop = F]
    bzk = as.vector(zk %*% b)
  }
  Lambda.xk = sum.I(xk, ">=", tt, dLambda)
  uk = exp( - Lambda.xk * exp(bzk))
  if (se.fit){
    xk.g.tt = 1 * (xk >= matrix(tt, byrow = T,
                                nrow = length(xk), ncol = length(tt)))
    Psi.theta = cbind(obj$Psi.theta$beta, obj$Psi.theta$dLambda)
    A = xk.g.tt * exp(bzk)
    if (!is.null(obj$beta)) A = cbind(zk * Lambda.xk * exp(bzk), A) 
    A = A * uk
    se.uk = sqrt(diag(A %*% obj$varcov$model %*% t(A)))
    Psi.uk = - Psi.theta %*% t(A)
    robust.se.uk = sqrt(apply(Psi.uk ^ 2, 2, sum)) / nrow(Psi.uk)
    
    list(surv = data.frame(fit = uk, model.se = se.uk, se = robust.se.uk),
         Psi.surv = Psi.uk)
  } else {
    list(surv = data.frame(fit = uk))
  }
}
  
