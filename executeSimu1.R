ncores = 2
cl = makeCluster(ncores)
registerDoSNOW(cl)
for (r in 1 : 2){
  cat("Rep", r, ":")
  file.nm = paste0(dir.nm, "n", nsamp, "_rep", r, ".RData")
  
  if (!file.exists(file.nm) & !(r %in% bad.ls)){
    set.seed(seeds[r])
    data = data.gen.fun(nsamp = nsamp, copula.fams = copula.fams, 
                        copula.links = copula.links, 
                        copula.pars = copula.pars, 
                        gammas = gammas, seed0 = seeds[r])
    
    cen.rate = data.frame(
      "T1" = 1 - mean(data$delta1),
      "T2" = 1 - mean(data$delta2),
      "TD" = 1 - mean(data$deltaD)
    ) %>% mutate(rep = r, nsamp = nsamp)
    
    user.time = system.time({
      out = try(
        cvine.trivic(
        data, times, deltas, fmlas, Gfuns, copula.fams, 
        copula.links, copula.fmlas, ncores = ncores, pc.method = "foreach"), 
        silent = TRUE
      )
    })
    
    if (class(out) != "try-error"){
      para.est = rbind(
        out$T1$beta %>% mutate(dist = "T1", true = betas$t1, para = rownames(.)),
        out$T2$beta %>% mutate(dist = "T2", true = betas$t2, para = rownames(.)),
        out$TD$beta %>% mutate(dist = "TD", true = betas$tD, para = rownames(.)),
        out$`Cop(1,3)` %>% mutate(dist = "C(1,D)", true = copula.pars$`(1,3)`,
                                  para = rownames(.)),
        out$`Cop(2,3)` %>% mutate(dist = "C(2,D)", true = copula.pars$`(2,3)`,
                                  para = rownames(.)),
        out$`Cop(1,2)` %>% mutate(dist = "C(1,2)", true = copula.pars$`(1,2)`,
                                  para = rownames(.))
      ) %>% mutate(rep = r, nsamp = nsamp)
      
      pp = tt.all$pp
      tt = tt.all$t1; Lambda.est = out$T1$cumhaz; 
      survT1 = Surv.baseline(Lambda.est = Lambda.est, tt = tt) %>%
        mutate(true = pp) %>% mutate(rep = r, nsamp = nsamp)
      tt = tt.all$t2; Lambda.est = out$T2$cumhaz; 
      survT2 = Surv.baseline(Lambda.est = Lambda.est, tt = tt) %>%
        mutate(true = pp) %>% mutate(rep = r, nsamp = nsamp)
      tt = tt.all$tD; Lambda.est = out$TD$cumhaz; 
      survTD = Surv.baseline(Lambda.est = Lambda.est, tt = tt) %>%
        mutate(true = pp) %>% mutate(rep = r, nsamp = nsamp)
      
      Psi.cop = out$Psi[c("Cop(1,3)", "Cop(2,3)", "Cop(1,2)")]
      
      save(cen.rate, para.est, Psi.cop, survT1, survT2, survTD, file = file.nm)
    } else {
      print("Error")
    }
    
    print(user.time)
    
    
  }
}
stopCluster(cl)
