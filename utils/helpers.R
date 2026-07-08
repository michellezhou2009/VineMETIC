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
           "none" = {
             out = sapply(1 : length(uu1), function(i){
               fun(u1 = uu1[i], u2 = uu2[i],
                   par12 = par12, par1 = par1, par2 = par2,
                   index12 = index12, index1 = index1, index2 = index2)
             })
           },
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
