

density_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopCDF(cu1, cu2, family = index12, par = par12)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar12_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    MyCopula(index12)$dC.para(cu1, cu2, par12)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) 
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) + 
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Du1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) * 
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar12.Dpar12_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    MyCopula(index12)$dC.para2(cu1, cu2, par12)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv2(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar1.Dpar1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) ^ 2 +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv2(u1, u3, deriv = "par", family = index1, par = par1) 
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopHfuncDeriv2(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) ^ 2  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv2(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) +
      # term 2
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "par", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv2(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    
    # term1
    BiCopDeriv2(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) ^ 2 *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv2(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "par",family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      # term 2
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar2.Dpar2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) ^ 2 +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv2(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) ^ 2* 
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv2(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1 
    BiCopHfuncDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) ^ 2 * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv2(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopDeriv(u2, u3, deriv = "par",family = index2, par = par2) +
      # term2
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2",family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) ^ 2 *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) + 
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv2(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv2(u2, u3, deriv = "par", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Du1.Du1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopHfuncDeriv2(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 3 +
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      2 * BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) +
      # term2
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) * 
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term 1
    BiCopDeriv2(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 3 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      2 * BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Du2.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1 
    BiCopHfuncDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 3 +
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      2 * BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      # term2
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u2, u3, deriv = "u1", family = index2, par = par2)
    
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 3 +
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      2 * BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv2(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}


density.Dpar12.Dpar1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) 
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu2, cu1, deriv = "par1u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopHfuncDeriv(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv2(cu1, cu2, deriv = "par1u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) + 
      BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar12.Dpar2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu1, cu2, deriv = "par1u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopHfuncDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv2(cu1, cu2, deriv = "par1u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar12.Du1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu2, cu1, deriv = "par1u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 +
      BiCopHfuncDeriv(cu2, cu1, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) * 
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv2(cu1, cu2, deriv = "par1u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar12.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv2(cu1, cu2, deriv = "par1u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopHfuncDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv2(cu1, cu2, deriv = "par1u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopDeriv(cu1, cu2, deriv = "par", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar1.Dpar2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2)  *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) 
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    MyCopula(index12)$dC.u1u1u2u2(cu1, cu2, par12) * 
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar1.Du1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3,  family = index1, par = par1) * 
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) 
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term 1
    BiCopHfuncDeriv2(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 +
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) +
      # term2
      BiCopHfuncDeriv(cu2, cu1, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) +
      BiCopHfunc1(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "par1u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv2(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1)  *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u1, u3, deriv = "par1u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar1.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1)  *
      BiCopPDF(u1, u3, family = index1, par = par1) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1 
    MyCopula(index12)$dC.u1u1u2u2(cu1, cu2, par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "par", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar2.Du1_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    MyCopula(index12)$dC.u1u1u2u2(cu1, cu2, par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (1, 1)
  return(out)
}

density.Dpar2.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)  *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2)  *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1) + 
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopHfuncDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2)  ^ 2 + 
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) * 
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      #term2
      BiCopHfuncDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopHfunc2(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv2(u2, u3, deriv = "par1u1", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    # term1
    BiCopDeriv2(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) +
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopHfuncDeriv(u2, u3, deriv = "par", family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      # term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "par", family = index2, par = par2) +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) *
      BiCopDeriv2(u2, u3, deriv = "par1u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}


density.Du1.Du2_trivic = function(index12){
  out = c();
  out[[1]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1)
  } # (0, 0)
  out[[2]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2  +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u2, u3, family = index2, par = par2) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1)
  } # (1, 0)
  out[[3]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) * 
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 +
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) * 
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (0, 1)
  out[[4]] = function(u1, u2, u3, par12, par1, par2, index12, index1, index2){
    cu1 = BiCopHfunc2(u1, u3, family = index1, par = par1)
    cu2 = BiCopHfunc2(u2, u3, family = index2, par = par2)
    #term1
    MyCopula(index12)$dC.u1u1u2u2(cu1, cu2, par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 + 
      BiCopDeriv(cu1, cu2, deriv = "u1", family = index12, par = par12) *
      BiCopPDF(u1, u3, family = index1, par = par1) ^ 2 *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2) +
      #term2
      BiCopDeriv(cu1, cu2, deriv = "u2", family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopPDF(u2, u3, family = index2, par = par2) ^ 2 + 
      BiCopPDF(cu1, cu2, family = index12, par = par12) *
      BiCopDeriv(u1, u3, deriv = "u1", family = index1, par = par1) *
      BiCopDeriv(u2, u3, deriv = "u1", family = index2, par = par2)
  } # (1, 1)
  return(out)
}


# funsData_trivic <- function(Funs, d1, d2, u1, u2, u3, alpha12, alpha1, alpha2,
#                             index12, index1, index2, integrate = TRUE, 
#                             pc.method = "foeach"){
#   
#   funs = Funs(index12)
#   grp = d1 + 2 * d2 + 1
#   out = as.list(1 : 4)
#   for (k in 1 : 4){
#     myfun = funs[[k]]
#     index = which(grp == k)
#     if (sum(index) > 0) {
#       if (integrate) {
#         if (pc.method == "foreach"){
#           ll = foreach(i = index, .packages = c("VineCopula", "dplyr")) %dopar% {
#             source("helpers.R")
#             myint = function(uu){
#               myfun(
#                 u1 = rep(u1[i], length(uu)), u2 = rep(u2[i], length(uu)),
#                 u3 = uu, par12 = rep(alpha12[i], length(uu)), 
#                 par1 = rep(alpha1[i], length(uu)), 
#                 par2 = rep(alpha2[i], length(uu)),
#                 index12 = index12, index1 = index1, index2 = index2)
#             }
#             pracma::integral(myint, xmin = 0, xmax = u3[i])
#           }
#           ll = unlist(ll)
#         } else {
#           ll =  parallel::mclapply(
#             index, function(i){
#               myint = function(uu){
#                 myfun(
#                   u1 = rep(u1[i], length(uu)), u2 = rep(u2[i], length(uu)),
#                   u3 = uu, par12 = rep(alpha12[i], length(uu)), 
#                   par1 = rep(alpha1[i], length(uu)), 
#                   par2 = rep(alpha2[i], length(uu)),
#                   index12 = index12, index1 = index1, index2 = index2)
#               }
#               pracma::integral(myint, xmin = 0, xmax = u3[i])
#             }
#           ) 
#           ll = unlist(ll)
#         }
#       } else {
#         ll = myfun(
#           u1 = u1[index], u2 = u2[index], u3 = u3[index],
#           par12 = alpha12[index], par1 = alpha1[index], par2 = alpha2[index],
#           index12 = index12, index1 = index1, index2 = index2)
#       }
#       out[[k]] = data.frame(index = index, ll = ll)
#     }
#   }
#   out = do.call(rbind, out) %>% arrange(index, )
#   out$ll
# }

funsData_trivic <- function(Funs, d1, d2, u1, u2, u3, alpha12, alpha1, alpha2,
                            index12, index1, index2, integrate = TRUE, 
                            pc.method = "foeach"){
  
  funs = Funs(index12)
  grp = d1 + 2 * d2 + 1
  out = as.list(1 : 4)
  for (k in 1 : 4){
    myfun = funs[[k]]
    index = which(grp == k)
    if (sum(index) > 0) {
      if (integrate) {
        if (pc.method == "foreach"){
          ll = foreach(i = index, .packages = c("VineCopula", "dplyr")) %dopar% {
            source("helpers.R")
            myint = function(uu){
              myfun(
                u1 = rep(u1[i], length(uu)), u2 = rep(u2[i], length(uu)),
                u3 = uu, par12 = rep(alpha12[i], length(uu)), 
                par1 = rep(alpha1[i], length(uu)), 
                par2 = rep(alpha2[i], length(uu)),
                index12 = index12, index1 = index1, index2 = index2)
            }
            pracma::integral(myint, xmin = 0, xmax = u3[i])
            # integrate(myint, lower = 0, upper = u3[i])$value
            }
          ll = unlist(ll)
        } else {
          ll =  parallel::mclapply(
            index, function(i){
              myint = function(uu){
                myfun(
                  u1 = rep(u1[i], length(uu)), u2 = rep(u2[i], length(uu)),
                  u3 = uu, par12 = rep(alpha12[i], length(uu)), 
                  par1 = rep(alpha1[i], length(uu)), 
                  par2 = rep(alpha2[i], length(uu)),
                  index12 = index12, index1 = index1, index2 = index2)
              }
              pracma::integral(myint, xmin = 0, xmax = u3[i])
              # integrate(myint, lower = 0, upper = u3[i])$value
            }
          ) 
          ll = unlist(ll)
        }
      } else {
        ll = myfun(
          u1 = u1[index], u2 = u2[index], u3 = u3[index],
          par12 = alpha12[index], par1 = alpha1[index], par2 = alpha2[index],
          index12 = index12, index1 = index1, index2 = index2)
      }
      out[[k]] = data.frame(index = index, ll = ll)
    }
  }
  out = do.call(rbind, out) %>% arrange(index, )
  out$ll
}


prepare_trivic = function(theta, index12, link12, Wmat12, control12, N, 
                          thetaT = FALSE, datT = NULL){
  
  theta0 = theta
  n.gamma = ncol(Wmat12); copula.para = theta0[1 : n.gamma]
  alphai = link12$h.fun(
    as.vector(Wmat12 %*% matrix(copula.para, byrow = F, ncol = 1))
  )
  yes.constraint = min(alphai) < control12$lwr | 
    max(alphai) > control12$upr
  
  if (!thetaT) {
    theta0 = theta[- c(1 : n.gamma12)]
    
    ## T1 ----
    Xi = datT$X1; deltaT = datT$d1; ZmatT = datT$Zmat1; 
    GfunT = datT$Gfun1; Wmat = datT$Wmat1; copula.link = datT$link1
    control = datT$control1
    
    n.bT = ncol(ZmatT); n.gamma = ncol(Wmat); 
    tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
    Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk)) 
    dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
    copula.para = theta0[1 : n.gamma]
    alphai = copula.link$h.fun(
      as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
    )
    yes.constraint =  yes.constraint |  
      min(alphai) < control$lwr | max(alphai) > control$upr
    if (n.bT == 0){
      dLambdaT = theta0[n.gamma + c(1 : n.tk)]; bzT = rep(0, N) 
    } else{
      bT = theta0[n.gamma + c(1 : n.bT)]; 
      dLambdaT = theta0[n.gamma + n.bT + c(1 : n.tk)] 
      bzT = as.vector(ZmatT %*% bT); 
    }
    dLambdaT.Xi = rep(0, N)
    dLambdaT.Xi[deltaT == 1] = dLambdaT[match(Xi[deltaT == 1], tk)]
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
    G.all = G.funs(GfunT)
    uiT = exp(- G.all$g.fun(LambdaT.Xi * exp(bzT)))
    yes.constraint = yes.constraint | any(dLambdaT < 0) 
    dat1 = list(
      alphai = alphai, ui = uiT, G.all = G.all, Lambda.Xi = LambdaT.Xi, 
      dLambda.Xi = dLambdaT.Xi, dLambda.tk = dLambdaT,
      bz = bzT, Xi.g.tk = Xi.g.tk, dN.tk = dN.tk, Zmat = ZmatT)
    
    ## T2 ----
    theta0 = theta0[-c(1: (n.gamma + n.bT + n.tk))]
    Xi = datT$X2; deltaT = datT$d2; ZmatT = datT$Zmat2; 
    GfunT = datT$Gfun2; Wmat = datT$Wmat2; copula.link = datT$link2
    control = datT$control1
    
    n.bT = ncol(ZmatT); n.gamma = ncol(Wmat); 
    tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
    Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk)) 
    dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
    copula.para = theta0[1 : n.gamma]
    alphai = copula.link$h.fun(
      as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
    )
    yes.constraint =  yes.constraint |  
      min(alphai) < control$lwr | max(alphai) > control$upr
    if (n.bT == 0){
      dLambdaT = theta0[n.gamma + c(1 : n.tk)]; bzT = rep(0, N) 
    } else{
      bT = theta0[n.gamma + c(1 : n.bT)]; 
      dLambdaT = theta0[n.gamma + n.bT + c(1 : n.tk)] 
      bzT = as.vector(ZmatT %*% bT); 
    }
    dLambdaT.Xi = rep(0, N)
    dLambdaT.Xi[deltaT == 1] = dLambdaT[match(Xi[deltaT == 1], tk)]
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
    G.all = G.funs(GfunT)
    uiT = exp(- G.all$g.fun(LambdaT.Xi * exp(bzT)))
    yes.constraint = yes.constraint | any(dLambdaT < 0) 
    dat2 = list(
      alphai = alphai, ui = uiT, G.all = G.all, Lambda.Xi = LambdaT.Xi, 
      dLambda.Xi = dLambdaT.Xi, dLambda.tk = dLambdaT,
      bz = bzT, Xi.g.tk = Xi.g.tk, dN.tk = dN.tk, Zmat = ZmatT)
    
    out = list(alphai = alphai, yes.constraint = yes.constraint,
               dat1 = dat1, dat2 = dat2)
  } else {
    out = list(alphai = alphai, yes.constraint = yes.constraint)
  } 
  return(out)
}



lln.trivic.PMLE = function(uu1, uu2, uu3, dd1, dd2, dd3, par12, par1, par2, 
                           index12, index1, index2, link12, Wmat12, 
                           yes.constraint, only.ll = FALSE, 
                           pc.method = "foeach"){
  
  if (yes.constraint) {
    if (only.ll) return(list(lln = -Inf)) else {
      return(list(lln = -Inf, dll = NA, ddlln = NA))
    }
  } else {
      N = length(uu1)
      # browser()
      cu1 = BiCopHfunc2(uu1, uu3, family = index1, par = par1)
      cu2 = BiCopHfunc2(uu2, uu3, family = index2, par = par2)
      if (index12 == 4) {
        cu1[cu1 == 1] = exp(- 1 / N)
        cu2[cu2 == 1] = exp(- 1 / N)
      }
      
      ll1.all = llfuns.bvic(
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1], d1 = dd1[dd3 == 1], 
        d2 = dd2[dd3 == 1], copula.index = index12, alphai = par12[dd3 == 1], 
        yes.constraint = yes.constraint, 
        only.ll = only.ll, theta1 = TRUE, dat_bvic =  NULL, 
        copula.link = link12, Wmat = Wmat12[dd3 == 1, ,drop = F], 
        yes.dllDu2 = FALSE)
      lln1.sum = ll1.all$lln * sum(dd3 == 1)
      
      dd = funsData_trivic(
        Funs = density_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      ll0 = log(dd)
      lln = (lln1.sum + sum(ll0)) / N
      
      if (only.ll) {
        out = list(lln = lln)
      } else{
        n.gamma12 = ncol(Wmat12)
        alpha12.lp = link12$hinv.fun(par12[dd3 == 0])
        ## ll.Dgamma12 ----
        ll1.Dgamma12 = ll1.all$dll
        dd.Dpar12 = funsData_trivic(
          Funs = density.Dpar12_trivic, 
          d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
          u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
          alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
          alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
          index2 = index2, integrate = TRUE, pc.method = pc.method)
        ll0.Dpar12 = dd.Dpar12 / dd
        ll0.Dgamma12 = ll0.Dpar12 * link12$dot.h.fun(alpha12.lp) * 
          Wmat12[dd3 == 0, drop = F]
        ll.Dgamma12 = matrix(NA, nrow = N, ncol = n.gamma12)
        ll.Dgamma12[dd3 == 1, ] = ll1.Dgamma12
        ll.Dgamma12[dd3 == 0, ] = ll0.Dgamma12
        
        ## ll.Dgamma12.Dgamma12 ----
        ll1.Dgamma12.Dgamma12 = ll1.all$ddlln * sum(dd3 == 1)
        dd.Dpar12.Dpar12 = funsData_trivic(
          Funs = density.Dpar12.Dpar12_trivic, 
          d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
          u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
          alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
          alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
          index2 = index2, integrate = TRUE, pc.method = pc.method)
        ll0.Dpar12.Dpar12 = dd.Dpar12.Dpar12 / dd - ll0.Dpar12 ^ 2
        ll0.Dgamma12.Dgamma12 = matrix(
          apply((ll0.Dpar12.Dpar12 * (link12$dot.h.fun(alpha12.lp) ^ 2) + 
                   ll0.Dpar12 * link12$ddot.h.fun(alpha12.lp)) * 
                  Wmat12[dd3 == 0, rep(1 : n.gamma12, n.gamma12), drop = F] * 
                  Wmat12[dd3 == 0, rep(1 : n.gamma12, each = n.gamma12), drop = F], 
                2, mean, na.rm = T),
          byrow = T, nrow = n.gamma12, ncol = n.gamma12) * sum(dd3 == 0)
        ll.Dgamma12.Dgamma12 = (ll1.Dgamma12.Dgamma12 + ll0.Dgamma12.Dgamma12) / 
          N
        out = list(lln = lln, dll = ll.Dgamma12, ddlln = ll.Dgamma12.Dgamma12)
      }
  }
  return(out)
}

dll.trivic.PMLE = function(uu1, uu2, uu3, dd1, dd2, dd3, par12, par1, par2, 
                           index12, index1, index2, yes.constraint, 
                           pc.method = "foreach"){
  if (yes.constraint) 
    stop("Some parameter estimates are not in their domain") else {
      N = length(uu1)
      cu1 = BiCopHfunc2(uu1, uu3, family = index1, par = par1)
      cu2 = BiCopHfunc2(uu2, uu3, family = index2, par = par2)
      if (index12 == 4) {
        cu1[cu1 == 1] = exp(- 1 / N)
        cu2[cu2 == 1] = exp(- 1 / N)
      }
      ## density ----
      out1 = funsData_bvic(
        Funs = density_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1])
      out0 = funsData_trivic(
        Funs = density_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd = rep(NA, N); dd[dd3 == 1] = out1; dd[dd3 == 0] = out0
      
      ## Dpar12 ----
      out1 = funsData_bvic(
        Funs = density.Dpar_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) 
      out0 = funsData_trivic(
        Funs = density.Dpar12_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar12 = rep(NA, N); dd.Dpar12[dd3 == 1] = out1; 
      dd.Dpar12[dd3 == 0] = out0
      
      ## Du1 ----
      out1 = funsData_bvic(
        Funs = density.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopPDF(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], 
                 family = index1, par = par1[dd3 == 1])
      out0= funsData_trivic(
        Funs = density.Du1_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Du1 = rep(NA, N); dd.Du1[dd3 == 1] = out1; dd.Du1[dd3 == 0] = out0
      
      
      ## Du2 ----
      out1 = funsData_bvic(
        Funs = density.Du2_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopPDF(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], 
                 family = index2, par = par2[dd3 == 1])
      out0 = funsData_trivic(
        Funs = density.Du2_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Du2 = rep(NA, N); dd.Du2[dd3 == 1] = out1; dd.Du2[dd3 == 0] = out0
      
      ## Dpar1
      out1 = funsData_bvic(
        Funs = density.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "par",
                        family = index1, par = par1[dd3 == 1])
      out0= funsData_trivic(
        Funs = density.Dpar1_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar1 = rep(NA, N); dd.Dpar1[dd3 == 1] = out1; dd.Dpar1[dd3 == 0] = out0
      
      ## Dpar2 ----
      out1 = funsData_bvic(
        Funs = density.Du2_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "par",
                        family = index2, par = par2[dd3 == 1])
      out0= funsData_trivic(
        Funs = density.Dpar2_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar2 = rep(NA, N); dd.Dpar2[dd3 == 1] = out1; dd.Dpar2[dd3 == 0] = out0
      
      ## Dpar12.Dparu1 ----
      out1 = funsData_bvic(
        Funs = density.Dpar.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopPDF(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], 
                 family = index1, par = par1[dd3 == 1])
      out0 = funsData_trivic(
        Funs = density.Dpar12.Du1_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar12.Du1 = rep(NA, N); dd.Dpar12.Du1[dd3 == 1] = out1; 
      dd.Dpar12.Du1[dd3 == 0] = out0
      
      ## Dpar12.Du2 ----
      out1 = funsData_bvic(
        Funs = density.Dpar.Du2_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopPDF(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], 
                 family = index2, par = par2[dd3 == 1])
      out0 = funsData_trivic(
        Funs = density.Dpar12.Du2_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar12.Du2 = rep(NA, N); dd.Dpar12.Du2[dd3 == 1] = out1; 
      dd.Dpar12.Du2[dd3 == 0] = out0
      
      ## Dpar12.Dpar1----
      out1 = funsData_bvic(
        Funs = density.Dpar.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "par",
                        family = index1, par = par1[dd3 == 1])
      out0= funsData_trivic(
        Funs = density.Dpar12.Dpar1_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar12.Dpar1 = rep(NA, N); dd.Dpar12.Dpar1[dd3 == 1] = out1; 
      dd.Dpar12.Dpar1[dd3 == 0] = out0
      
      ## Dpar12.Dpar2 ----
      out1 = funsData_bvic(
        Funs = density.Dpar.Du2_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "par",
                        family = index2, par = par2[dd3 == 1])
      out0 = funsData_trivic(
        Funs = density.Dpar12.Dpar2_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = TRUE, pc.method = pc.method)
      dd.Dpar12.Dpar2 = rep(NA, N); dd.Dpar12.Dpar2[dd3 == 1] = out1; 
      dd.Dpar12.Dpar2[dd3 == 0] = out0
      
      ll.Dpar12.Du1 = dd.Dpar12.Du1 / dd - (dd.Dpar12 / dd) * (dd.Du1 / dd)
      ll.Dpar12.Du2 = dd.Dpar12.Du2 / dd - (dd.Dpar12 / dd) * (dd.Du2 / dd)
      ll.Dpar12.Dpar1 = dd.Dpar12.Dpar1 / dd - (dd.Dpar12 / dd) * (dd.Dpar1 / dd)
      ll.Dpar12.Dpar2 = dd.Dpar12.Dpar2 / dd - (dd.Dpar12 / dd) * (dd.Dpar2 / dd)
      
      ## Du3 ----
      out1 = funsData_bvic(
        Funs = density.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "u2",
                        family = index1, par = par1[dd3 == 1]) +
        funsData_bvic(
          Funs = density.Du2_bvic, 
          d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
          u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
          copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "u2",
                        family = index2, par = par2[dd3 == 1]) 
      out0 = funsData_trivic(
        Funs = density_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = FALSE, pc.method = pc.method)
      dd.Du3 = rep(NA, N); dd.Du3[dd3 == 1] = out1; dd.Du3[dd3 == 0] = out0
      
      ## Dpar12.Du3 ----
      out1 = funsData_bvic(
        Funs = density.Dpar.Du1_bvic, 
        d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
        u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
        copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu1[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "u2",
                        family = index1, par = par1[dd3 == 1]) +
        funsData_bvic(
          Funs = density.Dpar.Du2_bvic, 
          d1 = dd1[dd3 == 1], d2 = dd2[dd3 == 1], 
          u1 = cu1[dd3 == 1], u2 = cu2[dd3 == 1],
          copula.index = index12, para = par12[dd3 == 1]) *
        BiCopHfuncDeriv(u1 = uu2[dd3 == 1], u2 = uu3[dd3 == 1], deriv = "u2",
                        family = index2, par = par2[dd3 == 1]) 
      out0 = funsData_trivic(
        Funs = density.Dpar12_trivic, 
        d1 = dd1[dd3 == 0], d2 = dd2[dd3 == 0], 
        u1 = uu1[dd3 == 0], u2 = uu2[dd3 == 0], u3 = uu3[dd3 == 0], 
        alpha12 = par12[dd3 == 0], alpha1 = par1[dd3 == 0], 
        alpha2 = par2[dd3 == 0], index12 = index12, index1 = index1, 
        index2 = index2, integrate = FALSE, pc.method = pc.method)
      dd.Dpar12.Du3 = rep(NA, N); dd.Dpar12.Du3[dd3 == 1] = out1; 
      dd.Dpar12.Du3[dd3 == 0] = out0
      ll.Dpar12.Du3 = dd.Dpar12.Du3 / dd - (dd.Dpar12 / dd) * (dd.Du3 / dd)
      
      return(list(ll.Dpar12.Du1 = ll.Dpar12.Du1, 
                  ll.Dpar12.Du2 = ll.Dpar12.Du2,
                  ll.Dpar12.Du3 = ll.Dpar12.Du3,
                  ll.Dpar12.Dpar1 = ll.Dpar12.Dpar1, 
                  ll.Dpar12.Dpar2 = ll.Dpar12.Dpar2))
      
      
    }
  
}


llfuns.trivic.MLE = function(uu1, uu2, uu3, dd1, dd2, dd3, par12, par1, par2, 
                             index12, index1, index2, link12, Wmat12, 
                             dat1, dat2, yes.constraint, pc.method = "foreach"){
  
  if (yes.constraint)
    return(list(lln = -Inf, dll = NA, ddlln = NA)) else {
      N = length(uu1)
      cu1 = BiCopHfunc2(uu1, uu3, family = index1, par = par1)
      cu2 = BiCopHfunc2(uu2, uu3, family = index2, par = par2)
      if (index12 == 4) {
        cu1[cu1 == 1] = exp(- 1 / N)
        cu2[cu2 == 1] = exp(- 1 / N)
      }
      
      
    }
}
