### In the SPAJoint_Null method, genetic loci are not considered. Instead, it calculates the moment generating function 
### and cumulative generating function for individual endpoints needed for subsequent saddle point computations.
JointSPA_Null <- function(mresidR = NULL, 
                          mresidS = NULL,
                          mu = NULL,
                          data=NULL,
                          pIDs=NULL,
                          gIDs=NULL,
                          range=c(-100,100),
                          length.out=10000)
{
  obj.check = check_input(pIDs, gIDs, range)
  p2g = obj.check$p2g
  pIDs = obj.check$pIDs
  Cova = data$X
  X1 = cbind(1, Cova)
  
  W = (1-mu)*mu
  WX = W*X1
  XWX_inv = solve(t(X1)%*%WX)
  XW = t(WX)
  y = data$R
  XXWX_inv = X1 %*% XWX_inv
  
  X.invXX = X1 %*% solve(t(X1)%*%X1)
  tX = t(X1)
  
  idx0 = qcauchy(1:length.out / (length.out + 1))
  idx1 = idx0 * max(range) / max(idx0)
  
  indep <- NULL
  print("Start calculating empirical CGF for martingale residuals...")
  c = 0 
  indep <- NULL
  for(i in idx1)
  {
    c = c+1
    u=i
    e_mresidS = exp(mresidS * u)
    M0u = mean(e_mresidS)
    M1u = mean(mresidS * e_mresidS)
    M2u = mean(mresidS^2 * e_mresidS)
    
    e_mresidS_2 = exp(2* mresidS * u)
    M0u_2 = mean(e_mresidS)
    M1u_2 = mean(2*mresidS * e_mresidS)
    M2u_2 = mean(4*mresidS^2 * e_mresidS)
    
    Ku0 = log(M0u)
    Ku1 = M1u/M0u
    Ku2 = (M0u*M2u - M1u^2)/M0u^2
    indep = rbind(indep, c(u, M0u, M1u, M2u, Ku0, Ku1, Ku2))
    if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }
  
  Ku0_emp = approxfun(indep[,1], indep[,5], rule=2)
  Ku1_emp = approxfun(indep[,1], indep[,6], rule=2)
  Ku2_emp = approxfun(indep[,1], indep[,7], rule=2)
  
  var.mresidR = var(mresidR)
  var.mresidS = var(mresidS)
  
  re=list(mresidR = mresidR,
          mresidS = mresidS,
          var.mresidR = var.mresidR,
          var.mresidS = var.mresidS,
          Ku0_emp = Ku0_emp,
          Ku1_emp = Ku1_emp,
          Ku2_emp = Ku2_emp,
          tX = tX,
          X.invXX = X.invXX,
          X1 = X1,
          X0 = Cova,
          mu = mu,
          W = W,
          XW = XW,
          XXWX_inv = XXWX_inv,
          y = y,
          p2g = p2g,
          gIDs = gIDs,
          pIDs = pIDs,
          data = data
  )
  return(re)
}
