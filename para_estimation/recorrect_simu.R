# Parameter estimation function
# logit(mu) = t(x) %*% betaR + Geno.mtx %*% gammaR + u
# lambda(t;x,G) = lambda0(t) * exp(t(x) %*% betaS + Geno.mtx %*% gammaS + v)

Coef_simu <- function(data, nSNP, theta1=theta1, theta2=theta2, ro=ro, itmax=500, epsilon=0.001)
{
  nb <- ncol(data)-3-nSNP
  n <- length(data[,1])
  
  # z is an identity matrix of order n
  z <-  rep(1, 1)
  for(i in 2:n){
    z <- c(z, rep(0, n), rep(1, 1))
  }
  z = matrix(z, ncol = n)
  r <- cbind(data,z)
  r <- r[sort.list(r[,1]),]
  t <- as.vector(r[,1])
  delta <- as.vector(r[,2]) 
  R <- as.vector(ifelse(r[,3] == 1, 1, 0)) 
  x <- as.matrix(r[,4:(3 + nb)]) 
  g <- as.matrix(r[, 3 + nb + 1:nSNP])
  z <- as.matrix(r[, 3 + nb + nSNP + 1:n])
  
  M <- matrix(rep(0,n^2),ncol=n)
  for(k1 in 1:n)
  {
    for(j1 in 1:k1)
    {
      M[k1,j1]<-1
    }
  } 
  one <- rep(1, n)
  tt_event <- cbind(t, delta)
  
  for(k2 in 2:n)
  {
    if(all(tt_event[k2,] == tt_event[k2-1,]))
    {
      M[,k2] <- M[,k2-1]
    }
  }
  
  #Sets the initial value of the parameter
  flag.para <- 0
  flag.var <- 0
  Diverge <- 0
  betaR <- as.matrix(rep(0.1, nb))
  gammaR <- as.matrix(rep(0.1, nSNP))
  betaS <- as.matrix(rep(0.1, nb))
  gammaS <- as.matrix(rep(0.1, nSNP)) 
  u <- as.matrix(rep(0, n))
  v <- as.matrix(rep(0, n))
  iter0.sum <- 0
  iter.sum <- 0
  count_ro <- 0 
  yetaS <- x %*% betaS + g %*% gammaS + z %*% u
  ksiR <- x %*% betaR + g %*% gammaR + z %*% v
 
  for(iter0 in 1:itmax)
  {
    iter0.sum <- iter0.sum + 1
    for(iter in 1:itmax)
    {
      iter.sum <- iter.sum + 1
      
      W1 <- diag(as.vector(exp(yetaS)))
      A <- diag(as.vector(delta / (t(M) %*% exp(yetaS))))
      B <- diag(as.vector(M %*% diag(A)))
      R11 <- W1 %*% B - W1 %*% M %*% A %*% A %*% t(M) %*% W1
      
      W2 <- as.vector(exp(ksiR) / (1 + exp(ksiR))^2)
      R22 <- diag(W2)
    
      x.g <- cbind(x,g)
      
      Info11 <- t(x.g) %*% R11 %*% x.g 
      Info12 <- matrix(rep(0,(nb+nSNP)*(nb+nSNP)), ncol=(nb+nSNP)) 
      Info13 <- t(x.g) %*% R11 %*% z 
      Info14 <- matrix(rep(0,(nb+nSNP)*n), ncol=n) 
      
      Info21 <- t(Info12) 
      Info22 <- t(x.g) %*% R22 %*% x.g 
      Info23 <- matrix(rep(0,(nb+nSNP)*n), ncol=n) 
      Info24 <- t(x.g) %*% R22 %*% z 
      
      Info31 <- t(Info13) 
      Info32 <- t(Info23) 
      Info33 <- t(z) %*% R11 %*% z + diag(rep((1/(theta1 * (1 - ro^2))),n))
      Info34 <- diag(rep((-ro/(sqrt(theta1 * theta2) * (1 - ro^2))),n))
      
      Info41 <- t(Info14) 
      Info42 <- t(Info24) 
      Info43 <- t(Info34)
      Info44 <- t(z) %*% R22 %*% z + diag(rep((1/(theta2 * (1 - ro^2))),n)) 
      
      Info1 <- cbind(Info11, Info12, Info13, Info14)
      Info2 <- cbind(Info21, Info22, Info23, Info24)
      Info3 <- cbind(Info31, Info32, Info33, Info34)
      Info4 <- cbind(Info41, Info42, Info43, Info44)
      Info0 <- rbind(Info1, Info2, Info3, Info4)
      # Info <- solve(Info0)
      Info <- ginv(Info0)
      
      First.ksiR <- as.matrix(R - (exp(ksiR) / (1 + exp(ksiR)))) 
      First.yetaS <- as.matrix(delta - W1 %*% M %*% A %*% one) 
      
      First.betaS <- t(x) %*% First.yetaS 
      First.gammaS <- t(g) %*% First.yetaS 
      First.betaR <- t(x) %*% First.ksiR
      First.gammaR <- t(g) %*% First.ksiR 
      First.u <- t(z) %*% First.yetaS - ((u*theta2 - v*ro*sqrt(theta1*theta2))/(theta1*theta2*(1 - ro^2)))
      First.v <- t(z) %*% First.ksiR - ((v*theta1 - u*ro*sqrt(theta1*theta2))/(theta1*theta2*(1 - ro^2)))
      First <- rbind(First.betaS,First.gammaS,First.betaR,First.gammaR,First.u,First.v)
      para0 <- rbind(betaS, gammaS, betaR, gammaR, u, v)
      para1 <- para0 + Info %*% First
      
      betaS <- as.matrix(para1[1:nb]) 
      gammaS <- as.matrix(para1[nb+1:nSNP]) 
      betaR <- as.matrix(para1[nb+nSNP+1:nb]) 
      gammaR <- as.matrix(para1[2*nb+nSNP+1:nSNP]) 
      u <- as.matrix(para1[2*nb+2*nSNP+1:n]) 
      v <- as.matrix(para1[2*nb+2*nSNP+n+1:n]) 
      
      yetaS <- x %*% betaS + g %*% gammaS + z %*% u
      ksiR <- x %*% betaR + g %*% gammaR + z %*% v
      # cat("para0=",para0,"\npara1=",para1)
      if(max(abs(para1 - para0)) < epsilon)
      {
        flag.para <- 1
        break
      }
      
      if(max(abs(para1 - para0)) > 30)
      {
        cat("break")
        Diverge <- 1
        non <- rep(NA,nSNP)
        result = list(Diverge=Diverge,flag.para=flag.para,flag.var=flag.var,
                      betaR=betaR,gammaR=gammaR,betaS=betaS,gammaS=gammaS,
                      theta1=theta1,theta2=theta2,ro=ro,
                      se.betaR=NA,se.gammaR=non,
                      se.betaS=NA,se.gammaS=non,
                      se.theta1=NA,se.theta2=NA,se.ro=NA,
                      ksiR=ksiR,yetaS=yetaS,u=u,v=v,
                      count_ro=NA)
        return(result)
        break
      }
    }
    
    if(Diverge == 1){
      cat("inner diverge leads to finish")
      break
    }
    
    # cat("betaS=",betaS,"betaR=",betaR,"gammaR=",gammaR,"gammaS=",gammaS)
    # cat("u=",u,"v=",v)
    
    d <- c(u, v)
    T33 <- as.matrix(Info[2*nb + 2*nSNP + 1:(2*n),2*nb + 2*nSNP + 1:(2*n)]) # ***
    J1 <- rbind(cbind(diag(rep(1,n)), diag(rep(0,n))),
                cbind(diag(rep(0,n)), diag(rep(0,n))))
    J2 <- rbind(cbind(diag(rep(0,n)), diag(rep(1,n))),
                cbind(diag(rep(1,n)), diag(rep(0,n))))
    J3 <- rbind(cbind(diag(rep(0,n)), diag(rep(0,n))),
                cbind(diag(rep(0,n)), diag(rep(1,n))))
    L1 <- sum(diag(J1 %*% (T33 + d %*% t(d))))
    L2 <- sum(diag(J2 %*% (T33 + d %*% t(d))))/2
    L3 <- sum(diag(J3 %*% (T33 + d %*% t(d))))
    
    old <- c(theta1,theta2,ro)
    theta1 <- L1 / n
    theta2 <- L3 / n
    ro <- L2 / sqrt(L1 * L3)
    update <- c(theta1,theta2,ro)
    
    if(max(abs(update-old)) > 100)
    {
      cat("outter Diverge\n")
      Diverge<-1
      non <- rep(NA,nSNP)
      result = list(Diverge=Diverge,flag.para=flag.para,flag.var=flag.var,
                    betaR=betaR,gammaR=gammaR,betaS=betaS,gammaS=gammaS,
                    theta1=theta1,theta2=theta2,ro=ro,
                    se.betaR=NA,se.gammaR=non,
                    se.betaS=NA,se.gammaS=non,
                    se.theta1=NA,se.theta2=NA,se.ro=NA,
                    ksiR=ksiR,yetaS=yetaS,u=u,v=v,
                    count_ro=NA)
      return(result)
      break
    }
    if(max(abs(update-old)) < epsilon)
    {
      flag.var <- 1
      break
    }
  }
  cat("iter0.sum=",iter0.sum,"iter.sum=",iter.sum,"\n")
  
  ########## standard error ###########
  if(flag.var == 1)
    {
      ksiR <- x %*% betaR + g %*% gammaR + z %*% v
      yetaS <- x %*% betaS + g %*% gammaS + z %*% u
    
      W1 <- diag(as.vector(exp(yetaS)))
      A <- diag(as.vector(delta / (t(M) %*% exp(yetaS))))
      B <- diag(as.vector(M %*% diag(A)))
      R11 <- W1 %*% B - W1 %*% M %*% A %*% A %*% t(M) %*% W1 
      W2 <- as.vector((exp(ksiR) / (1 + exp(ksiR))^2))
      R22 <- diag(W2) 
      x.g <- cbind(x,g)
    
      Info11 <- t(x.g) %*% R11 %*% x.g
      Info12 <- matrix(rep(0,(nb+nSNP)*(nb+nSNP)), ncol=(nb+nSNP))
      Info13 <- t(x.g) %*% R11 %*% z
      Info14 <- matrix(rep(0,(nb+nSNP)*n), ncol=n) 
      Info21 <- t(Info12)
      Info22 <- t(x.g) %*% R22 %*% x.g
      Info23 <- matrix(rep(0,(nb+nSNP)*n), ncol=n)
      Info24 <- t(x.g) %*% R22 %*% z
      Info31 <- t(Info13)
      Info32 <- t(Info23)
      Info33 <- t(z) %*% R11 %*% z + diag(rep((1/(theta1 * (1 - ro^2))),n))
      Info34 <- diag(rep((-ro/(sqrt(theta1 * theta2) * (1 - ro^2))),n))
      Info41 <- t(Info14)
      Info42 <- t(Info24)
      Info43 <- t(Info34)
      Info44 <- t(z) %*% R22 %*% z + diag(rep((1/(theta2 * (1 - ro^2))),n))
      Info1 <- cbind(Info11, Info12, Info13, Info14)
      Info2 <- cbind(Info21, Info22, Info23, Info24)
      Info3 <- cbind(Info31, Info32, Info33, Info34)
      Info4 <- cbind(Info41, Info42, Info43, Info44)
      Info0 <- rbind(Info1, Info2, Info3, Info4)
      Info <- ginv(Info0)
    
      se.betaS <- sqrt(diag(Info)[1:nb])
      se.gammaS <- sqrt(diag(Info)[nb+1:nSNP])
      se.betaR <- sqrt(diag(Info)[nb+nSNP+1:nb])
      se.gammaR <- sqrt(diag(Info)[2*nb+nSNP+1:nSNP])

    Info.random <- Info[2*nb + 2*nSNP + 1:(2*n),2*nb + 2*nSNP + 1:(2*n)]
    
    Omega <- rbind(cbind(diag(rep(theta1,n)),diag(rep(ro * sqrt(theta1 * theta2),n))),
                   cbind(diag(rep(ro * sqrt(theta1 * theta2),n)),diag(rep(theta2,n))))
    
    Omega.inverse <- (rbind(cbind(diag(rep(theta2,n)),diag(rep(-ro * sqrt(theta1 * theta2),n))),
                            cbind(diag(rep(-ro * sqrt(theta1 * theta2),n)),diag(rep(theta1,n)))) / 
                        (theta1 * theta2 * (1 - ro^2)))
    
    I_n <- diag(rep(1, n))
    I_0 <- diag(rep(0, n))
  
    ORo11 <- sqrt(theta1 * theta2)
    ORo12 <- 2 * theta1 * theta1 * theta2 * (1 - ro * ro)
    Omega1_theta1_11 <- diag(rep((-2 * theta2 / ORo12),n))
    Omega1_theta1_12 <- diag(rep(ro * ORo11 / ORo12,n))
    Omega1_theta1_21 <- diag(rep(ro * ORo11 / ORo12,n))
    Omega1_theta1_1 <- cbind(Omega1_theta1_11, Omega1_theta1_12)
    Omega1_theta1_2 <- cbind(Omega1_theta1_21, I_0)
    Omega1_theta1 <- rbind(Omega1_theta1_1, Omega1_theta1_2) 
    
    ORo21 <- sqrt(theta1 * theta2)
    ORo22 <- 2 * theta1 * theta2 * theta2 * (1 - ro * ro)
    Omega1_theta2_12 <- diag(rep(ro * ORo21 / ORo22,n))
    Omega1_theta2_21 <- diag(rep(ro * ORo21 / ORo22,n))
    Omega1_theta2_22 <- diag(rep(-2 * theta1 / ORo22,n))
    Omega1_theta2_1 <- cbind(I_0, Omega1_theta2_12)
    Omega1_theta2_2 <- cbind(Omega1_theta2_21, Omega1_theta2_22)
    Omega1_theta2 <- rbind(Omega1_theta2_1, Omega1_theta2_2)
    
    ORo31 <- 1 + ro^2
    ORo32 <- sqrt(theta1 * theta2)
    ORo33 <- theta1 * theta2 * (1 - ro * ro)^2
    Omega1_ro11 <- diag(rep((2 * ro * theta2)/ORo33,n))
    Omega1_ro12 <- diag(rep((-ORo31 * ORo32)/ORo33,n))
    Omega1_ro21 <- diag(rep((-ORo31 * ORo32)/ORo33,n))
    Omega1_ro22 <- diag(rep((2 * ro * theta1)/ORo33,n))
    Omega1_ro1 <- cbind(Omega1_ro11, Omega1_ro12)
    Omega1_ro2 <- cbind(Omega1_ro21, Omega1_ro22)
    Omega1_ro <- rbind(Omega1_ro1, Omega1_ro2) 
    
    K1 <- T33 %*% Omega1_theta1
    K2 <- Omega %*% Omega1_theta1
    K3 <- T33 %*% Omega1_theta2
    K4 <- Omega %*% Omega1_theta2
    K5 <- T33 %*% Omega1_ro
    K6 <- Omega %*% Omega1_ro
    
    a11 <- sum(diag((K1 - K2) %*% (K1 - K2)))
    a12 <- sum(diag(K1 %*% K3 + K2 %*% K4 - 2 * (K1 %*% K4)))
    a13 <- sum(diag(K1 %*% K5 + K2 %*% K6 - 2 * (K1 %*% K6)))
    a21 <- t(a12)
    a22 <- sum(diag((K3 - K4) %*% (K3 - K4)))
    a23 <- sum(diag(K3 %*% K5 + K4 %*% K6 - 2 * (K3 %*% K6)))
    a31 <- t(a13)
    a32 <- t(a23)
    a33 <- sum(diag((K5 - K6) %*% (K5 - K6)))
    a1 <- cbind(a11, a12, a13)
    a2 <- cbind(a21, a22, a23)
    a3 <- cbind(a31, a32, a33)
    a <- rbind(a1, a2, a3)
    
    var.d <- 2 * ginv(a)
    cat("var=",var.d,'\n')
    
    se.theta1 <- sqrt(var.d[1,1])
    se.theta2 <- sqrt(var.d[2,2])
    se.ro <- sqrt(var.d[3,3])
    
    cat('OVER','\n')
    cat("betaR=",betaR,"se=",se.betaR,"betaR/se=",betaR/se.betaR,"p-value=",2*(1-pnorm(abs(betaR/se.betaR))),'\n')
    cat("gammaR=",gammaR,"se=",se.gammaR,"gammaR/se=",gammaR/se.gammaR,"p-value=",2*(1-pnorm(abs(gammaR/se.gammaR))),'\n')
    cat("betaS=",betaS,"se=",se.betaS,"betaS/se=",betaS/se.betaS,"p-value=",2*(1-pnorm(abs(betaS/se.betaS))),'\n')
    cat("gammaS=",gammaS,"se=",se.gammaS,"gammaS/se=",gammaS/se.gammaS,"p-value=",2*(1-pnorm(abs(gammaS/se.gammaS))),'\n')
    cat("theta1=",theta1,"se=",se.theta1,"theta1/se=",theta1/se.theta1,"p-value=",2*(1-pnorm(abs(theta1/se.theta1))),'\n')
    cat("theta2=",theta2,"se=",se.theta2,"theta2/se=",theta2/se.theta2,"p-value=",2*(1-pnorm(abs(theta2/se.theta2))),'\n')
    cat("ro=",ro,"se=",se.ro,"ro/se=",ro/se.ro,"p-value=",2*(1-pnorm(abs(ro/se.ro))),'\n')
    result <- list(Diverge=Diverge,flag.para=flag.para,flag.var=flag.var,
                   betaR=betaR,gammaR=gammaR,betaS=betaS,gammaS=gammaS,
                   theta1=theta1,theta2=theta2,ro=ro,
                   se.betaR=se.betaR,se.gammaR=se.gammaR,
                   se.betaS=se.betaS,se.gammaS=se.gammaS,
                   se.theta1=se.theta1,se.theta2=se.theta2,se.ro=se.ro,
                   ksiR=ksiR,yetaS=yetaS,u=u,v=v,
                   count_ro=count_ro)
  }else{
    result <- list(flag.var=flag.var)
  }
  return(result)
} 
