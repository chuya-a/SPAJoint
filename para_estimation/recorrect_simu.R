# 参数估计函数
# 随机效应与个体有关，每个个体有不同的随机效应。
# 数据data的格式为时间，是否死亡(事件)，肿瘤响应，协变量，基因型
# logit(mu) = t(x) %*% betaR + Geno.mtx %*% gammaR + u
# lambda(t;x,G) = lambda0(t) * exp(t(x) %*% betaS + Geno.mtx %*% gammaS + v)

Coef_simu <- function(data, nSNP, theta1=theta1, theta2=theta2, ro=ro, itmax=500, epsilon=0.001)
{
  # theta1=1
  # theta2=1
  # ro=0.1
  # nSNP=nCG
  # itmax=500
  # epsilon=0.001
  # data <- as.matrix(data)
  nb <- ncol(data)-3-nSNP # 协变量的个数
  n <- length(data[,1])
  
  # 按照个体来区分随机效应
  # 每个人的随机效应是不同的，z是一个n阶的单位阵
  # z <- diag(rep(1, n))
  z <-  rep(1, 1)
  for(i in 2:n){
    z <- c(z, rep(0, n), rep(1, 1))
  }
  z = matrix(z, ncol = n)
  r <- cbind(data,z)
  # 按照时间重新排序
  r <- r[sort.list(r[,1]),]
  # 从重排后的数据中获取参与计算的向量和矩阵
  t <- as.vector(r[,1]) #time
  delta <- as.vector(r[,2]) #event
  R <- as.vector(ifelse(r[,3] == 1, 1, 0)) #纵向模型中的响应R
  x <- as.matrix(r[,4:(3 + nb)]) #协变量
  g <- as.matrix(r[, 3 + nb + 1:nSNP])
  z <- as.matrix(r[, 3 + nb + nSNP + 1:n])
  
  # 列出计算过程中各种不变的矩阵
  M <- matrix(rep(0,n^2),ncol=n)
  for(k1 in 1:n)
  {
    for(j1 in 1:k1)
    {
      M[k1,j1]<-1
    }
  } #得到下三角矩阵M
  one <- rep(1, n)
  tt_event <- cbind(t, delta)
  
  # 对tied数据调整M矩阵
  for(k2 in 2:n)
  {
    if(all(tt_event[k2,] == tt_event[k2-1,]))
    {
      M[,k2] <- M[,k2-1]
    }
  }
  
  #设置参数的初始值
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
  # cat("1")
  
  
  #迭代过程
  for(iter0 in 1:itmax)# 外循环：估计方差分量 theta1, theta2, ro
  {
    iter0.sum <- iter0.sum + 1
    # ro <- 0
    for(iter in 1:itmax)# 内循环：估计betaR,gammaR,betaS,gammaS,u,v
    {
      # cat("2")
      # ro <- 0
      iter.sum <- iter.sum + 1
      
      W1 <- diag(as.vector(exp(yetaS)))
      A <- diag(as.vector(delta / (t(M) %*% exp(yetaS))))
      B <- diag(as.vector(M %*% diag(A)))
      R11 <- W1 %*% B - W1 %*% M %*% A %*% A %*% t(M) %*% W1 #n*n
      
      W2 <- as.vector(exp(ksiR) / (1 + exp(ksiR))^2)
      R22 <- diag(W2) #n*n
    
      x.g <- cbind(x,g)
      #信息矩阵
      # x: n*nb
      # g: n*nSNP
      # z: n*n
      
      Info11 <- t(x.g) %*% R11 %*% x.g #(nb+nSNP)*(nb+nSNP)
      Info12 <- matrix(rep(0,(nb+nSNP)*(nb+nSNP)), ncol=(nb+nSNP)) #(nb+nSNP)*(nb+nSNP)
      Info13 <- t(x.g) %*% R11 %*% z #(nb+nSNP)*n
      Info14 <- matrix(rep(0,(nb+nSNP)*n), ncol=n) #(nb+nSNP)*n
      
      Info21 <- t(Info12) #(nb+nSNP)*(nb+nSNP)
      Info22 <- t(x.g) %*% R22 %*% x.g #(nb+nSNP)*(nb+nSNP)
      Info23 <- matrix(rep(0,(nb+nSNP)*n), ncol=n) #(nb+nSNP)*n
      Info24 <- t(x.g) %*% R22 %*% z #(nb+nSNP)*n
      
      Info31 <- t(Info13) #n*(nSNP+nb)
      Info32 <- t(Info23) #n*(nb+nSNP)
      Info33 <- t(z) %*% R11 %*% z + diag(rep((1/(theta1 * (1 - ro^2))),n))#n*n ***
      Info34 <- diag(rep((-ro/(sqrt(theta1 * theta2) * (1 - ro^2))),n)) #n*n  ***
      
      Info41 <- t(Info14) #n*(nb+nSNP)
      Info42 <- t(Info24) #n*(nb+nSNP)
      Info43 <- t(Info34) #n*n  ***
      Info44 <- t(z) %*% R22 %*% z + diag(rep((1/(theta2 * (1 - ro^2))),n)) #n*n  ***
      
      Info1 <- cbind(Info11, Info12, Info13, Info14)
      Info2 <- cbind(Info21, Info22, Info23, Info24)
      Info3 <- cbind(Info31, Info32, Info33, Info34)
      Info4 <- cbind(Info41, Info42, Info43, Info44)
      Info0 <- rbind(Info1, Info2, Info3, Info4)
      # cat("3, iter.sum = ", iter.sum)
      # Info <- solve(Info0) #信息矩阵的逆矩阵
      Info <- ginv(Info0)
      
      # (nb+nb+nSNP+nSNP+n+n)*(nb+nb+nSNP+nSNP+n+n)
      
      #线性预测的一阶导
      First.ksiR <- as.matrix(R - (exp(ksiR) / (1 + exp(ksiR)))) #矩阵:n*1
      First.yetaS <- as.matrix(delta - W1 %*% M %*% A %*% one) #矩阵:n*1
      #迭代得到的参数值
      First.betaS <- t(x) %*% First.yetaS #向量:nb*1
      First.gammaS <- t(g) %*% First.yetaS #向量:nSNP*1
      First.betaR <- t(x) %*% First.ksiR #向量:nb*1
      First.gammaR <- t(g) %*% First.ksiR #向量:nSNP*1
      First.u <- t(z) %*% First.yetaS - ((u*theta2 - v*ro*sqrt(theta1*theta2))/(theta1*theta2*(1 - ro^2)))# ***
      First.v <- t(z) %*% First.ksiR - ((v*theta1 - u*ro*sqrt(theta1*theta2))/(theta1*theta2*(1 - ro^2)))# ***
      First <- rbind(First.betaS,First.gammaS,First.betaR,First.gammaR,First.u,First.v)
      para0 <- rbind(betaS, gammaS, betaR, gammaR, u, v)
      para1 <- para0 + Info %*% First
      
      #重置参数的值
      betaS <- as.matrix(para1[1:nb]) #每一个协变量都有一个系数
      gammaS <- as.matrix(para1[nb+1:nSNP]) #每个SNP都有一个gamma
      betaR <- as.matrix(para1[nb+nSNP+1:nb]) #每一个协变量都有一个系数
      gammaR <- as.matrix(para1[2*nb+nSNP+1:nSNP]) #每个SNP都有一个gamma
      u <- as.matrix(para1[2*nb+2*nSNP+1:n]) #每个线性预测都有一个随机效应v  ***
      v <- as.matrix(para1[2*nb+2*nSNP+n+1:n]) # ***
      # nb+1:nb = (nb+1):(nb+nb)
      
      yetaS <- x %*% betaS + g %*% gammaS + z %*% u
      ksiR <- x %*% betaR + g %*% gammaR + z %*% v
      # cat("para0=",para0,"\npara1=",para1)
      # cat("\npara1 - para0", max(abs(para1 - para0)))
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
      # para0 <- para1
    }
    #内层的for循环结束,betaR,gammaR,betaS,gammaS,u,v已经求出来了,用这些求方差分量的参数
    
    if(Diverge == 1){
      cat("inner diverge leads to finish")
      break
    }
    
    # cat("betaS=",betaS,"betaR=",betaR,"gammaR=",gammaR,"gammaS=",gammaS)
    
    # cat("u=",u,"v=",v)
    # cat("4")
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
    # cat("L1=",L1,"L2=",L2,"L3=",L3,'\n')
    # cat("\nmax(abs(update-old))=",max(abs(update-old)))
    
    if(max(abs(update-old)) > 100)# 发散
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
    if(max(abs(update-old)) < epsilon) # 收敛
    {
      flag.var <- 1
      break
    }
    # theta1 <- theta10
    # theta2 <- theta20
    # ro <- ro0
    # cat("iter0=",iter0,"theta1=",theta1,"theta2=",theta2,"ro=",ro,'\n')
  }#外层的for循环结束,所有的参数的估计值都已经求出来了
  cat("iter0.sum=",iter0.sum,"iter.sum=",iter.sum,"\n")
  
  ########## standard error ###########
  if(flag.var == 1)
    {
      ksiR <- x %*% betaR + g %*% gammaR + z %*% v
      yetaS <- x %*% betaS + g %*% gammaS + z %*% u
    
      W1 <- diag(as.vector(exp(yetaS)))
      A <- diag(as.vector(delta / (t(M) %*% exp(yetaS))))
      B <- diag(as.vector(M %*% diag(A)))
      R11 <- W1 %*% B - W1 %*% M %*% A %*% A %*% t(M) %*% W1 #n*n
      W2 <- as.vector((exp(ksiR) / (1 + exp(ksiR))^2))
      R22 <- diag(W2) #n*n
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
    
      # cat("SE.matrix = ", as.numeric(sqrt(diag(Info))))
      # diag取出矩阵的对角线元素并保存成向量的格式
      # SD：各个参数的标准差
      # SE：各个参数的标准误，是多个样本平均数的标准差
      # 但是下面求的是各个参数的标准差？？？
      se.betaS <- sqrt(diag(Info)[1:nb])
      se.gammaS <- sqrt(diag(Info)[nb+1:nSNP])
      se.betaR <- sqrt(diag(Info)[nb+nSNP+1:nb])
      se.gammaR <- sqrt(diag(Info)[2*nb+nSNP+1:nSNP])
    
    #随机效应的部分
    # theta1 <- theta10
    # theta2 <- theta20
    # ro <- ro0
    Info.random <- Info[2*nb + 2*nSNP + 1:(2*n),2*nb + 2*nSNP + 1:(2*n)]
    
    Omega <- rbind(cbind(diag(rep(theta1,n)),diag(rep(ro * sqrt(theta1 * theta2),n))),
                   cbind(diag(rep(ro * sqrt(theta1 * theta2),n)),diag(rep(theta2,n))))
    
    Omega.inverse <- (rbind(cbind(diag(rep(theta2,n)),diag(rep(-ro * sqrt(theta1 * theta2),n))),
                            cbind(diag(rep(-ro * sqrt(theta1 * theta2),n)),diag(rep(theta1,n)))) / 
                        (theta1 * theta2 * (1 - ro^2)))
    
    #随机效应相关
    ##求方差的过程
    I_n <- diag(rep(1, n))
    I_0 <- diag(rep(0, n))
  
    #Omega.inverse对theta1求一阶偏导
    ORo11 <- sqrt(theta1 * theta2)
    ORo12 <- 2 * theta1 * theta1 * theta2 * (1 - ro * ro)
    Omega1_theta1_11 <- diag(rep((-2 * theta2 / ORo12),n))
    Omega1_theta1_12 <- diag(rep(ro * ORo11 / ORo12,n))
    Omega1_theta1_21 <- diag(rep(ro * ORo11 / ORo12,n))
    Omega1_theta1_1 <- cbind(Omega1_theta1_11, Omega1_theta1_12)
    Omega1_theta1_2 <- cbind(Omega1_theta1_21, I_0)
    Omega1_theta1 <- rbind(Omega1_theta1_1, Omega1_theta1_2) 
    
    #Omega.inverse对theta2求一阶偏导
    ORo21 <- sqrt(theta1 * theta2)
    ORo22 <- 2 * theta1 * theta2 * theta2 * (1 - ro * ro)
    Omega1_theta2_12 <- diag(rep(ro * ORo21 / ORo22,n))
    Omega1_theta2_21 <- diag(rep(ro * ORo21 / ORo22,n))
    Omega1_theta2_22 <- diag(rep(-2 * theta1 / ORo22,n))
    Omega1_theta2_1 <- cbind(I_0, Omega1_theta2_12)
    Omega1_theta2_2 <- cbind(Omega1_theta2_21, Omega1_theta2_22)
    Omega1_theta2 <- rbind(Omega1_theta2_1, Omega1_theta2_2)
    
    #Omega.inverse对ro求一阶偏导
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
    
    #估计方差分量的渐进方差
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
    
    # if(var.d[1] < 0 ){
    #   se.theta1 <- 0
    # }else{
      se.theta1 <- sqrt(var.d[1,1])
    # }
    # if(var.d[2] < 0 ){
    #   se.theta2 <- 0
    # }else{
      se.theta2 <- sqrt(var.d[2,2])
    # }
    # if(var.d[3] < 0 ){
    #   se.ro <- 0
    # }else{
      se.ro <- sqrt(var.d[3,3])
    # }
    
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
