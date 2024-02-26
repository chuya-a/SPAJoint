setwd("D:/JointSPA")
library("MASS")
library("survival")
library("cubature")
library("dplyr")
# source("JointSPA.R")
source("recorrect_simu.R")


###simulate coefficients
N <- 100
nSNP <- 1
MAF <- 0.3 #0.3, 0.01, 0.001
num <- 500 #循环的次数

flag.para0 <- flag.var0 <- rep(0,num)
betaR <- betaS <- rep(0,num)
theta1 <- theta2 <- ro <- rep(0,num)
gammaR <- gammaS<- matrix(rep(0,num*nSNP),ncol=nSNP)
se.betaR <- se.betaS <- rep(0,num)
se.theta1 <- se.theta2 <- se.ro <- rep(0,num)
se.gammaR <- se.gammaS <- matrix(rep(0,num*nSNP),ncol=nSNP)
ksiR <- yetaS <- matrix(rep(0,N*num),ncol=num)
# z1 <- diag(rep(1, N))
# 分为100个组，每组1人
z1 <-  rep(1, 1)
for(i in 2:N){
  z1 <- c(z1, rep(0, N), rep(1, 1))
}
z1 = matrix(z1, ncol = N)
count <- 0

betaR.simu <- as.matrix(rep(0.8))
gammaR.simu <- as.matrix(rep(-0.2,nSNP),ncol=nSNP)
betaS.simu <- as.matrix(rep(1.2))
gammaS.simu <- as.matrix(rep(-0.8,nSNP),ncol=nSNP)
theta1.simu <- 0.5
theta2.simu <- 0.5
ro.simu <- 0.2

mu <- rep(0, 2*N)
sigma <- rbind(cbind(diag(rep(theta1.simu,N)), diag(rep(ro.simu*sqrt(theta1.simu*theta2.simu),N))),
               cbind(diag(rep(ro.simu*sqrt(theta1.simu*theta2.simu),N)),diag(rep(theta2.simu,N))))

t1 <- Sys.time()
j <- 1
repeat
{
  # 基因型矩阵如何设置才能用于求线性预测
  # 把基因型矩阵中的AA，Cc转换为基因含量矩阵，这些参数用基因型中最小频率等位基因
  # 的个数表示，即0,1,2，该矩阵叫做MAF矩阵
  Geno.mtx <- matrix(rbinom(N * nSNP,2,MAF),N,nSNP)
  # rownames(Geno.mtx) = paste0("IID-",1:N)
  # colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
  
  x.simu <- runif(N,0,1)
  # x.simu = as.vector(sample(c(rep(0,N/2),rep(1,N/2))))
  d <- mvrnorm(1, mu, sigma)
  u1 <- as.matrix(d[1:N])
  v1 <- as.matrix(d[N+1:N])
  yetaS.simu <- x.simu %*% betaS.simu + Geno.mtx %*% gammaS.simu + z1 %*% u1 
  ksiR.simu <- x.simu %*% betaR.simu + Geno.mtx %*% gammaR.simu + z1 %*% v1
  
  T.simu <- rep(0,N)
  C.simu <- runif(N,50,500)
  for(i in 1:N)
  {
    T.simu[i] <- -500 * log(1 - runif(1,0,1)) / exp(yetaS.simu[i])
  }
  
  delta <- ifelse(T.simu <= C.simu,1,0) #1：观测到失败事件
  t.observe <- ifelse(delta == 1,T.simu,C.simu) #观测到的事件时间
  
  P0 <- 1 / (1 + exp(ksiR.simu)) #R=0的概率
  u <- runif(N,0,1) #logit(miu)中的miu
  R.simu <- rep(0,N) #tumor response
  
  for(i in 1:N)
  {
    if(u[i] > P0[i])
    {R.simu[i] <- 1}
    # if((u[i] <= P0[i]))
    # {R.simu[i] <- 0}
  }

  # ID <- paste0("IID-",1:N)
  data <- cbind(t.observe,delta,R.simu,x.simu,Geno.mtx)
  cat("j=",j,'\n')
  
  result <- Coef_simu(data=data,nSNP=nSNP,theta1=0.2,theta2=1,ro=0.3,itmax=1000,epsilon=0.001)
  
  if(result$flag.var == 0)
  {
    cat("flag.var=",result$flag.var,"\n")
    next
  }
  
  betaR[j] <- result$betaR
  gammaR[j,] <- result$gammaR
  betaS[j] <- result$betaS
  gammaS[j,] <- result$gammaS
  theta1[j] <- result$theta1
  theta2[j] <- result$theta2
  ro[j] <- result$ro
  # yetaR[,j] <- result$yetaR
  # yetaS[,j] <- result$yetaS

  se.betaR[j] <- result$se.betaR
  se.gammaS[j,] <- result$se.gammaS
  se.betaS[j] <- result$se.betaS
  se.gammaS[j,] <- result$se.gammaS
  se.theta1[j] <- result$se.theta1
  se.theta2[j] <- result$se.theta2
  se.ro[j] <- result$se.ro
  
  flag.para0[j] <- result$flag.para
  flag.var0[j] <- result$flag.var
  

  j <- j+1
  
  if(j == num+1) break

}#repeat循环结束

gammaR.var <- as.matrix(apply(gammaR,2,var))
gammaS.var <- as.matrix(apply(gammaS,2,var))

cat("betaR=",mean(betaR),"SE.betaR=",mean(se.betaR),"SD.betaR=",sqrt(var(betaR)),'\n')
cat("betaS=",mean(betaS),"SE.betaS=",mean(se.betaS),"SD.betaS=",sqrt(var(betaS)),'\n')
cat("gammaR=",apply(gammaR,2,mean),"SE.gammaR=",apply(se.gammaR,2,mean),"SD.gammaR=",apply(gammaR.var,2,sqrt),'\n')
cat("gammaS=",apply(gammaS,2,mean),"SE.gammaS=",apply(se.gammaS,2,mean),"SD.gammaS=",apply(gammaS.var,2,sqrt),'\n')
cat("theta1=",mean(theta1),"SE.theta1=",mean(se.theta1),"SD.theta1=",sqrt(var(theta1)),'\n')
cat("theta2=",mean(theta2),"SE.theta2=",mean(se.theta2),"SD.theta2=",sqrt(var(theta2)),'\n')
cat("ro=",mean(ro),"SE.ro=",mean(se.ro),"SD.ro=",sqrt(var(ro)),'\n')
cat("count=",result$count_ro,'\n')
# cat("yetaR=",apply(yetaR,1,mean),'\n')
# cat("yetaS=",apply(yetaS,1,mean),'\n')

t2 <- Sys.time()
t <- t2-t1
print(t)