setwd("D:/JointSPA")
library(MASS)
library(xlsx)
library(survival)
library(data.table)
library(cubature)
library(dplyr)
library(BB)
library(Matrix)
library(fMultivar)
library(SPAtest)
library(nleqslv)
source("logit&cox-Esti.R",encoding = "utf-8")
source("jointSPA_bivariate.R",encoding = "utf-8")

N <- 10000
nSNP <- 10000
maf <- 0.01 # MAF = 0.3, 0.1, 0.01
event_rate <- 0.002 # 0.2%, 1%, 10%, 20%, 50%
lamda <- c(0.168918,0.167564,0.152331,0.135406,0.084628)
num <- 1 
betaR <- as.matrix(c(0.5))
betaS <- as.matrix(c(0.5))
gammaR <- as.matrix(rep(-0.2,nSNP),ncol=nSNP)
gammaS <- as.matrix(rep(-0.8,nSNP),ncol=nSNP)
theta1 <- 0.5
theta2 <- 0.5
ro <- 0.6

count.spa <- 0
count.spa0 <- 0
count.spa1 <- 0
count.ro1 <- 0
count.ro2 <- 0
count.ro0 <- 0
countS <- 0
countR <- 0
pval <- matrix(rep(0,num*nSNP),ncol=num)
pvalS <- matrix(rep(0,num*nSNP),ncol=num)
pvalR <- matrix(rep(0,num*nSNP),ncol=num)
pval.ro0 <- matrix(rep(0,num*nSNP),ncol=num)

ri <- rep(0,N)
ci <- rep(0,N)
for(i in 1:N){
  ri[i] <- i
  ci[i] <- i
}
z1 <- sparseMatrix(ri,ci,x=1)
miu <- c(0, 0)
sigma <- matrix(c(theta1,ro*sqrt(theta1*theta2),
               ro*sqrt(theta1*theta2),theta2),nrow=2)

t1 <- Sys.time()
sink(file="C:/Users/Chuya/Desktop/output.txt")
for(k in 1:num)
{
  print(k)
  
  Geno.mtx = matrix(rbinom(N * nSNP,2,maf),N,nSNP)
  rownames(Geno.mtx) = paste0("IID-",1:N)
  colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
  
  X <- rnorm(N)
  d <- mvrnorm(N, miu, sigma)
  u1 <- as.matrix(d[,1])
  v1 <- as.matrix(d[,2])
  yetaS <- X %*% betaS + z1 %*% u1
  ksiR <- X %*% betaR + z1 %*% v1
  Ci <- rweibull(N,1,scale=0.15)
  Ui <- runif(N,0,1)
  Ti <- lamda[1]*sqrt((-log(Ui))/exp(yetaS))
  t <- as.vector(pmin(Ti,Ci))
  event <- rbinom(N,1,event_rate)
  
  P0 <- 1 / (1 + exp(ksiR))
  u <- runif(N,0,1)
  R <- rep(0,N)
  
  mu <- exp(ksiR) / (1+exp(ksiR))
  mu <- as.vector(exp(ksiR) / (1+exp(ksiR)))
  R = as.vector(rbinom(N,1,mu))
  
  ID <- paste0("IID-",1:N)
  Phen.mtx <- data.frame(ID=ID,
                         t=t,
                         event=event,
                         R=R,
                         X=X)
  cox.obj <- coxph(Surv(t,event) ~ X + u1, data = Phen.mtx, x=T)
  mresidS <- cox.obj$resid
  mresidR <- Phen.mtx$R - mu

  re <- JointSPA_Null(mresidS=mresidS,mresidR=mresidR,mu=mu,data=Phen.mtx,
                      pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
  obj.jointSPA <- JointSPA(re=re,Geno.mtx=Geno.mtx)
  pval[,k] <- obj.jointSPA[,3]
  print(head(obj.jointSPA))
  for(l in 1:nSNP){
    if(obj.jointSPA[l,3] < 0.05){
      count.spa = count.spa + 1
    }
    if(obj.jointSPA[l,3] == 0){
      count.spa0 = count.spa0 + 1
    }
    if(obj.jointSPA[l,3] == 1){
      count.spa1 = count.spa1 + 1
    }
  }
  
  re.cox <- Simu_cox0(mresid = mresidS, data=Phen.mtx,
                      pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
  P.cox <- Simu_cox1(re.cox, Geno.mtx=Geno.mtx)
  cat("Cox\n")
  print(head(P.cox))
  pvalS[,k] <- P.cox[,3]

  re.logit <- Simu_logit0(mu = mu,mresidR = mresidR, data=Phen.mtx)
  P.logit <- Simu_logit1(re.logit, Geno.mtx=Geno.mtx)
  cat("Logit\n")
  print(head(P.logit))
  pvalR[,k] <- P.logit[,3]
  
  pval.ro0[,k] <- pvalS[,k]*pvalR[,k]
  
  output = matrix(NA, nSNP, 11)
  output = cbind(P.cox[,1],P.cox[,2],P.cox[,3],P.logit[,3],pval.ro0[,k],P.cox[,4],
                  P.logit[,4],P.cox[,5],P.logit[,5],P.cox[,6],P.logit[,6])
  colnames(output) = c("MAF","missing.rate","pvalS","pvalR","pval.ro",
                       "StatS","StatR","VarS","VarR","zS","zR")
  rownames(output) = paste0("SNP-",1:nSNP)
  
  for(l in 1:nSNP){
    if(P.cox[l,3] < 0.05){
      countS = countS + 1
    }
    if(P.logit[l,3] < 0.05){
      countR = countR + 1
    }
    if((P.cox[l,3] < 0.025) & (P.logit[l,3] < 0.025)){
      count.ro1 = count.ro1 + 1}
    if((P.cox[l,3] < 0.025) | (P.logit[l,3] < 0.025)){
      count.ro2 = count.ro2 + 1}
    if(pval.ro0[l,k] < 0.05){
      count.ro0 = count.ro0 + 1}
  }
  gc()
}
sink()
t2 <- Sys.time()


res.p <- c(count.spa/(num*nSNP))
cat("\nP.spa=",res.p)
cat("\ncount.spa=", count.spa, "\ncount.spa0=",count.spa0)

res.p1 <- c(countS/(num*nSNP), countR/(num*nSNP))
cat("\nP.cox=",res.p1[1], ", P.logit=",res.p1[2])
cat("\ncountS=", countS, ", countR=", countR)

res.p.ro <- c(count.ro1,count.ro2)/(num*nSNP)
cat("P.ro1=",res.p.ro[1], ", count.ro1=", count.ro1)
cat("P.ro2=",res.p.ro[2], ", count.ro2=", count.ro2)
print(t2-t1)
