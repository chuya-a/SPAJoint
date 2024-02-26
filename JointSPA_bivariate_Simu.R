setwd("D:/JointSPA")
library(MASS)
library(xlsx)
# library(seqminer)
library(survival)
library(data.table)
library(cubature)
library(dplyr)
library(BB)
library(Matrix)
# library(rootSolve)
# library(akima)
# library(nlme)
# library(mgcv)
library(fMultivar)
library(SPAtest)
library(nleqslv)
source("logit&cox-Esti.R",encoding = "utf-8")
source("jointSPA_bivariate.mod.mod.R",encoding = "utf-8")

N <- 10000
nSNP <- 10000
maf <- 0.01 # 0.3, 0.1, 0.01
event_rate <- 0.002 # 0.2%, 1%, 10%, 20%, 50%
lamda <- c(0.168918,0.167564,0.152331,0.135406,0.084628)
num <- 1 #重复次数
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

# ri、ci中的坐标是一一对应的
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
# M <- lower.tri(matrix(rep(1, N*N), ncol=N), diag=TRUE)
# M <- ifelse(M == TRUE, 1, 0)
# M <- as.matrix(M) #下三角矩阵

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
  # event <- ifelse(Ti <= Ci, 1, 0)
  
  event <- rbinom(N,1,event_rate)
  # t <- runif(N)
  
  # 逻辑函数
  P0 <- 1 / (1 + exp(ksiR))
  u <- runif(N,0,1)
  R <- rep(0,N) #tumor response
  # for(l in 1:N)
  # {
  #   if(u[l] > P0[l])
  #   {R[l] <- 1}
  # }
  
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
  
  # I <- t(M) #上三角函数
  # pij <- matrix(rep(0,N*N),ncol=N)
  # for(i in 1:N){
  #   for(j in 1:N){
  #     pij[i,j] <- (M[i,j] * event[j] * exp(yetaS[i])) / (I[j,] %*% exp(yetaS))
  #   }
  # }
  # lambda <- as.matrix(apply(pij,1,sum))
  # mresidS <- as.vector(Phen.mtx$event - lambda)
  mresidR <- Phen.mtx$R - mu
  
  # glm.obj <- glm(y~X+v1,data=Phen.mtx,family = binomial)
  # mu1 <- glm.obj$fitted.values
  # R1 <- glm.obj$y
  
  # source("jointSPA_bivariate.mod.mod.R",encoding = "utf-8")
  # 不使用插值函数，直接采用公式计算联合CGFs
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
  
  # for(i in 1:nSNP){
  #   for(j in 1:10){
  #     if(pval[i,j] < 0.01)count.spa = count.spa+1
  #   }
  # }
  # count.spa = 0
  # 当ro=0时，计算独立时的联合概率值
  # source("jointSPA_ro0.mod.R",encoding = "utf-8")
  # re.ro0 <- JointSPA_Null.ro0(mresidS=mresidS,mresidR=mresidR,mu=mu,data=Phen.mtx,
  #                                 pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
  # obj.jointSPA.ro0 <- JointSPA.ro0(re=re.ro0,Geno.mtx=Geno.mtx)
  # pval.ro0[,k] <- obj.jointSPA.ro0[,3]
  # 
  # print(head(obj.jointSPA.ro0))
  # for(l in 1:nSNP){
  #   if(obj.jointSPA.ro0[l,3] < 0.05){
  #     count.ro0 = count.ro0 + 1
  #   }
  # }
  
  # 对cox和logit分别用SPA计算P值
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
  
  # 当X，Y相互独立时，不论是连续的还是离散的，
  # 二元的鞍点近似的概率值都是一元的两个概率值的乘积
  # 但我的这两个P值是有相关性的，X,Y不是相互独立的
  # 那么，Pval.indep <= pvalS * pvalR (当且仅当相互独立时等号成立)
  pval.ro0[,k] <- pvalS[,k]*pvalR[,k]
  # pval.ro[,k] <- min(1, 2*pvalS[,k]*pvalR[,k])
  
  output = matrix(NA, nSNP, 11)
  output = cbind(P.cox[,1],P.cox[,2],P.cox[,3],P.logit[,3],pval.ro0[,k],P.cox[,4],
                  P.logit[,4],P.cox[,5],P.logit[,5],P.cox[,6],P.logit[,6])
  colnames(output) = c("MAF","missing.rate","pvalS","pvalR","pval.ro",
                       "StatS","StatR","VarS","VarR","zS","zR")
  rownames(output) = paste0("SNP-",1:nSNP)
  
  for(l in 1:300){
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
  # l <- c('X','ksiR','yetaS','event','Ci','Ui','Ti','t','Phen.mtx','Geno.mtx',
  #        'pij','lambda','re','re.cox','re.logit')
  # rm(list = l)
}
sink()
t2 <- Sys.time()

# count.spa <- 0
# count.spa0 <- 0
# count.spa1 <- 0
# count.ro1 <- 0
# count.ro2 <- 0
# count.ro0 <- 0
# countS <- 0
# countR <- 0
# 
# for(i in 1:nSNP){
#   for(j in 1:num){
#     if(pval[i,j] < 0.01)count.spa = count.spa + 1
#     if(pval[i,j] == 0)count.spa0 = count.spa0 + 1
#     if(pval[i,j] == 1)count.spa1 = count.spa1 + 1
#     if(pvalR[i,j] < 0.05)countR = countR + 1
#     if(pvalS[i,j] < 0.05)countS = countS + 1
#     if((pvalR[i,j] < 0.025) & (pvalS[i,j] < 0.025)){
#       count.ro1 = count.ro1 + 1}
#     if((pvalR[i,j] < 0.025) | (pvalS[i,j] < 0.025)){
#       count.ro2 = count.ro2 + 1}
#     if(pval.ro0[i,j] < 0.05)count.ro0 = count.ro0 + 1    
#   }
# }

# 二元鞍点近似
# count.spa <- count.spa1 - count.spa2
# res.spa <- c(count.spa/(num*nSNP-count.spa2-count.spa3))
# cat("\nP.spa=",res.spa)
# cat("\ncount.spa=", count.spa)
res.p <- c(count.spa/(num*nSNP))
cat("\nP.spa=",res.p)
cat("\ncount.spa=", count.spa, "\ncount.spa0=",count.spa0)
# 单独的鞍点近似
res.p1 <- c(countS/(num*nSNP), countR/(num*nSNP))
cat("\nP.cox=",res.p1[1], ", P.logit=",res.p1[2])
cat("\ncountS=", countS, ", countR=", countR)

res.p.ro <- c(count.ro1,count.ro2)/(num*nSNP)
cat("P.ro1=",res.p.ro[1], ", count.ro1=", count.ro1)
cat("P.ro2=",res.p.ro[2], ", count.ro2=", count.ro2)

# # res.p.ro0 <- count.ro0/(num*nSNP)
# # cat("P.ro0=",res.p.ro0, ", count.ro0=", count.ro0)
# 
# write.xlsx(pval,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/ro=0.1.xlsx",sheetName = "pval")
# # write.xlsx(pval.ro0,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0/ro=0.8.xlsx",sheetName = "pval.ro0")
# write.xlsx(pvalR,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/ro=0.1.xlsx",sheetName = "pvalR",append = TRUE)
# write.xlsx(pvalS,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/ro=0.1.xlsx",sheetName = "pvalS",append = TRUE)
# write.xlsx(pval.ro,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/ro=0.8.xlsx",sheetName = "pval.ro",append = TRUE)
# 
# write.xlsx(obj.jointSPA,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.2/0.8.obj.xlsx",sheetName = "obj.jointSPA")
# # write.xlsx(obj.jointSPA.ro0,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0/0.8.obj.xlsx",sheetName = "obj.jointSPA.ro0")
# write.xlsx(P.cox,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/0.8.obj.xlsx",sheetName = "P.cox",append = TRUE)
# write.xlsx(P.logit,file= "D:/JointSPA/re <- JointSPA_Null(mresidS=mresidS,mresidR=mresidR,mu=mu,data=Phen.mtx,
pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
obj.jointSPA <- JointSPA(re=re,Geno.mtx=Geno.mtx)
对比实验/MAF=0.3,event=0.5/ro=0.8/0.8.obj.xlsx",sheetName = "P.logit",append = TRUE)
# write.xlsx(output,file= "D:/JointSPA/对比实验/MAF=0.3,event=0.5/ro=0.8/0.8.obj.xlsx",sheetName = "output",append = TRUE)

print(t2-t1)
