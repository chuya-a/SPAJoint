### Here is an example of a simple application.

N <- 10000
nSNP <- 1000
maf <- 0.3
betaR <- as.matrix(c(0.5))
betaS <- as.matrix(c(0.5))
gammaR <- as.matrix(rep(-0.2,nSNP),ncol=nSNP)
gammaS <- as.matrix(rep(-0.8,nSNP),ncol=nSNP)
theta1 <- 0.5
theta2 <- 0.5
ro <- 0.6

# NOTE: The row and column names of genotype matrix are required.
Geno.mtx = matrix(rbinom(N * nSNP,2,maf),N,nSNP)
rownames(Geno.mtx) = paste0("IID-",1:N)
colnames(Geno.mtx) = paste0("SNP-",1:nSNP)

d <- mvrnorm(N, miu, sigma)
u1 <- as.matrix(d[,1])
v1 <- as.matrix(d[,2])
yetaS <- X %*% betaS + z1 %*% u1
ksiR <- X %*% betaR + z1 %*% v1
mu <- as.vector(exp(ksiR) / (1+exp(ksiR)))
Phen.mtx <- data.frame(ID=paste0("IID-",1:N),
                       t=runif(N),
                       event=rbinom(N,1,0.5),
                       R=rbinom(N,1,mu),
                       X=rnorm(N))
I <- t(M)
pij <- matrix(rep(0,N*N),ncol=N)
for(i in 1:N){
  for(j in 1:N){
  pij[i,j] <- (M[i,j] * event[j] * exp(yetaS[i])) / (I[j,] %*% exp(yetaS))
  }
}
lambda <- as.matrix(apply(pij,1,sum))
mresidS <- as.vector(Phen.mtx$event - lambda)
mresidR <- Phen.mtx$R - mu

re <- JointSPA_Null(mresidS=mresidS,mresidR=mresidR,mu=mu,data=Phen.mtx,
                      pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
obj.jointSPA <- JointSPA(re=re,Geno.mtx=Geno.mtx)
head(obj.jointSPA)
