JointSPA_Null <- function(mresidR = NULL, 
                          mresidS = NULL,
                          mu = NULL,
                          data=NULL,
                          pIDs=NULL,
                          gIDs=NULL,
                          range=c(-100,100),
                          length.out=10000)
{
  # pIDs=Phen.mtx$ID
  # gIDs=rownames(Geno.mtx)
  # range=c(-100,100)
  # length.out=100
  # data = Phen.mtx
  
  ###check input arguments
  obj.check = check_input(pIDs, gIDs, range)
  p2g = obj.check$p2g
  pIDs = obj.check$pIDs
  Cova = data$X
  X1 = cbind(1, Cova)
  
  # logit模型
  W = (1-mu)*mu
  WX = W*X1
  XWX_inv = solve(t(X1)%*%WX)
  XW = t(WX)
  y = data$R
  XXWX_inv = X1 %*% XWX_inv
  
  # cox模型
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


JointSPA = function(re,
                    Geno.mtx=Geno.mtx,
                    G.model = "Add",
                    impute.method = "fixed",
                    missing.cutoff = 0.15,
                    min.maf = 0.0001,
                    CovAdj.cutoff = 5e-5)
{
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))
  
  ###Prepare the main output data frame
  nSNP = ncol(Geno.mtx)
  result = matrix(NA, nSNP, 15)
  colnames(result) = c("MAF","missing.rate","p1.value.spa","p2.value.spa",
                       "StatR","StatS","VarR","VarS","zR","zS",
                       "para.flag","zeta1","zeta2","ceta_t0","zeta_u0")
  rownames(result) = colnames(Geno.mtx)
  
  ###Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  for(i in 1:nSNP){
    # if(i %% 100 == 0) 
      cat("i = ", i)
    g <- Geno.mtx[,i]
    result.SNP = JointSPA_one_SNP(g,
                                  re,
                                  G.model,
                                  impute.method,
                                  missing.cutoff,
                                  min.maf,
                                  CovAdj.cutoff)
    result[i,] = result.SNP
  }
  
  print("\nAnalysis Complete.")
  print(Sys.time())
  return(result)
}

JointSPA_one_SNP = function(g,
                            re,
                            # S1=S1,
                            # S2=S2,
                            # ro1=ro1,
                            G.model = "Add",
                            impute.method = "fixed",
                            missing.cutoff = 0.15,
                            min.maf = 0.0001,
                            CovAdj.cutoff = 5e-5)
{
  # G.model = "Add"
  # impute.method = "fixed"
  # missing.cutoff = 0.15
  # min.maf = 0.0001
  # CovAdj.cutoff = 5e-5
  # i <- 1
  # g <- Geno.mtx[,2]

  # 根据基因型求得分统计量
  N = length(g)
  MAF = mean(g, na.rm=T)/2
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N
  
  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }
  
  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
    # g = 1-g # 取了个反
  }

  if(G.model=="Add"){}
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf || missing.rate > missing.cutoff)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA))

  if(!is.null(re$p2g))
    g = g[re$p2g]
  
  ###求p值
  ## 分数检验
  #logit和cox的Score统计量
  R = sum(g * re$mresidR)
  S = sum(g * re$mresidS)
  pos1 = which(g != 0)
  # 大于0.5的为pos1

  
  # 基因型矩阵未进行协变量的调整
  # 中心化的基因型矩阵，2*MAF为g的均值
  G11 = g - 2*MAF
  G12 = g - 2*MAF

  R.var1 = sum(G11 * re$W * G11)
  S.var1 = re$var.mresidS * sum(G12^2)
  # 相当于减去均值0再除以标准差，经过处理后数据符合标准正态分布
  # 将要求的概率分布转化为标准正态分布N~(0,1)范围,就可以使用标准正态分布函数求概率值
  R.Z1 = sum(G11 * re$mresidR)
  S.Z1 = sum(G12 * re$mresidS)
  z11 = R/sqrt(R.var1)
  z12 = S/sqrt(S.var1)

  N0 = N-length(pos1)
  G11norm = G11/sqrt(R.var1) # 标准化基因型矩阵(sd=1)
  G12norm = G12/sqrt(S.var1)
  G11N1 = G11norm[pos1]
  G12N1 = G12norm[pos1]
  G11N0 = -2*MAF/sqrt(R.var1)  # since N0=0, this value actually does not matter
  G12N0 = -2*MAF/sqrt(S.var1)

  # G1N1=G11N1
  # G1N0=G11N0
  # G2N1=G12N1
  # G2N0=G12N0
  # N1set=pos1
  # r=abs(z11)
  # s=abs(z12)
  # G1=G11
  # G2=G12
  # cat("未调整基因型G11=",G11, "\n未调整基因型G12=",G12, "N0=",N0,"\n")
  pval1 = GetProb_SPA1(re=re,G1N1=G11N1,G1N0=G11N0,G2N1=G12N1,G2N0=G12N0,
                       N0=N0,r=abs(z11),s=abs(z12),G1=G11,G2=G12,g=g,pos1=pos1)
  # return(c(MAF, missing.rate, pval1[1], pval1[2], R, S, R.var1, S.var1, z11, z12,
  #          pval1[3],pval1[4],pval1[5],pval1[6],pval1[7]))

  if(pval1[1] > CovAdj.cutoff){
    return(c(MAF, missing.rate, pval1[1], pval1[2], R, S, R.var1, S.var1, z11, z12,
             pval1[3],pval1[4],pval1[5],pval1[6],pval1[7]))
  }
  cat("\nP < 5e-5, 采用协变调整的基因型计算P值")

  # 基因型进行协变调整的情况
  G21 = as.vector(g - re$XXWX_inv %*% (re$XW[,pos1] %*% g[pos1]))
  G22 = as.vector(g - re$X.invXX %*% (re$tX[,pos1,drop=F] %*% g[pos1]))
  
  # G*E中的y.mu = y-mu = mresidR
  R.var2 = sum(G21 * re$W * G21)
  S.var2 = re$var.mresidS * sum(G22^2)
  R.Z2 = sum(G21 * re$mresidR)
  S.Z2 = sum(G22 * re$mresidS)
  z21 = R/sqrt(R.var2)
  z22 = S/sqrt(S.var2)

  G21norm = G21/sqrt(R.var2)
  G22norm = G22/sqrt(S.var2) # 标准化基因型矩阵(sd=1)
  N1set = 1:N
  N0 = 0
  G21N1 = G21norm
  G22N1 = G22norm
  G21N0 = 0 # since N0=0, this value actually does not matter
  G22N0 = 0

  # G1N1=G21N1
  # G1N0=G21N0
  # G2N1=G22N1
  # G2N0=G22N0
  # N1set=N1set
  # r=abs(z21)
  # s=abs(z22)
  # G1=G21
  # G2=G22
  # cat("\n协变调整基因型G21=",G21,"\n协变调整基因型G22=",G22,"N0=",N0,"\n")
  pval2 = GetProb_SPA2(re=re,G1N1=G21N1,G1N0=G21N0,G2N1=G22N1,G2N0=G22N0,
                       pos1=N1set,N0=N0,r=abs(z21),s=abs(z22),G1=G21,G2=G22,g=g)
  return(c(MAF, missing.rate, pval2[1], pval2[2], R, S, R.var2, S.var2, z21, z22,
           pval2[3],pval2[4],pval2[5],pval2[6],pval2[7]))
}

# 未经过协变调整的基因型求P值的函数
GetProb_SPA1 = function(re, G1N1, G1N0, G2N1, G2N0, N0, r, s, G1, G2, g, pos1)
{
  # sd1 = sqrt(var(S1))
  # sd2 = sqrt(var(S2))
  # mu1 = mean(S1)
  # mu2 = mean(S2)
  
  # zeta_t1 = getroot_t2(G1N1=G1N1, G1N0=G1N0, pos1=pos1, N0=N0, re=re, r=r)
  res.roott = getroot_t0(G1N1=G1N1, G1N0=G1N0, pos1=pos1, N0=N0, re=re, r=r)
  if(res.roott$flag.root == 1){
    zeta_t0 = res.roott$zeta
  }else{
    zeta_t0 = getroot_t1(init=0.1,G1N1=G1N1,G1N0=G1N0,G1=G1,re=re,pos1=pos1,N0=N0)$root
  }
  
  res.rootu = getroot_u0(G2N1=G2N1,G2N0=G2N0,pos1=pos1,N0=N0,re=re,s=s)
  if(res.rootu$flag.root == 1){
    zeta_u0 = res.rootu$zeta
  }else{
    zeta_u0 = NaN
    return(c(0,0,0,zeta_t=NaN,zeta_u=NaN,zeta_t0,zeta_u0))
  }
  # 用理论公式求解鞍点
  # init = getroot_K11(G1N1,G1N0,G2N1,G2N0,pos1,N0,re,ro,sd1,sd2,r,s)
  # cat("\n单独的鞍点为：", zeta_u0, zeta_t0, "\n")
  # cat("\n理论公式求得的鞍点为：", init,"\n")
  # 牛顿迭代求解鞍点
  # 检查迭代时的初始值是否满足条件M0>0
  # M01 = (mean(joint_CGF(t=zeta_t0,u=zeta_u0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1)$M0) + 
  #          joint_CGF0(t=zeta_t0,u=zeta_u0,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)$M0)
  # if(M01 != 0){
  #   out = getroot_K1(zeta1=0,zeta2=0,G1N1=G1N1,G2N1=G2N1,G1N0=G1N0,G2N0=G2N0,
  #                    re=re,N0=N0,pos1=pos1,r=r,s=s)
  # }
  mresid1 <- mresidR1(G1N1=G1N1,re,N0,pos1,flag=1)
  mresid2 <- mresidS2(G2N1=G2N1,re,N0,pos1,flag=1)
  mresid <- cbind(mresid1,mresid2)
  rm(list=c('mresid1','mresid2'))
  out = getroot_K1(zeta1=0,zeta2=0,G1N1=G1N1,G2N1=G2N1,G1N0=G1N0,G2N0=G2N0,
                   re=re,N0=N0,pos1=pos1,r=r,s=s,mresid=mresid)
  if(out$flag.para == 1)
  {
    # cat("\n联合鞍点第一组为：",out$root)
    zeta_t = out$root[1]
    zeta_u = out$root[2]
  }
  else{
    # zeta_t = zeta_t0
    # zeta_u = zeta_u0
    pval = c(0,0)
    return(c(pval,out$flag.para,zeta_t=NaN,zeta_u=NaN,zeta_t0,zeta_u0))
  }
  
  # out12 = getroot_K12(zeta1=init[1],zeta2=init[2],G1N1=G1N1,G2N1=G2N1,re=re,r=r,s=s)
  # zeta_t = out12[1]
  # zeta_u = out12[2]
  
  kt0 = (sum(joint_CGF(t=zeta_t0,u=0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)$K0) + 
           N0*joint_CGF0(t=zeta_t0,u=0,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)$K0)
  kt1 = (sum(joint_CGF(t=zeta_t,u=0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)$K0) + 
           N0*joint_CGF0(t=zeta_t,u=0,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)$K0)
  
  ku0 = (sum(joint_CGF(t=0,u=zeta_u0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)$K0) + 
           N0*joint_CGF0(t=0,u=zeta_u0,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)$K0)
  ku1 = (sum(joint_CGF(t=0,u=zeta_u,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)$K0) + 
           N0*joint_CGF0(t=0,u=zeta_u,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)$K0)
  
  joint.CGFs0 = joint_CGF0(t=zeta_t,u=zeta_u,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)
  joint.CGFs1 = joint_CGF(t=zeta_t,u=zeta_u,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)
  rm('mresid1')
  
  k = sum(joint.CGFs1$K0) + N0*(joint.CGFs0$K0)
  ktt = sum(joint.CGFs1$Ktt) + N0*(joint.CGFs0$Ktt)
  ktu = sum(joint.CGFs1$Ktu) + N0*(joint.CGFs0$Ktu)
  kuu = sum(joint.CGFs1$Kuu) + N0*(joint.CGFs0$Kuu)
  kut = sum(joint.CGFs1$Kut) + N0*(joint.CGFs0$Kut)
  k2 = rbind(cbind(ktt,ktu), cbind(kut,kuu))
  rm(list=c('joint.CGFs1','joint.CGFs0'))
  
  #用鞍点分布函数逼近法，也即直接近似得到s,r的分布直接求出概率
  # P(S<s,R<r)=I11=fai(s1,y1,ro)
  # v0 = sign(zeta2)*sqrt(-2*(k - kt0 - (zeta1-zeta_t)*R - zeta2*S))
  # 
  # 求得的根号内的值为负，使得v0=NaN
  # Re():取复数的实部，Im():取复数的虚部，i为虚部单位
  
  wu0_t = sqrt(abs(-2*(k - ku1 - zeta_t*r)))
  wu0_u = sqrt(abs(-2*(k - kt1 - zeta_u*s)))
  # 用较小的那个，大的会使得误差更显著
  if(wu0_t > wu0_u){
    wu0 = wu0_u
    s1 = sign(zeta_u0)*sqrt(abs(-2*(ku0 - zeta_u0*s)))
    v0 = sign(zeta_t)*sqrt(abs(-2*(k - ku0 - (zeta_u-zeta_u0)*s - zeta_t*r)))
    G = 1/sqrt((kuu - ktu^2/ktt))
    b = (wu0-s1)/v0
    r1 = (v0 - b*s1) / sqrt(1+b^2)
    ro1 = -b / sqrt(1+b^2)
    
    I11 = pnorm2d(r1, s1, rho=ro1)
    I12 = pnorm2d(wu0)*dnorm2d(v0)*(1/v0 - 1/(zeta_u*G))
    I21 = pnorm2d(v0)*dnorm2d(r1)*(1/wu0 - 1/(zeta_t * sqrt((ktt))))
    # I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_u * 1/sqrt(kuu)) * (1/v0 - 1/(zeta_t*G)) / (2*pi)
    I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_t * 1/sqrt((ktt))) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
    h2 = -(s1^2 - 2* ro1*s1*r1 + r1^2) / (2*(1-ro1^2))
    I31 = exp(h2)* (1/wu0 - 1/zeta_t * 1/sqrt((ktt))) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
  }else{
    wu0 = wu0_t
    r1 = sign(zeta_t0)*sqrt(abs(-2*(kt0 - zeta_t0*r)))
    v0 = sign(zeta_u)*sqrt(abs(-2*(k - kt0 - (zeta_t-zeta_t0)*r - zeta_u*s)))
    G = 1/sqrt((kuu - ktu^2/ktt))
    b = (wu0-r1)/v0
    s1 = (v0 - b*r1) / sqrt(1+b^2)
    ro1 = -b / sqrt(1+b^2)
    
    I11 = pnorm2d(r1, s1, rho=ro1)
    I12 = pnorm2d(wu0)*dnorm2d(v0)*(1/v0 - 1/(zeta_u*G))
    I21 = pnorm2d(v0)*dnorm2d(r1)*(1/wu0 - 1/(zeta_t * sqrt((ktt))))
    # I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_t * 1/sqrt(ktt)) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
    I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_t * 1/((ktt))) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
    h2 = -(r1^2 - 2* ro1*r1*s1 + s1^2) / (2*(1-ro1^2))
    I31 = exp(h2)* (1/wu0 - 1/zeta_t * 1/sqrt((ktt))) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
  }
  # wu0 = sign(zeta_t)*sqrt(abs(-2*(k - ku1 - zeta_t*r)))
  # cat("\nr1=",r1,",s1=",s1,",v0=",v0,",wu0=",wu0,",ro1=",ro1,",b=",b,",G=",G)
  
  # 用fMultivar包，计算二维正态分布的概率密度函数和分布函数
  # pnorm2d：求的是标准二元正态分布函数的概率值均值为0，方差为1
  # ktt_e = sqrt(abs(ktt))
  # cat("\nI11=",I11,",I12=",I12,",I21=",I21,",I22=",I22,",I31=",I31)
  pval1 = I11+I12+I21+I22
  pval2 = I11+I12+I21+I31
  if(is.na(pval1)){
    pval1 = 0
  }
  if(!is.na(pval1)){
    if(pval1 > 1) pval1 = 1
    if(pval1 < 0) pval1 = 0
    if(pval2 > 1) pval2 = 1
    if(pval2 < 0) pval2 = 0
  }
  if(pval1 >= 5e-5 & pval1 <= 1)pval1=pval1
  if(pval2 >= 5e-5 & pval2 <= 1)pval1=pval2
  
  # pval3 = pnorm2d(r/sd1, s/sd2, rho=ro) # 相当于采用了正态方法去计算联合概率值
  cat("\npval1=",pval1,",pval2=",pval2,"\n")
  
  #用鞍点密度积分法计算p值(由s,r的概率密度函数积分求得概率值)
  # f = function(x){
  #   1/(2*pi*sqrt(det(k2))) * exp(k - zeta_t*x[1] - zeta_u*x[2])
  # }
  # # 只有当联合鞍点的值小于0时，才能求出联合概率
  # pval3 = adaptIntegrate(f=f,c(0,0),c(r,s))$integral
  # cat("\npval3=",pval3)
  pval = c(pval1,pval2)
  return(c(pval,flag=2,zeta_t,zeta_u,zeta_t0,zeta_u0))
}

GetProb_SPA2 = function(re, G1N1, G1N0, G2N1, G2N0, pos1, N0, r, s, G1, G2, g)
{
  # N0=N0
  # g=g
  flag = 0
  # sd1 = sqrt(var(S1))
  # sd2 = sqrt(var(S2))
  # mu1 = mean(S1)
  # mu2 = mean(S2)
  
  # zeta_t1 = getroot_t2(G1N1=G1N1, G1N0=G1N0, pos1=pos1, N0=N0, re=re, r=r)
  res.roott = getroot_t0(G1N1=G1N1, G1N0=G1N0, pos1=pos1, N0=N0, re=re, r=r)
  if(res.roott$flag.root == 1){
    zeta_t0 = res.roott$zeta
  }else{
    zeta_t0 = getroot_t1(init=0.1,G1N1=G1N1,G1N0=G1N0,G1=G1,re=re,pos1=pos1,N0=N0)$root
  }
  
  res.rootu = getroot_u0(G2N1=G2N1,G2N0=G2N0,pos1=pos1,N0=N0,re=re,s=s)
  if(res.rootu$flag.root == 1){
    zeta_u0 = res.rootu$zeta
  }else{
    zeta_u0 = NaN
    return(c(0,0,0,zeta_t=NaN,zeta_u=NaN,zeta_t0,zeta_u0))
  }
  # 用理论公式求解鞍点
  # init = getroot_K11(G1N1,G1N0,G2N1,G2N0,pos1,N0,re,ro,sd1,sd2,r,s)
  # cat("\n单独的鞍点为：", zeta_u0, zeta_t0, "\n")
  # cat("\n理论公式求得的鞍点为：", init,"\n")
  # 牛顿迭代求解鞍点
  # 检查迭代时的初始值是否满足条件M0>0
  
  # M01 = joint_CGF(t=zeta_t0,u=zeta_u0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag)$M0
  # M01 = mean(M01)
  # if(M01 != 0){
  #   out = getroot_K2(zeta1=0,zeta2=0,G1N1=G1N1,G2N1=G2N1,G1N0=G1N0,G2N0=G2N0,
  #                    re=re,N0=N0,pos1=pos1,r=r,s=s)
  # }
  
  mresid1 <- mresidR1(G1N1=G1N1,re,N0,pos1,flag)
  mresid2 <- mresidS2(G2N1=G2N1,re,N0,pos1,flag)
  mresid <- cbind(mresid1,mresid2)
  rm(list=c('mresid1','mresid2'))
  out = getroot_K2(zeta1=0,zeta2=0,G1N1=G1N1,G2N1=G2N1,G1N0=G1N0,G2N0=G2N0,
                   re=re,N0=N0,pos1=pos1,r=r,s=s,mresid=mresid)
  if(out$flag.para == 1)
  {
    # cat("\n联合鞍点第一组为：",out$root)
    zeta_t = out$root[1]
    zeta_u = out$root[2]
  }
  else{
    # zeta_t = zeta_t1
    zeta_t = zeta_t0
    zeta_u = zeta_u0
  }
  
  # out12 = getroot_K12(zeta1=init[1],zeta2=init[2],G1N1=G1N1,G2N1=G2N1,re=re,r=r,s=s)
  # zeta_t = out12[1]
  # zeta_u = out12[2]
  
  kt0 = sum(joint_CGF(t=zeta_t0,u=0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag,mresid=mresid)$K0)
  kt1 = sum(joint_CGF(t=zeta_t,u=0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag,mresid=mresid)$K0)
  
  ku0 = sum(joint_CGF(t=0,u=zeta_u0,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag,mresid=mresid)$K0)
  ku1 = sum(joint_CGF(t=0,u=zeta_u,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag,mresid=mresid)$K0)
  
  joint.CGFs = joint_CGF(t=zeta_t,u=zeta_u,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag,mresid=mresid)
  rm('mresid')
  
  k = sum(joint.CGFs$K0)
  ktt = sum(joint.CGFs$Ktt)
  ktu = sum(joint.CGFs$Ktu)
  kuu = sum(joint.CGFs$Kuu)
  kut = sum(joint.CGFs$Kut)
  k2 = rbind(cbind(ktt,ktu), cbind(kut,kuu))
  rm('joint.CGFs')
  
  #用鞍点分布函数逼近法，也即直接近似得到s,r的分布直接求出概率
  # P(S<s,R<r)=I11=fai(s1,y1,ro)
  # Re():取复数的实部，Im():取复数的虚部，i为虚部单位
  
  wu0_t = sqrt(abs(-2*(k - ku1 - zeta_t*r)))
  wu0_u = sqrt(abs(-2*(k - kt1 - zeta_u*s)))
  # 用较小的那个，大的会使得误差更显著
  if(wu0_t > wu0_u){
    wu0 = wu0_u
    s1 = sign(zeta_u0)*sqrt(abs(-2*(ku0 - zeta_u0*s)))
    v0 = sign(zeta_t)*sqrt(abs(-2*(k - ku0 - (zeta_u-zeta_u0)*s - zeta_t*r)))
    G = 1/sqrt(abs(kuu - ktu^2/ktt))
    b = (wu0-s1)/v0
    r1 = (v0 - b*s1) / sqrt(1+b^2)
    ro1 = -b / sqrt(1+b^2)
    
    I11 = pnorm2d(r1, s1, rho=ro1)
    I12 = pnorm2d(wu0)*dnorm2d(v0)*(1/v0 - 1/(zeta_t*G))
    I21 = pnorm2d(v0)*dnorm2d(s1)*(1/wu0 - 1/(zeta_u * sqrt(abs(ktt))))
    # I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_u * 1/sqrt(kuu)) * (1/v0 - 1/(zeta_t*G)) / (2*pi)
    I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_u * 1/sqrt(abs(ktt))) * (1/v0 - 1/(zeta_t*G)) / (2*pi)
    h2 = -(s1^2 - 2* ro1*s1*r1 + r1^2) / (2*(1-ro1^2))
    I31 = exp(h2)* (1/wu0 - 1/zeta_u * 1/sqrt(abs(ktt))) * (1/v0 - 1/(zeta_t*G)) / (2*pi)
  }else{
    wu0 = wu0_t
    r1 = sign(zeta_t0)*sqrt(abs(-2*(kt0 - zeta_t0*r)))
    v0 = sign(zeta_u)*sqrt(abs(-2*(k - kt0 - (zeta_t-zeta_t0)*r - zeta_u*s)))
    G = 1/sqrt(abs(kuu - ktu^2/ktt))
    b = (wu0-r1)/v0
    s1 = (v0 - b*r1) / sqrt(1+b^2)
    ro1 = -b / sqrt(1+b^2)
    
    I11 = pnorm2d(r1, s1, rho=ro1)
    I12 = pnorm2d(wu0)*dnorm2d(v0)*(1/v0 - 1/(zeta_u*G))
    I21 = pnorm2d(v0)*dnorm2d(r1)*(1/wu0 - 1/(zeta_t * sqrt(ktt)))
    # I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_t * 1/sqrt(ktt)) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
    I22 = exp(k - zeta_t*r - zeta_u*s) * (1/wu0 - 1/zeta_t * 1/sqrt(ktt)) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
    h2 = -(r1^2 - 2* ro1*r1*s1 + s1^2) / (2*(1-ro1^2))
    I31 = exp(h2)* (1/wu0 - 1/zeta_t * 1/sqrt(ktt)) * (1/v0 - 1/(zeta_u*G)) / (2*pi)
  }
  # wu0 = sign(zeta_t)*sqrt(abs(-2*(k - ku1 - zeta_t*r)))
  # cat("\nr1=",r1,",s1=",s1,",v0=",v0,",wu0=",wu0,",ro1=",ro1,",b=",b,",G=",G)
  
  # 用fMultivar包，计算二维正态分布的概率密度函数和分布函数
  # pnorm2d：求的是标准二元正态分布函数的概率值均值为0，方差为1
  # ktt_e = sqrt(abs(ktt))
  # cat("\nI11=",I11,",I12=",I12,",I21=",I21,",I22=",I22,",I31=",I31)
  pval1 = I11+I12+I21+I22
  pval2 = I11+I12+I21+I31
  if(is.na(pval1)){
    pval1 = 0
  }
  if(!is.na(pval1) & !is.na(pval2)){
    if(pval1 > 1) pval1 = 1
    if(pval1 < 0) pval1 = 0
    if(pval2 > 1) pval2 = 1
    if(pval2 < 0) pval2 = 0
  }
  # if((pval1 >= 0.05) & (pval1 <= 1))pval1=pval1
  # if((pval2 >= 0.05) & (pval2 <= 1))pval1=pval2
  
  # pval3 = pnorm2d(r/sd1, s/sd2, rho=ro) # 相当于采用了正态方法去计算联合概率值
  cat("\npval1=",pval1,",pval2=",pval2,"\n")
  
  #用鞍点密度积分法计算p值(由s,r的概率密度函数积分求得概率值)
  # f = function(x){
  #   1/(2*pi*sqrt(det(k2))) * exp(k - zeta_t*x[1] - zeta_u*x[2])
  # }
  # # 只有当联合鞍点的值小于0时，才能求出联合概率
  # pval3 = adaptIntegrate(f=f,c(0,0),c(r,s))$integral
  # cat("\npval3=",pval3)
  pval = c(pval1,pval2)
  return(c(pval,out$flag.para,zeta_t,zeta_u,zeta_t0,zeta_u0))
}

# 求两个鞍点的值
# 1.用牛顿迭代的方法计算鞍点，cox和logit单独的鞍点作为初始值
getroot_K1 = function(zeta1, zeta2, G1N1, G1N0, G2N1, G2N0, re, N0, pos1, r, s, mresid)
{
  # zeta1 <- zeta_t0
  # zeta2 <- zeta_u0
  itmax=500
  epsilon=0.001
  iter.sum = 0
  flag.para <- 0
  Diverge <- 0
  
  # 采用牛顿迭代的方法求解，初始值用各自单独的鞍点
  for(iter in 1:itmax)# 内循环：估计betaR,gammaR,betaS,gammaS,u,v
  {
    iter.sum <- iter.sum+1
    CGFs0 = joint_CGF0(t=zeta1,u=zeta2,G1N1=G1N0,G2N1=G2N0,re=re,N0=N0,pos1=pos1,flag=0)
    CGFs1 = joint_CGF(t=zeta1,u=zeta2,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,pos1=pos1,flag=1,mresid=mresid)
    Ktt <- sum(CGFs1$Ktt) + N0*(CGFs0$Ktt)
    Ktu <- sum(CGFs1$Ktu) + N0*(CGFs0$Ktu)
    
    Kut <- sum(CGFs1$Kut) + N0*(CGFs0$Kut)
    Kuu <- sum(CGFs1$Kuu) + N0*(CGFs0$Kuu)
    
    Info1 <- Ktt+Kut
    Info2 <- Ktu+Kuu
    # cat("\nKtt=",Ktt, " Ktu=",Ktu," Kut=",Kut, " Kuu=",Kuu,"\n")
    Info0 <- rbind(Info1, Info2)
    
    # 原函数
    f1 <- sum(CGFs1$Kt) + N0*(CGFs0$Kt) - r
    f2 <- sum(CGFs1$Ku) + N0*(CGFs0$Ku) - s
    f <- rbind(f1,f2)
    # f <- f1+f2
    # cat("f=",f,"\n")
    para0 <- rbind(zeta1,zeta2)
    para1 <- para0 - f/Info0
    # cat("para1=",para1,"\n")
    
    rm(list=c('CGFs1','CGFs0'))
    if (is.na(para1[1]) | is.na(para1[2])) {
      Diverge = 1
      break
    }
    if(max(abs(para1 - para0)) > 100)
    {
      # cat("\nBreak, can't convergent!\n")
      Diverge <- 1
      break
    }
    if(max(abs(para1 - para0)) < epsilon)
    {
      # cat("max(abs(para1 - para0))=",max(abs(para1 - para0))," iter=",iter.sum,"\n")
      # cat("\nout= ",para1)
      flag.para <- 1
      break
    }
    if (iter.sum == itmax){
      # cat("\nCan't convergent!")
      Diverge <- 1
      break
    }
    #重置参数的值
    zeta1 <- para1[1]
    zeta2 <- para1[2]
  }
  return(list(root = para1, iter.sum = iter.sum, flag.para = flag.para))
}
getroot_K2 = function(zeta1, zeta2, G1N1, G1N0, G2N1, G2N0, re, N0, pos1, r, s,mresid)
{
  # zeta1 <- 0
  # zeta2 <- 0
  itmax=500
  epsilon=0.001
  iter.sum = 0
  flag.para <- 0
  Diverge <- 0
  
  # 采用牛顿迭代的方法求解，初始值用各自单独的鞍点
  for(iter in 1:itmax)
  {
    iter.sum <- iter.sum + 1
    CGFs = joint_CGF(t=zeta1,u=zeta2,G1N1=G1N1,G2N1=G2N1,re=re,N0=N0,
                     pos1=pos1,flag=0,mresid=mresid)
    Ktt <- sum(CGFs$Ktt)
    Ktu <- sum(CGFs$Ktu)
    
    Kut <- sum(CGFs$Kut)
    Kuu <- sum(CGFs$Kuu)
    
    Info1 <- Ktt+Kut
    Info2 <- Ktu+Kuu
    # cat("\nKtt=",Ktt, " Ktu=",Ktu," Kut=",Kut, " Kuu=",Kuu,"\n")
    Info0 <- rbind(Info1, Info2)
    
    # 原函数
    f1 <- sum(CGFs$Kt) - r
    f2 <- sum(CGFs$Ku) - s
    f <- rbind(f1,f2)
    # f <- f1+f2
    # cat("f=",f,"\n")
    para0 <- rbind(zeta1,zeta2)
    para1 <- para0 - f/Info0
    # cat("para0=",para0,"  para1=",para1,"\n")
    
    rm('CGFs')
    if (is.na(para1[1]) | is.na(para1[2])) {
      Diverge = 1
      break
    }
    if(max(abs(para1 - para0)) < epsilon)
    {
      # cat("max(abs(para1 - para0))=",max(abs(para1 - para0))," iter=",iter.sum,"\n")
      # cat("\nout= ",para1)
      flag.para <- 1
      break
    }
    if(max(abs(para1 - para0)) > 100)
    {
      # cat("\nBreak, can't convergent!\n")
      Diverge <- 1
      break
    }
    if (iter.sum == itmax){
      # cat("\nCan't convergent!")
      Diverge <- 1
      break
    }
    #重置参数的值
    
    zeta1 <- para1[1]
    zeta2 <- para1[2]
  }
  return(list(root = para1, iter.sum = iter.sum, flag.para = flag.para))
}

# 没有加入基因型的MGF和CGF
joint_CGF = function(t, u, G1N1, G2N1, re, N0, pos1, flag, mresid)
{
  # t = zeta1
  # u = zeta2
  # G1N1 = G1N1
  # G2N1 = G2N1
  # flag = 0
  n = length(G1N1)
  if (N0!=0 & flag == 1){
    mu = re$mu[pos1]
    mresidS = re$mresidS[pos1]
    mresidR = re$mresidR[pos1]
  }
  if(N0 == 0)
  {
    mresidR = re$mresidR
    mresidS = re$mresidS
  }
  
  cumul <- NULL
  e_mresid <- as.matrix(exp((mresid[,1:n]*t) + (mresid[,(n+1):(2*n)]*u)))
  
  M0 = apply(e_mresid, 2, mean)
  Mt = G1N1 * apply(mresidR*e_mresid, 2, mean)
  Mu = G2N1 * apply(mresidS*e_mresid, 2, mean)
  Mtt = G1N1^2 * apply(mresidR^2 * e_mresid, 2, mean)
  Mtu = G1N1*G2N1 * apply(mresidR*mresidS* e_mresid, 2, mean)
  Muu = G2N1^2 * apply(mresidS^2 * e_mresid, 2, mean)
  Mut = G1N1*G2N1 * apply(mresidR*mresidS* e_mresid, 2, mean)
  rm('e_mresid')
  
  K0 = log(M0)
  Kt = Mt / M0
  Ku = Mu / M0
  
  Ktt = (Mtt * M0 - Mt^2) / M0^2
  Ktu = (Mtu * M0 - Mt*Mu) / M0^2
  Kuu = (Muu * M0 - Mu^2) / M0^2
  Kut = (Mut * M0 - Mu*Mt) / M0^2
  cumul = data.frame(K0=K0, Kt=Kt, Ku=Ku, 
                     Ktt=Ktt, Ktu=Ktu, 
                     Kuu=Kuu, Kut=Kut, 
                     M0=M0)
  
  # for(i in 1:n){
  #   # 加入了基因型后的统计量的CGFs和MGFs
  #   # e_mresidR = exp(mresidR * t[i])
  #   # e_mresidS = exp(mresidS * u[i])
  #   e_mresid = exp(mresidR*t[i] + mresidS*u[i])
  #   M0 = mean(e_mresid)
  #   
  #   Mt = G1N1[i] * mean(mresidR*e_mresid)
  #   Mu = G2N1[i] * mean(mresidS*e_mresid)
  #   Mtt = G1N1[i]^2 * mean(mresidR^2 * e_mresid)
  #   Mtu = G1N1[i]*G2N1[i] * mean(mresidR*mresidS* e_mresid)
  #   Muu = G2N1[i]^2 * mean(mresidS^2 * e_mresid)
  #   Mut = G1N1[i]*G2N1[i] * mean(mresidR*mresidS* e_mresid)
  #   
  #   K0 = log(M0)
  #   Kt = Mt / M0
  #   Ku = Mu / M0
  #   
  #   Ktt = (Mtt * M0 - Mt^2) / M0^2
  #   Ktu = (Mtu * M0 - Mt*Mu) / M0^2
  #   Kuu = (Muu * M0 - Mu^2) / M0^2
  #   Kut = (Mut * M0 - Mu*Mt) / M0^2
  #   
  #   cumul = rbind(cumul, c(K0, Kt, Ku, Ktt, Ktu, Kuu, Kut, M0))
  # }
  # cumul = data.frame(K0=cumul[,1], Kt=cumul[,2], Ku=cumul[,3], 
  #                    Ktt=cumul[,4], Ktu=cumul[,5], 
  #                    Kuu=cumul[,6], Kut=cumul[,7], 
  #                    M0=cumul[,8])
  return(cumul)
}

joint_CGF0 = function(t, u, G1N1, G2N1, re, N0, pos1, flag)
{
  # t = init[1]
  # u = init[2]
  # G1N1 = G1N0
  # G2N1 = G2N0
  # flag = 0
  t = t*G1N1
  u = u*G2N1
  if(flag == 0){
    mu = re$mu[!(seq_along(re$mu) %in% pos1)]
    mresidS = re$mresidS[!(seq_along(re$mresidS) %in% pos1)]
    mresidR = re$mresidR[!(seq_along(re$mresidR) %in% pos1)]
  }
  n = length(G1N1)
  cumul <- NULL
  e_mresid <- NULL
  for(i in 1:n){
    e_mresid = cbind(exp(mresidR*t[i] + mresidS*u[i]))
  }
  M0 = apply(e_mresid, 2, mean)
  Mt = G1N1 * apply(mresidR*e_mresid, 2, mean)
  Mu = G2N1 * apply(mresidS*e_mresid, 2, mean)
  Mtt = G1N1^2 * apply(mresidR^2 * e_mresid, 2, mean)
  Mtu = G1N1*G2N1 * apply(mresidR*mresidS* e_mresid, 2, mean)
  Muu = G2N1^2 * apply(mresidS^2 * e_mresid, 2, mean)
  Mut = G1N1*G2N1 * apply(mresidR*mresidS* e_mresid, 2, mean)
  rm('e_mresid')
  
  K0 = log(M0)
  Kt = Mt / M0
  Ku = Mu / M0
  
  Ktt = (Mtt * M0 - Mt^2) / M0^2
  Ktu = (Mtu * M0 - Mt*Mu) / M0^2
  Kuu = (Muu * M0 - Mu^2) / M0^2
  Kut = (Mut * M0 - Mu*Mt) / M0^2
  cumul = data.frame(K0=K0, Kt=Kt, Ku=Ku, 
                     Ktt=Ktt, Ktu=Ktu, 
                     Kuu=Kuu, Kut=Kut, 
                     M0=M0)
  
  # for(i in 1:n){
  #   # 加入了基因型后的统计量的CGFs和MGFs
  #   # e_mresidR = exp(mresidR * t[i])
  #   # e_mresidS = exp(mresidS * u[i])
  #   e_mresid = exp(mresidR*t[i] + mresidS*u[i])
  #   M0 = mean(e_mresid)
  #   
  #   Mt = G1N1[i] * mean(mresidR*e_mresid)
  #   Mu = G2N1[i] * mean(mresidS*e_mresid)
  #   Mtt = G1N1[i]^2 * mean(mresidR^2 * e_mresid)
  #   Mtu = G1N1[i]*G2N1[i] * mean(mresidR*mresidS* e_mresid)
  #   Muu = G2N1[i]^2 * mean(mresidS^2 * e_mresid)
  #   Mut = G1N1[i]*G2N1[i] * mean(mresidR*mresidS* e_mresid)
  #   
  #   K0 = log(M0)
  #   Kt = Mt / M0
  #   Ku = Mu / M0
  #   
  #   Ktt = (Mtt * M0 - Mt^2) / M0^2
  #   Ktu = (Mtu * M0 - Mt*Mu) / M0^2
  #   Kuu = (Muu * M0 - Mu^2) / M0^2
  #   Kut = (Mut * M0 - Mu*Mt) / M0^2
  #   
  #   cumul = rbind(cumul, c(K0, Kt, Ku, Ktt, Ktu, Kuu, Kut, M0))
  # }
  # cumul = data.frame(K0=cumul[,1], Kt=cumul[,2], Ku=cumul[,3], 
  #                    Ktt=cumul[,4], Ktu=cumul[,5], 
  #                    Kuu=cumul[,6], Kut=cumul[,7], 
  #                    M0=cumul[,8])
  return(cumul)
}

mresidR1 = function(G1N1, re, N0, pos1, flag)
{
  n = length(G1N1)
  if(N0 == 0)
  {
    mresidR = re$mresidR
    mresidS = re$mresidS
  }
  if(N0!= 0 & flag == 1){
    mu = re$mu[pos1]
    mresidS = re$mresidS[pos1]
    mresidR = re$mresidR[pos1]
  }
  mresidR = as.matrix(mresidR,ncol=1)
  G1N1 = as.matrix(G1N1,nrow=1)
  mresid = mresidR %*% t(G1N1)
  return(mresid)
}
mresidS2 = function(G2N1, re, N0, pos1, flag)
{
  n = length(G2N1)
  blocksize = n/10
  
  if(N0 == 0)
  {
    mresidR = re$mresidR
    mresidS = re$mresidS
  }
  if(N0!= 0 & flag == 1){
    mu = re$mu[pos1]
    mresidS = re$mresidS[pos1]
    mresidR = re$mresidR[pos1]
  }
  mresidS = as.matrix(mresidS,ncol=1)
  # mresidS = split(mresidS, rep(1:10, each=blocksize))
  G2N1 = as.matrix(G2N1,nrow=1)
  # G2N1 = split(G2N1, rep(1:10), each=blocksize)
  mresid <- mresidS %*% t(G2N1)
  return(mresid)
}

# 求单独logit的鞍点：两种方法
getroot_t1 = function(init, re, G1N1, G1N0, G1, pos1, N0, tol = .Machine$double.eps^0.25, maxiter = 1000)
{
  q <- sum(G1 * re$data$R)
  mu = re$mu
  g.pos <- sum(G1N1[which(G1N1 > 0)])
  g.neg <- sum(G1N1[which(G1N1 < 0)])
  if (q >= g.pos || q <= g.neg) {
    return(list(root = Inf, n.iter = 0, Is.converge = TRUE))
  }
  else {
    t <- init
    K1_eval <- Kt1(t, G1N1, G1N0, pos1, N0, re, q)
    prevJump <- Inf
    rep <- 1
    repeat {
      K2_eval <- Kt2(t, G1N1, G1N0, pos1, N0, re)
      tnew <- t - K1_eval/K2_eval
      if (is.na(tnew)) {
        conv = FALSE
        break
      }
      if (abs(tnew - t) < tol) {
        conv <- TRUE
        break
      }
      if (rep == maxiter) {
        conv <- FALSE
        break
      }
      newK1 <- Kt1(tnew, G1N1, G1N0, pos1, N0, re, q)
      if (sign(K1_eval) != sign(newK1)) {
        if (abs(tnew - t) > prevJump - tol) {
          tnew <- t + sign(newK1 - K1_eval) * prevJump/2
          newK1 <- Kt1(tnew, G1N1, G1N0, pos1, N0, re, q)
          prevJump <- prevJump/2
        }
        else {
          prevJump <- abs(tnew - t)
        }
      }
      rep <- rep + 1
      t <- tnew
      K1_eval <- newK1
    }
    return(list(root = t, n.iter = rep, Is.converge = conv))
  }
}
getroot_t0 = function(t, G1N1, G1N0, pos1, N0, re, r)
{
  out <- tryCatch({
    result <- uniroot(Kt1, c(-50, 50), extendInt = "upX",
                      G1N1 = G1N1, G1N0 = G1N0, pos1 = pos1,
                      N0 = N0, re = re, r = r)
    zeta <- ifelse(is.na(result$root), 0, result$root)
    flag.root <- 1
    list(zeta = zeta, flag.root = flag.root)
  }, error = function(e) {
    # Code to execute when no root is found or an error occurs
    zeta <- 0
    flag.root <- 0  # Custom flag value
    list(zeta = zeta, flag.root = flag.root)
  })
  return(out)
}
# 原始方法
getroot_t2 = function(t, G1N1, G1N0, pos1, N0, re, r)
{
  out <- uniroot(Kt1, c(-50,50), extendInt = "upX",
            G1N1=G1N1, G1N0=G1N0, pos1=pos1, 
            N0=N0, re=re, r=r)
  zeta = out$root
  return(zeta)
}
# 求单独cox的鞍点，t=0
getroot_u0 = function(re, G2N1, G2N0, pos1, N0, s)
{
  out <- tryCatch({
    result <-uniroot(Ku1, c(-50,50), extendInt = "upX",
                G2N0=G2N0, G2N1=G2N1, pos1=pos1,
                N0=N0, s=s, re=re)
    zeta <- ifelse(is.na(result$root), 0, result$root)
    flag.root <- 1
    list(zeta = zeta, flag.root = flag.root)
  }, error = function(e) {
    # Code to execute when no root is found or an error occurs
    zeta <- 0
    flag.root <- 0  # Custom flag value
    list(zeta = zeta, flag.root = flag.root)
  })
  return(out)
}
# cox、logit单独的累积生成函数
Kt0 = function(t, G1N1, G1N0, pos1, N0, re)
{
  n.t = length(t)
  Kt0 = rep(0,n.t)
  if(N0 != 0){
    mu0 = re$mu[!seq_along(re$mu) %in% pos1]
  }else{
    mu0 = 0
  }
  for(i in 1:n.t){
    t1 = t[i]
    temp01 = log(1-mu0+mu0*exp(G1N0*t1))
    temp02 = G1N0*mu0
    temp11 = log(1-re$mu[pos1]+re$mu[pos1]*exp(G1N1*t1))
    temp12 = G1N1*re$mu[pos1]
    Kt0[i] = sum(temp01) - t1*sum(temp02) + sum(temp11) - t1*sum(temp12)
  }
  return(Kt0)
}
Ku0 = function(u, G2N1, G2N0, pos1, N0, re)
{
  n.u = length(u)
  Ku0 = rep(0,n.u)
  for(i in 1:n.u){
    u1 = u[i]
    t2N1 = u1*G2N1
    t2N0 = u1*G2N0
    Ku0[i] = N0*re$Ku0_emp(t2N0) + sum(re$Ku0_emp(t2N1))
  }
  return(Ku0)
}
# logit的累积生成函数一阶导
Kt1 = function(t, G1N1, G1N0, pos1, N0, re, r)
{
  n.t = length(t)
  Kt1 = rep(0,n.t)
  if(N0 != 0){
    mu0 = re$mu[!seq_along(re$mu) %in% pos1]
  }else{
    mu0 = 0
  }
  for(i in 1:n.t){
    t1 = t[i]
    temp01 <- (1 - mu0) * exp(-G1N0 * t1) + mu0
    temp02 <- mu0 * G1N0
    temp11 <- (1 - re$mu[pos1]) * exp(-G1N1 * t1) + re$mu[pos1]
    temp12 <- re$mu[pos1] * G1N1
    Kt1[i] <- sum(temp02/temp01) + sum(temp12/temp11) - r
  }
  return(Kt1)
}
# cox的累积生成函数一阶导
Ku1 = function(u, G2N1, G2N0, pos1, N0, re, s)
{
  n.u = length(u)
  Ku1 = rep(0,n.u)
  for(i in 1:n.u){
    u1 = u[i]
    t2N1 = u1*G2N1
    t2N0 = u1*G2N0
    Ku1[i] = N0*G2N0*re$Ku1_emp(t2N0) + sum(G2N1*re$Ku1_emp(t2N1)) - s
  }
  return(Ku1)
}
Kt2 = function(t, G1N1, G1N0, pos1, N0, re)
{
  n.t <- length(t)
  Kt2 <- rep(0, n.t)
  if(N0 != 0){
    mu0 = re$mu[!seq_along(re$mu) %in% pos1]
  }else{
    mu0 = 0
  }
  for (i in 1:n.t) {
    t1 <- t[i]
    temp01 <- (1 - mu0) * exp(-G1N0 * t1) + mu0^2
    temp02 <- (1 - mu0) * mu0 * G1N0^2 * exp(-G1N0 * t1)
    temp11 <- ((1 - re$mu[pos1]) * exp(-G1N1 * t1) + re$mu[pos1])^2
    temp12 <- (1 - re$mu[pos1]) * re$mu[pos1] * G1N1^2 * exp(-G1N1 * t1)
    
    Kt2[i] <- sum(temp02/temp01) + sum(temp12/temp11)
  }
  return(Kt2)
}


check_input = function(pIDs, gIDs, range)
{
  if(is.null(pIDs) & is.null(gIDs))
    stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  
  pIDs = as.character(pIDs)
  gIDs = as.character(gIDs)
  
  if(any(!is.element(pIDs, gIDs)))
    stop("All elements in pIDs should be also in gIDs.")
  
  if(anyDuplicated(gIDs)!=0)
    stop("Argument 'gIDs' should not have a duplicated element.")
  
  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")
  
  # if(length(mresid)!=length(pIDs))
  #   stop("length(mresid)!=length(pIDs) where mresid is the martingale residuals from coxph() in survival package.")
  
  p2g = NULL
  if(length(pIDs)!=length(gIDs)){
    p2g = match(pIDs, gIDs)
  }else{
    if(any(pIDs != gIDs))
      p2g = match(pIDs, gIDs)
  }
  return(list(p2g=p2g,pIDs=pIDs))
}