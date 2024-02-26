# getOption("timeout")
# options(timeout = 60)

# BiocManager::install("SummarizedExperiment",force = TRUE)
# BiocManager::install("ChAMP",force = TRUE)

# library(GenomeInfoDbData)
# library(GenomeInfoDb)
# library(GenomicRanges)
# library(DelayedArray)
# library(SummarizedExperiment)
# library(minfi)
# BiocManager::install("geneLenDataBase")
# install.packages("D:/broswer/geneLenDataBase_1.34.0.tar.gz", repos = NULL, type = "source")


# # 计算无进展生存期(PFS)：
# # 先计算进展时间TTP：
# # 1.如果患者死亡，TTP="days to death",pd[,4]
# # 2.如果患者未死亡但发生了最后随访，TTP="days to last follow up",pd[,5]
# # 3.如果患者既未死亡也未发生最后随访,但有新肿瘤时间不为空，
# # TTP="days_to_new_tumor_event_after_initial_treatment",pd[,9]
# # 4.三者都无，无法计算TTP
# # 然后根据ORR计算PFS
# # 1.如果肿瘤为部分响应/全响应，PFS=TTP-pd[,10]
# n <- length(pd$sampleID)
# # 先计算TTP
# for(i in 1:n){
#   if(pd$event[i] == 1) {pd$TTP[i] = pd$days_to_death[i]}
#   else if(pd$event[i] == 0 & pd$lost_follow_up[i] == "NO"){
#     pd$TTP[i] = pd$days_to_last_followup[i]
#   }
#   else if(pd$event[i] == 0 & pd$lost_follow_up[i] == "YES" & (pd$days_to_new_tumor_event_after_initial_treatment[i] != "NA")){
#     pd$TTP[i] = days_to_new_tumor_event_after_initial_treatment[i]
#   }
#   else {pd$TTP[i] = 0}
# }
# # 计算PFS

# # 将无法确定TTE数据的样本删除
# time_non0 <- NULL
# for(i in 1:n){
#   if(pd$time[i] != 0){
#     time_non0 = rbind(time_non0, rownames(pd)[i])
#   }
# }
# pd <- pd[time_non0,]
# pd["TCGA-GD-A76B-01",11] <- pd["TCGA-GD-A76B-01",5]


setwd("F:/real data/TCGA-BLCA")
library(BiocManager)
library(data.table)
library(dplyr)
library(tibble)
library(xlsx)

library(ChAMP)

pd.all <- read.delim("TCGA.BLCA.sampleMap_BLCA_clinicalMatrix",header = T)
data <- read.delim("F:/real data/TCGA-BLCA/BLCA_survival.txt")


identical(rownames(data),rownames(pd.all))

pd.all$PFS_Status <- data$PFI
pd.all$PFS <- data$PFI.time
colnames(pd.all)
pd <- pd.all[,c("sampleID","vital_status","days_to_death","days_to_last_followup","lost_follow_up",
                "PFS_Status","PFS","primary_therapy_outcome_success","gender","weight")]

# pd$patient <- substr(pd$sampleID,1,12)


# 删除异常值
# "TCGA-GV-A3QG-01","TCGA-ZF-AA4W-01"的死亡时间和上次随访时间为[Discrepancy]
pd <- pd[!grepl("TCGA-GV-A3QG-01|TCGA-ZF-AA4W-01",rownames(pd)),]

# 处理肿瘤响应中的空值，直接删除
N <- length(pd$primary_therapy_outcome_success)
pd$primary_therapy_outcome_success <- ifelse((pd$primary_therapy_outcome_success==""),0,pd$primary_therapy_outcome_success)

row_non0 <- NULL
for(i in 1:N){
  if(pd$primary_therapy_outcome_success[i] != "0")
    row_non0 = rbind(row_non0, rownames(pd)[i])
}
pd <- pd[row_non0,]


# 处理缺失值，将响应值符号化，确定event、协变量X1,X2、ORR
pd$days_to_death <- ifelse((pd$vital_status=="LIVING"), 0, pd$days_to_death)
pd$days_to_last_followup <- ifelse((pd$days_to_last_followup==""),0,pd$days_to_last_followup)
pd$primary_therapy_outcome_success <- ifelse((pd$primary_therapy_outcome_success=="Progressive Disease" | 
                                                pd$primary_therapy_outcome_success=="Stable Disease"), 0, 1)
pd$event <- ifelse(pd$vital_status=='LIVING',0,1)
pd$X1 <- ifelse(pd$gender=='MALE',1,0)
pd$X2 <- pd$weight

pd <- na.omit(pd)
rownames(pd) <- pd$sampleID
N <- length(pd$sampleID)

a <- fread("HumanMethylation450",data.table = F)
a <- na.omit(a)
rownames(a) <- a[,1]
a <- a[,-1]
a <- a[,!grepl("TCGA-GV-A3QG-01|TCGA-ZF-AA4W-01",colnames(a))]
suit <- NULL
x <- match(rownames(pd),colnames(a))
for(i in 1:N){
  if(!is.na(x[i])){suit=rbind(suit,rownames(pd)[i])}
}
a <- a[,suit]
pd <- pd[suit,]
betaData <- as.matrix(a)

# pd <- pd[match(colnames(a),rownames(pd)),]
identical(colnames(betaData),rownames(pd))

myLoad <- champ.filter(beta = betaData, pd = pd)
beta.values <- myLoad$beta

# 处理甲基化信号值
if (!is.matrix(beta.values)) {
  beta.values <- as.matrix(beta.values)
}
Problems <- which(beta.values < 0 | beta.values > 1)
if (is.character(beta.values[1, 1])) {
  stop("beta.values contains characters, matrix/data.frame must be all numeric\n")
}
if (length(Problems) != 0) {
  beta.values[Problems] <- NA
}
onevalues <- which(beta.values == 1)
zerovalues <- which(beta.values == 0)
if (length(onevalues) > 0 | length(zerovalues) > 0) {
  if (length(onevalues) > 0) {
    beta.values[onevalues] <- NA
    beta.values[onevalues] <- max(beta.values, na.rm = T)
  }
  if (length(zerovalues) > 0) {
    beta.values[zerovalues] <- NA
    beta.values[zerovalues] <- min(beta.values, 
                                   na.rm = T)
  }
}
beta.values = log(beta.values/(1 - beta.values))




pd <- myLoad$pd
betaData <- beta.values

# require(openxlsx)
# write.xlsx(pd,file= "D:/broswer/TCGA-BLCA/pd.xlsx",sheetName = "pd")
save(betaData, file = "D:/broswer/TCGA-BLCA/betaData.RData")
save(pd, file = "D:/broswer/TCGA-BLCA/pd.RData")











