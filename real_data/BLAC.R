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
pd <- pd[!grepl("TCGA-GV-A3QG-01|TCGA-ZF-AA4W-01",rownames(pd)),]

N <- length(pd$primary_therapy_outcome_success)
pd$primary_therapy_outcome_success <- ifelse((pd$primary_therapy_outcome_success==""),0,pd$primary_therapy_outcome_success)

row_non0 <- NULL
for(i in 1:N){
  if(pd$primary_therapy_outcome_success[i] != "0")
    row_non0 = rbind(row_non0, rownames(pd)[i])
}
pd <- pd[row_non0,]

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
identical(colnames(betaData),rownames(pd))

myLoad <- champ.filter(beta = betaData, pd = pd)
beta.values <- myLoad$beta

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
