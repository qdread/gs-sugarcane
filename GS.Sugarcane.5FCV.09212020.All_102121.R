rm(list=ls())
setwd("/homes/omo/public_html/Collaborations/Sugarcane/new_vcf_august/")
date="09212020" #08042020
#setwd("")
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer)
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sugarcane/new_vcf_august/Sugarcane.GS.RRBLUP.ADE.RKHS.SVR.RF.Algorithm.09212020.R")
source("https://raw.githubusercontent.com/ekfchan/evachan.org-Rscripts/master/rscripts/geno_to_allelecnt.R")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")

GD <- read.delim("http://people.beocat.ksu.edu/~omo/Collaborations/Sugarcane/new_vcf_august/Numerical.MAF2p.432.Individuals.txt", header=T, stringsAsFactors = F)
GD[1:6,1:6]
#dat=myGAPIT$GD
dat=GD
nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
#dat=dat-1
rownames(dat)=nan
dat[1:6,1:6]
dat[dat==1] <- 2
dat[dat==0.5] <- 1
dat[1:6,1:6]

df.nam.geno <- data.frame(rownames(dat))
names(df.nam.geno) <- "Taxa"

df2 <- as.matrix(dat)
#df2-1
df2[1:6,1:6]
any(is.na(df2))
sum(is.na(df2))


phen.taxa <- read.csv("http://people.beocat.ksu.edu/~omo/Collaborations/Sugarcane/new_vcf_august/Phenotypes.Rust.GS.09212020.csv", header=T)
head(phen.taxa)

colnames(phen.taxa)[1] <- "Taxa"
phen.taxa$Taxa <- as.character(phen.taxa$Taxa)
head(phen.taxa)

com.taxa <- merge(phen.taxa, df.nam.geno)
dim(com.taxa)
head(com.taxa)
com.taxa$Taxa <- as.character(com.taxa$Taxa)

com.taxa.df <- data.frame(com.taxa$Taxa)
colnames(com.taxa.df) <- "Taxa"

phen<-com.taxa



polyRAD.Geno <- dat[match(com.taxa$Taxa, rownames(dat)),]

head(com.taxa)
names(com.taxa)
com.taxa2 <- com.taxa[,c(1,2,4,6,8,10,11)]# -c(2,4,6,8,10,12,14,16,18,20)
head(com.taxa2)
#phenames <- as.character(colnames(com.taxa2[,-c(1:4)]))
phenames <- as.character(colnames(com.taxa2[,-1]))
phenames <- phenames
k = 5
cycles=25

number.of.folds=k

species <- "Sugarcane" #either SC or Msa
# m_train1



diploid.polyploid.5cfv(phenames=phenames, m_train.pheno=com.taxa2, diploid.Geno=polyRAD.Geno, species=species, k=k, cycles = cycles,
                                  number.of.folds=number.of.folds, date=date)



phenames=phenames; m_train.pheno=com.taxa2; diploid.Geno=polyRAD.Geno; species=species; k=k; cycles = cycles; number.of.folds=number.of.folds; date=date
