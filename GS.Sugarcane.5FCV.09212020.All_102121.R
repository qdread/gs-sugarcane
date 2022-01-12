# Paths modified by QDR, 05 Nov 2021
# Install version 4.1.0 of sommer, current as of 21 Sep 2020
# install.packages('~/R/sommer_4.1.0.tar.gz', repos = NULL, type = 'source')

date="09212020" #08042020
data_dir = 'C:/Users/qdread/onedrive_usda/ars_projects/islam/data'
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer)
source("Sugarcane.GS.RRBLUP.ADE.RKHS.SVR.RF.Algorithm.09212020.R")
source("geno_to_allelecnt.R")
source("CoincidenceIndex.R")

GD <- read.delim(file.path(data_dir, "Numerical.MAF2p.432.Individuals.txt"), header=T, stringsAsFactors = F)
GD[1:6,1:6]

dat=GD
nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)

rownames(dat)=nan
dat[1:6,1:6]
dat[dat==1] <- 2
dat[dat==0.5] <- 1
dat[1:6,1:6]

df.nam.geno <- data.frame(rownames(dat))
names(df.nam.geno) <- "Taxa"

df2 <- as.matrix(dat)

df2[1:6,1:6]
any(is.na(df2))
sum(is.na(df2))


phen.taxa <- read.csv(file.path(data_dir, "Phenotypes.Rust.GS.09212020.csv"), header=T)
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
