###Adapted from the data processing file accompanying:
###Lock, E. F., & Dunson, D. B. (2013). Bayesian consensus clustering. Bioinformatics (Oxford, England), 29(20), 2610–2616. http://doi.org/10.1093/bioinformatics/btt425


setwd("~/Dropbox/GwenDataIntegration/BCC")
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library(impute)
library("MCMCpack")

###Data can be downloaded from the TCGA Data Portal at: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/

GE      <- read.csv("BRCA.exp.348.med.csv", header = TRUE)
miRNA   <- read.csv("BRCA.348.precursor.txt", header = TRUE)
Protein <- read.table("rppaData-403Samp-171Ab-Trimmed.txt", header = TRUE)
Meth    <- read.table("BRCA.Methylation.574probes.802.txt", header = TRUE)

################Match columns (samples) between sources##############################
namesExp     <- names(GE)[2:349]
namesmiRNA   <- names(miRNA)[2:349]
namesProtein <- names(Protein)[2:404]
namesMeth    <- names(Meth)

namesExp     <- substr(namesExp,1,16)
namesmiRNA   <- substr(namesmiRNA,1,16)
namesProtein <- substr(namesProtein,1,16)

MatchProt    <- match(namesExp,namesProtein,nomatch=0)
MatchMeth    <- match(namesExp,namesMeth   ,nomatch=0)

miRNA        <- miRNA[,2:349]
miRNA.mat    <- as.matrix(miRNA[,order(namesmiRNA)])

Protein.mat  <- Protein[,2:404]
Protein.mat  <- as.matrix(Protein.mat[,MatchProt])

Meth.mat     <- as.matrix(Meth[,MatchMeth])

###################Data processing#############################
set.seed(1)
###Impute missing expression values
#load package 'impute' from CRAN
Exp.mat      <- as.matrix(GE[,2:349])
Exp.mat      <- impute.knn(Exp.mat) ##Impute missing values via KNN (K=10)
Exp.mat      <- Exp.mat$data

rowSums(miRNA.mat==0)
##Remove miRNAs with > 50% missing values
miRNA.mat    <- miRNA.mat[rowSums(miRNA.mat==0) < 348*0.5,]

processedExpression  <- Exp.mat[apply(Exp.mat,1,sd)>1.5,] ###Filter to select only most variable genes
processedMethylation <- sqrt(Meth.mat)    ##take square root of methylation data
processedmiRNA       <- log(1+miRNA.mat) ##take log of miRNA data
processedProtein     <- scale(Protein.mat,center=TRUE,scale=TRUE) #Column center/scale protein

save(processedExpression , file = 'processedExpression.RDa')
save(processedMethylation, file = 'processedMethylation.RDa')
save(processedmiRNA      , file = 'processedmiRNA.RDa')
save(processedProtein    , file = 'processedProtein.RDa')
