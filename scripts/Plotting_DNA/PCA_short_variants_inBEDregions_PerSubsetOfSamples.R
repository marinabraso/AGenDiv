#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")


######################################################################
# Libraries & functions
library(ade4)
library(RColorBrewer)

######################################################################
# Files & folders
args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
GENOTYPEMatricesTMPfile <- args[1]
SamplesOrderInVCF <- args[2]
PCA_PerSample <- args[3]
PCA_PerVariant <- args[4]
PCA_PropVariance <- args[5]
strsamples <- args[6]
script <- sub(".*=", "", commandArgs()[4])

system(paste0("mkdir -p $(dirname ", PCA_RData, ")"))
OutFolder <- system(paste0("echo $(dirname ", PCA_RData, ")"), intern=TRUE)
SelectedSamples <- unlist(strsplit(strsamples, " "))

######################################################################
# Functions

PCA_calc <- function(Data){
	pca <- dudi.pca(Data, center=T, scale=T, scannf=F, nf=5)
	propvar <- 100 * pca$eig/sum(pca$eig)
	co.df <- as.data.frame(pca$co)
	li.df <- as.data.frame(pca$li)
	return(list("pca"=pca, "propvar"=propvar, "co.df"=co.df, "li.df"=li.df))
}

######################################################################
# Read data

# Samples Order
SamplesOrder <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)

# Reading gentoype matrices
cat("----> Reading gentoype matrices\n")
Data <- read.table(GENOTYPEMatricesTMPfile, sep=" ", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(Data) <- c("Var", SamplesOrder)

# Select columns corresponding to the set of samples
Data <- Data[,c("Var", SelectedSamples)]
print(head(Data))

######################################################################
# Running PCA

pca <- PCA_calc(Data[,SelectedSamples])
write.table(pca$co.df, file=PCA_PerSample, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(cbind(Data[,1], pca$li.df), file=PCA_PerVariant, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(pca$propvar, file=PCA_PropVariance, sep = "\t", col.names = TRUE, row.names = TRUE)
















