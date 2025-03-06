#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)
library(ade4)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPCAValues_files <- args[1]
strPCAPropVar_files <- args[2]
strObsShPriv <- args[3]
strRandShPriv <- args[4]
strPSMC <- args[5]
SamplesOrderInVCF <- args[6]
Metadata <- args[7]
PDF <- args[8]
REPORT <- args[9]
strchrs <- args[10]
Rconfig <- args[11]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PCAValues_files <- unlist(strsplit(strPCAValues_files, " "))
PCAPropVar_files <- unlist(strsplit(strPCAPropVar_files, " "))
ObsShPriv <- unlist(strsplit(strObsShPriv, " "))
RandShPriv <- unlist(strsplit(strRandShPriv, " "))
PSMC <- unlist(strsplit(strPSMC, " "))
chrs <- unlist(strsplit(strchrs, " "))

######################################################################
# Functions

read_set_of_files <- function(files){
	data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	colnames(data) <- seq(1,length(colnames(data)),1)
	if(length(files)>1){
		for(f in c(2:length(files))){
			fdata <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
			colnames(fdata) <- seq(1,length(colnames(fdata)),1)
			data <- rbind(data, fdata)
		}
	}
	return(data)
}

######################################################################
# Read data

# Metadata
MetData <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
MetData <- unique(MetData[,c("Sample", "Population")])
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
MetData <- MetData[match(Samples, MetData$Sample),]
print(head(MetData)) 

# Reading PCA results for several goups of variants
cat("----> Reading PCA results for several goups of variants\n")
PCAresults <- data.frame(matrix(ncol = 0, nrow = length(Samples)))
rownames(PCAresults) <- Samples
TypesOfRegions <- c()
for(f in c(1:length(PCAValues_files))){
	TypeRegions <- unlist(strsplit(unlist(strsplit(PCAValues_files[f], "PerSample_in"))[2], "Regions_PerAllSamples"))[1]
	val <- read.table(PCAValues_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
	pvar <- read.table(PCAPropVar_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
	colnames(val) <- paste(TypeRegions, colnames(val), sep="_")
	colnames(pvar) <- paste(TypeRegions, "PropVar", sep="_")
	PCAresults <- cbind(PCAresults, val, pvar)
	TypesOfRegions <- c(TypesOfRegions, TypeRegions)
}
print(head(PCAresults))
print(TypesOfRegions)
quit()

######################################################################
# Plotting
pdf(PDF, width=10, height=15)

par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,4,4),nrow=3,ncol=2,byrow=F), widths=c(6,6), heights=c(2,2,2), TRUE)


# A
par(mar=c(0,8,4,5))
plot_ObsSim("Standard deviation\nof heterozygosity", ObsDataM$sd, ObsDataA$sd, SimHet, "sd", colNe, c(0,.1))
writePlotLabel("A")
# B
par(mar=c(0,8,4,5))
plot_ObsSim("Propotion of\nsingletons + doubletons", SFSAtl$PropSingleDoubletons, SFSMed$PropSingleDoubletons, SimSFS, "PropSingleDoubletons", colNe, c(0,1))
writePlotLabel("B")
# Legend
uNe <- unique(SimHet$Ne)
umu <- unique(SimHet$mu)
w <- .5
par(mar=c(0,7,0,5))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", xlim=c(0,length(uNe)+3), col=NA)
for(i in c(1:length(uNe))){
	text(i+2, 7, labels = NeValForm[which(NeVal==uNe[i])], xpd=NA)
	text(i+2, 5, labels = muValForm[which(muVal==umu[i])], xpd=NA)
	polygon(c(i+2-w/2, i+2+w/2, i+2+w/2, i+2-w/2), c(1, 1, 1+w*2, 1+w*2), col=colNe[i], border=NA, lwd=1)
}
text(1, 7, labels = expression("N"["e"]), xpd=NA)
text(1, 5, labels = expression(mu), xpd=NA)
text(1, 10, labels = "Atlantic", xpd=NA)
text(2, 10, labels = "Mediterranean", xpd=NA)
text(6, 10, labels = "Simulations", xpd=NA)

# C
par(mar=c(4,4,4,4))
plot_PCA("", 
		pca$co.df$Comp1, paste0("PC1 (", round(pca$propvar[1], digits = 2), "%)"), 
		pca$co.df$Comp2, paste0("PC2 (", round(pca$propvar[2], digits = 2), "%)"), 
		pcacolors, pcapch, pcacex, NA, c("black", "black", colNe), c("Atlantic", "Mediterranean", Nemu.combinations), c(0,1), c(-0.6,.3))
writePlotLabel("C")





dev.off()















