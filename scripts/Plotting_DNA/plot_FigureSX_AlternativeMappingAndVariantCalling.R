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
strHetTotal_files <- args[3]
SamplesOrderInVCF <- args[4]
Metadata <- args[5]
PDF <- args[6]
REPORT <- args[7]
Rconfig <- args[8]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PCAValues_files <- unlist(strsplit(strPCAValues_files, " "))
PCAPropVar_files <- unlist(strsplit(strPCAPropVar_files, " "))
HetTotal_files <- unlist(strsplit(strHetTotal_files, " "))

######################################################################
# Functions

plot_PCA <- function(vec1, lab1, vec2, lab2, main, cexplot, colors, shapes, pointlabels, legcolors, leglabels, legshapes){
	linelab <- 4
	cexlab <- 1.5
	cexpoints <- 3
	if(cexplot == 1){
		linelab <- 0
		cexlab <- 0.8
		cexpoints <- 1.5
	}
	vec1 <- as.numeric(vec1)
	vec2 <- as.numeric(vec2)
	ylim <- c(min(vec2), max(vec2))
	xlim <- c(min(vec1), max(vec1))
	yflank <- (ylim[2]-ylim[1])*0.05
	xflank <- (xlim[2]-xlim[1])*0.05
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-yflank,ylim[2]+yflank), xlim=c(xlim[1]-xflank,xlim[2]+xflank), col=NA)
	mtext(lab1, side = 1, line = linelab, cex=cexlab)
	mtext(lab2, side = 2, line = linelab, cex=cexlab)
	mtext(main, side = 3, line = linelab, cex=cexlab)
	points(vec1, vec2, pch=shapes, col=colors, bg=modif_alpha(colors), cex=cexpoints)
	text(vec1, vec2, labels=pointlabels, font=1, cex=.7)
	if(length(leglabels)>1 & !is.na(leglabels[1])){
		legend("right", as.character(leglabels), pch=legshapes, text.col="black", col=legcolors, bty = "n", cex=2, xjust = 0, yjust = 0)
	}
	#axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=2, cex.axis=1.5)
	#axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

######################################################################
# Read data

# Metadata
MetData <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
MetData <- unique(MetData[,c("Sample", "Population", "Sex")])
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
MetData <- MetData[match(Samples, MetData$Sample),]
MetData$ColorP <- population.colors[match(MetData$Population, c("Roscoff", "Banyuls"))]
MetData$ColorS <- sex.colors[match(MetData$Sex, unique(MetData$Sex))]
print(MetData)

# Reading PCA results for several goups of variants
cat("----> Reading PCA results for several goups of variants\n")
PCAresults <- data.frame(matrix(ncol = 0, nrow = length(Samples)))
rownames(PCAresults) <- Samples
PCApropvar <- data.frame(matrix(ncol = 0, nrow = 10))
TypesOfRegions <- c()
for(f in c(1:length(PCAValues_files))){
	TypeRegions <- unlist(strsplit(unlist(strsplit(PCAValues_files[f], "PerSample_in"))[2], "Regions_PerAllSamples"))[1]
	val <- read.table(PCAValues_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
	pvar <- as.data.frame(read.table(PCAPropVar_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)[c(1:10),])
	colnames(val) <- paste(TypeRegions, colnames(val), sep="_")
	colnames(pvar) <- TypeRegions
	PCAresults <- cbind(PCAresults, val)
	PCApropvar <- cbind(PCApropvar, pvar)
	TypesOfRegions <- c(TypesOfRegions, TypeRegions)
}
print(head(PCAresults))
print(head(PCApropvar))
print(TypesOfRegions)

######################################################################
## Start plotting Figure 2
pdf(PDF, width=15, height=10)

par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)

### A
## PCA
par(mar=c(7,7,2,2))
plot_PCA(	PCAresults$Callable_Comp1, paste("PC1", round(PCApropvar$Callable[1], digits = 2), "%"), 
			PCAresults$Callable_Comp2, paste("PC2", round(PCApropvar$Callable[2], digits = 2), "%"), 
			"", 3, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, c(population.colors, "black", "black"), c(populations, "Female", "Male"), c(15,15,19,17))
writePlotLabel("A")

### B
## Distance tree
par(mar=c(3,3,2,2))
tree <- readPNG(paste0(system("pwd", intern=TRUE), "/", TreeImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(tree, 1, 1, 10, 10)
writePlotLabel("B")


### C
## SHared Private variants analysis
par(mar=c(7,7,2,2))
plot_violin_SeveralRandomVsObserved(c("PrivateA", "PrivateM", "VariantShared"), OSharedPriv, RSharedPriv, c(9,20)*1000000, c(population.colors, "forestgreen"), 1000000, "Millions of variants", c("Private\nMediterranean", "Private\nAtlantic", "Shared"))
writePlotLabel("C")

### D
## PSMC
#plot_PSMC_samples_wo_mu(PSMCdata, MetData$Sample, c(PSMChighmu, PSMClowmu), Blangenerationtime, MetData$ColorP)
#plot.new()
#writePlotLabel("D")


### Supplementary figure 3
## PCA for different sets of variants
print("Plotting Supplementary figure 3")
par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)
par(mar=c(3,3,1,1))
plot_PCA(	PCAresults$SNPs_Comp1, paste("PC1", round(PCApropvar$SNPs[1], digits = 2), "%"), 
			PCAresults$SNPs_Comp2, paste("PC2", round(PCApropvar$SNPs[2], digits = 2), "%"), 
			"SNPs", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$INDELs_Comp1, paste("PC1", round(PCApropvar$INDELs[1], digits = 2), "%"), 
			PCAresults$INDELs_Comp2, paste("PC2", round(PCApropvar$INDELs[2], digits = 2), "%"), 
			"INDELs", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$Exons_Comp1, paste("PC1", round(PCApropvar$Exons[1], digits = 2), "%"), 
			PCAresults$Exons_Comp2, paste("PC2", round(PCApropvar$Exons[2], digits = 2), "%"), 
			"Exonic variants", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$Introns_Comp1, paste("PC1", round(PCApropvar$Introns[1], digits = 2), "%"), 
			PCAresults$Introns_Comp2, paste("PC2", round(PCApropvar$Introns[2], digits = 2), "%"), 
			"Intronic variants", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$Promoters_Comp1, paste("PC1", round(PCApropvar$Promoters[1], digits = 2), "%"), 
			PCAresults$Promoters_Comp2, paste("PC2", round(PCApropvar$Promoters[2], digits = 2), "%"), 
			"Promoter variants", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$Intergenic_Comp1, paste("PC1", round(PCApropvar$Intergenic[1], digits = 2), "%"), 
			PCAresults$Intergenic_Comp2, paste("PC2", round(PCApropvar$Intergenic[2], digits = 2), "%"), 
			"Intergenic variants", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$BiAllelic_Comp1, paste("PC1", round(PCApropvar$BiAllelic[1], digits = 2), "%"),
			PCAresults$BiAllelic_Comp2, paste("PC2", round(PCApropvar$BiAllelic[2], digits = 2), "%"), 
			"Biallelic variants", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot.new()
plot_PCA(	PCAresults[["AbsFreq_36-40_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_36-40_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_36-40_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_36-40_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [36,40)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_40-45_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_40-45_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_40-45_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_40-45_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [40,45)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_45-50_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_45-50_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_45-50_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_45-50_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [45,50)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_50-55_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_50-55_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_50-55_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_50-55_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [50,55)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_55-60_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_55-60_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_55-60_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_55-60_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [55,60)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_60-65_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_60-65_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_60-65_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_60-65_AllVariants"]][2], digits = 2), "%"), 
			"AAF in the interval [60,65)", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults[["AbsFreq_65-73_AllVariants_Comp1"]], paste("PC1", round(PCApropvar[["AbsFreq_65-73_AllVariants"]][1], digits = 2), "%"), 
			PCAresults[["AbsFreq_65-73_AllVariants_Comp2"]], paste("PC2", round(PCApropvar[["AbsFreq_65-73_AllVariants"]][2], digits = 2), "%"),
			"AAF in the interval [65,71]", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)


dev.off()
#############
## Report
print(OSharedPriv)
## Write Table S3 (Shared Private)
write("### Table S3 (Shared Private)", file = REPORT, append = FALSE)
tableS3 <- OSharedPriv
tableS3$PercentTotal <- tableS3$Observed/sum(tableS3$Observed)*100
tableS3$RMeans <- rowMeans(RSharedPriv)
pvals <- c()
sds <- c()
for(r in c(1:length(RSharedPriv[,1]))){
	print(rownames(tableS3)[r])
	moreextreme <- 0
	print(tableS3$Observed[r])
	print(mean(unlist(RSharedPriv[r,])))
	if(tableS3$Observed[r] > mean(unlist(RSharedPriv[r,]))){
		moreextreme <- length(which(tableS3$Observed[r] <= RSharedPriv[r,]))
	} else {
		moreextreme <- length(which(tableS3$Observed[r] >= RSharedPriv[r,]))
	}
	print(moreextreme)
	p_val <- moreextreme/length(RSharedPriv[r,])
	if(p_val==0){
		p_val <- paste0("<", 1/length(RSharedPriv[r,]))
	}
	pvals <- c(pvals, p_val)
	sd_val <- sd(as.numeric(RSharedPriv[r,]))
	sds <- c(sds, sd_val)
}
tableS3$RSds <- sds
tableS3$Rp_val <- pvals
print(tableS3)
write.table(tableS3, file = REPORT, append = TRUE, row.names=TRUE, sep="\t", quote = FALSE)










