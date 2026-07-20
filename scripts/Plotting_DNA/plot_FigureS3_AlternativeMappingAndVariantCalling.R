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
strAtlSamples <- args[8]
strMedSamples <- args[9]
Rconfig <- args[10]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PCAValues_files <- unlist(strsplit(strPCAValues_files, " "))
PCAPropVar_files <- unlist(strsplit(strPCAPropVar_files, " "))
HetTotal_files <- unlist(strsplit(strHetTotal_files, " "))
AllSamples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
AtlSamples <- unlist(strsplit(strAtlSamples, " "))
MedSamples <- unlist(strsplit(strMedSamples, " "))

######################################################################
# Functions

plot_PCA <- function(vec1, lab1, vec2, lab2, main, cexplot, colors, shapes, pointlabels, legcolors, leglabels, legshapes){
	linelab <- 4
	cexlab <- 1.5
	cexpoints <- 3
	#if(cexplot == 1){
	#	linelab <- 0
	#	cexlab <- 0.8
	#	cexpoints <- 1.5
	#}
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

plot_HetPerPop <- function(Hetdf, atls, meds, main, ylab, xlabs, ylim){
	w <- 0.05
	flanking <- 0.2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, flanking*3+(length(atls)+length(meds))*w), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1)
	mtext(main, side = 3, line = 4, cex=1)
	pos <- flanking
	hetsA <- Hetdf[1,atls]
	hetsM <- Hetdf[1,meds]
	print("P value of a t.test comparing Atlantic and Mediterranean populations:")
	print(t.test(hetsA, hetsM)$p.value)
	for(s in c(1:length(atls))){
		polygon(c(pos,pos+w,pos+w,pos), c(0,0,hetsA[s],hetsA[s]), col=population.colors[1], border=NA)
		pos <- pos+w
	}
	midpA <- flanking+(pos-flanking)/2
	pos <- pos+flanking
	stM <- pos
 	for(s in c(1:length(meds))){
		polygon(c(pos,pos+w,pos+w,pos), c(0,0,hetsM[s],hetsM[s]), col=population.colors[2], border=NA)
		pos <- pos+w
	}
	midpM <- stM+(pos-stM)/2
   	#violin(Hetdf[1,atls], 1, w, NA, modif_alpha(population.colors[1]), 2, 2)
	#points(jitter(rep(1, length(atls)), amount=w/3), Hetdf[,atls], pch=21, bg=modif_alpha(population.colors[1],0.2), col=modif_alpha(population.colors[1]), cex=1)
    #violin(Hetdf[1,meds], 2, w, NA, modif_alpha(population.colors[2]), 2, 2)
	#points(jitter(rep(2, length(meds)), amount=w/3), Hetdf[,meds], pch=21, bg=modif_alpha(population.colors[2],0.2), col=modif_alpha(population.colors[2]), cex=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	#axis(1, at = c(midpA, midpM), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(1, at = c(midpA, midpM), labels=xlabs, lwd=NA, las=1, line=1, cex.axis=1.5)
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
	TypeRegions <- paste0(unlist(strsplit(unlist(strsplit(PCAValues_files[f], "Plotting_DNA/Observed_"))[2], "/PCA_short_variants_inBEDregions"))[1], "_", unlist(strsplit(unlist(strsplit(PCAValues_files[f], "PerSample_in"))[2], "Regions_PerAllSamples"))[1])
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
## Start plotting Figure SX
pdf(PDF, width=15, height=10)

par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)


## Heterozygosity per sample
print("Heterozygosity per sample")
layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)
par(mar=c(4,4,2,2))
Het_bwa_gatk <- read.table(grep("Observed_bwa_mem2_GATK.*Callable", HetTotal_files, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(Het_bwa_gatk) <- AllSamples
print(Het_bwa_gatk)
plot_HetPerPop(Het_bwa_gatk*100, AtlSamples, MedSamples, "bwa mem2 + GATK", "Heterozygosity (%)", populations, c(2,3.5))
Het_minm_gatk <- read.table(grep("Observed_minimap2_GATK.*Callable", HetTotal_files, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(Het_minm_gatk) <- AllSamples
print(Het_minm_gatk)
plot_HetPerPop(Het_minm_gatk*100, AtlSamples, MedSamples, "minimap2 + GATK", "Heterozygosity (%)", populations, c(2,3.5))
Het_bwa_freeb <- read.table(grep("Observed_bwa_mem2_freebayes.*Callable", HetTotal_files, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(Het_bwa_freeb) <- AllSamples
print(Het_bwa_freeb)
plot_HetPerPop(Het_bwa_freeb*100, AtlSamples, MedSamples, "bwa mem2 + freebayes", "Heterozygosity (%)", populations, c(2,3.5))
Het_minm_freeb <- read.table(grep("Observed_minimap2_freebayes.*Callable", HetTotal_files, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(Het_minm_freeb) <- AllSamples
print(Het_minm_freeb)
plot_HetPerPop(Het_minm_freeb*100, AtlSamples, MedSamples, "minimap2 + freebayes", "Heterozygosity (%)", populations, c(2,3.5))

## PCA
print("PCA")
plot_PCA(	PCAresults$bwa_mem2_GATK_Callable_Comp1, paste("PC1", round(PCApropvar$bwa_mem2_GATK_Callable[1], digits = 2), "%"), 
			PCAresults$bwa_mem2_GATK_Callable_Comp2, paste("PC2", round(PCApropvar$bwa_mem2_GATK_Callable[2], digits = 2), "%"), 
			"bwa mem2 + GATK", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$minimap2_GATK_Callable_Comp1, paste("PC1", round(PCApropvar$minimap2_GATK_Callable[1], digits = 2), "%"), 
			PCAresults$minimap2_GATK_Callable_Comp2, paste("PC2", round(PCApropvar$minimap2_GATK_Callable[2], digits = 2), "%"), 
			"minimap2 + GATK", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$bwa_mem2_freebayes_Callable_Comp1, paste("PC1", round(PCApropvar$bwa_mem2_freebayes_Callable[1], digits = 2), "%"), 
			PCAresults$bwa_mem2_freebayes_Callable_Comp2, paste("PC2", round(PCApropvar$bwa_mem2_freebayes_Callable[2], digits = 2), "%"), 
			"bwa mem2 + freebayes", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)
plot_PCA(	PCAresults$minimap2_freebayes_Callable_Comp1, paste("PC1", round(PCApropvar$minimap2_freebayes_Callable[1], digits = 2), "%"), 
			PCAresults$minimap2_freebayes_Callable_Comp2, paste("PC2", round(PCApropvar$minimap2_freebayes_Callable[2], digits = 2), "%"), 
			"minimap2 + freebayes", 1, MetData$ColorP, unname(c(Female = 19, Male = 17)[MetData$Sex]), NA, NA, NA, NA)


dev.off()
#############
## Report
write("### Mapping and calling comparison", file = REPORT, append = FALSE)
#############










