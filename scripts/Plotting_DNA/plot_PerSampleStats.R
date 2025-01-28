#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
SampleMatrix <- args[2]
chrSampleMatrices <- args[3]
PCApropFile <- args[4]
PDF <- args[5]
strchrs <- args[6]
script <- sub(".*=", "", commandArgs()[4])

#####################
# develop / debug
#script <- "./scripts/Plotting_DNA/plot_PerSampleStats.R"
#MetadataFile <- "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#SampleMatrix <- "results/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSampleStats.tab.gz"
#chrSampleMatrices <- "results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr1.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr2.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr3.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr4.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr5.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr6.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr7.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr8.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr9.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr10.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr11.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr12.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr13.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr14.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr15.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr16.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr17.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr18.tab results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr19.tab"
#PCApropFile <- "results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCprop_PCA_all_all_all_all.txt.gz"
#PDF <- "results/Plotting_DNA/plot_PerSampleStats/plot_PerSampleStats.pdf"
#strchrs <- "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
#####################


source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
system(paste0("mkdir -p $(dirname ", PDF, ")"))
chrSampleMatrices <- unlist(strsplit(chrSampleMatrices, " "))
chrs <- unlist(strsplit(strchrs, " "))
PCApropBasename <- unlist(strsplit(PCApropFile, "_all_all_all_all"))[1]

######################################################################
# Read data
Metadata <- read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
samples <- unique(Metadata$Sample)
pcolors <- c("purple", "orange")[match(Metadata$Population, unique(Metadata$Population))]

Data <- read.table(SampleMatrix, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
colnames(Data) <- system(paste0("zcat ", SampleMatrix, " | head -1 | sed 's/(/_/g' | sed 's/)//g' | sed 's/\\t/\\n/g'"), intern = TRUE)
Data <- Data[match(samples, Data$Sample),]

Data$propminor <- Data$minalleles / Data$callsites *100
Data$propsingletons <- Data$singletons / Data$callsites *100
Data$propdoubletons <- Data$doubletons / Data$callsites *100
Data$proptripletons <- Data$tripletons / Data$callsites *100
WGpi <- Data$pi[1]


chrData <- data.frame(matrix(ncol = length(colnames(Data))+1, nrow = 0))
colnames(chrData) <- c(colnames(Data),"chr")
cWGpi <- c()
for(c in c(1:length(chrs))){
	cData <- read.table(chrSampleMatrices[c], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	cData$chr <- rep(chrs[c], length(cData[,1]))
	cData <- cData[match(samples, cData$Sample),]
	chrData <- rbind(chrData, cData)
	cWGpi <- c(cWGpi, cData$pi[1])
}
chrData$propminor <- chrData$minalleles / chrData$callsites *100
chrData$propsingletons <- chrData$singletons / chrData$callsites *100
chrData$propdoubletons <- chrData$doubletons / chrData$callsites *100
chrData$proptripletons <- chrData$tripletons / chrData$callsites *100


head(Metadata)
head(Data)
head(chrData)


PCAcolnames <- c("all", "polym", "bial", "exon", "0.5to0.6", "0.6to0.7", "0.7to0.8", "0.8to0.9", "0.9to1")
PCAfilenames <- c("all_all_all_all", "all_polymorphic_all_all", "all_all_biallelic_all", "exon_all_all_all", "all_all_all_0.5to0.6", "all_all_all_0.6to0.7", "all_all_all_0.7to0.8", "all_all_all_0.8to0.9", "all_all_all_0.9to1") 
for(s in c(1:length(PCAcolnames))){
	PCAprop.s <- read.table(paste0(PCApropBasename, "_", PCAfilenames[s], ".txt.gz"), sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	colnames(PCAprop.s) <- PCAcolnames[s]
	if(exists("PCAprop")){
		PCAprop <- cbind(PCAprop,PCAprop.s)
	}else{
		PCAprop <- PCAprop.s
	}

}
head(PCAprop)


######################################################################
# Plotting
pdf(PDF, width=20, height=10)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2,2), heights=c(1,1), TRUE)

#barplot_samples(Data$Sample, Data, "hetvar", "# of heterozyguous sites / # variant sites *100", c(0,20), pcolors)
barplot_samples(Data$Sample, Data, "hetcall", "# of heterozyguous sites / # callable sites *100", c(0,5), pcolors, c(WGpi*100, 0.156, 0.843))
#barplot_samples(Data$Sample, Data, "Fstat", "1- observed heterozyguous sites / expected heterozyguous sites", c(0,.3), pcolors)
barplot_samples(Data$Sample, Data, "propminor", "# of minor allales / # callable sites *100", c(0,5), pcolors)
barplot_samples(Data$Sample, Data, "propsingletons", "# of singletons / # callable sites *100", c(0,.5), pcolors)
barplot_samples(Data$Sample, Data, "propdoubletons", "# of doubletons / # callable sites *100", c(0,.5), pcolors)
barplot_samples(Data$Sample, Data, "proptripletons", "# of tripletons / # callable sites *100", c(0,.5), pcolors)

layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2,2), heights=c(1,1), TRUE)
#barplot_samples_chrs(samples, chrs, chrData, "hetvar", "# of heterozyguous sites / # variant sites *100", c(0,20), pcolors)
barplot_samples_chrs(samples, chrs, chrData, "hetcall", "# of heterozyguous sites / # callable sites *100", c(0,5), pcolors)
#barplot_samples_chrs(samples, chrs, chrData, "Fstat", "1-observed heterozyguous sites / expected heterozyguous sites", c(0,.3), pcolors)
barplot_samples_chrs(samples, chrs, chrData, "propminor", "# of minor allales / # callable sites *100", c(0,.5), pcolors)
barplot_samples_chrs(samples, chrs, chrData, "propsingletons", "# of singletons / # callable sites *100", c(0,.5), pcolors)
barplot_samples_chrs(samples, chrs, chrData, "propdoubletons", "# of doubletons / # callable sites *100", c(0,.5), pcolors)
barplot_samples_chrs(samples, chrs, chrData, "proptripletons", "# of tripletons / # callable sites *100", c(0,.5), pcolors)


layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2,2), heights=c(1,1), TRUE)
#points_samples_chrs(samples, chrs, chrData, "hetvar", "# of heterozyguous sites / # variant sites *100", c(5,20))
points_samples_chrs(samples, chrs, chrData, "hetcall", "# of heterozyguous sites / # callable sites *100", c(2,4))
#points_samples_chrs(samples, chrs, chrData, "Fstat", "1-observed heterozyguous sites / expected heterozyguous sites", c(0,.3))
points_samples_chrs(samples, chrs, chrData, "propminor", "# of minor allales / # callable sites *100", c(2,6))
points_samples_chrs(samples, chrs, chrData, "propsingletons", "# of singletons / # callable sites *100", c(0,2))
points_samples_chrs(samples, chrs, chrData, "propdoubletons", "# of doubletons / # callable sites *100", c(0,2))
points_samples_chrs(samples, chrs, chrData, "proptripletons", "# of tripletons / # callable sites *100", c(0,2))
points_samples_chrs_legend(chrs)

layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2,2), heights=c(1,1), TRUE)
pi_perChr(chrs, cWGpi, WGpi, "Average pairwise differences (pi)", c(0,.05))

PC1lim <- c(-1.2,-0.8)
PC2lim <- c(-0.15,0.15)
PC3lim <- c(-.3,.2)
PC4lim <- c(-.2,.2)
layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=F), widths=c(1,1,1,1), heights=c(1,1), TRUE)
par(mar=c(7,7,4,4),oma=c(1,1,1,1))
PCAtypes <- c("all", "polym", "bial", "exon")
PCAlabels <- c("All variants", "Polymorphisms", "Biallelic variants", "Exonic variants")
for(i in c(1:length(PCAtypes))){
	plot_scatter(PCAlabels[i], Data[, paste0("PC1_", PCAtypes[i])], paste0("PC1 (", format(round(PCAprop[1, PCAtypes[i]], 2), nsmall = 2), "%)"), Data[, paste0("PC2_", PCAtypes[i])], paste0("PC2 (", format(round(PCAprop[2, PCAtypes[i]], 2), nsmall = 2), "%)"), pcolors, c(min(Data[, paste0("PC1_", PCAtypes[i])]),max(Data[, paste0("PC1_", PCAtypes[i])])), c(min(Data[, paste0("PC2_", PCAtypes[i])]),max(Data[, paste0("PC2_", PCAtypes[i])])), 3, Data$Sample)
	plot_scatter(PCAlabels[i], Data[, paste0("PC3_", PCAtypes[i])], paste0("PC3 (", format(round(PCAprop[3, PCAtypes[i]], 2), nsmall = 2), "%)"), Data[, paste0("PC4_", PCAtypes[i])], paste0("PC4 (", format(round(PCAprop[4, PCAtypes[i]], 2), nsmall = 2), "%)"), pcolors, c(min(Data[, paste0("PC3_", PCAtypes[i])]),max(Data[, paste0("PC3_", PCAtypes[i])])), c(min(Data[, paste0("PC4_", PCAtypes[i])]),max(Data[, paste0("PC4_", PCAtypes[i])])), 3, Data$Sample)
}
layout(matrix(c(1,2,3,4,5,6,7,8,9,10),nrow=2,ncol=5,byrow=F), widths=c(1,1,1,1,1), heights=c(1,1), TRUE)
par(mar=c(7,7,4,4),oma=c(1,1,1,1))
PCAtypes <- c("0.5to0.6", "0.6to0.7", "0.7to0.8", "0.8to0.9", "0.9to1")
PCAlabels <- c("Major frequency from 0.5 to 0.6", "Major frequency from 0.6 to 0.7", "Major frequency from 0.7 to 0.8", "Major frequency from 0.8 to 0.9", "Major frequency from 0.9 to 1")
for(i in c(1:length(PCAtypes))){
	plot_scatter(PCAlabels[i], Data[, paste0("PC1_", PCAtypes[i])], paste0("PC1 (", format(round(PCAprop[1, PCAtypes[i]], 2), nsmall = 2), "%)"), Data[, paste0("PC2_", PCAtypes[i])], paste0("PC2 (", format(round(PCAprop[2, PCAtypes[i]], 2), nsmall = 2), "%)"), pcolors, c(min(Data[, paste0("PC1_", PCAtypes[i])]),max(Data[, paste0("PC1_", PCAtypes[i])])), c(min(Data[, paste0("PC2_", PCAtypes[i])]),max(Data[, paste0("PC2_", PCAtypes[i])])), 3, Data$Sample)
	plot_scatter(PCAlabels[i], Data[, paste0("PC3_", PCAtypes[i])], paste0("PC3 (", format(round(PCAprop[3, PCAtypes[i]], 2), nsmall = 2), "%)"), Data[, paste0("PC4_", PCAtypes[i])], paste0("PC4 (", format(round(PCAprop[4, PCAtypes[i]], 2), nsmall = 2), "%)"), pcolors, c(min(Data[, paste0("PC3_", PCAtypes[i])]),max(Data[, paste0("PC3_", PCAtypes[i])])), c(min(Data[, paste0("PC4_", PCAtypes[i])]),max(Data[, paste0("PC4_", PCAtypes[i])])), 3, Data$Sample)
}
layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=F), widths=c(1,1,1,1), heights=c(1,1), TRUE)
par(mar=c(7,7,4,4),oma=c(1,1,1,1))
plot_scatter("UMAP all variants", Data[, "UMAP1_all"], "UMAP1", Data[, "UMAP2_all"], "UMAP2", pcolors, c(min(Data[, "UMAP1_all"]),max(Data[, "UMAP1_all"])), c(min(Data[, "UMAP2_all"]),max(Data[, "UMAP2_all"])), 3, Data$Sample)

dev.off()


