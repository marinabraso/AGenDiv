#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
strMatrices <- args[2]
PDF <- args[3]
strchrs <- args[4]
Rconfig <- args[5]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#script <- "./scripts/Plotting_DNA/plot_PerSiteStats.R"
#MetadataFile <- "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt" 
#strMatrices <- "results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_10000_1000.chr17.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_10000_1000.chr18.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_10000_1000.chr19.tab.gz"
#PDF <- "results/Plotting_DNA/plot_PerWindowStats/plot_PerWindowStats.pdf"
#strchrs <- "chr17 chr18 chr19"
#Rconfig <- "config/AmphiHetDupExp_plot.R"
#####################







Matrices <- unlist(strsplit(strMatrices, " "))
chrs <- unlist(strsplit(strchrs, " "))
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
colfunc <- colorRampPalette(chr.colors)
chr.rampcolors <- colfunc(length(chrs))


######################################################################
# Read data
Metadata <- read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Metadata <- unique(Metadata[,c("Sample","Population","Sex")])
samples <- unique(Metadata$Sample)
pcolors <- population.colors[match(Metadata$Population, unique(Metadata$Population))]

header <- system(paste0("zcat ", Matrices[1], " | head -1 | sed 's/\\t/;/g'"), intern=TRUE)
header <- unlist(strsplit(header, ";"))
Data <- data.frame(matrix(ncol = length(header)+2, nrow = 0))
colnames(Data) <- header
for(c in c(1:length(chrs))){
	chrData<-read.table(Matrices[c], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	Data <- rbind(Data, chrData)
}
Data$propGC <- Data$GCcontent/(Data$ATcontent+Data$GCcontent)
Data$propExon <- Data$exoncontent/Data$callablelength
Data$proplowfreqVar <- (Data$NumVarFreq0-.05 + Data$NumVarFreq.05-.15)/Data$varsites


head(Metadata)
head(Data)
wlength <- Data$end[1]-Data$st[1]
Data <- Data[which(Data$valid==1),]


######################################################################
# Plotting
pdf(PDF, width=15, height=10)
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
par(mar=c(5,5,1,1),oma=c(1,1,1,1))

plot_scatter(NULL, Data$callablelength/wlength, "callable length / length", Data$PI, "PI", "black", c(0,1), c(0,.06))
plot_scatter(NULL, Data$callablelength/wlength, "callable length / length", Data$PC1all, "PC1", "black", c(0,1), c(-2,2))
plot_scatter(NULL, Data$callablelength/wlength, "callable length / length", Data$propGC, "GC content", "black", c(0,1), c(0.3,0.6))
plot_scatter(NULL, Data$callablelength/wlength, "callable length / length", Data$propExon, "Exon content", "black", c(0,1), c(0,1))
plot_scatter(NULL, Data$callablelength/wlength, "callable length / length", Data$proplowfreqVar, "Proportion of low frequency variants", "black", c(0,1), c(0,1))

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter(NULL, Data$PI, "Pi", Data$PC1all, "PC1", "black", c(0,.06), c(-2,2))
plot_scatter(NULL, Data$PI, "Pi", Data$propGC, "GC content", "black", c(0,.06), c(0.3,0.6))
plot_scatter(NULL, Data$PI, "Pi", Data$proplowfreqVar, "Proportion of low frequency variants", "black", c(0,.06), c(0,1))
plot_scatter(NULL, Data$propGC, "GC content", Data$PC1all, "PC1", "black", c(0.3,0.6), c(-2,2))

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter(NULL, Data$mean_pop_FST, "mean_pop_FST", Data$PC1all, "PC1", "black", c(0,0.1), c(-2,2))
plot_scatter(NULL, Data$mean_pop_FST, "mean_pop_FST", Data$PC2all, "PC2", "black", c(0,0.1), c(-.2,.2))
plot_scatter(NULL, Data$mean_sex_FST, "mean_sex_FST", Data$PC1all, "PC1", "black", c(0,0.02), c(-2,2))
plot_scatter(NULL, Data$mean_sex_FST, "mean_sex_FST", Data$PC2all, "PC2", "black", c(0,0.02), c(-.2,.2))


plot_Histogram(Data$PI*100, "Average pairwise differences (%)", "Relative frequency", c(0,0.05), c(0,7), 0.05)
plot_Histogram(Data$varsites/Data$callablelength*100, "Variant sites (%)", "Relative frequency", c(0,0.05), c(0,50), 0.5)
print(summary(Data$varsites/Data$callablelength*100))

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter(NULL, Data$PI*100, "Average pairwise differences (%)", Data$exoncontent/Data$callablelength*100, "Exon content (%)", "black", c(0,7), c(0,100))
plot_scatter(NULL, Data$varsites/Data$callablelength*100, "Variant sites (%)", Data$exoncontent/Data$callablelength*100, "Exon content (%)", "black", c(0,50), c(0,100))


dev.off()






