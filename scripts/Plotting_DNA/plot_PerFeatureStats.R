#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")


######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
strFeatMatrices <- args[2]
PDF <- args[3]
strchrs <- args[4]
script <- sub(".*=", "", commandArgs()[4])


#####################
# debug / develop
#MetadataFile = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#PDF="results/Plotting_DNA/plot_PerFeatureStats/plot_PerFeatureStats.pdf"
#strFeatMatrices= "results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr1.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr2.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr3.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr4.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr5.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr6.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr7.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr8.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr9.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr10.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr11.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr12.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr13.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr14.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr15.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr16.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr17.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr18.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr19.tab.gz"
#strchrs <- "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
#script <- "scripts/Plotting_DNA/plot_PerFeatureStats.R"
#####################

FeatMatrices <- unlist(strsplit(strFeatMatrices, " "))
chrs <- unlist(strsplit(strchrs, " "))
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
ppalcolor <- "royalblue4"
types.of.features <- c("exon", "intron", "promoter", "intergenic")
types.of.features.color <- c("seagreen4", "palegreen4", "olivedrab4","azure4")
types.of.CodingSites <- c("S", "N", "D")
pvalthresh <- c(0.05, 0.01, 0.001)





######################################################################
# Read data
Metadata<-read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Metadata <- unique(Metadata[,c("Sample","Population","Sex")])
Metadata <- Metadata[which(!(Metadata$Sample=="RU5D" & Metadata$Sex=="Male")),]
samples <- unique(Metadata$Sample)
pcolors <- c("purple", "orange")[match(Metadata$Population, unique(Metadata$Population))]

header <- system(paste0("zcat < ", FeatMatrices[1], " | head -1 | sed 's/\\t/;/g'"), intern=TRUE)
header <- unlist(strsplit(header, ";"))
Data <- data.frame(matrix(ncol = length(header)+2, nrow = 0))
colnames(Data) <- header
for(c in c(1:length(chrs))){
	chrData<-read.table(FeatMatrices[c], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	Data <- rbind(Data, chrData)
}
Data$len <- Data$end-Data$st
Data$totalhets <- rowSums(Data[,paste0("hets",samples)])
Data$propS <- Data$variantSynon/Data$totalSynon
Data$propN <- Data$variantnonSynon/Data$totalnonSynon
Data$propD <- Data$variantDeprec/Data$totalDeprec
Data <- Data[which(Data$callablesites>0),]

head(Metadata)
head(Data)

######################################################################
# Plotting
pdf(PDF, width=15, height=10)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))


layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_perGenomicFeature_boxplot(Data$callablesites/Data$len, Data$type, types.of.features, "Callable length / total length", types.of.features.color, c(0,1), pvalthresh, OutFolder)
plot_perGenomicFeature_boxplot(Data$varsites/Data$callablesites, Data$type, types.of.features, "Variant sites", types.of.features.color, c(0,1), pvalthresh, OutFolder)
plot_perGenomicFeature_boxplot(Data$polymsites/Data$callablesites, Data$type, types.of.features, "Plymorphic sites", types.of.features.color, c(0,0.5), pvalthresh, OutFolder)
plot_perGenomicFeature_boxplot(Data$totalhets/Data$callablesites/length(samples)*100, Data$type, types.of.features, "Mean heterozygosity (%)", types.of.features.color, c(0,8), pvalthresh, OutFolder)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_perGenomicFeature_violin(Data$callablesites/Data$len, Data$type, types.of.features, "Callable length / total length", types.of.features.color, c(0,1), pvalthresh, OutFolder)
plot_perGenomicFeature_violin(Data$varsites/Data$callablesites, Data$type, types.of.features, "Variant sites", types.of.features.color, c(0,1), pvalthresh, OutFolder)
plot_perGenomicFeature_violin(Data$polymsites/Data$callablesites, Data$type, types.of.features, "Plymorphic sites", types.of.features.color, c(0,0.5), pvalthresh, OutFolder)
plot_perGenomicFeature_violin(Data$totalhets/Data$callablesites/length(samples)*100, Data$type, types.of.features, "Mean heterozygosity (%)", types.of.features.color, c(0,8), pvalthresh, OutFolder)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_propSND_boxplot(Data$varsites/Data$callablesites, Data[,c("propS", "propN", "propD")], Data$type, types.of.features, types.of.CodingSites, "Variable sites / callable sites", types.of.features.color, c(0,1), pvalthresh, OutFolder)


plot_perGenomicFeature_perSample(Data, paste0("hets",samples), "type", "callablesites", types.of.features, "Heterozygosity", Metadata$Population, pcolors, c(0,0.05))




dev.off()






