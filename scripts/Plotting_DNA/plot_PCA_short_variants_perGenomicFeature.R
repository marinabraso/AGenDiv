#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(ade4)
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
PCAplotPDF <- args[2]
strchrs <- args[3]
strsamples <- args[4]
Matricesbasename <- args[5]
DivListsbasename <- args[6]
Rconfig <- args[7]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#script <- "scripts/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature.R"
#MetadataFile <- "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#Matricesbasename <- "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered" 
#DivListsbasename <- "results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.classified" 
#strsamples <- "F10D;F1D;F2D;F3D;F4D;F6D;F7D;F8D;F9D;M2D;M3D;M4D;M5D;M6D;M7D;M8D;M9D;RF10D;RF1D;RF2D;RF3D;RF4D;RF5D;RF6D;RF7D;RF8D;RM10D;RM1D;RM2D;RM3D;RM4D;RM5D;RM6D;RM7D;RM8D;RU5D"
#PCAplotPDF <- "results/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/plot_PCA_short_variants_perGenomicFeature.pdf"
#strchrs <- "chr1;chr2;chr3;chr4;chr5;chr6;chr7;chr8;chr9;chr10;chr11;chr12;chr13;chr14;chr15;chr16;chr17;chr18;chr19"
#Rconfig <- "config/AmphiHetDupExp_plot.R"
#####################




chrs <- unlist(strsplit(strchrs, ";"))
samples <- unlist(strsplit(strsamples, ";"))
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
OutFolder <- system(paste0("echo $(dirname ", PCAplotPDF, ")"), intern=TRUE)

###########################################################################
# Read data
Metadata<-read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Metadata <- Metadata[match(Metadata$Sample, samples),]

Metadata <- Metadata[,c("Sample", "Population", "Sex")]
Metadata <- Metadata[order(Metadata$Sample), ]
Metadata <- Metadata[!duplicated(Metadata), ]
head(Metadata)
Populations <- unique(Metadata$Population)
Sexes <- unique(Metadata$Sex)
Metadata$ColorP <- population.colors[match(Metadata$Population, Populations)]
Metadata$ColorS <- sex.colors[match(Metadata$Sex, Sexes)]

#########################################################
cat("----> Reading types of variants\n")
typesFeature <- c("_exon","_intron","_promoter","_genic","_intergenic")
typesVariant <- c("_biallelic","_nonbiallelic")
Variables <- data.frame(matrix(ncol = 3, nrow = 0))
for(c in chrs){
	cat(paste0("   ", c, "\n"))
	for(t in typesFeature){
		print(paste0(DivListsbasename, ".", c, t,".multiallelic.txt"))
		var <- read.table(paste0(DivListsbasename, ".", c, t,".multiallelic.txt"), h=FALSE, sep = " ", check.names = F, stringsAsFactors = F)
		Variables <- rbind(Variables, cbind(var, rep(c, length(var)), rep(gsub("_","",t), length(var))))
	}
	for(t in typesVariant){
		print(paste0(DivListsbasename, ".", c, ".", gsub("_","",t), ".txt"))
		var <- read.table(paste0(DivListsbasename, ".", c, ".", gsub("_","",t), ".txt"), h=FALSE, sep = " ", check.names = F, stringsAsFactors = F)
		Variables <- rbind(Variables, cbind(var, rep(c, length(var)), rep(gsub("_","",t), length(var))))
	}
}
colnames(Variables) <- c("num","chr","type")
print(head(Variables))

#########################################################
cat("----> Reading alleles\n")
Alleles <- data.frame(matrix(ncol = length(samples)+2, nrow = 0))
for(c in c(1:length(chrs))){
	cat(paste0("   ", chrs[c], "\n"))
	chralleles <- read.table(paste0(Matricesbasename, ".", chrs[c], ".genotype.numericmatrix"), h=FALSE, sep = " ", check.names = F, stringsAsFactors = F)
	chralleles$chr <- chrs[c]
	Alleles <- rbind(Alleles, chralleles)
}
colnames(Alleles) <- c("num", samples, "chr")
print(head(Alleles))


######################################################################
# Plotting
cat("----> Plotting\n")
pdf(PCAplotPDF, width=20, height=20)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

cat("    Multiallelic\n")
pca <- PCA_calc(Alleles[,samples])
plot_x4_PCA(pca, "Mulitallelic variants", Metadata, population.colors, Populations, sex.colors, Sexes, range(pca$co.df$Comp1), range(pca$co.df$Comp2), range(pca$co.df$Comp3), range(pca$co.df$Comp4))
write.table(pca$co.df, file=paste0(dirname(PCAplotPDF), "/PCvalues_PerSample.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(cbind(Alleles, pca$li.df), file=paste0(dirname(PCAplotPDF), "/Numericmatrix_PCvalues_PerVariant.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(pca$propvar, file=paste0(dirname(PCAplotPDF), "/PCpropVariance.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
#rm(pca); gc()


#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Biallelic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="biallelic")], Variables$num[which(Variables$type=="biallelic")])),samples])
#plot_x4_PCA(pca, "Biallelic variants", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Non-biallelic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="nonbiallelic")], Variables$num[which(Variables$type=="nonbiallelic")])),samples])
#plot_x4_PCA(pca, "Non-biallelic variants", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Exonic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="exon")], Variables$num[which(Variables$type=="exon")])),samples])
#plot_x4_PCA(pca, "Exonic variants", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Intronic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="intron")], Variables$num[which(Variables$type=="intron")])),samples])
#plot_x4_PCA(pca, "Intronic variants", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Promoter\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="promoter")], Variables$num[which(Variables$type=="promoter")])),samples])
#plot_x4_PCA(pca, "Promoter variants", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Genic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type %in% c("exon", "intron", "promoter"))], Variables$num[which(Variables$type %in% c("exon", "intron", "promoter"))])),samples])
#plot_x4_PCA(pca, "Genic variants (promoter, exon, intron)", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()
#layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#cat("    Intergenic\n")
#pca <- PCA_calc(Alleles[which(paste(Alleles$chr, Alleles$num) %in% paste(Variables$chr[which(Variables$type=="intergenic")], Variables$num[which(Variables$type=="intergenic")])),samples])
#plot_x4_PCA(pca, "Intergenic variants (outside promoter, exon or intron)", Metadata, population.colors, Populations, sex.colors, Sexes)
#rm(pca); gc()



dev.off()
















