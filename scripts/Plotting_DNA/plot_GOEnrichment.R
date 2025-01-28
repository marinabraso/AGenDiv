#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(fgsea)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strFeatMatrices <- args[1]
GOAnnotation <- args[2]
PDF <- args[3]
Rconfig <- args[4]
script <- sub(".*=", "", commandArgs()[4])


#strFeatMatrices <- "results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr1.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr2.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr3.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr4.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr5.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr6.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr7.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr8.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr9.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr10.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr11.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr12.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr13.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr14.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr15.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr16.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr17.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr18.tab.gz results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr19.tab.gz"
#GOAnnotation <- "data/BraLan3Gene_2_GOterm.txt"
#PDF <- "results/Plotting_DNA/plot_GOEnrichment/plot_GOEnrichment.pdf"
#Rconfig <- "config/AmphiHetDupExp_plot.R"


source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
FeatMatrices <- unlist(strsplit(strFeatMatrices, " "))


######################################################################
# Functions

read_set_of_Matrices <- function(files){
	data <- read.table(files[1], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	if(length(files)>1){
		for(f in c(2:length(files))){
			fdata <- read.table(files[f], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
			data <- rbind(data, fdata)
		}
	}
	return(data)
}

######################################################################
# Read data

FeatData <- read_set_of_Matrices(FeatMatrices)
FeatData <- FeatData[which(FeatData$type=="exon" & FeatData$callablesites > 10),]
FeatData$PropVar <- FeatData$varsites / FeatData$callablesites
FeatData$PropMeanHet <- rowMeans(FeatData[,grep("hets", colnames(FeatData))]) / FeatData$callablesites
FeatData$Gene <- unlist(lapply(FeatData$name, function(x){unlist(strsplit(x, "_"))[1]}))
print(head(FeatData))
GeneValues <- aggregate(FeatData[, c("PropVar","PropMeanHet")], list(FeatData$Gene), mean)
colnames(GeneValues) <- c("Gene", "PropVar", "PropMeanHet")
print(head(GeneValues))
ListPropVar <- 1:length(GeneValues[,1])
names(ListPropVar) <- GeneValues$Gene[order(GeneValues$PropVar, decreasing = TRUE)]
print(head(ListPropVar))
ListPropMeanHet <- 1:length(GeneValues[,1])
names(ListPropMeanHet) <- GeneValues$Gene[order(GeneValues$PropMeanHet, decreasing = TRUE)]
print(head(ListPropMeanHet))

GOAnn <- unique(read.table(GOAnnotation, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F, quote="\""))
colnames(GOAnn) <- c("Gene", "GO", "Description")
print(head(GOAnn))
GeneUniverse <- unique(GOAnn[,1])
GOAnnList <- split(GOAnn$Gene,GOAnn$GO)
print(head(GOAnnList))

fgseaPropVar <- fgsea(pathways = GOAnnList, 
	stats    = ListPropVar,
	minSize  = 15,
	maxSize  = 500,
	scoreType = "pos")
fgseaPropMeanHet <- fgsea(pathways = GOAnnList, 
	stats    = ListPropMeanHet,
	minSize  = 15,
	maxSize  = 500,
	scoreType = "pos")
SignifGOPropVar <- fgseaPropVar$pathway[which(fgseaPropVar$padj<0.01)]
SignifGOPropMeanHet <- fgseaPropMeanHet$pathway[which(fgseaPropMeanHet$padj<0.01)]
head(SignifGOPropVar)
head(SignifGOPropMeanHet)
print(unique(GOAnn[which(GOAnn$GO %in% SignifGOPropVar),c("GO", "Description")]))



######################################################################
# Plotting
pdf(PDF, width=10, height=5)
par(mar=c(7,3,2,2),oma=c(1,3,1,1))

#layout(matrix(seq(1,6,1),nrow=3,ncol=length(Types),byrow=T), TRUE)



dev.off()















