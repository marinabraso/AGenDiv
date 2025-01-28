#!/usr/bin/env Rscript



######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
Matricesstr <- args[1]
samplesorderFile <- args[2]
DividedSitesstr <- args[3]
OutSites <- args[4]
OutSamples <- args[5]
OutProp <- args[6]
strchrs <- args[7]
typeClust <- args[8] #PCA 
typeFeat <- args[9] #all/exons/introns/promoters/intergenic/genic 
AllOrPolym <- args[10] #all/polymorphic 
AllOrBiall <- args[11] #all/biallelic 
AllOrFreqint <- args[12] #all/specific major frequency intervals e.g. 0.5to0.6
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#script <- "scripts/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature.R"
#Matricesstr <- "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.numericmatrix.gz results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr17.genotype.numericmatrix.gz" 
#DividedSitesstr <- "results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.chr6_classified.genotype.gz results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.chr6_classified.genotype.gz"
#samplesorderFile <- "metadata/SamplesOrderInVCF.chr19.txt"
#OutSites <- "results/VariantAnalysis_DNA/Sample_clustering_PCA/Numericmatrix_PCvalues_PerVariant_UMAP_all_all_all.txt.gz"
#OutSamples <- "results/VariantAnalysis_DNA/Sample_clustering_PCA/Numericmatrix_PCvalues_PerSample_UMAP_all_all_all.txt.gz"
#OutProp <- "results/VariantAnalysis_DNA/Sample_clustering_PCA/Numericmatrix_PCprop_UMAP_all_all_all.txt.gz"
#strchrs <- "chr16 chr17 chr18 chr19"
#typeClust <- "UMAP"
#typeFeat <- "intergenic"
#AllOrPolym <- "polymorphic"
#AllOrBiall <- "biallelic"
#AllOrFreqint <- "0.5to0.6"
#####################


######################################################################
# Libraries & functions
if(typeClust == "PCA"){
	library(ade4)
}else if(typeClust == "UMAP"){
	library(umap)
}else{
	library(FactoMineR)	
}

print(sessionInfo())

Matricesbasename <- unlist(strsplit(unlist(strsplit(Matricesstr, " "))[1], ".chr\\d+."))[1]
Matricesextention <- unlist(strsplit(unlist(strsplit(Matricesstr, " "))[1], ".chr\\d+."))[2]
DividedSitesbasename <- unlist(strsplit(unlist(strsplit(DividedSitesstr, " "))[1], ".chr\\d+."))[1]
DividedSitesextention <- unlist(strsplit(unlist(strsplit(DividedSitesstr, " "))[1], ".chr\\d+."))[2]
OutSitesbasename <- unlist(strsplit(OutSites, ".gz"))[1]
OutSamplesbasename <- unlist(strsplit(OutSamples, ".gz"))[1]
OutPropbasename <- unlist(strsplit(OutProp, ".gz"))[1]
samples <- unlist(read.table(samplesorderFile, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F))
chrs <- unlist(strsplit(strchrs, " "))
#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
system(paste0("mkdir -p $(dirname ", OutSites, ")"))
print(typeClust)
print(typeFeat)
print(AllOrPolym)
print(AllOrBiall)
print(AllOrFreqint)


###########################################################################
# Functions

# Principal component analysis with ade4
PCA_calc <- function(Data){
	pca <- dudi.pca(Data, center=T, scale=T, scannf=F, nf=5)
	propvar <- 100 * pca$eig/sum(pca$eig)
	sampleCoord <- as.data.frame(pca$co)
	siteCoord <- as.data.frame(pca$li)
	return(list("pca"=pca, "propvar"=propvar, "sampleCoord"=sampleCoord, "siteCoord"=siteCoord))
}

# Principal component analysis with FactoMineR
PCA_calc_FactoMineR <- function(Data){
	pca <- PCA(Data, graph = F)
	propvar <- 100 * pca$eig/sum(pca$eig)
	sampleCoord <- as.data.frame(pca$ind)
	siteCoord <- as.data.frame(pca$var)
	return(list("pca"=pca, "propvar"=propvar, "sampleCoord"=sampleCoord, "siteCoord"=siteCoord))
}

# Uniform manifold approximation and projection (UMAP)
UMAP_calc <- function(Data){
	umap <- umap(t(Data), center=T, scale=T, scannf=F, nf=5)
	propvar <- NULL
	sampleCoord <- umap$layout
	colnames(sampleCoord) <- c("Comp1", "Comp2")
	siteCoord <- NULL
	return(list("umap"=umap, "propvar"=propvar, "sampleCoord"=sampleCoord, "siteCoord"=siteCoord))
}



###########################################################################
# Read data

cat("----> Getting list of selected variants\n")
header <- system(paste0("zcat ", DividedSitesbasename, ".", chrs[1], "_", DividedSitesextention, " | head -1 | sed 's/\\t/;/g'"), intern=TRUE)
header <- unlist(strsplit(header, ";"))
VarInfo <- data.frame(matrix(ncol = length(header)+1, nrow = 0))
for(c in c(1:length(chrs))){
	cat(paste0("   ", chrs[c], "\n"))
	chrinfo <- read.table(paste0(DividedSitesbasename, ".", chrs[c], "_", DividedSitesextention), h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
	chrinfo$chr <- chrs[c]
	print(head(chrinfo))
	VarInfo <- rbind(VarInfo, chrinfo)
}
colnames(VarInfo) <- c(header, "chr")
print(head(VarInfo))


cat("----> Reading alleles\n")
Alleles <- data.frame(matrix(ncol = length(samples)+2, nrow = 0))
for(c in c(1:length(chrs))){
	cat(paste0("   ", chrs[c], "\n"))
	chralleles <- read.table(paste0(Matricesbasename, ".", chrs[c], ".", Matricesextention), h=FALSE, sep = " ", check.names = F, stringsAsFactors = F)
	chralleles$chr <- chrs[c]
	Alleles <- rbind(Alleles, chralleles)
}
colnames(Alleles) <- c("num", samples, "chr")
print(head(Alleles))

cat("----> Filtering\n")

if(typeFeat == "all"){
	InFeat <- VarInfo[, c("Num","chr")]	
}else if (typeFeat == "exon" | typeFeat == "intron" | typeFeat == "promoter" | typeFeat == "intergenic"){
	InFeat <- VarInfo[which(VarInfo$Feature == typeFeat), c("Num","chr")]
}else if(typeFeat == "genic"){
	InFeat <- VarInfo[which(VarInfo$Feature != "intergenic"), c("Num","chr")]	
}else{
	cat(paste0("Error: unrecognized value of typeFeat: ", typeFeat))
	quit()
}

if(AllOrPolym == "all"){
	InPolym <- VarInfo[, c("Num","chr")]	
}else if(AllOrPolym == "polymorphic"){
	InPolym <- VarInfo[which(VarInfo$Polym == 1), c("Num","chr")]
}else{
	cat(paste0("Error: unrecognized value of AllOrPolym: ", AllOrPolym))
	quit()
}

if(AllOrBiall == "all"){
	InBiall <- VarInfo[, c("Num","chr")]	
}else if(AllOrBiall == "biallelic"){
	InBiall <- VarInfo[which(VarInfo$NumAlleles == 2), c("Num","chr")]
}else{
	cat(paste0("Error: unrecognized value of AllOrBiall: ", AllOrBiall))
	quit()
}

if(AllOrFreqint == "all"){
	InFreq <- VarInfo[, c("Num","chr")]	
}else if(length(unlist(strsplit(AllOrFreqint, "to")))==2){
	interval <- unlist(strsplit(AllOrFreqint, "to"))
	print(interval)
	InFreq <- VarInfo[which(VarInfo$MajorFreq/(length(samples)*2) >= interval[1] & VarInfo$MajorFreq/(length(samples)*2) < interval[2]), c("Num","chr")]
}else{
	cat(paste0("Error: unrecognized value of AllOrFreqint: ", AllOrFreqint))
	quit()
}


InAlleles <- Alleles[which(paste0(Alleles$num,Alleles$chr) %in% paste0(InFeat$Num,InFeat$chr) & paste0(Alleles$num,Alleles$chr) %in% paste0(InPolym$Num,InPolym$chr) & paste0(Alleles$num,Alleles$chr) %in% paste0(InBiall$Num,InBiall$chr) & paste0(Alleles$num,Alleles$chr) %in% paste0(InFreq$Num,InFreq$chr)),]
print(head(InAlleles))


######################################################################
# Clustering & printing

cat("----> Clustering\n")
if(typeClust == "PCA"){
	cat("         PCA\n")
	clust <- PCA_calc(InAlleles[,samples])
}else if(typeClust == "UMAP"){
	cat("         UMAP\n")
	clust <- UMAP_calc(InAlleles[,samples])
}else{
	cat(paste0("Error: unrecognized value of typeClust: ", typeClust))
	quit()
}

cat("----> Printing\n")
write.table(clust$sampleCoord, file=OutSamplesbasename, sep = "\t", col.names = TRUE, row.names = TRUE)
if(length(clust$siteCoord[,1])>0){
	write.table(cbind(InAlleles, clust$siteCoord), file=OutSitesbasename, sep = "\t", col.names = TRUE, row.names = TRUE)
}
if(length(clust$propvar)>0){
	write.table(clust$propvar, file=OutPropbasename, sep = "\t", col.names = TRUE, row.names = TRUE)
}

system(paste0("touch ", OutSamplesbasename))
system(paste0("touch ", OutSitesbasename))
system(paste0("touch ", OutPropbasename))
system(paste0("gzip ", OutSamplesbasename))
system(paste0("gzip ", OutSitesbasename))
system(paste0("gzip ", OutPropbasename))




