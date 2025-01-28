#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(RColorBrewer)
library(ash)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
strchrFilterStats <- args[2]
strSiteMatrices <- args[3]
GenomicFeaturesCallableBED <- args[4]
CallableRegions <- args[5]
PDF <- args[6]
strchrs <- args[7]
Rconfig <- args[8]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#script <- "./scripts/Plotting_DNA/plot_PerSiteStats.R"
#MetadataFile <- "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt" 
#strchrFilterStats <- "results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr1.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr2.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr3.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr4.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr5.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr6.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr7.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr8.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr9.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr10.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr11.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr12.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr13.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr14.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr15.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr16.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr17.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr18.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr19.tab.gz"
#strSiteMatrices <- "results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr1.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr2.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr3.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr4.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr5.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr6.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr7.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr8.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr9.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr10.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr11.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr12.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr13.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr14.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr15.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr16.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr17.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr18.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr19.tab.gz"
#GenomicFeaturesCallableBED <- "results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_ExtraCallableRegions.bed"
#CallableRegions <- "results/VariantCalling_DNA/7_join_extra_callable_regions/NonExtraCallableRegions.bed.gz"
#PDF <- "results/Plotting_DNA/plot_PerSiteStats/plot_PerSiteStats.pdf"
#strchrs <- "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
#Rconfig <- "config/AmphiHetDupExp_plot.R"
#####################





SiteMatrices <- unlist(strsplit(strSiteMatrices, " "))
chrFilterStats <- unlist(strsplit(strchrFilterStats, " "))
chrs <- unlist(strsplit(strchrs, " "))
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
heatmap.colors <- c("lightgoldenrodyellow", "gold", "orange", "red")
colfunc <- colorRampPalette(chr.colors)
chr.rampcolors <- colfunc(length(chrs))


######################################################################
# Read data

### Metadata
Metadata <- read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Metadata <- unique(Metadata[,c("Sample","Population","Sex")])
samples <- unique(Metadata$Sample)
pcolors <- population.colors[match(Metadata$Population, unique(Metadata$Population))]


### Site matrices
header <- colnames(read.table(SiteMatrices[length(chrs)], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F))
header <- unlist(strsplit(header, ";"))
Data <- data.frame(matrix(ncol = length(header)+2, nrow = 0))
colnames(Data) <- header
for(c in c(1:length(chrs))){
	chrData<-read.table(SiteMatrices[c], sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
	Data <- rbind(Data, chrData)
}

### Add columns to Data
Data$MeanCoverage <- rowMeans(Data[,paste0("cov", samples)])
Data$MeanCovFemales <- rowMeans(Data[,paste0("cov", Metadata$Sample[which(Metadata$Sex=="Female")])])
Data$MeanCovMales <- rowMeans(Data[,paste0("cov", Metadata$Sample[which(Metadata$Sex=="Male")])])
Data$AbsDiffCovMaleFemale <- abs(Data$MeanCovFemales-Data$MeanCovMales)
Data$AbsMajorFreq <- unlist(lapply(Data$freq, Extract_AbsMajorFreq))
Data$MajorAllele <- unlist(lapply(Data$freq, Extract_MajorAllele))
Data$EqGenoFemales <- unlist(apply(Data[,paste0("geno", Metadata$Sample[which(Metadata$Sex=="Female")])], 1, Check.EqGenotypes))
Data$EqGenoMales <- unlist(apply(Data[,paste0("geno", Metadata$Sample[which(Metadata$Sex=="Male")])], 1, Check.EqGenotypes))
Data$SampleSingleton <- unlist(apply(Data[,c(grep("AbsMajorFreq", colnames(Data)), grep("geno", colnames(Data)))], 1, Extract_SingletonSample))

head(Metadata)
head(Data)
polymData <- Data[which(Data$polymorhism==1),]
print(length(polymData[,1]))
print(length(Data[,1]))

### Summary for genomic feature type
WGlength <- as.numeric(system(paste("zcat", CallableRegions, "| awk -F'\\t' '{l=l+$3-$2+1;}END{print l}'"), intern=TRUE))
print(WGlength)
Sumlengths <- c()
NumVar <- c()
MeanHet <- c()
MeanPI <- c()
MeanFreq <- c()
NumPolym <- c()
MeanPopFst <- c()
MeanPIPolym <- c()
MeanFreqPolym <- c()
for(feat in types.of.features){
	system_out <- system(paste0("cat ", GenomicFeaturesCallableBED, " | grep '", feat, "' | sort -k1,1 -k2,2V | bedtools merge -i stdin | awk -F'\\t' '{l+=$3-$2+1;}END{print l}'"), intern=TRUE)
	Sumlengths <- c(Sumlengths,as.numeric(system_out))
	NumVar <- c(NumVar,length(grep(feat, Data$GenomicFeatureType)))
	MeanHet <- c(MeanHet,mean(Data$oHet[grep(feat, Data$GenomicFeatureType)]))
	MeanPI <- c(MeanPI,mean(Data$PI[grep(feat, Data$GenomicFeatureType)]))
	MeanFreq <- c(MeanFreq,mean(na.omit(Data$AbsMajorFreq[grep(feat, Data$GenomicFeatureType)]/(length(samples)*2))))
	NumPolym <- c(NumPolym,length(grep(feat, polymData$GenomicFeatureType)))
	MeanPopFst <- c(MeanPopFst,mean(na.omit(Data$pop_FST[grep(feat, Data$GenomicFeatureType)])))
}
Feat.info <- as.data.frame(cbind(types.of.features, as.numeric(Sumlengths), as.numeric(NumVar), as.numeric(MeanHet), as.numeric(MeanPI), as.numeric(MeanFreq), as.numeric(NumPolym), as.numeric(MeanPopFst)))
colnames(Feat.info) <- c("feature", "length","numvar", "meanHet", "meanPI", "meanFreq", "numpolym", "meanPopFst")
Feat.info$proplength <- as.numeric(Feat.info$length)/WGlength*100
Feat.info$propvar <- as.numeric(Feat.info$numvar)/length(Data[,1])*100
Feat.info$proppolym <- as.numeric(Feat.info$numpolym)/length(polymData[,1])*100
print(Feat.info)



HeadFilterStats <- read.table(text=system(paste0("zcat ", strchrFilterStats, " | head -1"), intern = TRUE), sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
FilterStats <- read.table(text=system(paste0("zcat ", strchrFilterStats, " | grep -v 'CHROM'"), intern = TRUE), sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
colnames(FilterStats) <- colnames(HeadFilterStats)
head(FilterStats)
FilterStats <- FilterStats[which(FilterStats$FILTER=="PASS"),]
Data <- cbind(Data, FilterStats)

# hetmap sample comparison
SComparison <- matrix(rep(NA,length(samples)*length(samples)), nrow = length(samples), ncol = length(samples))
for(s1 in c(1:length(samples))){
	for(s2 in c(1:length(samples))){
		if(s1!=s2){
			SComparison[s1,s2] <- length(Data[which(Data[,paste0("geno", samples[s1])]==Data[,paste0("geno", samples[s2])]),1])/length(Data[,1])
		}
	}
}
rownames(SComparison) <- samples
colnames(SComparison) <- samples
# Agglomeration method: '"ward.D"', '"ward.D2"', '"single"', '"complete"', '"average"' (= UPGMA), '"mcquitty"' (= WPGMA), '"median"' (= WPGMC) or '"centroid"' (= UPGMC).
hclu<-hclust(dist(SComparison), method="average")
SComparison.sorted <- SComparison[samples[hclu$order],samples[hclu$order]]



save(Data,file="results/Plotting_DNA/plot_PerSiteStats/Data.RData")
save(SComparison,file="results/Plotting_DNA/plot_PerSiteStats/SComparison.RData")

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

PC1lim <- c(-8,8)
PC1step <- (PC1lim[2]-PC1lim[1])/40
PC2lim <- c(-5,5)
PC2step <- (PC2lim[2]-PC2lim[1])/40
oHetlim <- c(0,36)
PIlim <- c(0,1)

layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter_heatmap("", Data$oHet, "Number of heterozygotes", Data$MinCoverage, "Minimum coverage", heatmap.colors, oHetlim, c(0,20), 1, 1)
plot_scatter_heatmap("", Data$PI, "Average pairwise differences", Data$MinCoverage, "Minimum coverage", heatmap.colors, PIlim, c(0,20), 1/20, 1)
plot_scatter_heatmap("", Data$PI, "Average pairwise differences", Data$AbsMajorFreq, "Absolute major frequency", heatmap.colors, PIlim, c(0,75), 1/20, 1)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$oHet, "Number of heterozygotes", heatmap.colors, PC1lim, oHetlim, PC1step, 1)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$PI, "Average pairwise differences", heatmap.colors, PC1lim, PIlim, PC1step, 1/50)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$AbsMajorFreq, "Absolute major frequency", heatmap.colors, PC1lim, c(0,75), PC1step, 1)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$MinCoverage, "Minimum coverage", heatmap.colors, PC1lim, c(5,25), PC1step, 1)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$MeanCoverage, "Mean coverage", heatmap.colors, PC1lim, c(0,30), PC1step, 1)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$pop_FST, "pop_FST", heatmap.colors, PC1lim, c(-.1,1), PC1step, 1/10)
plot_scatter_heatmap("", Data$PC1all, "PC1", Data$sex_FST, "sex_FST", heatmap.colors, PC1lim, c(-.1,1), PC1step, 1/10)
plot_scatter_heatmap("", Data$PC2all, "PC2", Data$pop_FST, "pop_FST", heatmap.colors, PC2lim, c(-.1,1), PC2step, 1/10)
plot_scatter_heatmap("", Data$PC2all, "PC2", Data$sex_FST, "sex_FST", heatmap.colors, PC2lim, c(-.1,1), PC2step, 1/10)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_scatter_heatmap("", Data$AbsMajorFreq, "Absolute major frequency", Data$QUAL, "QUAL", heatmap.colors, c(0,75), c(0,1000), 1, 50)
plot_scatter_heatmap("", Data$AbsMajorFreq, "Absolute major frequency", Data$DP, "DP", heatmap.colors, c(0,75), c(0,1000), 1, 50)
plot_scatter_heatmap("", Data$AbsMajorFreq, "Absolute major frequency", Data$numalleles, "Number of alleles", heatmap.colors, c(0,75), c(2,10), 1, 1)


layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_RelativeFreqDiscreteValues_histogram(Data$numalleles, length(Data$numalleles), seq(2,10,1), "Number of alleles", "Relative frequency", "", c(0,1), "black")

par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_histogram_absoluteSFS(Data$AbsMajorFreq, seq(72,12,-1), "Absolute major allele frequency", "All variants", c(0,.5), ppal.color)
plot_histogram_absoluteSFS(Data$AbsMajorFreq[which(Data$numalleles==2)], seq(72,12,-1), "Absolute major allele frequency", "Biallelic variants", c(0,.6), ppal.color)
plot_histogram_absoluteSFS(Data$AbsMajorFreq[which(Data$numalleles>2)], seq(72,12,-1), "Absolute major allele frequency", "Multiallelic variants (>2 alleles)", c(0,.5), ppal.color)
#plot_histogram_absoluteSFS(Data$AbsMajorFreq[grep("intergenic",Data$GenomicFeatureType)], seq(72,12,-1), "Absolute major allele frequency", "Intergenic variants", c(0,.5), ppal.color)
#plot_histogram_absoluteSFS(Data$AbsMajorFreq[grep("promoter",Data$GenomicFeatureType)], seq(72,12,-1), "Absolute major allele frequency", "Promoter variants", c(0,.5), ppal.color)
#plot_histogram_absoluteSFS(Data$AbsMajorFreq[grep("exon",Data$GenomicFeatureType)], seq(72,12,-1), "Absolute major allele frequency", "Exonic variants", c(0,.5), ppal.color)
#plot_histogram_absoluteSFS(Data$AbsMajorFreq[grep("intron",Data$GenomicFeatureType)], seq(72,12,-1), "Absolute major allele frequency", "Intronic variants", c(0,.5), ppal.color)

#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_joined_histogram_SFS(Data$AbsMajorFreq, Data$GenomicFeatureType, types.of.features, seq(1,72,1), "Absolute major allele frequency", NULL, c(0,.8), types.of.features.color)





hweData <- Data[which(!is.na(Data$P_hwe)),]
propDef <- length(hweData[which(hweData$P_HET_DEFICIT_hwe<=pvalthresh),1])/length(hweData[,1])
print(paste("propDef:", propDef))
propExc <- length(hweData[which(hweData$P_HET_EXCESS_hwe<=pvalthresh),1])/length(hweData[,1])
print(paste("propExc:", propExc))

par(mar=c(5,5,2,2),oma=c(1,1,1,1), xaxs="i", yaxs="i")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(5,1), heights=c(5,5), TRUE)

heatmap(SComparison, scale = "none", 
	labRow=samples, 
	labCol=samples, 
	RowSideColors = pcolors,
	ColSideColors = pcolors,
	hclustfun = function(d) hclust(d, method="average"))

#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_boxplot_PerFeat(Data$oHet/length(samples), Data$GenomicFeatureType, types.of.features, "Poportion of heterozygotes", "All variants", types.of.features.color, c(0,.5), pvalthresh, OutFolder)
#plot_boxplot_PerFeat(Data$PI, Data$GenomicFeatureType, types.of.features, "Average pairwise differences (pi)", "All variants", types.of.features.color, c(0,.5), pvalthresh, OutFolder)
#plot_boxplot_PerFeat(Data$AbsMajorFreq/(length(samples)*2), Data$GenomicFeatureType, types.of.features, "Major allele frequency", "All variants", types.of.features.color, c(.5,1), pvalthresh, OutFolder)
#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_boxplot_PerFeat(polymData$oHet/length(samples), polymData$GenomicFeatureType, types.of.features, "Poportion of heterozygotes", "Only polymorphism", types.of.features.color, c(0,.7), pvalthresh, OutFolder)
#plot_boxplot_PerFeat(polymData$PI, polymData$GenomicFeatureType, types.of.features, "Average pairwise differences (pi)", "Only polymorphism", types.of.features.color, c(0,.8), pvalthresh, OutFolder)
#plot_boxplot_PerFeat(polymData$AbsMajorFreq/(length(samples)*2), polymData$GenomicFeatureType, types.of.features, "Major allele frequency", "Only polymorphism", types.of.features.color, c(.5,1), pvalthresh, OutFolder)

#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_boxplot_PerFeat_SND(Data$oHet/length(samples), Data$GenomicFeatureType, Data$SynNonSynDeg, types.of.features, types.of.CodingSites, "Poportion of heterozygotes", "All variants", types.of.features.color, c(0,.5), pvalthresh, OutFolder)
#plot_boxplot_PerFeat_SND(Data$PI, Data$GenomicFeatureType, Data$SynNonSynDeg, types.of.features, types.of.CodingSites, "Average pairwise differences (pi)", "All variants", types.of.features.color, c(0,.5), pvalthresh, OutFolder)
#plot_boxplot_PerFeat_SND(1-Data$AbsMajorFreq/(length(samples)*2), Data$GenomicFeatureType, Data$SynNonSynDeg, types.of.features, types.of.CodingSites, "Minor allele frequency", "All variants", types.of.features.color, c(0,.3), pvalthresh, OutFolder)
#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_boxplot_PerFeat_SND(polymData$oHet/length(samples), polymData$GenomicFeatureType, polymData$SynNonSynDeg, types.of.features, types.of.CodingSites, "Poportion of heterozygotes", "Only polymorphism", types.of.features.color, c(0,.7), pvalthresh, OutFolder)
#plot_boxplot_PerFeat_SND(polymData$PI, polymData$GenomicFeatureType, polymData$SynNonSynDeg, types.of.features, types.of.CodingSites, "Average pairwise differences (pi)", "Only polymorphism", types.of.features.color, c(0,.8), pvalthresh, OutFolder)
#plot_boxplot_PerFeat_SND(1-polymData$AbsMajorFreq/(length(samples)*2), polymData$GenomicFeatureType, polymData$SynNonSynDeg, types.of.features, types.of.CodingSites, "Minor allele frequency", "Only polymorphism", types.of.features.color, c(0,.3), pvalthresh, OutFolder)

#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_propVar_propLength_PerFeat(Feat.info, types.of.features, "propvar", "proplength", "All variants", c(0,60))
#plot_propVar_propLength_PerFeat(Feat.info, types.of.features, "proppolym", "proplength", "Only polymorphism", c(0,60))

#layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
#plot_value_PerFeat(Feat.info, types.of.features, "meanHet", "Mean heterozygosity", "All variants", c(0,5), types.of.features.color)
#plot_value_PerFeat(Feat.info, types.of.features, "meanPI", "Mean PI", "All variants", c(0,.2), types.of.features.color)
#plot_value_PerFeat(Feat.info, types.of.features, "meanFreq", "Mean maximum frequency", "All variants", c(.9,1), types.of.features.color)
#plot_value_PerFeat(Feat.info, types.of.features, "meanPopFst", "meanPopFst", "All variants", c(.9,1), types.of.features.color)



#par(mar=c(1,5,1,2),oma=c(1,1,1,1))
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$sex_FST, "sex Fst", c(-.1,1), chrs, chr.rampcolors)
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$pop_FST, "pop Fst", c(-.1,1), chrs, chr.rampcolors)
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$MeanCovFemales, "MeanCovFemales", c(5,20), chrs, chr.rampcolors)
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$MeanCovMales, "MeanCovMales", c(5,20), chrs, chr.rampcolors)
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$AbsDiffCovMaleFemale, "AbsDiffCovMaleFemale", c(0,150), chrs, chr.rampcolors)
#AlongSequence_AllChr(Data[, c("chr","st","end"),], Data$PC1all, "PC1", c(-8,8), chrs, chr.rampcolors)


dev.off()






