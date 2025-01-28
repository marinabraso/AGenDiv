#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
ChrsFilterStats <- args[1]
SNPableRegions <- args[2]
CallableRegions <- args[3]
ExtraCallableRegions <- args[4]
ChrsCoverageStats <- args[5]
ChrLengths <- args[6]
outputtext <- args[7]
PDF <- args[8]
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
system(paste0("mkdir -p $(dirname ", PDF, ")"))


#####################################
## To debug
#ChrsFilterStats <- "results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr18.tab.gz results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.chr19.tab.gz"
#SNPableRegions <- "results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz"
#CallableRegions <- "results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz"
#ExtraCallableRegions <- "results/VariantCalling_DNA/7_join_extra_callable_regions/ExtraCallableRegions.bed.gz"
#ChrsCoverageStats <- "results/VariantAnalysis_DNA/ProcessingCoverageData/CoverageWindowStats_chr18.bed.gz results/VariantAnalysis_DNA/ProcessingCoverageData/CoverageWindowStats_chr19.bed.gz"
#ChrLengths <- "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#outputtext <- "results/Plotting_DNA/plot_GeneralVariantFilterStats/GeneralVariantFilterStats.txt"
#PDF <- "results/Plotting_DNA/plot_GeneralVariantFilterStats/GeneralVariantFilterStats.pdf"
#source("scripts/Plotting_DNA/plot_GeneralVariantFilterStats_functions.R")
#####################################

QD_threshold <- 2.0
QUAL_thresholdNovaSeq6000 <- 24.0
QUAL_thresholdHiSeq4000 <- 30.0
SOR_threshold <- 3.0
FS_threshold <- 60.0
MQ_threshold <- 40.0
MQRankSum_threshold <- -12.5
ReadPosRankSum_threshold <- -8.0
MinCov_threshold <- 5

######################################################################
# Read data
SRegions<-read.table(SNPableRegions, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(SRegions) <- c("chr", "st", "end","len")
CRegions<-read.table(CallableRegions, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(CRegions) <- c("chr", "st", "end","len")
ERegions<-read.table(ExtraCallableRegions, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ERegions) <- c("chr", "st", "end","len")
ChrLen<-read.table(ChrLengths, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ChrLen) <- c("Chr", "Length", "Order", "CummLength")
ChrLen <- ChrLen[grep("chr", ChrLen$Chr),]
ChrLen$Chr <- unlist(lapply(ChrLen$Chr, function(x){substr(x, 2, nchar(x))}))

head(SRegions)
head(CRegions)
head(ERegions)
head(ChrLen)

HeadFilterStats <- read.table(text=system(paste0("zcat ", ChrsFilterStats, " | head -1"), intern = TRUE), sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
FilterStats <- read.table(text=system(paste0("zcat ", ChrsFilterStats, " | grep -v 'CHROM'"), intern = TRUE), sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
colnames(FilterStats) <- colnames(HeadFilterStats)
head(FilterStats)
Filters <- unique(unlist(strsplit(paste0(unique(FilterStats$FILTER),collapse=","), ",")))
Types <- c("All", "SNP", "INDEL")
FiltersSummary <- data.frame(matrix(ncol = length(Types), nrow = 0))
colnames(FiltersSummary) <- Types

for(f in c(1:length(Filters))){
	frow <- length(grep(Filters[f], FilterStats$FILTER))
	frow <- c(frow, length(grep(Filters[f], FilterStats$FILTER[which(FilterStats$TYPE=="SNP")])))
	frow <- c(frow, length(grep(Filters[f], FilterStats$FILTER[which(FilterStats$TYPE=="INDEL")])))
	FiltersSummary <- rbind(FiltersSummary, frow)
}
frow <- length(FilterStats$FILTER)
frow <- c(frow, length(FilterStats$FILTER[which(FilterStats$TYPE=="SNP")]))
frow <- c(frow, length(FilterStats$FILTER[which(FilterStats$TYPE=="INDEL")]))
FiltersSummary <- rbind(FiltersSummary, frow)
colnames(FiltersSummary) <- Types
rownames(FiltersSummary) <- c(Filters, "Total")
head(FiltersSummary)
write.table(FiltersSummary, file=outputtext, sep = "\t", col.names = TRUE, row.names = TRUE)

CovWind <- read.table(text=system(paste0("zcat ", ChrsCoverageStats), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(CovWind) <- c("chr", "start", "end", "avMinCov", "avMaxCov", "avAvCov", "avSdCov")
head(CovWind)


sbins <- c(0,50,100,150,200,250,300)
Sebins <- c(50,100,150,200,250,300,max(SRegions$len))
Cebins <- c(50,100,150,200,250,300,max(CRegions$len))
Eebins <- c(50,100,150,200,250,300,max(ERegions$len))
colfunc <- colorRampPalette(c("brown2", "slateblue3"))
cbins <- colfunc(length(sbins))
print(cbins)

for(b in c(1:length(sbins))){
	vSumSNPable <- c()
	vSumCallable <- c()
	vSumECallable <- c()
	for(c in ChrLen$Chr){
		vSumSNPable <- c(vSumSNPable, sum(SRegions$len[which(SRegions$chr==c & SRegions$len>=sbins[b] & SRegions$len<Sebins[b])]))
		vSumCallable <- c(vSumCallable, sum(CRegions$len[which(CRegions$chr==c & CRegions$len>=sbins[b] & CRegions$len<Cebins[b])]))
		vSumECallable <- c(vSumECallable, sum(ERegions$len[which(ERegions$chr==c & ERegions$len>=sbins[b] & ERegions$len<Eebins[b])]))
	}
	ChrLen[,paste0("SumSNPable",b)] <- vSumSNPable
	ChrLen[,paste0("PercSNPable",b)] <- ChrLen[,paste0("SumSNPable",b)]*100/ChrLen$Length
	ChrLen[,paste0("SumCallable",b)] <- vSumCallable
	ChrLen[,paste0("PercCallable",b)] <- ChrLen[,paste0("SumCallable",b)]*100/ChrLen$Length
	ChrLen[,paste0("SumECallable",b)] <- vSumECallable
	ChrLen[,paste0("PercECallable",b)] <- ChrLen[,paste0("SumECallable",b)]*100/ChrLen$Length
}
head(ChrLen[,paste0("PercSNPable", c(1,2,3,4,5,6,7))])
head(ChrLen[,paste0("PercCallable", c(1,2,3,4,5,6,7))])
head(ChrLen[,paste0("PercECallable", c(1,2,3,4,5,6,7))])

######################################################################
# Plotting
pdf(PDF, width=15, height=10)
par(mar=c(0,5,1,1),oma=c(4,4,4,4))
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), TRUE)

plot_all_chr_rows_len_bins(ChrLen, CRegions, sbins, Cebins, cbins)
plot_all_chr_rows_len_bins(ChrLen, SRegions, sbins, Sebins, cbins)
plot_all_chr_rows_len_bins(ChrLen, ERegions, sbins, Eebins, cbins)

par(mar=c(5,5,2,2),oma=c(1,1,1,1))
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)

barplot_chromosomes(unlist(lapply(ChrLen$Chr, function(x){substr(x, 4, nchar(x))})), ChrLen, "SNPable", sbins, Sebins, cbins, "% of length", c(0,100))
barplot_chromosomes(unlist(lapply(ChrLen$Chr, function(x){substr(x, 4, nchar(x))})), ChrLen, "Callable", sbins, Cebins, cbins, "% of length", c(0,60))
barplot_chromosomes(unlist(lapply(ChrLen$Chr, function(x){substr(x, 4, nchar(x))})), ChrLen, "Extra Callable", sbins, Eebins, cbins, "% of length", c(0,60))

#To develop
#barplot_filters(FiltersSummary)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_density(SRegions$len, "SNPable regions length", "forestgreen", c(0,500), c(0,.02), sbins)
plot_density(CRegions$len, "Callable regions length", "forestgreen", c(0,500), c(0,.02), sbins)
plot_density(ERegions$len, "Extra callable regions length", "forestgreen", c(0,500), c(0,.02), sbins)







dev.off()








