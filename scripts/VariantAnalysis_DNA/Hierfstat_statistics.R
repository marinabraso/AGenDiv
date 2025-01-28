#!/usr/bin/env Rscript

Sys.setenv(LANGUAGE='en')

######################################################################
# Libraries & functions
library(hierfstat)
library(gaston)

######################################################################
# Files & folders
args <- commandArgs(trailingOnly=TRUE)
Metadata <- args[1]
SampleOrder <- args[2]
DataFile <- args[3]
OutFile <- args[4]
script <- sub(".*=", "", commandArgs()[4])


#####################
# develop / debug
#OutFile <- "results/VariantAnalysis_DNA/tmp_overallFst_hierfstat/Hierfstat_OverallStats.tbl"
#DataFile <- "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.bed.gz"
#SampleOrder <- "metadata/SamplesOrderInVCF.chr19.txt"
#Metadata <- "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#####################

ResultsFolder <- dirname(OutFile)
system(paste0("mkdir -p ", ResultsFolder))

#####################
cat("read data\n")
# Read the big genotype file
data <- read.table(DataFile, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)

# Read the sample order
colSamples <- read.table(SampleOrder, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)

# read and order metadata according to sample order
metadata <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
metadata <- unique(metadata[,c("Sample","Population","Sex")])
rownames(metadata) <- metadata$Sample
metadata <- metadata[unlist(colSamples),]
print(head(metadata))
print(dim(metadata))

# Transpose data & add population column
data <- as.data.frame(t(data))
print(head(data[,1:8]))
print(dim(data))
data <- cbind(metadata$Population, data)
print(head(data[,1:8]))
print(dim(data))


#####################
cat("hierfstat basic.stats\n")
## basic.stats hierfstat (see hierfstats manual basic.stats section)
# Ho: observed heterozygosity
# Hs: mean of all within population gene diversity
# Ht: overall gene diversity
# Dst: Ht − Hs
# Htp: Hs + Dstp
# Dstp: np/(np − 1)Dst
# Fst: Dst/Ht
# Fstp: Dstp/Htp
# Fis: 1 − Ho/Hs
# Dest: np/(np − 1)(Htp − Hs)/(1 − Hs)
bstats <- basic.stats(data, digits=6)

## per population and per locus observed heterozygosity
cat("write.table Ho\n")
write.table(bstats$Ho, file = paste(ResultsFolder, "/Hierfstat_ObsHeterozygosity_perlocus_perpop.tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
## per population and per locus expected heterozygosity
cat("write.table Hs\n")
write.table(bstats$Hs, file = paste(ResultsFolder, "/Hierfstat_GeneDiversity_perlocus_perpop.tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
## per pop and per locus Fis
cat("write.table Fis\n")
write.table(bstats$Fis, file = paste(ResultsFolder, "/Hierfstat_Fis_perlocus_perpop.tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
## per locus summary stats across all pops 
cat("write.table perloc\n")
write.table(bstats$perloc, file = paste(ResultsFolder, "/Hierfstat_SummaryStats_perlocus.tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
## Overall stats for the whole dataset
# Note: 
	# overall Fst = mean(Dst)/mean(Ht) # loci Fst NA values avoided
	# overall Fis = 1-mean(Ho)/mean(Hs) # idem
cat("write.table overall\n")
write.table(bstats$overall, file = OutFile, quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)

# bootstapping
nboot <- 100
conf <- 0.05 # confidence
ci <- c(conf/2,1-conf/2) # confidence interval
cat("betas\n")
betas <- betas(data, nboot=nboot, lim=ci)
write.table(betas, file = paste(ResultsFolder, "/Hierfstat_betas_", nboot,"_", conf,".tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
cat("boot.ppfis\n")
boot.ppfis <- boot.ppfis(data, nboot=nboot, quant=ci)
write.table(boot.ppfis, file = paste(ResultsFolder, "/Hierfstat_boot.ppfis_", nboot,"_", conf,".tbl", sep =""), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)






