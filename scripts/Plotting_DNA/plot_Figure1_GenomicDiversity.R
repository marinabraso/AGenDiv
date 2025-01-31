#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
######################################################################
# Libraries & functions
library(RColorBrewer)
library(png)

######################################################################
# Files & folders
args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPiPerRegFiles <- args[1]
strPiTotalFiles <- args[2]
strHetPerRegFiles <- args[3]
strHetTotalFiles <- args[4]
strFreqPerSiteFiles <- args[5]
SamplesOrderInVCF <- args[6]
MapImage <- args[7]
PDF <- args[8]
REPORT <- args[9]
strAtlSamples <- args[10]
strMedSamples <- args[11]
Rconfig <- args[12]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PiPerRegFiles <- unlist(strsplit(strPiPerRegFiles, " "))
PiTotalFiles <- unlist(strsplit(strPiTotalFiles, " "))
HetPerRegFiles <- unlist(strsplit(strHetPerRegFiles, " "))
HetTotalFiles <- unlist(strsplit(strHetTotalFiles, " "))
FreqPerSiteFiles <- unlist(strsplit(strFreqPerSiteFiles, " "))
AtlSamples <- unlist(strsplit(strAtlSamples, " "))
MedSamples <- unlist(strsplit(strMedSamples, " "))

print(AtlSamples)
print(MedSamples)

######################################################################
# Functions

read_SingleValueFiles <- function(files){
	Data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    if(length(files)>1){
        for(f in c(2:length(files))){
            Data <- c(Data, read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F))
        }
    }
	return(Data)
}

plot_HetPerPop <- function(Hetdf, atls, meds, pivalues, ylab, zlab, xlabs, ylim){
	w <- 0.7
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5, 2.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	#mtext(zlab, side = 4, line = 3, cex=1.5)
    violin(Hetdf[1,atls], 1, w, NA, modif_alpha(population.colors[1]), 2)
	points(jitter(rep(1, length(atls)), amount=w/3), Hetdf[,atls], pch=21, bg=modif_alpha(population.colors[1],0.2), col=modif_alpha(population.colors[1]), cex=2)
    violin(Hetdf[1,meds], 2, w, NA, modif_alpha(population.colors[2]), 2)
	points(jitter(rep(2, length(meds)), amount=w/3), Hetdf[,meds], pch=21, bg=modif_alpha(population.colors[2],0.2), col=modif_alpha(population.colors[2]), cex=2)
	#points(1, pivalues[2], pch=18, col=population.colors[1], cex=2)
	#points(2, pivalues[3], pch=18, col=population.colors[2], cex=2)
	#abline(h=pivalues[1], lty=2)
	#axis(4, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1,2), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1,2), labels=xlabs, lwd=NA, las=1, line=2, cex.axis=2)
	box()
}

plot_PiPerRegionPerPop <- function(RegFiles, TotalFiles, ylab, xlabs, ylim){
	w <- 0.2
    distPop <- 0.3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, 10), col=NA, xaxs = "i", yaxs = "i")
	mtext(ylab, side = 2, line = 4, cex=1.5)
    # Exon
    PiRegAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Exons*100), col = types.of.features.color[1], lwd=2, lty = 1)
    # Intron
    PiRegAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
   lines(ecdf(PiRegAll.Introns*100), col = types.of.features.color[2], lwd=2, lty = 1)
    # Promoter
    PiRegAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Promoters*100), col = types.of.features.color[3], lwd=2, lty = 1)
    # Intergenic
    PiRegAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Intergenic*100), col = types.of.features.color[4], lwd=2, lty = 1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1:length(xlabs)), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1:length(xlabs)), labels=xlabs, lwd=NA, las=1, line=2, cex.axis=2)
	box()
}

######################################################################
# Read data & plotting
VariantTypes = c("Observed_Data", "Observed_SNPs", "Observed_INDELs")
RegionTypes = c("Exons", "Introns", "Promoters", "Intergenic")
SampleGoups = c("AllSamples", "AtlSamples", "MedSamples")
AllSamples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)


## Start plot
pdf(PDF, width=10, height=15)
par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

### A
## Map of sampling
par(mar=c(3,3,2,2))
map <- readPNG(paste0(system("pwd", intern=TRUE), "/", MapImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(map, 1, 1, 10, 10)
legend("topleft", populations, pch=19, col=population.colors, bty = "n", cex=2, xjust = 0, yjust = 0)
writePlotLabel("A")

### B
## Plot heterozygosity per sample
par(mar=c(7,7,2,2))
HetAll <- read.table(grep("Observed_Data.*Callable", HetTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetAll) <- AllSamples
PiTAll <- read.table(grep("Observed_Data.*Callable.*AllSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTAtl <- read.table(grep("Observed_Data.*Callable.*AtlSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTMed <- read.table(grep("Observed_Data.*Callable.*MedSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
plot_HetPerPop(HetAll*100, AtlSamples, MedSamples, c(PiTAll,PiTAtl,PiTMed)*100, "Heterozygosity (%)", "Average pairwise differences (%)", populations, c(0,5))
writePlotLabel("B")
## Write heterozygosity values to report
write("### Heterozygosity per sample", file = REPORT, append = FALSE)
write(paste(AllSamples,collapse="\t"), file = REPORT, append = TRUE)
write(paste(HetAll,collapse="\t"), file = REPORT, append = TRUE)
## Write pi values to report
write("### Average pairwise differences total", file = REPORT, append = TRUE)
write("All\tAtlantic\tMediterranean", file = REPORT, append = TRUE)
write(paste(PiTAll, PiTAtl, PiTMed, sep="\t"), file = REPORT, append = TRUE)

### C
## Plots of pi on different regions (Exons, Introns, Promoters, Intergenic)
plot_PiPerRegionPerPop(PiPerRegFiles, PiTotalFiles, "Average pairwise differences (%)", RegionTypes, c(0,1))







dev.off()















