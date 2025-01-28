#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
HetAll <- args[1]
strPiTObs <- args[2]
SamplesOrder <- args[3]
PDF <- args[4]
strAtlSamples <- args[5]
strMedSamples <- args[6]
Rconfig <- args[7]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PiTObs <- unlist(strsplit(strPiTObs, " "))
AtlSamples <- unlist(strsplit(strAtlSamples, " "))
MedSamples <- unlist(strsplit(strMedSamples, " "))

print(AtlSamples)
print(MedSamples)
######################################################################
# Functions
plot_violinHetPerPop_pointsPI <- function(Hetdf, atls, meds, pivalues, ylab, zlab, xlabs, ylim){
	print("buu")
	w <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5, 2.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	#mtext(zlab, side = 4, line = 3, cex=1.2)

	md <- density(as.numeric(Hetdf[1,meds]), adjust = 2)
	mynorm <- md$y/max(md$y)*w/2
	polygon(c(1+mynorm, rev(1-mynorm)), c(md$x,rev(md$x)), col=modif_alpha(population.colors[1]), border=population.colors[1], lwd=1)
	points(jitter(rep(1, length(meds)), amount=w/2), Hetdf[,meds], pch=21, bg=modif_alpha(population.colors[1],0.2), col=modif_alpha(population.colors[1]))
	ad <- density(as.numeric(Hetdf[1,atls]), adjust = 2)
	aynorm <- ad$y/max(ad$y)*w/2
	polygon(c(2+aynorm, rev(2-aynorm)), c(ad$x,rev(ad$x)), col=modif_alpha(population.colors[2]), border=population.colors[2], lwd=1)
	points(jitter(rep(2, length(atls)), amount=w/2), Hetdf[,atls], pch=21, bg=modif_alpha(population.colors[2],0.2), col=modif_alpha(population.colors[2]))


	#points(1, pivalues[2], pch=18, col=population.colors[1], cex=2)
	#points(2, pivalues[3], pch=18, col=population.colors[2], cex=2)
	#abline(h=pivalues[1], lty=2)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(1, at = c(1,2), labels=xlabs, lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_HetBar <- function(Hetdf, atls, meds, ylab, xlabs, ylim){
	w <- 0.8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, length(Hetdf[1,])), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	colors <- c(rep(1,length(meds)), rep(2,length(meds)))
	for(s in c(1:length(Hetdf[1,]))){
		polygon(c(s-w/2, s+w/2, s+w/2, s-w/2), c(0,0,Hetdf[1,s],Hetdf[1,s]), col=population.colors[colors[s]], border=NA, lwd=1)
	}
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(1, at = c(length(meds)/2, length(meds)+length(atls)/2), labels=xlabs, lwd.ticks=0, las=1, cex.axis=1.5)
	box()
}

######################################################################
# Read data

HetAll <- read.table(HetAll, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
SamplesOrder <- read.table(SamplesOrder, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetAll) <- SamplesOrder[1,]
print(head(HetAll))

PiTAll <- read.table(PiTObs[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
print(head(PiTAll))
PiTAtl <- read.table(PiTObs[2], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
print(head(PiTAtl))
PiTMed <- read.table(PiTObs[3], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
print(head(PiTMed))


######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

plot_violinHetPerPop_pointsPI(HetAll*100, AtlSamples, MedSamples, c(PiTAll,PiTAtl,PiTMed)*100, "Heterozygosity (%)", "Average pairwise differences", c("Mediterranean", "Atlantic"), c(2,3))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot_HetBar(HetAll*100, AtlSamples, MedSamples, "Heterozygosity (%)", c("Mediterranean", "Atlantic"), c(0,10))
t <- t.test(as.numeric(HetAll[1,AtlSamples]), y = as.numeric(HetAll[1,MedSamples]), alternative = "two.sided")
print(t)


dev.off()















