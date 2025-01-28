#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strObsHetfiles <- args[1]
strSimHetfiles <- args[2]
SamplesOrderInVCF <- args[3]
Metadata <- args[4]
PDF <- args[5]
Rconfig <- args[6]
script <- sub(".*=", "", commandArgs()[4])
print(Rconfig)
source(Rconfig)
quit()
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
ObsHetfiles <- unlist(strsplit(strObsHetfiles, " "))
SimHetfiles <- unlist(strsplit(strSimHetfiles, " "))

print(ObsHetfiles[1])
print(SimHetfiles[1])
######################################################################
# Functions

plot_StandDevsHet_ObsSim <- function(main, mvec, avec, sdf, colNe, ylim){
	w <- 0.5
	adj <- 2
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	print(c(0,length(Nevalues)+2))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("", side = 1, line = 4, cex=1.2)
	mtext("Standard deviation of the heterozygosity", side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	#points(jitter(rep(1,length(mvec)), amount=0.25), mvec, pch=19, cex=1, col="black")
	d <- density(mvec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1+ynorm, rev(1-ynorm)), c(d$x,rev(d$x)), col=modif_alpha("black"), border="black", lwd=1)
	points(1, median(mvec), pch=19, cex=1, col="black")
	#points(jitter(rep(2,length(avec)), amount=0.25), avec, pch=19, cex=1, col="black")
	d <- density(avec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(2+ynorm, rev(2-ynorm)), c(d$x,rev(d$x)), col=modif_alpha("black"), border="black", lwd=1)
	points(2, median(avec), pch=19, cex=1, col="black")
	for(nm in c(1:length(uNemu))){
		vec <- c()
		for(rep in unique(sdf$Rep)){
			points(jitter(2+nm, amount=0.25), sdf$sd[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep)], pch=19, cex=1, col=colNe[nm])
			vec <- c(vec, sdf$sd[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep)])
		}
		d <- density(vec, adjust = adj)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(2+nm+ynorm, rev(2+nm-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(colNe[nm]), border=colNe[nm], lwd=1)
		#points(2+nm, median(vec), pch=19, cex=1, col=colNe[nm])
	}
	#axis(1, at = seq(1,length(Nevalues)+1), labels = NA, lwd.ticks=1, las=1, cex.axis=1)
	#axis(1, at = seq(1,length(Nevalues)+1), labels = c("Real", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1.5, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Ne_mu <- function(sdf, colNe, leglabs){
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	rNe <- c(min(log2(sdf$Ne)),max(log2(sdf$Ne)))
	print(rNe)
	rmu <- c(min(log2(sdf$mu)),max(log2(sdf$mu)))
	print(rmu)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-5,105), xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("Ne", side = 2, line = 4, cex=1.2)
	mtext("mu", side = 4, line = 4, cex=1.2, col="darkred")
	Neline <- c()
	muline <- c()
	for(nm in c(1:length(uNemu))){
		vecNemu <- as.numeric(as.character(unlist(strsplit(uNemu[nm], " "))))
		points(2+nm, (log2(vecNemu[1])-rNe[1])/(rNe[2]-rNe[1])*100, pch=19, cex=1, col="black")
		points(2+nm, (log2(vecNemu[2])-rmu[1])/(rmu[2]-rmu[1])*100, pch=19, cex=1, col="darkred")
		Neline <- c(Neline, vecNemu[1])
		muline <- c(muline, vecNemu[2])
	}
	lines(c(1:length(uNemu))+2, (log2(Neline)-rNe[1])/(rNe[2]-rNe[1])*100, col="black")
	lines(c(1:length(uNemu))+2, (log2(muline)-rmu[1])/(rmu[2]-rmu[1])*100, col="darkred")
	axis(1, at = seq(1,length(Nevalues)+2), labels = NA, lwd.ticks=1, las=1, cex.axis=.5)
	axis(1, at = seq(1,length(Nevalues)+2), labels = c("Mediterranean", "Atlantic", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1, las=1, cex.axis=.5)
	axis(2, at = (log2(Nevalues)-rNe[1])/(rNe[2]-rNe[1])*100, labels=Nevalues, lwd.ticks=1, las=1, cex.axis=.5)
	axis(4, at = (log2(muvalues)-rmu[1])/(rmu[2]-rmu[1])*100, labels=muvalues, lwd.ticks=1, las=1, cex.axis=.5, col="darkred", col.axis="darkred")
	box()
}

read_set_of_Hetfiles <- function(files){
	data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	colnames(data) <- seq(1,length(colnames(data)),1)
	if(length(files)>1){
		for(f in c(2:length(files))){
			fdata <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
			colnames(fdata) <- seq(1,length(colnames(fdata)),1)
			data <- rbind(data, fdata)
		}
	}
	return(data)
}

######################################################################
# Read data
# Metadata
MetData <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
print(head(MetData))
# Observed
ObsData <- read_set_of_Hetfiles(ObsHetfiles)
ObsData <- ObsData[,c(4:length(colnames(ObsData)))]
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ObsData) <- Samples[1,]
ObsDataM <- ObsData[,unique(MetData$Sample[which(MetData$Population=="Banyuls")])]
ObsDataA <- ObsData[,unique(MetData$Sample[which(MetData$Population=="Roscoff")])]
ObsDataM$sd <- apply(ObsDataM, 1, sd)
ObsDataA$sd <- apply(ObsDataA, 1, sd)
print(head(ObsDataM))
print(head(ObsDataA))

# Simulated
SimData <- read_set_of_Hetfiles(SimHetfiles)
SimData$sd <- apply(SimData, 1, sd)
Nes <- c()
mus <- c()
Reps <- c()
for(f in SimHetfiles){
	info <- unlist(strsplit(f, "/"))[5]
	Nes <- c(Nes, unlist(strsplit(info, "_"))[4])
	mus <- c(mus, unlist(strsplit(info, "_"))[3])
	Reps <- c(Reps, unlist(strsplit(unlist(strsplit(f, "/"))[6], "_"))[4])
}
SimData$Ne <- as.numeric(Nes)
SimData$mu <- as.numeric(mus)
SimData$Rep <- as.numeric(Reps)
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(SimData$Ne)))
print(head(SimData))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("", "Standard deviation of the heterozygosity", ObsDataM$sd, ObsDataA$sd, SimData, "sd", colNe, c(0,.1))
par(mar=c(7,7,0,7))
plot_Ne_mu(SimData, colNe)


dev.off()















