#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
ObsShPriv <- args[1]
strRandShPriv <- args[2]
PDF <- args[3]
Rconfig <- args[4]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
RandShPriv <- unlist(strsplit(strRandShPriv, " "))

print(RandShPriv[1])
print(ObsShPriv)

######################################################################
# Functions

plot_violin_RandomVsObserved <- function(t, xlab, ylab, odf, rdf, comp, ylim, scalefactor, col){
	print(t)
	w <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim/scalefactor, xlim=c(0.5,1.5), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	values <- as.numeric(rdf[which(rdf$Type==t),c(2:length(rdf[1,]))])/scalefactor
	d <- density(values, adjust = 2)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1-ynorm, rev(1+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(col), border=col, lwd=1)
	points(jitter(rep(1, length(rdf[1,])-1), amount=w/2), values, pch=21, bg=modif_alpha(col,0.2), col=modif_alpha(col))
	points(1, odf[which(odf$Type==t),2]/scalefactor, pch=18, col=col, cex=2)
	#xarr <- 1.4
	#arrows(x0=xarr, y0=odf[which(odf$Type==t),2]/scalefactor, x1=xarr, y1=comp$RandMeans[which(comp$Type==t)]/scalefactor, code=3, angle=90, length=0.05, col="black", lwd=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5)/scalefactor, lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_violin_AllRandomVsObserved <- function(tps, odf, rdf, ylim){
	w <- 0.5
	scalefactor <- 1000000
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim/scalefactor, xlim=c(1-.5,length(tps)+.5), col=NA)
	mtext("Millions of variants", side = 2, line = 3, cex=1.2)
	for(i in c(1:length(tps))){
		print(tps[i])
		values <- as.numeric(rdf[which(rdf$Type==tps[i]),c(2:length(rdf[1,]))])/scalefactor
		d <- density(values, adjust = 2)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(i+ynorm, rev(i-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(ppal.color), border=ppal.color, lwd=1)
		points(jitter(rep(i, length(rdf[1,])-1), amount=w/2), values, pch=21, bg=modif_alpha(ppal.color,0.2), col=modif_alpha(ppal.color))
		points(i, odf[which(odf$Type==tps[i]),2]/scalefactor, pch=18, col=ppal.color2, cex=2)
	}
	axis(1, at = c(1:length(tps)), labels=tps, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5)/scalefactor, lwd.ticks=1, las=1, cex.axis=1)
	box()
}

read_set_of_Numberfiles <- function(files){
	data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	data <- data[,c(2,1)]
	if(length(files)>1){
		for(f in c(2:length(files))){
			fdata <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
			fdata <- fdata[,c(2,1)]
			data <- cbind(data, fdata[,2])
		}
	}
	colnames(data) <- c("Type", paste0("Col", seq(1,length(files))))
	return(data)
}

######################################################################
# Read data

# Read different bootstrapps
RandData <- read_set_of_Numberfiles(RandShPriv)
ObsData <- read_set_of_Numberfiles(ObsShPriv)

Types <- c("PrivateA", "PrivateM", "VariantShared", "VariantDifferent", "FixedDifferent")
RandMeans <- c()
for(t in Types){
	print(t)
	if(length(grep(t, RandData$Type)) == 0){
		RandData <- rbind(RandData, c(t, as.numeric(rep(0, length(RandData[1,])-1))))
	}
	if(length(grep(t, ObsData$Type)) == 0){
		ObsData <- rbind(ObsData, c(t, rep(0, length(ObsData[1,])-1)))
	}
	RandMeans <- c(RandMeans, mean(as.numeric(RandData[which(RandData$Type==t),c(2:length(colnames(RandData)))])))
}
ObsData <- ObsData[match(Types, ObsData$Type),]
RandData <- RandData[match(Types, RandData$Type),]

print(RandData)
print(ObsData)
FC <- ObsData$Col1/RandMeans
Comparison <- cbind(ObsData, RandMeans, FC)
print(Comparison)


######################################################################
# Plotting
pdf(PDF, width=10, height=5)
par(mar=c(7,3,2,2),oma=c(1,3,1,1))

layout(matrix(seq(1,length(Types),1),nrow=3,ncol=length(Types),byrow=T), TRUE)

plot_violin_RandomVsObserved("PrivateM", "Private\nMediterranean", "Millions of variants", ObsData, RandData, Comparison, c(9,20)*1000000, 1000000, population.colors[1])
plot_violin_RandomVsObserved("PrivateA", "Private\nAtlantic", "", ObsData, RandData, Comparison, c(9,20)*1000000, 1000000, population.colors[2])
plot_violin_RandomVsObserved("VariantShared", "Shared", "", ObsData, RandData, Comparison, c(9,20)*1000000, 1000000, "forestgreen")
plot_violin_RandomVsObserved("VariantDifferent", "Variant\nDifferent", "", ObsData, RandData, Comparison, c(1,3)*1000000, 1000000, "darkred")
plot_violin_RandomVsObserved("FixedDifferent", "Fixed\nDifferent", "Variants", ObsData, RandData, Comparison, c(0,100), 1, "darkred")


dev.off()















