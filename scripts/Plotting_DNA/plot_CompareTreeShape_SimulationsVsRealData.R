#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
SimRatios <- args[1]
AtlRealRatios <- args[2]
MedRealRatios <- args[3]
PDF <- args[4]
TimesRandomTest <- args[5]
Rconfig <- args[6]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#####################

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)



######################################################################
# Functions

plot_ratioIntExt_Ne_mu_real <- function(main, rvec, sdf, colNe, ylim){
	w <- 0.5
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	print(c(0,length(Nevalues)+2))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,length(Nevalues)+2), col=NA)
	mtext("", side = 1, line = 4, cex=1.2)
	mtext("Internal / External branch length", side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	#points(jitter(rep(1,length(rvec)), amount=0.25), rvec, pch=19, cex=1, col="black")
	d <- density(rvec, adjust = 1)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1+ynorm, rev(1-ynorm)), c(d$x,rev(d$x)), col=modif_alpha("black"), border="black", lwd=1)
	points(1, mean(rvec), pch=19, cex=1, col="black")
	points(1, median(rvec), pch=21, cex=1, col="black")
	for(nm in c(1:length(uNemu))){
		vec <- c()
		for(rep in unique(sdf$Rep)){
			#points(jitter(1+nm, amount=0.25), sdf$ratio[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep)], pch=19, cex=1, col=colNe[nm])
			vec <- c(vec, sdf$ratio[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep)])
		}
		d <- density(vec, adjust = 1)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(1+nm+ynorm, rev(1+nm-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(colNe[nm]), border=colNe[nm], lwd=1)
		points(1+nm, mean(vec), pch=19, cex=1, col=colNe[nm])
		points(1+nm, median(vec), pch=21, cex=1, col=colNe[nm])
	}
	#axis(1, at = seq(1,length(Nevalues)+1), labels = NA, lwd.ticks=1, las=1, cex.axis=1)
	#axis(1, at = seq(1,length(Nevalues)+1), labels = c("Real", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1.5, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Ne_mu_old <- function(sdf, colNe){
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	maxNe <- max(sdf$Ne)
	maxmu <- max(sdf$mu)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0,length(Nevalues)+2), col=NA)
	mtext("Ne", side = 2, line = 4, cex=1.2)
	mtext("mu", side = 4, line = 4, cex=1.2)
	Neline <- c()
	muline <- c()
	for(nm in c(1:length(uNemu))){
		vecNemu <- as.numeric(as.character(unlist(strsplit(uNemu[nm], " "))))
		points(1+nm, log(vecNemu[1])/log(maxNe)*100, pch=19, cex=1, col="black")
		points(1+nm, (vecNemu[2])/log(maxmu)*100, pch=19, cex=1, col="darkred")
		Neline <- c(Neline, vecNemu[1])
		muline <- c(muline, vecNemu[2])
	}
	lines(c(1:length(uNemu))+1, Neline/maxNe*100, col="black")
	lines(c(1:length(uNemu))+1, muline/maxmu*100, col="darkred")
	axis(1, at = seq(1,length(Nevalues)+1), labels = NA, lwd.ticks=1, las=1, cex.axis=.5)
	axis(1, at = seq(1,length(Nevalues)+1), labels = c("Real", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1, las=1, cex.axis=.5)
	axis(2, at = log(Nevalues)/log(maxNe)*100, labels=Nevalues, lwd.ticks=1, las=1, cex.axis=.5)
	axis(4, at = log(muvalues)/log(maxmu)*100, labels=muvalues, lwd.ticks=1, las=1, cex.axis=.5, col="darkred")
	box()
}

#####################
# MAIN
#####################

######################################################################
# Read Simultions data
SData <- read.table(SimRatios, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
SData$ratiomeans <- SData$internalmean / SData$externalmean
SData$ratiosds <- SData$internalsd / SData$externalsd
SData$ratiomedians <- SData$internalmedian / SData$externalmedian
SData$ratiototals <- SData$internaltotal / SData$externaltotal
print(head(SData))
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(SData$Ne)))

# Read Real data
AltRData <- read.table(AtlRealRatios, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
AltRData$ratiomeans <- AltRData$internalmean / AltRData$externalmean
AltRData$ratiosds <- AltRData$internalsd / AltRData$externalsd
AltRData$ratiomedians <- AltRData$internalmedian / AltRData$externalmedian
AltRData$ratiototals <- RData$internaltotal / AltRData$externaltotal
print(head(AltRData))
MedRData <- read.table(RealRatios, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
MedRData$ratiomeans <- MedRData$internalmean / MedRData$externalmean
MedRData$ratiosds <- MedRData$internalsd / MedRData$externalsd
MedRData$ratiomedians <- MedRData$internalmedian / MedRData$externalmedian
MedRData$ratiototals <- MedRData$internaltotal / MedRData$externaltotal
print(head(MedRData))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("", RData$ratiomeans, -medRData$ratiomeans,  AltRData$ratiomeans, SData, "ratiomeans", colNe, c(0,35))
par(mar=c(7,7,0,7))
plot_Ne_mu(SData, colNe)



dev.off()






