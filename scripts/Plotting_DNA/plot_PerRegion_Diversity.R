#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")


######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPiTObs <- args[1]
strPiTObsExons <- args[2]
strPiTObsIntrons <- args[3]
strPiTObsPromoters <- args[4]
strPiTObsIntergenic <- args[5]
strPiPerRObsExons <- args[6]
strPiPerRObsIntrons <- args[7]
strPiPerRObsPromoters <- args[8]
strPiPerRObsIntergenic <- args[9]
PDF <- args[10]
Rconfig <- args[11]
script <- sub(".*=", "", commandArgs()[4])



#####################
# debug / develop
#####################

#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PiTObs <- unlist(strsplit(strPiTObs, " "))
PiTObsExons <- unlist(strsplit(strPiTObsExons, " "))
PiTObsIntrons <- unlist(strsplit(strPiTObsIntrons, " "))
PiTObsPromoters <- unlist(strsplit(strPiTObsPromoters, " "))
PiTObsIntergenic <- unlist(strsplit(strPiTObsIntergenic, " "))
PiPerRObsExons <- unlist(strsplit(strPiPerRObsExons, " "))
PiPerRObsIntrons <- unlist(strsplit(strPiPerRObsIntrons, " "))
PiPerRObsPromoters <- unlist(strsplit(strPiPerRObsPromoters, " "))
PiPerRObsIntergenic <- unlist(strsplit(strPiPerRObsIntergenic, " "))

print(PiPerRObsExons)
######################################################################
# Functions

plot_violin_AllPi <- function(TE, TIr, TP, TIe, RE, RIr, RP, RIe, ylim){
	width <- 0.3
	popdist <- 0.16
	scalefactor <- 1000000
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1-.5,4+.5), col=NA)
	mtext("Average pairwise differences ( )", side = 2, line = 3, cex=1.2)
	violin(RE[,3], 1-popdist, width, population.colors[1])
	points(1-popdist, TE[,3], pch=18, col=population.colors[1], cex=2)
	violin(RE[,2], 1+popdist, width, population.colors[2])
	points(1+popdist, TE[,2], pch=18, col=population.colors[2], cex=2)

	violin(RIr[,3], 2-popdist, width, population.colors[1])
	points(2-popdist, TIr[,3], pch=18, col=population.colors[1], cex=2)
	violin(RIr[,2], 2+popdist, width, population.colors[2])
	points(2+popdist, TIr[,2], pch=18, col=population.colors[2], cex=2)

	violin(RP[,3], 3-popdist, width, population.colors[1])
	points(3-popdist, TP[,3], pch=18, col=population.colors[1], cex=2)
	violin(RP[,2], 3+popdist, width, population.colors[2])
	points(3+popdist, TP[,2], pch=18, col=population.colors[2], cex=2)

	violin(RIe[,3], 4-popdist, width, population.colors[1])
	points(4-popdist, TIe[,3], pch=18, col=population.colors[1], cex=2)
	violin(RIe[,2], 4+popdist, width, population.colors[2])
	points(4+popdist, TIe[,2], pch=18, col=population.colors[2], cex=2)
	axis(1, at = c(1:4), labels=c("Exons","Introns","Promoters","Intergenic"), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

violin <- function(values, x, w, col){
	d <- density(values, adjust = 2)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(x+ynorm, rev(x-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(col), border=col, lwd=1)
	#points(jitter(rep(x, length(values)), amount=w/2), values, pch=21, bg=modif_alpha(col,0.2), col=modif_alpha(col))
}

plot_boxplot_AllPi <- function(TE, TIr, TP, TIe, RE, RIr, RP, RIe, ylim){
	width <- 0.25
	popdist <- 0.16
	scalefactor <- 1000000
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1-.5,4+.5), col=NA)
	mtext("Average pairwise differences ( )", side = 2, line = 3, cex=1.2)
	DrawBox(RE[,3], 1-popdist, population.colors[1], width)
	points(1-popdist, TE[,3], pch=18, col=population.colors[1], cex=2)
	DrawBox(RE[,2], 1+popdist, population.colors[2], width)
	points(1+popdist, TE[,2], pch=18, col=population.colors[2], cex=2)

	DrawBox(RIr[,3], 2-popdist, population.colors[1], width)
	points(2-popdist, TIr[,3], pch=18, col=population.colors[1], cex=2)
	DrawBox(RIr[,2], 2+popdist, population.colors[2], width)
	points(2+popdist, TIr[,2], pch=18, col=population.colors[2], cex=2)

	DrawBox(RP[,3], 3-popdist, population.colors[1], width)
	points(3-popdist, TP[,3], pch=18, col=population.colors[1], cex=2)
	DrawBox(RP[,2], 3+popdist, population.colors[2], width)
	points(3+popdist, TP[,2], pch=18, col=population.colors[2], cex=2)

	DrawBox(RIe[,3], 4-popdist, population.colors[1], width)
	points(4-popdist, TIe[,3], pch=18, col=population.colors[1], cex=2)
	DrawBox(RIe[,2], 4+popdist, population.colors[2], width)
	points(4+popdist, TIe[,2], pch=18, col=population.colors[2], cex=2)
	axis(1, at = c(1:4), labels=c("Exons","Introns","Promoters","Intergenic"), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


DrawBox <- function(values, pos, col, w=.8, den=NULL, text=FALSE, cextext=1){
	s <- boxplot(values, plot=FALSE)
	#points(rep(pos,length(s$out)), s$out, col=modif_alpha("black", 0.3), pch=16, cex=1.5, xpd = NA)
	arrows(x0=pos, y0=s$stats[1], x1=pos, y1=s$stats[5], angle=90, code=3, length=w/10, lwd=2, xpd = NA)
	#lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modif_alpha(col), border=modif_alpha(col), lwd=3)		
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=3, density=den)
	}else{
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modif_alpha(col), border=col, lwd=3, density=den)		
	}
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border="black", lwd=3)		
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}


read_pi <- function(files){
	print(files[1])
	Data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	for(f in c(2:length(files))){
		fData <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
		Data <- cbind(Data, fData[,1])
	}
	colnames(Data) <- c("piAll","piAtl","piMed")
	print(length(Data[which(Data$piAll==0 & Data$piAtl==0 & Data$piMed==0),1]))
	Data <- Data[which(!(Data$piAll==0 & Data$piAtl==0 & Data$piMed==0)),]
	print(head(Data))
	return(Data)
}


######################################################################
# Read data

PiT <- read_pi(PiTObs)
PiTExons <- read_pi(PiTObsExons)
PiTIntrons <- read_pi(PiTObsIntrons)
PiTPromoters <- read_pi(PiTObsPromoters)
PiTIntergenic <- read_pi(PiTObsIntergenic)
PiRExons <- read_pi(PiPerRObsExons)
PiRIntrons <- read_pi(PiPerRObsIntrons)
PiRPromoters <- read_pi(PiPerRObsPromoters)
PiRIntergenic <- read_pi(PiPerRObsIntergenic)

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))




layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot_violin_AllPi(PiTExons, PiTIntrons, PiTPromoters, PiTIntergenic, PiRExons, PiRIntrons, PiRPromoters, PiRIntergenic, c(0,.2))


plot_boxplot_AllPi(PiTExons, PiTIntrons, PiTPromoters, PiTIntergenic, PiRExons, PiRIntrons, PiRPromoters, PiRIntergenic, c(0,.2))



dev.off()















