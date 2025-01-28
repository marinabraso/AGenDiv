#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPiWObs <- args[1]
wBEDfile <- args[2]
NonExtraCallableRegions <- args[3]
ChrLenfile <- args[4]
PDF <- args[5]
Rconfig <- args[6]
script <- sub(".*=", "", commandArgs()[4])



#####################
# debug / develop
#####################

#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PiWObs <- unlist(strsplit(strPiWObs, " "))


######################################################################
# Functions


read_pi <- function(files){
	print(files[1])
	Data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	for(f in c(2:length(files))){
		fData <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
		Data <- as.data.frame(cbind(Data, fData))
	}
	colnames(Data) <- c("piAll","piAtl","piMed")
	print(head(Data))
	return(Data)
}


plot_Histogram <- function(values, xlab, ylab, ylim, xlim, vabline=NULL){
	step <- 0.05
	h <- hist(values, plot=F, breaks=seq(0,100,step))
	w <- h$mids[2]-h$mids[1]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col="black", border=NA)
	}
	abline(v=vabline)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_all_chr_rows_len <- function(lenDF, regDF, color){
	height <- 0.6
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,length(lenDF$Chr)), xlim=c(0,10), col=NA)
	for(c in c(1:length(lenDF$Chr))){
		print(lenDF$Chr[c])
		pos <- length(lenDF$Chr)-c
		sts <- regDF$st/max(lenDF$Length)*10
		ends <- regDF$end/max(lenDF$Length)*10
		draw_chr_row(height, pos, sts, ends, color)
		polygon(c(0,lenDF$Length[c]/max(lenDF$Length)*10,lenDF$Length[c]/max(lenDF$Length)*10,0), c(pos-height/2,pos-height/2,pos+height/2,pos+height/2), col=NA, border="black")
	}
	axis(2, at = c(1:length(lenDF$Chr))-1, labels=rev(lenDF$Chr), lwd=NA, lwd.ticks=NA, las=1, cex.axis=1)
}


draw_chr_row <- function(h, p, starts, ends, color){
	for(r in c(1:length(starts))){
		polygon(c(starts[r],ends[r],ends[r],starts[r]), c(p-h/2,p-h/2,p+h/2,p+h/2), col=color, border=NA)
	}
}

######################################################################
# Read data

PiW <- read_pi(PiWObs)
print(head(PiW))
wBED <- read.table(text=system(paste0("zcat < ", wBEDfile), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(wBED) <- c("chr", "st", "end")
print(head(wBED))

Data <- cbind(wBED, PiW)
print(head(Data))
print(dim(Data))

ChrLen <- read.table(text=system(paste0("cat ", ChrLenfile, " | sed 's/^>//g' | grep 'chr'"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
ChrLen <- ChrLen[,c(1,2)]
colnames(ChrLen) <- c("Chr", "Length")
print(head(ChrLen))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)

plot_Histogram(Data$piAll*100, "Average pairwise differences (%)", "Relative frequency", c(0,0.05), c(0,10))
plot_Histogram(Data$piAtl*100, "Average pairwise differences (%)", "Relative frequency", c(0,0.05), c(0,10))
plot_Histogram(Data$piMed*100, "Average pairwise differences (%)", "Relative frequency", c(0,0.05), c(0,10))

layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot_all_chr_rows_len(ChrLen, Data, "black")

dev.off()















