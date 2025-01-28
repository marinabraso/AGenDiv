#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")


######################################################################
# Libraries & functions
library(RColorBrewer)
library(ash)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
Freqfile <- args[1]
PDF <- args[2]
Rconfig <- args[3]
script <- sub(".*=", "", commandArgs()[4])



#####################
# debug / develop
#####################

#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)

######################################################################
# Functions

plot_SFS <- function(values, mids, xlab, main, ylim, col){
	w <- mids[2]-mids[1]
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	print(main)
	print(h$counts)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext("Proportion of sites", side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	## Draw expectation of SFS
	#eh <- h$counts[1]
	#for(i in c(2:length(h$counts))){
	#	eh <- c(eh, h$counts[1]/i)
	#}
	#for(i in c(1:length(h$mids))){
	#	polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,eh[i],eh[i]), col=modif_color(col, .3), border=NA)
	#}
	# Draw real values
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	abline(v=36)
	axis(1, at = seq(h$mids[1],h$mids[length(h$mids)],5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
	#return(h)
}



plot_RelativeFreq_Histogram <- function(values, mids, xlab, ylab, main, ylim, col){
	w <- mids[2]-mids[1]
	total <- length(values)
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	h$counts <- h$counts/total
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	axis(1, at = h$mids, labels=, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}





######################################################################
# Read data

MajorFreq <- read.table(text=system(paste0("zcat < ", Freqfile, " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
NumAlleles <- read.table(text=system(paste0("zcat < ", Freqfile, " | awk '{print NF}'"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
Data <- cbind(MajorFreq, NumAlleles)
colnames(Data) <- c("MajorFreq", "NumAlleles")
print(head(Data))


######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

# SFS
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_SFS(Data$MajorFreq, seq(72,12,-1), "Absolute major allele frequency", "All variants", c(0,.5), ppal.color)
plot_SFS(Data$MajorFreq[which(Data$NumAlleles==2)], seq(72,12,-1), "Absolute major allele frequency", "Biallelic variants", c(0,.6), ppal.color)
plot_SFS(Data$MajorFreq[which(Data$NumAlleles>2)], seq(72,12,-1), "Absolute major allele frequency", "Multiallelic variants (>2 alleles)", c(0,.5), ppal.color)
# SNPs vs INDELS
# Exons, introns...

# Histogram_NumAlleles
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_RelativeFreq_Histogram(Data$NumAlleles, seq(2,10,1), "Number of alleles", "Relative frequency", "", c(0,1), "black")
# SNPs vs INDELS
# Exons, introns...


dev.off()















