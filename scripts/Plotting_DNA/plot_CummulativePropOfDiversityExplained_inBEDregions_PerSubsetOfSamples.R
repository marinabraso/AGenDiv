#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strCumPropfiles <- args[1]
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
CumPropfiles <- unlist(strsplit(strCumPropfiles, " "))

######################################################################
# Functions

plot_lines_CumProp <- function(df, xlab, ylab, ylim, col){
	xaxisval <- as.numeric(colnames(df))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(min(xaxisval), max(xaxisval)), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	for(i in c(1:length(df[,1]))){
		lines(xaxisval, df[i,], lwd=0.5, col=col)
		points(xaxisval, df[i,], pch=23, cex=1.5, col="black", bg=col)
	}
	axis(1, at = seq(1,max(xaxisval),5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}





######################################################################
# Read data

# Read different bootstrapps
Data <- read.table(CumPropfiles[1], sep=" ", header=FALSE, check.names = F, stringsAsFactors = F)
for(f in c(2:length(CumPropfiles))){
	fData <- read.table(CumPropfiles[f], sep=" ", header=FALSE, check.names = F, stringsAsFactors = F)
	Data <- rbind(Data, fData)
}
colnames(Data) <- seq(length(colnames(Data)),1,-1)
print(Data)


######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_lines_CumProp(Data, "Number of individuals", "% of explained variance", c(0,1), ppal.color)

dev.off()















