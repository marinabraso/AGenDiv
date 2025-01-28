#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
blfile <- args[1]
PDF <- args[2]
Rconfig <- args[3]
samplesize <- args[4]
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

plot_ratioIntExt_Ne_mu <- function(main, df, colNe, ylim){
	colnames(df) <- c("Ne", "mu", "Rep", "ratio")
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,length(unique(df$Ne))+1), col=NA)
	mtext("", side = 1, line = 4, cex=1.2)
	mtext("Internal / External branch length", side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	uNemu <- unique(paste(df$Ne, df$mu, sep=" "))
	print(uNemu)
	for(nm in c(1:length(uNemu))){
		print(uNemu[nm])
		for(rep in unique(df$Rep)){
			points(nm, df$ratio[which(paste(df$Ne, df$mu, sep=" ")==uNemu[nm] & df$Rep==rep)], pch=19, cex=1, col="black")
		}
	}
	axis(1, at = seq(1,length(unique(df$Ne))), labels = NA, lwd.ticks=1, las=1, cex.axis=1)
	axis(1, at = seq(1,length(unique(df$Ne))), labels = gsub(" ", "\n", uNemu), lwd.ticks=NA, lwd=NA, line=2, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

#####################
# MAIN
#####################

######################################################################
# Read data
Data <- read.table(blfile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Data$ratiomean <- Data$internalmean / Data$externalmean
Data$ratiomedian <- Data$internalmedian / Data$externalmedian
Data$ratiototal <- Data$internaltotal / Data$externaltotal
print(head(Data))
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(Data$Ne)))
print(unique(Data$Ne))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot_ratioIntExt_Ne_mu("ratiomean", Data[,c("Ne", "mu", "Rep", "ratiomean")], colNe, c(0,20))
plot_ratioIntExt_Ne_mu("ratiomedian", Data[,c("Ne", "mu", "Rep", "ratiomedian")], colNe, c(0,20))
plot_ratioIntExt_Ne_mu("ratiototal", Data[,c("Ne", "mu", "Rep", "ratiototal")], colNe, c(0,20))



dev.off()






