#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strObsCumPropfiles <- args[1]
strSim1CumPropfiles <- args[2]
strSim2CumPropfiles <- args[3]
PDF <- args[4]
Rconfig <- args[5]
script <- sub(".*=", "", commandArgs()[4])



#####################
# debug / develop
#####################

#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
ObsCumPropfiles <- unlist(strsplit(strObsCumPropfiles, " "))
Sim1CumPropfiles <- unlist(strsplit(strSim1CumPropfiles, " "))
Sim2CumPropfiles <- unlist(strsplit(strSim2CumPropfiles, " "))


######################################################################
# Functions

plot_lines_CumProp <- function(odf, s1df, s2df, xlab, ylab, ylim, col, legendlist, ErrBarsOrNot){
	xaxmin = min(c(as.numeric(colnames(odf)),as.numeric(colnames(s1df)),as.numeric(colnames(s2df))))
	xaxmax = max(c(as.numeric(colnames(odf)),as.numeric(colnames(s1df)),as.numeric(colnames(s2df))))
	xaxisval = seq(xaxmin,xaxmax,1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(xaxmin, xaxmax), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	if(ErrBarsOrNot){
		draw_lines_wErrBars(odf, col)
		draw_lines_wErrBars(s1df, "gold")
		draw_lines_wErrBars(s2df, "forestgreen")
	}else{
		draw_lines(odf, col)
		draw_lines(s1df, "gold")
		draw_lines(s2df, "forestgreen")
	}
	axis(1, at = xaxisval, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	legend("bottomright", legendlist, pch=21, lwd=NA, text.col="black", col=c(col, "gold", "forestgreen"), bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

draw_lines <- function(df, c){
	class(df)
	xval <- as.numeric(colnames(df))
	for(i in c(1:length(df[,1]))){
		lines(xval, df[i,], lwd=0.5, col=c)
		points(xval, df[i,], pch=23, cex=1, col="black", bg=c)
	}
}

draw_lines_wErrBars <- function(df, c){
	class(df)
	xval <- as.numeric(colnames(df))
	means <- colMeans(df)
	sds <- apply(df, 2, sd)
	arrows(x0=xval, y0=means-sds, x1=xval, y1=means+sds, code=3, angle=90, length=0.05, col="black", lwd=1)
	lines(xval, means, lwd=0.5, col=c)
	points(xval, means, pch=23, cex=1, col="black", bg=c)
}

read_set_of_CumPropfiles <- function(files){
	data <- read.table(files[1], sep=" ", header=FALSE, check.names = F, stringsAsFactors = F)
	for(f in c(2:length(files))){
		fdata <- read.table(files[f], sep=" ", header=FALSE, check.names = F, stringsAsFactors = F)
		data <- rbind(data, fdata)
	}
	colnames(data) <- seq(1,length(colnames(data)),1)
	return(data)
}

######################################################################
# Read data

# Read different bootstrapps
ObsData <- read_set_of_CumPropfiles(ObsCumPropfiles)
Sim1Data <- read_set_of_CumPropfiles(Sim1CumPropfiles)
Sim2Data <- read_set_of_CumPropfiles(Sim2CumPropfiles)
print(ObsData)
print(Sim1Data)
print(Sim2Data)

legendList <- c("Observed data")
info <- unlist(strsplit(Sim1CumPropfiles[1], "/"))[5]
legendList <- c(legendList,  paste0("Ne = ", unlist(strsplit(info, "_"))[4], ", mu = ",unlist(strsplit(info, "_"))[3]))
info <- unlist(strsplit(Sim2CumPropfiles[1], "/"))[5]
legendList <- c(legendList,  paste0("Ne = ", unlist(strsplit(info, "_"))[4], ", mu = ",unlist(strsplit(info, "_"))[3]))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(3), heights=c(1), TRUE)
plot_lines_CumProp(ObsData, Sim1Data, Sim2Data, "Number of individuals", "% of explained variance", c(0,1), ppal.color, legendList, FALSE)
plot_lines_CumProp(ObsData, Sim1Data, Sim2Data, "Number of individuals", "% of explained variance", c(0,1), ppal.color, legendList, TRUE)

dev.off()















