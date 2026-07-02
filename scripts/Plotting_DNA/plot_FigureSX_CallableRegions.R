#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(RColorBrewer)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
ExtraCallableRegions <- args[1]
ChrLengths <- args[2]
FunctionalRegionsCallableSpan <- args[3]
FunctionalRegionsTotalSpan <- args[4]
PDF <- args[5]
REPORT <- args[6]
Rconfig <- args[7]
script <- sub(".*=", "", commandArgs()[4])
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))

######################################################################
# Functions


draw_chr_row <- function(h, p, starts, ends, color){
	for(r in c(1:length(starts))){
		polygon(c(starts[r],ends[r],ends[r],starts[r]), c(p-h/2,p-h/2,p+h/2,p+h/2), col=color, border=NA)
	}
}

plot_all_chr_rows_len_bins <- function(lenDF, regDF, ss, es, cs){
	height <- 0.6
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,length(lenDF$Chr)), xlim=c(0,10), col=NA)
	for(c in c(1:length(lenDF$Chr))){
		print(lenDF$Chr[c])
		pos <- length(lenDF$Chr)-c
		for(b in c(1:length(ss))){
			sts <- regDF$st[which(regDF$chr==lenDF$Chr[c] & regDF$len>=ss[b] & regDF$len<es[b])]/max(lenDF$Length)*10
			ends <- regDF$end[which(regDF$chr==lenDF$Chr[c] & regDF$len>=ss[b] & regDF$len<es[b])]/max(lenDF$Length)*10
			draw_chr_row(height, pos, sts, ends, cs[b])
		}
		polygon(c(0,lenDF$Length[c]/max(lenDF$Length)*10,lenDF$Length[c]/max(lenDF$Length)*10,0), c(pos-height/2,pos-height/2,pos+height/2,pos+height/2), col=NA, border="black")
	}
	axis(2, at = c(1:length(lenDF$Chr))-1, labels=rev(lenDF$Chr), lwd=NA, lwd.ticks=NA, las=1, cex.axis=1)
	legend("bottomright", paste(ss, "<= x <", es), pch=19, text.col="black", col=cs, bty = "n", cex=1, xjust = 0, yjust = 0)
}


barplot_chromosomes <- function(chrs, df, type, ss, es, cs, ylab, ylim){
	print(paste("barplot_chromosomes",type))
	w <- 0.8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(1-w/2,length(chrs)+w/2), col=NA)
	mtext(type, side = 3, line = 3, cex=1.5)
	mtext(ylab, side = 2, line = 3, cex=1.5)
	mtext("Chromosomes", side = 1, line = 3, cex=1.5)
	for(c in c(1:length(chrs))){
		print(chrs[c])
		cumst <- 0
		for(b in c(1:length(ss))){
			print(paste(b, df[c,paste0("Perc", type, b)], cs[b]))
			polygon(c(c-w/2,c+w/2,c+w/2,c-w/2), c(cumst,cumst,cumst+df[c,paste0("Perc", type, b)],cumst+df[c,paste0("Perc", type, b)]), col=cs[b], border=NA)
			cumst <- cumst+df[c,paste0("Perc", type, b)]
		}
	}
	axis(1, at = c(1:length(chrs)), labels=chrs, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	legend("topright", paste(ss, "<= x <", es), pch=19, text.col="black", col=cs, bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

plot_density <- function(vec, lab, color, xlim, ylim, thresh=NULL, thresh2=NULL, printlab=TRUE){
	print(lab)
	d <- density(na.omit(vec), adjust = 2)
	print(c(min(d$x), max(d$x)))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA, xaxs="i", yaxs="i")	
	if(printlab==TRUE){
		mtext(lab, side = 1, line = 3, cex=1.5)
		mtext("Density", side = 2, line = 3, cex=1.5)
	}	
	polygon(c(min(d$x)-10,d$x,max(d$x)+10), c(-10,d$y,10), col=color, border=NA)
	lines(d$x, d$y, lwd=2)
	abline(v=thresh, lty=2, lwd=2)
	abline(v=thresh2, lty=2, lwd=2)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_DiffInRelPropCallable_FunctionalRegions <- function(df, ylab, ylim){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, (length(df$Feature)-1)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	for(i in c(1:(length(df$Feature)-1))){
		print(i)
		print(df$Feature[i])
		print(df$Diff[i])
		if(df$Diff[i]>0){
			polygon(c(i-w,i+w,i+w,i-w), c(0,0,df$Diff[i],df$Diff[i]), col="forestgreen", border="forestgreen")
		}else{
			polygon(c(i-w,i+w,i+w,i-w), c(df$Diff[i],df$Diff[i],0,0), col="forestgreen", border="forestgreen")
		}
	}
	abline(h=0, col="black")
	axis(1, at = c(1:(length(df$Feature)-1)), labels=df$Feature[1:(length(df$Feature)-1)], lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/10), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

######################################################################
# Read data
ERegions<-read.table(ExtraCallableRegions, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ERegions) <- c("chr", "st", "end","len")
ChrLen<-read.table(ChrLengths, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ChrLen) <- c("Chr", "Length", "Order", "CummLength")
ChrLen <- ChrLen[grep("chr", ChrLen$Chr),]
ChrLen$Chr <- unlist(lapply(ChrLen$Chr, function(x){substr(x, 2, nchar(x))}))

head(ERegions)
head(ChrLen)


sbins <- c(0,50,100,150,200,250,300)
Eebins <- c(50,100,150,200,250,300,max(ERegions$len))
colfunc <- colorRampPalette(c("brown2", "slateblue3"))
cbins <- colfunc(length(sbins))
print(cbins)

for(b in c(1:length(sbins))){
	vSumECallable <- c()
	for(c in ChrLen$Chr){
		vSumECallable <- c(vSumECallable, sum(ERegions$len[which(ERegions$chr==c & ERegions$len>=sbins[b] & ERegions$len<Eebins[b])]))
	}
	ChrLen[,paste0("SumECallable",b)] <- vSumECallable
	ChrLen[,paste0("PercECallable",b)] <- ChrLen[,paste0("SumECallable",b)]*100/ChrLen$Length
}
head(ChrLen[,paste0("PercECallable", c(1,2,3,4,5,6,7))])


# % of callable length for the total length
CallableSpan <- read.table(FunctionalRegionsCallableSpan, sep="\t", header=FALSE)
FeatSpan <- read.table(FunctionalRegionsTotalSpan, sep="\t", header=FALSE)
FeatSpan <- rbind(FeatSpan, c("Total", sum(ChrLen$Length)))
FeatSpan <- cbind(FeatSpan, CallableSpan[,2])
colnames(FeatSpan) <- c("Feature", "Total", "Callable")
FeatSpan$Total <- as.numeric(FeatSpan$Total)
FeatSpan$Callable <- as.numeric(FeatSpan$Callable)
FeatSpan$PercCall <- FeatSpan$Callable*100/FeatSpan$Total
FeatSpan$RelPropTotal <- FeatSpan$Total/FeatSpan$Total[which(FeatSpan$Feature=="Total")]
FeatSpan$RelPropCallable <- FeatSpan$Callable/FeatSpan$Callable[which(FeatSpan$Feature=="Total")]
FeatSpan$Diff <- FeatSpan$RelPropCallable - FeatSpan$RelPropTotal
head(FeatSpan)

######################################################################
# Plotting
pdf(PDF, width=15, height=10)
par(mar=c(0,5,1,1),oma=c(4,4,4,4))
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), TRUE)

#plot_all_chr_rows_len_bins(ChrLen, ERegions, sbins, Eebins, cbins)

par(mar=c(5,5,2,2),oma=c(1,1,1,1))
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)

barplot_chromosomes(unlist(lapply(ChrLen$Chr, function(x){substr(x, 4, nchar(x))})), ChrLen, "ECallable", sbins, Eebins, cbins, "% of length", c(0,60))

plot_density(ERegions$len, "Extra callable regions length", "forestgreen", c(0,500), c(0,.02), sbins)

plot_DiffInRelPropCallable_FunctionalRegions(FeatSpan, "% of callable length", c(-1,1))






dev.off()

