#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)
library(ade4)
library(ggtree)
library(png)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPCAValues_files <- args[1]
strPCAPropVar_files <- args[2]
strObsShPriv <- args[3]
strRandShPriv <- args[4]
jpsmc <- args[5]
jpsmc_thetas <- args[6]
SamplesOrderInVCF <- args[7]
Metadata <- args[8]
DistTree <- args[9]
PDF <- args[10]
REPORT <- args[11]
strchrs <- args[12]
PSMChighmu <- as.numeric(args[13])
PSMClowmu <- as.numeric(args[14])
Blangenerationtime <- as.numeric(args[15])
Rconfig <- args[16]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
TreeFolder <- system(paste0("echo $(dirname ", DistTree, ")"), intern=TRUE)
PCAValues_files <- unlist(strsplit(strPCAValues_files, " "))
PCAPropVar_files <- unlist(strsplit(strPCAPropVar_files, " "))
ObsShPriv <- unlist(strsplit(strObsShPriv, " "))
RandShPriv <- unlist(strsplit(strRandShPriv, " "))
chrs <- unlist(strsplit(strchrs, " "))

######################################################################
# Functions

plot_PCA <- function(vec1, lab1, vec2, lab2, colors, samples, legcolors, leglabels, legpos){
	vec1 <- as.numeric(vec1)
	vec2 <- as.numeric(vec2)
	ylim <- c(min(vec2), max(vec2))
	xlim <- c(min(vec1), max(vec1))
	yflank <- (ylim[2]-ylim[1])*0.05
	xflank <- (xlim[2]-xlim[1])*0.05
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-yflank,ylim[2]+yflank), xlim=c(xlim[1]-xflank,xlim[2]+xflank), col=NA)
	mtext(lab1, side = 1, line = 4, cex=1.5)
	mtext(lab2, side = 2, line = 4, cex=1.5)
	points(vec1, vec2, pch=21, col=colors, bg=modif_alpha(colors), cex=4)
	text(vec1, vec2, labels=samples, font=1, cex=.7)
	if(length(leglabels)>1 & !is.na(leglabels[1])){
		legend("right", as.character(leglabels), pch=15, text.col="black", col=legcolors, bty = "n", cex=2, xjust = 0, yjust = 0)
	}
	#axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=2, cex.axis=1.5)
	#axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_violin_RandomVsObserved <- function(t, xlab, ylab, odf, rdf, ylim, scalefactor, col){
	print(t)
	print(odf)
	print(rdf[,c(1:5)])
	w <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim/scalefactor, xlim=c(0.5,1.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	mtext(xlab, side = 1, line = 4, cex=1.5)
	values <- as.numeric(rdf[t,c(1:length(rdf[1,]))])/scalefactor
	d <- density(values, adjust = 2)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1-ynorm, rev(1+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(col), border=col, lwd=1)
	points(jitter(rep(1, length(rdf[1,])), amount=w/2), values, pch=21, bg=modif_alpha(col,0.2), col=modif_alpha(col))
	points(1, odf[t,1]/scalefactor, pch=18, col=col, cex=2)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5)/scalefactor, lwd.ticks=1, las=1, cex.axis=1)
	box()
}


plot_violin_SeveralRandomVsObserved <- function(tps, odf, rdf, ylim, col, scalefactor, ylab, xlabels){
	w <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim/scalefactor, xlim=c(1-.5,length(tps)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	for(i in c(1:length(tps))){
		print(tps[i])
		values <- as.numeric(rdf[tps[i],c(1:length(rdf[1,]))])/scalefactor
		d <- density(values, adjust = 2)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(i+ynorm, rev(i-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(col[i]), border=col[i], lwd=1)
		#points(jitter(rep(i, length(rdf[1,])), amount=w/2), values, pch=21, bg=modif_alpha(col[i],0.2), col=modif_alpha(col[i]))
		points(i, odf[tps[i],1]/scalefactor, pch=19, col=col[i], cex=2)
	}
	axis(1, at = c(1:length(tps)), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(1, at = c(1:length(tps)), labels=xlabels, lwd=NA, line = 3, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/11)/scalefactor, lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_PSMC_samples_wo_mu <- function(df, samples, murange, g, colors){
	Xcolumn1 <- "T_k_lowmu"
	Xcolumn2 <- "T_k_highmu"
	Ycolumn1 <- "N_k_lowmu"
	Ycolumn2 <- "N_k_highmu"
	yscaling <- 100000
	ylim <- c(0, max(df[,Ycolumn1])/yscaling)
	print(ylim)
	xlim <- log(c(round(min(df[which(df[,Xcolumn1]!=0),Xcolumn1])/1000, digits=0)*1000, round(max(df[,Xcolumn1])/100000, digits=0)*100000))
	print(c(round(min(df[which(df[,Xcolumn1]!=0),Xcolumn1])/1000, digits=0)*1000, round(max(df[,Xcolumn1])/100000, digits=0)*100000))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)		
	mtext(paste0("Years ago (g=", g, ", ", expression(mu), "=", murange[1], "-", murange[2],")"), side = 1, line = 3, cex=1)
	mtext(paste0("Effective population size (x",yscaling,")"), side = 2, line = 3, cex=1)
	for(s in c(1:length(samples))){
		sdf <- df[which(df$sample==samples[s]),]
		x <- sdf[,Xcolumn1]
		y <- sdf[,Ycolumn1]
		stepsx <- c(rbind(x,x))[2:length(c(rbind(x,x)))]
		print(stepsx)
		stepsy <- c(rbind(y,y))[1:length(c(rbind(x,x)))-1]
		lines(log(stepsx), stepsy/yscaling, lwd=2, col=colors[s])
		print(stepsy)
		quit()
	}
	xaxlab1 <- c(10^3,10^4,10^5,10^6,10^7,10^8,10^9)
	xaxlab2 <- xaxlab1*murange[1]/murange[2]
	yaxlab1 <- seq(0,12,2)
	yaxlab2 <- yaxlab1*murange[1]/murange[2]

	axis(1, at = log(xaxlab1), labels=paste0(xaxlab1,"\n",xaxlab2), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = yaxlab1, labels=paste0(yaxlab1,"\n",yaxlab2), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}


######################################################################
# Read data

# Metadata
MetData <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
MetData <- unique(MetData[,c("Sample", "Population", "Sex")])
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
MetData <- MetData[match(Samples, MetData$Sample),]
MetData$ColorP <- population.colors[match(MetData$Population, c("Roscoff", "Banyuls"))]
MetData$ColorS <- sex.colors[match(MetData$Sex, unique(MetData$Sex))]
print(MetData)

# Reading PCA results for several goups of variants
cat("----> Reading PCA results for several goups of variants\n")
PCAresults <- data.frame(matrix(ncol = 0, nrow = length(Samples)))
rownames(PCAresults) <- Samples
PCApropvar <- data.frame(matrix(ncol = 0, nrow = 10))
TypesOfRegions <- c()
for(f in c(1:length(PCAValues_files))){
	TypeRegions <- unlist(strsplit(unlist(strsplit(PCAValues_files[f], "PerSample_in"))[2], "Regions_PerAllSamples"))[1]
	val <- read.table(PCAValues_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
	pvar <- as.data.frame(read.table(PCAPropVar_files[f], h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)[c(1:10),])
	colnames(val) <- paste(TypeRegions, colnames(val), sep="_")
	colnames(pvar) <- TypeRegions
	PCAresults <- cbind(PCAresults, val)
	PCApropvar <- cbind(PCApropvar, pvar)
	TypesOfRegions <- c(TypesOfRegions, TypeRegions)
}
print(head(PCAresults))
print(head(PCApropvar))
print(TypesOfRegions)

# Private and shared variants
OSharedPriv <- read.table(ObsShPriv, h=FALSE, sep = "\t", check.names = F, stringsAsFactors = F)
rownames <- OSharedPriv[,2]
OSharedPriv <- as.data.frame(OSharedPriv[,1])
colnames(OSharedPriv) <- "Observed"
rownames(OSharedPriv) <- rownames
print(OSharedPriv)
RSharedPriv <- data.frame(matrix(ncol = 0, nrow = length(rownames)))
rownames(RSharedPriv) <- rownames
for(r in c(0:(length(RandShPriv)-1))){
	rand <- read.table(RandShPriv[r+1], h=FALSE, sep = "\t", check.names = F, stringsAsFactors = F)
	tmpdf <- as.data.frame(rand$V1[match(rownames(RSharedPriv), rand$V2)])
	colnames(tmpdf) <- paste0("Rand", r)
	RSharedPriv <- cbind(RSharedPriv, tmpdf)
}
RSharedPriv[is.na(RSharedPriv)] <- 0
print(RSharedPriv[,c(1:5)])

# PSMC
PSMCdata <- read.table(jpsmc, h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
print(head(PSMCdata))
PSMCthetas <- read.table(jpsmc_thetas, h=TRUE, sep = "\t", check.names = F, stringsAsFactors = F)
print(head(PSMCthetas))

# If TreeImage verion of the treefile done with FigTree does not exist, create a substitute with ggtree
TreeImage <- paste0(OutFolder, "/SamplesDistanceTree_formated.png")
if(!file.exists(TreeImage)){
	png(filename = TreeImage, width = 1000, height = 1000)
	plot.new()
	tree <- read.tree(DistTree)
	p <- ggtree(tree, layout="daylight", size=1)
	p <- p + geom_tiplab(size=10, color="black", fontface="bold")
	#p <- p + geom_treescale(x = 0, y = 0, width=0.001)
	print(p)
	dev.off()
}

######################################################################
## Start plotting Figure 2
pdf(PDF, width=15, height=10)

par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)

### A
## PCA
par(mar=c(7,7,2,2))
plot_PCA(	PCAresults$Callable_Comp1, paste("PC1", round(PCApropvar$Callable[1], digits = 2), "%"), 
			PCAresults$Callable_Comp2, paste("PC2", round(PCApropvar$Callable[2], digits = 2), "%"), 
			MetData$ColorP, MetData$Sample, population.colors, populations, "topleft")
writePlotLabel("A")

### B
## Distance tree
par(mar=c(3,3,2,2))
tree <- readPNG(paste0(system("pwd", intern=TRUE), "/", TreeImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(tree, 1, 1, 10, 10)
writePlotLabel("B")


### C
## SHared Private variants analysis
par(mar=c(7,7,2,2))
plot_violin_SeveralRandomVsObserved(c("PrivateA", "PrivateM", "VariantShared"), OSharedPriv, RSharedPriv, c(9,20)*1000000, c(population.colors, "forestgreen"), 1000000, "Millions of variants", c("Private\nMediterranean", "Private\nAtlantic", "Shared"))
writePlotLabel("C")

### D
## PSMC
#plot_PSMC_samples_wo_mu(PSMCdata, MetData$Sample, c(PSMChighmu, PSMClowmu), Blangenerationtime, MetData$ColorP)
#plot.new()
#writePlotLabel("D")



dev.off()
#############
## Report
print(OSharedPriv)
## Write Table S3 (Shared Private)
write("### Table S3 (Shared Private)", file = REPORT, append = FALSE)
tableS3 <- OSharedPriv
tableS3$PercentTotal <- tableS3$Observed/sum(tableS3$Observed)*100
tableS3$RMeans <- rowMeans(RSharedPriv)
pvals <- c()
sds <- c()
for(r in c(1:length(RSharedPriv[,1]))){
	print(rownames(tableS3)[r])
	moreextreme <- 0
	print(tableS3$Observed[r])
	print(mean(unlist(RSharedPriv[r,])))
	if(tableS3$Observed[r] > mean(unlist(RSharedPriv[r,]))){
		moreextreme <- length(which(tableS3$Observed[r] <= RSharedPriv[r,]))
	} else {
		moreextreme <- length(which(tableS3$Observed[r] >= RSharedPriv[r,]))
	}
	print(moreextreme)
	p_val <- moreextreme/length(RSharedPriv[r,])
	if(p_val==0){
		p_val <- paste0("<", 1/length(RSharedPriv[r,]))
	}
	pvals <- c(pvals, p_val)
	sd_val <- sd(as.numeric(RSharedPriv[r,]))
	sds <- c(sds, sd_val)
}
tableS3$RSds <- sds
tableS3$Rp_val <- pvals
print(tableS3)
write.table(tableS3, file = REPORT, append = TRUE, row.names=TRUE, sep="\t", quote = FALSE)










