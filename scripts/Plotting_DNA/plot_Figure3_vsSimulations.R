#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
library(RColorBrewer)
library(ade4)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strObsHetfiles <- args[1]
strSimHetfiles <- args[2]
ObsSFSAtl <- args[3]
ObsSFSMed <- args[4]
strSimSFS <- args[5]
ObsNAllelesAtl <- args[6]
ObsNAllelesMed <- args[7]
strSimNAlleles <- args[8]
SamplesOrderInVCF <- args[9]
Metadata <- args[10]
PDF <- args[11]
REPORT <- args[12]
SimulationsSampleSize <- as.numeric(args[13])
Rconfig <- args[14]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
ObsHetfiles <- unlist(strsplit(strObsHetfiles, " "))
SimHetfiles <- unlist(strsplit(strSimHetfiles, " "))
SimSFSfiles <- unlist(strsplit(strSimSFS, " "))
SimNAllelesfiles <- unlist(strsplit(strSimNAlleles, " "))

######################################################################
# Functions

read_set_of_files <- function(files){
	data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	colnames(data) <- seq(1,length(colnames(data)),1)
	if(length(files)>1){
		for(f in c(2:length(files))){
			fdata <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
			colnames(fdata) <- seq(1,length(colnames(fdata)),1)
			data <- rbind(data, fdata)
		}
	}
	return(data)
}

PCA_calc <- function(Data){
	pca <- dudi.pca(Data, center=T, scale=T, scannf=F, nf=5)
	propvar <- 100 * pca$eig/sum(pca$eig)
	co.df <- as.data.frame(pca$co)
	li.df <- as.data.frame(pca$li)
	return(list("pca"=pca, "propvar"=propvar, "co.df"=co.df, "li.df"=li.df))
}


# Plots for Real data vs. Simulations comparison
plot_ObsSim <- function(ylab, mvec, avec, sdf, column, colNe, ylim){
	w <- 0.5
	adj <- 2
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	ylimw <- ylim[2]-ylim[1]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-ylimw/10,ylim[2]+ylimw/10), xlim=c(0,length(Nevalues)+3), col=NA)
	mtext(ylab, side = 2, line = 4, cex=.9)
	# Plot Atlantic windows
	d <- density(avec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1+ynorm, rev(1-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(population.colors[1]), border=population.colors[1], lwd=1)
	points(1, median(avec), pch=19, cex=1, col=population.colors[1])
	# Plot Mediterranean windows
	d <- density(mvec, adjust = adj)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(2+ynorm, rev(2-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(population.colors[2]), border=population.colors[2], lwd=1)
	points(2, median(mvec), pch=19, cex=1, col=population.colors[2])
	# Plot Simulations
	for(nm in c(1:length(uNemu))){
		vec <- c()
		for(rep in unique(sdf$Rep)){
			points(jitter(2+nm, amount=0.25), sdf[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep), column], pch=19, cex=.5, col=colNe[nm])
			vec <- c(vec, sdf[which(paste(sdf$Ne, sdf$mu, sep=" ")==uNemu[nm] & sdf$Rep==rep), column])
		}
		d <- density(vec, adjust = adj)
		ynorm <- d$y/max(d$y)*w/2
		polygon(c(2+nm+ynorm, rev(2+nm-ynorm)), c(d$x,rev(d$x)), col=modif_alpha(colNe[nm]), border=colNe[nm], lwd=1)
	}
	#axis(1, at = seq(1,length(Nevalues)+2), labels = NA, lwd.ticks=1, las=1, cex.axis=.5)
	#axis(1, at = seq(1,length(Nevalues)+2), labels = c("Atlantic", "Mediterranean", gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1, las=1, cex.axis=.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/2), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Ne_mu <- function(sdf, colNe, leglabs){
	uNemu <- unique(paste(sdf$Ne, sdf$mu, sep=" "))
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	rNe <- c(min(log2(sdf$Ne)),max(log2(sdf$Ne)))
	print(rNe)
	rmu <- c(min(log2(sdf$mu)),max(log2(sdf$mu)))
	print(rmu)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-5,105), xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("Ne", side = 2, line = 4, cex=1.2)
	mtext("mu", side = 4, line = 4, cex=1.2, col="darkred")
	Neline <- c()
	muline <- c()
	for(nm in c(1:length(uNemu))){
		vecNemu <- as.numeric(as.character(unlist(strsplit(uNemu[nm], " "))))
		points(2+nm, (log2(vecNemu[1])-rNe[1])/(rNe[2]-rNe[1])*100, pch=19, cex=1, col="black")
		points(2+nm, (log2(vecNemu[2])-rmu[1])/(rmu[2]-rmu[1])*100, pch=19, cex=1, col="darkred")
		Neline <- c(Neline, vecNemu[1])
		muline <- c(muline, vecNemu[2])
	}
	lines(c(1:length(uNemu))+2, (log2(Neline)-rNe[1])/(rNe[2]-rNe[1])*100, col="black")
	lines(c(1:length(uNemu))+2, (log2(muline)-rmu[1])/(rmu[2]-rmu[1])*100, col="darkred")
	axis(1, at = seq(1,length(Nevalues))+2, labels = NA, lwd.ticks=1, las=1, cex.axis=.5)
	axis(1, at = seq(1,length(Nevalues))+2, labels = c(gsub(" ", "\n", uNemu)), lwd.ticks=NA, lwd=NA, line=1, las=1, cex.axis=.5)
	axis(2, at = (log2(Nevalues)-rNe[1])/(rNe[2]-rNe[1])*100, labels=Nevalues, lwd.ticks=1, las=1, cex.axis=.5)
	axis(4, at = (log2(muvalues)-rmu[1])/(rmu[2]-rmu[1])*100, labels=muvalues, lwd.ticks=1, las=1, cex.axis=.5, col="darkred", col.axis="darkred")
	box()
}

plot_mu <- function(sdf, colNe, leglabs){
	umu <- unique(sdf$mu)
	muvalues <- unique(sdf$mu)[order(unique(sdf$mu))]
	rmu <- c(min(log2(sdf$mu)),max(log2(sdf$mu)))
	print(rmu)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-20,120), xlim=c(0,length(muvalues)+3), col=NA)
	mtext("mu", side = 2, line = 4, cex=1)
	for(nm in c(1:length(umu))){
		points(2+nm, (log2(umu[nm])-rmu[1])/(rmu[2]-rmu[1])*100, pch=19, cex=1.5, col=colNe[nm])
	}
	axis(2, at = (log2(muvalues)-rmu[1])/(rmu[2]-rmu[1])*100, labels=muvalues, lwd.ticks=1, las=1, cex.axis=1)
	#axis(1, at = seq(1,length(muvalues))+2, labels = c(gsub(" ", "\n", umu)), lwd.ticks=NA, lwd=NA, las=1, cex.axis=1)
	box()
}

plot_Ne <- function(sdf, colNe, leglabs){
	uNe <- unique(sdf$Ne)
	Nevalues <- unique(sdf$Ne)[order(unique(sdf$Ne))]
	rNe <- c(min(log2(sdf$Ne)),max(log2(sdf$Ne)))
	print(rNe)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-20,120), xlim=c(0,length(Nevalues)+3), col=NA)
	mtext("Ne", side = 2, line = 4, cex=1)
	for(nm in c(1:length(uNe))){
		points(2+nm, (log2(uNe[nm])-rNe[1])/(rNe[2]-rNe[1])*100, pch=19, cex=1.5, col=colNe[nm])
	}
	axis(2, at = (log2(Nevalues)-rNe[1])/(rNe[2]-rNe[1])*100, labels=Nevalues, lwd.ticks=1, las=1, cex.axis=1)
	#axis(1, at = seq(1,length(Nevalues))+2, labels = c(gsub(" ", "\n", uNe)), lwd.ticks=NA, lwd=NA, las=1, cex.axis=1)
	box()
}

plot_PCA <- function(main, vec1, lab1, vec2, lab2, colors, pch, cexpt, samples, legcolors, leglabels, xlim, ylim){
	vec1 <- as.numeric(vec1)
	vec2 <- as.numeric(vec2)
	ylim <- c(min(vec2), max(vec2))
	yflank <- (ylim[2]-ylim[1])*0.03
	xlim <- c(min(vec1), max(vec1))
	xflank <- (xlim[2]-xlim[1])*0.03
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-yflank,ylim[2]+yflank), xlim=c(xlim[1]-xflank,xlim[2]+xflank), col=NA)
	mtext(lab1, side = 1, line = 2, cex=.9)
	mtext(lab2, side = 2, line = 2, cex=.9)
	mtext(main, side = 3, line = 1, cex=2)
	points(vec1, vec2, pch=pch, col=colors, cex=cexpt)
	#text(vec1, vec2, labels=samples, pos=3, font=1, cex=1.5)
	if(length(which(leglabels %in% c("Atlantic","Mediterranean")))>0){
		legend("bottomleft", c("Atlantic","Mediterranean"), pch=c(3,4), text.col="black", col="black", bty = "n", cex=1, xjust = 0, yjust = 0)
	}
	#axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=2, cex.axis=1.5)
	#axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

######################################################################
# Read data
# Metadata
MetData <- read.table(Metadata, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
MetData <- unique(MetData[,c("Sample", "Population")])
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
MetData <- MetData[match(Samples, MetData$Sample),]
print(head(MetData)) 

# Observed Het
ObsData <- read_set_of_files(ObsHetfiles)
ObsData <- ObsData[,c(4:length(colnames(ObsData)))]
Samples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(ObsData) <- Samples[1,]
ObsDataM <- ObsData[,unique(MetData$Sample[which(MetData$Population=="Banyuls")])]
ObsDataA <- ObsData[,unique(MetData$Sample[which(MetData$Population=="Roscoff")])]
ObsDataM$sd <- apply(ObsDataM, 1, sd)
ObsDataA$sd <- apply(ObsDataA, 1, sd)
print(head(ObsDataM))
print(head(ObsDataA))

# Observed SFS
SFSAtl <- read_set_of_files(ObsSFSAtl)
SFSAtl <- SFSAtl[,c(4:length(colnames(SFSAtl)))]
colnames(SFSAtl) <- paste0("n", seq(length(MetData$Sample[which(MetData$Population=="Roscoff")])*2-1, length(MetData$Sample[which(MetData$Population=="Roscoff")])*2-length(colnames(SFSAtl)), -1))
SFSAtl$PropSingletons <- SFSAtl[,grep("n[0-9]", colnames(SFSAtl))[1]]/rowSums(SFSAtl[,grep("n[0-9]", colnames(SFSAtl))])
SFSAtl$PropSingleDoubletons <- rowSums(SFSAtl[,grep("n[0-9]", colnames(SFSAtl))[c(1,2)]])/rowSums(SFSAtl[,grep("n[0-9]", colnames(SFSAtl))])
print(head(SFSAtl))
SFSMed <- read_set_of_files(ObsSFSMed)
SFSMed <- SFSMed[,c(4:length(colnames(SFSMed)))]
colnames(SFSMed) <- paste0("n", seq(length(MetData$Sample[which(MetData$Population=="Banyuls")])*2-1, length(MetData$Sample[which(MetData$Population=="Banyuls")])*2-length(colnames(SFSMed)), -1))
SFSMed$PropSingletons <- SFSMed[,grep("n[0-9]", colnames(SFSMed))[1]]/rowSums(SFSMed[,grep("n[0-9]", colnames(SFSMed))])
SFSMed$PropSingleDoubletons <- rowSums(SFSMed[,grep("n[0-9]", colnames(SFSMed))[c(1,2)]])/rowSums(SFSMed[,grep("n[0-9]", colnames(SFSMed))])
print(head(SFSMed))

# Observed number of alleles
NAllelesAtl <- read_set_of_files(ObsNAllelesAtl)
NAllelesAtl <- NAllelesAtl[,c(4:length(colnames(NAllelesAtl)))]
colnames(NAllelesAtl) <- paste0("n", c(seq(2,2+length(colnames(NAllelesAtl))-2),"+10"))
NAllelesAtl$PropBial <- NAllelesAtl$n2 / rowSums(NAllelesAtl)
print(head(NAllelesAtl))
NAllelesMed <- read_set_of_files(ObsNAllelesMed)
NAllelesMed <- NAllelesMed[,c(4:length(colnames(NAllelesMed)))]
colnames(NAllelesMed) <- paste0("n", c(seq(2,2+length(colnames(NAllelesMed))-2),"+10"))
NAllelesMed$PropBial <- NAllelesMed$n2 / rowSums(NAllelesMed)
print(head(NAllelesMed))

# Simulated Het
SimHet <- read_set_of_files(SimHetfiles)
SimHet$sd <- apply(SimHet, 1, sd)

# Simulated SFS
SimSFS <- read_set_of_files(SimSFSfiles)
SimSFS <- SimSFS[,c(4:length(colnames(SimSFS)))]
colnames(SimSFS) <- paste0("n", seq(SimulationsSampleSize*2-1, SimulationsSampleSize*2-length(colnames(SimSFS)), -1))
SimSFS$PropSingletons <- SimSFS[,grep("n[0-9]", colnames(SimSFS))[1]]/rowSums(SimSFS[,grep("n[0-9]", colnames(SimSFS))])
SimSFS$PropSingleDoubletons <- rowSums(SimSFS[,grep("n[0-9]", colnames(SimSFS))[c(1,2)]])/rowSums(SimSFS[,grep("n[0-9]", colnames(SimSFS))])

# Simulated number of alleles
SimNAlleles <- read_set_of_files(SimNAllelesfiles)
SimNAlleles <- SimNAlleles[,c(4:length(colnames(SimNAlleles)))]
colnames(SimNAlleles) <- paste0("n", c(seq(2,2+length(colnames(SimNAlleles))-2),"+10"))
SimNAlleles$PropBial <- SimNAlleles$n2 / rowSums(SimNAlleles)

# Add Ne mu & Rep columns
Nes <- c()
mus <- c()
Reps <- c()
for(f in SimSFSfiles){
	info <- unlist(strsplit(f, "/"))[5]
	Nes <- c(Nes, unlist(strsplit(info, "_"))[4])
	mus <- c(mus, unlist(strsplit(info, "_"))[3])
	Reps <- c(Reps, unlist(strsplit(unlist(strsplit(f, "/"))[6], "_"))[2])
}
SimHet$Ne <- as.numeric(Nes)
SimHet$mu <- as.numeric(mus)
SimHet$Rep <- as.numeric(Reps)
SimSFS$Ne <- as.numeric(Nes)
SimSFS$mu <- as.numeric(mus)
SimSFS$Rep <- as.numeric(Reps)
SimNAlleles$Ne <- as.numeric(Nes)
SimNAlleles$mu <- as.numeric(mus)
SimNAlleles$Rep <- as.numeric(Reps)
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(SimSFS$Ne)))
Nemu.combinations <- unique(paste0(format(SimSFS$Ne, scientific = FALSE), "_", format(SimSFS$mu, scientific = FALSE)))
print(Nemu.combinations)
maxcolsSFS <- max(length(SimSFS[1,]), length(SFSAtl[1,]), length(SFSMed[1,]))
colfunc <- colorRampPalette(c("darkred", "gold"))
colsSFS <- colfunc(maxcolsSFS)

NeVal <- 		c(1500000, 					1000000, 				750000, 				  500000, 				  100000, 				  75000, 					50000, 					10000, 					7500, 				  5000)
NeValForm <- 	c(expression("1.5·10" ^ 6), expression("1·10" ^ 6), expression("7.5·10" ^ 5), expression("5·10" ^ 5), expression("10" ^ 5), expression("7.5·10" ^ 4), expression("5·10" ^ 4), expression("10" ^ 4), expression("7.5·10" ^ 3), expression("5·10" ^ 3))
muVal <- c(0.0000000050, 				0.0000000075, 				0.0000000100, 		  0.0000000150, 			 0.0000000750, 				0.0000001000, 		   0.0000001500, 			  0.0000007500, 			 0.0000010000, 			0.0000015000)
muValForm <- c(expression("5·10" ^ -9), expression("7.5·10" ^ -9), expression("10" ^ -8), expression("1.5·10" ^ -8), expression("7.5·10" ^ -8), expression("10" ^ -7), expression("1.5·10" ^ -7), expression("7.5·10" ^ -7), expression("10" ^ -6), expression("1.5·10" ^ -6))

print(head(SimHet))
print(head(SimSFS))
print(head(SimNAlleles))
######################################################################
# PCA analysis of SFS
rSimSFS <- SimSFS[,grep("n[0-9]", colnames(SimSFS))]/rowSums(SimSFS[,grep("n[0-9]", colnames(SimSFS))])
rSFSAtl <- SFSAtl[,grep("n[0-9]", colnames(SFSAtl))]/rowSums(SFSAtl[,grep("n[0-9]", colnames(SFSAtl))])
rSFSMed <- SFSMed[,grep("n[0-9]", colnames(SFSMed))]/rowSums(SFSMed[,grep("n[0-9]", colnames(SFSMed))])
rownames(rSimSFS) <- paste0("Sim_", format(SimSFS$Ne, scientific = FALSE), "_", format(SimSFS$mu, scientific = FALSE), "_", SimSFS$Rep)
rownames(rSFSAtl) <- paste0("Atl", seq(1, length(rSFSAtl[,1])))
rownames(rSFSMed) <- paste0("Med", seq(1, length(rSFSMed[,1])))
minnval <- min(length(grep("n[0-9]", colnames(SimSFS))), length(grep("n[0-9]", colnames(SFSAtl))), length(grep("n[0-9]", colnames(SFSMed))))
JoinedData <- cbind(t(rSFSAtl[,grep("n[0-9]", colnames(rSFSAtl))[seq(1,minnval)]]), t(rSFSMed[,grep("n[0-9]", colnames(rSFSMed))[seq(1,minnval)]]), t(rSimSFS[,grep("n[0-9]", colnames(rSimSFS))[seq(1,minnval)]]))
pca <- PCA_calc(JoinedData)

pcacolors <- rep("black", length(colnames(JoinedData)))
pcacolors[grep("Atl", colnames(JoinedData))] <- "black" # population.colors[1]
pcacolors[grep("Med", colnames(JoinedData))] <- "black" # population.colors[2]
pcapch <- rep(1, length(colnames(JoinedData)))
pcapch[grep("Atl", colnames(JoinedData))] <- 3
pcapch[grep("Med", colnames(JoinedData))] <- 4
pcacex <- rep(.5, length(colnames(JoinedData)))
for(c in c(1:length(Nemu.combinations))){
	pcacolors[grep(paste0("Sim_", Nemu.combinations[c]), colnames(JoinedData))] <- colNe[c]
	pcapch[grep(paste0("Sim_", Nemu.combinations[c]), colnames(JoinedData))] <- 19
	pcacex[grep(paste0("Sim_", Nemu.combinations[c]), colnames(JoinedData))] <- 1.5
}

######################################################################
# Plotting
pdf(PDF, width=10, height=15)

# Het
print("Plotting Het")
par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("Standard deviation of the heterozygosity", ObsDataM$sd, ObsDataA$sd, SimHet, "sd", colNe, c(0,.1))
par(mar=c(7,7,0,7))
plot_Ne_mu(SimHet, colNe)

# Number alleles 
print("Plotting Number alleles")
par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("Propotion of biallelic sites", NAllelesAtl$PropBial, NAllelesMed$PropBial, SimNAlleles, "PropBial", colNe, c(0.8,1))
par(mar=c(7,7,0,7))
plot_Ne_mu(SimNAlleles, colNe)

# SFS
print("Plotting SFS")
par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("Propotion of singletons", SFSAtl$PropSingletons, SFSMed$PropSingletons, SimSFS, "PropSingletons", colNe, c(0,1))
par(mar=c(7,7,0,7))
plot_Ne_mu(SimNAlleles, colNe)

par(mar=c(0,7,2,7),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(2), heights=c(1,.5), TRUE)
plot_ObsSim("Propotion of singletons + doubletons", SFSAtl$PropSingleDoubletons, SFSMed$PropSingleDoubletons, SimSFS, "PropSingleDoubletons", colNe, c(0,1))
par(mar=c(7,7,0,7))
plot_Ne_mu(SimNAlleles, colNe)

# PCA SFS
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(1), heights=c(1,1), TRUE)
plot_PCA("PCA based on Site Frequency Spectrum", 
		pca$co.df$Comp1, paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		pca$co.df$Comp2, paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		pcacolors, pcapch, pcacex, NA, c(population.colors, colNe), c("Atlantic", "Mediterranean", Nemu.combinations), c(0,1), c(-0.6,.3))
plot_PCA("PCA based on Site Frequency Spectrum", 
		pca$co.df$Comp1[grep("Sim_", colnames(JoinedData))], paste("PC1", round(pca$propvar[1], digits = 2), "%"), 
		pca$co.df$Comp2[grep("Sim_", colnames(JoinedData))], paste("PC2", round(pca$propvar[2], digits = 2), "%"), 
		pcacolors[grep("Sim_", colnames(JoinedData))], pcapch[grep("Sim_", colnames(JoinedData))], pcacex[grep("Sim_", colnames(JoinedData))], NA, c(colNe), c(Nemu.combinations), c(0,1), c(-0.6,.3))


par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,4,4),nrow=3,ncol=2,byrow=F), widths=c(6,6), heights=c(2,2,2), TRUE)


# A
par(mar=c(0,8,4,5))
plot_ObsSim("Standard deviation\nof heterozygosity", ObsDataM$sd, ObsDataA$sd, SimHet, "sd", colNe, c(0,.1))
writePlotLabel("A")
# B
par(mar=c(0,8,4,5))
plot_ObsSim("Propotion of\nsingletons + doubletons", SFSAtl$PropSingleDoubletons, SFSMed$PropSingleDoubletons, SimSFS, "PropSingleDoubletons", colNe, c(0,1))
writePlotLabel("B")
# Legend
uNe <- unique(SimHet$Ne)
umu <- unique(SimHet$mu)
w <- .5
par(mar=c(0,7,0,5))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", xlim=c(0,length(uNe)+3), col=NA)
for(i in c(1:length(uNe))){
	text(i+2, 7, labels = NeValForm[which(NeVal==uNe[i])], xpd=NA)
	text(i+2, 5, labels = muValForm[which(muVal==umu[i])], xpd=NA)
	polygon(c(i+2-w/2, i+2+w/2, i+2+w/2, i+2-w/2), c(1, 1, 1+w*2, 1+w*2), col=colNe[i], border=NA, lwd=1)
}

# C
par(mar=c(4,4,4,4))
plot_PCA("", 
		pca$co.df$Comp1, paste0("PC1 (", round(pca$propvar[1], digits = 2), "%)"), 
		pca$co.df$Comp2, paste0("PC2 (", round(pca$propvar[2], digits = 2), "%)"), 
		pcacolors, pcapch, pcacex, NA, c("black", "black", colNe), c("Atlantic", "Mediterranean", Nemu.combinations), c(0,1), c(-0.6,.3))
writePlotLabel("C")





dev.off()















