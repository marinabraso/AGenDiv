#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")


######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
FreqObs <- args[1]
TypesOfVariants <- args[2]
FunctionalFeatOfVariants <- args[3]
SharedPrivatePop <- args[4]
SynNonSynFile <- args[5]
SpanFuncRegions <- args[6]
PDF <- args[7]
Rconfig <- args[8]
script <- sub(".*=", "", commandArgs()[4])

#############
# Dev / debug
# FreqObs <- "results/Plotting_DNA/Observed_Data/FrequencyPerSite_inBEDregions_PerSubsetOfSamples/FrequencyPerSite_inCallableRegions_PerAllSamples.txt.gz"
# TypesOfVariants <- "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/PerSite_VariantType_PerAllSamples.txt.gz"
# FunctionalFeatOfVariants <- "results/Plotting_DNA/FunctionalFeatures_BEDs/PerSite_VariantType.txt.gz"
# SharedPrivatePop <- "results/Plotting_DNA/SharedPrivatePopulation_BEDs/AllVariants_BothPopulations_Observed.bed.gz"
# SynNonSynFile <- "results/Plotting_DNA/SynonymousNonSynonymous_BEDs/AllVar_SynNonSyn.bed.gz"
# SpanFuncRegions <- "results/Plotting_DNA/get_FunctionalRegionsCallableSpand/FunctionalRegionsCallableSpand.txt"
# PDF <- "results/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.pdf"
# Rconfig <- "config/AmphiHetDupExp_plot.R"




source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
GroupsOfSamples <- c("AllSamples", "AtlSamples", "MedSamples")

######################################################################
# Functions

plot_SFS <- function(values, mids, xlab, main, ylim, col){
	print(main)
	w <- mids[2]-mids[1]
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext("Proportion of sites", side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	#abline(v=36)
	axis(1, at = seq(h$mids[1],h$mids[length(h$mids)],5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}



plot_RelativeFreq_Histogram <- function(values, mids, xlab, ylab, main, ylim, col){
	print(main)
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

plot_ViolinPlot <- function(values, xlab, ylab, ylim){
	w <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5,1.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	d <- density(values, adjust = 2)
	ynorm <- d$y/max(d$y)*w/2
	polygon(c(1-ynorm, rev(1+ynorm)), c(d$x,rev(d$x)), col=modif_alpha(ppal.color), border=ppal.color, lwd=1)
	#points(jitter(rep(1, length(values)), amount=w/2), values, pch=21, bg=modif_alpha(ppal.color,0.2), col=modif_alpha(ppal.color))
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Histogram <- function(values, xlab, ylab, ylim, xlim, vabline){
	step <- 0.025
	h <- hist(values, plot=F, breaks=seq(-1+step/2,1,step))
	print(h$mids)
	print(h$breaks)
	w <- h$mids[2]-h$mids[1]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(xlab, side = 1, line = 3, cex=1.2)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	print(h$mids)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col="black", border=NA)
	}
	abline(v=vabline)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}


plot_SpanVsVariants_FunctionalRegions <- function(df, ylab, ylim){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, length(colnames(df))-.5), col=NA)
	mtext(ylab, side = 2, line = 3, cex=1.2)
	for(i in c(1:(length(colnames(df))-1))){
		polygon(c(i-w,i,i,i-w), c(0,0,df[1,i],df[1,i]), col="gold", border="gold")
		polygon(c(i,i+w,i+w,i), c(0,0,df[2,i],df[2,i]), col="forestgreen", border="forestgreen")
		text(i-w/2+.1, df[1,i]+8,  labels =paste0(format(round(df[1,i], 1), nsmall = 1),"%"), pos=3, srt = 90, cex=1.5)
		text(i+w/2+.1, df[2,i]+8,  labels =paste0(format(round(df[2,i], 1), nsmall = 1),"%"), pos=3, srt = 90, cex=1.5)
	}
	axis(1, at = c(1:(length(colnames(df))-1)), labels=colnames(df)[c(1:(length(colnames(df))-1))], lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	legend("topleft", c("% callable length", "% of variants"), pch=19, text.col="black", col=c("gold","forestgreen"), bty = "n", cex=2, xjust = 0, yjust = 0)
	box()
}

# with downsampling
find_generalmajorfreqAtl <- function(x){
	if(is.na(x[1])){return(NA)}
	else{
		t <- 38
		ss <- 34
		if(x[2]==x[3]){
			f <- as.numeric(x[1])
		}else{
			f <- t-as.numeric(x[1])
		}
		newfreq <- length(grep(1, sample(c(rep(1,f),rep(0,t-f)), ss)))
		if(newfreq == ss){
			return(NA)
		}else{
			return(newfreq)
		}
	}
}

find_generalmajorfreqMed <- function(x){
	if(is.na(x[1])){return(NA)}
	else{
		if(x[2]==x[3]){
			return(as.numeric(x[1]))
		}else{
			return(34-as.numeric(x[1]))
		}
	}
}


######################################################################
# Read Data


MajorFreq <- read.table(text=system(paste0("zcat < ", FreqObs, " | cut -f4 | rev | cut -f1 -d':' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
MajorAllele <- read.table(text=system(paste0("zcat < ", FreqObs, " | cut -f4 | rev | cut -f2 -d':' | cut -f1 -d',' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
NumAlleles <- read.table(text=system(paste0("zcat < ", FreqObs, " | awk '{split($4,a,\",\"); print length(a)}' "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
Type <- read.table(text=system(paste0("zcat < ", TypesOfVariants, " | cut -f4 "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FuncFeat <- read.table(text=system(paste0("zcat < ", FunctionalFeatOfVariants, " | cut -f4 "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
PopInfo <- read.table(text=system(paste0("zcat < ", SharedPrivatePop, " | cut -f7 "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
#SynNonSyn <- read.table(text=system(paste0("zcat < ", SynNonSynFile, " | cut -f4 "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
MajorFreqAtl <- read.table(text=system(paste0("zcat < ", SharedPrivatePop, " | sed 's/\t.\t/\tNA\t/g' | cut -f5 | rev | cut -f1 -d':' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
MajorFreqMed <- read.table(text=system(paste0("zcat < ", SharedPrivatePop, " | sed 's/\t.\t/\tNA\t/g' | sed 's/\t\t/\tNA\t/g' | cut -f6 | rev | cut -f1 -d':' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
MajorAlleleAtl <- read.table(text=system(paste0("zcat < ", SharedPrivatePop, " | sed 's/\t.\t/\tNA\t/g' | cut -f5 | rev | cut -f2 -d':' | cut -f1 -d',' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
MajorAlleleMed <- read.table(text=system(paste0("zcat < ", SharedPrivatePop, " | sed 's/\t.\t/\tNA\t/g' | sed 's/\t\t/\tNA\t/g' | cut -f6 | rev | cut -f2 -d':' | cut -f1 -d',' | rev "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
print(length(MajorAllele))
print(length(MajorAlleleMed))
print(length(MajorAlleleAtl))





#Data <- as.data.frame(cbind(MajorFreq, MajorFreqAtl, MajorFreqMed, NumAlleles, Type, FuncFeat, PopInfo, SynNonSyn))
Data <- as.data.frame(cbind(MajorFreq, MajorFreqAtl, MajorFreqMed, NumAlleles, Type, FuncFeat, PopInfo, MajorAllele, MajorAlleleAtl, MajorAlleleMed))
Data$MajorFreq <- as.numeric(Data$MajorFreq)
Data$MajorFreqAtl <- as.numeric(Data$MajorFreqAtl)
Data$MajorFreqMed <- as.numeric(Data$MajorFreqMed)
Data$NumAlleles <- as.numeric(Data$NumAlleles)
Data$GeneralMajorFreqAtl <- apply(Data[,c("MajorFreqAtl", "MajorAllele", "MajorAlleleAtl")], 1, find_generalmajorfreqAtl)
Data$GeneralMajorFreqMed <- apply(Data[,c("MajorFreqMed", "MajorAllele", "MajorAlleleMed")], 1, find_generalmajorfreqMed)
print(head(Data))

SpanObs <- read.table(file=SpanFuncRegions, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
Categ <- SpanObs[,1]
Span <- SpanObs[,2]
Variants <- c()
for(i in Categ[c(1:length(Categ)-1)]){
	Variants <- c(Variants, length(Data[which(Data$FuncFeat==tolower(i) | Data$FuncFeat==tolower(substring(i,1,nchar(i)-1))),1]))
}
Variants <- c(Variants, length(Data[,1]))
SpanObs <- as.data.frame(rbind(Span, Variants))
colnames(SpanObs) <- Categ
print(head(SpanObs))
RelSpanObs <- SpanObs/SpanObs$Callable

print("Number of variants")
print(length(Data[,1]))
print("% of variable sites among callable sites")
print(length(Data[,1])/SpanObs$Callable[1]*100)
print("Number of SNPs")
print(length(Data[which(Data$Type=="SNP"),1]))
print("% of SNPs")
print(length(Data[which(Data$Type=="SNP"),1])/length(Data[,1])*100)
print("Number of INDELs")
print(length(Data[which(Data$Type=="INDEL"),1]))
print("% of INDELs")
print(length(Data[which(Data$Type=="INDEL"),1])/length(Data[,1])*100)




SharBialData <- Data[which(Data$PopInfo=="VariantShared" & Data$NumAlleles==2),]
#SharBialData$GeneralMajorFreqAtl <- apply(SharBialData[,c("MajorFreqAtl", "MajorAllele", "MajorAlleleAtl")], 1, find_generalmajorfreqAtl)
#SharBialData$GeneralMajorFreqMed <- apply(SharBialData[,c("MajorFreqMed", "MajorAllele", "MajorAlleleMed")], 1, find_generalmajorfreqMed)
print(head(SharBialData))
DiffFreqAtlMed <- (SharBialData$GeneralMajorFreqAtl/34)-(SharBialData$GeneralMajorFreqMed/34)





######################################################################
# Plotting
pdf(PDF, width=15, height=10)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

# SFS
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_SFS(Data$MajorFreq, seq(72,12,-1), "Absolute major allele frequency", "All variants", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$NumAlleles==2)], seq(72,12,-1), "Absolute major allele frequency", "Biallelic variants", c(0,.6), "black")
plot_SFS(Data$MajorFreq[which(Data$NumAlleles>2)], seq(72,12,-1), "Absolute major allele frequency", "Multiallelic variants (>2 alleles)", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$Type=="SNP")], seq(72,12,-1), "Absolute major allele frequency", "SNPs", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$Type=="INDEL")], seq(72,12,-1), "Absolute major allele frequency", "INDELs", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$FuncFeat=="exon")], seq(72,12,-1), "Absolute major allele frequency", "Exons", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$FuncFeat=="intron")], seq(72,12,-1), "Absolute major allele frequency", "Introns", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$FuncFeat=="promoter")], seq(72,12,-1), "Absolute major allele frequency", "Promoters", c(0,.5), "black")
plot_SFS(Data$MajorFreq[which(Data$FuncFeat=="intergenic")], seq(72,12,-1), "Absolute major allele frequency", "Intergenic", c(0,.5), "black")



#
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_RelativeFreq_Histogram(Data$NumAlleles, seq(2,8,1), "Number of alleles", "Relative frequency", "All variants", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$Type=="SNP")], seq(2,8,1), "Number of alleles", "Relative frequency", "SNPs", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$Type=="INDEL")], seq(2,8,1), "Number of alleles", "Relative frequency", "INDELs", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$FuncFeat=="exon")], seq(2,,1), "Number of alleles", "Relative frequency", "Exons", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$FuncFeat=="intron")], seq(2,,1), "Number of alleles", "Relative frequency", "Introns", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$FuncFeat=="promoter")], seq(2,,1), "Number of alleles", "Relative frequency", "Promoters", c(0,1), "black")
plot_RelativeFreq_Histogram(Data$NumAlleles[which(Data$FuncFeat=="intergenic")], seq(2,,1), "Number of alleles", "Relative frequency", "Intergenic", c(0,1), "black")


#population.colors[1]

# Difference in MAF of shared variants in Atl and Med
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
plot_SFS(Data$MajorFreqAtl[which(Data$PopInfo=="VariantShared")], seq(38,1,-1), "Absolute major allele frequency", "Atl", c(0,.5), ppal.color)
plot_SFS(Data$MajorFreqMed[which(Data$PopInfo=="VariantShared")], seq(38,1,-1), "Absolute major allele frequency", "Med", c(0,.6), ppal.color)

plot_Histogram(DiffFreqAtlMed, "Freq. Atl - Freq. Med", "Proportion of shared variants", c(0,.5), c(-.5,.5), 0)
plot_Histogram(DiffFreqAtlMed[which(DiffFreqAtlMed!=0)], "Freq. Atl - Freq. Med", "Proportion of shared variants", c(0,.5), c(-.5,.5), 0)




plot_SFS(Data$GeneralMajorFreqAtl[which(Data$PopInfo=="PrivateA")], seq(34,1,-1), "Absolute major allele frequency", "Atl", c(0,.8), ppal.color)
plot_SFS(Data$GeneralMajorFreqMed[which(Data$PopInfo=="PrivateM")], seq(34,1,-1), "Absolute major allele frequency", "Med", c(0,.8), ppal.color)

plot_SpanVsVariants_FunctionalRegions(RelSpanObs*100, "%", c(0,100))


dev.off()















