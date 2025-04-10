#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
######################################################################
# Libraries & functions
library(RColorBrewer)
library(png)
library(ggtree)
library(ape)

######################################################################
# Files & folders
args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPiPerRegFiles <- args[1]
strPiTotalFiles <- args[2]
strHetPerRegFiles <- args[3]
strHetTotalFiles <- args[4]
strHetFixedWindowsFiles <- args[5]
strFreqPerSiteFiles <- args[6]
SamplesOrderInVCF <- args[7]
FunctionalRegionsCallableSpan <- args[8]
PerSiteFeatureType <- args[9]
SpeciesTree <- args[10]
MapImage <- args[11]
Lynch2023 <- args[12]
Leffler2012 <- args[13]
Romiguier2014 <- args[14]
CorbettDetig2015 <- args[15]
PDF <- args[16]
REPORT <- args[17]
strAtlSamples <- args[18]
strMedSamples <- args[19]
Rconfig <- args[20]
script <- sub(".*=", "", commandArgs()[4])

source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)
PiPerRegFiles <- unlist(strsplit(strPiPerRegFiles, " "))
PiTotalFiles <- unlist(strsplit(strPiTotalFiles, " "))
HetPerRegFiles <- unlist(strsplit(strHetPerRegFiles, " "))
HetTotalFiles <- unlist(strsplit(strHetTotalFiles, " "))
HetFixedWindowsFiles <- unlist(strsplit(strHetFixedWindowsFiles, " "))
FreqPerSiteFiles <- unlist(strsplit(strFreqPerSiteFiles, " "))
AtlSamples <- unlist(strsplit(strAtlSamples, " "))
MedSamples <- unlist(strsplit(strMedSamples, " "))


######################################################################
# Functions

plot_HetPerPop <- function(Hetdf, atls, meds, ylab, zlab, xlabs, ylim){
	w <- 0.05
	flanking <- 0.2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, flanking*3+(length(atls)+length(meds))*w), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	pos <- flanking
	hetsA <- Hetdf[1,atls]
	hetsM <- Hetdf[1,meds]
	for(s in c(1:length(atls))){
		polygon(c(pos,pos+w,pos+w,pos), c(0,0,hetsA[s],hetsA[s]), col=population.colors[1], border=NA)
		pos <- pos+w
	}
	midpA <- flanking+(pos-flanking)/2
	pos <- pos+flanking
	stM <- pos
 	for(s in c(1:length(meds))){
		polygon(c(pos,pos+w,pos+w,pos), c(0,0,hetsM[s],hetsM[s]), col=population.colors[2], border=NA)
		pos <- pos+w
	}
	midpM <- stM+(pos-stM)/2
   	#violin(Hetdf[1,atls], 1, w, NA, modif_alpha(population.colors[1]), 2, 2)
	#points(jitter(rep(1, length(atls)), amount=w/3), Hetdf[,atls], pch=21, bg=modif_alpha(population.colors[1],0.2), col=modif_alpha(population.colors[1]), cex=1)
    #violin(Hetdf[1,meds], 2, w, NA, modif_alpha(population.colors[2]), 2, 2)
	#points(jitter(rep(2, length(meds)), amount=w/3), Hetdf[,meds], pch=21, bg=modif_alpha(population.colors[2],0.2), col=modif_alpha(population.colors[2]), cex=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	#axis(1, at = c(midpA, midpM), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(1, at = c(midpA, midpM), labels=xlabs, lwd=NA, las=1, line=1, cex.axis=1.5)
	box()
}

plot_Comparison_diversity_estimates <- function(df, categ, mypivalues, aPi1, aPi2, aPi3, extralab){
	distbewcateg <- 40
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=c(0,length(df[,1])+distbewcateg*(length(categ)+1)+distbewcateg*3), xaxs="i", col=NA)
	mtext("Average pairwise differences (%)", side = 2, line = 3, cex=1.5)
	pos <- distbewcateg
	midp <- c()
	for(c in categ){
		print(c)
		points(pos+c(1:length(df[which(df$class==c),1])), df$Diversity[which(df$class==c)], pch=16, col=df$colors[which(df$class==c)], cex=1.5)
		midp <- c(midp, pos+(length(df[which(df$class==c),1]))/2)
		pos <- pos+length(df[which(df$class==c),1])+distbewcateg
	}
	abline(h=mypivalues[1], lty=2)
	extrapointspos <- c(pos, pos+distbewcateg/2, pos+distbewcateg/2*2, pos+distbewcateg/2*3, pos+distbewcateg/2*4)
	points(extrapointspos[1], mypivalues[2], pch=18, col=population.colors[1], cex=2)
	text(extrapointspos[1], mypivalues[2], labels=extralab[1], pos=1, cex=1.5)
	points(extrapointspos[2], mypivalues[3], pch=18, col=population.colors[2], cex=2)
	text(extrapointspos[2], mypivalues[3], labels=extralab[2], pos=1, cex=1.5)
	points(extrapointspos[3], aPi1, pch=18, col="black", cex=2)
	text(extrapointspos[3], aPi1, labels=extralab[3], pos=1, cex=1.5)
	points(extrapointspos[4], aPi2, pch=18, col="black", cex=2)
	text(extrapointspos[4], aPi2, labels=extralab[4], pos=1, cex=1.5)
	points(extrapointspos[5], aPi3, pch=18, col="black", cex=2)
	text(extrapointspos[5], aPi3, labels=extralab[5], pos=1, cex=1.5)
	axis(1, at = midp, labels=categ, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(0,20,1), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_DifferenceOfPercOfVarSites_FunctionalRegions <- function(df, ylab, ylim, reg){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, length(reg)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.2)
	meanA <- df$PercentVarAtl[which(df$Feature=="Callable")]
	meanM <- df$PercentVarMed[which(df$Feature=="Callable")]
	for(i in c(1:length(reg))){
		difperA <- df$PercentVarAtl[which(df$Feature==reg[i])]-meanA
		difperM <- df$PercentVarMed[which(df$Feature==reg[i])]-meanM
		polygon(c(i-w,i,i,i-w), c(0,0,difperA,difperA), col=population.colors[1], border=population.colors[1])
		polygon(c(i,i+w,i+w,i), c(0,0,difperM, difperM), col=population.colors[2], border=population.colors[2])
		textdirection <- 1
		textposition <- 3
		if(difperA<0){
			textdirection <- -1
			textposition <- 1
		}
		text(i-w/2, difperA+textdirection*.5,  labels =paste0(format(round(difperA, 1), nsmall = 1),"%"), pos=textposition, srt = 90, cex=1.5)
		text(i+w/2, difperM+textdirection*.5,  labels =paste0(format(round(difperM, 1), nsmall = 1),"%"), pos=textposition, srt = 90, cex=1.5)
	}
	abline(h=0, col="black")
	axis(1, at = c(1:(length(reg))), labels=reg, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_SFS <- function(values, main, ylab, xlab, max, min, ylim, col){
	w <- 1
	h <- hist(values, plot=F, breaks=seq(min(values)-w/2, max+w/2, w))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(max-w/2, min+w/2), xaxs="i", col=NA)
	mtext(main, side = 3, line = -3, cex=1.5)
	mtext(xlab, side = 1, line = 3, cex=1.5)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	h$counts <- h$counts/sum(h$counts)
	print(paste(max, min))
	print(h$mids)
	for(i in c(length(h$mids):(length(h$mids)-(max-min)))){
		print(h$mids[i])
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	axis(1, at = seq(min,max,1), lwd.ticks=1, las=1, cex.axis=.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_Heterozygous_Sites_in_Windows <- function(hetdf, samp, col, wsize, ymax){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ymax), xlim=c(0, 10), col=NA, xaxs = "i", yaxs = "i")
	mtext("Frequency (in millions)", side = 2, line = 4, cex=1.5)
	mtext(paste0("Heterozygous sites in ", wsize, "bp windows"), side = 1, line = 4, cex=1.5)
	for(s in c(1:length(samp))){
		h <- hist(as.numeric(hetdf[,s]), breaks=seq(-0.5,wsize+.5,1), plot=FALSE)
		lines(h$mids, h$counts/1000000, col=col[s])
	}
	lines(0:50, dgeom(x = 0:50, prob = 0.03))
	axis(1, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,ymax,ymax/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

Chi.squared.test_Feature_SpanVsVariants <- function(df, vartype, feat){
	#print(paste0("### Chi.squared.test for feature ", feat, " and ", vartype))
	Contingency.table <- matrix(c(df[which(df$Feature==feat), vartype], 
						df$Callable[which(df$Feature==feat)] - df[which(df$Feature==feat), vartype], 
						df[which(df$Feature=="Callable"), vartype] - df[which(df$Feature==feat), vartype], 
						df$Callable[which(df$Feature=="Callable")] - df$Callable[which(df$Feature==feat)] - df[which(df$Feature=="Callable"), vartype] + df[which(df$Feature==feat), vartype]), 
						nrow = 2, byrow = TRUE,
						dimnames = list(Feature = c("Present", "Absent"), Variants = c("Present", "Absent")))
	chisq <- chisq.test(Contingency.table)
	#print(chisq)
	#print(format.pval(chisq$p.value, digits = 20))
	#print("Observed")
	#print(chisq$observed)
	#print("Expected")
	#print(chisq$expected)
	a <- as.numeric(Contingency.table[1,1])  # A Present, B Present
	b <- as.numeric(Contingency.table[1,2])  # A Present, B Absent
	c <- as.numeric(Contingency.table[2,1])  # A Absent, B Present
	d <- as.numeric(Contingency.table[2,2])  # A Absent, B Absent
	# Compute the Odds Ratio
	OR <- (a * d) / (b * c)
	#print("Odds Ratio")
	#print(OR)
	return(list("pval" = format.pval(chisq$p.value, digits = 5), "OR" = OR))
}

######################################################################
# Read data
VariantTypes = c("Observed_Data", "Observed_SNPs", "Observed_INDELs")
RegionTypes = c("Exons", "Introns", "Promoters", "Intergenic")
SampleGoups = c("AllSamples", "AtlSamples", "MedSamples")
AllSamples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
PopColors <- rep(population.colors[1], length(AllSamples))
PopColors[match(AllSamples, MedSamples)] <- population.colors[2]

## Comparison with other diversity estimates
#class.colors <- c("#fde725", "#b5de2b", "black", "#6ece58", "#35b779", "#1f9e89", "#26828e", "#31688e", "#3e4989", "#482878", "#440154")
#class <- c("Vertebrata", "Tunicata", "Chordata", "Echinodermata", "Arthropoda", "Nematoda", "Mollusca", "Cnidaria", "Porifera", "Fungi", "Vascular plants")
class.colors <- c("#fde725", "#b5de2b", "black", "#6ece58", "#35b779", "#1f9e89", "#26828e", "#3e4989", "#482878", "#440154")
class <- c("Vertebrata", "Tunicata", "Chordata", "Echinodermata", "Arthropoda", "Nematoda", "Mollusca", "Cnidaria", "Fungi", "Vascular plants")
print(class)
# Read & prepare Lynch2023 data
Lyn23 <- read.table(Lynch2023, sep="\t", header=TRUE)
colnames(Lyn23) <- c("Species", "class", "Diversity")
print(unique(Lyn23$class))
Lyn23$class[which(Lyn23$class=="Vertebrates")] <- "Vertebrata"
Lyn23 <- Lyn23[which(Lyn23$class%in%class),]
uLyn23 <- unique(Lyn23[,c("Species","class")])
uLyn23$Diversity <- as.numeric(lapply(uLyn23$Species, function(x){mean(as.numeric(Lyn23$Diversity[which(Lyn23$Species == x)]))}))
Lyn23 <- uLyn23
Lyn23$source <- "Lyn23"
print(head(Lyn23))
# Read & prepare Romiguier2014 data
Rom14 <- read.table(Romiguier2014, sep="\t", header=TRUE)
Rom14 <- Rom14[,c(1,6,11)]
colnames(Rom14) <- c("Species", "class", "Diversity")
Rom14$Species <- gsub("_", " ", Rom14$Species)
Rom14$Diversity <- Rom14$Diversity*100 # diversity in %
Rom14 <- Rom14[which(Rom14$class %in% class),]
Rom14$class[which(Rom14$Species %in% c("Ciona intestinalis B", "Ciona intestinalis A", "Cystodytes dellechiajei purple", "Cystodytes dellechiajei blue"))] <- "Tunicata"
Rom14$class[which(Rom14$class=="Chordata")] <- "Vertebrata"
uRom14 <- unique(Rom14[,c("Species","class")])
uRom14$Diversity <- as.numeric(lapply(uRom14$Species, function(x){mean(as.numeric(Rom14$Diversity[which(Rom14$Species == x)]))}))
Rom14 <- uRom14
Rom14$source <- "Rom14"
print(head(Rom14))
# Read & prepare Leffler2012 data
Lef12 <- read.table(Leffler2012, sep="\t", header=TRUE)
Lef12 <- Lef12[which(Lef12$ThetaW.or.Pi=="Pi"),]
Lef12 <- Lef12[which(Lef12$Type.of.chromosome=="Autosome"),]
Lef12 <- Lef12[,c(1,3,5,6,9)]
colnames(Lef12) <- c("Species", "Phylum", "Total.loci", "Diversity", "Type.of.site")
Lef12$class <- Lef12$Phylum
Lef12$class[which(Lef12$class=="Magnoliophyta")] <- "Vascular plants"
Lef12$class[which(Lef12$class=="Pinophyta")] <- "Vascular plants"
Lef12$class[which(Lef12$class=="Basidiomycota")] <- "Fungi"
Lef12$class[which(Lef12$class=="Ascomycota")] <- "Fungi"
Lef12 <- Lef12[which(Lef12$class%in%class),]
Lef12 <- Lef12[order(Lef12$class), ]
Lef12$Diversity <- lapply(Lef12$Diversity, function(x){mean(as.numeric(unlist(strsplit(x, ";"))))})
Lef12$Species <- as.character(lapply(Lef12$Species, function(x){paste(unlist(strsplit(x, " "))[c(1,2)], collapse=" ")}))
uLef12 <- unique(Lef12[,c("Species","class")])
uLef12$Diversity <- as.numeric(lapply(uLef12$Species, function(x){mean(as.numeric(Lef12$Diversity[which(Lef12$Species == x)]))}))
Lef12 <- uLef12
Lef12$class[which(Lef12$Species %in% c("Ciona savignyi", "Ciona roulei"))] <- "Tunicata"
Lef12$class[which(Lef12$class=="Chordata")] <- "Vertebrata"
Lef12$source <- "Lef12"
print(head(Lef12))
# unification of datasets; the most recent prevails
unif <- rbind(Lyn23, Rom14[which(!Rom14$Species%in%Lyn23$Species),])
unif <- rbind(unif, Lef12[which(!Lef12$Species%in%unif$Species),])
unif$colors <- class.colors[match(unif$class,class)]
unif <- unif[order(unif$class, -unif$Diversity), ]
print(head(unif))
write.table(unif, file = REPORT, append = FALSE, row.names=FALSE, sep="\t", quote = FALSE)

# % of callable length & number of vairants of different features 
FeatSpan <- read.table(FunctionalRegionsCallableSpan, sep="\t", header=FALSE)
colnames(FeatSpan) <- c("Feature", "Callable")
# Atlantic variants
FeatSpan$VarAtl[which(FeatSpan$Feature=="Callable")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarAtl[which(FeatSpan$Feature=="Exons")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Exons.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarAtl[which(FeatSpan$Feature=="Introns")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Introns.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarAtl[which(FeatSpan$Feature=="Promoters")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Promoters.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarAtl[which(FeatSpan$Feature=="Intergenic")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Intergenic.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$PercentVarAtl <- FeatSpan$VarAtl/FeatSpan$Callable*100
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$ORAtl[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Exons")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Exons")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Exons")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Exons")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Introns")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Introns")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Introns")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Introns")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Promoters")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Promoters")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Promoters")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Promoters")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Intergenic")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Intergenic")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Intergenic")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Intergenic")$OR)
# Mediterranean variants
FeatSpan$VarMed[which(FeatSpan$Feature=="Callable")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Exons")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Exons.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Introns")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Introns.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Promoters")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Promoters.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Intergenic")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Intergenic.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$PercentVarMed <- FeatSpan$VarMed/FeatSpan$Callable*100
FeatSpan$pvalMed[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$ORMed[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$pvalMed[which(FeatSpan$Feature=="Exons")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Exons")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Exons")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Exons")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Introns")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Introns")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Introns")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Introns")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Promoters")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Promoters")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Promoters")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Promoters")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Intergenic")] <- Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Intergenic")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Intergenic")] <- as.numeric(Chi.squared.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Intergenic")$OR)
print(FeatSpan)


# If TreeImage verion of the treefile done with FigTree & Incksape does not exist, create a substitute with ggtree
TreeImage <- paste0(OutFolder, "/SpeciesTree_formated.png")
if(!file.exists(TreeImage)){
	tree <- read.tree(SpeciesTree)
	print(tree)
	tree <- ape::root(tree, "Drosophila_melanogaster")
	p <- ggtree(tree, linewidth=2)
	p <- p + geom_tiplab(size=10, color="black", fontface="bold")
	p <- p + xlim(0, max(p$data$x) + 0.2)
	print(p)
	ggsave(TreeImage, width = 15, height = 15)
}

######################################################################
## Figure 1
print("#### Figure 1")
pdf(PDF, width=15, height=10)
par(oma=c(2,2,2,2))
layout(matrix(c(1,2,3,4,4,5),nrow=2,ncol=3,byrow=T), widths=c(1,1,1), heights=c(1,1), TRUE)

### A
# One-to-one orthoologs phylogenetic tree 
print("A - One-to-one orthoologs phylogenetic tree")
par(mar=c(3,3,2,2),xpd=T)
tree <- readPNG(paste0(system("pwd", intern=TRUE), "/", TreeImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(tree, 0.5, 0.5, 10, 10)
par(mar=c(3,3,2,2),xpd=F)
writePlotLabel("A")

### B
## Map of sampling
print("B - Map of sampling")
par(mar=c(3,3,2,2))
map <- readPNG(paste0(system("pwd", intern=TRUE), "/", MapImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(map, 1, 1, 10, 10)
legend("topleft", populations, pch=19, col=population.colors, bty = "n", cex=1.5, xjust = 0, yjust = 0)
writePlotLabel("B")

### C
## Heterozygosity per sample
print("C - Heterozygosity per sample")
par(mar=c(7,7,2,2))
HetAll <- read.table(grep("Observed_Data.*Callable", HetTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetAll) <- AllSamples
plot_HetPerPop(HetAll*100, AtlSamples, MedSamples, "Heterozygosity (%)", "Average pairwise differences (%)", populations, c(2.5,3))
writePlotLabel("C")

### D
## Comparison with other diversity estimates
print("D - Comparison with other diversity estimates")
par(mar=c(7,7,2,2))
PutnamPi <- 4
HuangPi <- 5.37
BiPi <- 3.04
PiTAll <- read.table(grep("Observed_Data.*Callable.*AllSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTAtl <- read.table(grep("Observed_Data.*Callable.*AtlSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTMed <- read.table(grep("Observed_Data.*Callable.*MedSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
plot_Comparison_diversity_estimates(unif, class[which(class != "Chordata")], c(PiTAll,PiTAtl,PiTMed)*100, PutnamPi, HuangPi, BiPi, c("B. lanceolatum\nAtlantic", "B. lanceolatum\nMediterranean", "B. floridae\nPutnam 2008", "B. belcheri\nHuang 2014", "B. belcheri\nBi 2020"))
writePlotLabel("D")

### E
## % of variant sites on different regions (Exons, Introns, Promoters, Intergenic)
print("E - % of variant sites on different regions (Exons, Introns, Promoters, Intergenic)")
par(mar=c(7,7,2,2))
plot_DifferenceOfPercOfVarSites_FunctionalRegions(FeatSpan, "% of variants - total % of variants", c(-5,5), RegionTypes)
writePlotLabel("E")

quit()
######################################################################
## Figure S1?
print("#### Figure S1")
layout(matrix(c(1,4,
				2,4,
				2,5,
				3,5),nrow=4,ncol=2,byrow=T), widths=c(2,1), heights=c(2/3,2/6,2/6,2/3), TRUE)

### A
## Site frequency spectrum
print("A - Site frequency spectrum")
par(mar=c(3,7,2,2))
AllFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AllSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AllFreq, "All samples", "", "", max(AllFreq), max(AllFreq)/2, c(0,0.2), "black")
writePlotLabel("A")
AtlFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AtlFreq, "Atlantic samples", "Proportion of sites", "", max(AtlFreq), max(AtlFreq)/2, c(0,0.2), population.colors[1])
MedFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*MedSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(MedFreq, "Mediterranean samples", "", "Absolute major allele frequency", max(MedFreq), max(MedFreq)/2, c(0,0.2), population.colors[2])

### B
## Distribution of the distance between het sites
print("B - Distribution of the distance between het sites")
par(mar=c(7,7,2,2))
windowsize <- 50
maxyHetSites <- 2 # millions
HetRegAll.Fixed <- read.table(grep(paste0("Observed_Data.*",windowsize), HetFixedWindowsFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
HetRegAll.Fixed <-  HetRegAll.Fixed[,4:(3+length(AllSamples))]
colnames(HetRegAll.Fixed) <- AllSamples
HetRegAll.Fixed <- HetRegAll.Fixed*windowsize
print(head(HetRegAll.Fixed))
plot_Heterozygous_Sites_in_Windows(HetRegAll.Fixed[AtlSamples], AtlSamples, PopColors[1], windowsize, maxyHetSites)
plot_Heterozygous_Sites_in_Windows(HetRegAll.Fixed[MedSamples], MedSamples, PopColors[2], windowsize, maxyHetSites)
writePlotLabel("B")


dev.off()
#############
## Report

## Write heterozygosity values to report
write("### Heterozygosity per sample", file = REPORT, append = FALSE)
write(paste(AllSamples,collapse="\t"), file = REPORT, append = TRUE)
write(paste(HetAll,collapse="\t"), file = REPORT, append = TRUE)

## Prepare Table 1
# Only Pi values?
typesofdata <- c("Data", "SNPs", "INDELs", "Data", "Data", "Data", "Data")
typesofregions <- c("Callable", "Callable", "Callable", "Exons", "Introns", "Promoters", "Intergenic")
Table1 <- cbind(c("All","Atlantic", "Mediterranean"))
for(i in c(1:length(typesofdata))){
    All <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*AllSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Atl <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*AtlSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Med <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*MedSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Table1 <- cbind(Table1, c(All, Atl, Med))
}
Table1 <- as.data.frame(Table1)
colnames(Table1) <- c("GSamples", "Total", "SNPs", "INDELs", "Exons", "Introns", "Promoters", "Intergenic")
print(Table1)
write("### Table 1", file = REPORT, append = TRUE)
write.table(Table1, file = REPORT, append = TRUE, row.names=FALSE, sep="\t", quote = FALSE)

# Add unified version of Lynch2023, Romiguier2014, Leffler2012 diversity estimates
write.table(unif, file = REPORT, append = TRUE, row.names=FALSE, sep="\t", quote = FALSE)








