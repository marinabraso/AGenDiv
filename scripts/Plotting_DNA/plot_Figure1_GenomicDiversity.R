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
CoverageFile <- args[7]
SamplesOrderInVCF <- args[8]
FunctionalRegionsCallableSpan <- args[9]
PerSiteFeatureType <- args[10]
SpeciesTree <- args[11]
MapImage <- args[12]
Metadata <- args[13]
Lynch2023 <- args[14]
Leffler2012 <- args[15]
Romiguier2014 <- args[16]
CorbettDetig2015 <- args[17]
IndepEstimates <- args[18]
PDF <- args[19]
REPORT <- args[20]
strAtlSamples <- args[21]
strMedSamples <- args[22]
Rconfig <- args[23]

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
	print("P value of a t.test comparing Atlantic and Mediterranean populations:")
	print(t.test(hetsA, hetsM)$p.value)
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

plot_Comparison_diversity_estimates <- function(df, categ, meth, pch, mypivalues, extralab){
	distbewcateg <- 10
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,20), xlim=c(0,length(df[,1])+distbewcateg*(length(categ)+1)+distbewcateg*3), xaxs="i", col=NA)
	mtext("Average pairwise differences (%)", side = 2, line = 3, cex=1.5)
	pos <- distbewcateg
	midp <- c()
	for(c in categ){
		print(c)
		points(pos+c(1:length(df[which(df$class==c),1])), df$Diversity[which(df$class==c)], pch=df$pch[which(df$class==c)], col=df$colors[which(df$class==c)], cex=1.5)
		midp <- c(midp, pos+(length(df[which(df$class==c),1]))/2)
		pos <- pos+length(df[which(df$class==c),1])+distbewcateg
	}
	abline(h=mypivalues[1], lty=2) 
	points(pos, mypivalues[2], pch=16, col=population.colors[1], cex=2)
	text(pos, mypivalues[2], labels="B. lanceolatum\nAtlantic", pos=1, cex=1.5)
	points(pos+distbewcateg/2, mypivalues[3], pch=16, col=population.colors[2], cex=2)
	text(pos+distbewcateg/2, mypivalues[3], labels="B. lanceolatum\nMediterranean", pos=1, cex=1.5)
	for(i in extralab){
		text(0,20, labels=i, pos=1, cex=1.5)
	}
	legend("topright", meth, pch=pch, col="black", bty = "n", cex=1.5, xjust = 0, yjust = 0)
	axis(1, at = midp, labels=categ, lwd.ticks=1, las=2, cex.axis=1.5)
	axis(2, at = seq(0,20,1), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_DifferenceOfPercOfVarSites_FunctionalRegions <- function(df, ylab, ylim, reg){
	w <- .3
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, length(reg)+.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
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
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/10), lwd.ticks=1, las=1, cex.axis=1.5)
	legend("topleft", populations, pch=19, col=population.colors, bty = "n", cex=1.5, xjust = 0, yjust = 0)
	box()
}

plot_SFS <- function(values, main, ylab, xlab, max, min, ylim, col){
	w <- 1
	h <- hist(values, plot=F, breaks=seq(min(values)-w/2, max+w/2, w))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(max+w/2, min-w/2), xaxs="i", col=NA)
	mtext(main, side = 3, line = -3, cex=1.5)
	mtext(xlab, side = 1, line = 3, cex=1.5)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	h$counts <- h$counts/sum(h$counts)
	for(i in c(length(h$mids):(length(h$mids)-(max-min)))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
	}
	axis(1, at = h$mids, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plot_Nalleles <- function(values, main, ylab, xlab, max, min, ylim, col){
	w <- 1
	h <- hist(values, plot=F, breaks=seq(min(values)-w/2, max+w/2, w))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(min-w/2, max+w/2), xaxs="i", col=NA)
	mtext(main, side = 3, line = -3, cex=1.5)
	mtext(xlab, side = 1, line = 3, cex=1.5)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	h$counts <- h$counts/sum(h$counts)*100
	print(h$mids)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=col, border=NA)
		text(h$mids[i], h$counts[i]+.5,  labels =paste0(format(round(h$counts[i], 1), nsmall = 1),"%"), pos=3, cex=1)
	}
	axis(1, at = h$mids, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	box()
}

plot_Heterozygous_Sites_in_Windows <- function(x, samp, col, wsize){
	w <- 1
	ymax <- 0.5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ymax), xlim=c(-.5, 10.5), col=NA, xaxs = "i", yaxs = "i")
	mtext("Density", side = 2, line = 1.5, cex=.7)
	mtext(paste0("Heterozygous sites in ", wsize, "bp windows"), side = 1, line = 1.5, cex=.7)
	mtext(samp, side = 3, line = 0, cex=.7)
	p_hat <- 1 / (1 + mean(x))
	k <- 0:max(x)
	pmf <- dgeom(k, prob = p_hat)
	h <- hist(x, probability = TRUE, plot=F, breaks = seq(-0.5, max(x) + 0.5, by = 1))
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$density[i],h$density[i]), col=col, border=NA)
	}
	lines(0:max(x), pmf, col = "black", lwd=2)
	par(mgp = c(3, 0.4, 0))
	axis(1, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=.7)
	axis(2, at = seq(0,ymax,ymax/5), lwd.ticks=1, las=1, cex.axis=.7)
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

Fisher.test_Feature_SpanVsVariants <- function(df, vartype, feat) {
	spfeat_variants <- df[which(df$Feature==feat), vartype]
	spfeat_callable <- df[which(df$Feature==feat), "Callable"]
	intergenic_variants <- df[which(df$Feature=="Intergenic"), vartype]
	intergenic_callable <- df[which(df$Feature=="Intergenic"), "Callable"]
	a <- spfeat_variants
	b <- spfeat_callable - spfeat_variants
	c <- intergenic_variants
	d <- intergenic_callable - intergenic_variants
	Contingency.table <- matrix(c(	a, b,
									c, d),
		nrow = 2,
		byrow = TRUE,
		dimnames = list(
			Feature = c("SpFeat", "Intergenic"),
			Status  = c("Variant", "No_variant")
		)
	)
	print(Contingency.table)
	fisher <- fisher.test(Contingency.table)
	return(list(
		pval = format.pval(fisher$p.value, digits = 5),
		OR   = unname(fisher$estimate)
	))
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
class <- c("Fungi", "Vascular plants", "Porifera", "Cnidaria", "Nematoda", "Annelida", "Mollusca", "Arthropoda", "Echinodermata", "Tunicata", "Vertebrata", "Cephalochordata")
class.colors <- c("#440154", "#472878", "#3e4989", "#30678e", "#25818e", "#1f9e89", "#35b779", "#6dce58", "#b5de2b", "#fde725", "#fee825", "#000000")
print(class)
method.pch <- c(16, 17, 1, 2, 3, 4, 8, 7)
method <- c("Whole-genome, multiple individuals", "Whole-genome, single individual", "Non-coding", "Single chromosome", "Transcriptome", "Silent sites in protein coding genes", "Few loci (<1000)", "Other")
print(method)
# Read & prepare Lynch2023 data
Lyn23 <- read.table(Lynch2023, sep="\t", header=TRUE)
colnames(Lyn23) <- c("Species", "class", "Diversity")
Lyn23$class[which(Lyn23$class=="Vertebrates")] <- "Vertebrata"
Lyn23 <- Lyn23[which(Lyn23$class%in%class),]
uLyn23 <- unique(Lyn23[,c("Species","class")])
uLyn23$Diversity <- as.numeric(lapply(uLyn23$Species, function(x){mean(as.numeric(Lyn23$Diversity[which(Lyn23$Species == x)]))}))
Lyn23 <- uLyn23
Lyn23$source <- "Lyn23"
Lyn23$method <- "Silent sites in protein coding genes"
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
Rom14$method <- "Transcriptome"
print(head(Rom14))
# Read & prepare Leffler2012 data
Lef12 <- read.table(Leffler2012, sep="\t", header=TRUE)
Lef12 <- Lef12[which(Lef12$ThetaW.or.Pi=="Pi"),]
Lef12 <- Lef12[,c(1,3,4,5,6,8,9,10)]
colnames(Lef12) <- c("Species", "Phylum", "Classif", "Total.loci", "Diversity", "Type.of.chromosome", "Type.of.site", "Sampling.Strategy")
Lef12$class <- Lef12$Phylum
Lef12$class[which(Lef12$class=="Magnoliophyta")] <- "Vascular plants"
Lef12$class[which(Lef12$class=="Pinophyta")] <- "Vascular plants"
Lef12$class[which(Lef12$class=="Basidiomycota")] <- "Fungi"
Lef12$class[which(Lef12$class=="Ascomycota")] <- "Fungi"
Lef12$class[which(Lef12$Classif=="Mammalia")] <- "Vertebrata"
Lef12$class[which(Lef12$Classif=="Aves")] <- "Vertebrata"
Lef12$class[which(Lef12$Classif=="Actinopterygii")] <- "Vertebrata"
Lef12$class[which(Lef12$Classif=="Reptilia")] <- "Vertebrata"
Lef12$class[which(Lef12$Classif=="Ascidiacea")] <- "Tunicata"
Lef12 <- Lef12[which(Lef12$class%in%class),]
Lef12 <- Lef12[order(Lef12$class), ]
Lef12$method <- NA
Lef12$method[which(Lef12$Total.loci == "genome-wide" & Lef12$Type.of.site == "genome-wide" & Lef12$Sampling.Strategy == "One population")] <- "Whole-genome, single individual"
Lef12$method[which(Lef12$Total.loci == "genome-wide" & Lef12$Type.of.site == "genome-wide" & Lef12$Sampling.Strategy == "One population" & Lef12$class == "Porifera")] <- "Whole-genome, multiple individuals"
Lef12$method[which(Lef12$Total.loci == "genome-wide" & Lef12$Type.of.site == "genome-wide" & grepl("Multiple populations", Lef12$Sampling.Strategy))] <- "Whole-genome, multiple individuals"
Lef12$method[which((Lef12$Total.loci == "genome-wide" | Lef12$Total.loci == "chromosome-wide") & (Lef12$Type.of.site == "synonymous" | Lef12$Type.of.site == "intronic" | Lef12$Type.of.site == "4-fold degenerate"))] <- "Silent sites in protein coding genes"
Lef12$method[which((Lef12$Total.loci == "genome-wide" | Lef12$Total.loci == "chromosome-wide") & Lef12$Type.of.site == "non-coding")] <- "Non-coding"
Lef12$method[which((Lef12$Total.loci == "genome-wide" | Lef12$Total.loci == "chromosome-wide") & Lef12$Type.of.chromosome == "Sex")] <- "Single chromosome"
Lef12$method[which(grepl("bp", Lef12$Total.loci))] <- "Few loci (<1000)" # two cases with few loci that include "bp" in the column
wmethod.Lef12 <- Lef12[which(!is.na(Lef12$method)),]
Lef12 <- Lef12[which(is.na(Lef12$method)),]
Lef12$Total.loci <- unlist(lapply(Lef12$Total.loci, function(x){sum(as.numeric(unlist(strsplit(x, ";"))))}))
Lef12$method[which(as.numeric(Lef12$Total.loci) < 1000)] <- "Few loci (<1000)"
Lef12$method[which(is.na(Lef12$method) & Lef12$Type.of.site == "synonymous")] <- "Silent sites in protein coding genes"
Lef12$method[which(is.na(Lef12$method))] <- "Other"
Lef12 <- rbind(wmethod.Lef12, Lef12)
Lef12$Diversity <- as.numeric(lapply(Lef12$Diversity, function(x){mean(as.numeric(unlist(strsplit(x, ";"))))}))
print(class(Lef12$Diversity))
Lef12$Species <- as.character(lapply(Lef12$Species, function(x){paste(unlist(strsplit(x, " "))[c(1,2)], collapse=" ")}))
Lef12$Paste <- paste(Lef12$Species, Lef12$class, Lef12$method, sep="_")
uLef12 <- unique(Lef12$Paste)
uDiversity <- as.numeric(lapply(uLef12, function(x){ if(x == "Drosophila pseudoobscura_Arthropoda_Silent sites in protein coding genes"){print(Lef12[which(paste(Lef12$Species, Lef12$class, Lef12$method, sep="_") == x),])}
																										return(mean(as.numeric(Lef12$Diversity[which(paste(Lef12$Species, Lef12$class, Lef12$method, sep="_") == x)])))}))
Lef12 <- as.data.frame(cbind(unlist(strsplit(uLef12, "_"))[seq(1,length(unlist(strsplit(uLef12, "_"))),3)], unlist(strsplit(uLef12, "_"))[seq(2,length(unlist(strsplit(uLef12, "_"))),3)], unlist(strsplit(uLef12, "_"))[seq(3,length(unlist(strsplit(uLef12, "_"))),3)], as.numeric(uDiversity)))
colnames(Lef12) <- c("Species", "class", "method", "Diversity")
Lef12$Diversity <- as.numeric(Lef12$Diversity)
Lef12$class[which(Lef12$Species %in% c("Ciona savignyi", "Ciona roulei"))] <- "Tunicata"
Lef12$class[which(Lef12$class=="Chordata")] <- "Vertebrata"
Lef12$source <- "Lef12"
print(head(Lef12))
# Read & prepare the list of independent estimates collected from literature
Iestim <- read.table(IndepEstimates, sep="\t", header=TRUE)
Iestim <- Iestim[,c(1,2,3,4)]
colnames(Iestim) <- c("Diversity", "Species", "class", "method")
Iestim$Diversity <- as.numeric(Iestim$Diversity)*100 # diversity in %
Iestim$source <- "IndepEstimates"
print(head(Iestim))
# unification of datasets; the most recent prevails
unif <- rbind(Lyn23, Rom14[which(!Rom14$Species%in%Lyn23$Species),])
unif <- rbind(unif, Lef12[which(!Lef12$Species%in%unif$Species),])
unif <- rbind(unif, Iestim)
unif$colors <- class.colors[match(unif$class,class)]
unif$pch <- method.pch[match(unif$method, method)]
unif <- unif[order(unif$class, -unif$Diversity), ]
print(head(unif))

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
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Exons")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Exons")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Exons")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Exons")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Introns")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Introns")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Introns")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Introns")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Promoters")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Promoters")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Promoters")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Promoters")$OR)
FeatSpan$pvalAtl[which(FeatSpan$Feature=="Intergenic")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Intergenic")$pval
FeatSpan$ORAtl[which(FeatSpan$Feature=="Intergenic")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarAtl", "Intergenic")$OR)
# Mediterranean variants
FeatSpan$VarMed[which(FeatSpan$Feature=="Callable")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Exons")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Exons.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Introns")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Introns.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Promoters")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Promoters.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$VarMed[which(FeatSpan$Feature=="Intergenic")] <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Intergenic.*MedSamples", FreqPerSiteFiles, value = TRUE), " | wc -l"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
FeatSpan$PercentVarMed <- FeatSpan$VarMed/FeatSpan$Callable*100
FeatSpan$pvalMed[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$ORMed[which(FeatSpan$Feature=="Callable")] <- NA
FeatSpan$pvalMed[which(FeatSpan$Feature=="Exons")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Exons")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Exons")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Exons")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Introns")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Introns")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Introns")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Introns")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Promoters")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Promoters")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Promoters")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Promoters")$OR)
FeatSpan$pvalMed[which(FeatSpan$Feature=="Intergenic")] <- Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Intergenic")$pval
FeatSpan$ORMed[which(FeatSpan$Feature=="Intergenic")] <- as.numeric(Fisher.test_Feature_SpanVsVariants(FeatSpan, "VarMed", "Intergenic")$OR)
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
layout(matrix(c(1,2,3,4,4,4),nrow=2,ncol=3,byrow=T), widths=c(1,1,1), heights=c(1,1), TRUE)

### A
# One-to-one orthoologs phylogenetic tree 
# print("A - One-to-one orthoologs phylogenetic tree")
# par(mar=c(3,3,2,2),xpd=T)
# tree <- readPNG(paste0(system("pwd", intern=TRUE), "/", TreeImage))
# plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
# rasterImage(tree, 0.5, 0.5, 10, 10)
# par(mar=c(3,3,2,2),xpd=F)
# writePlotLabel("A")

### A
## Map of sampling
print("A - Map of sampling")
par(mar=c(3,3,2,2))
map <- readPNG(paste0(system("pwd", intern=TRUE), "/", MapImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(map, 1, 1, 10, 10)
legend("topleft", populations, pch=19, col=population.colors, bty = "n", cex=1.5, xjust = 0, yjust = 0)
writePlotLabel("A")

### B
## Heterozygosity per sample
print("B - Heterozygosity per sample")
par(mar=c(7,7,2,2))
HetAll <- read.table(grep("Observed_Data.*Callable", HetTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetAll) <- AllSamples
plot_HetPerPop(HetAll*100, AtlSamples, MedSamples, "Heterozygosity (%)", "Average pairwise differences (%)", populations, c(2.5,3))
writePlotLabel("B")

### C
## % of variant sites on different regions (Exons, Introns, Promoters, Intergenic)
print("D - % of variant sites on different regions (Exons, Introns, Promoters, Intergenic)")
par(mar=c(7,7,2,2))
plot_DifferenceOfPercOfVarSites_FunctionalRegions(FeatSpan, "% of variability - genome average", c(-5,5), RegionTypes)
writePlotLabel("C")
### D
## Comparison with other diversity estimates
print("C - Comparison with other diversity estimates")
par(mar=c(7,7,2,2))
PiTAll <- read.table(grep("Observed_Data.*Callable.*AllSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTAtl <- read.table(grep("Observed_Data.*Callable.*AtlSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTMed <- read.table(grep("Observed_Data.*Callable.*MedSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
plot_Comparison_diversity_estimates(unif, class, method, method.pch, c(PiTAll,PiTAtl,PiTMed)*100, c("Branchiostoma floridae", "Branchiostoma belcheri", "Branchiostoma belcheri", "Schizophyllum commune", "Caenorhabditis brenneri", "Galeolaria caespitosa", "Echinometra sp. EZ", "Ciona savignyi"))
writePlotLabel("D")

######################################################################
## Figure S3
print("#### Figure S3")
layout(matrix(c(1,4,
				2,5,
				3,6),nrow=3,ncol=2,byrow=T), widths=c(2,1), heights=c(1,1,1), TRUE)

## Site frequency spectrum
print("Site frequency spectrum")
par(mar=c(7,7,2,2))
### A
AllFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AllSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AllFreq, "All samples", "Proportion of sites", "Absolute major allele frequency", max(AllFreq), max(AllFreq)/2, c(0,0.5), "black")
writePlotLabel("A")
### 
AtlFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AtlFreq, "Atlantic samples", "Proportion of sites", "Absolute major allele frequency", max(AtlFreq), max(AtlFreq)/2, c(0,0.5), population.colors[1])
### 
MedFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*MedSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(MedFreq, "Mediterranean samples", "Proportion of sites", "Absolute major allele frequency", max(MedFreq), max(MedFreq)/2, c(0,0.5), population.colors[2])
# Histogram of number of alleles per variants site (SNPs)
print("Histogram of number of alleles per variant site")
par(mar=c(7,7,2,2))
### C
AllNalleles <- read.table(text = system(paste0("zcat ", grep("Observed_SNPs.*Callable.*AllSamples", FreqPerSiteFiles, value = TRUE), " | cut -f4 | awk '{n=split($0,a,\",\"); print n}' "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_Nalleles(AllNalleles, "SNPs", "% of sites", "Number of alleles", max(AllNalleles), 2, c(0,100), "grey60")
writePlotLabel("B")
## Histogram of number of alleles per variants site (INDELs)
AllNalleles <- read.table(text = system(paste0("zcat ", grep("Observed_INDELs.*Callable.*AllSamples", FreqPerSiteFiles, value = TRUE), " | cut -f4 | awk '{n=split($0,a,\",\"); print n}' "), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_Nalleles(AllNalleles, "INDELs", "% of sites", "Number of alleles", max(AllNalleles), 2, c(0,100), "grey60")


######################################################################
## Figure S4
print("#### Figure S4")
layout(matrix(seq(1, length(AllSamples), 1), nrow=6,ncol=6,byrow=T), widths=c(4,4,4,4,4,4), heights=c(3,3,3,3,3,3), TRUE)
## Distribution of the distance between het sites
print("Distribution of the distance between het sites")
par(mar=c(4,4,0,0))
windowsize <- 50
maxyHetSites <- 2 # millions
HetRegAll.Fixed <- read.table(grep(paste0("Observed_Data.*",windowsize), HetFixedWindowsFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
HetRegAll.Fixed <-  HetRegAll.Fixed[,4:(3+length(AllSamples))]
colnames(HetRegAll.Fixed) <- AllSamples
HetRegAll.Fixed <- HetRegAll.Fixed*windowsize
### 
for(s in c(1:length(AtlSamples))){
	plot_Heterozygous_Sites_in_Windows(HetRegAll.Fixed[,AtlSamples[s]], AtlSamples[s], population.colors[1], windowsize)
}
for(s in c(1:length(MedSamples))){
	plot_Heterozygous_Sites_in_Windows(HetRegAll.Fixed[,MedSamples[s]], MedSamples[s], population.colors[2], windowsize)
}

dev.off()
#############
## Report

## Write heterozygosity values to report
write(paste(AllSamples,collapse="\t"), file = REPORT, append = FALSE)
write("### Heterozygosity per sample", file = REPORT, append = TRUE)
write(paste(HetAll,collapse="\t"), file = REPORT, append = TRUE)

## Write average and standard deviation of coverage per sample
Coverage <- read.table(CoverageFile, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
write("### Mean and sd coverage per sample", file = REPORT, append = TRUE)
write(paste(Coverage[,1],collapse="\t"), file = REPORT, append = TRUE)
write(paste(Coverage[,2],collapse="\t"), file = REPORT, append = TRUE)

## Write proportion of variants in several regions
typesofdata <- c("Data", "SNPs", "INDELs", "Data", "Data", "Data", "Data")
typesofregions <- c("Callable", "Callable", "Callable", "Exons", "Introns", "Promoters", "Intergenic")
Table <- cbind(c("All","Atlantic", "Mediterranean"))
for(i in c(1:length(typesofdata))){
    All <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*AllSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Atl <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*AtlSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Med <- read.table(grep(paste0("Observed_", typesofdata[i], ".*", typesofregions[i], ".*MedSamples"), PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
    Table <- cbind(Table, c(All, Atl, Med))
}
Table <- as.data.frame(Table)
colnames(Table) <- c("GSamples", "Total", "SNPs", "INDELs", "Exons", "Introns", "Promoters", "Intergenic")
print(Table)
write("### Proportion of variants", file = REPORT, append = TRUE)
write.table(Table, file = REPORT, append = TRUE, row.names=FALSE, sep="\t", quote = FALSE)

# Add unified version of Lynch2023, Romiguier2014, Leffler2012 diversity estimates
write("### Unified version of compilation of diversity estimates", file = REPORT, append = TRUE)
write.table(unif, file = REPORT, append = TRUE, row.names=FALSE, sep="\t", quote = FALSE)








