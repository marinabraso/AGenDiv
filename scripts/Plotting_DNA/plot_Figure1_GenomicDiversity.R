#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
######################################################################
# Libraries & functions
library(RColorBrewer)
library(png)

######################################################################
# Files & folders
args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
strPiPerRegFiles <- args[1]
strPiTotalFiles <- args[2]
#strPiFixedWindowsFiles <- args[3]
strHetPerRegFiles <- args[3]
strHetTotalFiles <- args[4]
strHetFixedWindowsFiles <- args[5]
strFreqPerSiteFiles <- args[6]
SamplesOrderInVCF <- args[7]
MapImage <- args[8]
Lynch2023 <- args[9]
Leffler2012 <- args[10]
Romiguier2014 <- args[11]
CorbettDetig2015 <- args[12]
PDF <- args[13]
REPORT <- args[14]
strAtlSamples <- args[15]
strMedSamples <- args[16]
Rconfig <- args[17]
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

read_SingleValueFiles <- function(files){
	Data <- read.table(files[1], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    if(length(files)>1){
        for(f in c(2:length(files))){
            Data <- c(Data, read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F))
        }
    }
	return(Data)
}

plot_HetPerPop <- function(Hetdf, atls, meds, pivalues, ylab, zlab, xlabs, ylim){
	w <- 0.7
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5, 2.5), col=NA)
	mtext(ylab, side = 2, line = 4, cex=1.5)
	#mtext(zlab, side = 4, line = 3, cex=1.5)
    violin(Hetdf[1,atls], 1, w, NA, modif_alpha(population.colors[1]), 2)
	points(jitter(rep(1, length(atls)), amount=w/3), Hetdf[,atls], pch=21, bg=modif_alpha(population.colors[1],0.2), col=modif_alpha(population.colors[1]), cex=2)
    violin(Hetdf[1,meds], 2, w, NA, modif_alpha(population.colors[2]), 2)
	points(jitter(rep(2, length(meds)), amount=w/3), Hetdf[,meds], pch=21, bg=modif_alpha(population.colors[2],0.2), col=modif_alpha(population.colors[2]), cex=2)
	points(1, pivalues[2], pch=18, col=population.colors[1], cex=2)
	points(2, pivalues[3], pch=18, col=population.colors[2], cex=2)
	abline(h=pivalues[1], lty=2)
	#axis(4, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1,2), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1,2), labels=xlabs, lwd=NA, las=1, line=2, cex.axis=2)
	box()
}

plot_PiPerRegionPerPop <- function(RegFiles, TotalFiles, ylab, xlabs, ylim){
	w <- 0.2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, 10), col=NA, xaxs = "i", yaxs = "i")
	mtext(ylab, side = 2, line = 4, cex=1.5)
    # Exon
    PiRegAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Exons*100), col = types.of.features.color[1], lwd=2, lty = 1)
    # Intron
    PiRegAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
   lines(ecdf(PiRegAll.Introns*100), col = types.of.features.color[2], lwd=2, lty = 1)
    # Promoter
    PiRegAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Promoters*100), col = types.of.features.color[3], lwd=2, lty = 1)
    # Intergenic
    PiRegAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
    lines(ecdf(PiRegAll.Intergenic*100), col = types.of.features.color[4], lwd=2, lty = 1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1:length(xlabs)), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = c(1:length(xlabs)), labels=xlabs, lwd=NA, las=1, line=2, cex.axis=2)
	box()

    ##########
    plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0, 10), col=NA, xaxs = "i", yaxs = "i")
	mtext("Density", side = 2, line = 3, cex=1.5)
	mtext(ylab, side = 1, line = 4, cex=1.5)
    # Exon
    PiRegAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Exons <- read.table(grep("Observed_Data.*Exons.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	d.Exons <- density(as.numeric(PiRegAll.Exons*100), adjust = 1.5)
    # Intron
    PiRegAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Introns <- read.table(grep("Observed_Data.*Introns.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	d.Introns <- density(as.numeric(PiRegAll.Introns*100), adjust = 1.5)
    # Promoter
    PiRegAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Promoters <- read.table(grep("Observed_Data.*Promoters.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	d.Promoters <- density(as.numeric(PiRegAll.Promoters*100), adjust = 1.5)
    # Intergenic
    PiRegAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", RegFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
    PiTAll.Intergenic <- read.table(grep("Observed_Data.*Intergenic.*AllSamples", TotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	d.Intergenic <- density(as.numeric(PiRegAll.Intergenic*100), adjust = 1.5)


    ydist <- 0.2
    height <- 0.4
    maxy <- max(d.Exons$y, d.Introns$y, d.Promoters$y, d.Intergenic$y)
	polygon(d.Exons$x, ydist*3 + d.Exons$y/maxy*height, col=modif_alpha(types.of.features.color[1]), border=types.of.features.color[1], lwd=1)
 	points(PiTAll.Exons*100, ydist*3+height/10, pch=19, col=types.of.features.color[1], cex=1.5)   
	polygon(d.Introns$x, ydist*2 + d.Introns$y/maxy*height, col=modif_alpha(types.of.features.color[2]), border=types.of.features.color[2], lwd=1)
 	points(PiTAll.Introns*100, ydist*2+height/10, pch=19, col=types.of.features.color[2], cex=1.5)   
	polygon(d.Promoters$x, ydist*1 + d.Promoters$y/maxy*height, col=modif_alpha(types.of.features.color[3]), border=types.of.features.color[3], lwd=1)
 	points(PiTAll.Promoters*100, ydist*1+height/10, pch=19, col=types.of.features.color[3], cex=1.5)   
	polygon(d.Intergenic$x, d.Intergenic$y/maxy*height, col=modif_alpha(types.of.features.color[4]), border=types.of.features.color[4], lwd=1)
 	points(PiTAll.Intergenic*100, height/10, pch=19, col=types.of.features.color[4], cex=1.5)   
	axis(1, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=1.5)
    legend("topright", types.of.features, pch=19, text.col="black", col=types.of.features.color, bty = "n", cex=1, xjust = 0, yjust = 0)
	box()

}

plot_SFS <- function(values, main, mids, ylim){
	w <- mids[2]-mids[1]
	h <- hist(values, plot=F, breaks=seq(mids[1]-w/2, mids[length(mids)]+w/2, w))
	#print(h$counts)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(mids[1]-w/2, mids[length(mids)]+w/2), col=NA)
	mtext("", side = 1, line = 3, cex=1.2)
	mtext("Proportion of sites", side = 2, line = 3, cex=1.2)
	h$counts <- h$counts/sum(h$counts)
	for(i in c(1:length(h$mids))){
		polygon(c(h$mids[i]-w/2,h$mids[i]+w/2,h$mids[i]+w/2,h$mids[i]-w/2), c(0,0,h$counts[i],h$counts[i]), col=ppal.color, border=NA)
	}
	#abline(v=32)
	axis(1, at = seq(h$mids[1],h$mids[length(h$mids)],1), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
	#return(h)
}

plot_Comparison_diversity_estimates <- function(df, categ){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=c(1,length(unif[,1])+10*length(categ)), col=NA)
	mtext("Diversity (%)", side = 2, line = 3, cex=1.2)
	pos <- 0
	midp <- c()
	for(c in categ){
		print(c)
		points(pos+c(1:length(df[which(df$class==c),1])), df$Diversity[which(df$class==c)], pch=16, col=unif$colors[which(df$class==c)])
		midp <- c(midp, pos+(length(df[which(df$class==c),1]))/2)
		pos <- pos+length(df[which(df$class==c),1])+10
	}
	print(pos)
	abline(h=PiTAll*100, lty=2)
	axis(1, at = midp, labels=categ, lwd.ticks=1, las=2, cex.axis=1)
	axis(2, at = seq(0,20,1), lwd.ticks=1, las=1, cex.axis=1)
	box()

}

######################################################################
# Read data & plotting
VariantTypes = c("Observed_Data", "Observed_SNPs", "Observed_INDELs")
RegionTypes = c("Exons", "Introns", "Promoters", "Intergenic")
SampleGoups = c("AllSamples", "AtlSamples", "MedSamples")
AllSamples <- read.table(SamplesOrderInVCF, sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)


## Start plot
pdf(PDF, width=10, height=15)
par(oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

### A
## Map of sampling
par(mar=c(3,3,2,2))
map <- readPNG(paste0(system("pwd", intern=TRUE), "/", MapImage))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", col=NA, xaxs = "i", yaxs = "i")
rasterImage(map, 1, 1, 10, 10)
legend("topleft", populations, pch=19, col=population.colors, bty = "n", cex=2, xjust = 0, yjust = 0)
writePlotLabel("A")

### B
## Plot heterozygosity per sample
par(mar=c(7,7,2,2))
HetAll <- read.table(grep("Observed_Data.*Callable", HetTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetAll) <- AllSamples
PiTAll <- read.table(grep("Observed_Data.*Callable.*AllSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTAtl <- read.table(grep("Observed_Data.*Callable.*AtlSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
PiTMed <- read.table(grep("Observed_Data.*Callable.*MedSamples", PiTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[1,1]
plot_HetPerPop(HetAll*100, AtlSamples, MedSamples, c(PiTAll,PiTAtl,PiTMed)*100, "Heterozygosity (%)", "Average pairwise differences (%)", populations, c(0,5))
writePlotLabel("B")

### C
## Plots of pi on different regions (Exons, Introns, Promoters, Intergenic)
# Fix: take regions into account as a whole, not divided in its callable regions
#plot_PiPerRegionPerPop(PiPerRegFiles, PiTotalFiles, "Average pairwise differences (%)", RegionTypes, c(0,1))
plot.new()

### D
## Distribution of het in fixed-length windows
HetRegAll.Fixed <- read.table(grep("Observed_Data", HetFixedWindowsFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
HetRegAll.Fixed <-  HetRegAll.Fixed[,4:(3+length(AllSamples))]
colnames(HetRegAll.Fixed) <- AllSamples
HetTAll <- read.table(grep("Observed_Data.*Callable", HetTotalFiles, value = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
colnames(HetTAll) <- AllSamples
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.5), xlim=c(0, 10), col=NA, xaxs = "i", yaxs = "i")
mtext("Density", side = 2, line = 3, cex=1.5)
mtext("Heterozygosity (%)", side = 1, line = 4, cex=1.5)
for(s in AllSamples){
    d <- density(as.numeric(HetRegAll.Fixed[,s]*100), adjust = 1)
    polygon(d$x, d$y, col=NA, border=ppal.color, lwd=1)
    points(HetTAll[s]*100, max(d$y)/10, pch=19, col=ppal.color, cex=1.5)   
}
    d <- density(as.numeric(rowMeans(HetRegAll.Fixed)*100), adjust = 1)
    polygon(d$x, d$y, col=NA, border="black", lwd=1)
    points(rowMeans(HetTAll)*100, max(d$y)/10, pch=19, col="black", cex=1.5)   

axis(1, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=1.5)
box()
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=6,ncol=2,byrow=T), widths=c(2,2), heights=c(1), TRUE)

### E
## Site frequency spectrum
AllFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AllSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AllFreq, "All samples", seq(max(AllFreq),min(AllFreq),-1), c(0,0.5))
AtlFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*AtlSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(AtlFreq, "Atlantic samples", seq(max(AtlFreq),min(AtlFreq),-1), c(0,0.5))
MedFreq <- read.table(text = system(paste0("zcat ", grep("Observed_Data.*Callable.*MedSamples", FreqPerSiteFiles, value = TRUE), " | rev | cut -f1 -d':' | rev"), intern = TRUE), sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)[,1]
plot_SFS(MedFreq, "Mediterranean samples", seq(max(MedFreq),min(MedFreq),-1), c(0,0.5))



### F
## Comparison with other diversity estimates
class.colors <- c("#fde725", "#b5de2b", "black", "#6ece58", "#35b779", "#1f9e89", "#26828e", "#31688e", "#3e4989", "#482878", "#440154")
class <- c("Vertebrata", "Tunicata", "Chordata", "Echinodermata", "Arthropoda", "Nematoda", "Mollusca", "Cnidaria", "Porifera", "Fungi", "Vascular plants")
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
plot_Comparison_diversity_estimates(unif, class[which(class != "Chordata")])



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








