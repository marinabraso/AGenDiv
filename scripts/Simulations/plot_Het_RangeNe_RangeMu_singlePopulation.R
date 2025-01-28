#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
strfiles <- args[1]
PDF <- args[2]
Rconfig <- args[3]
samplesize <- args[4]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#####################

files <- unlist(strsplit(strfiles, " "))
#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)



######################################################################
# Functions

plot_sdHet_Ne_mu <- function(main, df, colNe, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,length(unique(df$Ne))+1), col=NA)
	mtext("", side = 1, line = 4, cex=1.2)
	mtext("Heterozygosity (%)", side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	uNemu <- unique(paste(df$Ne, df$mu, sep=" "))
	print(uNemu)
	leglist <- c()
	for(nm in c(1:length(uNemu))){
		print(uNemu[nm])
		for(rep in unique(df$Rep)){
			hetvalues <- df[which(paste(df$Ne, df$mu, sep=" ")==uNemu[nm] & df$Rep==rep), paste0("Het",c(1:samplesize))]
			#print(hetvalues[1,])
			print(sd(hetvalues[1,]))
			points(nm, sd(hetvalues[1,]), pch=21, cex=1, col="black")
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
header <- c("Generation", paste0("Het",c(1:samplesize)), "mu", "Ne", "Rep", "TypeOfAlleles")
Data <- data.frame(matrix(ncol = 25, nrow = 0))
colnames(Data) <- header
for(f in c(1:length(files))){
	data <- read.table(files[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	arr <- unlist(strsplit(files[f],"/"))
	arr <- unlist(strsplit(arr[5],"_"))
	chrsize <- as.numeric(arr[3])
	mu <- as.numeric(arr[4])
	Ne <- as.numeric(arr[5])
	Rep <- as.numeric(arr[6])
	arr <- unlist(strsplit(arr[7],".het"))
	TypeOfAlleles <- arr[1]
	data <- cbind(data, rep(mu, length(data[,1])), rep(Ne, length(data[,1])), rep(Rep, length(data[,1])), rep(TypeOfAlleles, length(data[,1])))
	colnames(data) <- header
	Data <- rbind(Data, data)
}
print(head(Data))
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(Data$Ne)))
print(unique(Data$Ne))
maxGen <- max(unique(Data$Generation))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
#plot_sdHet_Ne_mu("", Data[which(Data$Generation==maxGen & Data$TypeOfAlleles=="InfiniteAlleles"),], colNe, c(0,5))
plot_sdHet_Ne_mu("", Data[which(Data$Generation==maxGen),], colNe, c(0,.1))



dev.off()






