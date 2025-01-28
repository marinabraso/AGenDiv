#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
strPIfiles <- args[1]
PDF <- args[2]
Rconfig <- args[3]
script <- sub(".*=", "", commandArgs()[4])

#####################
# debug / develop
#####################

PIfiles <- unlist(strsplit(strPIfiles, " "))
#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)



######################################################################
# Functions

plot_lines_RangeNemu <- function(main, df, colNe, xlim, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(xlim[1],xlim[2]+(xlim[2]-xlim[1])/20), col=NA)
	mtext("Generations (x1000)", side = 1, line = 4, cex=1.2)
	mtext("Average pairwise differences", side = 2, line = 4, cex=1.2)
	mtext(main, side = 3, line = 1, cex=1.2)
	uNe <- unique(sort(df$Ne))
	umu <- unique(sort(df$mu))
	leglist <- c()
	for(ne in c(1:length(uNe))){
		for(mu in c(1:length(umu))){
			if(length(df[which(df$mu==umu[mu] & df$Ne==uNe[ne]),1])>0){
				for(rep in unique(df$Rep)){
					for(typ in "FiniteAlleles"){ #for(typ in unique(df$TypeOfAlleles)){
						subdf <- df[which(df$mu==umu[mu] & df$Ne==uNe[ne] & df$Rep==rep & df$TypeOfAlleles==typ),] 
						lines(subdf$Generation, subdf$PI, lwd=2, col=colNe[ne])
						points(subdf$Generation, subdf$PI, pch=c(21, 23)[grep(typ, c("InfiniteAlleles", "FiniteAlleles"))], cex=1, col="black", bg=colNe[ne])
					}
				}
				points(xlim[2]+(xlim[2]-xlim[1])/20, 4*umu[mu]*uNe[ne], pch=20, cex=2, col="black")
				leglist <- c(leglist, paste0("Ne = ", uNe[ne], ", mu = ",umu[mu]))
			}
		}
	}
	axis(1, at = unique(df$Generation)[seq(1,length(unique(df$Generation)), 2)], labels = unique(df$Generation)[seq(1,length(unique(df$Generation)), 2)]/1000, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/6), lwd.ticks=1, las=1, cex.axis=1)
	legend("topleft", c(leglist, "Infinite sites", "Finite sites"), pch=c(rep(NA,length(uNe)),21,23), lwd=c(2,2,NA,NA), text.col="black", col=c(colNe, "black", "black"), bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

plot_lines_EachNemu <- function(df, colType, xlim, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(xlim[1],xlim[2]+(xlim[2]-xlim[1])/20), col=NA)
	mtext("Generations (x1000)", side = 1, line = 4, cex=1.2)
	mtext("Average pairwise differences", side = 2, line = 4, cex=1.2)
	mtext(paste("Ne =", df$Ne[1], ", mu =", df$mu[1]), side = 3, line = 1, cex=1.2)
	for(rep in unique(df$Rep)){
		for(typ in "FiniteAlleles"){ #for(typ in unique(df$TypeOfAlleles)){
			subdf <- df[which(df$Rep==rep & df$TypeOfAlleles==typ),]
			subdf <- subdf[order(subdf$Generation),]
			print(subdf)
			lines(subdf$Generation, subdf$PI, lwd=2, col=colType[grep(typ, c("InfiniteAlleles", "FiniteAlleles"))])
			points(subdf$Generation, subdf$PI, pch=c(21, 23)[grep(typ, c("InfiniteAlleles", "FiniteAlleles"))], cex=1, col="black", bg=colType[grep(typ, c("InfiniteAlleles", "FiniteAlleles"))])
		}
	}
	points(xlim[2]+(xlim[2]-xlim[1])/20, 4*df$mu[1]*df$Ne[1], pch=20, cex=2, col="black")
	axis(1, at = unique(df$Generation)[seq(1,length(unique(df$Generation)), 2)], labels = unique(df$Generation)[seq(1,length(unique(df$Generation)), 2)]/1000, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/6), lwd.ticks=1, las=1, cex.axis=1)
	#legend("topleft", c("Infinite sites", "Finite sites"), pch=c(21,23), text.col="black", col=colType, bty = "n", cex=1, xjust = 0, yjust = 0)
	box()
}

#####################
# MAIN
#####################

######################################################################
# Read data
header <- c("Generation", "PI", "mu", "Ne", "Rep", "TypeOfAlleles")
Data <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(Data) <- header
for(f in c(1:length(PIfiles))){
	print(PIfiles[f])
	data <- read.table(PIfiles[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	arr <- unlist(strsplit(PIfiles[f],"/"))
	arr <- unlist(strsplit(arr[5],"_"))
	chrsize <- as.numeric(arr[3])
	mu <- as.numeric(arr[4])
	Ne <- as.numeric(arr[5])
	Rep <- as.numeric(arr[6])
	arr <- unlist(strsplit(arr[7],".pi"))
	print(arr)
	TypeOfAlleles <- arr[1]
	data <- cbind(data, rep(mu, length(data[,1])), rep(Ne, length(data[,1])), rep(Rep, length(data[,1])), rep(TypeOfAlleles, length(data[,1])))
	colnames(data) <- header
	Data <- rbind(Data, data)
}
print(Data)
print(class(Data$mu))
colfunc <- colorRampPalette(c("gold","forestgreen"))
colNe <- colfunc(length(unique(Data$Ne)))

######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
for(nm in unique(paste(Data$Ne, Data$mu, sep="_"))){
	print(nm)
	plot_lines_EachNemu(Data[which(paste(Data$Ne, Data$mu, sep="_")==nm),], c("black", "darkred"), c(min(Data$Generation),max(Data$Generation)), c(0,0.06))
}



dev.off()






