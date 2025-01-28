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
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)

######################################################################
# Read data
header <- c("Generation", "PI", "mu", "Ne", "Rep")
Data <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(Data) <- header
for(f in c(1:length(PIfiles))){
	print(PIfiles[f])
	data <- read.table(PIfiles[f], sep="\t", header=FALSE, check.names = F, stringsAsFactors = F)
	arr <- unlist(strsplit(PIfiles[f],"/"))
	arr <- unlist(strsplit(arr[5],"_"))
	chrsize <- as.numeric(arr[4])
	mu <- as.numeric(arr[5])
	Ne <- as.numeric(arr[6])
	Rep <- as.numeric(arr[7])
	data <- cbind(data, rep(mu, length(data[,1])), rep(Ne, length(data[,1])), rep(Rep, length(data[,1])))
	colnames(data) <- header
	Data <- rbind(Data, data)
}
print(Data)
print(class(Data$mu))
colfunc <- colorRampPalette(c("gold","red"))
colNe <- colfunc(length(unique(Data$Ne)))
print(colNe)
print(unique(Data$Ne))
######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

plot_lines_Ne_mu("", Data, colNe, seq(1,length(unique(Data$mu)),1), c(min(Data$Generation),max(Data$Generation)), c(0,max(Data$PI)))



dev.off()






