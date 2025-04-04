#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

######################################################################
# Files & general parameters

args <- commandArgs(trailingOnly=TRUE)
tmpfiles <- args[1]
samplesstr <- args[2]
thetas <- args[3]
outTXT <- args[4]
genertime <- as.numeric(args[5])
highmu  <- as.numeric(args[6])
lowmu  <- as.numeric(args[7])
step  <- as.numeric(args[8])

print(highmu)
print(lowmu)
print(step)
samples <- unlist(strsplit(samplesstr, ";"))
thetas <- unlist(strsplit(thetas, ";"))
files <- unlist(strsplit(tmpfiles, ";"))

###########################################################################
###########################################################################
### Read Data
columns = c("sample", "k", "t_k","lambda_k","T_k_highmu","T_k_lowmu","N_k_highmu","N_k_lowmu")
data = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(data) = columns

for(s in c(1:length(samples))){
	print(samples[s])
	d <- read.table(files[s], h=F, sep = "\t", check.names = F, stringsAsFactors = F)
	d <- cbind(rep(samples[s],length(d[,1])), d)
	colnames(d) <- c("sample", "k", "t_k","lambda_k")
	N0_highmu <- as.numeric(thetas[s])/(4*highmu*step)
	print(paste("N0_highmu", N0_highmu))
	N0_lowmu <- as.numeric(thetas[s])/(4*lowmu*step)
	print(paste("N0_lowmu", N0_lowmu))
	d$T_k_highmu <- d$t_k*2*N0_highmu*genertime
	d$T_k_lowmu <- d$t_k*2*N0_lowmu*genertime
	d$N_k_highmu <- d$lambda_k*N0_highmu
	d$N_k_lowmu <- d$lambda_k*N0_lowmu
	#d$four_N_mu <- d$lambda*thetas[s]/step
	data<-rbind(data, d)
}

# N0=theta/(4*mu*step)
# N=lambda*N0=lambda*theta/(4*mu*step)
# lambda*theta/step=4*N*mu
# T = t_k*2*N0*genertime = t_k*2*genertime*theta/(4*mu*step)
# 
head(data)
write.table(data, file = outTXT, quote = F, sep="\t", col.names = TRUE, row.names = FALSE)
