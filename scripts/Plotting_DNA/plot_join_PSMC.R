#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))

######################################################################
# Files & general parameters

args <- commandArgs(trailingOnly=TRUE)
MetadataFile <- args[1]
SeaLevelFile <- args[2]
tmpfiles <- args[3]
samplesstr <- args[4]
thetas <- args[5]
outTXT <- args[6]
outPDF <- args[7]
genertime <- as.numeric(args[8])

highmu <- 10^(-8)
lowmu <- 5*10^(-9)
# Cretaceous–Paleogene (K–Pg) extinction event
kpg <- 66*10^6  # Cohen, K.M., Finney, S.C., Gibbard, P.L. & Fan, J.-X. (2013; updated) The ICS International Chronostratigraphic Chart. Episodes 36: 199-204.
glaciarmax <- 20e6
step <- 100
samples <- unlist(strsplit(samplesstr, ";"))
thetas <- unlist(strsplit(thetas, ";"))
files <- unlist(strsplit(tmpfiles, ";"))
print(tmpfiles)

colfunc <- colorRampPalette(c("olivedrab1", "limegreen", "darkgreen"))
scolors <- colfunc(length(samples))

###########################################################################
###########################################################################
### Read Data
Metadata<-read.table(MetadataFile, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F)
Metadata <- unique(Metadata[,c("Sample","Population","Sex")])
Metadata <- Metadata[which(!(Metadata$Sample=="RU5D" & Metadata$Sex=="Male")),]
Metadata <- Metadata[order(Metadata$Sample, samples),]
pcolors <- c("purple", "orange")[match(Metadata$Population, unique(Metadata$Population))]
print(samples)
print(Metadata)

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
write.table(data, file = outTXT, quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

# N0=theta/(4*mu*step)
# N=lambda*N0=lambda*theta/(4*mu*step)
# lambda*theta/step=4*N*mu
# T = t_k*2*N0*genertime = t_k*2*genertime*theta/(4*mu*step)
# 
head(data)

# Miller KG, Browning JV, Schmelz WJ, Kopp RE, Mountain GS, Wright JD. Cenozoic sea-level and cryospheric evolution from deep-sea geochemical and continental margin records. Sci Adv. 2020 May 15;6(20):eaaz1346.
# Site Age_kya deltaO18 deltaO182 SeaLevel_m
SeaLevel<-read.table(SeaLevelFile, sep="\t", header=TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
SeaLevel$Age<-SeaLevel$Age_kya*1000
head(SeaLevel)


###########################################################################
###########################################################################
### Plotting
pdf(outPDF, width=20, height=10)
par(mar=c(10,10,5,5),oma=c(6,6,6,6))
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(5), heights=c(3), TRUE)

#plot_PSMC_samples_w_mu(data, samples, "highmu", highmu, genertime, pcolors)
#plot_PSMC_samples_w_mu(data, samples, "lowmu", lowmu, genertime, pcolors)
plot_PSMC_samples_wo_mu(data, samples, c(highmu,lowmu), genertime, pcolors, glaciarmax, SeaLevel[,c("SeaLevel_m","Age")])

dev.off()

###########################################################################
###########################################################################

