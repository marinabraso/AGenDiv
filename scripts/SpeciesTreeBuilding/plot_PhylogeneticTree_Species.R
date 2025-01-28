#!/usr/bin/env Rscript


######################################################################
# Libraries & functions
library(ape)
library(phangorn)

######################################################################
# Files & folders

args <- commandArgs(trailingOnly=TRUE)
script <- sub(".*=", "", commandArgs()[4])
SeqFA <- args[1]
PDF <- args[2]
Rconfig <- args[3]

#####################
# debug / develop
#script <- "./scripts/SpeciesTreeBuilding/plot_PhylogeneticTree_Species.R"
#SeqFA <- "./results/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR/Concatenate_MSA_One2OneOR.fa.gz"
#PDF <- "./results/SpeciesTreeBuilding/plot_PhylogeneticTree_Species/plot_PhylogeneticTree_Species.pdf"
#Rconfig <- "config/AmphiHetDupExp_plot.R"
#####################







#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(Rconfig)
system(paste0("mkdir -p $(dirname ", PDF, ")"))
OutFolder <- system(paste0("echo $(dirname ", PDF, ")"), intern=TRUE)


######################################################################
# Read data

Seq = read.FASTA(file=SeqFA, type = "DNA")

# Convert seq to a phyDat object
Seq_phyDat = phyDat(Seq, type = "DNA", levels = NULL)

DistancesMat <- dist.ml(Seq_phyDat, model="JC69")

tree_UPGMA  <- upgma(DistancesMat)
tree_NJ  <- NJ(DistancesMat)



######################################################################
# Plotting
pdf(PDF, width=10, height=15)
par(mar=c(7,7,2,2),oma=c(1,1,1,1))

plot(tree_NJ, "unrooted", main="NJ tree", cex = 0.4, type = "p")




dev.off()


# Bootstraping
#fit = pml(ConsensusSeq_NJ, ConsensusSeq_phyDat)
#fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
#logLik(fitJC)
#bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0), mc.cores=1)
#plotBS(midpoint(fitJC$tree), bs, p = 50, type="p", cex = 0.7)









