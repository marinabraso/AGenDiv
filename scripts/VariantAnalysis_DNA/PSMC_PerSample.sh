#!/bin/bash

# https://github.com/lh3/psmc

consensusseq=$1
psmcout=$2
sample=$3
fq2psmcfa=$4
psmc=$5
psmc2history=$6
history2ms=$7
psmc_plot=$8


mkdir -p $(dirname ${psmcout})


###################
# fq2psmcfa
# transforms the consensus sequence into a fasta-like format where the i-th character in the output sequence indicates whether there is at least one heterozygote in the bin [100i, 100i+100)
# reads the input with kseq_init(): a FASTA and FASTQ parser in C
# s = bin size [100]
echo "----> fq2psmcfa"; >&2 echo "----> fq2psmcfa"
${fq2psmcfa} -q20 ${consensusseq} > $(dirname ${psmcout})/${sample}.psmcfa
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"


###################
# PSMC
# Program: psmc (Pairwise SMC Model)
# Version: 0.6.5-r67
# -p STR      pattern of parameters [4+5*3+4] # 4+5*3+4=23 intervals of time but 7 free interval parameters that span 4-3-3-3-3-3-4 intervals of time each
# -t FLOAT    maximum 2N0 coalescent time [15]
# -N INT      maximum number of iterations [30]
# -r FLOAT    initial theta/rho ratio [4] # only initial, it is afterwards inferred
# -p & -t have to be manually tunned to ensure:
# Tried different t values, decided to go for the default (15) as it works well 
# Tried different combinations of evently spanned intervals for p parameters, got similar results. Decided to keep the "4+25*2+4+6" (sown to work in several species) 
echo "----> psmc"; >&2 echo "----> psmc"
${psmc} -N25 -t15 -r5 -p "4+25*2+4+6" -o $(dirname ${psmcout})/${sample}.psmc $(dirname ${psmcout})/${sample}.psmcfa
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"







