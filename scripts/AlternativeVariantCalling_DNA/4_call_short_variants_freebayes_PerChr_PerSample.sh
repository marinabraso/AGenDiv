#!/bin/bash

chr=$1
inputBAMs=$2
genome=$3
outputVCF=$4

mkdir -p $(dirname ${outputVCF})


echo freebayes \
  -f ${genome} \
  #--use-best-n-alleles 4 \ 
  --theta 0.02 \
  #--haplotype-length 0 \
  --min-alternate-fraction 0.05 \
  --min-alternate-count 2 \
  --min-base-quality 10 \
  --min-mapping-quality 20 \
  --genotype-qualities \
  ${inputBAMs} > ${outputVCF}
