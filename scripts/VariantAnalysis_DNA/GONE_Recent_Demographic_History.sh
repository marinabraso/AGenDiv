#!/bin/bash

# https://github.com/esrud/GONE downloaded on 07/09/2023

genome=$1
VCF=$2
out=$3
gonesc=$4
inputparams=$5
progfolder=$6

mkdir -p $(dirname ${out})
wd=$(pwd)
cp ${inputparams} $(dirname ${out})
cp -r $(dirname ${gonesc})/* $(dirname ${out})
cp -r ${progfolder} $(dirname ${out})

###################
# VCF to map&ped
#vcftools --gzvcf ${VCF} --plink --out $(dirname ${out})/$(basename ${VCF%.vcf*})

cd $(dirname ${out})
echo "./$(basename ${gonesc}) $(basename ${VCF%.vcf*})"
./$(basename ${gonesc}) $(basename ${VCF%.vcf*})
cd ${wd}


touch ${out}







