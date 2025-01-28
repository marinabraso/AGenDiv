#!/bin/bash

chrVCFs=$1
outVCF=$2
passoutVCF=$3
passSNPoutVCF=$4

mkdir -p $(dirname ${outVCF})
rm ${outVCF}* ${passoutVCF}* ${passSNPoutVCF}* 2> ~/null

# Concatenate vcf
echo "----> Concatenate vcf"; >&2 echo "----> Concatenate vcf"
bcftools concat -O z ${chrVCFs} > ${outVCF}
bcftools index ${outVCF}
tabix -p vcf ${outVCF}

# Build vcf with only PASS variants
echo "----> Build vcf with only PASS variants"; >&2 echo "----> Build vcf with only PASS variants"
bcftools view -i 'FILTER="PASS"' -O z ${outVCF} > ${passoutVCF}
bcftools index ${passoutVCF}
tabix -p vcf ${passoutVCF}


# Filter vcf with vcftools to keep only SNPs
# This file will be used to build the consensus cequenses per infividual and the resulting sample tree
echo "----> Filter vcf with vcftools and bcftools to keep only SNPs"; >&2 echo "----> Filter vcf with vcftools and bcftools to keep only SNPs"
vcftools --recode --gzvcf ${passoutVCF} --remove-indels --stdout | bcftools view -e 'ALT="*"' - > ${passSNPoutVCF%.gz}
bgzip ${passSNPoutVCF%.gz}
bcftools index ${passSNPoutVCF}
tabix -p vcf ${passSNPoutVCF}





