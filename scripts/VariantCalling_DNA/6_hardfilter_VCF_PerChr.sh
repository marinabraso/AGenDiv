#!/bin/bash

VCF=$1
outputVCF=$2

mkdir -p $(dirname ${outputVCF})

# Hard filtering
# SNPs
gatk SelectVariants -V ${VCF} -select-type SNP -O ${VCF%.vcf.gz}_snps.vcf.gz 
gatk VariantFiltration -V ${VCF%.vcf.gz}_snps.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${VCF%.vcf.gz}_snps_filtered.vcf.gz 

# Indels
gatk SelectVariants -V ${VCF} -select-type INDEL -O ${VCF%.vcf.gz}_indels.vcf.gz 
gatk VariantFiltration -V ${VCF%.vcf.gz}_indels.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O ${VCF%.vcf.gz}_indels_filtered.vcf.gz 

# Merge filtered SNPs and INDELs
gatk MergeVcfs -I ${VCF%.vcf.gz}_snps_filtered.vcf.gz -I ${VCF%.vcf.gz}_indels_filtered.vcf.gz -O ${outputVCF}
status=$?
[ $status -eq 0 ] && echo "Done" && rm ${VCF} 2> ~/null || exit "Failed ($status)"


rm ${VCF%.vcf.gz}_snps.vcf.gz* 2> ~/null
rm ${VCF%.vcf.gz}_indels.vcf.gz* 2> ~/null
rm ${VCF%.vcf.gz}_snps_filtered.vcf.gz* 2> ~/null
rm ${VCF%.vcf.gz}_indels_filtered.vcf.gz* 2> ~/null





