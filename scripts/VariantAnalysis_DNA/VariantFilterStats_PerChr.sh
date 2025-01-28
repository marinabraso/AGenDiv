#!/bin/bash


inputVCF=$1
outputTable=${2%.gz}

mkdir -p $(dirname ${outputTable})

# Extracts ALL variants, also the filtered ones
gatk VariantsToTable -V ${inputVCF} -O ${outputTable} --show-filtered \
   -F CHROM \
   -F POS \
   -F TYPE \
   -F ID \
   -F REF \
   -F ALT \
   -F QUAL \
   -F DP \
   -F QD \
   -F SB \
   -F FS \
   -F SOR \
   -F BQ \
   -F MQ \
   -F MQRankSum \
   -F ReadPosRankSum \
   -F FILTER


gzip ${outputTable}


# -F CHROM: Chromosome
# -F POS: Start coordinate of the variant
# -F TYPE: SNP/...
# -F ID: Variant ID
# -F REF: Reference allele
# -F ALT: ALternative allele
# -F QUAL: Quality score out of 100
# -F DP: Combined depth across samples
# -F QD: Quality by Depth (QUAL/DP??) --> filter out QD < 2
# -F SB: Strand bias (raw counts of reads supporting each allele on the forward and reverse strand)
# -F FS: Fisher's Exact Test on strand bias (of the alternate relative to the reference) 0 = no strand bias --> filter out FS > 60
# -F SOR: Strand odds ratio (strand bias estimate) --> filter out SOR > 3
# -F BQ: RMS base quality at this position
# -F MQ: RMS mapping quality (square root of the average of the squares of the mapping qualities at the site) --> filter out MQ < 40 (or 50)
# -F MQRankSum: MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference and alternative alleles) --> filter out MQRankSum < -12.5 (?)
# -F ReadPosRankSum: ReadPosRankSumTest (compares whether the positions of the reference and alternate alleles are different within the reads) --> filter out ReadPosRankSum < -8
# -F FILTER: either PASS or the reason (filter) for which the variant was filtered out


# not used here
# -GF AD: Allelic depth (for allele for sample)







