#!/bin/bash

genome=$1
VCF=$2
out=${3%.gz}
outconcatenated=${4%.gz}
sample=$5



#####
# dev / debug


mkdir -p $(dirname ${out})
rm ${out}* ${outconcatenated}* 2> ~/null


# Build consensus sequence
echo "----> Build consensus sequence"; >&2 echo "----> Build consensus sequence"
bcftools consensus -I -s ${sample} -f ${genome} ${VCF} -o ${out} # -I to print IUPAC codes to have the diploid consensus sequence in an haploid-like format
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo ">"${sample} > ${outconcatenated}
cat ${out} | grep -v '>' >> ${outconcatenated}

echo "----> gzip consensus sequence"; >&2 echo "----> gzip consensus sequence"
gzip ${out} ${outconcatenated}










