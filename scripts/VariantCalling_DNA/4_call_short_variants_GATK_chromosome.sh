#!/bin/bash

chr=$1
inputBAM=$2
genome=$3
outputGVCF=$4

mkdir -p $(dirname ${outputGVCF})

samtools index ${inputBAM}

gatk --java-options "-Xmx4g" HaplotypeCaller \
	--reference ${genome} \
	--input ${inputBAM} \
	--output ${outputGVCF} \
	--native-pair-hmm-threads 32 \
	--emit-ref-confidence GVCF \
	--intervals ${chr}
status=$?

if [[  status -eq 0 && -s ${outputGVCF} ]]
then
	echo "removing ${inputBAM}"
	rm ${inputBAM}* 2> ~/.null
fi

