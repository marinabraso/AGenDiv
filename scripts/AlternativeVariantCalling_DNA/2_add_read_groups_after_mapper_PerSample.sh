#!/bin/bash


inbam=$1
outbam=$2
sample=$3
lane=$4
platform=$5

mkdir -p $(dirname ${outbam})

gatk AddOrReplaceReadGroups \
    -I ${inbam} \
    -O ${outbam} \
    -RGID ${sample}_${lane}_${platform} \
    -RGLB "WGSBlan" \
    -RGPL ILLUMINA \
    -RGPU ${sample}_${lane}_${platform} \
    -RGSM ${sample}
status=$?

samtools index ${outbam}
status2=$?


if [[  status -eq 0 && status2 -eq 0 && -s ${outbam} ]]
then
	echo "removing ${inbam}"
	rm ${inbam} 2> ~/.null
fi


