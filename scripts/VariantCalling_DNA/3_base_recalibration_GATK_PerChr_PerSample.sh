#!/bin/bash

hfVCF=$1
insampleBAM=$2
genome=$3
outsampleTXT=$4
outsampleBAM=$5


mkdir -p $(dirname ${outsampleTXT}) $(dirname ${outsampleBAM})


#########################
# Base recalibration

echo ${insampleBAM} ${outsampleTXT} ${outsampleBAM}
gatk BaseRecalibrator -I ${insampleBAM} -R $genome \
--known-sites ${hfVCF} \
-O ${outsampleTXT%.txt}.bqsr
# Apply base recalibration
gatk ApplyBQSR -R $genome -I ${insampleBAM} --bqsr-recal-file ${outsampleTXT%.txt}.bqsr -O ${outsampleBAM}
status=$?
if [[  status -eq 0 && -s ${outsampleBAM} ]]
then
	echo "BQSR done, removing ${insampleBAM}"
	rm ${insampleBAM%.bam}.ba* 2> ~/.null
fi

# Analyse Covariates
gatk AnalyzeCovariates -bqsr ${outsampleTXT%.txt}.bqsr -csv ${outsampleTXT}







