#!/bin/bash


Rscript=$(echo $0 | sed 's/\.sh/.R/g')
MetadataFile=$1
strchrsGENOTYPEMatrices=$2
strchrsDividedLists=$3
SamplesOrderInVCF=$4
PDF=$5
strchrs=$6
Rconfig=$7

mkdir -p $(dirname ${PDF})

strchrs=( $(echo ${strchrs} | sed 's/ /;/g') )
chrsGENOTYPEMatrices=( $(echo ${strchrsGENOTYPEMatrices} | sed 's/ /\n/g') )
chrsDividedLists=( $(echo ${strchrsDividedLists} | sed 's/ /\n/g') )
strsamples=$(cat ${SamplesOrderInVCF} | sed 's/\s/;/g')

echo "----> Plotting"
${Rscript} ${MetadataFile} ${PDF} ${strchrs} ${strsamples} ${chrsGENOTYPEMatrices[0]%.chr*} ${chrsDividedLists[0]%.chr*} ${Rconfig}





















