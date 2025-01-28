#!/bin/bash

noncallable=$1
inVCF=$2
outVCF=$3
SamplesOrderInVCF=$4
chr=$5

mkdir -p $(dirname ${outVCF})
if [[  ! -s ${outVCF} ]]
then
	zcat ${noncallable} | cut -f1,2,3 > $(dirname ${outVCF})/$(basename ${noncallable%.gz})
	# index noncallable bed
	echo "### index noncallable bed"; >&2 echo "### index noncallable bed"
	gatk IndexFeatureFile -I $(dirname ${outVCF})/$(basename ${noncallable%.gz}) -O $(dirname ${outVCF})/$(basename ${noncallable%.gz}).idx
fi

# VariantFiltration
echo "### VariantFiltration"; >&2 echo "### VariantFiltration"
gatk VariantFiltration -V ${inVCF} \
-mask $(dirname ${outVCF})/$(basename ${noncallable%.gz}) \
--mask-name "not_callable" \
-O ${outVCF}
status=$?
#if [[  status -eq 0 && -s ${outVCF} ]]
#then
#	echo "${outVCF} done; removing residual files from recalibration rounds"
#	recalfolders=$(dirname ${inVCF})
#	rm ${recalfolders%3}1/*${chr}* 2> ~/null
#	rm ${recalfolders%3}2/*${chr}* 2> ~/null
#	rm ${recalfolders}/*${chr}* 2> ~/null
#fi


# create samples.txt file
echo "### samples.txt file"; >&2 echo "### samples.txt file"
zcat ${outVCF} | head -5000 | grep '^#CHROM' | cut -f 10- > ${SamplesOrderInVCF}

