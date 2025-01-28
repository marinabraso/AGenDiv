#!/bin/bash

NonExtraCallableRegions=$1
NonCallableRegions=$2
SNPableRegions=$3
inVCF=$4
outVCF=$5
SamplesOrderInVCF=$6
chr=$7

mkdir -p $(dirname ${outVCF})
if [[  ! -s $(dirname ${outVCF})/$(basename ${NonExtraCallableRegions%.gz}).idx ]]
then
	echo "### index NonExtraCallableRegions bed"; >&2 echo "### index NonExtraCallableRegions bed"
	zcat ${NonExtraCallableRegions} | cut -f1,2,3 > $(dirname ${outVCF})/$(basename ${NonExtraCallableRegions%.gz})
	gatk IndexFeatureFile -I $(dirname ${outVCF})/$(basename ${NonExtraCallableRegions%.gz}) -O $(dirname ${outVCF})/$(basename ${NonExtraCallableRegions%.gz}).idx
fi
if [[  ! -s $(dirname ${outVCF})/$(basename ${NonCallableRegions%.gz}).idx ]]
then
	echo "### index NonCallableRegions bed"; >&2 echo "### index NonCallableRegions bed"
	zcat ${NonCallableRegions} | cut -f1,2,3 > $(dirname ${outVCF})/$(basename ${NonCallableRegions%.gz})
	gatk IndexFeatureFile -I $(dirname ${outVCF})/$(basename ${NonCallableRegions%.gz}) -O $(dirname ${outVCF})/$(basename ${NonCallableRegions%.gz}).idx
fi
if [[  ! -s $(dirname ${outVCF})/$(basename ${SNPableRegions%.gz}).idx ]]
then
	echo "### index SNPableRegions bed"; >&2 echo "### index SNPableRegions bed"
	zcat ${SNPableRegions} | cut -f1,2,3 > $(dirname ${outVCF})/$(basename ${SNPableRegions%.gz})
	gatk IndexFeatureFile -I $(dirname ${outVCF})/$(basename ${SNPableRegions%.gz}) -O $(dirname ${outVCF})/$(basename ${SNPableRegions%.gz}).idx
fi


# VariantFiltration
echo "### VariantFiltration NonExtraCallableRegions"; >&2 echo "### VariantFiltration NonExtraCallableRegions"
gatk VariantFiltration -V ${inVCF} \
--mask $(dirname ${outVCF})/$(basename ${NonExtraCallableRegions%.gz}) \
--mask-name "not_extra_callable" \
-filter "AF == 1.00" --filter-name "reference_private_non_variant" \
-O ${outVCF}_1

echo "### VariantFiltration NonCallableRegions"; >&2 echo "### VariantFiltration NonCallableRegions"
gatk VariantFiltration -V ${outVCF}_1 \
--mask $(dirname ${outVCF})/$(basename ${NonCallableRegions%.gz}) \
--mask-name "not_callable" \
-O ${outVCF}_2

echo "### VariantFiltration SNPableRegions"; >&2 echo "### VariantFiltration SNPableRegions"
gatk VariantFiltration -V ${outVCF}_2 \
--mask $(dirname ${outVCF})/$(basename ${SNPableRegions%.gz}) \
--filter-not-in-mask true \
--mask-name "not_SNPable" \
-O ${outVCF}

rm ${outVCF}_1* ${outVCF}_2*  2> ~/.null

echo "### Filter variants with missing genotypes"; >&2 echo "### Filter variants with missing genotypes"
echo "         Note: done with vcftools, not gatk VariantFiltration. As aconsequence the filtered variants are removed from the vcf (not tagged as not PASS)"; >&2 echo "         Note: done with vcftools, not gatk VariantFiltration. As aconsequence the filtered variants are removed from the vcf (not tagged as not PASS)"
vcftools --recode --gzvcf ${outVCF} --max-missing 1 --out ${outVCF%.vcf.gz}_nomissinggeno
nomissinggenoFile=$(ls ${outVCF%.vcf.gz}_nomissinggeno*)
mv ${nomissinggenoFile} ${outVCF%.gz}
bgzip ${outVCF%.gz}
rm ${outVCF%.gz}
tabix -p vcf ${outVCF}

# create samples.txt file
echo "### samples.txt file"; >&2 echo "### samples.txt file"
zcat ${outVCF} | head -5000 | grep '^#CHROM' | cut -f 10- > ${SamplesOrderInVCF}

