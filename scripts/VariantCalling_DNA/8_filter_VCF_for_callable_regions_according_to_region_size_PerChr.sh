#!/bin/bash

callable=$1
inVCF=$2
ChrLengths=$3
outVCF=${4%.0_50.vcf.gz}

mkdir -p $(dirname ${outVCF})


mins=( 0 50 100 150 200 )
max=$(zcat ${callable} | awk 'BEGIN{max=0}{if(max<$4){max=$4}}END{print max}')
maxs=( 50 100 150 200 ${max} )

for b in ${!mins[@]}
do
	echo "min=${mins[$b]}; max=${maxs[$b]}"

	echo "### Select callable regions"; >&2 echo "### Select callable regions"
	zcat ${callable} | awk -v min=${mins[$b]} -v max=${maxs[$b]} '{if($4>=min && $4<max){print $1,$2,$3}}' > $(dirname ${outVCF})/$(basename ${callable%.bed.gz})_${mins[$b]}_${maxs[$b]}.bed

	echo "### Complementary bed to get non-callable regions"; >&2 echo "### Complementary bed to get non-callable regions"
	awk '{if(FNR==NR){a[$1]=$2;next}if(chr!=$1){if(chr!=""){print chr"\t"st"\t"a[">"chr];}chr=$1;st=0}end=$2-1; print chr"\t"st"\t"end; st=$3}END{print chr"\t"st"\t"a[">"chr];}' ${ChrLengths} $(dirname ${outVCF})/$(basename ${callable%.bed.gz})_${mins[$b]}_${maxs[$b]}.bed > $(dirname ${outVCF})/Non$(basename ${callable%.bed.gz})_${mins[$b]}_${maxs[$b]}.bed

	echo "### index noncallable bed"; >&2 echo "### index noncallable bed"
	gatk IndexFeatureFile -I $(dirname ${outVCF})/$(basename ${noncallable%.gz}) -O $(dirname ${outVCF})/$(basename ${noncallable%.gz}).idx

	echo "### VariantFiltration"; >&2 echo "### VariantFiltration"
	gatk VariantFiltration -V ${inVCF} \
	-mask $(dirname ${outVCF})/Non$(basename ${callable%.bed.gz})_${mins[$b]}_${maxs[$b]}.bed \
	--mask-name "not_callable" \
	-O ${outVCF}_${mins[$b]}_${maxs[$b]}.vcf.gz
done
