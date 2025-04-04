#!/bin/bash


Rscript=$(echo $0 | sed 's/\.sh/.R/g')
psmcfilesstr=$1
output=$2
output_thetas=$3
genertime=$4
samplesstr=$5
rounditerations=$6
highmu=$7
lowmu=$8
step=$9
mkdir -p $(dirname ${output})

rounditerations1=$(( rounditerations+1 ))
psmcfiles=( $(echo ${psmcfilesstr} | sed 's/ /\n/g') )
samples=( $(echo ${samplesstr} | sed 's/ /\n/g') )
samplesstr=( $(echo ${samplesstr} | sed 's/ /;/g') )

for s in ${!samples[@]}
do
	echo ${samples[$s]} ${psmcfiles[$s]}
	lineStIT=$(grep -n 'RD' ${psmcfiles[$s]} | grep -w "${rounditerations}$" | cut -f1 -d':')
	lengthIT=$((  $(grep -n 'RD' ${psmcfiles[$s]} | grep -w "${rounditerations1}$" | cut -f1 -d':') - ${lineStIT} - 2 ))
	# k, t_k & lambda_k in a tmp file
	tail -n +${lineStIT} ${psmcfiles[$s]} | head -${lengthIT} | grep '^RS' | cut -f2,3,4 > $(dirname ${output})"/"${samples[$s]}_IT${rounditerations}.tmp
	tmpfiles=$(echo ${tmpfiles}";"$(dirname ${output})"/"${samples[$s]}_IT${rounditerations}.tmp)
	theta=$(tail -n +${lineStIT} ${psmcfiles[$s]} | head -${lengthIT} | grep '^TR' | cut -f2)
	thetas=$(echo ${thetas}";"${theta})
done
tmpfiles=$(echo ${tmpfiles} | sed 's/^;//g')
thetas=$(echo ${thetas} | sed 's/^;//g')
echo ${samplesstr} | sed 's/;/\t/g' > ${output_thetas}
echo ${thetas} | sed 's/;/\t/g' >> ${output_thetas}


${Rscript} ${tmpfiles} ${samplesstr} ${thetas} ${output} ${genertime} ${highmu} ${lowmu} ${step}

rm $(dirname ${output})/*.tmp 2> ~/null



lineStIT=$(grep -n 'RD' results/VariantAnalysis_DNA/PSMC_PerSample/F10D.psmc | grep -w "20$" | cut -f1 -d':')
lengthIT=$((  $(grep -n 'RD' results/VariantAnalysis_DNA/PSMC_PerSample/F10D.psmc | grep -w "21$" | cut -f1 -d':') - ${lineStIT} - 2 ))
tail -n +${lineStIT} results/VariantAnalysis_DNA/PSMC_PerSample/F10D.psmc | head -${lengthIT} | grep '^RS' | cut -f2,3,4