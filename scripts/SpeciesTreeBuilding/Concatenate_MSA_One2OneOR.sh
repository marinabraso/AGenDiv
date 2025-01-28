#!/bin/bash


strOG_MSAs=$1
output=${2%.gz}
strSpecies=$3
MSAType=$4

OG_MSAs=( $(echo ${strOG_MSAs} | sed 's/ /\n/g') )
Species=( $(echo ${strSpecies} | sed 's/ /\n/g') )

rm -r $(dirname ${output})/*${MSAType}* 2> ~/null
mkdir -p $(dirname ${output})

for msa in ${strOG_MSAs}
do
	echo ${msa}
	for s in ${strSpecies}
	do
		MSAdone=$(grep 'Unable' ${msa} | wc -l)
		if [[ ${MSAdone} -eq 0 ]]
		then
			linenum=$(zcat ${msa} | grep -n ${s} | cut -f1 -d':' | awk '{num=$1+1; print num}')
			zcat ${msa} | tail -n +${linenum} | sed '/>/q' | grep -v '>' >> $(dirname ${output})"/Species_"${s}"_"${MSAType}"_Concatenated_One2OneOrthologs.tmp"
		fi
	done
done

for s in ${strSpecies}
do
	echo ">"${s} >> ${output%.fa}_tmp.fa
	cat $(dirname ${output})"/Species_"${s}"_"${MSAType}"_Concatenated_One2OneOrthologs.tmp" >> ${output%.fa}_tmp.fa
done

picard NormalizeFasta -I ${output%.fa}_tmp.fa -O ${output}
gzip ${output}

#rm ${output%.fa}_tmp.fa $(dirname ${output})"/Species_*_Concatenated_One2OneOrthologs.tmp" 2> ~/null





