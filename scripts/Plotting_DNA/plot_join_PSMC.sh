#!/bin/bash


Rscript=$(echo $0 | sed 's/\.sh/.R/g')
Matadata=$1
psmcfilesstr=$2
SeaLevel=$3
outTXT=$4
outPDF=$5
genertime=$6
samplesstr=$7
mkdir -p $(dirname ${outPDF})

rounditerations=20
rounditerations1=$(( rounditerations+1 ))
psmcfiles=( $(echo ${psmcfilesstr} | sed 's/ /\n/g') )
samples=( $(echo ${samplesstr} | sed 's/ /\n/g') )
samplesstr=( $(echo ${samplesstr} | sed 's/ /;/g') )

for s in ${!samples[@]}
do
	echo ${samples[$s]} ${psmcfiles[$s]}
	lineStIT=$(grep -n 'RD' ${psmcfiles[$s]} | grep "${rounditerations}$" | cut -f1 -d':')
	lengthIT=$((  $(grep -n 'RD' ${psmcfiles[$s]} | grep "${rounditerations1}$" | cut -f1 -d':') - ${lineStIT} - 2 ))
	# k, t_k & lambda_k in a tmp file
	tail -n +${lineStIT} ${psmcfiles[$s]} | head -${lengthIT} | grep '^RS' | cut -f2,3,4 > ${psmcfiles[$s]%.psmc}_IT${rounditerations}.tmp
	tmpfiles=$(echo ${tmpfiles}";"${psmcfiles[$s]%.psmc}_IT${rounditerations}.tmp)
	theta=$(tail -n +${lineStIT} ${psmcfiles[$s]} | head -${lengthIT} | grep '^TR' | cut -f2)
	thetas=$(echo ${thetas}";"${theta})
done
tmpfiles=$(echo ${tmpfiles} | sed 's/^;//g')
thetas=$(echo ${thetas} | sed 's/^;//g')
echo ${samplesstr} | sed 's/;/\t/g' > ${outTXT%.txt}_thetas.txt
echo ${thetas} | sed 's/;/\t/g' >> ${outTXT%.txt}_thetas.txt

## Plotting
echo "Plotting"
${Rscript} ${Matadata} ${SeaLevel} ${tmpfiles} ${samplesstr} ${thetas} ${outTXT} ${outPDF} ${genertime}



