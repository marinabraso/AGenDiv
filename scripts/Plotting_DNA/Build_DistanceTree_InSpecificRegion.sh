#!/bin/bash


inputsConsesusSeq=$1
inBED=$2
outtree=$3
outfa=$4
threads=$5

mkdir -p $(dirname ${outtree})
if [[ -s ${outtree} ]] || [[ -s ${outfa} ]]
then
	rm ${outtree} ${outfa} 2> ~/null
fi

for s in ${inputsConsesusSeq}
do
	sample=$(basename ${s%.*.*} | cut -f2 -d'_')
	echo ${sample} 2> {log.out}
	zcat ${s} > ${outtree%.*}_${sample}.fa
	bedtools getfasta -fi ${outtree%.*}_${sample}.fa -bed ${inBED} | sed "s/>.\\+/>${sample}/g" >> ${outfa}
	rm ${outtree%.*}_${sample}.f* 2> ~/null
done
iqtree -redo -nt ${threads} -s ${outfa} -m GTR --prefix ${outtree%.*} > {log.out} 2>> {log.err}
