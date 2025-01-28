#!/bin/bash

inbamsHiSeq=$1
inbamsNovaSeq1=$2
inbamsNovaSeq2=$3
outputbams=$4
samples=$5

inbamsHiSeq=( $(echo ${inbamsHiSeq} | sed 's/ /\n/g') )
inbamsNovaSeq1=( $(echo ${inbamsNovaSeq1} | sed 's/ /\n/g') )
inbamsNovaSeq2=( $(echo ${inbamsNovaSeq2} | sed 's/ /\n/g') )
outputbams=( $(echo ${outputbams} | sed 's/ /\n/g') )
samples=( $(echo ${samples} | sed 's/ /\n/g') )

for s in ${!samples[@]}
do
	samtools sort -n -o ${inbamsHiSeq[$s]%.bam}.readsorted.bam ${inbamsHiSeq[$s]}
	samtools sort -n -o ${inbamsNovaSeq1[$s]%.bam}.readsorted.bam ${inbamsNovaSeq1[$s]}
	samtools sort -n -o ${inbamsNovaSeq2[$s]%.bam}.readsorted.bam ${inbamsNovaSeq2[$s]}
	rm -r ${outputbams[$s]}.parts 2> ~/.null
	gatk MarkDuplicatesSpark -I ${inbamsHiSeq[$s]%.bam}.readsorted.bam -I ${inbamsNovaSeq1[$s]%.bam}.readsorted.bam -I ${inbamsNovaSeq2[$s]%.bam}.readsorted.bam -O ${outputbams[$s]}
	status=$?
	rm ${inbamsHiSeq[$s]%.bam}.readsorted.ba* ${inbamsNovaSeq1[$s]%.bam}.readsorted.ba* ${inbamsNovaSeq2[$s]%.bam}.readsorted.ba* 2> ~/null
	if [[  status -eq 0 && ${inbamsHiSeq[$s]%.bam} =~ .*round3.*round2.* ]]
	then
		echo "removing ${inbamsHiSeq[$s]} ${inbamsNovaSeq1[$s]} ${inbamsNovaSeq2[$s]}"
		rm ${inbamsHiSeq[$s]%.bam}.ba* ${inbamsNovaSeq1[$s]%.bam}.ba* ${inbamsNovaSeq2[$s]%.bam}.ba* 2> ~/null
	fi

done

