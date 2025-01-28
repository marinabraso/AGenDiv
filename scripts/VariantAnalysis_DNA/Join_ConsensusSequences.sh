#!/bin/bash


ConsensusSeqs=$1
out=${2%.gz}


#####
# dev / debug


mkdir -p $(dirname ${out})
rm ${out}* 2> ~/null

seqlength=0
fst=0
for i in ${ConsensusSeqs}
do
	echo $i
	sampleseqlength=$(zcat $i | grep -v '^>' | awk '{sum+=length($0)}END{print sum}')
	if [[ fst -eq 0 ]]
	then
		fst=1
		seqlength=${sampleseqlength}
	fi
	if [[ ${sampleseqlength} -ne ${seqlength} ]]
	then
		echo "$i does not have the same length as another sequence: ${sampleseqlength} -ne ${seqlength}"
		exit 1
	fi
done

zcat ${ConsensusSeqs} > ${out}
bgzip ${out}




