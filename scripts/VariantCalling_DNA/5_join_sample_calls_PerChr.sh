#!/bin/bash

chr=$1
genome=$2
inputVCFs=$3
mapfile=$4
outputVCF=$5
samples=$6

mkdir -p $(dirname ${outputVCF})
rm ${mapfile} 2> ~/.null

unzipoutputVCF=${outputVCF%.*}
samples=( $(echo ${samples} | sed 's/ /\n/g') )
inputVCFs=( $(echo ${inputVCFs} | sed 's/ /\n/g') )

for s in ${!samples[@]}
do
	echo "${samples[$s]}	${inputVCFs[$s]}" >> ${mapfile}
done

DBpath=$(dirname ${outputVCF})/GenomicsDBImport.${chr}
tmppath=$(dirname ${outputVCF})/GenomicsDBImporttmp.${chr}

rm -r ${DBpath} 2> ~/.null
rm -r ${tmppath} 2> ~/.null
mkdir -p ${tmppath}

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
       --genomicsdb-workspace-path ${DBpath} \
       --sample-name-map ${mapfile} \
       --tmp-dir ${tmppath} \
       --reader-threads 5 \
       --L ${chr}
status1=$?

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R ${genome} \
   -V gendb://${DBpath} \
   -O ${unzipoutputVCF}
status2=$?

bgzip ${unzipoutputVCF}
status3=$?
tabix -p vcf ${outputVCF}
status4=$?

if [[  status1 -eq 0 && status2 -eq 0 && status3 -eq 0 && status4 -eq 0 && -s ${outputVCF} ]]
then
		echo "all done, correct"
		for s in ${!samples[@]}
		do
			echo "removing ${inputVCFs[$s]}"
			rm ${inputVCFs[$s]}* 2> ~/.null
		done
		rm -r ${DBpath} 2> ~/.null
		rm -r ${tmppath} 2> ~/.null	
fi



