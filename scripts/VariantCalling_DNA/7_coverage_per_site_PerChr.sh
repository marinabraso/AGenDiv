#!/bin/bash

chrBAMs=$1
covgenomeTARGZ=$2
mincovTARGZ=$3
samples=$4
chr=$5

covgenome=${covgenomeTARGZ%.tar.gz}
echo ${covgenome}
mkdir -p $(dirname ${covgenome})

chrBAMs=( $(echo ${chrBAMs} | sed 's/ /\n/g') )
samples=( $(echo ${samples} | sed 's/ /\n/g') )
rm $(dirname ${covgenome})/Depth.allsamples.split* 2> ~/null

for s in ${!samples[@]}
do
	echo ${samples[$s]} ${chrBAMs[$s]}

	# samtools depth
	samtools depth -Q 24 -r ${chr} -aa ${chrBAMs[$s]} -o ${covgenome}.${samples[$s]}.tmp
	status=$?
	if [[  status -eq 0 && -s ${covgenome}.${samples[$s]}.tmp ]]
	then
		echo "removing ${chrBAMs[$s]}"
		rm ${chrBAMs[$s]} 2> ~/null
	fi

	# Join for all samples
	split -l 10000000 ${covgenome}.${samples[$s]}.tmp ${covgenome}.${samples[$s]}.split
	rm ${covgenome}.${samples[$s]}.tmp 2> ~/null
	for f in $(ls ${covgenome}.${samples[$s]}.split* | rev| cut -c1,2 | rev)
	do
		echo ${f}
		if [[ -s ${covgenome}.split${f} ]]
		then
			paste ${covgenome}.split${f} <(cut -f3 ${covgenome}.${samples[$s]}.split${f}) > ${covgenome}.tmp${f}
			mv ${covgenome}.tmp${f} ${covgenome}.split${f} 
		else
			mv ${covgenome}.${samples[$s]}.split${f} ${covgenome}.split${f}
		fi
	done
	rm ${covgenome}.${samples[$s]}.split* 2> ~/null
done

# Reunify all splitted files
cat ${covgenome}.split* | gzip > ${covgenomeTARGZ}
rm ${covgenome}.split* 2> ~/null

## Count number of covered (depth >=2 & >=3) samples per site
#zcat ${covgenomeTARGZ} | awk -F '\t' '{num2=0;num3=0; for(i=3;i<=NF;i++){if($i>=2){num2++}if($i>=3){num3++}} print $1"\t"$2"\t"num2"\t"num3}' | gzip > $(dirname ${covgenome})/Depth.NumCov23.tbl

# Calc minimum coverage per site
zcat ${covgenomeTARGZ} | awk -F '\t' '{min=1000000; for(i=3;i<=NF;i++){if($i<min){min=$i}} print $1"\t"$2"\t"min}' | gzip > ${mincovTARGZ}

## From min coverage, substitute all numbers bigger than 9 with a + and convert the format to fasta
#zcat ${mincovTARGZ} | awk '{if($3>9){$3="+"}print $0}' | awk '{if(c==$1){if(n<60){seq=seq$3;n++}else{print seq;seq=$3;n=1}}else{print seq;seq=$3;n=1;c=$1;print ">"c}}END{print seq}' | tail -n +2 > ${mincovTARGZ%.tar.gz}.fa









