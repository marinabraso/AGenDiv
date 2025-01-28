#!/bin/bash


ChrGenotypesFile=$1
outtree=$2
threads=$3

mkdir -p $(dirname ${outtree})
if [[ -s ${outtree} ]]
then
	rm ${outtree} 2> ~/null
fi
lastGenLine=$(zcat ${ChrGenotypesFile} | grep -n 'Generation' | tail -1 | cut -f1 -d ':')
echo ${lastGenLine}
zcat ${ChrGenotypesFile} | tail -n +${lastGenLine} | tail -n +2 | sed 's/0/A/g' | sed 's/2/C/g' | sed 's/1/G/g' | sed 's/3/T/g' | awk '{print ">"NR; print $0}' > ${outtree%.*}.fa.tmp
iqtree -redo -nt ${threads} -s ${outtree%.*}.fa.tmp -m GTR  -st DNA --prefix ${outtree%.*}
#iqtree -redo -nt ${threads} -s ${outtree%.*}.fa.tmp -m GTR  -st MORPH --prefix ${outtree%.*}
#rm ${outtree%.*}.fa.tmp 2> ~/null
