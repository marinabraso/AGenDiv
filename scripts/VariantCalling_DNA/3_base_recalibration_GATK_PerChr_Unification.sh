#!/bin/bash

BRecalR1HiSeq=$1
BRecalR2HiSeq=$2
BRecalR1NovaSeq1=$3
BRecalR2NovaSeq1=$4
BRecalR1NovaSeq2=$5
BRecalR2NovaSeq2=$6
BaseRecalTXT=$7
samples=$8
chr=$9

BRecalR1HiSeq=( $(echo ${BRecalR1HiSeq} | sed 's/ /\n/g') )
BRecalR2HiSeq=( $(echo ${BRecalR2HiSeq} | sed 's/ /\n/g') )
BRecalR1NovaSeq1=( $(echo ${BRecalR1NovaSeq1} | sed 's/ /\n/g') )
BRecalR2NovaSeq1=( $(echo ${BRecalR2NovaSeq1} | sed 's/ /\n/g') )
BRecalR1NovaSeq2=( $(echo ${BRecalR1NovaSeq2} | sed 's/ /\n/g') )
BRecalR2NovaSeq2=( $(echo ${BRecalR2NovaSeq2} | sed 's/ /\n/g') )
samples=( $(echo ${samples} | sed 's/ /\n/g') )
mkdir -p $(dirname ${BaseRecalTXT})
head -1 ${BRecalR1HiSeq[0]} | awk '{print "Round,Chr,Sample,"$0}' > ${BaseRecalTXT}

# Analyse Covariates
for s in ${!samples[@]}
do
	echo ${samples[$s]}
	tail -n +2 ${BRecalR1HiSeq[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "1,"chr","sample","$0}' >> ${BaseRecalTXT}
	tail -n +2 ${BRecalR2HiSeq[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "2,"chr","sample","$0}' >> ${BaseRecalTXT}
	tail -n +2 ${BRecalR1NovaSeq1[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "1,"chr","sample","$0}' >> ${BaseRecalTXT}
	tail -n +2 ${BRecalR2NovaSeq1[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "2,"chr","sample","$0}' >> ${BaseRecalTXT}
	tail -n +2 ${BRecalR1NovaSeq2[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "1,"chr","sample","$0}' >> ${BaseRecalTXT}
	tail -n +2 ${BRecalR2NovaSeq2[$s]} | awk -v sample=${samples[$s]} -v chr=${chr} '{print "2,"chr","sample","$0}' >> ${BaseRecalTXT}
done


