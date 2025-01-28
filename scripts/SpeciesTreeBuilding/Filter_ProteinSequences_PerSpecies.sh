#!/bin/bash


inputP=$1
inputT=$2
outputP=${3%.gz}
outputT=${4%.gz}
version=$5

mkdir -p $(dirname ${outputP})
rm ${outputP}* ${outputT}* 2> ~/.null
speciesInVersion=$(echo ${version} | cut -f1 -d'.')

if [[ $speciesInVersion == "Branchiostoma_lanceolatum" ]] || [[ $speciesInVersion == "Branchiostoma_floridae" ]] 
then
	zcat ${inputP} | sed 's/>.*\[protein_id=\([A-Za-z\.0-9_]\+\)\] .*/>\1/g' | sed 's/\.1//g' | sed 's/CAH/BLAN/g' | sed 's/>.*\[locus_tag=BLAG_LOCUS\([A-Za-z\.0-9_]\+\)\] .*/>BLAN\1/g' | sed 's/XP_/BFLO/g' | sed 's/NP_/BFLO/g' > ${outputP}.tmp
	zcat ${inputT} | sed 's/>.*\[protein_id=\([A-Za-z\.0-9_]\+\)\] .*/>\1/g' | sed 's/\.1//g' | sed 's/CAH/BLAN/g' | sed 's/>.*\[locus_tag=BLAG_LOCUS\([A-Za-z\.0-9_]\+\)\] .*/>BLAN\1/g' | sed 's/XP_/BFLO/g' | sed 's/NP_/BFLO/g' > ${outputT}.tmp
else
	zcat ${inputP} | sed 's/>.*gene:\([A-Za-z0-9_]*\).*transcript:\([A-Za-z0-9_]*\).*/>\1\t\2/g' | sed 's/GeneID_/LOC/g' > ${outputP}.tmp
	zcat ${inputT} | sed 's/>\([A-Za-z0-9_]*\).*gene:\([A-Za-z0-9_]*\).*/>\2\t\1/g' > ${outputT}.tmp
fi


# Choose the longuest isoform for each protein
cat ${outputP}.tmp | awk '{if($1 ~ />/){print seq; print $0; seq="";next}seq=seq""$0}END{print seq}' | tail -n +2 | \
awk '{if($1 ~ />/){gID=$1; tID=$2; next;} if(!a[gID] || length($0) > length(a[gID])){a[gID]=$0; b[gID]=tID}}END{for(i in a){print i"\t"b[i]; print a[i]}}' > ${outputP}.tmp2 
# Transcript sequences of the corresponding transcript
awk '{if(FNR==NR){a[$1"_"$2]=1;next}if($1 ~ />/){gt=$1"_"$2}if(a[gt]){print $0}}' <(cat ${outputP}.tmp2 | grep '>') ${outputT}.tmp | cut -f1 | gzip > ${outputT}.gz
# Check that all proteins had at least one transcript
awk '{if(FNR==NR){a[$1]=1;next}if($1 ~ />/){g=$1}if(a[g]){print $0}}' <(zcat ${outputT}.gz | grep '>') ${outputP}.tmp2 | cut -f1 | gzip >  ${outputP}.gz

rm ${outputT}.tmp* ${outputP}.tmp* 2> ~/null









numprot=$(zcat ${outputP}.gz | grep '>' | wc -l)
numtransc=$(zcat ${outputT}.gz | grep '>' | wc -l)
echo ${numprot}
echo ${numtransc}
if [[ ${numprot} -ne ${numtransc} ]]
then
	#rm ${outputP}.gz ${outputT}.gz
	echo "----> Not the same number of proteins and transcripts for $speciesInVersion"; >&2 echo "----> Not the same number of proteins and transcripts for $speciesInVersion"
	exit 1
fi
if [[ ${numprot} -eq 0 ]] || [[ ${numtransc} -eq 0 ]]
then
	#rm ${outputP}.gz ${outputT}.gz
	echo "----> No genes listed for $speciesInVersion"; >&2 echo "---->  No genes listed for $speciesInVersion"
	exit 1
fi



