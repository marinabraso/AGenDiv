#!/bin/bash

ProteinSequencesFiles=$1
TranscriptSequencesFiles=$2
OG=$3
outPSeqs=${4%.gz}
outTSeqs=${5%.gz}
outPSeqsAln=${6%.gz}
outTSeqsAln=${7%.gz}
outPSeqsAlnTrim=${8%.gz}
outTSeqsAlnTrim=${9%.gz}
species=${10}

############
# debug / develop
############




ProteinSequencesFilesArray=( $(echo ${ProteinSequencesFiles} | sed 's/ /\n/g') )
TranscriptSequencesFilesArray=( $(echo ${TranscriptSequencesFiles} | sed 's/ /\n/g') )
mkdir -p $(dirname ${outPSeqs})
rm ${outPSeqs}* ${outPSeqsAln}* ${outPSeqsAlnTrim}* ${outTSeqs}* ${outTSeqsAln}* ${outTSeqsAlnTrim}* 2> ~/null


#############################################
echo "----> Get the protein & transcript sequences of all genes"; >&2 echo "----> Get the protein & transcript sequences of all genes"
gIDs=( $(zcat ${OG} | cut -f1 -d'_') )
gSpecies=( $(zcat ${OG} | cut -f2- -d'_') )

for g in ${!gIDs[@]}
do 
	echo ${gIDs[$g]} ${gSpecies[$g]}
	spidx=$(echo ${species} | sed 's/ /\n/g' | grep -n ${gSpecies[$g]} | cut -f1 -d':' | awk '{idx=$1-1; print idx}')

	# Control for repeated gene IDs
	numfound=$(zcat ${ProteinSequencesFilesArray[$spidx]} | grep -w ${gIDs[$g]} | wc -l)	
	if [[ ${numfound} -gt 1 ]]
	then 
		echo "ERROR more than one annotation found for "${gIDs[$g]}" of OG file: "${OG}
		exit 1
	fi

	# Get protein sequence
	linenum=$(zcat ${ProteinSequencesFilesArray[$spidx]} | grep -nw ${gIDs[$g]} | cut -f1 -d':' | awk '{num=$1+1; print num}')
	echo ">"${gIDs[$g]}"_"${gSpecies[$g]} >> ${outPSeqs}
	zcat ${ProteinSequencesFilesArray[$spidx]} | tail -n +${linenum} | sed '/>/q' | grep -v '>' >> ${outPSeqs}
    # Get transcript sequence
	linenum=$(zcat ${TranscriptSequencesFilesArray[$spidx]} | grep -nw ${gIDs[$g]} | cut -f1 -d':' | awk '{num=$1+1; print num}')
	if [[ $linenum = "" ]]
	then
		echo "ERROR! Transcript not found for ${gIDs[$g]} in ${gSpecies[$g]}"; >&2 echo "ERROR! Transcript not found for ${gIDs[$g]} in ${gSpecies[$g]}"
		exit 1
	fi
	echo ">"${gIDs[$g]}"_"${gSpecies[$g]} >> ${outTSeqs}
    zcat ${TranscriptSequencesFilesArray[$spidx]} | tail -n +${linenum} | sed '/>/q' | grep -v '>' >> ${outTSeqs}
done


#############################################
echo "----> Do a multiple sequence alignment with MAFFT"; >&2 echo "----> Do a multiple sequence alignment with MAFFT"

# Using L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)
# Equivalent to : mafft --localpair --maxiterate 1000 input [> output]
linsi ${outPSeqs} > ${outPSeqsAln} 

failedMSA=0
if [[ ! -s ${outPSeqsAln} ]]
then
	echo "Unable to build MSA for this OG" > ${outPSeqsAln}
	failedMSA=1
fi



#############################################
echo "----> Backtranslate the aligned protein sequences"; >&2 echo "----> Backtranslate the aligned protein sequences"

treebest backtrans ${outPSeqsAln} ${outTSeqs} > ${outTSeqsAln}

failedBacktrans=0
if [[ ! -s ${outTSeqsAln} ]]
then
    echo "Unable to backtranslate MSA for ${OG}" > ${outTSeqsAln}
	failedBacktrans=1
fi



#############################################
echo "----> Trim MSA with ClipKIT"; >&2 echo "----> Trim MSA with ClipKIT"
#Steenwyk JL, Buida TJ 3rd, Li Y, Shen XX, Rokas A. ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS Biol. 2020 Dec 2;18(12):e3001007. doi: 10.1371/journal.pbio.3001007. PMID: 33264284; PMCID: PMC7735675.

if [[ ${failedMSA} -eq 0 ]]
then
	clipkit ${outPSeqsAln} -o ${outPSeqsAlnTrim}
	if [[ ${failedBacktrans} -eq 0 ]]
	then
		clipkit ${outTSeqsAln} -o ${outTSeqsAlnTrim}
	else
		echo "Unable to trim MSA for this OG" > ${outTSeqsAlnTrim}
	fi
else
	echo "Unable to trim MSA for this OG" > ${outPSeqsAlnTrim}
	echo "Unable to trim MSA for this OG" > ${outTSeqsAlnTrim}
fi


Plen=$(cat ${outPSeqsAlnTrim} | awk '{if($1 ~ />/){print seq; print $0; seq="";next}seq=seq""$0}END{print seq}' | tail -n +2 | awk '{if($1 ~ />/){next}print length($1)}' | head -1)
Tlen=$(cat ${outTSeqsAlnTrim} | awk '{if($1 ~ />/){print seq; print $0; seq="";next}seq=seq""$0}END{print seq}' | tail -n +2 | awk '{if($1 ~ />/){next}print length($1)}' | head -1)
echo ${Plen}
echo ${Tlen}

if [[ ${Plen}*3 -ne ${Tlen} ]]
then
	echo "----> Length of the trimmed transcript alignment not 3 times the trimmed protein alignment. Need to add gaps at the begining"; >&2 echo "----> Length of the trimmed transcript alignment not 3 times the trimmed protein alignment. Need to add gaps at the begining"
	dif=$((${Tlen}-${Plen} * 3 ))
	echo ${dif}
	if [[ $diff -gt 0 ]]
	then
		toadd=$((3 - $dif % 3))
	else
		toadd=$((- $dif % 3))
	fi
	echo "Added positions ${toadd}"; >&2 echo "Added positions ${toadd}"
	cat ${outTSeqsAlnTrim} | awk '{if($1 ~ />/){print seq; print $0; seq="";next}seq=seq""$0}END{print seq}' | tail -n +2 | awk -v toadd=${toadd} 'BEGIN{for(i=1;i<=toadd;i++){s=s"-"}}{if($1 ~ />/){print $0; next} print $0""s}' > ${outTSeqsAlnTrim}.tmp
	mv ${outTSeqsAlnTrim}.tmp ${outTSeqsAlnTrim}
fi


Tlen=$(cat ${outTSeqsAlnTrim} | awk '{if($1 ~ />/){print seq; print $0; seq="";next}seq=seq""$0}END{print seq}' | tail -n +2 | awk '{if($1 ~ />/){next}print length($1)}' | head -1)
echo ${Tlen}
if [[ ${Tlen}%3 -ne 0 ]]
then
	echo "Error final transcript trimmed alignment length (${Tlen}) not multiple of 3"; >&2 echo "Error final transcript trimmed alignment length (${Tlen}) not multiple of 3"
	exit 1
fi

#############################################
echo "----> Compressing output files"; >&2 echo "----> Compressing output files"
gzip ${outPSeqsAln} ${outPSeqs} ${outPSeqsAlnTrim} ${outTSeqsAln} ${outTSeqs} ${outTSeqsAlnTrim}





