#!/bin/bash

OGs=$1
OutDir=$2
species=$3
SpeciesGenePrefix=$4

############
# debug / develop
#OGs=results/OrthogroupsAndGeneTrees/BuildOGs_withOrthoLoger/OGs.tab.gz
#OutDir=results/OrthogroupsAndGeneTrees/FirstFilterOG
#species="Pan_troglodytes Homo_sapiens Gorilla_gorilla Macaca_fascicularis Macaca_mulatta Callithrix_jacchus Cavia_porcellus Heterocephalus_glaber Mus_musculus Rattus_norvegicus Oryctolagus_cuniculus Ovis_aries Bos_taurus Sus_scrofa Canis_lupus_familiaris Monodelphis_domestica Ornithorhynchus_anatinus Gallus_gallus Xenopus_tropicalis Danio_rerio"
#SpeciesGenePrefix="Pan_troglodytes.ENSPTR Homo_sapiens.ENS Gorilla_gorilla.ENSGGO Macaca_fascicularis.ENSMFA Macaca_mulatta.ENSMMU Callithrix_jacchus.ENSCJA Cavia_porcellus.ENSCPO Heterocephalus_glaber.ENSHGL Mus_musculus.ENSMUS Rattus_norvegicus.ENSRNO Oryctolagus_cuniculus.ENSOCU Ovis_aries.ENSOAR Bos_taurus.ENSBTA Sus_scrofa.ENSSSC Canis_lupus_familiaris.ENSCAF Monodelphis_domestica.ENSMOD Ornithorhynchus_anatinus.ENSOAN Gallus_gallus.ENSGAL Xenopus_tropicalis.ENSXET Danio_rerio.ENSDAR"
############




mkdir -p ${OutDir}
rm ${OutDir}/* 2> ~/null
numSpecies=$(echo ${species} | sed 's/ /\n/g' | wc -l)
OGsIDs=$(zcat ${OGs} | grep -v '^#' | cut -f1 -d' ' | sort | uniq)

for OG in ${OGsIDs}
do
	OutputFile=${OutDir}"/OG"${OG}".tbl"
	rm ${OutputFile} 2> ~/null
	numberGenes=$(zcat ${OGs} | grep -v '^#' | grep -w "^${OG}" | wc -l)
	IDs=$(zcat ${OGs} | grep -v '^#' | grep -w "^${OG}" | cut -f2 -d' ')
	ListSp=""
	if [[ ${numberGenes} -eq ${numSpecies} ]]
	then
		echo "Correct number of genes in "$OG
		for g in ${IDs}
		do 
			# Get species
			gprefix=$(echo ${g:0:10} | sed 's/G*[0-9]\+//g')
			sp=$(echo ${SpeciesGenePrefix} | sed 's/ /\n/g' | sed 's/\./\t/g' | grep -w ${gprefix} | cut -f1)
			ListSp=$(echo ${ListSp}" "${sp})
		done
		uniqueSpnumber=$(echo ${ListSp} | sed 's/ /\n/g' | sort | uniq | wc -l)
		if [[ ${uniqueSpnumber} -eq ${numSpecies} ]]
		then
			echo "Correct number of unique species in "$OG
			for g in ${IDs}
			do 
				echo $g
				# Get species
				gprefix=$(echo ${g:0:10} | sed 's/G*[0-9]\+//g')
				sp=$(echo ${SpeciesGenePrefix} | sed 's/ /\n/g' | sed 's/\./\t/g' | grep -w ${gprefix} | cut -f1)
				echo ${g}"_"${sp} >> ${OutputFile}
			done
			gzip ${OutputFile}
		fi
	fi
done




