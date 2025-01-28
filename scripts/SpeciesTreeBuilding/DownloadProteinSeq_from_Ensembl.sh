#!/bin/bash


proteomes=$1
transcriptomes=$2
species=$3
spVersions=$4
EnsRelease=$5
EnsMetazoaRelease=$6

############
# debug / develop
#proteomes="data/ProteinSequences/Homo_sapiens.GRCh38_proteins.fa.gz data/ProteinSequences/Mus_musculus.GRCm39_proteins.fa.gz data/ProteinSequences/Danio_rerio.GRCz11_proteins.fa.gz"
#species="Homo_sapiens Mus_musculus Danio_rerio"
#spVersions="Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11"
#EnsRelease=111
############

speciesArray=( $(echo ${species} | sed 's/ /\n/g') )
proteomesArray=( $(echo ${proteomes} | sed 's/ /\n/g') )
transcriptomesArray=( $(echo ${transcriptomes} | sed 's/ /\n/g') )
spVersionsArray=( $(echo ${spVersions} | sed 's/ /\n/g') )

pwd=$(pwd)
rm -r $(dirname ${proteomesArray[0]}) 2> ~/null
mkdir -p $(dirname ${proteomesArray[0]})
cd $(dirname ${proteomesArray[0]})

for s in ${!speciesArray[@]}
do
	echo "----> ${speciesArray[$s]} ${proteomesArray[$s]}"; >&2 echo "----> ${speciesArray[$s]} ${proteomesArray[$s]}"
	version=$(echo ${spVersionsArray[$s]} | cut -f2- -d'.')
	speciesInVersion=$(echo ${spVersionsArray[$s]} | cut -f1 -d'.')
	echo "         ${version} ${speciesInVersion}"; >&2 echo "         ${version} ${speciesInVersion}"
	if [ $speciesInVersion == "Strongylocentrotus_purpuratus" ]
	then
		wget "http://ftp.ensemblgenomes.org/pub/metazoa/release-${EnsMetazoaRelease}/fasta/${speciesInVersion,,}/pep/${speciesInVersion}.${version}.pep.all.fa.gz"
		wget "http://ftp.ensemblgenomes.org/pub/metazoa/release-${EnsMetazoaRelease}/fasta/${speciesInVersion,,}/cds/${speciesInVersion}.${version}.cds.all.fa.gz"
		mv ${speciesInVersion}.${version}.pep.all.fa.gz ${pwd}/${proteomesArray[$s]}
		mv ${speciesInVersion}.${version}.cds.all.fa.gz ${pwd}/${transcriptomesArray[$s]}
	elif [ $speciesInVersion == "Branchiostoma_lanceolatum" ]
	then
		wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/927/797/965/GCA_927797965.1_BraLan3/GCA_927797965.1_BraLan3_translated_cds.faa.gz"
		wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/927/797/965/GCA_927797965.1_BraLan3/GCA_927797965.1_BraLan3_cds_from_genomic.fna.gz"
		mv GCA_927797965.1_BraLan3_translated_cds.faa.gz ${pwd}/${proteomesArray[$s]}
		mv GCA_927797965.1_BraLan3_cds_from_genomic.fna.gz ${pwd}/${transcriptomesArray[$s]}
	elif [ $speciesInVersion == "Branchiostoma_floridae" ]
	then
		wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/815/GCF_000003815.2_Bfl_VNyyK/GCF_000003815.2_Bfl_VNyyK_translated_cds.faa.gz"
		wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/815/GCF_000003815.2_Bfl_VNyyK/GCF_000003815.2_Bfl_VNyyK_cds_from_genomic.fna.gz"
		mv GCF_000003815.2_Bfl_VNyyK_translated_cds.faa.gz ${pwd}/${proteomesArray[$s]}
		mv GCF_000003815.2_Bfl_VNyyK_cds_from_genomic.fna.gz ${pwd}/${transcriptomesArray[$s]}
	else
		wget "https://ftp.ensembl.org/pub/release-${EnsRelease}/fasta/${speciesInVersion,,}/pep/${speciesInVersion}.${version}.pep.all.fa.gz"
		wget "https://ftp.ensembl.org/pub/release-${EnsRelease}/fasta/${speciesInVersion,,}/cds/${speciesInVersion}.${version}.cds.all.fa.gz"
		mv ${speciesInVersion}.${version}.pep.all.fa.gz ${pwd}/${proteomesArray[$s]}
		mv ${speciesInVersion}.${version}.cds.all.fa.gz ${pwd}/${transcriptomesArray[$s]}
	fi
done

cd ${pwd}












