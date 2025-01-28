#!/bin/bash


proteomes=$1
OGsOutput=$2
orthologer_dir=$3
brhclus_dir=$4
script_prompt=$5
species=$6

############
# debug / develop
############

speciesArray=( $(echo ${species} | sed 's/ /\n/g') )
proteomesArray=( $(echo ${proteomes} | sed 's/ /\n/g') )
pwd=$(pwd)
OLinputFile=path2proteome.tsv ## needs to be named like this
rm -r $(dirname ${OGsOutput}) 2> ~/.null

############
# create working directory
echo "----> create working directory"; >&2 echo "----> create working directory"
mkdir -p $(dirname ${OGsOutput})
cd $(dirname ${OGsOutput})

############
# copy FASTA files into data directory
echo "----> copy FASTA files into data directory"; >&2 echo "----> copy FASTA files into data directory"
for s in ${!speciesArray[@]}
do
	echo ${speciesArray[$s]} ${pwd}/${proteomesArray[$s]}
	gzip -dc ${pwd}/${proteomesArray[$s]} > ${speciesArray[$s]^^}.fs
	echo "+"${speciesArray[$s]^^}" "${speciesArray[$s]^^}".fs" >> ${OLinputFile}
done

# Create the whole project structure
echo "----> Create the whole project structure"; >&2 echo "----> Create the whole project structure"
${pwd}/${orthologer_dir}/manage_project.sh -s -f ${OLinputFile}

# Prepare execution with run_loger_prompt.sh (${script_prompt})
echo "----> Prepare execution with run_loger_prompt.sh (${script_prompt})"; >&2 echo "----> Prepare execution with run_loger_prompt.sh (${script_prompt})"
${pwd}/${script_prompt} ${pwd}/${brhclus_dir} # runs sbin/setup.sh go and takes care of promt required inputs

# Setup project with setup exec created by OrthoLoger
echo "----> Setup project with setup exec created by OrthoLoger"; >&2 echo "----> Setup project with setup exec created by OrthoLoger"
#mv setup_project_*.sh setup_project.sh # change the defauld exec name created by OrthoLoger that contains the user name
./setup_project_*.sh

# Change config file
echo "----> Change config file"; >&2 echo "----> Change config file"
PYTHON3_PATH=$(which python3)
sed -i "s|PYTHON3=\".*\"|PYTHON3=\"$PYTHON3_PATH\"|" orthologer_conf.sh

sed -i "s|DIR_BRHCLUS=\"\"|DIR_BRHCLUS=\"${pwd}/${brhclus_dir}\"|" orthologer_conf.sh

CONDA_PATH=$(dirname ${PYTHON3_PATH})
sed -i "s|DIR_BLAST=\".*\"|DIR_BLAST=\"$CONDA_PATH\"|" orthologer_conf.sh
sed -i "s|DIR_SWIPE=\".*\"|DIR_SWIPE=\"$CONDA_PATH\"|" orthologer_conf.sh
sed -i "s|DIR_DIAMOND=\".*\"|DIR_DIAMOND=\"$CONDA_PATH\"|" orthologer_conf.sh
sed -i "s|DIR_MMSEQS=\".*\"|DIR_MMSEQS=\"$CONDA_PATH\"|" orthologer_conf.sh
sed -i "s|DIR_CDHIT=\".*\"|DIR_CDHIT=\"$CONDA_PATH\"|" orthologer_conf.sh

# Check the pipeline setup
echo "----> Check the pipeline setup"; >&2 echo "----> Check the pipeline setup"
./orthologer.sh -xp

# Run OrthoLoger
echo "----> Run OrthoLoger"; >&2 echo "----> Run OrthoLoger"
./orthologer.sh -r all

echo "Done"; >&2 echo "Done"
cd ${pwd}

echo "Copying output to correct path & zipping it"; >&2 echo "Copying output to correct path & zipping it"
cat $(dirname ${OGsOutput})/Results/path2proteome_orthogroups.txt | gzip > ${OGsOutput}


## Output format

# <cluster id> <gene label> <type> <seq len> <seq start> <seq end> <pid> <score> <e-value>
# dThe type is a flag that indicates how the gene was clustered.
# 0   - part of a triangle (BRHs)
# 2   - BRH with a type 0 without being a type 0 itself
# 4,5 - same as type 2 but with only a fraction of the sequence aligned
# 1   - a pair (single BRH)
# 3   - a chain of pairs
# 7   - inparalog
# 9   - 97% id - genes that were clustered via a single representative



