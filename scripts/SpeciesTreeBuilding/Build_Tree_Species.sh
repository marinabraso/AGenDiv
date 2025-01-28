#!/bin/bash

ConcatenatedMSA=$1
outtree=$2
SpTreeConstrainTree=$3
threads=$4

############
# debug / develop
############

rm -r $(dirname  $outtree)
mkdir -p $(dirname  $outtree)
echo "${SpTreeConstrainTree}" > $(dirname  ${outtree})/SpTreeConstrainTree.nwk
zcat $ConcatenatedMSA | sed 's/_Concatenated_One2OneOrthologs//g' >$(dirname  ${outtree})/Concatenate_MSA_One2OneOR.fa
iqtree -redo -nt $threads -s $(dirname  ${outtree})/Concatenate_MSA_One2OneOR.fa -m GTR20 -B 1000 -g $(dirname  $outtree)/SpTreeConstrainTree.nwk --prefix ${outtree%.*}
