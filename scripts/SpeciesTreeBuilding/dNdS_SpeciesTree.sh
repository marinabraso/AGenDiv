#!/bin/bash

inMSAfa=$1
outMSAphy=${2%.gz}
outCodemlOut=${3%.out}
Tree=$4


############
# debug / develop
############

rm -r $(dirname  $outCodemlOut) 2> ~/.null
mkdir -p $(dirname  $outCodemlOut)

# Modif species on names on tree for codeml
### Add header!
numspecies=$(echo ${Tree} | sed 's/[();, \#1]/ /g' | awk '{print NF}')
echo "${numspecies} 1" > $(dirname  $outCodemlOut)/SpeciesTree_base.tree
echo ${Tree} | sed 's/\([A-Z]\)[a-z]*_\([a-z][a-z][a-z]\)[a-z]*/\1\2/g' >> $(dirname  $outCodemlOut)/SpeciesTree_base.tree
#cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree
#cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree | sed 's/Drer/Drer\#1/g' | sed 's/Hsap/Hsap\#1/g' | sed 's/Mmus)/Mmus\#1)\#1/g' > $(dirname  $outCodemlOut)/SpeciesTree_foregroundVertebrates.tree
#cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree | sed 's/Blan/Blan\#1/g' | sed 's/Bflo/Bflo\#1/g' > $(dirname  $outCodemlOut)/SpeciesTree_foregroundAmphioxus.tree
#cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree | sed 's/Cint/Cint\#1/g' > $(dirname  $outCodemlOut)/SpeciesTree_foregroundTunicates.tree
speciesorderintree=$(cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree | sed 's/[();, \#1]/ /g')

# For MSA
# fasta to phylip
# Modif species on names on alignment for codeml
zcat $inMSAfa | sed 's/\./-/g' | sed 's/>\([A-Z]\)[a-z]*_\([a-z][a-z][a-z]\).*/>\1\2/g' | awk '{if($1 ~ />/){seq=""; num++}else{seq=seq$1}}END{print num" "length(seq)}' > $outMSAphy 
zcat $inMSAfa | sed 's/\./-/g' | sed 's/>\([A-Z]\)[a-z]*_\([a-z][a-z][a-z]\).*/>\1\2/g' | awk '{if($1 ~ />/){print name"  "seq; name=$1; seq=""}else{seq=seq$1}}END{print name"  "seq;}' | sed 's/>//g' | tail -n +2 > $outMSAphy.tmp
for s in $speciesorderintree
do
    cat $outMSAphy.tmp | awk -v s=$s '{if($1 == s){print $0}}' >> $outMSAphy
done
rm $outMSAphy.tmp

########################
# Prepare basic control file of codeml
########################

## Paths
echo "seqfile = ${outMSAphy}" > $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl
echo "treefile = $(dirname  ${outCodemlOut})/SpeciesTree_base.tree" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl
echo "outfile = ${outCodemlOut}.out" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl

## Output modes
echo "noisy = 9" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0,1,2,3,9: how much rubbish on the screen
echo "verbose = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 1: detailed output, 0: concise output
echo "runmode = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: user tree; 1: semi-automatic; 2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

## Specifics of the execution
echo "seqtype = 1" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl  # 1:codons; 2:AAs; 3:codons-->AAs
echo "CodonFreq = 2" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl  # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table; ndata = 10
echo "clock = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0:no clock, 1:clock; 2:local clock
echo "aaDist = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a; 7:AAClasses
jonespath=$(find $CONDA_DEFAULT_ENV -name jones.dat)
echo "aaRatefile = ${jonespath}" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # only used for aa seqs with model=empirical(_F), dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
echo "model = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl
# models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches
# models for AAs or codon-translated AAs: 0:poisson, 1:proportional,2:Empirical,3:Empirical+F,5:FromCodon0, 6:FromCodon1, 8:REVaa_0, 9:REVaa(nr=189)
echo "NSsites = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;13:3normal>0
echo "icode = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0:universal code; 1:mammalian mt; 2-11:see belowPAML M ANUAL
echo "Mgene = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0:rates, 1:separate;
echo "fix_kappa = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 1: kappa fixed, 0: kappa to be estimated
echo "kappa = 2" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # initial or fixed kappa
echo "fix_omega = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 1: omega or omega_1 fixed, 0: estimate
echo "omega = .4" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # initial or fixed omega, for codons or codon-based AAs
echo "fix_alpha = 1" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: estimate gamma shape parameter; 1: fix it at alpha
echo "alpha = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # initial or fixed alpha, 0:infinity (constant rate)
echo "Malpha = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # different alphas for genes
echo "ncatG = 3" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # # of categories in dG of NSsites models
echo "fix_rho = 1" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: estimate rho; 1: fix it at rho
echo "rho = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # initial or fixed rho, 0:no correlation
echo "getSE = 1" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: don't want them, 1: want S.E.s of estimates
echo "RateAncestor = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
echo "Small_Diff = .5e-6" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl
echo "cleandata = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # remove sites with ambiguity data (1:yes, 0:no)?
echo "fix_blength = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
echo "method = 0" >> $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl # 0: simultaneous; 1: one branch at a time

# Run codeml
#echo "----> Vertebrate branch"; >&2 echo "----> Vertebrate branch"
#echo "          ----> Tree"
#cat $(dirname  $outCodemlOut)/SpeciesTree_base.tree
#echo "          ----> Seq"
#cat ${outMSAphy} | cut -c1-99 > ${outMSAphy}.tmp
#cat ${outMSAphy}.tmp | sed 's/8 [0-9]\+/8 93/g' > ${outMSAphy}
#cat ${outMSAphy}
#echo "          ----> Control"
#cat $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl


echo -e "\n\n----> Run codeml model=0"; >&2 echo "----> Run codeml model=0"
codeml $(dirname  $outCodemlOut)/dNdS_SpeciesTree.ctl


gzip $outMSAphy
