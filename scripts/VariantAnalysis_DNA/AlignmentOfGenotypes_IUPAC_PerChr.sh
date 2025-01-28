#!/bin/bash


chrGenotypes=$1
outMSA=${2%.gz}



###################################
### To debug / develop
#chrGenotypes="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.bed.gz"
#outMSA="results/VariantAnalysis_DNA/AlignmentOfGenotypes_IUPAC_PerChr/AlignmentOfGenotypes_IUPAC_{chr}.fa.gz"
###################################

mkdir -p $(dirname ${outMSA})


#zcat ${chrGenotypes} | cut -f5-  | awk '{lenmax=0; for(i=1;i<=NF;i++){split($i, a, ":"); for(j in a){if(length(a[j])>lenmax){lenmax=length(a[j])}}} if(lenmax>1){rstr=""; for(i=1;i<=NF;i++){split($i, a, ":"); rstr=rstr"\t"; for(j in a){dif=lenmax-length(a[j]); astr=a[j]; for(k=1;k<=dif;k++){astr=astr"-"}; rstr=rstr":"astr}} print rstr}else{print $0}}' | sed 's/\t:/\t/g' | sed 's/^\t//g'

# Only SNPs 
zcat ${chrGenotypes} | cut -f5-  | awk '{lenmax=0; gap=0; for(i=1;i<=NF;i++){split($i, a, ":"); for(j in a){if(a[j]=="*"){gap=1} if(length(a[j])>lenmax){lenmax=length(a[j])}}} if(lenmax==1 && gap==0){print $0} delete(al)}' | \
sed 's/A:A/A/g' | \
sed 's/T:T/T/g' | \
sed 's/G:G/G/g' | \
sed 's/C:C/C/g' | \
sed 's/A:G/R/g' | sed 's/G:A/R/g' | \
sed 's/C:T/Y/g' | sed 's/T:C/Y/g' | \
sed 's/G:C/S/g' | sed 's/C:G/S/g' | \
sed 's/A:T/W/g' | sed 's/T:A/W/g' | \
sed 's/G:T/K/g' | sed 's/T:G/K/g' | \
sed 's/A:C/M/g' | sed 's/C:A/M/g' | grep ':' | more 



A	Adenine
C	Cytosine
G	Guanine
T (or U)	Thymine (or Uracil)
R	A or G
Y	C or T
S	G or C
W	A or T
K	G or T
M	A or C
B	C or G or T
D	A or G or T
H	A or C or T
V	A or C or G
N	any base
. or -	gap




