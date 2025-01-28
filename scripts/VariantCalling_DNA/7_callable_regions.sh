#!/bin/bash

SNPableQuality=$1
mincov=$2
CallableRegions=$3
NonCallableRegions=$4
lcovthreshold=$5
ucovthreshold=$6
chrs=$7


############################################
## To debug
#SNPableQuality="results/VariantCalling_DNA/7_SNPable_Regions/Branchiostoma_lanceolatum.BraLan3_mask_35_50.fa"
#CallableRegions="results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz"
#NonCallableRegions="results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz"
#mincov="results/VariantCalling_DNA/7_join_coverage_per_site/Depth.MinCov.tbl.gz"
#lcovthreshold=5
#ucovthreshold=20
#chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
############################################




chrs=( $(echo ${chrs} | sed 's/ /\n/g') )
mkdir -p $(dirname ${CallableRegions}) tmp


# SNPableQuality from fa to tbl
echo "### fa2tbl"; >&2 echo "### fa2tbl"
#cat ${SNPableQuality} | awk '{if($1 ~ />/){chr=$1;num=0;next}split($1, a, "");for(i in a){num++;print chr"\t"num"\t"a[i]}}' | sed 's/>//g' | grep '^chr' | sort -k1,2V > $(dirname ${CallableRegions})/$(basename ${SNPableQuality}).tbl
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

# Paste SNPableQuality and mincov
echo "### Checking number of lines before paste"; >&2 echo "### Checking number of lines before paste"
nmincov=$(zcat ${mincov} | wc -l)
nSNPable=$(cat $(dirname ${CallableRegions})/$(basename ${SNPableQuality}).tbl | wc -l)
if [[ ${nmincov} -ne ${nSNPable} ]]
then
	exit "Error in number of lines before paste (nmincov=${nmincov} & nSNPable=${nSNPable})"
fi
echo "### Paste"; >&2 echo "### Paste"
#paste $(dirname ${CallableRegions})/$(basename ${SNPableQuality}).tbl <(zcat ${mincov} | grep '^chr' | sort -k1,2V | cut -f3) | gzip > $(dirname ${CallableRegions})/SNPable_MinCov.tbl.gz
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### tbl2bed & apply SNPable & cov filters"; >&2 echo "### tbl2bed & apply SNPable & cov filters"
zcat $(dirname ${CallableRegions})/SNPable_MinCov.tbl.gz | \
	awk -v lct=$lcovthreshold -v uct=$ucovthreshold '{if($3==3 && $4>=lct && $4<uct){print $1"\t"$2"\t1"}else{print $1"\t"$2"\t0"}}' | \
	awk '{if(NR==1){chr=$1;st=$2;end=$2;len=1;io=$3;next}if($3!=io || $1!=chr){print chr"\t"st"\t"end"\t"io"\t"len;chr=$1;st=$2;end=$2;len=1;io=$3;next}end=$2;len++}END{print chr"\t"st"\t"end"\t"io"\t"len;}' | gzip > $(dirname ${CallableRegions})/AllWindows_SNPable_coverage.bed.gz
dstatus=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### tbl2bed applying only SNPable filter"; >&2 echo "### tbl2bed applying only SNPable filter"
zcat $(dirname ${CallableRegions})/SNPable_MinCov.tbl.gz | \
	awk '{if($3==3){print $1"\t"$2"\t1"}else{print $1"\t"$2"\t0"}}' | \
	awk '{if(NR==1){chr=$1;st=$2;end=$2;len=1;io=$3;next}if($3!=io || $1!=chr){print chr"\t"st"\t"end"\t"io"\t"len;chr=$1;st=$2;end=$2;len=1;io=$3;next}end=$2;len++}END{print chr"\t"st"\t"end"\t"io"\t"len;}' | awk '{if($4 == 1){print $1"\t"$2"\t"$3"\t"$5}}' | gzip > $(dirname ${CallableRegions})/SNPableRegions.bed.gz
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"


echo "### only callable"; >&2 echo "### only callable"
zcat $(dirname ${CallableRegions})/AllWindows_SNPable_coverage.bed.gz | awk '{if($4 == 1){print $1"\t"$2"\t"$3"\t"$5}}' | gzip > ${CallableRegions}
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"


echo "### only non-callable"; >&2 echo "### only non-callable"
# substracting one bp to each start to avoid gatk VariantFiltration to include variants in the first bp of the noncallable region
zcat $(dirname ${CallableRegions})/AllWindows_SNPable_coverage.bed.gz | awk '{if($4 == 0){st=$2-1;print $1"\t"st"\t"$3"\t"$5}}' | gzip > ${NonCallableRegions}
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

