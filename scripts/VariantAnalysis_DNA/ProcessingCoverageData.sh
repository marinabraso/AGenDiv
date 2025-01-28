#!/bin/bash

Coverage=$1
ChrLengths=$2
OutBasename=${3%.bed.gz}
window=$4
step=$5
chr=$6

###################################
### To debug
#Coverage="results/VariantCalling_DNA/7_join_coverage_per_site/Depth.allsamples.tbl.gz"
#ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#OutBasename="results/VariantAnalysis_DNA/ProcessingCoverageData/CoverageWindowStats_chr18"
#window=100000
#step=10000
#chr="chr18"
###################################

mkdir -p $(dirname ${OutBasename})

chrlen=$(grep -w "^>${chr}" ${ChrLengths} | cut -f2)

echo "----> Build window bed file"
if [[ ! -s $(dirname ${OutBasename})"/"${chr}"_windows.bed" ]]
then
	for ((w=0; w<=${chrlen}-step; w+=step))
	do
	end=$(( ${w}+${window} ))
	echo "${chr} ${w} ${end}" >> $(dirname ${OutBasename})"/"${chr}"_windows.bed"
	done
	sed -i 's/ /\t/g' $(dirname ${OutBasename})"/"${chr}"_windows.bed"
fi

echo "----> Calculate mincoverage, maxcoverage, average coverage and coverage standard deviation per site in the chr"
zcat ${Coverage} | grep ${chr} | awk '{min=1000;max=0;sum=0;for(i=3;i<=NF;i++){if($i<min){min=$i}if($i>max){max=$i}sum+=$i}av=sum/(NF-2);var=0; for(i=3;i<=NF;i++){var+=($i-av)*($i-av)} sd=sqrt(var/(NF-2)); print $1"\t"$2"\t"$2"\t"min"\t"max"\t"av"\t"sd}' | gzip > $(dirname ${OutBasename})"/"${chr}"_stats.bed.gz"

echo "----> Intersect per-base stats with windows and do the average for each window"
bedtools intersect -a $(dirname ${OutBasename})"/"${chr}"_windows.bed" -b <(zcat $(dirname ${OutBasename})"/"${chr}"_stats.bed.gz") -wao | cut -f1,2,3,7,8,9,10 | awk '{if(w == $1"\t"$2"\t"$3){a["min"]+=$4;a["max"]+=$5;a["av"]+=$6;a["sd"]+=$7;a["num"]+=1;next}if(NR!=1){print w"\t"a["min"]/a["num"]"\t"a["max"]/a["num"]"\t"a["av"]/a["num"]"\t"a["sd"]/a["num"];} w = $1"\t"$2"\t"$3; a["min"]=$4;a["max"]=$5;a["av"]=$6;a["sd"]=$7;a["num"]=1}END{print w"\t"a["min"]/a["num"]"\t"a["max"]/a["num"]"\t"a["av"]/a["num"]"\t"a["sd"]/a["num"];}' | gzip > ${OutBasename}.bed.gz



rm $(dirname ${OutBasename})"/"${chr}"_windows.bed"  2> ~/null
rm $(dirname ${OutBasename})"/"${chr}"_stats.bed.gz"  2> ~/null













