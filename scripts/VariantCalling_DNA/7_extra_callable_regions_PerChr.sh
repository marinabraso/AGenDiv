#!/bin/bash

CallableRegions=$1
ChrLengths=$2
ExtraCallableRegions=$3
NonExtraCallableRegions=$4
ValidWindows=$5
window=$6
step=$7
mincallable=$8
chr=$9


############################################
## To debug
#CallableRegions="results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz"
#ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#ExtraCallableRegions="results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ExtraCallableRegions_chr19.bed.gz"
#NonExtraCallableRegions="results/VariantCalling_DNA/7_extra_callable_regions_PerChr/NonExtraCallableRegions_chr19.bed.gz"
#ValidWindows="results/VariantCalling_DNA/7_join_extra_callable_regions/ValidWindows_chr19.bed.gz"
#window=100000 # window size for sliding window  
#step=10000 # step size for sliding window 
#mincallable=0.3 # minimum proportion of the sliding window being callable to be accepted as really callable
#chr="chr19"
############################################




chrs=( $(echo ${chrs} | sed 's/ /\n/g') )
chrlen=$(grep -w "^>${chr}" ${ChrLengths} | cut -f2)
mkdir -p $(dirname ${ExtraCallableRegions})
chrCallableRegions=$(dirname ${ExtraCallableRegions})/CallableRegions_${chr}.bed
zcat ${CallableRegions} | grep ${chr} > ${chrCallableRegions}

echo "----> Build window bed file"
if [[ ! -s $(dirname ${ExtraCallableRegions})"/Windows_"${chr}".bed" ]]
then
	for ((w=0; w<=${chrlen}-step; w+=step))
	do
	end=$(( ${w}+${window} ))
	echo "${chr} ${w} ${end}" >> $(dirname ${ExtraCallableRegions})"/Windows_"${chr}".bed"
	done
	sed -i 's/ /\t/g' $(dirname ${ExtraCallableRegions})"/Windows_"${chr}".bed"
fi

echo "----> Bedtools intersect window bed file with CallableRegions & keep only windows that fulfill the condition of a minimal proportion of their lengths being callable"
bedtools intersect -a $(dirname ${ExtraCallableRegions})"/Windows_"${chr}".bed" -b <(cat ${chrCallableRegions}) -wao | awk '{ov=$NF+1; print $1"\t"$2"\t"$3"\t"ov}' | awk '{if(NR==1){w=$1"\t"$2"\t"$3; sum=$4;next}if(w==$1"\t"$2"\t"$3){sum+=$4}else{print w"\t"sum; w=$1"\t"$2"\t"$3; sum=$4;}}END{print w"\t"sum;}' | awk -v minpcal=${mincallable} '{pcal=$4/($3-$2+1);if(pcal>=minpcal){print $1"\t"$2"\t"$3}}' > ${ValidWindows%.gz}

echo "----> Bedtools intersect valid windows with callable regions, keep all callable regions that overlap with at least one valid window"
bedtools intersect -a <(cat ${chrCallableRegions}) -b <(bedtools merge -i ${ValidWindows%.gz}) -wa > ${ExtraCallableRegions%.gz}

echo "----> build the complementary regions of ExtraCallableRegions to get NonExtraCallableRegions"
# substracting one bp to each end value to fit gatk VariantFiltration comprehension of non-valid regions (it does not consider the start bp but it does cosider the end pb of each non-valid region)
cat ${ExtraCallableRegions%.gz} | awk -v endchr=${chrlen} 'BEGIN{st=0}{end=$2-1;print $1"\t"st"\t"end; st=$3}END{print $1"\t"st"\t"endchr}' > ${NonExtraCallableRegions%.gz}

gzip ${ExtraCallableRegions%.gz} ${NonExtraCallableRegions%.gz} ${ValidWindows%.gz}
rm $(dirname ${ExtraCallableRegions})"/Windows_"${chr}".bed" ${chrCallableRegions} 2> ~/null








