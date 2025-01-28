#!/bin/bash


ChrLengths=$1
GTF=$2
CallableRegions=$3
OutGenomicFeatures=$4
OutGenomicFeaturesCallab=$5

##############
# To debug / develop
#ChrLenghts="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#GTF="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf"
#CallableRegions="results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz"
#OutGenomicFeatures="results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures.bed"
#OutGenomicFeaturesCallab="results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_CallableRegions.bed"
##############



mkdir -p $(dirname ${OutGenomicFeatures})

OutBasename=$(dirname ${OutGenomicFeatures})/$(basename ${OutGenomicFeatures%_*})

echo "----> Per gene merged exons"
cat ${GTF} | awk '{if($3 == "exon"){print $0}}' | cut -f1,4,5,9 | sed 's/gene_id "\(\S\+\)";.*/\1/g' | sort -k4,4 -k2,2V | uniq | awk '{if(NR==1){g=$4;chr=$1;st=$2;end=$3;next;}if($4!=g){print chr"\t"st"\t"end"\t"g; g=$4;chr=$1;st=$2;end=$3;next;} if($2<end){end=$3}else{print chr"\t"st"\t"end"\t"g;st=$2;end=$3;}}' > ${OutBasename}_exons.bed
echo "----> Introns"
cat ${OutBasename}_exons.bed | awk '{if(NR==1 || g!=$4){g=$4;st=$3+1;next;} end=$2-1; if(end>st){print $1"\t"st"\t"end"\t"g;} st=$3+1}' > ${OutBasename}_introns.bed
echo "----> Promoters"
cat ${OutBasename}_exons.bed | awk '{if(NR==1 || g!=$4){g=$4;st=$2-1000;end=$2-1; if(st<0){st=0};if(end>0){print $1"\t"st"\t"end"\t"g}}}' > ${OutBasename}_promoters.bed
echo "----> Merge all genic"
bedtools merge -i <(cat ${OutBasename}_exons.bed ${OutBasename}_introns.bed ${OutBasename}_promoters.bed | sort -k1,1 -k2,2V) > ${OutBasename}_genic.bed
echo "----> Ouput Intergenic from the negtive of the genic"
awk '{if(NR==FNR){a[$1]=$2;next}if(FNR==1){end=$2-1;if(end>0){print $1"\t0\t"end;}chr=$1;st=$3+1;next}if(chr!=$1){print chr"\t"st"\t"a[">"chr];end=$2-1;if(end>0){print $1"\t0\t"end;}chr=$1;st=$3+1;next}end=$2-1;if(end>st){print chr"\t"st"\t"end;}st=$3+1}END{print chr"\t"st"\t"a[">"chr];}' ${ChrLengths} ${OutBasename}_genic.bed > ${OutBasename}_intergenic.bed
echo "----> Join all types, select only chromosomes & asign an ID to all of them"
awk 'BEGIN{str="exon intron promoter intergenic"; split(str, a, " "); ninterg=0}{if(FNR==1){file++} if(file==4){ninterg++;id="Interg"ninterg}else{id=$4} print $1"\t"$2"\t"$3"\t"a[file]"\t"id}' ${OutBasename}_exons.bed ${OutBasename}_introns.bed ${OutBasename}_promoters.bed ${OutBasename}_intergenic.bed | grep '^chr' | sed 's/\t$/\t-/g' > ${OutGenomicFeatures%.bed}.tmp
cat ${OutGenomicFeatures%.bed}.tmp | sort -k5,5 -k2,2V | awk '{if($4=="exon"){if(ge==$5){idxe++;}else{ge=$5;idxe=1}print $0"_e"idxe}else if($4=="intron"){if(gi==$5){idxi++;}else{gi=$5;idxi=1}print $0"_i"idxi}else if($4=="promoter"){print $0"_p"}else{print $0}}' | sort -k1,1V -k2,2V > ${OutGenomicFeatures}
echo "----> Intersect with callable regions"
bedtools intersect -a ${OutGenomicFeatures} -b <(zcat ${CallableRegions} | cut -f1,2,3 | sort -k1,1 -k2,2V) > ${OutGenomicFeaturesCallab}




