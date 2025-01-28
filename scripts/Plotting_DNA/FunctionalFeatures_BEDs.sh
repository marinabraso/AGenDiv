#!/bin/bash


ChrLengths=$1
GTF=$2
varBED=$3
exons=$4
introns=$5
promoters=$6
intergenic=$7
PerSite_FeatureType=$8

##############
# To debug / develop
#ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#GTF="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf.gz"
#exons="results/Plotting_DNA/FunctionalFeatures_BEDs/exons.bed.gz"
#introns="results/Plotting_DNA/FunctionalFeatures_BEDs/introns.bed.gz"
#promoters="results/Plotting_DNA/FunctionalFeatures_BEDs/promoters.bed.gz"
#intergenic="results/Plotting_DNA/FunctionalFeatures_BEDs/intergenic.bed.gz"
#PerSite_FeatureType="results/Plotting_DNA/FunctionalFeatures_BEDs/PerSite_VariantType.txt.gz"
##############

mkdir -p $(dirname ${exons})

echo "----> Per gene merged exons"
zcat < ${GTF} | grep '^chr' | awk '{if($3 == "exon"){print $0}}' | cut -f1,4,5,9 | sed 's/gene_id "\([A-Z0-9]*\)";.*/\1/g' | sort -k4,4 -k2,2V | uniq | awk '{if(NR==1){g=$4;chr=$1;st=$2;end=$3;next;}if($4!=g){print chr"\t"st"\t"end"\t"g; g=$4;chr=$1;st=$2;end=$3;next;} if($2<end){end=$3}else{print chr"\t"st"\t"end"\t"g;st=$2;end=$3;}}' | gzip > ${exons}

echo "----> Introns"
zcat < ${exons} | awk '{if(NR==1 || g!=$4){g=$4;st=$3+1;next;} end=$2-1; if(end>st){print $1"\t"st"\t"end"\t"g;} st=$3+1}' | gzip > ${introns}

echo "----> Promoters"
zcat < ${exons} | awk '{if(NR==1 || g!=$4){g=$4;st=$2-1000;end=$2-1; if(st<0){st=0};if(end>0){print $1"\t"st"\t"end"\t"g}}}' | gzip > ${promoters}

echo "----> Merge all genic"
bedtools merge -i <(for i in ${exons} ${introns} ${promoters};do zcat < $i; done | sort -k1,1 -k2,2V) > $(dirname ${exons})/genic.tmp.bed

echo "----> Ouput Intergenic from the negtive of the genic"
awk '{if(NR==FNR){a[$1]=$2;next}if(FNR==1){end=$2-1;if($2>1){print $1"\t1\t"end;}chr=$1;st=$3+1;next}if(chr!=$1){print chr"\t"st"\t"a[">"chr];end=$2-1;if($2>1){print $1"\t1\t"end;}chr=$1;st=$3+1;next}end=$2-1;if(end>=st){print chr"\t"st"\t"end;}st=$3+1}END{print chr"\t"st"\t"a[">"chr];}' <(cat ${ChrLengths} | grep '^>chr') $(dirname ${exons})/genic.tmp.bed | gzip > ${intergenic}

echo "----> Join all types & asign an ID to all of them"
awk 'BEGIN{str="exon intron promoter intergenic"; split(str, a, " "); ninterg=0}{if(FNR==1){file++} if(file==4){ninterg++;id="Interg"ninterg}else{id=$4} print $1"\t"$2"\t"$3"\t"a[file]"\t"id}' <(zcat < ${exons}) <(zcat < ${introns}) <(zcat < ${promoters}) <(zcat < ${intergenic}) | sed 's/\t$/\t-/g' > $(dirname ${exons})/tmp.bed
cat $(dirname ${exons})/tmp.bed | sort -k5,5 -k2,2V | awk '{if($4=="exon"){if(ge==$5){idxe++;}else{ge=$5;idxe=1}print $0"_e"idxe}else if($4=="intron"){if(gi==$5){idxi++;}else{gi=$5;idxi=1}print $0"_i"idxi}else if($4=="promoter"){print $0"_p"}else{print $0}}' | sort -k1,1V -k2,2V > ${PerSite_FeatureType}.tmp

bedtools intersect -a <(zcat ${varBED} | cut -f1,2,3) -b ${PerSite_FeatureType}.tmp -wo | cut -f1,2,3,7,8 | awk '{if(var!=$1"\t"$2"\t"$3){print var"\t"type"\t"ids; var=$1"\t"$2"\t"$3; type=$4; ids=$5;next} type=type","$4; ids=ids","$5}END{print var"\t"type"\t"ids;}' | tail -n +2 | gzip > ${PerSite_FeatureType}

#rm $(dirname ${exons})/tmp.bed $(dirname ${exons})/genic.tmp.bed ${PerSite_FeatureType}.tmp 2> ~/null 











