#!/bin/bash


MetadataFile=$1
genotypeInput=$2
classifiedBED=$3
OutBasename=${4%.gz}
maxfreqploy=$5

##############################
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#genotypeInput="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.bed.gz"
#classifiedBED="results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_ExtraCallableRegions.bed"
#OutBasename="results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.chr19_classified.genotype"
#maxfreqploy="0.97"
##############################
nsamples=$(tail -n +2 ${MetadataFile} | cut -f3 | sort -u | wc -l)

echo "----> Intersect variant and features, keep variant ID, feature type and feature ID"
bedtools intersect -a <(zcat ${genotypeInput} | cut -f1,2,3,4 ) -b <(cat ${classifiedBED}) -wo 2> ~/.null | cut -f4,8,9 > ${OutBasename}_CompleteCorrespVariantFeature.tmp

echo "----> Merge same-vairant entries by listing the features & features IDs"
awk '{if(NR==FNR){if(a[$1]){a[$1]=a[$1]","$2; b[$1]=b[$1]","$3}else{a[$1]=$2;b[$1]=$3}next} print $4"\t"a[$4]"\t"b[$4]}' ${OutBasename}_CompleteCorrespVariantFeature.tmp <(zcat ${genotypeInput}) | awk -F'\t' '{n=split($2,a,",");if(n>1){for(i in a){b[a[i]]=1}str="";for(i in b){str=i","str}print $1"\t"str"\t"$3;delete a; delete b}else{print $0}}' | sed 's/,\t/\t/g' > ${OutBasename}.tmp

echo "----> Frequency of each allele in each position"
zcat ${genotypeInput} | cut -f5- | sed 's/:/\t/g' | perl -ne '@a = split(" ", $_); %h=(); $h{$_}++ for @a; $str=""; for(sort { $h{$a} <=> $h{$b} } keys(%h)){$str=$str."$_:$h{$_} "} print $str."\n"' 2> ~/.null  > ${OutBasename}.freqalleles

echo "----> Determine which sites are polymorfic according to ${maxfreqploy}"
cat ${OutBasename}.freqalleles | sed 's/\S\+://g' | awk -F' ' -v nsamples="${nsamples}" -v maxfreqploy="${maxfreqploy}" '{if($NF/(nsamples*2) <= maxfreqploy){print 1}else{print 0}}' > ${OutBasename}.polymsites


ntmp=$(cat ${OutBasename}.tmp | wc -l)
nfreq=$(cat ${OutBasename}.freqalleles | wc -l)
npolym=$(cat ${OutBasename}.polymsites | wc -l)
nmajor=$(cat ${OutBasename}.major | wc -l)
if [[ ${ntmp} -ne ${nfreq} || ${ntmp} -ne ${npolym} ]]
then
	echo "ERROR in number of lines"; exit 1
fi
echo "Num;Feature;ID;NumAlleles;Polym;MajorFreq" | sed 's/;/\t/g'  > ${OutBasename}
paste ${OutBasename}.tmp <(cat ${OutBasename}.freqalleles | awk -F ' ' '{print NF}') ${OutBasename}.polymsites <(cat ${OutBasename}.freqalleles | awk -F ':' '{print $NF}') >> ${OutBasename}
gzip ${OutBasename}

rm ${OutBasename}*tmp 2> ~/null
rm ${OutBasename}.freqalleles 2> ~/null
rm ${OutBasename}.polymsites 2> ~/null



