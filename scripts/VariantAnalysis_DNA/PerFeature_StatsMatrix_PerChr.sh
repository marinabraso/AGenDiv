#!/bin/bash


MetadataFile=$1
ChrLengths=$2
chrVCF=$3
genotypeFileBasename=${4%.chr*}
CallableRegions=$5
SNPableRegions=$6
GenomicFeaturesBED=$7
GenomicFeaturesCallableBED=$8
SamplesOrderInVCF=$9
SynNonSynBED=${10}
SiteMatrix=${11}
OutBasename=${12%.tab.gz}
maxfreqploy=${13}
chr=${14}


###################################
### To debug
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#chrVCF="results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.chr14.vcf.gz"
#genotypeFileBasename="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered"
#CallableRegions="results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz"
#SamplesOrderInVCF="metadata/SamplesOrderInVCF.chr14.txt"
#GenomicFeaturesBED="results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures.bed"
#SynNonSynBED="results/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/SND_positions_callable.bed"
#SiteMatrix="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr14.tab.gz"
#OutBasename="results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.chr14"
#SNPableRegions="results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz"
#maxfreqploy=0.97
#chr="chr14"
###################################


mkdir -p $(dirname ${OutBasename})

populations=( $(tail -n +2 ${MetadataFile} | cut -f9 | sort -u) )
echo ${populations[@]}
sexes=( $(tail -n +2 ${MetadataFile} | cut -f10 | sort -u | grep -v 'Unknown') )
echo ${sexes[@]}
samples=( $(tail -n +2 ${MetadataFile} | cut -f3 | sort -u) )
echo ${samples[@]}
nsamples=$(tail -n +2 ${MetadataFile} | cut -f3 | sort -u | wc -l)
echo ${nsamples}
chrlen=$(grep -w "^>${chr}" ${ChrLengths} | cut -f2)
echo ${chrlen}

VCFsortedSamples=$(cat ${SamplesOrderInVCF})


echo "----> Calc the callable length from all features"
awk '{if(NR==FNR){if(a[$4"\t"$5]){a[$4"\t"$5]+=($3-$2)}else{a[$4"\t"$5]=($3-$2)}next}if(a[$4"\t"$5]){print $0"\t"a[$4"\t"$5]}else{print $0"\t0"}}' <(bedtools intersect -a <(cat ${GenomicFeaturesBED} | grep -w "^${chr}") -b ${CallableRegions}) <(cat ${GenomicFeaturesBED} | grep -w "^${chr}") >  ${OutBasename}.bed

echo "----> Calc the SNPable length from all features"
awk '{if(NR==FNR){if(a[$4"\t"$5]){a[$4"\t"$5]++$3-$2}else{a[$4"\t"$5]=$3-$2}next}if(a[$4"\t"$5]){print a[$4"\t"$5]}else{print "0"}}' <(bedtools intersect -a <(cat ${GenomicFeaturesBED} | grep -w "^${chr}") -b <(zcat ${SNPableRegions} | grep -w "^${chr}")) <(cat ${GenomicFeaturesBED} | grep -w "^${chr}") >  ${OutBasename}.SNPable

echo "----> Calculate the frequency of each allele in each position"
zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f5- | sed 's/:/\t/g' | awk '{if(NR==1){print $0;} print $0}' | perl -ne 'while($a=<>){@a = split(" ", $a); %h=(); $h{$_}++ for @a; $str=""; for(sort { $h{$a} <=> $h{$b} } keys(%h)){$str=$str."$_:$h{$_} "} print $str."\n"}' 2> ~/.null  > ${OutBasename}_${maxfreqploy}.freqalleles
cat ${OutBasename}_${maxfreqploy}.freqalleles | sed 's/\S\+://g' | awk -F' ' -v nsamples="${nsamples}" -v maxfreqploy="${maxfreqploy}" '{if($NF/(nsamples*2) <= maxfreqploy){print 1}else{print 0}}' > ${OutBasename}_${maxfreqploy}.persite.polymsites
paste <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) ${OutBasename}_${maxfreqploy}.persite.polymsites | sed 's/\t$//g' > ${OutBasename}_${maxfreqploy}.persite.polymsites.bed

echo "----> Bedtools intersect feature bed file with 0hom1het & polymsites"
bedtools intersect -a <(cut -f1,2,3,4,5 ${OutBasename}.bed) -b <(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz) -wao | awk 'NF{NF--};1' > ${OutBasename}.0hom1het.intersect
bedtools intersect -a <(cut -f1,2,3,4,5 ${OutBasename}.bed) -b ${OutBasename}_${maxfreqploy}.persite.polymsites.bed -wao | awk 'NF{NF--};1' > ${OutBasename}_${maxfreqploy}.polymsites.intersect

echo "----> Heterozygous sites & Polymorphic sites in feature regions"
cat ${OutBasename}.0hom1het.intersect | awk '{if(w == $1" "$2" "$3" "$4" "$5){nvar++}else{print w" "nvar; w=$1" "$2" "$3" "$4" "$5; if($NF != "."){nvar=1}else{nvar=0}}}END{print w" "nvar;}' | tail -n +2 | sed 's/\s/\t/g' > ${OutBasename}.varsites
cat ${OutBasename}.0hom1het.intersect | awk '{if(w == $1" "$2" "$3" "$4" "$5){for(i=10;i<=NF;i++){a[i]+=$i;}}else{str=w; for(i in a){str=str" "a[i]} print str; delete a; w=$1" "$2" "$3" "$4" "$5; for(i=10;i<=NF;i++){a[i]+=$i;}}}END{str=w; for(i in a){str=str" "a[i]} print str;}' | tail -n +2 | sed 's/\s/\t/g' > ${OutBasename}.hetsites
cat ${OutBasename}_${maxfreqploy}.polymsites.intersect | awk '{if(w == $1" "$2" "$3" "$4" "$5){npoly+=$NF}else{print w" "npoly; w=$1" "$2" "$3" "$4" "$5; npoly=$NF}}END{print w" "npoly;}' | sed 's/-1/0/g' | tail -n +2 | sed 's/\s/\t/g' > ${OutBasename}_${maxfreqploy}.polymsites

echo "----> Synonymous and non-synonymous sites per feature"
bedtools intersect -wao -a <(cat ${OutBasename}.bed | cut -f1,2,3,4,5) -b <(cat ${SynNonSynBED}) | awk '{col=NF-1;if(feat==$1"\t"$2"\t"$3"\t"$4"\t"$5){a[$col]+=1}else{if(a["S"]){s=a["S"]}else{s=0}if(a["N"]){n=a["N"]}else{n=0}if(a["D"]){d=a["D"]}else{d=0}print feat"\t"s"\t"n"\t"d;feat=$1"\t"$2"\t"$3"\t"$4"\t"$5;delete a}}END{if(a["S"]){s=a["S"]}else{s=0}if(a["N"]){n=a["N"]}else{n=0}if(a["D"]){d=a["D"]}else{d=0}print feat"\t"s"\t"n"\t"d;}' | tail -n +2 | awk '{if($4=="exon"){print $0}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t0\t0"}}' > ${OutBasename}.SND


echo "----> Synonymous and non-synonymous sites per feature"
SynNonSynDegColumn=$(zcat ${SiteMatrix} | head -1 | awk '{for(i=1;i<=NF;i++){if($i == "SynNonSynDeg"){print i}}}')
bedtools intersect -wao -a <(cat ${OutBasename}.bed | cut -f1,2,3,4,5) -b <(zcat ${SiteMatrix} | tail -n +2 | cut -f1,2,3,4,5,${SynNonSynDegColumn}) | awk '{col=NF-1;if(feat==$1"\t"$2"\t"$3"\t"$4"\t"$5){a[$col]+=1}else{if(a["S"]){s=a["S"]}else{s=0}if(a["N"]){n=a["N"]}else{n=0}if(a["D"]){d=a["D"]}else{d=0}print feat"\t"s"\t"n"\t"d;feat=$1"\t"$2"\t"$3"\t"$4"\t"$5;delete a}}END{if(a["S"]){s=a["S"]}else{s=0}if(a["N"]){n=a["N"]}else{n=0}if(a["D"]){d=a["D"]}else{d=0}print feat"\t"s"\t"n"\t"d;}' | tail -n +2 | awk '{if($4=="exon"){print $0}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t0\t0"}}' > ${OutBasename}.variantSND






echo "----> Checking that number of lines coincide"
nbed=$(cat ${OutBasename}.bed | wc -l)
nSNPable=$(cat ${OutBasename}.SNPable | wc -l)
nvarsites=$(cat ${OutBasename}.varsites | wc -l)
nhetsites=$(cat ${OutBasename}.hetsites | wc -l)
npolymsites=$(cat ${OutBasename}_${maxfreqploy}.polymsites | wc -l)
nSND=$(cat ${OutBasename}.SND | wc -l)
nvarSND=$(cat ${OutBasename}.variantSND | wc -l)
if [[ ${nbed} -ne ${nvarsites} || ${nbed} -ne ${nhetsites} || ${nbed} -ne ${npolymsites} || ${nbed} -ne ${nSND} || ${nbed} -ne ${nvarSND} ]]
then
	exit "ERROR in number of lines"
fi

echo "----> Pasting in a single file"
str="chr;st;end;type;name;callablesites;SNPablesites;varsites"
for s in ${VCFsortedSamples[@]}
do
	str=${str}";"hets${s}
done
str=${str}";polymsites;totalSynon;totalnonSynon;totalDeprec;variantSynon;variantnonSynon;variantDeprec"
echo ${str} 
echo ${str} | sed 's/;/\t/g' > ${OutBasename}.tab

paste <(cat ${OutBasename}.bed) \
${OutBasename}.SNPable \
<(cat ${OutBasename}.varsites | cut -f6) \
<(cat ${OutBasename}.hetsites | cut -f6- | sed 's/ /\t/g') \
<(cat ${OutBasename}_${maxfreqploy}.polymsites | cut -f6) \
<(cat ${OutBasename}.SND | cut -f6,7,8 ) \
<(cat ${OutBasename}.variantSND | cut -f6,7,8 ) >> ${OutBasename}.tab
gzip ${OutBasename}.tab

rm ${OutBasename}.bed 2> ~/null
rm ${OutBasename}.SNPable 2> ~/null
rm ${OutBasename}.varsites 2> ~/null
rm ${OutBasename}.hetsites 2> ~/null
rm ${OutBasename}*polymsite* 2> ~/null
rm ${OutBasename}*freqalleles 2> ~/null
rm ${OutBasename}*intersect 2> ~/null
rm ${OutBasename}*SND 2> ~/null










