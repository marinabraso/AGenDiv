#!/bin/bash


MetadataFile=$1
ChrLengths=$2
chrVCF=$3
genotypeFileBasename=${4%.chr*}
ExtraCallableRegions=$5
ValidWindows=$6
GenomicFeaturesBED=$7
SamplesOrderInVCF=$8
SiteMatrix=$9
genomeFA=${10}
OutBasename=${11%.tab.gz}
window=${12}
step=${13}
maxfreqploy=${14}
chr=${15}


###################################
### To debug
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
#chrVCF="results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.chr18.vcf.gz"
#genotypeFileBasename="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered"
#ExtraCallableRegions="results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ExtraCallableRegions_chr18.bed.gz"
#ValidWindows="results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ValidWindows_chr18.bed.gz"
#GenomicFeaturesBED="results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_ExtraCallableRegions.bed"
#SamplesOrderInVCF="metadata/SamplesOrderInVCF.chr189.txt"
#SiteMatrix="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr18.tab.gz"
#genomeFA="results/VariantAnalysis_DNA/mask_callable_regions_from_genome/Branchiostoma_lanceolatum.BraLan3_genome_CallableRegionsMasked.fa"
#OutBasename="results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr18"
#window=100000
#step=10000
#maxfreqploy="0.97"
#chr="chr18"
###################################


mkdir -p $(dirname ${OutBasename})

populations=( $(tail -n +2 ${MetadataFile} | cut -f9 | sort -u) )
sexes=( $(tail -n +2 ${MetadataFile} | cut -f10 | sort -u | grep -v 'Unknown') )
nsamples=$(tail -n +2 ${MetadataFile} | cut -f3 | sort -u | wc -l)
chrlen=$(grep -w "^>${chr}" ${ChrLengths} | cut -f2)
VCFsortedSamples=$(cat ${SamplesOrderInVCF})

echo "----> Calculate the frequency of each allele in each position"
zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f5- | sed 's/:/\t/g' | perl -ne 'while($a=<>){@a = split(" ", $a); %h=(); $h{$_}++ for @a; $str=""; for(sort { $h{$a} <=> $h{$b} } keys(%h)){$str=$str."$_:$h{$_} "} print $str."\n"}' 2> ~/.null  > ${OutBasename}_${maxfreqploy}.freqalleles
cat ${OutBasename}_${maxfreqploy}.freqalleles | sed 's/\S\+://g' | awk -F' ' -v nsamples="${nsamples}" -v maxfreqploy="${maxfreqploy}" '{if($NF/(nsamples*2) <= maxfreqploy){print 1}else{print 0}}' > ${OutBasename}_${maxfreqploy}.polymsites
paste <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) ${OutBasename}_${maxfreqploy}.polymsites > ${OutBasename}_${maxfreqploy}.polymsites.bed

echo "----> Build window bed file"
if [[ ! -s ${OutBasename}.bed ]]
then
	for ((w=0; w<=${chrlen}-step; w+=step))
	do
	end=$(( ${w}+${window} ))
	echo "${chr} ${w} ${end}" >> ${OutBasename}.bed
	done
	sed -i 's/ /\t/g' ${OutBasename}.bed
fi
nbed=$(cat ${OutBasename}.bed | wc -l)

echo "----> Bedtools intersect window bed file with 0hom1het & polymsites"
bedtools intersect -a ${OutBasename}.bed -b <(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed) -wao | awk 'NF{NF--};1' > ${OutBasename}.0hom1het.intersect
bedtools intersect -a ${OutBasename}.bed -b ${OutBasename}_${maxfreqploy}.polymsites.bed -wao | awk 'NF{NF--};1' > ${OutBasename}_${maxfreqploy}.polymsites.intersect

echo "----> Heterozygous sites & Polymorphic sites in sliding windows"
cat ${OutBasename}.0hom1het.intersect | awk '{if(w == $1" "$2){nvar++}else{print w" "nvar; w=$1" "$2; nvar=1}}END{print w" "nvar;}' | tail -n +2 | sed 's/ /\t/g' > ${OutBasename}.varsites
cat ${OutBasename}.0hom1het.intersect | awk '{if(w == $1" "$2){for(i=8;i<=NF;i++){a[i]+=$i;}}else{str=w; for(i in a){str=str" "a[i]} print str; delete a; w=$1" "$2; for(i=8;i<=NF;i++){a[i]+=$i;}}}END{str=w; for(i in a){str=str" "a[i]} print str;}' | tail -n +2 | sed 's/ /\t/g' > ${OutBasename}.hetsites
cat ${OutBasename}_${maxfreqploy}.polymsites.intersect | awk '{if(w == $1" "$2){npoly+=$NF}else{print w" "npoly; w=$1" "$2; npoly=$NF}}END{print w" "npoly;}' | tail -n +2 | sed 's/ /\t/g' > ${OutBasename}_${maxfreqploy}.polymsites
nvarsites=$(cat ${OutBasename}.varsites | wc -l)
nhetsites=$(cat ${OutBasename}.hetsites | wc -l)
npolymsites=$(cat ${OutBasename}_${maxfreqploy}.polymsites | wc -l)
if [[ ${nbed} -ne ${nvarsites} || ${nbed} -ne ${nhetsites} || ${nbed} -ne ${npolymsites} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Pi (average pairwise differences) in sliding windows" # only considers PASS variants
vcftools --gzvcf ${chrVCF} --window-pi ${window} --window-pi-step ${step} --out ${OutBasename}_tmp
awk '{if(NR==FNR){a[$1"\t"$2]=$3;next}if(a[$1"\t"$2]){print a[$1"\t"$2]}else{print 0}}' <(tail -n +2 ${OutBasename}_tmp.windowed.pi | cut -f2,3,5 | awk '{st=$1-1; print st"\t"$2"\t"$3}') <(cat ${OutBasename}.bed | cut -f2,3) > ${OutBasename}.windowed.pi
npi=$(cat ${OutBasename}.windowed.pi | wc -l)
if [[ ${nbed} -ne ${npi} ]]
then
	echo "ERROR in number of lines"; exit 1
fi



echo "----> Fst"
if [[ ! -s $(dirname $OutBasename)/${populations[1]}.txt || ! -s $(dirname $OutBasename)/${sexes[1]}.txt || ! -s $(dirname $OutBasename)/${populations[1]}_${sexes[1]}.txt ]]
then
	for pop in ${populations[@]}
	do
		tail -n +2 ${MetadataFile} | grep -w ${pop} | cut -f3 > $(dirname $OutBasename)/${pop}.txt
		for sex in ${sexes[@]}
		do
			tail -n +2 ${MetadataFile} | grep -w ${pop} | grep -w ${sex} | cut -f3 > $(dirname $OutBasename)/${pop}_${sex}.txt
		done
	done
	for sex in ${sexes[@]}
	do
		tail -n +2 ${MetadataFile} | grep -w ${sex} | cut -f3 > $(dirname $OutBasename)/${sex}.txt
	done
fi

strp=""
for pop in ${populations[@]}
do
	strp=${strp}" --weir-fst-pop "$(dirname $OutBasename)"/"${pop}".txt"
done
strs=""
for sex in ${sexes[@]}
do
	strs=${strs}" --weir-fst-pop "$(dirname $OutBasename)"/"${sex}".txt"
done

echo "----> By population" # only considers PASS variants
vcftools --gzvcf ${chrVCF} ${strp} --fst-window-size ${window} --fst-window-step ${step} --out ${OutBasename}_pop_tmp

echo "----> By sex" # only considers PASS variants
vcftools --gzvcf ${chrVCF} ${strs} --fst-window-size ${window} --fst-window-step ${step} --out ${OutBasename}_sex_tmp

echo "----> By sex for each population"
for pop in ${populations[@]}
do
	echo "----> ${pop}"
	strps=""
	for sex in ${sexes[@]}
	do
		strps=${strps}" --weir-fst-pop "$(dirname $OutBasename)"/"${pop}_${sex}".txt"
	done
	echo "      ----> Fst in sliding windows" # only considers PASS variants
	vcftools --gzvcf ${chrVCF} ${strps} --fst-window-size ${window} --fst-window-step ${step} --out ${OutBasename}_${pop}_sex_tmp
done

awk '{if(NR==FNR){a[$1"\t"$2]=$3"\t"$4;next}if(a[$1"\t"$2]){print a[$1"\t"$2]}else{print NA"\t"NA}}' <(tail -n +2 ${OutBasename}_pop_tmp.windowed.weir.fst | cut -f2,3,5,6 | awk '{st=$1-1; print st"\t"$2"\t"$3"\t"$4}') <(cat ${OutBasename}.bed | cut -f2,3) > ${OutBasename}_pop.windowed.weir.fst
awk '{if(NR==FNR){a[$1"\t"$2]=$3"\t"$4;next}if(a[$1"\t"$2]){print a[$1"\t"$2]}else{print NA"\t"NA}}' <(tail -n +2 ${OutBasename}_sex_tmp.windowed.weir.fst | cut -f2,3,5,6 | awk '{st=$1-1; print st"\t"$2"\t"$3"\t"$4}') <(cat ${OutBasename}.bed | cut -f2,3) > ${OutBasename}_sex.windowed.weir.fst
npfst=$(cat ${OutBasename}_pop.windowed.weir.fst | wc -l)
nsfst=$(cat ${OutBasename}_sex.windowed.weir.fst | wc -l)
if [[ ${nbed} -ne ${npfst} || ${nbed} -ne ${nsfst} ]]
then
	exit "ERROR in number of lines"
fi
for pop in ${populations[@]}
do
	awk '{if(NR==FNR){a[$1"\t"$2]=$3"\t"$4;next}if(a[$1"\t"$2]){print a[$1"\t"$2]}else{print NA"\t"NA}}' <(tail -n +2 ${OutBasename}_${pop}_sex_tmp.windowed.weir.fst | cut -f2,3,5,6 | awk '{st=$1-1; print st"\t"$2"\t"$3"\t"$4}') <(cat ${OutBasename}.bed | cut -f2,3) > ${OutBasename}_${pop}_sex.windowed.weir.fst
	nfst=$(cat ${OutBasename}_${pop}_sex.windowed.weir.fst | wc -l)
	if [[ ${nbed} -ne ${nfst} ]]
	then
		exit "ERROR in number of lines"
	fi
done



echo "----> Genomic features"
# Genomic features that intersect with each window
bedtools intersect -a ${OutBasename}.bed -b <(grep -w ${chr} ${GenomicFeaturesBED} | grep -v "Interg") -wao | cut -f1,2,3,8 | sort | uniq | awk '{if(w == $1"\t"$2"\t"$3){list=list","$4}else{print w"\t"list; list=$4; w=$1"\t"$2"\t"$3}}END{print w"\t"list;}' | tail -n +2 > ${OutBasename}.genomicfeatures.bed
# Exon content
bedtools intersect -a ${OutBasename}.bed -b <(bedtools merge -i <(grep -w ${chr} ${GenomicFeaturesBED} | grep 'exon' | cut -f1,2,3 | sort -k2,2V)) -wao | awk '{print $1"\t"$2"\t"$3"\t"$NF}' | awk '{if(w == $1"\t"$2"\t"$3){val += $4}else{print w"\t"val; val = $4; w = $1"\t"$2"\t"$3}}END{print w"\t"val}' | tail -n +2 > ${OutBasename}.exoncontent.bed

nec=$(cat ${OutBasename}.exoncontent.bed | wc -l)
ngf=$(cat ${OutBasename}.genomicfeatures.bed | wc -l)
if [[ ${nbed} -ne ${ngf} || ${nbed} -ne ${nec} ]]
then
	exit "ERROR in number of lines"
fi

echo "----> Callable length per window"
bedtools intersect -a ${OutBasename}.bed -b <(zcat ${ExtraCallableRegions}) -wao | awk '{if(w==$1"\t"$2"\t"$3){cl+=$NF}else{print w"\t"cl; w=$1"\t"$2"\t"$3;cl=$NF}}END{print w"\t"cl;}' | tail -n +2 > ${OutBasename}.callablelength.bed
ncl=$(cat ${OutBasename}.callablelength.bed | wc -l)
if [[ ${nbed} -ne ${ncl} ]]
then
	exit "ERROR in number of lines"
fi

echo "----> Callable length per window"
awk '{if(NR==FNR){a[$1"\t"$2"\t"$3]=1;next}if(a[$1"\t"$2"\t"$3]){print 1}else{print 0}}' <(zcat ${ValidWindows}) ${OutBasename}.bed > ${OutBasename}.validwindows.bed
nvw=$(cat ${OutBasename}.validwindows.bed | wc -l)
if [[ ${nbed} -ne ${nvw} ]]
then
	exit "ERROR in number of lines"
fi

echo "----> PC average values"
PCColNums=$(zcat ${SiteMatrix} | head -1 | awk '{str=""; for(i=1;i<=NF;i++){if($i ~ /PC/){str=str","i}}print str}')
PCColNames=$(zcat ${SiteMatrix} | head -1 | awk '{str=""; for(i=1;i<=NF;i++){if($i ~ /PC/){str=str";"$i}}print str}')

bedtools intersect -a ${OutBasename}.bed -b <(zcat ${SiteMatrix} | cut -f1,2,3${PCColNums} | tail -n +2) -wao | awk '{if(w == $1"\t"$2"\t"$3){for(i=7;i<NF;i++){if($i != "NA"){v[i]+=$i;n[i]++}}}else{str=w; for(i=7;i<NF;i++){if(n[i]){val=v[i]/n[i]}else{val="NA"};str=str"\t"val}print str;delete(v);delete(n);w=$1"\t"$2"\t"$3; for(i=7;i<NF;i++){if($i != "NA"){v[i]+=$i;n[i]++}}}}END{str=w; for(i=7;i<NF;i++){if(n[i]){val=v[i]/n[i]}else{val="NA"};str=str"\t"val}print str;}' | tail -n +2 > ${OutBasename}.PCaverages.bed
nPC=$(cat ${OutBasename}.PCaverages.bed | wc -l)
if [[ ${nbed} -ne ${nPC} ]]
then
	exit "ERROR in number of lines"
fi


echo "----> Num alternative alleles per sample"
GenoFreqNalleleColumns=$(zcat ${SiteMatrix} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/ || $i=="freq" || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
bedtools intersect -a ${OutBasename}.bed -b <(zcat ${SiteMatrix} | tail -n +2 | cut -f1,2,3,${GenoFreqNalleleColumns} | awk '{if($NF==2){fc=NF-1; split($fc,a,";|:"); str=$1"\t"$2"\t"$3; for(i=4;i<fc;i++){split($i,b,":");alt=0; for(j in b){if(b[j]==a[1]){alt++}} str=str"\t"alt} print str}}') -wao | awk '{if(w == $1"\t"$2"\t"$3){for(i=7;i<NF;i++){if($i != "NA"){v[i]+=$i}}}else{str=w; for(i=7;i<NF;i++){str=str"\t"v[i]}print str;delete(v);w=$1"\t"$2"\t"$3; for(i=7;i<NF;i++){v[i]+=$i}}}END{str=w; for(i=7;i<NF;i++){str=str"\t"v[i]}print str;}' | tail -n +2 > ${OutBasename}.AlternAlleles.bed

nalt=$(cat ${OutBasename}.AlternAlleles.bed | wc -l)
if [[ ${nbed} -ne ${nalt} ]]
then
	exit "ERROR in number of lines"
fi


echo "----> GC content"
awk '{if(NR==FNR){a[$1"\t"$2"\t"$3]=$4"\t"$5;next}if(a[$1"\t"$2"\t"$3]){print $1"\t"$2"\t"$3"\t"a[$1"\t"$2"\t"$3]}else{print $1"\t"$2"\t"$3"\tNA\tNA"}}' <(bedtools nuc -fi ${genomeFA} -bed ${OutBasename}.bed | tail -n +2 | cut -f1,2,3,4,5) ${OutBasename}.bed > ${OutBasename}.GCcontent.bed
nGC=$(cat ${OutBasename}.GCcontent.bed | wc -l)
if [[ ${nbed} -ne ${nGC} ]]
then
	echo "ERROR in number of lines"; exit 1
fi


echo "----> Num of variants per bins of frequency"
FreqNalleleColumns=$(zcat ${SiteMatrix} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i=="freq" || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
bedtools intersect -a ${OutBasename}.bed -b <(zcat ${SiteMatrix} | tail -n +2 | cut -f1,2,3,${FreqNalleleColumns} | awk '{if($NF == 2){fc=NF-1; split($fc,a,";|:"); freq = a[2]/(a[2]+a[4]); print $1"\t"$2"\t"$3"\t"freq}}') -wao | awk '{printf "%s\t%s\t%s\t%.0f\n", $1, $2, $3, $7*10}' | awk '{if(w == $1"\t"$2"\t"$3){a[$4]++}else{str=w; for(i=0;i<=5;i++){if(a[i]){str=str"\t"a[i];}else{str=str"\t"0}} print str; w = $1"\t"$2"\t"$3; delete a; a[$4]++}}END{str=w; for(i=0;i<=5;i++){if(a[i]){str=str"\t"a[i];}else{str=str"\t"0}} print str;}' | tail -n +2 > ${OutBasename}.BinsOfMinFreq.bed
nbinfreq=$(cat ${OutBasename}.BinsOfMinFreq.bed | wc -l)
controlncol=$(cat ${OutBasename}.BinsOfMinFreq.bed | awk '{print NF}' | sort | uniq | wc -l)
if [[ ${nbed} -ne ${nbinfreq} || ${controlncol} -ne 1 ]]
then
	exit "ERROR in number of lines or number of columns"
fi





################################################
echo "----> Pasting in a single file"
str="chr;st;end;varsites;valid;callablelength"
for s in ${VCFsortedSamples[@]}
do
	str=${str}";"hets${s}
done
str=${str}";polymsites;PI;weighted_pop_FST;mean_pop_FST;weighted_sex_FST;mean_sex_FST"
for pop in ${populations[@]}
do
	str=${str}";weighted_"${pop}"_sex_FST;mean_"${pop}"_sex_FST"
done
str=${str}";features;exoncontent"${PCColNames}
for s in ${VCFsortedSamples[@]}
do
	str=${str}";"altAl${s}
done
str=${str}";ATcontent;GCcontent;NumVarFreq0-.05;NumVarFreq.05-.15;NumVarFreq.15-.25;NumVarFreq.25-.35;NumVarFreq.35-.45;NumVarFreq.45-.5"
echo ${str} 
echo ${str} | sed 's/;/\t/g' > ${OutBasename}.tab

paste ${OutBasename}.bed \
<(cat ${OutBasename}.varsites | cut -f3) \
<(cat ${OutBasename}.validwindows.bed) \
<(cat ${OutBasename}.callablelength.bed | cut -f4) \
<(cat ${OutBasename}.hetsites | cut -f3-) \
<(cat ${OutBasename}_${maxfreqploy}.polymsites | cut -f3) \
${OutBasename}.windowed.pi \
${OutBasename}_pop.windowed.weir.fst \
${OutBasename}_sex.windowed.weir.fst > ${OutBasename}.tmp
for pop in ${populations[@]}
do
	paste ${OutBasename}.tmp ${OutBasename}_${pop}_sex.windowed.weir.fst > ${OutBasename}.tmp2
	mv ${OutBasename}.tmp2 ${OutBasename}.tmp
done
paste ${OutBasename}.tmp \
<(cat ${OutBasename}.genomicfeatures.bed | cut -f4) \
<(cat ${OutBasename}.exoncontent.bed | cut -f4) \
<(cat ${OutBasename}.PCaverages.bed | cut -f4-) \
<(cat ${OutBasename}.AlternAlleles.bed | cut -f4-) \
<(cat ${OutBasename}.GCcontent.bed | cut -f4,5) \
<(cat ${OutBasename}.BinsOfMinFreq.bed | cut -f4-)  >> ${OutBasename}.tab

gzip ${OutBasename}.tab

rm ${OutBasename}.bed 2> ~/null
rm ${OutBasename}.varsites 2> ~/null
rm ${OutBasename}.hetsites 2> ~/null
rm ${OutBasename}_${maxfreqploy}.polymsite* 2> ~/null
rm ${OutBasename}_${maxfreqploy}.freqalleles 2> ~/null
rm ${OutBasename}.windowed.pi 2> ~/null
rm ${OutBasename}*.weir.fst 2> ~/null
rm ${OutBasename}*intersect 2> ~/null
rm ${OutBasename}.tmp 2> ~/null
#rm ${OutBasename}*genomicfeatures.bed 2> ~/null
rm ${OutBasename}.callablelength.bed 2> ~/null
rm ${OutBasename}*tmp* 2> ~/null
rm ${OutBasename}.PCaverages.bed 2> ~/null
rm ${OutBasename}.AlternAlleles.bed 2> ~/null
rm ${OutBasename}.BinsOfMinFreq.bed 2> ~/null
rm ${OutBasename}.GCcontent.bed 2> ~/null
#rm ${OutBasename}.exoncontent.bed 2> ~/null
rm ${OutBasename}.validwindows.bed 2> ~/null




