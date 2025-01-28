#!/bin/bash


MetadataFile=$1
chrVCF=$2
genotypeFileBasename=${3%.chr*}
Coverage=${4%.MinCov*}
SamplesOrderInVCF=$5
GenomicFeatureDivision=$6
SynNonSynBED=$7
PCAbasename=${8%PCA*}
OutBasename=${9%.tab.gz}
maxfreqploy=${10}
chr=${11}

###################################
### To debug
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#chrVCF="results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.chr19.vcf.gz"
#genotypeFileBasename="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered"
#Coverage="results/VariantCalling_DNA/7_join_coverage_per_site/Depth"
#SamplesOrderInVCF="metadata/SamplesOrderInVCF.chr19.txt"
#GenomicFeatureDivision="results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.chr19_classified.genotype.gz"
#OutBasename="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr19"
#SynNonSynBED="results/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/SND_positions_callable.bed"
#PCAbasename="results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCvalues_PerVariant_"
#window=10000
#step=1000
#maxfreqploy="0.97"
#chr="chr19"
###################################

mkdir -p $(dirname ${OutBasename})

populations=( $(tail -n +2 ${MetadataFile} | cut -f9 | sort -u) )
echo ${populations[@]}
sexes=( $(tail -n +2 ${MetadataFile} | cut -f10 | sort -u | grep -v 'Unknown') )
echo ${sexes[@]}
nsamples=$(tail -n +2 ${MetadataFile} | cut -f3 | sort -u | wc -l)
echo ${nsamples}
VCFsortedSamples=$(cat ${SamplesOrderInVCF})



echo "----> Frequency of each allele in each position"
zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f5- | sed 's/:/\t/g' | perl -ne '@a = split(" ", $_); %h=(); $h{$_}++ for @a; $str=""; for(sort { $h{$a} <=> $h{$b} } keys(%h)){$str=$str."$_:$h{$_} "} print $str."\n"' 2> ~/.null  > ${OutBasename}.freqalleles
nbed=$(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | wc -l)
nfreqalleles=$(cat ${OutBasename}.freqalleles | wc -l)
if [[ ${nbed} -ne ${nfreqalleles} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Number of alleles"
cat ${OutBasename}.freqalleles | awk -F '\t' '{n=split($1,a," "); print n}' > ${OutBasename}.numalleles
nnumalleles=$(cat ${OutBasename}.numalleles | wc -l)
if [[ ${nbed} -ne ${nnumalleles} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Min & Max length among alleles (to distinguish SNPs from INDELs)"
zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | awk '{for(i=5;i<=NF;i++){split($i,a,":");for(j in a){al[a[j]]=1}} min=999999999; max=0; for(i in al){n=split(i,a,"");if(n>max){max=n}if(n<min){min=n}if(i=="*"){min=0}}; if(min==1 && max==1){type="SNP"}else{type="INDEL"}print $1"\t"$2"\t"$3"\t"$4"\t"type; delete al}' > ${OutBasename}.typeVar
ntypeVar=$(cat ${OutBasename}.typeVar | wc -l)
if [[ ${nbed} -ne ${ntypeVar} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Determine which sites are polymorfic according to ${maxfreqploy}"
cat ${OutBasename}.freqalleles | sed 's/\S\+://g' | awk -F' ' -v nsamples="${nsamples}" -v maxfreqploy="${maxfreqploy}" '{if($NF/(nsamples*2) <= maxfreqploy){print 1}else{print 0}}' > ${OutBasename}.polymsites
npolymsites=$(cat ${OutBasename}.polymsites | wc -l)
if [[ ${nbed} -ne ${npolymsites} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Number of heterozygots per site"
zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz | cut -f5- | awk '{for(i=1;i<=NF;i++){sum+=$i}print sum; sum=0}' > ${OutBasename}.het
nhet=$(cat ${OutBasename}.het | wc -l)
if [[ ${nbed} -ne ${nhet} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Pi (average pairwise differences) per site" # only considers PASS variants
vcftools --gzvcf ${chrVCF} --remove-filtered-all --site-pi --out ${OutBasename}
npi=$(tail -n +2 ${OutBasename}.sites.pi | wc -l)
if [[ ${nbed} -ne ${npi} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Prepare files for Fst per site"
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

echo "----> Prepare strings for Fst per site"
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

echo "----> Fst by population" # only considers PASS variants
vcftools --gzvcf ${chrVCF} --remove-filtered-all ${strp} --out ${OutBasename}_pop
echo "----> Fst by sex" # only considers PASS variants
vcftools --gzvcf ${chrVCF} --remove-filtered-all ${strs} --out ${OutBasename}_sex

echo "----> Fst by sex for each population"
for pop in ${populations[@]}
do
	echo "    ----> ${pop}"
	strps=""
	for sex in ${sexes[@]}
	do
		strps=${strps}" --weir-fst-pop "$(dirname $OutBasename)"/"${pop}_${sex}".txt"
	done
	vcftools --gzvcf ${chrVCF} --remove-filtered-all ${strps} --out ${OutBasename}_${pop}_sex
done
npfst=$(tail -n +2 ${OutBasename}_pop.weir.fst | wc -l)
nsfst=$(tail -n +2 ${OutBasename}_sex.weir.fst | wc -l)
if [[ ${npfst} -ne ${npfst} || ${nsfst} -ne ${nsfst} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Hardy-Weinberg Equilibrium test"  # only considers PASS variants
# (as defined by Wigginton, Cutler and Abecasis (2005))"
vcftools --gzvcf ${chrVCF} --remove-filtered-all --hardy --out ${OutBasename}
awk '{if(FNR == NR){a[$1"."$2]=$0;next}if(a[$1"."$2]){print a[$1"."$2]}else{print $1"\t"$2"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"}}' <(tail -n +2 ${OutBasename}.hwe | sed 's/\//\t/g') <(tail -n +2 ${OutBasename}.sites.pi) > ${OutBasename}_complete.hwe
nhwe=$(cat ${OutBasename}_complete.hwe | wc -l)
if [[ ${nbed} -ne ${nhwe} ]]
then
	echo "ERROR in number of lines"; exit 1
fi


echo "----> Coverage per site"
bedtools intersect -a <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) -b <(zcat ${Coverage}.MinCov.tbl.gz | grep -w ${chr} | awk '{print $1"\t"$2"\t"$2"\t"$3}') -wao  | awk '{if($6>=$2 && $6<=$3){print $0}}' | cut -f8 > ${OutBasename}.mincov
bedtools intersect -a <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) -b <(zcat ${Coverage}.allsamples.tbl.gz | grep -w ${chr} | awk '{str=$1"\t"$2"\t"$2; for(i=3;i<=NF;i++){str=str"\t"$i} print str}') -wao | awk '{if($6>=$2 && $6<=$3){print $0}}' | cut -f8- | awk '{str=$1; for(i=2;i<NF;i++){str=str"\t"$i}print str}' > ${OutBasename}.allcov
nmincov=$(cat ${OutBasename}.mincov | wc -l)
nallcov=$(cat ${OutBasename}.allcov | wc -l)
if [[ ${nbed} -ne ${nmincov} || ${nsfst} -ne ${nallcov}  ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Synonymous and non-Synonymous positions"  # only SNPs are taken into account (INDELs == NA)
awk '{if(NR==FNR){a[$2]=$4;next}if($2==$3 && $NF=="SNP" && a[$2]){print a[$2]}else{print "NA"}}' <(cat ${SynNonSynBED} | grep -w ${chr}) <(cat ${OutBasename}.typeVar) > ${OutBasename}.SnS
nSND=$(cat ${OutBasename}.SnS | wc -l)
if [[ ${nbed} -ne ${nSND} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> PCA axis values"  
Columns=$(zcat ${PCAbasename}PCA_all_all_all_all.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | head -1 | awk '{str= ""; for(i=1;i<=NF;i++){if($i ~ /Axis/ || $i == "num"){str=str","i}}print str}' | sed 's/^,//g')
zcat ${PCAbasename}PCA_all_all_all_all.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PCall
zcat ${PCAbasename}PCA_all_polymorphic_all_all.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PCpolym
zcat ${PCAbasename}PCA_all_all_biallelic_all.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PCbial
zcat ${PCAbasename}PCA_exon_all_all_all.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PCexon
zcat ${PCAbasename}PCA_all_all_all_0.5to0.6.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PC0.5to0.6
zcat ${PCAbasename}PCA_all_all_all_0.6to0.7.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PC0.6to0.7
zcat ${PCAbasename}PCA_all_all_all_0.7to0.8.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PC0.7to0.8
zcat ${PCAbasename}PCA_all_all_all_0.8to0.9.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PC0.8to0.9
zcat ${PCAbasename}PCA_all_all_all_0.9to1.txt.gz | sed 's/"//g' | sed 's/^num/Rrow\tnum/' | grep -w ${chr} | cut -f${Columns} | awk '{if(NR==1){num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i};next}if(v==$1){for(i=2;i<=NF;i++){a[i]+=$i;}num+=1}else{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str; delete a; num=1;v=$1;for(i=2;i<=NF;i++){a[i]+=$i}}}END{str=v"\t"num; for(i=2;i<=NF;i++){mean=a[i]/num;str=str"\t"mean;}print str;}'  > ${OutBasename}.PCtmp
awk '{if(NR==FNR){a[$1]=$0;next}if(a[$4]){print a[$4]}else{print $4"\tNA\tNA\tNA\tNA\tNA\tNA"}}' ${OutBasename}.PCtmp <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz | cut -f1,2,3,4) >  ${OutBasename}.PC0.9to1
nPC=$(cat ${OutBasename}.PCall | wc -l)
if [[ ${nbed} -ne ${nPC} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

echo "----> Checking number of lines of GenomicFeatureDivision"
ngenfeat=$(zcat ${GenomicFeatureDivision} | tail -n +2 | wc -l)
if [[ ${nbed} -ne ${ngenfeat} ]]
then
	echo "ERROR in number of lines"; exit 1
fi

###############################################
echo "----> Pasting in a single file"
str="chr;st;end;numvar"
for s in ${VCFsortedSamples[@]}
do
	str=${str}";"geno${s}
done

str=${str}";freq;numalleles;typeVar;polymorhism;oHet;PI;pop_FST;sex_FST"
for pop in ${populations[@]}
do
	str=${str}";"${pop}"_sex_FST"
done
str=${str}";GenomicFeatureType;GenomicFeatureID;SynNonSynDeg;oHOM1_hwe;oHET_hwe;oHOM2_hwe;eHOM1_hwe;eHET_hwe;eHOM2_hwe;ChiSq_hwe;P_hwe;P_HET_DEFICIT_hwe;P_HET_EXCESS_hwe;MinCoverage"
for s in ${VCFsortedSamples[@]}
do
	str=${str}";cov"${s}
done
for s in all polym bial exon 0.5to0.6 0.6to0.7 0.7to0.8 0.8to0.9 0.9to1
do
	str=${str}";PC1"${s}";PC2"${s}";PC3"${s}";PC4"${s}";PC5"${s}
done
echo ${str} 
echo ${str} | sed 's/;/\t/g' > ${OutBasename}.tab

paste <(zcat ${genotypeFileBasename}.${chr}.genotype.bed.gz) \
<(cat ${OutBasename}.freqalleles | sed 's/ $//g' | sed 's/ /;/g') \
${OutBasename}.numalleles \
<(cat ${OutBasename}.typeVar | cut -f5) \
${OutBasename}.polymsites \
${OutBasename}.het \
<(tail -n +2 ${OutBasename}.sites.pi | cut -f3) \
<(tail -n +2 ${OutBasename}_pop.weir.fst | cut -f3) \
<(tail -n +2 ${OutBasename}_sex.weir.fst | cut -f3) > ${OutBasename}.tmp
for pop in ${populations[@]}
do
	paste ${OutBasename}.tmp <(tail -n +2 ${OutBasename}_${pop}_sex.weir.fst | cut -f3) > ${OutBasename}.tmp2
	mv ${OutBasename}.tmp2 ${OutBasename}.tmp
done
paste ${OutBasename}.tmp <(zcat ${GenomicFeatureDivision} | tail -n +2 | cut -f2,3) ${OutBasename}.SnS <(cat ${OutBasename}_complete.hwe | cut -f3-) ${OutBasename}.mincov ${OutBasename}.allcov >> ${OutBasename}.tmp2
mv ${OutBasename}.tmp2 ${OutBasename}.tmp
for s in all polym bial exon 0.5to0.6 0.6to0.7 0.7to0.8 0.8to0.9 0.9to1
do
	paste ${OutBasename}.tmp <(cat ${OutBasename}.PC${s} | cut -f3-) > ${OutBasename}.tmp2
	mv ${OutBasename}.tmp2 ${OutBasename}.tmp
done
cat ${OutBasename}.tmp >> ${OutBasename}.tab
gzip ${OutBasename}.tab

rm ${OutBasename}.tmp   2> ~/null
rm ${OutBasename}*hwe  2> ~/null
rm ${OutBasename}*.weir.fst  2> ~/null
rm ${OutBasename}.sites.pi  2> ~/null
rm ${OutBasename}.freqalleles  2> ~/null
rm ${OutBasename}.numalleles  2> ~/null
rm ${OutBasename}.typeVar  2> ~/null
rm ${OutBasename}.polymsites 2> ~/null
rm ${OutBasename}.het 2> ~/null
rm ${OutBasename}.mincov 2> ~/null
rm ${OutBasename}.allcov 2> ~/null
rm ${OutBasename}.SnS 2> ~/null
rm ${OutBasename}.PCtmp 2> ~/null
rm ${OutBasename}.PC* 2> ~/null










