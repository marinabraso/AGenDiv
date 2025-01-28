#!/bin/bash


MetadataFile=$1
VCF=$2
chrshomhetFiles=$3
CallableRegions=$4
SamplesOrderInVCF=$5
strchrsSiteMatrices=$6
PCAbasename=${7%PCA_*}
UMAPbasename=${8%UMAP*}
OutBasename=${9%.tab.gz}
chrs=${10}


############
# debug / develop
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#VCF="results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardExtraCallableFiltered.vcf.gz"
#chrshomhetFiles="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr18.0hom1het.bed.gz results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.0hom1het.bed.gz"
#CallableRegions="results/VariantCalling_DNA/7_join_extra_callable_regions/ExtraCallableRegions.bed.gz"
#SamplesOrderInVCF="metadata/SamplesOrderInVCF.chr19.txt"
#strchrsSiteMatrices="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr18.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr19.tab.gz"
#PCAbasename="results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCvalues_PerSample_"
#UMAPbasename="results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCvalues_PerSample_"
#OutBasename="results/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSampleStats"
#chrs="chr18 chr19"
############

mkdir -p $(dirname ${OutBasename})

VCFsortedSamples=( $(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') )
chrs=( $(echo ${chrs} | sed 's/ /\n/g') )
chrsSiteMatrices=( $(echo ${strchrsSiteMatrices} | sed 's/ /\n/g') )
rm ${OutBasename}.tmp ${OutBasename}.het 2> ~/.null

echo "----> vcftools heterozygosity"  # only considers PASS variants
vcftools --gzvcf ${VCF} --remove-filtered-all --het --out ${OutBasename}

echo "----> Expected heterozygous & homozygous sites per sample" 
# !!!! Only considering biallelic variants
# Column with allele frequency data
freqcol=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{for(i=0;i<=NF;i++){if($i=="freq"){print i}}}')
# total number of biallelic sites --> same for all samples
nBiallSites=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${freqcol} | awk '{n=split($1,a,";"); if(n==2){num++}}END{print num}')
# Expected number of heterozygous sites per sample:  --> same for all samples
# sum(2pq)=sum(2*numPalleles/(numsamples*2)*numQalleles/(numsamples*2)) for all biallelic sites
eHet=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${freqcol} | awk -v nsamp=${#VCFsortedSamples[@]} '{n=split($1,a,";"); if(n==2){ehet=2; for(i in a){split(a[i],b,":");ehet*=b[2]/nsamp/2;} printf "%.6f\n", ehet;}}' | awk '{sum+=$1;}END{print sprintf("%.9f", sum)}')
# Expected number of homozygous sites per sample: --> same for all samples
# sum(p^2+q^2)=sum(numPalleles*numPalleles/(nsamp*nsamp*4)+numQalleles*numQalleles/(nsamp*nsamp*4)) for all biallelic sites
eHom=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${freqcol} | awk -v nsamp=${#VCFsortedSamples[@]} '{n=split($1,a,";"); if(n==2){ehom=0; for(i in a){split(a[i],b,":");ehom+=b[2]*b[2]/(nsamp*nsamp*4);} printf "%.6f\n", ehom;}}' | awk '{sum+=$1;}END{print sprintf("%.9f", sum)}')

echo "----> Observed heterozygous & homozygous sites per sample"
# number of observed heterozygous sites per each sample
hetsites=( $(zcat ${chrshomhetFiles} | awk '{for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){str=str" "a[i]} print str}' | sed 's/^ //g' ) )
# number of observed homozygous sites per each sample
homsites=( $(zcat ${chrshomhetFiles} | awk '{num++;for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){a[i]=num-a[i]; str=str" "a[i]} print str}' | sed 's/^ //g' ) )
# number of variant sites (at least two alleles present) --> same for all samples
varsites=$(zcat ${chrshomhetFiles} | awk '{for(i=5;i<=NF;i++){s+=$i}if(s>0){print s}s=0}' | wc -l )
# number of callable sites (check r7_callable_regions for details) --> same for all samples
callsites=$(zcat ${CallableRegions} | awk '{len+=$3-$2}END{print len}' )
if [ ${callsites} -eq 0 ]; then
	echo "ERROR callsites = 0"; exit 1;
fi
# number of observed biallelic heterozygous sites per each sample
biallhetsites=( $(awk '{if(NR==FNR){n=split($1,a,";"); if(n==2){blines[FNR]=1;}next} if(blines[FNR]){print $0}}' <(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${freqcol}) <(zcat ${chrshomhetFiles}) | awk '{for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){str=str" "a[i]} print str}' | sed 's/^ //g' ) )

echo "----> Number of minor alleles per sample"
GenoFreqColumns=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/ || $i=="freq"){str=str","i}}print str}' | sed 's/^,//g')
minalleles=( $(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of singletons per sample"
singletons=( $(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-1)){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of doubletons per sample" # only biallelic
doubletons=( $(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-2) && length(a)==4){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of tripletons per sample" # only biallelic
tripletons=( $(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-3) && length(a)==4){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Average pairwise differences"
GenoColumns=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/){str=str","i}}print str}' | sed 's/^,//g')
GenoNalleleColumns=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/ || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
pi=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoColumns} | sed 's/:/\t/g' | awk -v callsites=${callsites} '{if(NR==1){NS=NF;c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;a[c]=0}}}comp=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){comp++;if($i!=$j){a[comp]++}}}}END{c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;sum+=a[c]/callsites}}print sum/c}')
piBi=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${GenoNalleleColumns} | sed 's/:/\t/g' | awk -v callsites=${callsites} '{if(NR==1){NS=NF;c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;a[c]=0}}}if($NF==2){comp=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){comp++;if($i!=$j){a[comp]++}}}}}END{c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;sum+=a[c]/callsites}}print sum/c}')
echo "----> population Fst"
popFstColumn=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{for(i=1;i<=NF;i++){if($i == "pop_FST"){print i}}}')
popFstNalleleColumns=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i == "pop_FST" || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
pFst=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${popFstColumn} | awk '{if($1!="-nan"){sum+=$1;c+=1}}END{fst=sum/c; print fst}')
pFstBi=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${popFstNalleleColumns} | awk '{if($2!="-nan" && $1==2){sum+=$2;c+=1}}END{fst=sum/c; print fst}')
echo "----> sex Fst"
popFstColumn=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{for(i=1;i<=NF;i++){if($i == "sex_FST"){print i}}}')
popFstNalleleColumns=$(zcat ${chrsSiteMatrices[0]} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i == "sex_FST" || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
sFst=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${popFstColumn} | awk '{if($1!="-nan"){sum+=$1;c+=1}}END{fst=sum/c; print fst}')
sFstBi=$(zcat ${strchrsSiteMatrices} | grep -Pv '^chr\tst' | cut -f${popFstNalleleColumns} | awk '{if($2!="-nan" && $1==2){sum+=$2;c+=1}}END{fst=sum/c; print fst}')


echo "----> PCA axis values"  
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_all.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PCall
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_polymorphic_all_all.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PCpolym
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_biallelic_all.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PCbial
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_exon_all_all_all.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PCexon
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_0.5to0.6.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PC0.5to0.6
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_0.6to0.7.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PC0.6to0.7
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_0.7to0.8.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PC0.7to0.8
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_0.8to0.9.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PC0.8to0.9
awk '{if(NR==FNR){a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6;next} print a[$1]}' <(zcat ${PCAbasename}PCA_all_all_all_0.9to1.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.PC0.9to1

echo "----> UMAP axis values"  
awk '{if(NR==FNR){a[$1]=$2"\t"$3;next} print a[$1]}' <(zcat ${UMAPbasename}UMAP_all_all_all_all.txt.gz | sed 's/"//g' | tail -n +2) <(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') >  ${OutBasename}.UMAPall

echo "----> Print output"
hetheader=$(head -1 ${OutBasename}.het | cut -f2- )
str="Sample hetsites homsites biallhetsites minalleles singletons doubletons tripletons hetvar hetcall varsites callsites nBiallSites pi piBiallSites eHet eHom Fstat meanpopFst meanpopFstBiallSites meansexFst meansexFstBiallSites "${hetheader}
for s in all polym bial exon 0.5to0.6 0.6to0.7 0.7to0.8 0.8to0.9 0.9to1
do
	str=${str}" PC1_"${s}" PC2_"${s}" PC3_"${s}" PC4_"${s}" PC5_"${s}
done
str=${str}" UMAP1_all UMAP2_all"
echo ${str} 
echo ${str} | sed 's/ /\t/g' > ${OutBasename}.tab
for i in ${!VCFsortedSamples[@]} # for each sample
do
	hetvar=$(echo "${hetsites[$i]}/${varsites}*100" | bc -l | sed 's/^\./0./g') # heterozygosity = het sites / het+hom sites (variable sites) *100
	hetcal=$(echo "${hetsites[$i]}/${callsites}*100" | bc -l | sed 's/^\./0./g') # heterozygosity = het sites / callable sites *100
	Fstat=$(echo "1-${biallhetsites[$i]}/${eHet}" | bc -l | sed 's/^\./0./g')
	echo ${VCFsortedSamples[$i]} ${hetsites[$i]} ${homsites[$i]} ${biallhetsites[$i]} ${minalleles[$i]} ${singletons[$i]} ${doubletons[$i]} ${tripletons[$i]} ${hetvar} ${hetcal} ${varsites} ${callsites} ${nBiallSites} ${pi} ${piBi} ${eHet} ${eHom} ${Fstat} ${pFst} ${pFstBi} ${sFst} ${sFstBi} | sed 's/ /\t/g' >> ${OutBasename}.tmp
done
paste ${OutBasename}.tmp <(cut -f2- ${OutBasename}.het | tail -n +2) > ${OutBasename}.tmp2
mv ${OutBasename}.tmp2 ${OutBasename}.tmp
for s in all polym bial exon 0.5to0.6 0.6to0.7 0.7to0.8 0.8to0.9 0.9to1
do
	paste ${OutBasename}.tmp <(cat ${OutBasename}.PC${s}) > ${OutBasename}.tmp2
	mv ${OutBasename}.tmp2 ${OutBasename}.tmp
done
paste ${OutBasename}.tmp <(cat ${OutBasename}.UMAPall) > ${OutBasename}.tmp2
cat ${OutBasename}.tmp2 >> ${OutBasename}.tab
gzip ${OutBasename}.tab

rm ${OutBasename}.tmp* ${OutBasename}.het ${OutBasename}.PC* ${OutBasename}.UMAP* 2> ~/null 

