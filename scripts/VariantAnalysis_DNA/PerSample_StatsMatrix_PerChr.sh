#!/bin/bash


MetadataFile=$1
chrVCF=$2
genotypeFileBasename=${3%.chr*}
CallableRegions=$4
SamplesOrderInVCF=$5
SiteMatrix=$6
OutBasename=${7%.tab}
chr=$8

###################################
### To debug
#MetadataFile="metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt"
#chrVCF="results/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/ShortVariants_HardCallableFiltered.chr19.vcf.gz"
#genotypeFileBasename="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered"
#CallableRegions="results/VariantCalling_DNA/7_join_extra_callable_regions/ExtraCallableRegions.bed.gz"
#SamplesOrderInVCF="metadata/SamplesOrderInVCF.chr19.txt" 
#SiteMatrix="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr19.tab.gz"
#OutBasename="results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.chr19"
#chr="chr19"
###################################

mkdir -p $(dirname ${OutBasename})

VCFsortedSamples=( $(cat ${SamplesOrderInVCF} | sed 's/\t/\n/g') )
rm ${OutBasename}.tmp ${OutBasename}.het


echo "----> vcftools heterozygosity"  # only considers PASS variants
vcftools --gzvcf ${chrVCF} --remove-filtered-all --het --out ${OutBasename}


echo "----> Expected heterozygous & homozygous sites per sample" 
# !!!! Only considering biallelic variants
# Column with allele frequency data
freqcol=$(zcat ${SiteMatrix} | head -1 | awk '{for(i=0;i<=NF;i++){if($i=="freq"){print i}}}')
# total number of biallelic sites --> same for all samples
nBiallSites=$(zcat ${SiteMatrix} | tail -n +2 | cut -f${freqcol} | awk '{n=split($1,a,";"); if(n==2){num++}}END{print num}')
# Expected number of heterozygous sites per sample:  --> same for all samples
# sum(2pq)=sum(2*numPalleles/(numsamples*2)*numQalleles/(numsamples*2)) for all biallelic sites
eHet=$(zcat ${SiteMatrix} | tail -n +2 | cut -f${freqcol} | awk -v nsamp=${#VCFsortedSamples[@]} '{n=split($1,a,";"); if(n==2){ehet=2; for(i in a){split(a[i],b,":");ehet*=b[2]/nsamp/2;} printf "%.6f\n", ehet;}}' | awk '{sum+=$1;}END{print sum}')
# Expected number of homozygous sites per sample: --> same for all samples
# sum(p^2+q^2)=sum(numPalleles*numPalleles/(nsamp*nsamp*4)+numQalleles*numQalleles/(nsamp*nsamp*4)) for all biallelic sites
eHom=$(zcat ${SiteMatrix} | tail -n +2 | cut -f${freqcol} | awk -v nsamp=${#VCFsortedSamples[@]} '{n=split($1,a,";"); if(n==2){ehom=0; for(i in a){split(a[i],b,":");ehom+=b[2]*b[2]/(nsamp*nsamp*4);} printf "%.6f\n", ehom;}}' | awk '{sum+=$1;}END{print sum}')

echo "----> Observed heterozygous & homozygous sites per sample"
# number of observed heterozygous sites per each sample
hetsites=( $(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz | awk '{for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){str=str" "a[i]} print str}' | sed 's/^ //g' ) )
# number of observed homozygous sites per each sample
homsites=( $(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz | awk '{num++;for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){a[i]=num-a[i]; str=str" "a[i]} print str}' | sed 's/^ //g' ) )
# number of variant sites (at least two alleles present) --> same for all samples
varsites=$(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz | cut -f5- | awk '{for(i=1;i<=NF;i++){s+=$i}if(s>0){print s}s=0}' | wc -l )
# number of callable sites (check r7_callable_regions for details) --> same for all samples
callsites=$(zcat ${CallableRegions} | awk -v chr=${chr} '{if($1==chr){len+=$3-$2;next}}END{print len}' )
# number of observed biallelic heterozygous sites per each sample
biallhetsites=( $(awk '{if(NR==FNR){n=split($1,a,";"); if(n==2){blines[FNR]=1;}next} if(blines[FNR]){print $0}}' <(zcat ${SiteMatrix} | cut -f${freqcol} | tail -n +2) <(zcat ${genotypeFileBasename}.${chr}.0hom1het.bed.gz) | awk '{for(i=5;i<=NF;i++){a[i]+=$i;}}END{str=""; for(i in a){str=str" "a[i]} print str}' | sed 's/^ //g' ) )

echo "----> Number of minor alleles per sample"
GenoFreqColumns=$(zcat ${SiteMatrix} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/ || $i=="freq"){str=str","i}}print str}' | sed 's/^,//g')
minalleles=( $(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of singletons per sample"
singletons=( $(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-1)){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of doubletons per sample" # only biallelic
doubletons=( $(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-2) && length(a)==4){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )
echo "----> Number of tripletons per sample" # only biallelic
tripletons=( $(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoFreqColumns} | awk '{split($NF,a,";|:"); nchr=(NF-1)*2; if(a[length(a)]==(nchr-3) && length(a)==4){str=""; for(i=1;i<NF;i++){split($i,b,":");alt=0; for(j in b){if(b[j]!=a[length(a)-1]){alt++}} str=str" "alt} print str}}' | sed 's/^ //g' | awk '{for(i=1;i<=NF;i++){a[i]+=$i}}END{for(i=1;i<=NF;i++){print a[i]}}') )

# average pairwise differences
GenoColumns=$(zcat ${SiteMatrix} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/){str=str","i}}print str}' | sed 's/^,//g')
pi=$(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoColumns} | sed 's/:/\t/g' | awk -v callsites=${callsites} '{if(NR==1){NS=NF;c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;a[c]=0}}}comp=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){comp++;if($i!=$j){a[comp]++}}}}END{c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;sum+=a[c]/callsites}}print sum/c}')
GenoNalleleColumns=$(zcat ${SiteMatrix} | head -1 | awk '{str="";for(i=1;i<=NF;i++){if($i ~ /geno/ || $i=="numalleles"){str=str","i}}print str}' | sed 's/^,//g')
piBi=$(zcat ${SiteMatrix} | grep -Pv '^chr\tst' | cut -f${GenoNalleleColumns} | sed 's/:/\t/g' | awk -v callsites=${callsites} '{if(NR==1){NS=NF;c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;a[c]=0}}}if($NF==2){comp=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){comp++;if($i!=$j){a[comp]++}}}}}END{c=0;for(i=1;i<NS;i++){for(j=i+1;j<=NS;j++){c++;sum+=a[c]/callsites}}print sum/c}')



echo "----> Print output"
echo "Sample hetsites homsites biallhetsites minalleles singletons doubletons tripletons hetvar hetcall varsites callsites nBiallSites pi piBiallSites eHet eHom Fstat" | sed 's/ /\t/g' > ${OutBasename}.tmp
for i in ${!VCFsortedSamples[@]} # for each sample
do
	hetvar=$(echo "${hetsites[$i]}/${varsites}*100" | bc -l | sed 's/^\./0./g') # heterozygosity = het sites / het+hom sites (variable sites) *100
	hetcal=$(echo "${hetsites[$i]}/${callsites}*100" | bc -l | sed 's/^\./0./g') # heterozygosity = het sites / callable sites *100
	Fstat=$(echo 1-${biallhetsites[$i]}/${eHet} | bc -l | sed 's/^\./0./g')
	echo ${VCFsortedSamples[$i]} ${hetsites[$i]} ${homsites[$i]} ${biallhetsites[$i]} ${minalleles[$i]} ${singletons[$i]} ${doubletons[$i]} ${tripletons[$i]} ${hetvar} ${hetcal} ${varsites} ${callsites} ${nBiallSites} ${pi} ${piBi} ${eHet} ${eHom} ${Fstat} | sed 's/ /\t/g' >> ${OutBasename}.tmp
done
echo "VCFsortedSamples ${VCFsortedSamples[$i]}"
echo "hetsites ${hetsites[$i]}"
echo "homsites ${homsites[$i]}"
echo "biallhetsites ${biallhetsites[$i]}"
echo "minalleles ${minalleles[$i]}"
echo "hetvar ${hetvar}"
echo "hetcal ${hetcal}"
echo "varsites ${varsites}"
echo "callsites ${callsites}"
echo "nBiallSites ${nBiallSites}"
echo "pi ${pi}"
echo "piBi ${piBi}"
echo "eHet ${eHet}"
echo "eHom ${eHom}"
echo "Fstat ${Fstat}"

paste ${OutBasename}.tmp <(cut -f2- ${OutBasename}.het) > ${OutBasename}.tab

rm ${OutBasename}.tmp ${OutBasename}.het