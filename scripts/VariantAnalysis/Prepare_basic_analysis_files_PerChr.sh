#!/bin/bash


InputChrVCF=${1%.gz}
OutBasename=${2%.genotype.bed.gz}

#####
# Dev & debug
#InputChrVCF="results/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/ShortVariants_HardCallableFiltered.chr19.vcf.gz"
#OutBasename="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.bed.gz"




mkdir -p $(dirname ${OutBasename})
rm ${OutBasename}* 2> ~/null

echo "----> Extrancting genotypes and variant info with gatk VariantsToTable"
echo "   ---> Extract vcf"
gzip -dk ${InputChrVCF}.gz

echo "   ---> Extract genotypes" # only PASS variants are printed
gatk VariantsToTable -V ${InputChrVCF} -O ${OutBasename}.genotype.tmp -GF GT
cat ${OutBasename}.genotype.tmp | sed 's/[|\/]/:/g' | sed '1 s/\.GT//g' > ${OutBasename}.genotype

echo "   ---> Extract variants" # only PASS variants are printed
gatk VariantsToTable -V ${InputChrVCF} -O ${OutBasename}.tmp.sites -F CHROM -F POS -F ID -F REF -F ALT
cat ${OutBasename}.tmp.sites | tail -n +2 | awk 'BEGIN{num=1}{print $1"\t"$2"\t"$2"\t"num"\t"$4"\t"$5;num++}' > ${OutBasename}.bed
echo "   ---> Paste genotype and vairant info"
paste <(cut -f1,2,3,4 ${OutBasename}.bed) <(tail -n +2 ${OutBasename}.genotype) > ${OutBasename}.genotype.bed

echo "----> Create genotype numeric matrix (#alleles-1 lines for each site) for PCA plotting"
cat ${OutBasename}.genotype.bed | sed 's/:/ /g' | awk -F' ' '{for(i=5;i<=NF;i+=2){a[$i][i]++; j=i+1; a[$j][i]++; b[$j]++;	b[$i]++;} maxal="XX"; maxval=0;	for(al in b){if(b[al]>maxval){ maxal = al; maxval = b[al];}}for(al in a){if(maxal!=al){	str="";	for(i=5;i<=NF;i+=2){if(a[al][i]){n=a[al][i]}else{n=0}str=str" "n} print NR""str}} delete a; delete b;}' > ${OutBasename}.genotype.numericmatrix

echo "----> Heterozygosity format conversion"
cat ${OutBasename}.genotype.bed | cut -f5- | awk -F '\t' '{str=""; for(i=1;i<=NF;i++){split($i,a,":"); if(a[1]==a[2]){str=str"\t0"}else{str=str"\t1"}} print str}' | sed 's/^\t//g' > ${OutBasename}.0hom1het
paste <(cut -f1,2,3,4 ${OutBasename}.bed) ${OutBasename}.0hom1het > ${OutBasename}.0hom1het.bed

gzip ${OutBasename}.genotype.bed
gzip ${OutBasename}.genotype.numericmatrix
gzip ${OutBasename}.0hom1het.bed

rm ${OutBasename}.genotype.tmp 2> ~/null
rm ${InputChrVCF}  2> ~/null
rm ${OutBasename}.tmp.sites  2> ~/null
rm ${OutBasename}.genotype 2> ~/null
rm ${OutBasename}.bed 2> ~/null
rm ${OutBasename}.0hom1het 2> ~/null

