#!/bin/bash



Rscript=$(echo $0 | sed 's/\.sh/.R/g')
MetadataFile=$1
chrInBasenames=$2
SamplesOrderInVCF=$3
FSTATFile=$4
chrs=$5
rm tmp/MaxNumAlleles.tbl 2> ~/null

### TMP
#chrs="18 19"
#chrInBasenames="results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr18.genotype.bed results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.chr19.genotype.bed"
###

mkdir -p $(dirname ${FSTATFile})

chrInBasenames=( $(echo ${chrInBasenames} | sed 's/.genotype.bed/\n/g') )
chrs=( $(echo ${chrs} | sed 's/ /\n/g') )


# Creating the FSTAT file for hierfstat to read
echo "----> Creating tmp file with the genotypes encoded in two digit numbers and a header for the popultion of each individual"
echo "      ----> header"
awk -F '\t' '{if(FNR==NR){a[$1]=$2;next;} str=""; for(i=1;i<=NF;i++){str=str"\t"a[$i]} print str}' <(tail -n +2 ${MetadataFile} | cut -f3,9 | sort | uniq | sed 's/Banyuls/1/g' | sed 's/Roscoff/2/g') ${SamplesOrderInVCF} | sed 's/^\t//g' > $(dirname ${FSTATFile})/Fstat_genotypes_transposed.tmp

echo "      ----> body"
for c in ${!chrs[@]}
do
	echo ${chrs[$c]}
	cat ${chrInBasenames[$c]}.genotype.bed | sed 's/:/\t/g' | awk -F' ' 'BEGIN{maxn=0}{num=1; for(i=6;i<=NF;i++){if(!a[$i]){a[$i]=num;num++}}if(num>maxn){maxn=num}str=""; for(i=6;i<=NF;i+=2){j=i+1; if(a[$i]>9){str=str"\t"a[$i]}else{str=str"\t0"a[$i]}if(a[$j]>9){str=str""a[$j]}else{str=str"0"a[$j]}} print str; delete a;}END{print maxn >> "tmp/MaxNumAlleles.tbl"}' | sed 's/^\t//g' >> $(dirname ${FSTATFile})/Fstat_genotypes_transposed.tmp
done

echo "----> First line of FSTAT"
echo "      ----> the number of samples, np"
np=$(awk -F '\t' '{print NF}' ${SamplesOrderInVCF})
echo "      ----> the number of loci, nl"
nl=$(wc -l $(dirname ${FSTATFile})/Fstat_genotypes_transposed.tmp | awk '{n=$1-1; print n}')
echo "      ----> the highest number used to label an allele, nu"
nu=$(awk -F '\t' '{if(NR==1){max=$1}if(max<$1){max=$1}}END{print max}' tmp/MaxNumAlleles.tbl)
if [ $nu -gt 99 ]
then
	echo "Maximum number of alleles > 99. FSTAT file needs to be changed to a 3-digit pattern"; 
	exit 1; 
fi
echo "      ----> the number of digits: 2"
ndigits=2 
echo "$np $nl $nu $ndigits" > ${FSTATFile%.gz}

echo "----> Loci names lines"
for ((l=1; l<=$nl; l++))
do
	echo "loc."$l >> ${FSTATFile%.gz}
done

echo "----> Transpose the previously build tmp file"
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' $(dirname ${FSTATFile})/Fstat_genotypes_transposed.tmp >> ${FSTATFile%.gz}


echo "----> gzip"
gzip ${FSTATFile%.gz}
#rm $(dirname ${FSTATFile})/Fstat_genotypes_transposed.tmp





