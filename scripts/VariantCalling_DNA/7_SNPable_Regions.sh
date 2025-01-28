#!/bin/bash

# Following http://lh3lh3.users.sourceforge.net/snpable.shtml as of 5th of July 2022


GenomeIn=$1
GenomeOut=$2
splitfa=$3
gen_raw_mask=$4
gen_mask=$5
apply_mask_s=$6

OutBaseName=$(echo ${GenomeOut%_genome_masked.fa})
mkdir -p $(dirname ${GenomeOut})
rm ${OutBaseName}_kmerseq_aln.sam.gz 2> ~/null

echo "### splitfa"; >&2 echo "### splitfa"
${splitfa} ${GenomeIn} 35 > ${OutBaseName}_kmerseq.fa
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### bwa aln"; >&2 echo "### bwa aln"
bwa aln -R 1000000 -O 3 -E 3 -t 4 ${GenomeIn} ${OutBaseName}_kmerseq.fa > ${OutBaseName}_kmerseq_aln.sai  
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### bwa samse"; >&2 echo "### bwa samse"
bwa samse ${GenomeIn} ${OutBaseName}_kmerseq_aln.sai ${OutBaseName}_kmerseq.fa > ${OutBaseName}_kmerseq_aln.sam
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### gzip"; >&2 echo "### gzip"
gzip ${OutBaseName}_kmerseq_aln.sam
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### gen_raw_mask"; >&2 echo "### gen_raw_mask"
gzip -dc ${OutBaseName}_kmerseq_aln.sam.gz | ${gen_raw_mask} > ${OutBaseName}_rawMask_35.fa  
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### gen_mask"; >&2 echo "### gen_mask"
${gen_mask} -l 35 -r 0.5 ${OutBaseName}_rawMask_35.fa > ${OutBaseName}_mask_35_50.fa
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### apply_mask_s"; >&2 echo "### apply_mask_s"
${apply_mask_s} ${OutBaseName}_mask_35_50.fa ${GenomeIn} > ${GenomeOut} 
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"

echo "### Extract masked regions"; >&2 echo "### Extract masked regions"
cat ${OutBaseName}_mask_35_50.fa | \
awk '{
	if($1 ~ />/){
		print seq
		print $0
		seq=""
	}else{
		seq=seq$0
	}
}END{print seq}' | tail -n +2 | \
awk -F ' ' '{
	if($1 ~ /^>/){
		chr=$1
		ins=0
		next
	}
	split($1, s, "")
	for(i in s){
		if(s[i]!=3 && ins==0){
			st=i
			ins=1
		}
		if(s[i]==3 && ins==1){
			end=i-1
			ins=0
			print chr"\t"st"\t"end
		}
	}
}' > ${OutBaseName}_unSNPableReg_35_50.bed
status=$?
[ $status -eq 0 ] && echo "Done" || exit "Failed ($status)"











