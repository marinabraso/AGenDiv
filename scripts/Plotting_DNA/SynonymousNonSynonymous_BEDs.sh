#!/bin/bash


GTF=$1
genome=$2
OutSyn=$3
OutNonSyn=$4
SynNonSyn=$5
scriptSnSfromSeq=$6
scriptProtfromTransc=$7

##############
# To debug / develop
#GTF="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf"
#genome="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
#scriptSnSfromSeq="scripts/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromSeq.pl"
#scriptProtfromTransc="scripts/VariantAnalysis_DNA/ProteinSeq_fromTranscriptSeq.pl"
##############


mkdir -p $(dirname ${OutSyn})


echo "----> Select longest transcript per gene"
cat ${GTF} | awk '{if($3 == "CDS"){print $0}}' | cut -f1,4,5,9 | sed 's/gene_id "\(\S\+\)";.*transcript_id "\(\S\+\)";.*/\1\t\2/g' | sort -k1,1 -k5,5 -k2,2V | awk -F'\t' '{if(t==$5){l+=$3-$2;next}else{print c,l,g,t;c=$1;l=$3-$2;g=$4;t=$5}}END{print c,l,g,t;}' | tail -n +2 | tr ' ' '\t' | sort -k1,1 -k3,3 -k2,2V | awk '{if(g==$3){t=$4}else{print g,t;g=$3;t=$4}}END{print g,t;}' | tail -n +2 | tr ' ' '\t' > $(dirname ${OutSyn})/LonguestTranscriptPerGene2.tab

echo "----> Extract CDS bed of longest transcript per gene"
awk -F'\t' '{if(NR==FNR){a[$1]=$2;next}if(a[$5]==$6){print $0}}' $(dirname ${OutSyn})/LonguestTranscriptPerGene2.tab <(cat ${GTF} | awk '{if($3=="CDS"){print $0}}' | cut -f1,4,5,7,9 | sed 's/gene_id "\(\S\+\)";.*transcript_id "\(\S\+\)";.*/\1\t\2/g') | sort -k1,1 -k5,5 -k2,2V | awk '{print $1"\t"$2"\t"$3"\t"$4":"$5":"$6}' > $(dirname ${OutSyn})/LonguestTranscript_CDS.bed

echo "----> Extract CDS sequence from genome per exon"
bedtools getfasta -name -fi ${genome} -bed <(cat $(dirname ${OutSyn})/LonguestTranscript_CDS.bed | awk '{st=$2-1; print $1"\t"st"\t"$3"\t"$4}') -bedOut | sed 's/:/\t/g' | awk '{st=$2+1; print $1"\t"st"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $(dirname ${OutSyn})/CDS_exon_sequences.bed

echo "----> Join transcript sequence"
cat $(dirname ${OutSyn})/CDS_exon_sequences.bed | awk '{if($5 == g){seq=seq$7}else{print g"\t"sense"\t"seq; g=$5; seq =$7; sense=$4}}END{print g"\t"sense"\t"seq;}' | tail -n +2 > $(dirname ${OutSyn})/CompleteCDS_sequence.tab

echo "----> Get Syn and nonSyn positions from transcript sequences"
# check 5 transcript with N!!!
cat $(dirname ${OutSyn})/CompleteCDS_sequence.tab | grep -v "N" | awk -v scr=${scriptSnSfromSeq} '{cmd=scr"\t"$3"\t"$2; cmd | getline nseq; close(cmd); print $1"\t"nseq}' > $(dirname ${OutSyn})/SnS_CompleteCDS_sequence.tab
#cat $(dirname ${OutSyn})/CompleteCDS_sequence.tab | grep -v "N" | awk -v scr=${scriptProtfromTransc} '{cmd=scr"\t"$3"\t"$2; cmd | getline nseq; close(cmd); print $1"\t"nseq}' > $(dirname ${OutSyn})/Prot_sequence.tab


echo "----> Filter out SnS sequences that are not the same lenght as the transcript"
## Step to be checked:
# differences in length come from errors in the transcript sequence
# N values within the transcript or len%3!=0
# only a few cases show problems so I'm skipping them for now
awk '{if(NR==FNR){a[$1]=$3;next}l1=split(a[$1],a1,"");l2=split($2,a2,"");if(l1==l2){print $0}}' $(dirname ${Out})/CompleteCDS_sequence.tab $(dirname ${Out})/SnS_CompleteCDS_sequence.tab > $(dirname ${Out})/SnS_CompleteCDS_sequence_filtered.tab

echo "----> Divide SnS sequences by exon"
awk '{if(NR==FNR){a[$1]=$2;next} el=split($7,eCDS,""); if(a[$5]){eSnS=substr(a[$5],1,el);a[$5]=substr(a[$5],el);print $0"\t"eSnS}}' <(cat $(dirname ${Out})/SnS_CompleteCDS_sequence_filtered.tab) <(cat $(dirname ${Out})/CDS_exon_sequences.bed) > $(dirname ${Out})/SnS_CDS_exon_sequences.bed

echo "----> Convert to per site bed"
cat $(dirname ${Out})/SnS_CDS_exon_sequences.bed | awk '{split($8,a,""); for(i=$2;i<=$3;i++){print $1"\t"i"\t"i"\t"a[i-$2+1]}}' > ${Out}

echo "----> Convert to per site bed"
bedtools intersect -wa -a <(cat ${Out} | grep '^chr') -b <(zcat ${CallableRegions}) > ${OutCallab}



