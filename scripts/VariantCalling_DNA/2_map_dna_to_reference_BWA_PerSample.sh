#!/bin/bash


sample=$1
lane=$2
platform=$3
genome=$4
fastq1=$5
fastq2=$6
readgroupsfile=$7
bam=$8
bammarkdup=$9

mkdir -p $(dirname ${bam})
rm -r ${bammarkdup}* 2> ~/null

lanewolettter=$(echo ${lane} | sed 's/^L//g')
readline=$(grep -P "@RG\tID:${lanewolettter}.${platform}\tSM:${sample}" ${readgroupsfile} | sed 's/\t$//g' | sed 's/\t/\\t/g')
echo ${readline}

bwa mem -R $(echo ${readline}) ${genome} ${fastq1} ${fastq2} | samtools sort -O bam - > ${bam}

gatk MarkDuplicatesSpark -I ${bam} -O ${bammarkdup}

