#!/bin/bash


genome=$1
fastq1=$2
fastq2=$3
bam=$4
threads=$5

mkdir -p $(dirname ${bam})

minimap2 -ax sr -B4 -O6,24 -E2,1 -t ${threads} ${genome} ${fastq1} ${fastq2} | samtools sort -O bam -o ${bam} -
		# -B4 → lowers mismatch penalty (more tolerant)
		# -O6,24 → gap open penalties (affects indels)
		# -E2,1 → gap extension penalties
