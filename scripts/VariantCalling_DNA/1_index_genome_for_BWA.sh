#!/bin/bash

genome=$1
output=$2

mkdir -p $(dirname ${output})

bwa index ${genome}
samtools faidx ${genome}
picard CreateSequenceDictionary REFERENCE=${genome} OUTPUT=${output}
