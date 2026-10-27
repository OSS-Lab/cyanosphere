#!/bin/bash
# Concatenate the lanes of each sample into one read pair under input_reads/<sample>/,
# as STRONG expects; edit the sample list before running.

mkdir input_reads
for SAMPLE in W_08 W_09
do
    mkdir input_reads/${SAMPLE}
    files=( raw_data/${SAMPLE}/*_1.fq.gz )
    if [[ ${#files[@]} -gt 1 ]]
    then
        zcat raw_data/${SAMPLE}/*_1.fq.gz | gzip > input_reads/${SAMPLE}/${SAMPLE}_R1.fq.gz
        zcat raw_data/${SAMPLE}/*_2.fq.gz | gzip > input_reads/${SAMPLE}/${SAMPLE}_R2.fq.gz
    else
        cp raw_data/${SAMPLE}/*_1.fq.gz input_reads/${SAMPLE}/${SAMPLE}_R1.fq.gz
        cp raw_data/${SAMPLE}/*_2.fq.gz input_reads/${SAMPLE}/${SAMPLE}_R2.fq.gz
    fi
done
