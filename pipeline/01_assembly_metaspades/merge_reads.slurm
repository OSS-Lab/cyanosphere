#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=4000
#SBATCH --time=48:00:00

# Define paths
INPUT_DIR="/springbrook/share/lifesci/lfsnhf/metagenomics/raw_reads/WAMOS"
OUTPUT_DIR="${INPUT_DIR}/merged_reads"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Process each sample
for SAMPLE in W_07
do
    echo "Processing sample: ${SAMPLE}"
    
    # Create sample output directory
    mkdir -p ${OUTPUT_DIR}/${SAMPLE}
    
    # Count how many read files we have
    read_count=$(ls ${INPUT_DIR}/${SAMPLE}/*_1.fq.gz 2>/dev/null | wc -l)
    
    if [[ ${read_count} -gt 1 ]]
    then
        echo "Multiple read files found for ${SAMPLE}, concatenating..."
        zcat ${INPUT_DIR}/${SAMPLE}/*_1.fq.gz | gzip > ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_R1.fq.gz
        zcat ${INPUT_DIR}/${SAMPLE}/*_2.fq.gz | gzip > ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_R2.fq.gz
    else
        echo "Single read file found for ${SAMPLE}, copying..."
        cp ${INPUT_DIR}/${SAMPLE}/*_1.fq.gz ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_R1.fq.gz
        cp ${INPUT_DIR}/${SAMPLE}/*_2.fq.gz ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_R2.fq.gz
    fi
    
    echo "Completed processing ${SAMPLE}"
done

echo "All samples processed successfully!"
