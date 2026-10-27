#!/bin/bash
# Write the anvi'o samples file (sample, r1, r2) for all paired read files under
# $METAG/raw_reads/, naming each sample after its folder.
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --job-name=generate_samples
#SBATCH --output=generate_samples.out
#SBATCH --error=generate_samples.err

# Define the input data directory in the shared folder
DATA_DIR="$METAG/raw_reads/"
# Make sure the output is created in the current working directory
OUTPUT_FILE="$(pwd)/samples.txt"

# Rename all .fq.gz files to .fastq.gz
echo "Renaming all .fq.gz files to .fastq.gz..."
# Look for files in the shared directory instead of current directory
find "$DATA_DIR" -type f -name "*.fq.gz" | while read -r file; do
    newfile="${file%.fq.gz}.fastq.gz"
    echo "Renaming $file to $newfile"
    mv "$file" "$newfile"
done

# Create samples.txt with header
# Use explicit output file path
echo -e "sample\tr1\tr2" > "$OUTPUT_FILE"

# Find all directories containing fastq.gz files
echo "Searching for fastq.gz files and pairing them..."
# Look for files in the shared directory instead of current directory
# Use absolute paths instead of relative paths
find "$DATA_DIR" -type f -name "*1.fastq.gz" | while read -r r1; do
    # Construct the expected R2 file by replacing *1 with *2
    r2="${r1%1.fastq.gz}2.fastq.gz"
    
    # Check if the R2 file exists
    if [ -f "$r2" ]; then
        # Get the sample name from the immediate parent directory of the read files
        sample=$(basename "$(dirname "$r1")")
        
        # Use absolute paths instead of relative paths
        r1_path="$r1"
        r2_path="$r2"
        
        # Write to samples.txt
        echo -e "$sample\t$r1_path\t$r2_path" >> "$OUTPUT_FILE"
    else
        echo "Warning: No matching R2 file found for $r1"
    fi
done

echo "Done! Check $OUTPUT_FILE for the results."
