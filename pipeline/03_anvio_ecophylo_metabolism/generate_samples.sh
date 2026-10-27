#!/bin/bash
# Write samples.txt (sample, r1, r2) for the paired read files below the current
# folder, naming each sample after its folder.

# Rename all .fq.gz files to .fastq.gz
echo "Renaming all .fq.gz files to .fastq.gz..."
find . -type f -name "*.fq.gz" | while read -r file; do
    newfile="${file%.fq.gz}.fastq.gz"
    echo "Renaming $file to $newfile"
    mv "$file" "$newfile"
done

# Create samples.txt with header
echo -e "sample\tr1\tr2" > samples.txt

# Find all directories containing fastq.gz files
echo "Searching for fastq.gz files and pairing them..."
# Use find to locate all _1.fastq.gz files recursively
find . -type f -name "*_1.fastq.gz" | while read -r r1; do
    # Construct the expected R2 file by replacing _1 with _2
    r2="${r1%_1.fastq.gz}_2.fastq.gz"
    
    # Check if the R2 file exists
    if [ -f "$r2" ]; then
        # Get the sample name from the immediate parent directory of the read files
        sample=$(basename "$(dirname "$r1")")
        
        # Clean up the paths by removing double slashes and ensuring relative path format
        r1_path=$(echo "./${r1#./}" | sed 's|//|/|g')
        r2_path=$(echo "./${r2#./}" | sed 's|//|/|g')
        
        # Write to samples.txt
        echo -e "$sample\t$r1_path\t$r2_path" >> samples.txt
    else
        echo "Warning: No matching R2 file found for $r1"
    fi
done

echo "Done! Check samples.txt for the results."
