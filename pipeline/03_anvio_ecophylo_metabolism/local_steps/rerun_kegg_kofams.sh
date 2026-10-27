#!/bin/bash

# Re-annotate the NVNN bin contigs databases with KOfam against one KEGG snapshot
# (TARGET_HASH), skipping databases already annotated with it.

# Set the base directory containing your bins/MAGs
BASE_DIR="METABAT_MAGs"

# Number of threads to use (adjust based on your system)
THREADS=4

# Target KEGG hash (the new version)
TARGET_HASH="66e53d49e65a"

# Log file
LOG_FILE="kegg_rerun_$(date +%Y%m%d_%H%M%S).log"

echo "Starting KEGG KOfam re-annotation at $(date)" | tee -a "$LOG_FILE"
echo "Target KEGG hash: $TARGET_HASH" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

# Counter for progress tracking
total=0
processed=0
skipped=0
failed=0

# Count total number of CONTIGS.db files
for dir in "$BASE_DIR"/NONANO_Bin_* "$BASE_DIR"/NONANO_MAG_*; do
    if [ -d "$dir" ]; then
        if [ -f "$dir/CONTIGS.db" ]; then
            ((total++))
        fi
    fi
done

echo "Found $total CONTIGS.db files to process" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Function to check KEGG hash in database
check_kegg_hash() {
    local db_path=$1
    # Extract the hash from anvi-db-info output
    local current_hash=$(anvi-db-info "$db_path" 2>/dev/null | grep -A 20 "gene_function_calls" | grep "modules_db_hash" | awk '{print $NF}')
    echo "$current_hash"
}

# Process each bin/MAG directory
for dir in "$BASE_DIR"/NONANO_Bin_* "$BASE_DIR"/NONANO_MAG_*; do
    if [ -d "$dir" ]; then
        contigs_db="$dir/CONTIGS.db"
        
        if [ -f "$contigs_db" ]; then
            ((processed++))
            bin_name=$(basename "$dir")
            
            # Check if this database already has the correct KEGG version
            current_hash=$(check_kegg_hash "$contigs_db")
            
            if [ "$current_hash" == "$TARGET_HASH" ]; then
                echo "[$processed/$total] SKIPPING: $bin_name (already has correct KEGG version)" | tee -a "$LOG_FILE"
                ((skipped++))
                continue
            fi
            
            echo "[$processed/$total] Processing: $bin_name" | tee -a "$LOG_FILE"
            
            # Run anvi-run-kegg-kofams
            if anvi-run-kegg-kofams -c "$contigs_db" -T "$THREADS" --just-do-it 2>&1 | tee -a "$LOG_FILE"; then
                echo "  done: $bin_name" | tee -a "$LOG_FILE"
            else
                echo "  FAILED: $bin_name" | tee -a "$LOG_FILE"
                ((failed++))
            fi
            
            echo "" | tee -a "$LOG_FILE"
        fi
    fi
done

echo "================================================" | tee -a "$LOG_FILE"
echo "Completed at $(date)" | tee -a "$LOG_FILE"
echo "Total checked: $processed" | tee -a "$LOG_FILE"
echo "Skipped (already done): $skipped" | tee -a "$LOG_FILE"
echo "Newly processed: $((processed - skipped - failed))" | tee -a "$LOG_FILE"
echo "Failed: $failed" | tee -a "$LOG_FILE"
echo "Log saved to: $LOG_FILE" | tee -a "$LOG_FILE"
