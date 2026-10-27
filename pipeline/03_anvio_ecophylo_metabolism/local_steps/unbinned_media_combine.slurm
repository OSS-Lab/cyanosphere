#!/bin/bash
# Concatenate the per-sample module tables of the media experiment. Note: tail on
# several files inserts '==> file <==' lines; data/raw_kegg was exported from the
# per-sample files by 07_export_raw_tables.py, not from this table.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

(head -1 NVNN_P10_R1_modules.txt; tail -n +2 *_modules.txt) > combined_modules.txt
