#!/bin/bash
# KEGG module completeness and copy number for all metagenomes in metagenomes.txt,
# as matrices.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

anvi-estimate-metabolism -M metagenomes.txt --matrix-format --add-copy-number
