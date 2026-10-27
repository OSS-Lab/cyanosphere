#!/bin/bash
# KEGG modules with coverage and copy number on the unbinned contigs of one media
# sample (NVNN_P10_R4 shown); run once per sample with the name changed.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

anvi-estimate-metabolism -c ../03_CONTIGS/NVNN_P10_R4-contigs.db -p ../05_ANVIO_PROFILE/NVNN_P10_R4/NVNN_P10_R4/PROFILE.db -O NVNN_P10_R4 --add-coverage --add-copy-number
