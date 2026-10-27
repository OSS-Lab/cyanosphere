#!/bin/bash
# Dereplicate the MAGs of all co-assemblies with fastANI (anvi-dereplicate-genomes),
# ANI threshold 0.90, representative chosen by quality score (Qscore).
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

anvi-dereplicate-genomes -e external_genomes.txt -o fastANI --program fastANI --similarity-threshold 0.90 --representative-method Qscore -T 30
