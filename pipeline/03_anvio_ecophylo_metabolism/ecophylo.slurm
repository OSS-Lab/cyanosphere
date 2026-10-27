#!/bin/bash
# EcoPhylo (anvi'o): ribosomal-protein clusters at 97% identity (MMseqs2 --min-seq-id 0.97) from contigs and MAGs,
# reads mapped back to the cluster representatives; config2.json, marker in hmm-list.txt.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

anvi-run-workflow -w ecophylo -c config2.json 
