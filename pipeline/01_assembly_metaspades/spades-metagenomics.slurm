#!/bin/bash
# Per-sample metaSPAdes assembly, contigs databases (HMMs, KOfams, COGs, SCG taxonomy)
# and read-recruitment profiles with the anvi'o metagenomics workflow; settings in
# spades-metagenomics-config.json, samples in all_samples.txt.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=168
#SBATCH --mem-per-cpu=4591
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file


anvi-run-workflow -w metagenomics \
                  -c spades-metagenomics-config.json 
