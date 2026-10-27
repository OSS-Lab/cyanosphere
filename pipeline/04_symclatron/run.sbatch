#!/bin/bash
# symclatron lifestyle classification (free-living / symbiont) of the 88 dereplicated
# MAGs from their predicted proteins; output data/symclatron_results.tsv.
#SBATCH --job-name=symclatron
#SBATCH --output=%x-%j.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
set -euo pipefail
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate symclatron
cd $HOME/symclatron_run
symclatron setup
symclatron classify --genome-dir proteins --output-dir results --threads 16
echo "DONE"
wc -l results/symclatron_results.tsv
