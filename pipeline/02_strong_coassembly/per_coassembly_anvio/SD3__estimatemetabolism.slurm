#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

# # Clean up loaded modules to avoid conflicts
# module purge

# # Load necessary modules (adjust as needed based on the system)
# module load Conda/22.9.0

# Activate your Conda environment
# conda activate anvio-dev

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run anvi-run-kegg-kofams for each .db file
for db in *.db; do
    anvi-run-kegg-kofams -c "$db" --log-bitscores
done

# Generate an external genomes file
anvi-script-gen-genomes-file --input-dir ./ \
                             --output-file external_genomes.txt

# Estimate metabolism using the external genomes file
anvi-estimate-metabolism -e external_genomes.txt \
                         -O metabolism_estimated \
                         --matrix-format
