#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=33
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file

## Clean up loaded modules to avoid conflicts
# module purge

## Load necessary modules (adjust as needed based on the system)
# module load Conda/22.9.0

## Activate your Conda environment
# conda activate anvio-dev

## Set the number of threads for OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## Generate contigs.db from .fasta files
#for file in *.fasta; do
#   # Extract the base name (without the .fa extension)
#   base_name=$(basename "$file" .fasta)
#   # Run the anvi-gen-contigs-database command for each .fasta file
#   anvi-gen-contigs-database -f "$file" -o "${base_name}.db" --num-threads $SLURM_CPUS_PER_TASK --project-name "M81_M85_M134_${base_name}"
#done

# Estimate taxonomy for this .db file
mkdir -p GTDB_Taxonomy
# anvi-run-hmms -c "$db" -T $SLURM_CPUS_PER_TASK
anvi-run-scg-taxonomy -c "$db" -T $SLURM_CPUS_PER_TASK

# Run the anvi-estimate-scg-taxonomy command with the new output file name
anvi-estimate-scg-taxonomy -c "$db" -o "GTDB_Taxonomy/${base_name}_taxonomy.txt"

## Run anvi-run-kegg-kofams for each .db file
#  for db in *.db; do
#   anvi-run-kegg-kofams -c "$db" --log-bitscores -T $SLURM_CPUS_PER_TASK
#  done

## Generate an external genomes file
# anvi-script-gen-genomes-file --input-dir ./ \
#                              --output-file external_genomes.txt

# Estimate metabolism using the external genomes file
anvi-estimate-metabolism -e external_genomes.txt \
                         -O metabolism_estimated \
                         --matrix-format
                         --include-metadata
