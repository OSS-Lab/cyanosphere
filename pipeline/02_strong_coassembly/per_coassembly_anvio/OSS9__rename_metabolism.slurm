#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=job_output_%j.out    # Output file
#SBATCH --error=job_error_%j.err      # Error file


## Set the number of threads for OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


## Rename projects without .db
for db in *.db; do
	base_name=$(basename "$db" .db)
    anvi-db-info --self-key project_name --self-value "M82_M86_M135_${base_name}" "${base_name}.db"  --just-do-it
done

#Step2: Generate an external genomes file

anvi-script-gen-genomes-file --input-dir ./ \
                             --output-file external_genomes.txt

# Estimate metabolism using the external genomes file
anvi-estimate-metabolism -e external_genomes.txt \
                         --matrix-format \
                         --include-metadata





