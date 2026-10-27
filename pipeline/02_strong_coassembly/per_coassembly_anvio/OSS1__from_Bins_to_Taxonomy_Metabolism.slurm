#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1          # Only one task per node
#SBATCH --cpus-per-task=24           # Number of CPU cores per task
#SBATCH --mem-per-cpu=3700           # Memory per CPU core
#SBATCH --time=48:00:00
#SBATCH --array=0-52                 # Array job, assuming you have 53 .db files
#SBATCH --output=job_output_%A_%a.out # Output file, %A is the job ID, %a is the array index
#SBATCH --error=job_error_%A_%a.err   # Error file

# Activate your Conda environment
# conda activate anvio-dev

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Get the list of .db files
db_files=(*.db)

# Determine which .db file this task should process based on the array index
db=${db_files[$SLURM_ARRAY_TASK_ID]}

# # Estimate taxonomy for this .db file
# mkdir -p ../GTDB_Taxonomy
# anvi-run-hmms -c "$db" -T $SLURM_CPUS_PER_TASK
# anvi-run-scg-taxonomy -c "$db" -T $SLURM_CPUS_PER_TASK

# # Extract the base name (without the .db extension)
# base_name=$(basename "$db" .db)

# # Run the anvi-estimate-scg-taxonomy command with the new output file name
# anvi-estimate-scg-taxonomy -c "$db" -o "../GTDB_Taxonomy/${base_name}_taxonomy.txt" -T $SLURM_CPUS_PER_TASK

# Run anvi-run-kegg-kofams for this .db file
anvi-run-kegg-kofams -c "$db" --just-do-it --log-bitscores -T $SLURM_CPUS_PER_TASK

# Capture the job ID for the current task for dependencies later
current_job_id=$SLURM_ARRAY_JOB_ID

# After all annotation tasks are complete, generate the external genomes file and estimate metabolism
if [ $SLURM_ARRAY_TASK_ID -eq 52 ]; then
    # Generate an external genomes file
    external_genomes_job_id=$(sbatch --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK --wrap="anvi-script-gen-genomes-file --input-dir ./ --output-file external_genomes.txt")

    # Estimate metabolism using the external genomes file
    sbatch --dependency=afterok:${external_genomes_job_id##* } --wrap="anvi-estimate-metabolism -e external_genomes.txt -O metabolism_estimated --matrix-format"
fi

