#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=31418
#SBATCH --time=48:00:00
#SBATCH --partition=hmem

# Load singularity module (if needed)
module load singularity 2>/dev/null || true

export OMP_NUM_THREADS=48
container=/home/shared/STRONG/containers/STRONG-b25b173.sif
outputdir=$(pwd)

echo "=== Checking remaining jobs ==="
singularity exec ${container} snakemake --directory ${outputdir} \
    --cores 1 \
    --config LOCAL_DIR=/STRONG \
    --configfile=${outputdir}/config.yaml \
    --snakefile /STRONG/SnakeNest/SCogSubGraph.snake \
    --dry-run --verbose

echo "=== Unlocking any stuck processes ==="
singularity exec ${container} snakemake --directory ${outputdir} \
    --unlock \
    --configfile=${outputdir}/config.yaml \
    --snakefile /STRONG/SnakeNest/SCogSubGraph.snake

echo "=== Resuming pipeline ==="
singularity run ${container} "/STRONG/bin/STRONG ${outputdir} --threads 48 --verbose --keep-going"
