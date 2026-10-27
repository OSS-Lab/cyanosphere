#!/bin/bash
# Resume a STRONG co-assembly with DESMAN strain resolution switched off
# (config_no_desman.yaml); used where the full run exceeded the memory or time limit.
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=31418
#SBATCH --time=24:00:00
#SBATCH --partition=hmem

module load singularity 2>/dev/null || true

export OMP_NUM_THREADS=48
container=/home/shared/STRONG/containers/STRONG-b25b173.sif
outputdir=$(pwd)

echo "=== Running STRONG without DESMAN strain resolution ==="
echo "This should complete much faster and use less memory"

singularity run ${container} "/STRONG/bin/STRONG ${outputdir} --threads 48 --verbose --config config_no_desman.yaml"
