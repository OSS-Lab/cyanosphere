#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00

# this value should be less or equal to --cpus-per-task
# larger --cpus-per-task, the more memory is allocated (useful).
# however, many threads can cause slow downs, so less threads sometimes is desireble
export OMP_NUM_THREADS=24

srun ./mergeReads.sh

