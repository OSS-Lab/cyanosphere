#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00

# this value should be less or equal to --cpus-per-task 
# larger --cpus-per-task, the more memory is allocated (useful).
# however, many threads can cause slow downs, so less threads sometimes is desireble
export OMP_NUM_THREADS=48

# specify the singularity container to launch
container=/home/shared/STRONG/containers/STRONG-b25b173.sif # <- keep this the same

# Set output directory:following directs to current directory
outputdir=$(pwd)
#set to what you want like this
#outputdir=/home/lifesci/lfsmbn/tempData

# run the container 
singularity run ${container} "/STRONG/bin/STRONG ${outputdir} --unlock --threads 48 --verbose"
#singularity run ${container} "/STRONG/bin/STRONG ${outputdir} --threads 48 --verbose"
