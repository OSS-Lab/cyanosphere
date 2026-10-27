#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00

MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$MY_NUM_THREADS

# Load mamba module, initiate shell, and load environment
#module load Mamba/4.14.0-0
source activate gtdbtk-2.3.2

# Create input directory and copy fastas
mkdir gtdbtk
mkdir gtdbtk/input
cp 01_FASTA/*.fasta gtdbtk/input

# Run gtdbtk
gtdbtk classify_wf --cpus 24 --genome_dir gtdbtk/input --extension fasta --out_dir gtdbtk_out --skip_ani_screen
