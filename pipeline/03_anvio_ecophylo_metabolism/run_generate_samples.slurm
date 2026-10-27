#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --job-name=generate_samples
#SBATCH --output=generate_samples.out
#SBATCH --error=generate_samples.err

# Change to the CYANOSPHERE directory
cd /home/lifesci/lfsnhf/ANVIO_SPACE/CYANOSPHERE

# Make sure the script is executable
chmod +x generate_samples.sh

# Run the script
./generate_samples.sh
