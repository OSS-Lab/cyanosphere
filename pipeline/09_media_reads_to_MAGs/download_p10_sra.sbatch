#!/bin/bash
#SBATCH --job-name=p10_sra
#SBATCH --output=%x-%j.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
# The 14 P10 media metagenomes from SRA (BioProject PRJNA1432497), into 01.RawData/<sample>/
# with the same naming as the Novogene delivery; MD5.txt is written from the downloaded files
# (SRA verifies its own transfers).
set -euo pipefail
module load Miniconda3/24.7.1-0
source activate /springbrook/share/lifesci/lfsnhf/conda_envs/media_mapping
cd /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs
mkdir -p sra_tmp
while IFS=$'\t' read -r s run name; do
    d=01.RawData/$s
    [ -s $d/${s}_SRA_1.fq.gz ] && [ -s $d/${s}_SRA_2.fq.gz ] && continue
    mkdir -p $d
    prefetch --max-size 60G -O sra_tmp $run
    fasterq-dump --split-files -e $SLURM_CPUS_PER_TASK -t sra_tmp -O sra_tmp sra_tmp/$run
    pigz -p $SLURM_CPUS_PER_TASK -c sra_tmp/${run}_1.fastq > $d/${s}_SRA_1.fq.gz
    pigz -p $SLURM_CPUS_PER_TASK -c sra_tmp/${run}_2.fastq > $d/${s}_SRA_2.fq.gz
    (cd $d && md5sum ${s}_SRA_1.fq.gz ${s}_SRA_2.fq.gz > MD5.txt)
    rm -rf sra_tmp/$run sra_tmp/${run}_*.fastq
    echo "$s ($run, $name) done"
done < p10_sra_runs.tsv
echo DONE
