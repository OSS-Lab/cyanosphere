#!/bin/bash
#SBATCH --job-name=media_to_MAGs
#SBATCH --output=%x-%j.out
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
# Dereplicate the reference (88 enrichment MAGs + media MAGs) at 95% ANI, then map all
# 35 media samples to the representatives and compute per-genome coverage.
set -euo pipefail
module load Miniconda3/24.7.1-0
ENV=/springbrook/share/lifesci/lfsnhf/conda_envs/media_mapping
[ -d $ENV ] || conda create -y -p $ENV -c conda-forge -c bioconda drep bowtie2 samtools coverm fastani
source activate $ENV
cd /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs
T=$SLURM_CPUS_PER_TASK

# completeness/contamination are already filtered; skip CheckM inside dRep
[ -d drep/dereplicated_genomes ] || \
    dRep dereplicate drep -g genomes/*.fasta -sa 0.95 --ignoreGenomeQuality -p $T

mkdir -p bam
cat drep/dereplicated_genomes/*.fasta > reference.fa
[ -f reference.fa.1.bt2 ] || [ -f reference.fa.1.bt2l ] || bowtie2-build --threads $T reference.fa reference.fa
# 01.RawData/<sample>/<sample>_*_1.fq.gz and _2 (Novogene delivery layout, copied from the NAS)
for d in 01.RawData/M_*; do
    s=$(basename $d)
    # verify the copy against Novogene's checksums before mapping
    if [ ! -f bam/$s.md5ok ]; then
        (cd $d && md5sum -c MD5.txt --ignore-missing --quiet) && touch bam/$s.md5ok \
            || { echo "CHECKSUM FAILED: $s" >&2; continue; }
    fi
    r1=$(ls $d/*_1.f*q.gz | head -1); r2=$(ls $d/*_2.f*q.gz | head -1)
    [ -s bam/$s.bam ] && continue
    bowtie2 -p $T -x reference.fa -1 "$r1" -2 "$r2" --no-unal 2> bam/$s.bowtie2.log \
        | samtools sort -@ 8 -o bam/$s.bam -
done

coverm genome -b bam/*.bam -d drep/dereplicated_genomes -x fasta -t $T \
    -m mean --min-read-percent-identity 95 --min-covered-fraction 10 > coverm_mean.tsv
coverm genome -b bam/*.bam -d drep/dereplicated_genomes -x fasta -t $T \
    -m relative_abundance --min-read-percent-identity 95 > coverm_relative_abundance.tsv
echo DONE
