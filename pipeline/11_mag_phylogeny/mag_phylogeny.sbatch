#!/bin/bash
#SBATCH --job-name=mag_phylogeny
#SBATCH --output=%x-%j.out
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# Phylogenomic tree of the 88 dereplicated enrichment MAGs (SI Fig S8): Bacteria_71
# single-copy genes found in >= 80% of MAGs, concatenated (anvi'o, muscle), trimmed
# (trimal -gt 0.5), maximum likelihood (IQ-TREE, ModelFinder, 1000 ultrafast bootstraps).
set -euo pipefail
E=/springbrook/share/lifesci/lfsnhf/conda_envs/anvio-dev/bin
export PATH=$E:$PATH
G=/springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs/genomes
W=/springbrook/share/lifesci/lfsnhf/metagenomics/runs/mag_phylogeny
mkdir -p $W/dbs && cd $W
ls $G/*.fasta | grep -v -E "/(VIT|NOVIT|NONANO)_" > mags.txt
echo "$(wc -l < mags.txt) MAGs"
one() {
    f=$1; n=$(basename $f .fasta | tr '.-' '__')
    [ -f dbs/$n.hmms.done ] && return 0
    anvi-script-reformat-fasta $f -o dbs/$n.fa --simplify-names --seq-type NT >/dev/null 2>&1
    anvi-gen-contigs-database -f dbs/$n.fa -o dbs/$n.db -n $n -T 2 --quiet --force-overwrite
    anvi-run-hmms -c dbs/$n.db -I Bacteria_71 -T 2 --quiet && touch dbs/$n.hmms.done
}
export -f one
cat mags.txt | xargs -P 16 -I{} bash -c 'one {}'
echo -e "name\tcontigs_db_path" > external_genomes.txt
for d in dbs/*.db; do echo -e "$(basename $d .db)\t$W/$d"; done >> external_genomes.txt
n=$(($(wc -l < external_genomes.txt) - 1)); min=$((n * 8 / 10))
anvi-get-sequences-for-hmm-hits -e external_genomes.txt --hmm-source Bacteria_71 --return-best-hit \
    --get-aa-sequences --concatenate-genes --min-num-bins-gene-occurs $min --align-with muscle \
    -o concatenated_Bacteria_71.faa
trimal -in concatenated_Bacteria_71.faa -out concatenated_Bacteria_71.trimmed.faa -gt 0.5
iqtree -s concatenated_Bacteria_71.trimmed.faa -m MFP -B 1000 -T AUTO --threads-max $SLURM_CPUS_PER_TASK \
    --prefix mags_bac71 -redo
echo DONE
