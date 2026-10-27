#!/bin/bash
#SBATCH --job-name=media_euk
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
# What replaced the cyanobacteria in the media cultures that lost them: V_P10_R3 (M_285),
# NV_P10_R2-R4 (M_293-M_295); NV_P10_R1 (M_292) as a culture that kept them.
# Subsample 10 M read pairs, assemble (MEGAHIT), classify contigs >= 3 kb into
# bacteria / eukarya / plastid / mitochondrion (Tiara), and extract 18S and 16S rRNA
# genes (anvi'o rRNA HMMs) for identification. Read-only on the inputs.
#   sbatch --output=$HOME/media_euk-%j.out media_eukaryote.sbatch M_293 M_294 ...   (sample folders in 01.RawData)
set -euo pipefail
A=/springbrook/share/lifesci/lfsnhf/conda_envs/anvio-dev/bin
T=/springbrook/share/lifesci/lfsnhf/conda_envs/tiara_env/bin
R=/springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs/01.RawData
W=/springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_eukaryote
mkdir -p $W && cd $W
for s in "$@"; do
  mkdir -p $s && cd $s
  # head closes the pipe early; zcat's SIGPIPE must not trip pipefail
  for i in 1 2; do
    [ -f r$i.fq ] || { (zcat $R/$s/${s}_SRA_$i.fq.gz || true) | head -n 40000000 > r$i.tmp && mv r$i.tmp r$i.fq; }
  done
  [ -f asm/final.contigs.fa ] || $A/megahit -1 r1.fq -2 r2.fq -o asm -t $SLURM_CPUS_PER_TASK --min-contig-len 1000
  $A/anvi-script-reformat-fasta asm/final.contigs.fa -o contigs.fa -l 3000 --simplify-names --seq-type NT >/dev/null
  $T/tiara -i contigs.fa -o tiara.txt --tf all -t $SLURM_CPUS_PER_TASK --probabilities >/dev/null 2>&1 || true
  $A/anvi-gen-contigs-database -f contigs.fa -o contigs.db -T $SLURM_CPUS_PER_TASK --quiet --force-overwrite
  $A/anvi-run-hmms -c contigs.db -I Ribosomal_RNA_18S,Ribosomal_RNA_16S,Ribosomal_RNA_28S -T $SLURM_CPUS_PER_TASK --quiet --just-do-it
  for h in Ribosomal_RNA_18S Ribosomal_RNA_16S; do
    $A/anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-source $h -o $h.fa --quiet 2>/dev/null || true
  done
  python3 - <<'PY'
import csv, collections
L = {}
for line in open("contigs.fa"):
    if line.startswith(">"): n = line[1:].strip(); L[n] = 0
    else: L[n] += len(line.strip())
cls = collections.Counter()
try:
    for r in csv.DictReader(open("tiara.txt"), delimiter="\t"):
        k = r["class_fst_stage"] + ("/" + r["class_snd_stage"] if r.get("class_snd_stage") not in (None, "", "n/a") else "")
        cls[k] += L.get(r["sequence_id"], 0)
except FileNotFoundError:
    pass
tot = sum(cls.values()) or 1
with open("tiara_summary.txt", "w") as f:
    for k, v in cls.most_common(): f.write(f"{k}\t{v}\t{100*v/tot:.1f}%\n")
PY
  cd $W
done
echo DONE
