#!/bin/bash
#SBATCH --job-name=fig1_lineages
#SBATCH --output=%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
# See README.md. Waits until every sample copied from the NAS passes its Novogene MD5
# check (the copy may still be running when this starts), then maps and profiles.
set -uo pipefail
R=/springbrook/share/lifesci/lfsnhf/metagenomics/runs/fig1_lineages
E=/springbrook/share/lifesci/lfsnhf/conda_envs
export PATH=$E/anvio-dev/bin:$PATH
cd $R
T=$SLURM_CPUS_PER_TASK
SAMPLES="M_79 M_80 M_81 M_82 M_83 M_84 M_85 M_86 M_62 M_70 M_76"
DB=03_CONTIGS/Ribosomal_L19-contigs.db
mkdir -p bam 05_ANVIO_PROFILE sra_tmp

# control end point from SRA, while the NAS copies arrive
if [ ! -s 01.RawData/M_132/M_132_SRA_2.fq.gz ]; then
    mkdir -p 01.RawData/M_132
    $E/media_mapping/bin/prefetch --max-size 60G -O sra_tmp SRR37528228 && \
    $E/media_mapping/bin/fasterq-dump --split-files -e $T -t sra_tmp -O sra_tmp sra_tmp/SRR37528228 && \
    for i in 1 2; do gzip -c sra_tmp/SRR37528228_$i.fastq > 01.RawData/M_132/M_132_SRA_$i.fq.gz; done && \
    rm -rf sra_tmp/SRR37528228* && touch bam/M_132.md5ok
fi

anvi-migrate --migrate-safely $DB >/dev/null 2>&1 || true
[ -f ref.fa ] || anvi-export-contigs -c $DB -o ref.fa
[ -f ref.1.bt2 ] || bowtie2-build -q ref.fa ref

profile() {
    s=$1; d=01.RawData/$s
    [ -f 05_ANVIO_PROFILE/$s/PROFILE.db ] && return 0
    # every lane: some samples (e.g. M_62) were delivered as several read pairs; bowtie2 takes
    # comma-separated lists, sorted identically so that R1 and R2 files pair up
    r1=$(ls $d/*_1.f*q.gz | sort | paste -sd, -); r2=$(ls $d/*_2.f*q.gz | sort | paste -sd, -)
    echo "$s lanes: $r1" >&2
    bowtie2 -p $T -x ref -1 "$r1" -2 "$r2" --no-unal 2> bam/$s.bowtie2.log | samtools view -F 4 -bS - > bam/$s-RAW.bam
    anvi-init-bam bam/$s-RAW.bam -o bam/$s.bam -T $T && rm bam/$s-RAW.bam
    anvi-profile -i bam/$s.bam -c $DB -o 05_ANVIO_PROFILE/$s -S $s --profile-SCVs --min-contig-length 0 -T $T
}
profile M_132

# wait for the NAS copies: a sample is ready when its files pass MD5.txt (checked only once
# the folder has stopped growing for 10 min)
deadline=$((SECONDS + 30 * 3600))
while :; do
    left=""
    for s in $SAMPLES; do
        [ -f bam/$s.md5ok ] && continue
        d=01.RawData/$s
        if [ -f $d/MD5.txt ] && ls $d/*_2.f*q.gz >/dev/null 2>&1; then
            a=$(du -sb $d | cut -f1); sleep 600; b=$(du -sb $d | cut -f1)
            if [ "$a" = "$b" ] && (cd $d && md5sum -c MD5.txt --quiet); then
                touch bam/$s.md5ok; echo "$(date) $s copied and verified"; profile $s; continue
            fi
        fi
        left="$left $s"
    done
    [ -z "$left" ] && break
    [ $SECONDS -gt $deadline ] && { echo "TIMEOUT waiting for:$left"; break; }
    echo "$(date) waiting for:$left"; sleep 600
done

done_profiles=$(ls -d 05_ANVIO_PROFILE/*/PROFILE.db)
anvi-merge $done_profiles -c $DB -o 06_MERGED -S Ribosomal_L19_lineages --overwrite-output-destinations
anvi-export-splits-and-coverages -p 06_MERGED/PROFILE.db -c $DB -o . -O coverage_mean --splits-mode --report-contigs 2>/dev/null \
  || anvi-export-splits-and-coverages -p 06_MERGED/PROFILE.db -c $DB -o . -O coverage_mean --splits-mode
anvi-export-splits-and-coverages -p 06_MERGED/PROFILE.db -c $DB -o . -O coverage_Q2Q3 --splits-mode --use-Q2Q3-coverages

# control: new M_132 against the existing end-point table
$E/media_mapping/bin/python - <<'PY' > control.txt
import pandas as pd, numpy as np
old = pd.read_csv("rpl19_cluster_coverage.csv", index_col=0)["M_132"]
for k in ("mean", "Q2Q3"):
    new = pd.read_csv(f"coverage_{k}-COVs.txt", sep="\t", index_col=0)["M_132"]
    j = pd.concat([old, new], axis=1, keys=["old", "new"]).fillna(0)
    print(k, "r =", round(np.corrcoef(j.old, j.new)[0, 1], 4), "| max abs diff", round((j.old - j.new).abs().max(), 3),
          "| detected old/new", int((j.old > 0).sum()), int((j.new > 0).sum()))
PY
cat control.txt
echo DONE
