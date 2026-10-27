#!/bin/bash
# Copy the reference genomes to blythe: 88 dereplicated enrichment MAGs + media MAGs
# that pass completion >= 50% and redundancy <= 10% (data/media_mag_quality.csv).
# Run from the folder that holds the project archive ($NAS). Streams with rsync, nothing extra written locally.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)          # this repository, wherever it is cloned
NAS=CyanoSphere_NAS_20260929
DEST=blythe:/springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs
MEDIA=$NAS/04_media_experiment/coassemblies_V_NV_NVNN
[ -d "$MEDIA" ] || MEDIA=$NAS/MEDIA

ssh blythe "mkdir -p /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs/genomes"
rsync -av $NAS/02_MAGs_dereplicated/MAG_contigs_dbs_and_metabolism/MAG_fasta/*.fasta "$DEST/genomes/"

# media MAGs: <bin>-contigs.fa from each co-assembly's anvi-summarize output
tail -n +2 "$REPO/data/media_mag_quality.csv" | cut -d, -f1 | while read -r b; do
    cond=${b%%_*}
    echo "$MEDIA/$cond/SUMMARY/bin_by_bin/$b/$b-contigs.fa"
done > /tmp/media_mag_paths.txt
rsync -av --files-from=<(sed "s|^$MEDIA/||" /tmp/media_mag_paths.txt) "$MEDIA/" "$DEST/media_bins/"
ssh blythe "cd /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs && \
    for f in \$(find media_bins -name '*-contigs.fa'); do b=\$(basename \$f -contigs.fa); cp \$f genomes/\$b.fasta; done && \
    ls genomes | wc -l"
rsync -av "$REPO"/pipeline/09_media_reads_to_MAGs/{map_media.sbatch,media_samples.tsv} "$DEST/"
