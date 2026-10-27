#!/bin/bash
# List the files that could go into the anvi'o data pack, with sizes, from the project archive.
# Read-only. Run on any machine that sees the archive (the NAS mounted on the laptop, or blythe):
#
#   bash inventory.sh /path/to/CyanoSphere_NAS_20260929 > datapack_inventory.tsv
#
# Output columns: size_bytes, kind, path (relative to the archive root).
set -euo pipefail
ROOT=${1:?give the archive root}
cd "$ROOT"

size() { if stat -f%z "$1" >/dev/null 2>&1; then stat -f%z "$1"; else stat -c%s "$1"; fi; }

printf "size_bytes\tkind\tpath\n"
find . \( -name '*.db' -o -name '*.nwk' -o -name '*.newick' -o -name '*-contigs.fa' \
          -o -name '*.fasta' -o -name '*.fa' -o -name '*_modules.txt' -o -name 'config*.json' \
          -o -name 'NAS_INDEX.tsv' \) -type f 2>/dev/null |
while read -r f; do
    case "$f" in
        *CONTIGS.db|*-contigs.db|*contigs.db) k=contigs-db ;;
        *PROFILE.db)                          k=profile-db ;;
        *AUXILIARY-DATA.db)                   k=auxiliary-db ;;
        *GENOMES.db|*PAN.db)                  k=pan-db ;;
        *.db)                                 k=other-db ;;
        *.nwk|*.newick)                       k=tree ;;
        *_modules.txt)                        k=metabolism ;;
        *.json)                               k=config ;;
        *NAS_INDEX.tsv)                       k=index ;;
        *)                                    k=fasta ;;
    esac
    printf "%s\t%s\t%s\n" "$(size "$f")" "$k" "${f#./}"
done
