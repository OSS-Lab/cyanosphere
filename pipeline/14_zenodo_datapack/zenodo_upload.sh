#!/bin/bash
# Upload the anvi'o data pack to a new Zenodo record, as a draft with a reserved DOI.
# Nothing is published: check the draft on zenodo.org and press Publish there.
#
#   export ZENODO_TOKEN=...        # zenodo.org > Account > Applications > Personal access token
#                                  # scopes: deposit:write and deposit:actions
#   bash zenodo_upload.sh metadata.json CYANOSPHERE_DATAPACK_*.tar.gz MD5SUMS.txt
#
# metadata.json holds the record's title, authors, description and licence (see metadata.json
# in this folder). Large files go to the record's bucket with a streaming PUT, so they never
# need to fit in memory; use tmux, as an upload of tens of GB takes hours.
set -euo pipefail
: "${ZENODO_TOKEN:?export ZENODO_TOKEN first}"
API=${ZENODO_API:-https://zenodo.org/api}     # https://sandbox.zenodo.org/api for a dry run
META=$1; shift

dep=$(curl -sS -f -X POST "$API/deposit/depositions" -H "Authorization: Bearer $ZENODO_TOKEN" \
           -H "Content-Type: application/json" -d '{}')
id=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$dep")
bucket=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["links"]["bucket"])' <<<"$dep")
echo "draft deposition $id"

for f in "$@"; do
    echo "uploading $f ($(du -h "$f" | cut -f1))"
    curl -sS -f --retry 3 -o /dev/null -X PUT "$bucket/$(basename "$f")" \
         -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file "$f"
done

out=$(curl -sS -f -X PUT "$API/deposit/depositions/$id" -H "Authorization: Bearer $ZENODO_TOKEN" \
           -H "Content-Type: application/json" --data @"$META")
python3 - "$out" <<'EOF'
import json, sys
m = json.loads(sys.argv[1])["metadata"]
print("reserved DOI:", m.get("prereserve_doi", {}).get("doi"))
EOF
echo "check the draft at https://zenodo.org/deposit/$id and press Publish when it is right"
