#!/usr/bin/env bash
# Download the Ren et al. (2022) "Snow single-scattering database" (II-TM/IGOM)
# from the Texas Data Repository. Run on a host that can reach dataverse.tdl.org
# (e.g. a NOAA HPC login node). Resumable: re-run to continue interrupted files.
#   usage: bash download_ren_snow.sh [output_dir]
set -euo pipefail
DOI="doi:10.18738/T8/LGJ9SA"
BASE="https://dataverse.tdl.org"
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
OUT="${1:-ren_snow}"

mkdir -p "$OUT"; cd "$OUT"
echo ">> fetching file listing"
curl -s -G "$BASE/api/datasets/:persistentId/" --data-urlencode "persistentId=$DOI" -A "$UA" -o ds.json

python3 - <<'PY' > files.txt
import json
for f in json.load(open('ds.json'))['data']['latestVersion']['files']:
    d = f['dataFile']
    print(d['id'], d.get('filename', f.get('label', '')))
PY

echo ">> files:"; awk '{printf "   %s  %s\n",$1,$2}' files.txt
while read -r id name; do
  echo ">> downloading $name (id $id)"
  wget -c --user-agent="$UA" -O "$name" "$BASE/api/access/datafile/$id"
done < files.txt

echo ">> done. verify sizes (ls -l), then unpack: gunzip -k *.gz"
