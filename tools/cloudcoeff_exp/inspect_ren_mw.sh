#!/usr/bin/env bash
# Inspect the extracted Ren snow MW data. Handles the SnowflakeModel_MW_partN/ layout.
# Run in the directory containing SnowflakeModel_MW_part*  (extracts any not yet done)
#   usage: bash inspect_ren_mw.sh [dir]
set -euo pipefail
cd "${1:-.}"

# extract only parts not already extracted
for f in SnowflakeModel_MW_part*.tar.gz; do
  [ -e "$f" ] || continue
  d="${f%.tar.gz}"
  [ -d "$d" ] || { echo ">> extracting $f"; tar xzf "$f"; }
done

echo; echo "=== directory tree under SnowflakeModel_MW_part* (depth 3) ==="
find SnowflakeModel_MW_part* -maxdepth 3 -type d 2>/dev/null | sort

echo; echo "=== ALL isca.dat paths (reveals Scheme/T layout & part->scheme mapping) ==="
find . -name isca.dat | sort

isca=$(find . -name isca.dat | head -1); leaf=$(dirname "$isca")
echo; echo "=== sample leaf: $leaf ==="; ls -l "$leaf"
echo; echo "=== isca.dat: first 3 rows + line count ==="; head -3 "$isca"; wc -l "$isca"

p="$leaf/P11.dat"
if [ -f "$p" ]; then
  echo; echo "=== P11.dat: angles(line1 NF), values(line2 NF), total lines ==="
  awk 'NR==1{print "angles:",NF} NR==2{print "values:",NF} END{print "lines:",NR}' "$p"
fi

echo; echo "=== MW_freq.txt (if present) ==="
mf=$(find . -name 'MW_freq.txt' | head -1)
[ -n "${mf:-}" ] && { echo "$mf"; wc -l "$mf"; head -3 "$mf"; } || echo "(none found)"
