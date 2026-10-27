#!/bin/bash
# Stage 2 - species-level screen (Kaiju) of samples that passed stage 1.
# Usage: bash scr/scr_21_species_screen.sh [sample ...]
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

BIN=$KAIJU_ENV
OUT=$WORK/21_species
mkdir -p "$OUT"

SAMPLES=("$@")
[ ${#SAMPLES[@]} -eq 0 ] && mapfile -t SAMPLES < "$WORK/01_passing_samples.tsv"

for s in "${SAMPLES[@]}"; do
  r1=$(ls "$PROJECT"/samples/storage2/${s}_1.fastq* 2>/dev/null | head -1)
  r2=$(ls "$PROJECT"/samples/storage2/${s}_2.fastq* 2>/dev/null | head -1)
  [ -z "$r1" ] && { echo "skip $s: no R1"; continue; }
  if [ -n "$r2" ]; then
    "$BIN/kaiju-multi" -t "$THREADS" -f "$KAIJU_DB" \
      -i "$r1" -j "$r2" -o "$OUT/$s.kaiju.out" -v
  else
    "$BIN/kaiju-multi" -t "$THREADS" -f "$KAIJU_DB" \
      -i "$r1" -o "$OUT/$s.kaiju.out" -v
  fi
  "$BIN/kaiju2table" -t "$KAIJU_TAXDUMP/nodes.dmp" \
    -n "$KAIJU_TAXDUMP/names.dmp" -r species \
    -o "$OUT/$s.species.tsv" "$OUT/$s.kaiju.out"
done
