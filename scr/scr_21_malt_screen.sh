#!/bin/bash
# Stage 2 - species-level screen with MALT (BLASTN + naive LCA).
# Runs only on samples that passed the gene-presence filter.
# Usage: bash scr/scr_21_malt_screen.sh [sample ...]
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

OUT=$WORK/21_malt
mkdir -p "$OUT"

SAMPLES=("$@")
[ ${#SAMPLES[@]} -eq 0 ] && mapfile -t SAMPLES < "$WORK/01_passing_samples.tsv"

for s in "${SAMPLES[@]}"; do
  r1=$(ls "$PROJECT"/samples/storage2/${s}_1.fastq* 2>/dev/null | head -1)
  r2=$(ls "$PROJECT"/samples/storage2/${s}_2.fastq* 2>/dev/null | head -1)
  [ -z "$r1" ] && { echo "skip $s: no R1"; continue; }
  d=$OUT/$s; mkdir -p "$d"
  input=("$r1")
  [ -n "$r2" ] && input+=("$r2")
  "$MALT_BIN/malt-run" -m BlastN -at "$MALT_ALIGN_TYPE" \
    -i "${input[@]}" -d "$MALT_INDEX" -t "$THREADS" \
    -b "$MALT_MIN_BITSCORE" -id "$MALT_MIN_PIDENT" -mq "$MALT_MAX_ALN" \
    -top "$MALT_TOP_PERCENT" -supp "$MALT_SUPPORT" \
    -o "$d" -a "$d" -f Tab
  echo "MALT done: $s"
done
