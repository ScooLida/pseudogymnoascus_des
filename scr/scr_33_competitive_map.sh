#!/bin/bash
# Stage 3c - competitive mapping with BWA aln.
# All selected reads are aligned simultaneously to the combined panel; up to
# MAX_ALT alternative hits are kept in the XA tag (-n), so that ambiguous reads
# can later be assigned to a group/genus instead of a species.
#
# Usage: bash scr/scr_33_competitive_map.sh SAMPLE [SAMPLE ...]
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

MAX_ALT=${MAX_ALT:-100}
P=$WORK/panel
OUT=$WORK/33_map
mkdir -p "$OUT"

for s in "$@"; do
  r1=$(ls "$PROJECT"/samples/storage2/${s}_1.fastq* 2>/dev/null | head -1)
  r2=$(ls "$PROJECT"/samples/storage2/${s}_2.fastq* 2>/dev/null | head -1)
  [ -z "$r1" ] && { echo "skip $s: no R1"; continue; }
  d=$OUT/$s; mkdir -p "$d"
  bwa aln -t "$THREADS" -n 0.01 -o 2 -e 0.01 -l 1024 "$P/panel.fasta" "$r1" > "$d/r1.sai"
  if [ -n "$r2" ]; then
    bwa aln -t "$THREADS" -n 0.01 -o 2 -e 0.01 -l 1024 "$P/panel.fasta" "$r2" > "$d/r2.sai"
    bwa sampe -n "$MAX_ALT" -N "$MAX_ALT" -P "$P/panel.fasta" \
      "$d/r1.sai" "$d/r2.sai" "$r1" "$r2" > "$d/$s.sam"
  else
    bwa samse -n "$MAX_ALT" "$P/panel.fasta" "$d/r1.sai" "$r1" > "$d/$s.sam"
  fi
  samtools view -bhS "$d/$s.sam" | samtools sort -@ "$THREADS" -o "$d/$s.bam" -
  samtools index "$d/$s.bam"
  rm -f "$d"/*.sai "$d/$s.sam"
  echo "mapped $s"
done
