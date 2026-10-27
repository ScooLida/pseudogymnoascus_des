#!/bin/bash
# Stage 3e - mapDamage authentication of reads assigned to a target species.
# Usage: bash scr/scr_35_mapdamage.sh TARGET_SPECIES SAMPLE [SAMPLE ...]
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

TARGET=$1; shift
P=$WORK/panel
OUT=$WORK/35_mapdamage
mkdir -p "$OUT"
export PATH=/mnt/Storage/ltursunova_data/miniconda3/envs/busco_env/bin:$PATH

samtools faidx "$P/panel.fasta"
awk -F"\t" -v t="$TARGET" '$2==t {print $1}' "$P/taxonomy.tsv" > "$OUT/genomes.txt"
awk -F"\t" 'NR==FNR{keep[$1]=1; next} {split($1,a,"|"); if (a[1] in keep) print $1"\t0\t"$2}' \
  "$OUT/genomes.txt" "$P/panel.fasta.fai" > "$OUT/regions.bed"
samtools faidx -r "$OUT/regions.bed" "$P/panel.fasta" > "$OUT/${TARGET// /_}.fasta"

for s in "$@"; do
  d=$OUT/$s; mkdir -p "$d"
  awk -F"\t" -v t="$TARGET" '$2==t {print $1}' "$WORK/34_assign/$s.assign.tsv" > "$d/qnames.txt"
  samtools view -bh -N "$d/qnames.txt" -L "$OUT/regions.bed" "$WORK/33_map/$s/$s.bam" | \
    samtools sort -@ "$THREADS" -o "$d/$s.target.bam" -
  samtools index "$d/$s.target.bam"
  mapDamage -i "$d/$s.target.bam" -r "$OUT/${TARGET// /_}.fasta" -d "$d/mapdamage" --no-progress
done
