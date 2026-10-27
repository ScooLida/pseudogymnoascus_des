#!/bin/bash
# Stage 3c - competitive mapping through Paleomix.
#
# Generates a Paleomix BAM makefile for the competitive panel and runs it.
# Default aligner is Bowtie2 in end-to-end mode (no 5'/3' soft-clipping); up to
# MAX_ALN alternative alignments per read are retained with -k.  Set
# COMPETITIVE_ALIGNER=BWA-MEM for BWA-MEM (XА-tag alternatives, soft-clips).
#
# Usage: bash scr/scr_33_competitive_map.sh
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

OUT=$WORK/33_map
mkdir -p "$OUT"
MK=$OUT/makefile_competitive.yaml
PANEL_FASTA=${PANEL_FASTA:-$WORK/panel/panel.fasta}
SAMPLES_DIR=${SAMPLES_DIR:-$PROJECT/samples/storage2}
ALIGNER=${COMPETITIVE_ALIGNER:-BWA-MEM}
MAX_ALN=${COMPETITIVE_MAX_ALN:-25}

{
  cat <<YAML
Options:
  Platform: Illumina
  QualityOffset: 33
  AdapterRemoval:
    --mm: 3
    --minlength: 25
    --collapse: yes
    --trimns: yes
    --trimqualities: yes
  Aligners:
    Program: $ALIGNER
YAML
  if [ "$ALIGNER" = "Bowtie2" ]; then
    cat <<YAML
    Bowtie2:
      MinQuality: 0
      FilterUnmappedReads: yes
      --end-to-end:
      "-k": $MAX_ALN
YAML
  else
    cat <<YAML
    BWA-MEM:
      MinQuality: 0
      FilterUnmappedReads: yes
      "-L": ${BWA_MEM_CLIP_PENALTY:-100}
YAML
  fi
  cat <<YAML
  Features:
    PCRDuplicates: filter
    mapDamage: rescale
    Coverage: yes
    Depths: yes
    Summary: yes
Prefixes:
  panel:
    Path: $PANEL_FASTA
Competitive:
YAML
  while read -r s; do
    [ -z "$s" ] && continue
    printf '  %s:\n    %s:\n      Path: %s/%s_[12].fastq*\n' "$s" "$s" "$SAMPLES_DIR" "$s"
  done < "$WORK/01_passing_samples.tsv"
} > "$MK"

paleomix bam run "$MK" --max-threads "$THREADS" --force
echo "Paleomix competitive mapping finished: $OUT"
