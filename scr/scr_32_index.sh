#!/bin/bash
# Stage 3b - indexes for the panel and destructans reference.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
P=$WORK/panel
bwa index "$P/panel.fasta"
samtools faidx "$P/panel.fasta"
if [ -f "$P/destructans_ref.fasta" ]; then
  bwa index "$P/destructans_ref.fasta"
  samtools faidx "$P/destructans_ref.fasta"
fi
echo "indexed"
