#!/bin/bash
# Paths and parameters. Copy to config/config.sh and edit.

export PROJECT=${PROJECT:-/home/users/ltursunova/Big_data/psy}
export WORK=${WORK:-/mnt/Storage2/ltursunova/destructans_panel}

# Stage 1 - gene presence
export GENES="Pd_18S Pd_ITS Pd_28S Pd_MCM7 Pd_TEF1alpha Pd_RPB2"
export MIN_DEPTH=${MIN_DEPTH:-3}
export BREADTH_MIN=${BREADTH_MIN:-50.5}
export MIN_GENES=${MIN_GENES:-2}
export BREADTH_TSV=${BREADTH_TSV:-$PROJECT/genes_breadth.tsv}
export GENE_BAM_DIR=${GENE_BAM_DIR:-$PROJECT/My_grib_genes}

# Stage 2 - MALT species screen (BLASTN + naive LCA)
export MALT_BIN=${MALT_BIN:-/mnt/Storage/ltursunova_data/miniconda3/envs/malt_env/bin}
export MALT_INDEX=${MALT_INDEX:-/mnt/Storage2/ltursunova/malt_db/index}
export MALT_ALIGN_TYPE=${MALT_ALIGN_TYPE:-SemiGlobal}
export MALT_MIN_BITSCORE=${MALT_MIN_BITSCORE:-50}
export MALT_MIN_PIDENT=${MALT_MIN_PIDENT:-0}
export MALT_MAX_ALN=${MALT_MAX_ALN:-25}
export MALT_TOP_PERCENT=${MALT_TOP_PERCENT:-10}
export MALT_SUPPORT=${MALT_SUPPORT:-0.001}

# Stage 3 - competitive mapping
export SELECTED=${SELECTED:-/mnt/Storage2/ltursunova/pseudogymnoascus/selection/selected_95.txt}
export GENOME_DIR=${GENOME_DIR:-/mnt/Storage2/ltursunova/pseudogymnoascus/fna}
export DESTRUCTANS_REF=${DESTRUCTANS_REF:-$PROJECT/grib_combine_named.fasta}
export THREADS=${THREADS:-32}
export PANEL_FASTA=${PANEL_FASTA:-$WORK/panel/panel.fasta}
export SAMPLES_DIR=${SAMPLES_DIR:-$PROJECT/samples/storage2}
export COMPETITIVE_ALIGNER=${COMPETITIVE_ALIGNER:-BWA-MEM}
export COMPETITIVE_MAX_ALN=${COMPETITIVE_MAX_ALN:-25}
export BWA_MEM_CLIP_PENALTY=${BWA_MEM_CLIP_PENALTY:-100}

# Ancient-DNA rule
export DELTAS_MIN=${DELTAS_MIN:-11.0}
export CTOT_MIN=${CTOT_MIN:-1.1}
