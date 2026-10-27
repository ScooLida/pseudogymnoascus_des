#!/bin/bash
# Paths and parameters. Copy to config/config.sh and edit.

export PROJECT=${PROJECT:-/home/users/ltursunova/Big_data/psy}
export WORK=${WORK:-/mnt/Storage2/ltursunova/destructans_panel}

# Stage 1
export GENES="Pd_18S Pd_ITS Pd_28S Pd_MCM7 Pd_TEF1alpha Pd_RPB2"
export MIN_DEPTH=${MIN_DEPTH:-3}
export BREADTH_MIN=${BREADTH_MIN:-50.0}
export MIN_GENES=${MIN_GENES:-2}
export BREADTH_TSV=${BREADTH_TSV:-$PROJECT/genes_breadth.tsv}
export GENE_BAM_DIR=${GENE_BAM_DIR:-$PROJECT/My_grib_genes}

# Stage 2
export KAIJU_ENV=${KAIJU_ENV:-/mnt/Storage/ltursunova_data/miniconda3/envs/kaiju_env/bin}
export KAIJU_DB=${KAIJU_DB:-/mnt/Storage2/kaiju/missing_human.fmi}
export KAIJU_TAXDUMP=${KAIJU_TAXDUMP:-/mnt/Storage2/kaiju/taxdump_2026_extracted}
export KAIJU_TAXON_LIST=${KAIJU_TAXON_LIST:-$WORK/config/fungi_taxa.txt}

# Stage 3
export SELECTED=${SELECTED:-/mnt/Storage2/ltursunova/pseudogymnoascus/selection/selected_99.9.txt}
export GENOME_DIR=${GENOME_DIR:-/mnt/Storage2/ltursunova/pseudogymnoascus/fna}
export DESTRUCTANS_REF=${DESTRUCTANS_REF:-$PROJECT/grib_combine_named.fasta}
export THREADS=${THREADS:-32}

# Ancient-DNA rule
export DELTAS_MIN=${DELTAS_MIN:-11.0}
export CTOT_MIN=${CTOT_MIN:-1.1}
