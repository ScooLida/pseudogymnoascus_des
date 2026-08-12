#!/bin/bash
# Calculate breadth of coverage for every BAM sample and target gene.
# Each row contains one sample, each column contains one gene, and each cell
# contains the percentage of positions with read depth at least MIN_DEPTH.

# Configuration
OUT_DIR="./genes_report"
BAM_DIR="./My_grib_genes/my_genome"
BAM_GLOB="*/*rescaled.bam"
BREADTH_FILE="$OUT_DIR/genes_breadth.tsv"
MIN_DEPTH=3
GENES=("Pd_18S" "Pd_ITS" "Pd_28S" "Pd_MCM7" "Pd_TEF1alpha" "Pd_RPB2")

shopt -s nullglob
mkdir -p "$OUT_DIR"
BAM_FILES=("$BAM_DIR"/*/*rescaled.bam)

if [ "${#BAM_FILES[@]}" -eq 0 ]; then
    echo "Error: no BAM files matching $BAM_GLOB were found in $BAM_DIR."
    exit 1
fi

{
    printf "Sample"
    printf "\t%s" "${GENES[@]}"
    printf "\n"
} > "$BREADTH_FILE"
total_files=${#BAM_FILES[@]}
echo "Samples found: $total_files"
echo "Coverage threshold: depth >= $MIN_DEPTH"

counter=0
for bam in "${BAM_FILES[@]}"; do
    counter=$((counter + 1))

    sample_name=$(basename "$bam")
    sample_name="${sample_name%.rescaled.bam}"
    sample_name="${sample_name%.my_genome.bam}"
    sample_name="${sample_name%.bam}"

    echo "[$counter/$total_files] Processing: $sample_name"

    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ] && [ ! -f "${bam}.csi" ]; then
        samtools index "$bam"
    fi

    breadth_values=()
    for gene in "${GENES[@]}"; do
        breadth=$(samtools depth -a -r "$gene" "$bam" | awk -v min_depth="$MIN_DEPTH" '
            {
                total++
                if ($3 >= min_depth) covered++
            }
            END {
                if (total == 0) print 0
                else print (covered / total) * 100
            }
        ')
        breadth_values+=("$breadth")
    done

    printf "%s" "$sample_name" >> "$BREADTH_FILE"
    printf "\t%s" "${breadth_values[@]}" >> "$BREADTH_FILE"
    printf "\n" >> "$BREADTH_FILE"
done

echo "Coverage calculation complete."
echo "Results saved to: $BREADTH_FILE"
