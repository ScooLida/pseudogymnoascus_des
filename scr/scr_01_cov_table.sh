#!/bin/bash
# Calculate per-gene breadth of coverage for every BAM sample.
# A reference position is counted as covered when its read depth is at least
# MIN_DEPTH. Results are written to a tab-separated summary table.

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

printf "Sample\tGene\tBreadthPercentage\n" > "$BREADTH_FILE"
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

        printf "%s\t%s\t%s\n" "$sample_name" "$gene" "$breadth" >> "$BREADTH_FILE"
    done
done

echo "Coverage calculation complete."
echo "Results saved to: $BREADTH_FILE"
