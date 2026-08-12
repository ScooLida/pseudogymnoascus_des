#!/bin/bash
# Select samples with sufficient breadth of coverage across several genes from
# the wide sample-by-gene coverage table. For selected samples, extract
# per-position depth tables for every target gene.

# Configuration
OUT_DIR="./genes_report"
INPUT_BREADTH="$OUT_DIR/genes_breadth.tsv"
TARGET_SAMPLES_FILE="$OUT_DIR/target_samples.txt"
BAM_DIR="./My_grib_genes/my_genome"
MIN_PERCENT=60
MIN_SUCCESS_GENES=2
GENES=("Pd_18S" "Pd_ITS" "Pd_28S" "Pd_MCM7" "Pd_TEF1alpha" "Pd_RPB2")

mkdir -p "$OUT_DIR"

if [ ! -f "$INPUT_BREADTH" ]; then
    echo "Error: input table $INPUT_BREADTH was not found."
    exit 1
fi

# Create a persistent list of samples passing the selection criteria.
awk -v threshold="$MIN_PERCENT" -v min_genes="$MIN_SUCCESS_GENES" '
BEGIN { FS = "\t" }
NR == 1 {
    for (column = 2; column <= NF; column++) gene_column[$column] = column
    next
}
{
    sample = $1
    successful_genes = 0
    for (gene in gene_column) {
        column = gene_column[gene]
        if (($column + 0) >= threshold) successful_genes++
    }
    if (successful_genes >= min_genes) print sample
}
' "$INPUT_BREADTH" | sort -u > "$TARGET_SAMPLES_FILE"

mapfile -t target_samples < "$TARGET_SAMPLES_FILE"
target_count=${#target_samples[@]}

if [ "$target_count" -eq 0 ]; then
    echo "No samples passed: at least $MIN_SUCCESS_GENES genes with breadth >= $MIN_PERCENT%."
    exit 0
fi

echo "Selected samples: $target_count"
echo "Sample list saved to: $TARGET_SAMPLES_FILE"
echo "Samples that passed the filter:"
cat "$TARGET_SAMPLES_FILE"

for gene in "${GENES[@]}"; do
    printf "Sample\tPosition\tDepth\n" > "$OUT_DIR/${gene}_temp.tsv"
done

counter=0
for target in "${target_samples[@]}"; do
    counter=$((counter + 1))
    echo "[$counter/$target_count] Processing sample: $target"

    mapfile -t candidate_bams < <(find "$BAM_DIR" -type f -name "${target}*rescaled.bam" -print)
    if [ "${#candidate_bams[@]}" -eq 0 ]; then
        echo "Warning: no BAM file found for $target; skipping."
        continue
    fi
    if [ "${#candidate_bams[@]}" -gt 1 ]; then
        echo "Warning: multiple BAM files found for $target; skipping to avoid ambiguity."
        continue
    fi
    target_bam="${candidate_bams[0]}"

    for gene in "${GENES[@]}"; do
        samtools depth -a -r "$gene" "$target_bam" | \
            awk -v sample="$target" '{ print sample "\t" $2 "\t" $3 }' >> "$OUT_DIR/${gene}_temp.tsv"
    done
done

for gene in "${GENES[@]}"; do
    output_file="$OUT_DIR/${gene}_fil_cov_${target_count}.tsv"
    mv "$OUT_DIR/${gene}_temp.tsv" "$output_file"
    echo "Created: $output_file"
done
