#!/bin/bash
# Convert per-gene long-format coverage tables into wide TSV matrices.
# Rows represent samples, columns represent genomic positions, and cell values
# contain read depth. Despite the historical name, this script does not create XLS files.

# Configuration
OUT_DIR="./genes_report"
INPUT_FILE_PATTERN="$OUT_DIR/*_fil_cov_*.tsv"
WIDE_SUFFIX="_wide.tsv"

# This file discovery intentionally follows the existing project convention.
files_to_process=$(ls $INPUT_FILE_PATTERN 2>/dev/null | grep -v "$WIDE_SUFFIX")

if [ -z "$files_to_process" ]; then
    echo "Error: no files matching $INPUT_FILE_PATTERN were found in $OUT_DIR."
    exit 1
fi

for input_file in $files_to_process; do
    output_file="${input_file%.tsv}${WIDE_SUFFIX}"
    echo "Converting $(basename "$input_file") -> $(basename "$output_file")"

    awk -F '\t' '
    BEGIN { OFS = "\t" }
    NR > 1 {
        sample = $1
        position = $2
        depth = $3

        if (!(position in position_seen)) {
            positions[++position_count] = position
            position_seen[position] = 1
        }
        if (!(sample in sample_seen)) {
            samples[++sample_count] = sample
            sample_seen[sample] = 1
        }
        data[sample, position] = depth
    }
    END {
        printf "Sample"
        for (i = 1; i <= position_count; i++) printf "\t%s", positions[i]
        printf "\n"

        for (j = 1; j <= sample_count; j++) {
            sample = samples[j]
            printf "%s", sample
            for (i = 1; i <= position_count; i++) {
                position = positions[i]
                value = (data[sample, position] != "") ? data[sample, position] : 0
                printf "\t%s", value
            }
            printf "\n"
        }
    }
    ' "$input_file" > "$output_file"
done

echo "Wide TSV matrices saved to: $OUT_DIR"
