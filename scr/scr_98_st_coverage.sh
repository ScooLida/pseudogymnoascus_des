#!/bin/bash
# Summarize read count, breadth of coverage, and mean depth for each read-group sample
# in one BAM file. The result is written as a comma-separated table.

# Configuration
INPUT_BAM="./My_grib3_18s.my_genome.bam"
OUTPUT_FILE="./merg_st.txt"

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: input BAM file $INPUT_BAM was not found."
    exit 1
fi

printf "Sample,numreads,coverage,meandepth\n" > "$OUTPUT_FILE"

SAMPLES=$(samtools view -H "$INPUT_BAM" | awk '
    /^@RG/ {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^SM:/) {
                sample = $i
                sub(/^SM:/, "", sample)
                sub(/[.-][0-9]$/, "", sample)
                print sample
            }
        }
    }
' | sort -u)

for sample in $SAMPLES; do
    echo "Processing sample group: $sample"

    summary=$(samtools view -h "$INPUT_BAM" | \
        awk -v sample="$sample" '
            /^@/ { print; next }
            $0 ~ ("RG:Z:" sample "([[:space:]]|$)") { print }
        ' | \
        samtools coverage - | \
        awk '
            /^#/ { next }
            NF >= 7 {
                region_length = $3 - $2 + 1
                total_length += region_length
                total_reads += $4
                covered_bases += $5
                depth_sum += $7 * region_length
            }
            END {
                if (total_length == 0 || total_reads == 0) exit 1
                coverage = 100 * covered_bases / total_length
                mean_depth = depth_sum / total_length
                printf "%d\t%.6f\t%.6f", total_reads, coverage, mean_depth
            }
        ') || continue

    IFS=$'\t' read -r numreads coverage meandepth <<< "$summary"
    printf "%s,%s,%s,%s\n" "$sample" "$numreads" "$coverage" "$meandepth" >> "$OUTPUT_FILE"
done

echo "Coverage summary saved to: $OUTPUT_FILE"
