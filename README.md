# Detecting Fungal Contamination in Museum Specimens

This repository contains scripts for screening museum specimens for fungal DNA,
with a focus on the genus *Pseudogymnoascus*. The project is under active
development and is intended for short-read Illumina sequencing data in FASTQ or
compressed FASTQ.GZ format.

## Coverage Workflow

The main workflow evaluates coverage of six target genes:

`Pd_18S`, `Pd_ITS`, `Pd_28S`, `Pd_MCM7`, `Pd_TEF1alpha`, and `Pd_RPB2`.

Run the scripts from the repository root in this order:

1. `scr/scr_01_cov_table.sh` scans BAM files in `./My_grib_genes/my_genome`,
   measures the percentage of positions with depth at least `MIN_DEPTH=3`,
   and writes `genes_report/genes_breadth.tsv`.
2. `scr/scr_02_target_samples.sh` selects samples with at least
   `MIN_SUCCESS_GENES=2` genes reaching `MIN_PERCENT=60` breadth. It writes
   `genes_report/target_samples.txt` and per-gene depth tables named
   `*_fil_cov_<number>.tsv`.
3. `scr/scr_03_cov_table_xls.sh` converts those long-format tables into wide
   TSV matrices named `*_wide.tsv`. Despite the historical `xls` name, no XLS
   files are produced.

The independent `scr/scr_98_st_coverage.sh` script summarizes read count,
coverage, and mean depth for read-group samples in the BAM file configured by
`INPUT_BAM`. Its output is `merg_st.txt`.

## Configuration

The main input and output paths, target genes, and filtering thresholds are
declared near the beginning of each script. Adjust them before running the
workflow if the project directory layout differs.

The scripts use `samtools`, `awk`, and standard Unix tools. BAM files must be
indexed or indexable by `samtools index`.

## Notes

`scr/scr_99_notes.sh` contains shell notes and examples. It is not part of the
coverage workflow.
