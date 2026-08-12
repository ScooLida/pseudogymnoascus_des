# Detecting Fungal Contamination in Museum Specimens
**#in progress**

This repository contains scripts for screening museum specimens for fungal DNA,
with a focus on the genus *Pseudogymnoascus*. The project is under active
development and is intended for short-read Illumina sequencing data in FASTQ or
compressed FASTQ.GZ format.

## Coverage Workflow

This workflow evaluates coverage of six target genes:

`Pd_18S`, `Pd_ITS`, `Pd_28S`, `Pd_MCM7`, `Pd_TEF1alpha`, and `Pd_RPB2`.

Run the scripts from the repository root in this order:

1. `scr/scr_01_cov_table.sh` scans BAM files in `./My_grib_genes/my_genome`,
   measures the percentage of positions with depth at least `MIN_DEPTH=3`,
   and writes `genes_report/genes_breadth.tsv` as a wide table: one row per
   sample, one column per target gene, and one coverage percentage at each
   sample-gene intersection.
2. `scr/scr_02_target_samples.sh` selects samples with at least
   `MIN_SUCCESS_GENES=2` genes reaching `MIN_PERCENT=60` breadth. It writes
   `genes_report/target_samples.txt` and per-gene depth tables named
   `*_fil_cov_<number>.tsv`.
3. `scr/scr_03_cov_table_xls.sh` converts those long-format tables into wide
   TSV matrices named `*_wide.tsv`. Despite the historical `xls` name, no XLS
   files are produced.
