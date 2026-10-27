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

11. `scr/scr_11_gene_presence.sh` scans BAM files in `./My_grib_genes/my_genome`,
   uses reads with `MAPQ >= MIN_MAPQ=20`, measures the percentage of positions
   with depth at least `MIN_DEPTH=3`,
   and writes `genes_report/genes_breadth.tsv` as a wide table: one row per
   sample, one column per target gene, and one coverage percentage at each
   sample-gene intersection.
12. `scr/scr_12_target_samples.sh` selects samples with at least
   `MIN_SUCCESS_GENES=2` genes reaching `MIN_PERCENT=60` breadth. It writes
   `genes_report/target_samples.txt` and per-gene depth tables named
   `*_fil_cov_<number>.tsv`. The depth tables use the same `MAPQ >= 20`
   filter as the breadth calculation.
13. `scr/scr_13_wide_tables.sh` converts those long-format tables into wide
   TSV matrices named `*_wide.tsv`.
14. `scr/scr_14_depth_plots.R` builds per-gene depth plots from the selected
   samples.

## Species / competitive-panel workflow

Stages after the coverage filter. Copy `config/config.example.sh` to
`config/config.sh` and edit the paths first.

21. `scr/scr_21_malt_screen.sh` runs a species-level screen with **MALT**
    (BLASTN, semi-global, naive LCA) **only** on the samples that passed the
    gene-presence filter, producing RMA6 files and Tab alignments in
    `21_malt/<sample>/` (open the RMA6 in MEGAN for the LCA taxonomy).
31. `scr/scr_31_build_panel.py` builds the competitive-mapping panel from the
    selected *Pseudogymnoascus* genomes plus close relatives. *P. destructans*
    is deliberately **excluded** from the panel and written separately to
    `destructans_ref.fasta` (used for reference/tree projection). Genome FASTAs
    are taken from `--genome-dir`, relatives from `--relatives`
    (species, FASTA, taxid).
32. `scr/scr_32_index.sh` builds the BWA and faidx indexes.
33. `scr/scr_33_competitive_map.sh SAMPLE ...` maps reads simultaneously to the
    combined panel with `bwa aln`, keeping alternative hits in the XA tag.
34. `scr/scr_34_assign_coverage.py` assigns each read to the best-scoring genome
    (ambiguous reads go to genus or are marked ambiguous), collapses PCR
    duplicates by 5' coordinate + strand, and reports breadth, mean depth and
    reads (including MAPQ >= 20).
35. `scr/scr_35_mapdamage.sh SPECIES SAMPLE ...` runs mapDamage on the reads
    assigned to a target species (fragment lengths and 5'-C>T / 3'-G>A damage).
36. `scr/scr_36_clip_qc.py` checks 5'/3' end clipping (soft/hard clips) in the
    gene BAMs and in the competitive-panel BAMs; run it before downstream
    coverage/damage interpretation.
41. `scr/scr_41_report.py` produces the final report: rows = samples, columns =
    species with % coverage, depth (MAPQ >= 20) and read counts.

### Ancient-DNA authentication

The project rule for "possibly ancient" is **DeltaS >= 11 % and
5'-C>T(pos1) >= 1.1 %** (from `mapDamage`).
