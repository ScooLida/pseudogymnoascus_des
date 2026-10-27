#!/usr/bin/env python3
"""Stage 4 - final species report.

Rows = samples, columns = species with three values each:
% coverage, mean depth (MAPQ >= 20) and number of reads.
"""
import argparse
import csv
import glob
import os
import re


def read_stats(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["taxon"]] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-reads", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.assign_dir, "*.species_stats.tsv")))
    per_sample = {}
    species = set()
    for f in files:
        s = os.path.basename(f).replace(".species_stats.tsv", "")
        st = read_stats(f)
        per_sample[s] = st
        species.update(st)

    species = sorted(species)
    cols = ["Sample"]
    for sp in species:
        cols += [f"{sp}|cov_pct", f"{sp}|depth_mq20", f"{sp}|reads"]

    long_path = os.path.join(args.out_dir, "species_report_long.tsv")
    with open(long_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample", "species", "breadth_pct", "depth_mq20", "reads", "breadth_mq20_pct"])
        for s in sorted(per_sample):
            for sp in species:
                r = per_sample[s].get(sp)
                if r and int(r["reads"]) >= args.min_reads:
                    w.writerow([s, sp, r["breadth_pct"], r["mean_depth_mq20"], r["reads"], r["breadth_mq20_pct"]])

    wide_path = os.path.join(args.out_dir, "species_report.tsv")
    with open(wide_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for s in sorted(per_sample):
            row = [s]
            for sp in species:
                r = per_sample[s].get(sp)
                if r and int(r["reads"]) >= args.min_reads:
                    row += [r["breadth_pct"], r["mean_depth_mq20"], r["reads"]]
                else:
                    row += ["0", "0", "0"]
            w.writerow(row)
    print("wrote", long_path)
    print("wrote", wide_path)


if __name__ == "__main__":
    main()
