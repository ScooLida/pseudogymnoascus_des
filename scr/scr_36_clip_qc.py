#!/usr/bin/env python3
"""Stage 36 - 5'/3' end-clipping QC for mapped reads.

Reports, per BAM, how often read ends are soft/hard clipped and by how many
bases.  Run it on both the gene BAMs and the competitive-panel BAMs.

A read is counted as clipped at an end when the first (5' end of the read) or
last (3' end of the read) CIGAR operation is S or H.  Because the 5'/3' ends
are defined relative to the read (not the reference), the first CIGAR op is
always the read's 5' end and the last op is always its 3' end.
"""
import argparse
import csv
import glob
import os
import re
import subprocess


def clip_stats(bam):
    p = subprocess.Popen(["samtools", "view", "-F", "4", bam],
                         stdout=subprocess.PIPE, text=True)
    n = 0
    c5 = c3 = 0
    s5 = s3 = 0
    for line in p.stdout:
        f = line.split("\t")
        if len(f) < 6:
            continue
        cigar = f[5]
        if cigar == "*":
            continue
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if not ops:
            continue
        n += 1
        l5, o5 = int(ops[0][0]), ops[0][1]
        l3, o3 = int(ops[-1][0]), ops[-1][1]
        if o5 in "SH":
            c5 += 1
            s5 += l5
        if o3 in "SH":
            c3 += 1
            s3 += l3
    p.wait()
    return {
        "reads": n,
        "pct_5p_clip": 100.0 * c5 / n if n else 0.0,
        "pct_3p_clip": 100.0 * c3 / n if n else 0.0,
        "mean_clip_5p": s5 / c5 if c5 else 0.0,
        "mean_clip_3p": s3 / c3 if c3 else 0.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bams", nargs="+", help="BAM files or globs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    files = []
    for b in args.bams:
        files.extend(glob.glob(b))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "reads", "pct_5p_clip", "pct_3p_clip",
                    "mean_clip_5p", "mean_clip_3p"])
        for bam in sorted(files):
            s = os.path.basename(bam)
            for suf in (".bam", ".rescaled"):
                s = s.replace(suf, "")
            st = clip_stats(bam)
            w.writerow([s, st["reads"], f"{st['pct_5p_clip']:.3f}",
                        f"{st['pct_3p_clip']:.3f}",
                        f"{st['mean_clip_5p']:.2f}", f"{st['mean_clip_3p']:.2f}"])
            print("done", bam)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
