#!/usr/bin/env python3
"""Stage 37 - drop reads with 5'/3' soft or hard clipping.

BWA-MEM cannot be configured to never soft-clip, so this script enforces an
end-to-end-only BAM: alignments whose first or last CIGAR operation is S or H
are removed.  Use it before coverage/damage analysis if clipping is not
acceptable.

Usage: python3 scr_37_filter_clipped.py --bam in.bam --out clean.bam
"""
import argparse
import re
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    view = subprocess.Popen(["samtools", "view", "-h", args.bam],
                            stdout=subprocess.PIPE, text=True)
    sam = subprocess.Popen(["samtools", "view", "-b", "-o", args.out, "-"],
                           stdin=subprocess.PIPE)
    kept = dropped = 0
    for line in view.stdout:
        if line.startswith("@"):
            sam.stdin.write(line.encode())
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        if ops and ops[0][1] in "SH" or ops and ops[-1][1] in "SH":
            dropped += 1
            continue
        kept += 1
        sam.stdin.write(line.encode())
    sam.stdin.close()
    view.wait(); sam.wait()
    subprocess.run(["samtools", "index", args.out], check=True)
    print(f"kept {kept}, dropped {dropped} ({100.0*dropped/(kept+dropped) if kept+dropped else 0:.2f}%)")


if __name__ == "__main__":
    main()
