#!/usr/bin/env python3
"""Stage 3d - read assignment, PCR de-duplication and per-species coverage.

Uses the competitive BWA alignments (primary hit + XA alternatives).  A read is
assigned to a species when its best alignment score is unique (within
--delta-as); otherwise it is assigned to the genus (if all candidates share a
genus) or marked ambiguous.

PCR duplicates are collapsed by the 5' mapping coordinate + strand of the
assigned genome, so every original fragment is counted once.
"""
import argparse
import collections
import csv
import os
import subprocess
import sys


def parse_taxonomy(path):
    tax = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = row["genome"]
            sp = row["species"]
            genus = sp.split()[0] if sp else g
            tax[g] = (sp, genus)
    return tax


def genome_of(ref):
    return ref.split("|")[0]


def parse_cigar(cigar, pos):
    """Return (ref_span, list_of_covered_ref_positions)."""
    import re
    span = 0
    cov = []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            cov.extend(range(pos + span, pos + span + n))
            span += n
        elif op in "DN":
            span += n
        elif op in "I":
            pass
    return span, cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--delta-as", type=int, default=10)
    ap.add_argument("--mapq-depth", type=int, default=20)
    args = ap.parse_args()

    tax = parse_taxonomy(args.taxonomy)

    # genome lengths from BAM header
    hdr = subprocess.run(["samtools", "view", "-H", args.bam],
                         check=True, capture_output=True, text=True).stdout
    lengths = collections.Counter()
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:])
            lengths[genome_of(f["SN"])] += int(f["LN"])

    per_species = collections.defaultdict(lambda: {
        "reads": 0, "cov": set(), "aligned_bases": 0, "mq_bases": 0, "mq_cov": set(),
    })
    seen = set()

    p = subprocess.Popen(["samtools", "view", "-F", "4", args.bam],
                         stdout=subprocess.PIPE, text=True)
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["read", "assigned", "level", "best_ref", "best_as", "n_candidates"])
        for line in p.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            qname, flag, rname, pos, mapq, cigar, seq = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[6], f[9]
            if rname == "*":
                continue
            tags = {}
            for t in f[11:]:
                if t.startswith("AS:i:"):
                    tags["AS"] = int(t[5:])
                elif t.startswith("XA:Z:"):
                    tags["XA"] = t[5:]

            cands = [(rname, tags.get("AS", 0))]
            if "XA" in tags:
                for alt in tags["XA"].split(";"):
                    if not alt:
                        continue
                    a = alt.split(",")
                    if len(a) >= 5:
                        try:
                            cands.append((a[0], int(a[4])))
                        except ValueError:
                            pass
            if not cands:
                continue
            best = max(x[1] for x in cands)
            winners = {genome_of(r) for r, s in cands if s >= best - args.delta_as}
            genomes = [g for g in winners if g in tax]
            if not genomes:
                continue
            species = {tax[g][0] for g in genomes}
            genera = {tax[g][1] for g in genomes}
            if len(species) == 1:
                assigned, level = species.pop(), "species"
            elif len(genera) == 1:
                assigned, level = genera.pop(), "genus"
            else:
                assigned, level = "ambiguous", "ambiguous"

            # de-duplicate by 5' coordinate + strand on the best reference
            _, cov = parse_cigar(cigar, pos)
            if not cov:
                continue
            strand = "-" if flag & 16 else "+"
            key = (assigned, rname, min(cov) if strand == "+" else max(cov), strand)
            if key in seen:
                w.writerow([qname, assigned, level, rname, best, len(cands)])
                continue
            seen.add(key)

            st = per_species[assigned]
            st["reads"] += 1
            st["cov"].update(cov)
            st["aligned_bases"] += len(cov)
            if mapq >= args.mapq_depth:
                st["mq_bases"] += len(cov)
                st["mq_cov"].update(cov)
            w.writerow([qname, assigned, level, rname, best, len(cands)])
    p.wait()

    stats = os.path.splitext(args.out)[0] + ".species_stats.tsv"
    with open(stats, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["taxon", "reads", "breadth_pct", "mean_depth", "mean_depth_mq20", "breadth_mq20_pct"])
        for name, st in sorted(per_species.items()):
            # genome length: sum over genomes assigned to this taxon
            glen = sum(v for g, v in lengths.items() if tax.get(g, ("", ""))[0] == name or tax.get(g, ("", ""))[1] == name)
            glen = glen or 1
            w.writerow([name, st["reads"],
                        f"{100.0 * len(st['cov']) / glen:.4f}",
                        f"{st['aligned_bases'] / glen:.4f}",
                        f"{st['mq_bases'] / glen:.4f}",
                        f"{100.0 * len(st['mq_cov']) / glen:.4f}"])
    print("wrote", args.out, "and", stats)


if __name__ == "__main__":
    main()
