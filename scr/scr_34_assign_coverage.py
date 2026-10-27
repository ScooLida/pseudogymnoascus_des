#!/usr/bin/env python3
"""Stage 3d - read assignment, PCR de-duplication and per-species coverage.

Works with competitive alignments from Paleomix (BWA-MEM with XA tags or
Bowtie2 end-to-end with multiple records per read).  The BAM is name-sorted
internally so that all alignments of one read are processed together.

A read is assigned to a species when its best alignment score is unique within
--delta-as; otherwise it goes to the genus (if all candidates share a genus) or
is marked ambiguous.  PCR duplicates are collapsed by 5' coordinate + strand.
"""
import argparse
import collections
import csv
import os
import re
import subprocess
import tempfile


def parse_taxonomy(path):
    tax = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sp = row["species"]
            tax[row["genome"]] = (sp, sp.split()[0] if sp else row["genome"])
    return tax


def genome_of(ref):
    return ref.split("|")[0]


def parse_cigar(cigar, pos):
    span, cov = 0, []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            cov.extend(range(pos + span, pos + span + n))
            span += n
        elif op in "DN":
            span += n
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

    hdr = subprocess.run(["samtools", "view", "-H", args.bam],
                         check=True, capture_output=True, text=True).stdout
    lengths = collections.Counter()
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:])
            lengths[genome_of(f["SN"])] += int(f["LN"])

    tmp = tempfile.NamedTemporaryFile(suffix=".nsort.bam", delete=False).name
    subprocess.run(["samtools", "sort", "-n", "-o", tmp, args.bam], check=True)

    per_species = collections.defaultdict(lambda: {
        "reads": 0, "cov": set(), "aligned_bases": 0, "mq_bases": 0, "mq_cov": set()})
    seen = set()

    out = open(args.out, "w", newline="")
    w = csv.writer(out, delimiter="\t")
    w.writerow(["read", "assigned", "level", "best_ref", "best_as", "n_candidates"])

    def process(qname, recs):
        cands = []
        for r in recs:
            cands.append((r["ref"], r["as"]))
            for alt in r["xa"]:
                cands.append(alt)
        if not cands:
            return
        best = max(a for _, a in cands)
        winners = {genome_of(r) for r, a in cands if a >= best - args.delta_as}
        genomes = [g for g in winners if g in tax]
        if not genomes:
            return
        species = {tax[g][0] for g in genomes}
        genera = {tax[g][1] for g in genomes}
        if len(species) == 1:
            assigned, level = species.pop(), "species"
        elif len(genera) == 1:
            assigned, level = genera.pop(), "genus"
        else:
            assigned, level = "ambiguous", "ambiguous"
        # best record among the winning genomes
        bestrec = max((r for r in recs if genome_of(r["ref"]) in winners),
                      key=lambda r: r["as"])
        _, cov = parse_cigar(bestrec["cigar"], bestrec["pos"])
        if not cov:
            return
        strand = "-" if bestrec["flag"] & 16 else "+"
        key = (assigned, bestrec["ref"], min(cov) if strand == "+" else max(cov), strand)
        if key in seen:
            w.writerow([qname, assigned, level, bestrec["ref"], bestrec["as"], len(cands)])
            return
        seen.add(key)
        st = per_species[assigned]
        st["reads"] += 1
        st["cov"].update(cov)
        st["aligned_bases"] += len(cov)
        if bestrec["mapq"] >= args.mapq_depth:
            st["mq_bases"] += len(cov)
            st["mq_cov"].update(cov)
        w.writerow([qname, assigned, level, bestrec["ref"], bestrec["as"], len(cands)])

    p = subprocess.Popen(["samtools", "view", "-F", "4", tmp],
                         stdout=subprocess.PIPE, text=True)
    cur, recs = None, []
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, flag, ref, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[6]
        if ref == "*":
            continue
        a = {"ref": ref, "pos": pos, "mapq": mapq, "cigar": cigar,
             "flag": flag, "as": 0, "xa": []}
        for t in f[11:]:
            if t.startswith("AS:i:"):
                a["as"] = int(t[5:])
            elif t.startswith("XA:Z:"):
                for alt in t[5:].split(";"):
                    p2 = alt.split(",")
                    if len(p2) >= 5:
                        try:
                            a["xa"].append((p2[0], int(p2[4])))
                        except ValueError:
                            pass
        if qname != cur:
            if cur is not None:
                process(cur, recs)
            cur, recs = qname, []
        recs.append(a)
    if cur is not None:
        process(cur, recs)
    p.wait()
    out.close()
    os.unlink(tmp)

    stats = os.path.splitext(args.out)[0] + ".species_stats.tsv"
    with open(stats, "w", newline="") as f:
        w2 = csv.writer(f, delimiter="\t")
        w2.writerow(["taxon", "reads", "breadth_pct", "mean_depth",
                     "mean_depth_mq20", "breadth_mq20_pct"])
        for name, st in sorted(per_species.items()):
            glen = sum(v for g, v in lengths.items()
                       if tax.get(g, ("", ""))[0] == name or tax.get(g, ("", ""))[1] == name) or 1
            w2.writerow([name, st["reads"],
                         f"{100.0 * len(st['cov']) / glen:.4f}",
                         f"{st['aligned_bases'] / glen:.4f}",
                         f"{st['mq_bases'] / glen:.4f}",
                         f"{100.0 * len(st['mq_cov']) / glen:.4f}"])
    print("wrote", args.out, "and", stats)


if __name__ == "__main__":
    main()
