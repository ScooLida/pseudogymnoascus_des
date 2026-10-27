#!/usr/bin/env python3
"""Stage 3a - build the competitive-mapping panel.

The panel contains the selected Pseudogymnoascus genomes (one representative per
ANI cluster) plus the listed close relatives (or their closest NCBI proxies).
P. destructans is deliberately EXCLUDED from the panel and written to a separate
reference used for read/tree projection.
"""
import argparse
import csv
import os
import re
import sys

DESTRUCTANS_RE = re.compile(r"Pseudogymnoascus\s+destructans", re.I)


def read_asm_summary(path):
    """assembly_summary*.txt -> {gca: {'organism','taxid'}}"""
    meta = {}
    if not path or not os.path.exists(path):
        return meta
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 8:
                continue
            meta[a[0]] = {"organism": a[7], "taxid": a[5]}
    return meta


def gca_of(name):
    m = re.match(r"(GC[AF]_\d+\.\d+)", name)
    return m.group(1) if m else name


def fasta_records(path):
    with open(path) as f:
        header, seq = None, []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", required=True, help="genome names, one per line")
    ap.add_argument("--genome-dir", required=True, help="dir with <name>.fna")
    ap.add_argument("--assembly-summary", default="")
    ap.add_argument("--relatives", default="", help="TSV: species<TAB>fasta<TAB>taxid")
    ap.add_argument("--destructans", default="", help="destructans reference FASTA")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    meta = read_asm_summary(args.assembly_summary)
    panel_path = os.path.join(args.out_dir, "panel.fasta")
    tax_path = os.path.join(args.out_dir, "taxonomy.tsv")

    selected = [l.strip() for l in open(args.selected) if l.strip()]
    records = []
    excluded = []
    for name in selected:
        gca = gca_of(name)
        info = meta.get(gca, {})
        species = info.get("organism", name)
        if DESTRUCTANS_RE.search(species):
            excluded.append((name, species))
            continue
        tid = info.get("taxid", "")
        fna = os.path.join(args.genome_dir, name if name.endswith(".fna") else name + ".fna")
        if not os.path.exists(fna):
            print("WARNING missing", fna, file=sys.stderr)
            continue
        records.append((gca, species, tid, fna))

    if args.relatives and os.path.exists(args.relatives):
        with open(args.relatives) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                species, fasta, tid = (line.rstrip("\n").split("\t") + ["", "", ""])[:3]
                if DESTRUCTANS_RE.search(species):
                    continue
                if not os.path.exists(fasta):
                    print("WARNING missing relative", fasta, file=sys.stderr)
                    continue
                records.append((os.path.basename(fasta).split(".")[0], species, tid, fasta))

    with open(panel_path, "w") as out, open(tax_path, "w", newline="") as tf:
        tw = csv.writer(tf, delimiter="\t")
        tw.writerow(["genome", "species", "taxid", "source"])
        for gca, species, tid, fna in records:
            for h, seq in fasta_records(fna):
                h = h.split()[0]
                out.write(f">{gca}|{h}|{species.replace(' ', '_')}|taxid={tid}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
            tw.writerow([gca, species, tid, fna])
    print(f"panel: {panel_path} ({len(records)} genomes)")
    for n, s in excluded:
        print(f"excluded destructans: {n} ({s})")

    if args.destructans and os.path.exists(args.destructans):
        dpath = os.path.join(args.out_dir, "destructans_ref.fasta")
        with open(dpath, "w") as out:
            for h, seq in fasta_records(args.destructans):
                out.write(f">destructans|{h.split()[0]}|Pseudogymnoascus_destructans\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
        print("destructans reference:", dpath)


if __name__ == "__main__":
    main()
