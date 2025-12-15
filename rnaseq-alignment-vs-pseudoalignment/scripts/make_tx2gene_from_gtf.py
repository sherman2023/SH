#!/usr/bin/env python3
import argparse, re
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gtf = Path(args.gtf)
    out = Path(args.out)

    gene_re = re.compile(r'gene_id "([^"]+)"')
    tx_re = re.compile(r'transcript_id "([^"]+)"')

    pairs = set()
    with gtf.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            g = gene_re.search(attrs)
            t = tx_re.search(attrs)
            if g and t:
                pairs.add((t.group(1), g.group(1)))

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as w:
        for tx, gene in sorted(pairs):
            w.write(f"{tx}\t{gene}\n")

    print(f"[tx2gene] wrote {len(pairs)} mappings to {out}")

if __name__ == "__main__":
    main()
