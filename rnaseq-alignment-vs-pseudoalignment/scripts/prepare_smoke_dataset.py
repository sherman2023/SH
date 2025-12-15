#!/usr/bin/env python3
import argparse, gzip, os, random
from pathlib import Path

DNA = "ACGT"

def rc(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]

def write_fasta(path: Path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_gtf(path: Path, features):
    with open(path, "w") as f:
        for row in features:
            f.write("\t".join(map(str, row)) + "\n")

def rand_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(DNA) for _ in range(n))

def simulate_paired_reads(transcript_seq: str, n_pairs: int, read_len: int, frag_len: int, rng: random.Random):
    # Generate Illumina-style paired-end reads from a transcript sequence.
    # read1 = fragment start -> read_len
    # read2 = reverse-complement of fragment end -> read_len
    L = len(transcript_seq)
    if frag_len < 2 * read_len:
        frag_len = 2 * read_len
    if L < frag_len:
        raise ValueError("Transcript shorter than fragment length.")
    reads = []
    for i in range(n_pairs):
        start = rng.randint(0, L - frag_len)
        frag = transcript_seq[start:start+frag_len]
        r1 = frag[:read_len]
        r2 = rc(frag[-read_len:])
        reads.append((r1, r2))
    return reads

def write_fastq_gz(path: Path, seqs, prefix: str):
    # seqs: list of sequences
    with gzip.open(path, "wt") as f:
        for i, s in enumerate(seqs, start=1):
            f.write(f"@{prefix}:{i}\n")
            f.write(s + "\n")
            f.write("+\n")
            f.write("I" * len(s) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="Output directory (e.g., smoke/)")
    ap.add_argument("--read-len", type=int, default=75)
    ap.add_argument("--frag-len", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--pairs-per-sample", type=int, default=4000)
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    data_dir = outdir / "data" / "fastq"
    ref_dir = outdir / "ref"
    cfg_dir = outdir / "config"
    results_dir = outdir / "results"
    for d in [data_dir, ref_dir, cfg_dir, results_dir]:
        d.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)

    # Toy genome: chr1 length 4000
    genome_len = 4000
    genome = list(rand_seq(genome_len, seed=args.seed))

    # Define 3 genes on chr1 (single-exon transcripts)
    genes = [
        ("GENE1", "TX1",  101,  700, "+"),
        ("GENE2", "TX2", 1001, 1600, "+"),
        ("GENE3", "TX3", 2001, 2600, "+"),
    ]

    # Fill gene regions with deterministic sequences (so transcripts are exact substrings)
    for gi, (gene, tx, start, end, strand) in enumerate(genes, start=1):
        seq = rand_seq(end - start + 1, seed=args.seed + 100 * gi)
        genome[start-1:end] = list(seq)

    genome_seq = "".join(genome)

    # Write genome fasta
    write_fasta(ref_dir / "genome.fa", [("chr1", genome_seq)])

    # Write transcripts fasta (extract from genome)
    tx_records = []
    tx2gene = []
    for gene, tx, start, end, strand in genes:
        tx_seq = genome_seq[start-1:end]
        tx_records.append((tx, tx_seq))
        tx2gene.append((tx, gene))
    write_fasta(ref_dir / "transcripts.fa", tx_records)

    # Write tx2gene.tsv
    with open(ref_dir / "tx2gene.tsv", "w") as f:
        for tx, gene in tx2gene:
            f.write(f"{tx}\t{gene}\n")

    # Write simple GTF (gene, transcript, exon)
    gtf_rows = []
    source = "smoke"
    for gene, tx, start, end, strand in genes:
        gene_attr = f'gene_id "{gene}"; gene_name "{gene}";'
        tx_attr = f'gene_id "{gene}"; transcript_id "{tx}";'
        exon_attr = f'gene_id "{gene}"; transcript_id "{tx}"; exon_number "1";'
        gtf_rows.append(["chr1", source, "gene", start, end, ".", strand, ".", gene_attr])
        gtf_rows.append(["chr1", source, "transcript", start, end, ".", strand, ".", tx_attr])
        gtf_rows.append(["chr1", source, "exon", start, end, ".", strand, ".", exon_attr])
    write_gtf(ref_dir / "annotation.gtf", gtf_rows)

    # Sample sheet: mimic airway structure with 4 cell lines and treated/untrt
    samples = []
    cell_lines = ["C1", "C2", "C3", "C4"]
    run_id = 1
    for cell in cell_lines:
        for dex in ["untrt", "trt"]:
            samples.append((f"SMK{run_id}", cell, dex))
            run_id += 1

    with open(cfg_dir / "samples.tsv", "w") as f:
        f.write("run\tcell\tdex\n")
        for r, cell, dex in samples:
            f.write(f"{r}\t{cell}\t{dex}\n")

    # Simulate reads from transcripts with clear differential expression:
    # - GENE1 up in treated
    # - GENE2 constant
    # - GENE3 down in treated
    tx_seq = {name: seq for name, seq in tx_records}

    for r, cell, dex in samples:
        # total library size varies by cell line slightly
        base_pairs = args.pairs_per_sample + (cell_lines.index(cell) * 300)
        n_pairs = base_pairs

        if dex == "untrt":
            frac = {"TX1": 0.20, "TX2": 0.60, "TX3": 0.20}
        else:
            frac = {"TX1": 0.65, "TX2": 0.25, "TX3": 0.10}

        reads1, reads2 = [], []
        for tx, f_tx in frac.items():
            n_tx = int(n_pairs * f_tx)
            pairs = simulate_paired_reads(tx_seq[tx], n_tx, args.read_len, args.frag_len, rng)
            for r1, r2 in pairs:
                reads1.append(r1)
                reads2.append(r2)

        # Shuffle (to avoid blocks by gene)
        combined = list(zip(reads1, reads2))
        rng.shuffle(combined)
        reads1 = [a for a, b in combined]
        reads2 = [b for a, b in combined]

        out1 = data_dir / f"{r}_1.fastq.gz"
        out2 = data_dir / f"{r}_2.fastq.gz"
        write_fastq_gz(out1, reads1, prefix=r)
        write_fastq_gz(out2, reads2, prefix=r)

    print(f"[smoke] wrote dataset to: {outdir}")
    print(f"[smoke] samples: {cfg_dir / 'samples.tsv'}")
    print(f"[smoke] fastq dir: {data_dir}")
    print(f"[smoke] ref dir:   {ref_dir}")

if __name__ == "__main__":
    main()
