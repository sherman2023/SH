# RNA-seq differential expression: alignment-based vs pseudoalignment/quasi-mapping

This repository provides a **fully scripted**, reproducible comparison of two RNA-seq DE workflows:

- **Approach A (alignment-based)**: STAR → featureCounts → DESeq2  
- **Approach B (pseudoalignment/quasi-mapping)**: Salmon → tximport → DESeq2  

## One-command run (recommended)

From a fresh clone:

```bash
bash run_all.sh
```

By default this runs a **SMOKE TEST** (small synthetic dataset) end-to-end to verify the toolchain and workflow wiring.

### Full airway dataset run (heavier)
This downloads airway FASTQs from SRA, downloads GENCODE/GRCh38 references, builds STAR+Salmon indices, and runs both workflows:

```bash
SMOKE_TEST=0 bash run_all.sh
```

Optional:
```bash
THREADS=16 SMOKE_TEST=0 GENCODE_RELEASE=49 bash run_all.sh
```

## Outputs

Smoke test outputs:
- `smoke/results/deseq2/AvsB_summary_stats.tsv`
- `smoke/results/deseq2/AvsB_log2FC_scatter.png`

Full run outputs:
- `results/deseq2/AvsB_summary_stats.tsv`
- `results/deseq2/AvsB_log2FC_scatter.png`

## Reproducibility notes

- The workflow is **fully scripted** (one-command entrypoint).
- Software dependencies are pinned in `envs/conda_env.yml`.
- The DE analysis writes `sessionInfo.txt` under the results directory.

## Repository layout

- `run_all.sh`: one-command entrypoint
- `scripts/`: workflow steps (download, indexing, approach A/B)
- `analysis/deseq2_compare.R`: DESeq2 + A-vs-B comparison
- `config/samples.tsv`: airway sample sheet (used in full run)
