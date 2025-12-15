#!/usr/bin/env bash
set -euo pipefail

# Downloads airway FASTQs from SRA using prefetch + fasterq-dump.
# Uses environment variables:
#   SAMPLES_TSV (default: config/samples.tsv)
#   DATA_DIR (default: data)
#   THREADS (default: 8)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
DATA_DIR="${DATA_DIR:-${ROOT_DIR}/data}"
THREADS="${THREADS:-8}"

FASTQ_DIR="${DATA_DIR}/fastq"
mkdir -p "${FASTQ_DIR}"

echo "[airway] using samples: ${SAMPLES_TSV}"
runs=$(cut -f1 "${SAMPLES_TSV}" | tail -n +2)

cd "${FASTQ_DIR}"

for s in ${runs}; do
  if [ -f "${s}_1.fastq.gz" ] && [ -f "${s}_2.fastq.gz" ]; then
    echo "[skip] ${s} FASTQs already exist"
    continue
  fi
  echo "[download] ${s}"
  prefetch "${s}"
  fasterq-dump --split-files --threads "${THREADS}" --outdir . "${s}"
  pigz -p "${THREADS}" "${s}_1.fastq" "${s}_2.fastq"
done

echo "[done] airway FASTQs in ${FASTQ_DIR}"
