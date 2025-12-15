#!/usr/bin/env bash
set -euo pipefail

# Approach B: Salmon quant
# Env variables:
#   SAMPLES_TSV, DATA_DIR, REF_DIR, RESULTS_DIR, THREADS

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
DATA_DIR="${DATA_DIR:-${ROOT_DIR}/data}"
REF_DIR="${REF_DIR:-${ROOT_DIR}/ref}"
RESULTS_DIR="${RESULTS_DIR:-${ROOT_DIR}/results}"
THREADS="${THREADS:-8}"

FASTQ_DIR="${DATA_DIR}/fastq"
SALMON_INDEX="${REF_DIR}/salmon_index"
OUT_QUANT="${RESULTS_DIR}/approachB_salmon"
mkdir -p "${OUT_QUANT}"

if [ ! -d "${SALMON_INDEX}" ]; then
  echo "ERROR: Salmon index missing at ${SALMON_INDEX}" >&2
  exit 1
fi

runs=$(cut -f1 "${SAMPLES_TSV}" | tail -n +2)

for s in ${runs}; do
  r1="${FASTQ_DIR}/${s}_1.fastq.gz"
  r2="${FASTQ_DIR}/${s}_2.fastq.gz"
  if [ ! -f "${r1}" ] || [ ! -f "${r2}" ]; then
    echo "ERROR: missing FASTQs for ${s} in ${FASTQ_DIR}" >&2
    exit 1
  fi
  if [ -f "${OUT_QUANT}/${s}/quant.sf" ]; then
    echo "[skip] Salmon already done: ${s}"
    continue
  fi
  echo "[Salmon] ${s}"
  salmon quant -i "${SALMON_INDEX}" -l A \
    -1 "${r1}" -2 "${r2}" \
    -p "${THREADS}" \
    --gcBias \
    -o "${OUT_QUANT}/${s}"
done

echo "[done] Approach B quants in: ${OUT_QUANT}/<sample>/quant.sf"
