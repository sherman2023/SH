#!/usr/bin/env bash
set -euo pipefail

# Approach A: STAR -> featureCounts
# Env variables:
#   SAMPLES_TSV, DATA_DIR, REF_DIR, RESULTS_DIR, THREADS

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${SAMPLES_TSV:-${ROOT_DIR}/config/samples.tsv}"
DATA_DIR="${DATA_DIR:-${ROOT_DIR}/data}"
REF_DIR="${REF_DIR:-${ROOT_DIR}/ref}"
RESULTS_DIR="${RESULTS_DIR:-${ROOT_DIR}/results}"
THREADS="${THREADS:-8}"

FASTQ_DIR="${DATA_DIR}/fastq"
STAR_INDEX="${REF_DIR}/star_index"
GTF="${REF_DIR}/annotation.gtf"

OUT_ALIGN="${RESULTS_DIR}/approachA_star"
OUT_COUNTS="${RESULTS_DIR}/approachA_featurecounts"
mkdir -p "${OUT_ALIGN}" "${OUT_COUNTS}"

if [ ! -d "${STAR_INDEX}" ]; then
  echo "ERROR: STAR index missing at ${STAR_INDEX}" >&2
  exit 1
fi
if [ ! -f "${GTF}" ]; then
  echo "ERROR: GTF missing at ${GTF}" >&2
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
  if [ -f "${OUT_ALIGN}/${s}.Aligned.sortedByCoord.out.bam" ]; then
    echo "[skip] STAR already done: ${s}"
  else
    echo "[STAR] ${s}"
    STAR --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${r1}" "${r2}" \
         --readFilesCommand zcat \
         --runThreadN "${THREADS}" \
         --outFileNamePrefix "${OUT_ALIGN}/${s}." \
         --outSAMtype BAM SortedByCoordinate
  fi
  samtools index -@ "${THREADS}" "${OUT_ALIGN}/${s}.Aligned.sortedByCoord.out.bam"
done

echo "[featureCounts] counting genes"
bams=$(ls "${OUT_ALIGN}"/*.Aligned.sortedByCoord.out.bam)
featureCounts -T "${THREADS}" -p -B -C \
  -a "${GTF}" \
  -o "${OUT_COUNTS}/featureCounts.counts.tsv" \
  ${bams}

echo "[done] Approach A counts: ${OUT_COUNTS}/featureCounts.counts.tsv"
