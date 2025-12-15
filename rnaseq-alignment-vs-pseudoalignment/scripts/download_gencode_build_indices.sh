#!/usr/bin/env bash
set -euo pipefail

# Downloads GENCODE GRCh38 reference and builds STAR + Salmon indices.
# Env variables:
#   REF_DIR (default: ref)
#   DATA_DIR (default: data)  [used to infer read length for STAR sjdbOverhang]
#   THREADS (default: 8)
#   GENCODE_RELEASE (default: 49)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${REF_DIR:-${ROOT_DIR}/ref}"
DATA_DIR="${DATA_DIR:-${ROOT_DIR}/data}"
THREADS="${THREADS:-8}"
GENCODE_RELEASE="${GENCODE_RELEASE:-49}"

mkdir -p "${REF_DIR}"

# Files
GENOME_FA_GZ="${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
GTF_GZ="${REF_DIR}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
TX_FA_GZ="${REF_DIR}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"

# Canonical GENCODE FTP base
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

download() {
  local url="$1"
  local out="$2"
  if [ -f "$out" ]; then
    echo "[skip] exists: $out"
    return
  fi
  echo "[download] $url"
  if command -v curl >/dev/null 2>&1; then
    curl -L --retry 5 --retry-delay 2 -o "$out" "$url"
  else
    wget -O "$out" "$url"
  fi
}

download "${BASE}/GRCh38.primary_assembly.genome.fa.gz" "${GENOME_FA_GZ}"
download "${BASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" "${GTF_GZ}"
download "${BASE}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz" "${TX_FA_GZ}"

# Decompress (keep gz as cache)
if [ ! -f "${REF_DIR}/genome.fa" ]; then
  echo "[gunzip] genome.fa"
  gunzip -c "${GENOME_FA_GZ}" > "${REF_DIR}/genome.fa"
fi
if [ ! -f "${REF_DIR}/annotation.gtf" ]; then
  echo "[gunzip] annotation.gtf"
  gunzip -c "${GTF_GZ}" > "${REF_DIR}/annotation.gtf"
fi
if [ ! -f "${REF_DIR}/transcripts.fa" ]; then
  echo "[gunzip] transcripts.fa"
  gunzip -c "${TX_FA_GZ}" > "${REF_DIR}/transcripts.fa"
fi

# Build tx2gene.tsv from GTF (python helper)
if [ ! -f "${REF_DIR}/tx2gene.tsv" ]; then
  echo "[build] tx2gene.tsv"
  python "${ROOT_DIR}/scripts/make_tx2gene_from_gtf.py" \
    --gtf "${REF_DIR}/annotation.gtf" \
    --out "${REF_DIR}/tx2gene.tsv"
fi

# Determine read length from first FASTQ (for sjdbOverhang = read_len - 1)
FASTQ_DIR="${DATA_DIR}/fastq"
first_fq=$(ls "${FASTQ_DIR}"/*_1.fastq.gz 2>/dev/null | head -n 1 || true)
if [ -z "${first_fq}" ]; then
  echo "ERROR: no FASTQs found in ${FASTQ_DIR}. Run download_airway_fastq.sh first." >&2
  exit 1
fi
read_len=$(zcat "${first_fq}" | awk 'NR==2 {print length($0); exit}')
sjdb_overhang=$((read_len - 1))
echo "[info] detected read length: ${read_len} => sjdbOverhang=${sjdb_overhang}"

# STAR index
STAR_DIR="${REF_DIR}/star_index"
if [ ! -d "${STAR_DIR}" ]; then
  mkdir -p "${STAR_DIR}"
  echo "[build] STAR genome index"
  STAR --runThreadN "${THREADS}" \
       --runMode genomeGenerate \
       --genomeDir "${STAR_DIR}" \
       --genomeFastaFiles "${REF_DIR}/genome.fa" \
       --sjdbGTFfile "${REF_DIR}/annotation.gtf" \
       --sjdbOverhang "${sjdb_overhang}"
else
  echo "[skip] STAR index exists: ${STAR_DIR}"
fi

# Salmon index
SALMON_DIR="${REF_DIR}/salmon_index"
if [ ! -d "${SALMON_DIR}" ]; then
  echo "[build] Salmon transcriptome index"
  salmon index -t "${REF_DIR}/transcripts.fa" -i "${SALMON_DIR}" -p "${THREADS}" --gencode
else
  echo "[skip] Salmon index exists: ${SALMON_DIR}"
fi

echo "[done] references + indices in ${REF_DIR}"
