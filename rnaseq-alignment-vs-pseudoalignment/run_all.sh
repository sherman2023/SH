#!/usr/bin/env bash
set -euo pipefail

# One-command entrypoint.
# Default: SMOKE_TEST=1 (fast synthetic dataset).
# Full run: SMOKE_TEST=0 (airway + GENCODE references).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${THREADS:-8}"
SMOKE_TEST="${SMOKE_TEST:-1}"
GENCODE_RELEASE="${GENCODE_RELEASE:-49}"

# -----------------------------
# Helpers
# -----------------------------
have_cmd() { command -v "$1" >/dev/null 2>&1; }

download_file() {
  local url="$1"
  local out="$2"
  if have_cmd curl; then
    curl -L --retry 5 --retry-delay 2 -o "$out" "$url"
  elif have_cmd wget; then
    wget -O "$out" "$url"
  else
    echo "ERROR: need curl or wget to download: $url" >&2
    exit 1
  fi
}

ensure_micromamba() {
  local mm_dir="${ROOT_DIR}/.micromamba/bin"
  local mm="${mm_dir}/micromamba"
  if [ -x "$mm" ]; then
    echo "[ok] micromamba exists: $mm"
    echo "$mm"
    return
  fi

  mkdir -p "$mm_dir"
  echo "[setup] downloading micromamba..."

  local uname_s
  uname_s="$(uname -s)"
  local uname_m
  uname_m="$(uname -m)"

  local platform=""
  case "${uname_s}_${uname_m}" in
    Linux_x86_64) platform="linux-64" ;;
    Linux_aarch64|Linux_arm64) platform="linux-aarch64" ;;
    Darwin_x86_64) platform="osx-64" ;;
    Darwin_arm64) platform="osx-arm64" ;;
    *)
      echo "ERROR: unsupported platform ${uname_s}_${uname_m} for auto micromamba install" >&2
      exit 1
      ;;
  esac

  local tarball="${ROOT_DIR}/.micromamba/micromamba.tar.bz2"
  download_file "https://micro.mamba.pm/api/micromamba/${platform}/latest" "$tarball"
  # The tarball contains bin/micromamba
  tar -xjf "$tarball" -C "${ROOT_DIR}/.micromamba"
  mv "${ROOT_DIR}/.micromamba/bin/micromamba" "$mm"
  chmod +x "$mm"
  rm -f "$tarball"
  echo "[ok] micromamba installed: $mm"
  echo "$mm"
}

ensure_env() {
  local micromamba="$1"
  local env_prefix="${ROOT_DIR}/.micromamba/env"
  if [ -d "$env_prefix" ]; then
    echo "[ok] env exists: $env_prefix"
    echo "$env_prefix"
    return
  fi
  echo "[setup] creating conda env at: $env_prefix"
  "${micromamba}" create -y -p "$env_prefix" -f "${ROOT_DIR}/envs/conda_env.yml"
  echo "$env_prefix"
}

mm_run() {
  local micromamba="$1"
  local env_prefix="$2"
  shift 2
  "${micromamba}" run -p "$env_prefix" "$@"
}

# -----------------------------
# Main
# -----------------------------
MICROMAMBA="$(ensure_micromamba)"
ENV_PREFIX="$(ensure_env "$MICROMAMBA")"

if [ "$SMOKE_TEST" -eq 1 ]; then
  echo "[mode] SMOKE_TEST=1 (synthetic toy dataset)"
  WORKDIR="${ROOT_DIR}/smoke"
  export SAMPLES_TSV="${WORKDIR}/config/samples.tsv"
  export DATA_DIR="${WORKDIR}/data"
  export REF_DIR="${WORKDIR}/ref"
  export RESULTS_DIR="${WORKDIR}/results"
  export THREADS

  mm_run "$MICROMAMBA" "$ENV_PREFIX" python "${ROOT_DIR}/scripts/prepare_smoke_dataset.py" --outdir "$WORKDIR"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/run_approachA_star_featurecounts.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/run_approachB_salmon.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" Rscript "${ROOT_DIR}/analysis/deseq2_compare.R"

  echo
  echo "[done] Smoke test results:"
  echo "  ${RESULTS_DIR}/deseq2/AvsB_summary_stats.tsv"
  echo "  ${RESULTS_DIR}/deseq2/AvsB_log2FC_scatter.png"
else
  echo "[mode] SMOKE_TEST=0 (airway + GENCODE references; large download/index build)"
  export SAMPLES_TSV="${ROOT_DIR}/config/samples.tsv"
  export DATA_DIR="${ROOT_DIR}/data"
  export REF_DIR="${ROOT_DIR}/ref"
  export RESULTS_DIR="${ROOT_DIR}/results"
  export THREADS
  export GENCODE_RELEASE

  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/download_airway_fastq.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/download_gencode_build_indices.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/run_approachA_star_featurecounts.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" bash "${ROOT_DIR}/scripts/run_approachB_salmon.sh"
  mm_run "$MICROMAMBA" "$ENV_PREFIX" Rscript "${ROOT_DIR}/analysis/deseq2_compare.R"

  echo
  echo "[done] Full run results:"
  echo "  ${RESULTS_DIR}/deseq2/AvsB_summary_stats.tsv"
  echo "  ${RESULTS_DIR}/deseq2/AvsB_log2FC_scatter.png"
fi
