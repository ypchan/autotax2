#!/usr/bin/env bash
set -euo pipefail

AUTOTAX2_CMD="${AUTOTAX2_CMD:-autotax2}"
BARRNAP_EXPECTED_VERSION="1.10.5"

SILVA_FASTA=""
DATASET_FASTA=""
OUTDIR=""
DOMAIN=""
DATASET_NAME=""
PREFIX=""
THREADS="4"
STRICT_TOOL_VERSION="0"

usage() {
  cat <<'USAGE'
Run an optional real-tool autotax2 integration workflow.

Required:
  --silva-fasta PATH
  --dataset-fasta PATH
  --outdir PATH
  --domain Archaea|Bacteria
  --dataset-name NAME
  --prefix PREFIX

Optional:
  --threads INT              default: 4
  --strict-tool-version      fail if barrnap is not 1.10.5

Environment:
  AUTOTAX2_CMD               command to run autotax2, default: autotax2
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --silva-fasta)
      SILVA_FASTA="${2:-}"; shift 2 ;;
    --dataset-fasta)
      DATASET_FASTA="${2:-}"; shift 2 ;;
    --outdir)
      OUTDIR="${2:-}"; shift 2 ;;
    --domain)
      DOMAIN="${2:-}"; shift 2 ;;
    --dataset-name)
      DATASET_NAME="${2:-}"; shift 2 ;;
    --prefix)
      PREFIX="${2:-}"; shift 2 ;;
    --threads)
      THREADS="${2:-}"; shift 2 ;;
    --strict-tool-version)
      STRICT_TOOL_VERSION="1"; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

require_value() {
  local name="$1"
  local value="$2"
  if [[ -z "$value" ]]; then
    echo "Missing required argument: $name" >&2
    usage >&2
    exit 2
  fi
}

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "Input file not found: $path" >&2
    exit 2
  fi
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found on PATH: $cmd" >&2
    exit 2
  fi
}

first_version() {
  sed -nE 's/.*([0-9]+(\.[0-9]+)+).*/\1/p' | head -n 1
}

require_nonempty() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "Expected non-empty output missing: $path" >&2
    exit 1
  fi
}

require_value "--silva-fasta" "$SILVA_FASTA"
require_value "--dataset-fasta" "$DATASET_FASTA"
require_value "--outdir" "$OUTDIR"
require_value "--domain" "$DOMAIN"
require_value "--dataset-name" "$DATASET_NAME"
require_value "--prefix" "$PREFIX"
require_file "$SILVA_FASTA"
require_file "$DATASET_FASTA"

require_cmd "$AUTOTAX2_CMD"
require_cmd barrnap
require_cmd sina
require_cmd vsearch
require_cmd gzip

if [[ -d "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Output directory exists and is not empty: $OUTDIR" >&2
  echo "Choose an empty/nonexistent --outdir for this integration run." >&2
  exit 2
fi
mkdir -p "$OUTDIR"

echo "[autotax2 integration] Checking external tool versions"
BARRNAP_VERSION="$(barrnap --version 2>&1 | first_version || true)"
sina --version >/dev/null 2>&1 || true
vsearch --version >/dev/null 2>&1 || true

if [[ "$BARRNAP_VERSION" != "$BARRNAP_EXPECTED_VERSION" ]]; then
  message="barrnap version is ${BARRNAP_VERSION:-unknown}; expected $BARRNAP_EXPECTED_VERSION"
  if [[ "$STRICT_TOOL_VERSION" == "1" ]]; then
    echo "ERROR: $message" >&2
    exit 1
  fi
  echo "WARNING: $message" >&2
fi

PREPARE_STRICT_ARGS=()
if [[ "$STRICT_TOOL_VERSION" == "1" || "$BARRNAP_VERSION" == "$BARRNAP_EXPECTED_VERSION" ]]; then
  PREPARE_STRICT_ARGS+=(--strict-tool-version)
fi

echo "[autotax2 integration] init"
"$AUTOTAX2_CMD" init \
  --silva-fasta "$SILVA_FASTA" \
  --outdir "$OUTDIR" \
  --domain "$DOMAIN"

echo "[autotax2 integration] resolve-silva"
"$AUTOTAX2_CMD" resolve-silva \
  --build "$OUTDIR" \
  --threads "$THREADS"

echo "[autotax2 integration] prepare-dataset"
"$AUTOTAX2_CMD" prepare-dataset \
  --build "$OUTDIR" \
  --name "$DATASET_NAME" \
  --prefix "$PREFIX" \
  --fasta "$DATASET_FASTA" \
  --domain "$DOMAIN" \
  --threads "$THREADS" \
  "${PREPARE_STRICT_ARGS[@]}"

echo "[autotax2 integration] orient-sina"
"$AUTOTAX2_CMD" orient-sina \
  --build "$OUTDIR" \
  --dataset "$DATASET_NAME" \
  --threads "$THREADS"

echo "[autotax2 integration] cluster-search"
"$AUTOTAX2_CMD" cluster-search \
  --build "$OUTDIR" \
  --dataset "$DATASET_NAME" \
  --threads "$THREADS"

echo "[autotax2 integration] place"
"$AUTOTAX2_CMD" place \
  --build "$OUTDIR" \
  --dataset "$DATASET_NAME"

echo "[autotax2 integration] export all"
"$AUTOTAX2_CMD" export all \
  --build "$OUTDIR" \
  --gzip

echo "[autotax2 integration] summarize"
"$AUTOTAX2_CMD" summarize \
  --build "$OUTDIR"

echo "[autotax2 integration] validate"
"$AUTOTAX2_CMD" validate \
  --build "$OUTDIR"

SINTAX="$OUTDIR/export/sintax/autotax2.sintax.fa.gz"
DADA2_GENUS="$OUTDIR/export/dada2/autotax2_toGenus_trainset.fa.gz"
DADA2_SPECIES="$OUTDIR/export/dada2/autotax2_assignSpecies.fa.gz"
QIIME2_FASTA="$OUTDIR/export/qiime2/reference_sequences.fasta.gz"
QIIME2_TAX="$OUTDIR/export/qiime2/reference_taxonomy.tsv"
GLOBAL_SUMMARY="$OUTDIR/reports/global_summary.tsv"
DATASET_DELTA="$OUTDIR/reports/dataset_delta_summary.tsv"
VALIDATION_MD="$OUTDIR/reports/validation_report.md"

echo "[autotax2 integration] post-run file checks"
for path in \
  "$SINTAX" \
  "$DADA2_GENUS" \
  "$DADA2_SPECIES" \
  "$QIIME2_FASTA" \
  "$QIIME2_TAX" \
  "$GLOBAL_SUMMARY" \
  "$DATASET_DELTA" \
  "$VALIDATION_MD"; do
  require_nonempty "$path"
done

echo "[autotax2 integration] SINTAX format checks"
gzip -cd "$SINTAX" | grep -m1 ';tax=d:' >/dev/null
if gzip -cd "$SINTAX" | awk -F';tax=' '/^>/{split($2, a, ";"); if (a[1] ~ /(g__|s__)/) bad=1} END{exit bad}'; then
  :
else
  echo "SINTAX tax values contain g__ or s__ prefixes." >&2
  exit 1
fi

echo "[autotax2 integration] DADA2 assignSpecies format checks"
gzip -cd "$DADA2_SPECIES" | grep -m1 '^>' >/dev/null
if gzip -cd "$DADA2_SPECIES" | grep '^>' | grep -q ';'; then
  echo "DADA2 assignSpecies headers contain semicolons." >&2
  exit 1
fi
if gzip -cd "$DADA2_SPECIES" | grep '^>' | grep -Eq 'g__|s__'; then
  echo "DADA2 assignSpecies headers contain g__ or s__ prefixes." >&2
  exit 1
fi

echo "[autotax2 integration] QIIME2 taxonomy format checks"
expected_header="$(printf 'Feature ID\tTaxon')"
actual_header="$(head -n 1 "$QIIME2_TAX")"
if [[ "$actual_header" != "$expected_header" ]]; then
  echo "QIIME2 taxonomy header mismatch: $actual_header" >&2
  exit 1
fi

echo "[autotax2 integration] OK: $OUTDIR"
