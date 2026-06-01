#!/usr/bin/env bash
set -euo pipefail

AUTOTAX2_CMD="${AUTOTAX2_CMD:-autotax2}"
SINA_BIN="${SINA_BIN:-sina}"
VSEARCH_BIN="${VSEARCH_BIN:-vsearch}"

SILVA_FASTA=""
TYPE_STRAIN_METADATA=""
GTDB_AR53_TAXONOMY=""
GTDB_BAC120_TAXONOMY=""
DATASET_FASTA=""
OUTDIR=""
DOMAIN=""
DATASET_NAME=""
PREFIX=""
THREADS="4"
SINA_REFERENCE=""
SINA_SEARCH_DB=""
SEARCH_CANDIDATES="1"
SEARCH_MIN_SIM="0.5"
SEARCH_MAX_RESULT="10"
REQUIRE_SINA_CANDIDATES="0"
STRICT_VALIDATE="0"

usage() {
  cat <<'USAGE'
Run an optional real-tool autotax2 integration workflow.

Required:
  --silva-fasta PATH              official SILVA NR99 taxonomy FASTA
  --type-strain-metadata PATH     official SILVA full_metadata TSV/TSV.gz
  --gtdb-ar53-taxonomy PATH       GTDB r232 ar53_taxonomy TSV/TSV.gz
  --gtdb-bac120-taxonomy PATH     GTDB r232 bac120_taxonomy TSV/TSV.gz
  --dataset-fasta PATH            externally extracted SSU/16S FASTA
  --outdir PATH                   empty or nonexistent output build directory
  --domain Archaea|Bacteria
  --dataset-name NAME
  --prefix PREFIX

Optional:
  --threads INT                   default: 4
  --sina-bin PATH                 default: $SINA_BIN or sina
  --vsearch-bin PATH              default: $VSEARCH_BIN or vsearch
  --sina-reference PATH           optional SINA reference/PTDB path
  --sina-search-db PATH           optional SINA search database
  --search-candidates             default; ask SINA for nearest_slv candidates
  --no-search-candidates          skip SINA candidate search
  --search-min-sim FLOAT          default: 0.5
  --search-max-result INT         default: 10
  --require-sina-candidates       fail cluster if no SINA targets match registry
  --strict-validate               run final validate with --strict

Environment:
  AUTOTAX2_CMD                    command to run autotax2, default: autotax2
  SINA_BIN                        SINA executable, default: sina
  VSEARCH_BIN                     VSEARCH executable, default: vsearch
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --silva-fasta)
      SILVA_FASTA="${2:-}"; shift 2 ;;
    --type-strain-metadata)
      TYPE_STRAIN_METADATA="${2:-}"; shift 2 ;;
    --gtdb-ar53-taxonomy)
      GTDB_AR53_TAXONOMY="${2:-}"; shift 2 ;;
    --gtdb-bac120-taxonomy)
      GTDB_BAC120_TAXONOMY="${2:-}"; shift 2 ;;
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
    --sina-bin)
      SINA_BIN="${2:-}"; shift 2 ;;
    --vsearch-bin)
      VSEARCH_BIN="${2:-}"; shift 2 ;;
    --sina-reference)
      SINA_REFERENCE="${2:-}"; shift 2 ;;
    --sina-search-db)
      SINA_SEARCH_DB="${2:-}"; shift 2 ;;
    --search-candidates)
      SEARCH_CANDIDATES="1"; shift ;;
    --no-search-candidates)
      SEARCH_CANDIDATES="0"; shift ;;
    --search-min-sim)
      SEARCH_MIN_SIM="${2:-}"; shift 2 ;;
    --search-max-result)
      SEARCH_MAX_RESULT="${2:-}"; shift 2 ;;
    --require-sina-candidates)
      REQUIRE_SINA_CANDIDATES="1"; shift ;;
    --strict-validate)
      STRICT_VALIDATE="1"; shift ;;
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

require_nonempty() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "Expected non-empty output missing: $path" >&2
    exit 1
  fi
}

require_value "--silva-fasta" "$SILVA_FASTA"
require_value "--type-strain-metadata" "$TYPE_STRAIN_METADATA"
require_value "--gtdb-ar53-taxonomy" "$GTDB_AR53_TAXONOMY"
require_value "--gtdb-bac120-taxonomy" "$GTDB_BAC120_TAXONOMY"
require_value "--dataset-fasta" "$DATASET_FASTA"
require_value "--outdir" "$OUTDIR"
require_value "--domain" "$DOMAIN"
require_value "--dataset-name" "$DATASET_NAME"
require_value "--prefix" "$PREFIX"
require_file "$SILVA_FASTA"
require_file "$TYPE_STRAIN_METADATA"
require_file "$GTDB_AR53_TAXONOMY"
require_file "$GTDB_BAC120_TAXONOMY"
require_file "$DATASET_FASTA"
if [[ -n "$SINA_REFERENCE" ]]; then
  require_file "$SINA_REFERENCE"
fi
if [[ -n "$SINA_SEARCH_DB" ]]; then
  require_file "$SINA_SEARCH_DB"
fi

require_cmd "$AUTOTAX2_CMD"
require_cmd "$SINA_BIN"
require_cmd "$VSEARCH_BIN"
require_cmd gzip

if [[ -d "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Output directory exists and is not empty: $OUTDIR" >&2
  echo "Choose an empty/nonexistent --outdir for this integration run." >&2
  exit 2
fi
mkdir -p "$OUTDIR"

echo "[autotax2 integration] Checking external tool availability"
"$SINA_BIN" --version >/dev/null 2>&1 || true
"$VSEARCH_BIN" --version >/dev/null 2>&1 || true

echo "[autotax2 integration] init"
"$AUTOTAX2_CMD" init \
  --silva-fasta "$SILVA_FASTA" \
  --type-strain-metadata "$TYPE_STRAIN_METADATA" \
  --gtdb-ar53-taxonomy "$GTDB_AR53_TAXONOMY" \
  --gtdb-bac120-taxonomy "$GTDB_BAC120_TAXONOMY" \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  --vsearch-bin "$VSEARCH_BIN"

echo "[autotax2 integration] validate after SILVA setup"
"$AUTOTAX2_CMD" validate \
  --build "$OUTDIR" \
  --no-check-exports

echo "[autotax2 integration] prepare"
"$AUTOTAX2_CMD" prepare \
  --build "$OUTDIR" \
  --name "$DATASET_NAME" \
  --prefix "$PREFIX" \
  --fasta "$DATASET_FASTA" \
  --domain "$DOMAIN"

CLUSTER_ARGS=(
  cluster
  --build "$OUTDIR"
  --dataset "$DATASET_NAME"
  --threads "$THREADS"
  --vsearch-bin "$VSEARCH_BIN"
  --sina-bin "$SINA_BIN"
)
if [[ -n "$SINA_REFERENCE" ]]; then
  CLUSTER_ARGS+=(--sina-reference "$SINA_REFERENCE")
fi
if [[ "$SEARCH_CANDIDATES" == "1" ]]; then
  CLUSTER_ARGS+=(
    --search-candidates
    --search-min-sim "$SEARCH_MIN_SIM"
    --search-max-result "$SEARCH_MAX_RESULT"
  )
  if [[ -n "$SINA_SEARCH_DB" ]]; then
    CLUSTER_ARGS+=(--search-db "$SINA_SEARCH_DB")
  fi
else
  CLUSTER_ARGS+=(--no-search-candidates)
fi
if [[ "$REQUIRE_SINA_CANDIDATES" == "1" ]]; then
  CLUSTER_ARGS+=(--require-sina-candidates)
fi

echo "[autotax2 integration] cluster with automatic orientation"
"$AUTOTAX2_CMD" "${CLUSTER_ARGS[@]}"

echo "[autotax2 integration] place"
"$AUTOTAX2_CMD" place \
  --build "$OUTDIR" \
  --dataset "$DATASET_NAME"

echo "[autotax2 integration] summarize"
"$AUTOTAX2_CMD" summarize \
  --build "$OUTDIR" \
  --overwrite

echo "[autotax2 integration] export all"
"$AUTOTAX2_CMD" export all \
  --build "$OUTDIR" \
  --gzip

VALIDATE_ARGS=(validate --build "$OUTDIR")
if [[ "$STRICT_VALIDATE" == "1" ]]; then
  VALIDATE_ARGS+=(--strict)
fi

echo "[autotax2 integration] validate final"
"$AUTOTAX2_CMD" "${VALIDATE_ARGS[@]}"

DATASET_DIR="$(find "$OUTDIR/datasets" -mindepth 1 -maxdepth 1 -type d -name "*_${DATASET_NAME}" | sort | tail -n 1)"
if [[ -z "$DATASET_DIR" ]]; then
  echo "Could not find dataset output directory for $DATASET_NAME under $OUTDIR/datasets" >&2
  exit 1
fi

SINTAX="$OUTDIR/export/sintax/autotax2.sintax.fa.gz"
DADA2_GENUS="$OUTDIR/export/dada2/autotax2_toGenus_trainset.fa.gz"
DADA2_SPECIES="$OUTDIR/export/dada2/autotax2_assignSpecies.fa.gz"
QIIME2_FASTA="$OUTDIR/export/qiime2/reference_sequences.fasta.gz"
QIIME2_TAX="$OUTDIR/export/qiime2/reference_taxonomy.tsv"
GLOBAL_SUMMARY="$OUTDIR/reports/global_summary.tsv"
DATASET_DELTA="$OUTDIR/reports/dataset_delta_summary.tsv"
VALIDATION_MD="$OUTDIR/reports/validation_report.md"
EXPORT_VALIDATION="$OUTDIR/export/export_validation.tsv"
SILVA_EVIDENCE="$OUTDIR/silva/silva_unresolved_evidence.tsv"
PLACEMENT_EVIDENCE="$DATASET_DIR/placement_evidence.tsv"
CLUSTER_SUMMARY="$DATASET_DIR/cluster_search_summary.tsv"
SINA_ORIENTED="$DATASET_DIR/sina.oriented.fa"
SINA_SUMMARY="$DATASET_DIR/sina.summary.tsv"

echo "[autotax2 integration] post-run file checks"
for path in \
  "$SINTAX" \
  "$DADA2_GENUS" \
  "$DADA2_SPECIES" \
  "$QIIME2_FASTA" \
  "$QIIME2_TAX" \
  "$GLOBAL_SUMMARY" \
  "$DATASET_DELTA" \
  "$VALIDATION_MD" \
  "$EXPORT_VALIDATION" \
  "$SILVA_EVIDENCE" \
  "$PLACEMENT_EVIDENCE" \
  "$CLUSTER_SUMMARY" \
  "$SINA_ORIENTED" \
  "$SINA_SUMMARY"; do
  require_nonempty "$path"
done

if [[ "$SEARCH_CANDIDATES" == "1" ]]; then
  require_nonempty "$DATASET_DIR/sina.candidates.tsv"
  require_nonempty "$DATASET_DIR/sina_candidate_diagnostics.tsv"
fi

echo "[autotax2 integration] evidence/report checks"
grep -q "$(printf 'silva_unresolved_evidence_rows\t')" "$GLOBAL_SUMMARY"
grep -q "$(printf 'placement_evidence_rows\t')" "$DATASET_DELTA"
if [[ "$SEARCH_CANDIDATES" == "1" ]]; then
  grep -q "$(printf 'sina_candidate_queries\t')" "$DATASET_DELTA"
fi

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
