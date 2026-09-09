#!/usr/bin/env bash
# ChromBERT-scAI imputation example
# Edit paths below before running.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="${INPUT:-${ROOT}/data/example.h5ad}"
OUTPUT="${OUTPUT:-${ROOT}/output/example_imputed_ALS.h5ad}"
EMB_DIR="${EMB_DIR:-/path/to/atac_embs_dim.zarr}"

mkdir -p "$(dirname "$OUTPUT")"

python -u "${ROOT}/chrombert_scai/impute_als.py" \
  --input "${INPUT}" \
  --output "${OUTPUT}" \
  --emb_dir "${EMB_DIR}" \
  --emb_dim 768 \
  --topics 128 \
  --lambdaA 1 \
  --lambdaB 0.1 \
  --als_iters 20 \
  --cg_iters 500 \
  --cg_tol 1e-5 \
  --binary \
  --seed 42

echo "Wrote ${OUTPUT}"
