#!/usr/bin/env bash
# Run every benchmark stage in sequence.  Stages share the cached features and
# write independent result tables, so a failure part-way through loses only the
# stage that was running.
set -euo pipefail
PY="${PYTHON:-/tmp/ptg-venv/bin/python}"
cd "$(dirname "$0")/.."
mkdir -p logs
for stage in "$@"; do
  echo "=== ${stage} $(date -Is) ==="
  "$PY" scripts/run_benchmark.py --stage "$stage" --n-repeats "${N_REPEATS:-10}" \
    > "logs/${stage}.log" 2>&1
  echo "=== ${stage} done $(date -Is) ==="
done
