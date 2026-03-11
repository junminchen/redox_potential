#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

usage() {
  cat <<'USAGE'
Usage:
  bash run_examples.sh test
  bash run_examples.sh batch
  bash run_examples.sh dft-validate
  bash run_examples.sh dft-run
  bash run_examples.sh sweep-help
  bash run_examples.sh all
USAGE
}

run_test() {
  python tests/test_redox_framework.py
}

run_batch() {
  python scripts/cli/batch_screen_formulations.py \
    --input configs/formulations/formulation_batch_example_extended.json \
    --output-dir results/formulation_batch_examples
}

run_dft_validate() {
  python - <<'PY'
import json
from pathlib import Path

workflow = Path('configs/dft_workflows/example_ec_dmc_reduction.json')
data = json.loads(workflow.read_text())
print(f'Workflow: {workflow}')
for molecule in data.get('molecules', []):
    geometry = (workflow.parent / molecule['geometry_file']).resolve()
    status = 'OK' if geometry.exists() else 'MISSING'
    print(f"[{status}] {molecule['name']}: {geometry}")
PY
}

run_dft_real() {
  python scripts/dft/generate_redox_config_from_dft.py \
    --workflow configs/dft_workflows/example_ec_dmc_reduction.json \
    --output-dir results/dft_interface_example
}

run_sweep_help() {
  python scripts/cli/run_voltage_sweep.py --help
}

cmd="${1:-}"
case "$cmd" in
  test)
    run_test
    ;;
  batch)
    run_batch
    ;;
  dft-validate)
    run_dft_validate
    ;;
  dft-run)
    run_dft_real
    ;;
  sweep-help)
    run_sweep_help
    ;;
  all)
    run_test
    run_batch
    run_dft_validate
    ;;
  *)
    usage
    exit 1
    ;;
esac
