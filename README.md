# redox_constantV

首先看手册：

- [docs/methods/USER_MANUAL_CN.md](/home/am3-peichenzhong-group/Documents/project/project_solv_structure/redox_potential/redox_constantV/docs/methods/USER_MANUAL_CN.md)

## Layout

- `configs/`
  - `md/`: MD/OpenMM configs
  - `redox/`: reduction-side screening configs
  - `oxidation/`: oxidation-side screening configs
  - `pathways/`: decomposition / pathway configs
  - `formulations/`: batch formulation inputs
  - `dft_workflows/`: PySCF workflow inputs
- `core/`: core Python modules used by the workflows
- `scripts/`
  - `cli/`: main runnable simulation/screening entry points
  - `dft/`: report generation and DFT-to-config utilities
  - `plotting/`: plotting utilities
  - `utils/`: one-off structure/build helpers
- `structures/`
  - `pdb_bank/`: single-molecule templates
  - `systems/`: assembled simulation systems
- `ff/`: force-field XML files
- `docs/`
  - `methods/`: current workflow docs
  - `provenance/`: parameter/source tracing
  - `notes/`: auxiliary notes
  - `reports/`: historical reports
- `tests/`: lightweight regression tests
- `results/`: generated outputs and archived debug runs
- `archive/legacy/`: older scripts and job files kept for reference

## Main entry points

- MD sweep:
  - `python scripts/cli/run_voltage_sweep.py --config configs/md/config_opls.json --redox-config configs/redox/config_redox.json --pdb structures/systems/start_with_electrodes.pdb`
- Batch screening:
  - `python scripts/cli/batch_screen_formulations.py --input configs/formulations/formulation_batch_example_extended.json --output-dir results/formulation_batch_test`
- DFT to screening config:
  - `python scripts/dft/generate_redox_config_from_dft.py --workflow configs/dft_workflows/example_ec_dmc_reduction.json --output-dir results/dft_interface_test`
