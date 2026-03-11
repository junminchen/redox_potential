# redox_potential

OpenMM-based constant-potential MD and redox screening workflow for battery electrolytes.

## Included

- Core simulation scripts and configs
- Force-field XML files in `ffdir/`
- Molecule PDB templates in `pdb_bank/`
- Input structure `start_with_electrodes.pdb`
- Workflow and reference docs
- Voltage-sweep summary in `results/voltage_sweep_20260311_224838/`
- Experiment-facing electrolyte screening report in `results/experimental_redox_report_20260311_234115/`

## Key Entry Points

- `run_voltage_sweep.py`: constant-potential sweep driver
- `generate_experimental_redox_report.py`: build LSV/CV-like plots and a report
- `config_redox_electrolyte_effective.json`: effective EC/DMC onset-screening model
- `config_redox_real_template.json`: template for real redox-active molecules
- `HIGH_VOLTAGE_STABILITY_WORKFLOW.md`: roadmap toward formulation-level 5 V stability screening

## Main Outputs

- `results/experimental_redox_report_20260311_234115/electrolyte_lsv_cv_like.png`
- `results/experimental_redox_report_20260311_234115/REPORT.md`
- `results/experimental_redox_report_20260311_234115/electrolyte_redox_summary.csv`

## Excluded from git

Large generated trajectories and transient run directories are ignored.
