# redox_potential

OpenMM-based constant-potential MD and voltage sweep workflow for redox-potential studies.

## Included

- Core simulation scripts and configs
- Force-field XML files in `ffdir/`
- Molecule PDB templates in `pdb_bank/`
- Input structure `start_with_electrodes.pdb`
- Workflow and reference docs
- Lightweight result summary in `results/voltage_sweep_20260311_224838/`

## Excluded from git

Large generated trajectories and transient run directories are ignored.

## Key result

See:
- `results/voltage_sweep_20260311_224838/voltage_sweep_results.csv`
- `results/voltage_sweep_20260311_224838/voltage_sweep_summary.png`
- `results/voltage_sweep_20260311_224838/RUN_SUMMARY.md`
