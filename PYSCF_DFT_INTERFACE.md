# PySCF DFT Interface

## Purpose

This interface connects molecule-level DFT calculations to the existing
`RedoxMC` screening configs used in this repository.

It does not replace formulation-level calibration. It provides a reproducible
bridge:

1. molecular geometry
2. DFT single-point energies for multiple charge states
3. relative free-energy offsets
4. `state_free_energy_offsets_kjmol` for screening

## What It Produces

The workflow generator writes two files:

- `dft_raw_results.json`
  - raw state energies
  - relative energies
  - mapped effective offsets
  - predicted half-wave voltages implied by the current mapping
- `redox_config_generated.json`
  - directly consumable by `RedoxMC`
  - uses the same `molecules -> state_free_energy_offsets_kjmol` structure as the existing hand-built configs

## Mapping Principle

For each non-reference state:

- `DeltaG_DFT = G(state) - G(reference)`
- `DeltaG_effective = DeltaG_DFT + empirical_shift`

The screening model then uses:

- `E_state(V) = DeltaG_effective - Delta z * F * V`

So the implied crossing voltage is:

- `V_half_wave_vs_vacuum = DeltaG_effective / (Delta z * F)`
- `V_half_wave_vs_li = V_half_wave_vs_vacuum + vacuum_to_li_offset`

where:

- `Delta z = z_state - z_reference`
- `F = 96.48533212331002 kJ mol^-1 V^-1`

## Why The Empirical Shift Exists

Raw isolated-molecule DFT values are usually not enough for electrolyte
screening because they miss:

- concentration effects
- ion pairing
- mixed-solvent competition
- interface polarization
- decomposition chemistry after charge transfer

So the interface keeps an explicit `effective_shift_kjmol` hook. That lets you:

- start from DFT instead of guessing
- then calibrate toward experiment or cluster/interface calculations

## Example Usage

Validation-only pass:

```bash
python generate_redox_config_from_dft.py \
  --workflow dft_workflows/example_ec_dmc_reduction.json \
  --output-dir results/dft_interface_validation \
  --allow-missing-pyscf
```

Real DFT pass after installing PySCF:

```bash
python generate_redox_config_from_dft.py \
  --workflow dft_workflows/example_ec_dmc_reduction.json \
  --output-dir results/dft_interface_ec_dmc
```

GPU4PySCF pass:

```bash
python generate_redox_config_from_dft.py \
  --workflow dft_workflows/example_ec_dmc_reduction_gpu.json \
  --output-dir results/dft_interface_ec_dmc_gpu
```

## Input Geometry

Currently supported:

- `*.pdb`
- `*.xyz`

This is enough to hook directly into the existing `pdb_bank/` directory.

## Current Limits

- single-point DFT only
- no conformer search yet
- no thermal free-energy correction yet
- no explicit cluster sampling yet
- no automatic formulation-level concentration correction yet
- GPU path currently expects a PCM-family solvent model, not `ddCOSMO`

## Recommended Workflow

1. Use this interface to get DFT-anchored offsets.
2. Compare the implied onset or `E1/2` against experiment.
3. Apply small empirical shifts where needed.
4. Feed the generated config into the batch-screening workflow.
