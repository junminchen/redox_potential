#!/usr/bin/env python
"""
Generate RedoxMC screening configs from PySCF DFT calculations.

Input: a workflow JSON describing molecules, geometry files, charge states,
DFT settings, and optional empirical shifts.

Output:
- `dft_raw_results.json`: per-state raw energies and mapped offsets
- `redox_config_generated.json`: screening config ready for RedoxMC
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_DIR = REPO_ROOT / "core"
if str(CORE_DIR) not in sys.path:
    sys.path.insert(0, str(CORE_DIR))

from pyscf_redox_interface import (
    ChargeStateSpec,
    build_redox_config_entry,
    build_screening_config_document,
    compute_state_profiles,
    load_geometry,
    load_workflow_spec,
    pyscf_available,
    serialize_state_results,
)


def molecule_charge_specs(spec: dict) -> list[ChargeStateSpec]:
    items = []
    for state_cfg in spec.get("charge_states", []):
        items.append(
            ChargeStateSpec(
                state=int(state_cfg["state"]),
                charge=int(state_cfg["charge"]),
                spin=int(state_cfg.get("spin", 0)),
                label=str(state_cfg.get("label", "")),
            )
        )
    if not items:
        raise ValueError(f"Molecule {spec.get('name', '<unknown>')} defines no charge states")
    return items


def resolve_input_path(workflow_path: Path, raw_path: str) -> Path:
    candidate = Path(raw_path)
    if candidate.is_absolute():
        return candidate.resolve()

    workflow_relative = (workflow_path.parent / candidate).resolve()
    if workflow_relative.exists():
        return workflow_relative

    cwd_relative = (Path.cwd() / candidate).resolve()
    if cwd_relative.exists():
        return cwd_relative

    repo_relative = (REPO_ROOT / candidate).resolve()
    if repo_relative.exists():
        return repo_relative

    return workflow_relative


def main():
    parser = argparse.ArgumentParser(description="Generate RedoxMC config from PySCF DFT data")
    parser.add_argument("--workflow", required=True, help="Workflow JSON for PySCF redox interface")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--allow-missing-pyscf",
        action="store_true",
        help="Write geometry/metadata validation errors only; do not fail if PySCF is not installed",
    )
    args = parser.parse_args()

    workflow_path = Path(args.workflow).resolve()
    workflow = load_workflow_spec(str(workflow_path))
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = workflow.get("metadata", {})
    vacuum_to_li_offset_ev = float(metadata.get("vacuum_to_li_offset_ev", 1.4))
    default_dft = workflow.get("default_dft", {})

    if not pyscf_available() and not args.allow_missing_pyscf:
        raise SystemExit(
            "PySCF is not installed in the active environment. "
            "Re-run after installing it, or pass --allow-missing-pyscf for a dry validation pass."
        )

    raw_output = {
        "metadata": metadata,
        "workflow_file": str(workflow_path),
        "pyscf_available": pyscf_available(),
        "molecules": {},
    }
    config_entries = {}

    for mol_cfg in workflow.get("molecules", []):
        name = mol_cfg["name"]
        geometry_file = resolve_input_path(workflow_path, mol_cfg["geometry_file"])
        geometry_format = mol_cfg.get("geometry_format")
        atoms = load_geometry(str(geometry_file), geometry_format)

        charge_specs = molecule_charge_specs(mol_cfg)
        reference_state = int(mol_cfg.get("reference_state", 0))
        dft_settings = dict(default_dft)
        dft_settings.update(mol_cfg.get("dft", {}))
        method = dft_settings.get("method", "b3lyp")
        basis = dft_settings.get("basis", "6-31+g(d,p)")
        solvent = dft_settings.get("solvent")
        backend = dft_settings.get("backend", "cpu")
        use_density_fit = bool(dft_settings.get("use_density_fit", backend == "gpu4pyscf"))
        effective_shift = {
            int(state): float(value)
            for state, value in mol_cfg.get("effective_shift_kjmol", {}).items()
        }

        if pyscf_available():
            state_results = compute_state_profiles(
                atoms=atoms,
                states=charge_specs,
                method=method,
                basis=basis,
                solvent=solvent,
                backend=backend,
                use_density_fit=use_density_fit,
                reference_state=reference_state,
                vacuum_to_li_offset_ev=vacuum_to_li_offset_ev,
                effective_shift_kjmol=effective_shift,
            )
            config_entries[name] = build_redox_config_entry(
                molecule_name=name,
                residue_names=mol_cfg.get("residue_names", [name]),
                states=charge_specs,
                state_results=state_results,
                partition_coefficient=float(mol_cfg.get("partition_coefficient", 0.0)),
            )
            raw_output["molecules"][name] = {
                "geometry_file": str(geometry_file),
                "geometry_format": geometry_format or geometry_file.suffix.lstrip(".").lower(),
                "dft": dft_settings,
                "reference_state": reference_state,
                "residue_names": mol_cfg.get("residue_names", [name]),
                "state_results": serialize_state_results(state_results),
            }
        else:
            raw_output["molecules"][name] = {
                "geometry_file": str(geometry_file),
                "geometry_format": geometry_format or geometry_file.suffix.lstrip(".").lower(),
                "dft": dft_settings,
                "reference_state": reference_state,
                "residue_names": mol_cfg.get("residue_names", [name]),
                "state_results": [],
                "note": "PySCF unavailable; validation-only pass completed.",
            }

    comments = [
        "Generated from PySCF workflow input.",
        "Raw DFT-derived offsets may still require empirical or cluster/interface calibration.",
    ]
    screening_config = build_screening_config_document(
        metadata=metadata,
        molecule_entries=config_entries,
        comments=comments,
    )

    with (output_dir / "dft_raw_results.json").open("w") as f:
        json.dump(raw_output, f, indent=2)
    with (output_dir / "redox_config_generated.json").open("w") as f:
        json.dump(screening_config, f, indent=2)

    print(output_dir / "dft_raw_results.json")
    print(output_dir / "redox_config_generated.json")


if __name__ == "__main__":
    main()
