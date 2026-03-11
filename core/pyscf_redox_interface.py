"""
PySCF-backed DFT interface for generating screening parameters used by RedoxMC.

This module does two jobs:
1. Run optional single-point DFT calculations for multiple molecular charge states.
2. Map the resulting state free-energy differences onto the project's
   `state_free_energy_offsets_kjmol` representation.

The intent is not to claim that raw gas-phase DFT is already a formulation-level
electrolyte prediction. Instead, this gives a reproducible bridge:

    geometry -> DFT state energies -> screening offsets -> RedoxMC config

Optional empirical shifts can then be applied on top of the DFT values to align
with experimental onset data or more complete cluster/interface calculations.
"""

from __future__ import annotations

import importlib.util
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


HARTREE_TO_KJMOL = 2625.499638
FARADAY_KJMOL_PER_V = 96.48533212331002


@dataclass
class ChargeStateSpec:
    state: int
    charge: int
    spin: int
    label: str = ""


@dataclass
class DFTStateResult:
    state: int
    charge: int
    spin: int
    total_energy_hartree: float
    relative_energy_kjmol: float
    effective_offset_kjmol: float
    predicted_half_wave_v_vs_vacuum: Optional[float]
    predicted_half_wave_v_vs_li: Optional[float]


def pyscf_available() -> bool:
    return importlib.util.find_spec("pyscf") is not None


def gpu4pyscf_available() -> bool:
    return importlib.util.find_spec("gpu4pyscf") is not None


def infer_geometry_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".xyz":
        return "xyz"
    if suffix == ".pdb":
        return "pdb"
    raise ValueError(f"Unsupported geometry format for {path}")


def load_xyz(path: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    with path.open() as f:
        lines = [line.rstrip() for line in f if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"XYZ file {path} is too short")

    try:
        natoms = int(lines[0])
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom count in {path}: {lines[0]!r}") from exc

    atoms = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom line in {path}: {line!r}")
        symbol = parts[0]
        coords = tuple(float(value) for value in parts[1:4])
        atoms.append((symbol, coords))
    if len(atoms) != natoms:
        raise ValueError(f"Expected {natoms} atoms in {path}, found {len(atoms)}")
    return atoms


def load_pdb(path: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    atoms = []
    with path.open() as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip()
            atom_name = line[12:16].strip()
            if not element:
                element = "".join(ch for ch in atom_name if ch.isalpha())[:2].strip().title()
            if not element:
                raise ValueError(f"Failed to infer element from PDB line: {line.rstrip()}")
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((element, (x, y, z)))
    if not atoms:
        raise ValueError(f"No atoms found in PDB file {path}")
    return atoms


def load_geometry(path: str, geometry_format: Optional[str] = None) -> List[Tuple[str, Tuple[float, float, float]]]:
    geometry_path = Path(path).resolve()
    fmt = geometry_format or infer_geometry_format(geometry_path)
    if fmt == "xyz":
        return load_xyz(geometry_path)
    if fmt == "pdb":
        return load_pdb(geometry_path)
    raise ValueError(f"Unsupported geometry format: {fmt}")


def format_atom_block(atoms: Iterable[Tuple[str, Tuple[float, float, float]]]) -> str:
    return "\n".join(
        f"{symbol} {x:.10f} {y:.10f} {z:.10f}"
        for symbol, (x, y, z) in atoms
    )


def _require_pyscf():
    if not pyscf_available():
        raise ImportError(
            "PySCF is not installed. Install it in the active environment, for example: "
            "`pip install pyscf` or `conda install -c conda-forge pyscf`."
        )
    from pyscf import dft, gto
    from pyscf import solvent as pyscf_solvent

    return gto, dft, pyscf_solvent


def apply_solvent_model(mf, solvent: Optional[Dict[str, Any]], backend: str):
    if not solvent:
        return mf

    model = solvent.get("model", "ddcosmo").lower()
    eps = float(solvent.get("eps", 1.0))

    if model == "ddcosmo":
        if backend == "gpu4pyscf":
            raise ValueError("gpu4pyscf workflow does not use ddCOSMO in this interface; use a PCM model instead")
        _, _, pyscf_solvent = _require_pyscf()
        mf = pyscf_solvent.ddCOSMO(mf)
        mf.with_solvent.eps = eps
        return mf

    if model in {"cpcm", "iefpcm", "cosmo", "ss(v)pe", "ssvpe"}:
        pcm_method = {
            "cpcm": "C-PCM",
            "iefpcm": "IEF-PCM",
            "cosmo": "COSMO",
            "ss(v)pe": "SS(V)PE",
            "ssvpe": "SS(V)PE",
        }[model]
        mf = mf.PCM()
        mf.with_solvent.method = pcm_method
        mf.with_solvent.eps = eps
        return mf

    raise ValueError(f"Unsupported solvent model: {model}")


def run_pyscf_single_point(
    atoms: List[Tuple[str, Tuple[float, float, float]]],
    charge: int,
    spin: int,
    method: str,
    basis: str,
    solvent: Optional[Dict[str, Any]] = None,
    backend: str = "cpu",
    use_density_fit: bool = False,
    max_cycle: int = 100,
    conv_tol: float = 1e-9,
) -> float:
    gto, dft, _ = _require_pyscf()

    mol = gto.M(
        atom=format_atom_block(atoms),
        basis=basis,
        charge=charge,
        spin=spin,
        unit="Angstrom",
        verbose=0,
    )

    if spin == 0:
        mf = dft.RKS(mol)
    else:
        mf = dft.UKS(mol)
    mf.xc = method
    mf.max_cycle = max_cycle
    mf.conv_tol = conv_tol

    backend = backend.lower()
    if use_density_fit and hasattr(mf, "density_fit"):
        mf = mf.density_fit()

    if backend == "gpu4pyscf":
        if not gpu4pyscf_available():
            raise ImportError(
                "gpu4pyscf is not installed. Install the CUDA-matched package, e.g. "
                "`pip install gpu4pyscf-cuda12x cutensor-cu12` for CUDA 12.x."
            )
        if not hasattr(mf, "to_gpu"):
            raise RuntimeError("Current PySCF object does not expose to_gpu(); upgrade PySCF/gpu4pyscf")
        mf = mf.to_gpu()
    elif backend != "cpu":
        raise ValueError(f"Unsupported backend: {backend}")

    mf = apply_solvent_model(mf, solvent, backend)

    energy = mf.kernel()
    if not mf.converged:
        raise RuntimeError(
            f"PySCF SCF did not converge for charge={charge}, spin={spin}, method={method}, basis={basis}"
        )
    return float(energy)


def compute_state_profiles(
    atoms: List[Tuple[str, Tuple[float, float, float]]],
    states: List[ChargeStateSpec],
    method: str,
    basis: str,
    solvent: Optional[Dict[str, Any]],
    backend: str,
    use_density_fit: bool,
    reference_state: int,
    vacuum_to_li_offset_ev: float,
    effective_shift_kjmol: Optional[Dict[int, float]] = None,
) -> List[DFTStateResult]:
    effective_shift_kjmol = effective_shift_kjmol or {}
    energies_hartree: Dict[int, float] = {}
    charge_by_state = {spec.state: spec.charge for spec in states}

    for spec in states:
        energies_hartree[spec.state] = run_pyscf_single_point(
            atoms=atoms,
            charge=spec.charge,
            spin=spec.spin,
            method=method,
            basis=basis,
            solvent=solvent,
            backend=backend,
            use_density_fit=use_density_fit,
        )

    if reference_state not in energies_hartree:
        raise ValueError(f"Reference state {reference_state} not present in charge-state spec")

    e_ref = energies_hartree[reference_state]
    q_ref = charge_by_state[reference_state]
    results: List[DFTStateResult] = []

    for spec in sorted(states, key=lambda item: item.state):
        rel_kjmol = (energies_hartree[spec.state] - e_ref) * HARTREE_TO_KJMOL
        effective_offset_kjmol = rel_kjmol + float(effective_shift_kjmol.get(spec.state, 0.0))
        delta_charge = spec.charge - q_ref
        if delta_charge == 0:
            v_vac = None
            v_li = None
        else:
            v_vac = effective_offset_kjmol / (delta_charge * FARADAY_KJMOL_PER_V)
            v_li = v_vac + vacuum_to_li_offset_ev
        results.append(
            DFTStateResult(
                state=spec.state,
                charge=spec.charge,
                spin=spec.spin,
                total_energy_hartree=energies_hartree[spec.state],
                relative_energy_kjmol=rel_kjmol,
                effective_offset_kjmol=effective_offset_kjmol,
                predicted_half_wave_v_vs_vacuum=v_vac,
                predicted_half_wave_v_vs_li=v_li,
            )
        )

    return results


def build_redox_config_entry(
    molecule_name: str,
    residue_names: List[str],
    states: List[ChargeStateSpec],
    state_results: List[DFTStateResult],
    partition_coefficient: float = 0.0,
) -> Dict[str, Any]:
    result_by_state = {result.state: result for result in state_results}
    charge_states = {str(spec.state): float(spec.charge) for spec in states}
    offsets = {
        str(spec.state): float(result_by_state[spec.state].effective_offset_kjmol)
        for spec in states
    }
    reference_state = min(
        (spec.state for spec in states if spec.charge == 0),
        default=min(spec.state for spec in states),
    )

    electron_affinities = {}
    ordered_states = sorted(states, key=lambda item: item.charge, reverse=True)
    for first, second in zip(ordered_states, ordered_states[1:]):
        rel_kjmol = (
            result_by_state[second.state].relative_energy_kjmol
            - result_by_state[first.state].relative_energy_kjmol
        )
        electron_affinities[f"{first.state},{second.state}"] = -rel_kjmol / FARADAY_KJMOL_PER_V

    if str(reference_state) in offsets:
        offsets[str(reference_state)] = 0.0

    return {
        "residue_names": residue_names,
        "charge_states": charge_states,
        "electron_affinities": electron_affinities,
        "state_free_energy_offsets_kjmol": offsets,
        "allowed_states": [spec.state for spec in sorted(states, key=lambda item: item.state)],
        "partition_coefficient": partition_coefficient,
    }


def build_screening_config_document(
    metadata: Dict[str, Any],
    molecule_entries: Dict[str, Dict[str, Any]],
    comments: Optional[List[str]] = None,
) -> Dict[str, Any]:
    doc: Dict[str, Any] = {}
    if comments:
        doc["_comment"] = comments
    doc["metadata"] = metadata
    doc["molecules"] = molecule_entries
    return doc


def serialize_state_results(state_results: List[DFTStateResult]) -> List[Dict[str, Any]]:
    return [
        {
            "state": item.state,
            "charge": item.charge,
            "spin": item.spin,
            "total_energy_hartree": item.total_energy_hartree,
            "relative_energy_kjmol": item.relative_energy_kjmol,
            "effective_offset_kjmol": item.effective_offset_kjmol,
            "predicted_half_wave_v_vs_vacuum": item.predicted_half_wave_v_vs_vacuum,
            "predicted_half_wave_v_vs_li": item.predicted_half_wave_v_vs_li,
        }
        for item in state_results
    ]


def load_workflow_spec(path: str) -> Dict[str, Any]:
    with Path(path).resolve().open() as f:
        return json.load(f)
