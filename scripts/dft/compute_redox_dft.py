#!/usr/bin/env python3
"""
Compute DFT redox energies for all available molecules.
B3LYP/6-31+g(d,p) + PCM(ε=40) single-point vertical EA.

Run with:
    /opt/homebrew/Caskroom/miniforge/base/bin/python3 compute_redox_dft.py
"""

import sys, json, time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
STRUCTURES = REPO_ROOT / "structures" / "pdb_bank"
OUT_JSON   = REPO_ROOT / "results" / "dft_redox_all_20260324.json"

sys.path.insert(0, str(REPO_ROOT / "core"))

import numpy as np
import pyscf
from pyscf import gto, dft
from pyscf.solvent import ddcosmo

print(f"PySCF {pyscf.__version__}  |  structures: {STRUCTURES}")

Ha_to_eV    = 27.211386
VACUUM_TO_LI = 1.4
F_KJMOL_PER_V = 96.485

# ── molecules to compute ──────────────────────────────────────────────────────
# neutral species: compute 0 (neutral) → -1 (radical anion)
# salt anions:     compute 0 (neutral anion) → -1 (radical anion, spin=1)
# VC: no PDB available — skip

MOLECULES = {
    # ── solvents (neutral) ──────────────────────────────────────────────────
    "EC":   {"pdb": "EC.pdb",   "cation": None},
    "DMC":  {"pdb": "DMC.pdb",  "cation": None},
    "DEC":  {"pdb": "DEC.pdb",  "cation": None},
    "EMC":  {"pdb": "EMC.pdb",  "cation": None},
    "DME":  {"pdb": "DME.pdb",  "cation": None},
    "PC":   {"pdb": "PC.pdb",   "cation": None},
    # ── salts (anion only; Li+ spectator) ─────────────────────────────────
    # Salts: neutral = anion at its natural charge (-1)
    #        radical  = reduced form (-2, spin=1)
    "LiFSI":   {"pdb": "FSI.pdb",  "cation": None, "neutral_chg": -1},   # FSI- → FSI•2-
    "LiPF6":   {"pdb": "PF6.pdb",  "cation": None, "neutral_chg": -1},   # PF6- → PF6•2-
    "LiTFSI":  {"pdb": "TFSI.pdb", "cation": None, "neutral_chg": -1},  # TFSI- → TFSI•2-
    "LiBOB":   {"pdb": "BOB.pdb",  "cation": None, "neutral_chg": -1},  # BOB- → BOB•2-
    "LiDFP":   {"pdb": "DFP.pdb",  "cation": None, "neutral_chg": -1},  # DFP- → DFP•2-
    # ── additives (neutral) ────────────────────────────────────────────────
    "FEC": {"pdb": "FEC.pdb", "cation": None},
    "PS":  {"pdb": "PS.pdb",  "cation": None},
    # ── others (may be useful) ──────────────────────────────────────────────
    "GBL":  {"pdb": "GBL.pdb",  "cation": None},  # gamma-butyrolactone
    "BF4":  {"pdb": "BF4.pdb",  "cation": None, "neutral_chg": -1},  # BF4- anion
    "NO3":  {"pdb": "NO3.pdb",  "cation": None, "neutral_chg": -1},  # NO3- anion
}


def clean_pdb(pdb_path: Path) -> str:
    """Extract ATOM/HETATM lines, return a simple element+xyz string."""
    lines = []
    for line in open(pdb_path):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        elem = line[76:78].strip()
        if not elem:
            elem = line[12:14].strip().capitalize()  # fallback
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        lines.append(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    return "".join(lines)


def sp(pdb_path: Path, charge: int, spin: int, eps: float = 40.0) -> float:
    """Single-point PCM-DFT energy in Hartree."""
    mol = gto.Mole()
    mol.atom = clean_pdb(pdb_path)
    mol.basis = "6-31+g(d,p)"
    mol.charge = charge
    mol.spin = spin
    mol.max_memory = 4000
    mol.verbose = 0
    mol.build()

    mf = dft.UKS(mol) if spin else dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.grids.level = 3
    mf.with_solvent = ddcosmo.DDCOSMO(mol)
    mf.with_solvent.eps = eps

    try:
        return float(mf.kernel())
    except Exception:
        # Fallback: RKS with smaller grid
        mf2 = dft.RKS(mol) if not spin else dft.UKS(mol)
        mf2.xc = "b3lyp"
        mf2.grids.level = 2
        mf2.with_solvent = ddcosmo.DDCOSMO(mol)
        mf2.with_solvent.eps = eps
        mf2.dm = mf2.init_guess_by_1e()
        return float(mf2.kernel())


def compute_ea(pdb_name: str, label: str, cation: str, neutral_chg: int = 0):
    """
    Compute vertical EA for a molecule.
    For salts (neutral_chg=-1): computes anion → radical-anion.
    For neutral molecules (neutral_chg=0): computes neutral → radical-anion.
    Returns dict with all computed values.
    """
    pdb = STRUCTURES / pdb_name
    if not pdb.exists():
        return {"status": "SKIP", "note": f"PDB not found: {pdb}"}

    radical_chg = neutral_chg - 1    # one extra electron
    radical_spin = 1 if neutral_chg == 0 else 1   # always spin=1 for radical

    t0 = time.time()
    try:
        E_neutral = sp(pdb, neutral_chg, 0)      # closed-shell neutral/anion
        E_anion   = sp(pdb, radical_chg, radical_spin)   # radical (spin=1)
    except Exception as ex:
        return {"status": "ERROR", "note": str(ex)}

    elapsed = time.time() - t0

    dG_ha   = E_anion - E_neutral          # Ha (positive = anion less stable)
    EA_eV   = -dG_ha * Ha_to_eV           # electron affinity (eV)
    V_vs_Li = EA_eV + VACUUM_TO_LI        # vs Li/Li+
    offset_kjmol = EA_eV * F_KJMOL_PER_V  # RedoxMC offset in kJ/mol

    result = {
        "status":         "OK",
        "pdb":            pdb_name,
        "E_neutral_ha":   E_neutral,
        "E_anion_ha":     E_anion,
        "dG_ha":          dG_ha,
        "EA_eV":          EA_eV,
        "E_vs_Li":        V_vs_Li,
        "offset_kjmol":   offset_kjmol,
        "elapsed_s":      round(elapsed, 1),
    }
    print(
        f"  {label:<10}  E_neutral={E_neutral:.6f}  E_anion={E_anion:.6f}  "
        f"EA={EA_eV:+.4f} eV  V_Li={V_vs_Li:.4f} V  [{elapsed:.0f}s]"
    )
    return result


# ── Main ──────────────────────────────────────────────────────────────────────
print(f"\nComputing vertical EA for {len(MOLECULES)} molecules...\n")

all_results = {}
for name, spec in MOLECULES.items():
    label = name
    print(f"[{name}]")
    res = compute_ea(spec["pdb"], label, spec["cation"], neutral_chg=spec.get("neutral_chg", 0))
    all_results[name] = res

# ── Save raw results ───────────────────────────────────────────────────────────
with open(OUT_JSON, "w") as f:
    json.dump(all_results, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
print(f"\nRaw results → {OUT_JSON}")

# ── Print summary table ────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("  DFT VERTICAL EA SUMMARY  (B3LYP/6-31+g(d,p) PCM eps=40)")
print("=" * 65)
print(f"  {'Molecule':<12} {'EA (eV)':>10} {'V vs Li/Li+':>14}  {'Status':>8}")
print("-" * 65)

ok_results = {}
for name, r in all_results.items():
    if r["status"] == "OK":
        print(f"  {name:<12} {r['EA_eV']:>+10.4f} {r['E_vs_Li']:>+14.4f}  {'OK':>8}")
        ok_results[name] = r
    else:
        print(f"  {name:<12} {'N/A':>10} {'N/A':>14}  {r['status']:>8}")

# ── Save corrected RedoxMC config ───────────────────────────────────────────────
REDOX_CONFIG = REPO_ROOT / "results" / "redox_config_dft_corrected.json"
config = {
    "metadata": {
        "method": "B3LYP/6-31+g(d,p)",
        "solvent": "PCM (ddCOSMO, eps=40)",
        "vacuum_to_li_offset_ev": VACUUM_TO_LI,
        "temperature_k": 300.0,
        "note": "Vertical EA. Geometry optimization NOT applied — use with caution for absolute values.",
    },
    "molecules": {},
}
for name, r in ok_results.items():
    config["molecules"][name] = {
        "residue_names": [name],
        "charge_states": {"0": 0.0, "-1": -1.0},
        "electron_affinities": {"0,-1": float(r["EA_eV"])},
        "state_free_energy_offsets_kjmol": {"0": 0.0, "-1": round(r["offset_kjmol"], 4)},
        "allowed_states": [0, -1],
        "partition_coefficient": 0.0,
        "_dft_E_vs_Li": round(r["E_vs_Li"], 4),
    }

with open(REDOX_CONFIG, "w") as f:
    json.dump(config, f, indent=2)
print(f"\nRedoxMC config → {REDOX_CONFIG}")
print("\nDone.")
