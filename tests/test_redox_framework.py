#!/usr/bin/env python
"""
test_redox_framework.py
-----------------------
Unit tests for redox MC and voltage sweep framework (Tasks 3, 4, 5).

Quick validation that all modules load and basic functionality works.

Usage:
    python tests/test_redox_framework.py
"""

import sys
import json
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = REPO_ROOT / "core"
if str(CORE_DIR) not in sys.path:
    sys.path.insert(0, str(CORE_DIR))

def test_imports():
    """Test that all modules can be imported."""
    print("[TEST 1] Testing imports...")
    try:
        from redox_mc import RedoxMC, RedoxParameters, load_redox_parameters_from_json
        print("  ✓ redox_mc module imported successfully")
    except Exception as e:
        print(f"  ✗ FAILED to import redox_mc: {e}")
        return False

    try:
        from voltage_sweep import VoltageSweepSimulation, RedoxPotentialAnalyzer
        print("  ✓ voltage_sweep module imported successfully")
    except Exception as e:
        print(f"  ✗ FAILED to import voltage_sweep: {e}")
        return False

    return True


def test_redox_parameters():
    """Test RedoxParameters creation."""
    print("\n[TEST 2] Testing RedoxParameters creation...")
    try:
        from redox_mc import RedoxParameters

        # Create a simple redox-active molecule (like benzoquinone)
        bq = RedoxParameters(
            name="BQ",
            residue_names=["BQ", "BQU"],
            charge_states={0: 0.0, -1: -1.0, -2: -2.0},
            electron_affinities={(0, -1): 1.53, (-1, -2): -0.2},
            allowed_states=[0, -1, -2],
            partition_coefficient=-2.5,
        )

        assert bq.name == "BQ"
        assert len(bq.allowed_states) == 3
        assert bq.electron_affinities[(0, -1)] == 1.53

        print(f"  ✓ Created RedoxParameters for {bq.name}")
        print(f"    - Allowed states: {bq.allowed_states}")
        print(f"    - EA(0→-1): {bq.electron_affinities[(0, -1)]} eV")
        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        return False


def test_redox_mc():
    """Test RedoxMC class initialization and basic calculations."""
    print("\n[TEST 3] Testing RedoxMC calculations...")
    try:
        from redox_mc import RedoxMC, RedoxParameters

        # Create a test molecule
        bq = RedoxParameters(
            name="BQ",
            residue_names=["BQ"],
            charge_states={0: 0.0, -1: -1.0},
            electron_affinities={(0, -1): 1.53},
            allowed_states=[0, -1],
            partition_coefficient=-1.0,
        )

        # Create MC sampler
        mc = RedoxMC(
            redox_params=bq,
            temperature_k=300.0,
            voltage_v=-3.5,
        )

        # Test energy calculations
        e_state_0 = mc.compute_state_energy(0, pe_nonbonded=0.0)
        e_state_m1 = mc.compute_state_energy(-1, pe_nonbonded=0.0)

        print(f"  ✓ RedoxMC initialized")
        print(f"    - Temperature: {mc.temperature_k} K")
        print(f"    - Voltage: {mc.voltage_v} V")
        print(f"    - kT: {mc.kT:.3f} kJ/mol")
        print(f"    - E(state 0): {e_state_0:.3f} kJ/mol")
        print(f"    - E(state -1): {e_state_m1:.3f} kJ/mol")
        print(f"    - ΔE: {(e_state_m1 - e_state_0):.3f} kJ/mol")
        assert e_state_m1 != e_state_0

        # Test MC sampling
        summary = mc.run_mc_sampling(n_steps=100)
        occupancy = mc.get_state_occupancy()

        acc_rate = mc.get_acceptance_rate()
        print(f"  ✓ Completed 100 MC moves")
        print(f"    - Acceptance rate: {acc_rate:.3f}")
        print(f"    - Occupancy: {occupancy}")
        assert summary["state_history_length"] == 100
        assert sum(occupancy.values()) > 0.99

        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_load_config():
    """Test loading redox configuration from JSON."""
    print("\n[TEST 4] Testing config JSON loading...")
    try:
        from redox_mc import load_redox_parameters_from_json

        config_file = REPO_ROOT / "configs" / "redox" / "config_redox.json"
        if not config_file.exists():
            print(f"  ⊘ Skipped (config file not found at {config_file})")
            return True

        redox_params = load_redox_parameters_from_json(str(config_file))

        print(f"  ✓ Loaded redox parameters from {config_file}")
        print(f"    - Molecules: {list(redox_params.keys())}")

        for mol_name, params in redox_params.items():
            print(f"      {mol_name}: states {params.allowed_states}, "
                  f"residues {params.residue_names}")

        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_redox_potential_analyzer():
    """Test RedoxPotentialAnalyzer."""
    print("\n[TEST 5] Testing RedoxPotentialAnalyzer...")
    try:
        import pandas as pd
        from voltage_sweep import RedoxPotentialAnalyzer

        # Create mock sweep results
        voltages = np.linspace(-4.5, -2.5, 21)
        # Sigmoidal occupancy curve (fraction reduced)
        fraction_reduced = 1.0 / (1.0 + np.exp(20 * (voltages + 3.5)))

        df = pd.DataFrame({
            'voltage_v': voltages,
            'EC_fraction_reduced': fraction_reduced,
        })

        # Find half-wave potential
        e_half = RedoxPotentialAnalyzer.find_half_wave_potential(
            df, 'EC', voltage_col='voltage_v'
        )

        print(f"  ✓ RedoxPotentialAnalyzer working")
        print(f"    - Voltage range: {voltages[0]:.2f} to {voltages[-1]:.2f} V")
        print(f"    - E_1/2 (half-wave): {e_half:.3f} V vs vacuum")

        # Convert to Li reference
        e_half_li = RedoxPotentialAnalyzer.convert_to_li_reference(e_half)
        print(f"    - E_1/2 (vs Li/Li+): {e_half_li:.3f} V")

        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("="*70)
    print("REDOX FRAMEWORK UNIT TESTS (Tasks 3, 4, 5)")
    print("="*70)

    tests = [
        test_imports,
        test_redox_parameters,
        test_redox_mc,
        test_load_config,
        test_redox_potential_analyzer,
    ]

    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"\n[FATAL ERROR in {test_func.__name__}]: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)

    # Summary
    print("\n" + "="*70)
    n_passed = sum(results)
    n_total = len(results)
    print(f"SUMMARY: {n_passed}/{n_total} tests passed")

    if all(results):
        print("✓ ALL TESTS PASSED - Framework is ready to use!")
        print("\nNext steps:")
        print("  1. Edit configs/redox/config_redox.json with your redox-active molecules")
        print("  2. Run: python scripts/cli/run_voltage_sweep.py ...")
        return 0
    else:
        print("✗ SOME TESTS FAILED - Please check errors above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
