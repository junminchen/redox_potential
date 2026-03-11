#!/usr/bin/env python
"""
run_voltage_sweep.py
--------------------
Driver script for voltage sweep simulations in 1M LiPF6 EC/DMC electrolyte.

Performs constant-potential MD across a range of electrode potentials,
extracting redox molecule charge state distributions and simulating CV/LSV curves.

Usage:
    python run_voltage_sweep.py --config config_opls.json \
                                 --redox-config config_redox.json \
                                 --pdb start_with_electrodes.pdb \
                                 --output-dir sweep_results_ec_dmc
"""

import argparse
import json
from datetime import datetime
import sys

from voltage_sweep import VoltageSweepSimulation, RedoxPotentialAnalyzer
from redox_mc import load_redox_parameters_from_json


def main():
    parser = argparse.ArgumentParser(
        description="Voltage sweep simulation for redox potential mapping"
    )
    parser.add_argument("--config", required=True, help="OPLS MD config JSON")
    parser.add_argument(
        "--redox-config", required=True, help="Redox parameters config JSON"
    )
    parser.add_argument("--pdb", required=True, help="Input PDB (electrodes + electrolyte)")
    parser.add_argument(
        "--output-dir", default="sweep_results", help="Output directory"
    )
    parser.add_argument(
        "--v-start", type=float, default=-4.5, help="Start voltage (V vs vacuum)"
    )
    parser.add_argument(
        "--v-end", type=float, default=-2.5, help="End voltage (V vs vacuum)"
    )
    parser.add_argument(
        "--v-step", type=float, default=0.1, help="Voltage step (V)"
    )
    parser.add_argument(
        "--equil-steps", type=int, default=5000, help="Equilibration steps per voltage"
    )
    parser.add_argument(
        "--sample-steps", type=int, default=10000, help="Sampling steps per voltage"
    )
    parser.add_argument(
        "--platform", help="OpenMM platform (CUDA, OpenCL, CPU). Defaults to config file."
    )
    parser.add_argument(
        "--lsv-peak-v", type=float, help="Experimental LSV peak (V vs Li/Li+) for comparison"
    )

    args = parser.parse_args()

    print("="*80)
    print("VOLTAGE SWEEP SIMULATION FOR 1M LiPF6 EC/DMC ELECTROLYTE")
    print("="*80)
    print(f"Start time: {datetime.now()}")
    print(f"PDB file: {args.pdb}")
    print(f"MD config: {args.config}")
    print(f"Redox config: {args.redox_config}")
    print(f"Output dir: {args.output_dir}")
    print(f"Voltage range: {args.v_start:.2f} to {args.v_end:.2f} V (step {args.v_step:.2f} V)")
    print(f"Platform override: {args.platform if args.platform is not None else 'config file default'}")
    print("="*80)

    # Load and display redox parameters
    try:
        redox_params = load_redox_parameters_from_json(args.redox_config)
        print(f"\n[INFO] Loaded redox parameters for {len(redox_params)} molecule(s):")
        for mol_name, params in redox_params.items():
            print(f"  - {mol_name}: states {params.allowed_states}")
            for (s1, s2), ea in params.electron_affinities.items():
                print(f"      EA({s1}→{s2}) = {ea:.2f} eV")
    except Exception as e:
        print(f"[ERROR] Failed to load redox config: {e}")
        sys.exit(1)

    # Initialize and run voltage sweep
    try:
        sweep = VoltageSweepSimulation(
            pdb_file=args.pdb,
            config_file=args.config,
            redox_config_file=args.redox_config,
            output_dir=args.output_dir,
        )

        print("\n[INFO] Starting voltage sweep simulation...")
        df_results = sweep.run_sweep(
            v_start=args.v_start,
            v_end=args.v_end,
            v_step=args.v_step,
            equil_steps_per_v=args.equil_steps,
            sample_steps_per_v=args.sample_steps,
            platform=args.platform,
        )

        print(f"\n[SUCCESS] Voltage sweep completed!")
        print(f"Collected {len(df_results)} voltage points")
        print(f"\nResults summary:")
        print(df_results[["voltage_v", "q_cathode_mean", "q_anode_mean"]].to_string(index=False))

        # Optional: compare with experimental LSV if provided
        if args.lsv_peak_v is not None:
            print(f"\n[INFO] Comparing with experimental LSV peak at {args.lsv_peak_v:.3f} V vs Li/Li+:")
            for mol_name in redox_params.keys():
                try:
                    comparison = RedoxPotentialAnalyzer.compare_with_lsv(
                        df_results,
                        mol_name,
                        args.lsv_peak_v,
                    )
                    print(f"\n  {mol_name}:")
                    print(f"    Simulated E_1/2 (vacuum):  {comparison['e_half_wave_v_vs_vacuum']:.3f} V")
                    print(f"    Simulated E_1/2 (Li/Li+):  {comparison['e_half_wave_v_vs_li']:.3f} V")
                    print(f"    LSV peak (Li/Li+):         {comparison['lsv_peak_v_vs_li']:.3f} V")
                    print(f"    Difference:                {comparison['difference_mv']:.1f} mV")
                    if comparison["agreement_excellent"]:
                        print(f"    Agreement: EXCELLENT (< 50 mV) ✓")
                    elif comparison["agreement_good"]:
                        print(f"    Agreement: GOOD (< 100 mV) ✓")
                    else:
                        print(f"    Agreement: POOR (> 100 mV)")
                except Exception as e:
                    print(f"    Warning: Could not analyze {mol_name}: {e}")

        print(f"\nEnd time: {datetime.now()}")
        print("="*80)

    except Exception as e:
        print(f"\n[FATAL ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
