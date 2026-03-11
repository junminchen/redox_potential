#!/usr/bin/env python
"""
Generate experiment-friendly LSV/CV-like plots and a markdown report from a
redox configuration.
"""

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_DIR = REPO_ROOT / "core"
if str(CORE_DIR) not in sys.path:
    sys.path.insert(0, str(CORE_DIR))

from redox_mc import RedoxMC, load_redox_parameters_from_json


def boltzmann_occupancy(redox_mc: RedoxMC) -> dict[int, float]:
    energies = {
        state: redox_mc.compute_state_energy(state, pe_nonbonded=0.0)
        for state in redox_mc.params.allowed_states
    }
    e_min = min(energies.values())
    weights = {
        state: np.exp(-(energy - e_min) / redox_mc.kT)
        for state, energy in energies.items()
    }
    total = sum(weights.values())
    return {state: weight / total for state, weight in weights.items()}


def crossing_voltage(x: np.ndarray, y: np.ndarray, target: float) -> float:
    for i in range(len(x) - 1):
        y1 = y[i] - target
        y2 = y[i + 1] - target
        if y1 == 0:
            return float(x[i])
        if y2 == 0:
            return float(x[i + 1])
        if y1 * y2 < 0:
            return float(x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return float("nan")


def main():
    parser = argparse.ArgumentParser(description="Generate experiment-friendly redox report")
    parser.add_argument("--redox-config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--v-start", type=float, default=-1.6, help="Start voltage vs vacuum")
    parser.add_argument("--v-end", type=float, default=-0.4, help="End voltage vs vacuum")
    parser.add_argument("--v-step", type=float, default=0.01, help="Voltage step vs vacuum")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args()

    params = load_redox_parameters_from_json(args.redox_config)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.redox_config) as f:
        cfg = json.load(f)
    offset = float(cfg.get("metadata", {}).get("vacuum_to_li_offset_ev", 1.4))

    voltages_vac = np.arange(args.v_start, args.v_end + args.v_step / 2, args.v_step)
    voltages_li = voltages_vac + offset

    curve_rows = []
    summary_rows = []

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    for mol_name, redox_params in params.items():
        fractions = []
        mean_charges = []

        for voltage_vac in voltages_vac:
            mc = RedoxMC(
                redox_params=redox_params,
                temperature_k=args.temperature_k,
                voltage_v=float(voltage_vac),
            )
            occ = boltzmann_occupancy(mc)
            reduced_fraction = sum(
                frac
                for state, frac in occ.items()
                if redox_params.charge_states[state] < mc.reference_charge
            )
            mean_charge = sum(
                occ[state] * redox_params.charge_states[state]
                for state in redox_params.allowed_states
            )
            fractions.append(reduced_fraction)
            mean_charges.append(mean_charge)
            curve_rows.append(
                {
                    "molecule": mol_name,
                    "voltage_vacuum_v": voltage_vac,
                    "voltage_li_v": voltage_vac + offset,
                    "fraction_reduced": reduced_fraction,
                    "mean_charge_e": mean_charge,
                }
            )

        fractions = np.array(fractions)
        mean_charges = np.array(mean_charges)
        reduction_response = -np.gradient(fractions, voltages_li)

        e10 = crossing_voltage(voltages_li[::-1], fractions[::-1], 0.10)
        e50 = crossing_voltage(voltages_li[::-1], fractions[::-1], 0.50)
        e90 = crossing_voltage(voltages_li[::-1], fractions[::-1], 0.90)
        peak_idx = int(np.argmax(reduction_response))
        peak_v = float(voltages_li[peak_idx])
        peak_height = float(reduction_response[peak_idx])

        summary_rows.append(
            {
                "molecule": mol_name,
                "onset_10pct_v_vs_li": e10,
                "half_wave_50pct_v_vs_li": e50,
                "onset_90pct_v_vs_li": e90,
                "lsv_proxy_peak_v_vs_li": peak_v,
                "lsv_proxy_peak_height": peak_height,
                "more_easily_reduced_rank_key": -e50 if not np.isnan(e50) else np.nan,
            }
        )

        axes[0].plot(voltages_li, reduction_response, label=f"{mol_name}")
        axes[1].plot(voltages_li, fractions, label=f"{mol_name}")

        if not np.isnan(e50):
            axes[0].axvline(e50, linestyle="--", alpha=0.25)
            axes[1].axvline(e50, linestyle="--", alpha=0.25)

    axes[0].set_title("LSV-like Reduction Response Proxy")
    axes[0].set_xlabel("Potential (V vs Li/Li+)")
    axes[0].set_ylabel(r"$-d(f_{reduced})/dV$ (arb. units)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].set_title("CV/Equilibrium Reduction Fraction")
    axes[1].set_xlabel("Potential (V vs Li/Li+)")
    axes[1].set_ylabel("Fraction Reduced")
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.suptitle("Electrolyte Redox Screening View")
    fig.tight_layout()
    png_path = output_dir / "electrolyte_lsv_cv_like.png"
    fig.savefig(png_path, dpi=220)

    curves_df = pd.DataFrame(curve_rows)
    curves_df.to_csv(output_dir / "electrolyte_redox_curves.csv", index=False)

    summary_df = pd.DataFrame(summary_rows).sort_values("half_wave_50pct_v_vs_li", ascending=False)
    summary_df.to_csv(output_dir / "electrolyte_redox_summary.csv", index=False)

    report_path = output_dir / "REPORT.md"
    with report_path.open("w") as f:
        f.write("# Electrolyte Redox Screening Report\n\n")
        f.write("## What This Figure Means\n\n")
        f.write("- Left panel: an LSV-like reduction-response proxy based on the voltage derivative of reduced-state occupancy.\n")
        f.write("- Right panel: the equilibrium reduced fraction, which is the clearest way to read reduction onset and ordering.\n")
        f.write("- Potentials are reported in the experimental convention `V vs Li/Li+`.\n\n")
        f.write("## Main Conclusions\n\n")
        for _, row in summary_df.iterrows():
            f.write(
                f"- `{row['molecule']}`: onset(10%) `{row['onset_10pct_v_vs_li']:.3f} V`, "
                f"E1/2 `{row['half_wave_50pct_v_vs_li']:.3f} V`, "
                f"LSV-like peak `{row['lsv_proxy_peak_v_vs_li']:.3f} V vs Li/Li+`.\n"
            )
        if len(summary_df) >= 2:
            ranking = " > ".join(summary_df["molecule"].tolist())
            f.write(f"- Easier to reduce ranking in this model: `{ranking}`.\n")
        f.write("\n## Method\n\n")
        f.write("1. Read the redox configuration for each molecule.\n")
        f.write("2. Convert the configured redox free-energy model into Boltzmann occupancies over a voltage grid.\n")
        f.write("3. Compute `fraction_reduced(V)` and use `-d(fraction_reduced)/dV` as an LSV-like proxy.\n")
        f.write("4. Report onset (10%), half-wave (50%), and near-complete reduction (90%) potentials.\n\n")
        f.write("## Principle\n\n")
        f.write("- The model is a constant-potential grand-canonical occupancy model.\n")
        f.write("- A more negative potential stabilizes reduced states.\n")
        f.write("- The crossing where oxidized and reduced populations are equal is the model half-wave potential.\n")
        f.write("- This is best interpreted as a thermodynamic reduction-onset screen, not a direct current prediction.\n\n")
        f.write("## Why Experimentalists Can Use It\n\n")
        f.write("- It outputs `V vs Li/Li+`, onset, and LSV-like peak positions directly.\n")
        f.write("- It makes molecular ranking obvious before expensive electrochemistry or full DFT workflows.\n")
        f.write("- It is well suited for formulation screening and sensitivity analysis.\n\n")
        f.write("## Cost and Throughput\n\n")
        f.write("- The equilibrium screening mode used here is cheap: seconds to minutes per formulation.\n")
        f.write("- The MD-coupled mode is much more expensive: minutes to hours per full sweep on one GPU.\n")
        f.write("- This supports rapid pre-screening followed by expensive refinement only for short-listed candidates.\n\n")
        f.write("## Accuracy and Limits\n\n")
        f.write("- Expected accuracy in the current effective model is qualitative to semi-quantitative.\n")
        f.write("- Ranking and onset order can be informative, but peak currents and exact LSV shapes are not predicted.\n")
        f.write("- Precision depends on how well the redox free-energy offsets reflect real solvation, ion pairing, and decomposition chemistry.\n")
        f.write("- EC/DMC reduction is often irreversible and tied to SEI-forming chemistry, so these outputs should be read as onset tendencies.\n\n")
        f.write("## Can It Be Used for Screening?\n\n")
        f.write("- Yes, for molecular ordering, identifying likely first-to-reduce species, and narrowing the voltage window of interest.\n")
        f.write("- Not yet as a direct replacement for experimental LSV/CV current curves.\n\n")
        f.write("## Next Optimizations\n\n")
        f.write("1. Couple RedoxMC to frame-resolved solvation energies from MD instead of the current effective offsets.\n")
        f.write("2. Add oxidation branches so high-voltage stability can be screened symmetrically.\n")
        f.write("3. Fit offsets against a training set of known experimental onsets.\n")
        f.write("4. Include reaction channels for irreversible solvent decomposition and SEI formation.\n")
        f.write("5. Build a formulation parser that maps composition -> molecules -> ranked redox windows.\n\n")
        f.write("## Future Product Direction\n\n")
        f.write("- Input: an electrolyte recipe (solvents, salt, additives, concentrations).\n")
        f.write("- Output: ranked molecular reduction/oxidation potentials, likely first-reduced species, and a high-voltage stability flag.\n")
        f.write("- For a target such as `5 V`, the tool would compare the oxidation side of each component against the applied window and report likely unstable species.\n")
        f.write("- The long-term value is fast formulation triage before synthesis and electrochemical testing.\n")

    print(png_path)
    print(report_path)


if __name__ == "__main__":
    main()
