#!/usr/bin/env python
"""
Batch-screen multiple electrolyte formulations and summarize reduction,
oxidation, and reactive high-voltage limits in one table/figure.
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


KBT_OVER_H_S = 6.251e12
KB_KJMOL_PER_K = 8.314462e-3


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


def threshold_voltage_from_rate(x: np.ndarray, y: np.ndarray, target: float) -> float:
    if len(x) == 0:
        return float("nan")
    if y[0] >= target:
        return float(x[0])
    if y[-1] < target:
        return float("nan")
    return crossing_voltage(x, y, target)


def build_branch_summary(params: dict, voltages_li: np.ndarray, offset: float, branch: str) -> pd.DataFrame:
    rows = []
    for mol_name, redox_params in params.items():
        metric = []
        for voltage_li in voltages_li:
            voltage_vac = voltage_li - offset
            mc = RedoxMC(redox_params=redox_params, temperature_k=300.0, voltage_v=float(voltage_vac))
            occ = boltzmann_occupancy(mc)
            if branch == "reduction":
                value = sum(
                    frac
                    for state, frac in occ.items()
                    if redox_params.charge_states[state] < mc.reference_charge
                )
            else:
                value = sum(
                    frac
                    for state, frac in occ.items()
                    if redox_params.charge_states[state] > mc.reference_charge
                )
            metric.append(value)

        metric = np.array(metric)
        rows.append(
            {
                "molecule": mol_name,
                f"{branch}_onset_10pct_v_vs_li": crossing_voltage(voltages_li, metric, 0.10),
                f"{branch}_half_wave_50pct_v_vs_li": crossing_voltage(voltages_li, metric, 0.50),
            }
        )
    return pd.DataFrame(rows)


def channel_barrier(channel: dict, voltage_li: float) -> float:
    onset = float(channel["onset_v_vs_li"])
    barrier = float(channel["barrier_at_onset_kjmol"])
    slope = float(channel.get("barrier_slope_kjmol_per_v", 0.0))
    floor = float(channel.get("floor_barrier_kjmol", 0.0))
    overdrive = max(0.0, voltage_li - onset)
    return max(floor, barrier - slope * overdrive)


def channel_rate(channel: dict, voltage_li: float, temperature_k: float) -> float:
    kT = KB_KJMOL_PER_K * temperature_k
    barrier = channel_barrier(channel, voltage_li)
    return KBT_OVER_H_S * np.exp(-barrier / kT)


def evaluate_formulation(entry: dict, base_dir: Path) -> tuple[dict, pd.DataFrame]:
    reduction_config = (base_dir / entry["reduction_config"]).resolve()
    oxidation_config = (base_dir / entry["oxidation_config"]).resolve()
    pathway_config = (base_dir / entry["pathway_config"]).resolve()
    target_v = float(entry.get("target_high_voltage_v", 5.0))
    rate_threshold = float(entry.get("rate_threshold_s^-1", 1.0e-4))

    red_params = load_redox_parameters_from_json(str(reduction_config))
    ox_params = load_redox_parameters_from_json(str(oxidation_config))
    with reduction_config.open() as f:
        red_cfg = json.load(f)
    with pathway_config.open() as f:
        pathway_cfg = json.load(f)

    offset = float(red_cfg.get("metadata", {}).get("vacuum_to_li_offset_ev", 1.4))
    temperature_k = float(pathway_cfg.get("metadata", {}).get("temperature_k", 300.0))
    voltages_li_red = np.arange(-0.2, 1.6 + 0.005, 0.005)
    voltages_li_ox = np.arange(3.5, 5.6 + 0.005, 0.005)

    red_df = build_branch_summary(red_params, voltages_li_red, offset, branch="reduction")
    ox_df = build_branch_summary(ox_params, voltages_li_ox, offset, branch="oxidation")
    mol_summary = red_df.merge(ox_df, on="molecule", how="outer")

    pathway_rows = []
    for mol_name, mol_cfg in pathway_cfg.get("molecules", {}).items():
        channels = mol_cfg.get("channels", [])
        if not channels:
            continue
        rate_matrix = np.array(
            [[channel_rate(channel, v, temperature_k) for v in voltages_li_ox] for channel in channels]
        )
        total_rate = rate_matrix.sum(axis=0)
        threshold_v = threshold_voltage_from_rate(
            voltages_li_ox, np.log10(total_rate + 1e-300), np.log10(rate_threshold)
        )
        target_idx = int(np.argmin(np.abs(voltages_li_ox - target_v)))
        fastest_idx = int(np.argmax(rate_matrix[:, target_idx]))
        pathway_rows.append(
            {
                "molecule": mol_name,
                "pathway_failure_threshold_v_vs_li": threshold_v,
                "pathway_rate_at_target_s^-1": float(total_rate[target_idx]),
                "fastest_channel_at_target": channels[fastest_idx]["name"],
            }
        )
    pathway_df = pd.DataFrame(pathway_rows)
    mol_summary = mol_summary.merge(pathway_df, on="molecule", how="left")

    reduction_limit = float(mol_summary["reduction_onset_10pct_v_vs_li"].max())
    oxidation_series = mol_summary["oxidation_onset_10pct_v_vs_li"].dropna()
    pathway_series = mol_summary["pathway_failure_threshold_v_vs_li"].dropna()
    oxidation_limit = float(oxidation_series.min()) if not oxidation_series.empty else float("nan")
    pathway_limit = float(pathway_series.min()) if not pathway_series.empty else float("nan")
    effective_high_limit = float(np.nanmin([oxidation_limit, pathway_limit]))

    limiting_reduction = mol_summary.loc[mol_summary["reduction_onset_10pct_v_vs_li"].idxmax(), "molecule"]
    limiting_oxidation = mol_summary.loc[oxidation_series.idxmin(), "molecule"] if not oxidation_series.empty else ""
    limiting_pathway = mol_summary.loc[pathway_series.idxmin(), "molecule"] if not pathway_series.empty else ""
    if np.isnan(pathway_limit) or oxidation_limit <= pathway_limit:
        limiting_factor = "reversible_oxidation"
        limiting_molecule = limiting_oxidation
    else:
        limiting_factor = "reactive_pathway"
        limiting_molecule = limiting_pathway

    formulation_row = {
        "formulation": entry["name"],
        "target_high_voltage_v": target_v,
        "reduction_low_limit_v_vs_li": reduction_limit,
        "oxidation_high_limit_v_vs_li": oxidation_limit,
        "pathway_high_limit_v_vs_li": pathway_limit,
        "effective_high_limit_v_vs_li": effective_high_limit,
        "stable_at_target": bool(target_v <= effective_high_limit),
        "limiting_molecule": limiting_molecule,
        "limiting_factor": limiting_factor,
        "limiting_reduction_molecule": limiting_reduction,
        "limiting_pathway_molecule": limiting_pathway,
    }
    return formulation_row, mol_summary


def main():
    parser = argparse.ArgumentParser(description="Batch screen multiple electrolyte formulations")
    parser.add_argument("--input", required=True, help="JSON file listing formulation entries")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    base_dir = input_path.parent
    with input_path.open() as f:
        formulations = json.load(f)

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    formulation_rows = []
    molecule_rows = []
    for entry in formulations:
        formulation_row, mol_summary = evaluate_formulation(entry, base_dir)
        formulation_rows.append(formulation_row)
        mol_summary = mol_summary.copy()
        mol_summary.insert(0, "formulation", entry["name"])
        molecule_rows.append(mol_summary)

    formulation_df = pd.DataFrame(formulation_rows).sort_values(
        ["stable_at_target", "effective_high_limit_v_vs_li"], ascending=[False, False]
    )
    molecule_df = pd.concat(molecule_rows, ignore_index=True) if molecule_rows else pd.DataFrame()

    formulation_df.to_csv(output_dir / "formulation_batch_summary.csv", index=False)
    molecule_df.to_csv(output_dir / "formulation_batch_molecule_breakdown.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    y = np.arange(len(formulation_df))

    axes[0].barh(y, formulation_df["effective_high_limit_v_vs_li"], color="#d95f0e")
    axes[0].set_yticks(y, labels=formulation_df["formulation"])
    for i, target_v in enumerate(formulation_df["target_high_voltage_v"]):
        axes[0].axvline(target_v, color="red", linestyle="--", alpha=0.15)
    axes[0].set_xlabel("Effective High-Voltage Limit (V vs Li/Li+)")
    axes[0].set_title("Formulation Ranking")
    axes[0].grid(True, axis="x", alpha=0.3)

    window_width = formulation_df["effective_high_limit_v_vs_li"] - formulation_df["reduction_low_limit_v_vs_li"]
    colors = ["#2ca25f" if stable else "#de2d26" for stable in formulation_df["stable_at_target"]]
    axes[1].barh(y, window_width, color=colors)
    axes[1].set_yticks(y, labels=formulation_df["formulation"])
    axes[1].set_xlabel("Predicted Stability Window Width (V)")
    axes[1].set_title("Target-Voltage Pass/Fail")
    axes[1].grid(True, axis="x", alpha=0.3)

    fig.suptitle("Batch Electrolyte Formulation Screening")
    fig.tight_layout()
    fig.savefig(output_dir / "formulation_batch_overview.png", dpi=220)

    with (output_dir / "FORMULATION_BATCH_REPORT_CN.md").open("w") as f:
        f.write("# 配方批处理筛选报告\n\n")
        f.write("## 结果摘要\n\n")
        for _, row in formulation_df.iterrows():
            status = "通过" if row["stable_at_target"] else "不通过"
            f.write(
                f"- `{row['formulation']}`: 目标 `{row['target_high_voltage_v']:.2f} V`, "
                f"有效高压上限 `{row['effective_high_limit_v_vs_li']:.3f} V`, "
                f"低压下限 `{row['reduction_low_limit_v_vs_li']:.3f} V`, "
                f"限制分子 `{row['limiting_molecule']}` ({row['limiting_factor']}), "
                f"结论 `{status}`。\n"
            )

        f.write("\n## 这套批处理在做什么\n\n")
        f.write("- 对每个配方自动读取还原、氧化、反应通道配置。\n")
        f.write("- 先给出分子级的还原/氧化 onset 与 E1/2。\n")
        f.write("- 再把反应通道势垒转成速率，给出更严格的高压失效上限。\n")
        f.write("- 最后把每个配方的稳定窗口和目标电压 pass/fail 汇总成一个总表。\n\n")

        f.write("## 为什么比单分子 redox potential 更适合配方初筛\n\n")
        f.write("- 单分子计算擅长回答“这个分子本身在理想化环境下大概何时得失电子”。\n")
        f.write("- 这套批处理更擅长回答“放进具体配方后，谁是短板、5 V 能不能过、先坏的是哪条通道”。\n")
        f.write("- 因此它适合做配方排序；单分子计算更适合做关键分子的精修和参数校准。\n")


if __name__ == "__main__":
    main()
