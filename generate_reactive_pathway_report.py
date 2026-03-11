#!/usr/bin/env python
"""
Generate a formulation stability report with explicit reaction channels and
voltage-dependent activation barriers.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def build_branch_summary(params: dict, voltages_li: np.ndarray, offset: float, branch: str):
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
                f"{branch}_onset_90pct_v_vs_li": crossing_voltage(voltages_li, metric, 0.90),
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


def build_pathway_summary(pathway_cfg: dict, voltages_li: np.ndarray, target_voltage_v: float, rate_threshold: float):
    temperature_k = float(pathway_cfg.get("metadata", {}).get("temperature_k", 300.0))
    rows = []
    pathway_curves = {}
    for mol_name, mol_cfg in pathway_cfg.get("molecules", {}).items():
        channels = mol_cfg.get("channels", [])
        if not channels:
            continue

        rate_matrix = []
        for channel in channels:
            rate_matrix.append([channel_rate(channel, v, temperature_k) for v in voltages_li])
        rate_matrix = np.array(rate_matrix)
        total_rate = rate_matrix.sum(axis=0)
        fastest_idx = int(np.argmax([channel_rate(ch, target_voltage_v, temperature_k) for ch in channels]))
        fastest_channel = channels[fastest_idx]
        threshold_voltage = threshold_voltage_from_rate(
            voltages_li, np.log10(total_rate + 1e-300), np.log10(rate_threshold)
        )

        pathway_curves[mol_name] = {
            "total_rate": total_rate,
            "channel_rates": rate_matrix,
            "channels": channels,
        }

        rows.append(
            {
                "molecule": mol_name,
                "pathway_failure_threshold_v_vs_li": threshold_voltage,
                "pathway_rate_at_target_s^-1": float(total_rate[np.argmin(np.abs(voltages_li - target_voltage_v))]),
                "fastest_channel_at_target": fastest_channel["name"],
                "fastest_channel_products": fastest_channel.get("dominant_products", ""),
            }
        )

    return pd.DataFrame(rows), pathway_curves


def main():
    parser = argparse.ArgumentParser(description="Generate reactive pathway stability report")
    parser.add_argument("--reduction-config", required=True)
    parser.add_argument("--oxidation-config", required=True)
    parser.add_argument("--pathway-config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-high-voltage-v", type=float, default=5.0)
    parser.add_argument("--rate-threshold", type=float, default=1.0e-4)
    args = parser.parse_args()

    red_params = load_redox_parameters_from_json(args.reduction_config)
    ox_params = load_redox_parameters_from_json(args.oxidation_config)
    with open(args.reduction_config) as f:
        red_cfg = json.load(f)
    with open(args.pathway_config) as f:
        pathway_cfg = json.load(f)
    offset = float(red_cfg.get("metadata", {}).get("vacuum_to_li_offset_ev", 1.4))

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    voltages_li_red = np.arange(-0.2, 1.6 + 0.005, 0.005)
    voltages_li_ox = np.arange(3.5, 5.6 + 0.005, 0.005)

    red_df = build_branch_summary(red_params, voltages_li_red, offset, branch="reduction")
    ox_df = build_branch_summary(ox_params, voltages_li_ox, offset, branch="oxidation")
    path_df, pathway_curves = build_pathway_summary(
        pathway_cfg, voltages_li_ox, args.target_high_voltage_v, args.rate_threshold
    )

    summary = red_df.merge(ox_df, on="molecule", how="outer").merge(path_df, on="molecule", how="outer")
    reduction_limit = summary["reduction_onset_10pct_v_vs_li"].max()
    oxidation_limit = summary["oxidation_onset_10pct_v_vs_li"].min()
    pathway_limit = summary["pathway_failure_threshold_v_vs_li"].min()
    effective_high_limit = np.nanmin([oxidation_limit, pathway_limit])
    stable_at_target = bool(args.target_high_voltage_v <= effective_high_limit)

    summary["reduction_rank"] = summary["reduction_half_wave_50pct_v_vs_li"].rank(ascending=False, method="dense")
    summary["oxidation_rank"] = summary["oxidation_half_wave_50pct_v_vs_li"].rank(ascending=False, method="dense")
    summary["pathway_rank"] = summary["pathway_failure_threshold_v_vs_li"].rank(ascending=False, method="dense")
    summary["likely_unstable_at_target"] = summary["pathway_rate_at_target_s^-1"] >= args.rate_threshold
    summary.to_csv(output_dir / "reactive_pathway_summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    y = np.arange(len(summary))

    axes[0].barh(y, summary["reduction_half_wave_50pct_v_vs_li"], color="#2c7fb8")
    axes[0].set_yticks(y, labels=summary["molecule"])
    axes[0].set_xlabel("Reduction E1/2 (V vs Li/Li+)")
    axes[0].set_title("Reduction Ranking")
    axes[0].grid(True, axis="x", alpha=0.3)

    axes[1].barh(y, summary["oxidation_half_wave_50pct_v_vs_li"], color="#d95f0e")
    axes[1].set_yticks(y, labels=summary["molecule"])
    axes[1].axvline(args.target_high_voltage_v, color="red", linestyle="--", alpha=0.7, label=f"Target {args.target_high_voltage_v:.1f} V")
    axes[1].set_xlabel("Oxidation E1/2 (V vs Li/Li+)")
    axes[1].set_title("Oxidation Ranking")
    axes[1].grid(True, axis="x", alpha=0.3)
    axes[1].legend()

    for mol_name, curves in pathway_curves.items():
        axes[2].plot(voltages_li_ox, np.log10(curves["total_rate"] + 1e-300), label=mol_name)
    axes[2].axvline(args.target_high_voltage_v, color="red", linestyle="--", alpha=0.7)
    axes[2].axhline(np.log10(args.rate_threshold), color="black", linestyle=":", alpha=0.7, label="Failure threshold")
    axes[2].set_xlabel("Potential (V vs Li/Li+)")
    axes[2].set_ylabel(r"log$_{10}$(estimated decomposition rate / s$^{-1}$)")
    axes[2].set_title("Explicit Pathway Kinetic Risk")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()

    fig.suptitle("Reactive High-Voltage Stability Screening")
    fig.tight_layout()
    fig.savefig(output_dir / "reactive_pathway_window.png", dpi=220)

    report = output_dir / "REACTIVE_PATHWAY_REPORT_CN.md"
    with report.open("w") as f:
        f.write("# 高压失效反应通道筛选报告\n\n")
        f.write("## 一句话结论\n\n")
        if stable_at_target:
            f.write(f"- 按当前反应通道+势垒模型，配方在 `{args.target_high_voltage_v:.1f} V vs Li/Li+` 下尚未超过最快失效通道的速率阈值，可继续验证。\n")
        else:
            f.write(f"- 按当前反应通道+势垒模型，配方在 `{args.target_high_voltage_v:.1f} V vs Li/Li+` 下已经跨过最快失效通道的速率阈值，不应视为高压稳定。\n")
        f.write(f"- 还原侧下限约 `{reduction_limit:.3f} V vs Li/Li+`。\n")
        f.write(f"- 可逆氧化上限约 `{oxidation_limit:.3f} V vs Li/Li+`。\n")
        f.write(f"- 反应通道失效上限约 `{pathway_limit:.3f} V vs Li/Li+`。\n")
        f.write(f"- 取更严格标准后的有效高压上限约 `{effective_high_limit:.3f} V vs Li/Li+`。\n\n")

        f.write("## 分子与通道判断\n\n")
        for _, row in summary.sort_values("pathway_failure_threshold_v_vs_li").iterrows():
            f.write(
                f"- `{row['molecule']}`: 还原 E1/2 `{row['reduction_half_wave_50pct_v_vs_li']:.3f} V`, "
                f"氧化 E1/2 `{row['oxidation_half_wave_50pct_v_vs_li']:.3f} V`, "
                f"反应失效阈值 `{row['pathway_failure_threshold_v_vs_li']:.3f} V vs Li/Li+`。\n"
            )
            f.write(
                f"  目标电压下最快通道: `{row['fastest_channel_at_target']}`, "
                f"估算总分解速率 `{row['pathway_rate_at_target_s^-1']:.3e} s^-1`。\n"
            )
            if isinstance(row["fastest_channel_products"], str) and row["fastest_channel_products"]:
                f.write(f"  主要产物/后果: {row['fastest_channel_products']}。\n")

        f.write("\n## 这一步比前一版多了什么\n\n")
        f.write("- 前一版只说“到这个电位附近可能不稳定”。\n")
        f.write("- 这一版进一步说“是哪一条反应通道先跑起来，以及在目标电压下它跑得有多快”。\n")
        f.write("- 因此它更接近实验上的失效，而不只是可逆氧化 onset。\n\n")

        f.write("## 原理\n\n")
        f.write("1. 还原/氧化部分仍然用分子态占有率模型给出 onset 和 E1/2。\n")
        f.write("2. 对每条失效通道，定义一个参考势垒 `DeltaG_dagger(onset)`。\n")
        f.write("3. 当电位继续升高时，势垒按设定斜率降低。\n")
        f.write("4. 用 `k = (kBT/h) * exp(-DeltaG_dagger/kT)` 把势垒转成分解速率。\n")
        f.write("5. 当总分解速率超过设定阈值时，把该电位视为可观测失效边界。\n\n")

        f.write("## 它现在算到什么程度\n\n")
        f.write("- 已经从“一个分解 onset”升级为“多条反应通道 + 电位相关势垒 + 速率阈值”。\n")
        f.write("- 已经能区分 EC 和 DMC 哪个更先坏，以及为什么坏。\n")
        f.write("- 仍然不是显式 DFT 反应路径，也不是显式产物网络模拟。\n\n")

        f.write("## 好处\n\n")
        f.write("- 比单个 oxidation onset 更接近真实失效机理。\n")
        f.write("- 能直接回答 5 V 下最先触发的是哪条化学通道。\n")
        f.write("- 对配方筛选比只看 HOMO/LUMO 或单个 redox potential 更实用。\n\n")

        f.write("## 仍然存在的边界\n\n")
        f.write("- 通道势垒目前仍是有效参数，不是第一性原理反应路径。\n")
        f.write("- 未显式包含表面催化、盐分解耦合、产物二次反应和 SEI/CEI 生长网络。\n")
        f.write("- 因此当前适合作为“机理感更强的筛选原型”，还不是终局预测器。\n\n")

        f.write("## 下一步如果要再往前走\n\n")
        f.write("1. 用 DFT/NEB 给每条关键通道提供真实势垒。\n")
        f.write("2. 把电极表面和 Li+ 配位态显式放进通道定义中。\n")
        f.write("3. 用常电位 MD 采样局域构型，再对代表构型做通道势垒统计。\n")
        f.write("4. 把单分子通道扩展成多步反应网络。\n")

    print(output_dir / "reactive_pathway_window.png")
    print(report)


if __name__ == "__main__":
    main()
