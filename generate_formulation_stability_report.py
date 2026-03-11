#!/usr/bin/env python
"""
Generate a formulation-level reduction/oxidation ranking and high-voltage
stability report.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def main():
    parser = argparse.ArgumentParser(description="Generate formulation stability report")
    parser.add_argument("--reduction-config", required=True)
    parser.add_argument("--oxidation-config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-high-voltage-v", type=float, default=5.0)
    args = parser.parse_args()

    red_params = load_redox_parameters_from_json(args.reduction_config)
    ox_params = load_redox_parameters_from_json(args.oxidation_config)
    with open(args.reduction_config) as f:
        red_cfg = json.load(f)
    offset = float(red_cfg.get("metadata", {}).get("vacuum_to_li_offset_ev", 1.4))

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    voltages_li_red = np.arange(-0.2, 1.6 + 0.005, 0.005)
    voltages_li_ox = np.arange(3.5, 5.6 + 0.005, 0.005)

    red_df = build_branch_summary(red_params, voltages_li_red, offset, branch="reduction")
    ox_df = build_branch_summary(ox_params, voltages_li_ox, offset, branch="oxidation")
    summary = red_df.merge(ox_df, on="molecule", how="outer")

    formulation_low_limit = summary["reduction_onset_10pct_v_vs_li"].max()
    formulation_high_limit = summary["oxidation_onset_10pct_v_vs_li"].min()
    window_width = formulation_high_limit - formulation_low_limit
    stable_at_target = bool(args.target_high_voltage_v <= formulation_high_limit)

    summary["reduction_rank"] = summary["reduction_half_wave_50pct_v_vs_li"].rank(ascending=False, method="dense")
    summary["oxidation_rank"] = summary["oxidation_half_wave_50pct_v_vs_li"].rank(ascending=False, method="dense")
    summary["likely_reduce_on_anode"] = summary["reduction_onset_10pct_v_vs_li"] > 0.2
    summary["likely_oxidize_by_target"] = summary["oxidation_onset_10pct_v_vs_li"] <= args.target_high_voltage_v
    summary.to_csv(output_dir / "formulation_stability_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    y = np.arange(len(summary))

    axes[0].barh(y, summary["reduction_half_wave_50pct_v_vs_li"], color="#2c7fb8")
    axes[0].set_yticks(y, labels=summary["molecule"])
    axes[0].axvline(0.0, color="black", linestyle="--", alpha=0.5, label="Li/graphite region")
    axes[0].set_xlabel("Reduction E1/2 (V vs Li/Li+)")
    axes[0].set_title("Reduction Liability")
    axes[0].grid(True, axis="x", alpha=0.3)

    axes[1].barh(y, summary["oxidation_half_wave_50pct_v_vs_li"], color="#d95f0e")
    axes[1].set_yticks(y, labels=summary["molecule"])
    axes[1].axvline(args.target_high_voltage_v, color="red", linestyle="--", alpha=0.7, label=f"Target {args.target_high_voltage_v:.1f} V")
    axes[1].set_xlabel("Oxidation E1/2 (V vs Li/Li+)")
    axes[1].set_title("High-Voltage Risk")
    axes[1].grid(True, axis="x", alpha=0.3)
    axes[1].legend()

    fig.suptitle("Formulation Stability Screening")
    fig.tight_layout()
    fig.savefig(output_dir / "formulation_stability_window.png", dpi=220)

    report = output_dir / "FORMULATION_STABILITY_REPORT_CN.md"
    with report.open("w") as f:
        f.write("# 电解液配方稳定性筛选报告\n\n")
        f.write("## 一句话结论\n\n")
        if stable_at_target:
            f.write(f"- 按当前有效模型，配方在 `{args.target_high_voltage_v:.1f} V vs Li/Li+` 下仍位于预测氧化起始电位之下，理论上可进入进一步验证。\n")
        else:
            f.write(f"- 按当前有效模型，配方在 `{args.target_high_voltage_v:.1f} V vs Li/Li+` 下已经超过预测氧化起始电位，不能视为高压稳定。\n")
        f.write(f"- 预测还原侧稳定下限约为 `{formulation_low_limit:.3f} V vs Li/Li+`。\n")
        f.write(f"- 预测氧化侧稳定上限约为 `{formulation_high_limit:.3f} V vs Li/Li+`。\n")
        f.write(f"- 预测电化学稳定窗口宽度约为 `{window_width:.3f} V`。\n\n")

        f.write("## 分子排序\n\n")
        for _, row in summary.sort_values("reduction_half_wave_50pct_v_vs_li", ascending=False).iterrows():
            f.write(
                f"- `{row['molecule']}`: 还原 E1/2 `{row['reduction_half_wave_50pct_v_vs_li']:.3f} V`, "
                f"氧化 E1/2 `{row['oxidation_half_wave_50pct_v_vs_li']:.3f} V vs Li/Li+`。\n"
            )
        f.write("\n## 对实验人员怎么读\n\n")
        f.write("- 左图越靠右，说明越容易在负极侧先被还原。\n")
        f.write("- 右图若落在目标高电压左边，说明在高压正极侧更容易先被氧化。\n")
        f.write("- 一种分子若既容易在低电位被还原，又容易在高电位被氧化，则它对宽电压窗口是不利的。\n\n")

        f.write("## 方法步骤\n\n")
        f.write("1. 为每个分子建立还原侧和氧化侧的有效自由能模型。\n")
        f.write("2. 用 Boltzmann 分布计算不同电位下的占有率。\n")
        f.write("3. 提取 onset(10%)、E1/2(50%) 和近完全转化点(90%)。\n")
        f.write("4. 用还原侧最高 onset 和氧化侧最低 onset 近似配方稳定窗口。\n")
        f.write("5. 将目标高压（如 5.0 V）与氧化侧窗口上限对比，给出风险判断。\n\n")

        f.write("## 原理\n\n")
        f.write("- 这是一个分子层面的热力学筛选模型。\n")
        f.write("- 它回答的是“哪一种分子更先进入电子转移不稳定区”。\n")
        f.write("- 它不是直接求实验电流，而是求 redox 起始与排序。\n\n")

        f.write("## 这样做的好处\n\n")
        f.write("- 能在实验之前快速筛掉明显不合适的配方。\n")
        f.write("- 能明确指出是哪一个分子限制了高压或低压稳定性。\n")
        f.write("- 能把实验 LSV/CV 的观察转成分子层面的解释和排序。\n\n")

        f.write("## 成本\n\n")
        f.write("- 有效模型筛选：秒级到分钟级。\n")
        f.write("- 加入常电位 MD 精修：单个 sweep 通常是分钟到小时级 GPU 计算。\n")
        f.write("- 加入反应路径和分解势垒：成本最高，但精度也最高。\n\n")

        f.write("## 能否用于筛选\n\n")
        f.write("- 可以，尤其适合溶剂、盐、添加剂的初筛和排序。\n")
        f.write("- 现阶段最适合用于“排队”和“缩小实验窗口”，不宜直接替代最终电化学测试。\n\n")

        f.write("## 精度判断\n\n")
        f.write("- 当前有效模型：定性到半定量。\n")
        f.write("- 如果用实验 onset 或 DFT/SCRF 数据校准，自由能参数可以更接近半定量预测。\n")
        f.write("- 如果要判断 5 V 真稳定，还必须补上氧化后的分解路径和界面动力学。\n\n")

        f.write("## 下一步优化\n\n")
        f.write("1. 用真实 DFT/溶剂化自由能替换有效 offset。\n")
        f.write("2. 加入氧化后的不可逆分解通道。\n")
        f.write("3. 用常电位 MD 提取界面溶剂化修正，降低纯参数模型误差。\n")
        f.write("4. 用实验 LSV/CV onset 数据回标定各类溶剂与添加剂。\n")
        f.write("5. 自动化输入配方 -> 输出稳定窗口、排序和风险标签。\n\n")

        f.write("## 面向未来的应用\n\n")
        f.write("- 目标形态是：给一个电解液配方，自动输出其中所有组分的还原/氧化电位、排序，以及在目标电压窗口下的稳定性判断。\n")
        f.write(f"- 对 `5.0 V` 这样的目标高压，系统将自动识别哪一个组分最先氧化，并给出是否需要替换溶剂/盐/添加剂的建议。\n")
        f.write("- 结合数据库和机器学习后，可以变成一个面向高通量电解液设计的前端筛选器。\n")

    print(output_dir / "formulation_stability_window.png")
    print(report)


if __name__ == "__main__":
    main()
