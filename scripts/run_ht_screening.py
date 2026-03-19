#!/usr/bin/env python3
"""
run_ht_screening.py
-------------------
高通量氧化还原电位筛选主程序。

给定分子数据库和配方生成脚本的结果，
批量跑 RedoxMC 筛选，输出结构化结果表。

Usage:
    python scripts/run_ht_screening.py \\
        --formulations db/formulations/ \\
        --output results/ht_screening/ \\
        --voltage-start -2.0 \\
        --voltage-end 0.0 \\
        --voltage-step 0.05

Output:
    results/ht_screening/
        ├── screening_summary.csv       # 所有配方的汇总表
        ├── screening_detailed.csv      # 详细数据
        ├── ranking_reduction.csv       # 还原稳定性排名
        ├── ranking_oxidation.csv       # 氧化稳定性排名
        └── screening_report.md         # Markdown 报告
"""

import argparse
import json
import sys
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = REPO_ROOT / "core"
if str(CORE_DIR) not in sys.path:
    sys.path.insert(0, str(CORE_DIR))

from redox_mc import RedoxMC

# 物理常数
FARADAY_KJMOL_PER_V = 96.48533212331002
KB_KJMOL_PER_K = 8.314462e-3
EV_TO_KJMOL = 96.4853
VACUUM_TO_LI_OFFSET_EV = 1.4  # 有机电解液的典型偏移


def boltzmann_occupancy(
    energies: Dict[int, float],
    temperature_k: float,
    charge_states: Dict[int, float],
) -> Dict[int, float]:
    """给定各态能量，计算 Boltzmann 占有率"""
    kT = KB_KJMOL_PER_K * temperature_k
    e_min = min(energies.values())
    weights = {state: np.exp(-(energy - e_min) / kT) for state, energy in energies.items()}
    total = sum(weights.values())
    return {state: weight / total for state, weight in weights.items()}


def crossing_voltage(
    voltages: np.ndarray,
    fractions: np.ndarray,
    target: float,
) -> float:
    """线性插值找 fraction = target 对应的电压 (V vs Li/Li+)"""
    for i in range(len(voltages) - 1):
        y1, y2 = fractions[i] - target, fractions[i + 1] - target
        if y1 == 0:
            return float(voltages[i])
        if y2 == 0:
            return float(voltages[i + 1])
        if y1 * y2 < 0:
            v1, v2 = voltages[i], voltages[i + 1]
            return float(v1 + (target - fractions[i]) * (v2 - v1) / (fractions[i + 1] - fractions[i]))
    return float("nan")


def screen_single_formulation(
    entry: Dict,
    voltage_start: float = -2.0,
    voltage_end: float = 0.0,
    voltage_step: float = 0.05,
    temperature_k: float = 300.0,
) -> Dict:
    """
    对单个配方跑 RedoxMC 筛选。

    电压扫描范围: voltage_start ~ voltage_end V vs vacuum
    转换为 Li/Li+ 需 + VACUUM_TO_LI_OFFSET_EV
    """
    redox_config = entry["redox_config"]
    molecules_cfg = redox_config.get("molecules", {})

    if not molecules_cfg:
        return {}

    # 构建 RedoxParameters 对象
    from redox_mc import RedoxParameters

    params_map = {}
    for mol_name, mol_cfg in molecules_cfg.items():
        charge_states = {int(k): float(v) for k, v in mol_cfg.get("charge_states", {}).items()}
        ea_dict = {}
        for key_str, ea_val in mol_cfg.get("electron_affinities", {}).items():
            s1, s2 = map(int, key_str.split(","))
            ea_dict[(s1, s2)] = float(ea_val)

        state_offsets_raw = mol_cfg.get("state_free_energy_offsets_kjmol", {})
        state_offsets = {int(k): float(v) for k, v in state_offsets_raw.items()} if state_offsets_raw else None

        params_map[mol_name] = RedoxParameters(
            name=mol_name,
            residue_names=mol_cfg.get("residue_names", [mol_name]),
            charge_states=charge_states,
            electron_affinities=ea_dict,
            allowed_states=list(charge_states.keys()),
            partition_coefficient=float(mol_cfg.get("partition_coefficient", 0.0)),
            state_free_energy_offsets_kjmol=state_offsets,
        )

    # 电压扫描
    voltages_vac = np.arange(voltage_start, voltage_end + voltage_step / 2, voltage_step)
    voltages_li = voltages_vac + VACUUM_TO_LI_OFFSET_EV

    results = {
        "formulation_name": entry.get("formulation_name", ""),
        "composition": entry.get("composition", ""),
        "salt": entry.get("salt_name", ""),
        "salt_conc_M": entry.get("salt_conc_M", 0.0),
        "additive": entry.get("additive_name", ""),
        "additive_wt_pct": entry.get("additive_wt_pct", 0.0),
    }

    for mol_name, params in params_map.items():
        # Reference charge
        ref_state = 0 if 0 in params.allowed_states else min(params.allowed_states)
        ref_charge = params.charge_states[ref_state]

        # 计算每个电压点下的还原占有率
        reduced_fractions = []
        oxidized_fractions = []

        for v_vac in voltages_vac:
            mc = RedoxMC(params, temperature_k=temperature_k, voltage_v=float(v_vac))
            occ = boltzmann_occupancy(
                {state: mc.compute_state_energy(state, pe_nonbonded=0.0) for state in mc.params.allowed_states},
                temperature_k,
                params.charge_states,
            )

            # 还原占有率（电荷比 ref 更负）
            reduced_frac = sum(
                frac for state, frac in occ.items() if params.charge_states[state] < ref_charge
            )
            oxidized_frac = sum(
                frac for state, frac in occ.items() if params.charge_states[state] > ref_charge
            )
            reduced_fractions.append(reduced_frac)
            oxidized_fractions.append(oxidized_frac)

        reduced_fractions = np.array(reduced_fractions)
        oxidized_fractions = np.array(oxidized_fractions)

        # 提取特征电位
        onset_10 = crossing_voltage(voltages_li, reduced_fractions, 0.10)
        onset_50 = crossing_voltage(voltages_li, reduced_fractions, 0.50)
        onset_90 = crossing_voltage(voltages_li, reduced_fractions, 0.90)

        results[f"{mol_name}_reduction_onset_10pct_V_vs_Li"] = onset_10
        results[f"{mol_name}_reduction_halfwave_V_vs_Li"] = onset_50
        results[f"{mol_name}_reduction_onset_90pct_V_vs_Li"] = onset_90

        # 氧化 onset
        ox_onset_10 = crossing_voltage(voltages_li, oxidized_fractions, 0.10)
        results[f"{mol_name}_oxidation_onset_10pct_V_vs_Li"] = ox_onset_10

    return results


def run_screening_batch(
    formulation_dir: Path,
    output_dir: Path,
    voltage_start: float = -2.0,
    voltage_end: float = 0.0,
    voltage_step: float = 0.05,
    n_workers: int = 4,
    max_formulations: int = None,
) -> pd.DataFrame:
    """批量筛选入口"""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 加载所有配方
    formulation_files = sorted(formulation_dir.glob("*.json"))
    formulation_files = [f for f in formulation_files if f.name != "formulation_index.json"]

    if max_formulations:
        formulation_files = formulation_files[:max_formulations]

    print(f"加载了 {len(formulation_files)} 个配方，开始筛选...")

    all_results = []
    start_time = time.time()

    for i, fpath in enumerate(formulation_files):
        if i % 50 == 0:
            elapsed = time.time() - start_time
            print(f"[{i}/{len(formulation_files)}] {fpath.stem} ... ", end="", flush=True)

        with open(fpath) as f:
            entry = json.load(f)

        result = screen_single_formulation(
            entry,
            voltage_start=voltage_start,
            voltage_end=voltage_end,
            voltage_step=voltage_step,
        )
        if result:
            all_results.append(result)

        if i % 50 == 0:
            print(f"✓ ({elapsed:.1f}s)")

    elapsed_total = time.time() - start_time
    print(f"\n筛选完成：{len(all_results)} 个配方，耗时 {elapsed_total:.1f}s")

    df = pd.DataFrame(all_results)
    return df


def generate_reports(df: pd.DataFrame, output_dir: Path):
    """生成汇总报告和排名"""
    # 保存 CSV
    summary_path = output_dir / "screening_summary.csv"
    df.to_csv(summary_path, index=False)
    print(f"汇总表已保存: {summary_path}")

    # 找还原电位最负的分子（最难还原 = 高压稳定性好）
    mol_cols = [c for c in df.columns if c.endswith("_halfwave_V_vs_Li")]
    if mol_cols:
        # 取所有分子中最正的 half-wave = 最难还原
        df["worst_reduction_V_vs_Li"] = df[mol_cols].max(axis=1)
        df["best_reduction_V_vs_Li"] = df[mol_cols].min(axis=1)

    # 排名
    if "worst_reduction_V_vs_Li" in df.columns:
        ranking = df.sort_values("worst_reduction_V_vs_Li", ascending=False)
        ranking_path = output_dir / "ranking_reduction_stability.csv"
        ranking.to_csv(ranking_path, index=False)
        print(f"还原稳定性排名已保存: {ranking_path}")

    # Markdown 报告
    report_path = output_dir / "screening_report.md"
    with open(report_path, "w") as f:
        f.write(f"# 高通量氧化还原筛选报告\n\n")
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"## 汇总\n\n")
        f.write(f"- 配方总数: {len(df)}\n")
        f.write(f"- 盐种类: {df['salt'].unique().tolist()}\n")
        if "solvents" in df.columns:
            try:
                f.write(f"- 溶剂组合: {df['solvents'].apply(lambda x: str(x)).unique().tolist()}\n")
            except Exception:
                pass
        f.write("\n")

        if "worst_reduction_V_vs_Li" in df.columns:
            f.write(f"## 还原稳定性排名（越高压越稳定）\n\n")
            f.write("| 排名 | 配方 | 最正还原电位 (V vs Li/Li+) |\n")
            f.write("|------|------|--------------------------|\n")
            top10 = df.sort_values("worst_reduction_V_vs_Li", ascending=False).head(10)
            for rank, (_, row) in enumerate(top10.iterrows(), 1):
                comp = row.get("composition", row.get("formulation_name", ""))
                f.write(f"| {rank} | {comp} | {row['worst_reduction_V_vs_Li']:.3f} |\n")

    print(f"报告已保存: {report_path}")


if __name__ == "__main__":
    import sys

    parser = argparse.ArgumentParser(description="高通量 RedoxMC 筛选")
    parser.add_argument("--formulations", type=str, default=str(REPO_ROOT / "db" / "formulations"), help="配方目录")
    parser.add_argument("--output", type=str, default=str(REPO_ROOT / "results" / "ht_screening"), help="输出目录")
    parser.add_argument("--voltage-start", type=float, default=-2.0)
    parser.add_argument("--voltage-end", type=float, default=0.0)
    parser.add_argument("--voltage-step", type=float, default=0.05)
    parser.add_argument("--max", type=int, default=None, help="最多跑多少个配方（用于测试）")
    args = parser.parse_args()

    formulation_dir = Path(args.formulations)
    output_dir = Path(args.output)

    if not formulation_dir.exists():
        print(f"配方目录不存在: {formulation_dir}")
        print("请先运行: python scripts/generate_formulations.py")
        sys.exit(1)

    df = run_screening_batch(
        formulation_dir,
        output_dir,
        voltage_start=args.voltage_start,
        voltage_end=args.voltage_end,
        voltage_step=args.voltage_step,
        max_formulations=args.max,
    )

    generate_reports(df, output_dir)
