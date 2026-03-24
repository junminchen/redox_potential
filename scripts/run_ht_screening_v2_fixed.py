#!/usr/bin/env python3
"""
run_ht_screening_v2.py
----------------------
高通量氧化还原双侧筛选：还原 + 氧化 + 稳定窗口。

核心改进：
- 同时计算还原侧和氧化侧 onset/halfwave
- 计算电化学稳定窗口（anodic_limit - cathodic_limit）
- 氧化参数从 molecule_library 直接读取

Usage:
    python scripts/run_ht_screening_v2.py \\
        --formulations db/formulations/ \\
        --molecule-db db/molecule_library.json \\
        --output results/ht_screening_v2/ \\
        --voltage-start -2.0 --voltage-end 2.0 --voltage-step 0.05
"""

import argparse
import json
import sys
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = REPO_ROOT / "core"
if str(CORE_DIR) not in sys.path:
    sys.path.insert(0, str(CORE_DIR))

from redox_mc import RedoxMC, RedoxParameters

# 物理常数
FARADAY_KJMOL_PER_V = 96.48533212331002
KB_KJMOL_PER_K = 8.314462e-3
EV_TO_KJMOL = 96.4853
VACUUM_TO_LI_OFFSET_EV = 1.4


def boltzmann_occupancy(
    energies: Dict[int, float],
    temperature_k: float,
    charge_states: Dict[int, float],
) -> Dict[int, float]:
    kT = KB_KJMOL_PER_K * temperature_k
    e_min = min(energies.values())
    weights = {state: np.exp(-(energy - e_min) / kT) for state, energy in energies.items()}
    total = sum(weights.values())
    return {state: w / total for state, w in weights.items()}


def crossing_voltage(voltages: np.ndarray, fractions: np.ndarray, target: float) -> float:
    for i in range(len(voltages) - 1):
        y1, y2 = fractions[i] - target, fractions[i + 1] - target
        if abs(y1) < 1e-10:
            return float(voltages[i])
        if abs(y2) < 1e-10:
            return float(voltages[i + 1])
        if y1 * y2 < 0:
            v1, v2 = voltages[i], voltages[i + 1]
            return float(v1 + (target - fractions[i]) * (v2 - v1) / (fractions[i + 1] - fractions[i]))
    return float("nan")


def build_redox_params_red(mol_cfg: dict, mol_name: str) -> RedoxParameters:
    charge_states = {int(k): float(v) for k, v in mol_cfg.get("charge_states", {}).items()}
    ea_dict = {}
    for key_str, ea_val in mol_cfg.get("electron_affinities", {}).items():
        s1, s2 = map(int, key_str.split(","))
        ea_dict[(s1, s2)] = float(ea_val)
    state_offsets_raw = mol_cfg.get("state_free_energy_offsets_kjmol", {})
    state_offsets = {int(k): float(v) for k, v in state_offsets_raw.items()} if state_offsets_raw else None
    return RedoxParameters(
        name=mol_name,
        residue_names=mol_cfg.get("residue_names", [mol_name]),
        charge_states=charge_states,
        electron_affinities=ea_dict,
        allowed_states=list(charge_states.keys()),
        partition_coefficient=float(mol_cfg.get("partition_coefficient", 0.0)),
        state_free_energy_offsets_kjmol=state_offsets,
    )


def build_redox_params_ox(mol_name: str, mol_data: dict) -> Optional[RedoxParameters]:
    """从分子库数据构建氧化侧 RedoxParameters"""
    ox = mol_data.get("oxidation", {})
    if not ox:
        return None

    E_vs_li = ox.get("E_vs_li", None)
    if E_vs_li is None:
        return None

    E_vs_vac = E_vs_li - VACUUM_TO_LI_OFFSET_EV

    # state 0 (neutral) → state +1 (oxidized)
    # energy offset = -E_vs_vac * EV_TO_KJMOL (positive E_vs_vac means oxidation is favorable)
    state_1_offset = E_vs_vac * EV_TO_KJMOL  # kJ/mol to shift from neutral to oxidized

    return RedoxParameters(
        name=mol_name,
        residue_names=[mol_name],
        charge_states={0: 0.0, 1: 1.0},
        electron_affinities={(0, 1): E_vs_vac + VACUUM_TO_LI_OFFSET_EV},
        allowed_states=[0, 1],
        partition_coefficient=0.0,
        state_free_energy_offsets_kjmol={0: 0.0, 1: state_1_offset},
    )


def compute_red_ox_for_molecule(
    params: RedoxParameters,
    voltages_vac: np.ndarray,
    voltages_li: np.ndarray,
    temperature_k: float = 300.0,
) -> dict:
    """对单个分子，计算还原和氧化在所有电压点的占有率"""
    all_states = params.allowed_states
    ref_state = 0 if 0 in all_states else min(all_states)
    ref_charge = params.charge_states[ref_state]

    # 确定还原态（更负电荷）和氧化态（更正电荷）
    neg_states = [s for s in all_states if params.charge_states[s] < ref_charge]
    pos_states = [s for s in all_states if params.charge_states[s] > ref_charge]

    red_fracs = []
    ox_fracs = []

    for v_vac in voltages_vac:
        mc = RedoxMC(params, temperature_k=temperature_k, voltage_v=float(v_vac))
        energies = {s: mc.compute_state_energy(s, pe_nonbonded=0.0) for s in all_states}
        occ = boltzmann_occupancy(energies, temperature_k, params.charge_states)

        if neg_states:
            red_frac = sum(occ.get(s, 0.0) for s in neg_states)
        else:
            red_frac = 0.0

        if pos_states:
            ox_frac = sum(occ.get(s, 0.0) for s in pos_states)
        else:
            ox_frac = 0.0

        red_fracs.append(red_frac)
        ox_fracs.append(ox_frac)

    red_fracs = np.array(red_fracs)
    ox_fracs = np.array(ox_fracs)

    return {
        "red_fracs": red_fracs,
        "ox_fracs": ox_fracs,
        "neg_states": neg_states,
        "pos_states": pos_states,
    }


def screen_single_formulation(
    entry: dict,
    mol_library: dict,
    voltage_start: float = -2.0,
    voltage_end: float = 5.0,   # extended to capture oxidation onsets up to ~6.4V vs Li
    voltage_step: float = 0.05,
    temperature_k: float = 300.0,
) -> dict:
    """
    对单个配方：
    1. 还原侧：用 formulation JSON 里的 redox_config（负电荷态）
    2. 氧化侧：从 molecule_library 查氧化参数（正电荷态）
    """
    red_config = entry.get("redox_config", {})
    mol_configs = red_config.get("molecules", {})

    voltages_vac = np.arange(voltage_start, voltage_end + voltage_step / 2, voltage_step)
    voltages_li = voltages_vac + VACUUM_TO_LI_OFFSET_EV

    results = {
        "formulation_name": entry.get("formulation_name", ""),
        "composition": entry.get("composition", ""),
        "salt": entry.get("salt_name", entry.get("salt", "")),
        "salt_conc_M": entry.get("salt_conc_M", 0.0),
        "additive": entry.get("additive_name", entry.get("additive", "")),
        "additive_wt_pct": entry.get("additive_wt_pct", 0.0),
    }

    if not mol_configs:
        return results

    all_red_onsets = []
    all_ox_onsets = []

    for mol_name, mol_cfg in mol_configs.items():
        # === 还原侧 ===
        try:
            red_params = build_redox_params_red(mol_cfg, mol_name)
            red_data = compute_red_ox_for_molecule(red_params, voltages_vac, voltages_li, temperature_k)
            red_fracs = red_data["red_fracs"]
            neg_states = red_data["neg_states"]

            if len(neg_states) > 0:
                # 还原占有率 = 1 - frac_in_neutral_or_better
                # 这里 frac(还原) = sum of negative charge states
                # f_red(V) = fraction in reduced state (neg charge)
                # At low V (very negative): red_frac → 1 (fully reduced)
                # At high V (positive): red_frac → 0 (oxidized/neutral)
                # We want onset where red_frac drops below 0.9 (i.e., 10% oxidized)
                red_onset = crossing_voltage(voltages_li, red_fracs, 0.10)
                red_halfwave = crossing_voltage(voltages_li, red_fracs, 0.50)
                all_red_onsets.append(red_onset)
            else:
                red_onset = float("nan")
                red_halfwave = float("nan")

            results[f"{mol_name}_red_onset_V_vs_Li"] = round(red_onset, 3) if not np.isnan(red_onset) else np.nan
            results[f"{mol_name}_red_halfwave_V_vs_Li"] = round(red_halfwave, 3) if not np.isnan(red_halfwave) else np.nan
        except Exception as e:
            results[f"{mol_name}_red_onset_V_vs_Li"] = np.nan
            results[f"{mol_name}_red_halfwave_V_vs_Li"] = np.nan

        # === 氧化侧 ===
        # 找该分子在 molecule_library 中的氧化数据
        mol_lib_entry = None
        for cat in ["solvent_library", "salt_library", "additive_library"]:
            if mol_name in mol_library.get(cat, {}):
                mol_lib_entry = mol_library[cat][mol_name]
                break

        if mol_lib_entry:
            try:
                ox_params = build_redox_params_ox(mol_name, mol_lib_entry)
                if ox_params:
                    ox_data = compute_red_ox_for_molecule(ox_params, voltages_vac, voltages_li, temperature_k)
                    ox_fracs = ox_data["ox_fracs"]
                    pos_states = ox_data["pos_states"]

                    if len(pos_states) > 0:
                        # 氧化占有率 = fraction in positive charge states
                        # At low V: ox_frac → 0 (no oxidation)
                        # At high V: ox_frac → 1 (fully oxidized)
                        # onset: when ox_frac = 0.10
                        ox_onset = crossing_voltage(voltages_li, ox_fracs, 0.10)
                        ox_halfwave = crossing_voltage(voltages_li, ox_fracs, 0.50)
                        all_ox_onsets.append(ox_onset)
                    else:
                        ox_onset = float("nan")
                        ox_halfwave = float("nan")

                    results[f"{mol_name}_ox_onset_V_vs_Li"] = round(ox_onset, 3) if not np.isnan(ox_onset) else np.nan
                    results[f"{mol_name}_ox_halfwave_V_vs_Li"] = round(ox_halfwave, 3) if not np.isnan(ox_halfwave) else np.nan
            except Exception as e:
                results[f"{mol_name}_ox_onset_V_vs_Li"] = np.nan
                results[f"{mol_name}_ox_halfwave_V_vs_Li"] = np.nan

    # 配方级稳定性窗口
    valid_red = [x for x in all_red_onsets if not np.isnan(x)]
    valid_ox = [x for x in all_ox_onsets if not np.isnan(x)]

    if valid_red:
        # cathodic limit = 最正的还原 onset（最难还原 = 负极最稳定）
        results["cathodic_limit_V_vs_Li"] = round(max(valid_red), 3)
    if valid_ox:
        # anodic limit = 最负的氧化 onset（最难氧化 = 正极最稳定）
        results["anodic_limit_V_vs_Li"] = round(min(valid_ox), 3)
    if valid_red and valid_ox:
        results["stability_window_V"] = round(min(valid_ox) - max(valid_red), 3)

    return results


def run_screening_batch(
    formulation_dir: Path,
    mol_library_path: Path,
    output_dir: Path,
    voltage_start: float = -2.0,
    voltage_end: float = 5.0,   # extended for oxidation side
    voltage_step: float = 0.05,
    max_formulations: int = None,
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(mol_library_path) as f:
        mol_library = json.load(f)

    formulation_files = sorted(formulation_dir.glob("*.json"))
    formulation_files = [f for f in formulation_files if f.name not in ("formulation_index.json",)]

    if max_formulations:
        formulation_files = formulation_files[:max_formulations]

    print(f"加载了 {len(formulation_files)} 个配方，开始双侧筛选（还原+氧化）...")
    all_results = []
    start_time = time.time()

    for i, fpath in enumerate(formulation_files):
        if i % 50 == 0:
            elapsed = time.time() - start_time
            if i > 0:
                eta = elapsed / i * (len(formulation_files) - i)
                print(f"[{i}/{len(formulation_files)}] {fpath.stem} ... ETA={eta:.0f}s ✓")
            else:
                print(f"[{i}/{len(formulation_files)}] {fpath.stem} ...", end="", flush=True)

        with open(fpath) as f:
            entry = json.load(f)

        result = screen_single_formulation(entry, mol_library, voltage_start, voltage_end, voltage_step)
        if result:
            all_results.append(result)

        if i % 50 == 0 and i > 0:
            pass  # already printed

    elapsed_total = time.time() - start_time
    print(f"\n筛选完成：{len(all_results)} 个配方，耗时 {elapsed_total:.1f}s")

    df = pd.DataFrame(all_results)
    return df


def generate_outputs(df: pd.DataFrame, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    # CSV
    df.to_csv(output_dir / "screening_summary.csv", index=False)
    print(f"汇总: {output_dir / 'screening_summary.csv'} ({len(df.columns)} 列)")

    # 窗口
    if "stability_window_V" in df.columns:
        win = df.sort_values("stability_window_V", ascending=False)
        win.to_csv(output_dir / "stability_windows.csv", index=False)

        # 还原排名
        if "cathodic_limit_V_vs_Li" in df.columns:
            red_rank = df.sort_values("cathodic_limit_V_vs_Li", ascending=False)
            red_rank.to_csv(output_dir / "ranking_reduction.csv", index=False)

        # 氧化排名
        if "anodic_limit_V_vs_Li" in df.columns:
            ox_rank = df.sort_values("anodic_limit_V_vs_Li", ascending=True)
            ox_rank.to_csv(output_dir / "ranking_oxidation.csv", index=False)

    # Markdown 报告
    report_path = output_dir / "screening_report.md"
    with open(report_path, "w") as f:
        f.write(f"# 高通量氧化还原双侧筛选报告\n\n")
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"## 汇总\n\n")
        f.write(f"- 配方总数: {len(df)}\n")
        f.write(f"- 盐种类: {df['salt'].dropna().unique().tolist()}\n\n")

        if "stability_window_V" in df.columns:
            f.write(f"## 电化学稳定窗口排名（越宽越好）\n\n")
            f.write("| 排名 | 配方 | 窗口 (V) | 阳极极限 (V vs Li/Li⁺) | 阴极极限 (V vs Li/Li⁺) |\n")
            f.write("|------|------|----------|------------------------|------------------------|\n")
            for rank, (_, row) in enumerate(df.sort_values("stability_window_V", ascending=False).head(15).iterrows(), 1):
                comp = row.get("composition", row.get("formulation_name", ""))
                win = row.get("stability_window_V", np.nan)
                an = row.get("anodic_limit_V_vs_Li", np.nan)
                cat = row.get("cathodic_limit_V_vs_Li", np.nan)
                f.write(f"| {rank} | {comp} | {win:.2f} | {an:.2f} | {cat:.2f} |\n")

        if "anodic_limit_V_vs_Li" in df.columns:
            f.write(f"\n## 氧化稳定性排名（越负越耐高压正极）\n\n")
            f.write("| 排名 | 配方 | 阳极极限 (V vs Li/Li⁺) |\n")
            f.write("|------|------|--------------------------|\n")
            for rank, (_, row) in enumerate(df.sort_values("anodic_limit_V_vs_Li", ascending=True).head(15).iterrows(), 1):
                comp = row.get("composition", row.get("formulation_name", ""))
                an = row.get("anodic_limit_V_vs_Li", np.nan)
                f.write(f"| {rank} | {comp} | {an:.2f} |\n")

        if "cathodic_limit_V_vs_Li" in df.columns:
            f.write(f"\n## 还原稳定性排名（越正越耐高压负极）\n\n")
            f.write("| 排名 | 配方 | 阴极极限 (V vs Li/Li⁺) |\n")
            f.write("|------|------|--------------------------|\n")
            for rank, (_, row) in enumerate(df.sort_values("cathodic_limit_V_vs_Li", ascending=False).head(15).iterrows(), 1):
                comp = row.get("composition", row.get("formulation_name", ""))
                cat = row.get("cathodic_limit_V_vs_Li", np.nan)
                f.write(f"| {rank} | {comp} | {cat:.2f} |\n")

    print(f"报告: {report_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="高通量 RedoxMC 双侧筛选")
    parser.add_argument("--formulations", type=str, default=str(REPO_ROOT / "db" / "formulations"))
    parser.add_argument("--molecule-db", type=str, default=str(REPO_ROOT / "db" / "molecule_library.json"))
    parser.add_argument("--output", type=str, default=str(REPO_ROOT / "results" / "ht_screening_v2"))
    parser.add_argument("--voltage-start", type=float, default=-2.0)
    parser.add_argument("--voltage-end", type=float, default=2.0)
    parser.add_argument("--voltage-step", type=float, default=0.05)
    parser.add_argument("--max", type=int, default=None)
    args = parser.parse_args()

    df = run_screening_batch(
        Path(args.formulations),
        Path(args.molecule_db),
        Path(args.output),
        args.voltage_start,
        args.voltage_end,
        args.voltage_step,
        args.max,
    )
    generate_outputs(df, Path(args.output))
