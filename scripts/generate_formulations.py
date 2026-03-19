#!/usr/bin/env python3
"""
generate_formulations.py
-----------------------
从分子数据库自动组合生成电解液配方配置文件。

支持：
- 溶剂组合（单一溶剂、二元混合、三元混合）
- 盐浓度梯度（0.5M - 5M）
- 添加剂筛选（无添加剂、低浓度、高浓度）

Usage:
    python scripts/generate_formulations.py --output db/formulations/
"""

import argparse
import json
import itertools
import sys
from pathlib import Path
from typing import Dict, List, Any

REPO_ROOT = Path(__file__).resolve().parents[1]
DB_DIR = REPO_ROOT / "db"
sys.path.insert(0, str(REPO_ROOT / "core"))


# 默认浓度规范
DEFAULT_CONCENTRATIONS_M = {
    "LiFSI": [1.0, 2.0, 4.0, 6.0],
    "LiPF6": [0.5, 1.0, 1.5],
    "LiTFSI": [1.0, 3.0, 5.0],
    "LiBOB": [0.1, 0.5, 1.0],
    "LiDFP": [0.1, 0.5, 1.0],
}

DEFAULT_ADDITIVE_WEIGHT_PCT = {
    "FEC": [0.0, 2.0, 5.0, 10.0],
    "VC": [0.0, 2.0, 5.0],
    "PS": [0.0, 1.0, 2.0],
    "LiBOB_additive": [0.0, 1.0, 2.0],
}

DEFAULT_SOLVENT_RATIOS = {
    "EC:DMC": [(3, 7), (5, 5), (7, 3)],
    "EC:EMC": [(3, 7), (1, 1)],
    "EC:DEC": [(3, 7), (5, 5)],
    "EC:PC": [(5, 5), (7, 3)],
    "EC:DME": [(7, 3), (5, 5), (3, 7)],
    "DMC:EMC": [(5, 5), (3, 7)],
}

SOLVENT_CATEGORIES = {
    "cyclic_carbonate": ["EC", "PC"],
    "linear_carbonate": ["DMC", "EMC", "DEC"],
    "ether": ["DME"],
}

SOLVENT_BLEND_COMBOS = [
    ["EC", "DMC"],
    ["EC", "EMC"],
    ["EC", "DEC"],
    ["EC", "PC"],
    ["EC", "DME"],
    ["DMC", "EMC"],
    ["DMC", "DME"],
    ["EMC", "DEC"],
    ["EC", "DMC", "EMC"],
]


def load_molecule_library() -> Dict:
    with open(DB_DIR / "molecule_library.json") as f:
        return json.load(f)


def build_solvent_blend_key(solvent_names: List[str], ratio: tuple) -> str:
    """生成溶剂混合的键名，如 EC:DMC_3:7"""
    return "_".join(solvent_names) + "_" + ":".join(str(x) for x in ratio)


def build_redox_config_from_molecules(
    mol_library: Dict,
    solvent_names: List[str],
    salt_name: str,
    additive_name: str = None,
    additive_wt_pct: float = 0.0,
    ratio: tuple = None,
) -> Dict:
    """
    从分子列表构建 redox config（只包含可还原组分）。
    溶剂混合时，按比例加权计算 state_free_energy_offset。
    """
    offset_vac_to_li = 1.4  # eV

    molecules = {}

    # 溶剂
    if len(solvent_names) == 1:
        solvents = {solvent_names[0]: 1.0}
    else:
        # ratio 是 (a, b) 表示 a:b 的重量比
        total = sum(ratio)
        solvents = {name: r / total for name, r in zip(solvent_names, ratio)}

    for solvent_name, wt_frac in solvents.items():
        mol_data = mol_library["solvent_library"][solvent_name]
        red = mol_data["reduction"]
        E_vac = red["E_vs_vacuum"]
        barrier = red["barrier_kjmol"]

        molecules[solvent_name] = {
            "residue_names": [solvent_name],
            "charge_states": {"0": 0.0, "-1": -1.0},
            "electron_affinities": {"0,-1": E_vac + offset_vac_to_li},
            "state_free_energy_offsets_kjmol": {
                "0": 0.0,
                "-1": -E_vac * 96.4853,
            },
            "allowed_states": [0, -1],
            "partition_coefficient": 0.0,
        }

    # 盐（作为还原组分处理，阴离子）
    if salt_name in mol_library["salt_library"]:
        salt_data = mol_library["salt_library"][salt_name]
        red = salt_data["reduction"]
        E_vac = red["E_vs_vacuum"]
        molecules[salt_name] = {
            "residue_names": [salt_name],
            "charge_states": {"0": 0.0, "-1": -1.0},
            "electron_affinities": {"0,-1": E_vac + offset_vac_to_li},
            "state_free_energy_offsets_kjmol": {
                "0": 0.0,
                "-1": -E_vac * 96.4853,
            },
            "allowed_states": [0, -1],
            "partition_coefficient": 0.0,
        }

    # 添加剂
    if additive_name and additive_wt_pct > 0:
        if additive_name in mol_library["additive_library"]:
            add_data = mol_library["additive_library"][additive_name]
            red = add_data["reduction"]
            E_vac = red["E_vs_vacuum"]
            molecules[additive_name] = {
                "residue_names": [additive_name],
                "charge_states": {"0": 0.0, "-1": -1.0},
                "electron_affinities": {"0,-1": E_vac + offset_vac_to_li},
                "state_free_energy_offsets_kjmol": {
                    "0": 0.0,
                    "-1": -E_vac * 96.4853,
                },
                "allowed_states": [0, -1],
                "partition_coefficient": 0.0,
            }

    composition_str = f"{salt_name} in "
    if len(solvent_names) == 1:
        composition_str += solvent_names[0]
    else:
        composition_str += ":".join(solvent_names) + f" ({':'.join(str(r) for r in ratio)})"
    if additive_name and additive_wt_pct > 0:
        composition_str += f" + {additive_wt_pct}wt% {additive_name}"

    config = {
        "metadata": {
            "electrolyte_composition": composition_str,
            "temperature_k": 300.0,
            "reference_electrode": "Li/Li+",
            "vacuum_to_li_offset_ev": offset_vac_to_li,
        },
        "molecules": molecules,
    }
    return config


def generate_formulation_entry(
    mol_library: Dict,
    solvent_names: List[str],
    salt_name: str,
    salt_conc_m: float,
    additive_name: str = None,
    additive_wt_pct: float = 0.0,
    ratio: tuple = None,
) -> Dict:
    """生成单个配方条目（对应一个 formulation JSON）"""
    composition = f"{salt_name} {salt_conc_m}M in "
    if len(solvent_names) == 1:
        composition += solvent_names[0]
        blend_key = solvent_names[0]
    else:
        ratio_str = ":".join(str(r) for r in ratio)
        composition += ":".join(solvent_names) + f" ({ratio_str})"
        blend_key = "_".join(solvent_names) + "_" + ratio_str

    additive_str = f"_{additive_name}_{additive_wt_pct}wt" if additive_name and additive_wt_pct > 0 else ""
    formulation_name = f"{salt_name}_{salt_conc_m}M_{blend_key}{additive_str}"

    redox_config = build_redox_config_from_molecules(
        mol_library, solvent_names, salt_name, additive_name, additive_wt_pct, ratio
    )

    return {
        "formulation_name": formulation_name,
        "composition": composition,
        "salt_name": salt_name,
        "salt_conc_M": salt_conc_m,
        "solvents": solvent_names,
        "solvent_ratio": list(ratio) if ratio else [1.0],
        "additive_name": additive_name,
        "additive_wt_pct": additive_wt_pct,
        "redox_config": redox_config,
    }


def generate_all_formulations(
    mol_library: Dict,
    solvent_combos: List = None,
    salt_concs: Dict = None,
    additive_combos: List = None,
    output_dir: Path = None,
) -> List[Dict]:
    """生成全套配方"""
    all_entries = []

    salt_concs = salt_concs or DEFAULT_CONCENTRATIONS_M
    additive_combos = additive_combos or list(DEFAULT_ADDITIVE_WEIGHT_PCT.keys())

    # 1. 纯溶剂 + 盐
    solvents = mol_library["solvent_library"]
    salts = mol_library["salt_library"]

    # 单一溶剂筛选
    for salt_name, concs in salt_concs.items():
        for conc in concs:
            for solvent_name in solvents:
                entry = generate_formulation_entry(
                    mol_library,
                    [solvent_name],
                    salt_name,
                    conc,
                )
                all_entries.append(entry)

    # 2. 溶剂混合 + 盐
    if solvent_combos is None:
        solvent_combos = SOLVENT_BLEND_COMBOS

    for combo in solvent_combos:
        if len(combo) == 2:
            solvent_key = "_".join(combo)
            if solvent_key in DEFAULT_SOLVENT_RATIOS:
                ratios = DEFAULT_SOLVENT_RATIOS[solvent_key]
            else:
                ratios = [(1, 1), (3, 7)]
        elif len(combo) == 3:
            ratios = [(5, 3, 2)]
        else:
            ratios = [(1,)]

        for salt_name, concs in salt_concs.items():
            for conc in concs:
                for ratio in ratios:
                    entry = generate_formulation_entry(
                        mol_library,
                        combo,
                        salt_name,
                        conc,
                        ratio=ratio,
                    )
                    all_entries.append(entry)

    # 3. 加添加剂
    additive_entries = []
    for entry in all_entries[:20]:  # 限制数量，只对top溶剂组合加添加剂
        for add_name in additive_combos:
            add_weights = DEFAULT_ADDITIVE_WEIGHT_PCT[add_name]
            for wt in add_weights:
                if wt == 0.0:
                    continue
                new_entry = generate_formulation_entry(
                    mol_library,
                    entry["solvents"],
                    entry["salt_name"],
                    entry["salt_conc_M"],
                    additive_name=add_name,
                    additive_wt_pct=wt,
                    ratio=tuple(entry["solvent_ratio"]) if len(entry["solvent_ratio"]) > 1 else None,
                )
                additive_entries.append(new_entry)

    all_entries.extend(additive_entries)

    # 保存
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        for entry in all_entries:
            fname = entry["formulation_name"].replace(" ", "_").replace(":", "_").replace(".", "p")
            with open(output_dir / f"{fname}.json", "w") as f:
                json.dump(entry, f, indent=2, ensure_ascii=False)

        # 保存汇总 index
        index_path = output_dir / "formulation_index.json"
        with open(index_path, "w") as f:
            # 不保存完整的 redox_config，只保留 metadata
            index = [{k: v for k, v in e.items() if k != "redox_config"} for e in all_entries]
            json.dump(index, f, indent=2, ensure_ascii=False)
        print(f"生成了 {len(all_entries)} 个配方，已保存到 {output_dir}")

    return all_entries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="生成电解液配方配置")
    parser.add_argument("--output", type=str, default=str(REPO_ROOT / "db" / "formulations"), help="输出目录")
    args = parser.parse_args()

    mol_library = load_molecule_library()
    output_dir = Path(args.output)
    generate_all_formulations(mol_library, output_dir=output_dir)
