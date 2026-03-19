# 高通量氧化还原电位筛选

## 概述

本模块基于 `redox_potential` 项目，实现电解液配方的自动化高通量氧化还原电位筛选。

**核心流程：**

```
分子数据库 (db/molecule_library.json)
    ↓
配方生成 (scripts/generate_formulations.py)
    ↓
RedoxMC 高通量筛选 (scripts/run_ht_screening.py)
    ↓
结果汇总 + 排名 (results/ht_screening/)
```

## 快速开始

### 1. 生成配方

```bash
python scripts/generate_formulations.py --output db/formulations/
```

生成 548+ 个配方，包含：
- 单一溶剂 + 各盐
- 二元/三元溶剂混合 + 各盐
- 含添加剂配方（FEC, VC, PS, LiBOB_additive）

### 2. 跑筛选

```bash
python scripts/run_ht_screening.py \
    --formulations db/formulations/ \
    --output results/ht_screening/ \
    --voltage-start -2.0 \
    --voltage-end 0.0 \
    --voltage-step 0.05
```

输出：
- `screening_summary.csv` — 所有配方汇总表
- `ranking_reduction_stability.csv` — 还原稳定性排名
- `screening_report.md` — Markdown 报告

### 3. 查看结果

```python
import pandas as pd
df = pd.read_csv("results/ht_screening/screening_summary.csv")

# 找最高压稳定的配方
mol_cols = [c for c in df.columns if c.endswith("_halfwave_V_vs_Li")]
df["worst_reduction"] = df[mol_cols].max(axis=1)
top = df.sort_values("worst_reduction", ascending=False).head(10)
print(top[["composition", "worst_reduction"]])
```

## 数据库结构

### 分子库 `db/molecule_library.json`

包含三大类分子：

| 类别 | 包含分子 |
|------|---------|
| 溶剂 (solvent_library) | EC, DMC, EMC, DEC, DME, PC |
| 锂盐 (salt_library) | LiFSI, LiPF6, LiTFSI, LiBOB, LiDFP |
| 添加剂 (additive_library) | FEC, VC, PS, LiBOB_additive |

每个分子包含：
- 基础物性（SMILES、密度、介电常数等）
- **还原参数**：E_vs_Li, E_vs_vacuum, barrier, irreversible, film_forming
- **氧化参数**：E_vs_Li, E_vs_vacuum, barrier

### 配方目录 `db/formulations/`

每个配方一个 JSON 文件，包含：
- 配方名称、成分、盐浓度、添加剂含量
- `redox_config`：可直接被 `RedoxMC` 使用的参数

## 输出字段说明

### CSV 关键列

| 列名 | 含义 |
|------|------|
| `formulation_name` | 配方唯一名称 |
| `composition` | 配方描述 |
| `salt`, `salt_conc_M` | 盐种类和浓度 |
| `additive`, `additive_wt_pct` | 添加剂种类和重量百分比 |
| `{mol}_reduction_halfwave_V_vs_Li` | 该组分还原 half-wave 电位 (V vs Li/Li⁺) |
| `{mol}_reduction_onset_10pct_V_vs_Li` | 该组分还原 10% 占有率电位 |
| `{mol}_oxidation_onset_10pct_V_vs_Li` | 该组分氧化 10% 占有率电位 |
| `worst_reduction_V_vs_Li` | 所有组分中最正的还原电位（越高越难还原=越适合高压正极） |

### 物理解释

**还原电位越正（越高）** → 组分越难被还原 → 更适合高压正极侧

**还原电位越负（越低）** → 组分越容易被还原 → 更适合负极成膜（SEI）

## 关键物理量

- 真空参考 vs Li/Li⁺ 参考偏移：**1.4 eV**（有机电解液典型值）
- 温度：**300 K**
- 计算方法：**Boltzmann 占有率模型**（无 MD 构型采样，速度快）

## 当前已知问题

1. **纯 Boltzmann 近似**：pe_solv = 0，未考虑真实溶剂化构型差异
2. **单分子近似**：多分子间相互作用通过平均场处理
3. **参数来源**：EA/IP 值来自文献，尚未用 DFT 批量计算校验

## 下一步改进方向

- [ ] 接入 PySCF DFT 计算批量生成 EA/IP 参数
- [ ] 接入真实溶剂化能校正（通过 MD 采样）
- [ ] 增加氧化侧筛选（不仅仅还原侧）
- [ ] 添加分解通道（pathway）筛选
- [ ] 可视化：电压-占有率曲线、雷达图
