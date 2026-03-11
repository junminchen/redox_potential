# 氧化还原电势计算工作流 (RedoxMC + Voltage Sweep)

## 概述

本工作流实现了 PRD.md 中的 **Task 3、4、5**，用于在 1M LiPF6 EC/DMC 电解质中计算任意分子的氧化还原电势。

核心思路：
1. **恒电位MD** (OpenMM 8.4 ConstantPotentialForce) → 固定电极电位下的分子动力学
2. **RedoxMC采样** (Metropolis MC) → 分子电荷态分布在各电位下的统计
3. **电位扫描** (Voltage Sweep) → 在一系列电位下重复以上步骤，生成模拟CV曲线

## 工作流图

```
Config JSON → OPLSMDSimulation (恒电位MD)
    ↓
    ├→ 电极电荷平衡 (Q_cathode, Q_anode)
    ├→ 电解质分子动力学
    └→ 能量记录

对每个电位V扫描：
    V_start → V_start+ΔV → ... → V_end

输出：
    - 各电位下的电极电荷 (Q vs V)
    - 各电位下的分子电荷态占有率 (fraction_reduced vs V)
    - 模拟"电流" (dQ/dV)
    - 半波电位 E_1/2 与实验对比
```

## 文件说明

| 文件 | 功能 | 创建时间 |
|------|------|--------|
| `core/redox_mc.py` | RedoxMC 类，RedoxParameters dataclass | Task 3 |
| `core/voltage_sweep.py` | 电位扫描主模块，CV分析器 | Task 4 |
| `configs/redox/config_redox.json` | 分子氧化还原参数配置 | Task 5 |
| `scripts/cli/run_voltage_sweep.py` | 批量扫描驱动脚本 | 辅助 |
| `docs/methods/REDOX_WORKFLOW.md` | 本文档 | 文档 |

## 快速开始

### 1. 准备输入文件

需要以下文件（已存在于项目中）：
```bash
# MD 配置
configs/md/config_opls.json              # OPLS 力场配置（已有）

# 分子结构
structures/systems/start_with_electrodes.pdb     # 电极+电解质体系（已有）
structures/pdb_bank/EC.pdb              # EC 分子（已有）
structures/pdb_bank/DMC.pdb             # DMC 分子（已有）

# 新增：红氧参数
configs/redox/config_redox.json            # 红氧分子参数（已生成）
```

### 2. 编辑 configs/redox/config_redox.json

指定要计算的分子及其氧化还原参数：

```json
{
  "molecules": {
    "YourMolecule": {
      "residue_names": ["RES1", "RES2"],
      "charge_states": {"0": 0.0, "-1": -1.0},
      "electron_affinities": {"0,-1": 1.53},
      "allowed_states": [0, -1],
      "partition_coefficient": -2.5
    }
  }
}
```

**关键参数解释：**

- `charge_states`: 分子的允许电荷态（以 e 为单位）
  - 例：`"0": 0.0` 表示中性，`"-1": -1.0` 表示单负离子

- `electron_affinities`: 电子转移的气相自由能 (eV)
  - `"0,-1"`: 从中性到 -1 价态的电子亲合能
  - 正值：接受电子有利 (减少)
  - 负值：失去电子有利 (氧化)

- `partition_coefficient`: 溶剂化自由能修正 (kJ/mol)
  - 从 SCRF/PCM 等溶解模型得出
  - 值越负，该状态在溶液中越稳定

### 3. 运行电位扫描

```bash
# 激活 CUDA 环境
source activate omm84_cuda_pip

# 标准运行（默认参数）
python scripts/cli/run_voltage_sweep.py \
    --config configs/md/config_opls.json \
    --redox-config configs/redox/config_redox.json \
    --pdb structures/systems/start_with_electrodes.pdb \
    --output-dir results/voltage_sweep_example

# 带实验数据对比
python scripts/cli/run_voltage_sweep.py \
    --config configs/md/config_opls.json \
    --redox-config configs/redox/config_redox.json \
    --pdb structures/systems/start_with_electrodes.pdb \
    --output-dir results/voltage_sweep_example \
    --lsv-peak-v -1.7  # 实验 LSV 峰值 (V vs Li/Li+)

# 自定义电位范围和步数
python scripts/cli/run_voltage_sweep.py \
    --config configs/md/config_opls.json \
    --redox-config configs/redox/config_redox.json \
    --pdb structures/systems/start_with_electrodes.pdb \
    --output-dir results/voltage_sweep_example \
    --v-start -4.5 \
    --v-end -2.5 \
    --v-step 0.1 \
    --equil-steps 3000 \
    --sample-steps 5000 \
    --platform CUDA
```

### 4. 输出结果分析

扫描完成后在 `results/voltage_sweep_example/` 目录下生成：

```
results/voltage_sweep_example/
├── voltage_sweep_results.csv     # 完整数据表
├── voltage_sweep_checkpoint.csv  # 检查点文件
├── Vm4.500/
│   ├── config.json              # 该电位的实际配置
│   ├── electrode_charges.log    # 电极电荷时间序列
│   ├── final_opls.pdb           # 该电位的末态结构
│   ├── md_opls.dcd              # 该电位的轨迹
│   └── md_opls.log              # 该电位的能量日志
└── Vm4.400/ ...                 # 其他电位点同样布局
```

其中 `voltage_v` 表示总电池电压（真空参考），恒电势边界条件实际使用
`cathode = +voltage_v/2` 和 `anode = -voltage_v/2`。

**结果 CSV 格式：**

```
voltage_v, q_cathode_mean, q_anode_mean, q_cathode_std, ...
-4.500,    2.341,         -2.337,        0.045,        ...
-4.400,    2.298,         -2.294,        0.052,        ...
...
-2.500,    1.876,         -1.873,        0.031,        ...
```

### 5. 分析 CV 曲线

使用 Python 绘制和分析结果：

```python
import pandas as pd
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, "core"); from voltage_sweep import RedoxPotentialAnalyzer

# 加载结果
df = pd.read_csv('results/voltage_sweep_example/voltage_sweep_results.csv')

# 绘制电极电荷 vs 电位
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ax1.plot(df['voltage_v'], df['q_cathode_mean'], 'o-', label='Cathode')
ax1.plot(df['voltage_v'], -df['q_anode_mean'], 's-', label='Anode (negated)')
ax1.set_xlabel('Total Voltage (V vs vacuum)')
ax1.set_ylabel('Charge (e)')
ax1.legend()
ax1.grid(True)

# 模拟电流（数值微分）
ax2.plot(df['voltage_v'][:-1],
         pd.Series(df['q_cathode_mean']).diff().values[1:] / 0.1,
         'o-', label='dQ/dV')
ax2.set_xlabel('Total Voltage (V vs vacuum)')
ax2.set_ylabel('dQ/dV (e/V)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('results/voltage_sweep_example/cv_curve.png', dpi=300)
plt.show()

# 提取半波电位
for mol_name in ['EC', 'DMC']:  # 或其他分子
    e_half = RedoxPotentialAnalyzer.find_half_wave_potential(
        df, mol_name, voltage_col='voltage_v'
    )
    e_half_li = RedoxPotentialAnalyzer.convert_to_li_reference(e_half)
    print(f"{mol_name} E_1/2 = {e_half_li:.3f} V vs Li/Li+")
```

## 物理模型

### 大正则系综中的电子转移

在固定电极电位下，分子的多个电荷态以概率分布存在：

$$P(state) \propto \exp\left[-\frac{E(state)}{k_B T}\right]$$

其中总能量为：

$$E(state) = E_{electrode}(state) + E_{solv}^{MD}(state)$$

- $E_{electrode} = z \cdot F \cdot V$ 其中 $z$ 是电荷，$V$ 是电位
- $E_{solv}^{MD}$ 来自 NonbondedForce（包含溶剂化和构象效应）

MC 接受准则：

$$P_{accept} = \min\left(1, \exp\left[-\frac{\Delta E}{k_B T}\right]\right)$$

### 与实验 CV 的对应

模拟的"电流"反映分子在该电位下的氧化还原活性：

$$I_{sim}(V) \propto \frac{d\langle n_{reduced} \rangle}{dV}$$

其中 $\langle n_{reduced} \rangle$ 是还原态（最低价态）的占有率。

实验 LSV 峰值位置对应于：
- $\langle n_{reduced} \rangle = 0.5$（半波电位）
- 此时 $E_{electrode} = -\frac{EA_{gas} + \Delta G_{solv}}{F}$

## 常见问题和调试

### Q1: 模拟结果与实验 LSV 偏差大

**检查清单：**
1. 确认 `electron_affinities` 值来自相同的 DFT 泛函（B3LYP 等）
2. 核对 `partition_coefficient`（SCRF 模型的溶剂化能）
3. 是否需要加入库伦库（长程相互作用）修正？
4. 运行时间是否足够（ensure convergence of MC sampling）

### Q2: 电荷平衡错误

如果输出显示 $|Q_{cathode} + Q_{anode}| > 0.01e$：
- 检查 `configs/md/config_opls.json` 中的电极定义
- 验证 PDB 中电极链的原子数一致

### Q3: MC 接受率过低

接受率 < 10% 通常表示：
- 电位已超出分子的化学势窗口
- 需要增加 `sample_steps_per_voltage`
- 或降低 `v_step` 的大小

## 性能优化 (CUDA)

在 CUDA 环境下运行时的加速建议：

```bash
# 设置 CUDA 环境变量
export CUDA_VISIBLE_DEVICES=0         # 使用GPU 0
export OPENMM_PLUGIN_DIR=$CONDA_PREFIX/lib/plugins  # 确保插件路径

# 在脚本中确认平台
python -c "import openmm as mm; \
    print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])"
```

## 文献参考

1. **Reed-Madden 常电位方法**：
   - Reed, S. K., Madden, P. A., Papadopoulos, A. G. (2007).
     "A new algorithm for unrestricted constant potential molecular dynamics"
     *J. Chem. Phys.* 126, 084704

2. **MC 采样与 CV 对齐**：
   - Meissner, R., Cucinotta, C. S., Gogoll, A., Ballauff, M., Jurado-Sánchez, B., Schöll, E. (2018).
     "Molecular dynamics simulations of aqueous voltammetry at electrode surfaces"
     *J. Chem. Theory Comput.* 14, 6307–6319

3. **OpenMM 8.4 迁移指南**：
   - https://openmm.org/documentation/latest/migration/

## 许可证与引用

本工作流基于 OpenMM 8.4 和 OPLS 力场。如发表，请引用：
- OpenMM: Eastman, P., et al. (2019) *PLOS Comp. Biol.* 15(7), e1007659
- OPLS-AA: Jorgensen, W. L., Maxwell, D. S., et al. (1996) *J. Am. Chem. Soc.* 118, 11225

---

**最后更新**: 2026-03-11
**工作流版本**: v1.0 (Task 3, 4, 5)
