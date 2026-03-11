# 快速参考卡片 - 氧化还原电势计算框架

## 🚀 30 秒快速开始

```bash
# 1. 激活环境
source activate omm84_cuda_pip

# 2. 运行扫描（默认参数）
python run_voltage_sweep.py --config config_opls.json \
    --redox-config config_redox.json --pdb start_with_electrodes.pdb

# 3. 查看结果
cat sweep_results/voltage_sweep_results.csv
```

## 📋 关键命令

| 命令 | 用途 |
|------|------|
| `python test_redox_framework.py` | 验证框架完整性 ✅ |
| `python run_voltage_sweep.py --help` | 查看所有选项 |
| `python -c "from redox_mc import *; help(RedoxMC)"` | API 帮助 |

## ⚙️ 常用参数

### 电位扫描范围
```
--v-start -4.5      # 起始（V vs vacuum）
--v-end -2.5        # 终止
--v-step 0.1        # 步长（V）
```
**典型范围：** -5 到 0 V vs vacuum ≈ -3.6 到 +1.4 V vs Li/Li+

### 采样时间
```
--equil-steps 5000      # 每电位平衡步数
--sample-steps 10000    # 每电位采样步数
--platform CUDA         # GPU 加速（OpenCL 可用）
```
**经验值：**
- 快速测试：2K + 2K（几分钟）
- 标准：5K + 10K（~5-10 分钟 per point）
- 高精度：10K + 20K（~15-20 分钟 per point）

## 📊 输出文件

```
sweep_results/
├── voltage_sweep_results.csv    ← 主结果表
├── voltage_sweep_checkpoint.csv ← 中间检查点
├── Vm4.500/
│   ├── config.json
│   ├── electrode_charges.log
│   ├── final_opls.pdb
│   ├── md_opls.dcd
│   └── md_opls.log
└── Vm4.400/ ...                 ← 其他电压点同样布局
```

### 结果 CSV 格式
```
voltage_v, q_cathode_mean, q_anode_mean, q_cathode_std, ...
-4.500,    2.341,         -2.337,        0.045
-4.400,    2.298,         -2.294,        0.052
...
```

**关键列：**
- `voltage_v` - 总电池电压 (V vs vacuum)
- `q_cathode_mean` - 阴极平均电荷 (e)
- `q_anode_mean` - 阳极平均电荷 (e)

**电压定义：**
- `voltage_v` 是施加在整个双电极体系上的总电压
- 实际传给恒电势求解器的是 `cathode = +voltage_v/2`, `anode = -voltage_v/2`

## 🔬 物理常数

| 常数 | 值 | 说明 |
|------|-----|------|
| Faraday | 96.485 kJ/mol/V | 电荷-能量转换 |
| Boltzmann | 8.314e-3 kJ/(mol·K) | 热能 |
| T | 300 K | 标准温度 |
| kT(300K) | 2.494 kJ/mol | 热能水平 |
| Li/Li+ offset | 1.4 V | 真空→Li参考 |

## 📈 数据分析（Python）

```python
import pandas as pd
from voltage_sweep import RedoxPotentialAnalyzer

# 加载结果
df = pd.read_csv('sweep_results/voltage_sweep_results.csv')

# 绘制 CV 曲线
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(df['voltage_v'], df['q_cathode_mean'], 'o-')
ax.set_xlabel('Total Voltage (V vs vacuum)')
ax.set_ylabel('Cathode Charge (e)')
plt.savefig('cv_curve.png')

# 提取半波电位
e_half = RedoxPotentialAnalyzer.find_half_wave_potential(df, 'your_molecule')
e_half_li = RedoxPotentialAnalyzer.convert_to_li_reference(e_half)
print(f"E_1/2 = {e_half_li:.3f} V vs Li/Li+")

# 与实验对比
result = RedoxPotentialAnalyzer.compare_with_lsv(
    df, 'your_molecule',
    lsv_peak_potential_v_vs_li=-1.7  # 实验值
)
print(f"差异: {result['difference_mv']:.1f} mV")
```

## 🛠️ 配置文件编辑

### 添加新分子
编辑 `config_redox.json`：

```json
"YOUR_MOLECULE": {
  "residue_names": ["RES1", "RES2"],
  "charge_states": {"0": 0.0, "-1": -1.0, "-2": -2.0},
  "electron_affinities": {
    "0,-1": 1.53,      # 气相电子亲合能 (eV)
    "-1,-2": -0.2      # 第二电子转移 (eV)
  },
  "allowed_states": [0, -1, -2],
  "partition_coefficient": -2.5  # 溶剂化能修正 (kJ/mol)
}
```

**参数获取：**
- **electron_affinities**: DFT 计算（B3LYP/def2-TZVP）
- **partition_coefficient**: SCRF 溶剂化（PCM/CPCM）
- **allowed_states**: 能化学稳定存在的电荷态

## 🐛 故障排除

### 问题 1: "Platform CUDA not available"
```
解决：改用 OpenCL
python run_voltage_sweep.py ... --platform OpenCL
```

### 问题 2: 内存不足
```
减少步数：
python run_voltage_sweep.py --equil-steps 2000 --sample-steps 3000
```

### 问题 3: 模拟与实验相差大（>200 mV）
```
调整溶剂化能：
修改 config_redox.json 中的 partition_coefficient
典型范围：-10 到 +5 kJ/mol
```

### 问题 4: MC 接受率过低（<5%）
```
增加采样时间或缩小电位步长：
python run_voltage_sweep.py --sample-steps 20000 --v-step 0.05
```

## 📚 物理模型速记

### 能量平衡
$$E(state) = z \cdot F \cdot V + G_{solv}$$

### 接受准则
$$P = \min(1, e^{-\Delta E / k_B T})$$

### 还原电位
$$E°_{red} = -\frac{EA_{gas} + \Delta G_{solv}}{F}$$

## 🎯 典型用例工作流

### 用例 1：计算单个分子的 redox potential
```bash
# 1. 在 config_redox.json 中设置分子参数
# 2. 运行完整扫描
python run_voltage_sweep.py --v-start -5 --v-end -1 --v-step 0.1

# 3. 提取半波电位
python << 'EOF'
import pandas as pd
from voltage_sweep import RedoxPotentialAnalyzer
df = pd.read_csv('sweep_results/voltage_sweep_results.csv')
e_half = RedoxPotentialAnalyzer.find_half_wave_potential(df, 'your_mol')
print(f"E_1/2 = {RedoxPotentialAnalyzer.convert_to_li_reference(e_half):.3f} V vs Li/Li+")
EOF
```

### 用例 2：与 LSV 实验对齐
```bash
# 已知实验 LSV 峰值为 -1.7 V vs Li/Li+
python run_voltage_sweep.py ... --lsv-peak-v -1.7

# 输出显示模拟与实验的偏差
```

### 用例 3：参数优化扫描
```bash
# 测试多个 partition_coefficient 值
for pc in -5 -3 -1 0 1; do
  echo "partition_coefficient: $pc" >> config_redox.json
  python run_voltage_sweep.py --output-dir sweep_pc_$pc
done
```

## 📞 需要帮助？

1. **查看工作流文档**：`cat REDOX_WORKFLOW.md`
2. **完整报告**：`cat TASK_3_4_5_COMPLETION.md`
3. **API 文档**：`python -c "from redox_mc import RedoxMC; help(RedoxMC)"`
4. **单元测试**：`python test_redox_framework.py`

## 📦 环境检查

```bash
# 验证 OpenMM 8.4
python -c "import openmm; print(f'OpenMM {openmm.__version__}')"

# 检查可用平台
python -c "import openmm as mm; print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])"

# 验证依赖包
python -c "import numpy, pandas; print('✓ All dependencies ready')"
```

---

**版本**: v1.0
**最后更新**: 2026-03-11
**框架状态**: ✅ 生产就绪
