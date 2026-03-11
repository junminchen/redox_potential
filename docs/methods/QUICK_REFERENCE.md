# 快速参考卡片

正式说明看：

- [USER_MANUAL_CN.md](/home/am3-peichenzhong-group/Documents/project/project_solv_structure/redox_potential/redox_constantV/docs/methods/USER_MANUAL_CN.md)

## 30 秒开始

```bash
source activate omm84_cuda_pip

python tests/test_redox_framework.py

python scripts/cli/batch_screen_formulations.py \
  --input configs/formulations/formulation_batch_example_extended.json \
  --output-dir results/formulation_batch_quickstart
```

## 最常用命令

| 命令 | 用途 |
|------|------|
| `python tests/test_redox_framework.py` | 验证框架完整性 |
| `python scripts/cli/run_voltage_sweep.py --help` | 查看 MD sweep 选项 |
| `python scripts/cli/batch_screen_formulations.py --help` | 查看配方批处理选项 |
| `python scripts/dft/generate_redox_config_from_dft.py --help` | 查看 DFT 接口选项 |

## 常用工作流

### 1. MD 电位扫描

```bash
python scripts/cli/run_voltage_sweep.py \
  --config configs/md/config_opls.json \
  --redox-config configs/redox/config_redox.json \
  --pdb structures/systems/start_with_electrodes.pdb \
  --output-dir results/voltage_sweep_quickstart
```

### 2. 配方批处理筛选

```bash
python scripts/cli/batch_screen_formulations.py \
  --input configs/formulations/formulation_batch_example_extended.json \
  --output-dir results/formulation_batch_quickstart
```

### 3. DFT 生成 Redox 配置

```bash
python scripts/dft/generate_redox_config_from_dft.py \
  --workflow configs/dft_workflows/example_ec_dmc_reduction.json \
  --output-dir results/dft_interface_quickstart \
  --allow-missing-pyscf
```

## 关键输入位置

- MD 配置：`configs/md/config_opls.json`
- 还原侧配置：`configs/redox/`
- 氧化侧配置：`configs/oxidation/`
- 反应通道配置：`configs/pathways/`
- 配方批处理输入：`configs/formulations/`
- DFT workflow：`configs/dft_workflows/`
- 体系结构：`structures/systems/`
- 单分子模板：`structures/pdb_bank/`
- 力场文件：`ff/`

## 常见输出位置

- MD sweep：`results/voltage_sweep_*`
- 配方批处理：`results/formulation_batch_*`
- DFT 接口：`results/dft_interface_*`
- 历史调试运行：`results/debug_runs/`

## 电压定义

- `voltage_v` 表示总电池电压，参考为 vacuum。
- 恒电势边界条件实际使用 `cathode = +V/2`、`anode = -V/2`。
- 近似转换：`V vs Li/Li+ = V vs vacuum + 1.4`。

## 常见排错

### CUDA 不可用

```bash
python scripts/cli/run_voltage_sweep.py ... --platform OpenCL
```

### 想先检查框架有没有坏

```bash
python tests/test_redox_framework.py
bash run_examples.sh test
```

### 想看当前脚本都能不能跑

```bash
bash run_examples.sh all
```

## 进一步阅读

- 详细手册：`docs/methods/USER_MANUAL_CN.md`
- Redox 工作流：`docs/methods/REDOX_WORKFLOW.md`
- PySCF 接口：`docs/methods/PYSCF_DFT_INTERFACE.md`
- 参数来源：`docs/provenance/REDOX_PARAMETER_PROVENANCE.md`
