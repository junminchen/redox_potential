# PRD：基于 OpenMM 8.4 的 Constant Potential Force + Redox MC 模拟升级

**文档版本：** v1.0  
**日期：** 2026-03-09  
**目标仓库：** `YiJung/redox_ConstantV`（或其fork）  
**目标读者：** 代码Agent / 开发工程师  
**优先级：** P0（核心功能迁移）

---

## 1. 背景与目标

### 1.1 现状问题

当前代码库（`run_openMM_1context.py` + `subroutines_1context.py`）使用的是旧版OpenMM API（`simtk.openmm`命名空间，对应OpenMM ≤7.x）。在OpenMM 8.0+中，`simtk`命名空间已被移除，所有接口迁移至`openmm`命名空间。具体破坏性变更包括：

| 旧API（≤7.x） | 新API（≥8.x） |
|---|---|
| `from simtk.openmm.app import *` | `from openmm.app import *` |
| `from simtk.openmm import *` | `from openmm import *` |
| `from simtk.unit import *` | `from openmm.unit import *` |
| `context.getState(...)` 部分参数 | 接口兼容，但行为有细微变化 |
| `NonbondedForce.updateParametersInContext` | API相同，但性能模型变化 |

此外，OpenMM 8.x 引入了新的 `CustomCVForce`、改进的 `CustomNonbondedForce`，以及对 `Platform` API 的重构，为实现**原生 ConstantPotentialForce 自定义插件**提供了更好的基础。

### 1.2 目标

1. **（P0）API迁移**：将现有代码从 `simtk.openmm` 完整迁移至 OpenMM 8.4 兼容的 `openmm` 命名空间，确保所有功能正常运行。
2. **（P0）实现 `ConstantPotentialForce` 类**：将目前分散在 `subroutines_1context.py` 中的恒电势电荷求解逻辑（`ConvergedCharge`、`Charge_solver`、`Scale_charge`等）封装为一个符合OpenMM 8.4接口规范的Force-like类，支持`CustomNonbondedForce`或`CustomExternalForce`驱动。
3. **（P1）Redox MC 模块化**：将 `MonteCarlo_redox` 逻辑解耦，支持任意redox分子体系的参数化配置，可对齐 LSV 实验 redox potential。
4. **（P1）电位扫描接口**：新增电位扫描（Voltage Sweep）功能，通过系列模拟计算各电极电位下的redox态占据概率，输出可直接与LSV实验数据对比的模拟"CV曲线"。
5. **（P2）性能优化**：利用OpenMM 8.4的CUDA/OpenCL性能提升，以及新的`VariableLangevinIntegrator`支持。

---

## 2. 技术规格

### 2.1 环境要求

```
Python >= 3.9
OpenMM == 8.4.*
numpy >= 1.24
```

PBS脚本（`2V.pbs`）中的 `cuda/9.0` 和 `python-3.6` 环境需更新为：
```bash
module load cuda/12.x
conda activate openmm84_env
```

---

### 2.2 任务一：API命名空间迁移（P0）

**文件：** `run_openMM_1context.py`、`subroutines_1context.py`

**变更要求：**

#### 2.2.1 Import语句替换

将所有文件头部的import替换如下：

```python
# 旧（需删除）
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

# 新（需添加）
from openmm.app import *
from openmm import *
from openmm.unit import *
```

#### 2.2.2 单位系统检查

OpenMM 8.x中`unit`模块路径改变，但`Quantity`语义不变。需全局搜索并验证以下用法仍正确：
- `nanometer`, `kilojoule_per_mole`, `kelvin`, `femtosecond`, `picosecond`
- `elementary_charge`（如有使用）

#### 2.2.3 Platform API更新

旧代码中若有 `Platform.getPlatformByName('OpenCL')` 或 `Platform.getPlatformByName('CUDA')`，需验证在OpenMM 8.4下仍可用（API未变，但属性命名可能变化）。

**验证标准：**
```bash
python -c "import openmm; print(openmm.__version__)"
python run_openMM_1context.py npt_mc_3pctbq.pdb md_2V.key
# 应正常启动，不报 ImportError 或 ModuleNotFoundError
```

---

### 2.3 任务二：ConstantPotentialForce 类实现（P0）

**新文件：** `constant_potential_force.py`

#### 2.3.1 类设计规格

```python
class ConstantPotentialForce:
    """
    Implements the constant potential electrode charge solver for OpenMM 8.4.
    
    The electrode charges are solved iteratively such that the electrostatic
    potential at each electrode matches the specified target voltage, using
    a method analogous to the Siepmann-Sprik / Reed-Madden approach:
    
        Φ_cathode = E_fermi - V/2
        Φ_anode   = E_fermi + V/2
    
    Charges are updated via an analytical induced charge correction:
        ΔQ_induced = ε₀ * A * ΔΦ / d_gap  (parallel plate approximation)
    
    with an iterative correction loop until |ΔQ| < tol.
    """

    def __init__(
        self,
        system: openmm.System,
        topology: openmm.app.Topology,
        cathode_atom_indices: list[int],
        anode_atom_indices: list[int],
        dummy_atom_indices: list[int],
        neutral_atom_indices: list[int],
        voltage: float,          # in Volts
        e_fermi: float,          # in eV (absolute vacuum reference)
        area_per_atom: float,    # nm^2 per graphene atom
        tol: float = 1e-4,       # convergence tolerance for charge (e)
        max_iter: int = 200,     # max charge solver iterations
    ):
        ...

    def update_charges(
        self,
        context: openmm.Context,
        positions,
        nonbonded_force: openmm.NonbondedForce,
        z_cathode: float,
        z_anode: float,
        cell_length_z: float,
    ) -> tuple[float, float]:
        """
        Solve for converged electrode charges at current positions.
        Returns (Q_cathode, Q_anode) in elementary charges.
        
        Algorithm:
        1. Compute electrostatic potential at cathode/anode planes from 
           current solvent charges using OpenMM's state forces.
        2. Calculate required induced charge via:
               q_i = q_i_prev + dq_analytical(Φ_target - Φ_current)
        3. Update NonbondedForce parameters and call 
               nonbonded_force.updateParametersInContext(context)
        4. Repeat until convergence.
        
        Note: OpenMM 8.4 supports updateParametersInContext() for
        NonbondedForce without full re-initialization.
        """
        ...

    def get_electrode_potential(
        self,
        is_cathode: bool,
    ) -> float:
        """Returns electrode electrochemical potential in kJ/mol."""
        ...
```

#### 2.3.2 与OpenMM 8.4集成方式

在 `subroutines_1context.py` 的 `MDsimulation` 类中：
- 将现有的 `ConvergedCharge()` 方法改为调用 `ConstantPotentialForce.update_charges()`
- 保留原方法作为backward-compatible wrapper（可选）

**关键OpenMM 8.4 API调用：**
```python
# 获取当前力（用于计算静电势）
state = context.getState(getForces=True, getPositions=True, getEnergy=True)

# 更新NonbondedForce中的粒子电荷（OpenMM 8.4仍支持此接口）
nonbonded_force.setParticleParameters(index, charge, sigma, epsilon)
nonbonded_force.updateParametersInContext(context)  # 不需要reinitialize

# 注意：OpenMM 8.4中reinitialize()调用成本较高，应避免在内循环中调用
# 仅在拓扑变化时才需要 context.reinitialize(preserveState=True)
```

---

### 2.4 任务三：Redox MC模块化（P1）

**文件：** `redox_mc.py`（新文件）

#### 2.4.1 RedoxMoleculeState 类

```python
@dataclass
class RedoxParameters:
    """Parameters for a redox-active molecule."""
    residue_name: str           # e.g., "qbzn"
    charge_states: dict         # {state_int: xml_filepath}
    # e.g., {0: "charge_qbzn.xml", -1: "charge_qsem.xml", -2: "charge_qhyd.xml"}
    electron_affinity: dict     # {(state_from, state_to): EA_in_eV}
    # e.g., {(0, -1): 1.53, (-1, -2): -3.5}  (EA: energy of neutralizing one electron)
    allowed_states: list[int]   # e.g., [0, -1, -2]
```

#### 2.4.2 RedoxMC 类

```python
class RedoxMC:
    """
    Metropolis Monte Carlo for electron transfer between electrode and
    redox-active molecules dissolved in electrolyte, at constant potential.
    
    Acceptance criterion (grand canonical for electrons):
    
        w = EA_gas(state_from → state_to) 
            ± E_electrode(in kJ/mol)       # + for reduction, - for oxidation
            + ΔPE_solv                      # from NonbondedForce energy difference
            + ΔE_intra                      # intramolecular contribution (default 0)
    
        P_accept = min(1, exp(-w / kT))    # kT in kJ/mol (= RT, T in Kelvin)
    
    This criterion ensures that at equilibrium, the fraction of molecules
    in state -1 vs 0 follows:
        ln(f_{-1}/f_0) = (EA_gas + ΔG_solv) / kT - e*E_electrode / kT
    which defines the thermodynamic reduction potential E°.
    """
    
    def __init__(
        self,
        redox_params: RedoxParameters,
        constant_potential_force: ConstantPotentialForce,
        temperature: float,  # Kelvin
    ):
        ...

    def attempt_electron_transfer(
        self,
        molecule_idx: int,
        context: openmm.Context,
        nonbonded_force: openmm.NonbondedForce,
        # ... other args
    ) -> bool:
        """
        Attempt one MC electron transfer move.
        Returns True if accepted.
        """
        ...

    def compute_redox_potential(
        self,
        n_samples: int = 1000,
    ) -> float:
        """
        Estimate the thermodynamic redox potential (V vs vacuum) from
        MC acceptance statistics.
        Method: find E_electrode at which acceptance rate for 
        reduction == acceptance rate for oxidation (detailed balance).
        """
        ...
```

---

### 2.5 任务四：电位扫描接口，对齐LSV实验数据（P1）

**新文件：** `voltage_sweep.py`

#### 2.5.1 功能描述

实现对 `E_electrode` 的扫描以生成模拟的循环伏安（CV）曲线，直接与LSV实验redox potential对齐：

```python
def voltage_sweep(
    pdb_file: str,
    input_file: str,
    v_start: float,       # 起始电压（V vs vacuum），如 -2.0
    v_end: float,         # 终止电压（V vs vacuum），如  0.0
    v_step: float,        # 步长（V），如 0.1
    n_equil_steps: int,   # 每个电位点的平衡步数
    n_sample_steps: int,  # 每个电位点的采样步数
    output_file: str = "simulated_cv.dat",
) -> pd.DataFrame:
    """
    Perform a potential sweep simulation analogous to LSV/CV experiment.
    
    At each electrode potential E_i:
    1. Run n_equil_steps MD steps at constant potential E_i
    2. Record mean redox state occupancy <n_redox>
    3. Compute "current" proxy: I_sim(E) ∝ d<n_redox>/dE
    
    Returns DataFrame with columns:
    - E_electrode_V_vs_vacuum: electrode potential (V)
    - E_electrode_V_vs_Li: electrode potential (V vs Li/Li+, using +1.4V offset for organic electrolyte)
    - mean_redox_state: mean charge state of redox molecules
    - fraction_reduced: fraction of molecules in most-reduced state
    - fraction_neutral: fraction of molecules in neutral state
    - simulated_current_proxy: d(fraction_reduced)/dE (normalized)
    - n_mc_accept: number of accepted MC moves at this potential
    - n_mc_total: total MC trial moves at this potential
    """
    ...
```

#### 2.5.2 与实验LSV对齐的校准

提供一个校准脚本 `calibrate_redox_potential.py`：

```python
def calibrate_to_lsv(
    sweep_results: pd.DataFrame,
    lsv_peak_potential_V_vs_Li: float,  # 实验LSV峰值电位
    reference_convention: str = "Li/Li+",  # 或 "SHE", "vacuum"
) -> dict:
    """
    Compare simulated redox potential to experimental LSV peak.
    
    The simulated E°_red is defined as E_electrode where 
    fraction_reduced = fraction_neutral = 0.5 (half-wave potential).
    
    Returns:
    - simulated_E_half_wave: simulated half-wave potential
    - experimental_E_peak: input LSV peak
    - offset_mV: systematic offset (useful for force field calibration)
    - alignment_quality: R² of simulated vs experimental CV shape
    """
    ...
```

---

### 2.6 任务五：输入文件格式更新（P0）

**文件：** `md_2V.key`（更新格式文档）

新增以下参数支持：

```ini
# ===== NEW in OpenMM 8.4 compatible version =====
openmm_version = 8.4          # 版本标记，用于API路由
integrator = LangevinMiddle   # 新增：支持 LangevinMiddle (更稳定) 或 Verlet
timestep(fs) = 2.0            # 新增：显式时间步长（旧版硬编码）
pressure(bar) = 1.0           # 新增：NPT气压（旧版硬编码）

# ===== Voltage Sweep (new) =====
sweep_mode = False            # 设为True启动电位扫描模式
sweep_v_start = -2.0          # V vs vacuum
sweep_v_end = 0.0
sweep_v_step = 0.1

# ===== Redox potential calibration =====
lsv_reference = Li/Li+        # 实验参比电极
vacuum_to_reference_offset = 1.4  # eV，真空→Li/Li+转换
```

---

### 2.7 任务六：PBS/SLURM脚本更新（P1）

**文件：** `2V.pbs` → 更新为现代集群环境

```bash
#!/bin/bash
#SBATCH -J redox_constV
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH -t 120:00:00

module load cuda/12.2
conda activate openmm84_env

export OPENMM_CUDA_COMPILER=$(which nvcc)
python run_openMM_1context.py npt_mc_3pctbq.pdb md_2V.key > energy_2V_3pctSq_300K.log
```

---

## 3. 文件变更清单

| 文件 | 操作 | 优先级 |
|------|------|--------|
| `run_openMM_1context.py` | 修改：更新import，适配新类接口 | P0 |
| `subroutines_1context.py` | 修改：更新import，将`ConvergedCharge`委托给`ConstantPotentialForce` | P0 |
| `constant_potential_force.py` | 新建：`ConstantPotentialForce`类 | P0 |
| `redox_mc.py` | 新建：`RedoxMC`类，`RedoxParameters` dataclass | P1 |
| `voltage_sweep.py` | 新建：电位扫描接口 | P1 |
| `calibrate_redox_potential.py` | 新建：实验对齐工具 | P1 |
| `md_2V.key` | 修改：新增参数字段（向后兼容） | P0 |
| `2V.pbs` | 修改：更新为SLURM+现代CUDA | P1 |
| `requirements.txt` | 新建：`openmm>=8.4, numpy>=1.24, pandas` | P0 |
| `tests/test_api_migration.py` | 新建：验证import和基础功能 | P0 |

---

## 4. 验收标准（Acceptance Criteria）

### AC-1：API迁移（P0）
```
GIVEN: OpenMM 8.4已安装
WHEN: python run_openMM_1context.py npt_mc_3pctbq.pdb md_2V.key
THEN: 
  - 无 ImportError / DeprecationError
  - 成功完成至少1个 charge_update 循环（200步MD + 1次电荷收敛）
  - 能量输出数值与旧版在±5%以内（相同随机种子）
```

### AC-2：ConstantPotentialForce（P0）
```
GIVEN: 简化测试体系（20个石墨烯原子 + 10个离子对）
WHEN: ConstantPotentialForce.update_charges() 被调用
THEN:
  - 阴极/阳极总电荷满足：|Q_cathode + Q_anode| < 0.01e（电荷中性）
  - 电荷收敛在 max_iter 内完成
  - updateParametersInContext 不触发 context.reinitialize()
```

### AC-3：Redox MC对齐LSV（P1）
```
GIVEN: voltage_sweep 在 -4.0 到 -3.0 V vs vacuum 范围扫描（对应苯醌体系）
WHEN: 计算 fraction_reduced(E)
THEN:
  - 半波电位 E_half（fraction=0.5处）落在 E_fermi + EA_neu ± 0.3V 范围内
  - 即：-4.6 + 1.53 = -3.07 V vs vacuum（≈ -1.67 V vs Li/Li+）
  - 与文献苯醌在有机电解质中的LSV峰值（~-1.7 V vs Li/Li+）一致
```

### AC-4：性能（P1）
```
GIVEN: 原始体系 npt_mc_3pctbq.pdb (约7000原子)
WHEN: 在同等GPU上运行100个charge_update周期
THEN: 总运行时间不超过旧版 OpenMM 7.x 的 120%
```

---

## 5. 关键技术注意事项

### 5.1 OpenMM 8.4 破坏性变更

1. **`simtk` 命名空间完全移除**（8.0起）：所有import必须改为`openmm`。
2. **`Integrator` 类变化**：`LangevinIntegrator` 仍可用，但推荐使用 `LangevinMiddleIntegrator` 以提高数值稳定性。
3. **PME精度提升**：OpenMM 8.x的Ewald误差控制更严格，可能导致初始能量略有差异。
4. **`NonbondedForce.updateParametersInContext()` 的线程安全性**：在多Context场景需加锁（当前代码是单Context，不受影响）。
5. **`CustomNonbondedForce` 新特性**：8.4支持 tabulated functions，可用于更高精度的电极-溶液静电势计算。

### 5.2 恒电势模拟的物理约束

- 电荷求解时必须满足**总电荷守恒**：$Q_\text{cathode} = -Q_\text{anode}$（对称体系）
- `area_per_atom` 的计算：`sheet_area / (N_graphene / 2)`，其中除以2是因为每层有两个sublattice点（现有代码已正确实现）
- `conv = 18.8973 / 2625.5`（bohr/nm × au per kJ/mol）单位转换因子需保留

### 5.3 LSV对齐的理论基础

热力学redox电位满足Nernst方程的零温极限：

$$E^\circ_{\text{red}} = \frac{1}{e}\left[EA_{\text{gas}} + \Delta G_{\text{solv}}(\text{red}) - \Delta G_{\text{solv}}(\text{ox})\right]$$

代码中 `EA_neu`、`EA_red1` 覆盖了气相项，而 `ΔPE_MD`（NonbondedForce能量差）覆盖了溶剂化自由能差（在MD时间平均意义上）。MC接受率统计即提供了对此自由能差的有效采样。

---

## 6. 不在范围内（Out of Scope）

- 量子化学层面的EA重算（EA_neu/EA_red1由外部DFT给出）
- 多步电子转移动力学（Marcus速率理论）
- 界面SEI膜的明确建模
- 锂离子溶剂化/嵌入过程（此版本仅针对液相redox分子）

---

## 7. 参考资料

- Reed, S. K. et al. "A new algorithm for unrestricted constant potential molecular dynamics" *J. Chem. Phys.* 126, 084704 (2007)
- Siepmann, J.I. & Sprik, M. "Influence of surface topology and electrostatic potential on water/electrode systems" *J. Chem. Phys.* 102, 511 (1995)  
- Wang, Z. et al. "Modelling electrochemical systems with finite field molecular dynamics" *Faraday Discuss.* 210, 1 (2018)
- OpenMM 8.4 Migration Guide: https://openmm.org/documentation/latest/migration/

