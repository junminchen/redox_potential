# 实验验证与模型推论说明

## 目的

这份说明专门回答两个问题：

1. 当前项目里哪些结论有文献或实验语境支持。
2. 当前项目里哪些数值只是为了让 screening 模型可用而做的有效映射或半经验推论。

建议把这份文件和下面两个文件一起看：

- [REDOX_PARAMETER_PROVENANCE.md](/home/am3-peichenzhong-group/Documents/project/project_solv_structure/redox_potential/redox_constantV/docs/provenance/REDOX_PARAMETER_PROVENANCE.md)
- [redox_parameter_provenance.csv](/home/am3-peichenzhong-group/Documents/project/project_solv_structure/redox_potential/redox_constantV/docs/provenance/redox_parameter_provenance.csv)

## 证据分级

本项目现在把证据分成三层。

### A. 文献直接数值锚点

- 文献中给出了可直接使用的数值。
- 例如某个溶剂在特定溶剂化或 Li+ 配位条件下的还原电位目标值。
- 这类数据是当前模型里最硬的数值来源。

### B. 文献趋势或实验语境锚点

- 文献没有给出一个可以直接抄到模型里的唯一常数。
- 但给出了足够清楚的实验趋势、排序、稳定性窗口或失效方向。
- 这类信息可以支持“谁更早还原”“谁在高压下更先失效”这类判断。

### C. 半经验映射或模型推论

- 这是为了让 screening 工作流能输出 onset、排序和风险边界而引入的有效参数。
- 它们通常不是严格第一性原理结果，也不是某篇实验文献里直接测出来的唯一值。
- 典型例子是：
  - `state_free_energy_offsets_kjmol`
  - 高压分解通道的 `barrier_at_onset`
  - 势垒随电位下降的 `slope`

## 最重要的阅读原则

- 如果一个结论属于 A 层，可以把它当作当前模型的主要标定依据。
- 如果一个结论属于 B 层，可以把它当作“实验一致的方向性约束”。
- 如果一个结论属于 C 层，只能把它当作 screening 假设，不能当作最终真值。

## 当前 EC/DMC 基线模型的来源拆分

| 主题 | 当前项目里的结论 | 文献来源 | 证据类型 | 我们在模型里怎么用 | 当前不能宣称什么 |
|---|---|---|---|---|---|
| EC 还原比 DMC 更早 | `EC E1/2 ~ 0.47 V vs Li/Li+`，`DMC E1/2 ~ 0.38 V vs Li/Li+` | PMC9052788: https://pmc.ncbi.nlm.nih.gov/articles/PMC9052788/ | A. 文献直接数值锚点 | 映射成 `configs/redox/config_redox_electrolyte_effective.json` 里的 `state_free_energy_offsets_kjmol` | 不能说这就是任意 EC/DMC 配方、任意电极表面的通用实验峰位 |
| EC 在高压下更容易先出问题 | `EC` 被放在比 `DMC` 更低的高压风险区 | Energy & Fuels review: https://pubs.acs.org/doi/10.1021/acs.energyfuels.4c02728 | B. 文献趋势锚点 | 用来约束氧化 onset 和分解通道排序：`EC` 先于 `DMC` | 不能说 `4.4 V`、`4.2 V` 这些值是文献里唯一直接测得的 EC 内禀常数 |
| 碳酸酯体系在 3.0-6.0 V 区间存在氧化响应语境 | 高压下碳酸酯基线电解液会出现明显氧化/失效问题 | PMC8695088: https://pmc.ncbi.nlm.nih.gov/articles/PMC8695088/ | B. 文献趋势锚点 | 用作 `LSV/CV-like` 图和高压筛选语境的实验背景 | 不能说我们现在的曲线已经等价于实验 LSV 电流曲线 |
| EC 环开裂/氧化诱导失效比 DMC 更早 | `EC radical-cation ring opening` 是当前模型里的优先失效通道 | Energy & Fuels review: https://pubs.acs.org/doi/10.1021/acs.energyfuels.4c02728 | B + C | 通道名字和排序来自文献趋势；通道 onset 与 barrier 是 screening 参数 | 不能说当前 barrier 已经等同于 DFT/NEB 势垒 |

## 为什么当前会出现 `4.343 V`、`4.200 V`、`3.500 V` 这类数

这些数需要分开理解。

### `4.343 V vs Li/Li+`

- 含义：当前“可逆氧化层”给出的较早氧化 onset。
- 来源：不是文献里直接给出的单一实验峰位。
- 它是基于“EC 在高压下更早成为短板”的文献趋势，映射到当前离散占有率模型后的结果。

### `4.200 V vs Li/Li+`

- 含义：当前“不可逆分解风险层”给出的 EC 分解起始边界。
- 来源：不是直接实验常数。
- 它是用来表达“EC 可在高压界面更早进入失效化学”的有效 threshold。

### `3.500 V vs Li/Li+`

- 含义：在加入速率阈值以后，当前反应通道模型给出的更严格 failure threshold。
- 来源：这是 C 层结果。
- 它强依赖：
  - 通道起始点
  - 势垒起点
  - 势垒随电位变化的斜率
  - 选用的速率阈值

所以 `3.500 V` 只能解释成：

- “按当前 screening 原型，如果采用这组通道参数和速率阈值，EC/DMC 基线体系会在 5 V 前进入明显失效区”

不能解释成：

- “实验上这套体系一定在 3.500 V 就绝对分解”

## 配方示例部分哪些有实验背景

下面这些示例配方都不是凭空杜撰，而是有文献背景的 screening 预设。

| 示例配方 | 文献锚点 | 文献可支持的结论 | 当前模型里怎么用 | 仍然属于推论的部分 |
|---|---|---|---|---|
| `4 m LiFSI in DME` | Scientific Reports: https://www.nature.com/articles/s41598-017-16268-7 | 高浓 LiFSI/DME 对锂金属侧有利，并能支持更高电压窗口 | 用作锂金属友好型示例配方 | 当前模型里的精确高压上限仍是有效参数结果 |
| `LiFSI:EC 1:4` | PMC8790720: https://pmc.ncbi.nlm.nih.gov/articles/PMC8790720/ | 浓盐 EC 体系可提升高压性能，相比普通碳酸酯更耐高压 | 用作高浓 EC 路线示例 | 当前 `4.893 V` 这种数值不是文献直接测得的普适常数 |
| `1:1.1 LiFSI in DMC` | PMC4931331: https://pmc.ncbi.nlm.nih.gov/articles/PMC4931331/ | 超浓 LiFSA/LiFSI-DMC 路线可支持 5.2 V 级高压正极 | 用作 5 V 候选示例 | 当前模型里 `5.093 V` 上限仍是 screening 拟合结果 |
| `1 M LiPF6 in EC:EMC (3:7)` | PMC11428382: https://pmc.ncbi.nlm.nih.gov/articles/PMC11428382/ | 作为常规碳酸酯基线是合理的实验语境 | 用作常规基线对照 | 当前模型对其高压上限的数值仍是半经验映射 |

## 给实验人员时应该怎么表述

推荐这样说：

- `EC` 比 `DMC` 更容易先还原，这个排序有文献和当前模型双重支持。
- 常规 `EC/DMC` 碳酸酯体系不适合作为 `5 V` 高压电解液，这和文献整体趋势一致。
- 当前模型给出的具体窗口边界用于筛选和排序，不应直接当作最终实验真值。
- 若要把某个具体电位边界作为结论提交，需要进一步做：
  - DFT/cluster 标定
  - 显式界面修正
  - 反应路径势垒计算
  - 与同体系 LSV/CV onset 对齐

不推荐这样说：

- “模型已经直接复现了实验 LSV 曲线。”
- “`3.500 V` 就是该电解液不可争议的真实分解电位。”
- “`4.343 V` 是 EC 的普适内禀氧化电位。”

## 当前最稳妥的结论模板

如果需要写在报告或汇报里，建议用下面这段话：

> 当前工作流已经把文献中的实验趋势和部分数值锚点映射到统一的 screening 框架中。  
> 其中，EC/DMC 的还原排序与高压风险顺序有明确文献支持；而具体的高压失效阈值和反应通道速率仍属于半经验模型结果，适合用于前筛、排序和实验窗口缩小，不宜直接替代最终电化学测量。

## 下一步怎样把“推论”升级成“更可信数值”

1. 用 `PySCF/DFT` 或 cluster 计算替换当前 `state_free_energy_offsets_kjmol` 的经验映射。
2. 用显式反应路径或 NEB 结果替换当前 `barrier_at_onset` 和 `slope`。
3. 用同一电解液体系的实验 `LSV/CV onset` 做回标定。
4. 把不同配方的文献实验窗口汇总成一张数据库表，统一约束 screening preset。
