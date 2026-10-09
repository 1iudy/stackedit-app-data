# icet 中文全模块教程

从团簇展开建模到构型采样与相图分析

阅读日期：2026 年 10 月 9 日。依据：icet 官方文档及其链接的官方 CE 扩展教程。适用读者：具备基础 Python、晶体结构和材料计算知识，希望建立可复现团簇展开研究流程的学生与研究人员。

本教程覆盖此次读取到的官方用户文档、进阶专题、公开 API 参考、FAQ、术语与引用信息。内容按学习顺序重新组织，用中文解释功能、输入、计算流程和判断方法；附录列出官方页面覆盖关系及公开接口。这里的“全模块”以公开文档为边界，并不意味着逐行解释全部 C++ 内部实现，也不把未发布功能算作已支持功能。

示例分为三类：**独立教学示例**使用人工设定的 CE 参数，不代表真实材料；**官方数据示例**需要下载官方数据；**上下文片段**接续前文的 `cs`、`sc`、`ce`、`A`、`y` 等变量。参数如截断距离、温度、步数和正则化强度均须针对研究体系重新验证。本文代码经过语法与文档接口核对；本次编写未在本机安装 icet，也未重新执行完整 DFT、拟合和相图计算。官方页面中的计时、误差和相变温度不作为本次运行结果。

Markdown 版适合编辑和复用代码；HTML 阅读版包含目录、搜索、代码复制和内嵌公式渲染资源，可离线阅读。

## 阅读路线

|目标|优先章节|应得到的结果|
|---|---|---|
|第一次使用|1 至 7，10，16 至 19|理解输入结构、训练 CE 并运行 MC|
|生成 DFT 候选结构|8，9，12，15|枚举、随机结构、SQS 或主动学习候选集|
|改善拟合|10 至 13|通过独立验证选择模型和约束|
|研究相图与有序化|14，17 至 25|自由能、序参量、相界及误差分析|
|研究表面、空位或应变|5，6，13，20，26|明确子晶格、归一化与模型适用范围|
|开发自定义采样|27，28，附录 A|正确维护计算器状态与并行任务|
|定位错误|30，附录 A 和 B|根据输入和接口找到具体原因|

## 1 软件定位与整体工作流

icet 用来建立和采样合金团簇展开模型。Python 接口便于与 ASE、第一性原理数据库、数值分析和绘图工具衔接，主要计算密集部分由 C++ 实现。随包提供的 `mchammer` 负责蒙特卡洛采样；拟合主要使用独立软件包 `trainstation`。

团簇展开把固定母晶格上的占位构型映射到能量或其他标量性质。常见问题包括有序相稳定性、混溶间隙、构型热力学和短程有序。它不是独立的 DFT 程序：训练标签要从外部计算或实验模型获得，也不会自动产生声子、电子或磁性自由能。

核心对象之间的关系如下。

```text
理想母晶格 + 允许占位 + 截断距离
                    |
               ClusterSpace
                    |
理想占位结构 + 性质标签 --> StructureContainer
                    |
             设计矩阵 A 与标签 y
                    |
       trainstation 拟合与交叉验证
                    |
             ClusterExpansion
              /             \
   新结构预测与凸包       超胞专用 Calculator
                                  |
                          Ensemble + Observer
                                  |
                            DataContainer
                                  |
                       统计量 自由能 相界
```

每个研究项目应先明确五件事：母晶格是什么；哪些位点允许哪些物种；拟合性质按什么单位归一化；训练集覆盖哪些构型；最终需要什么热力学量。后续大多数错误可以追溯到这五项中某项定义不清。

来源：[官方首页](https://icet.materialsmodeling.org/index.html)、[Workflow](https://icet.materialsmodeling.org/get_started/workflow.html)。

## 2 安装 环境与版本管理

### 2.1 安装稳定版本

此次官方安装页要求 Python 3.10 或更新版本。建议在独立环境内安装并固定版本。以下为环境配置建议，Python 3.11 是示例选择，不是 icet 的唯一可用版本。

```bash
conda create -n icet-study python=3.11
conda activate icet-study
conda install -c conda-forge icet
```

也可以在已激活的虚拟环境里使用：

```bash
python -m pip install icet
```

官方说明 PyPI 安装需要编译 C++ 核心，要求支持 C++17 的编译器；conda 包提供预编译二进制。不同平台和版本的包供应情况仍需以安装时的发行信息为准。Windows 编译失败时，先核查工具链和依赖，不要把错误简单归结为 Python 代码。

### 2.2 开发版本与测试

```bash
python -m pip install git+https://gitlab.com/materials-modeling/icet.git
```

研究生产环境应记录提交号或发行版本。运行官方测试时，测试源代码必须与安装版本对应：

```bash
git clone https://gitlab.com/materials-modeling/icet.git
cd icet
python -m pip install -e ".[test]"
pytest tests/
```

### 2.3 检查实际环境

```python
import sys
from importlib.metadata import version
from icet import ClusterSpace, StructureContainer, ClusterExpansion
from mchammer.ensembles import CanonicalEnsemble

print(sys.executable)
for name in ("icet", "ase", "trainstation", "numpy", "spglib"):
    print(name, version(name))
```

`ModuleNotFoundError` 常见原因是终端和 Jupyter 使用不同解释器。将 `sys.executable` 与安装命令所用解释器对应起来，再检查依赖。

主要依赖涉及 ASE、NumPy、SciPy、pandas、spglib、numba 和 trainstation；C++ 核心依赖 Eigen 和 pybind11。拟合算法还通过 trainstation 使用相关机器学习库。需要自行安装绘图和 Notebook 工具时，应与核心环境一起记录。

```bash
python -m pip freeze > requirements-frozen.txt
conda env export > environment.yml
```

来源：[Installation](https://icet.materialsmodeling.org/get_started/installation.html)、[FAQ](https://icet.materialsmodeling.org/backmatter/faq.html)。

## 3 团簇展开的数学基础

### 3.1 从构型到线性模型

设构型为 $\boldsymbol{\sigma}$，每个分量表示一个母晶格位点的占位。对某个强度性质 $q$，采用下列表示：

$$
q(\boldsymbol{\sigma}) = \sum_{\alpha} m_{\alpha} J_{\alpha}\,\Gamma_{\alpha}(\boldsymbol{\sigma})
=\boldsymbol{\Gamma}(\boldsymbol{\sigma})\cdot\boldsymbol{p}.
$$

$\Gamma_{\alpha}$ 是对对称等价团簇及对应点函数组合进行平均后的相关函数；$m_{\alpha}$ 是该项的重数；$J_{\alpha}$ 是未乘重数的有效团簇相互作用；icet 的拟合参数为 $p_{\alpha}=m_{\alpha}J_{\alpha}$。因此向 `ClusterExpansion` 传入 `opt.parameters` 时不能再额外乘一次重数。

几何团簇包括空团簇、单点、对、三体、四体等。轨道是晶体对称操作下相互等价的一组几何团簇。多组分体系中，一个几何轨道可对应多个点函数组合，所以“几何轨道数”与“参数数”不一定相等。

### 3.2 零阶项与多组分点函数

空团簇相关函数恒为 1，提供模型截距。二元体系可用两种离散占位值解释单点函数，但实际物种与符号的映射应以 icet 的约定为准。多组分体系使用一组完整的离散点函数，而不是简单地把原子序数当作连续回归变量。

不要跨不同的母晶格、物种集合、轨道合并规则或软件版本，直接混用参数数组。参数的意义取决于对应 `ClusterSpace` 的完整定义和索引顺序。

以 $M$ 种允许占位、整数占位编号 $\sigma$ 为例，官方采用的离散点函数可写为：

$$
\Theta_n(\sigma)=\begin{cases}
1,&n=0,\\
-\cos[\pi(n+1)\sigma/M],&n\text{ 为奇数},\\
-\sin[\pi n\sigma/M],&n>0\text{ 为偶数}.
\end{cases}
$$

团簇函数由各位点点函数相乘得到，再按轨道平均。实际使用无需自行重写这些函数；多子晶格、多物种排列和对称平均由 `ClusterSpace` 完成。

### 3.3 设计矩阵

对 $M$ 个结构和 $P$ 个参数，有：

$$
A_{i\alpha}=\Gamma_{\alpha}(\boldsymbol{\sigma}_i),\qquad
\boldsymbol{y}=A\boldsymbol{p}+\boldsymbol{\epsilon}.
$$

增加阶数和截断距离可以提高表达能力，也会增加参数、共线性和过拟合风险。用独立验证判断误差是否真的改善，比追求更低的训练误差更重要。模型不是对任意新晶体都有效，它描述的是给定母晶格上的构型空间。

来源：[Cluster expansions](https://icet.materialsmodeling.org/get_started/cluster_expansions.html)、[ClusterSpace](https://icet.materialsmodeling.org/moduleref/cluster_space.html)、[ClusterExpansion](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)。

## 4 结构 数据与能量约定

### 4.1 ASE 结构

icet 的结构输入主要是 `ase.Atoms`。必须核对晶胞、原子位置、周期边界条件、化学符号和母晶格关系。

```python
from ase.build import bulk
from ase.io import read, write

primitive = bulk("Ag", "fcc", a=4.09)
supercell = primitive.repeat((3, 3, 3))
supercell[0].symbol = "Pd"
supercell.wrap()
write("candidate.extxyz", supercell)
restored = read("candidate.extxyz")
```

`repeat((3,3,3))` 复制的是当前输入晶胞。FCC 原胞通常只有一个位点，不能把它误认为四原子的常规立方晶胞。

### 4.2 每原子 每位点与总能量

二元混合能的一种常用定义是：

$$
e_{\mathrm{mix}} = \frac{E_{\mathrm{tot}}}{N}
-(1-c)e_A^{\mathrm{ref}}-c e_B^{\mathrm{ref}}.
$$

先计算总能量除以归一化数，再减去组成加权的参考能。两个参考项都应减去。参考能必须使用一致的计算设置和明确的参考相；替换参考态会改变化学势和截距解释。

存在空位或固定背景子晶格时，`len(atoms)` 可能包含 `X` 位点，不能自动等同于实际原子数、金属原子数或化学式单位数。为每个标签明确写出“eV/母晶格位点”“eV/金属原子”或其他基准。`ce.predict` 返回拟合标签的单位，不会自动推断物理归一化。

### 4.3 数据记录

建议每条记录同时保存：理想占位结构、弛豫结构、总能量、归一化标签、参考能版本、计算参数、母晶格编号和结构来源。弛豫后的能量可以用于训练，但团簇向量应从相应的理想占位结构提取。

保留最终弛豫构型与初始占位的对应关系。发生跨位点迁移时，初始理想占位可能已不代表最终能量，需要重新映射。

来源：[Constructing a cluster expansion](https://icet.materialsmodeling.org/get_started/construct_cluster_expansion.html)、[Mapping structures](https://icet.materialsmodeling.org/advanced_topics/mapping_structures.html)、[官方 canonical 扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)。

## 5 创建 ClusterSpace

### 5.1 基本调用与截断距离

```python
from icet import ClusterSpace
from ase.build import bulk

primitive = bulk("Ag", "fcc", a=4.09)
cs = ClusterSpace(
    structure=primitive,
    cutoffs=[8.0, 5.0, 4.0],
    chemical_symbols=["Ag", "Pd"],
)
print(cs)
print("参数数", len(cs))
```

`cutoffs` 的第一项控制二体，第二项控制三体，第三项控制四体，继续追加对应更高阶。截断衡量团簇内部任意两位点之间的**最大距离**，单位与 ASE 坐标一致，通常为 Angstrom。输出表中的 `radius` 则是各位点距几何中心的平均距离，两者不能混用。

### 5.2 子晶格与固定占位

若所有位点具有相同允许物种，给出扁平列表；若不同位点有不同允许集合，外层列表长度必须等于输入结构的位点数。

```python
from ase import Atom

primitive2 = bulk("Pd", "fcc", a=4.0)
primitive2.append(Atom("H", position=(2.0, 2.0, 2.0)))
cs2 = ClusterSpace(
    primitive2,
    cutoffs=[4.5],
    chemical_symbols=[["Au", "Pd"], ["H", "X"]],
)
print(cs2)
sublattices = cs2.get_sublattices(primitive2.repeat((2, 2, 2)))
for sub in sublattices:
    print(sub.symbol, sub.chemical_symbols, sub.indices)
```

`X` 表示空位。只有一种允许物种的位点为非活跃位点，例如 `["O"]`；它们仍是结构的一部分，但没有占位自由度。子晶格字母要从实际打印结果读取，不应在复杂体系中凭原子顺序猜测。

### 5.3 对称性 容差与结构兼容性

`symprec` 用于识别晶体对称性；`position_tolerance` 和派生的分数坐标容差用于位置匹配。首先保证结构确实是母晶格的整数超胞，再考虑容差。任意放大容差可能把物理不同的位点错误合并。

```python
candidate = primitive.repeat((2, 3, 4))
cs.assert_structure_compatibility(candidate)
print(cs.space_group)
print(cs.number_of_orbits_by_order)
print(cs.to_dataframe().head())
cs.write("cluster_space.cs")
cs_reloaded = ClusterSpace.read("cluster_space.cs")
```

### 5.4 二维 表面与真空

官方要求 icet 结构在三个方向上都设置周期边界。二维模型或表面可通过在非物理周期方向加入充分真空，使团簇截断不会穿越真空连接相邻镜像。枚举和 SQS 工具的 `pbc` 设置与 `ClusterSpace` 的建模周期约定要分别理解。

母晶格对称性是晶体学对称性；任意非晶体点群，例如某些二十面体纳米颗粒的完整对称性，不属于当前文档支持范围。

来源：[ClusterSpace 参考](https://icet.materialsmodeling.org/moduleref/cluster_space.html)、[Core components](https://icet.materialsmodeling.org/moduleref/core.html)、[FAQ](https://icet.materialsmodeling.org/backmatter/faq.html)。

## 6 弛豫结构映射与诊断

`map_structure_to_reference` 将弛豫结构映射到参考晶格的理想超胞，返回理想结构和 `StructureMapping` 诊断对象。

```python
from icet.tools import map_structure_to_reference

reference = bulk("Au", "fcc", a=4.0)
relaxed = reference.repeat((3, 3, 3))
relaxed[0].symbol = "Pd"
relaxed.rattle(stdev=0.05, seed=42)
relaxed.set_cell(1.02 * relaxed.cell, scale_atoms=True)

mapped, info = map_structure_to_reference(relaxed, reference)
print(info)
print(info.drmax, info.dravg)
print(info.transformation_matrix)
print(info.warnings)
```

重点检查最大与平均位移、整数变换矩阵、移除的整体平移、应变、旋转、含糊位点和警告。当前接口通过属性访问诊断信息，不能假定所有版本都返回旧式字典。

有空位的体系可通过 `inert_species` 指定无空位子晶格上的物种，从而辅助确定超胞尺寸。若不存在这样的子晶格，文档给出 `assume_no_cell_relaxation=True` 的路径，但这意味着尺寸识别依赖无晶胞弛豫假设。间隙位点必须预先加入参考结构；参考中没有对应位点的额外原子不能靠映射凭空生成位点。

`tol_cell` 控制晶胞偏离容限，`find_translation` 控制整体平移搜索，`trigger_levels` 控制诊断警告阈值。映射支持旋转结构，但变换搜索可能更耗时。不要只依据“函数返回成功”接纳结构；跨位点迁移、过大应变或多解映射可能表明该母晶格不适合该数据。

```python
from icet.input_output.logging_tools import set_log_config
set_log_config(level="DEBUG")
```

来源：[Mapping structures](https://icet.materialsmodeling.org/advanced_topics/mapping_structures.html)、[映射 API](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.map_structure_to_reference)。

## 7 StructureContainer 与训练数据矩阵

```python
from icet import StructureContainer

sc = StructureContainer(cs)
structures = [primitive.copy(), primitive.repeat((2, 2, 2))]
structures[1][0].symbol = "Pd"
labels = [0.0, -0.01]  # 人工教学标签 不代表 DFT 结果
for i, (structure, value) in enumerate(zip(structures, labels)):
    sc.add_structure(
        structure,
        user_tag=f"example-{i}",
        properties={"mixing_energy": value},
    )
A, y = sc.get_fit_data(key="mixing_energy")
assert A.shape == (len(structures), len(cs))
assert y.shape == (len(structures),)
print(sc.available_properties)
```

这个两结构例子只展示数据接口，不足以训练有预测能力的模型。实际拟合须建立足够丰富的数据集。

`properties` 可以存多个目标字段；`get_fit_data(key=...)` 明确选择一个标签。默认键为 `energy`，与 `mixing_energy` 不同。`structure_indices` 可指定子集；`get_structure_indices(user_tag=...)` 可定位记录。容器可保存、恢复或转换为 DataFrame。

`add_structure` 的 `sanity_check` 负责基本检查，`allow_duplicate` 控制重复策略。重复团簇向量会改变有效样本权重，即使结构文件名不同也可能具有相同特征。重复计算有时有统计意义，应明确其目的。

```python
print(sc.get_condition_number())
sc.write("training.sc")
sc_reloaded = StructureContainer.read("training.sc")
```

设计矩阵列相关、训练组成过窄或候选构型太相似时，问题可能病态。样本数大于参数数也不能保证所有参数可辨识；稀疏拟合则不要求严格按“样本数超过参数数”的经验规则判断可用性。

来源：[StructureContainer](https://icet.materialsmodeling.org/moduleref/structure_container.html)、[Generating training structures](https://icet.materialsmodeling.org/advanced_topics/training_set_generation.html)。

## 8 枚举结构 空位与吸附体系

### 8.1 非等价占位枚举

```python
from icet.tools import enumerate_structures, enumerate_supercells

for i, structure in enumerate(enumerate_structures(
    primitive,
    sizes=range(1, 5),
    chemical_symbols=["Ag", "Pd"],
)):
    write(f"enumerated-{i:04d}.extxyz", structure)
```

`sizes` 是相对输入晶胞的超胞倍数，不是必然的原子数。函数返回迭代器，适合逐个处理，避免把大型枚举一次转换成列表。

```python
dilute = enumerate_structures(
    primitive,
    sizes=range(1, 9),
    chemical_symbols=["Ag", "Pd"],
    concentration_restrictions={"Pd": (0.0, 0.125)},
)
supercell_shapes = enumerate_supercells(primitive, sizes=range(1, 5))
```

浓度区间过滤候选，不能保证任意有限尺寸实现任意分数。`enumerate_supercells` 只枚举超胞形状，`enumerate_structures` 还枚举占位。`niggli_reduce` 控制晶胞化简行为；周期边界会影响可枚举的超胞方向。

### 8.2 空位与表面吸附

用 `X` 把空位当作离散占位物种。准备 DFT 输入时，可从理想结构副本中删除 `X`，但训练和采样所用理想结构必须保留位点。固定基底、可替换表层和可空置吸附位点可以分别定义允许物种。不能让结构读取程序把 `X` 错误当成真实化学元素参与势能计算。

枚举数量随尺寸、物种数和低对称性迅速增长。对大原胞改用随机生成、条件数筛选或主动学习；枚举用来寻找小周期有序结构也很有效。

来源：[Structure enumeration 进阶专题](https://icet.materialsmodeling.org/advanced_topics/structure_enumeration.html)、[Structure enumeration API](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.enumerate_structures)、[扩展枚举教程](https://ce-tutorials.materialsmodeling.org/part-1/enumeration.html)。

## 9 随机结构 条件数筛选与正交结构

### 9.1 随机占位

```python
from icet.tools.structure_generation import occupy_structure_randomly

pool = []
for seed in range(20):
    structure = primitive.repeat((2, 2, 2))
    occupy_structure_randomly(
        structure, cs,
        target_concentrations={"Ag": 0.5, "Pd": 0.5},
        random_seed=seed,
    )
    pool.append(structure)
```

此函数修改输入结构。组成乘以对应子晶格位点数应产生整数；不能要求 8 个位点恰好实现 1/3 浓度。研究用候选池应改变组成、超胞大小和形状，而不只是改变同一组成的随机种子。

### 9.2 最小化条件数

在较大的候选池中，`structure_selection_annealing` 通过交换候选集成员，降低设计矩阵病态程度。

```python
from icet.tools.training_set_generation import structure_selection_annealing

# 上下文片段 pool 应是覆盖组成和几何的充足候选池
indices, history = structure_selection_annealing(
    cluster_space=cs,
    monte_carlo_structures=pool,
    n_structures_to_add=min(10, len(pool)),
    n_steps=2000,
)
selected = [pool[i] for i in indices]
```

上面固定 50% 组成的演示候选池不适合识别全部组成依赖参数；正式使用本片段前应替换为丰富候选池。`base_structures` 可加入已有计算结构，返回的索引针对传入的候选池。

条件数可写为 $\operatorname{cond}(A)=s_{\max}/s_{\min}$，其中 $s$ 为奇异值。若物理约束造成严格列依赖，仅筛选结构不能消除该问题，需要减少冗余变量或使用适当约束。

### 9.3 正交结构生成

扩展教程展示了另一条路线：生成候选团簇向量，对已有向量正交化，根据单点分量确定组成，再用目标向量匹配工具寻找可实现结构。它提高特征空间覆盖，但任意数学向量并不一定对应实际构型；应检查目标与获得向量之间的距离。

条件数和正交性都是数据设计工具，不直接保证相图精度。还应保留低能有序构型、端点和最终 MC 会访问的构型。

来源：[Generating training structures](https://icet.materialsmodeling.org/advanced_topics/training_set_generation.html)、[随机结构教程](https://ce-tutorials.materialsmodeling.org/part-1/random-structures.html)、[条件数筛选教程](https://ce-tutorials.materialsmodeling.org/part-1/condition-number-selection.html)、[正交结构教程](https://ce-tutorials.materialsmodeling.org/part-1/orthogonal-structures.html)。

## 10 拟合 交叉验证与模型保存

### 10.1 官方 Ag Pd 数据实战

从 [官方 tutorial.zip](https://icet.materialsmodeling.org/tutorial.zip) 下载并解压数据，定位 `reference_data.db` 后在对应目录运行。PowerShell 可使用 `Invoke-WebRequest` 和 `Expand-Archive`，也可手工下载。

```python
from ase.db import connect
from icet import ClusterSpace, StructureContainer, ClusterExpansion
from trainstation import CrossValidationEstimator

db = connect("reference_data.db")
primitive = db.get(id=1).toatoms()  # 仅该官方数据库有此约定
cs = ClusterSpace(primitive, [13.5, 6.5, 6.0], ["Ag", "Pd"])
sc = StructureContainer(cs)
for row in db.select():
    sc.add_structure(
        row.toatoms(), user_tag=row.tag,
        properties={"mixing_energy": row.mixing_energy},
    )

A, y = sc.get_fit_data(key="mixing_energy")
opt = CrossValidationEstimator((A, y), fit_method="ardr")
opt.validate()
opt.train()
ce = ClusterExpansion(cs, opt.parameters, metadata=opt.summary)
ce.write("mixing_energy.ce")
sc.write("training.sc")
print(opt)
```

`validate()` 估计验证表现，`train()` 拟合最终模型，两者用途不同。某些官方说明文字与代码的算法名称不一致，执行时应以 `fit_method` 实际取值为准。不要直接把官方 Ag–Pd 的截断距离用于其他体系。

### 10.2 常见拟合方法

|方法名|基本机制|适用判断|
|---|---|---|
|`least-squares`|最小化残差平方和|数据丰富且矩阵良态时作基线|
|`ridge`|二次参数惩罚|希望稳定参数而不要求大量零参数|
|`lasso`|绝对值参数惩罚|生成稀疏模型，需要调 `alpha`|
|`adaptive-lasso`|迭代调整稀疏惩罚权重|比较稀疏解的准确性与稳定性|
|`ardr`|自动相关性判定|通过 `threshold_lambda` 等调稀疏程度|
|`rfe`|递归删除特征|通过 `n_features` 控制模型大小|
|`least-squares-with-reg-matrix`|用户定义二次正则矩阵|表达物理先验和参数耦合|

具体算法参数、标准化和返回属性以安装的 trainstation 文档为准。icet 负责特征与模型定义，trainstation 的全部独立功能不属于 icet 包内 API。

### 10.3 如何验证

均方根误差为：

$$
\mathrm{RMSE}=\sqrt{\frac{1}{M}\sum_i(\hat{y}_i-y_i)^2}.
$$

除训练和验证 RMSE，还应检查预测对参考的散点图、残差随组成和结构来源的变化、极端误差、端点、基态排序和数据量学习曲线。同一结构的变体或重复计算不能跨训练与验证集合造成信息泄漏。按来源、结构族或组成区间做额外独立测试，是本教程建议的验证补充。

低验证 RMSE 不等于相图可靠。两个竞争相的能差可能小于全局误差，也可能由训练范围外的低能构型决定。

### 10.4 模型文件与参数检查

```python
import numpy as np
from icet import ClusterExpansion

ce = ClusterExpansion.read("mixing_energy.ce")
cs_model = ce.get_cluster_space_copy()
test_structure = cs_model.primitive_structure.repeat((2, 2, 2))
prediction = ce.predict(test_structure)
direct = np.dot(cs_model.get_cluster_vector(test_structure), ce.parameters)
assert np.isclose(prediction, direct)
print(ce.to_dataframe().head())
```

表格可同时查看参数和未乘重数 ECI。`prune` 可裁掉零项或指定项，应保留原始模型并核对裁剪前后的预测。任何改变特征定义的操作都要和参数同步，不能仅删除参数数组中的几个数。

来源：[构建教程](https://icet.materialsmodeling.org/get_started/construct_cluster_expansion.html)、[参考数据对比](https://icet.materialsmodeling.org/get_started/compare_to_target_data.html)、[ECI 分析](https://icet.materialsmodeling.org/get_started/analyze_ecis.html)、[回归方法扩展教程](https://ce-tutorials.materialsmodeling.org/part-1/regression-methods.html)、[trainstation](https://trainstation.materialsmodeling.org/)。

## 11 截断选择 超参数扫描与学习曲线

先比较只含对的截断，再加入三体、四体；每次保存参数数、非零项数、训练和验证误差。截断和拟合算法相互影响，最终应联合比较，而不是认为逐阶扫描给出唯一最优解。

```python
# 上下文片段 training_records 是结构与混合能的列表
from trainstation import CrossValidationEstimator

results = []
for cutoffs in ([5.0], [7.0], [7.0, 4.0], [7.0, 4.0, 3.5]):
    cs_trial = ClusterSpace(primitive, cutoffs, ["Ag", "Pd"])
    sc_trial = StructureContainer(cs_trial)
    for structure, label in training_records:
        sc_trial.add_structure(structure, properties={"mixing_energy": label})
    opt_trial = CrossValidationEstimator(
        sc_trial.get_fit_data(key="mixing_energy"), fit_method="ardr",
    )
    opt_trial.validate()
    opt_trial.train()
    results.append((cutoffs, len(cs_trial), opt_trial.rmse_validation))
print(results)
```

超参数扫描复用相同数据、同一验证方案与随机种子。ARDR 扫描 `threshold_lambda`，RFE 扫描 `n_features`，LASSO 和 adaptive-LASSO 扫描 `alpha`。记录算法收敛状态，不能把未收敛模型与正常收敛模型的误差当成同等可靠结果。

```python
# 上下文片段 A y 已由完整训练集获得
records = []
for alpha in (1e-5, 1e-4, 1e-3):
    model = CrossValidationEstimator(
        (A, y), fit_method="lasso", alpha=alpha, max_iter=50000,
    )
    model.validate()
    model.train()
    records.append((alpha, model.rmse_validation, model.n_nonzero_parameters))
```

若多种模型的验证误差接近，优先检查模型稳定性和目标物理量，再考虑参数较少的模型。AIC、BIC 可提供模型复杂度信息，但不是相图正确性的证明。小数据集的验证分割会产生噪声，需要重复划分，并保留不参与调参的独立测试集。

学习曲线比较不同训练集大小下的验证误差。误差仍随数据量明显下降时，继续扩大数据可能有收益；误差进入平台时，应查标签噪声、遗漏的自由度或特征定义，而不是无限增加高阶项。

来源：[Selecting cutoffs](https://icet.materialsmodeling.org/advanced_topics/training_cutoffs_selection.html)、[Hyper-parameter scans](https://icet.materialsmodeling.org/advanced_topics/training_hyper_parameter_scans.html)、[扩展截断教程](https://ce-tutorials.materialsmodeling.org/part-1/cutoff-selection.html)、[扩展超参数教程](https://ce-tutorials.materialsmodeling.org/part-1/hyper-parameter-scans.html)。

## 12 模型集合 不确定性与主动学习

`EnsembleOptimizer` 通过不同训练子集或 bootstrap 产生多个参数向量。对候选结构，计算每个模型的预测，得到均值和分散程度。

```python
from trainstation import EnsembleOptimizer

eopt = EnsembleOptimizer((A, y), ensemble_size=50, fit_method="ardr")
eopt.train()
cv = cs.get_cluster_vector(candidate)
mean, std = eopt.predict(cv, return_std=True)
print("预测均值与模型分散", mean, std)
```

模型集合的不确定性可以写为：

$$
\overline{e}=\frac{1}{K}\sum_{k=1}^K e_k,\qquad
s_e^2=\frac{1}{K}\sum_{k=1}^K(e_k-\overline{e})^2.
$$

这里使用总体方差约定展示概念；实际返回值以 trainstation 实现为准。分散度衡量模型对训练样本变化的敏感性，不能保证捕捉所有系统性偏差。所有模型都漏掉同一种物理机制时，分散也可能很小。

主动学习的操作循环是：建立初始数据集；训练模型集合；生成覆盖组成和低能区的候选池；评估不确定性并去除重复候选；做新 DFT；映射和清洗新结构；重新训练；在独立测试集及目标性质上检查是否收敛。

```python
# 上下文片段 pool 是候选列表 eopt 是已经训练的模型集合
scores = []
for i, structure in enumerate(pool):
    cv = cs.get_cluster_vector(structure)
    mean, std = eopt.predict(cv, return_std=True)
    scores.append((float(std), i, float(mean)))
chosen_indices = [row[1] for row in sorted(scores, reverse=True)[:5]]
```

这个排序示例没有做相似性去重。实际批量选点应避免全选同一结构族，还要兼顾计算成本、低能相关性和空间覆盖。可分别用集合内各 CE 做 MC，观察相界随模型的变化；其成本通常远高于单点不确定性分析。

扩展教程部分旧数据用 `Be` 作为空位占位标签。复现这些数据时应保留其原始约定，新建模型可使用 `X`。不能在读取数据库后擅自替换符号而仍沿用原模型参数。

来源：[EnsembleOptimizer](https://icet.materialsmodeling.org/advanced_topics/training_ensemble_of_models.html)、[Active learning](https://ce-tutorials.materialsmodeling.org/part-1/active-learning.html)。

## 13 约束 加权拟合与低对称模型

### 13.1 精确齐次约束

`Constraints` 表达 $C\boldsymbol{p}=0$。例如让两个已含重数的参数相等：

```python
import numpy as np
from icet.tools import Constraints
from trainstation import Optimizer

C = np.zeros((1, A.shape[1]))
C[0, 2], C[0, 4] = 1.0, -1.0
constraint = Constraints(n_params=A.shape[1])
constraint.add_constraint(C)
opt_c = Optimizer((constraint.transform(A), y), fit_method="ridge")
opt_c.train()
p = constraint.inverse_transform(opt_c.parameters)
assert np.isclose(p[2], p[4])
```

若要约束未乘重数的 ECI 相等，应写成 $p_i/m_i-p_j/m_j=0$，不能无条件使用上面的系数。`get_mixing_energy_constraints(cs)` 可构造纯端点混合能为零的约束，需要确保标签本来就以这些端点作为零参考。

非齐次约束 $C\boldsymbol{p}=\boldsymbol{d}$ 不属于此齐次接口直接表示的形式。可通过变量平移处理，或使用软约束增广设计矩阵，但软约束并不保证精确满足。

### 13.2 权重与性质差约束

若希望最小化 $\sum_i w_i r_i^2$，应把每行矩阵和标签同时乘以 $\sqrt{w_i}$：

```python
weights = np.ones(len(y))
weights[0] = 5.0
row_scale = np.sqrt(weights)
A_weighted = A * row_scale[:, None]
y_weighted = y * row_scale
```

官方示例也直接以行乘数命名“权重”；如果行乘数是 5，对应平方损失的权重实际是 25。设置权重前必须区分这两种定义。

对结构 1 和 2 的能量差，可增广为：

$$
A_{\mathrm{aug}}=\begin{bmatrix}A\\s(\boldsymbol{\Gamma}_1-\boldsymbol{\Gamma}_2)\end{bmatrix},\qquad
\boldsymbol{y}_{\mathrm{aug}}=\begin{bmatrix}\boldsymbol{y}\\s\Delta e_{\mathrm{target}}\end{bmatrix}.
$$

表面偏析能若为总能差，而基础模型输出每位点能，应乘相同超胞的位点数后再定义该约束行。量纲与归一化必须一致。来自同一对结构的增广行应在验证时一起分组，避免性质差泄漏。

### 13.3 合并轨道

低对称表面会产生许多几何相近但晶体学不等价的轨道。`merge_orbits` 可以根据物理假设把这些轨道合并，例如把深层体相区域的同类团簇视作共享 ECI。

```python
# 上下文片段 合并索引必须先针对当前 cs 检查
for i, orbit in enumerate(cs.orbit_list):
    print(i, orbit)
# cs.merge_orbits({2: [3, 4]})
```

合并属于额外建模假设，不是发现了新的严格晶体对称性。修改后必须重新生成团簇向量并拟合。`prune_orbit_list` 则删除几何轨道，其索引与参数数组索引不能混为一谈。

### 13.4 贝叶斯先验与正则矩阵

保留所有轨道而允许相近 ECI 有小幅差异，可采用二次先验：

$$
\mathcal{L}=\|A_M\boldsymbol{J}-\boldsymbol{y}\|^2
 +\boldsymbol{J}^{\mathsf T}\Lambda\boldsymbol{J},\qquad
A_M=A\operatorname{diag}(\boldsymbol{m}).
$$

对两个 ECI 的耦合惩罚 $\lambda(J_i-J_j)^2$，矩阵中增加 $\lambda$ 到两个对角元，并向两个非对角元加入 $-\lambda$。矩阵应对称、半正定。若还给每个 ECI 施加独立惩罚，再添加对角项。

```python
m = np.asarray(cs.get_multiplicities(), dtype=float)
Lambda = np.eye(len(cs)) * 1e-4
i, j, coupling = 2, 3, 1.0  # 示例参数索引 不是任意体系的物理规则
Lambda[i, i] += coupling
Lambda[j, j] += coupling
Lambda[i, j] -= coupling
Lambda[j, i] -= coupling
opt_b = Optimizer(
    (A * m[None, :], y),
    fit_method="least-squares-with-reg-matrix",
    reg_matrix=Lambda,
    standardize=False,
)
opt_b.train()
ce_b = ClusterExpansion(cs, np.asarray(opt_b.parameters) * m)
```

这里把重数移到设计矩阵，使拟合变量是真实 ECI，最后再移回 icet 参数。关闭标准化是为了让本示例先验直接对应所写坐标；如果启用标准化，要核对正则矩阵是否在同一坐标中解释。先验强度必须通过验证选择。该方法得到正则化点估计，并不自动给出完整的后验采样。

来源：[Constraints and weights](https://icet.materialsmodeling.org/advanced_topics/training_with_constraints.html)、[Customizing cluster spaces](https://icet.materialsmodeling.org/advanced_topics/customizing_cluster_spaces.html)、[Bayesian cluster expansions](https://icet.materialsmodeling.org/advanced_topics/training_bayesian_cluster_expansions.html)、[低对称扩展教程](https://ce-tutorials.materialsmodeling.org/part-2/low-symmetry-ce.html)。

## 14 构型预测 凸包与基态

### 14.1 批量预测与零温凸包

```python
import numpy as np
from icet.tools import ConvexHull, enumerate_structures

cs_model = ce.get_cluster_space_copy()
structures = list(enumerate_structures(
    cs_model.primitive_structure, range(1, 7), ["Ag", "Pd"],
))
concentrations = np.array([s.symbols.count("Pd") / len(s) for s in structures])
energies = np.array([ce.predict(s) for s in structures])
hull = ConvexHull(concentrations, energies)
above = hull.get_energy_above_convex_hull(concentrations, energies)
on_hull = hull.is_on_convex_hull(concentrations, energies)
low_indices = hull.extract_low_energy_structures(
    concentrations, energies, energy_tolerance=0.005,
)
print(above, on_hull, low_indices)
```

这里假设模型输出 Ag–Pd 每原子混合能；若 `ce` 来自其他体系，须替换物种及浓度定义。凸包表示候选集合中相分离允许后的最低能包络。凸包以上能不是相对纯端点直线的高度，而是相对当前包络的高度。

小超胞枚举只覆盖有限周期结构。其凸包不证明不存在更大周期、更低能构型，也不能直接当作有限温度相图。

### 14.2 多组分与多子晶格凸包

多组分使用独立组成坐标。多子晶格应保留每个活跃子晶格的组成，而不是仅统计全局元素比例。`get_sublattice_concentrations`、`get_sublattice_site_fractions` 与 `ConvexHull.from_sublattice_concentrations` 用于组织这类数据。

`get_decomposition` 返回目标组成的相分解信息；`get_energy_at_convex_hull` 求包络能；`get_facets` 和 `get_facet_gradients` 给出凸包面及梯度；`get_chemical_potentials` 与 `get_species_chemical_potentials` 提供相应化学势信息。它们的具体结构和标签格式见附录接口链接。缺少端点、组成冗余或退化数据需要先处理，不能只看绘图形状。

### 14.3 混合整数规划基态搜索

当前文档的 `GroundStateFinder` 基于 HiGHS，只支持二元体系。给定超胞和占位数约束后，可优化该有限构型空间。

```python
from icet.tools.ground_state_finder import GroundStateFinder

cell = cs_model.primitive_structure.repeat((3, 3, 3))
finder = GroundStateFinder(ce, cell)
ground = finder.get_ground_state(
    species_count={"Pd": 10}, max_seconds=60, threads=1,
)
print(finder.optimization_status)
print(ce.predict(ground))
```

时间限制下返回候选不等于已证明全局最优。判断 `optimization_status`，并记录超胞、组成、求解器状态和运行限制；所谓全局性也是针对指定超胞和给定 CE。

来源：[Enumerating structures 入门教程](https://icet.materialsmodeling.org/get_started/enumerate_structures.html)、[ConvexHull 和 GroundStateFinder](https://icet.materialsmodeling.org/moduleref/tools.html)。

## 15 SQS 与目标团簇向量结构

SQS 用有限周期结构近似随机合金的团簇相关性。它优化的是指定相关函数集合，不能同时保证所有性质都等同于宏观随机材料。

```python
from ase.build import bulk
from icet import ClusterSpace
from icet.tools.structure_generation import (
    generate_sqs, generate_sqs_from_supercells,
    generate_sqs_by_enumeration, generate_target_structure,
)

prim_sqs = bulk("Au", "fcc", a=4.0)
cs_sqs = ClusterSpace(prim_sqs, [6.0, 4.0], ["Au", "Pd"])
sqs = generate_sqs(
    cs_sqs, max_size=8,
    target_concentrations={"Au": 0.5, "Pd": 0.5},
    include_smaller_cells=False,
    n_steps=10000, random_seed=42,
)
write("sqs.extxyz", sqs)
```

`max_size` 是输入原胞的最大倍数。若原胞有两个位点，倍数 16 对应 32 个位点。`include_smaller_cells=False` 才明确排除更小尺寸。对固定形状使用 `generate_sqs_from_supercells`；小构型空间可用 `generate_sqs_by_enumeration`，其最优性限于枚举范围和采用的目标函数。

目标函数比较团簇向量差，并奖励连续短程团簇的准确匹配；`optimality_weight`、`tol`、团簇截断和退火步数都会影响结果。至少比较多个随机种子，并逐壳检查相关性；最低总分不能替代特定关键壳层的检查。

```python
sqs_fixed = generate_sqs_from_supercells(
    cluster_space=cs_sqs,
    supercells=[prim_sqs.repeat((2, 2, 2))],
    target_concentrations={"Au": 0.5, "Pd": 0.5},
    n_steps=10000, random_seed=7,
)
```

多子晶格目标浓度写为 `{"A": {物种: 分数}, "B": {...}}`，子晶格字母从 `print(cs)` 获得。每个子晶格的分数分别归一化，且与各自位点数相容。

给定非随机目标向量时使用 `generate_target_structure` 或 `generate_target_structure_from_supercells`。向量长度应等于 `len(cs)`，零阶分量应符合恒定项约定，并与指定组成兼容。数学上给出的任意向量可能不可实现，函数只能寻找最接近者。由平均团簇向量生成的代表构型也不能替代完整热力学采样。

来源：[SQS generation](https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html)、[Special structures API](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_sqs)。

## 16 MC 计算器 单位与超胞

### 16.1 每位点性质转换为总量

```python
from mchammer.calculators import ClusterExpansionCalculator

cs_model = ce.get_cluster_space_copy()
structure = cs_model.primitive_structure.repeat((5, 5, 5))
calc = ClusterExpansionCalculator(structure, ce)
total = calc.calculate_total(occupations=structure.numbers)
```

默认 `scaling` 为超胞位点数，因此若 CE 输出 eV/位点，计算器输出 eV 总量，用于接受概率。若训练标签按金属原子、化学式单位或总量归一化，应设置对应缩放。错误缩放会改变能量与 $k_BT$ 的相对量级，从而改变采样温度的物理意义。

有固定子晶格和空位时，缩放不能只根据“看起来有多少原子”决定。例如每金属原子能量要乘金属位点数，而不是包含间隙位点后的总位点数。

### 16.2 自相互作用

```python
print(cs_model.is_supercell_self_interacting(structure))
print(cs_model.are_local_cluster_vectors_additive(structure))
```

超胞太小时，某个团簇可能通过周期镜像多次包含同一个位点。它会影响局部向量分解和局部更新算法。正式 MC 应核对该检查及警告，通常采用更大超胞；若改为缩短截断，必须重新验证和训练模型。

`use_local_energy_calculator=True` 使用局部环境加速性质变化计算。不要把关闭该选项当作普遍解决自相互作用的办法，必须检查所用版本支持行为以及局部与完整计算的一致性。

### 16.3 温度与玻尔兹曼常数

默认能量单位 eV，温度通常 K，玻尔兹曼常数约为 $8.6173\times10^{-5}$ eV/K。若使用无量纲 Ising 模型并设 `boltzmann_constant=1`，温度也成为相应能量标度中的量，不能把输出直接称作真实材料的开尔文温度。

来源：[Calculators](https://icet.materialsmodeling.org/moduleref/calculators.html)、[Cluster vectors](https://icet.materialsmodeling.org/advanced_topics/cluster_vectors.html)、[canonical 扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)。

## 17 正则系综与模拟退火

正则系综固定各物种数量，通过不同占位之间的交换探索构型。对称提议下接受概率为：

$$
P_{\mathrm{acc}}=\min\left(1,\exp[-\Delta E/(k_BT)]\right).
$$

以下独立教学示例使用人工参数，无需 DFT 数据。输出路径放在新目录中，避免与旧数据容器混用。

```python
from pathlib import Path
import numpy as np
from ase.build import bulk
from icet import ClusterSpace, ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import CanonicalEnsemble

run_dir = Path("canonical_demo")
run_dir.mkdir(exist_ok=False)
prim = bulk("Au", "fcc", a=4.0)
cs_demo = ClusterSpace(prim, [3.0], ["Au", "Pd"])
parameters = np.zeros(len(cs_demo))
pairs = [row["index"] for row in cs_demo.as_list if row["order"] == 2]
parameters[pairs[0]] = 0.02
ce_demo = ClusterExpansion(cs_demo, parameters)
supercell = prim.repeat((4, 4, 4))
supercell.symbols = ["Au"] * 32 + ["Pd"] * 32
calculator = ClusterExpansionCalculator(supercell, ce_demo)
mc = CanonicalEnsemble(
    structure=supercell, calculator=calculator,
    temperature=600, random_seed=42,
    dc_filename=str(run_dir / "run.dc"),
    ensemble_data_write_interval=len(supercell),
    trajectory_write_interval=10 * len(supercell),
)
mc.run(100 * len(supercell))
print(mc.data_container.data.tail())
```

一次 MC trial 是一次提议，不论接受还是拒绝都计数。常把 $N$ 次提议称为一个 cycle 或 sweep；多子晶格研究也会以活跃位点数定义 cycle，必须记录所用定义。MC 步数不能直接转换成真实动力学时间。

若全部活跃位点只有同一种占位，就无法交换，正则初始化可能报没有可用交换。固定组成需要确保活跃子晶格上存在可交换物种。

退火使用 `CanonicalAnnealing`，传入 `T_start`、`T_stop`、`n_steps` 和降温函数，然后调用不带步数的 `run()`。`estimated_ground_state` 只是退火找到的最低能候选，不是全局最优证明。多种初态和多个随机种子有助于判断是否陷入局部极小。

来源：[CanonicalEnsemble 和 CanonicalAnnealing](https://icet.materialsmodeling.org/moduleref/ensembles.html)、[正则采样扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)。

## 18 半巨正则系综与成分响应

SGC 固定总位点数和温度，让不同物种与化学势库交换。二元体系取 $\Delta\mu=\mu_B-\mu_A$，广义势为 $E-\Delta\mu N_B$，相应接受概率为：

$$
P_{\mathrm{acc}}=\min\left(1,\exp[-(\Delta E-\Delta\mu\Delta N_B)/(k_BT)]\right).
$$

```python
from mchammer.ensembles import SemiGrandCanonicalEnsemble

# 上下文片段 ce_demo supercell 来自第 17 章
mc_sgc = SemiGrandCanonicalEnsemble(
    supercell,
    ClusterExpansionCalculator(supercell, ce_demo),
    temperature=600,
    chemical_potentials={"Au": 0.0, "Pd": 0.05},
    random_seed=43,
)
mc_sgc.run(100 * len(supercell))
```

化学势用与模型一致的能量标度。若模型拟合的是减去端点参考能后的混合能，输入化学势对应同一参考约定。比较不同模型或实验化学势时须还原参考差。

扫描 `chemical_potentials` 可以获得组成响应。跨一级相变时会跳过某些平均组成，有限采样也可能出现滞后。建议双向扫描、不同初态与边界附近更密集采样；将浓度跳跃作为相界线索，再做收敛分析。

`SGCAnnealing` 结合化学势与降温，适合产生低能候选。与正则退火相同，降温路径中的瞬时统计不自动等同于每个温度的平衡统计。

来源：[Monte Carlo simulations](https://icet.materialsmodeling.org/get_started/run_monte_carlo.html)、[SGC and VCSGC 扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/sgc-vcsgc-simulations.html)、[SGC API](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble)。

## 19 VCSGC 方差约束与混溶间隙

VCSGC 通过 `phis` 控制平均组成、通过 `kappa` 抑制组成涨落，能够连续采样许多 SGC 会跳过的组成。二元单子晶格中，正二次偏置的能量形式为：

$$
\mathcal{H}_{\mathrm{bias}}=E+Nk_BT\kappa(c+\phi/2)^2,
\qquad \rho\propto\exp[-\mathcal{H}_{\mathrm{bias}}/(k_BT)].
$$

因此翻转导致的偏置变化为 $Nk_BT\kappa\Delta c(2c+\phi+\Delta c)$，与普通能量变化一起放入 Metropolis 的负指数。

**符号核对说明：**本次官方系综参考页中概率及接受概率的二次项印为正号，但官方实际 `do_vcsgc_flip` 实现将正二次偏置变化加到能量变化后执行常规接受判断；本文采用该实现对应的抑制方差约定。自由能导数也与官方实现一致。这个修正不能省略，否则会把惩罚偏离平均组成的项写成奖励项。

```python
from mchammer.ensembles import VCSGCEnsemble

mc_vcsgc = VCSGCEnsemble(
    supercell,
    ClusterExpansionCalculator(supercell, ce_demo),
    temperature=600,
    phis={"Pd": -1.0},
    kappa=200,
    random_seed=44,
)
mc_vcsgc.run(100 * len(supercell))
print(mc_vcsgc.data_container.observables)
```

对于每位点自由能 $f=F/N$，有：

$$
\frac{\partial f}{\partial c}=-k_BT\kappa(\phi+2\langle c\rangle).
$$

输出 `free_energy_derivative_Pd` 已是导数标度，在二元单子晶格例子中不要再除以总位点数。`potential` 和物种计数则需要依据定义转换为每位点能与浓度。多个子晶格时，浓度和导数要按相应子晶格定义处理。

当约束足够强时有 $\langle c\rangle\approx-\phi/2$，这是近似关系，不能把 `phi` 当作实际测量浓度。换成另一物种的二元等效参数是 $\phi_B=-2-\phi_A$。一个允许 $n$ 种物种的活跃子晶格应指定 $n-1$ 个 `phi`。

`kappa=200` 是常见保守起点，构造函数要求显式传入它，并无可省略的默认值。更大约束不必然更好：组成分布会变窄，接受率也可能降低，从而增大固定步数下导数误差。扫描少量 `kappa`，比较可达组成、接受率和误差，以选择合适值。

来源：[VCSGC API](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble)、[官方接受算法源码](https://gitlab.com/materials-modeling/icet/-/blob/master/mchammer/ensembles/thermodynamic_base_ensemble.py)、[variance constraint 扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/vcsgc-variance-constraint.html)。

## 20 HybridEnsemble 与位点排斥约束

### 20.1 不同子晶格使用不同系综

金属子晶格可以保持组成不变，氢与空位子晶格则与库交换。`HybridEnsemble` 的 `ensemble_specs` 为每种提议定义系综、子晶格及相关参数。

```python
from mchammer.ensembles import HybridEnsemble

# 上下文片段 cs2 来自第 5 章 ce2 是对应拟合模型
structure2 = cs2.primitive_structure.repeat((4, 4, 4))
subs = cs2.get_sublattices(structure2)
metal = next(i for i, sub in enumerate(subs) if "Au" in sub.chemical_symbols)
hydrogen = next(i for i, sub in enumerate(subs) if "H" in sub.chemical_symbols)
metal_sites = subs[metal].indices
for j, site in enumerate(metal_sites):
    structure2[site].symbol = "Au" if j < len(metal_sites) // 2 else "Pd"
for site in subs[hydrogen].indices:
    structure2[site].symbol = "X"

specs = [
    {"ensemble": "canonical", "sublattice_index": metal},
    {"ensemble": "semi-grand", "sublattice_index": hydrogen,
     "chemical_potentials": {"H": -0.1, "X": 0.0}},
]
mc_hybrid = HybridEnsemble(
    structure2, ClusterExpansionCalculator(structure2, ce2),
    temperature=300, ensemble_specs=specs,
    probabilities=[0.5, 0.5], random_seed=42,
)
mc_hybrid.run(100 * len(structure2))
```

`probabilities` 控制选择各类提议的概率，需要归一化；它是采样效率参数，不是不同物理过程的真实速率。可将某项改为 `ensemble="vcsgc"`，并给 `phis` 和 `kappa`。也可通过 `allowed_sites` 限制表面区域，通过 `allowed_symbols` 指定化学符号；底层提议方法使用的 `allowed_species` 则是原子序数，注意接口层次不同。

### 20.2 邻近位点不得同时占据

`neighbor_sites_to_avoid` 把位点索引映射到禁止同时占据的邻居索引。非 `X` 视为占据。映射必须对称，初态必须满足约束，周期邻居关系必须正确。

```python
from ase.neighborlist import neighbor_list

interstitial = set(subs[hydrogen].indices)
blocking = {int(site): set() for site in interstitial}
i_list, j_list = neighbor_list("ij", structure2, cutoff=3.0)
for i, j in zip(i_list, j_list):
    if i in interstitial and j in interstitial:
        blocking[int(i)].add(int(j))
        blocking[int(j)].add(int(i))
blocking = {site: sorted(neighbors) for site, neighbors in blocking.items()}
# 在 HybridEnsemble 初始化中加入 neighbor_sites_to_avoid=blocking
```

约束通过拒绝违规提议实现。接受准则的正确性并不保证允许构型空间连通；强约束可能让当前交换或翻转无法到达部分合法构型。重启时必须重新提供同一映射，数据容器只记录其摘要而不保存完整映射。

来源：[Hybrid ensembles](https://icet.materialsmodeling.org/advanced_topics/hybrid_ensembles.html)、[HybridEnsemble API](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble)、[ConfigurationManager](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager)。

## 21 Observer 观测器与有序参数

观测器将能量之外的量写入数据容器。应在正式采样前附加，并设置合理的观测间隔。

|观测器|测量对象|关键输入|
|---|---|---|
|`SiteOccupancyObserver`|指定位置组的物种占位率|`cluster_space`、`structure`、`sites`|
|`ShortRangeOrderObserver`|Warren–Cowley 短程有序|`radius`、可选 `pairs`|
|`StructureFactorObserver`|指定波矢的结构因子|`q_points`、可选形状因子|
|`ClusterCountObserver`|各轨道的占位团簇数量|可选 `orbit_indices`|
|`ClusterExpansionObserver`|另一 CE 所描述的性质|另一 `cluster_expansion`|
|`ConstituentStrainObserver`|构型应变贡献|`constituent_strain`|

```python
from mchammer.observers import ShortRangeOrderObserver, ClusterExpansionObserver

obs = ShortRangeOrderObserver(
    cs_demo, supercell, radius=4.1,
    pairs=[("Au", "Pd")], interval=len(supercell),
)
mc.attach_observer(obs)
mc.attach_observer(
    ClusterExpansionObserver(ce_demo, interval=len(supercell)),
    tag="ce_value",
)
mc.run(100 * len(supercell))
print(sorted(mc.data_container.observables))
```

第 17 章先前运行的数据不会自动补上新观测器结果，附加后的数据才有这些字段。另一 CE 观测器的值采用该 CE 的输出归一化，不应默认与 MC `potential` 总量相同。

### 21.1 短程有序

单子晶格中，异种配对的短程有序参数为：

$$
\alpha_{AB}^{(m)}=1-\frac{P_{B|A}^{(m)}}{c_B}.
$$

负值表示异种配对多于随机参考，正值表示异种配对较少。多子晶格实现使用按各子晶格组成构建的随机配对数作为基准，不能直接把总组成代入单子晶格公式。

当前输出键形如 `sro_Au_Pd_1`，最后的数字是对轨道壳层编号。相同距离不一定是同一壳层。某物种缺失时参数可能未定义，输出 `nan`；固定组成的有限单子晶格随机结构也有有限尺寸偏移，异种配对期望偏移为 $-1/(N-1)$，不应要求每个小随机超胞都精确给零。

### 21.2 结构因子与占位率

研究长程有序时选择与有序结构对应的倒空间波矢。`q_points` 应以观测器要求的坐标和单位输入，通常为笛卡尔倒空间坐标，不应直接把未转换的倒格矢分数坐标传入。`form_factors` 可以改变不同物种的散射权重。

位点占位率适合分析表面偏析或指定晶位的占位变化；`sites` 按观测器定义组织成带名称的位置组。初始化后结构的晶胞、位置和位点顺序应保持一致，仅改变占位。

来源：[Observers 全部参考](https://icet.materialsmodeling.org/moduleref/observers.html)、[正则扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)。

## 22 DataContainer 轨迹与统计误差

### 22.1 读取 提取与追加观测

```python
from mchammer import DataContainer

dc = DataContainer.read("canonical_demo/run.dc")
print(dc.ensemble_parameters)
print(dc.metadata)
print(dc.observables)
trial, energy = dc.get("mctrial", "potential", start=20 * len(supercell))
trajectory = dc.get("trajectory")
summary = dc.analyze_data("potential", start=20 * len(supercell))
print(summary)
```

`start` 是 MC trial 编号，而不是数组第几个元素，也不是 cycle 数。乘以所用 cycle 的位点数后再传入。`data` 提供 DataFrame，`get_average` 求均值，`analyze_data` 给均值、标准差、相关长度与均值误差。

`get("trajectory")` 返回 ASE 结构列表，前提是保存了占位轨迹。`apply_observer` 可事后对已存快照补算观测量，其数据只存在于快照对应步；它不会恢复未保存构型的信息。

```python
dc.apply_observer(obs)
dc.write("canonical_demo/run-with-sro.dc")
```

`append` 用于向容器添加步号和记录；手工构造数据时要保持标签类型、步号和间隔一致。不同观测器间隔的字段可能不能在所有步上同时取得，应检查共同采样点和缺失值。

### 22.2 标准差与均值误差不同

标准差描述构型波动；均值误差描述平均值估计精度。相关采样的有效独立观测数小于保存记录数。官方相关长度估计基于自相关函数首次降到 $\exp(-2)$ 以下的滞后；它不是所有统计学定义中的积分自相关时间。

```python
from mchammer.data_analysis import (
    get_autocorrelation_function, get_correlation_length, get_error_estimate,
)

acf = get_autocorrelation_function(energy)
length = get_correlation_length(energy)
error = get_error_estimate(energy)
```

等间隔标量序列才适合当前分析接口。序列太短、常量序列或自相关不衰减可能给出 `None`、`nan` 或不能稳定估计的结果。不能把这些结果改成零误差。

### 22.3 热容与收敛

正则系综中，每位点构型热容为：

$$
c_V=\frac{\langle E^2\rangle-\langle E\rangle^2}{Nk_BT^2}.
$$

此处 $E$ 为总能，$N$ 与报告每位点热容的定义一致。如果使用每位点能量的方差，分子前应乘 $N$ 而不是再除以 $N$。仅包含构型贡献时，不应将此热容直接等同于实验总热容。SGC 或 VCSGC 下浓度可变，不能无条件沿用正则波动公式。

检查能量、组成和序参量的时间序列；丢弃实际平衡期；比较前后分块均值；延长轨迹；换初态和种子；在相变附近提高采样量。保存得更密只增加相关数据，不必然提高有效信息量。

来源：[Data container 进阶专题](https://icet.materialsmodeling.org/advanced_topics/data_container.html)、[Data containers 与分析 API](https://icet.materialsmodeling.org/moduleref/data_containers.html)、[Analyzing Monte Carlo](https://icet.materialsmodeling.org/get_started/analyze_monte_carlo.html)。

## 23 从自由能导数到相图

### 23.1 数据汇总

先从每条 MC 轨迹提取平衡平均与误差，再写入汇总表。二元单子晶格示例如下。

```python
from glob import glob
import pandas as pd
from mchammer import DataContainer

records = []
for filename in sorted(glob("vcsgc-runs/*.dc")):
    dc_run = DataContainer.read(filename)
    nsites = dc_run.ensemble_parameters["n_atoms"]
    start = 50 * nsites  # 示例 需先用时间序列确认
    cstat = dc_run.analyze_data("Pd_count", start=start)
    estat = dc_run.analyze_data("potential", start=start)
    gstat = dc_run.analyze_data("free_energy_derivative_Pd", start=start)
    records.append({
        "filename": filename,
        "temperature": dc_run.ensemble_parameters["temperature"],
        "concentration": cstat["mean"] / nsites,
        "energy_per_site": estat["mean"] / nsites,
        "dfdc": gstat["mean"],
        "dfdc_error": gstat["error_estimate"],
    })
pd.DataFrame(records).to_csv("vcsgc-averages.csv", index=False)
```

这段假设 `Pd_count` 属于唯一活跃二元子晶格。金属加间隙体系必须使用相应子晶格位点数计算组成。若文件不含足够平衡样本，应先延长采样而不是强行输出平均。

### 23.2 数值积分

在每个温度分别按组成排序、处理重复和缺失点，积分 $g(c)=\partial f/\partial c$：

$$
f(c,T)=f(c_0,T)+\int_{c_0}^{c}g(c',T)\,dc'.
$$

```python
import numpy as np
from scipy.integrate import cumulative_trapezoid

# 上下文片段 df_T 是同一温度的数据表
df_T = df_T.sort_values("concentration")
c = df_T["concentration"].to_numpy()
g = df_T["dfdc"].to_numpy()
f_relative = cumulative_trapezoid(g, c, initial=0.0)
```

同一温度下一条自由能曲线的加性常数不影响共切线，但比较不同母晶格或不同相的绝对自由能必须统一参考。还要关注采样区间是否到达端点，以及积分误差与网格依赖。

### 23.3 共切线与相界

二元两相平衡组成 $c_a$ 和 $c_b$ 满足：

$$
f'(c_a)=f'(c_b)=\frac{f(c_b)-f(c_a)}{c_b-c_a}.
$$

自由能的下凸包或共切线可求分解区间。对多个温度重复此步骤得到相界。有限超胞的两相采样含界面贡献和有限尺寸效应，必须比较超胞大小和形状； noisy 导数积分后的细小非凸段不应立即解释为新相。

SGC 可通过组成跳跃和边界追踪定位溶解度，VCSGC 适合补齐跳过的组成及计算导数。有序无序相变还需结构因子、占位率或其他合适序参量，单看短程有序通常不足以准确给出相界。

来源：[MC 分析入门](https://icet.materialsmodeling.org/get_started/analyze_monte_carlo.html)、[SGC VCSGC 与相图扩展教程](https://ce-tutorials.materialsmodeling.org/part-3/sgc-vcsgc-simulations.html)、[正则系综与尺寸效应](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)。

## 24 Wang Landau 采样与态密度

### 24.1 原理与输入

Wang–Landau 估计能量态密度 $g(E)$，在能量空间迭代更新熵估计并使直方图足够平坦。得到 DOS 后，可重加权到不同温度：

$$
Z(T)=\sum_E g(E)\exp[-E/(k_BT)],\qquad
\langle E\rangle=\frac{\sum_E E\,g(E)\exp[-E/(k_BT)]}{Z(T)}.
$$

计算范围取决于提议类型：`swap` 保持组成，`flip` 允许组成变化。DOS 的归一化只确定到乘法常数时，热力学平均仍可获得，但绝对熵与绝对自由能需要状态数或其他参考固定该常数。

### 24.2 独立二维教学例子

```python
from ase import Atoms
import numpy as np
from icet import ClusterSpace, ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import WangLandauEnsemble

prim_wl = Atoms("Au", positions=[[0, 0, 0]], cell=[1, 1, 10], pbc=True)
cs_wl = ClusterSpace(prim_wl, [1.01], ["Au", "Pd"])
params = np.zeros(len(cs_wl))
pair = next(row["index"] for row in cs_wl.as_list if row["order"] == 2)
params[pair] = 2.0  # 人工 Ising 型参数
ce_wl = ClusterExpansion(cs_wl, params)
cell_wl = prim_wl.repeat((4, 4, 1))
cell_wl.symbols = ["Au"] * 8 + ["Pd"] * 8
wl = WangLandauEnsemble(
    cell_wl, ClusterExpansionCalculator(cell_wl, ce_wl),
    energy_spacing=1.0, trial_move="swap",
    fill_factor_limit=1e-4,
    ensemble_data_write_interval=160,
    trajectory_write_interval=1600,
    random_seed=42,
)
wl.run(200000)
print(wl.converged, wl.fill_factor)
```

`energy_spacing` 定义能量格点尺度，`fill_factor_limit` 定义修改因子的停止阈值，`flatness_limit` 和 `flatness_check_interval` 控制平坦性检查。步数只是演示预算，必须检查 `converged`，不能因为 `run()` 返回就宣布 DOS 收敛。

### 24.3 分析与并行窗口

```python
from mchammer.data_containers import (
    get_density_of_states_wl, get_average_observables_wl,
    get_average_cluster_vectors_wl,
)

dos, overlap_errors = get_density_of_states_wl(wl.data_container)
thermal = get_average_observables_wl(
    wl.data_container,
    temperatures=np.linspace(0.5, 4.0, 10),
    boltzmann_constant=1.0,
)
print(dos.head(), thermal.head())
```

`WangLandauDataContainer` 额外保存熵、直方图和修改因子历史；`fill_factor_limit` 可选择满足某阶段精度的历史信息。`get_average_cluster_vectors_wl` 可计算热平均团簇向量，再用于代表构型生成。

大型系统可用 `energy_limit_left/right` 划分有重叠的能量窗口，`get_bins_for_parallel_simulations` 辅助分窗。合并时，各窗口必须在相同能量网格、同一体系和连续重叠区上相容；输出的重叠误差用于检查拼接质量。只采低能区可以减少成本，但会限制可靠重加权的最高温度。

官方以二维 Ising 模型展示 DOS、热容、代表构型与尺寸效应。固定组成 `swap` 模型与允许磁化变化的 `flip` 模型对应不同状态空间，不应混用两者的精确结果作验证。

来源：[Wang-Landau simulations](https://icet.materialsmodeling.org/advanced_topics/wang_landau_simulations.html)、[WangLandauEnsemble](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble)、[WL 数据分析](https://icet.materialsmodeling.org/moduleref/data_containers.html)。

## 25 热力学积分与温度积分

### 25.1 热力学路径

从已知自由能参考哈密顿量 $H_A$ 到目标 $H_B$，定义 $H(\lambda)=(1-\lambda)H_A+\lambda H_B$，则：

$$
F_B-F_A=\int_0^1\left\langle\frac{\partial H}{\partial\lambda}\right\rangle_\lambda\,d\lambda.
$$

对于固定组成的无相互作用参考，其构型熵可从状态数计算：

$$
S_{\mathrm{ideal}}=k_B\ln\Omega,\qquad
\Omega=\frac{N!}{\prod_s N_s!},\qquad F_A=-TS_{\mathrm{ideal}}.
$$

多个独立且各自固定组成的子晶格，状态数为各子晶格状态数的乘积。使用位点阻塞等约束后，不受约束的理想混合状态数可能已不适用。

### 25.2 ThermodynamicIntegrationEnsemble

以下片段采用第 24 章人工模型并使用无量纲温度。

```python
from mchammer.ensembles import CanonicalEnsemble, ThermodynamicIntegrationEnsemble
from mchammer.free_energy_tools import get_free_energy_thermodynamic_integration

equil = CanonicalEnsemble(
    cell_wl, ClusterExpansionCalculator(cell_wl, ce_wl),
    temperature=20.0, boltzmann_constant=1.0,
)
equil.run(5000)
ti = ThermodynamicIntegrationEnsemble(
    equil.structure, ClusterExpansionCalculator(equil.structure, ce_wl),
    temperature=1.0, n_steps=100000, forward=True,
    boltzmann_constant=1.0, ensemble_data_write_interval=10,
)
ti.run()
temperatures, free_energy = get_free_energy_thermodynamic_integration(
    ti.data_container, cs_wl, forward=True,
    max_temperature=5.0, boltzmann_constant=1.0,
)
```

该类实现官方教程使用的参考路径，不是任意两个外部计算器插值的通用接口。计算器和后处理的玻尔兹曼常数、组成、子晶格及归一化应一致。正反向积分及更长路径可检查非平衡偏差；简单平均正反结果不能保证消除所有路径滞后。

### 25.3 温度积分

已知一个温度的自由能时：

$$
\frac{F(T_2)}{T_2}=\frac{F(T_1)}{T_1}
-\int_{T_1}^{T_2}\frac{U(T)}{T^2}\,dT.
$$

用 `CanonicalAnnealing` 得到缓慢温度路径，再调用 `get_free_energy_temperature_integration`，传入 `temperature_reference`，必要时传入 `free_energy_reference`。默认理想混合参考只适用于足够高温并满足相应状态计数假设的情况。

```python
from mchammer.ensembles import CanonicalAnnealing
from mchammer.free_energy_tools import get_free_energy_temperature_integration

anneal = CanonicalAnnealing(
    equil.structure, ClusterExpansionCalculator(equil.structure, ce_wl),
    T_start=20.0, T_stop=1.0, n_steps=100000,
    cooling_function="linear", boltzmann_constant=1.0,
    ensemble_data_write_interval=10,
)
anneal.run()
T_grid, F_grid = get_free_energy_temperature_integration(
    anneal.data_container, cs_wl, forward=True,
    temperature_reference=20.0, boltzmann_constant=1.0,
)
```

温度参考必须与路径端点对应。在相变附近增加路径长度或采用更充分平衡的温度网格，并用小体系枚举配分函数校验自由能定义与常数。

来源：[Thermodynamic and temperature integration](https://icet.materialsmodeling.org/advanced_topics/free_energy_canonical_ensemble.html)、[TI 及自由能函数 API](https://icet.materialsmodeling.org/moduleref/data_containers.html)。

## 26 ConstituentStrain 长程应变

普通有限截断 CE 难以表达相干相分离产生的长程弹性贡献。`ConstituentStrain` 通过倒空间构型结构因子和方向、组成相关的应变能函数引入该贡献。

总模型分为：

$$
e_{\mathrm{total}}(\boldsymbol{\sigma})
=e_{\mathrm{short}}(\boldsymbol{\sigma})+e_{\mathrm{CS}}(\boldsymbol{\sigma}).
$$

拟合短程 CE 前从参考标签中减去已定义的应变项；MC 再通过 `ConstituentStrainCalculator` 把它加回。若直接用已包含该应变贡献的标签拟合，又在采样时额外添加，会造成重复计数。

用户需要自己提供材料的 `strain_energy_function`，必要时提供 `k_to_parameter_function`。方向与组成依赖可以来自弹性理论或额外计算后的插值。包不会自动替用户计算刚度张量或完整方向应变数据库。

```python
# 上下文片段 strain_energy_function 是材料相关函数 必须自行提供
from icet.tools import ConstituentStrain
from mchammer.calculators import ConstituentStrainCalculator
from mchammer.observers import ConstituentStrainObserver

strain = ConstituentStrain(
    supercell=structure,
    primitive_structure=cs_model.primitive_structure,
    chemical_symbols=["Ag", "Cu"],
    concentration_symbol="Cu",
    strain_energy_function=strain_energy_function,
)
strain_value = strain.get_constituent_strain(structure.numbers)
# ce_short 须与上述 Ag Cu 结构相容 并已扣除应变贡献
calc_strain = ConstituentStrainCalculator(strain, ce_short)
```

`damping` 与材料应变函数的定义须联合检查，不能照抄其他材料的值。该模块官方只在 FCC 体系中测试，并提示其他晶体尤其非立方体系需单独验证。文档中的应变计算器仅适用于单点翻转，不能直接接正则交换系综。应变观测器可分离该贡献；倒空间非局部计算会增加 MC 成本。

底层 `KPoint`、Redlich–Kister 辅助函数及参数变换接口见附录。使用高级功能前，应先完成普通 CE 的数据和采样验证。

来源：[Constituent strain calculations](https://icet.materialsmodeling.org/advanced_topics/constituent_strain.html)、[ConstituentStrain API](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain)、[ConstituentStrainCalculator](https://icet.materialsmodeling.org/moduleref/calculators.html)。

## 27 局部团簇向量与自定义采样

### 27.1 局部向量

`cs.get_cluster_vector(structure)` 返回全局向量。编译计算器还提供局部向量，适合分析位点贡献和构造局部性质模型。

```python
# 上下文片段 calculator 是当前超胞的计算器
cpp_calc = calculator.cpp_calc
cpp_calc.set_occupations(supercell.numbers)
global_cv = cpp_calc.get_cluster_vector()
local_sum = sum(
    cpp_calc.get_local_cluster_vector(i) for i in range(len(supercell))
)
```

必须先检查 `cs.are_local_cluster_vectors_additive(supercell)`。同一团簇内重复周期位点会破坏简单可加关系。当前局部贡献按团簇阶数分摊，零阶项每位点贡献为 $1/N$。把局部向量用作原子能模型时，也要核对目标局部能量定义是否有唯一、可比较的分解。

### 27.2 高层计算器维护状态

计算器持有当前占位，手动调用时先同步，再评价提议，只有接受后才更新内部状态。

```python
import numpy as np

occ = supercell.numbers.copy()
calculator.set_occupations(occ)
a = 0
b = next(i for i, value in enumerate(occ) if value != occ[a])
sites = [a, b]
new = [int(occ[b]), int(occ[a])]
delta = calculator.calculate_change(
    sites=sites, current_occupations=occ,
    new_site_occupations=new,
)
new_occ = occ.copy()
new_occ[sites] = new
before = calculator.calculate_total(occupations=occ)
after = calculator.calculate_total(occupations=new_occ)
assert np.isclose(delta, after - before)
calculator.accept_change(sites=sites, species=new)
occ = new_occ
```

这段仅展示状态合同，不是完整 MC 算法。接受或拒绝必须按所研究系综的概率准则决定。调用 `calculate_change` 不会提交移动；只改变外部数组而忘记 `accept_change` 会使后续增量计算依据旧构型。

### 27.3 批量移动评价

底层一个 move 是按顺序执行的 `(位点索引, 新原子序数)` 列表。`get_cluster_vector_changes(moves)` 批量评价多个候选，候选之间独立且不提交状态；`apply_moves` 才提交选中移动。

```python
moves = [
    [(a, int(occ[b]))],
    [(a, int(occ[b])), (b, int(occ[a]))],
]
changes = cpp_calc.get_cluster_vector_changes(moves)
# 点乘得到的是当前计算器内部 CE 参数对应的每位点变化
p_internal = np.asarray(calculator.cluster_expansion.parameters)
delta_per_site = np.asarray(changes) @ p_internal
# 总能量还须按当前 scaling 转换 不能遗漏
```

内部计算器可能删除全零参数轨道，因此这里用 `calculator.cluster_expansion.parameters`，不无条件用原始 CE 参数。该属性返回副本，修改副本不会改变计算器。为了避免不一致，改变模型参数后重新创建计算器。

自定义带缓存的计算器必须同时实现绝对设置 `set_occupations` 和增量接受 `accept_change`；重启时只有后者不够。自定义观测器遵循 `BaseObserver` 的返回类型、间隔和 `get_observable` 合同。自定义采样中的提议若不对称，还需要相应的概率比修正。

底层状态读写的线程安全约定与是否实际并行不同。当前绑定调用持有 Python GIL，Python 线程不自动带来 C++ 计算加速；各进程创建独立计算器更适合并行采样。

来源：[Using a calculator directly](https://icet.materialsmodeling.org/advanced_topics/calculator_interface.html)、[Cluster vectors](https://icet.materialsmodeling.org/advanced_topics/cluster_vectors.html)、[Calculators](https://icet.materialsmodeling.org/moduleref/calculators.html)。

## 28 并行 断点续算与文件管理

不同温度、化学势、随机种子或模型通常可以独立运行。Windows 下必须把多进程入口放在 `if __name__ == "__main__":` 保护内；不要依赖 Linux 的进程复制行为。

```python
from multiprocessing import Pool
from pathlib import Path
from icet import ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import SemiGrandCanonicalEnsemble

def run_case(task):
    index, temperature, dmu = task
    ce_job = ClusterExpansion.read("mixing_energy.ce")
    cs_job = ce_job.get_cluster_space_copy()
    cell = cs_job.primitive_structure.repeat((5, 5, 5))
    calc_job = ClusterExpansionCalculator(cell, ce_job)
    mc_job = SemiGrandCanonicalEnsemble(
        cell, calc_job, temperature=temperature,
        chemical_potentials={"Ag": 0.0, "Pd": dmu},
        random_seed=1000 + index,
        dc_filename=f"parallel-runs/case-{index:04d}.dc",
    )
    mc_job.run(500 * len(cell))
    return index, mc_job.step

if __name__ == "__main__":
    Path("parallel-runs").mkdir(exist_ok=False)
    tasks = [(0, 600, -0.1), (1, 600, 0.1), (2, 900, -0.1)]
    with Pool(processes=2) as pool:
        result = pool.map_async(run_case, tasks)
        print(result.get())
```

每个进程在内部读取模型和创建计算器，写不同文件。官方建议使用 `map_async` 或 `starmap_async`，本示例沿用这一建议。进程数应同时考虑内存和底层数值库线程数，避免每个进程各自启动大量线程。

已有 `dc_filename` 可能触发读取旧容器并继续写入。开展新实验时使用新目录；真正续算时保持模型、结构、系综、参数、位点阻塞映射等一致，并核对步计数与旧记录。不要在改变温度后仍用原文件名假装是同一条轨迹。

`data_container_write_period` 的单位是秒，控制落盘周期；`ensemble_data_write_interval` 和 `trajectory_write_interval` 的单位是 MC trial，控制记录频率。它们作用不同。当前接口的 `None` 通常表示采用默认间隔，不能把它无条件解释为禁止保存轨迹。

需要连续扫描且前一个平衡构型有助于下一点平衡时，可按温度并行、每个温度内部顺序扫描化学势。该选择改变初始条件与效率，仍要检查滞后。

来源：[Parallel Monte Carlo simulations](https://icet.materialsmodeling.org/advanced_topics/parallel_monte_carlo.html)、[系综文件与记录参数](https://icet.materialsmodeling.org/moduleref/ensembles.html)、[Wang-Landau 重启说明](https://icet.materialsmodeling.org/advanced_topics/wang_landau_simulations.html)。

## 29 一套可复现的研究项目

以下结构是本教程的项目组织建议，可按计算集群调整。

```text
project/
  environment.yml
  config.json
  structures/primitive.extxyz
  structures/ideal/
  structures/relaxed/
  reference/reference_data.db
  models/training.sc
  models/mixing_energy.ce
  scripts/01_prepare.py
  scripts/02_fit.py
  scripts/03_validate.py
  scripts/04_run_mc.py
  scripts/05_analyze.py
  runs/canonical/
  runs/sgc/
  runs/vcsgc/
  results/metrics.csv
  results/free_energy.csv
  results/phase_boundaries.csv
```

完整实施按以下阶段推进。

1. **定义模型。** 固定母晶格、允许占位和参考态，记录能量归一化及目标研究范围。
2. **建立候选。** 结合端点、有序结构、随机结构和结构选择生成计算集。
3. **外部计算与映射。** 使用一致的 DFT 设置，剔除无有效晶格映射或标签不可信的记录。
4. **选择和训练。** 比较截断、算法与超参数，保存容器、模型和验证结果。
5. **物理验证。** 检查低能结构排序、凸包、组成外推、目标性质差和参数稳定性。
6. **采样设计。** 选合适系综，核对缩放与超胞，设置观测器、种子和独立输出路径。
7. **统计分析。** 按实际平衡期截取，检查自相关、误差、初态与尺寸依赖。
8. **相图或性质报告。** 统一自由能参考，报告收敛证据、模型敏感性和遗漏的自由度。

研究记录至少包括：icet 与依赖版本、原胞结构和物种定义、截断及容差、结构映射规则、数据筛选、拟合与验证方法、随机种子、MC 超胞、各子晶格组成、能量缩放、采样和记录间隔、平衡期及不确定性方法。

在官方实战中可顺序运行 `tutorial.zip` 的构建、参考对比、枚举、ECI 分析、SGC/VCSGC、数据汇总和绘图脚本；进阶资料在 [advanced_topics.zip](https://icet.materialsmodeling.org/advanced_topics.zip)。扩展 Notebook 与数据在 [Zenodo 官方教程资料](https://doi.org/10.5281/zenodo.10997197)。不同资料包的相对路径与空位编码可能不同，先检查 README 和文件树。

## 30 常见问题与排错

|现象|优先检查|处理方向|
|---|---|---|
|`Failed to find site by position`|是否为同一母晶格超胞、是否已弛豫|映射回理想晶格，核对晶胞，必要时 `wrap`|
|映射失败或诊断位移过大|变换矩阵、空位、间隙位点、应变|查看 DEBUG 日志，补齐参考位点，检查模型适用性|
|矩阵条件数很大|组成覆盖、冗余列、重复特征|改善候选集，删除冗余自由度或构造合理约束|
|训练误差低 验证误差高|过拟合或数据泄漏|独立测试，缩减基组，重新选择正则与数据|
|纯相混合能不为零|参考定义、端点数据、约束|核对两项参考能，必要时端点精确约束|
|参数数对不上|多组分组合、轨道合并、模型裁剪|使用匹配的 `ClusterSpace` 和参数，不手工猜索引|
|SQS 组成不可实现|浓度乘子晶格位点数是否为整数|改用相容尺寸或组成|
|SQS 晶胞形状奇怪|是否允许搜索多个超胞形状|改用固定超胞搜索接口|
|没有可用 canonical swaps|活跃子晶格是否只有一种占位|设置可交换混合占位或选择适合的系综|
|MC 温度行为不合理|每位点与总能混用、能量单位|核对 `scaling` 与 `boltzmann_constant`|
|MC 自相互作用警告|超胞与截断相对大小|扩大超胞，验证局部计算条件|
|SGC 扫描中浓度突跳|真实两相区或未平衡滞后|双向和长轨迹验证，必要时使用 VCSGC|
|VCSGC 导数数量级错|导数被再次除位点数、子晶格分母错|按对应导数定义处理|
|SRO 为 `nan`|某物种缺失或参考配对数为零|报告未定义，不替换为零|
|小随机超胞 SRO 不为零|固定组成有限尺寸偏移|比较同尺寸随机参考和热力学极限|
|误差估计 `None` 或 `nan`|序列太短、常量、ACF 未衰减|延长采样并检查分析对象|
|事后观测器没有足够记录|轨迹是否保存、记录间隔|重算已有快照，缺失轨迹需重新采样|
|WL DOS 不收敛|能量网格、不可达能量、窗口、初态|核对谱与平坦性，分窗并延长采样|
|并行脚本反复启动或挂起|Windows 主入口保护、共享状态|独立进程内初始化，使用主入口和异步映射|
|输出混入旧轨迹|文件存在触发续算|新实验新目录，续算核对原始参数|
|测试失败|源与测试版本、依赖版本|在一致环境运行对应版本测试|

官方 FAQ 说明主要充分测试的平台是 Linux 和 macOS；Windows 使用中要额外核对编译和进程入口。软件答疑入口是 [matsci.org icet 讨论区](https://matsci.org/c/icet/)。

来源：[FAQ](https://icet.materialsmodeling.org/backmatter/faq.html)、[结构映射](https://icet.materialsmodeling.org/advanced_topics/mapping_structures.html)、[SRO 参考](https://icet.materialsmodeling.org/moduleref/observers.html)、[MC 分析 API](https://icet.materialsmodeling.org/moduleref/data_containers.html)。

## 31 术语 功能边界与学术引用

|术语|含义|
|---|---|
|CE|Cluster expansion 团簇展开|
|cluster|母晶格上的有限位点集合|
|order|团簇包含的位点数|
|orbit|晶体对称操作下等价团簇的集合|
|cluster space|允许的团簇及多组分基函数定义|
|cluster vector|一个占位结构在该空间中的相关函数向量|
|zerolet|空团簇，对应恒为 1 的截距分量|
|singlet pair triplet|单点 二体 三体团簇|
|radius|团簇位点到几何中心的平均距离|
|cutoff|指定阶数团簇内允许的最大位点间距|
|multiplicity|轨道或参数项对应的重数|
|ECI|有效团簇相互作用，与拟合参数的重数约定须区分|
|sublattice|具有相应允许占位的子晶格位点集合|
|active inactive|允许多种占位的活跃位点与单一占位的固定位点|
|DFT|密度泛函理论，常用参考数据来源|
|FCC BCC|面心立方 体心立方晶格|
|SQS|特殊准随机结构|
|MC trial|一次蒙特卡洛提议|
|MCS sweep cycle|按约定数量位点对应的一组 MC trials|
|canonical|固定各物种数量的正则系综|
|SGC|固定化学势差的半巨正则系综|
|VCSGC|方差约束半巨正则系综|
|SRO|短程有序|
|DOS|构型能量态密度，此处不是电子态密度|
|WL|Wang–Landau 采样|
|TI TEI|热力学积分 温度积分|
|convex hull|组成与能量或自由能的下凸包|
|ground state|给定模型及构型空间内的最低能状态|

软件的固定晶格占位模型不自动包含位移自由度；由弛豫能建立 CE 也不意味着 MC 在每次提议中进行晶格弛豫。跨母晶格相竞争需要分别建模并统一参考；实验相图比较还需按研究目标考虑振动、磁性、电子、弹性和压力等贡献。ConstituentStrain 是特定长程应变模型，不等同于完整结构力学求解器。

使用 icet 应引用软件论文；以扩展教程为工作基础时同时引用教程论文，并针对枚举、SQS、VCSGC、MIP、贝叶斯或应变等功能引用对应方法。主要入口如下。

- 软件论文：M. Ångqvist 等，*icet – A Python Library for Constructing and Sampling Alloy Cluster Expansions*，Advanced Theory and Simulations，2019，[DOI 10.1002/adts.201900015](https://doi.org/10.1002/adts.201900015)。
- 教程论文：P. Ekborg-Tanner、P. Rosander、E. Fransson、P. Erhart，*Construction and sampling of alloy cluster expansions — A tutorial*，PRX Energy，2024，[DOI 10.1103/PRXEnergy.3.042001](https://doi.org/10.1103/PRXEnergy.3.042001)。
- 功能方法与依赖软件引用：[Credits](https://icet.materialsmodeling.org/credits.html)。
- 完整官方参考书目：[Bibliography](https://icet.materialsmodeling.org/backmatter/bibliography.html)。

正式投稿时从 DOI 或期刊导出文献条目，检查作者、卷页和所用方法。也应给 ASE、spglib、trainstation 和相关数值库适当引用。

来源：[Glossary](https://icet.materialsmodeling.org/backmatter/glossary.html)、[Credits](https://icet.materialsmodeling.org/credits.html)、[Bibliography](https://icet.materialsmodeling.org/backmatter/bibliography.html)。

## 附录 A 公开 API 完整检索

以下索引按本次主站 9 个模块参考页的公开声明生成，保留函数签名、默认值、属性与数据字段，并链接到相应官方锚点。继承方法会在多个类中重复出现，这些记录表示各类的公开可见接口。类型注解来自文档，不应把 `None` 或默认值解释为超出文档的行为。

本教程正文说明核心接口的操作与风险，索引用来查找全部声明及参数详情。内部 C++ 类与未正式列入参考页的接口不作为稳定公共 API 承诺；第 27 章另说明官方进阶文档公开的编译计算器调用方式。

### Calculators

[官方模块说明](https://icet.materialsmodeling.org/moduleref/calculators.html)

- [`class mchammer.calculators.ClusterExpansionCalculator(structure, cluster_expansion, name='Cluster Expansion Calculator', scaling=None, use_local_energy_calculator=True)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator)
- [`accept_change(*, sites=None, species=None)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.accept_change)
- [`calculate_change(*, sites, current_occupations, new_site_occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.calculate_change)
- [`calculate_total(*, occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.calculate_total)
- [`property cluster_expansion: ClusterExpansion`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.cluster_expansion)
- [`set_occupations(occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.set_occupations)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ClusterExpansionCalculator.sublattices)
- [`class mchammer.calculators.ConstituentStrainCalculator(constituent_strain, cluster_expansion, name='Constituent Strain Calculator', scaling=None)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator)
- [`accept_change(*, sites=None, species=None)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.accept_change)
- [`calculate_change(*, sites, current_occupations, new_site_occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.calculate_change)
- [`calculate_total(*, occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.calculate_total)
- [`property cluster_expansion: ClusterExpansion`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.cluster_expansion)
- [`set_occupations(occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.set_occupations)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.ConstituentStrainCalculator.sublattices)
- [`class mchammer.calculators.TargetVectorCalculator(structure, cluster_space, target_vector, weights=None, optimality_weight=1.0, optimality_tol=1e-05, name='Target vector calculator')`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.TargetVectorCalculator)
- [`accept_change(*, sites=None, species=None)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.TargetVectorCalculator.accept_change)
- [`calculate_total(occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.TargetVectorCalculator.calculate_total)
- [`set_occupations(occupations)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.TargetVectorCalculator.set_occupations)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.TargetVectorCalculator.sublattices)
- [`mchammer.calculators.compare_cluster_vectors(cv_1, cv_2, as_list, weights=None, optimality_weight=1.0, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/calculators.html#mchammer.calculators.compare_cluster_vectors)

### Cluster expansions

[官方模块说明](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)

- [`class icet.ClusterExpansion(cluster_space, parameters, metadata=None)`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`cluster_space`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`parameters`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`metadata`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property chemical_symbols: list[list[str]]`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`copy()`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property cutoffs: list[float]`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property fractional_position_tolerance: float`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`get_cluster_space_copy()`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property metadata: dict`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property orders: list[int]`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property parameters: list[float]`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property position_tolerance: float`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`predict(structure)`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property primitive_structure: Atoms`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`prune(indices=None, tol=0)`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`static read(filename)`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`property symprec: float`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`to_dataframe()`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)
- [`write(filename)`](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)

### Cluster space

[官方模块说明](https://icet.materialsmodeling.org/moduleref/cluster_space.html)

- [`class icet.ClusterSpace(structure, cutoffs, chemical_symbols, symprec=1e-05, position_tolerance=None)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace)
- [`are_local_cluster_vectors_additive(structure)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.are_local_cluster_vectors_additive)
- [`property as_list: list[dict]`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.as_list)
- [`assert_structure_compatibility(structure, vol_tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.assert_structure_compatibility)
- [`property chemical_symbols: list[list[str]]`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.chemical_symbols)
- [`copy()`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.copy)
- [`property cutoffs: list[float]`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.cutoffs)
- [`property fractional_position_tolerance: float`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.fractional_position_tolerance)
- [`get_cluster_vector(structure)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.get_cluster_vector)
- [`get_coordinates_of_representative_cluster(orbit_index)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.get_coordinates_of_representative_cluster)
- [`get_multiplicities()`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.get_multiplicities)
- [`get_possible_orbit_occupations(orbit_index)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.get_possible_orbit_occupations)
- [`get_sublattices(structure)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.get_sublattices)
- [`is_supercell_self_interacting(structure)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.is_supercell_self_interacting)
- [`merge_orbits(equivalent_orbits, ignore_permutations=False)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.merge_orbits)
- [`property number_of_orbits_by_order: dict`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.number_of_orbits_by_order)
- [`property orbit_list`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.orbit_list)
- [`property position_tolerance: float`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.position_tolerance)
- [`property primitive_structure: Atoms`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.primitive_structure)
- [`prune_orbit_list(indices)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.prune_orbit_list)
- [`static read(filename)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.read)
- [`property space_group: str`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.space_group)
- [`property symprec: float`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.symprec)
- [`to_dataframe()`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.to_dataframe)
- [`write(filename)`](https://icet.materialsmodeling.org/moduleref/cluster_space.html#icet.ClusterSpace.write)

### Core components

[官方模块说明](https://icet.materialsmodeling.org/moduleref/core.html)

- [`class icet.core.sublattices.Sublattices(allowed_species, primitive_structure, structure, fractional_position_tolerance)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices)
- [`property active_sublattices: list[Sublattice]`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.active_sublattices)
- [`property allowed_species: list[list[str]]`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.allowed_species)
- [`assert_occupation_is_allowed(chemical_symbols)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.assert_occupation_is_allowed)
- [`get_allowed_numbers_on_site(index)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.get_allowed_numbers_on_site)
- [`get_allowed_symbols_on_site(index)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.get_allowed_symbols_on_site)
- [`get_sublattice_index_from_site_index(index)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.get_sublattice_index_from_site_index)
- [`get_sublattice_sites(index)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.get_sublattice_sites)
- [`property inactive_sublattices: list[Sublattice]`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattices.inactive_sublattices)
- [`class icet.core.sublattices.Sublattice(chemical_symbols, indices, symbol)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattice)
- [`property symbol`](https://icet.materialsmodeling.org/moduleref/core.html#icet.core.sublattices.Sublattice.symbol)
- [`icet.tools.variable_transformation.get_transformation_matrix(structure, full_orbit_list)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.tools.variable_transformation.get_transformation_matrix)
- [`icet.tools.variable_transformation.transform_parameters(structure, full_orbit_list, parameters)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.tools.variable_transformation.transform_parameters)
- [`class icet.tools.constituent_strain.KPoint(kpt, multiplicity, structure_factor, strain_energy_function, damping)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.tools.constituent_strain.KPoint)
- [`icet.tools.constituent_strain_helper_functions.redlich_kister(x, *coeffs)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.tools.constituent_strain_helper_functions.redlich_kister)
- [`icet.tools.constituent_strain_helper_functions.redlich_kister_vector(x, *coeffs)`](https://icet.materialsmodeling.org/moduleref/core.html#icet.tools.constituent_strain_helper_functions.redlich_kister_vector)
- [`class mchammer.ConfigurationManager(structure, sublattices)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager)
- [`get_flip_state(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.get_flip_state)
- [`get_occupations_on_sublattice(sublattice_index)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.get_occupations_on_sublattice)
- [`get_swapped_state(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.get_swapped_state)
- [`is_constraint_violated(sites, species)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.is_constraint_violated)
- [`is_swap_possible(sublattice_index, allowed_species=None)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.is_swap_possible)
- [`property neighbor_sites_to_avoid: dict[int, list[int]] | None`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.neighbor_sites_to_avoid)
- [`property occupations: ndarray`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.occupations)
- [`set_neighbor_sites_to_avoid(neighbor_sites_to_avoid)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.set_neighbor_sites_to_avoid)
- [`set_occupations(occupations)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.set_occupations)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.sublattices)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.update_occupations)
- [`validate_constraint(neighbor_sites_to_avoid)`](https://icet.materialsmodeling.org/moduleref/core.html#mchammer.ConfigurationManager.validate_constraint)

### Data containers

[官方模块说明](https://icet.materialsmodeling.org/moduleref/data_containers.html)

- [`class mchammer.DataContainer(structure, ensemble_parameters, metadata={})`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer)
- [`analyze_data(tag, start=None, max_lag=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.analyze_data)
- [`append(mctrial, record)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.append)
- [`apply_observer(observer)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.apply_observer)
- [`property data: DataFrame`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.data)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.ensemble_parameters)
- [`get(*tags, start=0)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.get)
- [`get_average(tag, start=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.get_average)
- [`property metadata: dict`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.metadata)
- [`property observables: list[str]`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.observables)
- [`classmethod read(infile)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.read)
- [`write(outfile)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.DataContainer.write)
- [`class mchammer.WangLandauDataContainer(structure, ensemble_parameters, metadata={})`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer)
- [`append(mctrial, record)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.append)
- [`apply_observer(observer)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.apply_observer)
- [`property data: DataFrame`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.data)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.ensemble_parameters)
- [`property fill_factor: float`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.fill_factor)
- [`property fill_factor_history: DataFrame`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.fill_factor_history)
- [`get(*tags, fill_factor_limit=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.get)
- [`get_entropy(fill_factor_limit=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.get_entropy)
- [`get_histogram()`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.get_histogram)
- [`property metadata: dict`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.metadata)
- [`property observables: list[str]`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.observables)
- [`classmethod read(infile)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.read)
- [`write(outfile)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.WangLandauDataContainer.write)
- [`mchammer.data_containers.get_average_observables_wl(dcs, temperatures, observables=None, boltzmann_constant=8.617330337217213e-05, fill_factor_limit=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_containers.get_average_observables_wl)
- [`mchammer.data_containers.get_average_cluster_vectors_wl(dcs, cluster_space, temperatures, boltzmann_constant=8.617330337217213e-05, fill_factor_limit=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_containers.get_average_cluster_vectors_wl)
- [`mchammer.data_containers.get_density_of_states_wl(dcs, fill_factor_limit=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_containers.get_density_of_states_wl)
- [`mchammer.free_energy_tools.get_free_energy_thermodynamic_integration(dc, cluster_space, forward, max_temperature=inf, sublattice_probabilities=None, boltzmann_constant=8.617330337217213e-05)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.free_energy_tools.get_free_energy_thermodynamic_integration)
- [`mchammer.free_energy_tools.get_free_energy_temperature_integration(dc, cluster_space, forward, temperature_reference, free_energy_reference=None, sublattice_probabilities=None, max_temperature=inf, boltzmann_constant=8.617330337217213e-05)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.free_energy_tools.get_free_energy_temperature_integration)
- [`mchammer.data_analysis.analyze_data(data, max_lag=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_analysis.analyze_data)
- [`mchammer.data_analysis.get_autocorrelation_function(data, max_lag=None)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_analysis.get_autocorrelation_function)
- [`mchammer.data_analysis.get_correlation_length(data)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_analysis.get_correlation_length)
- [`mchammer.data_analysis.get_error_estimate(data, confidence=0.95)`](https://icet.materialsmodeling.org/moduleref/data_containers.html#mchammer.data_analysis.get_error_estimate)

### Ensembles

[官方模块说明](https://icet.materialsmodeling.org/moduleref/ensembles.html)

- [`class mchammer.ensembles.CanonicalEnsemble(structure, calculator, temperature, user_tag=None, boltzmann_constant=8.617330337217213e-05, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.calculator)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.ensemble_parameters)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.random_seed)
- [`run(number_of_trial_steps)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalEnsemble.write_data_container)
- [`class mchammer.ensembles.CanonicalAnnealing(structure, calculator, T_start, T_stop, n_steps, cooling_function='exponential', user_tag=None, boltzmann_constant=8.617330337217213e-05, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing)
- [`property T_start: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.T_start)
- [`property T_stop: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.T_stop)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.calculator)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.ensemble_parameters)
- [`property estimated_ground_state`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.estimated_ground_state)
- [`property estimated_ground_state_potential`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.estimated_ground_state_potential)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.get_random_sublattice_index)
- [`property n_steps: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.n_steps)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.random_seed)
- [`run()`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.CanonicalAnnealing.write_data_container)
- [`class mchammer.ensembles.SemiGrandCanonicalEnsemble(structure, calculator, temperature, chemical_potentials, boltzmann_constant=8.617330337217213e-05, user_tag=None, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.calculator)
- [`property chemical_potentials: dict[int, float]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.chemical_potentials)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.ensemble_parameters)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.random_seed)
- [`run(number_of_trial_steps)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SemiGrandCanonicalEnsemble.write_data_container)
- [`class mchammer.ensembles.SGCAnnealing(structure, calculator, T_start, T_stop, n_steps, chemical_potentials, cooling_function='exponential', boltzmann_constant=8.617330337217213e-05, user_tag=None, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.calculator)
- [`property chemical_potentials: dict[int, float]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.chemical_potentials)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.ensemble_parameters)
- [`property estimated_ground_state`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.estimated_ground_state)
- [`property estimated_ground_state_potential`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.estimated_ground_state_potential)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.random_seed)
- [`run()`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.SGCAnnealing.write_data_container)
- [`class mchammer.ensembles.VCSGCEnsemble(structure, calculator, temperature, phis, kappa, boltzmann_constant=8.617330337217213e-05, user_tag=None, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.calculator)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.ensemble_parameters)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.get_random_sublattice_index)
- [`property kappa: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.kappa)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.observers)
- [`property phis: dict[int, float]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.phis)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.random_seed)
- [`run(number_of_trial_steps)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.VCSGCEnsemble.write_data_container)
- [`class mchammer.ensembles.HybridEnsemble(structure, calculator, temperature, ensemble_specs, probabilities=None, boltzmann_constant=8.617330337217213e-05, user_tag=None, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, neighbor_sites_to_avoid=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.calculator)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.ensemble_parameters)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.observers)
- [`property probabilities: list[float]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.probabilities)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.random_seed)
- [`run(number_of_trial_steps)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.temperature)
- [`property trial_steps_per_ensemble: dict[str, int]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.trial_steps_per_ensemble)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.HybridEnsemble.write_data_container)
- [`class mchammer.ensembles.WangLandauEnsemble(structure, calculator, energy_spacing, energy_limit_left=None, energy_limit_right=None, trial_move='swap', fill_factor_limit=1e-06, flatness_check_interval=None, flatness_limit=0.8, window_search_penalty=2.0, user_tag=None, dc_filename=None, random_seed=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.attach_observer)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.calculator)
- [`property converged: bool | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.converged)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.data_container)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.ensemble_parameters)
- [`property fill_factor: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.fill_factor)
- [`property fill_factor_history: dict[int, float]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.fill_factor_history)
- [`property fill_factor_limit: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.fill_factor_limit)
- [`property flatness_check_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.flatness_check_interval)
- [`property flatness_limit: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.flatness_limit)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.random_seed)
- [`run(number_of_trial_steps)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.sublattices)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.WangLandauEnsemble.write_data_container)
- [`mchammer.ensembles.wang_landau_ensemble.get_bins_for_parallel_simulations(n_bins, energy_spacing, minimum_energy, maximum_energy, overlap=4, bin_size_exponent=1.0)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.wang_landau_ensemble.get_bins_for_parallel_simulations)
- [`class mchammer.ensembles.ThermodynamicIntegrationEnsemble(structure, calculator, temperature, n_steps, forward, user_tag=None, boltzmann_constant=8.617330337217213e-05, random_seed=None, dc_filename=None, data_container_write_period=600, ensemble_data_write_interval=None, trajectory_write_interval=None, sublattice_probabilities=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble)
- [`attach_observer(observer, tag=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.attach_observer)
- [`property boltzmann_constant: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.boltzmann_constant)
- [`property calculator: BaseCalculator`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.calculator)
- [`property data_container: BaseDataContainer`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.data_container)
- [`do_canonical_swap(sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.do_canonical_swap)
- [`do_sgc_flip(chemical_potentials, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.do_sgc_flip)
- [`do_thermodynamic_swap(sublattice_index, lambda_val, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.do_thermodynamic_swap)
- [`do_vcsgc_flip(phis, kappa, sublattice_index, allowed_species=None, allowed_sites=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.do_vcsgc_flip)
- [`property ensemble_parameters: dict`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.ensemble_parameters)
- [`get_random_sublattice_index(probability_distribution)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.get_random_sublattice_index)
- [`property observer_interval: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.observer_interval)
- [`property observers: dict[str, BaseObserver]`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.observers)
- [`property random_seed: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.random_seed)
- [`run()`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.run)
- [`property step: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.step)
- [`property structure: Atoms`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.structure)
- [`property sublattices: Sublattices`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.sublattices)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.temperature)
- [`update_occupations(sites, species)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.update_occupations)
- [`property user_tag: str | None`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.user_tag)
- [`write_data_container(outfile)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.ThermodynamicIntegrationEnsemble.write_data_container)
- [`class mchammer.ensembles.TargetClusterVectorAnnealing(structure, calculators, T_start=5.0, T_stop=0.001, random_seed=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing)
- [`property T_start: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.T_start)
- [`property T_stop: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.T_stop)
- [`property accepted_trials: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.accepted_trials)
- [`property best_score: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.best_score)
- [`property best_structure: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.best_structure)
- [`property current_score: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.current_score)
- [`generate_structure(number_of_trial_steps=None)`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.generate_structure)
- [`property n_steps: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.n_steps)
- [`property temperature: float`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.temperature)
- [`property total_trials: int`](https://icet.materialsmodeling.org/moduleref/ensembles.html#mchammer.ensembles.TargetClusterVectorAnnealing.total_trials)

### Observers

[官方模块说明](https://icet.materialsmodeling.org/moduleref/observers.html)

- [`class mchammer.observers.SiteOccupancyObserver(cluster_space, structure, sites, interval=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`class mchammer.observers.ShortRangeOrderObserver(cluster_space, structure, radius, pairs=None, interval=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property pairs: list[tuple[str, str]]`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property shells: list[dict]`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`class mchammer.observers.StructureFactorObserver(structure, q_points, symbol_pairs=None, form_factors=None, interval=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property form_factors: dict[str, float]`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property q_points: list[ndarray]`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`class mchammer.observers.ClusterCountObserver(cluster_space, structure, interval=None, orbit_indices=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_cluster_counts(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`class mchammer.observers.ClusterExpansionObserver(cluster_expansion, interval=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`class mchammer.observers.ConstituentStrainObserver(constituent_strain, interval=None)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`get_observable(structure)`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property interval: int`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property return_type: type`](https://icet.materialsmodeling.org/moduleref/observers.html)
- [`property tag: str`](https://icet.materialsmodeling.org/moduleref/observers.html)

### Structure containers

[官方模块说明](https://icet.materialsmodeling.org/moduleref/structure_container.html)

- [`class icet.StructureContainer(cluster_space)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`add_structure(structure, user_tag=None, properties=None, allow_duplicate=True, sanity_check=True)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`property available_properties: list[str]`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`property cluster_space: ClusterSpace`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`get_condition_number(structure_indices=None)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`get_fit_data(structure_indices=None, key='energy')`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`get_structure_indices(user_tag=None)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`static read(infile)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`to_dataframe()`](https://icet.materialsmodeling.org/moduleref/structure_container.html)
- [`write(outfile)`](https://icet.materialsmodeling.org/moduleref/structure_container.html)

### Tools

[官方模块说明](https://icet.materialsmodeling.org/moduleref/tools.html)

- [`icet.tools.map_structure_to_reference(structure, reference, *, inert_species=None, tol_positions=0.0001, tol_cell=0.25, symprec=1e-05, find_translation=True, trigger_levels=None, suppress_warnings=False, assume_no_cell_relaxation=False)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.map_structure_to_reference)
- [`class icet.tools.StructureMapping(drmax, dravg, transformation_matrix, translation, ambiguous_sites, warnings, strain_tensor, strain_tensor_eigenvalues, volumetric_strain, volume_dilation, isochoric_strain, isochoric_strain_eigenvalues, eigenstrain_norm, eigenstrain_rms, von_mises_strain, disorientation_angle, disorientation_axis, rotation)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping)
- [`drmax`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.drmax)
- [`dravg`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.dravg)
- [`transformation_matrix`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.transformation_matrix)
- [`translation`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.translation)
- [`ambiguous_sites`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.ambiguous_sites)
- [`warnings`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.warnings)
- [`strain_tensor`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.strain_tensor)
- [`strain_tensor_eigenvalues`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.strain_tensor_eigenvalues)
- [`volumetric_strain`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.volumetric_strain)
- [`volume_dilation`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.volume_dilation)
- [`isochoric_strain`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.isochoric_strain)
- [`isochoric_strain_eigenvalues`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.isochoric_strain_eigenvalues)
- [`eigenstrain_norm`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.eigenstrain_norm)
- [`eigenstrain_rms`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.eigenstrain_rms)
- [`von_mises_strain`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.von_mises_strain)
- [`disorientation_angle`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.disorientation_angle)
- [`disorientation_axis`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.disorientation_axis)
- [`rotation`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.StructureMapping.rotation)
- [`icet.tools.enumerate_structures(structure, sizes, chemical_symbols, concentration_restrictions=None, niggli_reduce=None, symprec=1e-05, position_tolerance=None)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.enumerate_structures)
- [`icet.tools.enumerate_supercells(structure, sizes, niggli_reduce=None, symprec=1e-05, position_tolerance=None)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.enumerate_supercells)
- [`icet.tools.training_set_generation.structure_selection_annealing(cluster_space, monte_carlo_structures, n_structures_to_add, n_steps, base_structures=None, cooling_start=5, cooling_stop=0.001, cooling_function='exponential', initial_indices=None)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.training_set_generation.structure_selection_annealing)
- [`icet.tools.structure_generation.generate_sqs(cluster_space, max_size, target_concentrations, include_smaller_cells=True, pbc=None, T_start=5.0, T_stop=0.001, n_steps=None, optimality_weight=1.0, random_seed=None, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_sqs)
- [`icet.tools.structure_generation.generate_sqs_by_enumeration(cluster_space, max_size, target_concentrations, include_smaller_cells=True, pbc=None, optimality_weight=1.0, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_sqs_by_enumeration)
- [`icet.tools.structure_generation.generate_sqs_from_supercells(cluster_space, supercells, target_concentrations, T_start=5.0, T_stop=0.001, n_steps=None, optimality_weight=1.0, random_seed=None, random_start=True, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_sqs_from_supercells)
- [`icet.tools.structure_generation.generate_target_structure(cluster_space, max_size, target_concentrations, target_cluster_vector, include_smaller_cells=True, pbc=None, T_start=5.0, T_stop=0.001, n_steps=None, optimality_weight=1.0, random_seed=None, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_target_structure)
- [`icet.tools.structure_generation.generate_target_structure_from_supercells(cluster_space, supercells, target_concentrations, target_cluster_vector, T_start=5.0, T_stop=0.001, n_steps=None, optimality_weight=1.0, random_seed=None, random_start=True, tol=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.generate_target_structure_from_supercells)
- [`icet.tools.structure_generation.occupy_structure_randomly(structure, cluster_space, target_concentrations, random_seed=None)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.structure_generation.occupy_structure_randomly)
- [`class icet.tools.ground_state_finder.GroundStateFinder(cluster_expansion, structure)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ground_state_finder.GroundStateFinder)
- [`property constraints: dict[str, highs_cons]`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ground_state_finder.GroundStateFinder.constraints)
- [`get_ground_state(species_count=None, max_seconds=inf, threads=0)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ground_state_finder.GroundStateFinder.get_ground_state)
- [`property model: Highs`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ground_state_finder.GroundStateFinder.model)
- [`property optimization_status: HighsModelStatus`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ground_state_finder.GroundStateFinder.optimization_status)
- [`class icet.tools.ConvexHull(concentrations, energies)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull)
- [`concentrations`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.concentrations)
- [`energies`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.energies)
- [`dimensions`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.dimensions)
- [`structures`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.structures)
- [`property concentration_labels: list[tuple[str, str]] | None`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.concentration_labels)
- [`extract_low_energy_structures(concentrations, energies, energy_tolerance)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.extract_low_energy_structures)
- [`classmethod from_sublattice_concentrations(concentrations, energies, site_fractions=None, cluster_space=None)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.from_sublattice_concentrations)
- [`get_chemical_potentials()`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_chemical_potentials)
- [`get_decomposition(target_concentrations)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_decomposition)
- [`get_energy_above_convex_hull(concentrations, energies)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_energy_above_convex_hull)
- [`get_energy_at_convex_hull(target_concentrations)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_energy_at_convex_hull)
- [`get_facet_gradients()`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_facet_gradients)
- [`get_facets()`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_facets)
- [`get_species_chemical_potentials(strict=True)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.get_species_chemical_potentials)
- [`is_on_convex_hull(concentrations, energies, energy_tolerance=1e-06)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConvexHull.is_on_convex_hull)
- [`icet.tools.get_sublattice_concentrations(structures, cluster_space)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.get_sublattice_concentrations)
- [`icet.tools.get_sublattice_site_fractions(structure, cluster_space)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.get_sublattice_site_fractions)
- [`class icet.tools.constraints.Constraints(n_params)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.constraints.Constraints)
- [`add_constraint(M)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.constraints.Constraints.add_constraint)
- [`inverse_transform(A)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.constraints.Constraints.inverse_transform)
- [`transform(A)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.constraints.Constraints.transform)
- [`icet.tools.constraints.get_mixing_energy_constraints(cluster_space)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.constraints.get_mixing_energy_constraints)
- [`class icet.tools.ConstituentStrain(supercell, primitive_structure, chemical_symbols, concentration_symbol, strain_energy_function, k_to_parameter_function=None, damping=1.0, tol=1e-06)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain)
- [`accept_change()`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain.accept_change)
- [`get_concentration(occupations)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain.get_concentration)
- [`get_constituent_strain(occupations, update_structure_factors=True)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain.get_constituent_strain)
- [`get_constituent_strain_change(occupations, atom_index)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.ConstituentStrain.get_constituent_strain_change)
- [`icet.tools.get_primitive_structure(structure, no_idealize=True, to_primitive=True, symprec=1e-05)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.get_primitive_structure)
- [`icet.tools.get_wyckoff_sites(structure, map_occupations=None, symprec=1e-05, include_representative_atom_index=False)`](https://icet.materialsmodeling.org/moduleref/tools.html#icet.tools.get_wyckoff_sites)


## 附录 B 官方页面覆盖清单

本清单逐页对应已读取的文档与教程章节。主站 47 个页面；扩展站读取了 17 个 URL，其中根地址与 `index.html` 是同一首页，因此合计为 63 个去重文档页面。站点首页、导航页和索引页单独计入清单；它们用于检查内容覆盖而不是独立教程章节。

|官方页面|对应章节|
|---|---|
|[icet — A pythonic approach to cluster expansions](https://icet.materialsmodeling.org/index.html)|1，阅读路线|
|[Advanced topics](https://icet.materialsmodeling.org/advanced_topics/index.html)|8 至 28|
|[Bibliography](https://icet.materialsmodeling.org/backmatter/bibliography.html)|31|
|[Glossary](https://icet.materialsmodeling.org/backmatter/glossary.html)|31|
|[Backmatter](https://icet.materialsmodeling.org/backmatter/index.html)|30，31|
|[Credits](https://icet.materialsmodeling.org/credits.html)|31|
|[Get started](https://icet.materialsmodeling.org/get_started/index.html)|1，10，29|
|[Reference](https://icet.materialsmodeling.org/moduleref/index.html)|附录 A|
|[Using a calculator directly](https://icet.materialsmodeling.org/advanced_topics/calculator_interface.html)|27|
|[Cluster vectors](https://icet.materialsmodeling.org/advanced_topics/cluster_vectors.html)|3，27|
|[Constituent strain calculations](https://icet.materialsmodeling.org/advanced_topics/constituent_strain.html)|26|
|[Customizing cluster spaces](https://icet.materialsmodeling.org/advanced_topics/customizing_cluster_spaces.html)|13|
|[Data container](https://icet.materialsmodeling.org/advanced_topics/data_container.html)|22|
|[Thermodynamic integration and temperature-integration simulations](https://icet.materialsmodeling.org/advanced_topics/free_energy_canonical_ensemble.html)|25|
|[Hybrid ensembles](https://icet.materialsmodeling.org/advanced_topics/hybrid_ensembles.html)|20|
|[Mapping structures](https://icet.materialsmodeling.org/advanced_topics/mapping_structures.html)|6|
|[Parallel Monte Carlo simulations](https://icet.materialsmodeling.org/advanced_topics/parallel_monte_carlo.html)|28|
|[Special quasirandom structures](https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html)|15|
|[Structure enumeration](https://icet.materialsmodeling.org/advanced_topics/structure_enumeration.html)|8|
|[Bayesian cluster expansions](https://icet.materialsmodeling.org/advanced_topics/training_bayesian_cluster_expansions.html)|13|
|[Selecting cutoffs](https://icet.materialsmodeling.org/advanced_topics/training_cutoffs_selection.html)|11|
|[EnsembleOptimizer](https://icet.materialsmodeling.org/advanced_topics/training_ensemble_of_models.html)|12|
|[Hyper-parameter scans](https://icet.materialsmodeling.org/advanced_topics/training_hyper_parameter_scans.html)|11|
|[Training of cluster expansions](https://icet.materialsmodeling.org/advanced_topics/training_of_ce.html)|10 至 12|
|[Generating training structures](https://icet.materialsmodeling.org/advanced_topics/training_set_generation.html)|9|
|[Training with constraints and weights](https://icet.materialsmodeling.org/advanced_topics/training_with_constraints.html)|13|
|[Wang-Landau simulations](https://icet.materialsmodeling.org/advanced_topics/wang_landau_simulations.html)|24|
|[FAQ](https://icet.materialsmodeling.org/backmatter/faq.html)|30|
|[Index](https://icet.materialsmodeling.org/backmatter/genindex.html)|附录 A|
|[Analyzing ECIs](https://icet.materialsmodeling.org/get_started/analyze_ecis.html)|3，10|
|[Analyzing Monte Carlo simulations](https://icet.materialsmodeling.org/get_started/analyze_monte_carlo.html)|22，23|
|[Cluster expansions](https://icet.materialsmodeling.org/get_started/cluster_expansions.html)|3|
|[Comparison with target data](https://icet.materialsmodeling.org/get_started/compare_to_target_data.html)|10|
|[Constructing a cluster expansion](https://icet.materialsmodeling.org/get_started/construct_cluster_expansion.html)|4，5，7，10|
|[Enumerating structures](https://icet.materialsmodeling.org/get_started/enumerate_structures.html)|8，14|
|[Installation](https://icet.materialsmodeling.org/get_started/installation.html)|2|
|[Monte Carlo simulations](https://icet.materialsmodeling.org/get_started/run_monte_carlo.html)|16，18，19|
|[Workflow](https://icet.materialsmodeling.org/get_started/workflow.html)|1，3，7，16|
|[Calculators](https://icet.materialsmodeling.org/moduleref/calculators.html)|16，26，27，附录 A|
|[Cluster expansions](https://icet.materialsmodeling.org/moduleref/cluster_expansion.html)|10，14，附录 A|
|[Cluster space](https://icet.materialsmodeling.org/moduleref/cluster_space.html)|3，5，27，附录 A|
|[Core components](https://icet.materialsmodeling.org/moduleref/core.html)|5，20，26，27，附录 A|
|[Data containers](https://icet.materialsmodeling.org/moduleref/data_containers.html)|22 至 25，附录 A|
|[Ensembles](https://icet.materialsmodeling.org/moduleref/ensembles.html)|17 至 20，24，25，附录 A|
|[Observers](https://icet.materialsmodeling.org/moduleref/observers.html)|21，附录 A|
|[Structure containers](https://icet.materialsmodeling.org/moduleref/structure_container.html)|7，附录 A|
|[Tools](https://icet.materialsmodeling.org/moduleref/tools.html)|6 至 9，13 至 15，26，附录 A|
|[扩展教程 Active learning](https://ce-tutorials.materialsmodeling.org/part-1/active-learning.html)|12|
|[扩展教程 Condition number selection of structures](https://ce-tutorials.materialsmodeling.org/part-1/condition-number-selection.html)|9|
|[扩展教程 Cutoff selection](https://ce-tutorials.materialsmodeling.org/part-1/cutoff-selection.html)|11|
|[扩展教程 Structure enumeration](https://ce-tutorials.materialsmodeling.org/part-1/enumeration.html)|8|
|[扩展教程 Hyper-parameter scans](https://ce-tutorials.materialsmodeling.org/part-1/hyper-parameter-scans.html)|11|
|[扩展教程 Part 1: CE construction](https://ce-tutorials.materialsmodeling.org/part-1/index.html)|8 至 13|
|[扩展教程 Orthogonal structure generation](https://ce-tutorials.materialsmodeling.org/part-1/orthogonal-structures.html)|9|
|[扩展教程 Random structure generation](https://ce-tutorials.materialsmodeling.org/part-1/random-structures.html)|9|
|[扩展教程 Regression methods](https://ce-tutorials.materialsmodeling.org/part-1/regression-methods.html)|10，11|
|[扩展教程 Part 2: Low-symmetry systems](https://ce-tutorials.materialsmodeling.org/part-2/index.html)|13|
|[扩展教程 Regular CE](https://ce-tutorials.materialsmodeling.org/part-2/low-symmetry-ce.html)|13|
|[扩展教程 Canonical ensemble](https://ce-tutorials.materialsmodeling.org/part-3/canonical-simulations.html)|17，21，22|
|[扩展教程 Part 3: MC sampling](https://ce-tutorials.materialsmodeling.org/part-3/index.html)|17 至 23|
|[扩展教程 SGC and VCSGC ensembles](https://ce-tutorials.materialsmodeling.org/part-3/sgc-vcsgc-simulations.html)|18，19，23|
|[扩展教程 The VCSCG variance constraint $\bar{\kappa}$](https://ce-tutorials.materialsmodeling.org/part-3/vcsgc-variance-constraint.html)|19|
|[扩展教程 Cluster expansion tutorials](https://ce-tutorials.materialsmodeling.org/index.html)|1，阅读路线|

## 附录 C 操作前后核对表

**建模前：**母晶格与允许物种明确；二维真空和周期约定正确；空位位点完整；参考能和标签归一化记录清楚。

**拟合前：**映射诊断通过；数据矩阵无无意重复和严重冗余；训练与验证划分符合研究目标；物理约束用正确参数坐标表达。

**MC 前：**模型与超胞匹配；位点占位合法；自相互作用已处理；能量缩放与温度单位一致；系综符合组成边界条件；初态满足排斥规则；输出目录、种子和观测间隔记录清楚。

**报告前：**实际平衡期已丢弃；轨迹长度足以得到可信统计；有限尺寸和模型敏感性已检查；自由能参考一致；软件与方法引用完整；教学参数未误写成真实材料结果。
