# MTP 势函数训练教程：从 DFT 数据准备到主动学习与 LAMMPS 验证

编写日期：2026-09-11  
适用范围：Shapeev 团队的 MLIP-2；补充 MLIP-3 的迁移边界和高熵合金实践。  
运行环境：下文命令使用 Linux/Bash，适用于 Linux 工作站、集群或 WSL；不是 Windows PowerShell 命令。  
核验说明：本文核对了开发者官方教程 Wiki、用户手册及源码；本次编写没有执行 DFT、编译 MLIP 或实际训练势函数。文中官方答案用于复核安装和工作流，不是本文计算结果。

## 阅读路线与来源

第一次接触 MTP，建议先读第 1—5 节理解数据和模型，再按第 6 节跑通官方 Mo 示例；随后完成第 7—9 节的性质验证和主动学习，最后用第 10 节设计自己的高熵合金训练任务。

本文主要依据以下开发者维护的材料，使用中文重新组织步骤，并加入数据质量、复现和高熵合金应用方面的实践说明：

- [MLIP 官方网站](https://mlip.skoltech.ru/)：本次读取主页超时，因此以开发者 GitLab 教程、源码和官方论文为核验依据。
- [MLIP-2 官方教程入口](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/home)及[示例文件仓库](https://gitlab.com/ashapeev/mlip-2-tutorials)。Tutorial 1 是 bcc-Mo 的训练、弛豫、能量—体积曲线与弹性常数；Tutorial 2 是 bcc-Nb 的分子动力学主动学习。
- [MLIP-2 用户手册](https://gitlab.com/ashapeev/mlip-2-paper-supp-info/-/blob/master/manual.pdf)：文件格式、命令和配置项的定义。
- [MLIP-2 源码](https://gitlab.com/ashapeev/mlip-2)与 [LAMMPS 接口](https://gitlab.com/ashapeev/interface-lammps-mlip-2)。
- [MLIP 方法论文](https://arxiv.org/abs/2007.08555)与 [MLIP-3 方法论文](https://arxiv.org/abs/2304.13144)。

本次核验版本：教程 Wiki 提交 `3357b75fb1c904005d92315b7f14a8485ab83f84`；MLIP-2 源码提交 `fc3b0721b4d0eb1befcc4eb48269dfc6812f37d4`。若未来教程页面步骤号调整，可按步骤标题及仓库目录定位。

> 三点需要先明确：本文所说的 MLIP 是实现 MTP 的软件包，不是其他同名 Python/JAX 项目；官方 Tutorial 2 用 EAM 作为教学参考模型，研究级 DFT 势需要替换这个标注环节；训练集误差低，只说明模型拟合了这些数据，不能单独证明目标材料性质可靠。

## 1. 训练 MTP 实际上在做什么

### 1.1 输入、输出和模型

一个训练样本是一个原子构型，包含晶胞、原子位置、元素种类，以及参考计算提供的能量、力和可选的应力信息。参考计算通常是 DFT。

MTP 对总能量进行局域分解：

$$
E(\mathcal R)=\sum_{i=1}^{N}V(\mathcal N_i),\qquad
V(\mathcal N_i)=\sum_\alpha \xi_\alpha B_\alpha(\mathcal N_i).
$$

其中 $\mathcal N_i$ 是截断半径内的原子环境，$B_\alpha$ 是由矩张量收缩得到的旋转不变量。矩张量可以写成：

$$
M_{\mu,\nu}(\mathcal N_i)
=\sum_{j\in\mathcal N_i}
f_\mu(r_{ij},z_i,z_j)\,\mathbf r_{ij}^{\otimes\nu}.
$$

对固定的径向函数，能量对外层系数 $\xi_\alpha$ 是线性的；实际多元素 MTP 还会优化径向系数，因此完整训练通常是非线性优化。不能把多元素 MTP 简化为“一次线性最小二乘就完成全部训练”。模型参数和算法背景见 [MLIP 方法论文](https://arxiv.org/abs/2007.08555)及[官方源码](https://gitlab.com/ashapeev/mlip-2)。

力由同一个能量函数求导得到：

$$
\mathbf F_i=-\frac{\partial E}{\partial\mathbf r_i}.
$$

这样能量和力共享一个势能面。DFT 并不需要提供唯一的“每原子能量标签”；训练使用构型总能量，局域原子能量是模型内部的分解。

### 1.2 一个完整训练项目的文件

| 文件 | 保存什么 | 后续用途 |
|---|---|---|
| `train.cfg` | 训练构型及参考标签 | 优化参数 |
| `valid.cfg` | 独立验证构型及参考标签 | 选择 level、权重和训练结果 |
| `test.cfg` | 最终保留测试集 | 最终评估，不用于反复调参 |
| `init.mtp` | 未拟合的完整基函数模板，或明确指定的已有势 | 训练起点 |
| `pot.mtp` | 拟合后的势函数 | 能量、力计算和 MD |
| `state.als` | 主动学习状态 | 计算外推等级与选择构型 |
| `mlip.ini` | 推理、弛豫或主动学习配置 | 被 MLIP/LAMMPS 接口读取 |
| `preselected.cfg` | 探索过程中超出选择阈值的候选构型 | 主动学习去冗余选择 |
| `selected.cfg` | 进一步选中的构型 | 送往参考计算 |
| `labeled.cfg` | 完成参考标注并检查后的新构型 | 加入下一轮训练 |

`selected.cfg` 中即使已经有能量或力，也可能来自当前 MTP；这些数值不能冒充 DFT 标签。

## 2. 安装与版本确认

### 2.1 获取代码和教程

先准备 C/C++ 编译器、make；MPI 版本还需要相应 MPI 工具链。BLAS 库选择和编译器要求以仓库的 `INSTALL.md` 为准。以下命令在一个新的学习目录中执行：

```bash
mkdir mtp-learning
cd mtp-learning
git clone https://gitlab.com/ashapeev/mlip-2.git
git clone https://gitlab.com/ashapeev/mlip-2-tutorials.git
cd mlip-2
./configure
make mlp
```

官方安装文档说明，默认配置编译 MPI 版本。若只想先使用串行版本，可在独立源码目录中选择：

```bash
./configure --no-mpi
make mlp
```

成功后可执行文件位于 `bin/mlp`。访问权限或许可证条件如有变化，以官方仓库当时的要求为准，不以旧教程中的历史说明代替当前条件。[官方安装说明](https://gitlab.com/ashapeev/mlip-2/-/blob/master/INSTALL.md)

### 2.2 记录实际使用的版本

在 MLIP 源码目录中执行：

```bash
git rev-parse HEAD
./bin/mlp list
./bin/mlp help train
./bin/mlp help convert-cfg
./bin/mlp help calc-grade
./bin/mlp help select-add
```

先读本机帮助，再套用本文命令。路径中最好避免空格。对于较长的研究流程，可以在 Bash 中设置一个绝对路径变量，例如：

```bash
export MTP_BIN=/absolute/path/to/mlip-2/bin/mlp
"$MTP_BIN" list
```

上面的路径必须替换为真实路径。下文官方示例使用 `../mlp`，自建项目使用 `$MTP_BIN`，两种写法指向的是同一个程序。

### 2.3 LAMMPS 需要匹配的接口

普通 LAMMPS 安装不一定认识 `pair_style mlip`。应按[官方接口 README](https://gitlab.com/ashapeev/interface-lammps-mlip-2/-/blob/master/README.md)编译配套版本，并记录 MLIP、接口和 LAMMPS 的提交号。不要假定旧接口能直接编译到任意最新 LAMMPS。

官方 Tutorial 2 对其对应版本建议使用串行 LAMMPS 执行主动选择，训练可以使用 MPI。这个限制是历史版本的具体条件，不能据此断言所有 MLIP 版本都不支持并行主动学习。[Tutorial 2：准备可执行文件](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Step-0)

## 3. 从零准备 DFT 数据集

以下是将官方工作流用于实际研究的数据准备建议。具体体系没有确定时，不存在一组普遍适用的 DFT 截断能、k 点、温度或数据量。

### 3.1 先写清势函数的适用范围

至少明确元素、组成区间、晶体结构、温度、压力、应变范围和目标性质。例如：

> 研究 CoCrFeMnNi 等原子比合金，在目标温度区间内的 fcc 固溶体、局部 hcp 堆垛和空位；目标是弹性响应与层错相关性质。

这个定义会指导训练集。若之后增加液态、裂纹或不同组成，就需要重新评估数据覆盖。

| 应用 | 建议重点采样 | 应保留的独立验证 |
|---|---|---|
| 平衡结构、弹性 | 体积变化、小应变、热扰动 | EOS、弹性常数、力 |
| 化学短程有序 | 随机排列、MC 交换构型、有序与偏聚环境 | 交换能差、SRO 参数随温度变化 |
| 层错和塑性 | fcc/hcp、剪切路径、层错附近局域化学环境 | 广义层错能曲线、关键能垒 |
| 扩散 | 空位、局部迁移路径、鞍点附近构型 | 迁移能垒、扩散趋势 |
| 熔化、凝固 | 高温晶体、液体、必要的界面环境 | 液态结构、固液行为 |
| 成分筛选 | 多个组成及相关子系统 | 留出组成的性质和能量差 |

纯元素、二元和三元结构有助于扩大适用范围，但并非所有固定组成 MTP 都必须穷尽全部低阶子系统。应根据拟应用的组成空间与性质做数据预算，并用留出测试检查是否充分。

### 3.2 构造候选结构

对高熵合金，应分别考虑几何变化与化学变化。几何采样可以是体积缩放、剪切、随机位移和热运动；化学采样可以是多个随机种子、多个 SQS、不同成分及元素交换轨迹。

SQS 近似指定随机合金的相关函数。一个 SQS 不能代表所有局域化学环境，也不等价于有限温度平衡短程有序。若研究 SRO，仅延长固定化学占位的普通 MD 轨迹通常不足以充分采样化学排列。[CrCoNi 短程有序研究](https://www.nature.com/articles/s41524-025-01722-2)

建议先用少量代表结构建立初始数据，之后通过主动学习补充。结构生成阶段可以借助已有经验势，但最终作为 DFT 训练标签的数据，仍需由规定的 DFT 方法重新计算。

### 3.3 固定 DFT 标注协议

建立一份计算说明，至少记录：

1. 泛函、赝势/PAW 数据集及版本、截断能。
2. k 点密度与网格设置规则，电子展宽方案。
3. 电子收敛阈值、是否自旋极化、初始磁矩和磁性状态处理。
4. 能量字段的选择，以及与力、应力的对应关系。
5. 是否固定晶胞、是否移动原子、是否采用对称性。

应以少量代表构型做截断能和 k 点收敛测试，使标签数值噪声低于研究所需的精度。大量未弛豫的热快照需要非零力标签；若全部先弛豫到极小值，再拿这些结构训练，模型将缺少离开极小值后的力信息。

含 Fe、Co、Cr、Mn 的体系还应注意：普通 MTP 没有显式自旋输入。同一个原子构型如果混入多个磁性分支的能量和力，会使监督目标不再单值。需要选择一致、可复现的磁性标注方案；普通原子势不会自动学出磁性自由度。

### 3.4 AIMD 采样和单点回标

从多个独立初始结构开展 AIMD，覆盖目标温度和体积。去掉未平衡阶段，按相关时间抽样；抽样间隔应由能量、力或结构量的相关性判断，不能固定认为“每 10 帧”或“每 100 帧”就独立。

为了统一精度，可以先用适合采样成本的设置生成轨迹，再对选中快照以统一的高精度设置做静态单点。回标时保留快照原子位置和晶胞，只求电子自洽；不能先弛豫掉主动学习想补充的异常环境。

每条数据应追溯到母结构、轨迹编号、温度、时间步、化学组成、参考计算目录及收敛状态。数据清洗应排除未收敛和格式错误；高能、大力构型如果属于真实目标过程，不应仅因“数值大”就删掉。

## 4. 转换为 CFG：最容易影响整个项目的一步

### 4.1 CFG 的关键约定

MLIP-2 的构型以 `BEGIN_CFG` 开始、`END_CFG` 结束；一个文件可串联多个完整构型。`Size` 给出原子数，`Supercell` 给出晶格向量，`AtomData` 指明原子字段。官方手册规定，`Energy` 是总能量，`type` 从 0 开始；`PlusStress` 是乘以体积的应力，单位为 eV，压缩为正。[官方 CFG 格式](https://gitlab.com/ashapeev/mlip-2-paper-supp-info/-/blob/master/manual.pdf)

下面只是一个无标签的 bcc 单元素几何示例，用来理解格式；没有能量和力，不能当作完整监督训练样本：

```text
BEGIN_CFG
 Size
    2
 Supercell
    3.2 0.0 0.0
    0.0 3.2 0.0
    0.0 0.0 3.2
 AtomData: id type cartes_x cartes_y cartes_z
    1 0 0.0 0.0 0.0
    2 0 1.6 1.6 1.6
END_CFG
```

带标签的数据还需要加入 `Energy` 及 `fx fy fz`；若拟合应力则加入 `PlusStress: xx yy zz yz xz xy`。不要为缺失标签填零，零代表一个真实的监督数值。

### 4.2 VASP OUTCAR 转换

先对一个确认正常的 OUTCAR 做转换：

```bash
"$MTP_BIN" convert-cfg OUTCAR trajectory.cfg --input-format=vasp-outcar
```

若只需要弛豫过程最后一个构型，可以使用：

```bash
"$MTP_BIN" convert-cfg OUTCAR final.cfg --input-format=vasp-outcar --last
```

`--last` 不适用于希望保留多个热快照的 AIMD 数据。转换工具不会自动替你判断电子收敛、数据独立性或训练适用性。核验版本的帮助中只列出了历史上测试的 VASP 版本，因此使用新版本 OUTCAR 时应人工比对至少一个样本。[转换命令源码](https://gitlab.com/ashapeev/mlip-2/-/blob/fc3b0721b4d0eb1befcc4eb48269dfc6812f37d4/src/mlp/mlp_commands.cpp)

核对总原子数、晶胞、原子顺序、总能量、选定原子的力和全部应力分量。不同 VASP 能量字段，例如自由能与去除熵项的能量，不应在数据集中随意混用；核验版本的 OUTCAR 读取代码读取 `TOTEN` 相关字段，若使用其他解析器要检查其能量约定。[OUTCAR 读取实现](https://gitlab.com/ashapeev/mlip-2/-/blob/fc3b0721b4d0eb1befcc4eb48269dfc6812f37d4/src/configuration.cpp)

### 4.3 高熵合金必须固定全局元素编号

例如始终采用：

| 元素 | CFG/MTP 的 type | 按顺序映射的 LAMMPS atom type |
|---|---:|---:|
| Co | 0 | 1 |
| Cr | 1 | 2 |
| Fe | 2 | 3 |
| Mn | 3 | 4 |
| Ni | 4 | 5 |

危险情况是：纯 Ni 的 OUTCAR 转换后可能把 Ni 编成局部 `type=0`，但五元训练集中的 `type=0` 是 Co。直接追加文件不会让程序自动理解两者差别。官方仓库也记录过这种多文件转换问题。[元素编号问题记录](https://gitlab.com/ashapeev/mlip-2/-/work_items/37)

正确做法是先单独转换每个计算，依据该计算的真实元素顺序建立“局部编号到全局编号”的映射，再合并。若重新排序原子，位置、元素、力和原子编号必须一起重排。不能只替换元素标签，也不要假设所有二元/三元 OUTCAR 的第一个元素都相同。

### 4.4 应力单位与符号的核查

设 CFG 中的量为 $W_{\alpha\beta}$。对于拉伸为正、单位为 $\mathrm{eV}/\mathrm{\AA}^3$ 的应力 $\sigma^{\mathrm{tension}}$，它与本文采用的压缩为正约定满足：

$$
W_{\alpha\beta}=-V\sigma^{\mathrm{tension}}_{\alpha\beta}.
$$

若参考输出已经是压缩为正的应力，就不应再额外改变符号。单位换算为：

$$
1\ \mathrm{eV}/\mathrm{\AA}^3
=160.21766208\ \mathrm{GPa},\qquad
1\ \mathrm{kbar}=0.1\ \mathrm{GPa}.
$$

因此把 GPa 数字直接写入 `PlusStress` 会错一个体积和单位因子。六分量还必须按字段名对应，不能把 `xy`、`yz`、`xz` 机械地按某个解析器的默认顺序粘贴。

一个有效的小检查是构造轻微均匀压缩，确认压力符号和能量—体积导数一致。还应在同一构型上对比 MLIP 的单点结果与 LAMMPS `run 0`，排除接口映射错误。

### 4.5 清洗、去重和分组划分

建议在切分前完成格式检查和近重复检查，但只能用训练部分确定数据驱动的拟合预处理参数。保留每类结构的构型数、原子数、最小距离、每原子能量和力分布。

可从 80%/10%/10% 的训练/验证/测试比例开始。这是工作建议，不是 MLIP 的硬性要求。更重要的是按母结构、轨迹或化学排列分组划分，避免同一 AIMD 轨迹的邻帧分别进入训练和测试。主动学习产生的新训练构型也不能与最终测试集近重复。

最终形成 `train.cfg`、`valid.cfg`、`test.cfg` 和记录划分来源的清单。不同组成、温度或缺陷的外推能力，应设置额外的留出测试，单独报告。

## 5. 设置 MTP 模板与训练目标

### 5.1 选择完整的未训练模板

从 `untrained_mtps/` 复制如 `08.mtp`、`12.mtp`、`16.mtp` 的完整文件。level 决定基函数集合，不是神经网络层数；不能仅修改文件名或随便增加一个 `level=20` 字段来提高复杂度。

高熵合金可以先比较较低、中等 level，在同一训练和验证集上评估精度、速度和稳定性，再决定是否增加复杂度。数据不够时，增大 level 不会自动补足未采样的化学或几何环境。

对于五元素模型，以下字段说明应该怎样设置；它们只是完整模板的一部分，不能单独保存为可训练的势文件：

```text
species_count = 5
radial_basis_type = RBChebyshev
    min_dist = <根据训练集确定>
    max_dist = 5.0
    radial_basis_size = 8
```

尖括号必须替换为数值。`max_dist=5.0` 和径向基大小 8 可作为试验起点，但不是对所有高熵合金验证过的最优设置。`radial_basis_size=8` 是基函数数量，没有长度单位。其余 `alpha_index_*` 等字段定义具体基函数组合，应保留模板内容。[官方模板](https://gitlab.com/ashapeev/mlip-2/-/tree/master/untrained_mtps)

### 5.2 `min_dist` 应怎样选

先运行：

```bash
"$MTP_BIN" mindist train.cfg
```

检查最短距离对应的构型是否合理，然后把 `min_dist` 设置为略低于有效训练集最短距离的数值。官方 Mo 示例是 `1.95588`，设置为 `1.9`；Nb 示例是 `2.85455`，设置为 `2.8`。这些数值仅属于各自数据集。[Mo 距离检查](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-1)、[Nb 初始训练](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Step-1)

`min_dist` 是径向基的相关设置，不是阻止原子靠近的硬墙。模拟走到训练集未覆盖的更短距离时，仍可能发生不可靠外推。使用 `--update-mindist` 也不等于获得了准确的短程排斥。

### 5.3 损失函数与权重

可用下面的示意表达式理解训练：

$$
L(\theta)=\sum_s\left[
w_E a_s(\Delta E_s)^2
w_F b_s\sum_i\|\Delta\mathbf F_{si}\|^2
w_W c_s\sum_{\alpha\beta}(\Delta W_{s,\alpha\beta})^2
\right].
$$

$a_s,b_s,c_s$ 表示构型大小、归一化等因素。这个式子用于解释能量、力、应力如何共同约束模型，不代表所有版本的逐项实现完全相同；实际权重还受 `--weighting` 等设置影响。

本次核验的 MLIP-2 命令帮助给出能量、力、应力默认权重分别为 `1`、`0.01`、`0.001`，最大迭代数为 `1000`，`--bfgs-conv-tol` 默认 `1e-3`，并支持 `vibrations`、`molecules`、`structures` 等构型大小加权方式。[训练命令实现](https://gitlab.com/ashapeev/mlip-2/-/blob/fc3b0721b4d0eb1befcc4eb48269dfc6812f37d4/src/mlp/mlp_commands.cpp)

“力权重 0.01”不能解释为力只占总损失的 1%；每个构型有多个力分量，各量的单位和归一化也不同。应根据验证误差和目标性质调权重。若没有应力标签，先把应力权重设为零，并明确该模型未通过应力训练得到约束。

## 6. 跑通官方 Tutorial 1：bcc-Mo 被动训练

### 6.1 准备目录

假设已经在同一个 `mtp-learning/` 目录下克隆两个官方仓库，并完成编译。回到该目录，执行：

```bash
cp mlip-2/bin/mlp mlip-2-tutorials/tutorial-1/mlp
cd mlip-2-tutorials/tutorial-1
```

官方 `tutorial-1/` 存放输入文件，`tutorial-1-answers/` 提供参考答案。[Tutorial 1 概览](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Home)

### 6.2 检查训练集和模板

从 `tutorial-1/` 进入：

```bash
cd 1.training_test_set
../mlp mindist train.cfg
```

对官方未修改的数据，预期 `Global mindist` 为 `1.95588`。然后进入训练目录，并复制模板：

```bash
cd ../2.training_validation_error
cp ../../../mlip-2/untrained_mtps/16.mtp init.mtp
```

上述相对路径对应本文前面两个仓库并列的布局。编辑完整的 `init.mtp`，设 `min_dist=1.9`，保留 `max_dist=5` 和默认径向基大小。检查该目录已有 `train.cfg` 和 `test.cfg`。[官方拟合步骤](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-2)

### 6.3 执行拟合

官方核心命令是：

```bash
../mlp train init.mtp train.cfg \
  --trained-pot-name=pot.mtp \
  --valid-cfgs=test.cfg
```

MPI 编译的程序可以在已申请到的计算资源内训练，例如：

```bash
mpirun -n 8 ../mlp train init.mtp train.cfg \
  --trained-pot-name=pot.mtp \
  --valid-cfgs=test.cfg > train.log 2>&1
```

这里的 8 只是资源示例，应与调度系统实际分配匹配。不要对串行二进制启动多个 MPI 副本，否则可能造成多个进程写同一输出文件。

官方演示将 `test.cfg` 传给 `--valid-cfgs`。在正式研究中，建议把调参集命名为 `valid.cfg`，最终测试集另行保留；参数名不会自动保证数据用途正确。

### 6.4 检查训练输出

应看到损失与能量、力、应力误差信息，并得到带拟合系数的 `pot.mtp`。训练结束后独立计算一次误差：

```bash
../mlp calc-errors pot.mtp train.cfg
../mlp calc-errors pot.mtp test.cfg
```

重点判断训练是否正常退出、误差是否存在极端离群值、验证误差是否明显高于训练误差。达到最大迭代数并不证明完全收敛，损失停止下降也不证明目标性质已达到精度要求。

### 6.5 自建研究项目的训练命令模板

已经准备好完整 `init.mtp`、`train.cfg`、`valid.cfg` 后，可以显式记录参数：

```bash
"$MTP_BIN" train init.mtp train.cfg \
  --trained-pot-name=pot.mtp \
  --curr-pot-name=checkpoint.mtp \
  --valid-cfgs=valid.cfg \
  --energy-weight=1 \
  --force-weight=0.01 \
  --stress-weight=0.001 \
  --weighting=vibrations \
  --max-iter=1000 \
  --bfgs-conv-tol=1e-3 > train.log 2>&1
```

这一块是基于已核验命令选项组织的研究起点，不是任何高熵合金的最优配方。`--curr-pot-name` 保存过程中的势，`--trained-pot-name` 保存最终势。按验证结果调整超参数时，应保留每次日志和对应输入。

多个随机初始化应从完整未训练模板出发，在不同目录中训练。对已经有系数的势继续拟合属于热启动，不能把它当作独立随机初始化；`--init-params=same` 也不能简单理解成“沿用已有模型的随机种子”。

## 7. 按官方流程做性质验证

### 7.1 弛豫结构

在官方 `tutorial-1/` 布局中，把训练后的 `pot.mtp` 复制到 `3.relaxation/`，使用该目录提供的 `init.cfg` 和 `relax.ini`：

```bash
../mlp relax relax.ini \
  --cfg-filename=init.cfg \
  --save-relaxed=relaxed.cfg \
  --force-tolerance=1e-6 \
  --stress-tolerance=1e-5
```

这是官方示例的收敛设置，不宜未经检查地用作所有体系的统一标准。检查得到的晶胞、能量和残余力，而不只是看文件是否生成。[官方弛豫步骤](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-3)

对高熵合金，应在相同边界条件下比较 DFT 和 MTP：例如都固定晶胞只弛豫原子，或都放开晶胞。比较“DFT 固定体积结果”和“MTP 完全弛豫结果”会混入额外差异。

### 7.2 能量—体积曲线

把 `pot.mtp` 和 `relaxed.cfg` 放入 `4.energy_volume_curve/`，执行：

```bash
../mlp compress-extend relaxed.cfg deformed.cfg
../mlp calc-efs pot.mtp deformed.cfg deformed_efs.cfg
python cfg-to-energy-volume.py deformed_efs.cfg E_V.txt
```

`compress-extend` 和提取脚本的用法依据当前官方示例，若旧版本不支持前者，就使用示例已有的变形构型或自行生成等比例缩放结构。将 MTP 曲线与目录中的 DFT 参考数据比较。[官方 EOS 步骤](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-4)

观察平衡体积、曲率和压缩/膨胀两侧的趋势。体积模量与最低点附近的曲率有关：

$$
B=V_0\left.\frac{\partial^2E}{\partial V^2}\right|_{V_0}.
$$

在自建数据上，可用多个应变幅度重复检查曲线是否平滑，排除异常极小值。若采用每原子能量，体积也应使用每原子体积以保持定义一致。

### 7.3 弹性常数

官方示例通过变形构型的应力有限差分与 LAMMPS 两种方式验证。定义拉伸为正的应力后：

$$
C_{ij}\approx\frac{\sigma_i(+\delta\epsilon_j)-\sigma_i(-\delta\epsilon_j)}{2\delta\epsilon_j}.
$$

计算前要把 CFG 中的 `PlusStress` 转成实际应力；剪切分量还应明确采用工程剪切还是张量剪切，避免差一个因子 2。对高熵合金的有限随机超胞，不应强行假定其完整弹性张量严格具有立方对称性；可以先计算完整张量，再做有说明的对称投影或构型平均。

官方答案给出的 bcc-Mo 结果如下，使用的是随教程提供的势：

| 方法 | C11 / GPa | C12 / GPa | C44 / GPa |
|---|---:|---:|---:|
| LAMMPS-MLIP | 473 | 155 | 106 |
| 有限差分 MLIP | 470 | 156 | 106 |
| 有限差分 DFT | 472 | 154 | 95 |

重新训练的模型受优化等因素影响，不必逐位复现这些数值；重点是量纲、趋势和偏差合理。表中 C44 也说明整体拟合良好不代表每个物性都与 DFT 一致。[官方弹性常数步骤](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-5)

### 7.4 真正的独立测试应怎样报告

对于构型大小不同的数据，推荐分别报告每原子能量 RMSE 和逐分量力 RMSE：

$$
\mathrm{RMSE}_{E/N}=\sqrt{\frac1S\sum_s\left(\frac{E_s^{\mathrm{MTP}}-E_s^{\mathrm{DFT}}}{N_s}\right)^2},
$$

$$
\mathrm{RMSE}_F=\sqrt{\frac{\sum_s\sum_{i=1}^{N_s}\sum_{a=x,y,z}(F_{sia}^{\mathrm{MTP}}-F_{sia}^{\mathrm{DFT}})^2}{3\sum_sN_s}}.
$$

应力若以 GPa 报告，应先逐构型除以其体积并换单位，再计算误差。不同版本程序日志中的总量、每原子量或分量误差定义可能不同，论文中应写明自己的口径。

还应按组成、温度、相和缺陷分别列出误差及极端值，并报告与目标直接相关的性质。模型没有被用来训练过某个构型，不代表该构型就与训练分布独立；近重复与同轨迹相关性都需要检查。

## 8. 将势函数接入 LAMMPS

在包含 `pot.mtp` 的工作目录建立 `mlip.ini`：

```text
mtp-filename pot.mtp
select FALSE
```

在与该 MLIP-2 接口匹配的 LAMMPS 输入中使用：

```text
pair_style mlip mlip.ini
pair_coeff * *
```

该写法已经由官方 Tutorial 1/2 使用。`mlip.ini` 在这里指定模型和选择模式；离线拟合参数通过 `mlp train` 传入。[官方接口示例](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Step-5)

下面是自建体系的 NVT 检查模板。运行前必须准备含正确质量和原子类型的 `input.data`，并把 `lmp` 换成实际配套可执行文件：

```text
units metal
atom_style atomic
boundary p p p
read_data input.data

pair_style mlip mlip.ini
pair_coeff * *

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
thermo 100
thermo_style custom step temp pe ke etotal press vol

run 0
velocity all create 300.0 314159 mom yes rot no dist gaussian
timestep 0.001
fix thermostat all nvt temp 300.0 300.0 0.1
run 10000
```

在 `units metal` 下时间单位是 ps，因此示例时间步为 1 fs，总长度 10 ps。这只是初步检查，不代表充分验证。应先确认 `run 0` 的能量与 MLIP 单点计算一致，再检查温度、压力、结构和外推等级。测试能量守恒时使用合适的 NVE 条件，不能要求带恒温器的 NVT 总能量像孤立体系一样守恒。

## 9. 主动学习：复现官方 Nb 示例，再替换为 DFT

### 9.1 先分清参考模型

官方 Tutorial 2 对 16 原子的 bcc-Nb 开展 300 K、0.1 ns 的 MD，使用 Farkas EAM 作为廉价参考，方便完整演示迭代过程。页面虽然在若干位置沿用了 “ab initio” 字样，实际 EAM 计算并不是第一性原理计算。[Tutorial 2 概览](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Home)

如果照原样跑完，你得到的是逼近该 EAM 参考的模型。要训练 DFT-MTP，应把初始标签和新增标签都按一致的 DFT 协议生成，或另行建立明确的多保真方法；不能无标识混合 EAM 和 DFT 标签。

### 9.2 外推等级意味着什么

主动学习用当前模型的参数敏感性和活动集评估新构型。当采用线性化表达时，可把候选行记为 $b$、活动矩阵记为 $A$，形式上用：

$$
c=bA^{-1},\qquad\gamma=\max_k|c_k|
$$

理解新信息是否超出已有集合。实际使用的行和选择方式以实现为准。maxvol/D-optimal 选择的作用是在模型特征或导数空间中减少冗余，优先补充信息不足的方向。[主动学习原始方法](https://arxiv.org/abs/1611.09346)

$\gamma$ 不是 eV 单位的误差，也不是“正确概率”。较小的 $\gamma$ 不能保证 DFT 标签无噪声、模型容量充分或磁性处理正确；应在代表性样本上用真实 DFT 误差校验它的诊断意义。

### 9.3 官方初始势训练

在 `tutorial-2/1.training_initial_potential/`，使用官方初始构型，复制完整 `08.mtp`，设置 `min_dist=2.8`、`max_dist=5`，再执行：

```bash
../mlp mindist train_init.cfg
../mlp train 08.mtp train_init.cfg --trained-pot-name=init.mtp
```

将 `init.mtp`、`train_init.cfg` 放入 `2.active_learning_md_lammps/`。按官方目录准备 `mlp` 和 `lmp`，然后进入该目录执行下面步骤。[Nb 初始势](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Step-1)

### 9.4 初始化当前势和活动集

```bash
cp init.mtp curr.mtp
cp train_init.cfg train.cfg
../mlp calc-grade curr.mtp train.cfg train.cfg out.cfg \
  --als-filename=state.als
```

此处两个 `train.cfg` 处于不同参数位置：一个用于建立活动集，另一个是要评估的构型集合。`out.cfg` 是输出，`state.als` 是随后 MD 使用的状态文件。[官方主动学习步骤](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Step-2)

### 9.5 开启在线监测

`mlip.ini` 设置为：

```text
mtp-filename curr.mtp
select TRUE
select:threshold 2.1
select:threshold-break 10.0
select:save-selected preselected.cfg
select:load-state state.als
```

然后运行：

```bash
../lmp -in in.nb_md
```

官方示例用 2.1 作为保存阈值、10.0 作为终止阈值。超过前者的构型进入候选文件；超过后者时停止 MD，以便补充数据后再探索。这是示例参数，不是跨体系通用的误差标准。

每轮最好使用独立目录保存输入、日志、候选集与状态文件，避免旧的 `preselected.cfg` 残留导致重复标注或错误判断是否收敛。还要区分“达到外推终止阈值”和时间步、邻居表、原子类型等普通程序错误。

### 9.6 从候选集中选择代表结构

```bash
../mlp select-add curr.mtp train.cfg preselected.cfg selected.cfg \
  --als-filename=state.als
```

选完仍建议用构型 ID 或结构匹配检查是否已经标注；不要仅根据 `selected` 或 `diff` 文件名认定其中每条都是全新数据。具体版本可能有不同输出行为，官方仓库已有相关使用报告。[选择命令说明](https://gitlab.com/ashapeev/mlip-2-paper-supp-info/-/blob/master/manual.pdf)、[重复选择问题记录](https://gitlab.com/ashapeev/mlip-2/-/issues/15)

### 9.7 将 EAM 标注替换为 VASP 单点计算

对正式 DFT 项目，执行以下顺序：

1. 导出选择的几何结构；确认每种元素的数量、顺序和全局 type。
2. 为每个构型建立单独 VASP 目录，准备与初始数据一致的 INCAR、KPOINTS 和匹配顺序的 POTCAR。
3. 固定构型做静态单点计算。典型设置涉及 `NSW=0` 和静态计算模式，但完整参数必须使用本体系已经收敛验证的协议。
4. 检查电子收敛、能量、力和应力，失败计算先修复，不能直接入库。
5. 转换结果，恢复全局元素编号，确认几何与送算构型相同。
6. 把通过检查的新标签存入 `labeled.cfg`。

MLIP-2 导出命令为：

```bash
../mlp convert-cfg selected.cfg POSCAR --output-format=vasp-poscar
```

单构型可能输出 `POSCAR`，多构型则输出 `POSCAR0`、`POSCAR1` 等。CFG 的数字类型不自动决定 POTCAR 的化学身份，因此必须根据全局映射核对导出文件和赝势顺序。

每个 DFT 任务完成后，用前面的 `convert-cfg ... --input-format=vasp-outcar` 导回；多组分时先重映射再合并。这里最重要的检查是：新能量和力来自 DFT，而不是候选文件残留的 MTP 预测。

### 9.8 合并、回训和重建状态

在一个新的轮次目录中准备旧的 `train.cfg`、`curr.mtp`、新 `labeled.cfg` 和固定 `valid.cfg`。下面用 MLIP 的转换/追加功能生成新训练文件，避免覆盖旧数据：

```bash
../mlp convert-cfg train.cfg train_next.cfg
../mlp convert-cfg labeled.cfg train_next.cfg --append

../mlp train curr.mtp train_next.cfg \
  --trained-pot-name=next.mtp \
  --valid-cfgs=valid.cfg \
  --update-mindist

../mlp calc-errors next.mtp valid.cfg
../mlp calc-grade next.mtp train_next.cfg train_next.cfg graded_next.cfg \
  --als-filename=next.als
```

`--append` 只负责格式追加，不负责去重和元素映射。重新训练时保留上一轮实际使用的权重和优化选项；上面的命令突出新增数据流程，默认参数未必适合你的模型。

通过本轮检查后，在下一轮目录把 `next.mtp`、`train_next.cfg`、`next.als` 分别作为当前势、训练集和状态，更新 `mlip.ini`，再开始新探索。模型参数改变后应重建与之对应的主动学习状态，避免混用旧 `.als`。

### 9.9 何时结束主动学习

官方演示在指定 MD 不再产生预选构型时停止。在研究项目中，还应确认多个随机初态、目标温度/压力/应变条件和相关缺陷都已探索，并且新增数据数量、验证误差和目标物性逐轮稳定。

固定最终测试集应在模型选择完成后评估。探索覆盖不足时，“本轮没有新构型”只能说明本轮轨迹没有触发阈值，不能证明所有适用范围都已覆盖。

## 10. 迁移到高熵合金的具体设计

### 10.1 一份初始数据设计表

以五元素 fcc 合金为例，下面是研究设计建议，构型数和温度点由预算及验证结果确定：

| 数据组 | 构型来源 | 主要用途 | 实施时的检查 |
|---|---|---|---|
| 随机/SQS 晶体 | 多个独立化学排列和超胞 | 无序固溶体基线 | 元素比例、局域相关、有限尺寸 |
| 热快照 | 多个温度/体积的 AIMD | 热振动与非平衡力 | 去相关、独立轨迹、DFT 收敛 |
| 变形晶体 | 体积与剪切扰动 | EOS、弹性和应力 | 避免把晶格线性应变当体积应变 |
| 化学交换构型 | MC、受控交换、必要的 DFT-MC | SRO 与偏聚 | 交换能差及有序倾向 |
| fcc/hcp 与层错 | 不同堆垛及滑移路径 | 层错能、竞争相 | 同组成、同边界条件的能量差 |
| 目标缺陷 | 空位、表面、晶界等 | 对应缺陷物性 | 低配位、局部畸变及尺寸效应 |
| 主动学习新增组 | MTP-MD/MC 的外推构型 | 补足实际模拟缺口 | 真实 DFT 回标、去重和来源追踪 |

如果目标包含多成分优化，就在这些组中增加组成维度；如果目标只限固态弹性，就不必无条件把所有液态、表面和裂纹类型都纳入初始训练。

### 10.2 三项论文经验及其应用边界

**随机合金局域环境与层错能。** Hodapp 和 Shapeev 的 MoNbTa 工作先构造大量小超胞候选，再通过主动选择减少 DFT 计算。适合借鉴的是“根据目标大尺度构型设计可计算的代表环境”；MoNbTa 是三元中熵合金，不能把其案例的数据规模当作五元体系所需数量。[原始论文](https://doi.org/10.1103/PhysRevMaterials.5.113802)

**化学与振动自由度都需要检查。** Ta-V-Cr-W 合金族的 MTP/GM-NN 对比研究同时考虑低温化学构型与高温 MD 数据。可借鉴的是分组数据和分组误差分析，以及区分能量与力的学习效果；主动学习对所有指标都一定优于合理随机采样，并不是该类研究能普遍保证的结论。[原始论文](https://www.nature.com/articles/s41524-023-01073-w)

**SRO 需要化学采样。** CrCoNi 工作系统比较训练集的化学采样、竞争相和液相内容，并评估多个独立训练的 MTP。对研究 SRO 的项目，直接启示是验证元素交换能差与平衡化学关联，而不能只依赖总体 RMSE。CrCoNi 也是三元中熵合金；其训练策略可为五元 HEA 提供参考，但需在目标体系重新验证。[原始论文](https://www.nature.com/articles/s41524-025-01722-2)

### 10.3 大尺寸模拟中的 DFT 成本

若直接把十万原子的外推构型送去 DFT，主动学习并没有解决标注成本问题。可以先在能做 DFT 的超胞中开展针对性探索，或采用经过验证的局域环境提取与重建方法。

MLIP-3 的重要扩展是对原子局域环境开展主动学习，更适合大体系中稀少的局域异常环境。但局域提取不意味着随便切一个小球就能得到无边界影响的 DFT 标签；边界重建、邻域完整性和电子结构仍需处理。[MLIP-3 论文](https://arxiv.org/abs/2304.13144)

使用 MLIP-3 时应从[其官方仓库](https://gitlab.com/ashapeev/mlip-3)对应手册和示例起步，重新核验训练命令、主动学习状态文件和 LAMMPS 接口。本文 MLIP-2 的配置块不承诺与 MLIP-3 直接兼容。

## 11. 常见问题与排查顺序

| 现象 | 优先排查 | 处理思路 |
|---|---|---|
| CFG 解析失败 | 字段名称、列数、构型边界、编码 | 先用单构型验证转换，再处理批量文件 |
| 多元素训练异常 | 全局 type 是否一致 | 核对每个 OUTCAR 的真实元素顺序 |
| 能量系统偏移 | 总能/每原子能量、DFT 能量字段、参考协议 | 统一标签定义；每原子形成能不要当总能输入 |
| 力误差高 | 标签错位、热扰动覆盖不足、权重不合适 | 核对原子顺序并补充非平衡构型 |
| 应力误差特别大 | GPa/kbar/eV、体积、符号、分量顺序 | 对一个压缩构型逐项手算核验 |
| 训练误差低而验证差 | 数据覆盖或模型复杂度 | 分组检查，增加相关数据并控制 level |
| 验证很好而 SRO 错 | 化学交换采样和小能量差 | 加入有序/交换构型并验证目标性质 |
| LAMMPS 不识别 mlip | 接口是否编译、版本是否匹配 | 按配套接口 README 构建 |
| MD 原子重叠或飞散 | 单位、质量、初始结构、时间步、外推 | 先做单点一致性与短程受控测试 |
| 主动学习一直选同一批 | 旧候选残留、状态未更新、去重失效 | 每轮独立目录、重建 ALS、按 ID 查重 |
| 外推等级低但物性错 | 模型能力、标签质量、评估构型覆盖 | 用真实 DFT 误差和物性测试诊断 |
| 增加迭代数没有改善 | 标签噪声、缺失环境或优化停滞 | 检查数据与独立初始化，不只增加 max-iter |

## 12. 最终归档与验收

建议保留如下研究目录。下列是组织示例，不要求 MLIP 必须使用这些文件名：

```text
hea-mtp/
  README.md
  metadata/
    versions.txt
    element-map.txt
    dft-settings.txt
    split-manifest.txt
  datasets/
    train.cfg
    valid.cfg
    test.cfg
  fits/
    level16-run01/
      init.mtp
      pot.mtp
      train.log
  active-learning/
    round-000/
    round-001/
  validation/
    errors/
    eos/
    elastic/
    target-properties/
  release/
    pot.mtp
    mlip.ini
    model-card.md
```

训练完成后，模型说明至少写清：适用元素和组成、相与温压范围、DFT 标签方法、数据量和切分方式、level 与径向设置、权重、独立测试误差、目标性质验证，以及尚未验证的使用范围。

验收时应同时满足以下条件：

- 元素类型、原子顺序、能量和应力单位已逐样本抽查。
- 训练/验证/测试没有轨迹或近重复泄漏。
- 所有新增训练标签都能追溯到实际参考计算。
- MLIP 与 LAMMPS 的同构型单点计算一致。
- 目标性质通过独立 DFT 或合适的实验对照，并记录偏差。
- 多初态的目标模拟稳定，主动学习监测覆盖已声明的工况。
- 势文件、对应数据、软件版本及训练日志已归档。

模型交付不应只有一个 `pot.mtp`。缺少元素映射、训练范围和验证记录时，后续使用者无法可靠判断它是否适合自己的模拟。

## 参考资料索引

1. [MLIP-2 教程首页](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/home)。
2. [MLIP-2 安装教程](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/installation-tutorial)与[源码安装说明](https://gitlab.com/ashapeev/mlip-2/-/blob/master/INSTALL.md)。
3. [Tutorial 1：Mo 被动训练](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-1/Home)。
4. [Tutorial 2：Nb 主动学习](https://gitlab.com/ashapeev/mlip-2-tutorials/-/wikis/Tutorial-2/Home)。
5. [官方用户手册与论文补充材料](https://gitlab.com/ashapeev/mlip-2-paper-supp-info)。
6. [本次核验的 MLIP-2 命令源码](https://gitlab.com/ashapeev/mlip-2/-/blob/fc3b0721b4d0eb1befcc4eb48269dfc6812f37d4/src/mlp/mlp_commands.cpp)。
7. [LAMMPS-MLIP-2 官方接口](https://gitlab.com/ashapeev/interface-lammps-mlip-2)。
8. Novikov 等：[The MLIP package: Moment Tensor Potentials with MPI and Active Learning](https://arxiv.org/abs/2007.08555)。
9. Podryabinkin、Shapeev：[Active learning of linearly parametrized interatomic potentials](https://arxiv.org/abs/1611.09346)。
10. Podryabinkin 等：[MLIP-3: Active learning on atomic environments with Moment Tensor Potentials](https://arxiv.org/abs/2304.13144)。
11. Hodapp、Shapeev：[Machine-learning potentials enable predictive and tractable high-throughput screening of random alloys](https://doi.org/10.1103/PhysRevMaterials.5.113802)。
12. [Performance of two complementary machine-learned potentials in modelling chemically complex systems](https://www.nature.com/articles/s41524-023-01073-w)。
13. [Capturing short-range order in high-entropy alloys with machine learning potentials](https://www.nature.com/articles/s41524-025-01722-2)。
