# ATAT（Alloy Theoretic Automated Toolkit）全面教程

> 基于 ATAT 官方文档（Axel van de Walle, Brown University）整理
> 文档版本：ATAT 3.x | 整理日期：2026-07-29

---

## 目录

1. [ATAT 概述](#1-atat-概述)
2. [理论基础](#2-理论基础)
3. [安装与配置](#3-安装与配置)
4. [核心文件格式](#4-核心文件格式)
5. [MAPS：二元合金团簇展开构建](#5-maps二元合金团簇展开构建)
6. [MMAPS：多组分团簇展开构建](#6-mmaps多组分团簇展开构建)
7. [EMC2：蒙特卡洛模拟](#7-emc2蒙特卡洛模拟)
8. [PHB：相界计算](#8-phb相界计算)
9. [MEMC2：多组分蒙特卡洛模拟](#9-memc2多组分蒙特卡洛模拟)
10. [MCSQS：特殊准随机结构生成](#10-mcsqs特殊准随机结构生成)
11. [CORRDUMP：关联函数计算](#11-corrdump关联函数计算)
12. [GCE：广义（张量）团簇展开](#12-gce广义张量团簇展开)
13. [FITFC：晶格动力学与声子计算](#13-fitfc晶格动力学与声子计算)
14. [FELEC：电子自由能](#14-felec电子自由能)
15. [MKTECI：合并自由能贡献](#15-mkteci合并自由能贡献)
16. [GENSTR：结构枚举](#16-genstr结构枚举)
17. [SQS2TDB：CALPHAD 数据库生成](#17-sqs2tdbcalphad-数据库生成)
18. [辅助工具](#18-辅助工具)
19. [典型工作流程](#19-典型工作流程)
20. [常见问题与技巧](#20-常见问题与技巧)

---

## 1. ATAT 概述

### 1.1 什么是 ATAT

ATAT（Alloy Theoretic Automated Toolkit）是由 Brown 大学 Axel van de Walle 教授开发的一套开源计算材料科学工具包。它的核心目标是通过**团簇展开（Cluster Expansion, CE）**方法，将第一性原理计算（如 VASP）得到的少量结构能量数据，外推为整个成分-构型空间的能量描述，从而实现：

- 合金相图的计算
- 有限温度热力学性质的预测
- 特殊准随机结构（SQS）的生成
- 声子谱和振动自由能的计算
- 弹性常数的预测
- 多组分体系的相稳定性分析

### 1.2 核心架构

ATAT 的设计哲学是**两步法**：

1. **构建团簇展开**（MAPS/MMAPS）：用第一性原理计算少量有序结构的能量，拟合出团簇展开的交互作用参数（ECI）。
2. **蒙特卡洛模拟**（EMC2/MEMC2/PHB）：利用已建立的团簇展开作为能量模型，进行大规模统计力学模拟，获取有限温度热力学性质。

### 1.3 工具列表

| 工具 | 功能 |
|------|------|
| `maps` | 二元合金团簇展开自动构建 |
| `mmaps` | 多组分（三元及以上）团簇展开构建 |
| `emc2` | 单相蒙特卡洛模拟（自由能面） |
| `memc2` | 多组分蒙特卡洛模拟 |
| `phb` | 两相平衡相界追踪 |
| `mcsqs` | 特殊准随机结构（SQS）生成 |
| `corrdump` | 关联函数计算、空间群确定、团簇枚举 |
| `gce` | 广义/张量团簇展开（应变、弹性常数等） |
| `fitfc` | 力常数拟合、声子态密度、色散关系 |
| `fitsvsl` / `svsl` | 声子频率拟合（替代方法） |
| `felec` | 电子自由能（从 DOS 计算） |
| `fmag` / `fempmag` | 磁性自由能 |
| `mkteci` | 合并多种自由能贡献的 ECI |
| `genstr` | 超结构枚举 |
| `gensqs` | SQS 生成（旧版） |
| `sqs2tdb` | 将 SQS 结果转换为 CALPHAD TDB 格式 |
| `calcelas` | 弹性常数计算 |
| `cellcvrt` | 结构文件格式转换 |
| `wycked` | Wyckoff 位置编辑 |
| `pollmach` | 并行任务调度 |
| `runstruct_vasp` | VASP 接口脚本 |

---

## 2. 理论基础

### 2.1 团簇展开（Cluster Expansion）

#### 2.1.1 基本思想

团簇展开的核心思想是：将合金的构型能量表示为**关联函数**的线性组合：

$$E(\sigma) = \sum_{\alpha} m_{\alpha} J_{\alpha} \Phi_{\alpha}(\sigma)$$

其中：
- $\sigma$ 表示一种原子构型（occupation configuration）
- $\alpha$ 遍历所有对称不等价的团簇（点、对、三体、四体等）
- $m_{\alpha}$ 是团簇 $\alpha$ 的多重度（multiplicity）
- $J_{\alpha}$ 是有效团簇交互作用（Effective Cluster Interaction, ECI）
- $\Phi_{\alpha}(\sigma)$ 是团簇 $\alpha$ 的关联函数（correlation function）

#### 2.1.2 占据变量与关联函数

对于二元体系，每个格点 $i$ 定义一个占据变量：

$$\sigma_i = \begin{cases} +1 & \text{if site } i \text{ is occupied by species A} \\ -1 & \text{if site } i \text{ is occupied by species B} \end{cases}$$

对于 $n$ 元体系（$n > 2$），使用 $(n-1)$ 维向量表示每个格点的占据状态，关联函数采用 cos/sin 函数：

$$\Phi_{\alpha}^{(k)} = \prod_{i \in \alpha} \cos\left(\frac{2\pi k \sigma_i}{n-1}\right) \quad \text{或} \quad \prod_{i \in \alpha} \sin\left(\frac{2\pi k \sigma_i}{n-1}\right)$$

#### 2.1.3 ECI 的拟合：SIM 方法

ATAT 使用 **Structure Inversion Method (SIM)**（也称为 Connolly-Williams 方法）来拟合 ECI。给定 $N$ 个已知能量的结构，求解超定方程组：

$$\mathbf{E} = \mathbf{\Pi} \mathbf{J}$$

其中 $\mathbf{\Pi}$ 是关联函数矩阵，$\mathbf{E}$ 是能量向量，$\mathbf{J}$ 是待求的 ECI 向量。

由于方程通常是欠定的（结构数少于团簇数），ATAT 采用**压缩感知（compressed sensing）**思想，通过最小化 ECI 的 $L_1$ 范数来选择最稀疏的团簇集合。

#### 2.1.4 交叉验证（Cross-Validation）

拟合质量通过留一法交叉验证（leave-one-out CV）评估：

$$CV = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (E_i - \hat{E}_{(-i)})^2}$$

其中 $\hat{E}_{(-i)}$ 是去掉第 $i$ 个结构后重新拟合得到的预测能量。CV 分数越低，CE 的预测能力越强。一般认为 CV < 25 meV/atom 表示良好的拟合。

#### 2.1.5 基态搜索

MAPS 在每次拟合后，通过枚举所有可能的超结构（使用 `genstr`），在 CE 能量面上搜索基态。如果预测的基态尚未被第一性原理验证，则自动提交计算。

### 2.2 蒙特卡洛方法

#### 2.2.1 半大正则系综

ATAT 的 MC 模拟采用**半大正则系综（semi-grand canonical ensemble）**：
- 总原子数 $N$ 固定
- 成分（浓度）可以变化
- 通过化学势 $\mu$ 控制成分

MC 移动为**原子交换**（swap）：随机选取两个不同种类的原子并交换位置。接受概率为：

$$P_{acc} = \min\left(1, \exp\left(-\frac{\Delta E - \mu \Delta N_A}{k_B T}\right)\right)$$

#### 2.2.2 热力学积分

自由能通过热力学积分获得。从高温（$T = \infty$，此时自由能可解析计算）出发，逐步降温：

$$G(T) = G(T_0) - \int_{T_0}^{T} S(T') \, dT'$$

ATAT 提供三种近似级别：
- **LTE（Low Temperature Expansion）**：低温展开，适用于有序相
- **MF（Mean Field）**：平均场近似
- **HTE（High Temperature Expansion）**：高温展开，适用于无序相

#### 2.2.3 相界确定

`phb` 代码同时运行两个相的 MC 模拟，通过追踪两相化学势相等的条件来确定相界：

$$\mu_1(T, x_1) = \mu_2(T, x_2)$$

### 2.3 晶格动力学

#### 2.3.1 力常数方法

`fitfc` 采用**实空间力常数**方法。通过对原子施加小位移，计算恢复力，拟合出原子间"弹簧"常数：

$$F_i = -\sum_j \Phi_{ij} \cdot u_j$$

其中 $\Phi_{ij}$ 是 $3 \times 3$ 力常数张量，$u_j$ 是原子 $j$ 的位移。

#### 2.3.2 准谐近似

通过在多个体积（应变）下计算声子谱，实现**准谐近似（Quasiharmonic Approximation, QHA）**：

$$F(V, T) = E_{static}(V) + F_{vib}(V, T)$$

对每个温度，最小化 $F(V, T)$ 得到平衡体积，从而获得热膨胀。

### 2.4 电子自由能

从电子态密度（DOS）计算电子激发对自由能的贡献：

$$F_{elec}(T) = -k_B T \int g(\epsilon) \ln\left(1 + e^{-(\epsilon - \mu_e)/k_B T}\right) d\epsilon$$

其中 $g(\epsilon)$ 是电子态密度，$\mu_e$ 是电子化学势（由电子数守恒确定）。

---

## 3. 安装与配置

### 3.1 编译安装

ATAT 使用 C++ 编写，依赖 LAPACK/BLAS 库。典型安装流程：

```bash
# 下载源码
git clone https://github.com/axelvandewalle/atat.git
cd atat

# 配置（编辑 makefile 设置编译器和库路径）
# 默认使用 g++，需要链接 -llapack -lblas

make
make install  # 安装到 /usr/local/bin 或指定目录
```

### 3.2 VASP 接口配置

ATAT 通过 `runstruct_vasp` 脚本与 VASP 交互。需要：

1. 确保 `vasp` 可执行文件在 PATH 中（或修改脚本中的路径）
2. 准备 `vasp.wrap` 模板文件（见第 4 节）
3. 配置 `pollmach` 的并行调度参数

### 3.3 pollmach 配置

`pollmach` 是 ATAT 的任务调度器，负责：
- 监控计算负载
- 在资源可用时触发新的计算
- 管理多个并行 VASP 任务

基本用法：
```bash
pollmach runstruct_vasp &
```

使用 `-lu`（lookup）选项可复用先前计算的 WAVECAR/CHGCAR 加速收敛：
```bash
pollmach -lu runstruct_vasp &
```

---

## 4. 核心文件格式

### 4.1 lat.in（晶格定义文件）

这是 ATAT 最核心的输入文件，定义了母体晶格（parent lattice）：

```
# 第一行：坐标系定义（二选一）
# 方式一：晶格参数
[a] [b] [c] [alpha] [beta] [gamma]
# 方式二：笛卡尔坐标
[ax] [ay] [az]
[bx] [by] [bz]
[cx] [cy] [cz]

# 接下来三行：晶格向量（在上述坐标系中）
[ua] [ub] [uc]
[va] [vb] [vc]
[wa] [wb] [wc]

# 原子位置及种类
[x1] [y1] [z1] Species1,Species2
[x2] [y2] [z2] Species3
...
```

**规则：**
- 逗号分隔的多种原子表示该位点可被多种原子占据（活性位点）
- 只有一种原子的位点不参与关联函数计算，但参与对称性确定
- `Vac` 表示空位

**示例：Cu-Au fcc 体系**
```
3.8 3.8 3.8 90 90 90
0   0.5 0.5
0.5 0   0.5
0.5 0.5 0
0 0 0 Cu,Au
```

**示例：Li-Co-O2 层状体系**
```
0.707 0.707 6.928 90 90 120
 0.3333  0.6667 0.3333
-0.6667 -0.3333 0.3333
 0.3333 -0.3333 0.3333
 0       0      0       Li,Vac
 0.3333  0.6667 0.0833  O
 0.6667  0.3333 0.1667  Co,Al
 0       0      0.25    O
```

### 4.2 vasp.wrap（VASP 参数模板）

```
[INCAR]
PREC = high
ISMEAR = -1
SIGMA = 0.1
NSW = 41
IBRION = 2
ISIF = 3
KPPRA = 1000
DOSTATIC
```

- `KPPRA`：k 点密度（每倒易原子约 1000）
- `DOSTATIC`：计算静态能量
- `ISIF = 3`：全自由度弛豫（包括体积）
- `NSW = 0`：不弛豫（用于力计算）

### 4.3 clusters.out（团簇文件）

```
# 对每个团簇：
[multiplicity]
[diameter]
[number of points]
[x1] [y1] [z1] [n_species-2] [cluster_function]
[x2] [y2] [z2] [n_species-2] [cluster_function]
...

# 空行分隔不同团簇
```

### 4.4 eci.out（ECI 文件）

每行一个 ECI 值，顺序与 `clusters.out` 中的团簇一一对应。ECI 已除以多重度。

### 4.5 str.out（结构文件）

格式与 `lat.in` 类似，但：
- 坐标系总是以 3×3 矩阵形式写出
- 每个位点只列出一种原子

### 4.6 rndstr.in（SQS 输入文件）

格式与 `lat.in` 类似，但原子位点使用**部分占据**表示：

```
3.8 3.8 3.8 90 90 90
0   0.5 0.5
0.5 0   0.5
0.5 0.5 0
0 0 0 Cu=0.5,Au=0.5
```

---

## 5. MAPS：二元合金团簇展开构建

### 5.1 原理

MAPS（Multicomponent Ab-initio Phase Stability）是 ATAT 的核心自动化引擎。它实现了一个**自适应迭代循环**：

1. 从初始结构集开始，拟合一个初步的 CE
2. 在 CE 能量面上搜索基态（通过枚举超结构）
3. 如果预测的基态尚未被 DFT 验证，提交新的 DFT 计算
4. 用新数据更新 CE
5. 重复直到收敛

收敛判据：
- 交叉验证分数 < 25 meV/atom
- 预测基态与真实基态一致
- ECI 随团簇尺寸衰减

### 5.2 输入文件

**必需文件：**
- `lat.in`：晶格定义（见 4.1 节）
- `xxxx.wrap`：第一性原理代码参数模板（如 `vasp.wrap`）

**可选文件：**
- `makelat` 工具可辅助生成 `lat.in`

### 5.3 运行方法

```bash
# 在工作目录中准备好 lat.in 和 vasp.wrap 后：
maps -d &

# -d 表示后台运行
# maps 会自动创建子目录并提交计算

# 当计算资源可用时：
touch ready

# 使用 pollmach 自动调度：
pollmach runstruct_vasp &
```

**工作流程：**
1. `maps` 启动后，创建初始结构目录
2. 用户（或 `pollmach`）创建 `ready` 文件通知 maps 可以提交新计算
3. `maps` 创建 `n/str.out` 和 `n/wait` 文件
4. `runstruct_vasp` 检测到 `wait` 文件，执行 VASP 计算
5. 计算完成后，能量写入 `n/energy`
6. `maps` 读取能量，更新 CE

### 5.4 输出文件

| 文件 | 内容 |
|------|------|
| `log.out` | 运行日志，包含基态检查结果和 CV 分数 |
| `clusters.out` | 选定的团簇集合 |
| `eci.out` | 拟合的 ECI |
| `gs_str.out` | 基态结构列表 |
| `fit.out` | 所有结构的拟合结果（能量、CV 误差等） |
| `predstr.out` | 预测能量（未计算的结构） |

### 5.5 收敛判断

运行 `mapsrep` 命令生成 5 个诊断图：
1. ECI vs 团簇尺寸
2. 拟合能量 vs DFT 能量
3. CV 分数随结构数变化
4. 基态凸包
5. ECI 衰减行为

**停止命令：**
```bash
touch stop       # 停止 maps
touch stoppoll   # 停止 pollmach
```

### 5.6 校准第一性原理参数

在正式运行前，建议校准：
- 截断能（ENCUT）
- k 点密度（KPPRA）
- 确保纯元素端点的能量收敛

---

## 6. MMAPS：多组分团簇展开构建

### 6.1 原理

`mmaps` 是 `maps` 的多组分扩展版本，支持三元、四元及更复杂的体系。核心算法相同（SIM + 基态搜索），但需要处理多维成分空间。

新增特性：
- 多维化学势空间
- 多子晶格支持
- 贝叶斯算法选项（`-fa=bayesian`）

### 6.2 输入文件

**必需：**
- `lat.in`：多组分晶格定义

**可选：**
- `ref_energy.in`：参考能量（每活性位点），用于计算形成能
- `nbclusters.in`：手动指定团簇数量（可在运行中修改，修改后 `touch refresh`）
- `crange.in`：限定拟合的成分范围
- `weights.in`：手动设置结构权重

**crange.in 示例：**
```
1.0*Al -1.0*Li +0.1*Co >= 0.5
```

### 6.3 运行方法

```bash
mmaps &
pollmach runstruct_vasp &
```

### 6.4 输出文件

| 文件 | 内容 |
|------|------|
| `maps.log` | 日志（警告、CV 分数） |
| `atoms.out` | 所有原子种类列表 |
| `fit.out` | 拟合结果（浓度向量、能量、拟合能量、权重） |
| `predstr.out` | 预测能量及状态（b=busy, e=error, u=unknown, g=ground state） |
| `gs.out` | 基态能量 |
| `gs_connect.out` | 基态凸包连接关系 |
| `gs_str.out` | 基态结构 |
| `chempot.out` | 稳定各相的化学势值 |
| `eci.out` | ECI |
| `clusters.out` | 团簇集合 |
| `ref_energy.out` | 使用的参考能量 |

### 6.5 通信协议

`mmaps` 与外部脚本的通信遵循严格的文件协议：

```
脚本创建 ready → mmaps 创建 n/str.out + n/wait → 脚本删除 wait，开始计算
→ 计算完成写入 n/energy（或 n/error）→ mmaps 更新 CE
```

### 6.6 贝叶斯算法

```bash
mmaps -fa=bayesian
```

贝叶斯方法提供更好的不确定性量化，特别适合数据稀疏的早期阶段。详见：
https://doi.org/10.1016/j.commatsci.2023.112571

---

## 7. EMC2：蒙特卡洛模拟

### 7.1 原理

`emc2` 在半大正则系综中执行蒙特卡洛模拟，计算单相的自由能面 $G(T, \mu)$。

**关键特点：**
- 固定总原子数，允许成分波动
- 化学势 $\mu$ 控制平衡成分
- 自动检测相变并停止
- 通过热力学积分计算自由能

### 7.2 输入文件

- `lat.in`：晶格定义
- `clusters.out`：团簇集合（由 maps 生成）
- `eci.out`：ECI（由 maps 生成）
- `gs_str.out`：基态结构（提供初始构型）

### 7.3 命令行参数

```bash
emc2 -gs=<n> -er=<R> -dx=<prec> -T0=<T_start> -T1=<T_end> -dT=<step> \
     -mu0=<mu_start> -mu1=<mu_end> -dmu=<mu_step> -k=8.617e-5
```

| 参数 | 含义 |
|------|------|
| `-gs=<n>` | 选择相（0=无序相，1~N=各基态） |
| `-er=<R>` | 模拟胞半径（Å），决定超胞大小 |
| `-dx=<prec>` | 成分精度目标（自动确定模拟时长） |
| `-T0, -T1, -dT` | 温度范围和步长（K） |
| `-T0, -T1, -db` | 用倒数温度步长（从 $T=\infty$ 开始） |
| `-mu0, -mu1, -dmu` | 化学势范围和步长（无量纲） |
| `-k=8.617e-5` | 玻尔兹曼常数（eV/K） |
| `-abs` | 化学势使用绝对值（eV） |
| `-eq=<n>` | 平衡步数（手动设置） |
| `-n=<n>` | 采样步数（手动设置） |
| `-can` | 正则系综模式 |
| `-innerT` | 内部温度循环 |

### 7.4 无量纲化学势

ATAT 默认使用**无量纲化学势**，其物理含义为：
- $\mu = i$（整数）时，在 $T = 0$ 稳定第 $i$ 和第 $i+1$ 个基态之间的两相平衡
- $\mu$ 在 $i$ 和 $i+1$ 之间时，在 $T = 0$ 稳定第 $i+1$ 个基态

这使得用户无需预先知道化学势的绝对值即可设置模拟。

### 7.5 输出

标准输出包含统计平均值：
- 平均成分
- 平均能量
- 比热
- 短程序参数

积分量（通过热力学积分获得）：
- Gibbs 自由能
- 熵

### 7.6 示例

```bash
# 计算第 2 个基态相的自由能面
emc2 -gs=2 -er=8 -dx=0.001 -T0=0 -T1=2000 -dT=50 \
     -mu0=1.5 -mu1=2.5 -dmu=0.1 -k=8.617e-5

# 从无穷温度开始（使用倒数温度步长）
emc2 -gs=0 -er=10 -dx=0.001 -T0=100000 -T1=300 -db=0.0001 \
     -mu0=1.0 -mu1=2.0 -dmu=0.05
```

---

## 8. PHB：相界计算

### 8.1 原理

`phb`（Phase Boundary）同时运行两个相的 MC 模拟，追踪两相共存条件。在每个温度步，调整化学势使两相化学势相等，从而确定相界上的 $(T, x)$ 点。

### 8.2 命令行参数

```bash
phb -gs1=<n1> -gs2=<n2> -T=<T_start> -dT=<step> -er=<R> -dx=<prec>
```

| 参数 | 含义 |
|------|------|
| `-gs1, -gs2` | 两个相的编号 |
| `-d1, -d2` | 不同母体晶格的 CE 目录 |
| `-T` | 起始温度（省略则从 0 K 开始） |
| `-dT` | 温度步长 |
| `-mu` | 起始化学势（有限温度起点时需要） |
| `-er` | 模拟胞半径 |
| `-dx` | 精度 |
| `-ltep` | 使用低温展开 |
| `-dmu` | 化学势步长 |

### 8.3 输出

输出列：温度 $T$、化学势 $\mu$、相 1 成分 $x_1$、相 2 成分 $x_2$、两相能量。

### 8.4 不同晶格间的相界

```bash
phb -gs1=1 -gs2=1 -d1=../fcc_CE -d2=../bcc_CE -T=0 -dT=50 -er=8
```

### 8.5 第三相出现

如果在追踪过程中第三相变得稳定，`phb` 会报告并停止。用户需要：
1. 确定新的两相平衡对
2. 从三相点重新开始

---

## 9. MEMC2：多组分蒙特卡洛模拟

### 9.1 原理

`memc2` 是 `emc2` 的多组分版本，处理多维化学势空间。对于 $n$ 元体系，需要 $(n-1)$ 个独立化学势。

### 9.2 使用方法

```bash
memc2 -gs=<n> -er=<R> -dx=<prec> -T0=<T> -T1=<T> -dT=<step> \
      -mu0=<mu1>,<mu2>,... -dmu=<dmu1>,<dmu2>,...
```

化学势以逗号分隔的向量形式给出。

### 9.3 与 mmaps 的配合

`mmaps` 生成的 `chempot.out` 文件提供了稳定各相的化学势值，可直接用于设置 `memc2` 的模拟参数。

---

## 10. MCSQS：特殊准随机结构生成

### 10.1 原理

特殊准随机结构（Special Quasirandom Structure, SQS）是一种有限大小的超胞，其短程序参数（关联函数）尽可能接近完全随机固溶体。

`mcsqs` 使用**蒙特卡洛模拟退火**算法：
1. 从随机构型开始
2. 随机交换原子
3. 以目标函数为"能量"进行 Metropolis 接受/拒绝
4. 逐步降低"温度"（模拟退火）

**目标函数：**

$$f = \frac{\sum_{\alpha} |\Phi_{\alpha}^{SQS} - \Phi_{\alpha}^{random}| \cdot e^{-w_d \cdot d_{\alpha}} \cdot w_n^{n_{\alpha}-2}}{\text{normalization}} - \sum_p w_r \cdot w_n^{p-2} \cdot d_{min}^{(p)}$$

其中：
- $d_{\alpha}$ 是团簇直径（归一化到最近邻距离）
- $n_{\alpha}$ 是团簇中的点数
- $d_{min}^{(p)}$ 是 $p$ 点及以下最小不匹配团簇的直径
- $w_r, w_n, w_d$ 是可调权重参数（默认 $w_n=1, w_d=0$）

### 10.2 输入文件

**rndstr.in**（随机态定义）：
```
3.8 3.8 3.8 90 90 90
0   0.5 0.5
0.5 0   0.5
0.5 0.5 0
0 0 0 Cu=0.5,Au=0.5
```

多子晶格示例：
```
0.707 0.707 6.928 90 90 120
 0.3333  0.6667 0.3333
-0.6667 -0.3333 0.3333
 0.3333 -0.3333 0.3333
 0       0      0       Li=0.75,Vac=0.25
 0.3333  0.6667 0.0833  O
 0.6667  0.3333 0.1667  Co=0.25,Ni=0.25,Al=0.5
 0       0      0.25    O
```

**注意：** 对称等价的位点必须具有相同的占据。使用 `rndstrgrp.out` 查看等价位点分组。如需覆盖，可使用不同标签（如 `Fe_a,Fe_b`）。

### 10.3 运行方法

```bash
# 第一步：生成团簇文件
mcsqs -2=8 -3=5 -4=4.2

# 第二步：运行 SQS 搜索
mcsqs -n=32 -2=8 -3=5 -4=4.2
```

| 参数 | 含义 |
|------|------|
| `-n=<N>` | 超胞中的原子数 |
| `-2=<R>` | 对团簇最大范围 |
| `-3=<R>` | 三体团簇最大范围 |
| `-4=<R>` | 四体团簇最大范围 |
| `-rc` | 使用自定义超胞（`sqscell.out`） |
| `-wr, -wn, -wd` | 目标函数权重 |
| `-T` | MC 温度 |
| `-pf=<file>` | 参数文件（运行时可调） |

### 10.4 输出文件

| 文件 | 内容 |
|------|------|
| `bestsqs.out` | 最佳 SQS 结构（ATAT 格式） |
| `bestcorr.out` | 关联函数对比（SQS vs 随机态） |
| `rndstrgrp.out` | 等价位点分组信息 |
| `mcsqs.log` | 日志 |

**bestcorr.out 格式：**
```
[点数] [直径] [SQS关联] [随机态关联] [差值]
```

### 10.5 停止与后处理

```bash
touch stopsqs   # 优雅停止

# 后处理：检查未包含在目标函数中的关联
corrdump -l=rndstr.in -ro -noe -2=10 -3=6 -s=bestsqs.out
corrdump -l=rndstr.in -ro -noe -2=10 -3=6 -s=bestsqs.out -rnd
```

### 10.6 自定义超胞

创建 `sqscell.out`：
```
[超胞数量]
[u1x] [u1y] [u1z]
[u2x] [u2y] [u2z]
[u3x] [u3y] [u3z]
...
```

超胞向量以 `rndstr.in` 中定义的轴的倍数表示（通常为整数）。

### 10.7 引用

> A. van de Walle, P. Tiwary, M. de Jong, D.L. Olmsted, M. Asta, A. Dick, D. Shin, Y. Wang, L.-Q. Chen, Z.-K. Liu, "Efficient stochastic generation of special quasirandom structures," Calphad Journal 42, pp. 13-18 (2013)

---

## 11. CORRDUMP：关联函数计算

### 11.1 原理

`corrdump` 是 ATAT 的底层工具，负责：
- 确定晶格的空间群
- 枚举对称不等价团簇
- 计算给定结构的关联函数

### 11.2 基本用法

```bash
# 枚举团簇
corrdump -l=lat.in -2=5 -3=4 -4=3.5

# 计算结构的关联函数
corrdump -l=lat.in -s=str.out -2=5 -3=4

# 使用 rndstr.in 格式（-ro 选项）
corrdump -l=rndstr.in -ro -noe -nop -clus -2=8 -3=5
```

### 11.3 选项

| 选项 | 含义 |
|------|------|
| `-l=<file>` | 晶格文件 |
| `-s=<file>` | 结构文件 |
| `-2=<R>` | 对团簇范围 |
| `-3=<R>` | 三体范围 |
| `-4=<R>` | 四体范围 |
| `-5=<R>` | 五体范围 |
| `-6=<R>` | 六体范围 |
| `-ro` | 读取 rndstr.in 格式 |
| `-noe` | 跳过空团簇 |
| `-nop` | 跳过点团簇 |
| `-clus` | 仅生成团簇 |
| `-rnd` | 输出随机态关联 |

### 11.4 输出文件

- `sym.out`：空间群对称操作
- `clusters.out`：团簇列表
- 标准输出：关联函数值（一行，按 clusters.out 顺序）
- `corrdump.log`：结构到理想晶格的映射调整

---

## 12. GCE：广义（张量）团簇展开

### 12.1 原理

广义团簇展开（Generalized Cluster Expansion）将 CE 的概念从标量（能量）推广到**张量**性质。例如：
- 应变（rank 2 张量）
- 弹性常数（rank 4 张量）

$$P_{ij}(\sigma) = \sum_{\alpha} m_{\alpha} J_{\alpha}^{ij} \Phi_{\alpha}(\sigma)$$

### 12.2 gcetensor.in 格式

```
[rank]
[对称等价索引对列表]
...
```

**应变（rank 2）：**
```
2
0 1
```

**弹性常数（rank 4）：**
```
4
0 1
2 3
0 2 1 3
```

### 12.3 运行方法

```bash
gce -l=lat.in -s=str.out -2=5 -3=4
```

### 12.4 clusters.out 中的张量数据

当"对象"是张量时，clusters.out 中每个团簇后附加：
```
[rank]
[3 3 3 ...]  (rank 个 3)
[3^rank]     (张量元素总数)
[张量元素...]
```

### 12.5 应用

- **成分应变（Constituent Strain）**：`csfit` 工具利用 GCE 拟合成分依赖的应变
- **弹性常数 CE**：预测不同构型的弹性张量

---

## 13. FITFC：晶格动力学与声子计算

### 13.1 原理

`fitfc` 通过以下步骤计算振动性质：
1. 对平衡结构中的原子施加小位移
2. 用 DFT 计算恢复力
3. 拟实力常数（"弹簧模型"）
4. 对角化动力学矩阵获得声子谱

支持**准谐近似**：在多个体积下重复上述过程，获得体积依赖的声子频率。

### 13.2 工作流程

#### 步骤 1：结构弛豫

```bash
# 准备 str.out（未弛豫结构）
runstruct_vasp   # 获得 str_relax.out（弛豫后结构）
```

确保 `vasp.wrap` 中设置全自由度弛豫（`ISIF = 3`）。

#### 步骤 2：生成扰动

```bash
fitfc -er=11.5 -ns=3 -ms=0.02 -dr=0.1
```

| 参数 | 含义 |
|------|------|
| `-er=<R>` | 位移原子的周期像之间的最小距离 |
| `-ns=<N>` | 应变级别数（1=纯谐，>1=准谐） |
| `-ms=<s>` | 最大应变（0.02 = 2%） |
| `-dr=<d>` | 位移幅度 |
| `-nrr` | 不重新弛豫（立方对称或各向同性膨胀） |

此命令创建 `vol_*` 子目录（每个应变级别一个）。

#### 步骤 2b：重新弛豫（准谐模式）

```bash
pollmach runstruct_vasp &
# 修改 vasp.wrap：允许除体积外所有自由度弛豫
touch stoppoll   # 完成后
```

#### 步骤 2c：重新生成扰动

```bash
fitfc -er=11.5 -ns=3 -ms=0.02 -dr=0.1
# 相同命令，但代码检测到新文件后继续
```

#### 步骤 3：计算力

```bash
pollmach -lu runstruct_vasp &
# vasp.wrap：NSW=0（不弛豫），使用 smearing，ICHARG=1
# -lu 复用 WAVECAR/CHGCAR 加速
```

#### 步骤 4：拟合力常数

```bash
fitfc -f -fr=5.0
# -f：拟合模式
# -fr：弹簧范围（不超过 -er 的一半）
```

### 13.3 不稳定模式处理

如果出现 "Unstable modes found. Aborting."：

```bash
# 查看不稳定模式
fitfc -fu   # 输出 vol_*/unstable.out

# 生成特定模式的扰动
fitfc -gu=<index>   # 正值=cos相位，负值=sin相位

# 在生成的子目录中运行 DFT
cd vol_*/p_uns_*
runstruct_vasp

# 重新拟合
fitfc -f -fr=5.0

# 强制生成 DOS（即使有不稳定模式）
fitfc -f -fr=5.0 -fn
```

### 13.4 声子色散曲线

创建输入文件（如 `disp.in`）：
```
20 0 0 0 0.5 0 0
20 0.5 0 0 0.5 0.5 0
20 0.5 0.5 0 0 0 0
```

每行定义一个 k 路径段：`[点数] [kx1] [ky1] [kz1] [kx2] [ky2] [kz2]`

```bash
fitfc -df=disp.in -si=str.out
```

### 13.5 输出文件

| 文件 | 内容 |
|------|------|
| `fitfc.log` | 日志 |
| `vol_*/vdos.out` | 声子态密度 |
| `vol_*/fc.out` | 力常数（拉伸和弯曲项） |
| `vol_*/svib_ht` | 高温振动熵极限（$k_B$/atom） |
| `vol_*/fvib` | 谐振动自由能 vs 温度 |
| `fitfc.out` | 温度、自由能、线性热膨胀 |
| `fvib` | 自由能（准谐/谐） |
| `svib` | 熵（eV/K） |
| `eos0.out` | 0K 状态方程 |
| `eigenfreq.out` | 声子色散频率 |

### 13.6 注意事项

- 在 CE 中使用时，加 `-me0` 选项避免静态能重复计算
- `-fr` 应从最近邻距离开始，逐步增大检查收敛
- 使用 `-dr=0.05`（较小位移）配合 `-lu` 选项效果更好

---

## 14. FELEC：电子自由能

### 14.1 原理

从电子态密度（DOS）计算电子激发自由能：

$$F_{elec}(T) = E_{band}(T) - T \cdot S_{elec}(T)$$

其中电子熵：

$$S_{elec} = -k_B \int g(\epsilon) [f \ln f + (1-f) \ln(1-f)] \, d\epsilon$$

$f(\epsilon)$ 是 Fermi-Dirac 分布函数。

### 14.2 输入文件

**dos.out**（电子态密度）：
```
[能量1] [DOS1]
[能量2] [DOS2]
...
```

### 14.3 使用方法

```bash
# 计算电子自由能
felec -T0=0 -T1=2000 -dT=10

# 将结果转换为 ECI
clusterexpand felec

# 合并到总 ECI
mkteci felec.eci
```

### 14.4 与振动自由能的合并

```bash
# 分别计算振动和电子 ECI 后：
mkteci fvib.eci felec.eci
# 生成合并的 eci.out
```

---

## 15. MKTECI：合并自由能贡献

### 15.1 功能

`mkteci` 将来自不同物理贡献的 ECI 合并为总的 ECI：

$$J_{\alpha}^{total} = J_{\alpha}^{static} + J_{\alpha}^{vib} + J_{\alpha}^{elec} + J_{\alpha}^{mag} + \cdots$$

### 15.2 使用方法

```bash
mkteci felec.eci fvib.eci fmag.eci
```

输出合并后的 `eci.out`，可直接用于 MC 模拟。

---

## 16. GENSTR：结构枚举

### 16.1 功能

`genstr` 枚举给定母体晶格的所有对称不等价超结构，用于：
- MAPS 中的基态搜索
- 系统性探索构型空间

### 16.2 使用方法

```bash
genstr -l=lat.in -n=4
# -n：超胞中的最大原子数
```

---

## 17. SQS2TDB：CALPHAD 数据库生成

### 17.1 功能

将 SQS 计算结果（形成能、振动自由能等）转换为 CALPHAD 格式的 TDB 文件，用于与商业热力学软件（如 Thermo-Calc）对接。

### 17.2 使用方法

```bash
sqs2tdb [options]
```

具体参数取决于体系，建议参考 `sqs2tdb -h`。

---

## 18. 辅助工具

### 18.1 pollmach

并行任务调度器：
```bash
pollmach runstruct_vasp &       # 基本用法
pollmach -lu runstruct_vasp &   # 复用波函数加速
touch stoppoll                   # 停止
```

### 18.2 runstruct_vasp

VASP 接口脚本：
```bash
runstruct_vasp          # 在当前目录运行 VASP
runstruct_vasp -h       # 帮助
```

功能：
1. 将 `vasp.wrap` + `str.out` 转换为 VASP 输入文件（POSCAR, INCAR, KPOINTS）
2. 运行 VASP
3. 提取能量写入 `energy` 文件

### 18.3 cellcvrt

结构文件格式转换（ATAT ↔ POSCAR ↔ CIF 等）。

### 18.4 wycked

Wyckoff 位置编辑器，用于修改结构中的原子占位。

### 18.5 calcelas

弹性常数计算，结合 GCE 框架。

### 18.6 makelat

辅助生成 `lat.in` 文件的交互工具。

---

## 19. 典型工作流程

### 19.1 二元合金相图计算（以 Cu-Au 为例）

```
步骤 1：准备输入
├── 创建 lat.in（fcc Cu-Au）
├── 创建 vasp.wrap（DFT 参数）
└── 校准 DFT 参数（截断能、k 点）

步骤 2：构建团簇展开
├── maps -d &
├── pollmach runstruct_vasp &
├── 等待收敛（CV < 25 meV/atom）
├── mapsrep  # 检查诊断图
└── touch stop

步骤 3：蒙特卡洛模拟
├── emc2 -gs=0 -er=10 -dx=0.001 -T0=0 -T1=1500 -dT=25 \
│        -mu0=-1 -mu1=3 -dmu=0.1
├── emc2 -gs=1 ...  # 各有序相
└── emc2 -gs=2 ...

步骤 4：相界计算
├── phb -gs1=0 -gs2=1 -T=0 -dT=25 -er=10 -dx=0.001
├── phb -gs1=0 -gs2=2 ...
└── phb -gs1=1 -gs2=2 ...

步骤 5：后处理
└── 绘制 T-x 相图
```

### 19.2 SQS 生成（以 Ni-Cr 合金为例）

```
步骤 1：准备 rndstr.in
└── 定义 fcc 晶格，Ni=0.8,Cr=0.2

步骤 2：生成 SQS
├── mcsqs -n=108 -2=8 -3=5 -4=4.2
├── 等待收敛或满意结果
└── touch stopsqs

步骤 3：验证
├── corrdump -l=rndstr.in -ro -noe -2=10 -3=6 -s=bestsqs.out
└── 对比 bestcorr.out 中的关联差异

步骤 4：DFT 计算
├── 将 bestsqs.out 转换为 POSCAR
└── 运行 VASP 获得 SQS 能量
```

### 19.3 声子计算（以 Ni3Al 为例）

```
步骤 1：结构弛豫
├── 准备 str.out（L12 Ni3Al）
├── runstruct_vasp  # 全弛豫
└── 获得 str_relax.out

步骤 2：准谐计算
├── fitfc -er=12 -ns=5 -ms=0.03 -dr=0.1
├── pollmach runstruct_vasp &  # 各体积弛豫
├── touch stoppoll
├── fitfc -er=12 -ns=5 -ms=0.03 -dr=0.1  # 重新生成扰动
├── pollmach -lu runstruct_vasp &  # 力计算
└── touch stoppoll

步骤 3：拟合与输出
├── fitfc -f -fr=5.5
├── 检查 vdos.out, fvib, svib
└── fitfc -df=disp.in  # 色散曲线
```

### 19.4 多组分体系（以 Li-Co-O2 为例）

```
步骤 1：准备多组分 lat.in
└── 包含 Li/Vac 和 Co/Al 子晶格

步骤 2：构建 CE
├── mmaps &
├── pollmach runstruct_vasp &
├── 监控 maps.log 中的 CV 和基态信息
└── touch stop

步骤 3：多组分 MC
├── 查看 chempot.out 确定化学势范围
├── memc2 -gs=1 -er=8 -dx=0.001 -T0=300 -T1=1000 -dT=50 \
│         -mu0=... -dmu=...
└── 分析相稳定性
```

---

## 20. 常见问题与技巧

### 20.1 MAPS 相关

**Q: CV 分数一直不收敛怎么办？**
- 检查 DFT 计算是否充分收敛（ENCUT, KPPRA）
- 确认结构正确映射到母体晶格
- 尝试增加 `-er` 范围以包含更多团簇
- 检查是否有结构计算失败（查看 error 文件）

**Q: 预测基态总是与真实不一致？**
- 可能需要更多结构数据
- 检查是否遗漏了重要的有序相
- 考虑增大枚举的超胞尺寸

**Q: 如何导入外部结构？**
- 将结构文件（`str.out` 格式）放入子目录
- MAPS 会自动扫描并尝试映射
- 注意：导入的结构必须是未弛豫的

### 20.2 MC 模拟相关

**Q: 模拟胞多大才够？**
- 逐步增大 `-er` 直到结果不再变化
- 一般 8-12 Å 对大多数体系足够
- 接近相变温度时可能需要更大

**Q: 如何判断相变？**
- `emc2` 会自动检测并停止
- 观察成分-化学势曲线的突变
- 比热发散是二级相变的标志

**Q: 无量纲化学势如何理解？**
- $\mu = 1.5$ 表示在 $T=0$ 时稳定第 2 个基态
- 范围 $[i, i+1]$ 对应第 $i+1$ 个基态的稳定区
- 使用 `-abs` 切换到绝对化学势（eV）

### 20.3 SQS 相关

**Q: 如何选择合适的超胞大小？**
- 原子数应为母体晶格原子数的整数倍
- 越大越好，但计算成本增加
- 典型：32, 64, 108, 128 原子

**Q: SQS 质量如何评估？**
- 查看 `bestcorr.out` 中差值列
- 所有包含在目标函数中的关联应接近零
- 用 `corrdump` 检查更大范围的关联

**Q: 对称等价位点必须相同占据吗？**
- 是的，这是默认要求
- 如需不同占据，使用不同标签（如 `Fe_a`, `Fe_b`）

### 20.4 声子计算相关

**Q: -fr 如何选取？**
- 从最近邻距离开始
- 逐步增大，检查声子 DOS 收敛
- 不超过 `-er` 的一半

**Q: 出现不稳定模式怎么办？**
- 可能是真实不稳定性（结构不在能量极小点）
- 也可能是拟合伪影
- 使用 `-fu` 查看，`-gu` 生成扰动验证
- 使用 `-fn` 强制输出（负频率表示不稳定）

### 20.5 一般技巧

- **并行效率**：使用 `pollmach -lu` 复用波函数
- **断点续算**：MAPS/MMAPS 的所有状态保存在文件中，可随时重启
- **手动干预**：可以手动创建 `n/energy` 或 `n/error` 文件
- **成分限制**：使用 `crange.in` 限定感兴趣的成分范围
- **权重调整**：后期微调时使用 `weights.in`
- **刷新**：修改 `nbclusters.in` 后执行 `touch refresh`

---

## 参考文献

1. A. van de Walle and M. Asta, "Self-driven lattice-model Monte Carlo simulations of alloy thermodynamic properties and phase diagrams," Modelling Simul. Mater. Sci. Eng. 10, 521 (2002).

2. A. van de Walle, G. Ceder, "Automating first-principles phase diagram calculations," J. Phase Equilib. 23, 348 (2002).

3. A. van de Walle, P. Tiwary, M. de Jong, et al., "Efficient stochastic generation of special quasirandom structures," Calphad 42, 13-18 (2013).

4. A. van de Walle, "Multicomponent multisublattice alloys, nonconfigurational entropy and other additions to the Alloy Theoretic Automated Toolkit," Calphad 33, 266 (2009).

5. A. van de Walle, G. Ceder, U. Waghmare, "First-principles theory of short-range order, size effects, and nanostructure formation in metallic alloys," Phys. Rev. Lett. 81, 3864 (1998).

6. ATAT 官方文档：https://axelvandewalle.github.io/www-avdw//atat/manual/manual.html

---

*本教程基于 ATAT 官方手册整理，涵盖 ATAT 3.x 版本的主要功能。具体命令选项请以 `command -h` 输出为准。*
