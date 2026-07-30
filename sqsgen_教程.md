# sqsgenerator（sqsgen）全面教程

> 基于 sqsgenerator 官方文档（https://sqsgenerator.readthedocs.io）整理
> 版本：sqsgenerator 0.5.6 | 整理日期：2026-07-30

---

## 目录

1. [概述](#1-概述)
2. [理论基础：Warren-Cowley 短程序参数](#2-理论基础warren-cowley-短程序参数)
3. [安装与配置](#3-安装与配置)
4. [核心概念：sqs.json 配置文件](#4-核心概念sqsjson-配置文件)
5. [参数详解](#5-参数详解)
6. [CLI 使用教程](#6-cli-使用教程)
7. [Python API 使用教程](#7-python-api-使用教程)
8. [子晶格优化](#8-子晶格优化)
9. [输出文件与结果分析](#9-输出文件与结果分析)
10. [实战示例](#10-实战示例)
11. [模板系统](#11-模板系统)
12. [高级用法与技巧](#12-高级用法与技巧)
13. [常见问题](#13-常见问题)

---

## 1. 概述

### 1.1 什么是 sqsgenerator

sqsgenerator（命令行工具名为 `sqsgen`）是一个用于生成**特殊准随机结构（Special Quasirandom Structure, SQS）** 的高性能计算工具。它通过最小化 Warren-Cowley 短程序（Short-Range Order, SRO）参数来寻找最接近理想随机固溶体的有限超胞构型。

### 1.2 核心特性

- **极快的计算速度**：核心算法用 C++ 编写，性能优异
- **两种搜索策略**：Monte-Carlo 随机搜索和系统性（字典序）穷举
- **多线程并行**：默认多线程运行，可选 MPI 支持（HPC 环境）
- **多子晶格同时优化**：单次运行中可同时优化多个子晶格
- **浏览器可用**：提供 WebAssembly 版本，无需安装即可在线使用
- **框架集成**：与 ASE、pymatgen、pyiron 无缝对接

### 1.3 三种使用方式

| 方式 | 适用场景 | 特点 |
|------|----------|------|
| Web App | 快速测试、教学 | 打开 sqsgen.gehringer.tech 即用 |
| Python 包 | 日常研究、脚本化 | CLI + Python API |
| 原生应用 | HPC 大规模计算 | MPI 并行，无需 Python |

### 1.4 引用

> D. Gehringer, M. Friak, D. Holec, "Models of configurationally-complex alloys made simple," Computer Physics Communications, 108664 (2023). DOI: 10.1016/j.cpc.2023.108664

---

## 2. 理论基础：Warren-Cowley 短程序参数

### 2.1 二元体系的 SRO 参数

sqsgenerator 的核心是计算和优化 **Warren-Cowley 短程序（WC-SRO）参数**。对于二元合金（A-B 体系），SRO 参数定义为：

$$\alpha_{AB} = 1 - \frac{N_{AB}}{NMx_Ax_B}$$

其中：
- $N_{AB}$：A-B 键（对）的数量
- $N$：胞内总原子数
- $x_A$（$x_B$）：A（B）的摩尔分数
- $M$：配位数（如 fcc 的 $M = 12$）

**物理含义：**
- $\alpha_{AB} < 0$：有序倾向（异类原子偏好相邻）
- $\alpha_{AB} > 0$：聚团倾向（同类原子偏好相邻）
- $\alpha_{AB} \approx 0$：完全无序（理想随机固溶体）

**sqsgenerator 的默认目标就是使 $\alpha \approx 0$。**

### 2.2 多配位壳层推广

实际晶体中不仅有最近邻，还有次近邻、第三近邻等。SRO 参数推广为：

$$\alpha^k_{AB} = 1 - \frac{N^k_{AB}}{NM^kx_Ax_B}$$

其中 $k$ 是配位壳层索引，$M^k$ 是第 $k$ 壳层的配位数。例如 fcc 晶格：
- 第 1 壳层：$M_1 = 12$
- 第 2 壳层：$M_2 = 6$
- 第 3 壳层：$M_3 = 24$

### 2.3 多组分推广

对于多组分体系（如 {A, B, C}），需要考虑所有原子对组合，SRO 参数成为三维张量：

$$\alpha^i_{\xi\eta} = 1 - \frac{N^i_{\xi\eta}}{NM^ix_{\xi}x_{\eta}}$$

其中 $\xi, \eta$ 遍历所有原子种类。

### 2.4 目标函数

将所有 SRO 参数压缩为单一标量目标函数：

$$\mathcal{O}(\sigma) = \sum_{i,\xi,\eta} \tilde{p}^i_{\xi\eta} \left| \tilde{\alpha}_{\xi\eta} - \alpha^i_{\xi\eta} \right|$$

其中：
- $\tilde{p}^i_{\xi\eta} = w^i \cdot p_{\xi\eta}$：合并的配对权重（shell_weights * pair_weights）
- $\tilde{\alpha}_{\xi\eta}$：目标 SRO 值（默认为 0，即完全随机）
- $\alpha^i_{\xi\eta}$：当前构型的实际 SRO 值

**计算优化：** 由于 $\frac{1}{NM^ix_\xi x_\eta}$ 不依赖于构型，可预计算为"预因子"：

$$f^i_{\xi\eta} = (NM^ix_\xi x_\eta)^{-1}$$

使得：

$$\alpha^i_{\xi\eta} = 1 - f^i_{\xi\eta} N^i_{\xi\eta}$$

### 2.5 搜索算法

sqsgenerator 通过两种方式探索构型空间：
- **random 模式**：随机打乱原子排列（Monte-Carlo），适合大体系
- **systematic 模式**：按字典序穷举所有排列，保证找到全局最优（仅适合小体系）

总排列数为：

$$N_{\text{iter}} = \frac{N!}{\prod_m^M N_m!}$$

---

## 3. 安装与配置

### 3.1 Web App（零安装）

直接访问 https://sqsgen.gehringer.tech 即可在浏览器中使用（多线程 WebAssembly）。

### 3.2 pip 安装

```bash
pip install sqsgenerator
```

预编译 wheel 支持：Linux (x86-64, aarch64)、macOS (universal2)、Windows (amd64)。
Python 版本：3.9 - 3.13。

### 3.3 conda 安装

```bash
conda install -c conda-forge sqsgenerator
```

### 3.4 从源码编译

**前置要求：**
- Python 3.9+
- C++20 兼容编译器（g++ 13+、clang++ 19+、Apple clang++ 15+、MSVC 19.29+）
- CMake 3.25+
- Git
- vcpkg（管理 C++ 依赖）

```bash
git clone --recursive https://github.com/dgehringer/sqsgenerator.git
cd sqsgenerator
pip install python/ -v
```

### 3.5 原生应用（HPC/MPI）

```bash
git clone --recursive https://github.com/dgehringer/sqsgenerator.git
cd sqsgenerator
./scripts/build-cli-mpi.sh build
```

编译后在 `build` 目录中生成 `sqsgen` 可执行文件。

### 3.6 验证安装

```bash
sqsgen --version
```

---

## 4. 核心概念：sqs.json 配置文件

sqsgenerator 使用 JSON（或 YAML）格式的配置文件来定义优化问题。默认文件名为 `sqs.json`。

### 4.1 最小配置示例

```json
{
    "structure": {
        "lattice": [
            [3.165, 0.0, 0.0],
            [0.0, 3.165, 0.0],
            [0.0, 0.0, 3.165]
        ],
        "coords": [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ],
        "species": ["W", "Re"],
        "supercell": [3, 3, 3]
    },
    "iterations": 50000000,
    "shell_weights": {
        "1": 1.0
    },
    "composition": {
        "W": 27,
        "Re": 27
    },
    "max_results_per_objective": 10
}
```

### 4.2 配置文件结构

| 键 | 必需 | 说明 |
|----|------|------|
| `structure` | 是 | 输入结构定义 |
| `composition` | 是 | 目标成分 |
| `iterations` | 否 | 迭代次数（默认 10^5） |
| `iteration_mode` | 否 | `random` 或 `systematic` |
| `sublattice_mode` | 否 | `interact` 或 `split` |
| `shell_weights` | 否 | 配位壳层权重 |
| `target_objective` | 否 | 目标 SRO 值（默认 0） |
| `pair_weights` | 否 | 原子对权重 |
| `max_results_per_objective` | 否 | 每个目标值保存的最大结构数 |
| `seed` | 否 | 随机种子（可重复性） |
| `thread_config` | 否 | 线程配置 |

---

## 5. 参数详解

### 5.1 structure（结构定义）

定义输入晶体结构。两种方式：

**方式一：直接在 JSON 中定义**
```json
{
  "structure": {
    "lattice": [[4.123, 0.0, 0.0], [0.0, 4.123, 0.0], [0.0, 0.0, 4.123]],
    "coords": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    "species": ["Cs", "Cl"],
    "supercell": [3, 3, 3]
  }
}
```

**方式二：从文件读取**
```json
{
  "structure": {
    "file": "ti-n.vasp",
    "supercell": [2, 2, 2]
  }
}
```

支持的文件格式（通过 ASE 或 pymatgen）：POSCAR、CIF、JSON 等。

**可选子键：**
- `reader`：指定读取后端（`ase` 或 `pymatgen`）
- `args`：传递给读取函数的额外参数
- `supercell`：超胞倍数 `[nx, ny, nz]`

### 5.2 composition（成分）

定义目标原子分布。**原子总数必须等于格点数。**

**简单模式（所有位点）：**
```json
{
  "composition": {
    "Ti": 18,
    "Al": 18,
    "Mo": 18
  }
}
```

**子晶格模式（指定 sites）：**
```json
{
  "composition": [
    {"sites": "Ti", "Ti": 16, "Al": 16},
    {"sites": "N", "N": 16, "B": 16}
  ]
}
```

**空位：** 使用 `"0"` 表示空位：
```json
{
  "composition": {"Al": 56, "0": 8}
}
```

### 5.3 iteration_mode（迭代模式）

| 模式 | 说明 |
|------|------|
| `random` | 随机打乱（Monte-Carlo），默认 |
| `systematic` | 字典序穷举所有排列（忽略 iterations 参数） |

`systematic` 仅适用于 `interact` 模式的小体系。

### 5.4 sublattice_mode（子晶格模式）

| 模式 | 说明 |
|------|------|
| `interact` | 所有原子相互交互，整体优化（默认） |
| `split` | 各子晶格独立优化 |

- `interact`：适用于固定某些原子位置、全局优化
- `split`：适用于多子晶格独立优化（如 TiN 中只优化 N 子晶格）

### 5.5 iterations（迭代次数）

- 默认：10^5（random 模式）
- `systematic` 模式下被忽略
- 典型值：10^6 ~ 10^8

### 5.6 shell_weights（壳层权重）

控制哪些配位壳层参与优化及其权重 $w^i$。

- 默认：$w^i = 1/i$（所有壳层）
- 仅第一壳层：`{"1": 1.0}`
- 前两个壳层：`{"1": 1.0, "2": 0.5}`

**性能提示：** 减少壳层数可显著提升速度。

### 5.7 target_objective（目标 SRO）

目标 SRO 值 $\tilde{\alpha}_{\xi\eta}$。默认为 0（完全随机）。

可设置为非零值来生成具有特定短程序的结构（如模拟有序化倾向）。

### 5.8 pair_weights（配对权重）

$\tilde{p}^i_{\xi\eta}$ 矩阵，用于差异化不同原子对的权重。

默认为空心矩阵（对角线为 0，非对角线为 $w^i/2$），即只优化异类原子对。

### 5.9 seed（随机种子）

```json
{"seed": 42}
```

- 设置后可重复结果
- **要求 `thread_config` 为 1**（多线程下不保证可重复性）
- 支持列表形式（每个子晶格一个种子）：`[42, null, 123]`

### 5.10 配位壳层自动检测

sqsgenerator 通过距离直方图自动检测配位壳层：
- `bin_width`：直方图 bin 宽度（默认 0.05 Angstrom）
- `peak_isolation`：峰隔离阈值（默认 0.25）
- `shell_radii`：手动指定壳层半径（覆盖自动检测）

---

## 6. CLI 使用教程

### 6.1 基本运行

```bash
# 使用默认 sqs.json
sqsgen run

# 指定输入文件
sqsgen run -i my_config.json
```

原生 CLI：
```bash
sqsgen -i my_config.json
```

### 6.2 查看结果

```bash
# 列出所有目标值和对应结构数
sqsgen output list

# 指定输出文件
sqsgen output -o re-w.first.sqs.mpack list
```

输出示例：
```
Mode: interact
min(O(sigma)): 0.00000
Num. objectives: 5

INDEX  OBJ.     N
0      0.00000  13
1      0.01852  12
2      0.03704  2
3      0.07407  1
4      0.25926  1
```

### 6.3 导出结构

```bash
# 导出最佳结构为 CIF
sqsgen output structure -f cif

# 导出指定目标值和索引的结构
sqsgen output structure -f cif --objective 4 --index 0

# 使用 pymatgen 后端
sqsgen output structure -f pymatgen.cif
```

文件命名规则：`sqs-{objective_index}-{structure_index}.{format}`

### 6.4 支持的输出格式

原生 CLI：CIF、POSCAR、PDB、JSON、sqsgen JSON
Python CLI：额外支持 ASE 和 pymatgen 的所有格式

### 6.5 模板命令

```bash
# 列出可用模板
sqsgen template list

# 使用模板生成输入文件
sqsgen template use re-w.first
```

### 6.6 计算辅助

```bash
# 计算总排列数
sqsgen compute total-permutations

# 估算计算时间
sqsgen compute estimated-time

# 查看解析后的参数
sqsgen params show sqs.json -p shell_weights
```

---

## 7. Python API 使用教程

### 7.1 核心函数

```python
from sqsgenerator import parse_config, optimize, load_result_pack
```

### 7.2 运行优化

```python
from sqsgenerator import parse_config, optimize

with open("re-w.first.json") as f:
    config = parse_config(f.read())

pack = optimize(config)
```

也接受 Python dict 和 numpy 数组作为输入。

### 7.3 加载已有结果

```python
from sqsgenerator import load_result_pack

with open("sqs.mpack", "rb") as f:
    pack = load_result_pack(f.read())
```

### 7.4 分析结果

```python
# 获取最佳结构
best = pack.best()

# 遍历所有目标值
for obj, solutions in pack:
    print(f"Objective: {obj}, Num. solutions: {len(solutions)}")
```

### 7.5 导出结构

```python
from sqsgenerator import write

for oi, (obj, solutions) in enumerate(pack):
    for si, solution in enumerate(solutions):
        write(solution.structure(), f"sqs-{oi}-{si}.pymatgen.cif")
```

### 7.6 分析 SRO 参数

**interact 模式：**
```python
solution = pack.best()

# 获取完整 SRO 数组 (num_shells, num_species, num_species)
sro_params = solution.sro()

# 获取特定原子对的 SRO（所有壳层）
sro_re_w = solution.sro("Re", "W")

# 获取特定壳层的所有 SRO
sro_shell1 = solution.sro(1)

# 获取单个值
alpha = float(solution.sro(1, "Re", "W"))
```

**split 模式：**
```python
solution = pack.best()
# 获取各子晶格结果
for sublattice_result in solution.sublattices():
    print(sublattice_result.sro())
```

### 7.7 SqsResult 属性

每个 `SqsResult` 对象包含：
- `structure()`：获取 Structure 对象
- `objective`：目标函数值
- `rank`：排列编号
- `permutation`：物种数组（字符串）

---

## 8. 子晶格优化

### 8.1 split 模式

将结构拆分为多个子晶格，独立优化每个子晶格。

**示例：TiN 中只在 N 子晶格上分布 B 和 N**
```json
{
  "iterations": 100000000,
  "sublattice_mode": "split",
  "structure": {
    "file": "ti-n.vasp",
    "supercell": [2, 2, 2]
  },
  "composition": [{
    "sites": "N",
    "B": 16,
    "N": 16
  }]
}
```

此配置中：
- Ti 子晶格保持不变
- 只优化 N-B、N-N、B-B 对
- 不涉及 Ti-Ti、Ti-N、Ti-B 对

### 8.2 interact 模式

所有原子参与交互。适用于：
- 固定某些原子位置（pinning）
- 全局优化所有原子对

### 8.3 多子晶格同时优化

```json
{
  "composition": [
    {"sites": "Ti", "Ti": 16, "Al": 16},
    {"sites": "N", "N": 16, "B": 16}
  ]
}
```

在 split 模式下，两个子晶格独立优化，一次运行完成。

---

## 9. 输出文件与结果分析

### 9.1 sqs.mpack 文件

输出为 MessagePack 二进制格式，包含：
- 输入配置
- 性能指标
- 所有计算的结构（压缩格式）
- 所有计算的 SRO 参数

结果按目标函数值 $\mathcal{O}(\sigma)$ 升序排列。

### 9.2 结果解读

```
INDEX  OBJ.     N
0      0.00000  13    <- 13 个结构达到完美随机（O=0）
1      0.01852  12
2      0.03704  2
```

- INDEX 0 对应最佳结果
- N 是该目标值下的结构数量
- 目标值越小越接近理想随机

### 9.3 与 DFT 工作流对接

```bash
# 导出为 VASP POSCAR 格式
sqsgen output structure -f poscar

# 导出为 CIF（用于可视化）
sqsgen output structure -f cif
```

---

## 10. 实战示例

### 10.1 BCC Re-W 等摩尔合金（54 原子）

```json
{
    "structure": {
        "lattice": [[3.165, 0.0, 0.0], [0.0, 3.165, 0.0], [0.0, 0.0, 3.165]],
        "coords": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        "species": ["W", "Re"],
        "supercell": [3, 3, 3]
    },
    "iterations": 50000000,
    "shell_weights": {"1": 1.0},
    "composition": {"W": 27, "Re": 27},
    "max_results_per_objective": 10
}
```

**解读：**
- B2 结构，$a = 3.165$ Angstrom
- $3 \times 3 \times 3$ 超胞 = 54 原子
- 5000 万次随机尝试
- 仅优化第一配位壳层
- 27 W + 27 Re = 50% 等摩尔

### 10.2 高熵氧化物（含空位）

```json
{
  "iterations": 1000000,
  "sublattice_mode": "split",
  "structure": {
    "file": "Co3O4.pymatgen.json"
  },
  "composition": [
    {"sites": "Co", "Cr": 4, "Mn": 5, "Co": 5, "Ni": 5, "Fe": 5},
    {"sites": "O", "0": 8, "O": 24}
  ]
}
```

**解读：**
- 从 Co3O4 尖晶石结构出发
- Co 子晶格：5 种元素（高熵）
- O 子晶格：引入 8 个空位（O_0.75）
- 两个子晶格独立优化

### 10.3 三元合金 Ti-Al-Mo（54 原子）

```json
{
  "structure": {
    "lattice": [[3.2, 0.0, 0.0], [0.0, 3.2, 0.0], [0.0, 0.0, 3.2]],
    "coords": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    "species": ["Ti", "Al"],
    "supercell": [3, 3, 3]
  },
  "composition": {"Ti": 18, "Al": 18, "Mo": 18}
}
```

### 10.4 系统性搜索（小体系）

```json
{
  "iteration_mode": "systematic",
  "structure": {
    "lattice": [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
    "coords": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    "species": ["A", "B"],
    "supercell": [2, 2, 2]
  },
  "composition": {"A": 8, "B": 8}
}
```

16 原子、2 种元素 -> 16!/(8!*8!) = 12870 种排列，可完全穷举。

---

## 11. 模板系统

sqsgenerator 内置了多个示例模板：

```bash
# 列出所有模板
sqsgen template list

# 生成模板输入文件
sqsgen template use re-w.first
# 生成 re-w.first.json

# 直接运行模板
sqsgen run -i re-w.first.json
```

模板是学习和测试的良好起点。

---

## 12. 高级用法与技巧

### 12.1 可重复性

```json
{
  "seed": 42,
  "thread_config": 1
}
```

**注意：** 设置 seed 时必须使用单线程。

### 12.2 性能优化

- 减少 `shell_weights` 中的壳层数（只用第一壳层最快）
- 增大 `iterations` 提高找到最优解的概率
- 使用多线程（默认启用）
- HPC 环境使用原生 MPI 版本

### 12.3 自定义配位壳层

```json
{
  "shell_radii": [3.2, 4.5, 5.8]
}
```

手动指定壳层半径（Angstrom），覆盖自动检测。

### 12.4 非零目标 SRO（部分有序）

设置 `target_objective` 为非零值，生成具有特定短程序的结构：
- 负值 -> 模拟有序化倾向
- 正值 -> 模拟聚团倾向

### 12.5 与 ASE/pymatgen 集成

```python
from sqsgenerator import parse_config, optimize, write
from ase.io import read

# 从 ASE 读取结构
atoms = read("my_structure.cif")

# 构建配置
config = parse_config({
    "structure": {
        "file": "my_structure.cif",
        "reader": "ase",
        "supercell": [2, 2, 2]
    },
    "composition": {"Cu": 32, "Au": 32},
    "iterations": 10000000
})

pack = optimize(config)
best = pack.best()

# 导出为 ASE Atoms 对象
write(best.structure(), "best_sqs.ase.extxyz")
```

### 12.6 从 MD 轨迹读取结构

```json
{
  "structure": {
    "file": "md.traj",
    "reader": "ase",
    "args": {"index": -1}
  }
}
```

---

## 13. 常见问题

### Q: 如何判断 SQS 质量？

查看输出中的最小目标函数值：
- $\mathcal{O} = 0$：完美 SQS（所有优化壳层的 SRO 均为 0）
- $\mathcal{O} > 0$：近似 SQS，值越小越好
- 增大 `iterations` 或增大超胞可改善质量

### Q: random vs systematic 如何选择？

- 原子数 > 20：必须用 `random`（排列数太大）
- 原子数 < 16：可用 `systematic` 保证全局最优
- 用 `sqsgen compute total-permutations` 检查排列数

### Q: 如何只优化特定子晶格？

使用 `split` 模式 + `composition` 中的 `sites` 键：
```json
{
  "sublattice_mode": "split",
  "composition": [{"sites": "N", "B": 16, "N": 16}]
}
```

### Q: 如何表示空位？

使用 `"0"` 作为物种名：
```json
{"composition": {"Al": 56, "0": 8}}
```

### Q: 输出文件太大怎么办？

减小 `max_results_per_objective`（默认保存所有等价结构）。

### Q: 如何与 VASP 配合使用？

```bash
# 导出为 POSCAR
sqsgen output structure -f poscar
# 文件名为 sqs-0-0.poscar，可直接用于 VASP 计算
```

### Q: 多线程下结果不可重复？

这是预期行为。需要可重复性时设置 `seed` 并将 `thread_config` 设为 1。

---

## 参考文献

1. D. Gehringer, M. Friak, D. Holec, "Models of configurationally-complex alloys made simple," Computer Physics Communications, 108664 (2023).

2. A. Zunger, S.-H. Wei, L.G. Ferreira, J.E. Bernard, "Special quasirandom structures," Physical Review Letters 65, 353 (1990).

3. J.M. Cowley, "An approximate theory of order in alloys," Physical Review 77, 669 (1950).

4. sqsgenerator 官方文档：https://sqsgenerator.readthedocs.io/en/latest/

5. GitHub 仓库：https://github.com/dgehringer/sqsgenerator

---

*本教程基于 sqsgenerator 0.5.6 官方文档整理。具体参数和命令请以 `sqsgen --help` 输出为准。*