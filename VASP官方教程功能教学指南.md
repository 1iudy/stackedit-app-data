# VASP 官方教程功能教学指南

> 本文档整理自 VASP 官方教程 [https://vasp.at/tutorials/latest/](https://vasp.at/tutorials/latest/)（抓取日期：2026-08-05，共 58 个教程页面、16 个主题分类、37 个 Part、80 余个示例），按官网顺序逐节讲解每个教程所教授的 VASP 功能。
>
> 说明：
> - 按用户要求**保留全部图片**（教程输出图以本地 `images/` 相对路径或 vasp.at 原始 URL 嵌入）；
> - **示例输入文件（POSCAR/INCAR/KPOINTS/POTCAR 等）的具体内容已省略**，仅保留输入要点与关键标签说明；
> - 每个 Part 均附官方教程链接与输入文件 zip 下载地址；
> - 数学公式采用 LaTeX 记法，行内 $...$、行间 $$...$$。

---

## 目录

- 入门指南（Get started）
- 1. 原子与分子（Atoms and Molecules）
- 2. 块体体系（Bulk Systems）
- 3. 磁性（Magnetism）
- 4. 分子动力学（Molecular Dynamics）
- 5. 机器学习力场（Machine-Learned Force Fields）
- 6. 表面（Surfaces）
- 7. 过渡态（Transition States）
- 8. 杂化泛函（Hybrid Functionals）
- 9. 线性响应（Linear Response）
- 10. GW 近似（GW Approximation）
- 11. Bethe–Salpeter 方程（BSE）
- 12. X 射线吸收谱（XAS）
- 13. 强关联体系（Strongly Correlated Systems）
- 14. 核磁共振（NMR）
- 15. 声子（Phonons）
- 16. 电子–声子相互作用（Electron-Phonon Interactions）
- 附录 A：教程输入文件下载清单（41 个 zip）
- 附录 B：全部教程页面索引（58 页）

---

## 入门指南（Get started）

> 来源：<https://vasp.at/tutorials/latest/>

### 1. 快速开始

- 把 VASP 可执行文件加入系统路径，并设定教程工作目录：`export PATH=$PATH:path/to/vasp.X.X.X/bin`、`export TUTORIALS="path/to/tutorials/working/directory"`；
- 新手建议从**原子与分子**（Atoms and molecules）分类开始。

### 2. 搭建工作环境

- **安装 VASP**：按 VASP Wiki 的《Installing VASP 6.X.X》指南安装；建议尽量启用可选特性（hdf5 支持、链接 Wannier90 等）；
- **下载教程**：每个教程开头提供 zip 下载（本页对应 `get_started.zip`）；远程机器可用 `curl -O` 或 `wget`；
- **JupyterLab**：教程以 Jupyter notebook 形式给出，推荐用 JupyterLab 打开（终端 + notebook + 浏览器一体）；建议水平分屏并排查看终端与 notebook；
- **py4vasp**：读取 VASP 输出的便捷 Python 接口，多数教程至少部分依赖它，推荐安装。

### 3. 入门技巧

- **学习顺序**：按左侧导航自上而下；初学者从《Atoms and molecules - Part 1: Introduction to VASP》开始，涵盖 DFT 基础与离子弛豫；理论进阶可看 VASP YouTube 频道或参加 workshop；
- **终端命令**：熟悉 `ls`、`cd`、`rm`、`cp`、`pwd`、`bash`，用 `man cmd` 查手册，`mpirun --help` 看并行运行帮助；`tail -f`、`cat`、`vim` 也很常用；
- **vasp_rm**：建议写一个 bash 脚本 `vasp_rm` 加入路径，用于清理 VASP 输出文件（CHG、CHGCAR、CONTCAR、DOSCAR、EIGENVAL、IBZKPT、OSZICAR、OUTCAR、PROCAR、WAVECAR、XDATCAR、vasprun.xml、wannier90.* 系列、vaspout.h5、ML_* 系列等），随计算类型增删清单。

### 4. 获取帮助

- [VASP Wiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)：完整手册与标签参考；
- [VASP 论坛](https://www.vasp.at/forum/)：问题讨论；
- [VASP YouTube 频道](https://www.youtube.com/channel/UCBATkNZ7pkAXU9tx7GVhlaw)：理论讲解视频。

---

## 1. 原子与分子（Atoms and Molecules）

孤立的原子与分子是运行第一批 VASP 计算、熟悉基本输入输出的良好起点。每个示例开头列出学习目标，结尾附有自测问题。

## 1.1 Part 1：VASP 入门（孤立氧原子）

> 来源：https://vasp.at/tutorials/latest/molecules/part1/ ｜ 练习文件：`molecules-part1.zip` ｜ 示例目录：e01_O-DFT、e02_O-SDFT、e03_O-SDFT-symm

本部分以"盒子中的孤立氧原子"为载体，完成 VASP 的第一次完整计算。三个示例层层递进：先跑通一次最基本的 DFT 计算并读懂输出；再打开自旋极化得到氧原子的磁矩；最后改变盒子对称性，理解对称性与展宽参数对能量可比性的影响。

### 示例 1：孤立氧原子的 DFT 计算（e01_O-DFT）

**教学目标**：能用伪代码级别的描述解释 DFT 自洽流程；会为孤立原子准备四个输入文件；能识别 stdout 与 OUTCAR 的基本结构；会提取分子/原子计算的相关能量；会从已有 Kohn-Sham 轨道重启计算。

**DFT 自洽循环（伪代码）**：
1. 给定电子电荷密度，构造哈密顿量；
2. 求解哈密顿量的本征函数与本征值；
3. 由本征态更新电荷密度；
4. 重复 1–3 直至收敛。

**输入文件要点**：
- `POSCAR`：8 Å 立方盒子（足够大以忽略周期镜像间相互作用），1 个 O 原子置于原点，`cart` 笛卡尔坐标。
- `INCAR`：仅两个标签——`SYSTEM`（作业名注释）与 `ISMEAR = 0`（Gaussian 展宽；配合 `SIGMA` 控制展宽宽度，用于处理能级占据）。
- `KPOINTS`：Monkhorst-Pack `1 1 1`，仅 Γ 点。孤立体系的波函数在实空间局域，不存在色散，一个 k 点即足够。
- `POTCAR`：O 的 PAW 赝势（由 POTCAR 库拼接而来，包含赝势、原子波函数、电荷密度等预计算信息）。

**运行与输出解读**：`mpirun -np N vasp_std` 运行。stdout 中逐行打印电子迭代（`DAV:` 表示 Davidson 算法、`RMM:` 表示 RMM-DIIS），每步给出总能 E、能量变化 dE、本征值变化 d eps、电荷密度残差 rms 等；收敛后给出 `F=`（自由能）、`E0=`（无熵总能）、`d E=`。OUTCAR 中可找到：本征值表（band energies 与 occupation）、费米能、以及 `energy without entropy`（总能）。

**核心结论（教程问答）**：
- 孤立氧原子的总能 `energy without entropy` 接近零并非巧合：所有 PAW 赝势都以"孤立、非自旋极化原子"为参考态生成，此计算相当于把参考系与自身比较；盒子更大、精度更高时能量会严格趋于 0。牢记：绝对总能无物理意义，只有相对能量（两相关体系之差、相对费米能的带能量等）可解释。
- 重启计算：默认从 `WAVECAR` 读取上一步 KS 轨道继续，收敛显著加快；全新计算应先删除 WAVECAR。

**复习问题**：`ISMEAR = 0` 的含义及其与 `SIGMA` 的关系；孤立原子为何只需 1 个 k 点；本例电荷密度更新了多少次；总能、本征值、费米能在 OUTCAR 何处；重启需要哪个文件。

### 示例 2：自旋极化氧原子（e02_O-SDFT）

**教学目标**：开启自旋极化（SDFT）；从两个自旋分量的本征值与占据数推断自旋磁化强度；知道何时使用 `vasp_gam`。

**原理**：SDFT 中哈密顿量是 2×2 矩阵，KS 轨道为双分量旋量，两个分量对应自旋上/下，允许有不同的本征值；两自旋分量通过同时包含两种自旋密度的有效哈密顿量分别求解、自洽耦合。

**输入要点**：与示例 1 相同的 POSCAR/KPOINTS，INCAR 增加 `ISPIN = 2`（自旋极化）。因 KPOINTS 只含 Γ 点，可用 `vasp_gam` 可执行文件：它把部分数组按实数而非复数处理，计算更省；但仅当确实只用 Γ 点时才成立。

**输出解读与物理分析**：
- stdout 出现关于 `MAGMOM` 的警告：磁性/非共线计算未显式给 `MAGMOM` 时默认每原子 1 μB，这种铁磁初设可能排除反铁磁解，教程建议手动设置。警告不代表计算错误，而是提醒检查设置是否符合预期。
- OUTCAR 本征值：自旋分量 1 的 2s 能级约 −10.08 eV、三个 2p 能级全占据（occupation = 1）；自旋分量 2 的 2p 能级约 −7.05 eV、每个占据 1/3（即第 4 个 2p 电子以自旋向下分摊在三个简并 p 轨道上）。两分量本征值不同，说明基态是自旋极化的。
- 这与洪特规则一致：先以同一自旋占据每个 2p 轨道，余下电子反向。自旋磁化强度 $\langle \mu_{s,z} \rangle = g_s \langle \sigma_z/2 \rangle \mu_B$（$g_s \approx 2$）；自旋向上的电子比向下多 2 个，故磁矩为 2 μB，与 stdout 收敛行 `mag = 1.9998` 一致。

**复习问题**：`ISPIN` 取何值关闭自旋极化；两自旋分量的带能量与费米能是否相同；如何由占据数计算磁化强度。

### 示例 3：低对称性盒子中的自旋极化氧原子（e03_O-SDFT-symm）

**教学目标**：找到施加在电荷密度上的对称性；理解对称性如何从输入中推断；判断两次计算的总能能否相互比较。

**原理与背景**：晶体有 32 种晶体学点群。VASP 默认（`ISYM` 默认值）搜索并利用体系的全部对称性以降低计算量，但必须警惕：不要强加基态本不具备的对称性。周期边界条件本身就可能对解强加对称性，必要时需手动关闭（参见 `ISYM`）。

**输入要点**：POSCAR 改为三个互不相等的晶格矢量（7.0/7.5/8.0 Å，夹角均 90°），即正交晶系（D$_{2h}$）；其余与示例 2 相同。对比：示例 2 的立方盒子给出 O$_h$ 对称性。若要做四方对称（D$_{4h}$），令两个晶格矢量等长即可。注意此计算不含自旋轨道耦合——SOC 是相对论修正，需用 `LSORBIT` 显式打开；`ISPIN = 2` 只是非相对论哈密顿量下分自旋分量自洽。

**输出解读与物理分析**：
- OUTCAR 的对称性分析：`Found a simple orthorhombic cell`，8 个空间群操作（均为纯点群操作），`point symmetry D_2h`。
- 能量对比：正交（低对称）计算 `energy without entropy` ≈ −1.515 eV，低于立方计算的≈0，因此低对称解更接近真实多体基态。GGA 交换关联泛函对多数原子偏好对称性破缺的解。
- `SIGMA` 实验：把 `SIGMA` 从默认 0.2 改为 0.01 后，收敛所需迭代次数增多，且 `energy without entropy` 变为 ≈ −1.891 eV，数值显著变化。结论：要比较两次计算的能量，必须使用相同的 `SIGMA`（展宽方案一致），否则熵项贡献不同、能量不可比。

**复习问题**：正交情形 OUTCAR 报告多少空间群操作；如何改成四方对称；改变对称性的同时改变 `SIGMA` 时能量还能否比较。

**本部分涉及的 VASP 功能与标签汇总**：四输入文件（POSCAR/INCAR/KPOINTS/POTCAR）、DFT 自洽循环与 Davidson（`DAV`）/RMM-DIIS（`RMM`）电子最小化、`ISMEAR = 0` Gaussian 展宽与 `SIGMA`、孤立体系 Γ 点采样、`vasp_gam` 与 `vasp_std` 的区别、stdout/OUTCAR 结构（本征值、费米能、`energy without entropy`）、WAVECAR 重启、`ISPIN = 2` 自旋极化 SDFT、`MAGMOM` 初设与警告、洪特规则与自旋磁化强度、`ISYM` 对称性查找与 D$_{2h}$/O$_h$ 点群、`LSORBIT`（自旋轨道耦合）、能量比较的前提（同 `SIGMA`、同泛函、同赝势）。

---

## 1.2 Part 2：VASP 中的分子（O₂、CO 键长、振动与 DOS）

> 来源：https://vasp.at/tutorials/latest/molecules/part2/ ｜ 练习文件：`molecules-part2.zip` ｜ 示例目录：e04_O2-bond、e05_CO-bond、e06_CO-vibration、e07_CO-partial-dos

本部分从原子走向双原子分子，覆盖几何弛豫（共轭梯度算法）、多元素赝势拼接、有限差分振动频率计算与分波态密度（PDOS）四个主题。

### 示例 4：O₂ 键长的共轭梯度弛豫（e04_O2-bond）

**教学目标**：用共轭梯度（CG）算法做几何弛豫求二聚体键长；在伪代码层面理解几何弛豫（Hellmann-Feynman 力与应力）与 CG 算法。

**原理**：描述两个氧原子成键的完整多体哈密顿量包含电子与离子自由度；离子远重于电子、运动更慢，二者可解耦（Born-Oppenheimer 近似）。几何弛豫即寻找力与应力为零的离子位置，伪代码为：
1. 给定离子位置，用 DFT 计算电子基态；
2. 用 Hellmann-Feynman 定理计算作用在离子上的力与应力（电子哈密顿量梯度的期望值）；
3. 把离子弛豫到瞬时基态；
4. 重复 1–3 直至满足收敛判据。

**输入要点**：POSCAR 为 8 Å 盒子中相距 1.22 Å 的两个 O 原子；INCAR 在 SDFT 设置（`ISMEAR = 0`、`ISPIN = 2`）基础上增加 `NSW = 5`（最多 5 个离子步）与 `IBRION = 2`（CG 算法）。

**CG 算法（伪代码）**：
1. 最陡下降步：取搜索方向为最大梯度方向，沿该方向做线搜索最小化（Brent 方法，步长由 `POTIM` 控制）直至该方向力足够小，更新离子位置；
2. 共轭梯度步：把新的最大梯度对前一步搜索方向做正交化得到共轭方向，再做线搜索最小化并更新位置；
3. 重复以上直至梯度足够小或达到 `NSW`。

CG 适用于自由度少（≲4）且初猜接近基态的体系；否则应换用 `IBRION` 的其他算法。

**输出与分析**：从 stdout 读出离子步与电子步数；最终结构在 `CONTCAR`（弛豫后的键长约 1.23 Å，可由 Direct 坐标差×晶格常数验证）。也可用 py4vasp：`py4vasp.Calculation.from_path(...)` 载入后 `structure.print()` 查看。

**复习问题**：CG 方法如何工作、何时失效；`POTIM` 与 `NSW` 各控制什么。

### 示例 5：CO 键长——多元素体系与赝势选择（e05_CO-bond）

**教学目标**：运行含多种原子物种的 VASP 计算；选择合适的赝势；调节 CG 步长；用 py4vasp 可视化结构。

**赝势背景**：赝势是对离子实真实库仑势的近似，避免核处发散势与核间平缓势的强烈对比。赝势按元素区分，还按价电子数目与"软/硬"程度区分：含短键二聚体（O₂、CO、N₂、F₂、P₂、S₂、Cl₂）的体系推荐**硬赝势**（如 C_h、O_h），其代价是需要更大的平面波基组、计算更贵；低精度场合可用标准赝势。

**输入要点**：
- POSCAR：两种物种各 1 个原子（`1 1`），开启选择性自由度（`sel`）行，每个原子后跟 `F F T`——只允许沿 z 方向（键轴）弛豫，防止分子整体漂移。
- POTCAR 缺失，需自行拼接：`cat pot/C_h/POTCAR pot/O_h/POTCAR > POTCAR`。拼接顺序必须与 POSCAR 中物种顺序一致。
- INCAR：`ISMEAR = 0`、`NSW = 5`、`IBRION = 2`。

**练习与要点**：删除 WAVECAR 后在 INCAR 追加 `POTIM = 0.2`（小于默认 0.5）重算，观察更小步长对线搜索的影响；用 `my_calc.structure[:].plot()` 动画演示各离子步结构。

**复习问题**：POTCAR 拼接顺序；改变 CG 步长的标签；何时应选硬赝势、为什么。

### 示例 6：CO 振动频率——有限差分与 Hessian 矩阵（e06_CO-vibration）

**教学目标**：用 VASP 计算分子振动频率；解释 Hessian 矩阵与声子频率的联系。

**原理**：真实分子键长在约 10–100 THz 频率下振动。分子振动模式的计算与 Γ 点声子同理：需要对 Born-Oppenheimer 能量面 $E(\mathbf{R})$（以离子位置为参数的电子基态能量）求一阶与二阶导数。振动频率由离子 Hessian 的质量加权本征值问题给出：

$$\det \left| \frac{1}{\sqrt{M_\mu M_\nu}} \frac{\partial^2 E(\mathbf{R})}{\partial R_{\mu i} \partial R_{\nu j}} - \omega^2 \right| = 0$$

其中 $M_\mu$ 为离子质量，二阶导数即 Hessian 矩阵。

**输入要点**：INCAR 设 `IBRION = 5`（有限差分法计算二阶导数/Hessian/声子频率）、`NFREE = 2`（中心差分）、`POTIM = 0.02`（位移步长 0.02 Å）、`NSW = 1`（触发离子步）、`ISMEAR = 0`；POTCAR 用 C、O 硬赝势；POSCAR 沿用 CO 构型（键长 1.143 Å，选择性自由度只放开 z）。

**输出与分析**：OUTCAR 末尾 `Eigenvectors and eigenvalues of the dynamical matrix` 给出振动本征值（THz、2PiTHz、cm$^{-1}$、meV 四种单位）与本征矢。CO 伸缩模计算值约 63.89 THz，与实验值 2143 cm$^{-1}$（≈64.25 THz，Mina-Camilde 等, J. Chem. Educ. 1996）高度接近。

**复习问题**：Hessian 矩阵的定义与有限差分构造方式；`POTIM` 的作用；计算振动频率必须设置哪些标签；频率以哪些单位输出。

### 示例 7：CO 的分波态密度（e07_CO-partial-dos）

**教学目标**：用 py4vasp 绘制 DOS；区分分波态密度（partial DOS）、局域态密度（local DOS）与总态密度（total DOS）。

**原理**：DOS 描述某能量区间内电子态的数目 $N = \int f(\epsilon) D(\epsilon)\, d\epsilon$（$f$ 为 Fermi 函数）；对有限温度的相互作用电子，它对应谱函数，可与光电子能谱（PES）等实验对比。三者区分：
- **total DOS**：整个体系（此处为整个 CO 分子）的态密度；
- **partial DOS**：投影到特定轨道特征（如 O 的 2p 轨道）；
- **local DOS**：对晶体（波函数依赖波矢 k）而言，把 k 依赖的 DOS 对 $\mathrm{d}\mathbf{k}$ 积分得到。

**输入要点**：INCAR 仅加 `LORBIT = 11`——该设置使 VASP 把投影信息写入 DOSCAR/PROCAR，是获得 PDOS 的关键标签。注意：纯静态 DFT 计算中 POSCAR 的选择性自由度标记不起作用。

**输出与分析**：用 py4vasp 绘图：`my_calc.dos.plot()` 画总 DOS，`my_calc.dos.plot("p")` 按 p 轨道特征画分波 DOS。

**复习问题**：`LORBIT = 11` 产生什么输出；如何绘制分波 DOS。

**本部分涉及的 VASP 功能与标签汇总**：Born-Oppenheimer 近似与几何弛豫框架、Hellmann-Feynman 力、`IBRION = 2` 共轭梯度弛豫（线搜索、`POTIM` 步长）、`NSW` 离子步数、CONTCAR 末态结构、选择性动力学（`sel` + `F F T`）、多物种 POTCAR 拼接顺序、硬赝势（C_h/O_h）适用场景（短键二聚体）、`IBRION = 5/6` 有限差分 Hessian 与振动频率（`NFREE`、`POTIM`）、振动频率四种单位（THz/2PiTHz/cm$^{-1}$/meV）、`LORBIT = 11` 与 DOSCAR/PROCAR 投影、total/partial/local DOS 概念、py4vasp 结构与 DOS 可视化。

---

## 1.3 Part 3：水分子（RMM-DIIS、截断能收敛、DFPT 振动与从头算分子动力学）

> 来源：https://vasp.at/tutorials/latest/molecules/part3/ ｜ 练习文件：`molecules-part3.zip` ｜ 示例目录：e08_H2O-relaxation、e09_H2O-energy-cutoff、e10_H2O-vibration、e11_H2O-MD

本部分以水分子为对象，覆盖四个进阶主题：RMM-DIIS 离子弛豫、平面波截断能收敛研究、有限差分与 DFPT 两种振动频率算法的对比、以及从头算分子动力学与对关联函数。

### 示例 8：H₂O 键长的 RMM-DIIS 弛豫（e08_H2O-relaxation）

**教学目标**：在伪代码层面解释 RMM-DIIS；会判断何时用 RMM-DIIS、何时用 CG；从零写 POSCAR（含缩放参数）；做两自由度的几何弛豫；设置离子弛豫收敛判据。

**原理**：RMM-DIIS（残差最小化方法 + 迭代子空间直接求逆，又称广义最小残差法 GMRES）是一种准牛顿算法。与 CG 不同，它最小化的对象是离子受力模方 $|\mathbf{F}|^2$ 而非总能。对可微势 $V(\mathbf{r};\mathbf{R})$，力与 Hessian 由 Hellmann-Feynman 定理给出：

$$F_{\mu i} = \int \frac{\partial V(\mathbf{r};\mathbf{R})}{\partial R_{\mu i}} n(\mathbf{r};\mathbf{R})\, d\mathbf{r}, \qquad H_{\mu i \nu j} = \frac{\partial^2 E}{\partial R_{\mu i} \partial R_{\nu j}}$$

RMM-DIIS 按 Hessian 谱把最小化分解到子空间中求解，因此非常快，但对初始化敏感、偶有需要人工干预的情况。适用场景：离子位置已接近极小值、自由度 ≳3 的体系；自由度 ≳20 且振动谱很宽时建议改用分子动力学（如阻尼 MD/模拟退火）。

**输入要点**：
- POSCAR 从零构建：O 置于原点，H–O 键长约 0.96 Å、键角 105°，盒子足够大以模拟孤立分子。本例第二行用缩放参数 `0.52918`（Bohr→Å 换算因子），晶格矩阵与坐标均乘以该因子（12×0.52918 ≈ 6.35 Å 盒子）；开启 `select` 选择性动力学，O 固定（`F F F`），每个 H 放开两个自由度（`T T F`）。
- INCAR：`IBRION = 1`（RMM-DIIS），`EDIFFG` 设定离子弛豫收敛判据（力的阈值）。

**结果与分析**：弛豫后 H–O 键长约 0.971 Å、键角约 104.36°，与实验值接近。教程演示两种验证：由 CONTCAR 的 Direct 坐标用 $\cos\alpha = \mathbf{b}_1\cdot\mathbf{b}_2/(|\mathbf{b}_1||\mathbf{b}_2|)$ 手算；或用 py4vasp 画 2×2×2 超胞后在图上右键测量距离与角度。

**复习问题**：`EDIFFG` 的作用；RMM-DIIS 伪代码；POSCAR 中如何设置离子自由度；缩放参数的作用。

### 示例 9：H₂O 的平面波截断能收敛研究（e09_H2O-energy-cutoff）

**教学目标**：理解并设置平面波基组的能量截断 `ENCUT`；理解 `PREC` 精度档位。

**原理**：PAW 方法中 KS 轨道由离域的赝轨道加原子中心增补项构成。赝轨道用平面波展开：

$$\langle \mathbf{r} | \tilde{\psi}_{i\mathbf{k}\sigma} \rangle = \frac{1}{\sqrt{\Omega_r}} \sum_{\mathbf{G}} C_{i\mathbf{G}\mathbf{k}}^\sigma \, e^{i(\mathbf{G}+\mathbf{k})\cdot\mathbf{r}}$$

无限多平面波才是完备基，实际按能量截断舍去高频分量：$\frac{1}{2}|\mathbf{G}+\mathbf{k}|^2 < E_\text{cutoff}$。所需截断能主要取决于赝势（硬赝势需要更高截断），也依赖体系与目标性质，必须做收敛研究。

**输入与操作**：教程在多个子目录 `ENCUTXXX` 中放置相同输入，仅 `ENCUT` 从 400 eV 递增到 820 eV；INCAR 设 `PREC = Normal`、`ISMEAR = 0`、`SIGMA = 0.1`。用 bash 循环逐目录运行 `vasp_std`，再用提供的 `convergence_study.sh` + gnuplot 脚本提取各截断下的总能并绘图。

![H₂O 总能随平面波截断能的收敛曲线](/imgs/2026-08-05/fig1.png)

**核心结论**：总能随截断增大收敛向 PAW 框架内的完备基组（CBS）极限。但总能本身是最难收敛的量，实际工作应针对感兴趣的量（振动频率、化学位移、带隙等）做截断能收敛研究，而不是只看总能。精度从根本上受基组完备性限制，没有直接公式预测某性质所需的截断能——唯一可靠做法是对该性质逐点收敛。

**复习问题**：什么是赝轨道；PAW 基组是否完备、为什么；`PREC` 控制什么。

### 示例 10：H₂O 振动频率——有限差分 vs DFPT（e10_H2O-vibration）

**教学目标**：解释有限差分法与密度泛函微扰理论（DFPT）的区别；会设置能带数目 `NBANDS`。

**原理对比**：
- **有限差分**：把每个离子沿各方向做微小位移（步长 `POTIM`），用差分算符 $\Delta^\mathbf{R}$ 数值地构造电荷密度变化 $\Delta^\mathbf{R} n = 2\,\mathrm{Re}\sum_i^{occ} \psi_i^* \Delta^\mathbf{R}\psi_i$，进而得到力与 Hessian。
- **DFPT**：KS 轨道对离子位移的一阶响应由标准一阶微扰理论给出 Sternheimer 方程：

$$(H - \epsilon_{i\sigma}) |\Delta^\mathbf{R} \psi_{i\sigma}\rangle = -(\Delta^\mathbf{R} V - \Delta^\mathbf{R}\epsilon_{i\sigma}) |\psi_{i\sigma}\rangle$$

自洽求解后得到 $\Delta^\mathbf{R} n$，再算力与 Hessian。DFPT 的优势是只需要占据带（有限差分需要空带参与、`NBANDS` 要足够）；理论上 DFPT 还可解耦不同波长的响应、在任意 q 点算声子而无需超胞（这一点 VASP 未实现，但对孤立分子振动无影响）。

**输入要点**：两个子目录对比。共同设置：`PREC = Accurate`（算力推荐）、`ENCUT = 520`（依据示例 9 的收敛研究）、`ISMEAR = 0`、`EDIFF = 1E-8`（严格电子收敛）、`NSW = 1`。IBRION6 目录：`IBRION = 6`（利用对称性的有限差分）、`NFREE = 2`、`POTIM = 0.015`、`NBANDS = 6`；IBRION8 目录：`IBRION = 8`（DFPT），无需 NBANDS。

**结果与分析**：N 原子分子有 3N = 9 个模式：3 个平动 + 3 个转动（频率≈0）+ 3 个振动。以 IBRION=8、EDIFF=1E-8 为参考：对称伸缩约 115.5 THz（3853 cm$^{-1}$）、反对称伸缩约 112.3 THz（3744 cm$^{-1}$）、弯曲约 47.6 THz（1589 cm$^{-1}$）；平动/转动模式接近 0，个别出现微小虚频属数值噪声。练习：以 DFPT 为基准，改变 `POTIM` 考察有限差分频率对步长的敏感性。

**复习问题**：两种方法的本质区别；为什么 DFPT 不需要空带；`NBANDS` 的作用。

### 示例 11：H₂O 的从头算 MD 与对关联函数（e11_H2O-MD）

**教学目标**：关闭对称性做标准从头算 MD；理解 Nosé-Hoover 恒温器的基本概念；会分析对关联函数。

**背景**：MD 用牛顿第二定律模拟给定温度下的原子运动，每个离子步即一个时间步；力由量子力学计算时即"从头算 MD"。MD 也可解离子弛豫问题——缓慢降温的模拟退火（如 `IBRION = 3` 阻尼 MD 引入摩擦），尤其适合自由度 ≳20、振动谱很宽的体系；但水分子这种小体系用 RMM-DIIS 更高效。MD 相对弛豫算法的独特优势：轨迹本身具有物理意义，可用来计算对关联函数 g(r)——在距某粒子中心给定距离处找到另一粒子的概率。

**输入要点**：INCAR：`PREC = Normal`、`ENCUT = 520`、`ISMEAR = 0`、`SIGMA = 0.1`、`IBRION = 0`（MD）、`NSW = 500`、`POTIM = 0.5`（0.5 fs 时间步）、`ISYM = 0`（MD 必须关闭对称性）、`SMASS = 0`（Nosé 恒温器质量参数，控制热浴惯性/温度振荡周期）、`TEBEG = TEEND = 2000`（恒温 2000 K）、`NBANDS = 8`。

**输出与分析**：
- py4vasp `energy[:].plot()`：总能随时间振荡，振幅来自温度——Nosé-Hoover 恒温器按温度与时间步微扰离子位置，每次电子弛豫后给出新的总能；起止温度相同，故振荡分布全程一致。
- `structure[:].plot()` 动画显示分子以对称伸缩为主的振动。
- 对关联函数写入 PCDAT，用脚本绘图：大距离处 g(r) 趋于体系密度，长程行为约从 6 Å 开始（盒子边长 0.52918×12 ≈ 6.35 Å；6 Å 以下仍有信号是因为相邻水分子的 H–O 距离小于盒边）；短距离处由分子内成对相互作用决定——约 1 Å 的峰对应 H–O 键（基态约 0.97 Å），约 1.53 Å 处对应 H–H 距离；振动使峰展宽，对称伸缩模（3853 cm$^{-1}$ → 1/3.853 ≈ 0.26 Å）导致 1.53 Å 峰出现约 0.26 Å 的劈裂，与图示一致。

![水分子 MD 的对关联函数](/imgs/2026-08-05/fig57.png)

**复习问题**：什么是对关联函数；MD 能否解离子弛豫问题、何时选择 MD；`SMASS` 控制什么；`TEBEG` 的单位。

**本部分涉及的 VASP 功能与标签汇总**：`IBRION = 1` RMM-DIIS（准牛顿、最小化 |F|²、按 Hessian 谱分子空间）、`EDIFFG` 离子收敛判据、POSCAR 缩放参数（Bohr/Å 换算）与选择性动力学、PAW 赝轨道平面波展开与 `ENCUT` 截断、`PREC` 精度档、收敛研究方法学（对目标性质收敛）、`IBRION = 6`（对称性有限差分）与 `IBRION = 8`（DFPT）、Sternheimer 方程、`NBANDS` 与空带需求、`EDIFF = 1E-8` 严格收敛、`IBRION = 0` + `MDALGO`/`SMASS` Nosé-Hoover 恒温 MD、`ISYM = 0`、`TEBEG`/`TEEND`、`POTIM` 时间步、PCDAT 对关联函数、py4vasp 能量/结构/DOS 分析。

---

## 2. 块体体系（Bulk Systems）

块体体系是不受表面影响的物质部分的模型。本教程聚焦具有小原胞、周期性边界条件的晶体建模，学习密度泛函理论计算与结构优化。

## 2.1 Part 1：硅——典型块体材料（晶格常数、DOS、能带）

> 来源：https://vasp.at/tutorials/latest/bulk/part1/ ｜ 练习文件：`bulk-part1.zip` ｜ 示例目录：e01_fcc-Si、e02_fcc-Si-DOS、e03_fcc-Si-band

本部分以 fcc 硅为载体，完成块体材料的三大基础任务：手动扫描晶格常数求平衡体积、计算态密度（DOS）、计算高对称路径能带结构。

**结构学背景**：三维有 14 种 Bravais 格子（7 个晶系 × 4 种心化类型）。原胞是周期体系的最小重复单元，由晶格矢量 $\mathbf{a}$、$\mathbf{b}$、$\mathbf{c}$ 张成。实空间格子的傅里叶变换给出倒空间（k 空间），第一布里渊区是倒空间中唯一定义的原胞；利用点群对称性可进一步约化为不可约布里渊区。高对称点/线（如 Γ 点）因对称性强制的简并而特别重要。

### 示例 1：fcc 硅的晶格常数（e01_fcc-Si）

**教学目标**：从晶格矩阵认出 fcc 结构；构造 Γ 中心 k 网格；通过在不同体积下手动跑 DFT 找晶格常数；会用 POSCAR 的通用缩放因子。

**结构要点**：fcc 常规单胞展现立方全部对称性（$a=b=c$，面心附加位点），但存在更小的原胞表示。对称选择的原胞晶格矢量为 $\mathbf{a} = \frac{a}{2}(\hat{x}+\hat{y})$、$\mathbf{b} = \frac{a}{2}(\hat{y}+\hat{z})$、$\mathbf{c} = \frac{a}{2}(\hat{z}+\hat{x})$。常规单胞含 4 个原子（8 顶点×1/8 + 6 面心×1/2），原胞只含 1 个原子。

![fcc 原胞示意图](/imgs/2026-08-05/fig2.png)

**输入要点**：
- POSCAR：第二行用占位符 `a` 作为通用缩放参数，三个晶格矢量写为 $(0.5,0.5,0)$、$(0,0.5,0.5)$、$(0.5,0,0.5)$（即 $\mathbf{b}/a$、$\mathbf{c}/a$、$\mathbf{a}/a$），把缩放参数解释为晶格常数即恢复上式；1 个 Si 原子在原点。
- INCAR：`ISTART = 0`（从头开始）、`ICHARG = 2`（初始电荷密度取原子电荷密度叠加）、`ENCUT = 240`（若不在 INCAR 中设置，VASP 默认取 POTCAR 中的 ENMAX，本例恰为 240 eV）、`ISMEAR = 0`、`SIGMA = 0.1`。
- KPOINTS：Monkhorst-Pack `11 11 11`——奇数网格保证 Γ 中心。

**计算流程**：把 `a` 依次替换为 3.5–4.3（步长 0.1）逐个运行，用脚本从 OSZICAR 提取自由能 `F=` 写入数据文件，再用 gnuplot 绘 E(a) 曲线。

![不同晶格常数下的总能曲线](/imgs/2026-08-05/fig3.png)

**结果**：在测试值中总能在 a = 3.9 Å 处最低；把最优 POSCAR 拷贝到下一示例继续使用。

**复习问题**：Γ 点的定义；写出另一种 fcc 原胞晶格矩阵；常规 fcc 单胞是否含非初基平移。

### 示例 2：fcc 硅的态密度（e02_fcc-Si-DOS）

**教学目标**：从输出中提取不可约 k 点数；估计点群对称性节省的计算量；理解 Monkhorst-Pack 方法的任务；会为 DOS 设置 `ISMEAR`；会用 py4vasp 画 DOS 并设置 `NEDOS`/`EMIN`/`EMAX`。

**对称性与 k 网格**：体系局域上服从 32 种晶体学点群之一，倒空间则由周期性边界条件隐含地赋予倒格子的点群对称性。VASP 默认利用这些对称性，但 KPOINTS 网格必须与 POSCAR 结构定义的倒格子点群对称性一致，才能正确约化。Monkhorst-Pack 方法据此找出等价 k 点、定义不可约 k 点集与权重（权重=该不可约点在全布里渊区中的等价点数），把对全部 k 点的求和换成带权求和。本例 15×15×15 = 3375 个 k 点约化为 **120 个不可约 k 点**（OUTCAR 的 IBZKPT 输出），节省约 28 倍；四面体积分还报告 20250 个四面体中 484 个不等价。

**输入要点**：INCAR 设 `ISMEAR = -5`（带 Blöchl 修正的四面体方法，DOS 计算首选；注意它对部分占据不满足变分性，会给出错误的力与应力，**不可用于金属的几何优化**）与 `LORBIT = 11`（写出 DOSCAR）。KPOINTS 加密到 `15 15 15`。

**高分辨 DOS 练习**：先用 py4vasp 画全范围 DOS（费米能移到 0），再从 OUTCAR 读费米能（本例约 9.91 eV），加 `ICHARG = 11`（读入 CHGCAR 固定密度，只做非自洽谱计算）、`NEDOS = 401`、`EMIN = -5`、`EMAX = 12`（注意 EMIN/EMAX 是绝对能量，不是相对费米能）重算，得到费米能附近的高分辨 DOS。

**复习问题**：为什么 `ISMEAR = -5` 不能用于几何弛豫；如何提高 DOS 分辨率；本例全布里渊区与不可约 k 点数各是多少、为何不同。

### 示例 3：fcc 硅的能带结构（e03_fcc-Si-band）

**教学目标**：沿倒空间高对称线获得能带结构；用 py4vasp 绘图。

**k 路径选择**：高对称点取决于空间群，可用 SeeK-path 等工具（上传 POSCAR 即得布里渊区与建议路径）。本例路径为 L–Γ–X–U 与 K–Γ。

![SeeK-path 给出的 fcc 布里渊区与高对称点](/imgs/2026-08-05/fig4.png)

**输入要点**：能带计算需要 5 个输入——常规四件套加上自洽计算的 CHGCAR。INCAR：`ICHARG = 11`（读 CHGCAR 并固定电荷密度，沿路径做非自洽计算）、`ENCUT = 240`、`ISMEAR = 0`、`SIGMA = 0.1`、`LORBIT = 11`。KPOINTS 用 line 模式：第三行首字母 `l`（line）开启路径模式，`10` 表示每段插入点数，`reciprocal` 表示用分数坐标，随后逐段给出端点与标签（L/G/G/X/X/U/K/G）。

**流程**：`cp ../e02_*/CHGCAR .` 后运行，`mycalc.band.plot()` 绘图。

**复习问题**：能带计算需要哪些输入文件；如何选择与设置 k 路径。

**本部分涉及的 VASP 功能与标签汇总**：Bravais 格子/原胞/布里渊区/不可约区概念、POSCAR 通用缩放因子与占位符技巧、`ISTART`/`ICHARG` 启动控制、`ENCUT` 与 POTCAR ENMAX 默认值、Monkhorst-Pack 网格与 Γ 中心化、对称性约化与 k 点权重、`ISMEAR = -5` 四面体+Blöchl 修正（仅限 DOS/能带，不用于弛豫）、`ISMEAR = 0` Gaussian 展宽、`LORBIT = 11` 与 DOSCAR、`NEDOS`/`EMIN`/`EMAX`、`ICHARG = 11` 非自洽计算、KPOINTS line 模式能带路径、SeeK-path 工具、OSZICAR 能量提取与 bash/gnuplot 批处理。

---

## 2.2 Part 2：更多硅——金刚石结构、全弛豫与 β-Sn（收敛研究与作业控制）

> 来源：https://vasp.at/tutorials/latest/bulk/part2/ ｜ 练习文件：`bulk-part2.zip` ｜ 示例目录：e04_cd-Si、e05_cd-Si-relaxation、e06_perturbed-cd-Si、e07_beta-tin-Si

本部分在 Part 1 基础上进阶：独立完成金刚石结构硅的全流程（体积弛豫+DOS+能带）、单次运行全弛豫并规避 Pulay 应力、受扰结构的离子弛豫与 STOPCAR 作业控制、β-Sn 相的 k 网格收敛研究。

### 示例 4：立方金刚石硅的晶格常数、DOS 与能带（e04_cd-Si）

**教学目标**：手动完成体积弛豫；以最少指导独立完成 DOS 与能带计算。

**方法学提醒**：本例利用"算法收敛到能量面附近局域极小"的特性。但要警惕：初参离基态太远时，计算可能收敛到亚稳态而非基态——宣称"第一性原理基态"前应与实验对照。这一特性的好处是可以比较不同晶型/磁构型的相对稳定性。

**输入要点**：POSCAR 用 fcc 原胞晶格矢量 + Direct 坐标的两个 Si 原子（±1/8 体对角线偏移构成金刚石双原子基元）；INCAR：`ENCUT = 540`、`ISMEAR = 0`、`SIGMA = 0.1`、`EDIFF = 1e-6`、`ISTART = 0`、`ICHARG = 2`；KPOINTS `11 11 11`。

**三步流程**：
1. **晶格常数**：在 [5.4, 5.6] Å 区间用脚本逐点计算，取能量最低的两个值缩小区间、取中点迭代，直到两位小数精度——最优值 5.46(8) Å。最后用最优值再跑一次写出 CHGCAR 并备份（`cp CHGCAR CHGCAR.bk`）。
2. **DOS**：换 `ISMEAR = -5` + `LORBIT = 11`，KPOINTS 加密到 15×15×15；DOS 对截断能不敏感，`ENCUT` 可降为 400。本例 cd-Si 无部分占据带，四面体法用于体积弛豫也无妨（金属则不行）。建议把 `vaspout.h5` 备份为 `backup.h5` 供 py4vasp 回溯。
3. **能带**：`ICHARG = 11` 读入 Step 1 备份的 CHGCAR，KPOINTS 改 line 模式自选高对称路径（本例 L–Γ–X–U、K–Γ），py4vasp `band.plot()` 绘图。

**复习问题**：`LORBIT = 11` 控制什么；比较 fcc 与金刚石结构硅的总能，哪个更稳定。

### 示例 5：单次运行弛豫体积、晶胞形状与离子位置（e05_cd-Si-relaxation）

**教学目标**：理解并规避 Pulay 应力；用一次 VASP 运行同时弛豫体积、晶胞形状与离子位置；掌握 `ISIF`。

**Pulay 应力原理**：平面波基组按 $\frac{1}{2}|\mathbf{G}+\mathbf{k}|^2 < E_\text{cutoff}$ 截断，而倒格矢 $\mathbf{G}$ 相对实空间晶胞定义——同样 ENCUT 下不同体积的两次计算实际上使用不同的基组。体积越小等效截断越高、能量越低，因此不完备基组对体积弛豫引入偏向小体积的非物理应力，即 Pulay 应力。规避方法：① 像示例 4 那样在每个体积从头开始计算（每个体积独立基组）；② 单次运行时把 ENCUT 设得足够高，并在体积显著变化后重启计算。

**输入要点**：POSCAR 故意把第二个 Si 原子的 z 坐标微扰到 0.130（偏离 0.125 的理想位置）；INCAR：`ENCUT = 500`（偏高以抑制 Pulay 应力）、`IBRION = 2`（CG）、`ISIF = 3`（弛豫离子位置+晶胞形状+体积）、`NSW = 5`（限制离子步数，防止体积变化过大后 Pulay 应力显著）、`EDIFFG = -0.001`（按力收敛，负值表示以力为判据，单位 eV/Å）、`EDIFF = 1E-06`。

**流程**：运行后 `cp CONTCAR POSCAR` 重启继续，直至体积变化在小数点后若干位才出现，即达收敛。

**复习问题**：Pulay 应力的来源；`ISIF` 各取值含义；为何限制 NSW。

### 示例 6：受扰金刚石硅的离子弛豫与作业控制（e06_perturbed-cd-Si）

**教学目标**：只弛豫离子位置（固定晶胞）；理解停止判据；会用 STOPCAR 优雅中止作业并续算。

**输入要点**：与示例 5 相同的受扰 POSCAR；INCAR：`ISTART = 0`、`ICHARG = 2`、`ENCUT = 360`、`IBRION = 2`、`ISIF = 2`（只弛豫离子位置）、`NSW = 10`、`EDIFFG = -0.0001`、`EDIFF = 1E-06`。

**要点与练习**：
- 首次运行因 `NSW = 10` 耗尽而停止，并未真正收敛——要时刻清楚计算为什么停下。
- 把 NSW 加到 15 重算，并在另一个终端写入 `echo LSTOP = .TRUE. > STOPCAR`，让 VASP 在完成当前离子步后优雅退出（输出文件完整可用，可直接续算）。
- 续算：删除 STOPCAR，`cp CONTCAR POSCAR`，把停止判据设为略好于已达到的水平（例如已跑 3 步、目标共 10 步，则设 NSW = 8）。
- 收敛后观察晶体对称性在多大程度上被恢复。

**复习问题**：`EDIFFG` 取负值的含义；STOPCAR 停止后的波函数/电荷密度能否使用；续算离子弛豫要用哪个文件替换 POSCAR。

### 示例 7：β-Sn 结构硅——k 网格收敛研究（e07_beta-tin-Si）

**教学目标**：选择合适的 k 网格密度；为单元素计算创建 POTCAR。

**任务**：教程提供的 β-Sn POSCAR（四方，a = b = 4.8 Å、c/a = 0.26 内参数，2 原子）体积不可信，需弛豫。用 `ISIF = 7`（固定形状只变体积）在不同 k 网格密度下做体积弛豫收敛研究，判据：第一晶格矢量长度精确到约 0.005 Å、总能精确到 2 meV；最后计算 DOS。

**流程要点**：
- POTCAR：`cat pot/Si/POTCAR > POTCAR`。
- k 循环脚本：对 k = 3–13，每个网格跑两轮弛豫（`cp CONTCAR POSCAR` 再跑），记录 OSZICAR 能量与 CONTCAR 第三行（缩放后晶格矢量），gnuplot 分别画能量-网格数与晶格长度-网格数曲线。
- 观察：增大 k 点数并不单调趋近某值，但方差随 k 增多而减小；Γ 中心与否会引入显著差异。
- DOS：从收敛结构出发，把网格加密到 18×18×18，分别用 Gaussian 展宽（`ISMEAR = 0`）与四面体法（`ISMEAR = -5`）计算，py4vasp 把两条 DOS 画在同一图（`gaussian.dos.plot().label("Gaussian") + tetrahedron.dos.plot().label("Tetrahedron")`）对比。

**复习问题**：总能随 k 点数是否单调；如何为你的计算确定合适的 k 点数。

**本部分涉及的 VASP 功能与标签汇总**：局域极小与亚稳态的方法学警示、体积弛豫二分迭代法、CHGCAR 备份与非自洽计算、Pulay 应力机理与规避（高 ENCUT、分体积独立计算、体积大变后重启）、`ISIF = 2/3/7`（离子/晶胞形状/体积的组合控制）、`EDIFFG` 负值力判据、STOPCAR（`LSTOP`）优雅中止与 CONTCAR→POSCAR 续算、k 网格收敛研究方法论、Γ 中心网格的影响、`vaspout.h5` 备份与 py4vasp 回溯、Gaussian vs 四面体 DOS 对比。

---

## 2.3 Part 3：自旋极化块体——fcc Ni 与反铁磁 NiO（泛函对比）

> 来源：https://vasp.at/tutorials/latest/bulk/part3/ ｜ 练习文件：`bulk-part3.zip` ｜ 示例目录：e08_fcc-Ni、e09_AFM-NiO

本部分进入磁性块体：先对铁磁 fcc Ni 做自旋极化的体积弛豫、DOS 与能带；再以反铁磁 NiO 为案例，横向对比 PBE、DFT+U、SCAN、MBJ、HSE06 五种泛函对强关联体系带隙与磁矩的描述能力。

### 示例 8：自旋极化 fcc 镍（e08_fcc-Ni）

**教学目标**：做共线自旋密度泛函（SDFT）计算；弛豫共线磁体的体积；解读磁性体系的 DOS；用 py4vasp 画磁性 DOS 与能带。

**原理**：SDFT 中哈密顿量是 2×2 矩阵，KS 轨道是自旋空间的双分量向量，两分量可有不同本征值；两自旋分量通过同时依赖两种自旋密度的交换关联泛函耦合。实用上 `ISPIN = 2` 给出共线磁矩；若自旋方向重要（各向异性等），需加 `LSORBIT = True` 并用 `vasp_ncl` 可执行文件。

**输入要点**（三套 INCAR/KPOINTS 分别对应体积弛豫、DOS、能带）：
- 体积弛豫：`ICHARG = 2`、`ISMEAR = 1`（Methfessel-Paxton，适合金属）、`SIGMA = 0.2`、`ENCUT = 540`、`IBRION = 2`、`ISIF = 7`、`NSW = 6`、`EDIFFG = -0.001`、`ISPIN = 2`、`MAGMOM = 1`（初始磁矩）；KPOINTS `9 9 9`。
- DOS：`ENCUT = 350`、`ISMEAR = -5`、`LORBIT = 11`、`ISPIN = 2`、`MAGMOM = 1`；KPOINTS `16 16 16`。
- 能带：`ICHARG = 11`（读备份 CHGCAR）、`ISMEAR = 1`、`SIGMA = 0.2`、`ENCUT = 540`；KPOINTS 为 L–Γ–X–U、K–Γ 线模式。

**流程与结果**：
1. 体积弛豫只跑最多 6 个离子步，防止体积变化过大后 Pulay 应力偏向小体积；反复 `cp CONTCAR POSCAR` 重启并用 `vimdiff CONTCAR POSCAR` 观察变化出现在哪一位小数。收敛晶格常数：CONTCAR 中晶格矢量分量 0.5025，故 a = 3.5×0.5025/0.5 ≈ 3.517 Å。
2. 备份 CHGCAR（`cp CHGCAR CHGCAR.bk`）供能带使用；用优化结构算 DOS，py4vasp 绘图——磁化使自旋上/下 DOS 分裂，可观察 3d 带的交换劈裂。
3. 能带计算读回备份 CHGCAR，`band.plot()` 绘图。

**复习问题**：自旋向上是否一定指向 +z；磁化对 DOS 的定性影响；Methfessel-Paxton 方法适用于什么任务、哪些体系不能用。

### 示例 9：NiO 的带隙与磁矩——泛函对比（e09_AFM-NiO）

**教学目标**：用 `XC` 标签选择交换关联泛函；用 `BANDGAP` 标签控制带隙输出详细度；设置与原子磁矩相关的 `MAGMOM`/`LORBIT`；用 py4vasp 提取带隙与原子磁矩。

**背景**：反铁磁 NiO 是强关联体系（Ni 3d 电子高度局域），标准 GGA（如 PBE）通常给出定性错误的电子与磁性。处理此类体系需要 DFT+U、meta-GGA（SCAN、MBJ）或杂化泛函（HSE06）等更高级方法（DMFT 等 beyond-DFT 方法不在本教程范围）。实验参考：带隙 4.0–4.3 eV，Ni 磁矩 1.9–2.2 μB。

**结构要点**：NiO 为岩盐结构，但为容纳 AFM II 反铁磁序（沿立方 [111] 方向堆叠的铁磁平面、磁矩逐层交替），采用菱方单胞。立方常规单胞含 8 个原子（1 + 8×1/8 + 6×1/2 + 12×1/4），菱方单胞含 4 个原子（岩盐原胞只有 2 个原子）。

![NiO 立方常规单胞与菱方原胞中的 AFM II 磁结构](/imgs/2026-08-05/fig5.png)![NiO 菱方原胞磁结构](/imgs/2026-08-05/fig6.png)

**输入要点**：POSCAR 为菱方单胞（缩放因子 4.17 Å，2 Ni + 2 O）。各泛函 INCAR 的公共设置：`LASPH = .TRUE.`（计入梯度校正的非球面贡献，对磁性/关联体系重要）、`ISMEAR = -5`、`ENCUT = 350`、`EDIFF = 1E-4`、`ISPIN = 2`、`MAGMOM = 2 -2 0 0`（指定 AFM II 初始磁构型——磁性体系中初猜强烈影响最终磁态，想要特定磁序必须显式设置）、`LORBIT = 11`（把原子磁矩写入 OUTCAR 与 vaspout.h5）、`BANDGAP = KPOINT`（打印更详细带隙信息）。泛函差异：PBE 用 `XC = PE`；DFT+U 用 `XC = CA`（LDA 基底）+ `LDAU = .TRUE.`、`LDAUTYPE = 2`、`LDAUL = 2 -1`（对 Ni 的 d 轨道加 U，O 不加）、`LDAUU = 8.00 0.00`、`LDAUJ = 0.95 0.00`；SCAN 用 `XC = SCAN`；MBJ 用 `XC = MBJ`；HSE06 用对应杂化设置。`XC` 标签自 VASP 6.4.3 引入，还可用于混合不同泛函。杂化泛函贵一到两个数量级，故 HSE06 用较粗 k 网格（带来约 0.2 eV 带隙误差）。

**结果与分析**：
- PBE：基本带隙 0.96 eV（py4vasp `mycalc.bandgap.fundamental()`），远小于实验 4.0–4.3 eV；最小直接带隙 1.13 eV 大于基本带隙 0.17 eV，故为间接带隙。`grep "fundamental gap" OUTCAR` 可见两自旋通道带隙相同——反铁磁体中自旋上/下带是简并的。原子磁矩（`mycalc.local_moment.projected_magnetic()`）：Ni ±1.37 μB、O 严格为 0（对称性要求），符号与 MAGMOM 设定一致；PBE 磁矩远小于实验值。
- DFT+U：带隙 3.44 eV，Ni 磁矩 1.74 μB。DFT+U 对局域 d/f 电子加入在位库仑/交换项，概念简单、成本与基底泛函相当，但 U 常需经验调节、削弱预测力。
- SCAN：带隙 2.58 eV，磁矩 1.61 μB。SCAN 是固态物理最流行的 meta-GGA，通过动能密度等引入密度二阶导数信息，满足精确泛函的多个解析约束且无经验参数。
- MBJ：带隙 4.61 eV，磁矩 1.75 μB。MBJ 势通常能相当准确地给出固体带隙，但它通过附加势而非能量泛函实现大带隙，没有对应的总能贡献，**不能用于结构优化或稳定性判断**。
- HSE06：带隙 4.79 eV，磁矩 1.68 μB。HSE06 是屏蔽杂化泛函（局域泛函+非局域 HF 交换的混合，屏蔽掉长程部分），k 网格收敛比 PBE0 快，但 HF 交换使计算贵一到两个数量级。
- 结论：带隙以 DFT+U、MBJ、HSE06 最接近实验（DFT+U 仍偏低，MBJ/HSE06 略偏高）；Ni 磁矩以 DFT+U 与 MBJ 最好（约 1.73 μB，仍低于实验区间）。注意计算磁矩只含自旋贡献，不含轨道贡献——NiO 的轨道磁矩估计为 0.3–0.45 μB，可弥合差距；磁矩数值还依赖投影球半径等数值细节。

**复习问题**：哪些泛函含可调参数；HSE06 与 MBJ 结果相近，哪个计算更高效。

**本部分涉及的 VASP 功能与标签汇总**：`ISPIN = 2` 共线 SDFT、`MAGMOM` 初始磁构型设定（铁磁 vs 反铁磁）、`LSORBIT`/`vasp_ncl`（自旋方向相关计算）、`ISMEAR = 1` Methfessel-Paxton（金属）、`ISIF = 7` 磁性体积弛豫与 Pulay 应力控制、`XC` 标签泛函选择（PE/CA/SCAN/MBJ/HSE06）、`LASPH` 非球面梯度贡献、DFT+U 全套标签（`LDAU`/`LDAUTYPE`/`LDAUL`/`LDAUU`/`LDAUJ`）、`BANDGAP = KPOINT` 带隙输出、`LORBIT = 11` 原子磁矩投影、py4vasp bandgap（fundamental/direct）与 local_moment 分析、间接/直接带隙判别、强关联体系泛函选择策略与实验对照。

---

## 2.4 Part 4：范德华力——石墨层间结合（Tkatchenko-Scheffler 方法）

> 来源：https://vasp.at/tutorials/latest/bulk/part4/ ｜ 练习文件：`bulk-part4.zip` ｜ 示例目录：e10_VdW-TS-binding、e11_interlayer-distance

本部分处理半局域 DFT 的著名短板——长程色散（范德华）相互作用，用 Tkatchenko-Scheffler（TS）方法计算石墨层间结合能与层间距。

### 示例 10：TS 方法计算石墨层间结合能（e10_VdW-TS-binding）

**教学目标**：在 DFT 中加入范德华相互作用；设置 `LWAVE`/`LCHARG`；计算两个结构的总能差。

**背景**：GGA 级半局域 DFT 低估长程相互作用。石墨的 PBE 层间结合能仅约 1 meV/原子，远低于 RPA 参考值 0.048 eV/原子（Lebègue 等, PRL 105, 196401, 2010）。TS 方法基于 Becke 1997 的幂级数形式，显式加入带阻尼的原子对色散修正 $C_6 R^{-6}$（实现与测试见 Bucko 等, PRB 87, 064110, 2013）。

**输入要点**：两个子目录 graphene（单层，c = 20 Å 真空）与 graphite（双层堆叠，c = 6.71 Å）：
- POSCAR：同样的面内晶格（a ≈ 2.456 Å 六角），graphene 2 原子、graphite 4 原子（两层各 2 原子，AB 堆叠坐标）。
- INCAR：`LWAVE = False`、`LCHARG = False`（只用总能，不写 WAVECAR/CHGCAR 省空间）、`ISMEAR = -5`、`SIGMA = 0.01`、`ALGO = Fast`、`PREC = Accurate`、`EDIFF = 1e-6`；范德华设置：`IVDW = 20`（TS 方法）、`LVDW_EWALD = True`（可选：对成对相互作用做 Ewald 求和，短程在实空间、长程在倒空间处理以避免奇异）。
- KPOINTS：graphene `15 15 1`（z 方向层间无相互作用，只需 1 个 k 点）；graphite `15 15 7`。

**计算与结果**：分别运行两个体系，结合能 = (石墨烯×层数 − 石墨) 总能/原子数。PBE 结果接近 0（几乎不结合），加 TS 后显著增大——但教程指出 TS 结合能相对 RPA 参考仍偏大。py4vasp 可用 `structure.plot(3)` 可视化 3×3×3 超胞。

**复习问题**：`LWAVE`/`LCHARG` 的作用；两结构能量差的物理意义。

### 示例 11：石墨的层间距（e11_interlayer-distance）

**教学目标**：用 TS 方法沿堆叠方向弛豫晶格常数（c 轴），理解 Ewald 求和的动机。

**背景**：GGA 低估长程色散导致沿堆叠方向晶格常数被高估：GGA-PBE 给出 c = 8.84 Å，而实验为 6.71 Å。TS 方法中的长程库仑部分可选用 Ewald 求和（实空间处理短程、倒空间处理长程，避免奇异）。

**输入要点**：POSCAR 中 c 用占位符 `c` 待替换；INCAR：`ISMEAR = 0`、`SIGMA = 0.1`、`ALGO = Fast`、`PREC = Accurate`、`ENCUT = 620`、`EDIFF = 1e-6`、`IVDW = 20`、`LVDW_EWALD = True`；KPOINTS `15 15 7`。

**计算流程**：固定其余自由度为实验值，在 [6.6, 6.8] Å 区间以 0.01 Å 步长扫描 c，用脚本从 OUTCAR 提取 `free energy` 写入 loop.dat，gnuplot 绘 E(c) 曲线找极小。

![石墨总能随层间距变化曲线](/imgs/2026-08-05/fig7.png)

**结果与讨论**：TS 方法给出层间距 6.65 Å，与实验 6.71 Å 吻合良好——尽管示例 10 中 TS 结合能相对 RPA 偏大：结合能（能量尺度）与平衡距离（几何）对色散修正的敏感度不同，几何量往往更稳健。思考题：如何改 KPOINTS 完全忽略层间相互作用（z 方向取 1 个 k 点+足够真空）；Ewald 求和的动机与基本思想。

**本部分涉及的 VASP 功能与标签汇总**：GGA 对长程色散的失效（结合能与晶格常数）、`IVDW = 20` Tkatchenko-Scheffler 方法、$C_6 R^{-6}$ 阻尼对势、`LVDW_EWALD` Ewald 求和（实空间短程/倒空间长程）、`LWAVE`/`LCHARG` 控制输出文件、RPA 参考值对比方法学、低维体系 k 网格各向异性设置（真空方向 1 点）、结合能与平衡几何对泛函的不同敏感度。

---

## 3. 磁性（Magnetism）

磁性源于电子自旋、量子统计与电子–电子相互作用，是一种集体的量子电动力学现象。教程介绍如何用自旋极化建模铁磁/反铁磁材料、DFT+U、非共线计算与磁各向异性。

## 3.1 Part 1：自旋极化计算（hcp Co、bcc Cr、NiO 的 eg/t2g、DFT+U 与赝势选择）

> 来源：https://vasp.at/tutorials/latest/magnetism/part1/ ｜ 练习文件：`magnetism-part1.zip` ｜ 示例目录：e01_hcp-Co、e02_bcc_Cr、e03_NiO、e04_NiO_LSDA+U

本部分覆盖自旋极化计算的四个典型场景：铁磁金属的交换劈裂、反铁磁体的磁单胞与传播矢量、晶体场分裂的轨道分辨 DOS、以及 DFT+U 参数与 PAW 赝势选择对结果的影响。

### 示例 1：铁磁 hcp Co 的 DOS 与能带（e01_hcp-Co）

**教学目标**：做自旋极化计算；获得总磁化与局域磁矩；画自旋极化 PDOS；自洽后第二步算能带。

**原理**：`ISPIN = 2` 时 KS 轨道是自旋空间双分量向量，对应全局自旋量子化轴的上/下分量，允许不同本征值，从而描述巡游电子的共线磁结构；两分量仅通过依赖两种自旋密度的交换关联泛函耦合。**不含 SOC 时自旋在实空间的方向是任意的**，只有相邻磁矩的相对方向有意义（可区分铁磁/反铁磁）。共线计算用 `vasp_std`；含 SOC 或非共线（局域自旋量子化轴）才需要 `vasp_ncl`。

**输入要点**：hcp Co 原胞（a ≈ 2.472 Å、c ≈ 4.021 Å，2 原子）。三套 INCAR：
- SCF：`ALGO = Normal`、`PREC = Accurate`、`EDIFF = 1e-5`、`ENCUT = 350`、`ISPIN = 2`、`MAGMOM = 2 2`（铁磁初设）、`LMAXMIX = 4`（磁性体系必须——让 d 通道的 PAW 占居矩阵进入混合，否则收敛慢甚至失败）、`LASPH = T`（非球面梯度校正，磁性体系推荐）、`ISMEAR = 2`（Methfessel-Paxton，金属）、`SIGMA = 0.1`。
- DOS：加 `ICHARG = 11`（非自洽）、`ISMEAR = -5`（四面体）、`NEDOS = 5001`、`LORBIT = 11`（局域磁矩+PDOS）。
- 能带（band 子目录）：`ICHARG = 11` + line 模式 KPOINTS（Γ–M–K–Γ）。

**结果与分析**：
- DOS 显示自旋上/下带在能量上相对移动——多数自旋与少数自旋之间的相对位移即交换劈裂（Stoner 模型语境）。`mycalc.dos.to_dict("d")` 可得能量网格、上/下分量与 d 投影、费米能。
- 局域磁矩：`mycalc.local_moment.to_dict()` 按 s/p/d 轨道给出每原子磁矩（本例 d 通道约 1.61 μB，合计约 1.53 μB/原子）。注意这只含自旋贡献不含轨道角动量贡献；d 轨道的轨道矩较小，但 f 轨道可能可观（轨道矩在 Part 2 处理）。计算值对应"有效磁矩"（μ子散射/居里定律可测），而非中子散射测的 $m = gm\mu_B$。
- 无 SOC 时自旋量子化轴不与实空间耦合，计算不含磁矩方向信息；两 Co 位点分磁化同号，说明铁磁耦合（与 MAGMOM 初设一致）。
- 能带：少数自旋带整体上移交换劈裂量——例如 M–K 间 −0.5 eV 附近的三条自旋上带，在自旋下通道约 1.3 eV 处可见相同色散。

**复习问题**：`ISPIN` 的作用；获得局域磁矩必须设哪个标签（`LORBIT`）；什么是交换劈裂、在哪里观察。

### 示例 2：反铁磁 bcc Cr 的磁单胞（e02_bcc_Cr）

**教学目标**：对反铁磁材料做自旋极化计算；重启自旋极化计算。

**原理——传播矢量与磁单胞**：传播矢量 $\mathbf{q}$ 描述磁矩在块体中的空间变化周期，可据以定义磁单胞（可能大于晶体学单胞）。bcc Cr 顺磁相原胞只有 1 个原子；AFM 相原胞含 2 个原子（(0,0,0) 与 (0.5,0.5,0.5)），磁矩方向交替，磁单胞是晶体学单胞的两倍。$\mathbf{q}=(1/2,1/2,1/2)$（倒格矢单位）表示磁矩以两个晶格常数为周期振荡，对应体对角线方向的单胞加倍。

**输入要点**（三套 INCAR）：
- SCF：`ISPIN = 2`、`MAGMOM = 2 -2`（AFM 初设）、`LMAXMIX = 4`、`LASPH = T`、`ISMEAR = 0`、`SIGMA = 0.1`、`LORBIT = 11`、`LKPOINTS_OPT = F`；KPOINTS 4×4×4。
- 细网格续算：`MAGMOM = -0.1 0.1`（只为破缺对称性）、`LKPOINTS_OPT = T`（同时按 KPOINTS_OPT 线模式算能带）、`ICHARG = 1`（从 WAVECAR/CHGCAR 重启读电荷密度）；KPOINTS 8×8×8。
- DOS 后处理：`ALGO = None`（纯后处理，不做电子最小化）、`ISMEAR = -5`、`ICHARG = 11`、`NEDOS = 5001`、`LH5 = T`。

**流程与要点**：
- 粗网格收敛后读 OUTCAR 的 `Total CPU time used` 记录成本，备份 OUTCAR；换细网格时**删除 WAVECAR**——k 点数改变后读它会失败，但失败的读取尝试仍有副作用：没有 WAVECAR 时 KS 轨道随机初始化、前几步固定电荷密度，而读取失败后电荷密度不固定，就吃不到粗网格已收敛的红利。
- 对比两个 OUTCAR（`grep LOOP`）：细网格明显更贵；**在位磁矩在粗网格上并未收敛**（1.64 → 1.87 μB），而在位电荷几乎不变——必须对目标量本身做收敛研究。
- `mycalc.band.plot("kpoints_opt")` 画 KPOINTS_OPT 路径能带（图例可点击开关）。
- AFM 的 DOS：位点 1 多数自旋 DOS 更大但**没有铁磁体那样的交换劈裂**；位点 2 因 AFM 序上/下互换。
- 平面波网格上的磁化密度是 SDFT 的直接副产品（局域磁矩则依赖 `LORBIT` 投影），可从 vaspwave.h5 读取并用 `density.to_contour("magnetization")` + `to_quiver()` 画箭头图。

**复习问题**：`MAGMOM` 是否影响电荷密度的对称化；`ALGO = None` 是什么；严格来说必须对总能做 k 网格收敛吗；k 网格改变后如何重启。

### 示例 3：反铁磁 NiO 中 eg 与 t2g 轨道的 DOS（e03_NiO）

**教学目标**：把费米能设在带隙中央；调节密度混合器加速收敛；在 py4vasp 中组合轨道画 PDOS。

**晶体场背景**：NiO 中 Ni²⁺ 被 6 个等价 O²⁻ 配位成 NiO₆ 八面体。按 O_h 点群特征标表，Ni d 轨道分裂为 $t_{2g}$（d_xy、d_xz、d_yz，三重简并，g=gerade 中心对称，2 表示对垂直主轴的镜面反对称）与 $e_g$（d_z²、d_x²−y²，二重简并）。$e_g$ 瓣指向 O²⁻ 配体、受负点电荷排斥更强，故八面体晶体场下 $e_g$ 能量更高。离子极限下 Ni²⁺ 为 [Ar]3d⁸：6 个电子填满 $t_{2g}$，2 个电子在 $e_g$ 上按洪特规则形成总自旋 1，预期最大磁化 2 μB/分子式。

![Ni d 轨道在八面体场中的 eg/t2g 分裂](/imgs/2026-08-05/fig58.jpg)

**输入要点**：AFM NiO 菱方单胞（同 bulk 教程 e09 结构）。INCAR 关键点：
- `ISPIN = 2`、`MAGMOM = 2.0 -2.0 2*0`、`LASPH = T`、`EDIFF = 1E-8`。
- `ISMEAR = -5` + `EFERMI = MIDGAP`：把费米能放在带隙中央（绝缘体 DOS 计算技巧）。
- 密度混合：`IMIX = 4`（Pulay 混合），电荷与磁化各有线性混合系数（`AMIX`/`BMIX`）与 Kerker 波矢截断参数（`AMIX_MAG`/`BMIX_MAG`）；本例用大 AMIX（0.8）+极小 BMIX（0.00001）的组合加速这个绝缘体的收敛。
- `LMAXMIX = 4`（d 通道进入混合）、`LORBIT = 11`、`EMIN = 3.5`、`EMAX = 10`、`NEDOS = 5001`。

**结果与分析**：VASP 直接投影的是五个 d 轨道（dz²、dxz、dyz、dx²−y²、dxy），$t_{2g}$/$e_g$ 是后处理组合：py4vasp 用 `dos.to_dict()` 取出各轨道分量后自行加减绘制（如 Ni1 eg up = dz2_up + dx2y2_up）。可见少数自旋在两个 Ni 位点间交替；每个位点上多数自旋几乎全占、少数自旋几乎全空。Ni 在位磁矩约 1.34 μB，低于实验约 1.9 μB（可参考 Kwon & Min, PRB 62, 73, 2000 的综述）。

**复习问题**：如何选线性密度混合；哪个标签把费米能设带隙中央；$t_{2g}$ 的 DOS 是 VASP 直接算的吗（否，需组合投影）。

### 示例 4：NiO 弛豫中的在位库仑作用与赝势选择（e04_NiO_LSDA+U）

**教学目标**：用状态方程拟合做体积弛豫；选择并对比不同 PAW 势；运行 Dudarev 等（PRB 57, 1505, 1998）方法的 LSDA+U 计算。

**背景**：DFT+U 针对强关联电子（过渡金属氧化物、f 电子体系）补足 LDA/GGA 对在位库仑排斥（Hubbard U）的低估——否则会把 Mott 绝缘体错判为金属。习惯称 LDA+U/LSDA+U，但可与 GGA、meta-GGA 等任意泛函组合。

**任务设计**：对 U = 0/3/8 eV（J = 0.95 eV）× 两种 PAW 组合（Ni+O 标准势 vs Ni_sv_GW+O_s）× 5 个体积缩放（0.95–1.05），共 30 个单点计算；用 ASE 批量生成以经验平衡晶格常数为中心的 POSCAR，用状态方程拟合得平衡晶格常数、带隙与在位磁矩。

**赝势选择论证**（教程核心方法学）：
- O 有 7 种可选势，原子组态相同，差别在软硬度（`_s` 特软低截断、`_h` 特硬适合短键分子）与 `_GW` 变体（适合涉及未占据态的计算如光学）。本例是块体、关注磁性（Ni 3d 基态性质），选便宜的 O_s 即可。
- Ni 有 4 种势，部分原子组态不同但都含 3d 壳层；额外把 p/s 纳入价电子对磁性影响不大；磁性是基态性质，无需 `_GW`。

**输入要点**：INCAR 在示例 3 基础上加 `ENCUT = 500`、`PREC = Acc`、`BANDGAP = WEIGHT`、`XC = CA`（LDA），U 变体加 `LDAU = .TRUE.`、`LDAUTYPE = 2`（Dudarev 形式，只依赖 U−J）、`LDAUL = 2 -1`、`LDAUU = 8.00 0.00`、`LDAUJ = 0.95 0.00`。

**结果与分析**：

![平衡晶格常数随 U 值变化（两种 PAW 势对比）](/imgs/2026-08-05/fig8.png)

- 平衡晶格常数对 U 的响应在两种 PAW 势间**截然不同**：Ni_sv_GW 把 Ni 的 p、s 轨道也放进价带，与 d 轨道杂化；U 增大时 d 态被推离费米能并收缩，纠缠的 sp 轨道随之膨胀，晶格常数增大。标准 Ni 势把 sp 冻结在芯里，无法膨胀，晶格常数反而因 d 收缩而略降。

![Ni 在位磁矩随 U 值变化](/imgs/2026-08-05/fig9.png)

- 增大 U 对在位磁矩影响显著，但磁矩源于 Ni d 轨道，PAW 势选择的影响很小。

![基带隙随 U 值变化](/imgs/2026-08-05/fig10.png)

- 带隙几乎不受 PAW 势选择影响。结论：选 PAW 势要综合考虑方法（此处 DFT+U 是核心）与目标量——最准的总是更贵的势（更多价电子或更高截断），但若目标量对势不敏感，就用便宜的。另注意：目标量依赖体积，许多文献直接用实验晶格常数，但平衡晶格常数本身依赖泛函选择（教程演示了带隙随晶格常数的变化）。

**复习问题**：什么是状态方程；选 PAW 势应依据什么量；`LDAUTYPE` 有哪些选项。

**本部分涉及的 VASP 功能与标签汇总**：`ISPIN = 2` 共线 SDFT、`vasp_std` vs `vasp_ncl`、`MAGMOM` 初设与简写语法（`2*0`、`4*5`）、`LMAXMIX`（磁性必设）、`LASPH`、`ISMEAR = 2` Methfessel-Paxton、交换劈裂与 Stoner 模型、`LORBIT = 11` 局域磁矩/PDOS、传播矢量 q 与磁单胞构造、`LKPOINTS_OPT` 一次计算同时出 SCF 与能带、`ICHARG = 1/11` 重启模式、`ALGO = None` 纯后处理、WAVECAR 与 k 网格变更的陷阱、`EFERMI = MIDGAP`、密度混合全套（`IMIX`/`AMIX`/`BMIX`/`AMIX_MAG`/`BMIX_MAG`）、晶体场理论与 eg/t2g 组合投影、DFT+U（`LDAU`/`LDAUTYPE`/`LDAUL`/`LDAUU`/`LDAUJ`）、状态方程拟合体积弛豫、PAW 势选择方法学（软硬度/_GW/价电子组态）、`BANDGAP = WEIGHT`。

---

## 3.2 Part 2：共线磁结构能量比较（Heisenberg 参数、SOC 与磁各向异性）

> 来源：https://vasp.at/tutorials/latest/magnetism/part2/ ｜ 练习文件：`magnetism-part2.zip` ｜ 示例目录：e05_NiO_Heisenberg、e06_bcc_Fe_SOC、e07_FeO_bulk_MAE

本部分从共线磁结构的总能比较出发提取 Heisenberg 模型参数，再进入自旋轨道耦合：先看 SOC 如何打破空间群对称性（bcc Fe 能带反交叉），最后完整走一遍 FeO 磁各向异性能（MAE）的计算流程。

### 示例 5：NiO 的 Heisenberg 模型参数（e05_NiO_Heisenberg）

**教学目标**：为 NiO 确定最近邻与次近邻 Heisenberg 耦合；按单胞大小缩放能量与 k 网格。

**模型与推导**：计入最近邻 $J_1$ 与次近邻 $J_2$ 的 Heisenberg 哈密顿量为

$$H = -J_1 \sum_{\langle i,j\rangle} \vec{e}_i\cdot\vec{e}_j - J_2 \sum_{\langle\langle i,j\rangle\rangle} \vec{e}_i\cdot\vec{e}_j$$

其中 $\vec{e}_i$ 为在位磁矩方向的单位矢量。考虑铁磁（FM）与两种反铁磁序——AF1（(001)/(110) 面内排序）与 AF2（(111) 面内排序，即 NiO 实际基态），每分子式能量为

$$E_{FM} = E_0 - 6J_1 - 3J_2, \quad E_{AF1} = E_0 + 2J_1 - 3J_2, \quad E_{AF2} = E_0 + 3J_2$$

反解得

$$J_1 = \frac{1}{8}(E_{AF1} - E_{FM}), \qquad J_2 = \frac{1}{24}(4E_{AF2} - 3E_{AF1} - E_{FM})$$

用共线磁性第一性原理计算得到三个总能即可定出模型参数，进而用模型哈密顿量求解临界温度等其他性质（详见 Archer 等, PRB 84, 115114, 2011）。

**输入要点**：
- AF2 用 4 原子菱方单胞（同 e03/e09），`MAGMOM = 5.0 -5.0 2*0`；AF1 与 FM 用 8 原子立方单胞（4 Ni + 4 O），AF1 按面交替设 MAGMOM，FM 用 `MAGMOM = 4*5.0 4*0`（`4*5` 是 `5 5 5 5` 的简写，不等价于单个 20）。
- 公共设置：`XC = CA` + LDA+U（U = 8、J = 0.95 eV）、`ENCUT = 450`、`PREC = Acc`、`EDIFF = 1E-7`、`ISMEAR = -5`（精确总能推荐四面体法）、`EFERMI = MIDGAP`、Pulay 混合参数、`LMAXMIX = 4`、`BANDGAP = WEIGHT`、`LORBIT = 11`、`LDAUPRINT = 1`。教程有意未充分收敛截断能与 k 网格以节省机时，实际应用必须做收敛测试。
- k 网格按晶格矢量长度反比选取（AF2 用 5×5×5、AF1/FM 用 4×4×4），保持可比网格密度；但精确总能更重要的是每个计算各自对网格密度收敛——极端情形下若一个磁结构是金属、另一个带隙，金属对长程相互作用更敏感、需要更细网格。

**计算与结果**：三个子目录依次运行后，用 py4vasp 读 `energy.read()["free energy TOTEN"]`，注意按单胞内分子式数归一（FM/AF1 除 4、AF2 除 2）再代入上式。教程结果：$J_1 \approx +1.4$ meV，$J_2 \approx -18.2$ meV——次近邻反铁磁耦合主导，与 NiO 的 AF2 基态一致。

**复习问题**：`MAGMOM = 2*3` 是否等于 `MAGMOM = 6`；精确总能最佳 ISMEAR；比较总能时可比 k 网格密度重要吗。

### 示例 6：SOC 打破空间群对称性——bcc Fe（e06_bcc_Fe_SOC）

**教学目标**：打开自旋轨道耦合；判断何时必须用 `vasp_ncl`。

**原理**：SOC 需在非共线磁性的自旋密度泛函框架下处理（Hobbs 等, PRB 62, 11556, 2000），体系用 2×2 自旋密度矩阵描述（含非对角 $n_{\uparrow\downarrow}$ 分量）。VASP 中 SOC 项采用零阶正则近似（ZORA）：

$$\hat{H}_{SO} = \frac{1}{(2c)^2}\frac{K^2(r)}{r}\frac{dV(r)}{dr}\,\vec{\sigma}\cdot\vec{L}, \qquad K(r) = \frac{1}{1 - V(r)/2c^2}$$

SOC 能量尺度很小，后果却可以显著：不含 SOC 时 bcc Fe 的 Γ–H₃ 与 Γ–H 两条高对称路径等价，含 SOC 后出现特征性的能带**反交叉**（anticrossing）。

![bcc Fe 含 SOC 的能带劈裂示意](/imgs/2026-08-05/fig11.png)

**输入要点**：bcc Fe 原胞（a ≈ 2.86 Å）。INCAR：`ENCUT = 300`、`EDIFF = 1e-8`、`ALGO = Normal`、`ISMEAR = 0`、`SIGMA = 0.05`、`LNONCOLLINEAR = T`、`LSORBIT = T`、`MAGMOM = 0.0 0.0 2.0`（注意格式变为每原子三分量矢量，z 方向磁化）、`LASPH = T`、`LMAXMIX = 4`、混合参数（`BMIX` 取极小但非零，0 会使某些版本崩溃）、`LORBIT = 11`。KPOINTS 为 12×12×12 均匀网格，另备 KPOINTS_OPT 线模式路径 H₃–Γ–H（141 个分段点）。

**流程**：均匀网格与路径在同一文件中给出，一次运行（无需重启）即同时完成 SCF 与能带，用 `vasp_ncl` 运行（ncl = noncollinear，处理完整 2×2 自旋密度矩阵，SOC 计算必需）。注意费米能要在均匀网格上计算，py4vasp 绘图时以 SCF 结果为准。

**复习问题**：何时必须用 vasp_ncl；SOC 为何要求非共线框架；反交叉的物理来源。

### 示例 7：FeO 的磁各向异性能（e07_FeO_bulk_MAE）

**教学目标**：完整掌握 MAE 计算工作流——构造原初磁单胞、确定磁化方向扫描路径、用 `SAXIS` 旋转自旋量子化轴、非自洽 SOC 计算并拟合各向异性能。

**方法概述**：MAE 是磁化方向依赖的微小能量差。教程采用"固定电荷密度、旋转 SAXIS"方案：先做共线 LDA+U SCF，再以 `ICHARG = 11` 固定密度做一系列含 SOC 的非共线非自洽计算，通过 `SAXIS` 改变自旋量子化轴间接旋转磁化方向。优点是不必担心计算收敛到别的磁化方向；前提是 SOC 是小微扰（密度与磁化几乎不因 SOC 改变），并假定在位轨道角动量由晶体场钉扎而非跟随磁矩方向。

**路径定义**：以极角 θ、φ 表征 MAE 路径（取 Schrön 等, PRB 86, 115134 的约定）：θ = 0 对应磁化平行立方 [111] 方向，φ = 0 对应磁化位于 (1−10) 面。难点在于 FeO 有轻微单斜畸变（不是严格立方），教程一步步演示：

![FeO 2x2x2 超胞](/imgs/2026-08-05/fig12.png)

1. **立方结构起步**：从 Materials Project 取 FeO 岩盐结构 CIF，用 ASE 构造能容纳 q = (1/2,1/2,1/2) 自旋波的 2×2×2 超胞。
2. **施加自旋波与畸变**：按 Schrön 等 Table II 的参数（FeO：r = 0.0014、e = 0.0253、t = 0.0251）构造应变张量施加单斜畸变，Spglib 确认空间群降为 C2/m (12)；自旋波通过把自旋向下位点的 Fe 临时替换为 Mn 作标记来可视化。

![施加单斜畸变后的超胞](/imgs/2026-08-05/fig13.png)

![带自旋波标记的超胞](/imgs/2026-08-05/fig14.png)

3. **原初磁单胞**：超胞便于可视化，但计算用原胞更便宜。用与 NiO 相同的变换 $A_{prim} = P\cdot A_{conv}$（坐标按 $(P^{-1})^T$ 变换并过滤到 [0,1)），得到 4 原子原初磁单胞（2 Fe + 2 O，MAGMOM 共线形式 `0 0 4 -4`、非共线形式 `0 0 0 0 0 0 0 0 4 0 0 -4`），Spglib 复核空间群仍为 C2/m。

![原初磁单胞](/imgs/2026-08-05/fig15.png)

4. **在畸变结构中辨认立方方向与面**：[111] 方向即原点指向第 7 个原子的方向（畸变不改变原子顺序）；(1−10) 面按 Miller 指数定义（截 x 轴于 1、截 y 轴于 −1、平行 z 轴）。
5. **生成 SAXIS 列表**：φ = 0 时把 e111（[111] 单位矢量）绕 e_110（[−110] 单位矢量）以 5° 步进旋转，θ 从 0 到 180° 共 37 个方向；φ = 90° 时先把旋转轴绕 e111 转 90° 再做同样的 θ 扫描。用 Scipy Rotation 生成列表，打印成 bash 数组。

**VASP 输入要点**：
- 共线 SCF（INCAR.ispin2）：`ENCUT = 520`、`NBANDS = 32`、`ISPIN = 2`、`LNONCOLLINEAR = F`、`LSORBIT = F`、`MAGMOM = 0 0 4 -4`、`GGA_COMPAT = F`（SOC 计算需关闭对称性兼容的梯度校正）、`EDIFF = 1E-8`、`ISMEAR = 0`、`SIGMA = 0.01`、`EFERMI = MIDGAP`、`XC = CA` + `LDAU`（Fe 的 d 加 U = 4 eV，`LDAUL = -1 2`）、混合参数与 `LMAXMIX = 4`。
- 路径计算（INCAR.path）：`NBANDS = 64`、`ISPIN = 1`、`LNONCOLLINEAR = T`、`LSORBIT = T`、MAGMOM 为每原子三分量形式、`ISMEAR = -5`、`LORBMOM = T`（输出轨道矩）、`LCHARG = F`、`ICHARG = 11`，末尾逐点追加 `SAXIS = x y z`。
- 重启规则：从共线计算重启非共线计算靠读入电荷密度+磁化，因此 ENCUT 必须一致（FFT 网格不变）；共线的 NBANDS 按每个自旋通道计，要张成同样的空间需设为非共线运行的两倍；`ICHARG = 11` 固定密度时建议保持同一 k 网格（4×4×4）。

**计算与分析**：run 脚本先跑共线 SCF 并备份 CHGCAR.ispin2，然后对 φ = 0 与 φ = 90 两条路径的每个 SAXIS 点复制 CHGCAR、跑 `vasp_ncl`、收集 OSZICAR 的 `F=`。以 θ 为横轴、相对 θ = 0 的能量为 MAE 作图，用 $K\sin^2\theta$ 拟合 φ = 90 路径提取各向异性常数 K；对比参考数据与自己计算的结果，并讨论两条路径是否等价——答案是否定的，因为单斜畸变降低了立方对称性。

**复习问题**：轨道矩定义在哪个坐标系；哪个标签控制旋量空间相对笛卡尔空间的取向（`SAXIS`）；MAE 计算工作流包含哪些步骤。

**本部分涉及的 VASP 功能与标签汇总**：Heisenberg 模型参数提取（FM/AF1/AF2 能量映射、按分子式归一）、磁结构比较的 k 网格策略、`LNONCOLLINEAR`/`LSORBIT`/`vasp_ncl`、MAGMOM 三分量格式、ZORA 形式的 SOC、能带反交叉、KPOINTS+KPOINTS_OPT 单次计算、`GGA_COMPAT = F`、`ICHARG = 11` 固定密度 SOC 路径、`SAXIS` 自旋量子化轴旋转、共线→非共线重启规则（ENCUT/NBANDS 关系）、`LORBMOM` 轨道矩、`LDAUPRINT`、ASE/Spglib/Scipy 构造磁单胞与扫描路径、MAE 的 $\sin^2\theta$ 拟合。

---

## 4. 分子动力学（Molecular Dynamics）

分子动力学（MD）按每个时间步作用在各粒子上的力模拟原子（与分子）的运动。学习如何在 VASP 中进行 MD 模拟。

## 4.1 Part 1：硅的熔化（从头算分子动力学基础）

> 来源：https://vasp.at/tutorials/latest/md/part1/ ｜ 练习文件：`md-part1.zip` ｜ 示例目录：e01_solid-cd-Si、e02_melting-Si、e03_monitoring

本部分是 AIMD 的完整入门：跑通一次 NVT 系综的从头算 MD、用对关联函数观察硅的熔化、学习恒温器家族与 ICONST 内坐标监控，并体会 AIMD 的机时代价。

### 示例 1：固态金刚石硅的 NVT 从头算 MD（e01_solid-cd-Si）

**教学目标**：说明什么是从头算 MD；用 pymatgen 从 CIF 创建超胞；区分 `vasp_gam` 与 `vasp_std`；用 py4vasp 绘制 MD 中各类能量随时间步的演化与轨迹动画。

**概念**：MD 用经典运动方程模拟给定温度下的原子运动，每个离子步是一个时间步；力由量子力学（第一性原理）计算时称为从头算 MD。正则系综（NVT）要求粒子数、体积、温度恒定；温度效应由恒温器引入——通过在运动方程中加入随机项或确定性附加动力学变量来模拟热浴。Nosé-Hoover 属于确定性方案，在小/刚性体系中可能缺乏遍历性（经典反例是单丁烷分子模拟），但对大体系完全适用。

**结构准备**：用 pymatgen 读取 Materials Project 的 cd-Si CIF（CIF 是国际晶体学联合会规定的晶体数据交换标准格式），`my_struc.make_supercell(2)` 构造 2×2×2 常规单胞超胞（64 原子）并写为 POSCAR。为什么用超胞：超胞尺寸限制了可容纳的晶格振动最大波长，MD 需要捕获影响动力学的晶格振动；实际工作应对原胞大小做收敛检查。

**输入要点**：
- 电子部分：`ISMEAR = 0`、`SIGMA = 0.1`、`LREAL = Auto`（投影算符放实空间，大体系提速）、`ALGO = VeryFast`（RMM-DIIS 电子弛豫）、`PREC = Low`（MD 可用低精度提速）、`ISYM = 0`（不施加对称性）。
- MD 部分：`IBRION = 0`（开启 MD）、`NSW = 30`（离子步数）、`POTIM = 3.0`（时间步 3 fs——MD 中 POTIM 是时间单位）、`MDALGO = 2`（Nosé-Hoover）、`SMASS = 1.0`（Nosé 质量/热惯性，调节体系与热库间的能量流）、`TEBEG = TEEND = 2000`（恒温 2000 K 即 NVT）、`ISIF = 2`（固定晶胞形状与体积）。
- KPOINTS：大超胞只用 Γ 点，用 `vasp_gam` 运行（该可执行文件只能做 Γ 点计算，内部把 KS 轨道按实数处理，比 `vasp_std` 快）。

**输出解读**：stdout 每个离子步末行给出 `T`（瞬时温度）、`E`（总能）、`F`（自由能）、`E0`、`EK`（离子动能）、`SP`（Nosé 恒温器势能）、`SK`（恒温器动能）等；py4vasp `mycalc.energy[:].plot()` 可分别绘出 TOTEN（离子-电子能）、EKIN、TEIN（瞬时温度）、ES（Nosé 势）、EPS（Nosé 动）、ETOTAL（守恒量）随时间步的曲线；`mycalc.structure[:].plot()` 动画显示 2000 K 下 90 fs 内金刚石硅开始熔化的过程。

**复习问题**：AIMD 的轨迹是量子粒子的路径吗（否，离子经典处理）；晶体结构信息的标准存储格式（CIF）；`vasp_gam` 适用什么计算；NVT 系综中总能是否守恒（ETOTAL 守恒，离子势能+动能+恒温器项之和）；MD 为何用大超胞。

### 示例 2：硅的熔化——对关联函数分析（e02_melting-Si）

**教学目标**：计算并解释对关联函数；重启 MD 模拟；指定模拟时长；在 XDATCAR 中找轨迹；用 `MAXMIX` 降低成本。

**概念**：液体结构用径向分布函数（对关联函数）研究。严格地说，它是距参考粒子中心距离 r 处的相对局域粒子密度（相对平均密度）：

$$g(r) = \frac{\Omega_r}{N^2}\left\langle \sum_i \sum_{j\neq i} \delta(\mathbf{r} - \mathbf{R}_{ij}) \right\rangle$$

其中系综平均按遍历性假设用时间平均实现。晶体有长程序，g(r) 出现尖锐的远距离峰；液体粒子结合牢固但不刚性，熔化时 g(r) 特征显著变化。对各向同性液体，任何对函数 $\mathcal{B}$ 的系综平均都可经 g(r) 计算：$\mathcal{B}_{obs} = \frac{1}{2}N\rho\int_0^\infty b(r)g(r)4\pi r^2 dr$（Allen & Tildesley《Computer Simulation of Liquids》2.6 节）。

**流程**：把示例 1 的 CONTCAR 拷贝为本例 POSCAR 续算 90 fs（总 180 fs）；对关联函数写入 PCDAT，用教程提供的 awk+gnuplot 脚本解析绘图。

![90 fs 与 180 fs 的对关联函数对比](/imgs/2026-08-05/fig16.png)

**要点**：90 fs 时体系仍接近晶体，长程序使 g(r) 在 4 Å 以外仍有明显峰；180 fs 后趋于液体特征（短程有序、长程平坦）。注意 g(r) 只在 r 小于超胞最短边长一半时有意义（此处 < 5.4 Å）。

![XDATCAR 轨迹与中心盒表示示意](/imgs/2026-08-05/fig17.png)

**轨迹表示**：XDATCAR 记录粒子穿越周期边界的真实轨迹；"中心盒表示"则对每步位置施加周期性折叠（把跑出盒子的原子替换为其周期镜像）。

**重启与调参**：`cp CONTCAR POSCAR` 续算；`NPACO`/`APACO` 控制 g(r) 的最大距离与网格点数，`NBLOCK`/`KBLOCK` 控制写出频率。练习把 `APACO = 5.4`、`NSW = 5`（再跑 15 fs）重启：新 PCDAT 只对应这最后 15 fs 的采样（总模拟时间 195 fs），统计质量明显下降。

![90/180/195 fs 对关联函数对比](/imgs/2026-08-05/fig18.png)

**复习问题**：长程序如何反映在 g(r) 中；CONTCAR 与 POSCAR 格式是否相同、CONTCAR 含什么信息；`APACO`/`NPACO` 的作用；轨迹在哪个输出文件（XDATCAR）。

### 示例 3：监控分子几何——恒温器对比与 ICONST（e03_monitoring）

**教学目标**：用 ICONST 文件指定几何坐标；通过 REPORT 文件监控坐标；估算模拟耗时；用 Langevin 恒温器模拟 NpT 系综；独立构造超胞 POSCAR。

**恒温器家族**（`MDALGO` 选择）：`MDALGO = 1` Andersen（完全随机：用随机碰撞模拟热浴）、`MDALGO = 2` Nosé-Hoover（确定性）、`MDALGO = 3` Langevin（随机+确定性）。对分子内强作用力（如键伸缩）的体系，Andersen 与 Langevin 无遍历性问题且平衡化高效；Nosé-Hoover 在小/刚性体系缺遍历性。实际选择主要取决于 VASP 中恒温器与热力学系综的可用组合（见 Wiki 的 Ensembles 分类）。

**ICONST**：用 ICONST 文件指定内部几何坐标（键长、角度、二面角、到原点距离等）。第一列定义坐标类型（如 LR 晶格矢量长度、LA 夹角、LV 体积、R 原子间距），随后整数按 POSCAR 顺序引用原子，最后一个整数指定动作（7 = 监控，0 = 约束）。监控量逐离子步写入 REPORT 文件。

**实践**：构造 cd-Si 原胞的 2×2×2 超胞（16 原子），本例 INCAR 还演示了 `ISIF = 3` + `MDALGO = 3`（Langevin + 晶格自由度 → NpT 系综）、`LANGEVIN_GAMMA`（原子摩擦系数）、`IVDW = 10`（DFT-D2 色散校正）等；KPOINTS 用 2×2×2 均匀网格（因此必须用 `vasp_std` 而非 `vasp_gam`）。运行后从 REPORT 提取 `mc> LV` 行绘制体积随 MD 步数的弛豫曲线。

**机时估算**：OUTCAR 中 `Total CPU time used` 显示 10 个 MD 步约 74 s，外推 10000 步需要 20 小时以上——10000 步是从 MD 得出可靠结论的典型量级，这直观说明 AIMD 的昂贵，也是引入机器学习力场的动机。

**复习问题**：ICONST 中 `R 1 6 0` 定义什么（原子 1 与 6 的距离约束）；如何估算长 MD 的机时；`LANGEVIN_GAMMA` 指定什么。

**本部分涉及的 VASP 功能与标签汇总**：`IBRION = 0` MD、`MDALGO`（1=Andersen、2=Nosé-Hoover、3=Langevin）、`SMASS` Nosé 质量、`TEBEG`/`TEEND`、`POTIM`（fs 时间步）、`NSW`、NVT/NpT 系综实现（`ISIF = 2/3` + 恒温器组合）、`LREAL = Auto`、`ALGO = VeryFast`、`PREC = Low`、`ISYM = 0`、`vasp_gam` vs `vasp_std`、stdout 能量分量（T/E/F/EK/SP/SK）、XDATCAR/PCDAT/REPORT 输出、`NPACO`/`APACO`/`NBLOCK`/`KBLOCK`、CONTCAR→POSCAR 重启 MD、ICONST 内坐标（LR/LA/LV/R/S + 动作码）、pymatgen CIF→超胞→POSCAR、py4vasp 能量与轨迹分析、AIMD 机时估算方法学。

---

## 4.2 Part 2：机器学习力场（在线训练、验证、系综平均与可迁移性）

> 来源：https://vasp.at/tutorials/latest/md/part2/ ｜ 练习文件：`md-part2.zip` ｜ 示例目录：e04_MLFF、e05_MLFF-test、e06_MLFF-production、e07_transferability

本部分展示 MLFF 的完整生命周期：在 AIMD 中在线训练 → 用物理性质测试 → 纯力场生产运行 → 检验可迁移性并计算热膨胀系数。

### 示例 4：在从头算 MD 中在线训练力场（e04_MLFF）

**教学目标**：在 ab-initio MD 过程中在线（on-the-fly）训练机器学习力场。

**动机**：良好的热力学平均需要长模拟时间，而每个 AIMD 步昂贵；经典力场便宜但依赖经验参数；MLFF 从第一性原理数据学习物理，兼顾精度与长时间模拟。

**输入要点**（16 原子 Si，NpT + Langevin）：`ML_LMLFF = T` 开启 MLFF、`ML_ISTART = 0` 从 AIMD 训练新力场；MD：`IBRION = 0`、`NSW = 10000`、`POTIM = 2.0`、`MDALGO = 3`、`LANGEVIN_GAMMA = 1`（原子摩擦）、`LANGEVIN_GAMMA_L = 10`（晶格摩擦）、`PMASS = 10`（晶格质量，NpT 所需）、`TEBEG = 400`、`ISIF = 3`；`ICONST` 用 LR/LA/LV 监控晶格矢量长度、夹角与体积；`RANDOM_SEED` 保证可重复性（学到的局域参考构型数强烈依赖随机种子，训练示例尽量短）；训练用 2×2×2 k 网格——力必须对 k 点收敛，虽然力场可以在小单胞训练后用于大体系。

**在线训练算法（三种触发）**：① 无力场时，把当前结构的局域参考构型加入训练集并构建力场；② 任意 MD 步若任一原子的贝叶斯力误差超过严格阈值 `ML_CDOUB × ML_CTIFOR`，加入训练集并重建；③ 构建后进入 `ML_NMDINT` 步热身期（误差判据放宽），之后若误差超过低阈值 `ML_CTIFOR`，结构进入待处理列表，列表达到 `ML_MCONF_NEW` 长度后批量加入训练集。注意贝叶斯误差只是样本外误差的估计，并非力场质量的真实指标。

**输出文件**：训练数据写入 ML_ABN，最终力场参数存于 ML_FFN。

### 示例 5：用离子弛豫测试力场（e05_MLFF-test）

**教学目标**：通过对比弛豫晶格参数与参考数据测试力场；在力场下使用共轭梯度算法。

**测试策略**：① 独立测试集（在随机构型上对比 DFT 与 MLFF 的力、应力、能量差）；② 对比物理性质（弛豫晶格参数、声子、相相对能量、弹性常数、缺陷形成能等）。本例采用后者。

**输入要点**：`ML_LMLFF = T` + `ML_ISTART = 2`（只读已训练力场、不做 DFT）；`IBRION = 2`、`ISIF = 3`；把 ML_FFN 软链接为 ML_FF 读入。bash 循环 13 个体积缩放因子（0.96–1.03）做静态单点，提取 OUTCAR 的 volume 与 free energy 绘制 E-V 曲线，与相同 DFT 设置的参考数据对比——MLFF 曲线与 DFT 高度吻合。

**概念**：训练只能得到样本内误差；贝叶斯误差是样本外误差的估计而非真实度量。

### 示例 6：系综平均的晶格常数与体积（e06_MLFF-production）

**教学目标**：把晶格常数/体积作为系综平均计算；用预训练力场运行纯 MLFF 的 MD 生产运行。

**概念**：晶格参数常被视为 T = 0 K 基态性质，但许多材料随温度发生结构转变，且实验测量的就是热力学系综平均——用 MLFF 在 400 K（训练温度）做生产运行。

**操作**：INCAR 不再含任何 ab-initio 标签（无 ENCUT 等），仅 MD + `ML_LMLFF = T` + `ML_ISTART = 2`；从 400 K 热化构型出发跑 3×10000 步。耗时对比：同样模拟 AIMD 需 20 小时以上，MLFF 不到 10 分钟。REPORT 文件记录每步 ICONST 监控量，用脚本提取 LV/LR 并对时间取平均：体积约 319.3 Å³、平均晶格矢量约 7.673 Å；三个晶格矢量长度的微小不对称源于模拟时间不够长。

### 示例 7：力场可迁移性与热膨胀系数（e07_transferability）

**教学目标**：检验力场是否适用于新的参数集/体系；独立从 REPORT 提取监控坐标；计算热膨胀系数。

**可迁移性判据**：训练后若 MD 中贝叶斯误差突增，说明进入了未知相空间（如相变/结构转变），需要继续训练——把 ML_ABN 拷贝为 ML_AB 并设 `ML_ISTART = 1` 即可续训（引入缺陷、表面等新相空间时同理）。但改变温度等参数未必需要重训：只要已有局域参考构型充分覆盖目标系综的相空间即可。检验方法（无相变前提下）：在训练参数集与目标参数集分别用 MLFF 跑足够长的 MD，对比监控量（如晶胞体积）的概率分布——200 K 遇到的所有局域构型都落在 400 K 相空间之内，则力场可迁移。

![200 K 与 400 K 体积概率密度对比（可迁移性验证）](/imgs/2026-08-05/fig19.png)

**计算与结果**：在 T200/T300 子目录软链接 e04 的 ML_FFN 为 ML_FF，各跑 3×10000 步 Langevin NpT MD；用 get_lattice_parameters.sh 从 REPORT 提取并平均：200 K 体积 319.065 Å³、平均晶格矢量 7.6709 Å；300 K 体积 319.192 Å³、7.6723 Å。联合示例 6 的 400 K 结果（319.287 Å³）画 V(T) 曲线，在此区间近似线性，体膨胀系数

$$\alpha_V = \frac{1}{V}\frac{dV}{dT} \approx 3.47\times10^{-6}\ \mathrm{K}^{-1}$$

实验通常测线膨胀系数 $\alpha_L = \frac{1}{L}\frac{dL}{dT}$（Si 实验值 2.6×10⁻⁶ K⁻¹，Okada & Tokumaru, JAP 56, 314, 1984）；立方晶格有 $\alpha_L = \alpha_V/3$，得 1.16×10⁻⁶ K⁻¹。也可用三个等价晶格矢量长度合并统计直接算 $\alpha_L$（得 1.38×10⁻⁶ K⁻¹），两种途径的差异反映有限采样统计误差。

**复习问题**：力场可迁移的条件；热膨胀系数的定义；监控几何坐标写入哪个文件（REPORT）。

**本部分涉及的 VASP 功能与标签汇总**：`ML_LMLFF`、`ML_ISTART`（0 从头训/1 续训/2 只读运行）、在线训练三触发机制（`ML_CTIFOR`/`ML_CDOUB`/`ML_NMDINT`/`ML_MCONF_NEW`）、贝叶斯力误差的解读边界、ML_AB/ML_ABN（训练数据）与 ML_FF/ML_FFN（力场参数）、`RANDOM_SEED` 可重复性、`PMASS`/`LANGEVIN_GAMMA_L` NpT 晶格动力学、E-V 曲线测试法、系综平均晶格常数、可迁移性的相空间覆盖判据、热膨胀系数（$\alpha_V$ 与 $\alpha_L$）。

---

## 4.3 Part 3：氯甲烷的氯离子取代反应（SN2 自由能面与反应速率）

> 来源：https://vasp.at/tutorials/latest/md/part3/ ｜ 练习文件：`md-part3.zip` ｜ 示例目录：e08_slow-growth-SN2、e09_free-MD-SN2、e10_reaction-rate-SN2

本部分以对称 SN2 反应（Cl⁻ + CH₃Cl → ClCH₃ + Cl⁻）为载体，演示约束 MD 与热力学积分计算自由能面、反应物态概率分布、过渡态广义速度，最终组装出反应速率常数与表观活化自由能。全部 MD 用 MLFF 完成。

### 示例 8：slow-growth 模拟自由能剖面（e08_slow-growth-SN2）

**教学目标**：用热力学积分获得自由能剖面；执行并评判 slow-growth 模拟质量；解读 ICONST 与 REPORT。

**背景与反应**：自由能剖面给出活化自由能——从稳定反应物态到过渡态所需的功。热力学积分沿反应路径（初末热力学态之间的路径，反应坐标 ξ 是其一维表示）积分自由能梯度。slow-growth AIMD 能捕获反应中可能发生的转变（如转动等软自由度）：变换足够慢则过程近似可逆，对体系做的功等于自由能。对照方案：fast-switching 用纯静态多次快速变换+Jarzynski 恒等式；静态 PES 拓扑分析（只需少数驻点的弛豫与振动分析）便宜但无法计入过渡中激活的软模，只给定性估计。

![SN2 反应示意](/imgs/2026-08-05/fig59.png)

本反应是对称 SN2：Cl⁻ 攻击 CH₃Cl 形成可检测的中间态，再释放原来的 Cl⁻；因进攻与离去基团相同，自由能剖面对称。

**输入要点**：
- POSCAR：12 Å 盒子中 C + 3H + 2Cl；POTCAR 顺序对应 C、H、Cl。
- `POMASS = 12.011 3.000 35.453`：把 H 换成更重的氚同位素质量，减缓 H 键振荡从而允许更大时间步。**同位素效应讨论**：实验中的动力学同位素效应是离子自由度的量子效应，经典 MD（含 MLFF-AIMD）中质量在热力学平均的方程中消去，观察不到同位素效应。
- `KSPACING = 1000`：设得极大以致只初始化 Γ 点（等价于单 k 点 KPOINTS）。
- `MDALGO = 1` + `ANDERSEN_PROB = 0.05` + `TEBEG = 300` + `ISIF = 2` → NVT 系综（Andersen 用随机碰撞模拟热浴）。
- ICONST 定义反应坐标：原子 1 为 C、5/6 为两个 Cl；`R 1 5 0`、`R 1 6 0` 监控两个 C-Cl 键长，`S 1. -1. 5` 以组合坐标（两键长之差）作约束。例如分子内 C-Cl 1.8 Å、远端 Cl⁻ 距 C 3 Å 时 ξ = −1.2 Å；过渡态两个 Cl 等距，ξ = 0。找到合适的反应坐标是面对新反应的主要挑战。
- `INCREM`：约束 MD 沿反应路径演化的变换速度（单位是每步=每 POTIM，而非每时间）。`LBLUEOUT = T`：计算自由能梯度并写入 REPORT。

**计算流程**：先用 INCAR.refit 对 ML_AB 中的 ab-initio 能量/力/应力重新拟合出力场（ML_FFN→ML_FF），再跑 run.sh：两段 MD（各 1500 步，总 3000 步），段间 `cp CONTCAR POSCAR` 续结构，并从 REPORT 读取上一段末尾的 RANDOM_SEED 写入新 INCAR 以保证速度初始化连续可重复。

**热力学积分**：自由能梯度 ∂A/∂ξ 在各 ξ′ 处采样后，

$$\Delta A_{A\to B} = \int_{\xi_A}^{\xi_B} \left.\frac{\partial A}{\partial \xi}\right|_{\xi'} d\xi'$$

slow-growth 中 $\mathrm{d}\xi' = \dot{\xi}'\,\mathrm{d}\tau$（$\dot{\xi}'$ 由 INCREM 设定）。若变换速度与时间步不够小则方法失效——质量检查：把过程反向跑，自由能剖面出现滞后（hysteresis）即说明参数过大；对称反应中剖面偏离对称也提示同样问题。slow-growth 只适合快速预览自由能面；高质量热力学积分用 Blue Moon 方法（收敛可控但更贵）。

**REPORT 解读**：每个 MD 步含 SHAKE 收敛信息（SHAKE 算法约束快自由度）、`cc> S ξ 值...`（约束坐标当前值）、`b_m>` 行（Blue Moon 量：lambda、|z|^(−1/2)、GkT 及其组合——第一个数即自由能梯度）、能量与温度块。用脚本把 `cc>` 第 3 列与 `b_m>` 第 2 列拼成数据文件，`scipy.integrate.cumulative_trapezoid` 数值积分得 ΔA(ξ)。ξ < 0 为反应物态、ξ > 0 为产物态、ξ = 0 为过渡态；代表构型约在步 0（反应物）、750（过渡态）、1499（产物）。

**复习问题**：哪个文件指定约束 MD 的约束（ICONST）；INCREM 沿哪条路径作用；自由能剖面总是对称的吗（否，仅对称反应）；slow-growth 的典型用途与失效条件。

### 示例 9：反应物态的概率分布（e09_free-MD-SN2）

**教学目标**：计算反应坐标的概率密度；监控几何参数。

**原理**：若反应在给定温度下不被激活，转变不会自发发生——可通过反应坐标概率密度确认：

$$P(\xi') = \langle \delta(\xi' - \xi(\mathbf{R})) \rangle_R, \quad \xi' < 0$$

从反应物系综出发做自由（无约束）MD 采样。

**输入要点**：软链接示例 8 的 ML_FF；INCAR 与示例 8 相比**去掉** `LBLUEOUT` 与 `INCREM`——不写自由能梯度、不强制沿反应路径，纯自由 MLFF-MD（`ML_MODE = run`，只读不学）；`NSW = 5000`、`POTIM = 2`、Andersen NVT 300 K。ICONST 中两个 C-Cl 键长动作码为 0（约束但不驱动），组合坐标 `S 1. -1. 7` 的动作码 7 表示监控。从 REPORT/轨迹统计 ξ 的直方图得 P(ξ)，可见反应物态在 ξ 负值处形成峰（可用谐函数拟合得到力常数与参考态概率 P(ξ_ref,R)）。

### 示例 10：过渡态平均速度与速率常数（e10_reaction-rate-SN2）

**教学目标**：对 MD 步数做收敛研究；计算反应速率与表观活化能；运行约束 MD。

**原理**：速率常数可写为

$$k_{R\to P} = \frac{\langle |\dot{\xi}^*| \rangle}{2} P(\xi_{ref,R}) \exp\left(-\frac{\Delta A_{\xi_{ref,R}\to\xi^*}}{k_B T}\right)$$

表观活化自由能与可逆功的关系为

$$\Delta A^\ddagger = \Delta A_{\xi_{ref,R}\to\xi^*} - k_B T \ln\left(\frac{h}{k_B T}\frac{\langle|\dot{\xi}^*|\rangle}{2} P(\xi_{ref,R})\right)$$

示例 8、9 已给出除 $\langle|\dot{\xi}^*|\rangle$ 外的全部量。广义速度定义为 $\dot{\xi} = \sum_i \sum_\mu \frac{\partial \xi}{\partial R_{i,\mu}} \dot{R}_{i,\mu}$。

**输入要点**：POSCAR 直接取过渡态构型；INCAR 同示例 9 但不设 INCREM，ICONST 组合坐标动作码为 0（约束在 ξ = 0）——即过渡态上的约束 MD；`LBLUEOUT = T` 写出相关量；`NSW = 600` 起步，逐步加长做步数收敛研究。从 REPORT 提取广义速度取绝对值平均，代入上式得速率常数与表观活化自由能。

**复习问题**：约束 MD 与自由 MD 的输入差异；速率公式中各量的来源；收敛研究如何做。

**本部分涉及的 VASP 功能与标签汇总**：约束 MD 全套（ICONST 的 R/S 坐标与动作码 0/5/7、`INCREM` 变换速度、`LBLUEOUT` 自由能梯度输出、REPORT 的 cc/b_m 行）、slow-growth 方法（可逆性、滞后检验、对称性检验）与 BlueMoon 方法定位、热力学积分数值实现、`POMASS` 同位素质量与经典 MD 无同位素效应的论证、`KSPACING` 极简 k 采样、`MDALGO = 1` + `ANDERSEN_PROB` NVT、`ML_MODE = run`、RANDOM_SEED 段间传递、ML_AB→ML_FFN 重拟合（refit）、SN2 反应速率理论组装（P(ξ_ref)、ΔA、⟨|ξ̇*|⟩ → k 与 ΔA‡）。

---

## 5. 机器学习力场（Machine-Learned Force Fields）

机器学习力场与从头算分子动力学结合，既能从第一性原理捕捉物理本质，又能以较低成本达到较长模拟时间。教程讲解 MLFF 的误差分析与超参数优化。

## 5.1 Part 1：机器学习力场的误差分析与超参数优化（LiH 案例）

> 来源：https://vasp.at/tutorials/latest/mlff/part1/ ｜ 练习文件：`mlff-part1.zip` ｜ 示例目录：`$TUTORIALS/MLFF/e01_error_analysis`、`e02_hyperparameter_scan`、`e03_sparsification`、`e04_production`
> 注：教程声明——本教程使用的力场与测试集是专门裁剪过的教学数据，仅用于在合理计算代价内演示完整工作流，不能直接用于科研应用。

本部分围绕 VASP 的机器学习力场（Machine-Learned Force Field, MLFF）展开，主线是"训练好一个力场之后，如何科学地评估它、改进它、并最终投入使用"。全部工作基于同一体系 LiH（锂氢化物），共四个层层递进的示例：误差分析 → 精度调优 → 性能调优 → 生产运行。

### 示例 1：训练集与测试集误差分析

**教学目标**：对给定参考数据计算 MLFF 模型的训练集误差与测试集误差，并学会解读两者的关系。

**原理与概念**：

- 误差估计回答两个问题：训练好的模型有多准？它对训练集之外的结构泛化能力如何？
- **训练集误差**（training-set error）：在训练集上评估，即拟合所用 DFT 数据（能量、力、应力）与 MLFF 预测值之间的误差。
- **测试集误差**（test-set error）：在训练过程中从未使用过的外部测试集上评估。良好实践要求测试集构型在与生产运行相同的条件下采样——例如测试结构应与生产运行有相同的原子数（即使训练是在更小体系上做的），且必须采自与生产运行相同的热力学相。这样测出的误差才能反映模型从训练集到生产条件的真实泛化能力。
- 三种典型诊断情形：
  1. **训练误差低 + 测试误差高 → 过拟合**：力场外推能力差，需要补充更多训练结构或调整超参数；
  2. **两者大致同量级且都足够低 → 力场可用**；
  3. **训练误差高 + 测试误差低 → 测试集有偏**（不够一般化），不能说明力场好。
- 延伸阅读：Tokita & Behler 关于 MLFF 验证的论文（arXiv:2308.08859，DOI: 10.1063/5.0160326）第 IV 节。

**输入要点**：本示例不做任何 DFT 计算，而是直接利用在线学习（on-the-fly learning）已经收集好的第一性原理数据 `ML_AB` 文件（含 LiH 的晶格矢量、原子位置、能量、力、应力张量）做**重新拟合（refit）**：

- `INCAR.refit`：仅两行——`ML_LMLFF = TRUE` 打开机器学习模式，`ML_MODE = refit` 选择重拟合模式。
- `INCAR.run`：`IBRION = -1`（不做离子更新，即单点）+ `ML_LMLFF = TRUE` + `ML_MODE = run`（纯力场预测模式）。
- `POSCAR`：LiH 原胞（立方，晶格常数约 4.03 Å，4 Li + 4 H）。
- `KPOINTS`：Gamma 单点（1×1×1）。
- `POTCAR`：只需原子种类数正确的"哑文件"，会被读取但不参与计算。

三个关键理解：

1. **为什么 refit 总应在 on-the-fly 训练获得 ML_AB 之后进行？** 因为 refit 用奇异值分解（SVD）求解线性回归（贝叶斯线性回归框架），比在线训练使用的 evidence approximation 更精确；且 refit 之后 `ML_MODE = run` 的生产运行更快（高性能算法依赖 refit 产生的模型格式）。
2. **为什么 POSCAR 不影响 refit 结果？** 拟合完全基于 ML_AB 中的数据；POSCAR 只用于做一次预测步（其能量/力会出现在 OUTCAR 中，但本练习不关心）。
3. **POTCAR 中的原子类型不重要**：机器学习算法从 ML_AB 独立地确定原子类型。

测试集位于 `$TUTORIALS/MLFF/test_set`，含 50 个独立结构的 POSCAR 及对应 DFT 参考数据（`vaspout.h5` 格式）。实际应用中测试集应更大，约 500 个结构，视力场最终用途而定。

**流程与结果**：

1. 训练集误差：运行 refit 后 `grep ERR ML_LOGFILE`。ML_LOGFILE 的 ERR 段给出对 ab initio 数据的预测 RMSE：第 3 列为能量误差（eV/atom），第 4 列为力误差（eV/Å），末列为应力误差（kB）。本例结果：能量 RMSE ≈ $1.21\times10^{-4}$ eV/atom，力 RMSE ≈ $8.13\times10^{-3}$ eV/Å，应力 RMSE ≈ 0.61 kbar。
2. 测试集误差：`cp ML_FFN ML_FF`（把 refit 生成的新力场文件作为运行用力场），换用 `INCAR.run`，然后用脚本 `test_set_analysis.sh` 对 50 个测试结构逐个做 `ML_MODE = run` 单点（`mpirun -np 4 vasp_gam`），把每个构型的 `vaspout.h5` 收集到 `MLFF_data/`。
3. 用 py4vasp 做误差统计：`MLFFErrorAnalysis.from_files(dft_data="./test_set/DFTdata/*.h5", mlff_data="./MLFF_data/*.h5")`，再调用 `get_energy_error_per_atom()`、`get_force_rmse()`、`get_stress_rmse()` 得到逐构型误差并绘图；加 `normalize_by_configurations=True` 得到对整个测试集平均的总误差。本例总误差：能量 $1.26\times10^{-4}$ eV/atom，力 0.0085 eV/Å，应力 0.39 kbar。
4. **结论判读**：力的训练误差略低于测试误差，但两者都小于 10 meV/Å 且同量级 → 模型可靠，能够泛化到训练集之外的结构。

**复习问题**：(1) 若误差分析显示训练误差低而测试误差高，说明了什么？(2) `ML_MODE = refit` 时需要 POSCAR 吗？(3) refit 模式下线性回归问题是如何求解的？

### 示例 2：超参数扫描提升精度（ML_RCUT1）

**教学目标**：理解超参数（hyperparameter）的概念，掌握系统化超参数扫描的方法，通过调优进一步提升力场精度。

**原理与概念**：

- **超参数**是 MLFF 模型中由用户设定、在训练过程中不被优化的参数。恰当的超参数优化可以**同时**提升力场的精度与性能。典型例子是计算局域原子环境描述符所用的截断半径 `ML_RCUT1`（二体）与 `ML_RCUT2`（三体）。
- 教程给出最重要的超参数子集（完整列表见 VASP Wiki）：

| 参数 | 功能 | 类别 |
|---|---|---|
| `ML_RCUT1` | 二体截断半径 | 描述符 |
| `ML_RCUT2` | 三体截断半径 | 描述符 |
| `ML_WTIFOR` | 原子力的拟合权重 | 拟合 |
| `ML_WTOTEN` | 能量的拟合权重 | 拟合 |
| `ML_WTSIF` | 应力张量的拟合权重 | 拟合 |
| `ML_EPS_LOW` | 局域参考构型稀疏因子 | 稀疏 |
| `ML_RDES_SPARSDES` | 描述符稀疏因子 | 稀疏 |

- **方法学（交叉验证思想）**：先确认没有过拟合（训练误差不显著低于测试误差）；然后**系统地**（而非随机地）改变超参数——随机改动几乎不可能改进力场；对不同超参数组合反复比较训练/测试误差。为节省计算时间与磁盘，扫描阶段只考察训练集误差（从 ML_LOGFILE 读取），最后仅对训练误差最低的模型做完整测试集分析。

**输入要点**：`INCAR.refit` 在示例 1 的基础上增加 `ML_RCUT1 = <待扫描值>`，并设 `KSPACING = 300`（极大的 k 点间距，等效 Gamma 单点；教程原文个别处拼作 `KPACING`，系笔误）。本例对多个 `ML_RCUT1` 取值（教程扫描约 9–15 Å 范围，蓝色为教程预存数据、红色为现场计算点）分别 refit，得到一系列 `ML_FF_x.y` 与 `ML_LOGFILE_x.y`。

**流程与结果**：

1. 对每个 `ML_RCUT1` 值 refit，从各 ML_LOGFILE 提取 ERR 行，用 plotly 绘制能量/力/应力训练误差随截断半径的变化曲线。
2. **误差通道的取舍**：原则上三类误差越小越好，但实际应用中可能需要偏重其一——计算缺陷形成能时应最小化能量误差；计算声子频率需要最小的力误差；要精确预测体积则关注应力误差。偏重方式是在拟合中设置权重：`ML_WTOTEN`（能量）、`ML_WTIFOR`（力）、`ML_WTSIF`（应力）。
3. 从图中选出最优值 `ML_RCUT1 = 13.0`，保留对应的 `ML_AB_13.0` 与 `ML_FFN_13.0`，重复示例 1 的测试集分析流程（`ln -s ML_FF_13.0 ML_FF` + `test_set_analysis.sh` + py4vasp）。
4. 优化后测试集误差：能量 $7.56\times10^{-5}$ eV/atom，力 0.0079 eV/Å，应力 0.33 kbar——相比示例 1 的标准模型全面提升。
5. 实践建议：理论上应对 VASP 所有超参数重复此流程，但实际不可行；VASP 默认值已在多种块体材料与分子数据库上优化过，通常合理。若你的力场要支撑大量生产计算，则值得对上面表格中的参数逐一扫描。

**复习问题**：(1) 什么是超参数？(2) 相比示例 1，超参数优化使力误差改善了多少？(3) 超参数优化能否同时提升力场性能？

### 示例 3：稀疏化与性能优化（Pareto 前沿）

**教学目标**：基于现有 ML_AB 文件优化力场的运行性能——用尽可能小的精度损失换取尽可能大的加速。

**原理与概念**：

- 高效的**稀疏化（sparsification）**是提升力场性能的关键手段。VASP 可用 CUR 算法削减角向（三体）描述符的数量（开关 `ML_LSPARSDES`）。
- **CUR 算法**：把描述符协方差矩阵 $A$ 分解为三个矩阵的乘积 $A \approx CUR$，其中 $C$ 由 $A$ 的少数真实列组成、$R$ 由 $A$ 的少数真实行组成、$U$ 是保证乘积逼近 $A$ 的小幺正矩阵。原始 CUR 算法用于挑选矩阵中少数重要列，而 VASP 反其道而行之——基于 **leverage scores** 与稀疏因子 `ML_RDES_SPARSDES` 找出**不重要**的角向描述符并将其移除。参考文献：Mahoney 等（PNAS 106, 697-702, 2009）与 Jinnouchi 等（J. Chem. Phys. 152, 234102, 2020）。
- 移除描述符以精度换性能。寻找最优折中可借助 **Pareto 前沿**：x 轴为代码计时、y 轴为测试集精度，每个点对应一个 `ML_RDES_SPARSDES` 取值，从前沿上挑选最优点。对每个取值需要三步：refit → 测试集精度分析（同示例 1）→ 计时。
- **计时的两个要点**：(a) 力评估例程的计时必须对多个 MD 步取平均才有意义；(b) 文件写出通常很慢，必须在计时运行中消除——通过 `ML_OUTBLOCK` 与 `ML_OUTPUT_MODE` 控制。

**输入要点**：

- `INCAR.sparse`：`ML_MODE = refit`，保留上一步确定的 `ML_RCUT1 = 13.0`，打开 `ML_LSPARSDES = True`，设 `ML_RDES_SPARSDES = 0.5`。
- `INCAR.timing`：一段 100 步的 NVT MD（`IBRION = 0`、`NSW = 100`、`POTIM = 2.0` fs、`MDALGO = 3` Langevin 恒温器、`TEBEG = 300`、`LANGEVIN_GAMMA = 3.0 3.0`、`ISIF = 2` 恒容、`RANDOM_SEED = 10 0 0`）+ `ML_MODE = RUN`，并设 `ML_OUTBLOCK = 100`、`ML_OUTPUT_MODE = 0` 以抑制输出。

**流程与结果**：

1. 描述符数量对比：`grep NDESC ML_LOGFILE` 的最后一行给出每种元素的径向/角向描述符数。稀疏化前：径向 24、角向 452（Li 与 H 各自）；稀疏化后：径向 24 不变、角向降为 225——**角向描述符减半（还差 1 个到恰好一半）**，径向描述符不受影响。
2. 精度代价：稀疏化模型测试集误差为能量 $7.94\times10^{-5}$ eV/atom、力 0.0082 eV/Å、应力 0.37 kbar。三个模型力误差对比：标准 MLFF 8.5 meV/Å → RCUT 扫描后 7.9 meV/Å → 稀疏化后 8.2 meV/Å。结论：稀疏化使测试集误差**轻微**上升，但角向描述符减半。
3. 性能收益：对想要对比的各力场分别跑计时 MD（抑制 I/O），保存 OUTCAR，用 `grep LOOP+ OUTCAR | awk '{sum+=$7;n++} END {print sum/n}'` 对每步 CPU 时间取平均。教程预置结果：稀疏化 0.0111 s/步 vs 未稀疏化 0.0127 s/步。注意：本教程数据库很小；**数据库越大，稀疏化对计时的影响越显著**，理论最大加速因子即为 `ML_RDES_SPARSDES` 的取值。
4. 计时结果受同一计算节点上并行任务数的影响很大，做性能对比时要控制变量。

**复习问题**：(1) 如何找到性能与精度的最优折中？(2) `ML_RDES_SPARSDES` 的物理含义是什么？(3) `ML_OUTPUT_MODE = 1` 会写出哪些内容？

### 示例 4：用优化后的 MLFF 模型做生产运行

**教学目标**：把经过误差分析、超参数优化与稀疏化的 MLFF 投入实际生产——运行分子动力学并计算可与实验对比的物理量（对分布函数）。

**原理与概念**：

- 模型训练完毕后，所有参数都存储在 `ML_FF` 文件中，生产运行的 INCAR 只需 `ML_LMLFF = TRUE` 与 `ML_MODE = run` 两个 ML 标签。训练时使用的超参数记录在 ML_FF 的**第一行**（其余为二进制），可用 `head -n 1 ML_FF` 检查，例如能看到 `"ML_RCUT1" : 1.3000E+01`、`"descriptors" : [249, 249]`（Li、H 的描述符总数）、训练结构数 122、局域参考构型数 [737, 752] 等元数据。
- 生产任务：在 $T = 450$ K、$p = 1$ bar 条件下对 LiH 做 MD，计算**对分布函数（pair-distribution function, g(r)）**。对分布函数与静态结构因子互为傅里叶变换，而静态结构因子可直接从 X 射线/中子衍射等散射实验获得（参见 L. Van Hove, Phys. Rev. 95, 249, 1954）——因此 g(r) 提供了用实验数据验证力场的途径。
- 温控采用 **Langevin 恒温器**（`MDALGO = 3`）。

**输入要点**：

- `INCAR`：MD 参数 `IBRION = 0`、`NSW = 2000`、`POTIM = 2.0` fs；NVT 系综 `MDALGO = 3`、`TEBEG = 450`、`LANGEVIN_GAMMA = 0.1 0.1`、`ISIF = 2`；`RANDOM_SEED = 10 0 0`；`ML_LMLFF = .TRUE.`、`ML_MODE = RUN`、`ML_OUTPUT_MODE = 1`（确保写出 PCDAT）；`ISYM = -1`。
- `POSCAR`：已预平衡的 64 原子超胞（32 Li + 32 H），并包含初始速度。
- `ML_FF`：用符号链接复用示例 3 的稀疏化力场 `ln -s ../e03_*/ML_FF.sparse ML_FF`，避免复制大文件。
- `POTCAR`：纯 MLFF 运行中的哑文件，必须存在但内容不重要。

**流程与结果**：

1. `mpirun -np 4 vasp_std` 运行 2000 步 MD。
2. 用 py4vasp 验证运行质量：`Calculation.from_path("./e04_production")` 后 `my_calc.energy[:].plot("TOTEN temperature")`——能量应围绕均值涨落（能量守恒的表征）、温度应在 450 K 附近涨落，两者同时满足说明 MD 运行正常。
3. 画对分布函数：`my_calc.pair_correlation.plot()`。
4. **不用 py4vasp 的替代方案**：(a) 能量与温度可从 OSZICAR 提取——`grep T= OSZICAR | awk '{print $7}' > energy.out`（第 7 列为能量）、`awk '{print $3}' > temperature.out`（第 3 列为温度），再用 gnuplot 画图；(b) 教程提供 awk 脚本 `extract_PCDAT.sh` 解析 PCDAT 文件头部的标度参数（pcskal、pcfein、npaco）并把各时间片累加平均，输出 `PCDAT.xy`（r [Å] 与四组 g(r)），用 gnuplot `set xlabel "r [Å]"; set ylabel "g(r)"; p 'PCDAT.xy' w l` 绘图。

**本部分小结**（教程给出的完成清单）：计算并解读训练/测试集误差；微调超参数提升精度；对力场力评估计时；通过角向描述符稀疏化优化性能；用最终力场计算可与实验对比的物理量。

**本部分涉及的 VASP 功能与标签汇总**：`ML_LMLFF`、`ML_MODE`（refit/run）、文件体系 ML_AB（训练数据）/ML_FFN（refit 输出）/ML_FF（运行用力场）/ML_LOGFILE（ERR、NDESC、NDESC_SIC 段）；超参数 `ML_RCUT1`/`ML_RCUT2`/`ML_WTIFOR`/`ML_WTOTEN`/`ML_WTSIF`/`ML_EPS_LOW`/`ML_RDES_SPARSDES`；稀疏化 `ML_LSPARSDES`（CUR 算法、leverage scores）；输出控制 `ML_OUTBLOCK`/`ML_OUTPUT_MODE`；`KSPACING`；MD 相关 `IBRION = 0`/`MDALGO = 3`/`LANGEVIN_GAMMA`/`RANDOM_SEED`/`ISIF`；PCDAT 对分布函数；py4vasp 的 `MLFFErrorAnalysis` 与 `Calculation`；过拟合诊断、交叉验证与 Pareto 前沿方法学。

---

## 6. 表面（Surfaces）

真实块体总有边界——表面。教程通过沿一个方向把原胞延伸到真空中构建表面/板（slab）模型：表面结构弛豫、态密度与能带、分子吸附以及 STM 图像模拟。

## 6.1 Part 1：镍表面（表面能、局域态密度、表面能带）

> 来源：https://vasp.at/tutorials/latest/surface/part1/ ｜ 练习文件：`surface-part1.zip` ｜ 示例目录：`$TUTORIALS/surface/e01_Ni100-relaxation`、`e02_Ni100-DOS`、`e03_Ni100-band`、`e04_Ni111-relaxation`

本部分以 fcc 镍为对象，系统讲授 VASP 中的表面（slab）建模与表征：表面弛豫与表面能、逐层局域态密度与局域磁矩、表面能带与带字符分析，以及密排面 (111) 的对比研究和 Methfessel-Paxton 占据展宽方法。四个示例环环相扣——示例 1 的 CONTCAR 作为示例 2、3 的结构输入，示例 1 的 CHGCAR 作为示例 3 的电荷密度输入。

### 示例 1：Ni(100) 表面的弛豫与表面能

**教学目标**：(1) 用选择性动力学（selective dynamics）模式配合 RMM-DIIS 算法弛豫指定原子层；(2) 掌握表面计算的 k 点与单胞设置；(3) 计算单原子体系的表面能。

**原理：slab 建模**：表面是任何块体材料的边界。VASP 中通过沿某一方向拉长单胞来建模表面——单胞的一部分填材料、另一部分留真空，构成表面平板（slab）。两个关键约束：

1. 拉长方向只取 **1 个 k 点**：该方向上任何额外的 k 点都等价于引入与周期性镜像单胞的相互作用，破坏"孤立表面"的物理图像。
2. 材料层与真空层都必须**足够厚**：使 slab 中部能呈现块体性质、相邻镜像 slab 之间无相互作用。换言之，该方向上不能依赖周期性边界条件，必须主动消除其影响。

**输入要点**（`e01_Ni100-relaxation`）：

- `POSCAR`：fcc Ni（晶格常数 3.53 Å）沿 (100) 切出的 5 层 slab，面内基矢 $(\pm0.5, 0.5, 0)\times a$，z 向基矢拉长至 $5a$。每层 1 个原子，共 5 个原子；真空厚度 = 3 个晶格周期 = 10.59 Å。启用 `Selective dynamics`：底部 3 层标记 `F F F`（固定，模拟块体环境），顶部 2 层标记 `T T T`（允许弛豫，模拟表面重构）。
- 块体参考计算（`bulk/` 子目录）：常规 fcc 单胞、1 个原子。
- `INCAR`：`ISPIN = 2` + `MAGMOM = 5*1`——自旋极化的共线自旋密度泛函（SDFT）计算，两个自旋分量在包含双方电荷密度的有效哈密顿量下分别求解 Kohn-Sham 轨道；`ICHARG = 2` 初始电荷密度取原子电荷叠加；`ENCUT = 270`；`ISMEAR = 2`、`SIGMA = 0.2`——Methfessel-Paxton 二阶占据展宽（适合金属弛豫）；`ALGO = Fast`；`EDIFF = 1E-6`；几何弛豫 `IBRION = 1`（RMM-DIIS）、`NSW = 100`、`POTIM = 0.8`、`ISIF = 2`（弛豫离子位置、不变晶胞形状与体积）。
- `KPOINTS`：slab 用 Monkhorst-Pack 9×9×1（奇数网格自动 Γ 中心），块体参考用 9×9×9。

弛豫的算法循环（伪代码）：① 给定离子位置，在 SDFT 下求电子基态；② 由 Hellmann-Feynman 定理计算作用在离子上的力与应力（电子哈密顿量梯度的期望值）；③ 用 RMM-DIIS 把离子向瞬时基态移动；④ 重复直到满足收敛判据。

**流程与结果**：

1. 分别运行 slab 与块体计算（`mpirun -np 2 vasp_std`）。用 py4vasp `mycalc.structure.plot(3)` 可视化结构（传 3 绘制 3×3×3 超胞网格；右键两次可测量原子间距）。
2. 弛豫结果（CONTCAR）：顶层双层发生明显**向内弛豫**——第 1、2 层间距相对块体值变化 $\Delta d_{12} = -4.70\%$，第 2、3 层间距变化 $\Delta d_{23} = +0.21\%$。OUTCAR 的 TOTAL-FORCE 段显示第 4、5 层受力已降至约 $10^{-3}$ eV/Å 量级；未显式设置时 `EDIFFG = 10*EDIFF = 1e-5`。
3. **表面能计算**。公式

$$\sigma = \frac{E_\mathrm{slab} - N E_\mathrm{bulk}}{2A}$$

仅对**对称且化学计量**的 slab 成立：因子 $1/2$ 是因为切开块体总是产生两个表面，只有两者等价（slab 对称）时公式才有效；$N$ 表示 slab 包含块体化学式单元的个数（要求 slab 化学计量比与块体一致）。
4. 教程采用"部分弛豫三步法"提取弛豫表面能：未弛豫表面能 $\sigma^\mathrm{unrel} = \frac{1}{2A}(E_\mathrm{slab}^\mathrm{unrel} - 5E_\mathrm{bulk}) = 143.0$ meV/Å$^2$（其中 $A = 3.53^2/2$ Å$^2$，$E_\mathrm{slab}^\mathrm{unrel} = -25.5567$ eV）；弛豫顶两层后的平均表面能 $\sigma^\mathrm{avg} = 141.5$ meV/Å$^2$；近似的完全弛豫表面能 $\sigma^\mathrm{rel} = 2\sigma^\mathrm{avg} - \sigma^\mathrm{unrel} = 140.1$ meV/Å$^2$。
5. **解读**：弛豫能仅约 3 meV/Å$^2$——对 fcc Ni 这样的密排金属很小，但其他材料（如半导体、氧化物表面）的弛豫可能大得多，方法学不变。文献 PBE 值为 138 meV/Å$^2$（Sci. Data 2016, 3, 160080），说明部分弛豫法以极小代价给出了好结果。
6. 物理图像：表面能为正——原子在能量上偏好块体配位环境，形成表面需要付出能量；半局域交换关联泛函通常使表面比实验更"稳定"（模拟表面能偏小）；若 (111) 表面能低于 (211)，说明沿 (211) 解理需要更多能量、该晶面更难暴露。

**复习问题**：(1) VASP 中如何建模表面？(2) 垂直表面方向取几个 k 点合适？为什么？(3) 什么是表面能？单原子 slab 如何计算？(4) `EDIFFG` 控制什么？

### 示例 2：Ni(100) 表面的局域态密度（LDOS）

**教学目标**：计算局域量（逐原子电荷、局域自旋磁化强度）；计算并绘制局域态密度。

**原理：PAW 框架下获取局域量的三条路线**——要回答"每个 Ni 原子上有多少电荷"，必须先回答"哪部分电荷密度属于哪个原子"：

1. **Wigner-Seitz 半径法**：用 `RWIGS` 指定属于每个原子的球区（POTCAR 提供默认值），`LORBIT < 10` 时 DOS 按此划分。缺点：体积划分含糊，球外电荷无归属。
2. **Wannier 函数法**：对（部分）Kohn-Sham 轨道做幺正变换构造局域态，如最大局域 Wannier 函数（MLWF）。缺点：当不同字符的能带强烈纠缠时，幺正变换难以构造。
3. **投影算符法**：用投影算符把 KS 态投影到由局域量子数（角动量 $l$、磁量子数 $m$、自旋 $\sigma$）标记的局域基上，即投影局域轨道（Phys. Rev. B 77, 205112, 2008）；`LORBIT >= 10` 时 DOS 采用此方案，投影权重与相位因子写入 `PROCAR` 文件。

**输入要点**（`e02_Ni100-DOS`）：

- `POSCAR` 直接取示例 1 的 CONTCAR（弛豫后结构）。
- `INCAR`：`ISMEAR = -5`——带 Blöchl 修正的四面体方法。**注意**：该方法对部分占据不是变分的，会给出错误的力与应力，因此**绝不能用于金属的几何优化**，只适合静态 DOS/能带计算；`ALGO = Normal`；`ISPIN = 2`、`MAGMOM = 5*1`；`LORBIT = 11`——使 DOSCAR 与 PROCAR 写出逐原子、逐 $(l, m)$ 的分波信息。
- `KPOINTS`：13×13×1（DOS 计算加密面内网格）。

**流程与结果**：

1. 运行后用 py4vasp 绘制逐层 DOS：`mycalc.dos.plot(selection="3, 5")` 对比最"块体化"的第 3 层与弛豫表面层（第 5 层）。核心结论：**表面层的投影 DOS 显示更大的交换劈裂**（exchange splitting）——表面 d 带变窄、交换作用增强。
2. OUTCAR 末尾给出按层的 total charge 与 magnetization (x)（s/p/d 分解）：表面层（1、5）的总电荷（~9.11/9.15）低于内层（~9.28-9.32），但 d 磁矩更高（0.751/0.744 vs 内层 ~0.67）——**表面层电荷密度降低而磁化增强**，这是 Ni 表面的经典结论（配位数减少 → d 带变窄 → Stoner 判据更易满足）。
3. 用 `mycalc.local_moment.plot()` 可视化每个位点的自旋磁矩箭头大小。
4. 对比实验：把 INCAR 改为 `LORBIT = 1` + `RWIGS = 1.286` 重新计算（Wigner-Seitz 球法），OUTCAR 会报告球区占单胞体积的比例（40.5%），逐层电荷/磁矩数值发生变化——说明局域量的绝对值依赖于划分方案，比较时应固定同一方案。

**复习问题**：(1) PAW 方法中计算局域量有哪些策略？(2) `RWIGS` 控制什么？(3) `PROCAR` 文件写入哪些信息？

### 示例 3：Ni(100) 表面能带（带字符分析）

**教学目标**：对比块体与表面计算的 Wigner-Seitz 原胞（布里渊区）；绘制带字符（band character）。

**原理**：slab 在实空间沿 z 拉长单胞，倒空间布里渊区随之沿对应方向**被压扁**，高对称点标签也随之改变。教程用 SeeK-path 工具（上传 POSCAR 即可）分别生成 fcc 块体的第一布里渊区（截角八面体）与 slab 的布里渊区：

![fcc 块体与 Ni(100) slab 的布里渊区对比](/imgs/2026-08-05/fig20.png)

由于定义布拉维格子的基矢变了，把块体与表面的高对称点相互对应时必须格外谨慎（图中已给出两者的 Γ-X-M-Γ 对应关系）。

**输入要点**（`e03_Ni100-band`）：

- 结构仍用示例 1 弛豫后的 CONTCAR；另需把示例 1 的 `CHGCAR` 复制到本目录。
- `INCAR`：关键是 `ICHARG = 11`——从 CHGCAR 读取电荷密度并**在整个运行中保持不变**（非自洽能带计算），OUTCAR 会打印 "Static calculation / charge density remains constant during run"；`ISMEAR = 2`、`SIGMA = 0.2`；`ISPIN = 2`、`MAGMOM = 5*1`；`LORBIT = 11`（为带字符投影提供 PROCAR 信息）。
- `KPOINTS`：line-mode 沿表面二维布里渊区路径 Γ(0,0,0) → X(0,0.5,0) → M(0.5,0.5,0) → Γ，每段 15 个点；z 方向无色散，故所有 k 点的第三分量为 0。

**流程与结果**：`cp ../e01_*/CHGCAR .` 后运行 `vasp_std`；用 py4vasp 画能带并按轨道字符着色：`mycalc.band.plot(selection="Ni(s, p, d)")`——每条能带按其 s/p/d 投影权重着色，可直观识别表面附近的 d 带特征与可能的表面态。

**复习问题**：(1) 与块体计算相比，表面计算的 Wigner-Seitz 原胞如何变化？(2) 什么是能带的"字符"（character）？

### 示例 4：Ni(111) 表面弛豫（Methfessel-Paxton 方法深入）

**教学目标**：解释 Methfessel-Paxton 方法的目的与应用；用 py4vasp 测量结构中的距离与键角。

**原理：金属的布里渊区积分难题**——电荷密度等量需要对第一布里渊区做数值积分，而金属在 $T = 0$ 的占据数是阶跃函数：同一能带在某些 k 点被占据、在另一些 k 点空置（部分占据的导带），被积函数不连续，需要极细的 k 网格才能收敛，代价巨大。`ISMEAR` 控制的各类展宽方法即为此而生。**Methfessel-Paxton 方法**（Phys. Rev. B 40, 3616, 1989）用阶跃函数的光滑近似函数替代，这些近似函数被构造为对指定阶数的多项式积分给出精确结果；阶数由 `ISMEAR = N > 0` 设定（本教程用二阶）。

**输入要点**（`e04_Ni111-relaxation`）：

- `POSCAR`：fcc Ni 沿 (111) 面切出的 5 层 slab（面内六角基矢，z 向基矢长度 5.774 Å 量级），真空约 $(1-0.4)\times 5.774\times 3.53 \approx 12.23$ Å；选择性动力学同样固定底 3 层、放开顶 2 层。
- `INCAR`：与示例 1 基本相同，但**关闭自旋极化**以节省时间，并新增：`NBANDS = 34`（高于默认值以改善收敛）、`LMAXMIX = 6`、`SIGMA = 0.15`、`NSW = 20`。
- `KPOINTS`：11×11×1（块体参考 11×11×11）。

**流程与结果**：

1. 运行 slab 与 bulk 计算；用 py4vasp `mycalc.structure.plot()` 可视化并测量（右键一个原子一次 + 另一个原子两次测距离；右键三个原子且最后一个两次测键角）。
2. 表面面积：(111) 面内为边长 $a' = 2.5$ Å、夹角 120° 的平行四边形，面积 $A = (a')^2 \sin(120°) = 5.41$ Å$^2$（等价于两个等边三角形面积 $a'^2\sqrt{3}/4$）——比 Ni(100) 更致密。
3. 表面能（同样的三步法）：$\sigma^\mathrm{unrel} = 122.7$ meV/Å$^2$；$\sigma^\mathrm{avg} = 122.1$ meV/Å$^2$；$\sigma^\mathrm{rel} = 2\sigma^\mathrm{avg} - \sigma^\mathrm{unrel} = 121.6$ meV/Å$^2$。
4. **对比结论**：(111) 面的初始受力就很小，弛豫后结构更接近块体；其弛豫能增益比 (100) 更小（层间距更小、表面更密排）；(111) 是 fcc Ni 所有表面中表面能最低的晶面。尽管本例忽略了磁性，与文献值 ~120 meV/Å$^2$ 的偏差依然很小。
5. **教程的三点告诫**：(a) 务必对 k 网格密度、ENCUT、以及表面计算特有的**真空厚度与 slab 厚度**做收敛性研究；(b) OUTCAR 中较大的 external pressure 是 **Pulay 应力**，可通过提高 ENCUT 消除——且比较两个计算的总能时必须使用相同 ENCUT，因此严格计算应在更高截断能下进行；(c) 教程直接把块体晶格常数当作已知，实际工作中晶格常数依赖交换关联泛函，应先对块体做结构弛豫。

**复习问题**：(1) Methfessel-Paxton 方法解决什么问题？(2) 如何在 py4vasp 结构图中测量两原子间距？(3) 如何验证 slab 厚度对表面能计算已经足够？

**本部分涉及的 VASP 功能与标签汇总**：slab 建模（真空层、垂直方向单 k 点、板中部块体化）、`Selective dynamics` 分层弛豫、`IBRION = 1`（RMM-DIIS）、`EDIFFG` 默认值（10×EDIFF）、表面能公式与"部分弛豫三步法"、`ISPIN`/`MAGMOM` 共线自旋极化、`ISMEAR = 2`（MP 展宽）与 `ISMEAR = -5`（四面体法，禁用于金属弛豫）、局域量计算（`RWIGS`、`LORBIT = 1/11`、PROCAR、投影算符法）、`ICHARG = 11` 非自洽能带计算、slab 布里渊区与高对称点（SeeK-path）、带字符绘图、`NBANDS`/`LMAXMIX`、Pulay 应力与 ENCUT 一致性、py4vasp 结构/DOS/能带/局域磁矩分析与几何测量。

---

## 6.2 Part 2：CO 吸附在 Ni(111) 表面（偶极修正、吸附能、功函数、PDOS 与振动）

> 来源：https://vasp.at/tutorials/latest/surface/part2/ ｜ 练习文件：`surface-part2.zip` ｜ 示例目录：`$TUTORIALS/surface/e05_CO_on_Ni111_relaxation`、`e06_CO_on_Ni111_adsorption`、`e07_CO_on_Ni111_pDOS`、`e08_CO_on_Ni111_Freqs`

本部分是表面催化计算的标准全流程演练：吸附结构弛豫（含偶极修正）→ 吸附能与功函数 → 分波态密度与 Blyholder 成键模型验证 → 吸附分子振动频率。体系为 CO 分子吸附在 Ni(111) 表面的顶位（ontop），覆盖度 100%（每个表面原胞一个 CO）。

### 示例 5：CO/Ni(111) 表面弛豫（100% 覆盖度与偶极修正）

**教学目标**：(1) 弛豫 100% 覆盖度的表面吸附分子；(2) 理解并应用偶极修正（dipole corrections）；(3) 评估吸附引起的 slab-吸附物几何变化。

**原理**：

- 分子在表面的吸附是表面科学的核心问题：即使对简单分子，构型空间也很大，高维势能面往往很平坦；分子可以不同取向吸附，平衡键角也会因与表面的相互作用和电荷转移而改变。CO 吸附是催化反应性的经典测试案例（有毒 CO 常需氧化为 CO₂）。
- **偶极修正的必要性**：吸附物会显著改变 slab 的电荷分布，使体系产生有限偶极矩。在周期性边界条件下，这个偶极矩会在真空区产生电势梯度（人为电场），而我们想模拟的是孤立 slab-吸附物体系。偶极修正（`LDIPOL`）对势、能量和力做修正，**推荐用于一切破坏反演对称、允许偶极矩存在的体系**。
- 为什么分两步做（先无修正电子最小化，再开修正弛豫）：修正通过在势中引入一个不连续跳跃形成"锯齿"来屏蔽偶极的长程相互作用，跳跃大小依赖偶极矩本身，过早引入会拖慢 KS 本征值收敛；但若不加修正，总能量对真空厚度 $L$ 呈 $1/L^2$ 依赖，真空收敛被破坏，等效于模拟"层状 slab 堆成的块体"。

**输入要点**（`e05_CO_on_Ni111_relaxation`）：

- 结构用 **ASE**（Python 原子模拟环境）构建：`fcc111('Ni', (1,1,5), a=3.53)` 生成 5 层 Ni(111) slab，`add_adsorbate(slab, adsorbate, 1.8, 'ontop')` 把 CO（初始键长 1.142 Å）放到顶位、初始高度 1.8 Å，z 向留 4 Å 真空。fcc(111) 表面有多种高对称吸附位（顶位 ontop、桥位 bridge、fcc/hcp 空位），本例选顶位——高覆盖度下最稳定的吸附位（Surf. Sci. 526, 332-340, 2003）；模拟低覆盖度则需要超胞。

![Ni(111) 表面高对称吸附位示意](/imgs/2026-08-05/fig21.png)

- `POSCAR`：7 原子（5 Ni + C + O），`Selective dynamics` 固定底部 3 层 Ni，放开顶部 2 层 Ni 与 C、O（z 向）。
- `DIPOL` 应设在电荷中心；事先未知时可用质心近似：ASE `slab.get_center_of_mass(scaled=True)` 给出 `DIPOL = 0.5 0.5 0.453`。
- 两阶段 INCAR：`INCAR.1`（无偶极修正的预备电子最小化：`ENCUT = 400`、`ISMEAR = 0`、`SIGMA = 0.15`、`ALGO = Fast`、`EDIFF = 1E-4`、`AMIX = 0.4`、`LMAXMIX = 6`）；`INCAR.relax`（`ISTART = 1` 从上一步波函数续算，`EDIFF = 1E-6`，弛豫 `IBRION = 1` RMM-DIIS、`NSW = 20`、`POTIM = 0.8`、`EDIFFG = -0.05`（力判据 eV/Å），偶极修正 `LDIPOL = .TRUE.`、`IDIPOL = 3`（仅沿 z 方向）、`DIPOL = 0.5 0.5 0.453`）。
- 孤立 CO 分子（`CO/` 子目录）：8×8 Å 截面、z 向约 9.1 Å 的盒子，`IBRION = 2`（共轭梯度）、`POTIM = 0.1`、`NSW = 10`，同样开偶极修正（`DIPOL = 0.5 0.5 0.509`）；Gamma 单点。弱极性分子从一开始就引入修正通常也能收敛。
- 细节讨论：(a) `ENCUT = 400` 主要由 O 的 POTCAR 中 ENMAX 决定；小分子建议测试硬势（C_h、O_h）是否改变结果——短键长用硬势描述更好但更贵；(b) Ni(111) 有磁性，严格应做 `ISPIN = 2` 计算，但 CO 无磁矩，非磁计算是良好近似（教程为演示流程采用非磁）；(c) 打破反演对称的长单胞容易出现 charge sloshing，收敛困难时调整电荷密度混合参数（AMIX 等）。

**流程与结果**：

1. 两阶段运行：`cp INCAR.1 INCAR; vasp_std`，然后 `cp INCAR.relax INCAR; vasp_std`。OUTCAR 中搜 `dipolmoment` 与 `dipol+quadrupol energy correction`——本例二者都很小（金属 slab + 弱极性吸附物）；若吸附物是带电离子（如 K）或自由基氧，偶极矩会大得多。
2. 几何分析（py4vasp）：CO-Ni 距离 1.762 Å；顶部两层 Ni 相对干净表面**向外弛豫 3.17%**（与 Part 1 示例 4 干净 (111) 面的微弱弛豫对比）。注意：示例 4 与本例 ENCUT 不同——比较**总能**必须同 ENCUT，但比较几何量（层间距）只要求各计算自身对该量收敛即可。
3. CO 键长变化：吸附态 1.154 Å vs 孤立分子 1.142 Å，**键长增大约 1%**——CO 与表面成键削弱了 C-O 键（与后文 back-donation 图像一致）。

**复习问题**：(1) fcc(111) 面有多少个对称不等价的高对称吸附位？(2) 电荷密度的偶极矩对真空势有什么影响？(3) 列举吸附计算需要收敛性研究的四个量。(4) 比较干净 slab 与吸附体系的层间距是否必须同 ENCUT？

### 示例 6：吸附能与 Ni(111) 功函数

**教学目标**：计算分子在表面的吸附能；计算表面功函数；可视化平均局域势。

**原理**：

- 吸附能定义：

$$E_\mathrm{ads} = E_\mathrm{total} - E_\mathrm{slab} - E_\mathrm{molecule}$$

其中三项分别是吸附体系、干净 slab、孤立分子的总能；各体系结构都应先弛豫。$E_\mathrm{ads}$ 依赖吸附位、覆盖度、晶面、化学组成等。
- **ENCUT 一致性**：Part 1 示例 4 的干净 slab 与示例 5 的吸附体系用了不同 ENCUT。各自性质可能都已收敛，但**比较总能必须把所有计算的 ENCUT 手动设为同一值**——因此本例把示例 4 的 CONTCAR 复制为 POSCAR，用示例 5 的 ENCUT=400 重算干净 slab 的总能。
- **功函数** $\phi$：真空能级与费米能级之差。计算依赖真空区定义良好的势：沿选定方向做线平均，通过寻找势为常数的区域识别真空区。干净金属表面本身几乎不积累电荷，但（小但有限的）偶极矩会干扰该搜索算法，故本例仍开偶极修正。输入中新增 `WRT_POTENTIAL = hartree ionic` 与 `LVHAR = .TRUE.`——把 Hartree 势（含离子贡献）写入文件供分析。

**输入要点**（`e06_CO_on_Ni111_adsorption`）：POSCAR = 示例 4 弛豫后的干净 Ni(111) slab；INCAR：`ENCUT = 400`、`ISMEAR = 0`、`NBANDS = 34`、`LDIPOL/IDIPOL/DIPOL`（0.5 0.5 0.199）、`WRT_POTENTIAL = hartree ionic`、`LVHAR = .TRUE.`；KPOINTS 11×11×1。

**流程与结果**：

1. 吸附能：用 py4vasp 读取三个计算最后一步的 `energy.to_numpy('without_entropy')` 相减——$E_\mathrm{ads} \approx -292$ meV（即约 0.29 eV 的吸附强度）。**泛函表现**：LDA 严重高估吸附能；半局域 GGA 略有改善但仍高估（Phys. Rev. B 59, 7413, 1999）。一个耐人寻味的对比：半局域泛函**低估**表面形成能、却**高估**吸附能。
2. 功函数：`mycalc.workfunction` 给出沿晶格矢量 3 的分析——两侧真空势 7.165/7.181 eV、费米能 1.997 eV，$\phi \approx 5.18$ eV；生产级计算（更厚 slab、更精确设置）给出 5.10 eV（Surf. Sci. 687, 2019）。`mycalc.workfunction.plot()` 可视化势线平均。
3. 两个关键问答：(a) 为什么用 `LVHAR`（仅 Hartree 势）而不是 `LVTOT`（含 xc 势）求真空能级——xc 势在真空区的数值评估有问题，物理上其对真空势的贡献应趋于零；所有 xc 效应已通过自洽场与费米能完全计入，故不含 xc 更准确。(b) slab 两侧真空势为何略有差异——slab 不对称：下表面未弛豫、上表面弛豫。

**复习问题**：(1) 表面计算何时应使用偶极修正相关标签？(2) 确定真空能级时应否包含交换关联势？为什么？

### 示例 7：CO/Ni(111) 的分波态密度（Blyholder 模型验证）

**教学目标**：计算吸附体系的分波 DOS（PDOS），在电子结构层面观察 CO 与金属表面的电荷捐赠/反馈（donation/back-donation）。

**原理：Blyholder 模型**（J. Phys. Chem. 68, 2772, 1964）——CO 的 $5\sigma$ 是最高占据分子轨道（HOMO）、$2\pi^*$ 是最低未占分子轨道（LUMO）。二者与金属表面相互作用形成成键/反键杂化轨道：反键的 CO $5\sigma$-金属杂化轨道位于费米能级之上 → 电子从 CO **捐赠（donation）**给表面；成键的 $2\pi^*$-金属杂化轨道被占据 → 表面向 CO **反馈（back-donation）**。本示例的目标就是在 PDOS 中观察这一电荷转移。

**输入要点**（`e07_CO_on_Ni111_pDOS`）：从示例 5 的 CONTCAR + CHGCAR 续算（`ICHARG = 1` 从 CHGCAR 读电荷密度），k 网格加密为 13×13×1（Γ 中心，DOS 更平滑）；INCAR 的 DOS 关键标签：`ISMEAR = -5`（四面体法+Blöchl 修正，推荐用于高精度总能与 DOS）、`NEDOS = 5001`、`EMIN = -10`、`EMAX = 10`、`LORBIT = 11`（用 PAW 投影把 KS 轨道分解到局域量子数，给出位点投影 DOS）、`NBANDS = 48`；偶极修正与 `LVHAR` 保留。

**流程与结果**（py4vasp 分析）：

1. **$5\sigma$ 通道**：读 `C(pz) O(pz) C(s) 5(d)` 的投影 DOS，组合 C 的 s+pz 与 O 的 pz、表面第 5 层 Ni d。可见 CO $5\sigma$ 态发生劈裂、与表面形成成键与反键态；费米能级上方有小峰——$5\sigma$ 吸附后未完全占据，对应向表面的电荷捐赠。
2. **平滑技巧**：表面 Ni d 的 DOS 噪声大（要真正光滑需极细 k 网格），教程用 Hann 窗卷积做平滑帮助辨认特征——但卷积可能引入假象或抹掉真实峰，**不应作为标准后处理，须逐案决定**。
3. **$\pi$ 通道**：读 `C(px) C(py) O(px) O(py) 3(d)`。$1\pi$ 态以 O 为主、$2\pi^*$ 以 C 为主。对比表面层（第 5 层）与"块体层"（第 3 层）Ni d：表面态在 -8.5 eV 与 -6.2 eV 附近形成额外权重、在 ~1 eV 出现未占态，使费米能级附近的 Ni d 表面态耗尽。费米能级以上以 C 为主的 CO 态即 CO-$2\pi^*$（back-donation 通道）；费米能级以下以 O 为主的是 CO-$1\sigma$ 态，它们形成完全占据的成键/反键组合，把 Ni d 态推离费米能级、使表面"钝化"。
4. **轨道分辨**：读 `1(dz2), 1(dxy), 3(dz2), 3(dxy), 5(dz2), 5(dxy)`——Ni $d_{3z^2-r^2}$（沿表面法线取向）与 CO 分子杂化明显；$d_{xy}$ 位于 slab 平面内，几乎不受影响。
5. 顺带重算功函数，对比示例 6 的干净表面：CO 吸附物通过改变局域势（表面偶极）使功函数发生变化——这正是吸附改变表面电子性质的直观证据。

**复习问题**：教程要求复习 `LORBIT` 与 `ISMEAR` 的取值理由；思考 donation/back-donation 在 PDOS 中的谱学指纹。

### 示例 8：CO/Ni(111) 的振动频率（有限差分法）

**教学目标**：用有限差分法计算表面分子的 Γ 点振动频率；通过选择性动力学裁剪分子的自由度。

**原理**：晶格振动是表面热激发最重要的贡献，振动谱对理解表面热力学性质至关重要；与高分辨电子能量损失谱（HREELS）、拉曼光谱等实验手段结合，表面声子计算能揭示吸附物的成键情况（综述：Wöll, Appl. Phys. A 53, 377, 1991）。VASP 的 `IBRION = 5` 用**有限差分法**：把每个允许自由度上的原子位移 ±POTIM，用 Hellmann-Feynman 力构建动力学矩阵并对角化。

**输入要点**（`e08_CO_on_Ni111_Freqs`）：

- `POSCAR`：全部 Ni 固定（`F F F`），C 与 O 标记 `F F T`——**只允许 z 方向运动**，即只计算垂直表面的振动模式。
- `INCAR`：`ICHARG = 1`（从示例 5 的 CHGCAR 续算）；频率计算精度要求：`EDIFF = 1E-7`（比默认 $10^{-4}$ 严格得多，确保力精确）、`PREC = Accurate`、`SIGMA = 0.05`；有限差分核心标签：`IBRION = 5`、`POTIM = 0.015`（位移幅度，单位 Å）、`NFREE = 2`（每个自由度正负两个位移）；偶极修正保留。
- 两个关键问答：(a) 自由度计数——C、O 各只能沿 z 动：同向运动 = CO-金属伸缩振动，反向运动 = C-O 伸缩模式，共 **2 个频率**；`NFREE = 2` 时每原子每方向 2 个位移，加上平衡位置的 1 次计算，共 **5 次静态计算**。(b) POTIM 为何这么小——频率计算假设谐振近似成立，位移太大时恢复力不再与位移成正比（非谐效应污染结果）。注意 `IBRION = 5` 会考虑选择性动力学设置（只对标记 T 的自由度位移）。

**流程与结果**：

1. `cp ../e05_*/CHGCAR .` 后运行。在 OUTCAR 中搜 `Eigenvectors and eigenvalues of the dynamical matrix`。
2. 两个模式：(a) **C-O 伸缩模式** $f = 64.0$ THz = 2135 cm$^{-1}$ = 264.7 meV，本征矢量显示 C（dz = -0.78）与 O（dz = +0.62）反向运动；(b) **CO-金属伸缩模式** $f = 12.1$ THz = 404.5 cm$^{-1}$ = 50.2 meV，C、O 同向运动（dz 同号）。Ni 原子的位移分量为零（被冻结）。
3. 拓展练习：把 C、O 的选择性动力学改为 `T T T`（放开全部方向）重算，会多出 4 个模式——两个是 C、O 平行表面反向运动（分子在 x、y 方向来回"倾斜"，频率略高于垂直集体运动）；另两个是 x、y 方向的集体平动，频率不到前者的 10%。

**复习问题**：(1) `IBRION = 5` 是否考虑选择性动力学？(2) 有限差分步长过大会有什么问题？

**本部分涉及的 VASP 功能与标签汇总**：ASE 结构构建与吸附位选择、`LDIPOL`/`IDIPOL`/`DIPOL` 偶极修正（原理与两步法策略）、吸附能公式与 ENCUT 一致性、功函数计算（`LVHAR`、`WRT_POTENTIAL`、真空能级与 xc 势取舍）、`ICHARG = 1` 续算、`ISMEAR = -5` 四面体法 DOS、`NEDOS`/`EMIN`/`EMAX`、`LORBIT = 11` 位点与轨道分辨 PDOS、Blyholder donation/back-donation 模型、DOS 卷积平滑的利弊、`IBRION = 5`/`NFREE`/`POTIM` 有限差分频率计算、选择性动力学裁剪振动自由度、OUTCAR 动力学矩阵本征值输出。

---

## 6.3 表面科学（三）：STM 扫描隧道显微镜模拟（恒定高度 / 恒定电流模式）

> 来源：https://vasp.at/tutorials/latest/surface/part3/ ｜ 练习文件：`surface-part3.zip` ｜ 示例目录：`$TUTORIALS/surface/e09_STM_of_graphite`、`e10_STM_of_graphene`

本部分演示如何用 VASP 计算**能带分解的部分电荷密度（band-decomposed partial charge density）**，再结合 py4vasp 内置的 Tersoff-Hamann 模型模拟 STM（扫描隧道显微镜）图像。两个示例分别覆盖 STM 的两种实验模式：石墨 (0001) 表面的恒定高度模式（示例 9）与石墨烯的恒定电流模式（示例 10）。

**STM 模拟的物理基础**：STM 利用量子隧穿在原子尺度上成像干净的（半）导电表面——可以在固定针尖高度下测电流，也可以在保持隧穿电流恒定的同时追踪针尖位置。测量信号是针尖位置、偏压与局域态密度的函数，因此图像中亮暗区域对应的表面结构并不总是直观的；模拟 STM 图像是解释实验结果、研究表面结构与重构乃至动力学的有效手段。VASP+py4vasp 采用 **Tersoff-Hamann 方法**（Phys. Rev. Lett. 50, 1998, 1983，基于 Bardeen 隧穿理论，Phys. Rev. Lett. 6, 57, 1961）：把针尖建模为无结构的点（单个 s 类波函数），从而把表面对隧穿电流的贡献与针尖的贡献解耦（后者被忽略）。该近似对许多表面（尤其是较大针尖-样品距离）成立；需要更精确结果时可用 Chen 导数规则（Phys. Rev. B 42, 8841, 1990），但必须显式处理针尖。

### 示例 9：石墨 (0001) 表面的恒定高度 STM 模拟

**教学目标**：(1) 用 VASP 后处理计算石墨 slab 的部分电荷密度；(2) 在不同针尖高度下模拟并绘制恒定高度 STM 图像。

**任务流程**：先做静态计算收敛电荷密度，把得到的 WAVECAR 作为后处理步的输入，计算费米能级附近能带的部分电荷密度，最后绘制模拟 STM 图像。本例模拟偏压 -250 mV 的恒定高度图像——即在表面上方固定距离处绘制平滑化的部分电荷密度。

**输入要点**（`e09_STM_of_graphite`）：

- `POSCAR`：5 层石墨 slab（10 个 C 原子），采用实验晶格常数（a = 2.46 Å）与层间距，z 向单胞约 26.7 Å。
- 静态计算 `INCAR`：非自旋极化静态计算，`ENCUT = 400`、`ISMEAR = 0`、`SIGMA = 0.2`、`ALGO = Fast`、`EDIFF = 1E-6`，另设 `NEDOS = 3001`、`LORBIT = 11` 便于 DOS 分析。KPOINTS：Γ 中心 17×17×1。
- **PREC 细节**：STM 模拟推荐 `PREC = Single`——此时精细 FFT 网格与粗网格一致（NGX = NGXF），电荷在真空中不做插值，避免真空区电荷密度出现高频振荡；教程保留默认 `PREC = Normal`（分辨率翻倍），由 py4vasp 在模拟时用高斯平滑去除潜在振荡。
- 部分电荷后处理 `partial_charges/INCAR`：
  - `LPARD = .TRUE.`：打开部分电荷后处理；
  - `LPARDH5 = .TRUE.`：把输出从一个或多个 PARCHG 文件改写到 vaspout.h5（便于 py4vasp 读取）；
  - `NBMOD = -3`：把 `EINT` 解释为**相对费米能级**的能量区间（而非绝对能量或带指标）；
  - `EINT = -0.25 0.0`：费米能以下 0.25 eV 的窗口，对应约 -250 mV 的样品偏压，只积分占据态；
  - `LSEPK/LSEPB = .FALSE.`：对所有 k 点与能带**求和**输出，而不是按 k 点/能带分开写文件。py4vasp 的 STM 模拟要求二者均为 `.FALSE.`。
- 后处理运行：`cp WAVECAR partial_charges/` 后运行 vasp_std（该步是单核后处理、未并行化，无需 mpi）。终端输出会报告落入能量窗口的 k 点与能带——本例只有**第 18、19 两条能带在 k 点 33 处**贡献。

**流程与结果**：用 py4vasp 模拟并绘图：`calc.partial_density.to_stm(supercell=7, tip_height=3)`（7×7 超胞重复、针尖高度 3 Å），并自定义色标。将模拟图与实验 STM 图对比——后者由 LMU/CeNS 在常温常压空气中测得（低温真空图像通常更清晰），但仍能看到原子分辨与清晰的三角格子：

![石墨 0001 表面的实验 STM 图像（常温常压）](https://vasp.at/tutorials/latest/_images/Graphite_ambient_STM.jpg)

模拟图像给出三角格子，与实验图（相差任意旋转角）吻合良好：直坐标 `0.0 0.0 0.5` 的表面原子呈亮斑，而 `0.333 0.666 0.5` 的另一个表面原子信号很弱。可自由修改针尖高度与偏压（改后处理步的 `EINT`）观察图像变化：**若负偏压幅度增大，更多能带参与成像**——次表面层有碳原子垫底的那个表面原子将变得可见，三角格子会转变为六角格子。

**复习问题**：(1) py4vasp 的 STM 实现能否模拟针尖性质？(2) 为什么 `NBMOD = -3` 对 STM 模拟很方便？它与实验偏压如何对应？(3) 为什么第二个表面原子只在偏压幅度增大后才出现在模拟图像中？

### 示例 10：石墨烯的恒定电流 STM 模拟

**教学目标**：(1) 计算石墨烯的部分电荷密度；(2) 模拟并绘制恒定电流模式的 STM 图像。

**原理**：恒流模式是 STM 最常用的数据采集模式——横向扫描时实时调节针尖高度使隧穿电流保持恒定，每个横向点的高度构成图像。py4vasp 中用 `to_stm(selection="constant_current")` 开启（默认为 `"constant_height"`），并可设定隧穿电流（nA，默认 1 nA）。算法：对每个横向网格点，从真空中部向表面方向扫描，直到平滑化部分电荷密度的值达到与所选电流成正比的目标值，记录该高度成图。

**输入要点**（`e10_STM_of_graphene`）：POSCAR 为单层石墨烯（2 个 C，面内晶格常数与示例 9 相同）；静态 INCAR、KPOINTS（17×17×1 Γ 中心）与示例 9 相同；部分电荷 INCAR 也相同（`LPARD/LPARDH5/NBMOD = -3/EINT = -0.25 0.0`；`LSEPK/LSEPB` 默认即 `.FALSE.`，不写也可以，但写明可避免歧义——再次强调 py4vasp 要求二者均非 `.TRUE.`）。

**流程与结果**：

1. 静态计算后可用 py4vasp 对比石墨烯与石墨顶层的 DOS：石墨烯两个原子等价、DOS 介于石墨第 9、10 号原子之间。注意石墨烯著名的零带隙电子结构在本计算中未正确再现——k 网格不够密，出现约 1 eV 的假隙；但这对 STM 图像模拟影响不大。
2. 后处理：复制 WAVECAR 到 `partial_charges/` 再运行，计算费米能 ±0.1 eV 附近的部分电荷。
3. 恒流模式成像：`calc.partial_density.to_stm(selection="constant_current", current=8, supercell=[5,4])`。偏压仍取 -250 meV 以便与示例 9 对比——对自由悬浮石墨烯，小偏压是合适的，否则针尖会对它产生显著吸引（J. Vac. Sci. Technol. B 31, 04D103, 2013）。py4vasp 支持选择自旋通道（"up"/"down"/"total"），磁性表面可能得到不同图像；本例非自旋极化，总是显示总电荷密度。
4. 结果：图像呈现石墨烯的六角对称性；两个原子化学环境相同，两个顶位上方的针尖高度等价（与石墨的三角格子形成鲜明对比——后者源于 AB 堆垛导致两个亚晶格位环境不同）。

**复习问题**：(1) py4vasp STM 模拟的默认模式是什么？(2) 若把电流从默认 1 nA 调大，针尖会更靠近还是更远离表面？(3) 能否用 py4vasp 分辨自旋相关的表面特征（类比自旋极化 STM）？

**本部分涉及的 VASP 功能与标签汇总**：部分电荷后处理全家桶（`LPARD`、`LPARDH5`、`NBMOD = -3`、`EINT`、`LSEPK`/`LSEPB`）、WAVECAR 续算工作流、`PREC = Single` 与 FFT 网格插值伪影、Bardeen 隧穿理论与 Tersoff-Hamann 近似、Chen 导数规则、py4vasp `partial_density.to_stm`（恒定高度/恒定电流两种模式、针尖高度、电流设定、超胞重复、自旋通道选择）、偏压极性与占据/未占态成像、石墨三角格子 vs 石墨烯六角格子的物理起源。

---

## 7. 过渡态（Transition States）

化学反应在反应物与产物两个极小值之间进行，势能面上连接二者的能量极大值即过渡态。学习用静态与动态方法建模反应过渡态。

## 7.1 过渡态（一）：基础过渡态方法（NEB / IDM / 肘形图 / IRC）

> 来源：https://vasp.at/tutorials/latest/transition_states/part1/ ｜ 练习文件：`transition_states-part1.zip` ｜ 示例目录：`$TUTORIALS/transition_states/e01_NEB`、`e02_IDM`、`e03_Elbow_plot`、`e04_IRC`
> 注：本教程使用的预训练机器学习力场（MLFF）可从教程页面给出的 Google Drive 链接下载，放入 `transition_states/mlff` 目录。

本部分是过渡态搜索方法的"四大件"入门：NEB 求最小能量路径（MEP）、改进二聚体方法（IDM）优化鞍点、双原子分子解离的肘形图（elbow plot）、以及内禀反应坐标（IRC）追踪反应全路径。四个示例分别以 Si 自扩散、H₂/Cu(111) 解离吸附、SN2 取代反应为载体，并穿插了 MLFF 加速振动分析与 IRC 的实战用法。

### 示例 1：NEB 求反应路径（Si 空位自扩散）

**教学目标**：(1) 准备并执行 NEB 计算；(2) 优化 image 结构与反应路径；(3) 可视化 image 结构与能量路径。

**原理**：所有化学反应都经由过渡态从反应物到达产物。**nudging elastic band（NEB）方法**（J. Chem. Phys. 113, 9978, 2000）在反应物与产物结构之间插值出若干**image**（通常 4–20 个）来建模最小能量路径（MEP）。image 之间以及 image 与端点之间用弹簧相连，构成一条连续的"弹性带"。带的受力分为两部分：垂直于带的"真实力"（$F_{i,\perp}$，来自势能面）与平行于带的弹簧力（$F_{i,\mathrm{spring}}$）——这种分解即所谓 "nudging"，它防止 image 沿路径方向滑向能量最低点。按合力 $F_i$ 弛豫整条带即得反应路径。image 越多路径越精细，但计算更贵、收敛更难。

![NEB 方法示意图](/imgs/2026-08-05/fig22.png)

**体系**：15 原子 Si 原胞超胞中的空位自扩散——金刚石立方 Si 中相邻位点间的扩散（Phys. Rev. B 59, 3969, 1999）。蓝色 Si 原子在 NEB 过程中移向红色空位位置：

![Si 空位超胞：蓝色 Si 原子向红色空位扩散](/imgs/2026-08-05/fig23.png)

**输入要点**（`e01_NEB`）：

- 目录结构：`00`（初态 V1 空位）→ `01`–`04`（四个 image）→ `05`（末态 V2 空位），NEB 必须从父目录运行，VASP 在各子目录中并行计算各 image（**MPI 进程数必须能被 image 数整除**）。教程提供的 POSCAR 已预弛豫，故只需几步 NEB 即收敛。
- 结构构建与插值（教程附脚本）：ASE `make_supercell` 建 16 原子超胞，分别删除倒数第 1/第 2 个原子得到两个相邻空位构型 POSCAR.V1/V2，先各自弛豫（`IBRION = 2` 共轭梯度、`EDIFFG = -0.01`、`ISIF = 0`），再用 ASE `NEB(...).interpolate(method="idpp")` 做 IDPP 插值（避免插值 image 中原子过于靠近；`apply_constraint` 会读取选择性动力学）。
- `INCAR.neb`：电子最小化 `ENCUT = 250`、`PREC = Normal`、`EDIFF = 1e-6`；`ISMEAR = 0`、`SIGMA = 0.05`；离子弛豫 `NSW = 100`、`EDIFFG = -0.04`（负值 = 力阈值 eV/Å，教程推荐力判据而非能量差判据）、`IBRION = 1`（准牛顿）、`POTIM = 0.8`、`ISIF = 2`；NEB 专有标签 `IMAGES = 4`、`SPRING = -5`（弹簧常数）。若输出出现 "BRIONS problems: POTIM should be increased"，说明收敛困难，可增大 POTIM 或换 IBRION 算法。KPOINTS：Γ 中心 2×2×2。
- 教程设置以速度优先；要高质量结果需加密 k 网格、增大超胞直到形成能等收敛，并把 EDIFFG 收紧到 -0.01 eV/Å。

**流程与结果**：

1. 父目录下运行 vasp_std 完成 NEB 弛豫。
2. 结构可视化：把六个目录的 CONTCAR 投影到指定平面叠画（ASE plot_atoms），可见绝大多数 Si 原子位置不变，扩散原子从初态平滑移动到空位，相邻 image 间距大致均匀：

![各 NEB image 投影结构](/imgs/2026-08-05/fig24.png)

3. 能量路径：用 py4vasp 逐 image 读 `energy.read()["free energy    TOTEN"]`，以扩散原子相对初态的位移为横坐标画相对能量曲线——能量剖面大致对称，峰值约 0.7 eV（即扩散势垒）。
4. 方法学：若想更精细地刻画该过渡过程，可在四个 image 之间再插入更多 image；但单次计算中 image 越少收敛越好，过多 image 会使优化变得非常昂贵。也可以取中间 image 作为下一步过渡态搜索的试探结构。

**复习问题**：(1) NEB 计算需要哪些标签？(2) NEB 找到的是什么路径？(3) 为什么不宜一次优化太多 image？

### 示例 2：改进二聚体方法（IDM）优化 H₂/Cu(111) 解离过渡态

**教学目标**：(1) 用振动计算求虚频模式；(2) 用 IDM 优化 H₂ 解离的过渡态；(3) 用虚频检验过渡态。

**原理**：NEB 给出反应路径，但**过渡态本身**（路径上的能量极大值、鞍点）需要专门方法优化。过渡态的判据是**只有一个虚频**。**改进二聚体方法（improved dimer method, IDM）**（J. Chem. Phys. 111, 7010, 1999）是基于力（而非 Hessian）的局域鞍点搜索算法，四步循环如下：

![IDM 四步弛豫示意图](/imgs/2026-08-05/fig25.png)

(a) 从试探结构（第一个点 $q$）的**最负振动模式**取初始方向 $u_\xi$；(b) 沿该方向在势能面上取第二个点 $q + \delta u_\xi$，两点构成"二聚体"；(c) 二聚体绕 $q$ 在势能面上旋转角度 $\phi_1$；(d) 把 $u_\xi$ 旋转 $\phi_{min}$ 以最小化负曲率 $\lambda$，用最小化算法定义搜索方向 $\bar{N}$，沿 $\bar{N}$ 平移 $\epsilon$ 步到 $q + \epsilon\bar{N}$。新结构重新定义二聚体，重复"旋转→平移"直至鞍点——即只剩一个负振动模式的过渡态。

**体系**：H₂ 在 Cu(111) 上的解离吸附（3 层 2×2 超胞，12 Cu + 2 H），该反应有大量实验与理论数据可对照（J. Chem. Theory Comput. 13, 3208, 2017）。

**输入要点**（`e02_IDM`，三个子目录 freq1/idm/freq2）：

- `POSCAR.idm`：底部 8 个 Cu 固定（F F F），顶层 4 个 Cu 与 2 个 H 放开（T T T）；坐标之后附加的 "! Trial unstable direction" 块就是试探不稳定方向（每原子一个三维矢量，从最负虚频的本征矢量提取）。注意 POSCAR 坐标段之后的附加行通常被解释为速度，这里是借该位置存放不稳定方向。
- `freq1/INCAR`（与 freq2 相同）：`IBRION = 5`（有限差分频率计算）、`NFREE = 2`、`NSW = 1`；`IDIPOL = 3` 消除 z 向周期性镜像的偶极影响；`ISMEAR = 1`、`SIGMA = 0.2`；**`ML_LMLFF = T` + `ML_MODE = run`——用预训练 MLFF 替代 DFT 力，大幅加速振动分析**。
- `idm/INCAR`：**`IBRION = 44`——启用 IDM 作为优化引擎**；`NSW = 20`（为演示速度只跑 20 步，真收敛需数百上千步）、`EDIFFG = -0.01`、`NELMIN = 5`（每个离子步的最少电子自洽步数）。
- KPOINTS：Γ 中心 3×3×1（表面体系，垂直方向单 k 点 + 大真空避免镜像相互作用）。

**流程与结果**：

1. **试探结构振动分析**（freq1）：链接 MLFF（`ln -s ../../mlff/MLFF_e02/ML_FFN ML_FF`）后运行。OUTCAR 中 `grep "cm-1"` 列出 18 个模式——16 个实频 + **2 个虚频**（`f/i=` 标记）：17 号 $i83.6$ cm$^{-1}$、18 号 $i964.4$ cm$^{-1}$。过渡态只应有一个虚频，IDM 将沿较大的模式弛豫、消除较小的那个。
2. 模式可视化：`vibFreq.py` 脚本生成 molden 文件，或 py4vasp `calculation.phonon.mode.print()` + `calculation.force_constant.to_molden()`。最负模式是 **H-H 键伸缩振动**（解离时正是要断裂的键），H₂ 位于桥位、箭头略偏离垂直方向：

![H2 吸附在 Cu(111) 桥位，橙色箭头为 H-H 伸缩振动模式](/imgs/2026-08-05/fig26.png)

3. **IDM 弛豫**（idm）：从 OUTCAR 提取 18 号模式本征矢量追加到 POSCAR（`grep -A 15 "18 f/i=" OUTCAR | awk ...`），运行 20 步 IDM。
4. **过渡态验证**（freq2）：再做一次振动分析，确认弛豫后结构只剩一个虚频——这是过渡态的必要判据。

**复习问题**：教程要求掌握 `IBRION = 5`/`NFREE`（频率）与 `IBRION = 44`（IDM）的设置逻辑，以及为什么表面计算垂直方向只取 Γ 点。

### 示例 3：H₂/Cu(111) 肘形图（2D 势能面切片）

**教学目标**：(1) 为肘形图准备结构网格；(2) 可视化势能面的 2D 切片。

**原理**：IDM 给出过渡态这一"点"，NEB 给出路径这条"线"；对双原子分子还可以直接画势能面（PES）的 2D 切片。双原子分子有 6 个自由度，固定表面后取分子质心高度 $z$ 与键长 $r$ 两个自由度做网格，逐点计算能量画等值线图，即**肘形图（elbow plot）**（Phys. Rev. Lett. 73, 1400, 1994; Phys. Rev. B 62, 8295）。本例沿袭示例 2 的 H₂/Cu(111) 解离吸附体系，H₂ 垂直置于桥位正上方：

![肘形图定义：r 与 z 网格上逐点计算能量](/imgs/2026-08-05/fig27.png)

**输入要点**（`e03_Elbow_plot`）：

- `POSCAR.surf`：3 层 Cu(111) 2×2 超胞（底部两层固定、顶层放开）。
- `INCAR`：`NSW = 0` 单点计算；`ISMEAR = 1`、`SIGMA = 0.2`；**`LCHARG = .FALSE.`、`LWAVE = .FALSE.`——不写 CHGCAR/WAVECAR 以节省磁盘**（本例要跑几十次单点）。KPOINTS：Γ 单点。
- 网格生成脚本：$z$ 取 1.0–3.0 Å、$r$ 取 0.5–2.5 Å，各 7 个点（含端点），共 **7×7 = 49 个结构**；每个结构建独立目录并复制 INCAR/POTCAR/KPOINTS；批量运行脚本用 **vasp_gam**——只用 Γ 点时波函数为实数，gamma 版 VASP 更高效。

**流程与结果**：

1. 用 py4vasp 逐目录读 TOTEN，在 z-r 网格上画填充等值线图（0.5 eV 间隔标黑色等高线）：

![肘形图（49 点，$\Gamma$ 点）](/imgs/2026-08-05/fig28.png)

2. 与更精细设置的参考图对比（121 点网格、3×3×1 k 网格）：两图外观相似，但精细版等高线更光滑（更密的点更好描述 PES），且能量明显更负（更密的 k 网格）；完全收敛还需更细网格：

![肘形图（121 点，3x3x1 k 网格）](/imgs/2026-08-05/fig29.png)

3. 物理解读：图中有**两个极小值**——一个对应桥位上的 H₂ 分子（分子态吸附），另一个对应解离后分别吸附在 fcc 与 hcp 位的两个 H 原子（解离吸附态）。肘形图直观展示了从分子态越过势垒到解离态的全过程。

**复习问题**：(1) 两张图为何差异明显？(2) 肘形图适用于什么体系？(3) 两个极小值对应什么吸附结构？

### 示例 4：内禀反应坐标（IRC）追踪 SN2 取代反应

**教学目标**：(1) 从过渡态出发沿 IRC 正向追踪到产物、反向追踪到反应物；(2) 可视化反应剖面与反应物/产物结构。

**原理**：**IRC**（J. Phys. Chem. A 106, 165, 2002）定义为过渡态在质量加权笛卡尔坐标下的**最陡下降路径**。从已充分优化（例如用 IDM）的过渡态出发，沿其虚频模式 $u_\xi$ 用**阻尼速度 Verlet（DVV）算法**正向推进到产物 P、反向推进到反应物 R，得到完整反应剖面。

![SN2 反应示意](/imgs/2026-08-05/fig30.png)

**体系**：Cl⁻ 对氯甲烷的 SN2 取代（12×12×12 Å³ 盒子，C + 3H + 2Cl，体系总电荷 -1）。过渡态中 CH₃ 呈平面构型、两个 Cl 与 C 等距；一个 Cl 进攻 C 的同时另一个 Cl 被弹出。为把计算时间压缩到几分钟，**全程使用 MLFF**。

**输入要点**（`e04_IRC`）：

- `POSCAR`：过渡态结构，坐标段后附 "! unstable direction optimized by the dimer method"（IDM 优化的不稳定方向，每原子一个矢量）。
- `INCAR` 的 IRC 标签：**`IBRION = 40`——启用 IRC 计算**；`IRC_DIRECTION = -1/+1` 初始运动方向（反向/正向）；`IRC_STOP = 20`（连续 20 步能量上升则终止）；`IRC_VNORM0 = 0.0005`（速度范数，越小越精确但所需离子步越多）；`NSW = 5000`；`ISIF = 2`；`NELECT = 22`（电荷 -1 → 电子数加 1；MLFF 层面并非必需）；`ML_LMLFF = .TRUE.` + `ML_MODE = run`。
- KPOINTS：Γ 单点（盒子中的分子）。ML_FF 用符号链接到 MLFF_e04 目录以节省磁盘。

**流程与结果**：

1. 正、反两个方向各建一个目录（`p` 与 `m`），唯一区别是 `IRC_DIRECTION` 的符号（用 `sed` 替换）。分别运行 vasp_gam——MLFF 使约 450 步的 IRC 在几秒内完成（纯 ab initio 需数小时）。
2. 分析：`ircShift.sh` 脚本从两个 OUTCAR 的 `IRC (A):` 行提取 IRC 坐标与能量（相对过渡态能量）并合并排序为 `ircShift.dat`；py4vasp 绘图得到一条**漂亮的对称曲线**——Cl⁻ 从 Cl-C 键的正对面进攻碳原子，使甲基翻转（Walden inversion），把原来的 Cl 以 Cl⁻ 形式弹出。由于反应物与产物化学上等价（输入=输出），曲线关于过渡态对称。
3. 结构可视化：py4vasp `calc.structure.plot([1,1,1])` 查看反应物/产物；`calc.structure[:].plot()` 播放 IRC 轨迹；教程还提供脚本把 XDATCAR 每 10 步抽帧转成 xyz 动画（JMol 播放），完整展示从反应物经过渡态到产物的平滑演变。

**复习问题**：(1) 改变哪个标签可以得到正向/反向轨迹？(2) 等价的 ab initio IRC 计算大约要多久？(3) 为什么本反应的 IRC 曲线是对称的？

**本部分涉及的 VASP 功能与标签汇总**：NEB（`IMAGES`、`SPRING`、多目录并行运行、MPI 进程整除约束、IDPP 插值）、IDM（`IBRION = 44`、试探不稳定方向写入 POSCAR、单虚频判据）、有限差分频率计算（`IBRION = 5`、`NFREE`、虚频 `f/i=` 输出、molden 导出）、肘形图方法学（`NSW = 0` 单点扫描、`LCHARG`/`LWAVE` 省盘、vasp_gam 加速）、IRC（`IBRION = 40`、`IRC_DIRECTION`、`IRC_STOP`、`IRC_VNORM0`、DVV 算法、`NELECT` 带电体系）、MLFF 加速（`ML_LMLFF`/`ML_MODE = run`）、`EDIFFG` 力判据、py4vasp 能量提取与轨迹可视化。

---

## 7.2 过渡态（二）：静态方法实战（沸石催化羰基化反应全流程）

> 来源：https://vasp.at/tutorials/latest/transition_states/part2/ ｜ 练习文件：`transition_states-part2.zip` ｜ 示例目录：`$TUTORIALS/transition_states/e05_Static_Vib_analysis` ～ `e09_Static_Thermodynamics`
> 注：本页使用的 MLFF 需从教程页面给出的 Google Drive 链接下载，放入 `transition_states/mlff` 目录。

Part 1 分别介绍了 NEB、IDM、IRC；真实体系（如沸石催化反应）要把这些方法**组合**使用。本部分以菱沸石（chabazite）催化的合成气制乙醇中的**羰基化速率控制步**为对象，完整走一遍静态过渡态工作流：试探结构振动分析 → IDM 优化过渡态 → 鞍点验证与自由能 → IRC 反应路径 → 反应热力学与速率常数。五个示例全程用预训练 MLFF 加速（示例 5 的振动分析若用 ab initio 需数小时，MLFF 只需几秒）。

**化学背景**：合成气（CO + H₂）在菱沸石催化下转化为乙醇 CH₃CH₂OH 经历六步基元反应：(1) CO + 2H₂ → CH₃OH；(2) CH₃OH + H-Z → H₂O + CH₃-Z；(3) CO + CH₃-Z → CH₃CO⁺ + Z⁻；(4) CH₃CO⁺ + Z⁻ → CH₃CO-Z；(5) CH₃OH + CH₃CO-Z → CH₃COOCH₃ + H-Z；(6) CH₃COOCH₃ + 2H₂ → CH₃CH₂OH + CH₃OH。**速率控制步（RDS）是第 3 步——甲基化菱沸石的羰基化**（PCCP 13, 2603, 2011），产物为自由的乙酰阳离子 CH₃CO⁺；视热力学条件它或接回骨架成 CH₃CO-Z，或经质子化变成 CH₂CO⁺（J. Catalysis 396, 166, 2021）。本教程只计算正向 RDS 的活化自由能 $\Delta A_{R\to P}^{\ddagger}$（$\ddagger$ 标记过渡态），故羰基化之后的产物状态不相关。

![菱沸石结构：红 O、蓝 Si、青 Al；沸石空腔提供巨大内表面](/imgs/2026-08-05/fig31.png)

### 示例 5：试探过渡态的振动分析

**教学目标**：(1) 对试探过渡态做振动分析；(2) 找出用于过渡态弛豫的试探不稳定方向。

**试探结构的构造逻辑**：根据反应物与预期产物的几何猜测过渡态——平面 CH₃ 基团放在沸石 O 原子与 CO 的 C 原子中间，两个距离都任意取 2.1 Å。若猜测接近鞍点，结构应只有一个虚频；若有多个虚频，取**最强虚频**作为 IDM 的试探方向。

**输入要点**（`e05_Static_Vib_analysis`）：POSCAR 为 42 原子（Al C₃ H₃ O₂₅ Si₁₁）菱沸石超胞内的试探 TS；`INCAR`：`IBRION = 5`（有限差分）、`NFREE = 2`、`POTIM = 0.01`、`NSW = 1`（对 IBRION=5，NSW 只需为正，VASP 会把它重置为单胞原子数的三倍）；`ML_LMLFF = T` + `ML_MODE = run`——用 MLFF 时无需任何电子结构标签。KPOINTS：Γ 单点（大超胞分子类体系）。

**流程与结果**：链接 MLFF（`ln -s mlff/MLFF_e05-12/ML_FFN ML_FF`）后运行 vasp_gam，几秒钟完成。`vibFreq.py` 脚本（或 py4vasp `calculation.phonon.mode.print()`）提取 126 个频率模式，前八个中有 **4 个虚频**：$i504$、$i184$、$i159$、$i34$ cm$^{-1}$，随后是三个近乎消失的模式（数值噪声/平动，频率 $\sim 10^{-5}$ cm$^{-1}$）。结论：试探结构还不是过渡态，取最强虚频（$i504$ cm$^{-1}$）的本征矢量作为 IDM 的试探方向。

### 示例 6：IDM 优化过渡态

**教学目标**：用 IDM 沿最负虚频弛豫到过渡态，并理解收敛行为（梯度 vs 曲率）。

**输入要点**（`e06_Static_TS_Opt`）：`POSCAR.init` = 示例 5 的 POSCAR + 末尾追加最强虚频对应的位移矢量（从 `vibFreq.py` 生成的 `hesseModes.dat` 取前 43 行，或用 py4vasp `calculation.force_constant.eigenvectors()[0]` 提取）；`INCAR`：`IBRION = 44`（IDM）、`NSW = 5000`、`EDIFFG = -0.005` eV/Å（比 Part 1 更严格的力判据）、MLFF 开启；Γ 单点。

**流程与结果**：

1. **梯度收敛**：`grep grad OUTCAR | awk '{print $4}'` 提取最大梯度画曲线——前 100 步快速下降，200 步后力基本不变。
2. **曲率收敛**：`grep curvature OUTCAR` 提取 PES 曲率 $\lambda$（每 4 个离子步更新一次）。与梯度相反，曲率初期上升，其振荡与梯度振荡同步但收敛慢得多——**曲率才是限制因素**，这解释了为什么力 200 步就收敛、计算却跑到约 600–800 步。过渡态要求曲率最小化（唯一负曲率方向）。
3. **结构变化**：py4vasp 对比初末结构，肉眼差别很小；用 ASE `geometry.get_distances`（含 pbc）定量测量——C-C 与 C-O 距离都从 2.10 Å 收缩到约 2.0 Å：过渡态把反应物拉得更近。
4. **效率对比**：若用 Hessian 方法，每步都要做一次频率计算，代价巨大；IDM 每步只需几次单点计算，高效得多。

**复习问题**：(1) 曲率每隔几个离子步更新一次？(2) IDM 相比逐步重算 Hessian 的优势是什么？

### 示例 7：鞍点验证与过渡态自由能

**教学目标**：(1) 验证弛豫后的结构是过渡态（一阶鞍点）；(2) 计算 IRC 所需的不稳定试探方向；(3) 计算过渡态的自由能。

**原理**：一阶鞍点应当**只有一个虚频**（扣除消失的平动模式）。该虚频模式即 IRC 的试探方向。此外由振动频率可计算关键热力学量：零点振动能

$$E_{ZPV} = \sum_i \frac{\hbar}{2} \nu_i$$

以及温度 $T$ 下的振动焓 $H_{vib}$ 与振动自由能 $A_{vib}$。注意这里的自由能是 **Helmholtz** 自由能（不是 Gibbs）——因为采用正则（NVT）系综，$N$、$V$、$T$ 固定。$A_{vib}$ 加到总能 $E_{tot}$ 上得到过渡态的总 Helmholtz 自由能 $A_{tot}$，后续算速率常数要用。

**输入要点**（`e07_Static_Saddle_point`）：POSCAR = 示例 6 的 CONTCAR；INCAR 与示例 5 相同（`IBRION = 5` 振动分析 + MLFF）。

**流程与结果**：

1. 频率对比（IDM 前 vs 后）：IDM 前有 3 个 >100 cm$^{-1}$ 的虚频加 2 个小的；IDM 后只剩 **1 个显著虚频 $i522$ cm$^{-1}$**，另有几个 <10 cm$^{-1}$ 的近零模式——这些是**消失的平动模式**，计算振动热力学量前必须手动剔除。鞍点验证通过。
2. 用 molden（`hesseModes_shifted.molden`）或 py4vasp 可视化虚频模式的运动方向。
3. 热力学量：用 ASE `HarmonicThermo`（滤掉虚频与消失模式，126 个自由度中保留 122 个）计算 $T = 440$ K 下的各项：ZPE = 4.387 eV，$U = 6.387$ eV，$S = 0.00890$ eV/K（$TS = 3.916$ eV），振动自由能 $A_{vib} = 2.471$ eV。
4. 总 Helmholtz 自由能：$E_{tot} = -319.263$ eV（OUTCAR TOTEN），$A_{tot}^{TS} = -316.792$ eV。孤立的自由能不能直接与实验对比——它的意义在于与反应物态作差得到可测量的活化自由能。

**复习问题**：(1) 剩余的近零模式是什么类型？(2) 零点振动能如何计算？

### 示例 8：IRC 反应路径（连接反应物与产物）

**教学目标**：(1) 沿不稳定试探方向追踪 IRC；(2) 找出并可视化反应物态与产物态。

**原理**：过渡态位于反应物与产物之间，只有一个不稳定振动模式；沿该模式正向传播到产物、反向到反应物（阻尼速度 Verlet 算法，J. Chem. Phys. A 106, 165, 2002），同时验证过渡态确实连接预期的反应物与产物。本例验证羰基化 TS 连接"乙酰阳离子 + Z⁻"与"CO + CH₃-Z"。

**输入要点**（`e08_Static_Reaction_pathway`）：POSCAR.init = 示例 6/7 的 TS 结构 + 虚频本征矢量；INCAR：`IBRION = 40`（IRC）、`IRC_VNORM0 = 0.0005`（速度矢量重标定的常数值，提供阻尼；越小越精确、步数越多）、`IRC_DIRECTION = -1`（反向→反应物；+1 为正向→产物）、MLFF 开启；Γ 单点。

**流程与结果**：

1. 在 `m`（反向）与 `p`（正向，`sed` 改 IRC_DIRECTION）两个目录分别运行 IRC。
2. `ircShift.sh` 合并两个 OUTCAR 的 `IRC (A):` 行得到完整反应剖面；py4vasp 绘图：**约 -3 Å 处是反应物**（初始 CO 与 CH₃-Z）；CO 靠近、CH₃ 翻转，能量升至 **0 Å 处的过渡态极大值**；随后自由乙酰阳离子形成，能量快速降至**约 +4 Å 处的产物**。
3. py4vasp 分别可视化反应物（m/）、过渡态（e06 的 CONTCAR）、产物（p/）结构；`pos2xyzAnim_rev.py`/`pos2xyzAnim_for.py` 把 XDATCAR 每 10 步抽帧拼成动画（JMol 播放）。
4. IRC 在连续 `IRC_STOP`（默认 20）步能量单调上升时终止，故**倒数第 20 个结构最接近势能极小**；`pickMin_irc.sh 20` 从 XDATCAR 提取它写入 `POSCAR.pick`。由于 IRC_VNORM0 很小（0.0005），该结构已非常接近极小点，无需再弛豫。

**复习问题**：`IRC_DIRECTION = 1` 与 `-1` 的端点分别给出什么结构？

### 示例 9：反应热力学（静态）——活化自由能与速率常数

**教学目标**：(1) 对 IRC 得到的反应物结构做振动分析；(2) 对比谐振子量子处理与半经典处理；(3) 计算反应活化自由能；(4) 计算速率常数。

**原理**：ab initio 给出总能，实验测量的是自由能等热力学量，两者桥接需要振动频率提供的热力学修正（零点能、热贡献、熵）。反应物与过渡态的自由能差即活化自由能 $\Delta A_{R\to P}^{\ddagger}$；代入 **Eyring 方程**得速率常数：

$$k_{R\to P} = \frac{k_B T}{h} e^{-\Delta A_{R\to P}^{\ddagger} / k_B T}$$

**输入要点**（`e09_Static_Thermodynamics`）：POSCAR = 示例 8 的 `m/POSCAR.pick`（IRC 反向端点）；INCAR 与示例 5 相同（振动分析 + MLFF）。

**流程与结果**：

1. 反应物验证：对端点结构做频率计算——**无虚频**，确认是能量极小点，即真正的反应物态。
2. 量子（谐振子）处理：与示例 7 同样的 ASE HarmonicThermo 流程，分别得到反应物与 TS 在 440 K 的振动自由能；$A_{tot}^{TS} - A_{tot}^{R}$ 给出 $\Delta A^{\ddagger}_{QM} = 1.051$ eV。
3. 半经典处理：`vibpartitionClassic_eV.py` 脚本（不含零点能）给出 TS 振动自由能 1.505 eV、反应物 1.474 eV，$\Delta A^{\ddagger}_{SC} = 1.069$ eV。
4. **两种近似的对比**：半经典与量子的 Helmholtz 自由能各自相差约 1 eV（零点能贡献），但**活化自由能之差仅 0.018 eV**——量子效应在能垒上几乎抵消。这说明在本条件下半经典近似有效——正是对分子动力学（MD 力计算本质上半经典、不含零点能）的重要验证。
5. **速率常数**（Eyring，440 K）：半经典 $k = 5.259$ s$^{-1}$，静态量子 $k = 8.486$ s$^{-1}$——振动能量的小差异（0.018 eV）经指数放大造成速率约 40% 的差别。微小的能量差经指数化后会在可测量上产生大变化，这就是热力学量必须仔细收敛的原因。

**复习问题**：(1) 振动量子效应对 $\Delta A^{\ddagger}$ 有多重要？(2) 这对速率常数的计算造成多大差别？

**本部分涉及的 VASP 功能与标签汇总**：沸石催化反应建模（菱沸石超胞、RDS 识别）、`IBRION = 5`/`NFREE`/`POTIM` 有限差分频率、虚频谱解读（试探结构多虚频 → 取最强者）、`IBRION = 44` IDM 过渡态优化（梯度 vs 曲率收敛行为）、一阶鞍点判据（单虚频、剔除平动模式）、`IBRION = 40` IRC（`IRC_VNORM0`、`IRC_DIRECTION`、`IRC_STOP`、`pickMin` 提取极小点结构）、振动热力学（ZPE、$H_{vib}$、$A_{vib}$、Helmholtz vs Gibbs、ASE HarmonicThermo、半经典对比）、Eyring 方程速率常数、MLFF 全程加速（`ML_LMLFF`/`ML_MODE = run`、vasp_gam）。

---

## 7.3 过渡态（三）：动态方法（反应坐标、Blue Moon 系综、动态热力学）

> 来源：https://vasp.at/tutorials/latest/transition_states/part3/ ｜ 练习文件：`transition_states-part3.zip` ｜ 示例目录：`$TUTORIALS/transition_states/e10_Dynamic_Coordinates`、`e11_Dynamic_Validation`、`e12_Dynamic_Thermodynamics`
> 注：本页 MLFF 与 Part 2 相同，需从教程页面链接下载放入 `transition_states/mlff`。

静态方法用单个驻点结构代表反应物/过渡态/产物，但真实体系要用**系综**描述。本部分继续 Part 2 的菱沸石羰基化反应，用分子动力学（MD）走完整的动态过渡态工作流：定义并测试反应坐标 → Blue Moon 系综验证反应坐标 → 动态热力学（热力学积分）计算活化自由能与速率常数。三个示例均在 440 K、CSVR 恒温器下进行，全程 MLFF 加速。

**为什么需要动态方法**：PES 上布满局域极小，对大多数体系没有哪个单一驻点能代表完整的反应物、过渡态或产物——需要一组结构的系综。系综效应可以很大，例如酸性菱沸石上烯烃异构化势垒的系综修正可达约 20 kJ/mol（J. Catal. 373, 361, 2019）；温度升高时熵甚至可能改变反应机理本身（J. Chem. Phys. 131, 214508, 2009）——这些都是静态方法无法捕捉的。

### 示例 10：定义反应坐标（MD 测试与概率密度）

**教学目标**：(1) 用 MD 模拟测试试探反应坐标；(2) 确定找到试探反应物结构的概率密度。

**原理**：要计算正向活化自由能 $\Delta A_{R\to P}^{\ddagger}$，必须先定义**反应坐标** $\xi$——它应是原子位置 $\mathbf{q}$ 的连续函数，且无穷缓慢地改变 $\xi$ 应是可逆的（Annu. Rev. Phys. Chem. 67, 669, 2016）。本反应取

$$\xi = r_2 - r_1$$

其中 $r_1$ 是 CO 的 C 与甲基 C 的距离，$r_2$ 是甲基 C 与其所连沸石 O 的距离。目标过渡态处的反应坐标记为 $\xi^*$。MD 产生的系综中，概率密度 $P(\xi)$ 的最可几结构被选作候选反应物参考态 $\xi_{ref,R}$，作为后续约束 MD 的起点。

![反应坐标定义示意](/imgs/2026-08-05/fig32.png)

**输入要点**（`e10_Dynamic_Coordinates`）：

- `POSCAR` = Part 2 示例 9 的 CONTCAR（静态反应物）。
- `INCAR.init`：MD 设置 `IBRION = 0`、`NSW = 10000`、`POTIM = 1.0` fs、`MDALGO = 5`（**CSVR 恒温器**——正则系综速度重标度采样）、`TEBEG = 440`、`CSVR_PERIOD = 20`（恒温器时间尺度，以 MD 步数计）；`ISYM = 0`、`NPAR = 2`；MLFF：`ML_LMLFF = T`、`ML_MODE = run`、`ML_OUTPUT_MODE = 0`（减少文件输出提升 MD 性能）、`ML_OUTBLOCK = 10`（每 10 步输出一次）；`RANDOM_SEED` 由运行脚本追加，定义初始速度并保证可重复。
- **`ICONST` 文件**——约束/监测坐标的定义：`R 3 25 0`（原子 3 与 25 的距离 r₂）、`R 3 2 0`（原子 3 与 2 的距离 r₁）、`S -1 1 7`（两个 R 的组合：ξ = r₂ − r₁）。末位标志 `0` = 约束、`7` = **仅监测**——本例只监测不约束。

**流程与结果**：

1. 运行脚本把 50 ps（5 段 × 10000 步 × 1 fs）的 MD 拆成 5 次连续运行：每段结束从 REPORT 提取随机数种子写入下一段的 INCAR（保证随机序列连续），CONTCAR 滚动为下一段 POSCAR。
2. REPORT 文件解读：开头给出 MD 初始条件（MDALGO、CSVR_PERIOD、RANDOM_SEED、自由度数 123）；每 10 步输出 `mc>`（被监测的反应坐标值）、`e_b>`（E_tot/E_pot/E_kin 等）、`tmprt>`（T_sim 目标温度与 T_inst 瞬时温度）。
3. 从各 report 提取 `mc>` 列绘制 ξ 随时间演化，做直方图得到概率密度 $P(\xi)$：**峰位 ξ ≈ 2.4 Å 即最佳采样区，也就是动态意义上的反应物态**。选最佳采样区的 ξ 作参考态可压低 $P(\xi)$ 的相对误差。
4. **静态 vs 动态反应坐标**：静态反应物（POSCAR.1）的 ξ = 1.95 Å，动态最可几值 ξ ≈ 2.4 Å，相差 0.45 Å——静态用单一结构、动态用分布描述反应物态，差异显著。
5. 从保存的 POSCAR.2–6 中找出 ξ 最接近峰位的结构（参考数据为 POSCAR.2，ξ ≈ 2.32 Å、$P \approx 0.87$ Å$^{-1}$），作为下一示例约束 MD 的起点。

**复习问题**：(1) 概率密度峰对应什么结构？(2) 静态与动态反应坐标有何差异？

### 示例 11：反应坐标的验证（约束 MD + Blue Moon 系综）

**教学目标**：(1) 用约束 MD 把反应物（R）转变到过渡态（TS）；(2) 检验自由能梯度的连续性；(3)（可选）检验 TS→R 的可逆性。

**原理**：反应坐标 $\xi$ 是否合适，可通过追踪从 R 到 TS 的自由能梯度 $(\partial A/\partial \xi)_{\xi^*}$ 来检验——它决定活化自由能。用**约束 MD** 追踪会产生有偏的统计平均，**Blue Moon 系综**正是为此设计：实现简单、便于统计误差估计、可诊断 $\xi$ 选择的问题。采用**慢生长（slow-growth）方法**：以转变速度 $\dot{\xi}$（每步 ξ 的改变量）沿 $\xi$ 从反应物扫描到产物，TS 位于自由能图的极大值处。$\dot{\xi}$ 必须足够小，使自由能梯度从 R 到 TS 连续——连续性成立即证明 $\xi$ 选择合理；最终检验是可逆性（从 TS 反向重复过程回 R）。

**输入要点**（`e11_Dynamic_Validation`）：POSCAR、ICONST 与示例 10 相同（注意此示例中 ICONST 的组合坐标处于约束状态）；INCAR 在示例 10 基础上新增两个标签：`INCREM = -1.1e-4`——慢生长协议的转变速度增量，**负值使 ξ 从反应物向过渡态方向推进**（$\dot{\xi}$ = 总距离/总步数）；`LBLUEOUT = T`——执行 Blue Moon 约束 MD，把自由能梯度写入 REPORT。`NSW = 4000`（6 段 × 4000 步 ≈ 24 ps）。

**流程与结果**：

1. 从示例 10 选定的 POSCAR 出发运行 6 段慢生长 MD。若出现 `Error: SHAKE algorithm did not converge!`，说明所选起点离反应物态太远，需重选更接近 $P(\xi)$ 峰位的结构（约束 MD 用 SHAKE 算法维持约束）。
2. 得到 7 个沿 ξ 均匀分布的结构（POSCAR.1–6 + 最终 CONTCAR）；`pos2xyzStatic.py` 拼成动画检查机理：第 1 帧 CH₃ 仍连在沸石 O（ZO）上、CO 在盒子别处；随后 CO 与 CH₃ 靠近、ZO-CH₃ 键断裂、CH₃ 翻转，最终形成乙酰阳离子 H₃C-CO⁺——与预期机理一致。注意这 7 帧各自是所在系综的代表结构。
3. REPORT 中的关键输出：`cc>`（约束坐标的当前值与目标值）、`b_m>`（Blue Moon 量：Lagrange 乘子 $\lambda$、$|z|^{-1/2}$、$Gk_BT$ 与自由能梯度 $|z|^{-1/2}(\lambda + Gk_BT)$）。
4. `fgrad.sh` 脚本从各 report 提取 ξ 与自由能梯度到 `grad.dat`，绘图检查梯度从 R 到 TS 的连续性——连续则 $\xi$ 选择得到确认。

**复习问题**：教程要求理解 INCREM 符号与扫描方向的关系、SHAKE 收敛失败的排查。

### 示例 12：动态热力学（活化自由能与速率常数）

**教学目标**：(1) 从 R 到 TS 沿反应坐标求自由能梯度；(2) 导出反应活化自由能；(3) 在动态近似下计算速率常数；(4) 对比静态与动态结果。

**原理**：R→TS 的转变以可逆功的形式给出活化自由能。具体流程是对一系列不同 ξ 的约束 MD（Blue Moon 系综）做热力学积分：

$$\Delta A_{R\to TS}^{\ddagger} = \Delta A_{\xi_{ref,R}\to\xi^*} - k_B T \ln\left[\frac{h}{k_B T}\frac{\langle|\dot{\xi}^*|\rangle}{2}P(\xi_{ref,R})\right]$$

其中可逆功由自由能梯度积分给出：

$$\Delta A_{\xi_{ref,R}\to\xi^*} = \int_{\xi_{ref,R}}^{\xi^*} \left(\frac{\partial A}{\partial \xi}\right)_{\xi^*} d\xi$$

自由能梯度本身是约束系综平均：

$$\left(\frac{\partial A}{\partial \xi}\right)_{\xi^*} = \frac{1}{\langle|Z|^{-1/2}\rangle_{\xi^*}} \langle |Z|^{-1/2}[\lambda + G k_B T]\rangle_{\xi^*}$$

$Z$ 是质量度规张量，$\lambda$ 是 SHAKE 约束的 Lagrange 乘子，$G$ 由 $|Z|$、$Z^{-1}$、$\nabla\xi$、$\nabla|Z|$ 导出（详见 VASP Wiki "Blue moon ensemble"）。最后用 Eyring 方程把 $\Delta A^{\ddagger}$ 换算成速率常数 $k$。

**输入要点**（`e12_Dynamic_Thermodynamics`）：7 个初始 POSCAR 取自示例 11 生成的结构；INCAR 与示例 11 类似，但**去掉 INCREM（不再慢生长，改为逐点约束 MD）**、保留 `LBLUEOUT = T`，并把 `NSW` 从 4000 提高到 **10000**（每个 ξ 点跑 10 ps 以保证统计收敛，7 个点共约 70–140 ps）。ICONST 中组合坐标 S 改为约束模式。

**流程与结果**：

1. 对 7 个 ξ 点分别运行约束 MD，从 REPORT 的 `b_m>` 提取各点自由能梯度；对梯度做样条积分得到自由能曲线 `dA_splines.dat`——极小点对应反应物、极大点对应 TS。
2. 补齐活化自由能公式的其余项：(a) $\langle|\dot{\xi}^*|\rangle$——用 `aux_vel.sh 440` 从 TS 点轨迹计算反应坐标速度（得 $1.05\times10^{13}$ Å/s）；(b) $P(\xi_{ref,R})$——示例 10 的概率密度直方图在 ξ_ref,R 处的取值。
3. **三种方法的最终对比**（440 K）：半经典 $\Delta A^{\ddagger} = 1.067$ eV；静态量子 $1.048$ eV；**动态量子 $1.127$ eV**。静态与半经典只差 0.018 eV（见 §7.2），而**动态比静态高约 0.08 eV**——这个更显著的差别来源于系综平均（而非单一结构）的引入。
4. **速率常数**（Eyring）：半经典 5.563 s$^{-1}$、静态量子 8.984 s$^{-1}$、**动态量子 1.142 s$^{-1}$**——自由能上约 0.1 eV 的小差别使速率常数变化近一个数量级。

**复习问题**：(1) 静态与动态方法的活化自由能差多大？(2) 反应自由能的小差别对速率常数影响多大？

**本部分涉及的 VASP 功能与标签汇总**：反应坐标定义与 `ICONST` 文件（R/S 组合坐标、flag 0 约束 / 7 监测）、`IBRION = 0` MD、`MDALGO = 5` CSVR 恒温器与 `CSVR_PERIOD`、`RANDOM_SEED` 续接多段 MD、`ML_OUTPUT_MODE = 0`/`ML_OUTBLOCK` MD 输出控制、REPORT 文件解读（`mc>`/`e_b>`/`tmprt>`/`cc>`/`b_m>`）、概率密度 $P(\xi)$ 分析、约束 MD 与 SHAKE、Blue Moon 系综（`LBLUEOUT = T`）、慢生长方法（`INCREM`）、自由能梯度与热力学积分、活化自由能公式（可逆功 + 速度/概率修正项）、Eyring 速率常数、静态 vs 半经典 vs 动态系综方法的定量对比。

---

## 8. 杂化泛函（Hybrid Functionals）

杂化泛函是交换关联泛函的特殊类别，把一定比例的 Fock 交换混入密度泛函，灵感来自对无关联体系有精确解的 Hartree–Fock 方法。

## 8.1 杂化泛函（一）：可用泛函总览（PBE0 / B3LYP / HF / 屏蔽杂化 / HSE06 能带）

> 来源：<https://vasp.at/tutorials/latest/hybrids/part1/> ｜ 练习文件：[hybrids-part1.zip](https://vasp.at/tutorials/latest/hybrids-part1.zip) ｜ 示例目录：`$TUTORIALS/hybrids/e01_Si-gap`、`e02_Ar-gap`、`e03_MgO-gap-optimization`、`e04_CaS-band`

计算效率最高的泛函是 LDA、GGA 与 meta-GGA（半局域类），但它们无法正确描述带隙。杂化泛函（hybrid functional）把一定比例的 Fock 交换混入半局域泛函，灵感来自对无关联体系有精确解的 Hartree–Fock 方法，目标是矫正 LDA/GGA 著名的带隙低估。非局域 Fock 势为

$$V_{\mathrm{x}}^{\mathrm{HF}}(\mathbf{r},\mathbf{r}') = -\frac{e^2}{2}\sum_{n\mathbf{k}} f_{n\mathbf{k}} \frac{\psi_{n\mathbf{k}}^{*}(\mathbf{r}')\,\psi_{n\mathbf{k}}(\mathbf{r})}{|\mathbf{r}-\mathbf{r}'|}$$

其中 $f_{n\mathbf{k}}$ 是能带占据数。由于它依赖 KS 轨道而非密度，杂化泛函严格说**不是**密度泛函。本页通过四个例子总览 VASP 中可用的杂化泛函：例 1 用 PBE0 改善 Si 带隙并学习杂化泛函的收敛设置；例 2 用 HF 与 B3LYP 计算 Ar 带隙并认识无屏蔽杂化泛函的根本缺陷；例 3 讲屏蔽杂化泛函（HSE）的原理并实战优化 MgO 带隙；例 4 给出杂化泛函能带结构的完整做法（CaS，PBE vs HSE06）。

### 示例 1：Si 带隙（PBE 与 PBE0）

**教学目标**：(1) 解释什么是杂化泛函；(2) 设置一个 PBE0 计算；(3) 计算带隙。

**任务与原理**：计算立方金刚石 Si 在 PBE 与 PBE0 下的带隙。PBE0 是 VASP 支持的杂化泛函之一，其构成为

$$E_{\mathrm{xc}}^{\mathrm{PBE0}} = \frac{1}{4} E_{\mathrm{x}}^{\mathrm{HF}} + \frac{3}{4} E_{\mathrm{x}}^{\mathrm{PBE}} + E_{\mathrm{c}}^{\mathrm{PBE}}$$

即 25% Fock 交换 + 75% PBE 交换 + 全部 PBE 关联。与许多其他泛函不同，PBE0 的系数**不含**任何对实验数据拟合的参数。KS 轨道反映的物理是否正确强烈依赖于交换关联泛函的类型，对带隙计算而言 meta-GGA 与杂化泛函比 LDA/GGA 更合适。

**输入要点**（目录 `e01_Si-gap`，示例文件内容略）：
- POSCAR：金刚石结构 Si 的 fcc 原胞（$a = 5.46873$ Å，2 原子）；KPOINTS 为 7×7×7 Γ 心网格。
- PBE 计算：`GGA = PE` 显式指定泛函（教程提问：查 `GGA` 的默认值，显式设置是否必要？——VASP 会按 POTCAR 自动选择，显式写出是好习惯但非必需）；`ISMEAR = -5`（四面体方法 + Blöchl 修正）给出最平滑的 DOS；`LORBIT = 11` 把投影 DOS 写入 DOSCAR 与 vaspout.h5。
- PBE0 计算（子目录 `PBE0`）：`LHFCALC = T` 开启 Fock 交换即激活 PBE0——只要 `AEXX`、`AGGAC`、`ALDAC` 保持默认（`AEXX = 0.25`、`AGGAX = 1-AEXX = 0.75`、`AGGAC = ALDAC = 1`），`GGA = PE` 就自动对应上式的系数。这些前置因子标签正是各交换/关联分量的权重。
- **算法与展宽的讲究**：杂化泛函电子弛豫用阻尼直接优化算法 `ALGO = Damped` + `TIME = 0.8`——它总能变分、比迭代对角化稳健，杂化泛函常需要这种稳健性。由于 Blöchl 修正使能量非变分，弛豫阶段不能用 `ISMEAR = -5`，改用高斯展宽 `ISMEAR = 0`、`SIGMA = 0.01`；收敛后再做单发计算拿 DOS。
- 教程声明本例带隙精度在 0.1 eV 以内；若要证实并提高精度，需对截断能 `ENCUT` 与 k 点数做收敛研究。

**流程与结果**：
1. 在 `e01_Si-gap` 运行 PBE（`vasp_std`），用 py4vasp 画 DOS：带隙对应费米能附近没有态的能量区间（能量轴已对齐费米能）。
2. 把 `WAVECAR` 复制到 `PBE0` 目录再运行（`mpirun -np 2 vasp_std`）——复用 PBE 波函数可大幅加速杂化泛函收敛；运行时 stdout 每步会多一行阻尼算法的信息。
3. **重要概念**：不同交换关联泛函的总能**不可直接比较**，只能比较带隙这类物理量。
4. 用 py4vasp 从占据数判断 HOMO/LUMO 求带隙：PBE 为 0.616 eV，PBE0 为 1.842 eV。普遍趋势是杂化泛函给出更大带隙，部分矫正 GGA 的低估。
5. 单发 DOS：把 PBE0 的 INCAR 改为 `ISMEAR = -5`、`ALGO = None`、`LORBIT = 11`（删去不再使用的 `SIGMA`、`TIME`），从已有 WAVECAR 重启。对比发现 PBE0 的 DOS 不只是带隙变大——能带间的相对位置与带宽都不同，这些细微变化会实质性影响材料性质。

**复习问题**：
1. 什么是杂化泛函？
2. 如何计算一个体系的带隙？
3. `LHFCALC` 标签控制什么？
4. 为什么 **ISMEAR = -5** 不能与 **ALGO = Damped** 搭配使用？

### 示例 2：Ar 带隙（PBE / HF / B3LYP）

**教学目标**：(1) 使用 B3LYP 泛函；(2) 使用 HF 方法；(3) 解释 HF 与 B3LYP 的局限性。

**任务与原理**：计算 fcc Ar 晶体的带隙，实验值高达 14.4 eV。B3LYP 最初为描述分子性质而设计，是**半经验**杂化泛函——其参数通过拟合测试集中原子/分子的解离能、电子-质子亲和能与电离能等实验数据得到。经验参数正是乘在 B3LYP 各分量前的系数：

$$E_{\mathrm{x}}^{\mathrm{B3LYP}} = 0.8\,E_{\mathrm{x}}^{\mathrm{LDA}} + 0.2\,E_{\mathrm{x}}^{\mathrm{HF}} + 0.72\,\Delta E_{\mathrm{x}}^{\mathrm{B88}}$$

$$E_{\mathrm{c}}^{\mathrm{B3LYP}} = 0.19\,E_{\mathrm{c}}^{\mathrm{VWN3}} + 0.81\,E_{\mathrm{c}}^{\mathrm{LYP}}$$

HF 计算则是 100% Fock 交换（`AEXX = 1`）且不含任何关联（$E_{\mathrm{c}} = 0$）。B3LYP 等**无屏蔽**泛函继承了 HF 方法的根本缺陷——在均匀电子气极限下失效。根源是 Fock 势傅里叶变换在 $\mathbf{q} = 0$ 处的奇点：

$$V_{\mathrm{x}}^{\mathrm{HF}}(\mathbf{q}) = -\frac{4\pi e^2}{2}\sum_{n\mathbf{k}} f_{n\mathbf{k}} \frac{\psi_{n\mathbf{k}}^{*}(\mathbf{r}')\,\psi_{n\mathbf{k}}(\mathbf{r})}{\mathbf{q}}$$

这会导致巡游态被完全压制。因此 HF 方法与无屏蔽杂化泛函**不推荐**用于金属、小带隙半导体等巡游性显著的系统（Paier et al., 2007）。Ar 是范德华结合的宽隙分子晶体，正好适合检验这些泛函。

**输入要点**（目录 `e02_Ar-gap`，示例文件内容略）：
- POSCAR：fcc Ar 原胞（晶格常数约 5.32 Å，1 原子）；KPOINTS 为 4×4×4；`ENCUT = 300`。
- PBE：`ISMEAR = -5`、`LORBIT = 11`。
- HF（子目录 `HF`）：`LHFCALC = T` 开启 Fock 交换、`AEXX = 1` 设定 HF 比例为满额，于是只用 Fock 交换、无半局域交换关联。教程追问其余标签的默认值：`GGA = PE`，且因 `AEXX = 1`，`AGGAX = 1 - AEXX = 0`，`AGGAC` 与 `ALDAC` 被自动置零——即 HF 只需两个标签即可实现。
- B3LYP（子目录 `B3LYP`）：用 `GGA = B3` 切换到 B3LYP 对应的半局域分量，再显式设置各分量权重 `AEXX = 0.2`、`AGGAX = 0.72`、`AGGAC = 0.81`、`ALDAC = 0.19`，与上面两式一一对应。
- HF 与 B3LYP 的电子弛豫同样用 `ALGO = Damped` + `TIME = 0.4`，弛豫阶段 `ISMEAR = 0`、`SIGMA = 0.01`。

**流程与结果**：先跑 PBE 并用 py4vasp 查看 DOS；把 WAVECAR 分别复制到 `HF` 与 `B3LYP` 目录运行，最后比较带隙：

| 方法 | 带隙 (eV) |
| --- | --- |
| PBE | 8.504 |
| HF | 17.612 |
| B3LYP | 10.313 |
| 实验 | 14.4 |

结论：PBE 严重低估（差约 6 eV），HF 与 B3LYP 也都偏离实验值 3 eV 以上——**没有一个方法接近实验值**。对范德华结合的分子晶体，这些泛函既缺色散作用又各有系统性偏差，说明"杂化"本身不是万能药，选泛函必须结合体系的成键特征。

**复习问题**：
1. PBE0、B3LYP 与 HF 对哪一类体系描述失败？为什么？（巡游/金属体系，源于 Fock 势 $\mathbf{q}=0$ 奇点对金属态的压制）
2. `AGGAC` 标签控制什么？

### 示例 3：MgO 屏蔽杂化泛函带隙优化

**教学目标**：(1) 理解屏蔽杂化泛函的基本原理；(2) 解释 HSE06 泛函的基本假设；(3) 针对带隙优化 Fock 交换比例与屏蔽长度。

**任务与原理**：优化屏蔽杂化泛函中 Fock 交换的比例，以复现 MgO 的实验带隙 7.7 eV。屏蔽杂化泛函把 Fock 交换的 $1/r$ 依赖分解为短程（SR）与长程（LR）两部分，引入屏蔽参数 $\mu$：

$$E_{\mathrm{xc}}^{\mathrm{screened}} = \frac{1}{4} E_{\mathrm{x}}^{\mathrm{HF,SR}}(\mu) + \frac{3}{4} E_{\mathrm{x}}^{\mathrm{PBE,SR}}(\mu) + E_{\mathrm{x}}^{\mathrm{PBE,LR}}(\mu) + E_{\mathrm{c}}^{\mathrm{PBE}}$$

短程极限即 PBE0；截去长程的 Fock 贡献后：(1) 避免了 Fock 交换的奇点——正是它压制所有金属态，因此屏蔽杂化泛函可用于金属；(2) 可通过对 Hartree–Fock 算符**降采样**（downsampling）进一步降低计算成本。固体物理常用的是 HSE06：屏蔽参数经经验优化取 $\mu = 0.2$ Å$^{-1}$，对应作用范围 $2/\mu \approx 10$ Å，通常覆盖几个近邻。VASP 中用 `HFSCREEN` 标签控制屏蔽参数（相关标签还有 `HFRCUT` 与 `HFALPHA`）。

带隙大小仍强烈依赖于 Fock 交换的比例。本例采取实用主义路线：调节混合比例去复现实验带隙。教程特别提醒——这样做**牺牲了第一性原理方法的预测能力**，更接近"模型物理"（定性理解体系）；且文献中偶尔把这种拟合泛函称为 HSE 是不对的：**HSE06 仅对应 `AEXX = 0.25` 且 `HFSCREEN = 0.2`**。

**输入要点**（目录 `e03_MgO-gap-optimization`，示例文件内容略）：
- POSCAR：岩盐结构 MgO 的 fcc 原胞（$a \approx 4.21$ Å，Mg、O 各 1 个）；KPOINTS 为 4×4×4 Γ 心网格；`ENCUT = 500`。
- PBE：`ISMEAR = -5`、`LORBIT = 11`。
- 屏蔽杂化（子目录 `screened_hybrid`）：`LHFCALC = T`、`HFSCREEN = 0.2`、`AEXX = 0.25`（即标准 HSE06），弛豫阶段 `ISMEAR = 0`、`SIGMA = 0.01`、`ALGO = Damped`、`TIME = 0.4`。

**流程与结果**：
1. 先跑 PBE 看 DOS 中的带隙。
2. **AEXX 循环**：在 `screened_hybrid` 目录运行 `loop_aexx.sh`——对 `AEXX = 0.40/0.45/0.50/0.55` 逐个重写 INCAR、从 PBE 的 WAVECAR 重启、算完把 `vaspout.h5` 存入 `aexx$值` 子目录。用 py4vasp 提取各带隙：

| AEXX | 带隙 (eV) |
| --- | --- |
| 0.40 | 7.219 |
| 0.45 | 7.580 |
| 0.50 | 7.946 |
| 0.55 | 8.313 |

带隙随 Fock 交换比例单调增大；要复现 7.7 eV，最优 `AEXX` 在 0.45–0.50 之间。

3. **HFSCREEN 循环**：固定接近最优的 `AEXX = 0.5`，用 `loop_hfscreen.sh` 扫描 `HFSCREEN = 0.15/0.20/0.25/0.30`：

| HFSCREEN (Å$^{-1}$) | 带隙 (eV) |
| --- | --- |
| 0.15 | 8.291 |
| 0.20 | 7.946 |
| 0.25 | 7.643 |
| 0.30 | 7.372 |

结果对屏蔽参数相当敏感：`HFSCREEN` 越小，Fock 交换的作用范围越长、保留的 Fock 交换越多，带隙越大。两个循环合起来展示了"材料专属杂化泛函"的完整拟合流程及其代价。

**复习问题**：
1. `HFSCREEN` 标签控制什么？
2. PBE0 是否依赖作用范围（range）？（否，PBE0 是无屏蔽的全程杂化泛函）
3. HSE06 泛函包含哪些贡献？（短程 Fock 交换 ×0.25 + 短程 PBE 交换 ×0.75 + 长程 PBE 交换 + PBE 关联）

### 示例 4：CaS 能带结构（PBE 与 HSE06）

**教学目标**：用杂化泛函计算能带结构——这是本页技术上最有含量的一例。

**任务与原理**：对 CaS 做 HSE06 计算，求能带结构与带隙。能带沿第一布里渊区中连接高对称点的**高对称线**计算——这些点因对称性简并而特别重要，能带极值（VBM/CBM，即 HOMO/LUMO）常出现在其上。若 VBM 与 CBM 位于同一 k 点则带隙为**直接**带隙，否则为**间接**带隙。生成 k 路径可用 SeeK-path 工具（接受 POSCAR，有网页界面）：

![SeeK-path 生成的 CaS 第一布里渊区](/imgs/2026-08-05/fig33.png)

*图：SeeK-path 生成的 CaS 第一布里渊区与高对称点路径。*

**输入要点**（目录 `e04_CaS-band`，示例文件内容略）：
- POSCAR：fcc CaS 原胞（$a \approx 5.66$ Å，Ca、S 各 1 个）；SCF 统一用 6×6×6 Γ 心网格；`ENCUT = 500`。
- PBE SCF：`ISMEAR = -5`、`LORBIT = 11`；`band/` 子目录做非自洽能带：`ICHARG = 11`（读取 CHGCAR 固定电荷密度）+ `ISMEAR = 0`，KPOINTS 用**线模式**（line-mode）沿 L–W–X–Γ–L 路径生成 k 点——这是产生路径 k 点的简便方式。
- HSE06 SCF：`LHFCALC = T`、`HFSCREEN = 0.2`（`AEXX` 默认 0.25），`ALGO = Damped`、`TIME = 0.4`、`ISMEAR = 0`、`SIGMA = 0.01`。
- **杂化泛函能带的关键技巧**：杂化泛函计算能带时既要高对称路径上的 k 点（画能带），又要均匀 k 网格（评估 Fock 交换能）。因此 `HSE06/band/KPOINTS` 用**显式 k 点列表**写入了 53 个 k 点 = 16 个带权重的不可约布里渊区 k 点（来自 IBZKPT，张成第一布里渊区用于 Fock 交换）+ 37 个**零权重** k 点（定义与 PBE 相同的路径）。零权重保证这些点只贡献本征值、不参与自洽电荷密度与 Fock 能的积分。
- `HSE06/band/INCAR` 额外设置 `NELMIN = 4`：从 WAVECAR 重启时体系总能已收敛，VASP 在 `NELMIN` 次迭代后就会停止更新轨道、只在该轨道基底下对角化 KS 哈密顿量；但新增 k 点上的轨道可能尚未充分收敛，基底不够好会导致本征值不准、能带出现跳变——`NELMIN = 4` 强制至少更新四次以避免这一问题。

**流程与结果**：
1. PBE：SCF → `cd band; cp ../CHGCAR .; vasp_std`，用 py4vasp 画能带（教程还引导查阅文档、画出能带的轨道特征 character）。
2. HSE06 SCF：复制 PBE 的 WAVECAR 后 `mpirun -np 2`；运行间隙的练习题是比较 `HSE06/band/KPOINTS`、`IBZKPT` 与 PBE OUTCAR 中的 k 点——结论：能带 k 点两者相同，杂化泛函额外需要均匀网格（权重按不可约布里渊区取），路径 k 点权重取零。
3. HSE06 能带：`cd HSE06/band; cp ../WAVECAR .; mpirun -np 4 vasp_std`（耗时较长）。作图时会看到**全部** k 点（含规则网格上的点），需缩放至两个 L 点标记之间查看路径能带。导带上若仍有轻微扭折（kink），可加大 `NELMIN` 或增加 `NBANDS` 改善。
4. 带隙对比：PBE 为 2.382 eV，HSE06 为 3.457 eV。带隙是**间接**的——LUMO 在 X 点、HOMO 在 Γ 点。与实验值 5.4 eV 相比，HSE06 显著改善但仍未完全达到（标准 HSE06 对 CaS 这类体系仍偏低，这正是例 3 参数优化思路的动机）。

**复习问题**：
1. 用杂化泛函计算能带结构需要哪些步骤？
2. 为什么本例中高对称路径上的 k 点权重为零？
3. 什么是间接带隙？

**本部分涉及的 VASP 功能与标签汇总**：
- 杂化泛函开关与配方：`LHFCALC`（开启 Fock 交换）、`AEXX`/`AGGAX`/`AGGAC`/`ALDAC`（四个交换/关联分量的前置系数）、`GGA`（`PE`→PBE0，`B3`→B3LYP 分量基）；HF = `LHFCALC = T` + `AEXX = 1`；PBE0 = 默认系数；HSE06 = `HFSCREEN = 0.2` + `AEXX = 0.25`。
- 屏蔽杂化：`HFSCREEN`（屏蔽参数 μ）、`HFRCUT`、`HFALPHA`；屏蔽带来的好处——避免 Fock 奇点、可用于金属、支持 HF 算符降采样（downsampling）。
- 电子算法：`ALGO = Damped` + `TIME`（杂化泛函稳健的阻尼直接优化）、`ALGO = None`（单发 DOS/能带）。
- 展宽与 DOS：弛豫阶段 `ISMEAR = 0` + 小 `SIGMA`（Blöchl 修正非变分），后处理单发 `ISMEAR = -5`；`LORBIT = 11` 输出投影 DOS。
- 重启与加速：WAVECAR 复用（PBE→杂化）、`ICHARG = 11` 非自洽能带、`NELMIN` 保证新增 k 点轨道更新、`NBANDS` 改善能带扭折。
- KPOINTS 技巧：线模式（半局域能带）vs 显式列表 = 带权 IBZKPT 网格 + 零权重路径点（杂化能带）。
- 分析：py4vasp 画 DOS/能带、按占据数求 HOMO–LUMO 带隙；不同泛函的总能不可比较。

## 9. 线性响应（Linear Response）

除基态性质外，VASP 还能计算体系对外界微扰的响应。微扰足够小时响应处于线性区，可应用线性响应理论。

## 9.1 响应函数（一）：静态与频率依赖的介电响应

> 来源：<https://vasp.at/tutorials/latest/response/part1/> ｜ 练习文件：[response-part1.zip](https://vasp.at/tutorials/latest/response-part1.zip) ｜ 示例目录：`$TUTORIALS/response/e01_SiC-static-finite-differences`、`e02_AlP-static-DFPT`、`e03_NaCl-optics`

除基态性质外，VASP 还能计算体系对外界微扰的响应。微扰足够小时响应处于线性区，可用线性响应理论处理。本页三个例子覆盖响应计算的两大路线：**有限差分法**（显式施加微扰）与 **DFPT 密度泛函微扰理论**（解析求解响应方程），对象分别是电场、应变与离子位移三类微扰。例 1 用有限差分计算 SiC 的全套静态响应张量；例 2 用 DFPT 计算 AlP 的静态介电常数并对比两种方法；例 3 计算 NaCl 频率依赖的介电函数（电子 + 离子贡献），并顺带给出声子频率与频率单位换算。

### 示例 1：有限差分法计算静态响应（SiC）

**教学目标**：(1) 说出并计算对小电场、应变、离子位移的静态响应；(2) 掌握离子位移与介电响应有限差分法相关的 INCAR 标签，使用 PEAD 方法；(3) 区分离子冻结（ion-clamped）与离子弛豫（relaxed-ion）版本的极化率、弹性模量与压电张量。

**任务与原理**：对 SiC 施加小的静态电场、应变与离子位移，用有限差分法计算响应。总能 $E$ 可在无微扰体系 $E^0$ 附近按电场 $\mathcal{E}$、应变 $\eta$、离子位移 $u$ 展开到二阶：

$$E(\mathbf{u}, \mathcal{E}, \eta) = E^{0} + \frac{\partial E}{\partial u_i} u_i + \frac{\partial E}{\partial \mathcal{E}_i} \mathcal{E}_i + \frac{\partial E}{\partial \eta_{ij}} \eta_{ij} + \frac{1}{2}\frac{\partial^2 E}{\partial u_i\partial u_j} u_i u_j + \frac{1}{2}\frac{\partial^2 E}{\partial \mathcal{E}_i\partial \mathcal{E}_j} \mathcal{E}_i \mathcal{E}_j + \frac{1}{2}\frac{\partial^2 E}{\partial \eta_{ij}\partial \eta_{km}} \eta_{ij} \eta_{km} + \frac{\partial^2 E}{\partial u_i\partial \mathcal{E}_j} u_i \mathcal{E}_j + \frac{\partial^2 E}{\partial u_i\partial \eta_{jk}} u_i \eta_{jk} + \frac{\partial^2 E}{\partial \mathcal{E}_i\partial \eta_{jk}} \mathcal{E}_i \eta_{jk}$$

（重复指标求和；注意对电场的线性响应也会与其他微扰组合出现在二阶项中）。一阶系数按物理意义命名（$\Omega_0$ 为原胞体积，力与应力由 Hellmann–Feynman 定理计算）：

$$F_i = -\Omega_0 \left.\frac{\partial E}{\partial u_i}\right|_{\mathcal{E},\eta}, \qquad P_i = -\left.\frac{\partial E}{\partial \mathcal{E}_i}\right|_{u,\eta}, \qquad \sigma_{ij} = \left.\frac{\partial E}{\partial \eta_{ij}}\right|_{u,\mathcal{E}}$$

其中极化 $P$ 按现代极化理论（Berry 相位）写成对全部占据态 Berry 相的求和。二阶系数即各类响应张量：

- **力常数（Hessian）矩阵** $K_{ij} = \Omega_0\, \partial^2 E/\partial u_i \partial u_j$；
- **离子冻结极化率** $\chi_{ij} = -\partial^2 E/\partial \mathcal{E}_i \partial \mathcal{E}_j = \partial P_i/\partial \mathcal{E}_j$；
- **零场离子冻结弹性模量** $C_{ijkm} = \partial^2 E/\partial \eta_{ij}\partial \eta_{km}$；
- **Born 有效电荷** $Z^*_{ij} = -(\Omega_0/e)\, \partial^2 E/\partial u_i \partial \mathcal{E}_j$（位移-电场交叉响应）；
- **内应变张量** $\Lambda_{ijk} = -\Omega_0\, \partial^2 E/\partial u_i \partial \eta_{jk}$（力-应变交叉响应）；
- **离子冻结压电张量** $e^{(0)}_{ijk} = -\partial^2 E/\partial \mathcal{E}_i \partial \eta_{jk}$（电场-应变交叉响应）。

有限差分可在两个独立层面实施：**离子位移与应变**用 `IBRION = 5/6`；**外加电场**用 `LCALCEPS`（其内部用 PEAD——离散化后的微扰表达式方法）。

**输入要点**（目录 `e01_SiC-static-finite-differences`，示例文件内容略）：
- POSCAR：fcc 闪锌矿 SiC 原胞（2 原子）；KPOINTS 为 4×4×4 Γ 心网格；另有 `KPOINTS_OPT` 提供能带高对称路径（Γ–X–U/K–Γ–L–W）。
- 基态 INCAR：默认 GGA-PBE，高斯展宽 `ISMEAR`/`SIGMA`；`ENCUT` 设得较高以**避免 Pulay 应力**（弹性模量计算的关键细节）；`EDIFF` 控制电子收敛。
- `electric/`：在 DFT 之上施加小静态电场，`LCALCPOL = T` 计算电子极化（Berry 相位）。`LCALCEPS = T` 启用 PEAD：用有限差分表达式求场极化 KS 轨道胞周期部分对 k 的导数，差分算符由重叠矩阵 $S_{mn}(\mathbf{k}_j, \mathbf{k}_{j+1})$ 构造。相关标签默认值：`LPEAD = T` 时 `EFIELD_PEAD = 0.01 0.01 0.01`（步长）、`IPEAD = 4`（步数）；设置 `EFIELD_PEAD` 必须考虑**带隙大小与 k 点数**。
- `electric-ionic/`：在电场基础上用 `IBRION` 开启离子位移有限差分，步长与步数由 `POTIM` 与 `NFREE` 设定；`ISIF = 3` 允许体积与晶胞形状变化——这是计算弹性模量的必要条件。

**流程与结果**：
1. 基态计算后，用 py4vasp 画 `KPOINTS_OPT` 路径能带、求带隙：SiC 为 1.285 eV（PBE）。
2. `electric/`：复制 WAVECAR 后运行。教程提问：该计算为何给出警告？——**电场过强可能发生 Zener 隧穿**；解决办法是减少 k 点和/或降低 `EFIELD_PEAD`（后续计算将其减小一个数量级）。
3. `electric-ionic/`：同样从基态 WAVECAR 出发运行，即可一次性获得前述全部静态响应量。用 py4vasp 依次打印：`force`（两原子受力大小相等方向相反，约 ±0.294 eV/Å）、`polarization`（离子偶极 4.35 |e|Å 三方向相等，电子偶极极小）、`stress`、`force_constant`（对角约 ±19.6 eV/Å²）、`elastic_modulus`、`born_effective_charge`（Si +2.757，C −2.757，各向同性）、`internal_strain`、`piezoelectric_tensor`。
4. **冻结 vs 弛豫**：实验可比的是**离子弛豫**版本，VASP 会同时输出两者，其关系为

$$\chi_{ij}|_{\eta} = \chi_{ij}|_{u,\eta} + \Omega_0^{-1} Z^*_{mi} (K^{-1})_{mn} Z^*_{nj}, \qquad C_{ijkm}|_{\mathcal{E}} = C_{ijkm}|_{u,\mathcal{E}} - \Omega_0^{-1} \Lambda_{nij} (K^{-1})_{no} \Lambda_{okm}$$

$$e^{(0)}_{ijk} = e^{(0)}_{ijk}|_{u} - \Omega_0^{-1} Z^*_{mi} (K^{-1})_{mn} \Lambda_{njk}$$

即弛豫版本 = 冻结版本 + 通过力常数逆矩阵传递的离子弛豫贡献。输出中可见：弹性模量剪切分量 clamped 2836 → relaxed 2575 kBar；压电张量非零分量 −1.060 → −0.173 C/m²——离子弛豫显著削弱 SiC 的压电响应。极化率这里不单独打印，因为 VASP 直接用它构造介电张量（下一示例展开）。同样信息也可在 OUTCAR 末尾找到。

**复习问题**：
1. Born 有效电荷是什么？VASP 用什么单位报告它？（$|e|$，含局域场效应）
2. **LCALCEPS = T** 时极化如何计算？（PEAD/Berry 相位）
3. `ISIF` 标签控制什么？
4. `LPEAD` 能用于金属体系吗？（不能——PEAD 需要带隙）

### 示例 2：DFPT 计算静态介电响应（AlP）

**教学目标**：(1) 把介电常数与极化率联系起来；(2) 用 DFPT 计算静态介电性质；(3) 解释长波近似以及宏观与微观响应的区别；(4) 说明 `LEPSILON` 与 `LCALCEPS` 的优缺点。

**任务与原理**：在长波近似、忽略局域场效应的条件下，用 DFPT 计算 AlP 的静态介电常数。介电常数把材料内电场 $\mathcal{E}_i(\mathbf{r},\omega)$ 与外加电场 $\mathcal{E}_{\mathrm{ext}}$ 联系起来；动量空间中为介电矩阵 $\epsilon_{\mathbf{G}\mathbf{G}'\,ij}(\mathbf{q},\omega)$。若材料均匀或电场在远大于原子间距的尺度上变化，可用**长波近似**：$\epsilon$ 只依赖 $\mathbf{r}-\mathbf{r}'$，且当**局域场效应**可忽略时 $\mathbf{G}\neq\mathbf{G}'$ 的非对角元消失。离子冻结电子介电常数与极化率的关系为

$$\epsilon^{\infty}_{ij} = \delta_{ij} + \frac{4\pi}{\epsilon_0} \chi_{ij}$$

上标 $\infty$ 指对远高于声子频率的交流（AC）/光频场的响应，即 $\mathbf{q}\to 0$ 极限（极化激元领域的惯用记号）。

DFPT 路线：电场对 KS 轨道胞周期部分的线性响应 $|\xi_{n\mathbf{k}}\rangle$ 通过线性 Sternheimer 方程求得

$$\left[ H(\mathbf{k}) - \epsilon_{n\mathbf{k}} S(\mathbf{k}) \right] |\xi_{n\mathbf{k}}\rangle = -\Delta H_{\mathrm{SCF}}(\mathbf{k}) |u_{n\mathbf{k}}\rangle - \hat{\mathbf{q}} |\nabla_{\mathbf{k}} u_{n\mathbf{k}}\rangle$$

其中 $\Delta H_{\mathrm{SCF}}$ 是局域场效应（KS 轨道变化引起的哈密顿量微观变化），`LRPA = F` 时取零；轨道对 k 的导数再由第二个 Sternheimer 方程求出。最后宏观静态介电矩阵为

$$\hat{\mathbf{q}} \cdot \epsilon^{\infty} \cdot \hat{\mathbf{q}} = 1 - \frac{f 8\pi e^2}{\Omega_0} \sum_{n\mathbf{k}} w_{\mathbf{k}} \langle \hat{\mathbf{q}} \nabla_{\mathbf{k}} u_{n\mathbf{k}} | \xi_{n\mathbf{k}} \rangle$$

**LEPSILON 与 LCALCEPS 对比**（核心考点）：两者做同一种二阶展开，但 `LEPSILON`（DFPT）**可用于金属**，而 `LCALCEPS`（有限差分/PEAD）不行；两者都**不需要对未占据态求和**。但 `LEPSILON` 目前不能用于显式依赖轨道的交换关联泛函（HF 类与杂化泛函）。

**输入要点**（目录 `e02_AlP-static-DFPT`，示例文件内容略）：
- POSCAR：闪锌矿 AlP 原胞（$a \approx 5.45$ Å）；KPOINTS 4×4×4；基态 INCAR：`ISMEAR = 0`、`SIGMA = 0.01`、`EDIFF = 1E-6`、`ENCUT = 400`。
- `electric-ionic/INCAR`：`LASPH = T`；电场微扰用有限差分 `LPEAD = T`；极化响应用 DFPT `LEPSILON = T`；离子位移用 DFPT `IBRION = 8`（微扰论版本的 Hessian/声子，对应有限差分的 `IBRION = 5/6`）。

**流程与结果**：基态 → 复制 WAVECAR 到 `electric-ionic` 运行。长波近似下宏观 AC 介电响应 $\epsilon^{\infty}_{ij}(\hat{\mathbf{q}},\omega) \approx \lim_{\mathbf{q}\to 0} \epsilon_{00\,ij}$；若局域场效应不可忽略，则取其逆：$1/\epsilon^{\infty}_{ij} \equiv \lim_{\mathbf{q}\to 0} \epsilon^{-1}_{00\,ij}$。py4vasp 打印结果：
- **介电张量**（无量纲，含 DFT 局域场效应）：离子冻结 $\epsilon^{\infty} = 7.553$（各向同性），离子弛豫 9.782——离子贡献使静态介电常数增大约 30%。
- **Born 有效电荷**：Al +2.230，P −2.230。
- **压电张量**：闪锌矿结构无对称中心故压电张量非零；clamped 分量如 0.555/−0.392 C/m²，relaxed 版本大幅减小（约 0.038/−0.027）——与 SiC 例的规律一致。

**复习问题**：
1. 什么是 AC 介电常数？它如何与极化率联系？
2. DFPT 与有限差分法计算介电常数各有什么优缺点？分别由哪些 INCAR 标签控制？
3. 什么时候局域场效应可以忽略？
4. 什么是静态响应？

### 示例 3：频率依赖的介电函数（NaCl）

**教学目标**：(1) 计算频率依赖的介电响应；(2) 用 Kramers–Kronig 关系联系介电函数的实部与虚部；(3) 把介电函数与光学性质联系起来；(4) 考虑离子对频率依赖介电响应的贡献；(5) 计算声子频率并换算频率单位。

**任务与原理**：在独立粒子近似下计算 NaCl 频率依赖的介电函数。完整介电函数包含电子与离子两部分，需在 DFT 基态之上分别计算再相加：

$$\varepsilon(\omega) = \varepsilon_{\mathrm{elec}}(\omega) + \varepsilon_{\mathrm{ion}}(\omega)$$

实部与虚部由 Kramers–Kronig 关系联系；从介电函数出发可进一步导出反射率、吸收与光导率等光学性质（详见 Wiki 的 `LOPTICS` 文档）。

**电子部分**：极化中轨道对 k 的导数可用二阶微扰论展开

$$\nabla_{\mathbf{k}} |u^{(\mathcal{E}_i)}_{n\mathbf{k}}\rangle = \sum_{n'\neq n} \frac{|u^{(\mathcal{E}_i)}_{n'\mathbf{k}}\rangle \left\langle u^{(\mathcal{E}_i)}_{n'\mathbf{k}} \left| \frac{\partial [H(\mathbf{k}) - \epsilon_{n\mathbf{k}} S(\mathbf{k})]}{\partial \mathbf{k}} \right| u^{(\mathcal{E}_i)}_{n\mathbf{k}} \right\rangle}{\epsilon_{n\mathbf{k}} - \epsilon_{n'\mathbf{k}}}$$

这正是 `LOPTICS = T` 所用公式——**引入了对未占据态的求和**，因此必须设置足够的 `NBANDS` 并对 NBANDS 做收敛测试。

**输入要点**（目录 `e03_NaCl-optics`，示例文件内容略）：
- POSCAR：岩盐 NaCl 的 fcc 原胞（$a \approx 5.56$ Å）；KPOINTS 4×4×4 Γ 心；基态 INCAR：`ALGO = Normal`、`ISMEAR = -5`、`EDIFF = 1E-6`、`LORBIT = T`。
- `electron/INCAR`：`ALGO = Exact`（光学计算需精确对角化）、`LASPH = T`、`LOPTICS = T`（频率依赖介电矩阵）、`CSHIFT = 0.1`（复数展宽，默认值）、`NBANDS = 32`（未占据带求和）、`NEDOS = 2000`（频率网格点数）。
- `ion-dfpt/INCAR`：纯 DFPT 方案——`IBRION = 8` + `LEPSILON = T`，`EDIFF = 1E-8`、`PREC = High`。
- `ion-finite-differences/INCAR`：纯有限差分方案——`IBRION = 6`、`NFREE = 2`、`POTIM = 0.015`、`LPEAD = T`、`LCALCEPS = T`。
- **赝势选择**：POTCAR 用的是 **GW 推荐的 Na、Cl 势**而非普通 DFT 势。原因：标准 DFT 中未占据带不贡献基态性质，而 GW 势对费米能级以上的性质也经过测试——凡是需要增大 `NBANDS` 的计算（光学、GW）都应优先选 GW 推荐势。

**流程与结果**：
1. 基态：py4vasp 给出 NaCl 带隙 5.292 eV（PBE），画 DOS。
2. `electron/`：复制 WAVECAR 运行。教程提供 `get_dielectric_function.sh` 脚本（awk 解析 vasprun.xml 中 `dielectricfunction` 段，输出实部/虚部到 `dielectric_function.dat`，取三个对角分量的平均），也可直接用 `mycalc.dielectric_function.plot()`。虚部在带隙以上出现吸收边，实部由 Kramers–Kronig 变换得到。
3. **离子部分**：DFPT（`ion-dfpt`）与有限差分（`ion-finite-differences`）两条路线都要跑并对比——**有限差分方案耗时得多**（每个原子位移一次独立计算），DFPT 一次完成。两者给出的声子频率与离子介电响应应一致。
4. **声子与介电函数的联系**：离子位移给出振动方程 $m_{\mathrm{ion}} \omega^2 \mathbf{u} = (\partial^2 E/\partial u_i \partial u_j)\,\mathbf{u}$，Hessian 即例 1 的力常数矩阵；OUTCAR 中 "Eigenvectors and eigenvalues of the dynamical matrix" 段给出声子频率与振型。离子介电函数的峰位即声子频率（Lyddane–Sachs–Teller 型响应），其最大计算频率自动取为 $1.2 \times \omega_{\max}$。
5. **质量加权细节**：OUTCAR 中两种原子的位移本征矢量平行但长度不同——这是质量加权的结果，两矢量长度之比恰为 $\sqrt{m_{\mathrm{Na}}/m_{\mathrm{Cl}}}$（质量取自 POTCAR 的 `POMASS`）。
6. **频率单位换算**（本例重要知识点）：由振动方程可得计算内禀单位 $[\omega] = \sqrt{(\mathrm{eV}/\text{Å}^2)/\mathrm{a.m.u.}} = 9.822517\times 10^{13}\ \mathrm{s}^{-1}$，换算得

$$1\ \mathrm{THz} = 4.1357\ \mathrm{meV} = 33.356\ \mathrm{cm}^{-1}, \qquad 1\ \mathrm{meV} = 0.242\ \mathrm{THz} = 8.066\ \mathrm{cm}^{-1}, \qquad 1\ \mathrm{cm}^{-1} = 0.030\ \mathrm{THz} = 0.124\ \mathrm{meV}$$

注意电子介电函数的频率轴单位是 eV 而离子部分是 2π·THz，两者不同，作图比较前须先换算（如 `w*4.135667/(2*np.pi)` 把 2π·THz 转为 meV）。最终能量尺度对比：离子贡献集中在几十 meV（声子能区），电子贡献在 eV 能区（带间跃迁），两者频率尺度相差约两个数量级。

**复习问题**：
1. 声子如何进入频率依赖的介电函数？
2. **LOPTICS = T** 时计算的是什么？
3. 如何设置电子介电函数的频率点数？（`NEDOS`）
4. 频率单位如何从 cm$^{-1}$ 换算到 eV？

**本部分涉及的 VASP 功能与标签汇总**：
- 静态响应（有限差分路线）：`LCALCPOL`（Berry 相位极化）、`LCALCEPS` + `LPEAD`（PEAD 电场响应）、`EFIELD_PEAD`/`IPEAD`（电场差分步长/步数，需按带隙与 k 点设置）、`IBRION = 5/6` + `POTIM`/`NFREE`（离子位移/应变差分）、`ISIF = 3`（弹性模量需变胞形）。
- 静态响应（DFPT 路线）：`LEPSILON = T`（介电张量 + Born 电荷 + 压电张量，可用于金属、无需未占据态，但不支持杂化泛函）、`IBRION = 7/8`（微扰论二阶导数）、`LRPA`（局域场开关）、`LASPH`。
- 频率依赖响应：`LOPTICS = T`（独立粒子光学介电函数，需 `NBANDS` 收敛 + GW 推荐势）、`CSHIFT`（复数展宽）、`NEDOS`（频率网格）、`ALGO = Exact`；离子贡献由 `IBRION = 8`+`LEPSILON`（DFPT）或 `IBRION = 6`+`LCALCEPS`（有限差分，更慢）给出。
- 输出与分析：py4vasp 的 `force`/`polarization`/`stress`/`force_constant`/`elastic_modulus`/`born_effective_charge`/`internal_strain`/`piezoelectric_tensor`/`dielectric_tensor`/`dielectric_function`；vasprun.xml 解析脚本；OUTCAR 声子本征段；clamped vs relaxed 张量关系；频率单位换算链。

## 10. GW 近似（GW Approximation）

GW 近似是多体微扰理论中的一种近似，自能取 Green 函数 G 与屏蔽库仑相互作用 W 之积。GW 通过准粒子能量给出体系的谱性质，是目前计算带隙最精确的多体方法之一。

## 10.1 GW 近似（一）：入门（$G_0W_0$、$GW_0$ 与 Wannier90 能带）

> 来源：<https://vasp.at/tutorials/latest/gw/part1/> ｜ 练习文件：[GW-part1.zip](https://vasp.at/tutorials/latest/GW-part1.zip) ｜ 示例目录：`$TUTORIALS/GW/e01_Si-G0W0`、`e02_Si-GW0-band`

GW 近似是多体微扰理论中处理电子关联的高级方法：自能取 Green 函数 $G$ 与屏蔽库仑相互作用 $W$ 之积，通过准粒子能量给出体系的谱性质，是目前计算带隙最精确的多体方法之一。本页含 GW 理论简介（第 0 节）、例 1 Si 的 $G_0W_0$ 单发计算（三步流程）、例 2 部分自洽 $GW_0$ 并用 Wannier90 插值出 GW 能带结构。

**GW 简介**：在 GW 中电子不再是独立粒子，而是在材料中运动时与其他电子相互作用、极化其环境——被自身诱导的极化云"着装"（dressed）的电子称为**准粒子**（quasiparticle）。准粒子的本征态与本征能由准粒子方程求解：

$$(T+V_{\mathrm{ext}}+V_{\mathrm{h}})\psi_{n\mathbf{k}}(\mathbf{r}) + \int \mathrm{d}\mathbf{r}'\,\Sigma(\mathbf{r},\mathbf{r}',E_{n\mathbf{k}})\psi_{n\mathbf{k}}(\mathbf{r}') = E_{n\mathbf{k}}\psi_{n\mathbf{k}}(\mathbf{r})$$

与 KS 方程相比，交换关联效应由**自能** $\Sigma(\mathbf{r},\mathbf{r}',E_{n\mathbf{k}})$ 引入——它是非局域的、依赖轨道且依赖能量的函数。GW 中自能近似为 $\Sigma = GW$（Green 函数与屏蔽库仑作用之积，方法因此得名）；库仑作用的介电屏蔽通常在 RPA 随机相位近似下计算，即电子与环境的相互作用严格限于极化事件。原则上 $G$ 与 $W$ 都依赖准粒子本征态与本征能，求解是自洽问题，但实践中几乎从不完全自洽迭代，而是分为两类：

- **$G_0W_0$（单发 GW）**：$G$ 与 $W$ 均取自 KS 本征态/本征能，准粒子方程只解一次；
- **$GW_0$（部分自洽）**：只对 Green 函数自洽——首次单发后，用**准粒子**本征能（但保留初始 KS 轨道）重构 $G$，再解新的准粒子方程，如此重复若干次。

### 示例 1：Si 的 $G_0W_0$ 带隙

**教学目标**：(1) 运行单发 GW（$G_0W_0$）计算；(2) 获得准粒子能量与重正化因子 $Z$。

**任务**：用 $G_0W_0$ 计算 Si 的准粒子能量、重正化因子与带隙。单发 GW 分三步：(1) DFT 基态计算；(2) 从 WAVECAR 重启，做一次精确对角化计算额外的未占据 KS 态（用于介电屏蔽），同时计算 KS 轨道胞周期部分对 Bloch 矢量 k 的导数（写入 WAVEDER）；(3) 真正的 $G_0W_0$ 计算，以第 (2) 步的 WAVECAR + WAVEDER 为起点。

**输入要点**（目录 `$TUTORIALS/GW/e01_Si-G0W0`，示例文件内容略）：
- POSCAR：金刚石 Si（$a = 5.430$ Å）；KPOINTS 为 6×6×6 Γ 心网格；`KPOINTS_OPT` 给出能带高对称路径（Γ–X–U/K–Γ–L–W），自洽后做非自洽计算得到沿这些线的色散。
- 基态 INCAR：PBE，`ISMEAR = -5`（四面体法 + Blöchl 修正，平滑 DOS）、`EDIFF = 1E-8`、`ENCUT = 400`。
- **POTCAR 必须用 `_GW` 后缀的势**（如 `PAW_PBE Si_GW`）：GW 势在更宽的能量范围内再现原子散射性质——因为 GW 需要大量未占据态。可用 OUTCAR/POTCAR 的 `TITEL` 字段检查。

**流程与结果**：
1. **基态**：运行后用 py4vasp 求带隙得 PBE 0.708 eV，并画 DOS 与 `KPOINTS_OPT` 路径能带。教程提问：两者带隙为何不一致？——DOS/HOMO-LUMO 差基于规则 k 网格上的本征值，而能带图沿 `KPOINTS_OPT` 的线取样。Si 是**间接带隙**：HOMO 在 Γ，LUMO 在 Γ–X 之间略偏离 X 的点，该点不在规则网格上，所以网格带隙偏大。
2. **未占据态**（`unoccupied-states/`）：从基态复制 WAVECAR。INCAR 要点：`NBANDS = 64`——GW 需要约 **3–4 倍于默认值**的 KS 轨道；`ALGO = Exact` + `NELM = 1`（已收敛，只做一次精确对角化）；`ISMEAR = 0`、`SIGMA = 0.05`（小展宽避免部分占据）；`LOPTICS = T` 使 Bloch 态对 k 的导数写入 WAVEDER。副产品：可用 py4vasp 画独立粒子近似（IPA）的电子介电函数（不含局域场效应）。
3. **$G_0W_0$**（`single-shot/`）：复制上一步的 WAVECAR 与 WAVEDER。INCAR 要点：`ALGO = EVGW0`（"本征值 GW"；VASP.5.X 写作 `GW0`）、`NELMGW = 1`（单次电子更新即 $G_0W_0$；VASP 6.3 之前用 `NELM`）、`NBANDS` 与上一步一致、`NOMEGA = 50`（频率网格，默认值）。
4. **结果读取**：OUTCAR 中 "QP shifts" 段给出每个 k 点每条带的 KS 能量、**QP 能量**、$\sigma(\mathrm{KS})$、$V_{xc}(\mathrm{KS})$、精确交换 $V_x^{\mathrm{pw}}$、**重正化因子 $Z$**（约 0.63–0.78，衡量准粒子权重）、占据数与 $\mathrm{Im}(\Sigma)$。HOMO 在 Γ（k 点 1，QP 5.1235 eV），LUMO 在 X（k 点 13，QP 6.2826 eV）：**$G_0W_0$ 带隙 = 6.2826 − 5.1235 = 1.1591 eV**，相比 PBE 的 0.708 eV 显著改善（实验约 1.17 eV）。
5. **IPA vs RPA 介电函数**：教程提供 `get_dielectric_function.sh` 脚本分别提取 `dielectric_function_IPA.dat` 与 `dielectric_function_RPA.dat`（也可 py4vasp `dielectric_function.plot("Re(IPA, RPA)")`）。RPA 计入行进电子极化环境：从 IPA 到 RPA，$\omega = 0$ 处实部（静态介电常数）增大、虚部谱形变化。VASP 在 GW 中默认用 **RPA** 计算屏蔽库仑作用。

**复习问题**：
1. $GW$ 近似的是什么物理量？（自能 $\Sigma$）
2. VASP 把准粒子能量与重正化因子写在哪里？（OUTCAR）
3. VASP 默认在什么近似级别计算 GW 的屏蔽库仑作用？（RPA）

### 示例 2：Si 的 $GW_0$ 能带结构（含 Wannier90 插值）

**教学目标**：(1) 执行部分自洽的 $GW_0$ 计算；(2) 用 Wannier 化近似得到 GW 能带结构。

**$GW_0$ 计算**（目录 `$TUTORIALS/GW/e02_Si-GW0-band`）：KPOINTS、POSCAR、POTCAR 与例 1 相同。INCAR 与 $G_0W_0$ 的差异只有两点：`NELMGW = 3`（多次电子更新 → 对 $G$ 部分自洽）；`NBANDSGW = 16`——增大实际计算准粒子修正的态数（默认只算费米面附近），为画能带做准备。需要从例 1 第 (2) 步复制 WAVECAR + WAVEDER 后运行。

OUTCAR 中每轮迭代都写 "QP shifts" 段：第一轮第二列是 KS 本征能，后续各列是上一轮的 QP 能量。三次迭代后准粒子带隙：HOMO（Γ）5.0779 eV，LUMO（X）6.2999 eV → **$GW_0$ 带隙 = 1.2220 eV**（对比 $G_0W_0$ 的 1.1591 eV：部分自洽使带隙进一步增大）。

**Wannier90 插值能带**（子目录 `band-structure`）：GW 中 VASP **只能**计算规则 k 网格上的准粒子能量，不能沿任意高对称线计算——GW 能带结构只能近似获得：用 Wannier90 对网格上的 $GW_0$ 准粒子能量做 Fourier 插值（QP 能量存在 GW 计算的 WAVECAR 中）。关键输入：
- `ALGO = None` + `NELM = 1`：VASP 读入文件后直接进入后处理阶段，不做实际计算；
- `LWANNIER90 = T`：后处理阶段写出 Wannier90 输入文件；
- `NUM_WANN = 8`：构造 8 个 Wannier 函数（Si 的 4 条价带 × 2 原子 → sp³ 成键/反键描述）；
- `WANNIER90_WIN`：内容原样复制到 `wannier90.win`——`exclude_bands 17-64`（只用前 16 条带）、投影 `Si:sp3`、解缠结窗口 `dis_win_min = -7`、`dis_win_max = 16`、`dis_num_iter = 100`、`guiding_centres = true`；
- `NBANDS` 须与 GW 计算一致；从 $GW_0$ 的 WAVECAR 重启。

流程：`vasp_std` 生成 `wannier90.*` → `wannier90.x wannier90.win` 构造**最大局域 Wannier 函数** → 检查 `wannier90.wout` 的 Final State（8 个 WF，总展宽约 21.98 Å²，中心位于键中点/原子位，表明 sp³ 成键图像成功）。接着把画图段（`restart = plot`、`bands_plot = true`、`kpoint_path` Γ–X–U/K–Γ–L–W、`bands_num_points 40`、gnuplot 格式）追加到 win 文件再跑一次，用 gnuplot 生成 `band.png`。最后可把 DFT 能带（来自 `e01_Si-G0W0`，经 py4vasp 导出并按 wannier90_band.dat 的横轴标度归一化，能量平移 DFT 费米能 5.628 eV）叠加进同一张图，生成 `band_DFT_vs_GW.png`——直观看到 GW 对带隙与色散的整体修正。

**复习问题**：
1. 为什么 GW 计算的能带结构要用 Wannier 化？（VASP 只能在规则网格上算 QP 能量，任意 k 点靠 Fourier 插值）
2. GW 计算应包含多少条带？（NBANDS 约为默认 3–4 倍，并做收敛测试）

**本部分涉及的 VASP 功能与标签汇总**：
- GW 算法：`ALGO = EVGW0`（本征值 GW）、`NELMGW`（=1 单发 $G_0W_0$，>1 部分自洽 $GW_0$）、`NOMEGA`（频率网格）、`NBANDSGW`（计算 QP 修正的带范围）。
- 三步流程：DFT 基态 → `ALGO = Exact` + `NELM = 1` + 大 `NBANDS` + `LOPTICS = T`（未占据态 + WAVEDER）→ GW（WAVECAR+WAVEDER 起点）。
- 赝势：GW 计算必须用 `_GW` 后缀 PAW 势。
- 输出：OUTCAR "QP shifts" 段（QP 能量、$\sigma$、$V_{xc}$、$Z$、$\mathrm{Im}\,\Sigma$）；IPA/RPA 介电函数脚本与 py4vasp。
- Wannier90 接口：`LWANNIER90 = T`、`NUM_WANN`、`WANNIER90_WIN`（投影、exclude_bands、解缠结窗口）；`ALGO = None` + `NELM = 1` 后处理模式；`wannier90.x` 最大局域化与能带插值（gnuplot 出图）。

## 11. Bethe–Salpeter 方程（BSE）

VASP 可基于 Bethe–Salpeter 方程计算含激子效应的频率依赖介电函数，可在 DFT 基态、杂化泛函或 GW 近似之上进行。

## 11.1 BSE（一）：金刚石碳的光吸收（GW+BSE、TDDDH）

> 来源：<https://vasp.at/tutorials/latest/bse/part1/> ｜ 练习文件：[BSE-part1.zip](https://vasp.at/tutorials/latest/BSE-part1.zip) ｜ 示例目录：`$TUTORIALS/BSE/e01_C_GW`、`e02_C_BSE`、`e03_C_TDHF`

VASP 可基于 Bethe–Salpeter 方程（BSE）计算含激子效应的频率依赖介电函数。本页含 BSE 形式理论简介（第 0 节）、例 1 金刚石 $G_0W_0$ 准备计算、例 2 TDA 下的 BSE 光吸收、例 3 用介电依赖杂化泛函（TDDDH）绕过 GW 的等价计算。体系为金刚石碳——介电屏蔽较强（$\epsilon_\infty \approx 5.7$）、激子弱束缚（Wannier–Mott 型），为 Part 2 的强束缚 Frenkel 激子（LiF）做铺垫。

**BSE 形式理论**：激发能是如下线性问题的本征值 $\omega_\lambda$：

$$\left(\begin{array}{cc} \mathbf{A} & \mathbf{B} \\ \mathbf{B}^* & \mathbf{A}^* \end{array}\right)\left(\begin{array}{l} \mathbf{X}_\lambda \\ \mathbf{Y}_\lambda \end{array}\right)=\omega_\lambda\left(\begin{array}{cc} -\mathbf{1} & \mathbf{0} \\ \mathbf{0} & \mathbf{1} \end{array}\right)\left(\begin{array}{l} \mathbf{X}_\lambda \\ \mathbf{Y}_\lambda \end{array}\right)$$

矩阵 $A_{vc}^{v'c'} = (\varepsilon_v-\varepsilon_c)\delta_{vv'}\delta_{cc'} + \langle cv'|V|vc'\rangle - \langle cv'|W|c'v\rangle$ 与 $B_{vc}^{v'c'} = \langle vv'|V|cc'\rangle - \langle vv'|W|c'c\rangle$ 描述占据态 $v$ 与未占据态 $c$ 之间的跃迁；相互作用含裸库仑 $V$ 与屏蔽库仑 $W$ 两项，能量与轨道通常取自 GW 计算。$A$ 的三项分别记为对角项 $D$、交换项 $K^x$（来自 $V$，排斥）与直接项 $K^D$（来自 $W$，吸引，形成激子）。$W$ 由 RPA 介电函数屏蔽：$\epsilon^{-1}_{\mathbf{G},\mathbf{G}'}(\mathbf{q},\omega) = \delta_{\mathbf{G},\mathbf{G}'} + \frac{4\pi}{|\mathbf{q}+\mathbf{G}|^2}\chi^{\mathrm{RPA}}_{\mathbf{G},\mathbf{G}'}(\mathbf{q},\omega)$；虽然 $W$ 依赖频率，实际 BSE 计算的标准做法是**静态近似** $W(\mathbf{q},\omega=0)$。**Tamm–Dancoff 近似（TDA）**忽略激发与退激发的耦合（$B$、$B^*$），退化为 $A X_\lambda = \omega_\lambda X_\lambda$。宏观介电函数由 BSE 本征值/本征矢构造：

$$\epsilon_M(\mathbf{q},\omega) = 1 + \lim_{\mathbf{q}\to 0} v(q)\sum_\lambda \left|\sum_{cvk}\langle c\mathbf{k}|e^{i\mathbf{qr}}|v\mathbf{k}\rangle X_\lambda^{cv\mathbf{k}}\right|^2 \left(\frac{1}{\omega_\lambda-\omega-i\delta} + \frac{1}{\omega_\lambda+\omega+i\delta}\right)$$

光学极限 $q\to 0$ 下 Taylor 展开给出 $\langle c\mathbf{k}|e^{i\mathbf{qr}}|v\mathbf{k}\rangle \propto \mathbf{q}\langle c\mathbf{k}|\nabla|v\mathbf{k}\rangle$（偶极近似）。TDDFT 的 Casida 形式给出同样的线性问题，但电子-空穴相互作用由交换关联核 $f_{\mathrm{xc}} = \delta v_{\mathrm{xc}}/\delta\rho$ 描述；本教程采用的核为 $\langle cv'|f_{\mathrm{xc}}|c'v\rangle = \gamma\langle cv'|V|c'v\rangle + (1-\gamma)\langle cv'|f^{\mathrm{ALDA}}_{\mathrm{x}}|vc'\rangle + \langle cv'|f^{\mathrm{ALDA}}_{\mathrm{c}}|vc'\rangle$，其中 $\gamma$ 描述吸引作用的屏蔽、由杂化泛函的 Fock 交换比例决定。BSE/TDHF 计算需要前置基态提供准粒子轨道、轨道对 k 的导数与能量（VASP 存于 WAVECAR/WAVEDER）；BSE 额外需要 GW 产生的静态屏蔽势（存于 `WFULLxxxx.tmp`/`Wxxxx.tmp`）。

### 示例 1：金刚石的前置 $G_0W_0$ 计算

**教学目标**：(1) 设置单发 $G_0W_0$ 计算；(2) 求准粒子（QP）带隙；(3) 画出 GW 计算得到的 RPA 介电函数。

**任务**：基于 PBE 基态对金刚石做 $G_0W_0$。GW 收敛需仔细检查未占据态数、k 点数与频率点数；小原胞时用精确对角化更快。三步流程与 §10 相同：DFT 基态 → 精确对角化算未占据带 → $G_0W_0$。

**输入要点**（目录 `e01_C_GW`，示例文件内容略）：
- POSCAR：金刚石原胞，取**实验晶格常数** 3.567 Å；KPOINTS 6×6×6（k 网格虽小，但声称带隙已收敛到 0.1 eV 以内）。
- `ground-state/INCAR`：PBE，`ISMEAR = 0`、`SIGMA = 0.01`、`ENCUT = 450`。
- `unoccupied-states/INCAR`：`NBANDS = 64`（含 60 条空带）、`ALGO = Exact` + `NELM = 1`（单次精确对角化）、`LOPTICS = T`（轨道导数写入 WAVEDER）、`CSHIFT = 0.6`（介电函数展宽）。
- $G_0W_0$ INCAR：`ALGO = G0W0`（本教程 VASP 6 写法，等价于本征值 GW）、`NELM = 1`（单发）、`KPAR = 4`（k 点并行）、`NOMEGA = 50`（默认 100，此处 50 已满足 0.1 eV 收敛标准）、`NBANDS` 与上一步一致。
- POTCAR：`_GW` 推荐势。

**流程与结果**：
1. 基态后区分两个带隙概念（教程考点）：**直接带隙**为同一 k 点最高占据与最低未占据带之差；**基本（fundamental）带隙**为 VBM 与 CBM 之差（可在不同 k 点）。金刚石：PBE 基本带隙 4.16 eV、直接带隙 5.60 eV；实验基本带隙 5.85 eV（已做零声子重正化 ZPR 修正——ZPR 使实验带隙"打开"，计算不含此效应，故比较时须修正实验值）。
2. `unoccupied-states/`：复制 WAVECAR 运行，检查 OUTCAR 带数正确、WAVEDER 存在；`LOPTICS` 副产品是 IPA 介电函数（py4vasp `dielectric_function.plot("Im")`）。
3. $G_0W_0$：复制 WAVECAR + WAVEDER 运行。OUTCAR "QP shifts" 段给出 QP 本征值：**QP 直接带隙 7.40 eV、基本带隙 5.58 eV**——与实验 5.85 eV 相当接近。再画 IPA 与 RPA 介电函数对比（RPA 来自 GW 计算本身，含屏蔽）。

**复习问题**：
1. 直接带隙与基本带隙的区别？
2. `LOPTICS` 算的介电函数与 GW 计算中得到的介电函数有何区别？（IPA vs RPA）

### 示例 2：BSE 光吸收（TDA）

**教学目标**：执行含激子效应的光吸收计算。

**任务**：用上一步的 QP 能量与屏蔽势计算光吸收谱。BSE 计算流程：用 QP 轨道与 RPA 屏蔽势构造 BSE 哈密顿量 → 解本征问题（式 8）→ 由式 (9) 得 BSE 介电函数。

**输入要点**（目录 `e02_C_BSE`，示例文件内容略）：核心是 `ALGO = BSE`；`ANTIRES = 0` 启用 Tamm–Dancoff 近似；`NBANDSO = 4`（BSE 中的价带数）与 `NBANDSV = 4`（BSE 中的导带数）——对全部 4 条占据带与 4 条空带构造哈密顿量，对应最高跃迁能约 37 eV；`NBANDS = 64` 必须与前置 GW 一致；`CSHIFT = 0.6`。

**流程与结果**：从 `e01_C_GW` 复制 **WAVECAR、WAVEDER 与 W\*.tmp 文件**（屏蔽势！缺一不可）后运行。vasprun.xml 给出：
- `opticaltransitions`：各激发能及振子强度——首个亮激发在 6.864 eV（振子强度 8.456，三重简并），8.706 eV 处还有一个强峰（15.968）；若干跃迁振子强度为零（暗态）。激子效应使吸收边相对 QP 带隙（直接 7.40 eV）红移。
- `dielectricfunction`：BSE 介电函数（含全部张量分量），可用 awk 脚本取 xx/yy/zz 平均输出 `optics.dat`。
- OUTCAR 检查点：`Bands included in the BSE ... VB(min)=1 VB(max)=4 CB(min)=5 CB(max)=8` 确认带范围；`BSE (scaLAPACK) attempting allocation ... rank=3456` 给出 BSE 矩阵秩。**技巧**：若 BSE 介电函数全为零，检查 stdout 中 "the WAVECAR file was read successfully"——WAVECAR/WAVEDER 读取失败是常见事故。

最后用 py4vasp 把 IPA、RPA、BSE 三条虚部谱叠在一张图上：BSE 相对 IPA/RPA 出现激子峰与谱形重排。

**复习问题**：
1. 电子-空穴相互作用如何影响介电常数？

### 示例 3：TDDDH 光吸收（绕过 GW）

**教学目标**：(1) 为模型介电函数确定参数；(2) 用含激子效应的 TDDFT（TDDDH）计算光吸收谱。

**原理**：GW 步计算量大；用杂化泛函路线可达到相近精度并绕开 GW。**介电依赖杂化泛函（DDH）**用模型函数描述 Fock 交换的屏蔽：

$$\epsilon^{-1}(|\mathbf{q}+\mathbf{G}|) = 1 - (1-\epsilon_\infty^{-1})\,e^{-|\mathbf{q}+\mathbf{G}|^2/4\mu^2}$$

同一函数也可用于屏蔽电子-空穴间的吸引库仑作用，即 **TDDDH**（Phys. Rev. Research 2, 032019 (2020)）。本例用 DDH 做基态，两步：DDH 基态 → DDH 精确对角化算未占据带。

**输入要点**（目录 `e03_C_TDHF`，示例文件内容略）：
- **参数拟合**：用 awk 从 `e01_C_GW/vasprun.xml` 提取 RPA 介电函数对角元 `eps_diag.dat`，将模型函数拟合到 RPA 数据得 $\epsilon_\infty = 5.7$、$\mu \approx 1.70$–1.79 Å$^{-1}$（教程用 `eps_model` 函数绘图验证拟合质量）。
- `ground-state/INCAR`：`LHFCALC = T`、`LMODELHF = T`（范围分离杂化泛函）、`HFSCREEN = 1.70`（范围分离参数，即 $\mu$）、`AEXX = 0.175`（精确交换比例 ≈ $1/\epsilon_\infty$）。
- `unoccupied-states/INCAR`：同上 + `NBANDS = 8`（4 条未占据带）、`ALGO = Exact`、`NELM = 1`、`LOPTICS = T`。
- TDDDH INCAR：`LMODELHF = T`、`HFSCREEN = 1.70`、`AEXX = 0.175`、`LFXC = T`（启用局域 xc 核）、`ALGO = TDHF`、`ANTIRES = 0`（TDA）、`NBANDSO = NBANDSV = 4`、`CSHIFT = 0.6`、`NBANDS = 8` 与基态一致。

**流程与结果**：DDH 基态（复制 PBE WAVECAR 加速杂化收敛）→ OUTCAR 给出 DDH 直接带隙 7.5 eV、基本带隙 5.6 eV（与 $G_0W_0$ 的 7.40/5.58 几乎一致——这就是 DDH 的价值）；DDH 未占据带得到 WAVECAR+WAVEDER；TDDDH 光吸收（检查 OUTCAR 中 BSE 哈密顿量秩与例 2 相同）。最终把 IPA/RPA/BSE/TDDDH 四条谱叠图：**TDDDH 与 BSE 谱高度吻合**。本例 6×6×6 稀疏网格不足以显示 RPA 与含 e-h 相互作用谱的差别，但激子效应的重要性可见于 Phys. Rev. Lett. 107, 186401 (2011) 图 2。金刚石介电屏蔽强（~5.7），激子弱束缚，属 Wannier–Mott 激子。

**复习问题**：
1. 如何确定式 (14) 模型介电函数的参数？（拟合 GW 计算的 RPA 介电函数）
2. TDDDH 相对 BSE 的主要优势？（省去昂贵的 GW 步，精度相近）

**本部分涉及的 VASP 功能与标签汇总**：
- GW 准备链（BSE 版）：PBE 基态 → `ALGO = Exact` + `NELM = 1` + `NBANDS = 64` + `LOPTICS = T` + `CSHIFT` → `ALGO = G0W0` + `NELM = 1` + `KPAR` + `NOMEGA = 50`（`_GW` 势）。
- BSE：`ALGO = BSE`、`ANTIRES = 0`（TDA）、`NBANDSO`/`NBANDSV`（BSE 活性价/导带）、静态 $W$（W\*.tmp 文件必须随 WAVECAR/WAVEDER 一并复制）；输出：vasprun.xml 的 opticaltransitions（激发能+振子强度）与 dielectricfunction；OUTCAR 的 BSE 带范围与矩阵秩。
- TDDDH：`LMODELHF = T`（模型介电屏蔽的范围分离杂化）、`HFSCREEN`（$\mu$）、`AEXX`（≈$1/\epsilon_\infty$）、`LFXC = T`（局域 xc 核）、`ALGO = TDHF`。
- 概念：直接 vs 基本带隙、ZPR 修正、TDA、光学极限偶极近似、Wannier–Mott vs Frenkel 激子、IPA→RPA→BSE 谱的逐层对比。

## 11.2 BSE（二）：LiF 光吸收——逐层近似对比（IP / RPA / TDA / full BSE）

> 来源：<https://vasp.at/tutorials/latest/bse/part2/> ｜ 练习文件：[BSE-part2.zip](https://vasp.at/tutorials/latest/BSE-part2.zip) ｜ 示例目录：`$TUTORIALS/BSE/e04_LiF_GW`、`e05_LiF_IP`、`e06_LiF_RPA`、`e07_LiF_BSE`、`e08_LiF_FBSE`

与金刚石不同，LiF 的电子间相互作用屏蔽很弱（$\epsilon_\infty = 1.9$），导致强束缚的 **Frenkel 激子**。本页用同一 $G_0W_0$ 起点，通过 `LHARTREE`/`LADDER`/`ANTIRES` 三个开关逐级开启 BSE 哈密顿量的各项，对比四个近似层级的光吸收谱并与实验对照——这是理解 BSE 中每一项物理贡献的最佳教学案例。

### 示例 4：LiF 的前置 $G_0W_0$ 计算

**教学目标**：(1) 设置单发 $G_0W_0$ 计算；(2) 求 QP 带隙。

**任务**：基于 PBE 电子结构计算 LiF 原胞的 $G_0W_0$ 能带结构与屏蔽库仑势。流程与 Part 1 金刚石完全相同（三步：DFT 基态 → 精确对角化未占据带 → $G_0W_0$），详见 Part 1。

**输入要点**（目录 `e04_LiF_GW`，示例文件内容略）：POSCAR 为岩盐 LiF 原胞（实验晶格常数 4.026 Å）；KPOINTS 6×6×6；`ENCUT = 500`；`ground-state`：PBE；`unoccupied-states`：`NBANDS = 64`、`ALGO = Exact`、`NELM = 1`、`LOPTICS = T`；$G_0W_0$：`ALGO = EVGW0`（VASP.5.X 写作 `GW0`）、`NELM = 1`、`KPAR = 4`、`NOMEGA = 50`（默认 100）、`NBANDS` 一致；`_GW` 推荐势。

**结果**：OUTCAR 中 LiF 的**直接带隙约 13.44 eV**。教程提问：LiF 是直接还是间接带隙？其 RPA 介电常数是多少？（弱屏蔽，$\epsilon_\infty \approx 1.9$）

### 示例 5：独立粒子近似（BSE-IP）

**教学目标**：计算独立粒子近似下的光吸收。

**原理**：IP 近似下 `LHARTREE = .FALSE.`、`LADDER = .FALSE.`，BSE 哈密顿量只剩对角项（$A = D$），粒子间相互作用完全忽略。对屏蔽强、电子-空穴束缚能小的体系，这一近似对光吸收的描述尚算合理——但对 LiF 不行。

**输入要点**（目录 `e05_LiF_IP`，示例文件内容略）：`ALGO = BSE`、`ANTIRES = 0`（TDA）、`NBANDSO = NBANDSV = 5`（各 5 条价/导带进入 BSE）、`CSHIFT = 0.2`、`OMEGAMAX = 40`（谱的最大能量），加上述两个 `.FALSE.` 开关。从 `e04_LiF_GW` 复制 WAVECAR、WAVEDER、W\*.tmp 后运行。

**结果**：与实验谱对比，符合很差——吸收边在 13.4 eV 处，**吸收边上的强激子峰完全缺失**，说明电子与空穴未形成束缚态、跃迁能未被拉低。该谱应等价于用同一套轨道与能量做 `LOPTICS = .TRUE.` 的结果。

**复习问题**：(1) IP 近似是否包含激子效应？(2) 对什么体系 IP 可以是好近似？

### 示例 6：BSE-RPA

**教学目标**：基于 $G_0W_0$ 准粒子计算 LiF 的 RPA 谱。

**原理**：RPA 层次不算梯形图（`LADDER = .FALSE.`），只含交换相互作用（`LHARTREE = .TRUE.`），即 $A = D + K^x$。

**输入要点**（目录 `e06_LiF_RPA`，示例文件内容略）：与例 5 相同，仅 `LHARTREE = .TRUE.`。

**结果**：与实验相比**并未改善**，谱反而略**蓝移**——这是排斥性交换相互作用 $K^x$ 的效果。注意这里的 RPA 只是 TDA 下的 bubble 图（$B$、$B^*$ 被忽略），因此不等价于 `ALGO = CHI` 的结果（后者含激发-退激发 bubble 的耦合）。RPA 与 IP 有同样的缺陷——没有激子；但 RPA 计入了与等离子激元激发的相互作用，对 **EELS 谱**往往给出好的描述。

**复习问题**：(1) RPA 是否包含激子效应？(2) 对什么体系 RPA 是好近似？

### 示例 7：BSE-TDA（首次出现激子）

**教学目标**：在 Tamm–Dancoff 近似下计算光吸收。

**原理**：LiF 是宽隙半导体且介电常数小（$\epsilon_\infty = 1.9$），屏蔽弱 → 电子-空穴吸引强 → 激子强局域、束缚能大（0.1–1 eV）。强 e-h 相互作用是碱卤化物与有机分子体系的典型特征。TDA 下 bubble 图与梯形图都计入（`LHARTREE = .TRUE.`、`LADDER = .TRUE.`），但 $B$、$B^*$ 被忽略，哈密顿量为 $A = D + K^x + K^D$——吸引的直接项 $K^D$（梯形图）正是激子的来源。

**输入要点**（目录 `e07_LiF_BSE`，示例文件内容略）：与例 6 相同，再加 `LADDER = .TRUE.`。

**结果**：谱与实验**显著改善**——首个强激子峰被正确捕获（吸收边以下出现尖锐激子线）。但部分特征的强度仍与实验不符，Part 3 将用更密的 k 网格解决这一问题。

**复习问题**：(1) TDA 是否包含激子效应？(2) TDA 忽略了式 (1) 中的哪些项？

### 示例 8：超越 TDA 的完整 BSE

**教学目标**：计算超越 TDA 的光吸收谱。

**原理**：完整 BSE 包含共振-反共振耦合 $H^{\mathrm{BSE}} = \begin{pmatrix} A & B \\ B^* & A^* \end{pmatrix}$。`ANTIRES = 2` 时 $B$、$B^*$ 进入计算，原则上需对角化秩为 $2N_cN_vN_k$ 的矩阵；但 VASP 利用**时间反演对称性**把问题化为两个秩 $N_cN_vN_k$ 的方程——因此完整 BSE 的计算量只约为 TDA 的**两倍**。

**输入要点**（目录 `e08_LiF_FBSE`，示例文件内容略）：与例 7 相同，仅 `ANTIRES = 2`。

**结果**：把实验谱与 BSE-IPA、BSE-RPA、BSE-TDA、BSE-FULL 五条谱叠图，完整呈现逐层近似的效应链：IP 只有吸收边（13.4 eV）→ RPA 蓝移、无激子 → TDA 出现强激子峰 → full BSE 进一步精修。对 LiF 这类宽隙绝缘体，TDA 与 full BSE 差异有限；反共振项在金属性/小隙体系或高能区更重要。

**复习问题**：(1) TDA 与完整 BSE 的区别？(2) 什么时候需要超越 TDA？

**本部分涉及的 VASP 功能与标签汇总**：
- BSE 哈密顿量开关：`LHARTREE`（bubble/交换 $K^x$ 图）、`LADDER`（梯形/直接 $K^D$ 图）、`ANTIRES`（0=TDA，2=含 $B$/$B^*$ 的完整 BSE，时间反演对称性使成本仅约翻倍）。
- 四层近似配方：IP = 两者 FALSE；RPA = LHARTREE 开、LADDER 关；TDA = 两者开 + ANTIRES=0；full BSE = 两者开 + ANTIRES=2。
- 谱控制：`OMEGAMAX`（能量上限）、`CSHIFT`（展宽）、`NBANDSO`/`NBANDSV`（活性带）。
- 物理对照：排斥 $K^x$ 导致蓝移；吸引 $K^D$ 产生激子峰（Frenkel 型，束缚能 0.1–1 eV）；RPA 适合 EELS；`LOPTICS` 谱与 BSE-IP 等价；TDA-RPA 不等价于 `ALGO = CHI`。

## 11.3 BSE（三）：高效布里渊区采样与激子分析

> 来源：<https://vasp.at/tutorials/latest/bse/part3/> ｜ 练习文件：[BSE-part3.zip](https://vasp.at/tutorials/latest/BSE-part3.zip) ｜ 示例目录：`$TUTORIALS/BSE/e09_LiF_TDHF_shifted_grid`、`e10_C_FATBAND`

BSE/TDHF 的矩阵秩为 $N_k \times N_v \times N_c$，k 点一多计算量急剧上升。本页教两个进阶技能：例 9 用**平移网格**以 8 个 4×4×4 粗网格计算加权平均，逼近 16×16×16 密网格的 TDDDH 吸收谱（LiF）；例 10 用 `NBSEEIG`/`BSEFATBAND` 文件对单个激子态做 fatband 分析（金刚石）。

### 示例 9：高效布里渊区采样（平移网格平均）

**教学目标**：用平移网格加速收敛的光吸收计算。

**方法**（三步）：
1. 用 4×4×4 网格的 DFT 计算确定不可约 k 点 $\mathbf{k}_n^{ir}$ 与权重 $W_n$；
2. 把 TDDDH 计算分别平移到每个 $\mathbf{k}_n^{ir}$，各自提取介电函数；
3. 按 $W_n$ 加权求和得到最终介电函数。

即：用 8 个平移的 4×4×4 网格计算逼近一个 16×16×16 网格的 LiF TDDDH 谱。

**输入要点**（目录 `e09_LiF_TDHF_shifted_grid`，示例文件内容略）：三份 INCAR 串联——`INCAR.DFT`（PBE，`ISYM = -1` 关闭布里渊区对称性——平移网格需要显式看到全部 k 点；`NBANDS = 12`、`KPAR = 4`）；`INCAR.DDH`（DDH 基态：`LHFCALC`、`LMODELHF = T`、`HFSCREEN = 1.47`、`AEXX = 0.53`——LiF 弱屏蔽故精确交换比例高达 0.53，对比金刚石 0.175；`LOPTICS = T`、`CSHIFT = 0.2`）；`INCAR.TDHF`（`ALGO = TDHF`、`ANTIRES = 0`、`LFXC = T`、`NBANDSO = NBANDSV = 5`、`OMEGAMAX = 25`、`NEDOS = 2001`）。KPOINTS 为 4×4×4 Γ 心。

**流程与结果**：
1. DFT 运行后，用 `grep -A11 "irreducible k-points:" OUTCAR | tail -8 > POINTS` 收集不可约 k 点与权重。
2. `run.sh` 脚本循环 8 个平移：每个子目录把平移量写入 KPOINTS 的 shift 行，依次跑 DFT→DDH→TDHF，`extract_optics.sh` 提取 optics.dat 并乘以该 k 点权重，最后 `paste` 合并、按总权重（24）归一化。
3. 与实验谱对比时须做**零声子重正化（ZPR）修正**：LiF 的 ZPR 高达 **1.15 eV**，计算不含该效应，故把 BSE 谱整体平移 −1.15 eV 后与实验吻合良好。
4. **两个重要提醒**：(a) 激子束缚能随相邻 k 点间距 $\mathbf{q} = \mathbf{k}-\mathbf{k}'$ **线性**收敛，多网格平均对谱形是合理近似，但会**高估激子束缚能**——精确束缚能必须用显式布里渊区采样；(b) 另一种改善 k 收敛的常用手段是小非对称平移 $\mathbf{s} = (1/12, 3/12, 5/12)$——打破规则网格对称性、纳入最多不等价 k 点。

**复习问题**：多平移网格近似能否用于精确的激子束缚能？（不能，会高估）

### 示例 10：激子分析（BSEFATBAND 与 fatband 图）

**教学目标**：检查单个激子态的性质——把 BSE/TDHF 本征矢写成 fatband 图（耦合系数以圆圈画在对应能带与 k 点上）。

**原理与输出格式**：VASP 可把最低的 `NBSEEIG` 个本征矢写入 `BSEFATBAND` 文件，结构为：首行 BSE 矩阵秩与 NBSEEIG；每个本征态一行头（`nBSE eigenvalue E_BSE IP-eigenvalue E_IP`），随后每行给出 k 点坐标、空穴本征值 $E_h$、电子本征值 $E_e$、$|X_{\mathrm{BSE}}|/W_k$（本征矢分量模除以 k 点权重）、空穴/电子带号 $NB_h$/$NB_e$、以及 $X_{\mathrm{BSE}}$ 的实部与虚部。

**任务 10.1 DDH 基态**（目录 `e10_C_FATBAND/DFT`、`DDH`）：金刚石（3.567 Å），KPOINTS 加密到 **8×8×8**；DFT 为 PBE；DDH 步 `HFSCREEN = 1.26`、`AEXX = 0.175`、`NBANDS = 8`（4 条空带）、`LOPTICS = T` 写 WAVEDER。先跑 DFT，复制 WAVECAR 后跑 DDH，并求金刚石的直接/基本带隙。

**任务 10.2 TDDDH 与 fatband**（目录 `e10_C_FATBAND/`）：INCAR 与 Part 1 的 TDDDH 类似（`ALGO = TDHF`、`ANTIRES = 0`、`LFXC = T`、`NBANDSO = NBANDSV = 4`），关键新增 `NBSEEIG = 10`——写出前 10 个激子本征矢。复制 DDH 的 WAVECAR+WAVEDER 运行。

**激子态定位**：先在 vasprun.xml 的 opticaltransitions 中找首个振子强度非零的跃迁——本例是第 4 个态（7.4057 eV，振子强度 8.077；前 3 个是暗态）。再到 BSEFATBAND 中查 `4BSE`：其 BSE 本征值 7.4057 eV、IP 本征值 7.8231 eV → 首个亮激子的束缚能 = 7.823 − 7.406 ≈ **0.42 eV**（因稀疏 k 网格而被显著高估——呼应例 9 的警告）。

**fatband 作图**：`fatband.sh` 提取沿 Γ–L 与 Γ–X 方向的数据（坐标：|k 点| × 空穴/电子本征值，圆圈半径 = $|X_{\mathrm{BSE}}|$），gnuplot（`plot_fatbands.gnuplot`）或教程给出的 plotly 脚本（空心圆环，蓝/红两色）出图。结论：最大贡献来自 Γ 点直接带隙附近的带顶，邻近 k 点贡献很小——说明该激子态在 **k 空间高度局域**（Wannier–Mott 激子但受强 e-h 吸引约束于 Γ 附近跃迁）。

**复习问题**：BSEFATBAND 文件包含哪些信息？

**本部分涉及的 VASP 功能与标签汇总**：
- 平移网格法：`ISYM = -1`（关对称性）+ 从 OUTCAR 提取不可约 k 点与权重（`irreducible k-points` 段）+ KPOINTS shift 行 + 加权平均脚本；非对称平移 $\{1/12, 3/12, 5/12\}$ 技巧。
- 激子分析：`NBSEEIG`（写出本征矢个数）、`BSEFATBAND` 文件（$E_{\mathrm{BSE}}$/$E_{\mathrm{IP}}$、$E_h$/$E_e$、$X_{\mathrm{BSE}}$、带指标）、fatband 可视化（gnuplot/plotly）；束缚能 = $E_{\mathrm{IP}} - E_{\mathrm{BSE}}$（k 网格须足够密）。
- 谱对比修正：LiF 零声子重正化 1.15 eV；`NEDOS` 控制谱点数。

## 12. X 射线吸收谱（XAS）

X 射线吸收谱涉及芯电子激发到导带，留下芯空穴与可形成束缚态（激子）的激发电子。VASP 有两种建模途径：基于 DFT 的超胞核空穴（SCH）方法，以及基于多体微扰理论、需求解 BSE 的方法。

## 12.1 XAS（一）：LiCl 的 X 射线吸收谱（超胞核空穴 / 全核空穴）

> 来源：<https://vasp.at/tutorials/latest/xas/part1/> ｜ 练习文件：[XAS-part1.zip](https://vasp.at/tutorials/latest/XAS-part1.zip) ｜ 示例目录：`$TUTORIALS/XAS/e01_SCH`、`e02_FCH`

X 射线吸收谱（XAS）中，X 射线把芯电子激发到导带，留下芯空穴与可形成束缚态（激子）的激发电子。VASP 的两种建模途径中，本页讲基于 DFT 的 ΔSCF 类方法：例 1 超胞核空穴（SCH）方法计算 LiCl 中 Li 的 K 边 XAS（含超胞收敛与芯空穴效应可视化），例 2 全核空穴（FCH）与激发电子核空穴（XCH）方法及与实验对比的规范做法。

**理论基础**：由于 X 射线波长远大于固体中的特征动量，从长波极限（$\mathbf{q}=0$）的介电函数虚部出发：

$$\epsilon^{(2)}_{\alpha\beta}(\omega, \mathbf{q}=0) = \frac{4\pi^2 e^2\hbar^2}{\Omega\omega^2 m_e^2}\sum_{c,v,\mathbf{k}} 2w_{\mathbf{k}}\,\delta(\varepsilon_{c\mathbf{k}}-\varepsilon_{v\mathbf{k}}-\omega)\,M^{v\to c}_\alpha {M^{v\to c}_\beta}^{*}$$

它正比于吸收谱。PAW 形式下的动量矩阵元含全电子部分波 $\phi_i$、赝部分波 $\tilde\phi_i$ 与投影符 $\tilde p_i$。XAS 只考虑单个原子位点的芯空穴，利用完备性关系 $\sum_i|\tilde p_i\rangle\langle\tilde\phi_i| = 1$ 可把矩阵元限制在单站点：

$$M^{\mathrm{core}\to c\mathbf{k}}_\alpha = \sum_i \langle\tilde\psi_{c\mathbf{k}}|\tilde p_i\rangle \langle\phi_i|\nabla_\alpha|\phi_{\mathrm{core}}\rangle$$

且求和可限制到导带 $c$（芯能级作为唯一的初态）。

### 示例 1：超胞核空穴（SCH）方法

**教学目标**：(1) 设置超胞核空穴计算；(2) 掌握超胞尺寸收敛；(3) 理解激子效应在 XAS 中的重要性。

**任务**：用 SCH 方法计算 LiCl 中 Li 的 K 边 XAS。采用终态近似（final-state approximation，`ICORELEVEL = 2`，Phys. Rev. B 63, 205419 (2001)）。

**输入要点**（目录 `e01_SCH`，示例文件内容略）：
- POSCAR 提供 1×1×1、2×2×2、3×3×3、4×4×4 四套超胞（LiCl 岩盐，$a = 5.106$ Å）。**关键技巧**：芯激发原子须作为**独立化学种**处理——例如 8 Li + 8 Cl 的 2×2×2 超胞写成 `Li Li Cl` = 1 7 8，因此 POTCAR 需要三份（`cat POTCAR_Li POTCAR_Li POTCAR_Cl > POTCAR`）。
- INCAR 核心标签：`ICORELEVEL = 2`（终态近似）、`CLNT = 1`（对哪个化学种打芯空穴——永远只选该种的一个原子，保证超胞内只有一个芯空穴）、`CLN = 1`/`CLL = 0`（主量子数 n、角量子数 l；K 边即 1s）、`CLZ = 1`（激发芯原子的电子数——SCH/XCH 中从芯能级拿走 1 个电子并自动放进导带）、`CH_LSPEC = .TRUE.`（按角动量分解谱）、`CH_SIGMA`（谱展宽）、`CH_NEDOS = 5000`（谱能量网格点数）、`CH_AMPLIFICATION`（谱放大倍数，取超胞内 Li 原子数——因为谱被超胞体积除，不放大则强度随超胞增大骤减）、`NBANDS` 取占据带数的两倍以保证足够空带。
- KPOINTS：Γ 单点（`vasp_gam` 运行）——k 点数也是需要收敛的量。

**流程与结果**：
1. `run_conv.sh` 对 2×2×2（16 原子）、3×3×3（54 原子）、4×4×4（128 原子）逐个计算，py4vasp 用 `dielectric_function.plot("Im(XAS)")` 叠图：**峰的劈裂与相对强度随超胞增大而变化**。为什么必须用超胞？——最小化芯空穴与其周期镜像的相互作用。增大超胞顺带改善了布里渊区的 k 点表示，但这不是主因：芯空穴间相互作用足够小后，直接增加 k 点数比继续加大超胞更便宜。
2. **激子效应可视化**：`ICORELEVEL = 1` 是初态近似——价电子不随芯空穴弛豫、不含激子效应。对比 4×4×4 超胞的 CH（终态）与 noCH（初态）谱：加入芯空穴后谱向低能移动、态更局域（首个峰尤甚）——全是电子-空穴相互作用（激子效应）的结果。
3. **部分电荷密度对照**：两步流程（step1 自洽 → step2 用 `LPARD = T` + `IBAND = 257` 提取最低导带的部分电荷；CH 版 step1 用 `ICORELEVEL = 2` + 88 条带做完全弛豫的芯空穴基态）。py4vasp `partial_density.to_ngl(isolevel=1)` 三维可视化：**有芯空穴时**最低未占据态局域在持芯空穴的原子周围；**无芯空穴时**完全离域（降低等值面可看得更清楚）。

**复习问题**：
1. 为什么 XAS 的相互作用可以限制在单一位点？
2. 为什么要增大 NBANDS？
3. 为什么含芯空穴的谱计算必须用超胞？
4. INCAR 中 CLNT、CLN、CLL、CLZ 各做什么？
5. 为什么激子效应在 XAS 中如此重要？

### 示例 2：全核空穴（FCH）与激发电子核空穴（XCH）

**教学目标**：(1) 执行全核空穴计算；(2) 掌握计算谱与实验对比的呈现规范。

**原理**：上一例的 SCH（`ICORELEVEL = 2` + `CLZ = 1`）把芯电子拿走 1 个并自动放入导带，即**激发电子-芯空穴（XCH）方法**（Phys. Rev. Lett. 96, 215502 (2006)）。另一种做法是把激发电子放回均匀背景电荷、不在自洽中显式处理激发电子分布，即**全核空穴（FCH）方法**（J. Chem. Phys. 120, 8632 (2004)）。许多情况下 **FCH 比 XCH 更接近实验**（原因详见 Unzog et al., Phys. Rev. B 106, 155133 (2022) 第 V 部分）。

**输入要点**（目录 `e02_FCH`，示例文件内容略）：4×4×4 超胞（128 原子），`INCAR.XCH` 与例 1 类似（`ICORELEVEL = 2`、`CLNT/CLN/CLL/CLZ`、`CH_SIGMA = 0.1`、`CH_NEDOS = 5000`、`CH_AMPLIFICATION = 64`、`NBANDS = 648`、`NCORE = 2`）；`INCAR.FCH` 唯一差别是加 `NELECT = 512`（比中性少 1 个电子）。

**流程与结果**：
1. **FCH 的 NELECT 确定法**：VASP 会自动用均匀背景电荷补偿带电晶胞，所以只需把电子数减一。先 `vasp_std --dry-run`，`grep NELECT OUTCAR` 得到 513，再在 INCAR.FCH 中写 `NELECT = 512`。
2. `run_job.sh` 在 4×4×4 超胞 Γ 点分别跑 FCH 与 XCH。
3. **与实验对比的规范**（本例核心教学点）：芯电子能量的计算值总是显著偏离实验，因此谱需要**平移**——常见做法是把实验第一峰移到零、把计算谱平移使第一峰极大值与实验重合（选哪个峰是任意的）；强度也常按任意单位处理（谱仪展宽等参数未知），把计算谱缩放使第一峰等高。实质上 XAS 比较的是**指纹**（峰位模式与相对强度），而非绝对位置和绝对强度。
4. **展宽匹配**：初算用 `CH_SIGMA = 0.1` eV，峰比实验窄得多；把 `CH_SIGMA` 改为 0.5 重跑（复用已有 WAVECAR，不重做整个 SCF）。理论上很难定量再现实验展宽（谱仪的能量依赖展宽等因素未显式模拟），展宽线型（高斯、洛伦兹、常数、能量依赖）可自由选择——文献中各种都有。
5. 平移缩放后（FCH −58.347 eV、XCH −59.369 eV，强度分别乘 0.276/0.166、0.276/0.089），两条谱与实验形状吻合；相对峰位与强度还可通过更彻底的参数收敛进一步改善（收敛计算见 Unzog et al.）。

**复习问题**：
1. FCH 与 XCH 的区别？
2. 计算 XAS 谱与实验对比需要做什么？（能量平移 + 强度缩放 + 展宽匹配）

**本部分涉及的 VASP 功能与标签汇总**：
- 芯空穴框架：`ICORELEVEL`（1=初态近似/无激子，2=终态近似）、`CLNT`（芯空穴化学种）、`CLN`/`CLL`（n、l 量子数，K 边 = 1s）、`CLZ`（激发电子数；XCH=1）。
- FCH 技巧：`NELECT` 减一 + VASP 自动背景电荷补偿；`--dry-run` 查 NELECT。
- 谱输出：`CH_LSPEC`（角动量分解）、`CH_SIGMA`（展宽）、`CH_NEDOS`（网格）、`CH_AMPLIFICATION`（按超胞原子数放大抵消体积归一）。
- 超胞规范：芯激发原子单独成种（POTCAR 三份）；超胞收敛以压制芯空穴周期镜像相互作用；`NBANDS` ≈ 2×占据带数；`vasp_gam` + Γ 点。
- 分析：py4vasp `dielectric_function.plot("Im(XAS)")`、`partial_density.to_ngl`（`LPARD` + `IBAND` 部分电荷）；谱的平移/缩放/展宽对比规范。

## 12.2 XAS（二）：BSE 方法计算 LiCl K 边 XAS

> 来源：<https://vasp.at/tutorials/latest/xas/part2/> ｜ 练习文件：[XAS-part2.zip](https://vasp.at/tutorials/latest/XAS-part2.zip) ｜ 示例目录：`$TUTORIALS/XAS/e03_LiCl_G0W0`、`e04_LiCl_BSE`

本页用多体微扰理论（MBPT）模拟岩盐 LiCl 的 X 射线吸收：例 3 先做 $G_0W_0$ 得到准粒子能量与 RPA 介电张量，例 4 用 BSE 计算含激子效应的 K 边 XAS、对比 RPA 层次、并用 `NBSEEIG`/`BSEHOLE` 计算与可视化激子波函数、区分暗/亮激子。

### 示例 3：LiCl 的 $G_0W_0$ 带隙

**教学目标**：(1) 运行单发 $G_0W_0$；(2) 获得 $G_0W_0$ 带隙；(3) 计算 RPA 介电函数。

**任务与流程**：与 §10 相同的三步——DFT 基态（PBE，6×6×6 Γ 心网格，`ENCUT = 300`，`EDIFF = 1E-8`）→ 精确对角化未占据带（`NBANDS = 64`，约默认值的 3–4 倍，`ALGO = Exact` + `NELM = 1` + `LOPTICS = T` 写 WAVEDER，副产品是 IPA 介电函数）→ $G_0W_0$。POTCAR 用 `_GW` 势（Li_GW 等）。

**输入要点**：`gw/INCAR` 中 `ALGO = EVGW0`、`NELMGW = 1`（单发）、`NBANDS = 64` 与上步一致、`NBANDSGW = 32`（为 32 条带计算 QP 修正）、`NOMEGA = 50`。

**结果**：DFT 基本带隙 6.4 eV（OUTCAR `fundamental gap:` 行）；$G_0W_0$ 后 OUTCAR "QP shifts" 段给出 **QP 基本带隙 8.9 eV**（Γ 点：VBM QP −2.4378 eV，CBM QP 6.4198 eV）。对比 IPA 与 RPA 介电函数（py4vasp `dielectric_function.plot("Re(IPA, RPA)")`）：实验介电常数 2.7，考察哪个近似更接近。GW 与 DFT 一样需要彻底收敛分析，且额外要对空带数 `NBANDS` 与频率点数 `NOMEGA` 收敛——本教程只做到较低收敛水平。

**复习问题**：
1. $GW$ 近似的是什么量？（自能）
2. VASP 把准粒子能量写在哪里？（OUTCAR）
3. VASP 默认在什么级别计算 GW 的屏蔽库仑作用？（RPA）

### 示例 4：BSE 方法计算 Li K 边 XAS

**教学目标**：(1) 模拟含/不含激子效应的 XAS；(2) 计算并可视化激子电荷密度；(3) 区分暗激子与亮激子。

**原理**：BSE 写成双粒子极化率 $L$ 的 Dyson 方程 $L = L_0 + L_0(v-W)L$，TDA 下化为 $A X_\lambda = \omega_\lambda X_\lambda$，双粒子哈密顿量

$$A_{vc}^{v'c'} = (\varepsilon_v - \varepsilon_c)\delta_{vv'}\delta_{cc'} + \langle cv'|V|vc'\rangle - \langle cv'|W|c'v\rangle$$

**XAS 与光学区 BSE 的关键区别**：光学区的占据态 $v,v'$ 只含价带，而 XAS 中 $v,v'$ 必须包含被激发的**芯态**。能量与轨道取自 GW 准粒子。屏蔽势取静态近似 $W(\mathbf{q}, \omega = 0)$。

**输入要点**（目录 `e04_LiCl_BSE/bse`，示例文件内容略）：
- 需要的文件：GW 步的 **WAVECAR**（QP 能量与轨道）、**WFULLxxxx.tmp**（RPA 介电张量）、WAVEDER（仅当包含价带轨道时需要）；POSCAR/KPOINTS/POTCAR 与 GW 步一致。
- INCAR：`ALGO = BSE`、`ANTIRES = 0`（TDA）、`NBANDS = 64`（与 GW 一致）、`NBANDSV = 16`（BSE 导带数）、**`NBANDSO = 0`**（XAS 不含价带——与光学区 BSE 的最大差异）、芯态标签与 ΔSCF 相同：`ICORELEVEL = 2`、`CLNT = 1`、`CLN = 1`、`CLL = 0`（1s）、`CSHIFT = 0.25`。

**流程与结果**：
1. 运行后检查 OUTCAR 的 BSE 哈密顿量秩 `rank = 3456`，应等于 $N_c \times N_v \times N_k$（芯态数 × 导带数 × 全布里渊区 k 点数）。
2. vasprun.xml 的 opticaltransitions：**第一个跃迁（56.93 eV）是暗激子**（振子强度 0），其后的 57.368 eV 三重简并态是**亮激子**（振子强度 1.35）。振子强度在偶极近似下按选择定则计算：违反选择定则的跃迁无光学活性（暗），遵守的有活性（亮）。
3. 与实验（Shirley 2004）对比：谱整体平移 −55.3 eV 后主要特征与相对位置吻合良好。**绝对位置再现不佳的原因**：PAW 中芯态被冻结、直接取自 POTCAR；GW 的 QP 修正只施加于价带与导带，芯态固定在 DFT 级别。
4. **收敛提醒**：6×6×6 网格远未收敛（需密得多的采样，见 Unzog et al. 2022 与 Olovsson et al. 2009）；介电函数还需对导带数彻底收敛——XAS 激子强局域，可能需要大量导带。**技巧**：若虚部全为零，检查 stdout 中 "the WAVECAR file was read successfully"。ΔSCF 与 BSE 的对比见 Unzog et al.：两者吻合良好，BSE 谱总体上更接近实验。
5. **RPA 对照**（`rpa/` 子目录）：`LADDER = .FALSE.`（关梯形图）+ `LHARTREE = .TRUE.` 即 RPA。对比 RPA 与 BSE 谱可见激子效应对 XAS 的决定性影响。一个值得注意的现象：BSE 算法中得到的"RPA 介电函数"实部接近 1，与 GW 步的 RPA 结果差异巨大——因为 BSE 中只考虑了芯电子，1s 芯态对该体系的屏蔽贡献极小；另一概念性差别是 BSE 的 RPA 在 TDA 下计算而 GW 步不是（本例数值上影响不大）。

**复习问题**：
1. RPA 与 BSE 介电函数的区别？
2. 如何估计激子束缚能？
3. 为什么 BSE 算法给出的 RPA 介电常数远小于 $G_0W_0$ 步的值？

**激子波函数分析**（`wf/` 子目录）：激子波函数是双粒子坐标函数 $\psi_\lambda(\mathbf{r}_e, \mathbf{r}_h) = \sum_{vc} X_{vc}^\lambda \phi_v^*(\mathbf{r}_h)\phi_c(\mathbf{r}_e)$；要在 3D 中可视化须固定其中一个坐标——XAS 中芯空穴强局域，通常**固定空穴坐标**。输入：`NBSEEIG = 3`（分析的本征矢数）、`BSEHOLE = 0.0 0.0 0.0`（把空穴固定在激发原子坐标）、`LCHARGH5 = .TRUE.`（激子电荷密度写入 vaspout.h5）、`PRECFOCK = low`（较小 FFT 网格以减小文件尺寸）。结果用 py4vasp 可视化：首个亮激子三重简并，须把三个分量求和（`bse.exciton.density.to_ngl("2+3+4", isolevel=1, center=True)`）——电荷密度很好地局域在两个最近邻 Li 原子之间、以 Li 原子为中心；首个暗激子可用 `"1"` 单独画出。

**复习问题**：
1. 激子局域在哪个原子？
2. 首个亮激子的局域程度如何？

**本部分涉及的 VASP 功能与标签汇总**：
- XAS-BSE 配方：`ALGO = BSE` + `NBANDSO = 0`（只含芯初态与导带）+ `NBANDSV`（导带数）+ `ICORELEVEL = 2`/`CLNT`/`CLN`/`CLL`（选芯态，与 ΔSCF 通用）；起点文件 WAVECAR + WFULLxxxx.tmp（+WAVEDER）。
- 近似开关：`LADDER`/`LHARTREE`（RPA vs BSE）、`ANTIRES`（TDA）。
- 激子分析：`NBSEEIG`（写本征矢）、`BSEHOLE`（固定空穴位置）、`LCHARGH5`、`PRECFOCK`；py4vasp `exciton.density.to_ngl`（分量求和语法 "2+3+4"）；暗/亮激子由 opticaltransitions 振子强度判定。
- 物理要点：芯态在 PAW 中冻结 → 绝对能量须平移；XAS 激子强局域 → 导带数与 k 网格收敛要求高。

## 13. 强关联体系（Strongly Correlated Systems）

强关联材料中电子–电子相互作用占主导，标准 DFT 的独立粒子近似不足以描述（通常含部分填充的局域轨道），可出现金属–绝缘体转变、磁性与非常规超导。教程覆盖 DFT+U、约束随机相位近似（cRPA）与动力学平均场理论（DMFT）等扩展。

## 13.1 强关联（一）：约束随机相位近似（cRPA）求 Hubbard 参数（NiO）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part1/> ｜ 练习文件：[strongly_correlated-part1.zip](https://vasp.at/tutorials/latest/strongly_correlated-part1.zip) ｜ 示例目录：`$TUTORIALS/strongly_correlated/e01_bands`、`e02_wannier_projection`、`e03_CRPA`

**背景**：强关联材料中电子-电子相互作用占主导，标准 DFT 的独立粒子近似无法充分描述（典型为局域的部分占据 $d$/$f$ 轨道元素），关联效应导致金属-绝缘体转变、磁性、非常规超导等现象。VASP 的扩展方法包括 cRPA（本部分）、DFT+U 与 DFT+DMFT（Part 2）。cRPA 用于计算模型哈密顿量的有效相互作用参数 $U$、$J$、$U'$（Phys. Rev. B 77, 085122 (2008)；Phys. Rev. B 112, 245102 (2025)）：核心思想是在 GW 的屏蔽库仑作用 $W$ 中**剔除目标态自身的屏蔽贡献**，得到的部分屏蔽库仑作用在张成目标空间的局域基下评估；目标空间通常低维（$\leq 5$ 态），便于用 DMFT 等高级理论处理。

cRPA 工作流五步：① DFT 基态（例 1）；② 未占据 KS 态（例 1）；③ DOS 投影与 fatband（例 1）；④ Wannier 投影（例 2）；⑤ cRPA 计算（例 3）。体系为原胞 NiO，12×12×12 k 网格。

### 示例 1：NiO 能带结构（基态、未占据态、投影 DOS 与 fatband）

**教学目标**：(1) 计算基态波函数；(2) 生成 cRPA 所需的未占据态；(3) 投影 DOS（含 fatband）并确定目标空间。

**① DFT 基态**（目录 `e01_bands`）：POSCAR 为岩盐 NiO 原胞（$a = 4.171$ Å）；INCAR 用 PBE，`ISMEAR = -1`（Fermi 展宽——后续 `LFINITE_TEMPERATURE` 必需）、`SIGMA = 0.1`、`EDIFF = 1E-8`、`NBANDS = 48`、`ALGO = Normal`；KPOINTS 12×12×12；POTCAR 用 `_GW` 势（Ni_GW、O_GW_new——cRPA 需要大量未占据态）。运行后 py4vasp 给出 **PBE 带隙仅 0.05 eV**——看似金属，但 NiO 因强关联实际是 **Mott 绝缘体**（这正是需要 DFT+U/DMFT 的动机）。

**② 未占据态**（`unoccupied/`）：`NBANDS = 96`（96 通常不够，但教程够用）、`ALGO = Exact`（精确对角化保证未占据轨道与能量质量）、`LOPTICS = T`（频率依赖介电矩阵）、`LFINITE_TEMPERATURE = T`（计算全部光学矩阵元，含未占据态——cRPA 后续要用）。产物 WAVECAR+WAVEDER 是 cRPA 的起点。

**③ 投影 DOS 与 fatband**（`fatbands/`）：`ICHARG = 11`（冻结电荷密度插值能带）、`LORBIT = 11`（d/p 投影）、`LWAVE = F`、`EMIN = -5`/`EMAX = 15`/`NEDOS = 2000`。复制 CHGCAR+WAVECAR 运行后，用 py4vasp 叠画 d/p/s 特征能带与 DOS：费米能附近有 **8 条带（第 2–9 带，约 −8 至 2 eV）**承载几乎全部的 d 与 p 特征。

**目标空间的选择逻辑**（本例核心概念）：这 8 条 Bloch 带混合了 s/p/d 特征（可查 PROCAR）。只选 Ni 的 5 条 d 态理论上可行，但整体局域性反而变差；纳入 O 的 p 态（离域到 O 原子上）允许混合，改善费米能附近 Wannier 轨道的描述。s 态不宜纳入：它们与更高的 p 态纠缠（约 9 eV 处的交叉），难以单独选取，且对费米面附近 5 条带的贡献可忽略。DOS 佐证：约 −5 eV 的 p 主导态含 d 贡献、费米能附近 d 主导态含 p 贡献——p 必须进基组；约 7 eV 处几乎无 d 贡献、0 eV 附近几乎无 s 贡献——s 可排除。

**复习问题**：零带隙或纠缠体系必须用什么展宽？（Fermi 展宽 `ISMEAR = -1` + `LFINITE_TEMPERATURE = T`）等。

### 示例 2：Wannier 投影（目标空间基组）

**教学目标**：用 Wannier90 投影 Wannier 轨道并作图可视化其成分。

**原理**：Wannier 函数由 Bloch 态经变换矩阵 $T$ 构造：

$$|w_{i{\bf R}}^\sigma\rangle = \frac{1}{N_k}\sum_{n{\bf k}} e^{-i{\bf kR}} T_{in}^{\sigma({\bf k})} |\psi_{n{\bf k}}^\sigma\rangle$$

（教程用 $T$ 记号以避免与 DFT+U 的 $U$ 混淆）。基组通常足够局域，使得周期镜像间相互作用可忽略，于是只需 ${\bf R}=0$ 原胞内的 Wannier 函数。模型哈密顿量（如 cRPA）只用 Bloch 函数的小子集——这些**目标态**通常位于化学势附近、强局域于离子周围；只有目标态能被 Wannier 基良好表示时，模型哈密顿量才能成功求解。

**输入要点**（目录 `e02_wannier_projection`，示例文件内容略）：
- `LWRITE_WANPROJ = .TRUE.`：写出 WANPROJ 文件（每个 k 点与带的 Wannier 投影矩阵 $T$），后续步骤可跳过重复投影；
- `NUM_WANN = 8`；`LWANNIER90_RUN = T`（VASP 以库模式运行 wannier90）；
- `WANNIER90_WIN` 内容：`num_bands = num_wann = 8`、`exclude_bands = 1, 10-96`（只保留第 2–9 带）、投影 `Ni:d` + `O:p`、`bands_plot = true` 沿 L–Γ–X–W–L–K–Γ 路径插值能带验证投影质量、`bands_plot_project : 1-5`（fatband 显示各带 d 特征）；
- `ALGO = None` + `LWAVE = F`：后处理步，绝不改动波函数。

**流程与结果**：复制 `unoccupied/WAVECAR` 运行，得到 8 个 Wannier 轨道（写入 WANPROJ）。gnuplot 出 fatband 图：

![Wannier 插值能带与 d/p 轨道特征（Ni d 红、O p 蓝）](/imgs/2026-08-05/fig60.png)

若 Wannier 态选得好，Wannier 基的能带应与 Bloch 基能带一致。图中费米能附近的 5 条 Wannier 带以 d 特征为主、含重要 p 贡献；约 −7 eV 的 3 条带以 p 特征为主、含显著 d 贡献。`wannier90.wout` 的 Final State 给出各 WF 中心与展宽：前 5 个中心在 Ni 原子 (0,0,0)、展宽约 0.45 Å²，后 3 个中心在 O 原子、展宽约 0.87 Å²——**cRPA 目标态应选展宽最小的 Ni d 轨道**（Phys. Rev. B 56, 12847 (1997) 式 14）。

**可选对照（仅 d）**：`d-only/` 目录把 Wannier 数降为 5、排除 3 条 p 带、去掉 `O:p` 投影。后续 cRPA 用仅 d 基组时，裸 Hubbard U 从 26.37 降到 22.85、屏蔽 U 从 7.08 降到 6.12（4×4×4 网格）——**大量相互作用被遗漏**，印证纳入 p 态的必要性。

### 示例 3：cRPA 计算 Hubbard $U$、$U'$、$J$

**教学目标**：(1) 用粗/细 k 网格计算 cRPA 的 $U$、$U'$、$J$；(2) 计算非中心（off-center）项；(3) 用 Ohno 势拟合 Hubbard U 势并研究等离激元的频率依赖。

**任务**：用 8 个 Wannier 轨道对 NiO 做 cRPA，对比 k 网格（4×4×4 vs 6×6×6）与频率点数（12 vs 24）对库仑势的影响。Wannier 基中的有效相互作用矩阵为

$$U^{\sigma\sigma'}_{ijkl} = \int \mathrm{d}{\bf r}\int \mathrm{d}{\bf r}'\, w_i^{*\sigma}({\bf r}) w_j^{\sigma}({\bf r})\, U({\bf r},{\bf r}',\omega)\, w_k^{*\sigma'}({\bf r}') w_l^{\sigma'}({\bf r}')$$

实践中以 Hubbard–Kanamori 参数为指引（单胞 $R=0$，对应 `Vijkl`/`Uijkl` 数据集）：

$${\cal U} = \frac{1}{N}\sum_{i\in{\cal T}} U_{iiii}, \qquad {\cal U}' = \frac{1}{N(N-1)}\sum_{i\neq j} U_{ijji}, \qquad {\cal J} = \frac{1}{N(N-1)}\sum_{i\neq j} U_{ijij}$$

扩展体系（$R>0$）对应 `VRijkl`/`URijkl`。实际 DFT+U 计算（Anisimov 1991；Liechtenstein 1998）推荐球平均，Part 2 讨论。

**输入要点**（目录 `e03_CRPA`，示例文件内容略）：
- 主 INCAR 与 `fine_grid/INCAR`：`ALGO = CRPA`；`LOCALIZED_BASIS` 由工作目录中的 WANPROJ 文件提供（复用上一步）；`NTARGET_STATES = 1 2 3 4 5`——**从 cRPA 屏蔽中剔除这 5 个 Wannier 态**（即目标态）；`LSCRPA = T`；`LTWO_CENTER = T`（多中心项）；`NBANDS` 与 WAVECAR 一致；`LFINITE_TEMPERATURE` 开启有限温度形式（需 Fermi 展宽），对零带隙/纠缠体系必需（单频点计算可关以提速）。
- `plasmons/INCAR`：`ALGO = CRPAR`（欧几里得时空算法版）、`NOMEGA = 12`、`LFINITE_TEMPERATURE = T`（Matsubara 形式恒用）、`MAXMEM = 64000`。
- KPOINTS：粗网格 4×4×4；`fine_grid` 为 6×6×6——**必须与 DFT 步的 12×12×12 可公度**（8×8×8 不行）。

**流程与结果**：
1. **粗网格**：复制 WAVECAR/WAVEDER/WANPROJ 运行。OUTCAR 中 `bare Hubbard U = 26.37`、`screened Hubbard U = 7.08`。py4vasp `effective_coulomb` 给出全套：裸 U/u/J = 26.37/24.61/0.88，屏蔽 U/u/J = 7.08/5.47/0.81。
2. **细网格**：裸相互作用几乎不变（已收敛）；屏蔽值 U 7.08→7.47、u 5.47→5.86（未收敛到 0.1 eV，需更密网格），J 基本不变（已收敛）。
3. **非中心项与 Ohno 拟合**：`LTWO_CENTER = T` 使 vaspout.h5 的 `uijkl` 多出 Wannier 中心 $R$ 维度，存储 $U_{ijkl}(R)$ 描述库仑势在胞间的衰减。Hubbard–Kanamori $U$ 的空间衰减常用 Ohno 势描述 $U(R)/U(0) = \sqrt{\delta/(R+\delta)}$。py4vasp `effective_coulomb.plot("U", radius_max=10)` 自动拟合 $\delta$；远距离数据因 k 采样不足出现噪声，教程演示了过滤离群点（>0.74 eV 且 R>10 Å）后用 `scipy.optimize.curve_fit` 重拟合（最优 $\delta \approx 0.0975$）。**实践提醒**：教程为演示而剔除噪声点；实际应加密 k 网格与提高截断能（噪声源于 FFT 网格过疏）。
4. **等离激元频率依赖**（`plasmons/`）：`ALGO = CRPAR` 用欧几里得时空算法（虚时间轴计算，度规 $(-+++)\to(++++)$，允许库仑核的稀疏表示）——小体系比 Minkowski 算法慢且吃内存（4 进程 4×4×4 约 8 GB），但大体系显著更快；注意 `estimated memory requirement` 打印后有长时间静默属正常。结果存于 `uijkl`（NOMEGA 个虚频点）。实频依赖通过 **AAA 解析延拓**（py4vasp `interpolate.AAAConfig(clean_up_tol=...)`）获得：24 个虚频点比 12 个能分辨更多结构；50 eV 以上的波纹状结构对 `clean_up_tol` 敏感，提示 Froissart doublets（数值噪声）；增加虚频点无用（Minimax 等距选点已最优，Phys. Rev. B 101, 205145 (2020)），**加密 k 点**才能稳定高频结构。
5. **物理结论**：第一个等离激元峰在 5.5 eV（大致对应 d 与 p 中心的间距）对 k 网格与频率点数不敏感——可靠，且低频区正是模型哈密顿量物理的关键；主等离激元峰在约 50 eV（集体激发，GW 类近似的典型特征）难收敛、物理意义有限；高频处屏蔽势 $U$ 逐渐趋近裸势 $V$。做计算时应专注收敛有清晰物理意义的低频峰。

**复习问题**：
1. 5.5 eV 的第一个等离激元峰对应什么？
2. 50 eV 的主等离激元峰对应什么？
3. 高频下屏蔽库仑势如何变化？
4. 零带隙/纠缠体系的 cRPA 需要哪个标签？（`LFINITE_TEMPERATURE = T` + Fermi 展宽）

**本部分涉及的 VASP 功能与标签汇总**：
- cRPA 流程：PBE 基态（`ISMEAR = -1`）→ `ALGO = Exact` + `NBANDS = 96` + `LOPTICS = T` + `LFINITE_TEMPERATURE = T` → Wannier 投影（`LWRITE_WANPROJ`、`NUM_WANN`、`LWANNIER90_RUN`、`WANNIER90_WIN` 的 exclude_bands/projections）→ `ALGO = CRPA`（+WANPROJ、`NTARGET_STATES`、`LSCRPA`、`LTWO_CENTER`）。
- 等离激元：`ALGO = CRPAR`（欧几里得时空/虚时间轴）、`NOMEGA`、`MAXMEM`；AAA 解析延拓（py4vasp interpolate）。
- 输出：OUTCAR bare/screened Hubbard U；vaspout.h5 `/results/crpa`（uijkl/vijkl/cijkl）；py4vasp `effective_coulomb`（含 Ohno 拟合、频率依赖绘图）。
- 概念：Mott 绝缘体、目标空间选择（d+p vs 仅 d）、Wannier 展宽判据、k 网格可公度性、Hubbard–Kanamori 平均。

## 13.2 强关联（二）：DFT+U 与 DFT+DMFT（NiO 全流程）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part2/> ｜ 练习文件：[strongly_correlated-part2.zip](https://vasp.at/tutorials/latest/strongly_correlated-part2.zip) ｜ 示例目录：`$TUTORIALS/strongly_correlated/e04_scf`、`e05_dftu`、`e06_dmft_os`、`e07_dmft_csc`、`e08_dmft_os_qmc`

本页先跑 VASP（DFT 与 DFT+U），再用 `solid_dmft`（TRIQS）以 VASP 的投影局域轨道（PLO）做 DFT+DMFT。五个例子构成完整链条：例 4 非磁 DFT 与 DMFT 准备（PLO 转换 + 从 cRPA 提取 $U/J$）、例 5 反铁磁 DFT+U 与 DOS、例 6 单发 DFT+DMFT（Hartree 求解器）、例 7 电荷自洽（CSC）DFT+DMFT、例 8 反磁 vs 顺磁、Hartree vs CT-QMC 对比。

### 示例 4：DFT SCF 与 DMFT 准备（PLO 与 U/J 提取）

**教学目标**：(1) 设置并运行 AFM 单胞 NiO 的 DFT 计算；(2) 控制 DMFT 所需 Ni-$d$ 投影子的生成；(3) 了解 DMFT 相关输出文件；(4) 从 cRPA 提取 DMFT 所需的 $U$、$J$。

**DFT+U 原理**：DFT+U 是对标准 DFT 的低成本修正，针对局域轨道（典型为 $d$/$f$ 壳层）中的强关联电子。半局域泛函存在虚假自相互作用，使电子人为离域、低估 Mott 绝缘体带隙。DFT+U 通过在目标子空间的双重占据上加 Hubbard 型在位库仑惩罚来矫正。Dudarev 旋转不变简化形式下总能读作

$$E = E_{\mathrm{DFT}} + \sum_I \left[ \frac{U^I}{2} \sum_{m\sigma \neq m'\sigma'} n_m^{I\sigma} n_{m'}^{I\sigma'} - \frac{U^I}{2} n^I(n^I - 1) \right]$$

第二项是**双计数修正**——扣除 $E_{\mathrm{DFT}}$ 中已包含的相互作用估计。括号项对整数占据消失（原子极限下不改变 $E_{\mathrm{DFT}}$），但惩罚分数占据、驱使体系走向局域的整数化解。DFT+U 便宜且常定性正确，但 $U$ 须谨慎选取（cRPA 或线性响应），且无法捕捉动力学关联与 Kondo 物理——这正是需要 DFT+DMFT 之处。

**输入要点**（目录 `e04_scf`，示例文件内容略）：POSCAR 为 AFM 单胞 NiO（4 原子，沿 [111] 方向反铁磁排列的常规单胞）；INCAR 收紧收敛（`EDIFF = 1E-8`、增大 `NELMIN`/`NBANDS`），`ISMEAR = -5`（四面体法 + Blöchl 修正）、`KSPACING` 生成 k 网格、`NCORE`/`KPAR` 并行控制；关键投影标签 `LORBIT = 14` + `LOCPROJ = "1 2 : d : Pr"`（对两个 Ni 做 d 投影，`Pr` 即 projector 类型）。

**PLO 转换**：VASP 的局域投影需转为 TRIQS/solid_dmft 使用的投影局域轨道（PLO）。`plo.cfg` 的关键条目：`EWINDOW = -10 10`（寻找相关 KS 态的能量窗）、`BANDS = 10 28`（投影的 VASP 带指标窗）、`LSHELL = 2`（d 壳层）、`IONS = 1 2`（两个 Ni 承载关联壳层）、`CORR = True`（标记为 DMFT 处理的关联壳层）。运行 `python converter.py` 生成 `vasp.h5` 与 PLO 步的投影 DOS（`pdos_*.dat`）。converter 日志给出壳层信息（l=2、2 个离子、维数 5、关联）与 5×5 占据矩阵 $n_m^{I\sigma}$（迹约 8.44，即杂质密度 16.88——两 Ni 合计）。

![VASP 局域投影 DOS 与 PLO 投影 DOS 对比（Ni d）](/imgs/2026-08-05/fig34.png)

对照 vaspout.h5 的 VASP Ni-d DOS 与 PLO 的 pdos：两者吻合良好，说明 PLO 充分捕获了计算所需的 Ni d 特征。（约 −5 eV 处高度差异源于 VASP 与 TRIQS 正交化过程的不同。）

**从 cRPA 提取 U/J**：Part 1 的部分屏蔽库仑作用存于 UIJKL 文件（四指标张量，$\omega = 0$）。教程提供更接近收敛的参考文件 `UIJKL_cRPA_ref`，用 TRIQS 工具 `fit_slater_fulld` 把四指标张量拟合到 Slater 参数化的相互作用哈密顿量：先用试探 Slater 参数（$F^0, F^2, F^4$）构造旋转不变库仑张量，约化为密度-密度（$U_{iijj}$）与交换（$U_{ijij}$）两指标矩阵，再最小二乘拟合。结果：**$U = 6.181$ eV、$J = 1.134$ eV**——后续所有 DFT+U 与 DFT+DMFT 计算一致使用这对参数。

**复习问题**：哪两个 VASP 文件含局域投影？TRIQS 的哪个文件含 PLO？什么是 PLO、如何用它定义 DFT+DMFT 的关联子空间？屏蔽库仑作用与裸库仑作用的区别？

### 示例 5：反铁磁 DFT+U 与 DOS

**教学目标**：用 cRPA 导出的相互作用参数运行 AFM DFT+U，并理解带隙打开机制。

**输入要点**（目录 `e05_dftu`，示例文件内容略）：核心 LDAU 标签——`LDAU = .TRUE.`、`LDAUTYPE = 1`（Liechtenstein 旋转不变形式）、`LDAUL = 2 -1`（Ni 加 d 壳层 U，O 不加）、`LDAUU = 6.181 0.00`、`LDAUJ = 1.134 0.00`（即 cRPA 拟合值）、`LDAUPRINT = 1`；磁性：`ISPIN = 2`、`MAGMOM = 2.0 -2.0 2*0`（两 Ni 反平行，O 为零）；另需 `LMAXMIX = 4`、`LASPH = .TRUE.`、`NEDOS = 5001`。POTCAR 用 Ni_pv 与 O 常规势。

**结果与物理**：

![DFT 与 DFT+U 的 NiO DOS 对比](/imgs/2026-08-05/fig35.png)

DFT（灰虚线）无带隙；DFT+U（紫）打开略大于 3 eV 的绝缘带隙（VBM/CBM 绿线标出），导带约 2.5 eV 出现强峰、价带整体下移。机理：NiO 是典型 Mott 绝缘体，Ni$^{2+}$ 为 $d^8$；八面体晶场把 Ni d 分裂为低能 $t_{2g}$ 三重态与高能 $e_g$ 二重态，8 个 d 电子使 $t_{2g}$ 全满（6 电子），其余 2 个占据 $e_g$。允许磁有序后交换分裂驱动高自旋 $S=1$ 组态：多数自旋 $e_g$ 全占、少数自旋 $e_g$ 空——DFT+U 把占据的多数自旋 $e_g$ 推向低能、空的少数自旋 $e_g$ 推向高能，在两者间打开带隙，DOS 在 2.5 eV 出现强 $e_g$ 峰；全满的 $t_{2g}$ 基本不受影响。带隙具**电荷转移**特征（O 2p → Ni $e_g$）而非纯 Mott–Hubbard 型，但 DFT+U 恢复了接近实验（约 4 eV）的带隙。

**复习问题**：哪些 INCAR 标签定义 DFT+U 的 $U$、$J$？为什么必须做自旋极化计算——只加 U 为何不开带隙？

### 示例 6：单发 DFT+DMFT（Hartree 求解器）

**教学目标**：(1) 用静态 Hartree 杂质求解器经 `solid_dmft` 设置并运行单发 DFT+DMFT；(2) 从 `observables_*.dat` 与 h5 存档监控并解读 DMFT 收敛诊断；(3) 后处理得到关联谱函数并与 DFT+U DOS 对比。

**原理**：DFT 收敛能量与密度，DFT+DMFT 收敛 **Green 函数与自能**。单粒子 Green 函数 $G(\mathbf{k},\omega)$ 编码向相互作用体系添加/移除电子的概率幅，是 DMFT 的自然语言：谱函数 $A(\mathbf{k},\omega) = -\frac{1}{\pi}\mathrm{Im}\,G(\mathbf{k},\omega+i0^+)$ 直接对应 ARPES 强度；k 积分谱函数在无相互作用 0 K 极限退化为 DOS，一般情况下还携带准粒子寿命与 Hubbard 卫星信息。非相互作用 Green 函数由 PLO 与 KS 本征值构造，k 积分得到 PLO 基下的局域 Green 函数 $G_{\mathrm{loc}}$；与局域相互作用哈密顿量一起定义 DMFT 的**杂质问题**，求解产生新自能并进入自洽循环。相互作用经 Dyson 方程 $G^{-1} = G_0^{-1} - \Sigma$ 进入：$\mathrm{Re}\,\Sigma$ 移动与重正化准粒子能量，$\mathrm{Im}\,\Sigma$ 给出有限寿命——两者都是 DFT+U 看不见而 DMFT 自然捕捉的。DFT+U 可视为 Hubbard 相互作用的静态平均场（自能无频率依赖）；DMFT 保留 $\Sigma(\omega)$ 的完整频率依赖，捕捉准粒子重正化、Hubbard 卫星、Kondo 效应。形式上 DMFT 把晶格问题映射到单站点 Anderson 杂质模型——关联杂质耦合到自洽确定的无相互作用浴（Rev. Mod. Phys. 78, 865 (2006)）。本例用 triqs/hartree_fock 静态求解器：对密度-密度相互作用应复现 DFT+U 类静态平均场物理，无需额外参数。

**流程**：复制 `e04_scf/vasp.h5`，`mpirun -np 4 solid_dmft`（4 核 6–10 分钟）。`dmft_config.toml` 关键设置：`seedname = "vasp"`、`jobname = "beta20-hartree-afm"`、`beta = 20`、`n_iw = 401`、`h_int_type = "density_density"`、`h_int_basis = "vasp"`、`U = 6.181`/`J = 1.134`（与 cRPA 提取值一致——必须核对）、`magnetic = true`、`magmom = [1.5, -1.5]`、`afm_order = true`、`n_iter_dmft = 5`（资源受限只做 5 次）、双计数 `dc_type = 0` + `dc = true` + `dc_dmft = true`、`[solver] type = "hartree"`、`[advanced] dc_fixed_occ = 8.0`。

**收敛监控**：`observables_imp0_up.dat` 每迭代一行——化学势 `mu`（0.066 → 0.987 eV 后稳定）、每轨道占据、杂质占据（up 通道 4.22 → 4.95）。收敛后求解器前后的杂质密度矩阵应一致。四个收敛判据：$\mu$、杂质占据、Weiss 场变化 $\delta\mathcal{G}^0(i\omega_n)$、DMFT 自洽条件 $\delta G_{\mathrm{imp}} = \|G^{\mathrm{loc}} - G^{\mathrm{imp}}\|$（后两者也存于 `conv_imp0.dat`）：

![DMFT 单发计算收敛诊断（化学势、杂质占据、Weiss 场与自洽条件）](/imgs/2026-08-05/fig36.png)

**后处理谱函数**：`calc_Aw_hartree.py` 从 h5 存档读取 Hartree 自能、$\mu$、双计数势与能，创建同网格 SumK 对象，计算 k 积分谱函数 $A(\omega) = -\frac{1}{\pi}\sum_k \mathrm{Im}\,G(k,\omega)$（含展宽 0.1 eV）。注意两条 UserWarning（Sigma 网格 −10 至 10 eV 未完全覆盖带能范围 −6.98 至 14.87 eV）：

![DFT+U DOS 与单发 DMFT（Hartree）谱函数对比](/imgs/2026-08-05/fig37.png)

DFT+U 与 DFT+DMFT@Hartree 结果高度相似——验证了两方法学上的等价性（静态极限下 DMFT@Hartree 还原 DFT+U），带隙均约 3 eV。

**复习问题**：谱函数可与什么实验对比？为何 DMFT 以 Green 函数为中心量？DMFT 如何检查收敛？

### 示例 7：电荷自洽（CSC）DFT+DMFT

**教学目标**：(1) 解释 CSC 条件及其对完全自洽 DFT+DMFT 的重要性；(2) 以收敛的单发结果为起点配置并运行 CSC DFT+DMFT；(3) 评估 CSC 对谱函数与带隙的影响。

**原理**：除 DMFT 自洽条件 $|G^{\mathrm{imp}} - G^{\mathrm{loc}}| \to 0$ 外，DFT+DMFT 还有第二个自洽条件。由晶格 Green 函数可算密度矩阵 $N_{\nu\nu'}(\mathbf{k}) = \frac{1}{\beta}\sum_{i\omega_n} G_{\nu\nu'}(\mathbf{k}, i\omega_n)$，它与 DFT 密度的关系为 $\rho(\mathbf{r}) = \sum_{\mathbf{k}}\sum_{\nu\nu'} \langle\mathbf{r}|\Psi_{\nu\mathbf{k}}\rangle N_{\nu\nu'}(\mathbf{k}) \langle\Psi_{\nu'\mathbf{k}}|\mathbf{r}\rangle$。因此 DFT 电荷密度与 DMFT 密度矩阵必须一致——即**电荷自洽条件**，要求每个 DMFT 步后做电荷密度更新（`ICHARG = 5`），DMFT 自能被反馈进 DFT 循环。

**输入要点**（目录 `e07_dmft_csc`，示例文件内容略）：VASP 侧 INCAR——CSC 标志 `ICHARG = 5`、`NELM = NELMIN = 1000`（设得极大，防止 VASP 提前终止——必须由 TRIQS 终止）、`NELMDL = -2`、密度混合 `IMIX = 4`（Broyden 2 + Pulay）、`AMIX = 0.02`、`BMIX = 0.2`、`LSYNCH5 = True`（运行中实时同步 vaspout.h5，供 py4vasp/solid_dmft 监控）、`LMAXMIX = 6`、投影 `LORBIT = 14` + `LOCPROJ = "1 2 : d : Pr"`、`EMIN = -6`/`EMAX = 18`（优化投影通道的能量窗）。dmft_config.toml 新增：从单发结果加载自能作为起点（`load_sigma = true`、`path_to_sigma = "../e06_dmft_os/..."`）；CSC 迭代结构 `n_iter_dmft_first = 2`（先跑 2 次 DMFT）、`n_iter_dmft_per = 2`（每次 CSC 更新后跑 2 次）、`n_iter_dmft = 10`（总计）；新增 `[dft]` 段：`plo_cfg`、`dft_code = "vasp"`、`projector_type = "plo"`、`n_iter = 4`（每次 CSC 更新的 DFT 迭代数）、`n_cores`、`mpi_exe`、`dft_exec`。

**流程与结果**：复制 vasp.h5、符号链接 e04 的 CHGCAR/WAVECAR 后运行（10–15 分钟）；`out` 文件中 VASP 与 solid_dmft 交替输出，直观展示两者如何循环；出现 `solid_dmft: Stopping VASP` 即完成（VASP 不停可 `pkill -9 vasp_std`）。再跑一次后处理脚本：

![单发与 CSC DMFT 谱函数对比](/imgs/2026-08-05/fig38.png)

**CSC 使带隙相对单发 DMFT 略微缩小**——符合预期：单发 DMFT 通常高估关联效应，CSC 会削减之。相对 DFT+U 带隙也略小，仍约 3 eV，低于实验的 **4.3 eV**（Phys. Rev. Lett. 53, 2339 (1984)）。

**复习问题**：CSC 全称是什么？为什么 DFT+DMFT 需要更新 DFT 电荷密度？CSC 相对 OS DMFT 的典型影响？

### 示例 8：AFM vs PM、Hartree vs QMC 对比

**教学目标**：(1) 在 Matsubara 轴与实频轴读取并解读 QMC 自能、识别强关联特征；(2) 区分 AFM 对称破缺带隙与顺磁动力学 Mott 带隙的不同物理起源；(3) 批判性对比 DFT+U、单发 Hartree DMFT 与顺磁 QMC 的谱函数。

**任务**：对顺磁 NiO 做连续时间量子蒙特卡洛（CT-QMC），与 AFM DFT+U/DMFT 对比。核心区别：Hartree 计算经 AFM 长程序破缺自旋对称（两 Ni 带相反静态磁矩）；顺磁（PM）解在平均上强制自旋对称（$\langle n_\uparrow\rangle = \langle n_\downarrow\rangle$），但保留静态平均场看不见的**动力学局域磁矩涨落**。CT-HYB 通过 $\Sigma(\omega)$ 的完整频率依赖捕捉这些涨落：即使无长程序，强在位库仑排斥 $U$ 压制双重占据，费米能级谱权重转移到上下 Hubbard 带，打开 **Mott 带隙**。这正是 Néel 温度以上（$T_N \approx 520$ K）的关键物理——NiO 在 PM 相仍是绝缘体，不是因为磁有序，而是关联驱动的局域化本身足以打开谱隙。预计算结果在 `e08_dmft_os_qmc/UcRPA-beta20-qmc1e+7/vasp.h5`（自行运行需数小时）。

**求解器要点**：`beta = 20` 对应 $T \approx 580$ K（$\beta = 1/k_BT$）——QMC 在高温更便宜，因为虚时间区间 $[0,\beta]$ 收缩。精度由 MC 移动总数与测量次数决定：`[solver]` 段 `n_warmup_cycles = 1e+4`（每个 MPI walker 热化，简单高温问题 5000–10000 足够）、`n_cycles_tot = 1e+7`（测量总数）、`length_cycle = 60`。Hartree 是静态自能（捕捉对称破缺与静态能移）；CT-QMC 是数值精确（统计误差内）的杂质求解器（捕捉 Hubbard 带、寿命、频率依赖自能）。

**自能分析**：CT-QMC 在 Matsubara 轴给出 $\Sigma(i\omega_n)$，两个特征最有信息量：低频行为（Fermi 液体中 $-\mathrm{Im}\,\Sigma(i\omega_n) \propto \omega_n \to 0$；大值或发散提示强散射或 Mott 倾向）与轨道分化（立方晶场分裂 Ni d 为 $t_{2g}$（三重简并，块指标 `up_0`）与 $e_g$（二重简并，`up_2`）；$e_g$ 与氧杂化更弱，预期关联更强）。

![Matsubara 轴自能：t2g 与 eg 对比](/imgs/2026-08-05/fig39.png)

实频自能需**解析延拓**（MaxEnt 最大熵法，对含噪 QMC 数据的标准做法）：solid_dmft 后处理脚本 `sigma_maxent.py`（`omega_min/max = ±30`、`maxent_error = 0.01`、`LineFitAnalyzer`、`continuator_type = "inversion_dc"`）把 `Sigma_freq_*` 延拓为 `Sigma_maxent_*`。实频自能给出准粒子权重 $Z = (1 - \partial\,\mathrm{Re}\,\Sigma/\partial\omega|_{\omega=0})^{-1}$ 与散射率 $\Gamma \propto -\mathrm{Im}\,\Sigma(\omega=0)$。注意解析延拓是不适定逆问题，高频细节须谨慎解读。

![实频自能 Re Σ 与 −Im Σ（t2g 与 eg）](/imgs/2026-08-05/fig40.png)

**顺磁谱函数总对比**：

![DFT+U、单发 Hartree DMFT 与顺磁 QMC 谱函数对比](/imgs/2026-08-05/fig41.png)

三层理论都打开带隙，但物理原因与谱形定性不同：DFT+U 与 DMFT@Hartree 谱几乎相同（都是依赖 AFM 对称破缺的静态平均场——多数自旋 $e_g$ 压到费米能以下、少数自旋 $e_g$ 推上去，给出约 3 eV 的锐利带隙与窄准粒子峰；两者吻合验证了 DMFT@Hartree 实现在静态极限正确还原 DFT+U）。PM CT-HYB 定性不同：未施加任何长程磁序却清晰开隙——**顺磁 Mott 绝缘体**的标志，关联驱动局域化本身（编码在 $\Sigma(\omega)$ 的完整频率依赖中）足以在 Néel 温度以上打开谱隙（计算温度约 580 K > $T_N$）。谱特征更宽（$\mathrm{Im}\,\Sigma$ 给出有限准粒子寿命），且高束缚能（约 6–8 eV）出现额外非相干谱权重——**下 Hubbard 带**，静态 Hartree 谱中完全缺失，是动力学关联的直接指纹。

**复习问题**：
1. DFT+U 与 DMFT@Hartree 都在 NiO 打开带隙——各自的物理机制是什么？为何谱函数几乎相同？
2. PM CT-HYB 在 Néel 温度以上计算，NiO 仍是绝缘体——这说明磁有序对绝缘态的作用如何？
3. 哪个轨道流形（$t_{2g}$ vs $e_g$）关联更强（更大 $-\mathrm{Im}\,\Sigma(i\omega_0)$）？基于晶场分裂与 $d^8$ 计数解释。
4. 什么是解析延拓？为什么从 CT-QMC 得到实频谱函数必须做它？其局限是什么？
5. （进阶）PM 解平均上强制自旋对称，DMFT 如何不破缺对称仍捕捉局域磁矩形成？

**本部分涉及的 VASP 功能与标签汇总**：
- DFT+U：`LDAU`/`LDAUTYPE = 1`（Liechtenstein）/`LDAUL`/`LDAUU`/`LDAUJ`/`LDAUPRINT`；`ISPIN = 2` + `MAGMOM`（AFM）；`LMAXMIX = 4`、`LASPH`。
- DMFT 接口：`LORBIT = 14` + `LOCPROJ`（局域投影）；`plo.cfg` → `converter.py` → `vasp.h5`（TRIQS）；CSC 用 `ICHARG = 5`、`NELM/NELMIN` 极大、`IMIX/AMIX/BMIX`、`NELMDL = -2`、`LSYNCH5`、`EMIN/EMAX` 投影能量窗。
- 参数流：cRPA UIJKL → Slater 拟合（TRIQS `fit_slater_fulld`）→ $U/J$ 一致贯穿 LDAU 与 dmft_config.toml。
- solid_dmft：`dmft_config.toml`（general/solver/dft/advanced 四段）、Hartree vs CT-HYB 求解器、`load_sigma`、`n_iter_dmft_first/per`、收敛观测量（mu、杂质占据、$\delta\mathcal{G}^0$、$\|G^{\mathrm{imp}}-G^{\mathrm{loc}}\|$）；后处理 k 积分谱函数（SumkDFTTools）与 MaxEnt 解析延拓。
- 物理图景：双计数修正、$d^8$ 晶场分裂（$t_{2g}^6 e_g^2$）、交换分裂开隙、电荷转移 vs Mott–Hubbard 带隙、AFM 对称破缺隙 vs 顺磁 Mott 隙、Hubbard 卫星。

## 13.3 强关联（三）：PBE+U 起点的 $G_0W_0$ 与 BSE 光学性质（NiO）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part3/> ｜ 练习文件：[strongly_correlated-part3.zip](https://vasp.at/tutorials/latest/strongly_correlated-part3.zip) ｜ 示例目录：`$TUTORIALS/strongly_correlated/e09_G0W0`（`dft/`、`unoccupied/`、`g0w0/`）、`e10_BSE/`（`absorption/`、`eels/`、`spin_flip/`）

本部分把前两部分铺垫的强关联工具链推进到多体微扰理论（MBPT）：例 9 以第一部分 cRPA 确定的 $U$、$J$ 做 PBE+U 基态，再算 $G_0W_0$ 准粒子能带与屏蔽作用 $W$；例 10 以 $G_0W_0$ 产出的 $W$ 为输入求解 Bethe-Salpeter 方程（BSE），计算 AFM NiO 的光吸收谱、电子能量损失谱（EELS）与自旋翻转激发（磁振子）。整条链路展示了"强关联起点选择 - 准粒子修正 - 激子效应"的完整第一性原理光谱学工作流，也是第 10、11 章 GW/BSE 通用流程在强关联体系上的深化版。

### 示例 9：$G_0W_0$ 电子结构（PBE+U 起点）

**教学目标**：(1) 以 PBE+U 波函数为起点运行 $G_0W_0$ 计算；(2) 计算大量未占据 Kohn-Sham 轨道及其对 $\mathbf k$ 的导数；(3) 在独立粒子近似（IPA）下得到介电函数与光吸收；(4) 确定 NiO 的 $G_0W_0$ 带隙；(5) 计算随机相位近似（RPA）介电常数并与实验对比。

**原理：为什么起点很重要**

$GW$ 近似中自能 $\Sigma = iGW$ 所用的屏蔽作用 $W$ 由 RPA 极化率函数给出：

$$\chi(1,2) = -i\int d(3,4)\, G(1,3)\, G(4,1)\, \Gamma(3,4;2)$$

其中顶点函数取 RPA 值 $\Gamma(3,4;2) = \delta(3-2)\,\delta(4-2)$，$G$ 为 Green 函数。

以 **PBE 电子结构为起点**的 RPA 屏蔽依赖一种"部分误差抵消"：

- PBE 低估带隙，低能电子激发偏容易，导致**高估**屏蔽；
- RPA 缺顶点修正（vertex corrections），导致**低估**屏蔽。

对 $s$/$p$ 离域带主导的普通半导体，两种误差方向相反、大小相当，抵消效果良好（Phys. Rev. Lett. 99, 246403 (2007)）。但对**局域 $d$ 态主导**的体系，抵消失效：AFM NiO 的 PBE 带隙被严重低估，RPA 把屏蔽推得过高，介电常数算到约 15，而实验值约 5.7（Phys. Rev. Materials 2, 073803 (2018)）；Mott 绝缘体更极端，PBE 直接给出金属态，RPA 屏蔽完全失真。

修复途径有两条：

1. 自洽求解 Hedin 方程并在极化率中加入顶点修正（$QSG\tilde{W}$ 近似，Phys. Rev. B 108, 165104 (2023)）——最严格但计算极其昂贵，超出本教程范围；
2. **换一个更好的起点电子结构**——本教程采用 PBE+U，其中 Hubbard $U$ 由第一部分 cRPA 计算确定（$U = 6.181$ eV、$J = 1.134$ eV）。带隙改善后，RPA 极化率不再被过小的带隙污染，屏蔽与介电常数随之合理。

**计算流程总览**（目录 `e09_G0W0`，三个子目录串接，示例文件内容略）。一般可用标准 PBE 计算的波函数重启来改善 PBE+U 的收敛，本例为简洁省略该预备步：

1. `dft/`：AFM NiO 的 PBE+U 基态，产出 WAVECAR；
2. `unoccupied/`：大 `NBANDS` 精确对角化 + 轨道导数，更新 WAVECAR、产出 WAVEDER 与 IPA 介电函数；
3. `g0w0/`：$G_0W_0$ 准粒子能量与屏蔽作用 $W$，产出 `W*` 文件供后续 BSE 使用。

**第一步：PBE+U 基态（`dft/`）**

输入要点：

- POSCAR：AFM NiO 扩为 4 原子原胞（2 Ni + 2 O，晶格常数 4.171 Angstrom 的 FCC 常规单胞），以容纳反铁磁序——两 Ni 沿 [111] 方向磁矩反平行。
- INCAR：`ISPIN = 2` 开启自旋极化；`MAGMOM = 2.0 -2.0 2*0` 设定 AFM 初猜（两 Ni 反平行，O 为零）；`LASPH = T` 计入 PAW 球内梯度 corrections 的非球面贡献，磁性计算应视为标准设置；`LMAXMIX = 4`——对 $d$ 电子必须设置（$f$ 电子需设 6），否则写入 CHGCAR 的电荷密度不完整、无法用于自旋极化重启；同时电荷密度混合器只考虑到 `LMAXMIX` 的角动量分量，设足也加速收敛。DFT+U 标签：`LDAU = .TRUE.`、`LDAUTYPE = 2`（**注意本部分改用 Dudarev 形式**，与第二部分示例 5 的 Liechtenstein `LDAUTYPE = 1` 不同，两种形式在只给 $U$、$J$ 时数值上相近）、`LDAUL = 2 -1`（Ni 的 d 壳层加 U，O 不加）、`LDAUU = 6.181 0.00`、`LDAUJ = 1.134 0.00`（cRPA 拟合值，三部分教程一以贯之）；`EDIFF = 1E-8` 收紧电子收敛。
- KPOINTS：3x3x3 Gamma 中心网格——教学够用，正式计算必须做 k 点收敛分析。
- POTCAR：**GW 计算应使用带 `_GW` 后缀的势**（本例 `Ni_GW`、`O_GW_new`）——它们在更宽的能量范围内拟合原子散射性质，以支撑大量未占据态的计算；可检查 POTCAR 的 `TITEL` 字段确认。追求更高精度时还可能需要把半芯 $s$/$p$ 态纳入价带的 `_sv`/`_pv` 势。

运行与结果：`mpirun -np 4 vasp_std`。在 OUTCAR 中搜索 `fundamental gap`，或用 py4vasp（`calculation.bandgap.print()`）读取：价带顶 3.8266 eV、导带底 7.3925 eV，**基本带隙 3.566 eV**（VBM 位于 $(1/3, 1/3, 1/3)$、CBM 位于 $(-1/3, 1/3, 1/3)$，两自旋分量完全相同），直接带隙 3.862 eV，费米能 4.057 eV。对照三个层次：AFM NiO 实验基本带隙 4.3 eV（Phys. Rev. Lett. 53, 2339 (1984)）；纯 PBE 严重低估至约 0.95 eV（可自行关掉 U 验证，参见 J. Phys. Chem. A 121, 3318 (2017)）；加 Hubbard 修正后带隙改善到约 3.6 eV——仍低估，但已足以让后续 RPA 屏蔽合理化。本步产出的 WAVECAR 是生成 GW 所需虚轨道的基础。

**第二步：未占据 Kohn-Sham 轨道与轨道导数（`unoccupied/`）**

目的：为 $G_0W_0$ 提供大量空态与轨道对 $\mathbf k$ 的导数。POSCAR/POTCAR/KPOINTS 与第一步相同，WAVECAR 从 `dft/` 复制；INCAR 在第一步基础上新增三个标签：

- `ALGO = Exact`：精确对角化 KS 哈密顿量，是小原胞中计算大量空带的优选方式；
- `NBANDS = 256`：GW 准粒子能量的精度要求显著增加带数（空态参与极化率与自能的求和）；
- `LOPTICS = .TRUE.`：基态收敛后计算频率依赖介电矩阵，同时写出轨道导数文件 **WAVEDER**（GW/BSE 所需）。

运行后除确认带隙外，可用 py4vasp 的 `dielectric_function.plot()` 分析 IPA 介电函数：悬停读 $\omega = 0$ 处实部即静态介电常数——**相对实验值 5.7 被高估**，这正是带隙低估的后果。反之，若带隙已很好再现，IP 近似会因缺少电子-空穴相互作用（激子效应）而**低估**介电常数。这两条规律是判断介电常数偏差来源的实用判据。

**第三步：$G_0W_0$ 准粒子与屏蔽势（`g0w0/`）**

输入要点：从 `unoccupied/` 复制 WAVECAR 与 WAVEDER；POSCAR/POTCAR/KPOINTS 同前。INCAR 关键设置：

- `NBANDS` 与上一步一致（教程示例为 128，需与 WAVECAR 中实际参与 GW 的带数匹配）；
- **$G_0W_0$ 步不再加 Hubbard 修正**——INCAR 中没有任何 LDAU 标签：$U$ 的作用已体现在起点波函数里，准粒子修正由 GW 自能本身负责，重复加 U 会双重计数；
- `ALGO = EVGW0` + `NELM = 1`：单次 $G_0W_0$，不迭代更新本征值；
- `NOMEGA = 50`：极化率与自能的实频率网格点数；
- `PRECFOCK = low`、`ENCUTGW = 200`：为教学目的降低的取值，正式计算必须对二者以及 k 网格、NBANDS、NOMEGA 做系统收敛分析。

运行与结果：准粒子能量写入 OUTCAR——在 `spin component 1` 行之后的 `QP-energies` 列取值，得 $G_0W_0$ 带隙 **4.6 eV**，合理接近实验值 4.3 eV（本教学计算未对 k 网格、频率数、带数完全收敛，故与严格收敛值略有偏差）。用 py4vasp 画 RPA 介电函数（`dielectric_function.plot("RPA(Re,Im)")`）得 RPA 介电常数 **5.76**，略高于实验；若以更严格的 k 点/NBANDS/NOMEGA 收敛，得介电常数 5.2、带隙 4.25 eV——RPA 相对实验略低估，主要因为缺激子效应（下一例的 BSE 将补上）。

**结论**：PBE+U 起点同时给出好的电子结构（基本带隙）与合理的屏蔽作用 $W$，后者正是后续 BSE 计算所需的输入——这就是本例同时输出"准粒子能量 + 屏蔽势"两个成果的原因。

**复习问题**：
1. RPA 介电函数缺失哪些重要效应？
2. `LOPTICS = .TRUE.` 计算介电函数用的是什么近似？

### 示例 10：NiO 的 BSE 光学性质

**教学目标**：(1) 用 BSE 计算 NiO 的光吸收谱；(2) 计算电子能量损失谱（EELS）；(3) 找出 NiO 的自旋翻转激发。

**原理：Bethe-Salpeter 方程**

多体微扰理论中，电子-空穴相互作用（激子效应）通过 BSE 计入光谱。完整 BSE 是伪本征值问题：

$$\left(\begin{array}{cc} \mathbf{A} & \mathbf{B} \\ \mathbf{B}^* & \mathbf{A}^* \end{array}\right) \left(\begin{array}{l} \mathbf{X}_\lambda \\ \mathbf{Y}_\lambda \end{array}\right) = \omega_\lambda \left(\begin{array}{cc} -\mathbf{1} & \mathbf{0} \\ \mathbf{0} & \mathbf{1} \end{array}\right) \left(\begin{array}{l} \mathbf{X}_\lambda \\ \mathbf{Y}_\lambda \end{array}\right)$$

矩阵元描述占据态 $v, v'$ 与未占据态 $c, c'$ 之间经裸库仑作用 $V$ 与屏蔽势 $W$ 耦合的跃迁：

$$A_{vc}^{v'c'} = (\varepsilon_v^{QP} - \varepsilon_c^{QP})\,\delta_{vv'}\,\delta_{cc'} + \langle cv'|V|vc'\rangle - \langle cv'|W|c'v\rangle$$

$$B_{vc}^{v'c'} = \langle vv'|V|cc'\rangle - \langle vv'|W|c'c\rangle$$

方程的解 $\omega_\lambda$ 与 $(X_\lambda, Y_\lambda)$ 即**激子态**。Tamm-Dancoff 近似（TDA）忽略共振与反共振部分的耦合（令 $B = 0$），BSE 退化为标准本征值问题 $A(\mathbf q)|X_\lambda\rangle = E_\lambda |X_\lambda\rangle$；有限动量 $\mathbf q$ 下哈密顿量为

$$A_{\mathbf k c v}^{\mathbf k' c' v'}(\mathbf q) = \delta_{\mathbf k\mathbf k'}\,\delta_{vv'}\,\delta_{cc'}\left(\varepsilon_{\mathbf k+\mathbf q, v}^{QP} - \varepsilon_{\mathbf k c}^{QP}\right) - \left(f_{\mathbf k+\mathbf q, v} - f_{\mathbf k c}\right) K_{\mathbf k c,\, \mathbf k+\mathbf q\, v,\, \mathbf k'+\mathbf q\, v',\, \mathbf k' c'}$$

核 $K$ 含裸库仑势 $V$ 与静态屏蔽势 $W$。**共线自旋情形**下，在基 $|\uparrow\uparrow\rangle, |\uparrow\downarrow\rangle, |\downarrow\uparrow\rangle, |\downarrow\downarrow\rangle$ 中核的自旋结构为

$$K = \left[\begin{array}{cccc} V-W & 0 & 0 & V \\ 0 & 0 & -W & 0 \\ 0 & -W & 0 & 0 \\ V & 0 & 0 & V-W \end{array}\right]$$

因此 BSE 哈密顿量解耦为两类解：**含裸库仑势 $V$ 的自旋守恒跃迁**，与**不含 $V$ 的自旋翻转跃迁**。自旋守恒激发可以是光学允许的（亮激子）或被偶极选择定则禁戒的（暗激子）；共线情形下自旋翻转跃迁**总是暗的**。若计入自旋轨道耦合使自旋非共线，上述分解不再成立，必须求解完整核 $K$。本例在 TDA 下求解自旋守恒解（Phys. Rev. B 92, 045209 (2015)）。

#### 10.1 光吸收谱

输入要点（目录 `e10_BSE/absorption`，POSCAR/KPOINTS/POTCAR 沿用 e09，示例文件内容略）：

- 运行前从 `e09_G0W0/g0w0/` 复制全部 `W*` 文件（屏蔽作用等 GW 步产物），BSE 核直接在其上构建；
- `ALGO = BSE` + `IBSE = 2`：三种 BSE 求解算法中，**精确对角化**最精确，直接给出本征值与本征矢——后者是计算振子强度所必需的；代价是按 BSE 矩阵规模 $N^3$ 标度；
- `NBANDSO = 16` / `NBANDSV = 16`：纳入 BSE 核的价带/导带数目，决定激子 Hilbert 空间大小；
- `OMEGAMAX = 16`：输出谱的能量上限（eV）；
- `CSHIFT = 0.3`：Lorentz 展宽宽度（见下）；
- `NBANDS = 128`、`PRECFOCK = low`、`ENCUTGW = 200`、`ISPIN = 2` + `MAGMOM` 沿用 GW 步设置。

运行与结果：`mpirun -np 4 vasp_std` 后用 py4vasp 的 `bse.dielectric_function.plot("BSE(Im)")`，与实验谱及教程提供的收敛参考数据（`optics_converged.dat`）对比：计算谱的峰位与整体形状与实验合理吻合；谱在 16 eV 处截断正是 `OMEGAMAX` 的限制。收敛参考结果表明 3x3x3 k 网格已足以正确再现主峰位置与谱形。

**寿命展宽**：静态近似下 BSE 的解是实值的、寿命无限。真实激子因与其他激发耦合（声子贡献最主导）而具有有限寿命；本方法完全忽略这些过程，因此通过给本征值加一个小虚位移来模拟寿命展宽，得到 Lorentz 线型，宽度由 `CSHIFT` 设定（此处 0.3 eV，通常调节到与实验主要特征宽度匹配）。

**振子强度与暗/亮激子**：vasprun.xml 中 `<varray name="opticaltransitions">` 块（可用 `grep -A 16 '<varray name="opticaltransitions" >' vasprun.xml` 查看）第一列是激子态能量、第二列是振子强度。本例前 12 个跃迁强度为零——违反偶极选择定则的**暗激子**（成组出现在约 1.61 eV 与 2.65 eV）；第一个**亮激子**出现在约 3.93 eV（振子强度 18.57），正对应吸收谱中可见的起始峰。

**复习问题**：为什么有些跃迁的振子强度为零？共线自旋情形下自旋翻转激发可以是亮的吗？Lorentz 展宽如何设定？

#### 10.2 EELS（电子能量损失谱）

原理：EELS 中入射电子携带有限动量，可以激发**有限动量转移** $q$ 的态——光子只能竖直激发（$q$ 近似为 0）。$q$ 通过 `KPOINT_BSE` 标签指定，且一次计算只能给出该特定 $q$ 矢量的谱；对应 k 点的指标可在 OUTCAR 的 `irreducible k-points` 列表中找到。EELS 信号由 BSE 宏观介电函数取逆虚部得到：

$$\mathrm{EELS}(\mathbf q, \omega) = -\mathrm{Im}\left[\epsilon^{-1}(\mathbf q, \omega)\right]$$

输入要点（目录 `e10_BSE/eels`，示例文件内容略）：

- 实验上通常在小动量转移处测谱，计算可能需要稠密 k 采样才能包含对应 $q$ 点；本例从 e09 的 OUTCAR 不可约 k 点列表（共 13 个）选第二个——$(0.25, 0, 0)$（权重 6），故 `KPOINT_BSE = 2`；
- **`ANTIRES = 2`：光学区 TDA 足够好，但有限动量计算必须包含共振-反共振耦合**（Phys. Rev. B 91, 045209 (2015)）；
- 观察等离激元需要宽能量范围（`OMEGAMAX = 45`），BSE 矩阵显著增大，精确对角化（`IBSE = 2`）太贵，改用迭代算法：`IBSE = 1`（时间演化）或 `IBSE = 3`（Lanczos）。Lanczos 通常是计算介电函数最快的算法，但**当前 Lanczos 实现不支持 full BSE（ANTIRES = 2）**，所以本例用时间演化算法 `IBSE = 1`；
- 注意 `IBSE = 1` 只把介电函数写入 vasprun.xml，**没有** `opticaltransitions` 块；
- 其余：`NBANDSO/NBANDSV = 12`、`CSHIFT = 0.3`、`NBANDS = 128`。

运行与结果：复制 `W*` 文件后运行，用 py4vasp 的 `dielectric_function.read("BSE")` 提取介电函数，按 $\mathrm{Im}\,\varepsilon / (\mathrm{Im}\,\varepsilon^2 + \mathrm{Re}\,\varepsilon^2)$（即 $-\mathrm{Im}\,\varepsilon^{-1}$）画损失函数，与实验谱及收敛参考对比：计算谱与实验合理一致；最显著特征是损失函数在约 **20 eV** 处陡升——对应**纵向等离子体频率**，即价电子的集体电荷振荡（纵向等离激元）。

**复习问题**：Tamm-Dancoff 近似对有限 $q$ 计算是好近似吗？20 eV 峰对应什么类型的激发？

#### 10.3 自旋翻转激发

输入要点（目录 `e10_BSE/spin_flip`，示例文件内容略）：在 10.1 的 INCAR 基础上仅改一处——`LHARTREE = .FALSE.`，关掉核中的裸库仑（Hartree）项；按上面共线核的自旋结构，去掉 $V$ 后剩下的就是**自旋翻转**部分。其余标签（`IBSE = 2`、`NBANDSO/NBANDSV = 16` 等）与 10.1 相同。

运行与结果：复制 `W*` 文件后运行，检查 vasprun.xml 的 `opticaltransitions`：

- 第一个跃迁能量 **1.5514 eV**，比吸收计算（10.1）的对应跃迁**低约 50 meV**——这一移动正来自排斥性裸库仑贡献的缺失（`LHARTREE = .FALSE.`）；
- 所有自旋翻转跃迁的振子强度均为零：它们违反光学偶极跃迁的自旋守恒，**无光学活性**（暗激发）；
- 这些激发包含集体自旋激发即**磁振子（magnon）**，可在 BSE 框架内分析——这是把 BSE 从电荷激发推广到自旋激发的典型应用。

**复习问题**：自旋翻转跃迁有光学活性吗？共线自旋翻转跃迁使 BSE 哈密顿量中哪一项变为零？

**本部分涉及的 VASP 功能与标签汇总**：
- 起点策略：PBE+U（cRPA 的 $U/J$）作为 $G_0W_0$ 起点；RPA 屏蔽的误差抵消逻辑（带隙低估导致高估屏蔽 vs 缺顶点修正导致低估屏蔽）；GW 步**不加** LDAU 标签（避免双重计数）。
- DFT+U（本部分版本）：`LDAUTYPE = 2`（Dudarev）、`LDAUL/LDAUU/LDAUJ`；`ISPIN = 2` + `MAGMOM`（AFM）；`LMAXMIX = 4`（自旋极化重启必需）、`LASPH = T`。
- 空态与导数：`ALGO = Exact`、`NBANDS`、`LOPTICS = .TRUE.`（产出 WAVEDER）；IPA 介电函数诊断规律（带隙低估则介电常数高估；带隙准确则 IPA 因缺激子效应低估介电常数）。
- $G_0W_0$：`ALGO = EVGW0` + `NELM = 1`、`NOMEGA`、`PRECFOCK`、`ENCUTGW`；`_GW` POTCAR（必要时配 `_sv`/`_pv`）；OUTCAR 的 `QP-energies`；本例结果：QP 带隙 4.6 eV、RPA 介电常数 5.76（严格收敛后 5.2 与 4.25 eV）。
- BSE 求解：`ALGO = BSE`；`IBSE = 2`（精确对角化，$N^3$ 标度，可取振子强度）/ `IBSE = 1`（时间演化，支持 full BSE，不写 opticaltransitions）/ `IBSE = 3`（Lanczos，最快但不支持 `ANTIRES = 2`）；`NBANDSO`/`NBANDSV`（核的价/导带窗）；`CSHIFT`（Lorentz 寿命展宽）；`OMEGAMAX`（谱能量上限）。
- EELS：`KPOINT_BSE`（有限动量转移 $q$，指标取自 OUTCAR 不可约 k 点表）、损失函数 $-\mathrm{Im}\,\epsilon^{-1}$、`ANTIRES = 2`（有限 $q$ 必需）、约 20 eV 纵向等离激元。
- 自旋翻转：`LHARTREE = .FALSE.` 去除裸库仑项、共线核解耦为含 $V$（自旋守恒）与不含 $V$（自旋翻转）两类、暗激发与磁振子。
- 数据流：`dft/WAVECAR` - `unoccupied/WAVECAR + WAVEDER` - `g0w0/W*` - BSE 各目录复制 `W*`。

## 14. 核磁共振（NMR）

核磁共振（NMR）谱学是研究分子、液体与固体化学环境的高灵敏技术。第一性原理模拟可计算化学屏蔽、四极耦合常数与超精细常数，帮助理解 NMR 谱。

## 14.1 NMR（一）：化学屏蔽（收敛、诱导电流、与实验对比、预测）

> 来源：<https://vasp.at/tutorials/latest/nmr/part1/> ｜ 练习文件：[nmr-part1.zip](https://vasp.at/tutorials/latest/nmr-part1.zip) ｜ 示例目录：`$TUTORIALS/NMR/e01_shielding`、`e02_induced_current`、`e03_experiment`、`e04_comparison`

本页四个例子构成 NMR 方法论的完整闭环：例 1 在金刚石上学习化学屏蔽的计算与参数收敛；例 2 把屏蔽的微观来源——诱导电流——可视化（苯的环电流）；例 3 学会把计算屏蔽与实验化学位移正确地对比（多化合物线性拟合）；例 4 反过来用拟合预测未知化合物的实验位移并评估误差来源。

**NMR 物理背景**：具有非零核自旋 $\mathbf I$ 的原子核拥有磁矩 $\boldsymbol\mu$，在外加恒定磁场 $\mathbf B_{ext}$ 中绕场方向以 **Larmor 频率** $\omega_L$ 进动（进动频率由 $\mathbf B_{ext}$ 强度与核的旋磁比 $\gamma$ 决定）。若在垂直于恒定外场的方向施加射频脉冲（交变磁场），当 $\omega_{rf}$ 接近 $\omega_L$ 时发生**磁共振**：核磁矩从参考系翻转到横向系综，随后弛豫回参考系的过程产生信号，即 NMR 的测量对象。

![核磁矩绕外磁场进动示意](/imgs/2026-08-05/fig42.png)

化学环境相同的同位素核（$\gamma$ 相同）未必以同一频率振荡：核外电子壳层的差异会**屏蔽**原子核——例如溴乙烷 Br-CH2-CH3 中亚甲基 -CH2- 的 1H 在空间上比甲基 -CH3 更靠近 Br，感受到不同的电子屏蔽。屏蔽改变核处的有效磁场，从而改变共振频率。相对标准参考物（四甲基硅烷 TMS）的频率移动即**化学位移** $\delta$，它使 NMR 成为探测分子与固体化学结构的利器：

$$\delta_{ij}(\mathbf R) = \sigma^{ref}_{ij} - \sigma_{ij}(\mathbf R)$$

其中 $ij$ 遍历笛卡尔方向，$\mathbf R$ 为核位置。**化学屏蔽张量**通过线性响应计算：

$$\sigma_{ij}(\mathbf R) = -\frac{\partial B^{ind}_i(\mathbf R)}{\partial B^{ext}_j}$$

其中 $B^{ind}$ 是感应磁场。

### 示例 1：金刚石中化学屏蔽的收敛

**教学目标**：(1) 计算化学屏蔽；(2) 对平面波截断能与 k 网格收敛化学屏蔽；(3) 对比化学屏蔽与总能量收敛行为的差异。

**任务**：计算 2 原子金刚石原胞中 C 的各向同性化学屏蔽，并对截断能与 k 点收敛。

**输入要点**（目录 `e01_shielding`，示例文件内容略）：

- POSCAR：金刚石常规原胞（晶格常数 3.567 Angstrom，2 个 C）。
- INCAR（模板 `INCAR.nmr`）：普通 DFT 部分——`ENCUT = 400`（待扫描）、`ISMEAR = 0`/`SIGMA = 0.01`、`EDIFF = 1E-8`、`PREC = Accurate`、`LASPH = .TRUE.`（PAW 球内密度梯度的非球面贡献）。**NMR 要求电子态高精度收敛**：`PREC = Accurate`、`EDIFF` 不高于 1E-8、足够高的 ENCUT。NMR 专属标签：
  - `LCHIMAG = .TRUE.`：开启化学屏蔽的线性响应计算；
  - `LNMR_SYM_RED = .TRUE.`：丢弃与屏蔽线性响应中 k 空间导数计算方式不一致的对称操作（保持对称性处理自洽）；
  - `NLSPLINE = .TRUE.`：用样条插值保证倒空间 PAW 投影算符可微；
  - `KPAR = 4`：并行处理的 k 点数，加速计算。
- KPOINTS：4x4x4 Gamma 中心网格（待扫描）；POTCAR：标准 C 势。

**流程与结果**：

1. **截断能扫描**：`bash encut_shielding.sh` 把 ENCUT 从 400 eV 以 100 eV 步长扫到 900 eV，每次运行后用 `grep -A 5 "iso_shield" OUTCAR | tail -1 | awk '{print $6}'` 提取 13C 屏蔽存入 `encut_shielding.dat`。stdout 显示计算分两阶段：先常规 DFT 收敛 CHGCAR/WAVECAR，随后是沿磁场轴的线性响应——输出中 `BDIR` 是外加磁场方向、`IDIR` 是感应磁场方向（1:x、2:y、3:z），`IQ` 是 `ICHIBARE` 定义的有限差分模板网格点，取值从 `-ICHIBARE` 到 `+ICHIBARE`（默认 `ICHIBARE = 1`，故有 3 个 IQ 点）。
2. **OUTCAR 屏蔽输出解读**（`CSA tensor` 表，定义见 J. Mason, Solid State Nucl. Magn. Reson. 2, 285 (1993)）：

| 标签 | 含义 |
|---|---|
| `EXCLUDING G=0 CONTRIBUTION` | 不含宏观长程贡献的屏蔽，适用于分子 |
| `INCLUDING G=0 CONTRIBUTION` | 含宏观长程贡献（尤其表面电流）的屏蔽，适用于晶体 |
| `ion` | 按 POSCAR 排序的原子编号 |
| `BDIR` | 外加磁场（B0）方向（BDIR=1 即沿 x 轴）；每行再分 X/Y/Z 分量，即逐离子打印完整屏蔽张量 |
| `iso_shield` | 各向同性屏蔽 $\sigma_{iso} = \mathrm{Tr}(\sigma)/3 = (\sigma_{11}+\sigma_{22}+\sigma_{33})/3$（括号内同时给出 span 与 skew） |
| `span` | 张幅 $\Omega = \sigma_{33} - \sigma_{22}$ |
| `skew` | 偏斜 $\kappa = 3(\sigma_{iso} - \sigma_{22})/(\sigma_{33} - \sigma_{11}) = 3(\sigma_{iso} - \sigma_{22})/\Omega$（span 为 0 时无定义） |

  注意这些屏蔽**已包含芯电子贡献**，是与实验比较的最终关键项（例 3 会用到）。

3. **收敛结论**：总能量 500 eV 即收敛到 10 meV 以内，但屏蔽收敛慢得多——500 eV 只能到约 1 ppm，0.1 ppm 需要不低于 800 eV。目标精度要视元素而定：1H 化学位移范围约 0-14 ppm、13C 约 0-250 ppm，而 195Pt 达正负 7000 ppm，此时 0.1 ppm 精度不现实，应放宽到 1 ppm。
4. **收敛方法论**：教程**不建议对总能收敛**（此处仅为对比演示）；应对能量差（如形成能）或目标性质本身收敛。多个技术参数应**顺序收敛**：先截断能，再用收敛的截断能（此处应用 800 eV）测 k 网格；若还有第三个参数，用前两者收敛值测之。分子体系则把 k 网格换成真空尺寸收敛。
5. **k 网格扫描**：`bash kpoints_shielding.sh` 把 KPOINTS 从 4x4x4 以步长 2 扫到 16x16x16（脚本用 sed 替换 `4 4 4`，12/14/16 的结果已预先算好写入数据文件）。能量在 6x6x6 收敛到约 10 meV，而屏蔽到 0.1 ppm 需要 14x14x14。**屏蔽对 k 点的依赖远大于对截断能的依赖**：400 到 900 eV 只差约 1 ppm，而 4x4x4 到 16x16x16 差近 90 ppm，近两个数量级——对化学屏蔽而言，布里渊区的精细描述（密 k 网格）比平面波基组完备性（高截断能）更重要。

**复习问题**：化学位移与化学屏蔽的区别？需要密 k 网格说明布里渊区描述的什么重要性？需要高截断能说明基组完备性的什么重要性？

### 示例 2：诱导电流成像（苯的环电流）

**教学目标**：(1) 计算诱导电流；(2) 可视化电流分布，理解苯环电流的物理后果。

**原理**：均匀外磁场 $\mathbf B_{ext}$ 施加到材料上会感应出电流，该感应电流又经 Biot-Savart 定律产生感应磁场。化学屏蔽张量 $\sigma$ 描述感应磁场 $\mathbf B_{in}^{(1)}$ 与外场的关系：

$$\mathbf B_{in}^{(1)}(\mathbf R) = -\sigma(\mathbf R)\,\mathbf B_{ext}$$

感应磁场可由一阶诱导电流 $\mathbf j^{(1)}(\mathbf r')$ 算出（Phys. Rev. B 63, 245101 (2001)）：

$$\mathbf B_{in}^{(1)}(\mathbf R) = \frac{1}{c}\int \mathbf j^{(1)}(\mathbf r') \times \frac{\mathbf R - \mathbf r'}{|\mathbf R - \mathbf r'|^3}\, d^3r'$$

PAW 框架下诱导电流是三项之和（Phys. Rev. B 76, 024401 (2007)）：

$$\mathbf j^{(1)}(\mathbf r') = \mathbf j^{(1)}_{bare}(\mathbf r') + \mathbf j^{(1)}_{\Delta d}(\mathbf r') + \mathbf j^{(1)}_{\Delta p}(\mathbf r')$$

即赝化价波函数的贡献（bare）+ 抗磁修正项 + 顺磁修正项，后两项补偿核附近赝波函数与真实全电子波函数的偏差。

**输入要点**（目录 `e02_induced_current`，示例文件内容略）：

- POSCAR：气相苯 C6H6 置于 8x8x8 Angstrom 立方胞（分子平面垂直 x 轴）。
- INCAR：与例 1 相同，另加 `WRT_NMRCUR = 1`——把诱导电流写出为 `NMRCURBX`（磁场沿 x 轴；y/z 轴对应 `NMRCURBY`/`NMRCURBZ`）。
- KPOINTS：Gamma 点单点网格——孤立分子需要尽可能大的胞 + 仅 Gamma 采样；**真空大小**是必须收敛的参数（增大直到目标性质如化学屏蔽收敛）。

**流程与结果**：

1. `mpirun -np 4 vasp_std` 后可用 py4vasp 的 `calc.structure.plot()` 查看苯——高度对称的平面六元环（CH 单元），可思考这种对称性如何影响电流图案。
2. 化学屏蔽仍在 OUTCAR 中；诱导电流在 NMRCURBX/Y/Z 文件里。
3. 用 `plot_nmrcur_slice.py` 可视化，四个参数：电流文件（NMRCURBX/Y/Z）、网格取样间隔 n（1=每点、2=隔点）、yz 截面在 x 轴的直坐标位置（0.5=胞中点）、箭头缩放指数（箭头长度 $L$ 按 $(L/10^{4.5})\times W$ 缩放，$W$ 为轴宽）。本例用带原子位置标注的版本：`python3 plot_nmrcur_slice_c6h6.py NMRCURBX 2 0.5 0.9`，输出 `fig_nmrcurbx-slice.png`（箭头缩放 0.9 需自行微调）。

![苯分子诱导电流切片图](/imgs/2026-08-05/fig43.png)

也可用 py4vasp 画 NMR 电流：`calc.current_density.to_quiver("NMR(bx)", a=0.5) + calc.current_density.to_contour("NMR(bx)", a=0.5)`，并用 `fig[0].max_number_arrows` 控制箭头数。

4. **物理**：苯分子平面上下存在离域 $\pi$ 电子环（芳香性）；外磁场下电子沿环流动形成**环电流**，在环外产生与外场同向、环内反向的强磁场。环上的质子（H）因此经历强烈**去屏蔽**——这是芳香化合物 NMR 谱的标志性特征。

**复习问题**：哪个标签开启诱导电流的写出？写到哪些文件？苯的 H 原子是被屏蔽还是去屏蔽？

### 示例 3：与实验化学位移对比

**教学目标**：(1) 计算一组含 O 化合物的化学屏蔽；(2) 把实验化学位移对计算屏蔽作图并线性拟合；(3) 理解多化合物拟合相对单一参考的优势。

**原理**：化学屏蔽不是可观测量——NMR 测的是进动频率 $\omega_L = -\gamma B_z$。现代谱仪磁场约 1-20 T：场越强，自旋态能量差越大、稳定自旋态的布居偏置越大，信噪比越高；但强场也改变样品频率 $\omega_{sample}$。为跨谱仪比较，取参考化合物频率 $\omega_{ref}$（Pure Appl. Chem. 80, 59 (2008)）定义相对化学位移：

$$\delta = \frac{\omega_{ref} - \omega_{sample}}{\omega_{sample}} \times 10^6$$

![自旋能级在外磁场下的分裂](/imgs/2026-08-05/fig44.png)

计算上由线性响应得到的是屏蔽 $\sigma$：诱导电流（例 2）经 Biot-Savart 定律给出感应磁场（bare、抗磁、顺磁三项分别对应 $\mathbf B^{(1)}_{bare}$、$\mathbf B^{(1)}_{\Delta d}$、$\mathbf B^{(1)}_{\Delta p}$），换算成对应屏蔽后加上**化学不变的芯贡献** $\sigma_{core}$ 得总屏蔽：

$$\sigma_{tot} = \sigma_{bare} + \sigma_{\Delta d} + \sigma_{\Delta p} + \sigma_{core}$$

**输入要点**（目录 `e03_experiment`，含 BaSnO3、MgO、SrO、SrTiO3 四个子目录，示例文件内容略）：INCAR 与例 1 相同（ENCUT=500），KPOINTS 为 6x6x6 Gamma 中心网格；另附 `experiment.dat`（实验位移：MgO 47、BaSnO3 143、SrO 390、SrTiO3 465 ppm）。

**流程与结果**：

1. 全跑四个体系很耗时，**建议只跑 MgO**；可选的 `O_shielding.sh` 会遍历四个目录运行 VASP 并把各向同性屏蔽提取到 `O_shielding.dat`（对含 O3 的钙钛矿用 `grep -A 15`、二元氧化物用 `grep -A 10` 取最后一个 O 的 iso_shield）。
2. 在 OUTCAR 的 `CSA tensor` 表中：前三行是屏蔽张量，末行括号内是各向同性屏蔽。固体应取 **INCLUDING G=0 CONTRIBUTION**（该行第 4 个数）。教程给出 MgO 中 O 的参考数值（更紧设置的示例）：张量对角元 198.9458（含 G=0），各向同性屏蔽 198.9458 ppm（不含 G=0 为 187.2125）。
3. 其余三个化合物的屏蔽由教程直接提供：BaSnO3 74.0956、SrO -219.9072、SrTiO3 -281.9145 ppm。
4. **屏蔽到位移的换算**：对元素 $N$，各向同性位移与屏蔽的关系为 $\delta_{iso}[N] = \delta_{ref}[N] - \sigma_{iso}[N]$。实验中参考物是明确选定的（如 1H 用 TMS，因其 NMR 峰独特、化学性质合适）；**计算中没有这样的标准参考**——任何单个参考计算都可能因方法局限引入误差并扭曲由其导出的位移。因此把**多个实验位移对计算屏蔽做线性拟合**（J. Chem. Phys. 146, 064115 (2017)），平均掉个别计算的误差：

$$\delta^{exp}_{iso}[N] = \delta_{ref}[N] + m\,\sigma^{calc}_{iso}[N]$$

5. 本教程四化合物拟合结果：**$\sigma_{ref} = 210.94$ ppm、$m = -0.8639$、Pearson $r = -0.981$**。若 $|r|$ 低于约 0.95，先检查 `O_shielding.dat` 与 `experiment.dat` 中化合物顺序是否一致（点要对得上）。

**复习问题**：化学屏蔽与化学位移的关系？参考物在实验数据分析中的作用？能否拿单个化合物的计算屏蔽直接与实验比较？

### 示例 4：用理论预测实验

**教学目标**：(1) 计算 BaZrO3 与 CaO 两个新化合物的 O 屏蔽；(2) 用例 3 的线性拟合预测其化学位移并与实验对比；(3) 把自己的系列与更大的文献系列（Laskowski）比较，评估误差来源。

**原理**：例 3 已展示实验位移与计算屏蔽的线性关系。BaSnO3/MgO/SrO/SrTiO3 是 Laskowski 系列（J. Chem. Phys. 146, 064115 (2017)）的子集；该文献还验证了 F 的屏蔽，并通过与全电子（AE）基组对比验证了 PAW 方法本身。由拟合方程即可**从计算屏蔽预测实验位移**。

**输入要点**（目录 `e04_comparison`，含 BaZrO3/CaO 两个子目录，示例文件内容略）：INCAR 与例 1 相同（ENCUT=500）、6x6x6 网格；`experiment.dat` 含六个化合物（四化合物 + BaZrO3 376、CaO 294 ppm）。建议只跑 CaO（BaZrO3 很耗时），教程提供 BaZrO3 的计算屏蔽。

**流程与结果**：

1. 预测：把新化合物的计算屏蔽代入例 3 的拟合直线，与实验值对比：

| | BaZrO3 | CaO |
|---|---|---|
| 计算 $\sigma$ / ppm | -161.697 | -148.444 |
| 预测 $\delta$ / ppm | 350.639 | 339.188 |
| 实验 $\delta$ / ppm | 376 | 294 |
| 差值 / ppm | 25.0 | -45.0 |

2. **偏差来源**：BaZrO3 偏离源于 5s/5p **半芯态未冻结**——浅芯态在晶场中被允许极化；CaO 偏离是已知效应——空的 Ca 3d 态在 DFT 中离价带顶太近，3d 态的这种靠近导致 O 位移偏差（J. Chem. Phys. 146, 064115 (2017)）。
3. **拟合质量对比**：

| | $\sigma_{ref}$ / ppm | $m$ | $r$ |
|---|---|---|---|
| 四化合物拟合 | 210.94 | -0.8639 | -0.981 |
| 六化合物拟合 | 208.07 | -0.8592 | -0.9248 |
| Laskowski（标准） | 220.65 | -0.8724 | -0.995 |
| Laskowski（优化） | 216.67 | -0.8558 | -0.996 |

4. 解读：四化合物系列截距与斜率与 Laskowski 优化系列相差不到 1 ppm 与 0.02、Pearson 系数接近 1——是好的拟合（四化合物是 Laskowski 系列中"行为良好"的子集）。加入 BaZrO3/CaO 组成六化合物系列后，截距变小、斜率偏离文献值、Pearson 系数下降——更难的化合物引入了额外复杂性。Laskowski O 系列共 10 个化合物："标准"组用 700 eV 截断能与更硬的 `_h` POTCAR（更小芯半径），"优化"组用 900 eV 与特制改进 POTCAR——两者都优于教程的快速设置。
5. **理想斜率 $m = -1$**：意味着计算与实验只差一个参考位移 $\sigma_{ref}$；$m$ 偏离 -1 源于 DFT 对交换关联效应捕捉的局限。**赝势选择同样关键**：多数情况普通 POTCAR 足够；少数情况（如 BaZrO3）需要芯半径更小的硬赝势。

**复习问题**：加入 BaZrO3 和 CaO 对系列有何影响？PAW 赝势的选择影响化学屏蔽吗？斜率 $m$ 为什么偏离理想的 -1？

**本部分涉及的 VASP 功能与标签汇总**：
- NMR 核心标签：`LCHIMAG = .TRUE.`（线性响应屏蔽）、`LNMR_SYM_RED`（与 k 空间导数自洽的对称约化）、`NLSPLINE`（可微 PAW 投影算符）、`ICHIBARE`（有限差分模板半径）、`KPAR`（k 点并行）。
- 精度要求：`PREC = Accurate`、`EDIFF` 不高于 1E-8、高 ENCUT、`LASPH = .TRUE.`。
- 电流成像：`WRT_NMRCUR = 1` 写出 NMRCURBX/Y/Z；PAW 三项电流分解（bare/抗磁/顺磁）+ 化学不变的芯贡献 $\sigma_{core}$；Biot-Savart 由电流得感应磁场。
- 输出解读：OUTCAR `CSA tensor` 表；EXCLUDING/INCLUDING G=0（分子/固体各取其一）；`iso_shield`/`span`/`skew`（Mason 约定）；`BDIR`/`IDIR`/`IQ`（stdout 线性响应网格）。
- 方法论：对目标性质（而非总能）顺序收敛（ENCUT - k 网格 - 真空）；屏蔽对 k 网格的敏感度远高于截断能（90 ppm vs 1 ppm）；多化合物线性拟合（$\sigma_{ref}$、$m$、Pearson $r$）替代单一参考；理想斜率 $m = -1$。
- 误差来源与赝势：半芯态极化（BaZrO3）、空 d 态位置（CaO）；必要时改用硬赝势（`_h`，更小芯半径）。

## 14.2 NMR（二）：耦合常数与双中心修正

> 来源：<https://vasp.at/tutorials/latest/nmr/part2/> ｜ 练习文件：[NMR-part2.zip](https://vasp.at/tutorials/latest/NMR-part2.zip) ｜ 示例目录：`$TUTORIALS/NMR/e05_hyperfine`、`e06_efg_ediff`、`e07_efg_converge`、`e08_two_center`

本页从化学屏蔽扩展到其他 NMR 相关谱学参数：例 5 超精细耦合张量（CH3 自由基，PBE vs HSE06）；例 6 电场梯度与四极耦合常数对 EDIFF 的敏感（MAPbI3）；例 7 四极耦合常数的截断能与 k 网格收敛及与 NQR 实验对比的教训；例 8 化学屏蔽的双中心修正（LiH）。

### 示例 5：超精细常数

**教学目标**：(1) 运行超精细耦合计算；(2) 把超精细耦合参数与实验和文献对比；(3) 比较 GGA 与杂化泛函对耦合常数的影响。

**原理**：电子也有自旋。电子产生的内禀磁场与核磁偶极矩的相互作用使原本简并的能级分裂——**超精细分裂**。描述该分裂的哈密顿量为核磁偶极矩 $\boldsymbol\mu_I$ 在磁场 $\boldsymbol B$ 中：

$$\hat H_D = -\boldsymbol\mu_I \cdot \boldsymbol B$$

此处磁场是内部的，由电子轨道与自旋角动量共同贡献：$\boldsymbol B = \boldsymbol B^l_{el} + \boldsymbol B^s_{el}$。核自旋 $S^I$ 与电子自旋 $S^e$ 经超精细张量 $A^I$ 耦合（Phys. Rev. B 88, 075202 (2013)）：

$$E = \sum_{ij} S^e_i\, A^I_{ij}\, S^I_j$$

**任务**：在 10x10x10 Angstrom 立方胞中用 PBE 与 HSE06 计算气相 CH3 自由基（含一个未配对电子）的超精细耦合张量。

**输入要点**（目录 `e05_hyperfine`，示例文件内容略）：

- INCAR.pbe：`ENCUT = 500`、`ISMEAR = 0`/`SIGMA = 0.01`、`EDIFF = 1E-6`、`PREC = Accurate`；NMR 专属——`LHYPERFINE = .TRUE.`（计算超精细耦合张量）、`NGYROMAG = 10.7084 42.577478461`（按离子种类给出核旋磁比，此处选 13C 与 1H）、`ISPIN = 2`（有未配对电子，必须自旋极化）、`LASPH = .TRUE.`。
- INCAR.hse06：在 PBE 基础上加 `LHFCALC = .TRUE.` + `GGA = PE` + `HFSCREEN = 0.2`（即 HSE06）；`ALGO = Damped`（阻尼速度摩擦算法，`LHFCALC` 的推荐设置）。
- KPOINTS：Gamma 单点（孤立分子）。

**输出解读**（OUTCAR）：超精细参数分为各向同性的 **Fermi 接触**项与各向异性的**偶极**贡献。`Fermi contact (isotropic) hyperfine coupling parameter (MHz)` 之后按 POSCAR 离子顺序列出（单位 MHz）：

| 标签 | 含义 |
|---|---|
| `A_pw` | Fermi 接触项的平面波贡献 |
| `A_1PS` | 赝单中心贡献 |
| `A_1AE` | 全电子单中心贡献 |
| `A_1c` | 单中心芯贡献（Phys. Rev. B 71, 115110 (2005)） |
| `A_tot` | Fermi 接触项总贡献（**不含** `A_1c`，需显式相加） |

随后是偶极贡献张量 $A_{ij}$（`A_xx A_yy A_zz A_xy A_xz A_yz`），最后对角化写出总超精细张量（约定 $|A_{zz}| > |A_{xx}| > |A_{yy}|$，并给不对称度 $(A_{yy}-A_{xx})/A_{zz}$）。注意三个 1H 核（离子 2-4）数值几乎相同——CH3 自由基有三重旋转对称，三个质子等价。

**流程与结果**：

1. `cp INCAR.pbe INCAR` 后运行，用 grep 提取 13C 的 A_tot 与 A_1c 存入 `hyperfine.dat`。
2. **收敛行为**：与化学屏蔽形成鲜明对比——超精细参数对截断能**收敛很快**（500 eV 已收敛，并非普遍规律）。固体应用还需收敛 k 网格（如金刚石 NV- 缺陷，Phys. Rev. B 88, 075202 (2013)）。
3. **与文献/实验对比（13C）**：文献 PBE 的 $A_{tot} = 182.5$ MHz；教程值与文献之差来自 `LASPH = .TRUE.`（文献未含 PAW 球内梯度的非球贡献；关掉 LASPH 即与文献一致，只剩胞大小 10 vs 20 Angstrom 的小差异）。实验值（EPR 电子顺磁共振，即未配对电子版的 NMR）为 107.4 MHz：不含芯贡献的 $A_{tot}$（约 169 MHz）与实验差很多；加上芯贡献（$A_{tot} + A_{1c}$）后文献值变为 86.1 MHz（接近实验），教程值约 71.9 MHz，仍差 20-30 MHz。
4. **HSE06 显著改善**：教程值 $A_{tot} + A_{1c} = 92.1$ MHz（$A_{tot} = 188.8$、$A_{1c} = -96.8$），文献 101.9 MHz，明显更接近实验 107.4 MHz。物理原因：杂化泛函含精确交换，有助于**局域化孤电子或缺陷态**（J. Chem. Phys. 125, 224106 (2006)）——孤电子越局域，核处自旋密度越准，Fermi 接触项越可靠。

**复习问题**：$A_{tot}$ 上要加哪一项才能改善与实验的对比？核旋磁比由哪个标签定义？为什么杂化泛函算耦合常数优于 GGA？

### 示例 6：电场梯度（EFG）

**教学目标**：(1) 计算电场梯度；(2) 通过收紧 SCF 能量收敛判据 EDIFF 使四极耦合常数收敛，理解 EFG 对 EDIFF 的强敏感。

**原理**：四极核（自旋大于 1/2）因核形变非球而具有四极电矩，与核处的**电场梯度**耦合。EFG 是静电势 $V$ 在核位置的二阶导数：

$$V_{ij} = \frac{\partial^2 V}{\partial x_i\,\partial x_j}$$

EFG 本身不可测，但核四极耦合常数 $C_q$ 可通过 NMR、EPR 或核四极共振（NQR）谱学测量：

$$C_q = \frac{eQV_{zz}}{h}$$

其中 $e$ 为电子电荷、$Q$ 为元素与同位素特定的核四极矩、$h$ 为普朗克常量。

**任务**：对四方相钙钛矿太阳能电池材料 MAPbI3（甲基铵铅碘），用仅 Gamma 的 k 网格计算 I 与 N 的核电四极矩，考察收紧 EDIFF 的影响。

**输入要点**（目录 `e06_efg_ediff`，48 原子 Pnma 胞：4 Pb + 12 I + 4 N + 4 C + 24 H，示例文件内容略）：INCAR 同例 1（14.1）另加两个标签——`LEFG` 开启 EFG 计算；`QUAD_EFG` 按 POTCAR 原子类型定义核四极矩（单位毫靶恩 millibarn）。KPOINTS 为 Gamma 单点。

**流程与结果**：

1. `bash ediff_efg.sh` 把 EDIFF 从 1E-4 扫到 1E-8 eV（sed 替换），用 `grep -A 30 'Cq(MHz)' OUTCAR` 提取耦合常数。OUTCAR 中先有 EFG 张量（`V_ii` 为沿 i 轴分量），再往下是 NMR 四极参数表：`Cq(MHz)`、`eta`（不对称参数 $(V_{yy}-V_{xx})/V_{zz}$）、`Q (mb)`。
2. **结构要点**：I 在 Pb 周围构成八面体，分**轴向**（与 Pb 同层）与**赤道**（与甲基铵 MA 分子同层）两类位点；NQR 可区分这两种 I 环境，而 N 的环境非常相似、实验上难以分辨。
3. **Cq 随 EDIFF 的变化**：

| log(EDIFF/eV) | Cq(N) / MHz | Cq(轴向 I) / MHz | Cq(赤道 I) / MHz |
|---|---|---|---|
| -4 | 0.471 | -274.522 | -308.696 |
| -5 | 0.471 | -274.084 | -308.270 |
| -6 | 0.471 | -274.081 | -308.276 |
| -7 | 0.471 | -274.075 | -308.271 |
| -8 | 0.471 | -274.075 | -308.270 |

  轴向/赤道 I 从 1E-4 到 1E-5 eV 变化很大，之后变化小于 0.1 MHz 视为收敛；N 在 1E-4 即已收敛（0.471 MHz 不变）。**EFG 常受 EDIFF 影响强烈**；性质对技术参数的敏感度因体系而异，每组新计算都必须做收敛测试。
4. **另注**：POTCAR 选择对四极耦合常数影响可能很大（GW 赝势、`_pv`/`_sv` 半芯电子），开始收敛之前应先检查。

**复习问题**：为什么算四极耦合常数要用紧的 EDIFF（如 1E-8 eV）？

### 示例 7：四极耦合常数的收敛

**教学目标**：(1) 对 k 点与截断能收敛 EFG；(2) 为不同元素确定合适的收敛限；(3) 学会与实验对比时警惕动力学效应。

**任务**：MAPbI3 中两种高丰度四极核——14N（丰度 99.9(6)%）与 127I（100%）——可用 NQR 在固体中观测。对四方相 MAPbI3 收敛其四极耦合常数，并尽量对比实验。

**输入要点**（目录 `e07_efg_converge`，示例文件内容略）：INCAR 为例 6 的 EFG 设置（初始 `ENCUT = 200` 供扫描），KPOINTS 初始 2x2x2 Gamma 中心网格。

**流程与结果**：

1. **截断能扫描**：`bash encut_efg.sh` 把 ENCUT 从 200 eV 以 100 eV 步长扫到 600 eV（500/600 eV 结果已提供）。**警告**：实际计算绝不应使用低于 POTCAR 中 ENMAX 的 ENCUT，此处仅为演示收敛。结果：

| ENCUT / eV | Cq(N) / MHz | Cq(轴向 I) / MHz | Cq(赤道 I) / MHz |
|---|---|---|---|
| 200 | 0.468 | -543.980 | -528.387 |
| 300 | 0.446 | -546.078 | -530.393 |
| 400 | 0.425 | -546.386 | -531.625 |
| 500 | 0.423 | -548.255 | -530.985 |
| 600 | 0.422 | -548.273 | -531.025 |

  赤道 I 收敛快（300 eV 内到约 0.1 MHz），轴向 I 类似；N 因数值小需要更严的限（0.01 MHz），400-500 eV 才完全收敛——**推荐至少 400 eV、最好 500 eV**。
2. **k 网格扫描**（400 eV 固定；顺序收敛原则下正式计算应使用 500 eV）：`bash kpoints_efg.sh` 扫 1x1x1 与 2x2x2（4x4x4、6x6x6 预制）：

| k 网格 | Cq(N) / MHz | Cq(轴向 I) / MHz | Cq(赤道 I) / MHz |
|---|---|---|---|
| 1x1x1 | 0.472 | -273.800 | -308.100 |
| 2x2x2 | 0.435 | -548.155 | -531.451 |
| 4x4x4 | 0.435 | -528.959 | -549.491 |
| 6x6x6 | 0.435 | -528.579 | -549.714 |

  **Gamma 点远未收敛**（轴向 I 从 -274 跳到 2x2x2 的 -548 MHz）；N 在 2x2x2 即收敛到 0.001 MHz 内，而 I 需要 6x6x6 才到 1 MHz 内——**不同元素的收敛严格度不同，四极耦合常数本身的量级是关键因素**（N 要求 0.001 MHz、I 要求 1 MHz）。
3. **与实验对比的教训**：Franssen et al.（J. Phys. Chem. Lett. 8, 61 (2017)）用 NQR 测得 14N 的 Cq 从 -100 到 75 摄氏度由 0.11 降到 0.00 MHz，与静态计算值（约 0.44 MHz）很不同；但该文用实验结构单独计算甲胺分子得 Cq = 0.45 MHz（与教程值接近）。差异根源：甲胺在毫秒量级重新取向，在 NMR 时间尺度上很快——对四种构型取平均后，计算 Cq 随温度从 0.09 降到 0.00 MHz，与实验趋势一致。**收敛的计算只是与实验对比的一部分，还需考虑未预期的动力学平均**。

**复习问题**：四极耦合常数对布里渊区描述的依赖有多强？对 MAPbI3，平面波基组与布里渊区描述哪个影响更大？

### 示例 8：双中心屏蔽贡献

**教学目标**：(1) 计算 LiH 有/无双中心修正的化学屏蔽；(2) 对截断能收敛；(3) 比较双中心贡献对 Li 与 H 的影响量级。

**背景**：化学屏蔽通常被当作**单中心性质**——忽略相邻 PAW 球中增补电流的贡献，这在一般情况下成立。但这些**双中心贡献**有时很大（J. Chem. Phys. 139, 014109 (2013)）：尤其当对周期表顶部行元素（H、B、C、N、O、F）使用硬赝势（`*_h`）且键短时。

**输入要点**（目录 `e08_two_center`，15x15x15 Angstrom 胞中的 LiH 分子，示例文件内容略）：INCAR 同例 1（ENCUT=400 起），关键标签 `LLRAUG`——开启屏蔽张量的双中心贡献；脚本先跑 `LLRAUG = .FALSE.` 再向 INCAR 追加 `LLRAUG = .TRUE.` 重跑。

**流程与结果**：

1. `bash llraug.sh` 在 400 eV 分别跑开/关双中心修正（500/600 eV 结果后补），屏蔽写入 `llraug.dat`；Li 与 H 的各向同性屏蔽分别用 `grep -A 5/-A 10 "iso_shield"` 提取。
2. 结果（含差值）：

| ENCUT / eV | $\sigma$(Li) / ppm | $\sigma$(Li)+2c / ppm | $\sigma$(H) / ppm | $\sigma$(H)+2c / ppm | Li 差值 | H 差值 |
|---|---|---|---|---|---|---|
| 400 | 87.2927 | 87.2965 | 24.6484 | 25.8085 | 0.0038 | 1.1601 |
| 500 | 88.6826 | 88.6863 | 24.6228 | 26.1131 | 0.0037 | 1.4903 |
| 600 | 88.8883 | 88.8918 | 24.6947 | 26.2643 | 0.0035 | 1.5696 |

  Li 与 H 的屏蔽本身都在 600 eV 基本收敛；关键区别在双中心贡献的量级：Li 仅约 0.004 ppm（几乎可忽略），H 却达约 1.1-1.6 ppm——**大三个数量级**，凸显双中心贡献的重要性。
3. **方法学验证**：双中心修正只源于 PAW 球的使用；完整原子中心基组下屏蔽为 Li 88.12 ppm、H 26.34 ppm（Phys. Chem. Chem. Phys. 18, 21145 (2016)）——与含双中心修正的 PAW 结果吻合良好，说明 `LLRAUG` 正确恢复了被 PAW 分区截断的物理。

**复习问题**：周期表哪些行的元素双中心电流贡献重要？双中心贡献对 Li 还是 H 更显著？

**本部分涉及的 VASP 功能与标签汇总**：
- 超精细：`LHYPERFINE = .TRUE.`、`NGYROMAG`（核旋磁比）、`ISPIN = 2`；OUTCAR 输出分解——Fermi 接触（A_pw/A_1PS/A_1AE/A_1c/A_tot）+ 偶极张量 + 对角化总张量；芯贡献 $A_{1c}$ 必须显式加回；`LASPH` 影响与文献对比；杂化泛函（`LHFCALC`/`GGA = PE`/`HFSCREEN`/`ALGO = Damped`）局域化孤电子改善耦合常数。
- EFG/四极：`LEFG`、`QUAD_EFG`（毫靶恩）；OUTCAR 的 $V_{ii}$ 与 Cq/eta/Q 表；EFG 对 EDIFF 极度敏感（需 1E-5 eV 以下）；Cq 量级决定各元素收敛限；ENCUT 不应低于 ENMAX。
- 收敛认识：性质收敛限按量级定（N 0.001-0.01 MHz vs I 0.1-1 MHz）；Gamma 点对 EFG 远不收敛；与实验对比需考虑动力学平均（MA 重新取向）。
- 双中心：`LLRAUG`；顶部行元素 + 硬赝势（`*_h`）+ 短键时重要；H 的双中心贡献比 Li 大三个数量级；与全电子基组结果互验。

## 14.3 Part 3：芳香性（Aromaticity）

> 来源：<https://vasp.at/tutorials/latest/nmr/part3/> ｜ 练习文件：[NMR-part3.zip](https://vasp.at/tutorials/latest/NMR-part3.zip) ｜ 示例目录：`$TUTORIALS/NMR/e09_nics_benzene`、`e10_current_c4h4`、`e11_nics_c4h4`

本部分把 NMR 屏蔽计算推广到**芳香性的磁学判据**：通过核独立化学位移（NICS, nucleus-independent chemical shift）与诱导电流的方向，从磁响应角度区分芳香（苯，六元环）与反芳香（环丁二烯，四元环）分子。三个例子层层递进：例 9 苯的 NICS 等高线图（环内屏蔽），例 10 环丁二烯的诱导电流（顺磁环流的可见性问题），例 11 环丁二烯的 NICS（环内去屏蔽）并与苯并排对比。

### 示例 9：苯的 NICS 等高线图

**教学目标**：(1) 计算核独立化学位移（NICS）；(2) 绘制苯分子平面内的 NICS 等高线图；(3) 用磁学判据判断芳香性。

**原理**：外磁场 $\mathbf B_{ext}$ 感应出电流 $\mathbf j^{(1)}$，该电流又产生诱导磁场 $\mathbf B_{in}^{(1)}$。芳香分子中共轭双键形成的电子环在能量上稳定分子，电子沿环流动形成环电流：环外质子（H 核）因此强烈去屏蔽——这是芳香化合物的标志性 NMR 特征；环内则被强烈屏蔽。**空间中任意一点（不必在原子核上）的化学屏蔽称为 NICS**（Chem. Rev. 105, 3842 (2005)）。常用环中心单点 NICS 判断芳香性，但 NICS 等高线图能更直观地展示各向同性诱导磁场在分子中的分布（J. Phys. Chem. A 117, 518 (2013)）。

**任务**：在 8x8x8 Angstrom 胞中计算穿过苯分子平面的 NICS 等高线图。

**输入要点**（目录 `e09_nics_benzene`，示例文件内容略）：

- INCAR：与 14.1 例 1 的 NMR 设置相同（`LCHIMAG`/`LNMR_SYM_RED`/`NLSPLINE`、ENCUT=400、EDIFF=1E-8、PREC=Accurate、LASPH），另加两个标签：
  - `NUCIND = .TRUE.`：开启 NICS 计算；
  - `LPOSNICS = .TRUE.`：在 **POSNICS** 文件指定的位置计算 NICS。
- POSNICS 文件：第一行是 NICS 点数，随后逐行给出直接坐标（x y z）。本例取 x 固定在 0.5（分子平面）、y/z 从 0.00 到 0.99 步长 0.01 的 100x100=10000 个网格点（教程提供 Python 脚本用双循环自动生成）。
- KPOINTS：Gamma 点网格（孤立分子）。

**流程与结果**：

1. `mpirun -np 4 vasp_std` 后，NICS 值在 OUTCAR 的 `Nucleus Independent Chemical Shielding, NICS (absolute)` 与 `SYMMETRIZED TENSORS` 之后逐点给出（每点 3x3 张量）。用脚本（按 POSNICS 行数计算 `grep -A 3 nics OUTCAR` 的提取行数）提取到 `nics` 文件，再读入 POSNICS 坐标，把每点张量取迹平均得各向同性 NICS，在 y-z 平面画等高线（配色红=去屏蔽、白=0、蓝=屏蔽，等值线间隔 5 ppm，并叠加分子骨架）。

![苯分子平面的 NICS 等高线图](/imgs/2026-08-05/fig45.png)

2. **解读**：图中信息丰富——C 核周围呈红色去屏蔽区，H 核周围呈蓝色屏蔽区；**最清晰的芳香性证据在环中心：环内呈现强屏蔽**，这源于 14.1 例 2 中计算过的抗磁环流（diamagnetic/diatropic current）。

**复习问题**：用哪个文件定义 NICS 的计算位置？NICS 结果在哪个输出文件中？

### 示例 10：反芳香电流（环丁二烯）

**教学目标**：(1) 计算并可视化环丁二烯的诱导电流；(2) 用电流方向判断环丁二烯的（反）芳香性；(3) 理解平面波电流与 PAW 单中心电流的可见性差异。

**原理**：芳香性的一眼判据是数环上的 $\pi$ 电子——芳香环有 $4n+2$ 个 $\pi$ 电子（苯 6 个，$n=1$）；反芳香环有 $4n$ 个。$n=1$ 的反芳香体系即环丁二烯 C4H4（四元环）。芳香性使苯稳定为平面正六边形；反芳香性使环丁二烯能量上不稳定，从正方形畸变为矩形——$\pi$ 电子不再离域于整个环，而是定域在两个双键上。

**输入要点**（目录 `e10_current_c4h4`，示例文件内容略）：INCAR 与 14.1 例 2 完全相同（`LCHIMAG` 系列 + `WRT_NMRCUR = 1` 写出 NMRCURBX）；POSCAR 为矩形环丁二烯（4 C + 4 H，8x8x8 Angstrom 胞）；Gamma 点网格。

**流程与结果**：

1. 运行后用 py4vasp 查看结构，可与苯的六元环明确区分（四元环 + 矩形畸变）。用 `plot_nmrcur_slice_c4h4.py NMRCURBX 2 0.5 0.7`（与 14.1 例 2 的脚本相同，仅标注原子位置不同）画电流切片；也可用 py4vasp 的 `current_density.to_quiver("NMR(bx)") + to_contour` 画平面波贡献电流。

![环丁二烯电流密度切片](/imgs/2026-08-05/fig46.png)

2. **观察**：黄色区域清楚显示双键对电流密度的强烈影响（电流定域化，可用 `plot_chgcar_slice.py` 画 CHGCAR 切片佐证）。苯与环丁二烯的平面波电流看起来都绕环流动——但关键区别在于：**环丁二烯实际产生顺时针的顺磁（paramagnetic/paratropic）电流，苯产生逆时针的抗磁（diamagnetic/diatropic）电流**（Chem. Commun. 21, 2220 (2001)）。环流方向是芳香/反芳香的指示器：逆时针=芳香，顺时针=反芳香；两电流方向相反，诱导磁场也相反。
3. **为什么图中看不到顺时针电流**？因为绘图只画出了**平面波贡献**的电流（易于在精细 FFT 网格上表达）。PAW 方法把电流拆分为平面波部分与单中心（one-center, 1c）部分：抗磁 1c 贡献很小，但**顺磁 1c 贡献很大**——而 1c 贡献处于"局域"基组中、不在 FFT 网格上，因此未被画出（自旋贡献在本例很小，不构成差异来源）。
4. **数值验证**：设 `LNMRLEG = .TRUE.`，可在 OUTCAR 中查找 `plane wave contribution`、`one center paramagnetic contribution`、`one center diamagnetic contribution` 三行，直接比较各贡献的量级——磁贡献的计算相对简单，只是电流本身的可视化受网格表示限制。

**复习问题**：电流在哪个输出文件？电流在什么网格上计算（精细 FFT 网格）？芳香环的电流方向与性质（逆时针、抗磁）？反芳香环的电流方向为何在图中看不到？

### 示例 11：环丁二烯的 NICS

**教学目标**：(1) 计算环丁二烯的 NICS 并绘制等高线图；(2) 用磁学判据区分芳香性与反芳香性；(3) 对比两种 NICS 计算方式（POSNICS vs 全 FFT 网格）。

**原理**：芳香环中电流逆时针（抗磁），反芳香环中电流顺时针（顺磁）；由 Biot-Savart 定律，两者诱导的磁场方向相反。

![芳香与反芳香环流产生的诱导磁场示意](/imgs/2026-08-05/fig47.png)

图中紫色为外磁场、蓝色为电子电流、红色为诱导磁场：左侧芳香苯与右侧反芳香环丁二烯的电流方向和诱导磁场方向恰好互换。

**输入要点**（目录 `e11_nics_c4h4`，示例文件内容略）：INCAR 与例 9 基本相同，但**不使用 POSNICS/LPOSNICS**——只设 `NUCIND = .TRUE.` + `LNICSALL = .TRUE.`（`LPOSNICS = .FALSE.`）。没有 POSNICS 文件且 `NUCIND = .TRUE.` 时，NICS 在**精细 FFT 网格的所有点**上计算，结果写入专门的 **NICS 输出文件**。

**流程与结果**：

1. 用 `plot_nics_slice_c4h4.py 1 0.5 5.0` 从 NICS 文件提取数值并画等高线图；由于算了全网格点，还可以用 py4vasp 的 `calc.nics.to_contour(a=0.5)` 取任意方向的切片（如 yz 平面）。

![环丁二烯分子平面的 NICS 等高线图](/imgs/2026-08-05/fig48.png)

2. **与苯的 NICS 图并排对比**：H 核周围的蓝色屏蔽、C 核周围的红色去屏蔽两者相似（仅陡峭程度略异）；主要区别在化学键与环内——苯的 C=C 双键区域（去屏蔽晕之间的屏蔽区）因芳香稳定化而明显更屏蔽（J. Phys. Chem. A 117, 518 (2013)）；**环内差异最显著：芳香的苯环内为屏蔽，反芳香的环丁二烯环内为去屏蔽**。
3. **物理图像**（Biot-Savart 定律）：环丁二烯的顺时针顺磁电流在环内诱导与外场同向的磁场，环内去屏蔽、环外屏蔽；苯则相反，逆时针抗磁电流在环内诱导与外场反向的磁场，环内屏蔽、环外去屏蔽。（注意：所绘电流仅为平面波贡献，见例 10。）
4. **两种 NICS 计算方式的取舍**：全 FFT 网格（`LNICSALL`）计算更快且可并行；POSNICS 方式允许在感兴趣的局部区域（如化学键、氢键）布置比精细 FFT 更密的网格，适合研究局部化学环境对屏蔽的影响，例如聚合物中的氢键（Macromolecules 49, 5548 (2016)）。

**复习问题**：`LNICSALL = .TRUE.` 在哪些位置计算 NICS（精细 FFT 网格全部点）？为什么比 POSNICS 方式更快（可并行）？POSNICS 自定义位置的优势（关注区域可用更密网格）？

**本部分涉及的 VASP 功能与标签汇总**：
- NICS 计算：`NUCIND = .TRUE.` 开启核独立化学位移。
- 指定位置 NICS：`LPOSNICS = .TRUE.` + POSNICS 文件（首行点数 + 直接坐标），结果在 OUTCAR；可在关注区域布置超密网格。
- 全网格 NICS：`LNICSALL = .TRUE.`（有 `NUCIND`、无 POSNICS 时的行为），在精细 FFT 网格所有点计算，写入 NICS 输出文件，更快、可并行、可任意切片。
- 电流分解：`LNMRLEG = .TRUE.` 在 OUTCAR 输出平面波、单中心顺磁、单中心抗磁三项贡献；NMRCURBX/Y/Z 只含平面波部分（1c 顺磁贡献大但不在 FFT 网格上）。
- 芳香性磁判据：环内屏蔽 + 逆时针抗磁环流 = 芳香（苯）；环内去屏蔽 + 顺时针顺磁环流 = 反芳香（环丁二烯）；Hückel 规则（$4n+2$ vs $4n$ 个 $\pi$ 电子）与几何畸变（矩形化）对应。

## 15. 声子（Phonons）

声子是扩展周期体系中原子核的集体激发。

## 15.1 Part 1：石墨烯声子（Graphene）

> 来源：<https://vasp.at/tutorials/latest/phonon/part1/> ｜ 练习文件：[phonon-part1.zip](https://vasp.at/tutorials/latest/phonon-part1.zip) ｜ 示例目录：`$TUTORIALS/phonon/e01_static-lattice`、`e02_fin-diff-force-constants`、`e03_dispersion`

本部分以石墨烯为例，走通 VASP 声子计算的完整流程：例 1 结构弛豫与静态晶格近似的局限；例 2 有限差分法计算力常数（超胞构建、对称性利用、LREAL 讨论）；例 3 由力常数傅里叶插值得到声子色散（可公度 q 点精确、其余插值），并讨论超胞尺寸收敛与虚频的物理含义。

### 示例 1：静态晶格近似（Static-lattice approximation）

**教学目标**：(1) 列举静态晶格近似的局限；(2) 对二维材料在固定体积下弛豫离子位置与晶胞形状；(3) 设置 OUTCAR 的详细程度（verbosity）；(4) 用 py4vasp 可视化晶体结构。

**原理**：标准第一性原理计算采用静态晶格模型——离子构成固定、刚性的周期结构，电子与离子自由度经 Born-Oppenheimer 近似分离。找稳定结构的典型流程：选取初始结构 - 量子力学处理电子 - 依据 Hellmann-Feynman 定理得到的力更新离子位置。静态晶格模型能描述很多现象，但在晶格动力学主导的问题上有明显局限：许多材料中**声子对比热的贡献大于电子**；静态模型无法描述热膨胀、热导、压电性、熔化、声波传播、X 射线衍射的温度依赖强度、中子散射的某些方面，以及部分光学/介电性质（如晶体在远低于带隙处出现反射率共振）。

**任务**：在固定体积下弛豫石墨烯的离子位置与晶胞形状。

**输入要点**（目录 `e01_static-lattice`，示例文件内容略）：

- POSCAR：单层碳原子蜂窝晶格（2 原子原胞），z 方向留约 8 Angstrom 真空层抑制周期性镜像间相互作用；使用 selective dynamics（`T T F`）固定 z 方向弛豫。
- INCAR：`SYSTEM = graphene`、`NWRITE = 3`（控制 OUTCAR 输出详细程度，本例借短计算观察额外输出）、`ENCUT = 400`；电子部分 `PREC = Accurate`、`EDIFF = 1e-8`、`ISMEAR = -1`、`SIGMA = 0.2`；离子部分 `IBRION = 2`（共轭梯度结构优化）、`NSW = 100`、`ISIF = 4`。
- 电子与离子自由度迭代更新：(i) 给定结构下 Kohn-Sham 方程定电子基态；(ii) Hellmann-Feynman 定理给出力与应力；(iii) `IBRION` 决定的优化算法产生新结构；(iv) 重复直到满足 `EDIFFG` 收敛判据。
- **`EDIFF` 必须足够小**：否则电子密度收敛不足，力/应力不准，结构优化可能虚假收敛或收敛到错误解——建议用更紧的 EDIFF 复核力与应力是否真的消失。
- **`ISIF = 4` 的含义与 Pulay 应力**：固定体积、允许改变晶胞形状与离子位置。虽然体积固定，晶格参数仍可变化，满足 `ENCUT` 的倒格矢数目可能随之改变，产生 **Pulay 应力**；对策是若干离子步后重启计算，或提高 `ENCUT`。
- KPOINTS：12x12x1 Gamma 中心网格——单层材料无需在表面法向细分 k 网格（法向 k 点只描述镜像层间相互作用）。POTCAR：`PAW C_s`（半芯 s 态纳入价带）。

**流程与结果**：`mpirun -np 2 vasp_std` 后用 py4vasp 的 `structure.plot()` 可视化（`supercell=(3,3,1)` 可看到典型石墨烯层状图案）。`EDIFFG` 默认为 10 倍 `EDIFF`（此处即 1e-7）；在 OUTCAR 中搜索 `TOTAL-FORCE (eV/Angst)` 可见每步力均为零——初始结构已是（亚）稳态。后续示例沿用 POSCAR 中对称化的结构。

**复习问题**：静态晶格近似有哪些局限？垂直表面方向 k 网格应取多少细分、为什么？`NWRITE` 控制什么？要变体积计算需改哪个标签（`ISIF`）？

### 示例 2：有限差分法计算力常数

**教学目标**：(1) 为声子计算构建超胞；(2) 给出力常数定义并用有限差分法计算；(3) 判断何时设置 `LREAL = Auto`；(4) 解释对称性在有限差分算法中的作用。

**原理**：Born-Oppenheimer 能量面 $E(\{\mathbf R\})$ 依赖离子位置 $\mathbf R_I = \mathbf R^0_I + \mathbf u_I$（$\mathbf u_I$ 为相对平衡位置的位移）。围绕平衡展开：

$$E(\{\mathbf R\}) = E(\{\mathbf R^0\}) - \sum_{I\alpha} F_{I\alpha}(\{\mathbf R^0\})\, u_{I\alpha} + \sum_{I\alpha J\beta} \Phi_{I\alpha J\beta}(\{\mathbf R^0\})\, u_{I\alpha} u_{J\beta} + \mathcal{O}(\mathbf R^3)$$

谐波近似在位移小时截断到二阶。展开系数的名称与符号约定对应物理意义：一阶系数是**力**，二阶系数是 **Hessian / 力常数矩阵**：

$$F_{I\alpha} = -\left.\frac{\partial E}{\partial R_{I\alpha}}\right|_{\mathbf R = \mathbf R^0}, \qquad \Phi_{I\alpha J\beta} = \left.\frac{\partial^2 E}{\partial R_{I\alpha}\,\partial R_{J\beta}}\right|_{\mathbf R = \mathbf R^0} = -\left.\frac{\partial F_{I\alpha}}{\partial R_{J\beta}}\right|_{\mathbf R = \mathbf R^0}$$

力由 Hellmann-Feynman 定理计算。要捕捉波长大于原胞的声子模式必须使用足够大的超胞——**所有声子计算都应对超胞尺寸收敛**；超胞放大 n 倍时，k 网格按同倍数缩小，保持采样密度不变。

**任务**：计算石墨烯的力常数。

**输入要点**（目录 `e02_fin-diff-force-constants`，示例文件内容略）：

- 用 pymatgen 的 `make_supercell` 从弛豫原胞生成面内 3x3、4x4、5x5、6x6 超胞（各 18/32/50/72 个 C，z 方向 8 Angstrom 真空不变）。
- INCAR：`ENCUT = 400`、`PREC = Accurate`、`NELMIN = 5`、`EDIFF = 1e-8`、`ISMEAR = -1`、`LREAL = .FALSE.`（精确的倒空间投影）、`LWAVE/LCHARG = .FALSE.`；离子部分 **`IBRION = 6`**（触发有限差分法计算力常数与声子模式）+ **`POTIM = 0.015`**（位移幅度，Angstrom）。
- KPOINTS 等密度缩小：3x3 超胞用 4x4x1、4x4 用 3x3x1、5x5 与 6x6 用 2x2x1。

**流程与结果**：

1. 运行后 VASP 打印 ADVICE：大超胞建议改用 `LREAL = Auto`（实空间投影算符更高效）；但**追求高精度时可保留倒空间投影方案 `LREAL = .FALSE.`**。教程借机引导理解投影算符与 aliasing：赝势中需要投影的部分是把平面波投影到原子局域分量的部分；`LREAL = Auto` 需配合自动确定的投影网格（ROPT）。
2. py4vasp 可视化 `calculation.structure[:]`：基态结构 + 有限差分用的微扰结构——3x3 超胞只有 5 个结构（1 平衡 + 4 个位移）。**为什么不是 3x 离子数个自由度**？对称性把碳原子互相映射：平移把所有原子映射回 2 原子原胞；六方结构面内格矢等长、夹角 60 度，绕第三格矢的 6 重旋转对称把自由度从 6 降到 4（每个位点只剩 1 个面内自由度）；反演对称再减 2（只需位移一个位点）。`IBRION = 5` 则不利用对称性、逐原子位移。
3. 力常数矩阵经 py4vasp 的 `force_constant.print()` 读取（单位 eV/Angstrom 平方，9 分量块）。石墨烯等二维材料中面内（x、y）与面外（z）力常数不耦合，对角块外元素为零。
4. **力常数随距离的行为**：把面外力常数（2,2 分量）对原子间距作图，随近邻距离振荡于正负之间；面内力常数旋转到局域参考系——纵向 $\hat{\mathbf r}_L$（两原子连线）、横向 $\hat{\mathbf r}_T$（面内垂直）、法向 $\hat{\mathbf z}$，即 $\tilde\Phi = R^T \Phi R$，$R = (\hat{\mathbf r}_L, \hat{\mathbf r}_T, \hat{\mathbf z})$——再对距离作图。增大超胞尺寸可获得更多原子对的力常数。
5. **负力常数的含义**：二阶力常数是总能对两个离子位置的二阶导数，符号决定这对离子相对移动时总能升降。但**出现少量负力常数不足以判定亚稳**——某对离子位移可能降低能量，其他正力常数可抵消之，整体仍稳定。
6. 远距离原子对的力常数幅度衰减，声子色散最终随超胞尺寸收敛。两个著名例外：(i) 带不同电荷的离子位移形成偶极，相互作用衰减缓慢（非解析修正，Part 2 处理）；(ii) 金属/半金属中核运动改变费米能级进而改变总能，形成长程的 **Kohn 反常**。

**复习问题**：`IBRION = 5` 和 `6` 的区别？力常数的定义？`LREAL` 控制什么？超胞尺寸为何对精确声子色散重要？

### 示例 3：声子色散（Phonon dispersion）

**教学目标**：计算并绘制石墨烯的声子色散。

**原理**：谐波近似下离子哈密顿量为

$$H = \frac{1}{2}\sum_{I\alpha} M_I \dot u^2_{I\alpha} + \frac{1}{2}\sum_{I\alpha J\beta} \Phi_{I\alpha J\beta}\, u_{I\alpha} u_{J\beta}$$

运动方程 $M_I \ddot u_{I\alpha} = -\Phi_{I\alpha J\beta} u_{J\beta}$。设沿波矢 $\mathbf q$ 传播的平面波形式解，求声子模式 $\varepsilon^\mu_{I\alpha}(\mathbf q)$ 与频率 $\omega^\mu(\mathbf q)$ 归结为动力学矩阵本征值问题：

$$\sum_{J\beta} \frac{1}{\sqrt{M_I M_J}}\, \Phi_{I\alpha J\beta}\, e^{i\mathbf q\cdot(\mathbf R_J - \mathbf R_I)}\, \varepsilon^\mu_{J\beta}(\mathbf q) = \omega^\mu(\mathbf q)^2\, \varepsilon^\mu_{I\alpha}(\mathbf q)$$

超胞含 N 个离子时动力学矩阵维度为 3N；把超胞原子位置写成 $\mathbf R_I = \mathbf L_i + \mathbf r_i$（原胞格矢整数倍 + 原胞内位置）后，可在原胞基下求解，维度降为 3n（n 为原胞原子数）。与超胞**可公度**的 q 点（$\mathbf q = \frac{s_1}{n_1}\mathbf b_1 + \frac{s_2}{n_2}\mathbf b_2 + \frac{s_3}{n_3}\mathbf b_3$，$s_i = 1...n_i$，限于第一布里渊区）上声子频率是精确计算的，其余 q 点由实空间力常数插值得到。

**输入要点**（目录 `e03_dispersion`，示例文件内容略）：

- 把示例 2 各超胞的 `vaspout.h5` 与 POSCAR 复制为 `vaspin.h5`（VASP 可读的力常数格式）；
- INCAR：`PHON_NWRITE = -2`、`LPHON_DISPERSION = .TRUE.`（计算色散）、`LPHON_POLAR = .FALSE.`（非极性材料，不需要偶极修正）、`LPHON_READ_FORCE_CONSTANTS = .TRUE.`（从 vaspin.h5 读入力常数）；
- QPOINTS：线模式（`line` + `reciprocal`）定义色散路径 GAMMA - M - K - GAMMA，每段 50 个点（六角布里渊区高对称路径）。

**流程与结果**：

1. 对四个超胞分别运行 `vasp_std`（只需 1 个核，纯插值），用 py4vasp 的 `phonon.band.plot()` 绘制并叠图比较。
2. **超胞尺寸的影响**：超胞过小时 Gamma 点附近出现鞍点（轻微虚频）；这类二维材料的二次型"软"声子模式收敛特别困难，需格外小心。**负声子频率表示结构不稳定**（0 K、谐波近似下）：平衡结构力为零，但沿某声子模式畸变会降低总能（亚稳）。虚频有时是计算不准（超胞不够、`ENCUT` 过低）或截断力常数插值的假象，有时有物理意义（沿该模式弛豫得到更低对称性、更低能量的结构）；本例 Gamma 附近的负模式通常指示计算问题。
3. **5x5 超胞的 Gamma 点光学支频率与其他超胞不同**：3x3/4x4/6x6 保持了等密度 k 采样，而 5x5 无法整除、只能用更粗的 k 网格；k 采样足够密时该差异会变小。
4. **倒空间采样图像**：对给定超胞，仅在 Gamma 中心网格对应的采样点上频率精确，其余为插值：

![超胞尺寸对应的布里渊区采样点](/imgs/2026-08-05/fig49.png)

  计算 M 点频率的最小超胞为 **2x2**，计算 K 点的最小超胞为 **3x3**。
5. **模式可视化**：用 [phononwebsite](https://henriquemiranda.github.io/phononwebsite) 交互式查看——从 `vaspout.h5` 生成 `graphene.json`（`phononweb.vaspphonon.VaspPhonon(...).write_json()`），网页中可选择原胞重复次数、调节振荡幅度、显示原子运动方向矢量；点击色散曲线上任一点即可在中央面板播放对应声子模式的动画。

**复习问题**：$\mathbf q = 0$ 声子的波长在晶体中延伸多远？什么是声子色散？力常数与声子频率如何联系？

**本部分涉及的 VASP 功能与标签汇总**：
- 结构弛豫：`IBRION = 2`（共轭梯度）+ `ISIF = 4`（固定体积变形状，注意 Pulay 应力）+ `EDIFFG`（默认 10xEDIFF）；`NWRITE` 控制 OUTCAR 详细度；二维材料法向 k 点取 1。
- 有限差分力常数：`IBRION = 6`（利用对称性）/ `IBRION = 5`（不用对称性）；`POTIM`（位移幅度）；等密度 k 网格原则；`LREAL = .FALSE.`（精确）vs `Auto`（大胞高效）。
- 数据流：`vaspout.h5` - `vaspin.h5` 传递力常数；`LPHON_READ_FORCE_CONSTANTS`、`LPHON_DISPERSION`、`PHON_NWRITE`、`LPHON_POLAR`。
- 色散插值：QPOINTS 线模式；可公度 q 点精确、其余傅里叶插值；M/K 点最小超胞 2x2/3x3。
- 物理判读：虚频的含义（计算假象 vs 真实不稳定性）；力常数衰减与收敛；Kohn 反常、偶极长程作用两个例外；负力常数不等于亚稳。

## 15.2 Part 2：MgO 声子与 LO-TO 分裂

> 来源：<https://vasp.at/tutorials/latest/phonon/part2/> ｜ 练习文件：[phonon-part2.zip](https://vasp.at/tutorials/latest/phonon-part2.zip) ｜ 示例目录：`$TUTORIALS/phonon/e04_mgo-relax`、`e05_mgo-fin-diff-force-constants`、`e06_mgo-dispersion`、`e07_mgo_dielectric_becs`、`e08_mgo-dispersion-polar`

本部分以离子晶体 MgO 为例，在 Part 1 流程之上补充体材料特有的内容：例 4 晶格常数的 k 点收敛；例 5 体材料超胞力常数；例 6 声子色散与声子 DOS（q 网格收敛、原子分解）；例 7 计算介电张量与 Born 有效电荷；例 8 长程偶极-偶极相互作用的非解析处理，得到含 LO-TO 分裂的色散。

### 示例 4：晶格常数（Lattice parameter）

**教学目标**：(1) 执行晶胞体积弛豫；(2) 监测晶格常数随 k 点网格密度的收敛。

**原理与策略**：精确晶格常数需要足够大的基组（`ENCUT`）与足够的布里渊区采样。推荐做法：观察目标量（本例为弛豫晶格常数）随 k 点密度的变化，选取能给出精确结果的**最小**密度，节省计算资源——因为后续超胞有限差分计算昂贵得多。

**输入要点**（目录 `e04_mgo-relax`，示例文件内容略）：POSCAR 为 MgO 常规原胞（FCC，晶格常数初值 4.211 Angstrom，Mg 在顶点、O 在体心）；INCAR：`PREC = Accurate`、`EDIFF = 1.E-08`、`ENCUT = 500`、`ISMEAR = 0`/`SIGMA = 0.05`、`LREAL = .FALSE.`、`LWAVE/LCHARG = .FALSE.`；弛豫部分 `IBRION = 2`、**`ISIF = 3`**（体积、形状、离子全弛豫，与 Part 1 的 ISIF=4 对照）、`NSW = 100`。

**流程与结果**：

1. 脚本对 N = 3、4、6、8、10、12 生成 NxNxN KPOINTS，逐个运行 VASP 并把 `vaspout.h5` 存为 `vaspout_$N.h5`。
2. 用 py4vasp 读取各文件的体积——常规胞含 4 个原胞，晶格常数取 $(4V)^{1/3}$——对 k 网格作图，判断收敛网格。
3. 可直接从收敛计算的结果生成后续 POSCAR：`mycalc.structure.to_POSCAR()` 写入 `e05_*/POSCAR-relaxed`（教程示例用 6x6x6）。

**复习问题**：为什么要对 k 点密度收敛晶格常数？

### 示例 5：有限差分法计算力常数（体材料）

**教学目标**：(1) 为体材料声子计算构建超胞；(2) 用有限差分法计算力常数。

**输入要点**（目录 `e05_mgo-fin-diff-force-constants`，示例文件内容略）：

- 与二维材料不同，体材料沿**三个**笛卡尔方向扩大原胞；立方对称允许三方向共用同一缩放因子。脚本生成 2x2x2、3x3x3、4x4x4 超胞（16/54/128 原子）。
- INCAR：`PREC = Accurate`、`EDIFF = 1.E-08`、`ISMEAR = 0`、`LREAL = .FALSE.`、`IBRION = 6`、`POTIM = 0.015`。
- KPOINTS 等密度对应原胞 6x6x6 采样：2x2x2 超胞用 3x3x3、3x3x3 与 4x4x4 用 2x2x2。

**流程与结果**：对各超胞分别运行 `vasp_std`。3x3x3 与 4x4x4 耗时较长，教程建议先用 2x2x2 继续后续章节（但要记住它对精确色散明显不足），把大超胞挂在后台运行，最后可用更大的力常数重跑示例 6、8。

**复习问题**：体材料与二维材料构建力常数超胞的主要区别是什么？

### 示例 6：声子色散与态密度（Phonon dispersion and DOS）

**教学目标**：计算 MgO 的声子色散与声子 DOS，并对 DOS 的 q 网格收敛。

**输入要点**（目录 `e06_mgo-dispersion`，示例文件内容略）：

- 把示例 5 的 `supercell_2x2x2/vaspout.h5` 复制为 `dispersion/vaspin.h5`；INCAR：`PHON_NWRITE = -2`、`LPHON_DISPERSION = .TRUE.`、`LPHON_POLAR = .FALSE.`、**`PHON_DOS = 1`**（计算声子 DOS）、`LPHON_READ_FORCE_CONSTANTS = .TRUE.`。设置 `LPHON_READ_FORCE_CONSTANTS` 后，VASP 从 `vaspin.h5` 读入力常数，在 QPOINTS 指定的 q 点上通过傅里叶插值构建动力学矩阵、绘出色散后退出。
- QPOINTS：seekpath 给出的 FCC 标准路径 GAMMA-X-L-W-W_2-K-U（线模式，每段 50 点）。

**流程与结果**：

1. 运行后用 py4vasp 的 `phonon.band.plot(selection="Mg O")` 画色散。
2. **声子 DOS 的 q 点收敛**：在 `dos/` 目录用脚本生成 10x10x10、20x20x20、30x30x30、40x40x40 的 Gamma 中心 QPOINTS 网格，逐个插值计算 DOS（结果存 `vaspout_$i.h5`）并叠图比较：10 与 20 之间仍有明显差别；30 与 40 差别可忽略，DOS 已收敛。**插值非常快**：2x2x2 力常数在 40x40x40 q 点上仅需约 5 秒。
3. **原子分解的声子 DOS**（`phonon.dos.plot(selection="Mg O")`）：Mg 原子在低频声子模式中权重最大，O 原子在高频部分权重最大。原因：O 质量小于 Mg，其参与振动的模式频率更高；频率大小还取决于力常数强度（由声子振动模式所形变化学键的强度决定）。

**复习问题**：与超胞不可公度的 q 点用什么插值得到声子模式（力常数傅里叶插值）？决定某声子模式频率大小的两个主要量是什么（原子质量与力常数）？

### 示例 7：长程偶极-偶极相互作用——介电张量与 Born 有效电荷

**教学目标**：计算 MgO 的（离子钳制）静态介电张量与 Born 有效电荷。

**原理**：半导体/绝缘体中电子对离子的屏蔽不完全，产生**长程（LR）原子间力常数**。显式计算需要无穷大超胞；有限尺寸截断会在声子色散中引起 Gibbs 振荡。解决办法（Gonze and Lee, Phys. Rev. B 55, 10355 (1997)）：把二阶力常数分解为短程与长程两部分：

$$\Phi_{I\alpha J\beta} = \Phi_{I\alpha J\beta}^\text{SR} + \Phi_{I\alpha J\beta}^\text{LR}$$

长程部分来自总能量中离子-离子贡献的长程部分的解析导数，用 Ewald 求和技术评估——实空间项对应短程、倒空间项对应长程（Ewald 参数 $\lambda$ 控制分离，代表截断长度）：

$$\Phi_{I\alpha J\beta}^\text{LR} = \frac{4\pi e^2}{\Omega_0} \sum_\mathbf{G} \frac{(\mathbf G \cdot \mathbf Z^*_{I\alpha})(\mathbf G \cdot \mathbf Z^*_{J\beta})}{\mathbf G \cdot \epsilon^\infty \cdot \mathbf G}\, e^{i\mathbf G \cdot (\mathbf R_J - \mathbf R_I)}\, \exp\left[\frac{-\mathbf G \cdot \epsilon^\infty \cdot \mathbf G}{4\lambda^2}\right]$$

其中 $\epsilon^\infty$ 是钳制离子介电张量、$\mathbf Z^*_{I\alpha}$ 是 Born 有效电荷；$\lambda$ 的选取使指数因子在 `PHON_G_CUTOFF` 定义的 G 矢量截断球内可忽略。**长程力常数由 Born 有效电荷张量与静态介电张量完全描述。**

**输入要点**（目录 `e07_mgo_dielectric_becs`，示例文件内容略）：INCAR 在标准设置上加 **`LEPSILON = .TRUE.`**——计算钳制离子静态介电张量与 Born 有效电荷；POSCAR 为弛豫后的 MgO 原胞；KPOINTS 为 6x6x6 Gamma 中心网格。

**流程与结果**：运行后用 py4vasp 从 `vaspout.h5` 读取：`calc.born_effective_charge` 给出 Born 有效电荷（含局域场效应，单位 |e|）——Mg 为 +1.979 的对角张量、O 为 -1.979（立方对称下各向同性，且两离子互为相反数，反映 MgO 的离子性略低于形式电荷 2）；介电张量同理读出。教程提供 Python 代码把这些量写成 INCAR 可用的 `PHON_DIELECTRIC` 与 `PHON_BORN_CHARGES` 标签，供示例 8 使用。

**复习问题**：力常数长程相互作用的起源是什么（离子屏蔽不完全导致的偶极-偶极作用）？描述长程力常数需要哪些量（介电张量与 Born 有效电荷）？

### 示例 8：含长程处理的声子色散与 DOS（LO-TO 分裂）

**教学目标**：利用静态介电张量与 Born 有效电荷计算极性材料力常数的长程部分，得到含 LO-TO 分裂的声子色散。

**原理**：力常数的分解对应动力学矩阵的分解：

$$D_{I\alpha J\beta}(\mathbf q) = D^\text{SR}_{I\alpha J\beta}(\mathbf q) + D^\text{LR}_{I\alpha J\beta}(\mathbf q)$$

长程部分在倒空间的表达式（含 Born 有效电荷 $\mathbf Z^*$ 与高频介电张量 $\epsilon^\infty$，用 Ewald 参数 $\lambda$ 平滑收敛）：

$$D^\text{LR}_{i\alpha j\beta}(\mathbf q) = \frac{4\pi e^2}{\Omega_0} \sum_\mathbf{G} \sum_l \frac{\big[(\mathbf G + \mathbf q)\cdot\mathbf Z^*_{i\alpha}\big]\big[(\mathbf G + \mathbf q)\cdot\mathbf Z^*_{j\beta}\big]}{(\mathbf G + \mathbf q)\cdot\epsilon^\infty\cdot(\mathbf G + \mathbf q)}\, e^{i(\mathbf q + \mathbf G)\cdot(\mathbf l + \mathbf r_i - \mathbf r_j)}\, e^{-\frac{(\mathbf G + \mathbf q)\cdot\epsilon^\infty\cdot(\mathbf G + \mathbf q)}{4\lambda^2}}$$

**实用流程**：(1) 用有限尺寸超胞计算 $\Phi_{I\alpha J\beta}$；(2) 减去解析长程部分得 $\Phi^\text{SR} = \Phi - \Phi^\text{LR}$；(3) 用 $\Phi^\text{SR}$ 傅里叶插值得 $D^\text{SR}(\mathbf q)$；(4) 在原胞中计算 $D^\text{LR}(\mathbf q)$ 并相加。

**输入要点**（目录 `e08_mgo-dispersion-polar`，示例文件内容略）：VASP 通过 **`LPHON_POLAR = .TRUE.`** 自动完成上述处理，并用 `PHON_DIELECTRIC` 提供介电张量、`PHON_BORN_CHARGES` 提供 Born 有效电荷；把示例 5 的 `vaspout.h5` 复制为 `dispersion/vaspin.h5` 与 `dos/vaspin.h5` 后运行。

**流程与结果**：

1. 对比示例 6 的色散图：加入长程处理后，极性材料在 Gamma 点附近出现 **LO-TO 分裂**（纵向光学支与横向光学支分开），且色散插值更加平滑。
2. 同样在不同 q 网格上计算 DOS（`bash run.sh`）并检查最密网格（`vaspout_40.h5`）的原子分解结果。
3. **加分任务**：用 3x3x3、4x4x4 超胞力常数重跑示例 6、8 并对比——验证长程处理如何让插值更平滑、大超胞如何改善短程部分。

**复习问题**：加入长程力常数特殊处理后，声子色散最主要的变化是什么（LO-TO 分裂 + 更平滑插值）？

**本部分涉及的 VASP 功能与标签汇总**：
- 体积弛豫：`IBRION = 2` + `ISIF = 3`；晶格常数对 k 网格收敛（目标量收敛策略）。
- 有限差分：`IBRION = 6` + `POTIM = 0.015`；体材料三方向等比扩胞；等密度 k 网格。
- 色散/DOS 插值：`LPHON_READ_FORCE_CONSTANTS` + `vaspin.h5` + QPOINTS（seekpath FCC 路径）；`PHON_DOS = 1`；DOS 对 q 网格收敛（30x30x30 足够）；插值极快（40x40x40 约 5 秒）。
- 介电与 Born 电荷：`LEPSILON = .TRUE.`（钳制离子介电张量 + Born 有效电荷，含局域场效应）；MgO 的 $Z^* = \pm 1.979$。
- LO-TO 分裂：`LPHON_POLAR = .TRUE.` + `PHON_DIELECTRIC` + `PHON_BORN_CHARGES`（+ `PHON_G_CUTOFF` 控制 G 截断）；Gonze-Lee 短程/长程分解与 Ewald 求和；Gamma 点 LO-TO 分裂与平滑插值。

## 16. 电子–声子相互作用（Electron-Phonon Interactions）

许多体系中电子与振动（声子）自由度可分别处理（电子远快于核运动）。该处理是近似的，可通过电子–声子耦合修正。电子–声子散射在诸多应用中占主导，如半导体迁移率、室温金属电导率。

## 16.1 Part 1：微扰理论计算带隙重正化（Bandgap renormalization）

> 来源：<https://vasp.at/tutorials/latest/electron-phonon/part1/> ｜ 练习文件：[electron-phonon-part1.zip](https://vasp.at/tutorials/latest/electron-phonon-part1.zip) ｜ 示例目录：`$TUTORIALS/electron_phonon/e01_diamond`、`e02_diamond_convergence`、`e03_MgO`

电子-声子相互作用计算至少分两步：第一步在**超胞**中用有限差分计算电子-声子势，第二步在**原胞**中计算电子-声子矩阵元 $g_{mn\nu}(\mathbf k, \mathbf q)$ 与带隙重正化。该流程把昂贵的超胞计算与原胞中的性质计算分离，兼顾效率与精度。本部分三个例子：例 1 金刚石的温度依赖带隙重正化（完整两步流程），例 2 零点重正化的三参数收敛研究（q 网格、中间态数、展宽 $\delta$），例 3 极性材料 MgO 的特殊处理与 ZPR 外推。

### 示例 1：金刚石的温度依赖带隙重正化

**教学目标**：(1) 从超胞计算电子-声子势；(2) 计算声子诱导的带隙重正化；(3) 比较零点与有限温度下的带隙。

**任务与体系**：对 16 原子金刚石超胞计算零点与温度依赖的声子诱导带隙重正化。金刚石是宽间接带隙的共价绝缘体，强共价键带来高频声子与显著电子-声子耦合，是研究电子-声子物理的理想模型体系。

**第一步：超胞计算电子-声子势**（目录 `e01_diamond/supercell`，示例文件内容略）：

- POSCAR：金刚石原胞的 2x2x2 重复（16 原子；教程因算力限制用小超胞，实际应对目标性质做 2x2x2 - 3x3x3 - 4x4x4 的超胞收敛测试）。
- INCAR 两个关键标签：`IBRION = 6`（激活有限差分驱动器，`POTIM = 0.015` 位移幅度）+ **`ELPH_POT_GENERATE = True`**（计算电子-声子势）；其余 `ismear = 0`/`sigma = 0.1`/`efermi = midgap`（绝缘体把费米能钉在带隙中央）、`prec = accurate`、`ediff = 1e-7`、`encut = 400`。
- KPOINTS 只有 1 个 k 点，可用 gamma 版 VASP（`vasp_gam`）；POTCAR 用 `C_s`。
- 运行流程：VASP 先做超胞电子基态，然后逐原子沿三个笛卡尔方向位移、计算 Kohn-Sham 势的变化，结果写入 **`phelel_params.hdf5`**。
- `ELPH_POT_GENERATE` 还会生成 **`CONTCAR_ELPH`**：与 CONTCAR 类似但包含**原胞**结构，应作为下一步的 POSCAR；不满意自动选的原胞可用 `ELPH_POT_LATTICE` 自定义。
- **FFT 网格匹配**：电子-声子势在超胞精细 FFT 网格上计算，按 `ELPH_POT_FFT_MESH` 定义的原胞 FFT 网格输出；默认要求该网格与原胞电子-声子计算所用 FFT 网格**完全一致**。网格尺寸由 `ENCUT`、`PREC` 决定（也可直接设 `NGX/NGY/NGZ`）：两步保持相同的 `ENCUT` 与 `PREC`（推荐）则网格自动匹配，无需额外输入。

**第二步：原胞计算带隙重正化**（目录 `e01_diamond/el-ph`，示例文件内容略）：

- 把 `phelel_params.hdf5` 与 `CONTCAR_ELPH`（改名为 POSCAR）复制到原胞目录。
- INCAR 与超胞计算基本相同（`ENCUT`、`PREC` 完全一致），只加 `ELPH_MODE = renorm`（元标签，自动设置合理默认值）与温度列表：
  - `ELPH_RUN`：启动电子-声子计算；
  - `ELPH_SELFEN_GAPS`：自动计算带隙重正化；
  - `ELPH_SELFEN_FAN` 与 `ELPH_SELFEN_DW`：自能中包含 Fan-Migdal（FM）与 Debye-Waller（DW）两项贡献；
  - `ELPH_SELFEN_TEMPS = 0 100 200 300 400 500 600`：以 K 为单位列出温度列表，单次运行计算所有温度下的重正化。
- KPOINTS 只含少量 k 点（4x4x4）用于自洽电子最小化；电子-声子物理量（矩阵元、电子能量、声子频率）在更密的 **`KPOINTS_ELPH`** 文件指定网格（10x10x10）上计算。

**流程与结果**：

1. 输出三处：标准输出（最简摘要，含直接/间接带隙的温度依赖重正化）、OUTCAR 末尾（详细，`grep -A 8 "Direct gap"` / `"Fundamental gap"` 可查）、`phelel_params.hdf5`（机器可读，供 py4vasp 后处理：`calc.electron_phonon.bandgap[0]`，数组下标对应同一运行中的多组收敛设置）。
2. 本例数值：直接带隙 KS 值 5.586 eV，0 K 准粒子带隙 5.195 eV（ZPR -391 meV），600 K 降至 5.079 eV（-508 meV）；基本（间接）带隙 KS 值 4.179 eV，0 K 准粒子带隙 3.862 eV（ZPR -316 meV），600 K 降至 3.754 eV（-424 meV）。
3. **物理**：零温下由于量子涨落仍存在晶格振动（零点运动），产生**零点重正化（ZPR）**——即便基态计算也需要电子-声子物理才能得到正确带隙。重正化符号为负：准粒子（重正化）带隙小于 Kohn-Sham 带隙，ZPR 通常使带隙收缩。温度升高，可参与散射的声子态增多，电子态被额外散射过程重正化，带隙进一步收缩（用 py4vasp 的 `bandgap.plot("direct_renorm")` 可绘直接带隙随温度曲线）。

**复习问题**：为什么零温也有重正化（零点运动）？`CONTCAR_ELPH` 包含什么（原胞结构）？

### 示例 2：金刚石零点重正化的收敛研究

**教学目标**：(1) 对 q 点数 $N_q$、中间态数 $N_b$、展宽参数 $\delta$ 做收敛研究；(2) 理解各收敛参数如何影响自能。

**原理**：只看 ZPR 即可——温度只通过费米/玻色分布函数进入自能，ZPR 收敛通常足以保证有限温度精度。$T=0$ 的 Fan-Migdal 自能：

$$\Sigma^{\text{FM}}_{n\mathbf k}(T=0) = \frac{1}{N_q} \sum_q^{N_q} \sum_m^{N_b} \sum_\nu \frac{|g_{mn\nu}(\mathbf k, \mathbf q)|^2}{\varepsilon_{n\mathbf k} - \varepsilon_{m\mathbf k+\mathbf q} \pm \omega_{\nu\mathbf q} + i\delta}$$

**三个收敛参数**：

1. **q 点网格 $N_q$**（`KPOINTS_ELPH`）：目录 `e02_diamond_convergence/q-points` 的各子目录给出 4x4x4、6x6x6、8x8x8、10x10x10 网格（POSCAR 与 `phelel_params.hdf5` 已提供，无需重跑超胞；`ELPH_SELFEN_TEMPS = 0` 只算 0 K）。逐个运行后画 ZPR-网格图：本例基本带隙在 10x10x10 即收敛，但直接带隙尚未收敛——实际中常需更高密度，每种材料、每个带隙都不同，务必复查。
2. **中间态数 $N_b$**（`ELPH_NBANDS_SUM`）：态求和收敛慢（高能中间态的散射过程物理相关性递减）。该标签可给列表（如 `20 40 ... 200`），**单次运行全部评估**。本例 200 条带仍勉强收敛，且曲线非单调、难以判断收敛；视体系与赝势可能需要数百至上千条中间带。`ELPH_MODE = renorm` 默认设 `ELPH_NBANDS = -2`（用满所有可用能带，数量等于平面波系数数）。要进一步增多只能提高 `ENCUT`——但这会改变原胞 FFT 网格尺寸，与 `phelel_params.hdf5` 中的网格不再匹配，**必须用更高 ENCUT 重跑超胞计算**（也可用 `NGX/NGY/NGZ` 直接控制网格）。因此推荐超胞计算一开始就用较高截断能，避免多次重跑昂贵的超胞步。
3. **展宽参数 $\delta$**（`ELPH_SELFEN_DELTA`）：源自多体微扰理论中保证传播子正确极点结构的无穷小虚位移，数值上可视为处理（近）简并态的展宽参数。应尽量小——过大引入非物理贡献，过小则受浮点精度限制引入数值噪声；默认 0.01 eV 通常是好的起点。目录 `delta/` 用列表 `0.0001 0.001 0.01 0.1 1` 单次运行评估：金刚石从 0.01 eV 降到 0.1 meV 结果无显著变化，默认值足够。但 $\delta$ 应与 q 网格配套选择——更密的 q 网格引入更多近能态、分母出现更小能量，**推荐做 q 网格与 $\delta$ 的双收敛研究**。这对极性材料尤其重要：LO 模式对应的电子-声子矩阵元在长波极限发散（J. Chem. Phys. 143, 102813 (2015)），所需 q 网格可能非常密。

**复习问题**：什么是 q 点？电子-声子计算能带不够收敛怎么办（提高 ENCUT 并重跑超胞）？展宽参数为何影响 ZPR？

### 示例 3：极性材料 MgO 的带隙重正化

**教学目标**：(1) 计算 Born 有效电荷与介电张量；(2) 计算电子-声子势；(3) 计算声子诱导带隙重正化；(4) 用外推法确定零点重正化。

**原理**：MgO 是极性材料，长程静电作用使某些声子模式（LO）的电子-声子矩阵元发散（类似声子的 LO-TO 分裂），需要密得多的 q 网格与谨慎的收敛测试；与声子类似，存在把长程相互作用并入矩阵元的特殊修正方案。长程（仅偶极）电子-声子矩阵元：

$$g_{mn\nu}^{\text{LR}}(\mathbf k, \mathbf q) = i \frac{4\pi}{V_{\text{pc}}} \sum_\kappa \sqrt{\frac{\hbar}{2 m_\kappa \omega_{\nu\mathbf q}}} \sum_\mathbf G \frac{(\mathbf q + \mathbf G)\cdot \mathbf Z^\star_\kappa \cdot \mathbf e_{\kappa,\nu\mathbf q}}{(\mathbf q + \mathbf G)\cdot \boldsymbol\epsilon_\infty \cdot (\mathbf q + \mathbf G)}\, \langle \psi_{m\mathbf k+\mathbf q} | e^{i(\mathbf q + \mathbf G)\cdot(\mathbf r + \boldsymbol\tau_\kappa)} | \psi_{n\mathbf k} \rangle$$

其中 $\mathbf Z^\star_\kappa$ 是 Born 有效电荷张量、$\boldsymbol\epsilon_\infty$ 是高频介电张量、$m_\kappa$ 是原子质量、$\omega_{\nu\mathbf q}$ 与 $\mathbf e_{\kappa,\nu\mathbf q}$ 是声子频率与本征矢。声子频率/本征矢与电子 Bloch 轨道可由 VASP 在电子-声子计算中得出；**Born 有效电荷与介电张量需事先计算并作为输入**。（"NAC"即 non-analytic-term correction，原指布里渊区中心动力学矩阵非解析行为的处理（Phys. Rev. B 55, 10355 (1997)），电子-声子中有类似概念，故沿用该词。）

**三步流程**：

1. **介电性质**（`e03_MgO/nac`）：INCAR 设 `LEPSILON = .TRUE.`（配合 `EFERMI = MIDGAP`、`ENCUT = 400`、6x6x6 网格），用密度泛函微扰理论计算电场响应，得到 Born 有效电荷张量与离子钳制静态介电张量；结果在 OUTCAR 末尾（`MACROSCOPIC STATIC DIELECTRIC TENSOR ...` 与 `BORN EFFECTIVE CHARGES ...`）或用 py4vasp 从 `vaspout.h5` 读取。
2. **超胞电子-声子势**（`e03_MgO/supercell`）：与金刚石类似（2x2x2 MgO 超胞、`elph_pot_generate`），但 INCAR 中必须填入上一步得到的 `PHON_DIELECTRIC` 与 `PHON_BORN_CHARGES`（教程提供脚本自动填充）；POTCAR 用 Mg 与 `O_s`。
3. **ZPR 与 q 点收敛**（`el-ph/` 下 4x4x4 至 10x10x10 四个子目录）：把 `phelel_params.hdf5` 与 `CONTCAR_ELPH` 复制到各子目录；因介电张量已在超胞步设置，此步 INCAR 与金刚石例子相同（`elph_mode = renorm`、`elph_selfen_temps = 0`）。MgO 只有直接带隙，输出中 "direct" 与 "fundamental" 相同。

**流程与结果**：

1. 关注量是 ZPR，输出记为 `KS-QP gap (meV)`（准粒子带隙减 Kohn-Sham 带隙）。ZPR 随 q 网格明显不收敛，而且**即使继续加密 q 网格也不收敛**：LO 声子模式的电子-声子矩阵元在长波极限（$\mathbf q \to 0$）发散。好在发散可积——由长程矩阵元公式，$\mathbf q \to 0$ 时分母行为 $g^{\text{LR}} \propto |\mathbf q|/|\mathbf q|^2 = 1/|\mathbf q|$，自能积分 $\int d^3q\, |g^{\text{LR}}|^2 \propto \int d^3q\, 1/|\mathbf q|^2$ 在原点附近可积，（非绝热）自能仍有限、ZPR 良定义。
2. **外推做法**：逐步加密 q 网格，把 ZPR 对**网格间距的倒数**（$1/q$）作图，取线性区外推到 $1/q \to 0$（教程用最后两点线性外推演示）。**注意**：教程仅在线性区演示外推思想；实际需要密得多（如 32x32x32 甚至更高）的 q 网格才能进入可靠的线性区。

**复习问题**：为什么要计算 Born 有效电荷与介电张量（描述极性材料电子-声子相互作用的长程部分）？为什么仅靠密 q 网格不足以收敛极性材料的 ZPR（LO 模式矩阵元长波发散，需外推）？

**本部分涉及的 VASP 功能与标签汇总**：
- 电子-声子势生成：`IBRION = 6` + `ELPH_POT_GENERATE`（超胞有限差分，写 `phelel_params.hdf5` 与 `CONTCAR_ELPH`）；`ELPH_POT_LATTICE` 自定义原胞；`vasp_gam` 可加速单 k 点超胞步。
- FFT 网格匹配：`ELPH_POT_FFT_MESH`、`ENCUT`/`PREC`、`NGX/NGY/NGZ`——两步须一致，改 ENCUT 必须重跑超胞。
- 带隙重正化：`ELPH_MODE = renorm`（元标签）、`ELPH_RUN`、`ELPH_SELFEN_GAPS`；自能贡献 `ELPH_SELFEN_FAN` + `ELPH_SELFEN_DW`；温度列表 `ELPH_SELFEN_TEMPS`；绝缘体用 `EFERMI = MIDGAP`。
- 双网格：KPOINTS（自洽）vs `KPOINTS_ELPH`（电子-声子物理量的密网格，兼作自能积分 q 网格）。
- 收敛控制：`ELPH_NBANDS_SUM`（列表单次评估）、`ELPH_NBANDS = -2`（默认用满平面波数）、`ELPH_SELFEN_DELTA`（默认 0.01 eV，与 q 网格双收敛）。
- 极性修正：`LEPSILON`（Born 电荷 + 介电张量）+ 超胞步的 `PHON_DIELECTRIC`/`PHON_BORN_CHARGES`；LO 矩阵元 $1/|\mathbf q|$ 发散可积；ZPR 对网格间距倒数外推。
- 物理图景：零点运动与 ZPR、带隙随温度收缩、FM/DW 两项自能、金刚石 ZPR 约 -316 meV（间接）/-391 meV（直接）。

## 16.2 Part 2：统计采样法计算电子–声子相互作用

> 来源：<https://vasp.at/tutorials/latest/electron-phonon/part2/> ｜ 练习文件：[electron-phonon-part2.zip](https://vasp.at/tutorials/latest/electron-phonon-part2.zip) ｜ 示例目录：`$TUTORIALS/electron-phonon/e04_SIMPLE_ONE-SHOT_C`、`e05_C_supercell_convergence`、`e06_C_temperature_dependence`、`e07_C_Monte_Carlo_sampling`

与 Part 1 的微扰理论路线不同，本部分用**统计采样**方法：温度 $T$ 下的可观测量 $O(T)$ 是对不同坐标组采样后的统计平均。完整 Monte Carlo（MC）采样中每组位移由平衡位置加随机位移得到：

$$x_T^{\text{MC},i} = x_{\text{eq}} + \Delta\tau^{\text{MC},i}$$

位移由笛卡尔坐标、原子质量 $M_\kappa$、原子 $\kappa$ 上的声子本征模 $\nu$、本征频率 $\omega_\nu$ 与正态分布随机数决定（声子本征模与频率在谐波近似下计算）。MC 需要大量结构、计算昂贵。Zacharias 与 Giustino 基于"超胞越大所需 MC 结构数越少"的经验观察提出 **one-shot（OS）方法**（ZG 构型）：只用一组位移（Phys. Rev. B 94, 075125 (2016)）：

$$\Delta\tau^{\text{OS}} = \sqrt{\frac{1}{M_\kappa}} \sum_\nu^{3(N-1)} (-1)^{\nu-1}\, \varepsilon_{\kappa,\nu}\, \sigma_{\nu,T}$$

本征模按本征频率升序求和（交替符号），位移幅度为

$$\sigma_{\nu,T} = \sqrt{(2n_{\nu,T}+1)\frac{\hbar}{2\omega_\nu}}, \qquad n_{\nu,T} = \left[\exp\left(\frac{\hbar\omega_\nu}{k_B T}\right) - 1\right]^{-1}$$

其中 $n_{\nu,T}$ 为 Bose-Einstein 占据数。one-shot 方法把 MC 求和约化为单次计算（New J. Phys. 20, 123008 (2018)）。

### 示例 4：金刚石的 one-shot 计算

**教学目标**：(1) 用 one-shot 方法创建（含特殊构型的）位移结构；(2) 解读电子-声子相互作用对带隙附近本征态的影响。

**输入要点**（目录 `e04_SIMPLE_ONE-SHOT_C`，16 原子金刚石常规胞，示例文件内容略）：

- 两套 INCAR：PBE 基态（`INCAR.dft`：`ALGO = Normal`、`ISMEAR = 0`/`SIGMA = 0.1`/`EFERMI = MIDGAP`、`LREAL = A`、`NCORE = 2`）与 one-shot 构型生成（`INCAR.os_disp`），后者关键标签：
  - `IBRION = 6`：获得动力学矩阵及其本征矢（one-shot 位移的基础）；
  - `PHON_LMC = .TRUE.`：启用 one-shot 与 Monte Carlo 计算；
  - `PHON_NSTRUCT`：大于 0 时设定 MC 采样的结构数；**等于 0 切换为 one-shot**；
  - `PHON_TLIST = 0.0`：温度列表——一次计算可同时生成多个温度的 one-shot POSCAR（动力学矩阵只需算一次）；
  - `PHON_NWRITE = -2`、`PHON_DOS = 1`：绘制声子态密度；
- KPOINTS 为 Gamma 单点（可用 `vasp_gam`）；QPOINTS 仅声子 DOS 计算需要，one-shot/随机采样的电子-声子计算不需要。

**流程与结果**：

1. 构型生成步产出 0 K 位移结构 `POSCAR.T=0.`（只有原子位置移动）；声子 DOS 的峰对应本征模能量，正是 one-shot 位移幅度所用的 $\omega_\nu$。
2. 分别对原胞与 one-shot 结构做标准 DFT（`vasp_gam`，OUTCAR 分别存为 `OUTCAR.2x2x2` 与 `OUTCAR.2x2x2_T0`），检查 Gamma 点本征值：
   - **完整超胞**：带隙附近因结构对称性本征能级强简并——价带顶三条带（band 30-32）简并于 10.136 eV，导带底六条带（band 33-38）简并于 14.677 eV，带隙约 **4.54 eV**；
   - **one-shot 结构**（占据数按 Bose-Einstein 分布定义）：简并全部解除（对称性破缺），价带 9.751-10.513 eV、导带 14.080-14.981 eV，带隙约 **3.57 eV**，比完整结构低近 1 eV——直观展示电子-声子耦合对带隙的影响。
3. 文献中常不用最高价带与最低导带之差作带隙，而是取原先简并的价带（band 30-32）与导带（band 33-38）的**平均 band energy** 之差，以改善与实验的可比性和收敛性（New J. Phys. 20, 123008 (2018)）。

**复习问题**：one-shot 方法为何要求解动力学矩阵？占据数用什么分布计算（Bose-Einstein）？电子-声子相互作用对带隙附近本征态有什么影响（解除简并、带隙收缩）？QPOINTS 文件为何使用、one-shot 计算是否必需？

### 示例 5：one-shot 方法的超胞收敛

**教学目标**：(1) 学会做简单的收敛测试；(2) 体会超胞尺寸对 one-shot 与随机采样方法的重要性；(3) 理解收敛设置的连锁影响。

**要点与流程**（目录 `e05_C_supercell_convergence`，示例文件内容略）：

- 声子 q 网格必须大到覆盖胞内所有声子；one-shot/MC 方法用**超胞尺寸**控制 q 网格大小。大超胞计算量增长快，因此做收敛测试、取最小够用尺寸。
- 对 2x2x2、3x3x3、4x4x4 超胞各做三步计算：完整超胞带隙 - 生成 one-shot POSCAR（`POSCAR.NxNxN_T0`）- 在位移结构上算重正化带隙；重正化 = 两者之差。
- **收敛设置的注意事项（重要）**：测试对象应是目标性质（带隙、力、声子模式）。求导会放大误差：力（能量一阶导）、动力学矩阵（二阶导）及其依赖量（声子模式）的误差层层累积。教程用很松的设置求快，初始小误差（核数/机器差异）可使 `POSCAR.NxNxN_T0` 产生约 **10 mAngstrom** 的结构差异，进而使 ZPR 差约 **0.05 eV**。因此只收敛带隙不够，**还必须收敛声子模式**，确保位移结构本身不变。改进方向：提高 `ENCUT`（如 600 eV）与 `PREC`、把 `ALGO` 从 Fast 改为 Normal、收紧 `EDIFF`。小 k 网格 + 低截断会引入大 **Pulay 应力**（基组不完备，力不收敛，声子模式也不收敛），可用 `grep "external pressure" OUTCAR.*` 检查。
- 带隙按文献标准做法取平均简并价带/导带能量之差（脚本用 awk 取末 3 条占据带与末 6 条空带的均值再相减），`grep 'fundamental gap' OUTCAR.xxx` 可提取带隙；结果存入 `gap_vs_supercell_size.dat` 并与参考数据 `gap_vs_supercell_size_ref.dat` 对比作图。
- 可视化位移：红色为原超胞、青色为位移后超胞：

![one-shot 方法位移前后 2x2x2 超胞对比](/imgs/2026-08-05/fig50.png)

- 带隙重正化随超胞尺寸收敛缓慢且**不单调**；用更严格设置后金刚石间接带隙可收敛（见 New J. Phys. 20, 123008 (2018) 表 2，及其他化合物表 3）。

**复习问题**：为什么要做收敛测试？得到带隙重正化需要几次计算（完整超胞 + one-shot 结构各一次，每个超胞尺寸共两次）？

### 示例 6：金刚石带隙的温度依赖

**教学目标**：(1) 计算温度依赖的带隙重正化；(2) 把温度依赖拟合到合适的经验曲线。

**流程与结果**（目录 `e06_C_temperature_dependence`，示例文件内容略）：

- 以 4x4x4 超胞为起点，用 `PHON_TLIST` 在单次计算中生成多个温度下的 one-shot POSCAR（`POSCAR.T=xxx`）；`run_calcs.sh` 对每个结构做基态计算（0-700 K 范围，教程为省时只保留 0/200/400 K，注释掉其余温度），`evaluate_gap.sh` 用平均简并带能量差计算各温度的带隙重正化，存入 `gap_vs_temp.dat`。
- 半导体/绝缘体带隙温度依赖常用 **modified-Varshni 关系**（Physica 34, 149 (1967)）：

$$E_{\text{GAP}}(T) = E_{\text{GAP}}(0) - \frac{A T^{4}}{(B+T)^2}$$

  其中 $A$、$B$ 为材料拟合参数。拟合显示计算带隙能被 Varshni 关系很好描述。
- **缺失的效应**：体积膨胀（热膨胀）贡献未包含，可通过准谐波（quasi-harmonic）计算补充。

**复习问题**：半导体/绝缘体带隙温度依赖用什么关系描述（modified-Varshni）？完整描述带隙温度依赖还需要什么（体积膨胀/准谐波效应）？

### 示例 7：Monte Carlo 采样 vs one-shot

**教学目标**：(1) 用 Monte Carlo 采样计算带隙；(2) 与 one-shot 方法比较。

**输入要点**（目录 `e07_C_Monte_Carlo_sampling`，示例文件内容略）：三套 INCAR——基态（`INCAR.dft`）、one-shot（`INCAR.elphon_oneshot`）、MC（`INCAR.elphon_mc`）；`PHON_NSTRUCT` 设定 MC 结构数；与 one-shot 用 `PHON_TLIST` 不同，**MC 的温度由 `TEBEG` 选择**。

**流程与结果**：

1. 先生成 0 K one-shot 结构并计算（结果存为 `OUTCAR.oneshot_T0` 参考）；再用谐波分布的密度矩阵生成 30 个 0 K 随机结构（`POSCAR.T=0.X`，Phys. Status Solidi A 188, 1209 (2001)），脚本批量做基态计算；`evaluate_gap.sh` 对 MC 结构带隙取算术平均，并与 one-shot 比较。
2. 教程参考数值：MC 平均 ZPR 随结构数——10 个 -0.366 eV、15 个 -0.366 eV、20 个 -0.370 eV、25 个 -0.369 eV、30 个 -0.369 eV；one-shot ZPR 为 -0.347 eV（读者自算结果会略有差异，如 -0.358/-0.351 eV）。**20 个结构时 MC 已相对 30 结构结果收敛到第二位小数**；one-shot 带隙与 MC 平均带隙非常接近，且超胞越大两者越趋于一致（New J. Phys. 20, 123008 (2018)）。
3. **实践建议**：超胞足够大时两者结果相近，推荐**始终使用 one-shot**；MC 采样仅用于特殊场合，例如需要多个有限温度谐波采样结构作为分子动力学起点。

**复习问题**：MC 采样得到含电子-声子作用带隙的计算步骤是什么？one-shot 与 MC 结果何时互相收敛（超胞足够大时）？

**本部分涉及的 VASP 功能与标签汇总**：
- 动力学矩阵：`IBRION = 6`（one-shot/MC 位移的基础）；声子 DOS：`PHON_DOS`、`PHON_NWRITE`、QPOINTS（仅 DOS 需要）。
- 统计采样开关：`PHON_LMC = .TRUE.`；`PHON_NSTRUCT`（0 = one-shot，大于 0 = MC 结构数）；`PHON_TLIST`（one-shot 多温度 POSCAR 一次生成）；`TEBEG`（MC 采样温度）。
- 带隙提取：`grep 'fundamental gap'`；平均简并价带/导带能量差（文献标准，改善与实验的可比性）。
- 收敛认识：超胞尺寸控制 q 网格；求导放大误差（力 - 动力学矩阵 - 声子模式）；位移结构对松设置敏感（10 mAngstrom - 0.05 eV）；必须同时收敛声子模式；Pulay 应力检查（`grep "external pressure"`）。
- 温度拟合：modified-Varshni 关系；缺热膨胀项，需准谐波补充。
- 方法选择：one-shot 优先；MC 用于需要多构型的场合（如 MD 起点）。

## 16.3 Part 3：电子–声子矩阵元与 VASP+phelel 工作流

> 来源：<https://vasp.at/tutorials/latest/electron-phonon/part3/> ｜ 练习文件：[electron-phonon-part3.zip](https://vasp.at/tutorials/latest/electron-phonon-part3.zip) ｜ 示例目录：`$TUTORIALS/electron_phonon/e08_matrix_elements`、`e09_phelel_band`、`e10_phelel_MgO`

本部分有两条主线：(1) 把电子–声子计算的核心原料——**电子–声子矩阵元**——直接输出到磁盘并可视化，学会把矩阵元图与电子、声子能带对照解读；(2) 引入 **phelel/velph 自动化工作流**，从仅有原胞 POSCAR 与 POTCAR 出发，自动完成 Born 有效电荷/介电张量、超胞电子–声子势、原胞自能计算的全流程，并用它复现 Part 1 手动流程的结果（金刚石能带、MgO 带隙重正化），体会两种工作流的等价性与取舍。

电子与晶格振动的相互作用修正电子能量，由电子自能描述。二阶微扰理论把自能拆为 Fan–Migdal（FM）与 Debye–Waller（DW）两项，二者都通过电子–声子矩阵元计算：

$$g_{mn\mathbf{k},\nu\mathbf{q}} = \left\langle \psi_{m\mathbf{k}+\mathbf{q}} \left| \partial_{\nu\mathbf{q}} V_{\text{KS}} \right| \psi_{n\mathbf{k}} \right\rangle$$

其中 $(\nu,\mathbf{q})$ 为声子模式的支指标与波矢，$m,n$ 是波矢分别为 $\mathbf{k}+\mathbf{q}$、$\mathbf{k}$ 的电子能带，$\partial_{\nu\mathbf{q}} V_{\text{KS}}$ 是 Kohn-Sham 势对声子模式 $(\nu,\mathbf{q})$ 的导数（电子–声子势）。常规自能计算中矩阵元即算即弃以保证效率；只有需要单独分析或绘制矩阵元时才显式输出。

### 示例 8：获取电子–声子矩阵元（ELPH_DRIVER = mels）

**教学目标**：(1) 计算电子–声子矩阵元；(2) 沿 q 路径绘制矩阵元；(3) 把矩阵元图与电子能带、声子能带结构联系起来。

**任务**：用 4x4x4 k 网格，沿 $\Gamma$-X-U-K-$\Gamma$-L-W-X 的 q 路径绘制金刚石（原胞）的平均电子–声子矩阵元。

**输入要点**（目录 `e08_matrix_elements`，示例文件内容略）：

- POSCAR：2 原子金刚石原胞（FCC 原胞基矢，晶格常数约 3.56 Angstrom），POTCAR 为 `PAW_PBE C_s 06Sep2000`；
- INCAR 分两组。常规设置：`ISMEAR = 0`、`SIGMA = 0.1`、`EFERMI = MIDGAP`（把费米能级钉在带隙中央，便于指认价/导带）、`PREC = ACCURATE`、`EDIFF = 1e-7`、`ENCUT = 400`。电子–声子设置：
  - `ELPH_RUN = .TRUE.`：启动电子–声子计算；
  - `ELPH_DRIVER = mels`：不计算自能，而是累积所有满足选择定则的矩阵元并写入 `vaspelph.h5`（HDF5 文件，布局见 VASP wiki）；
  - `ELPH_SELFEN_IKPT = 1`：按索引指定计算自能所用的 k 点，本例只取 $\Gamma$ 点；
  - `ELPH_SELFEN_BAND_START = 1` / `ELPH_SELFEN_BAND_STOP = 4`：自能计算涉及的电子能带范围（本例 1–4 带）；
  - `ELPH_NBANDS = 4`：电子–声子驱动器在密 k 网格上计算的能带数；
- KPOINTS_ELPH：line 模式、每段 21 个点、reciprocal 坐标，依次连接 $\Gamma$(0,0,0)→X(0.5,0,0.5)→U(0.625,0.25,0.625)→K(0.375,0.375,0.75)→$\Gamma$→L(0.5,0.5,0.5)→W(0.5,0.25,0.75)→X，构成面心立方布里渊区的标准高对称路径；
- 教程提供现成的 `phelel_params.hdf5`（含电子–声子势，可跳过超胞计算）与绘图脚本 `elph_mels.py`。

**流程与结果**：

1. 运行 `mpirun -np 4 vasp_std`，很快完成。`ELPH_DRIVER = mels` 模式下 VASP 不组装自能，只把沿 q 路径的矩阵元写入 `vaspelph.h5`；
2. 用 `elph_mels.py` 中的 `VaspElphG.from_file("vaspelph.h5")` 读取数据，配合 plotly 画四张图：
   - 图 1：电子能带 $\varepsilon_m(\mathbf{k}+\mathbf{q})$（数据字段 `eig_kp`，固定 $\mathbf{k}=\Gamma$）；
   - 图 2：声子能带 $\omega_{\nu,\mathbf{q}}$（数据字段 `pheigs`，第 5 支高亮）；
   - 图 3：平均矩阵元（`elph.get_averaged_g()`），取最高电子带 $m=n=4$；由于简并，前 9 个 q 点上 band 3 与 band 4 相连，脚本手工把这段换为 band 3 的数据；平均是对 $m=n=4$ 的全部简并态（含 $n=2,3,4$）进行的，平均方式见 `elph_mels.py` 或 Phys. Rev. B 100, 174304 (2019) 的式 (97)；
   - 图 4：固定声子支 $\nu=6$ 的矩阵元沿整条路径的散点图；
3. 物理图像：$\Gamma$ 点的电子态 $|n,\mathbf{k}\rangle$ 通过发射或吸收一支波矢 $\mathbf{q}$、支指标 $\nu$ 的声子，与 $|m,\mathbf{k}+\mathbf{q}\rangle$ 态耦合，耦合强度即 $g_{mn\mathbf{k},\nu\mathbf{q}}$；
4. 沿声子能带追踪 $\nu=6$ 可发现：若按**带指标**追踪，声子支交叉处能带导数突变，矩阵元曲线随之出现不连续；若按**物理态**追踪（在交叉处切换支指标 $\nu$），曲线连续。即不连续来自指标排序而非物理。

**复习问题**：(1) 如何沿声子能带与矩阵元图追踪同一个声子态？(2) 按声子带指标追踪为何出现不连续？这些不连续是否物理？

### 示例 9：金刚石电子能带（velph 入门）

**教学目标**：(1) 用 `velph` 生成 `velph.toml` 配置文件并自动准备能带、DOS 计算目录；(2) 绘制金刚石能带与态密度；(3) 认识 phelel/velph 自动化工作流的结构。

**背景**：Part 1 的手动流程依赖预先准备好的分目录输入文件，`ELPH_POT_GENERATE` 路线虽然直接但也有局限。[phelel](https://github.com/phonopy/phelel) 提供 VASP 电子–声子计算的高层接口，其命令行工具 `velph` 可自动完成：

- NAC（非解析项修正）所需的 Born 有效电荷与介电张量计算设置；
- 超胞生成与电子–声子势计算；
- 原胞电子–声子计算的准备；
- 文件管理并维持各步计算间设置的一致性。

这套集成工作流降低手工搭建与数据传递出错的可能；`phelel` 在原子位移选择上也更灵活，构建电子–声子势时使用更多对称操作，可缩短运行时间。本例演示 C 金刚石的完整 VASP+phelel 工作流，起点只有原胞 POSCAR 与 POTCAR。

**任务**：用 velph 以 k-spacing 0.6 搭建金刚石原胞的计算，并对 2x2x2 超胞计算能带结构与 DOS。

![velph 工作流示意](/imgs/2026-08-05/fig51.png)

**输入要点**（目录 `e09_phelel_band`）：velph 工作流只需要三样东西——

- POSCAR：金刚石原胞（与 Part 1 练习 2 介电计算同一结构）；
- POTCAR：`PAW_PBE C_s 06Sep2000`（与 Part 1 相同）；
- `template.toml`：为 velph 提供默认值之外的输入设置，本例在 `[vasp.incar]` 设 `encut = 400`，并在 `[vasp.phelel.incar]` 设 `elph_prepare = true`、`lwap = false`。

**流程与结果**：

1. **熟悉命令**：`velph --help` 与 `velph hints` 查看全部子命令与提示，`velph init --help` 查看初始化选项；
2. **初始化**：在示例目录执行 `velph init POSCAR . --max-num-atoms=100 --kspacing=0.6 --force --template-toml template.toml`。`velph init` 需要 POSCAR、项目目录与超胞尺寸（`--max-num-atoms` 限制超胞最大原子数）；`--kspacing=0.6` 对应 INCAR 的 `KSPACING`（倒空间 k 点间距，自动生成 KPOINTS），0.6 恰好给出与前面 MgO 例子相同的 2x2x2 超胞采样。命令产出 `velph.toml`，包含后续所有计算的全部设置，且大多与 INCAR 标签一一对应：`[vasp.*.incar]`、`[vasp.*.kpoints]` 小节，文件末尾还有 `[symmetry]`、`[unitcell.*]`、`[lattice]` 等结构信息。例如 `[vasp.el_bands.incar]` 给出 `ibrion = -1`、`nsw = 0`、`lorbit = 11`、`nedos = 5001`、`ismear = 0`、`sigma = 0.05`、`ediff = 1e-08`、`encut = 400`、`prec = "accurate"`、`lreal = false`、`lwave = false`、`lcharg = false`；`[vasp.el_bands.kpoints]` 为 5x5x5 网格、`[vasp.el_bands.kpoints_opt]` 为 51 点 line 路径、`[vasp.el_bands.kpoints_dense]` 为 61x61x61 密网格。这些设置都可手工修改，后续任何 `velph ... generate` 都会使用它们；
3. **生成并运行能带计算**：`velph el_bands generate` 产出 `el_bands/bands` 与 `el_bands/dos` 两套输入文件；先在 `el_bands/bands` 运行 `mpirun -np 4 vasp_std`，用 py4vasp 查看：`py4vasp.Calculation.from_path(".../el_bands/bands").band.plot("kpoints_opt")`；
4. **DOS 与并排绘图**：在 `el_bands/dos` 运行 VASP 后，用 `velph.cli.el_bands.plot_el_bandstructures(window=(-23, 18), vaspout_filename_bands=..., vaspout_filename_dos=...)` 把两个 `vaspout.h5` 的能带与 DOS 并排绘制：

![金刚石电子能带与态密度](/imgs/2026-08-05/fig52.png)

本例完成的 velph 工作流三步走：(1) `velph init` 初始化 `velph.toml`；(2) `velph el_bands generate` 生成 `el_bands/bands`、`el_bands/dos` 目录；(3) 分别运行能带与 DOS 计算——仅凭 POSCAR 与 POTCAR 两个文件就得到两个完整计算。下一例将用同样的工作流计算 MgO 的带隙重正化。

**复习问题**：(1) 使用 velph 时去哪里找计算的 INCAR 设置？（`velph.toml` 的 `[vasp.*.incar]` 小节）(2) 这种"单一配置文件驱动全流程"的结构有什么优点？

### 示例 10：极性材料 MgO 的介电性质与带隙重正化（velph 全流程）

**教学目标**：(1) 用 velph 计算极性材料的介电性质（Born 有效电荷与介电张量）；(2) 用集成工作流计算带隙重正化；(3) 与 Part 1 手动流程的结果对比，确认两条路线完全等价。

**任务**：用 velph 生成 `velph.toml`（kspacing=0.6），计算 Born 有效电荷与介电张量，自动生成 2x2x2 超胞、计算电子–声子势，再在原胞做电子–声子计算，得到 MgO 在 0–700 K 范围内的带隙重正化。

**输入要点**（目录 `e10_phelel_MgO`）：与示例 9 一样只需原胞 POSCAR（岩盐结构 MgO 原胞，Mg、O 各 1 个原子，与 Part 1 练习 3 同一结构）与 POTCAR（`PAW_PBE Mg 13Apr2007`、`PAW_PBE O_s 07Sep2000`，与 Part 1 相同），外加提供默认设置的 `template.toml`。

**流程与结果**：

1. **初始化**：`velph init POSCAR . --max-num-atoms=100 --kspacing=0.6 --force --template-toml template.toml` 生成 `velph.toml`；
2. **NAC/介电性质**：用 `grep -A 14 "\[vasp.nac\]" velph.toml` 查看 `vasp.nac` 小节——`vasp.nac.incar` 含基本收敛标签与 `LEPSILON = True`（DFPT 计算介电张量与 Born 有效电荷），`vasp.nac.kpoints` 定义 k 网格。教程初始化给出的是较粗网格，需手工把它改为 6x6x6 以与 Part 1 一致。然后 `velph nac generate` 生成 `nac` 目录（核对 INCAR 标签与 KPOINTS 确为 6x6x6），运行 VASP；用 py4vasp 打印结果：
   - clamped-ion 宏观静态介电张量：对角元均为 **3.167509**（无量纲，含局域场效应）；
   - Born 有效电荷：Mg 为 $+1.97612 \times$ 单位阵，O 为 $-1.97612 \times$ 单位阵（|e| 单位）；
   - 与 Part 1 示例 3（`e03_MgO/nac`）的手动结果逐项对比**完全一致**——两条工作流跑的是同一个计算，唯一区别是组织方式。velph 随后自动知道 Born 电荷与介电张量在哪里，并把它们带入后续步骤（NAC 修正）；
3. **超胞电子–声子势**：`velph phelel --help` 了解子命令（该步用 phelel 做超胞有限差分）。`velph phelel init` 找出所需原子位移，生成 `phelel` 目录与配置文件 `phelel_disp.yml`（写明位移哪些原子、应用哪些对称性）；`velph phelel generate` 据此为每个位移生成一个子目录（另加一个未微扰几何的子目录），每个子目录含完整 VASP 输入，彼此仅 POSCAR（原子位移）不同；每个子目录的 INCAR 都含新标签 `ELPH_PREPARE`——把局域势等必要信息写入 `vaspout.h5`。子目录个数取决于初始化时的对称性设置。用循环逐个运行：`for d in disp-*; do cd $d; mpirun -np 4 vasp_std; cd ..; done`。全部完成后回到项目根目录执行 `velph phelel differentiate`：它以正确设置调用 phelel、定位各 VASP 计算，执行有限差分并生成 `phelel_params.hdf5`。VASP+phelel 在位移方案选择上比 VASP-only 更灵活；
4. **原胞自能计算**：查看 `velph.toml` 的 `[vasp.selfenergy.incar]` 小节（最终电子–声子计算的 INCAR 标签），按需修改：把 `[vasp.selfenergy.kpoints_dense]` 的密 k 网格改为 8x8x8（未完全收敛但足以给出合理数值）；把 `ELPH_SELFEN_TEMPS` 设为 `[0, 100, 200, 300, 400, 500, 600, 700]`（0–700 K）；删除 `ELPH_NBANDS_SUM` 并设 `ELPH_NBANDS = -2`（不研究中间态数收敛，直接用全部能带）。这些修改可以手工编辑，也可以用 tomlkit 脚本批量完成（读入 `velph.toml`、修改 `doc["vasp"]["selfenergy"]["incar"]` 与 `kpoints_dense.mesh` 后写回）。然后 `velph selfenergy generate` 生成 `selfenergy` 目录——核对 INCAR 含上述标签、KPOINTS_ELPH 为 8x8x8 网格——运行 `mpirun -np 4 vasp_std`；
5. **结果对比**：用 py4vasp 读取 `calc.electron_phonon.bandgap[0]`，把本例（e10）与 Part 1 示例 3（`e03_MgO/el-ph/8x8x8`）的 `direct_renorm` 曲线叠加绘图：示例 3 只算了温度曲线的第一个点（把精力放在 k 网格收敛上），其 0 K 值与本例 0 K 点**完全一致**；本例进一步展示 ZPR 随温度的演化——温度升高时直接带隙重正化的绝对值增大，带隙持续收缩。这证明 velph 工作流用另一种组织方式搭建了完全相同的计算；
6. **工作流选择建议**：VASP-only 流程适合熟悉 VASP 手工计算的用户；VASP+phelel 适合熟悉 phonopy 风格工作流、刚接触 VASP 的用户，且目前灵活性更高；选哪条路线主要看个人偏好。

**复习问题**：(1) VASP+phelel 工作流在哪些环节辅助了你？(2) 什么情况下你会选 VASP-only，什么情况下选 VASP+phelel？

完成本部分后你应当能够：绘制电子–声子矩阵元并追踪一条声子支；用 velph 搭建并执行完整 VASP+phelel 工作流；自动计算 Born 有效电荷与介电张量；以最少手工干预生成超胞电子–声子势；用集成工作流计算带隙重正化。

### 功能与标签汇总（Electron-phonon Part 3）

| 功能 | 关键标签 / 工具 | 说明 |
|---|---|---|
| 启动电子–声子计算 | `ELPH_RUN` | 总开关 |
| 矩阵元输出 | `ELPH_DRIVER = mels` | 不算自能，矩阵元按选择定则累积写入 `vaspelph.h5` |
| 自能 k 点/能带选择 | `ELPH_SELFEN_IKPT`、`ELPH_SELFEN_BAND_START/STOP` | 按索引挑选 k 点与能带区间 |
| 密网格能带数 | `ELPH_NBANDS` | `-2` = 使用全部能带 |
| 温度列表 | `ELPH_SELFEN_TEMPS` | 一次计算多个温度的自能/带隙重正化 |
| q 路径定义 | KPOINTS_ELPH（line 模式） | 矩阵元沿 q 路径绘制 |
| 工作流管理 | `velph init/generate/differentiate` + `velph.toml` | 单配置文件驱动全流程，`[vasp.*.incar/kpoints]` 小节对应 INCAR/KPOINTS |
| 默认值模板 | `template.toml` | 覆盖 velph 默认（如 `encut`、`elph_prepare`、`lwap`） |
| 超胞位移 | `phelel_disp.yml`、`ELPH_PREPARE` | 每位移一个子目录，局域势等写入 `vaspout.h5` |
| 有限差分 | `velph phelel differentiate` | 自动调用 phelel 生成 `phelel_params.hdf5` |
| NAC 集成 | `velph nac` + `LEPSILON` | 自动衔接 Born 电荷/介电张量到后续步骤 |
| k 采样控制 | `--kspacing`（对应 `KSPACING`） | 倒空间间距自动决定各步 k 网格 |
| 结果读取/绘图 | py4vasp、`velph.cli.el_bands.plot_el_bandstructures` | `vaspout.h5` 为统一数据载体 |

---

## 16.4 Part 4：铁的电导率（Conductivity of iron）

> 来源：<https://vasp.at/tutorials/latest/electron-phonon/part4/> ｜ 练习文件：[electron-phonon-part4.zip](https://vasp.at/tutorials/latest/electron-phonon-part4.zip) ｜ 示例目录：`$TUTORIALS/electron-phonon/e11_conductivity_iron_crta`、`e12_conductivity_iron_serta`、`e13_conductivity_iron_exp`

本部分以 bcc 铁为载体，系统学习**金属电导率的第一性原理计算**：先用常数弛豫时间近似（CRTA）快速得到电导率并与实验对比，再做 k 网格收敛测试、加入铁磁性；然后切换到自能弛豫时间近似（SERTA），从第一性原理计算散射寿命，理解寿命分布如何改变电导率；最后考察晶格常数、密度泛函与赝势三个"非物理"因素对结果的巨大影响，体会"与实验吻合可能出于错误的原因"。

金属中电输运源于载流子（电子与空穴）在外电场下的运动。电导率 $\sigma$ 量化电流通过材料的难易程度，是电阻率 $\rho$ 的倒数 $\sigma = 1/\rho$（不要与密度混淆）。由线性化 Boltzmann 方程，输运分布函数为：

$$\mathcal{T}(\epsilon) = \frac{e^2}{N_\mathbf{k}\Omega} \sum_{n\mathbf{k}} \tau_{n\mathbf{k}}\, \mathbf{v}_{n\mathbf{k}} \otimes \mathbf{v}_{n\mathbf{k}}\, \delta(\epsilon_{n\mathbf{k}}-\epsilon)$$

其中 $e$ 为元电荷，$N_\mathbf{k}$ 为 k 点数，$\Omega$ 为体积，$\tau_{n\mathbf{k}}$ 是态 $(n,\mathbf{k})$（能带 $n$、k 点）的弛豫时间，$\mathbf{v}_{n\mathbf{k}}$ 为群速度分量，$\epsilon_{n\mathbf{k}}$ 为电子能量。Onsager 输运系数 $\mathcal{L}_{ij}$（$i,j \in \{1,2\}$，[J.M. Ziman, Electrons and Phonons (2001)](https://global.oup.com/academic/product/electrons-and-phonons-9780198507796)）为：

$$\mathcal{L}_{ij} = \int d\epsilon\, \mathcal{T}(\epsilon)\, (\epsilon-\mu)^{i+j-2} \left(-\frac{\partial f^0}{\partial \epsilon}\right)$$

由此得电导率 $\sigma = \mathcal{L}_{11}$（式 11.3）、Seebeck 系数 $S = \mathcal{L}_{11}/\mathcal{L}_{12}$（式 11.4）、电子热导 $\kappa_e = \frac{1}{T}(\mathcal{L}_{22} - \mathcal{L}_{21}\mathcal{L}_{11}^{-1}\mathcal{L}_{12})$（式 11.5）。热导还有晶格贡献，此处省略（Phys. Rev. B 85, 024102 (2012)）。

### 示例 11：常数弛豫时间近似（CRTA）

**教学目标**：(1) 在 CRTA 下计算简单金属的电导率并与实验对比；(2) 测试电导率对输运 k 网格的收敛；(3) 做自旋极化计算以反映铁的磁性，并与非磁结果比较。

**任务**：在 CRTA 下计算 bcc 铁的电导率（含自旋极化与否两套），并对输运系数所用 k 网格做收敛。

**原理**：CRTA 假设费米面上所有电子共享同一个平均散射寿命 $\tau$（弛豫时间）。$\tau$ 只是费米面电子速度平均值的整体缩放因子，因此 CRTA 是一种快速而优雅的方式，可以看到**能带几何本身**如何控制输运。金属费米面上的电子本就可以自由移动（没有带隙），CRTA 是合理的起点。

**输入要点**（目录 `e11_conductivity_iron_crta/crta`，示例文件内容略）：

- POSCAR：bcc 铁原胞（PBE 晶格常数），POTCAR 为 `PAW_PBE Fe 06Sep2000`；
- INCAR：常规部分 `PREC = Accurate`、`EDIFF = 1e-8`、`ISMEAR = -15`、`SIGMA = 0.1`；输运部分：
  - `ELPH_RUN`：启动电子–声子计算；
  - `ELPH_MODE = TRANSPORT`：一键布置输运计算，自动设定一组合理默认值（具体见 wiki）；
  - `ELPH_SELFEN_CARRIER_DEN = 0.0`：额外载流子（电子/空穴）密度，本例不掺杂；
  - `ELPH_SCATTERING_APPROX = CRTA`：弛豫时间近似方式（CRTA 不需要声子计算）；
  - `TRANSPORT_RELAXATION_TIME = 1e-14`：输运弛豫时间（秒），CRTA 的唯象寿命；
- KPOINTS：12x12x12 Gamma 网格（自洽基态）；KPOINTS_ELPH：24x24x24（非自洽 + 输运，比基态网格更密）。

**流程与结果**：单次 `mpirun -np 4 vasp_std` 运行内部完成五步：

1. 电子最小化获得基态解；
2. 在 KPOINTS_ELPH（或 `ELPH_KSPACING`）指定的另一套 k 网格上做非自洽计算；
3. 用位置算符与哈密顿量的对易子计算电子群速度；
4. 在 Gauss–Legendre 网格上评估输运函数、计算 Onsager 系数；
5. 计算输运系数。

OUTCAR 中可逐步核对细节：输运参数块（`elph_transport = T`、`elph_transport_driver = 2` 即 Gauss–Legendre 积分、`elph_transport_nedos = 501`、`transport_relaxation_time = .100E-13`、`elph_velocity_mode = 2`）；输运能量窗口 `[6.421 : 7.461] eV`（宽 1.040 eV）；413 个 k 点中选出 412 个计算自能，涉及能带 3–6；Transport calculator 块报告：CRTA、501 个积分点、平均弛豫时间 1.0000E-14 s、电子数与空穴数均为 3.8104E-01。最终输运表（找 `sigma S/m`）按温度列出 `mu eV`、`sigma S/m`、`mob cm^2/(V.s)`、`seebeck uV/K`、`peltier uV`、`kappa_e W/(m.K)`：

- 300 K：$\sigma = 9.25\times10^6$ S/m，迁移率 4.27 cm2/(V s)，Seebeck 12.9 uV/K，$\kappa_e = 70.37$ W/(m K)；0–500 K 范围内电导率从 $9.11\times10^6$ 缓升至 $9.47\times10^6$ S/m；
- 对比实验：铁的电阻率 97.1 nOhm m（CRC Handbook）对应 $\sigma = 1.03\times10^7$ S/m——第一次计算就相当接近！电子热导参考值 80.4 W/(m K)，计算值 70.4 W/(m K) 也合理；
- 但教程明确警告：这种吻合是**巧合**——CRTA 恰好在 300 K 附近与实验一致，后面会看到为什么。

**k 点收敛**（`kpoints/NxNxN` 目录，N = 10, 20, 30, 40, 50, 60）：

- 两套网格分工：KPOINTS 用于自洽基态 DFT，KPOINTS_ELPH 用于非自洽 DFT 与输运计算（也可用 `ELPH_KSPACING` 替代 KPOINTS_ELPH，类似 KSPACING 替代 KPOINTS）。收敛测试只变 KPOINTS_ELPH；INCAR 新增 `ELPH_ISMEAR = -24`：电子–声子计算前确定费米能级与化学势所用的展宽方法；
- 电导率对 k 采样极其敏感，原因有二：(1) 它依赖费米面几何；(2) Dirac delta 函数 $\delta(\varepsilon_{n\mathbf{k}}-\varepsilon_F)$ 只挑选恰在费米能级上的态。加密网格相当于给模糊的费米面图像"提高分辨率"；
- 运行方式：基态网格固定 12x12x12，`kpoints.sh` 脚本把基础计算的 WAVECAR 拷入每个 NxNxN 目录再逐个运行；用 `grep 300.00000 kpoints/*/OUTCAR | grep Gauss ...` 可快速提取各网格 300 K 电导率，或用 py4vasp 的 `calc.electron_phonon.transport[0].electronic_conductivity()` 循环读取并绘图（叠加实验线 1.04e7 S/m）；
- 结果：收敛并不完全平滑（60x60x60 还略有回升）。收敛慢的典型情形：**复杂费米面**拓扑、**各向异性材料**（不同晶向收敛速度不同，bcc 铁不属此类）、**轻有效质量**（高曲率能带）。生产计算必须监控电导率随 k 网格的变化；判据应以"接近收敛值"而非"接近实验值"为准。教程取与收敛值偏差 5% 以内为足够——对 bcc 铁即 KPOINTS_ELPH 用 30x30x30 或更密。

**磁性**（`kpoints_mag/NxNxN` 目录）：铁是磁性材料，在位磁矩约 2.2 $\mu_B$——前面的非磁计算漏掉了关键成分，属于"因错误的原因得到正确结果"。本节重做共线自旋极化（铁磁）计算：

- INCAR 增加 `ISPIN = 2`、`MAGMOM = 6`（初始磁矩），并行标签 `KPAR`（同时处理的 k 点数）；同一目录还演示了声子相关标签：`IBRION = 6`（有限原子位移）、`ELPH_POT_GENERATE`（由有限位移计算电子–声子势）、`LPHON_DISPERSION`（沿 QPOINTS 文件给定路径计算声子色散）；
- 运行 10x10x10、20x20x20、30x30x30，与预计算的 40x40x40 一起，和非磁结果同图对比；
- 结果：非磁曲线似乎更贴近实验值，但对铁的物理准确描述**必须**包含磁性；金属费米面复杂，两种情形都需要很密的 k 网格才能收敛。

**复习问题**：(1) 输运计算应把 `ELPH_MODE` 设为什么？(2) 铁的电导率为何对 k 网格密度收敛缓慢？(3) 若非磁解更接近实验，能否忽略铁的磁性？

### 示例 12：自能弛豫时间近似（SERTA）

**教学目标**：(1) 在 SERTA 下计算简单金属的电导率；(2) 测试 k 网格收敛并与 CRTA 对比；(3) 把散射寿命对 Kohn–Sham 能量作图，理解寿命分布如何影响电导率。

**任务**：对铁的 2x2x2 超胞，在 SERTA 下计算电导率，并对输运 k 网格做收敛。

**原理**：已知所需 k 网格密度后，可进阶到**从第一性原理计算弛豫时间**：在超胞中计算 Kohn–Sham 势对离子位移的导数，从而涵盖尽可能多的声子。输运计算改用 `ELPH_SCATTERING_APPROX = SERTA`（Phys. Rev. B 97 (2018) 121201(R)）：

$$\frac{1}{\tau^{\text{SERTA}}_{n\mathbf{k}}} = \sum_{n'\nu\mathbf{k}'} \frac{2\pi}{\hbar} w_{n\mathbf{k},n'\mathbf{k}'} |g^{\nu}_{n\mathbf{k},n'\mathbf{k}'}|^2 \left[ (n_{\nu\mathbf{q}}+1-f_{n'\mathbf{k}'})\, \delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{n'\mathbf{k}'}-\hbar\omega_{\nu\mathbf{q}}) + (n_{\nu\mathbf{q}}+f_{n'\mathbf{k}'})\, \delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{n'\mathbf{k}'}+\hbar\omega_{\nu\mathbf{k}'}) \right]$$

其中 $\tau^{\text{SERTA}}_{n\mathbf{k}}$ 是态 $(n,\mathbf{k})$ 的弛豫时间（散射时间/寿命），$w_{n\mathbf{k},n'\mathbf{k}'}$ 为权重（SERTA 取 1），$g^\nu_{n\mathbf{k},n'\mathbf{k}'}$ 是电子–声子耦合矩阵元，$f_{n\mathbf{k}}$ 是电子占据（Fermi–Dirac 分布），$n_{\nu\mathbf{q}}$ 是声子占据（Bose–Einstein 分布），两个 delta 函数分别对应声子发射与吸收过程。寿命进入输运分布函数（式 11.1）→ Onsager 系数（式 11.2）→ 电导率（式 11.3），因此 SERTA 中每条能带的能量 $\varepsilon_{n\mathbf{k}}$ 都会影响电导率。

**流程与结果**：

1. **超胞电子–声子势**（`e12_conductivity_iron_serta/supercell`）：POSCAR 为常规胞的 2x2x2 超胞（立方，$a = 5.4819$ Angstrom，16 个 Fe 原子）；INCAR 要点：金属展宽 `ISMEAR = 1`、`SIGMA = 0.1`，自旋极化 `ISPIN = 2`、`MAGMOM = 6.0*16`，`NBANDS = 96`，`KPAR = 4`，声子部分 `LPHON_DISPERSION = .TRUE.`、`IBRION = 6`、`POTIM = 0.015`、`ELPH_POT_GENERATE = TRUE`；KPOINTS 仅 2x2x2（超胞）。先用 `python3 kpath.py -o QPOINTS` 生成声子色散路径文件，再运行 `vasp_std`，产出后续电子–声子/输运计算必需的 `phelel_params.hdf5`；用 py4vasp 的 `mycalc.phonon.band.plot()` 可查看铁的声子色散；
2. **声子限制输运**（`kpoints/NxNxN`，N = 10, 15, 20, 30, 40）：计算形态与 CRTA 完全相同，区别只有两点——传入超胞的 `phelel_params.hdf5`，并设 `ELPH_SCATTERING_APPROX = SERTA`（不再需要 `TRANSPORT_RELAXATION_TIME`，因为寿命由第一性原理给出）。POSCAR 换回 1 原子 bcc 原胞，`MAGMOM = 2.0`（原胞 1 个原子），KPOINTS 12x12x12；`ELPH_SELFEN_CARRIER_DEN` 定义额外载流子浓度（cm-3 单位）。`kpoints.sh` 把 `phelel_params.hdf5` 拷入各目录并运行。教程只要求跑 10x10x10 与 15x15x15；20x20x20、30x30x30 约需 15 分钟，40x40x40 约 1 小时，教程提供预计算的 `vaspout.h5`；
3. **结果**：把 300 K 电导率对 k 网格作图（叠加示例 11 的 CRTA 磁性曲线与实验线）：SERTA 电导率比实验值**高约 40%**——计入声子散射反而恶化了与实验的一致性，"更多物理"有时只是"更多问题"；
4. **寿命分析**：联系式 11.1–11.3 可知寿命的变化必然改变电导率。用 py4vasp 读取 40x40x40 计算的 Fan 自能（`calc.electron_phonon.self_energy[0].to_dict()` 的 `fan` 字段），按 $\tau = \hbar/(2|\mathrm{Im}\,\Sigma^{\text{FM}}|)$ 换算寿命（$\hbar = 6.582119569\times10^{-16}$ eV s），对 KS 本征能量作图，并标出化学势（费米能级）竖线与 CRTA 的 $\tau = 10^{-14}$ s 横线：
   - 紫线即 CRTA 假设：所有载流子同一寿命 $10^{-14}$ s；
   - 实际载流子寿命分布很宽，高能量处离散尤其大——这正是 SERTA 与 CRTA 的本质区别：SERTA 的 $\tau$ 是一组随态变化的值，CRTA 固定为单一常数；
   - 多数寿命小于 $2\times10^{-14}$ s，但图右上区域有相当多的载流子寿命高达约 $9\times10^{-14}$ s。寿命越长电导率越高：SERTA 纳入了这些高能、高迁移率载流子，电导率因此大幅升高，与实验偏差更大；CRTA 忽略了它们，与实验吻合纯属巧合而非物理描述更好。可换用其他 k 网格观察收敛对寿命图的影响；
5. **网格说明**：KPOINTS 定义电子结构自洽所用网格；电子–声子计算的 k 与 q 网格由同一个 KPOINTS_ELPH 文件决定——网格越密，$\varepsilon_{n\mathbf{k}}$ 与对应寿命越多，散射通道越多，相当于提高图像分辨率。

**复习问题**：(1) CRTA 与 SERTA 哪个使用平均寿命？(2) 纳入一个寿命范围如何影响电导率？(3) k 与 q 网格如何定义？

### 示例 13：晶格常数、密度泛函与赝势的影响

**教学目标**：(1) 在实验晶格常数下计算铁的电导率；(2) 理解不同密度泛函如何影响电子结构与输运性质；(3) 比较标准与 GW 两种铁赝势的电导率结果。

**任务**：生成电子–声子势后，在实验晶格常数下用 SERTA（10x10x10、20x20x20、30x30x30 网格）计算 bcc 铁的电导率，同时变换赝势（standard vs GW）与密度泛函（CA/LDA vs PBE）。

**为什么要考察这些因素**：

- 金属电导率与能带结构内在关联，而能带结构依赖晶格。把晶格常数从 PBE 计算值换成实验值（CRC Handbook）会同时引起：(1) **带宽变化**——轨道重叠改变，影响导带带宽；(2) **费米面改变**——费米面大小与形状随体积变化，直接影响载流子数目与速度；(3) **态密度变化**——费米能级态密度 $N(\varepsilon_F)$ 随体积变化，影响电导率；
- 赝势选择同样显著影响电子性质与输运系数。不同赝势的差异在于：(1) **核–价电子划分**——多少电子显式处理、多少冻结在核内；(2) **截断能需求**——更"硬"的赝势需要更高平面波截断；(3) **能带精度**——对输运而言费米能级附近的精度尤其重要。本例将评估：更精确赝势能给最终结果带来多大差别？换交换关联泛函会引入多大相对误差？

**输入要点**（目录 `e13_conductivity_iron_exp/`）：目录结构为 4 种组合 × (超胞 + 3 套 SERTA 网格)：

- `exp_std_pbe`、`exp_std_lda`、`exp_gw_pbe`、`exp_gw_lda`，命名规则 `exp_PP_DF` = 实验晶格常数 + 赝势（std/GW）+ 泛函（LDA/PBE）；每个组合含 `supercell_exp_PP_DF` 与 `serta_exp_PP_DF_NxNxN`（N = 10, 20, 30）；
- 四种 POTCAR：LDA+GW 用 `PAW Fe_sv_GW 05Dec2013`，PBE+GW 用 `PAW_PBE Fe_sv_GW 05Dec2013`，LDA+std 用 `PAW Fe 03Mar1998`，PBE+std 用 `PAW_PBE Fe 06Sep2000`（GW 版本均为 Fe_sv，半芯态作价电子）；
- 注意：密度泛函默认由 POTCAR 文件定义，本例 INCAR 中无需再指定；`phelel_params.hdf5` 已预先提供，无需重跑超胞（如需重跑，`conductivity.sh` 中取消相应注释，约 5 分钟/个）；
- SERTA 输入与示例 12 相同（`ELPH_MODE = TRANSPORT`、`ELPH_SCATTERING_APPROX = SERTA`、`ISPIN = 2`、`MAGMOM = 2.0`）；超胞 KPOINTS 2x2x2，SERTA KPOINTS 12x12x12，KPOINTS_ELPH 分别为 10/20/30；
- 共 16 个计算（12 个 SERTA + 4 个超胞）。教程只要求运行 4 个 10x10x10 计算（`bash conductivity.sh` 循环 4 种组合：拷入 `phelel_params.hdf5`、运行 `vasp_std`）；20x20x20 每个约 3–5 分钟、30x30x30 约 10–30 分钟，全部跑完约 20–30 分钟，其结果以预计算 `vaspout.h5` 提供。

**流程与结果**：用 py4vasp 循环读取 4 种组合在 10/20/30 网格下的 300 K 电导率，连同示例 11 的 CRTA（磁性）与示例 12 的 SERTA（磁性、计算晶格常数）曲线，画成三并排子图：(1) PBE 下 std 与 GW 赝势对比；(2) std 赝势下 LDA 与 PBE 对比；(3) 全部计算 + "无穷 k 点极限"外推点（对 1/k 做线性外推，拟合优度 R2：PBE+std 0.9967、PBE+gw 0.9598、LDA+std 0.9999、LDA+gw 0.9859、CRTA 0.9776，线性度很好；40x40x40 已收敛的 SERTA 磁性曲线不再外推，直接取末点）。主要结论：

- **赝势是主导因素**：std 赝势随 k 网格加密逐渐**低估**实验电导率，GW 赝势**高估**；无穷 k 极限下这一格局不变，但 GW 更接近实验。GW 赝势更精确，制备时考虑了未占据态，能在高能区复现全电子散射性质；
- **泛函影响较小**：从 LDA（CA）换到 GGA（PBE）差别不大，两者都相对实验低估，且越接近收敛越明显；无穷极限下几乎无差别。LDA 对金属常表现良好，这一致性符合预期——但若用 LDA 弛豫结构会得到较差的晶格常数，泛函选择仍须谨慎；
- **晶格常数影响巨大**：对比同一设置（PBE + std 赝势）在计算晶格常数（示例 12）与实验晶格常数（本例）下的 SERTA 结果，差别约 **50%**，足以把预测从低估翻转为高估。理想情况下用实验晶格常数应最接近实验，但这里并非如此；用实验值还是泛函自洽值（不同泛函还不一样），目前并无成熟惯例；
- **综合结论**：示例 11 的 CRTA 曲线贴近实验完全是"因错误的原因"——常数弛豫时间（物理描述不如 SERTA）、k 网格未收敛（收敛后一致性反而变差）、赝势（std）不适合本问题、晶格常数非实验值，多重误差恰好抵消。与实验吻合很容易出于错误原因，必须通过系统的收敛测试逐一确定每个因素的误差量级。如有时间，可用收敛设置重算 CRTA，看它离实验有多远；
- 其他可考虑因素：超胞尺寸、不同磁相、缺陷影响、自洽步 k 网格、电子–电子散射（post-DFT 方法）。

**复习问题**：(1) 如何获得电导率的"无穷 k 点极限"？(2) GW POTCAR 考虑了标准 POTCAR 缺失的哪些态？

**Part 4 小结**：输运性质需要比基态计算密得多的 k 采样；CRTA 是电导率计算的良好起点；晶格常数通过能带结构对电导率产生显著影响。下一部分（Part 5）转向半导体中声子限制的迁移率。

### 功能与标签汇总（Electron-phonon Part 4）

| 功能 | 关键标签 | 说明 |
|---|---|---|
| 输运计算 | `ELPH_MODE = TRANSPORT` | 一键布置输运计算并设定默认值 |
| 散射近似 | `ELPH_SCATTERING_APPROX` | CRTA（唯象常数寿命）/ SERTA（第一性原理寿命） |
| 常数弛豫时间 | `TRANSPORT_RELAXATION_TIME` | CRTA 的全局寿命（s），SERTA 不需要 |
| 载流子密度 | `ELPH_SELFEN_CARRIER_DEN` | 额外电子/空穴浓度（cm-3） |
| 输运展宽 | `ELPH_ISMEAR` | 电子–声子计算前确定费米能级/化学势的展宽 |
| 输运 k 网格 | KPOINTS_ELPH / `ELPH_KSPACING` | 同时决定电子–声子的 k 与 q 网格，需远密于基态 |
| 磁性输运 | `ISPIN` + `MAGMOM`、`KPAR` | 铁磁 bcc 铁必须自旋极化；KPAR 为 k 点并行 |
| 电子–声子势 | `IBRION = 6` + `POTIM` + `ELPH_POT_GENERATE` | 超胞有限位移法，产出 `phelel_params.hdf5` |
| 声子色散 | `LPHON_DISPERSION` + QPOINTS | 沿给定 q 路径输出声子色散 |
| 寿命分析 | Fan 自能虚部 | $\tau = \hbar/(2|\mathrm{Im}\,\Sigma^{\text{FM}}|)$，经 py4vasp 读取 |
| 赝势/泛函变量 | PAW_PBE Fe / Fe_sv_GW 等 | 泛函默认由 POTCAR 定义 |

---

## 16.5 Part 5：半导体声子限制迁移率与 ZT 热电优值

> 来源：<https://vasp.at/tutorials/latest/electron-phonon/part5/> ｜ 练习文件：[electron-phonon-part5.zip](https://vasp.at/tutorials/latest/electron-phonon-part5.zip) ｜ 示例目录：`$TUTORIALS/electron-phonon/e14_mobility_si_crta`、`e15_mobility_si_serta`、`e16_thermoelectric_bi2te3`

本部分从金属转向**半导体与热电材料**：先在 CRTA 下计算硅的电子/空穴迁移率，系统比较输运分布函数的两种能量积分方案（线性网格 + Simpson vs Gauss–Legendre）并做能量网格与 k 网格收敛；再用 SERTA 从第一性原理声子散射重算迁移率；最后以 Bi2Te3 为例计算热电输运系数与 ZT 优值（含自旋轨道耦合），学会把 ZT 对化学势、载流子浓度与温度作图。

与零带隙金属（天然存在大量载流子）不同，半导体在 $T \neq 0$ 时本征载流子极少，需要通过掺杂提高电荷密度（n 型加电子、p 型加空穴）。载流子在材料中运动的难易程度用**迁移率** $\mu_e$（电子）或 $\mu_h$（空穴）描述，它是半导体器件性能的关键指标。把能带求和限制在导带或价带，可把电导率拆成电子与空穴贡献 $\sigma_e$、$\sigma_h$，并定义：

- 电子迁移率 $\mu_e = \sigma_{n \in \text{CB}} / n_e$，载流子密度 $n_e = \frac{1}{\Omega N_\mathbf{k}} \sum_{\mathbf{k}, n \in \text{CB}} f(\varepsilon, T, \mu)$；
- 空穴迁移率 $\mu_h = \sigma_{n \in \text{VB}} / n_h$，载流子密度 $n_h = \frac{1}{\Omega N_\mathbf{k}} \sum_{\mathbf{k}, n \in \text{VB}} [1 - f(\varepsilon, T, \mu)]$；

其中 $f$ 为 Fermi–Dirac 分布，$\mu$ 为给定温度下的化学势。**符号约定**：化学势恒记作 $\mu$；迁移率恒带下标（$\mu_e$、$\mu_h$），二者不要混淆。本部分用线性化 Boltzmann 输运方程（BTE）计算迁移率，重点关注电子–声子相互作用的影响（Rep. Prog. Phys. 83 036501 (2020)）。CRTA 迁移率只受能带几何影响（载流子轻重、能谷各向异性、应变引起的能带重排），是一种"能带结构指纹"，反映晶体在**没有散射**时的表现。

### 示例 14：CRTA 计算迁移率（硅）

**教学目标**：(1) 在 CRTA 下计算硅的迁移率；(2) 对 Onsager 系数所用能量积分网格（线性与 Gauss–Legendre）和 k 网格做数值收敛测试；(3) 比较输运分布函数两种积分方案。

**任务**：在 CRTA 下计算金刚石结构 Si 的电子迁移率，对能量网格（线性与 Gauss–Legendre）与 k 网格做收敛。

**输入要点**（目录 `e14_mobility_si_crta/crta`，示例文件内容略）：

- POSCAR：2 原子 Si 原胞（FCC，晶格常数约 5.43 Angstrom）；POTCAR 为 `PAW_PBE Si 05Jan2001`；
- INCAR：基态部分 `EDIFF = 1e-8`、`ISMEAR = 0`、`SIGMA = 0.05`、`PREC = Accurate`、`KPAR = 4`；输运部分：
  - `ELPH_MODE = TRANSPORT`：布置输运计算；CRTA 不需要事先在超胞中计算电子–声子势；
  - `ELPH_SELFEN_CARRIER_DEN = -1e20 ... -1e13 1e13 ... 1e20`：额外载流子浓度列表（cm-3）。**注意符号约定**：正值加电子（n 型），负值加空穴（p 型），与文献惯例相反；
  - `ELPH_SELFEN_TEMPS = 0 100 200 300`：计算电子–声子散射自能的温度列表；
  - `ELPH_ISMEAR = -24`：电子–声子计算前确定费米能级/化学势的展宽方法——是四面体法更省内存与时间的实现（可试 `-15`，约慢 6–7 倍）；
  - `ELPH_SCATTERING_APPROX = CRTA` + `TRANSPORT_RELAXATION_TIME = 1e-14`：常数弛豫时间 1e-14 s；
  - `ELPH_TRANSPORT_DRIVER = 2`：用 Gauss–Legendre 积分做输运分布函数的能量积分（默认值）；
- KPOINTS 8x8x8（自洽基态）；KPOINTS_ELPH 24x24x24（输运）。

**流程与结果**：

1. 运行 `vasp_std` 后，用 py4vasp 遍历 `calculation.electron_phonon.transport`，从 `to_dict()` 取 `mobility[:,0,0]` 与 metadata 的 `selfen_carrier_den`，把各温度下迁移率对载流子密度作图。可见迁移率依赖载流子密度：**高密度下依赖是物理的，低密度下迁移率应与密度无关**——低密度区出现的"鼓包"正提示需要对 `ELPH_TRANSPORT_NEDOS` 做收敛；
2. 与实验对比：硅的实验迁移率约为电子 1400 cm2/Vs、空穴 500 cm2/Vs（Phys. Rev. B 97 121201(R)），而 CRTA 只给出约 99 与 72 cm2/Vs。这并不意外——弛豫时间被固定为 1e-14 s，调节该参数可以复现实验值。CRTA 对半导体尤其粗糙，因为态密度在带边行为尖锐。

**线性网格与 Simpson 积分**（`linear_grid/{51,501,5001,10001}`）：

- Onsager 系数 $L_{ij} = \int d\epsilon\, \mathcal{T}(\epsilon)\, (\epsilon-\mu)^{i+j-2} (-\partial f^0/\partial \epsilon)$（式 14.1）在线性能量网格上用 Simpson 法积分。新标签：
  - `ELPH_TRANSPORT_DRIVER = 1`：选择线性网格；
  - `ELPH_FERMI_NEDOS` / `ELPH_TRANSPORT_NEDOS`：线性网格点数（本例两者同取 51/501/5001/10001）；
  - `ELPH_TRANSPORT_DFERMI_TOL = 1e-6`：定义计算 Onsager 系数的能量窗口——窗口这样选取，使化学势附近 Fermi–Dirac 分布导数积分中被排除的部分不超过该容差；
- 用 CRTA 测试网格差异（计算便宜）；**CRTA 收敛的网格是 SERTA 应尝试的下限**。`linear_grid.sh` 循环运行四个 nedos；
- 结果：501、5001、10001 三条曲线重合，各载流子密度下能量网格已收敛。但线性网格方案的主要缺点是：能量范围由 `ELPH_TRANSPORT_DFERMI_TOL` 决定，**只加密点数不足以保证迁移率收敛，还必须对 DFERMI_TOL 做收敛**（双重收敛）。选做练习：把容差从 1e-6 改到 1e-8 重算——容差更严时低密度下迁移率不再依赖载流子密度，1e16 之前的鼓包正是数值收敛问题的信号。

**Gauss–Legendre 积分**（`gauss_legendre/{51,501,5001,10001}`）：

- 为避免上述双重收敛的不便，VASP 实现了 Gauss–Legendre 积分：`ELPH_TRANSPORT_DRIVER = 2`（默认）。INCAR 与线性网格版相同但去掉 DRIVER=1 与 DFERMI_TOL；
- OUTCAR 的电子–声子累加器输出行末会标明用了哪种网格（`Linear grids` 或 `Gauss-Legendre grid`）。例如同一累加器在 0 K：线性网格给出全零，Gauss–Legendre 网格给出 sigma = 115.9 S/m、mob = 72.35 cm2/(V s)；
- 结果：按能量点数看收敛比线性网格**慢**，但不需要调节第二个参数（能量范围），把待收敛参数从两个减为一个，收敛更简单。Gauss–Legendre 网格**自动适应 Fermi 窗口宽度**，数值高效，无需用 `ELPH_TRANSPORT_DFERMI_TOL`、`ELPH_TRANSPORT_EMIN/EMAX` 手工定义能量窗口，只需经 `TRANSPORT_NEDOS` 给定求和点数；
- 这一点对 SERTA 尤其重要：第一性原理寿命计算昂贵（每个 $\tau_{n\mathbf{k}}$ 都要对 q 求和），能量范围用来筛选哪些 $\tau_{n\mathbf{k}}$ 需要计算，可大幅减少计算量。结论：线性网格算得快但需要经验，Gauss–Legendre 更省心。

**k 点收敛**（`kpoints/{12,18,24,32,48,56}^3`）：

- INCAR 把载流子密度简化为 `ELPH_SELFEN_CARRIER_DEN = -1e16 1e16`（一个 p 型、一个 n 型），能量网格固定 5001（Gauss–Legendre）；KPOINTS_ELPH 依次取 12–56 的立方网格；
- 用 py4vasp 分别取电子（+1e16）与空穴（-1e16）在 300 K 的迁移率对网格作图：从 12x12x12 起电子与空穴迁移率都开始收敛，但真正收敛所需网格远比教程所取更密；
- 更便捷的做法是把迁移率对 **1/k** 作图并线性外推（取最后 3 点拟合）到"无穷 k 网格"极限，得到收敛值；
- 真实体系补充：完美晶体之外还有电离杂质（缺陷）影响散射率；缺陷密度高时迁移率下降。掺杂同时增加载流子与缺陷，超过某一点后收益递减——缺陷会反过来限制迁移率。

**复习问题**：(1) 用哪些 INCAR 标签定义 Gauss–Legendre 能量网格？(2) 如何把迁移率外推到无穷 k 采样？

### 示例 15：SERTA 计算迁移率（硅）

**教学目标**：(1) 在超胞上计算电子–声子势；(2) 在 SERTA 下让迁移率对 k、q 网格**同步**收敛；(3) 观察电子与空穴迁移率随 k 收敛的变化。

**任务**：计算 2x2x2 Si 超胞的电子–声子势，然后在 SERTA 下对金刚石结构 Si 原胞的迁移率做 k/q 网格收敛。

**原理**：CRTA 的电子寿命是常数；SERTA 改为第一性原理寿命，纳入电子–声子耦合，迁移率随之改变。第一步仍是准备超胞电子–声子势（与 Part 4 示例 12 相同），然后做收敛研究。

**输入要点与流程**：

1. **超胞电子–声子势**（`e15_mobility_si_serta/supercell`）：POSCAR 为常规胞的 2x2x2 超胞（立方 10.86 Angstrom，64 个 Si 原子）；INCAR：`EDIFF = 1e-8`、`ISMEAR = 0`、`SIGMA = 0.05`、`PREC = Accurate`，超胞部分 `IBRION = 6`、`ELPH_POT_GENERATE = TRUE`；KPOINTS 仅 Gamma 1x1x1。每个子目录用 KPOINTS_ELPH 指定不同 q 网格。超胞计算约 5 分钟，教程直接提供 `phelel_params.hdf5`；
2. **k/q 网格收敛**（`kpoints/{12,18,24,32,48}^3`）：POSCAR 换回 2 原子原胞；INCAR 与示例 14 类似（同样的载流子密度列表与 0/100/200/300 K 温度），改 `ELPH_SCATTERING_APPROX = SERTA`（不再需要 `TRANSPORT_RELAXATION_TIME`），能量网格 `ELPH_FERMI_NEDOS = ELPH_TRANSPORT_NEDOS = 5001`；KPOINTS 8x8x8；
3. 关键机制：SERTA 模式下电子–声子势从 `phelel_params.hdf5` 读出，并**插值**到连接 KPOINTS_ELPH 中所有 k 点的 q 网格上——因此修改 KPOINTS_ELPH 会**同时**调整 k 与 q 网格，可以系统研究迁移率对两者密度的收敛；
4. `qpoints.sh` 把 `phelel_params.hdf5` 拷入各目录并运行 12x12x12 与 18x18x18（其余网格每个 10 分钟以上，用预计算结果）。等待时可查阅 `ELPH_SCATTERING_APPROX` 支持的各种散射率近似；
5. **结果**：
   - 与示例 14 线性网格的曲线对比，数值困难已解决——低密度下迁移率不再依赖载流子密度；
   - 随 k 网格加密，电子迁移率持续上升直至收敛；
   - 把电子与空穴迁移率分开对网格作图：**电子迁移率远大于空穴迁移率**，源于有效质量与平均弛豫时间的差异（Phys. Rev. B 68, 125210 (2003)）；空穴比电子所需收敛动量更低；
   - 动量为 q 的声子散射电子，可散射的声子由 q 网格上限决定：q 网格越大，参与散射的声子模式越多，计算越贵。应在保证 q 网格收敛的前提下使用尽可能小的 q 网格。

**复习问题**：(1) 用哪个 INCAR 标签定义载流子密度？(2) 电子–声子计算中的 q 网格是哪种粒子的动量？（声子）

### 示例 16：热电材料与 ZT 优值（Bi2Te3）

**教学目标**：(1) 计算 Bi2Te3 含与不含自旋轨道耦合（SOC）的能带；(2) 在 CRTA 下计算 Bi2Te3 的各种输运系数；(3) 把态密度、电导率、热电势（Seebeck）、功率因子与 ZT 对化学势作图（300 K）；(4) 绘制 ZT 随载流子浓度与温度变化的热图。

**任务**：计算 Bi2Te3 原胞的能带（含与不含自旋极化/SOC），再在 CRTA 下计算其输运系数，把电导率、DOS、Seebeck 系数、功率因子与 ZT 对化学势作图（300 K），并产出 ZT 随载流子浓度与温度的热图。

**原理**：热电材料把热转化为电，其效率由几何与热电优值 $Z$ 决定，常与温度 $T$ 合写为 ZT：

$$ZT = \frac{\sigma S^2 T}{\kappa_e + \kappa_l}$$

其中 $\sigma$ 为电导率，$S$ 为 Seebeck 系数，$\kappa_e$ 为电子热导（以上三者 VASP 可直接计算），$\kappa_l$ 为晶格热导，可用 Müller–Plathe 方法或 phon3py 计算。典型热电材料 ZT 在 1–1.5 之间；ZT > 2 的新材料对废热回收、无运动部件的紧凑制冷器等能源技术意义重大。ZT 公式把多个输运系数汇聚在一起，它们都可追溯到 Onsager 系数与输运分布函数。改善 ZT 是经典优化问题（PNAS 93, 7436 (1996)）；本例仿照 Sofo 等人（Phys. Rev. B 68, 125210 (2003)）在 CRTA 下计算 Bi2Te3 的电子输运系数。

**输入要点与流程**（目录 `e16_thermoelectric_bi2te3/`）：

1. **能带：无/有 SOC**（`bands` 与 `bands_soc` 子目录）：POSCAR 为 Bi2Te3 原胞（2 Bi + 3 Te，六方/菱方层状结构，c 轴约 10.17 Angstrom）；INCAR 基态设置 `EDIFF = 1e-08`、`ISMEAR = -15`、`SIGMA = 0.01`、`LCHARG/LWAVE = .FALSE.`、`KPAR = 4`；`bands_soc` 版额外打开 `LSORBIT = .TRUE.`（自旋轨道耦合）并设 `MAGMOM = 15*0.0`（各原子初始磁矩为零）。POTCAR：`PAW Bi_d 09Feb1998` + `PAW Te 03Oct2001`（Bi 用含半芯 d 电子的赝势）。分别在两个目录运行：无 SOC 用 `vasp_std`，有 SOC 用 `vasp_ncl`（非共线）；
2. 用 py4vasp 把两套能带（`band.plot('kpoints_opt')`）叠加绘制（y 轴限 [-1, 1] eV）。py4vasp 的典型工作流：`Calculation.from_path` 加载 → 提取能带/DOS/输运量 → 内置绘图或导出分析。图因绘图 k 点有限而较粗糙，可增加点数改善。**SOC 的效果**：带隙改变；价带升到 0 处的"口袋"（pocket）从不含 SOC 时的 $\Gamma$ 点移到含 SOC 时的 $T$ 点附近——这改变了能谷简并度（multiplicity）与有效电荷，进而改变电导率。Bi2Te3 这类含重元素的材料必须考虑 SOC；
3. **输运系数对化学势**（`crta` 目录）：POSCAR 同上；INCAR 在 SOC 基态之上加电子–声子/输运标签：`ELPH_RUN = .TRUE.`、`ELPH_DRIVER = EL`、`ELPH_MODE = TRANSPORT`、`ELPH_ISMEAR = -24`、`EFERMI_NEDOS = 51`、`ELPH_SELFEN_TEMPS = 300`、`ELPH_SCATTERING_APPROX = CRTA`、`TRANSPORT_RELAXATION_TIME = 2.2e-14`，以及**新标签** `ELPH_SELFEN_MU_RANGE = -1.0 1.0 101`——把化学势设为相对费米能级的能量平移（起点、终点、点数；费米能级即 $T \to 0$、无额外载流子时化学势的极限）。KPOINTS 由 kspacing=0.3 生成（div=[6,6,1]）；KPOINTS_ELPH 由 kspacing=0.08 生成（div=[21,21,3]，面内密、c 轴疏，反映层状结构）。用 `vasp_ncl` 运行；
4. 后处理用 `py4vasp` 的 `raw.access("electron_phonon_transport")` 与 `calculation.electron_phonon.transport.from_data`，从每个 transport 的 metadata 读 `selfen_mu`，`to_dict()` 取 `electronic_conductivity`、`seebeck`、`electronic_thermal_conductivity`；DOS 直接从 `vaspout.h5` 的 `/results/electron_phonon/electrons/dos/` 读取。画五联图（300 K，xx 与 zz 两个张量分量）：电导率、DOS、Seebeck、功率因子 $\sigma S^2$、ZT（取晶格热导 $\kappa_l$ 为 xx 方向 1.5、zz 方向 0.7 W/(m K)）：

![Bi2Te3 输运系数对化学势（300 K）](/imgs/2026-08-05/fig53.png)

5. **图的解读**：所有量都由 Onsager 系数联系，图形与 Sofo 等人 Fig. 6 非常相似（加密 k 网格、改用实验晶格常数可进一步改善）。第一张图中电导率 $\sigma$ 实质上是被温度展宽的 DOS（第二张图）：离费米能级越远载流子越多、电导率越大。回忆 $\sigma = \mathcal{L}_{11}$、$S$ 与 $\kappa_e$ 由 $\mathcal{L}_{ij}$ 组合而成：若把电导率在带边附近近似为线性 $\sigma = a(\varepsilon-\mu) + b$（式 16.2），则 $\sigma S^2 \propto a^2/b^2$（式 16.3）——Seebeck 系数的峰与功率因子的峰并不重合。最后一张图中 ZT 有两个峰：化学势相对费米能为正与负各一个，分别对应电子与空穴载流子峰；
6. **化学势 ↔ 载流子密度**：图中化学势是相对费米能级的平移。实验上通过掺杂把化学势移进导带/价带以增加载流子。用 py4vasp 的 `electron_phonon.chemical_potential.to_dict()`（`chemical_potential`、`carrier_density`、`fermi_energy`；载流子密度 = (电子数 − 48) / 体积）作图：

![化学势与载流子密度关系](/imgs/2026-08-05/fig54.png)

载流子密度与化学势直接相关但**非线性**。前面所有图都可以改以载流子密度为横轴，只是 x 轴标度不同。由于二者不独立，计算中只能二选一：设电荷密度（`ELPH_SELFEN_CARRIER_DEN`）或化学势（`ELPH_SELFEN_MU`），按想测量的对象选择；

7. **ZT 热图：化学势 × 温度**（`crta_zt_map/0.2_mu`）：KPOINTS_ELPH 可用 `kmesh_poscar_spglib.py` 生成（`python3 ../../kmesh_poscar_spglib.py -kspacing 0.4 -o KPOINTS_ELPH`）。INCAR 与上一步的关键区别：不再用单温度 `ELPH_SELFEN_TEMPS = 300`，而用 **`ELPH_SELFEN_TEMPS_RANGE = 0 700 41`**（最低温、最高温、网格点数）；化学势用 `ELPH_SELFEN_MU_RANGE = -0.5 0.5 51`；另设 `elph_transport_nedos = 1001`、`TRANSPORT_RELAXATION_TIME = 1e-14`。KPOINTS kspacing=0.3（div=[6,6,1]），KPOINTS_ELPH kspacing=0.2（div=[9,9,2]）。用 `vasp_ncl` 运行后，以 `transport.figure_of_merit(kappa_lattice=1.0)` 计算 ZT（热图中 $\kappa_l$ 取 1），对化学势与温度画 imshow 热图：**max ZT ≈ 1.496**，300 K 对应图中红线（即上一节的曲线）：

![ZT 随化学势与温度热图](/imgs/2026-08-05/fig55.png)

8. **ZT 热图：载流子密度 × 温度**（`crta_zt_map/0.2_carrier`）：按化学势作图是理论计算的常用做法，但实验上直接控制的是载流子浓度（掺杂），按载流子密度作图便于与实验对比、指导材料设计。把 `ELPH_SELFEN_MU_RANGE` 换成 **`ELPH_SELFEN_CARRIER_DEN_RANGE = -1e22 -1e17 51 1e17 1e22 51`**：前一组（-1e22 到 -1e17，51 点）是空穴的对数网格，后一组（1e17 到 1e22，51 点）是电子的对数网格（也可用 `ELPH_SELFEN_CARRIER_DEN` 逐点定义）。其余设置同 `0.2_mu`（`EFERMI_NEDOS = 51`、`ELPH_SELFEN_TEMPS_RANGE = 0 700 41`、CRTA）。绘图时把密度取符号对数、按正负号拆成空穴/电子两支，画三并排热图（p 掺杂 / 合并 / n 掺杂）：**max ZT（p 型，空穴）≈ 1.505，max ZT（n 型，电子）≈ 0.771**：

![ZT 随载流子密度与温度热图](/imgs/2026-08-05/fig56.png)

9. **讨论**：换横坐标后输运性质的总体趋势与特征不变——ZT、Seebeck、电导率随 $\mu$ 的峰谷在载流子密度坐标下有对应结构，两种表示包含相同物理，但载流子密度图对材料优化与实验对比更实用；化学势与载流子密度的关系在带边附近（低密度）高度非线性，解读时需留意。选做：用 `-kspacing 0.2` 生成更密的 KPOINTS_ELPH 重算，观察图分辨率的变化。

**复习问题**：(1) 用什么 INCAR 标签定义化学势的范围？(2) 化学势–载流子密度关系如何影响 ZT 图的解读，尤其在带边附近？

### 功能与标签汇总（Electron-phonon Part 5）

| 功能 | 关键标签 | 说明 |
|---|---|---|
| 迁移率/输运 | `ELPH_MODE = TRANSPORT` | 输出 sigma、mob、seebeck、peltier、kappa_e |
| 散射近似 | `ELPH_SCATTERING_APPROX = CRTA / SERTA` | CRTA 免声子计算；SERTA 需 `phelel_params.hdf5` |
| 输运积分方案 | `ELPH_TRANSPORT_DRIVER = 1 / 2` | 线性网格+Simpson（需配 DFERMI_TOL）/ Gauss–Legendre（默认，自适应 Fermi 窗） |
| 能量网格 | `ELPH_TRANSPORT_NEDOS`、`ELPH_FERMI_NEDOS` | 积分网格点数；CRTA 收敛值是 SERTA 下限 |
| 能量窗 | `ELPH_TRANSPORT_DFERMI_TOL`、`ELPH_TRANSPORT_EMIN/EMAX` | 仅线性网格需要手动收敛 |
| 费米能级方法 | `ELPH_ISMEAR = -24` | 高效四面体法（比 -15 快约 6–7 倍） |
| 载流子设定 | `ELPH_SELFEN_CARRIER_DEN`、`ELPH_SELFEN_CARRIER_DEN_RANGE`、`ELPH_SELFEN_CARRIER_PER_CELL` | 正值为 n 型、负值为 p 型（与文献惯例相反）；RANGE 按对数网格生成 |
| 化学势设定 | `ELPH_SELFEN_MU`、`ELPH_SELFEN_MU_RANGE` | 相对费米能级的平移；与载流子密度二选一 |
| 温度范围 | `ELPH_SELFEN_TEMPS`、`ELPH_SELFEN_TEMPS_RANGE` | 单列表 / 起止+点数 |
| SOC 能带与输运 | `LSORBIT` + `vasp_ncl` | Bi2Te3 等重元素材料必须；SOC 改变带隙与能谷位置 |
| 晶格热导 | Müller–Plathe 方法 / phon3py | ZT 分母中的 $\kappa_l$ 需另行计算 |
| ZT 计算 | py4vasp `transport.figure_of_merit(kappa_lattice=...)` | 从 `vaspout.h5` 读取输运累加器 |

---

## 附录 A：教程输入文件下载清单

每个 Part 的完整输入文件（POSCAR/INCAR/KPOINTS/POTCAR 及脚本）可从官网下载：

| 分类 / Part | 下载链接 |
|---|---|
| 入门页（Get started） | [get_started.zip](https://vasp.at/tutorials/latest/get_started.zip) |
| 原子与分子 Part 1 | [molecules-part1.zip](https://vasp.at/tutorials/latest/molecules-part1.zip) |
| 原子与分子 Part 2 | [molecules-part2.zip](https://vasp.at/tutorials/latest/molecules-part2.zip) |
| 原子与分子 Part 3 | [molecules-part3.zip](https://vasp.at/tutorials/latest/molecules-part3.zip) |
| 块体 Part 1 | [bulk-part1.zip](https://vasp.at/tutorials/latest/bulk-part1.zip) |
| 块体 Part 2 | [bulk-part2.zip](https://vasp.at/tutorials/latest/bulk-part2.zip) |
| 块体 Part 3 | [bulk-part3.zip](https://vasp.at/tutorials/latest/bulk-part3.zip) |
| 块体 Part 4 | [bulk-part4.zip](https://vasp.at/tutorials/latest/bulk-part4.zip) |
| 磁性 Part 1 | [magnetism-part1.zip](https://vasp.at/tutorials/latest/magnetism-part1.zip) |
| 磁性 Part 2 | [magnetism-part2.zip](https://vasp.at/tutorials/latest/magnetism-part2.zip) |
| 分子动力学 Part 1 | [md-part1.zip](https://vasp.at/tutorials/latest/md-part1.zip) |
| 分子动力学 Part 2 | [md-part2.zip](https://vasp.at/tutorials/latest/md-part2.zip) |
| 分子动力学 Part 3 | [md-part3.zip](https://vasp.at/tutorials/latest/md-part3.zip) |
| 机器学习力场 Part 1 | [mlff-part1.zip](https://vasp.at/tutorials/latest/mlff-part1.zip) |
| 表面 Part 1 | [surface-part1.zip](https://vasp.at/tutorials/latest/surface-part1.zip) |
| 表面 Part 2 | [surface-part2.zip](https://vasp.at/tutorials/latest/surface-part2.zip) |
| 表面 Part 3 | [surface-part3.zip](https://vasp.at/tutorials/latest/surface-part3.zip) |
| 过渡态 Part 1 | [transition_states-part1.zip](https://vasp.at/tutorials/latest/transition_states-part1.zip) |
| 过渡态 Part 2 | [transition_states-part2.zip](https://vasp.at/tutorials/latest/transition_states-part2.zip) |
| 过渡态 Part 3 | [transition_states-part3.zip](https://vasp.at/tutorials/latest/transition_states-part3.zip) |
| 杂化泛函 Part 1 | [hybrids-part1.zip](https://vasp.at/tutorials/latest/hybrids-part1.zip) |
| 线性响应 Part 1 | [response-part1.zip](https://vasp.at/tutorials/latest/response-part1.zip) |
| GW Part 1 | [GW-part1.zip](https://vasp.at/tutorials/latest/GW-part1.zip) |
| BSE Part 1 | [BSE-part1.zip](https://vasp.at/tutorials/latest/BSE-part1.zip) |
| BSE Part 2 | [BSE-part2.zip](https://vasp.at/tutorials/latest/BSE-part2.zip) |
| BSE Part 3 | [BSE-part3.zip](https://vasp.at/tutorials/latest/BSE-part3.zip) |
| XAS Part 1 | [XAS-part1.zip](https://vasp.at/tutorials/latest/XAS-part1.zip) |
| XAS Part 2 | [XAS-part2.zip](https://vasp.at/tutorials/latest/XAS-part2.zip) |
| 强关联 Part 1 | [strongly_correlated-part1.zip](https://vasp.at/tutorials/latest/strongly_correlated-part1.zip) |
| 强关联 Part 2 | [strongly_correlated-part2.zip](https://vasp.at/tutorials/latest/strongly_correlated-part2.zip) |
| 强关联 Part 3 | [strongly_correlated-part3.zip](https://vasp.at/tutorials/latest/strongly_correlated-part3.zip) |
| NMR Part 1 | [nmr-part1.zip](https://vasp.at/tutorials/latest/nmr-part1.zip) |
| NMR Part 2 | [NMR-part2.zip](https://vasp.at/tutorials/latest/NMR-part2.zip) |
| NMR Part 3 | [NMR-part3.zip](https://vasp.at/tutorials/latest/NMR-part3.zip) |
| 声子 Part 1 | [phonon-part1.zip](https://vasp.at/tutorials/latest/phonon-part1.zip) |
| 声子 Part 2 | [phonon-part2.zip](https://vasp.at/tutorials/latest/phonon-part2.zip) |
| 电子–声子 Part 1 | [electron-phonon-part1.zip](https://vasp.at/tutorials/latest/electron-phonon-part1.zip) |
| 电子–声子 Part 2 | [electron-phonon-part2.zip](https://vasp.at/tutorials/latest/electron-phonon-part2.zip) |
| 电子–声子 Part 3 | [electron-phonon-part3.zip](https://vasp.at/tutorials/latest/electron-phonon-part3.zip) |
| 电子–声子 Part 4 | [electron-phonon-part4.zip](https://vasp.at/tutorials/latest/electron-phonon-part4.zip) |
| 电子–声子 Part 5 | [electron-phonon-part5.zip](https://vasp.at/tutorials/latest/electron-phonon-part5.zip) |

---

## 附录 B：全部教程页面索引

| 页面 | 标题 | URL |
|---|---|---|
| `index` | Tutorials | <https://vasp.at/tutorials/latest/> |
| `molecules` | Atoms and Molecules | <https://vasp.at/tutorials/latest/molecules/> |
| `molecules__part1` | Part 1: Introduction to VASP | <https://vasp.at/tutorials/latest/molecules/part1/> |
| `molecules__part2` | Part 2: Molecules in VASP | <https://vasp.at/tutorials/latest/molecules/part2/> |
| `molecules__part3` | Part 3: Water | <https://vasp.at/tutorials/latest/molecules/part3/> |
| `bulk` | Bulk systems | <https://vasp.at/tutorials/latest/bulk/> |
| `bulk__part1` | Part 1: Silicon as a typical bulk material | <https://vasp.at/tutorials/latest/bulk/part1/> |
| `bulk__part2` | Part 2: More silicon | <https://vasp.at/tutorials/latest/bulk/part2/> |
| `bulk__part3` | Part 3: Spin-polarized bulk systems | <https://vasp.at/tutorials/latest/bulk/part3/> |
| `bulk__part4` | Part 4: Van der Waals forces | <https://vasp.at/tutorials/latest/bulk/part4/> |
| `magnetism` | Magnetism | <https://vasp.at/tutorials/latest/magnetism/> |
| `magnetism__part1` | Part 1: Spin-polarized calculations | <https://vasp.at/tutorials/latest/magnetism/part1/> |
| `magnetism__part2` | Part 2: Energy differences comparing collinear magnetic structures | <https://vasp.at/tutorials/latest/magnetism/part2/> |
| `md` | Molecular dynamics | <https://vasp.at/tutorials/latest/md/> |
| `md__part1` | Part 1: Melting silicon | <https://vasp.at/tutorials/latest/md/part1/> |
| `md__part2` | Part 2: Machine learning force fields | <https://vasp.at/tutorials/latest/md/part2/> |
| `md__part3` | Part 3: Substitution reaction of chloromethane by a chloride ion | <https://vasp.at/tutorials/latest/md/part3/> |
| `mlff` | Machine-learned force fields | <https://vasp.at/tutorials/latest/mlff/> |
| `mlff__part1` | Part 1: Error analysis and hyperparameter optimization of machine-learned force fields | <https://vasp.at/tutorials/latest/mlff/part1/> |
| `surface` | Surfaces | <https://vasp.at/tutorials/latest/surface/> |
| `surface__part1` | Part 1: A nickel surface | <https://vasp.at/tutorials/latest/surface/part1/> |
| `surface__part2` | Part 2: A CO molecule adsorbed on a nickel surface | <https://vasp.at/tutorials/latest/surface/part2/> |
| `surface__part3` | Part 3: STM simulations | <https://vasp.at/tutorials/latest/surface/part3/> |
| `transition_states` | Transition states | <https://vasp.at/tutorials/latest/transition_states/> |
| `transition_states__part1` | Part 1: Basic transition states | <https://vasp.at/tutorials/latest/transition_states/part1/> |
| `transition_states__part2` | Part 2: Static approaches | <https://vasp.at/tutorials/latest/transition_states/part2/> |
| `transition_states__part3` | Part 3: Dynamic approaches | <https://vasp.at/tutorials/latest/transition_states/part3/> |
| `hybrids` | Hybrid functionals | <https://vasp.at/tutorials/latest/hybrids/> |
| `hybrids__part1` | Part 1: An overview of available functionals | <https://vasp.at/tutorials/latest/hybrids/part1/> |
| `response` | Linear response | <https://vasp.at/tutorials/latest/response/> |
| `response__part1` | Part 1: Static and frequency-dependent response | <https://vasp.at/tutorials/latest/response/part1/> |
| `gw` | GW approximation | <https://vasp.at/tutorials/latest/gw/> |
| `gw__part1` | Part 1: Introduction | <https://vasp.at/tutorials/latest/gw/part1/> |
| `bse` | Bethe-Salpeter equation | <https://vasp.at/tutorials/latest/bse/> |
| `bse__part1` | Part 1: Optical absorption of diamond carbon | <https://vasp.at/tutorials/latest/bse/part1/> |
| `bse__part2` | Part 2: Optical absorption of LiF | <https://vasp.at/tutorials/latest/bse/part2/> |
| `bse__part3` | Part 3: Efficient Brillouin zone sampling and analysis of the excitons | <https://vasp.at/tutorials/latest/bse/part3/> |
| `xas` | X-ray absorption spectroscopy | <https://vasp.at/tutorials/latest/xas/> |
| `xas__part1` | Part 1: X-ray absorption spectrum of LiCl | <https://vasp.at/tutorials/latest/xas/part1/> |
| `xas__part2` | Part 2: XAS K-edge of LiCl via Bethe-Salpeter equation | <https://vasp.at/tutorials/latest/xas/part2/> |
| `strongly_correlated` | Strongly correlated systems | <https://vasp.at/tutorials/latest/strongly_correlated/> |
| `strongly_correlated__part1` | Part 1: Constrained random-phase approximation | <https://vasp.at/tutorials/latest/strongly_correlated/part1/> |
| `strongly_correlated__part2` | Part 2: DFT+U and DFT+DMFT | <https://vasp.at/tutorials/latest/strongly_correlated/part2/> |
| `strongly_correlated__part3` | Part 3: Bethe-Salpeter equation calculations | <https://vasp.at/tutorials/latest/strongly_correlated/part3/> |
| `nmr` | Nuclear magnetic resonance | <https://vasp.at/tutorials/latest/nmr/> |
| `nmr__part1` | Part 1: NMR - chemical shielding | <https://vasp.at/tutorials/latest/nmr/part1/> |
| `nmr__part2` | Part 2: Coupling constants and two-center corrections | <https://vasp.at/tutorials/latest/nmr/part2/> |
| `nmr__part3` | Part 3: Aromaticity | <https://vasp.at/tutorials/latest/nmr/part3/> |
| `phonon` | Phonons | <https://vasp.at/tutorials/latest/phonon/> |
| `phonon__part1` | Part 1: Graphene | <https://vasp.at/tutorials/latest/phonon/part1/> |
| `phonon__part2` | Part 2: MgO | <https://vasp.at/tutorials/latest/phonon/part2/> |
| `electron-phonon` | Electron-phonon interactions | <https://vasp.at/tutorials/latest/electron-phonon/> |
| `electron-phonon__part1` | Part 1: Bandgap renormalization from perturbation theory | <https://vasp.at/tutorials/latest/electron-phonon/part1/> |
| `electron-phonon__part2` | Part 2: Electron-phonon interactions from statistical sampling | <https://vasp.at/tutorials/latest/electron-phonon/part2/> |
| `electron-phonon__part3` | Part 3: Electron-phonon matrix elements and using VASP with phelel | <https://vasp.at/tutorials/latest/electron-phonon/part3/> |
| `electron-phonon__part4` | Part 4: Conductivity of iron | <https://vasp.at/tutorials/latest/electron-phonon/part4/> |
| `electron-phonon__part5` | Part 5: Phonon-limited mobility of semiconductors and the ZT figure of merit | <https://vasp.at/tutorials/latest/electron-phonon/part5/> |

---

*本文档由 VASP 官方教程（vasp.at/tutorials/latest）整理而成，教程内容版权归 VASP 团队所有。*
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTI0NTcwNTAxNl19
-->