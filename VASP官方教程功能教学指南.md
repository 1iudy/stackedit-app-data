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

> 来源：https://vasp.at/tutorials/latest/molecules/part1/ ｜ 练习文件：`molecules-part1.zip`

本部分通过"大盒子中的单个氧原子"这一最简单体系，带领学习者完成第一次 VASP 计算，覆盖以下教学目标：

- 用"伪代码"层次理解 DFT 自洽循环：给定电荷密度 → 构造哈密顿量 → 求解本征函数/本征值 → 更新电荷密度 → 迭代至收敛；
- 学会为孤立原子创建四个输入文件（POSCAR/INCAR/KPOINTS/POTCAR）；
- 认识 stdout 与 OUTCAR 的基本结构，提取分子/原子的相关能量；
- 学会从上一次的 Kohn-Sham 轨道（WAVECAR）重启计算。

**示例 1：孤立氧原子的 DFT 计算**

- POSCAR：8 Å 立方大盒子（避免周期性镜像间相互作用）+ 1 个 O 原子（笛卡尔坐标）。
- INCAR：仅设 `SYSTEM` 与 `ISMEAR = 0`（高斯展宽）；其余用默认值即可发起一次 DFT 计算。教程强调 INCAR 标签都有默认值，未设置时按默认执行。
- KPOINTS：孤立原子只需 1 个 $\Gamma$ 点（Monkhorst-Pack `1 1 1`）。
- POTCAR：PAW_PBE O 赝势。教程讲解如何读懂 POTCAR 头部信息：价电子组态（O 为 s2p4）、交换关联（LEXCH = PE 即 PBE）、原子化能 EATOM、ZVAL = 6（6 个价电子）等。
- 运行：`mpirun -np 2 vasp_std`。教程逐列讲解 stdout 中 Davidson 迭代表格：`N`（迭代数）、`E`（总能）、`dE`（总能变化）、`d eps`（本征值变化）、`ncg`（迭代对角化步数）、`rms`（KS 轨道残差）、`rms(c)`（电荷密度残差）；收敛末行给出 `F`（自由能）、`E0`（收敛总能）、`d E`（熵项 $TS$ 引起的能量差）。
- 关键知识点：初始电荷密度取自 POTCAR，前 4 步保持不变（热身期），之后每轮更新；默认收敛判据为 dE < 1e-4 且 d eps < 1e-4。
- OUTCAR 讲解：版本与执行信息、INCAR/POSCAR/KPOINTS 读取、对称性分析、FFT 网格、赝势信息、本征值与费米能、最终总能等分区。孤立原子的总能接近 0，因为能量以赝势的原子参考能为基准。
- 重启计算：保留 WAVECAR 再运行可作为重启，从已收敛波函数出发显著减少迭代步数；全新计算应先删除 WAVECAR。

**示例 2：自旋极化氧原子（SDFT）**

- 教学目标：开启自旋极化（`ISPIN = 2`）；从两个自旋分量的本征值与占据数推断自旋磁矩；知道何时用 `vasp_gam`。
- 输入仅比示例 1 多一行 `ISPIN = 2`。理论要点：SDFT 中哈密顿量是 2x2 矩阵，KS 轨道为二分量的自旋向上/向下分量，两者可拥有不同本征值；两个自旋分量通过同时考虑双方电荷密度的有效哈密顿量分别求解、自洽耦合。
- 运行用 `vasp_gam`：当 KPOINTS 只含 $\Gamma$ 点时可用该可执行文件，数组按实数处理，比 `vasp_std`（复数）更省算力，但仅对纯 $\Gamma$ 点计算有效。
- stdout 中出现 WARNING：磁性/非共线计算未设置 `MAGMOM`，默认每原子 1 $\mu_B$ 的铁磁初设可能破坏对称性、排除反铁磁解；教程借此强调应手动设置初始磁矩并理解警告含义（警告不代表计算错误，而是提醒检查设置）。
- OUTCAR 中分别查看 spin component 1/2 的本征值与占据数：氧原子 2p 壳层在自旋向上占 3 个电子、自旋向下占 1 个电子，得到磁矩约 2 $\mu_B$（stdout 末行 `mag= 1.9998`）。
- 本示例还引导学习者找出费米能级位置，理解孤立原子中"带隙"与占据的关系。

**示例 3：低对称性下的自旋极化氧原子（对称性的作用）**

- 教学目标：理解 VASP 默认利用对称性（`ISYM`）降低计算量；比较立方（O_h）与正交（D_2h）盒子下的结果。
- 方法：把立方盒子改为三边不等长的正交盒子（7.0/7.5/8.0 Å），其余输入不变。OUTCAR 的对称性分析显示 VASP 找到 8 个空间群操作、点群 D_2h（示例 2 中为 O_h）。
- 结论：低对称性计算能量更低，更接近真实多体基态；GGA 交换关联泛函对多数原子倾向于对称性破缺的解。
- 重要教训：比较两个计算的能量时必须使用相同的 `SIGMA`。把 `SIGMA` 从默认 0.2 改为 0.01 会使收敛所需迭代显著增多，且 `energy without entropy` 数值明显变化——高斯展宽宽度影响能量绝对值，故能量比较必须在相同展宽设置下进行。
- 补充要点：该计算不含自旋轨道耦合（SOC）；ISPIN=2 只是非相对论哈密顿量的两个独立自旋分量，SOC 是相对论修正，需用 `LSORBIT` 开启。
- 课后问题示例：如何在 OUTCAR 中查看空间群操作数？如何构造四方（D_4h）对称的盒子？改变对称性的同时改变 SIGMA 能否比较能量？

**本部分涉及的 VASP 功能与标签汇总**：四输入文件体系（POSCAR/INCAR/KPOINTS/POTCAR）、PAW 赝势与 POTCAR 元信息、`ISMEAR`/`SIGMA`（展宽）、`ISPIN`（自旋极化）、`MAGMOM`（初始磁矩）、`ISYM`（对称性利用）、`SYSTEM`、`vasp_std`/`vasp_gam` 可执行文件、Davidson 电子迭代与收敛判据（EDIFF 相关默认值）、WAVECAR 重启、OUTCAR/stdout 输出解读。

---

## 1.2 Part 2：VASP 中的分子（O₂、CO 键长、振动与 DOS）

> 来源：https://vasp.at/tutorials/latest/molecules/part2/ ｜ 练习文件：`molecules-part2.zip`

**示例 4：O₂ 分子键长（共轭梯度几何弛豫）**

- 教学目标：用共轭梯度（CG）算法做几何弛豫求双聚体键长；以伪代码层次理解几何弛豫（Hellmann-Feynman 力与应力）与 CG 算法；学会设置 CG 步长与离子步数。
- 理论框架：Born-Oppenheimer 近似下电子与离子自由度解耦；几何弛豫即寻找力与应力为零的离子构型。伪代码：给定离子位置 → DFT 求电子基态 → 由 Hellmann-Feynman 定理计算力/应力 → 离子向瞬时基态弛豫 → 重复至收敛。
- 输入要点：POSCAR 为 8 Å 盒子中相距 1.22 Å 的两个 O 原子；INCAR 在 SDFT 设置基础上新增 `NSW = 5`（离子步数）与 `IBRION = 2`（CG 算法）；仍是 $\Gamma$ 点单 k 点。
- CG 算法伪代码：先做最速下降步（搜索方向=最大梯度方向，沿该方向用 Brent 法做线搜索直至力分量足够小，步长由 `POTIM` 控制）；再做共轭梯度步（把最大梯度对上一步搜索方向正交化得到共轭方向，再线搜索）；重复至梯度足够小或达到 `NSW`。CG 适合自由度少（≲4）且初猜接近基态的情形，否则应换用其它 IBRION 方法。
- 输出分析：stdout 中区分离子步与电子步；CONTCAR 保存弛豫后结构；用 py4vasp 的 `Calculation.from_path(...).structure` 查看/可视化结构。

**示例 5：CO 键长（多原子种类与赝势选择）**

- 教学目标：对多元素体系运行 VASP；选择合适的赝势；调节 CG 步长；用 py4vasp 可视化。
- 赝势知识：赝势是对离子实库仑势的近似，避免核处发散势与离子间缓变势的巨大反差；PAW 势按元素、价电子数目以及"软/硬"区分。教程强调：含短键双聚体（O₂、CO、N₂、F₂、P₂、S₂、Cl₂）时建议使用硬（hard）赝势。
- POSCAR 新要素：两种元素分别给出原子数（`1 1`）；引入 `sel`（选择性自由度）行，每个原子后跟 F/F/T 只允许沿 z 方向移动（对纯 DFT 静态计算该设置无实际影响，只在离子弛豫中起作用）。
- 动手环节：本例不提供 POTCAR，需要学习者用 `cat` 按 POSCAR 中元素顺序拼接单元素 POTCAR（顺序必须与 POSCAR 一致）。
- 调参练习：默认 `POTIM = 0.5` 步长过大导致弛豫不理想，删除 WAVECAR 后追加 `POTIM = 0.2` 重新计算，并用 py4vasp 逐步可视化结构演化（`my_calc.structure[:].plot()`）。

**示例 6：CO 分子振动频率（有限差分 Hessian）**

- 教学目标：用 VASP 计算分子振动频率；理解 Hessian 矩阵与声子频率的关系。
- 理论：分子振动即 $\Gamma$ 点声子。需求解 Born-Oppenheimer 能面 E(R) 的一、二阶导数；频率来自方程 $\det\left|\frac{1}{\sqrt{M_\mu M_\nu}}\frac{\partial^2 E}{\partial R_{\mu i}\partial R_{\nu j}} - \omega^2\right| = 0$，其中二阶导矩阵即 Hessian。
- 输入要点：`IBRION = 5`（有限差分求二阶导/Hessian/声子频率）、`NFREE = 2`（中心差分）、`POTIM = 0.02`（位移步长 0.02 Å）、`NSW = 1`；POSCAR 用选择性自由度只允许沿键轴（z）方向位移，从而只激发拉伸模式；使用 C、O 硬赝势。
- 结果解读：OUTCAR 末尾"Eigenvectors and eigenvalues of the dynamical matrix"部分给出本征频率；计算值 63.89 THz 与实验值 2143 cm$^{-1}$（64.25 THz）吻合良好。频率以多种单位（THz、cm$^{-1}$、meV）输出。
- 说明：IBRION=5 与 6 的区别在于是否利用对称性减少位移构型。

**示例 7：CO 的分波态密度（PDOS）**

- 教学目标：用 py4vasp 绘制 DOS；区分分波 DOS（partial）、局域 DOS（local）与总 DOS（total）。
- 概念：DOS 描述某能量区间内电子态数目 $N = \int f(\varepsilon)D(\varepsilon)\,d\varepsilon$；分波 DOS 是向特定轨道特征（如 O 2p）的投影；总 DOS 针对整个体系；晶体中波函数依赖波矢 k，局域 DOS 是对 $d^3k$ 积分的结果。有限温度下相互作用电子的 DOS 对应谱函数，可与光电子能谱（PES）实验对照。
- 输入要点：INCAR 设 `LORBIT = 11` 以产生投影信息（写入 DOSCAR/PROCAR）；其余与静态 CO 计算相同。
- 分析：py4vasp 中 `my_calc.dos.plot()` 绘总 DOS，`my_calc.dos.plot("p")` 按轨道字符绘分波 DOS。

**本部分涉及的 VASP 功能与标签汇总**：几何弛豫（`IBRION = 2` CG、`NSW`、`POTIM`、力收敛判据 `EDIFFG`）、CONTCAR 输出、多元素 POSCAR 与 POTCAR 拼接顺序、软/硬 PAW 赝势选择、选择性自由度（Selective dynamics）、有限差分声子（`IBRION = 5/6`、`NFREE`、`POTIM`）、Hessian 与动力学矩阵、态密度与分波态密度（`LORBIT = 11`、DOSCAR/PROCAR）、py4vasp 结构与 DOS 分析。

---

## 1.3 Part 3：水分子（RMM-DIIS、截断能收敛、DFPT 振动与从头算分子动力学）

> 来源：https://vasp.at/tutorials/latest/molecules/part3/ ｜ 练习文件：`molecules-part3.zip`

**示例 8：H₂O 键长（RMM-DIIS 几何弛豫）**

- 教学目标：以伪代码层次解释 RMM-DIIS（即广义最小残差法 GMRES）；判断何时用 RMM-DIIS、何时用 CG；从零编写 POSCAR；使用 POSCAR 缩放参数；完成双自由度几何弛豫；设置离子弛豫收敛判据。
- 算法要点：RMM-DIIS 是准牛顿算法，最小化对象是离子受力的模 $|F|^2$（而非总能）。DIIS 在若干步迭代构造的子空间内用直接法（RMM）求解；力与 Hessian 由 Hellmann-Feynman 定理给出（教程给出力 $F_{\mu i}$ 与 Hessian $H_{\mu i\nu j}$ 对电荷密度 $n(r;R)$ 的积分表达式）。RMM-DIIS 按 Hessian 谱把最小化分解到子空间，因此很快，但对初值敏感；适合离子位置接近极小、自由度不少于约 3 个的体系；自由度约 20 以上且振动谱很宽时建议改用分子动力学。
- 动手环节：学习者需自己编写 POSCAR——O 在原点、H-O 键长约 0.96 Å、键角 105°、足够大的盒子，并用选择性动力学允许每个 H 原子两个自由度。POSCAR 缩放参数取 0.52918（Bohr→Å 换算技巧）。POTCAR 需按 POSCAR 元素顺序拼接（O 在前 H 在后），此处用标准赝势而非硬赝势以降低数值成本。
- 结果：CONTCAR 给出弛豫后的键长与键角，可用 py4vasp 打印/可视化。

**示例 9：H₂O 的平面波截断能（ENCUT 收敛研究）**

- 教学目标：设置平面波基组的能量截断；用 `PREC` 设置计算精度。
- 理论：PAW 方法中 KS 轨道 = 离域的赝轨道 + 原子中心校正项；赝轨道用平面波展开 $\langle r|\tilde{\psi}\rangle = \frac{1}{\sqrt{\Omega_r}} \sum_G C^\sigma_{iGk} e^{i(G+k)\cdot r}$；无穷多平面波才构成完备基，实际用能量截断 $\frac{1}{2}|G+k|^2 < E_\mathrm{cutoff}$ 截去强振荡分量。截断能的选择主要受赝势影响，也依赖体系与目标性质，必须做收敛研究。
- 实践：教程提供一组子目录（ENCUT 从 400 到 820 eV），用 bash 循环逐个运行，再用脚本+gnuplot 绘制总能-截断能曲线。

![H₂O 总能随平面波截断能的收敛曲线](/imgs/2026-08-05/fig1.png)

- 结论：总能在 PAW 框架内向完备基组（CBS）极限收敛；平面波基始终不完备，由 ENCUT 定义。总能难以收敛，更好的做法是针对目标性质（振动频率、化学位移、带隙等）做截断能收敛研究；`PREC` 控制精度档位（影响 FFT 网格等数值设置）。

**示例 10：H₂O 振动频率——有限差分 vs 密度泛函微扰理论（DFPT）**

- 教学目标：解释有限差分法与 DFPT 的区别；设置 KS 轨道/电子能带数目（`NBANDS`）。
- 理论：两种方法都能计算二阶导（Hessian 与声子频率）。有限差分显式位移离子、用中心差分求力常数；DFPT 则用一阶微扰论求 KS 轨道的变化 $\Delta^R \psi$，得到 Sternheimer 方程 $(H - \varepsilon_{i\sigma})|\Delta^R \psi_{i\sigma}\rangle = -(\Delta^R V - \Delta^R \varepsilon_{i\sigma})|\psi_{i\sigma}\rangle$，无需显式位移离子。
- 输入要点：有限差分用 `IBRION = 6`（利用对称性）+ `NFREE` + `POTIM`；DFPT 用 `IBRION = 8`（教程以 `IBRION = 8`、`EDIFF = 1E-8` 作为参考结果）。
- 结果分析：OUTCAR 动力学矩阵部分给出 9 个模式——3 个高频内禀振动（对称/反对称伸缩 ~3744–3853 cm$^{-1}$、弯曲 ~1589 cm$^{-1}$）与 6 个接近零的平动/转动模式（数值上不为零，体现力的漂移）。教程演示：增大 ENCUT（520→800 eV）可减小力漂移，但物理有意义的频率对截断能不敏感；并让学习者改变 POTIM 考察有限差分对步长的敏感性。
- 对比结论：DFPT 无需多次位移计算、可直接处理金属与任意 q 点，但实现更复杂、某些泛函/设置受限；有限差分直观、适用广，但成本随原子数增长且依赖步长选择。

**示例 11：H₂O 对关联函数（从头算分子动力学）**

- 教学目标：关闭对称性（`ISYM = 0`）；运行标准从头算 MD（`IBRION = 0`）；解释 Nosé-Hoover 恒温器的基本概念。
- 概念：MD 用牛顿第二定律模拟给定温度下原子运动，每个离子步是一个时间步；力由量子力学（ab initio）计算即从头算 MD。缓慢降温的 MD 可实现模拟退火（如引入摩擦的阻尼分子动力学 IBRION=3），适合自由度约 20 以上且振动谱宽的体系。MD 的优势是轨迹具有物理意义，可统计计算对关联函数（在给定距离处找到粒子的概率）。
- 输入要点：`IBRION = 0`（MD）、`NSW = 500`、`POTIM = 0.5`（0.5 fs 时间步）、`ISYM = 0`（MD 中不施加对称性）、`SMASS = 0` + `TEBEG = TEEND = 2000`（Nosé-Hoover 恒温器控温 2000 K）、`NBANDS = 8`。
- 分析：用 py4vasp 逐步绘制总能量（`my_calc.energy[:].plot()`），绘制对关联函数；讨论恒温器如何把温度效应引入 MD。

**本部分涉及的 VASP 功能与标签汇总**：RMM-DIIS 离子弛豫（`IBRION = 1`）与 CG 的适用场景对比、POSCAR 缩放参数与选择性动力学、离子收敛判据（`EDIFFG`）、平面波截断能（`ENCUT`）收敛研究与 `PREC` 精度档位、PAW 赝轨道与基组完备性、有限差分声子（`IBRION = 5/6`）与 DFPT（`IBRION = 7/8`、Sternheimer 方程）、`NBANDS` 设置、从头算 MD（`IBRION = 0`、`POTIM`、`NSW`）、Nosé-Hoover 恒温器（`SMASS`、`TEBEG`/`TEEND`）、MD 中关闭对称性（`ISYM = 0`）、对关联函数分析（py4vasp）。

---

## 2. 块体体系（Bulk Systems）

块体体系是不受表面影响的物质部分的模型。本教程聚焦具有小原胞、周期性边界条件的晶体建模，学习密度泛函理论计算与结构优化。

## 2.1 Part 1：硅——典型块体材料（晶格常数、DOS、能带）

> 来源：https://vasp.at/tutorials/latest/bulk/part1/ ｜ 练习文件：`bulk-part1.zip`

**示例 1：fcc 硅的晶格常数（体积扫描 + 状态方程）**

- 教学目标：从晶格矩阵认出 fcc 结构；创建 $\Gamma$ 中心 k 网格；在不同原胞体积下手动执行多次 DFT 计算找总能极小，从而确定晶格常数；使用 POSCAR 的普适缩放因子。
- 晶体学背景：14 种布拉维格子 = 7 个晶系 x 心化类型；原胞是最小重复单元；实空间格子的傅里叶变换给出倒空间（k 空间），第一布里渊区是倒空间中原胞，不可约布里渊区再用点群对称性约化；高对称点/线因对称性强制的简并而特别重要（如 $\Gamma$ 点）。
- 输入要点：POSCAR 用 fcc 原胞（1 个原子），对称原胞基矢取 $a/2(\hat{x}+\hat{y})$、$a/2(\hat{y}+\hat{z})$、$a/2(\hat{z}+\hat{x})$，第二行缩放因子位置放占位符 `a`，把它解释为晶格常数；教程对比常规单胞（4 原子）与原胞（1 原子）。INCAR 显式给出 `ISTART = 0`（从头开始）、`ICHARG = 2`（原子电荷密度叠加初始化）、`ENCUT = 240`、`ISMEAR = 0`、`SIGMA = 0.1`；若不设 ENCUT 则取 POTCAR 中的默认值（ENMAX）。KPOINTS 用 Monkhorst-Pack 11x11x11，奇数网格保证 $\Gamma$ 中心。
- 工作流：bash 脚本循环把 `a` 替换为 3.5–4.3（步长 0.1），逐个运行 `vasp_std`，从 OSZICAR 提取自由能 F 写入数据文件，再用 gnuplot 绘制能量-晶格常数曲线求极小（教程配图展示抛物线型曲线与极小点）。
- 关键概念：OSZICAR 提供逐步摘要；总能极小对应平衡晶格常数；更严谨的做法是用状态方程（如 Birch-Murnaghan）拟合 E(V)。

![fcc 原胞示意图](/imgs/2026-08-05/fig2.png)

![不同晶格常数下的总能曲线](/imgs/2026-08-05/fig3.png)

**示例 2：fcc 硅的态密度（四面体方法与不可约 k 点）**

- 教学目标：理解 k 网格对称性约化；用四面体方法计算 DOS；提高 DOS 能量分辨率。
- 输入要点：`ISMEAR = -5`（带 Blöchl 修正的四面体方法）——教程强调它对部分占据不是变分的，会给出错误的力与应力，因此金属的几何优化不能用它；`LORBIT = 11` 使 DOSCAR 写出投影信息；k 网格加密到 15x15x15 以获得平滑 DOS（几何弛豫则不需要这么密）。
- 对称性约化：OUTCAR 显示 3375 个 k 点约化为 120 个不可约 k 点，每个带权重（等效 k 点数目），布里渊区求和用不可约点加权重完成；k 网格必须保持与倒易点阵相同的点群对称性。
- 提高分辨率：从 OUTCAR 找到费米能（9.91 eV），在 INCAR 加 `ICHARG = 11`（读 CHGCAR 固定密度做非自洽）、`NEDOS = 401`、`EMIN = -5`、`EMAX = 12` 重新运行；注意 EMIN/EMAX 是绝对能量，不是相对费米能。
- 分析：py4vasp 的 `mycalc.dos.plot()` 绘图（自动把费米能移到 0）。

**示例 3：fcc 硅的能带结构（高对称线 + 非自洽计算）**

- 教学目标：沿倒空间高对称线计算能带；用 py4vasp 绘图。
- 流程：能带计算需要 5 个输入——常规 4 个文件外加自洽计算的 CHGCAR。INCAR 设 `ICHARG = 11`（读 CHGCAR 并保持电荷密度固定）、`LORBIT = 11`；KPOINTS 用 line 模式（第 3 行首字母 `l` 开启线模式），在 reciprocal 坐标下逐段给出 L–$\Gamma$–X–U 与 K–$\Gamma$ 路径及每段插值点数，高对称点可附标签。
- 工具：用 SeeK-path（Materials Cloud 在线工具）从 POSCAR 确定空间群与推荐 k 路径，对照本例 KPOINTS。

![SeeK-path 给出的 fcc 布里渊区与高对称点](/imgs/2026-08-05/fig4.png)

- 分析：`mycalc.band.plot()` 绘制能带；讨论 fcc Si 的间接带隙特征。

**本部分涉及的 VASP 功能与标签汇总**：`ISTART`/`ICHARG`（启动模式与电荷初始化：2=原子叠加、11=读 CHGCAR 非自洽）、`ENCUT` 与 POTCAR 的 ENMAX 默认值、Monkhorst-Pack k 网格与 $\Gamma$ 中心、对称性约化与不可约 k 点权重、`ISMEAR = -5` 四面体方法（及其不能用于弛豫的原因）、`LORBIT`、DOS 参数（`NEDOS`/`EMIN`/`EMAX`）、OSZICAR/CHGCAR/DOSCAR 输出、KPOINTS line 模式能带计算、py4vasp 的 dos/band 分析。

---

## 2.2 Part 2：更多硅——金刚石结构、全弛豫与 β-Sn（收敛研究与作业控制）

> 来源：https://vasp.at/tutorials/latest/bulk/part2/ ｜ 练习文件：`bulk-part2.zip`

**示例 4：立方金刚石硅的晶格常数、DOS 与能带（自主练习）**

- 教学目标：以最少指导独立完成体积弛豫 + DOS + 能带全流程。
- 方法：对金刚石结构 Si（fcc 原胞 + 基元 2 原子，Direct 坐标 ±1/8）在 [5.4, 5.6] Å 区间二分式扫描晶格常数（每次算 3 个点、取最低两个缩小区间），最终得 a 约为 5.468 Å；备份最优结构的 CHGCAR。
- 重要告诫：算法收敛到的是能量景观中附近的局部极小——ab initio 结果可能对应亚稳态而非基态；正因如此不同晶体结构/磁构型可以互相比较，且应与实验对照验证。
- DOS 步骤：ENCUT 可降低（DOS 对截断能不如体积弛豫敏感）、`ISMEAR = -5`、`LORBIT = 11`、15x15x15 网格；教程解释 cd-Si 非金属、无部分占据带，故四面体方法可用于弛豫。推荐把 `vaspout.h5` 备份后用 `py4vasp.Calculation.from_file()` 读取。
- 能带步骤：`ICHARG = 11` + 线模式 KPOINTS + 复用 CHGCAR。
- 比较问题：fcc Si 与 cd Si 哪个更稳定？（比较总能，金刚石结构更低。）

**示例 5：一次运行弛豫体积、晶胞形状与离子位置（ISIF 与 Pulay 应力）**

- 教学目标：避免 Pulay 应力与 Pulay 力；单次 VASP 运行弛豫体积+形状+离子位置；设置 `ISIF`。
- 理论：平面波基不完备导致体积弛豫偏向更小晶胞，即 Pulay 应力——体积越小等效截断能越高，而总能随截断能增大而降低；不同体积的两个计算即使用相同 ENCUT 基组也不同。逐个体积重新计算（示例 4 的做法）可避免 Pulay 应力；单次运行优化体积则需把 ENCUT 设得足够高，并在体积显著变化后重启。
- 输入要点：`IBRION = 2` + `ISIF = 4`（弛豫晶胞形状+离子位置）或更高 ISIF 值同时变体积；教程要求逐行给 INCAR 标签写注释。
- 结果分析：CONTCAR 中最优晶格常数约为 5.469 Å；计算因 NSW=5 停止，需要 CONTCAR→POSCAR 续算；Pulay 应力表现为 OUTCAR 中负的 external pressure，本例因截断能足够高而可忽略。
- 稳健性建议：同时改变太多自由度会使过程不稳定；更稳健的是先手动体积扫描，再用 ISIF=4 + IBRION=2 弛豫形状与位置；ISIF 弛豫体积是在"恒定基组"而非恒定截断能下进行的，体积/形状变化大时要重启以更新基组。

**示例 6：受扰金刚石硅的离子位置弛豫（STOPCAR 与重启）**

- 教学目标：用 STOPCAR 文件中止计算；重启离子弛豫。
- 场景：HPC 上长时间作业时，希望优雅中止而非强杀，保留中间结果；重启也用于提高精度或在收敛结果上追加计算。
- 输入要点：POSCAR 中轻微扰动一个 Si 原子坐标；INCAR：`IBRION = 2`、`ISIF = 2`（只弛豫离子位置）、`NSW = 10`、`EDIFFG = -0.0001`（负值=按力收敛，单位 eV/Å）、`EDIFF = 1E-06`。
- 操作：第一次运行因 NSW=10 耗尽而停止（未收敛）；第二次把 NSW 增大到 15，在另一终端 `echo "LSTOP = .TRUE." > STOPCAR` 使 VASP 在完成当前离子步后退出；然后删除 STOPCAR、`cp CONTCAR POSCAR`、把 NSW 设为略大于剩余所需步数（如 8）续算，直至收敛并观察晶体对称性的恢复。
- 关键点：STOPCAR 支持 `LSTOP`（立即在离子步边界停止）与 `LABORT`（中止）；被 STOPCAR 停止后波函数/电荷密度仍可复用。

**示例 7：β-Sn 结构硅（k 网格收敛研究）**

- 教学目标：选择合适的 k 网格密度；为单元素计算创建 POTCAR。
- 任务：对 β-Sn Si（POSCAR 提供的体积不可信）做体积弛豫的 k 网格收敛研究，并计算 DOS。
- 流程：bash 循环不同 k 网格（如 6–14），每个网格用 `ISIF = 7`（固定形状只变体积）弛豫，且每个 k 点连续弛豫两次（CONTCAR→POSCAR）；从 CONTCAR 提取第一晶格矢量长度、从 OSZICAR 提取总能，作图分析收敛。
- 结论与细节：增大 k 点数并非单调趋近某值，但方差随 k 点增多而减小；$\Gamma$ 中心网格与否会引入显著差异；收敛后加 `LORBIT = 11` 计算 DOS；用 py4vasp 的 `structure.plot(2)` 可视化多个原胞。

**本部分涉及的 VASP 功能与标签汇总**：手动体积扫描与二分法确定晶格常数、局部极小/亚稳态的辨别、`ISIF` 全家族（2=只弛豫离子、4=形状+离子、7=体积、3=全弛豫等）、Pulay 应力机理与规避（高 ENCUT、逐体积重启、恒定基组）、`EDIFF`/`EDIFFG`（负值按力判据）、STOPCAR（`LSTOP`/`LABORT`）优雅中止、CONTCAR→POSCAR 续算流程、k 网格密度收敛研究、vaspout.h5 与 py4vasp。

---

## 2.3 Part 3：自旋极化块体——fcc Ni 与反铁磁 NiO（泛函对比）

> 来源：https://vasp.at/tutorials/latest/bulk/part3/ ｜ 练习文件：`bulk-part3.zip`

**示例 8：自旋极化 fcc 镍（共线 SDFT 全流程）**

- 教学目标：执行共线自旋密度泛函（SDFT）计算；弛豫共线磁体的体积；解读磁性体系的 DOS；用 py4vasp 绘 DOS 与能带。
- 理论：SDFT 中哈密顿量为 2x2 矩阵，KS 轨道是自旋空间二分量矢量，两分量可有不同的本征值；自旋分量间的耦合来自依赖两种自旋密度的交换关联泛函。`ISPIN = 2` 给出共线磁矩；若自旋方向重要（磁各向异性等），需加自旋轨道耦合 `LSORBIT = .TRUE.` 并使用 `vasp_ncl`。
- 输入要点：教程给出三套 INCAR/KPOINTS（体积弛豫、DOS、能带）：
  - 弛豫：`ISMEAR = 1`（Methfessel-Paxton，适合金属）、`SIGMA = 0.2`、`ENCUT = 540`、`IBRION = 2`、`ISIF = 7`（体积弛豫）、`NSW = 6`、`EDIFFG = -0.001`、`ISPIN = 2`、`MAGMOM = 1`；
  - DOS：`ENCUT = 350`、`ISMEAR = -5`、`LORBIT = 11`、16x16x16 网格；
  - 能带：`ICHARG = 11` + 线模式 KPOINTS + 备份的 CHGCAR。
- 操作细节：金属体积弛豫中每次最多 6 个离子步，防止体积变化过大引入 Pulay 应力偏差；反复 CONTCAR→POSCAR 续算，用 `vimdiff CONTCAR POSCAR` 观察变化发生在第几位小数；最优晶格常数约为 3.517 Å。
- 结果讨论：磁化使自旋向上/向下能带发生交换劈裂，DOS 在费米面附近呈现自旋不对称（Ni 为弱铁磁体，3d 带极化）；Methfessel-Paxton 展宽适合金属总能/弛豫，但不可用于绝缘体/分子（占据非物理）。

**示例 9：NiO 的带隙与磁矩（强关联体系与泛函对比）**

- 教学目标：用 `XC` 标签选择交换关联泛函；用 `BANDGAP` 标签控制带隙输出细节；设置与原子磁矩相关的 `MAGMOM`、`LORBIT`；用 py4vasp 提取带隙与原子磁矩。
- 背景：反铁磁 NiO 是强关联体系（Ni 3d 电子局域），标准 GGA（PBE）给出定性错误的电子与磁性质；需要 DFT+U、meta-GGA（SCAN、MBJ）、杂化泛函（HSE06）等。实验参考：带隙 4.0–4.3 eV，Ni 磁矩 1.9–2.2 $\mu_B$。
- 结构：NiO 为岩盐结构，但为容纳 AFM II 反铁磁序（沿立方 [111] 方向交替取向的铁磁面）采用菱方原胞；`MAGMOM = 2 -2 0 0` 设定两个 Ni 初始磁矩反平行。

![NiO 立方常规单胞与菱方原胞中的 AFM II 磁结构](/imgs/2026-08-05/fig5.png)

![NiO 菱方原胞磁结构](/imgs/2026-08-05/fig6.png)

- 共同设置：`LASPH = .TRUE.`（梯度的非球面贡献，对磁性与带隙精度重要）、`ISMEAR = -5`、`ENCUT = 350`、`EDIFF = 1E-4`、`ISPIN = 2`、`LORBIT = 11`（打印原子磁矩）、`BANDGAP = KPOINT`（打印更详细带隙信息）；KPOINTS 用 $\Gamma$ 中心 4x4x4。
- 分析手段：py4vasp 的 `bandgap.fundamental()` 与 `bandgap.direct()` 比较判断直接/间接带隙；OUTCAR 中 `grep "fundamental gap"` 查看每个自旋通道的带隙及 VBM/CBM 的 k 点位置（反铁磁体中两自旋通道带隙相同，因自旋上下能带简并）；`local_moment.projected_magnetic()` 提取原子磁矩（O 原子磁矩因对称为零；磁矩数值还依赖投影球半径等数值细节）。
- 泛函对比结果：
  - PBE（XC = PE）：带隙 0.96 eV、Ni 磁矩 1.40 $\mu_B$，显著偏离实验；
  - DFT+U（XC = CA 即 LDA + `LDAU = .TRUE.`、`LDAUTYPE = 2`、`LDAUL = 2 -1`（对 Ni 的 d 轨道加 U）、`LDAUU = 8.0 0.0`、`LDAUJ = 0.95 0.0`）：带隙 3.44 eV、磁矩 1.74 $\mu_B$；概念简单、成本与底层泛函相当，但 U 常需经验调节，削弱预测能力；
  - SCAN（XC = SCAN，meta-GGA，含动能能量密度、满足精确泛函的解析约束、无经验参数）：带隙 2.58 eV、磁矩 1.61 $\mu_B$；
  - MBJ（XC = MBJ，改进的 Becke-Johnson 势）：带隙 4.61 eV、磁矩 1.75 $\mu_B$；带隙通常很准，但额外势不对应总能贡献，因此不能用于结构优化或稳定性比较；
  - HSE06（XC = PE + `LHFCALC = .TRUE.` + `HFSCREEN = 0.2`，屏蔽杂化，25% HF 交换、屏蔽长程部分）：带隙 4.79 eV、磁矩 1.68 $\mu_B$；比 PBE0 对 k 网格收敛更快，但 HF 交换使计算贵 1–2 个数量级（用更多核运行）。

**本部分涉及的 VASP 功能与标签汇总**：共线 SDFT（`ISPIN = 2`、`MAGMOM`）、SOC 与非共线计算入口（`LSORBIT`、`vasp_ncl`）、Methfessel-Paxton 展宽（`ISMEAR = 1`）适用场景、磁性金属体积弛豫与 Pulay 应力控制、`XC` 标签切换泛函（PE/CA/SCAN/MBJ）、DFT+U 全套标签（`LDAU`/`LDAUTYPE`/`LDAUL`/`LDAUU`/`LDAUJ`）、杂化泛函标签（`LHFCALC`/`HFSCREEN`）、`LASPH`、`BANDGAP` 输出控制、py4vasp 的 bandgap 与 local_moment 分析。

---

## 2.4 Part 4：范德华力——石墨层间结合（Tkatchenko-Scheffler 方法）

> 来源：https://vasp.at/tutorials/latest/bulk/part4/ ｜ 练习文件：`bulk-part4.zip`

**示例 10：石墨层间结合能（TS 色散校正）**

- 教学目标：在 DFT 中引入范德华相互作用；设置 `LWAVE` 与 `LCHARG`；计算两个结构的总能之差。
- 背景：半局域 GGA 低估长程相互作用——PBE 预测石墨层间结合能仅 ~1 meV/atom，远小于 RPA 参考值 0.048 eV/atom（Lebègue et al., PRL 105, 196401）。Tkatchenko-Scheffler（TS）方法基于 Becke 的幂级数 Ansatz，显式引入带阻尼的原子对色散修正 $C_6 R^{-6}$（实现见 Bucko et al., PRB 87, 064110）。
- 输入要点：分别对石墨烯（2 原子、z 向 20 Å 真空、k 网格 15x15x1——层间无相互作用故 z 向只需 1 个 k 点）与石墨（4 原子、c = 6.71 Å、15x15x7）计算；INCAR：`IVDW = 20`（TS 方法）、`LVDW_EWALD = .TRUE.`（可选：对成对相互作用做 Ewald 求和）、`ALGO = Fast`、`PREC = Accurate`、`EDIFF = 1e-6`、`ISMEAR = -5`、`SIGMA = 0.01`；`LWAVE = .FALSE.`、`LCHARG = .FALSE.` 避免写出 WAVECAR/CHGCAR（节省磁盘与时间）。
- 结合能计算：E_bind = E(graphite)/4 - E(graphene)/2，从 OUTCAR 的 `free energy` 行提取总能；教程强调比较能量应按原子数归一而非按单胞。
- 进阶（IVDW = 202）：多体色散展开——把 INCAR 改为 `IVDW = 202` + `LVDWEXPANSION = .TRUE.` 可显示 2 体至 6 体对色散能的贡献；得到结合能 -0.051 eV/atom，与 RPA 参考 0.048 eV/atom 相当接近。
- 讨论：GGA 能否正确描述范德华作用？TS 如何引入色散？能否避免写出 CHGCAR/WAVECAR？

**示例 11：石墨层间距（Ewald 求和与一维扫描）**

- 教学目标：计算层状化合物的层间距；解释 Ewald 求和的基本概念。
- 背景：TS 中需要积分的长程库仑势可选用 Ewald 求和处理——短程部分在实空间、长程部分在倒空间求和，以避免奇异性。GGA-PBE 因低估色散会显著高估石墨堆叠方向晶格常数（8.84 Å vs 实验 6.71 Å）。
- 方法：POSCAR 中 c 方向用占位符 `c`，bash 循环在 [6.6, 6.8] Å 以 0.01 Å 步长扫描，每次运行提取 `free energy` 写入 loop.dat，gnuplot 作图求极小。
- 结果：TS 方法给出层间距 6.65 Å，与实验 6.71 Å 吻合良好（尽管示例 10 中结合能偏大）；配图展示能量-层间距曲线。

![石墨总能随层间距变化曲线](/imgs/2026-08-05/fig7.png)

- 思考题：如何修改 KPOINTS 完全忽略层间作用？Ewald 求和的动机与基本思想？

**本部分涉及的 VASP 功能与标签汇总**：范德华/色散校正家族（`IVDW`：20=TS、202=TS+多体展开等）、`LVDW_EWALD`、`LVDWEXPANSION`、Ewald 求和概念、`LWAVE`/`LCHARG` 控制输出文件、两结构能量差（按原子归一）、真空层建模与 k 网格维度匹配（表面/层状体系 z 向 1 个 k 点）。

---

## 3. 磁性（Magnetism）

磁性源于电子自旋、量子统计与电子–电子相互作用，是一种集体的量子电动力学现象。教程介绍如何用自旋极化建模铁磁/反铁磁材料、DFT+U、非共线计算与磁各向异性。

## 3.1 Part 1：自旋极化计算（hcp Co、bcc Cr、NiO 的 eg/t2g、DFT+U 与赝势选择）

> 来源：https://vasp.at/tutorials/latest/magnetism/part1/ ｜ 练习文件：`magnetism-part1.zip`

**示例 1：铁磁 hcp Co 的 DOS 与能带**

- 教学目标：执行自旋极化计算；获得总磁化强度与局域磁矩；绘制自旋极化体系的分波 DOS；自洽后第二步计算能带。
- 理论：`ISPIN = 2` 下 KS 轨道是自旋空间二分量矢量（全局自旋量子化轴），两分量可拥有不同本征值，可描述巡游电子的共线磁结构；自旋分量间仅通过交换关联泛函耦合。不加 SOC 时自旋在实空间的方向任意，只有相邻磁矩的相对取向有意义（可区分铁磁/反铁磁）。共线计算用 `vasp_std`；含 SOC 或非共线结构（局域量子化轴）才需 `vasp_ncl`。
- 输入要点：INCAR.scf：`ALGO = Normal`、`PREC = Accurate`、`EDIFF = 1e-5`、`ENCUT = 350`、`ISPIN = 2`、`MAGMOM = 2 2`、`LMAXMIX = 4`（d 电子体系需把 l 混合通道提高到 4 以保证磁信息正确混合）、`LASPH = T`、`ISMEAR = 2`（Methfessel-Paxton）、`SIGMA = 0.1`。k 网格各方向细分要与倒格矢长度匹配（实空间 c/a 约为 1.6，故 z 向细分约为面内的 1/1.6）。
- 流程与结果：stdout 末行 `mag=` 给出总磁化（除以 2 个 Co 原子后与实验 1.7 $\mu_B$/Co 接近）；DOS 步改用四面体方法（自洽用 MP 更可靠，DOS 用四面体更平滑），`LORBIT = 11`、`NEDOS` 设置 DOS 网格点数；py4vasp 绘 `dos.plot("up, down")` 与 `dos.plot("d")`——磁化主要来自 d 轨道，自旋上下 DOS 形状相似但能量平移，即 Stoner 模型中的交换劈裂。
- 局域磁矩：OUTCAR 中搜 `magnetization`，或 py4vasp `local_moment.to_dict()` 得每个原子 s/p/d 轨道磁矩（求和得单点总矩 ~1.53 $\mu_B$）；讨论"磁矩"与"有效磁矩"（$g\sqrt{s(s+1)}\,\mu_B$）的实验对应关系；无 SOC 时磁矩方向只能给出示意性全局方向。
- 能带：把 CHGCAR/POSCAR/POTCAR 拷入 band 子目录做非自洽计算；少数自旋能带因交换劈裂整体上移（如 M-K 间 -0.5 eV 的自旋向上带与 1.3 eV 的自旋向下带具有相同色散）。

**示例 2：反铁磁 bcc Cr 的磁单胞（k 收敛、KPOINTS_OPT 与重启）**

- 教学目标：对反铁磁材料做自旋极化计算；重启自旋极化计算。
- 体系：bcc Cr，AFM 序传播矢量 q=(1/2,1/2,1/2)，常规单胞即磁单胞；`MAGMOM = 2 -2`（或仅用 ±0.1 破缺对称）设定反平行初始磁矩，使两个 Cr 位点不再等价。
- 输入要点：`LMAXMIX = 4`、`LASPH = T`；`LKPOINTS_OPT = F` 先忽略 KPOINTS_OPT；细网格步 `ICHARG = 1` 从 CHGCAR 读自旋上下电荷密度（此时 MAGMOM 只用于破缺对称而不初始化磁矩）；DOS 步用 `ALGO = None`（关闭电子最小化的纯后处理，需要已收敛的 WAVECAR/vaspwave.h5）、`ICHARG = 11`、`NEDOS = 5001`、`LH5 = T`。
- k 收敛教学：先 4x4x4 再 8x8x8；重启细网格时要删除旧 WAVECAR（k 点数不匹配读取会失败，且失败尝试会使前几步不再冻结电荷密度，丧失粗网格已收敛的优势）；对比 OUTCAR 的 LOOP 数与 `total charge`/`magnetization (x)`——磁矩在粗网格未收敛而电荷几乎不变，强调必须对目标物理量做收敛研究（4x4x4 得 ±1.64 $\mu_B$，8x8x8 得 ±1.87 $\mu_B$）。
- 能带与 DOS：`LKPOINTS_OPT = .TRUE.`（默认）时用 KPOINTS_OPT（线模式 $\Gamma$–X–M–$\Gamma$）计算能带，py4vasp `band.plot("kpoints_opt")`；DOS 显示 AFM 无交换劈裂但多数自旋 DOS 更大，位点 2 与位点 1 的多数/少数自旋互换；`density.to_contour("magnetization") + density.to_quiver()` 在平面波网格上可视化磁化强度密度（磁化密度是自旋极化 DFT 的直接副产品，不依赖投影）。

**示例 3：反铁磁 NiO 中 e_g 与 t_2g 轨道的 DOS（密度混合与 MIDGAP）**

- 教学目标：分析晶体场劈裂后的轨道分辨 DOS；掌握磁性绝缘体的收敛技巧。
- 设置要点：教程引导阅读密度混合文档——`IMIX = 4` 为 Pulay 混合；电荷与磁自由度各有混合系数（`AMIX`/`BMIX` 与 `AMIX_MAG`/`BMIX_MAG`）；磁性计算可改用线性混合（把 BMIX、BMIX_MAG 设很小但非零）；`EFERMI = MIDGAP` 把费米能设在带隙中央。
- 分析：py4vasp 读取 dxy/dyz/dz2/dxz/dx2y2 的投影 DOS，组合成 t2g 与 eg 不可约表示；t2g 全占据，eg 在每个 Ni 位点上多数自旋近全占、少数自旋近全空，少数自旋在两个 Ni 位点间交替——AFM II 序。PBE 带隙仅 ~1.0 eV，远小于实验 4.2 eV（Dudarev et al., PRB 57, 1505），源于 DFT 对交换相互作用的低估；Ni 磁矩 1.34 $\mu_B$ 也低于实验 ~1.9 $\mu_B$。
- 概念辨析：VASP 直接计算的是 m 量子数分辨的投影 DOS，t2g/eg 是后处理的线性组合。

**示例 4：NiO 的 DFT+U 与 PAW 赝势选择（状态方程拟合）**

- 教学目标：用状态方程拟合做体积弛豫；选择并对比不同 PAW 势；按 Dudarev 方法运行 LSDA+U。
- 背景：DFT+U 通过在 LDA/GGA 上添加在位库仑/交换项改善强关联电子（过渡金属氧化物、f 电子体系）；LDA/GGA 低估 Hubbard U 导致 Mott 绝缘体被误判为金属。Dudarev 方案中有效参数为 U-J。
- 设计：两种 PAW 势（PAW_Ni_sv_GWO_s：Ni 的 s、p 半芯态作价电子，GWO 精度；普通 PAW_NiO）x 三个 U 值（0、3、8 eV，J = 0.95 eV，`LDAUTYPE = 2`、`LDAUL = 2 -1`）x 五个晶格常数缩放（0.95–1.05），共 30 组体积扫描计算；`XC = CA`（LDA）。
- 分析：用 ASE 的 Birch-Murnaghan 状态方程拟合 E(V) 得平衡晶格常数；py4vasp 提取带隙（`bandgap.print()` 显示 fundamental gap、direct gap、VBM/CBM 的 k 点位置与费米能）与局域磁矩；用 pandas 汇总作图。

![平衡晶格常数随 U 值变化（两种 PAW 势对比）](/imgs/2026-08-05/fig8.png)

![Ni 在位磁矩随 U 值变化](/imgs/2026-08-05/fig9.png)

![基带隙随 U 值变化](/imgs/2026-08-05/fig10.png)

- 物理结论：Ni_sv_GW 势中 sp 轨道在价带、与 d 轨道杂化——U 增大把 d 态推离费米能并收缩，纠缠的 sp 轨道膨胀，晶格常数增大；普通 Ni 势把 sp 冻结在芯里无法膨胀，晶格常数反而因 d 收缩略降。磁矩源于 d 轨道，对 PAW 势选择不敏感但随 U 显著变化；带隙对 PAW 势几乎不敏感——选赝势要综合考虑方法（如 DFT+U）、目标物理量与计算成本；还要注意目标量对体积的依赖（平衡晶格常数本身依赖泛函）。

**本部分涉及的 VASP 功能与标签汇总**：`ISPIN`/`MAGMOM`/`LMAXMIX`/`LASPH`、`vasp_std` 与 `vasp_ncl` 的适用边界、`ISMEAR` 各模式（MP=1/2、四面体=-5）的分工、`LORBIT`/`NEDOS` 投影 DOS、KPOINTS_OPT 与 `LKPOINTS_OPT`（自洽网格+能带路径二合一）、`ALGO = None` 后处理模式、`ICHARG = 1/11` 重启方式、WAVECAR 与 k 网格不匹配的陷阱、密度混合参数（`IMIX`/`AMIX`/`BMIX`/`AMIX_MAG`/`BMIX_MAG`）、`EFERMI = MIDGAP`、DFT+U（`LDAU`/`LDAUTYPE`/`LDAUL`/`LDAUU`/`LDAUJ`）、`BANDGAP` 输出、PAW 势选择（价电子数、GWO 精度）、Birch-Murnaghan 状态方程拟合、py4vasp 的 local_moment/density/bandgap 可视化。

---

## 3.2 Part 2：共线磁结构能量比较（Heisenberg 参数、SOC 与磁各向异性）

> 来源：https://vasp.at/tutorials/latest/magnetism/part2/ ｜ 练习文件：`magnetism-part2.zip`

**示例 5：用 DFT+U 为 NiO 建立 Heisenberg 模型**

- 教学目标：从第一性原理总能提取 Heisenberg 模型的最近邻 $J_1$ 与次近邻 $J_2$ 交换耦合参数；学会按单胞大小归一能量与匹配 k 网格。
- 方法：Heisenberg 哈密顿量 $H = -J_1\sum_{\langle ij\rangle} e_i\cdot e_j - J_2\sum_{\langle\langle ij\rangle\rangle} e_i\cdot e_j$；对 FM、AF1（(001)/(110) 面有序）、AF2（(111) 面有序）三种磁序，每化学式单元能量满足 $E_\mathrm{FM} = E_0-6J_1-3J_2$、$E_\mathrm{AF1} = E_0+2J_1-3J_2$、$E_\mathrm{AF2} = E_0+3J_2$，反解得 $J_1 = (E_\mathrm{AF1}-E_\mathrm{FM})/8$、$J_2 = (4E_\mathrm{AF2}-3E_\mathrm{AF1}-E_\mathrm{FM})/24$（参考 Archer et al., PRB 84, 115114）。模型参数可用于求解临界温度等性质。
- 输入要点：三种磁构型用不同超胞/原胞与 `MAGMOM` 初猜（如 AFM1 用 4Ni+4O 常规单胞 `MAGMOM = 2*-5.0 2*5.0 4*0`；AFM2 用菱方原胞 2Ni+2O；`4*5` 是 `5 5 5 5` 的简写）；共同设置：LDA+U（XC = CA、U=8、J=0.95、`LDAUPRINT = 1`）、`ISMEAR = -5`（精确总能推荐四面体方法）、`EDIFF = 1E-7`（收紧）、`EFERMI = MIDGAP`、线性化密度混合（AMIX/BMIX 系列）、`LMAXMIX = 4`、`LASPH = T`、`BANDGAP = WEIGHT`。
- k 网格：细分与晶格矢量长度成反比以保持可比密度（AFM2 用 5x5x5，FM/AFM1 用 4x4x4）；但精确总能更重要的是各自对 k 密度收敛——若一构型金属、另一构型有带隙，金属体系需要更细网格。
- 结果：能量按化学式单元归一（/4 或 /2）后计算得 $J_1 \approx 1.43$ meV、$J_2 \approx -18.2$ meV；MAGMOM 决定对称操作，因此不同子目录收敛到不同磁构型（实际中即使不加对称约束，终态也常接近初猜——这既是便利也说明磁基态判定本身很难）。

**示例 6：bcc Fe 中的自旋轨道耦合（SOC 打破空间群对称性）**

- 教学目标：开启 SOC；判断何时用 `vasp_ncl` 而非 `vasp_std`。
- 理论：SOC 需使用非共线磁性的自旋密度泛函框架（Hobbs et al., PRB 62, 11556），体系用 2x2 自旋密度矩阵 n(r) 描述；SOC 项采用零阶正则近似（ZORA）：$H_\mathrm{SO} = \frac{1}{(2c)^2}\frac{K^2(r)}{r}\frac{dV}{dr}\,\sigma\cdot L$，$K(r) = 1/(1-V(r)/2c^2)$。SOC 能标很小但后果显著：无 SOC 时 bcc Fe 的 $\Gamma$–$H_3$ 与 $\Gamma$–H 高对称路径等价，加入 SOC 后出现典型的能带反交叉（anticrossing）。
- 输入要点：`LNONCOLLINEAR = T` + `LSORBIT = T`；`MAGMOM = 0.0 0.0 2.0`（三分量矢量，指定磁矩沿 z）；`EDIFF = 1e-8`；混合参数 AMIX/BMIX 系列（BMIX 设很小但非零以免某些版本崩溃）；`LMAXMIX = 4`、`LASPH = T`、`LORBIT = 11`；自洽 12x12x12 $\Gamma$ 网格 + KPOINTS_OPT 线模式（$H_3$–$\Gamma$–H，141 个交点）；用 `vasp_ncl` 运行。
- 结果：能带图中可见 SOC 引起的劈裂与反交叉，教程配图展示。

![bcc Fe 含 SOC 的能带劈裂示意](/imgs/2026-08-05/fig11.png)

**示例 7：FeO 的磁各向异性能（MAE）**

- 教学目标：完整的 MAE 计算工作流；SAXIS 旋转路径的构造；结合 ASE/Spglib 的结构处理。
- 背景：磁各向异性是 SOC 引起的相对论效应。FeO 为轻微单斜畸变的岩盐结构（空间群 C2/m），基态磁结构为 (111) 面 AFM 有序；教程复现 Schrön, Rödl & Bechstedt, PRB 86, 115134 的 Fig. 5(a)。
- 标准 MAE 工作流：① 用 `vasp_std` 做共线 SCF（`ISPIN = 2`、`ISYM = -1` 关闭对称性）；② 用 `ICHARG = 11` 固定共线电荷密度重启非共线+SOC 计算（`vasp_ncl`），沿 MAE 路径旋转 `SAXIS`。优点：改变磁化方向时不必担心收敛到别的方向；前提：SOC 是小微扰、电荷密度与磁化强度几乎不受 SOC 影响，且在位轨道角动量由晶体场而非磁矩方向决定。
- 路径构造（$\theta=0$ 对应磁化平行立方 [111]，$\varphi=0$ 对应磁化在 (1-10) 面内）：用 ASE 读取 FeO 岩盐 CIF → 构造 2x2x2 超胞 → 施加 $q=(1/2,1/2,1/2)$ 自旋波（用 Fe/Mn 占位标记自旋上下）与单斜畸变（Schrön 表 II 的应变参数 r=0.0014, e=0.0253, t=0.0251）→ 用 Spglib 确认空间群 C2/m → 求原初磁单胞 → 在畸变超胞中识别立方 [111] 方向与 (1-10) 面（旋转轴为 [-110]）→ 计算 $\theta\in[0,\pi]$、$\varphi=0$ 与 $\varphi=90^\circ$ 的 SAXIS 列表。
- 结果：沿路径逐点计算总能，得到 MAE 随角度的变化并复现文献曲线；讨论磁各向异性能的量级（$\mu$eV 级）与收敛要求（极紧 EDIFF、密集 k 网格、力定理或 torque 方法等进阶话题）。

![FeO 2x2x2 超胞](/imgs/2026-08-05/fig12.png)

![施加单斜畸变后的超胞](/imgs/2026-08-05/fig13.png)

![带自旋波标记的超胞](/imgs/2026-08-05/fig14.png)

![原初磁单胞](/imgs/2026-08-05/fig15.png)

**本部分涉及的 VASP 功能与标签汇总**：多磁构型总能比较与 Heisenberg 参数提取、`MAGMOM` 简写语法与对称性作用、`ISMEAR = -5` 精确总能、`EDIFF` 收紧、`EFERMI = MIDGAP`、`LDAUPRINT`、SOC 与非共线（`LSORBIT`/`LNONCOLLINEAR`/`SAXIS`/`MAGMOM` 三分量）、ZORA 形式的 SOC、`ISYM = -1`、MAE 两步工作流（共线 SCF → ICHARG=11 非共线 SOC）、ASE+Spglib 结构处理流水线。

---

## 4. 分子动力学（Molecular Dynamics）

分子动力学（MD）按每个时间步作用在各粒子上的力模拟原子（与分子）的运动。学习如何在 VASP 中进行 MD 模拟。

## 4.1 Part 1：硅的熔化（从头算分子动力学基础）

> 来源：https://vasp.at/tutorials/latest/md/part1/ ｜ 练习文件：`md-part1.zip`

**示例 1：固态金刚石硅的 NVT 从头算 MD**

- 教学目标：说明什么是从头算分子动力学；用 pymatgen 从 CIF 创建超胞；区分 `vasp_gam` 与 `vasp_std`；用 py4vasp 绘制 MD 中各类能量随时间的演化与轨迹。
- 概念：MD 用经典运动方程模拟给定温度下的原子运动，每个离子步是一个时间步；力由量子力学计算即"从头算 MD"。正则系综（NVT）要求粒子数、体积、温度恒定；温度效应由恒温器引入——通过随机项或确定性附加动力学变量模拟热浴，Nosé-Hoover 属于确定性方案（在小/刚性体系中可能缺乏遍历性，但对大体系适用）。
- 结构准备：用 pymatgen 读取 Materials Project 的 cd-Si CIF，构造 2x2x2 常规单胞超胞（64 原子）写为 POSCAR——超胞用于捕获影响动力学的晶格振动，实际工作应对原胞大小做收敛检查。
- 输入要点：`ISMEAR = 0`、`SIGMA = 0.1`、`LREAL = Auto`（实空间投影算符，大体系提速）、`ALGO = VeryFast`（RMM-DIIS 电子弛豫）、`PREC = Low`（MD 可用低精度提速）、`ISYM = 0`；MD 部分：`IBRION = 0`、`NSW = 30`、`POTIM = 3.0`（时间步 3 fs）、`MDALGO = 2`（Nosé-Hoover）、`SMASS = 1.0`（Nosé 质量/热惯性）、`TEBEG = TEEND = 2000`、`ISIF = 2`（固定晶胞）；大超胞用 $\Gamma$ 点、`vasp_gam` 运行。
- 输出解读：stdout 每个离子步末行给出 `T`（瞬时温度）、`E`（总能=离子势能 F + Nosé 恒温器势能 SP + 动能 SK + 电子…）、`EK` 等；轨迹写入 XDATCAR；用 py4vasp `energy[:].plot()` 绘各能量分量，`structure[:].plot()` 动画轨迹。

**示例 2：硅的熔化（对关联函数 PCDAT 分析）**

- 任务：从示例 1 的末态结构重启，继续 90 fs 模拟并分析对关联函数 g(r)（径向分布函数）。
- 要点：对关联函数写入 PCDAT 文件；教程提供 awk+gnuplot 脚本解析绘图。90 fs 时体系仍接近晶体，长程序使 g(r) 在 4 Å 以外仍有明显峰；继续模拟到 180 fs 后峰形趋于液体特征（短程有序、长程平坦）。注意 g(r) 只在 r 小于超胞最短边长一半时有意义（此处 < 5.4 Å）。
- 轨迹表示：XDATCAR 记录真实穿越边界的轨迹；"中心盒表示"则对每步位置施加周期性折叠。
- 重启与调参：`cp CONTCAR POSCAR` 续算；`NPACO`/`APACO` 控制 g(r) 的最大距离与网格点数、`NBLOCK`/`KBLOCK` 控制写出频率；示例把 APACO 设为 5.4 Å 再跑 15 fs，并讨论新 PCDAT 只对本段采样、统计质量下降的问题。

![90 fs 与 180 fs 的对关联函数对比](/imgs/2026-08-05/fig16.png)

![XDATCAR 轨迹与中心盒表示示意](/imgs/2026-08-05/fig17.png)

![90/180/195 fs 对关联函数对比](/imgs/2026-08-05/fig18.png)

**示例 3：监控分子几何（恒温器对比与 ICONST）**

- 教学目标：了解 VASP 中各恒温器（`MDALGO`）与系综的组合；用 ICONST 文件监控 MD 中的内坐标（键长、键角等）。
- 恒温器家族：`MDALGO = 1` Andersen（完全随机引入温度）、`MDALGO = 2` Nosé-Hoover（确定性）、`MDALGO = 3` Langevin（随机+确定性）。对分子内强作用力（键伸缩）体系，Andersen 与 Langevin 无遍历性问题且平衡化高效；Nosé-Hoover 在小/刚性体系缺遍历性（如单丁烷分子模拟）。实际选择主要取决于 VASP 中恒温器与热力学系综的可用组合（参见 Wiki 的 Ensembles 分类）。
- ICONST：VASP 用 ICONST 文件指定内部几何坐标（键长、角度、二面角），可在 MD 中监控这些量的演化（也可用于约束/伞形采样）。
- 实践：构造 cd-Si 原胞的 2x2x2 超胞（16 原子），本例 INCAR 还演示了 `IVDW = 10`（DFT-D2 色散校正）、`ISMEAR = -1`（Fermi 展宽）+ 小 SIGMA 的设置，运行 MD 并用 py4vasp 监控几何量随时间的变化。

**本部分涉及的 VASP 功能与标签汇总**：`IBRION = 0` MD、`MDALGO`（1=Andersen、2=Nosé-Hoover、3=Langevin）、`SMASS`、`TEBEG`/`TEEND`、`POTIM`（fs 时间步）、`NSW`、NVT 系综实现、`LREAL = Auto`、`ALGO = VeryFast`、`PREC = Low`、`ISYM = 0`、`vasp_gam`、XDATCAR/PCDAT 输出、`NPACO`/`APACO`/`NBLOCK`/`KBLOCK`、CONTCAR→POSCAR 重启 MD、ICONST 内坐标监控、pymatgen 超胞构建、py4vasp 轨迹与能量分析。

---

## 4.2 Part 2：机器学习力场（在线训练、验证、系综平均与可迁移性）

> 来源：https://vasp.at/tutorials/latest/md/part2/ ｜ 练习文件：`md-part2.zip`

**示例 4：在从头算 MD 中在线训练力场（on-the-fly MLFF）**

- 教学目标：在 ab-initio MD 过程中在线（on-the-fly）训练机器学习力场。
- 动机：良好的热力学平均需要长模拟时间，而每个 ab-initio MD 步昂贵；经典力场便宜但依赖经验参数；MLFF 从第一性原理数据学习物理，兼顾精度与长时间模拟。
- 输入要点（16 原子 Si，NpT + Langevin）：`ML_LMLFF = T` 开启 MLFF、`ML_ISTART = 0` 表示从 ab-initio MD 训练新力场；MD 设置：`IBRION = 0`、`NSW = 10000`、`POTIM = 2.0`、`MDALGO = 3`（Langevin）、`LANGEVIN_GAMMA = 1`（原子摩擦）、`LANGEVIN_GAMMA_L = 10`（晶格摩擦）、`PMASS = 10`（晶格质量，NpT 所需）、`TEBEG = 400`、`ISIF = 3`（位置+形状+体积都更新）；`ICONST` 用 LR/LA/LV 监控晶格矢量长度、夹角与体积；`RANDOM_SEED` 保证可重复性（训练示例尽量短，学到的局域参考构型数强烈依赖随机种子）；训练用 2x2x2 k 网格——力必须对 k 点收敛，虽然力场可以在小单胞训练后用于大体系。
- 在线训练算法（三种触发）：① 无力场时，把结构的局域参考构型加入训练集并构建力场；② 任意 MD 步若任一原子的贝叶斯力误差超过严格阈值 `ML_CDOUB x ML_CTIFOR`，加入训练集并重建；③ 构建后进入 `ML_NMDINT` 步热身期（对误差宽松），之后若误差超过低阈值 `ML_CTIFOR`，结构进入待处理列表，列表达到 `ML_MCONF_NEW` 长度后批量加入训练集。注意贝叶斯误差只是样本外误差的估计，并非力场质量的真实指标。
- 输出文件：训练数据写入 ML_ABN，最终力场参数存于 ML_FFN。

**示例 5：用离子弛豫测试力场（E-V 曲线对比）**

- 教学目标：通过对比弛豫晶格参数与参考数据测试力场；在力场下使用共轭梯度算法。
- 测试策略：① 独立测试集（随机构型上对比 DFT 与 MLFF 的力、应力、能量差）；② 对比物理性质（弛豫晶格参数、声子、相相对能量、弹性常数、缺陷形成能等）。本例采用后者。
- 输入要点：`ML_LMLFF = T` + `ML_ISTART = 2`（只读已训练力场、不做 DFT）；`IBRION = 2`、`ISIF = 3` 做单点/弛豫；ML_FFN 软链接为 ML_FF 读入。bash 循环 13 个体积缩放因子（0.96–1.03）做静态单点，提取 OUTCAR 的 volume 与 free energy 绘制 E-V 曲线，与相同 DFT 设置的参考数据对比——MLFF 曲线与 DFT 高度吻合。
- 概念：训练只能得到样本内误差；贝叶斯误差是样本外误差的估计而非真实度量。

**示例 6：系综平均的晶格常数与体积（纯力场生产运行）**

- 教学目标：把晶格常数/体积作为系综平均计算；用预训练力场运行 MD。
- 概念：晶格参数常被视为 T=0 K 基态性质，但许多材料随温度发生结构转变，且实验测量的就是热力学系综平均——用 MLFF 在 400 K（训练温度）做生产运行。
- 操作：INCAR 不再含任何 ab-initio 标签（无 ENCUT 等），仅 MD + `ML_LMLFF = T` + `ML_ISTART = 2`；从 400 K 热化构型出发跑 3x10000 步；对比耗时——同样模拟 ab-initio MD 需 >20 小时，MLFF 不到 10 分钟；REPORT 文件记录每步 ICONST 监控量，用脚本提取 LV/LR 并对时间取平均（得体积 ~319.3 Å$^3$、平均晶格矢量 ~7.673 Å）；三个晶格矢量长度的微小不对称源于模拟时间不够长。

**示例 7：力场可迁移性与热膨胀系数**

- 教学目标：检验力场对特定体系/参数集的适用性；独立从 REPORT 提取监控量；计算热膨胀系数。
- 可迁移性判据：纯力场 MD 中贝叶斯误差突增表明进入未知相空间区域（如相变），需要继续训练——把 ML_ABN 复制为 ML_AB 并设 `ML_ISTART = 1` 继续学习；但改变温度等参数不一定需要重训：只要训练时的局域参考构型已充分覆盖目标系综的相空间（无相变前提下）。验证方法：在训练参数与目标参数下分别用 MLFF 监控单胞体积分布，比较概率密度——教程图示表明 200 K 遇到的所有局域构型都属于 400 K 相空间的子集。
- 实践：在 200 K 与 300 K 各跑 3x10000 步 MLFF MD（沿用 400 K 训练的 ML_FF），提取系综平均体积，绘制 V(T) 曲线并拟合热膨胀系数。

![200 K 与 400 K 体积概率密度对比（可迁移性验证）](/imgs/2026-08-05/fig19.png)

**本部分涉及的 VASP 功能与标签汇总**：MLFF 全套标签（`ML_LMLFF`、`ML_ISTART`：0=在线训练、1=继续训练、2=仅使用、`ML_CTIFOR`、`ML_CDOUB`、`ML_NMDINT`、`ML_MCONF_NEW`）、ML_FFN/ML_AB/ML_ABN/ML_LOGFILE 文件体系、贝叶斯误差与样本内/外误差、NpT 系综（Langevin + `PMASS` + `LANGEVIN_GAMMA_L`）、`RANDOM_SEED`、ICONST 监控（LR/LA/LV）、REPORT 文件解析、力场验证方法学、热膨胀系数计算。

---

## 4.3 Part 3：氯甲烷的氯离子取代反应（SN2 自由能面与反应速率）

> 来源：https://vasp.at/tutorials/latest/md/part3/ ｜ 练习文件：`md-part3.zip`

**示例 8：慢增长模拟自由能面（热力学积分）**

- 教学目标：用热力学积分获得自由能面；执行慢增长（slow-growth）模拟并判断其质量；解读 ICONST 文件内容；从 REPORT 文件提取信息。
- 理论：自由能面给出反应的活化自由能（从反应物到过渡态所需可逆功）。热力学积分沿反应路径对自由能梯度 $\partial A/\partial\xi$ 积分：$\Delta A_{A\to B} = \int (\partial A/\partial\xi)|_{\xi'} d\xi'$。慢增长法在每个时间步以速度 $\dot{\xi}'$（`INCREM`）缓慢改变反应坐标，使过程近似可逆（功=自由能）；对比快切换模拟（多次快速转变 + Jarzynski 恒等式）与静态方法（PES 驻点的弛豫+振动分析，便宜但无法计入软模）。
- 体系：对称 SN2 反应 Cl$^-$ + CH$_3$Cl $\to$ [Cl$\cdots$CH$_3\cdots$Cl]$^-$ 中间态 $\to$ 产物，自由能面对称。
- 输入要点：反应坐标用 ICONST 定义——`R 1 5 0`、`R 1 6 0` 监控两个 C–Cl 键长，`S 1. -1. 0` 定义 $\xi = r(\mathrm{C-Cl_1}) - r(\mathrm{C-Cl_2})$（差值坐标）；INCAR：`IBRION = 0`、`MDALGO = 1`（Andersen 恒温器）+ `ANDERSEN_PROB = 0.05`、`TEBEG = 300`、`ISIF = 2`、`LBLUEOUT = T`（把自由能梯度写入 REPORT）、`INCREM = 8e-4`（慢增长转变速度）；MLFF：`ML_LMLFF = T` + `ML_MODE = RUN`（读 ML_AB 与 ML_FF、不学习），力场来自事先准备的 CH$_3$Cl/Cl$^-$ ab-initio 数据；`ML_MODE = refit` 可从现有 ML_AB 生成力场；`POMASS = 12.011 3.000 35.453` 把 H 换成氚（更重→氢振动更慢→可用更大 POTIM 时间步）；`KSPACING = 1000` 强制 $\Gamma$ 点。
- 同位素效应讨论：经典 MD（含 MLFF）中质量在热力学平均中消去，观察不到动力学同位素效应——那是离子自由度的量子效应。
- 运行与质量判断：分两段（run1/run2）共 2750 步；续段要 `cp CONTCAR POSCAR` 并从 REPORT 提取末行 RANDOM_SEED 写入 INCAR 保证速度连续。质量判据：反向运行若自由能面出现滞后（hysteresis），说明转变速度/时间步不够小。慢增长只用于快速浏览自由能面；高质量积分用 blue-moon 方法（VASP 也实现，收敛可控但更贵）。
- REPORT 解读：每个 MD 步记录 SHAKE 参数与收敛、`cc>`（约束坐标值 $\xi'$）、`b_m>`（blue-moon 量：lambda、|z|^{-1/2}、GkT 等）、能量、温度、RANDOM_SEED；后处理时 `cc>` 第一个数是 $\xi'$，`b_m>` 第一个数是 $\partial A/\partial\xi$，用脚本提取后梯形积分得 $\Delta A(\xi)$。

**示例 9：反应物态的概率分布**

- 教学目标：计算反应坐标的概率密度；监控几何参数。
- 方法：在反应物态做自由 MD（无约束驱动），统计 $P(\xi') = \langle\delta(\xi' - \xi(R))\rangle$；若反应有势垒，$\xi$ 不会自发越过过渡态，$P(\xi')$ 在反应物阱内归一。
- 输入：沿用 MLFF（软链接 ML_FF）、Andersen 恒温器、5000 步、ICONST 监控；从轨迹统计 $\xi$ 的直方图得 $P(\xi_\mathrm{ref},R)$（如 $\xi_\mathrm{ref} = -1.5$ Å 处 $P \approx 1.54$ Å$^{-1}$）。

**示例 10：过渡态平均速度、反应速率与表观活化能（约束 MD）**

- 教学目标：对 MD 步数做收敛研究；计算反应速率与表观活化自由能；运行约束分子动力学。
- 理论：速率常数 $k_{R\to P} = (\langle|\dot{\xi}^*|\rangle/2)\,P(\xi_\mathrm{ref},R)\,\exp(-\Delta A_{\xi_\mathrm{ref},R\to\xi^*}/k_BT)$；表观活化自由能 $\Delta A^\ddagger = \Delta A_{\xi_\mathrm{ref},R\to\xi^*} - k_BT\ln[(h/k_BT)(\langle|\dot{\xi}^*|\rangle/2)P(\xi_\mathrm{ref},R)]$；广义速度 $\dot{\xi} = \sum_{i\mu} (\partial\xi/\partial R_{i\mu})\dot{R}_{i\mu}$，在 SHAKE 框架下经逆质量度规张量 $Z$ 计算 $\langle|\dot{\xi}^*|\rangle = \sqrt{2k_BT/\pi}/\langle Z^{-1/2}\rangle_{\xi^*}$。
- 输入：POSCAR 取恰在过渡态的结构；INCAR 与示例 9 相同但不设 INCREM——配合 ICONST 的 S 约束即成为过渡态处的约束 MD；run.sh 循环 10 段 x600 步，每段从 REPORT 提取 RANDOM_SEED 传递、CONTCAR→POSCAR 续接，保证可重复性。
- 结果：$\Delta A(\xi_\mathrm{ref}\to\xi^*) \approx 0.41$ eV（$\xi_\mathrm{ref} = -1.5$ Å）；对广义速度按 MD 步数做收敛研究后代入公式得速率与表观活化能。

**本部分涉及的 VASP 功能与标签汇总**：热力学积分与慢增长（`INCREM`、`LBLUEOUT`）、blue-moon 约束 MD、ICONST 约束语法（R=键长、S=线性组合）、REPORT 文件（cc>/b_m>/SHAKE/RANDOM_SEED）、`ML_MODE`（RUN/refit）、`POMASS` 同位素质量、`KSPACING` 控制 k 点密度、Andersen 恒温器（`ANDERSEN_PROB`）、Jarzynski 恒等式与静态方法对比、SHAKE 算法、约束 MD 与反应速率理论。

---

## 5. 机器学习力场（Machine-Learned Force Fields）

机器学习力场与从头算分子动力学结合，既能从第一性原理捕捉物理本质，又能以较低成本达到较长模拟时间。教程讲解 MLFF 的误差分析与超参数优化。

## 5.1 Part 1：机器学习力场的误差分析与超参数优化（LiH 案例）

> 来源：https://vasp.at/tutorials/latest/mlff/part1/ ｜ 练习文件：`mlff-part1.zip`
> 注：教程声明所用数据为教学裁剪，仅演示工作流，不用于科研。

**示例 1：训练集与测试集误差分析**

- 教学目标：对给定参考数据计算 MLFF 的训练集与测试集误差。
- 概念：训练集误差 = 拟合所用 DFT 值与 MLFF 预测之差；测试集误差在训练未用的外部测试集上评估。良好实践：测试集构型应与生产运行同条件（同原子数、同热力学相）。三种诊断情形：训练误差低+测试误差高 = 过拟合（需更多训练结构或调超参数）；两者相近且足够低 = 可用；训练误差高+测试误差低 = 测试集有偏、不够一般。参考 Tokita & Behler（arXiv:2308.08859）第 IV 节。
- 工作流：用 `ML_LMLFF = TRUE` + `ML_MODE = refit` 从现有 ML_AB（在线学习收集的能量/力/应力数据）重新拟合模型；ML_LOGFILE 的 ERR 段给出训练集 RMSE（能量 eV/atom、力 eV/Å、应力 kbar）；把 ML_FFN 复制为 ML_FF，用 `ML_MODE = run` + `IBRION = -1` 对 50 个测试构型逐个单点计算，再用 py4vasp 的 `MLFFErrorAnalysis.from_files(dft_data=..., mlff_data=...)` 计算逐构型与总平均的能量/力/应力误差并绘图。结论：训练与测试力误差均 <10 meV/Å 且同量级 → 模型可靠、可外推。

**示例 2：超参数扫描提升精度（ML_RCUT1）**

- 教学目标：理解超参数（训练过程中不优化、由用户设定的参数）并系统调优。
- 关键超参数表：`ML_RCUT1`（二体描述符截断半径）、`ML_RCUT2`（三体）、`ML_WTIFOR`/`ML_WTOTEN`/`ML_WTSIF`（力/能量/应力拟合权重）、`ML_EPS_LOW`（局域参考构型稀疏因子）、`ML_RDES_SPARSDES`（描述符稀疏因子）。
- 方法：交叉验证思想——先确认无过拟合，再系统（而非随机）扫描超参数；对每个 ML_RCUT1 值 refit 并记录训练误差（ML_LOGFILE）与力场（ML_FF）；为省时省盘，扫描阶段只看训练误差，最后只对训练误差最低的模型做测试集分析。结果：最优截断半径下力误差进一步下降。

**示例 3：稀疏化与性能优化（Pareto 前沿）**

- 教学目标：基于现有 ML_AB 优化力场性能。
- 理论：稀疏化提升性能——VASP 用 CUR 算法（`ML_LSPARSDES = True`）约化角描述符：把描述符协方差矩阵 A 分解为 C*U*R，基于 leverage scores 与 `ML_RDES_SPARSDES` 找出不重要的角描述符并移除；精度换性能。
- 方法：对多个 ML_RDES_SPARSDES 值分别 refit → 测试集精度 → 计时，绘制"时间 vs 精度"的 Pareto 前沿选最优点。计时要点：对多个 MD 步的力评估例程取平均；用 `ML_OUTBLOCK` 与 `ML_OUTPUT_MODE = 0` 关闭文件写出（写盘很慢会污染计时）。
- 输出解读：ML_LOGFILE 的 NDESC 段给出每元素径向/角向描述符数（稀疏化后角描述符 452→225，减半）；NDESC_SIC 段给出角描述符自相互作用校正项数；稀疏化模型测试集误差与精确模型相当（本例甚至略好），但性能显著提升。

**示例 4：用优化后的模型做生产运行**

- 输入：MD（`IBRION = 0`、`NSW = 100`、`POTIM = 2.0`）+ NVT（`MDALGO = 3` Langevin、`LANGEVIN_GAMMA = 3.0 3.0`、`TEBEG = 300`、`ISIF = 2`）+ `ML_LMLFF = .TRUE.` + `ML_MODE = RUN` + `ML_OUTBLOCK = 100`、`ML_OUTPUT_MODE = 0`；链接前面选定的 ML_FF。
- 要点：refit 不需要 POSCAR；run 模式按 POSCAR 计算能量/力/应力；POTCAR 在纯 MLFF 计算中不使用。

**本部分涉及的 VASP 功能与标签汇总**：`ML_MODE` 全家族（refit/run/…）、ML_AB/ML_FFN/ML_FF/ML_LOGFILE 文件体系、ML_LOGFILE 的 ERR/NDESC/NDESC_SIC 段、超参数（`ML_RCUT1`/`ML_RCUT2`/`ML_WTIFOR`/`ML_WTOTEN`/`ML_WTSIF`/`ML_EPS_LOW`/`ML_RDES_SPARSDES`/`ML_LSPARSDES`）、CUR 稀疏化与 leverage scores、`ML_OUTBLOCK`/`ML_OUTPUT_MODE` 计时设置、py4vasp MLFFErrorAnalysis、过拟合诊断与 Pareto 前沿方法学。

---

## 6. 表面（Surfaces）

真实块体总有边界——表面。教程通过沿一个方向把原胞延伸到真空中构建表面/板（slab）模型：表面结构弛豫、态密度与能带、分子吸附以及 STM 图像模拟。

## 6.1 Part 1：镍表面（表面能、局域态密度、表面能带）

> 来源：https://vasp.at/tutorials/latest/surface/part1/ ｜ 练习文件：`surface-part1.zip`

**示例 1：Ni(100) 表面的弛豫与表面能**

- 教学目标：用选择性动力学 + RMM-DIIS 弛豫特定层；为表面建模设置 k 点与单胞；计算单原子体系的表面能。
- 平板（slab）建模：沿一个方向拉长单胞，一部分填材料、一部分留真空；拉长方向只取 1 个 k 点（更多 k 点会引入与镜像单胞的虚假相互作用）；材料与真空区都要足够厚，使板中部呈现块体性质——该方向不能依赖周期性边界条件，必须消除其影响。
- 输入要点：POSCAR 为 5 层 Ni(100)（a=3.53 Å，z 向拉长 5x），Selective Dynamics 固定底部 3 层（F F F）、放开顶部 2 层（T T T）；INCAR：`ISPIN = 2`、`MAGMOM = 5*1`、`ISMEAR = 2`（MP 二阶）、`ALGO = Fast`、`EDIFF = 1E-6`、`IBRION = 1`（RMM-DIIS）、`NSW = 100`、`POTIM = 0.8`、`ISIF = 2`；k 网格 9x9x1（块体参考用 9x9x9）。
- 力与收敛：OUTCAR 的 FORCES 段逐步给出各层受力（electron-ion、ewald、non-local 分量）；弛豫后顶层力降至 ~$10^{-3}$ eV/Å；未显式设置时 EDIFFG = 10xEDIFF。
- 表面能计算（对称、化学计量比 slab 才适用）：$\sigma = (E_\mathrm{slab} - N E_\mathrm{bulk})/(2A)$，因子 1/2 因为切开块体产生两个等价表面。教程给出三步法：未弛豫 $\sigma^\mathrm{unrel} = 143.0$ meV/Å$^2$；弛豫顶两层后的平均 $\sigma^\mathrm{avg} = 141.5$ meV/Å$^2$；近似弛豫表面能 $\sigma^\mathrm{rel} = 2\sigma^\mathrm{avg} - \sigma^\mathrm{unrel} = 140.1$ meV/Å$^2$，弛豫能仅 ~3 meV/Å$^2$（密排金属小，其他材料可能大得多），与文献 PBE 值 138 meV/Å$^2$ 接近。
- 物理解读：表面能为正说明原子偏好块体环境；半局域泛函往往低估表面能；不同晶面表面能大小决定解理难易（如 (111) < (211)）。

**示例 2：Ni(100) 表面的局域态密度（LDOS）**

- 教学目标：计算局域电荷与局域自旋磁化；计算并绘制局域 DOS。
- 方法学（PAW 下获取局域量的三条路线）：① Wigner-Seitz 半径法（`RWIGS`，POTCAR 有默认值，LORBIT<10 时用，体积划分含糊）；② Wannier 函数类幺正变换（能带纠缠时难构造）；③ 投影算符法——把 KS 态投影到局域量子数 (l, m, $\sigma$) 标记的局域基（VASP 的 LORBIT 大于等于 10/PROCAR 路线）。
- 输入要点：`LORBIT = 11` + `RWIGS`；计算后 OUTCAR 末尾给出按层分辨的 charge 与 magnetization (x)（s/p/d 分解）：表面层（1、5）d 磁矩 ~0.76 $\mu_B$ 高于内层 ~0.68 $\mu_B$——表面增强磁性的经典结论。
- 分析：py4vasp 绘制逐层 LDOS，讨论表面态与 d 带变窄；PROCAR 存储 k、带、原子、轨道分辨的投影权重。

**示例 3：Ni(100) 表面能带（带字符分析）**

- 教学目标：对比块体与表面计算的 Wigner-Seitz 原胞；绘制带字符（band character）。
- 概念：slab 使实空间单胞沿 z 拉长，倒空间布里渊区沿对应方向被压扁，高对称点标签改变——用 SeeK-path 分别检查块体与 slab 的布里渊区，谨慎对应两者的 $\Gamma$-X-M-$\Gamma$ 路径。

![fcc 块体与 Ni(100) slab 的布里渊区对比](/imgs/2026-08-05/fig20.png)

- 输入要点：`ICHARG = 11` 读示例 1 的 CHGCAR 做非自洽能带（OUTCAR 显示 "charge density remains constant during run"）；KPOINTS 线模式 $\Gamma$-X-M-$\Gamma$（表面二维布里渊区路径，z 方向无色散）；`LORBIT = 11`。
- 分析：py4vasp `band.plot(selection="Ni(s, p, d)")` 按轨道字符着色能带，识别表面态与 d 带特征。

**示例 4：Ni(111) 表面弛豫（MP 方法深入与几何测量）**

- 教学目标：解释 Methfessel-Paxton 方法的目的与应用；用 py4vasp 测量键角/距离。
- 理论：金属的部分占据带使布里渊区积分面对阶跃函数，需极细 k 网格；MP 方法用对阶跃函数的光滑近似（对指定阶数多项式积分给出精确结果，阶数由 ISMEAR>0 的数值设定）缓解该问题；教程以非自旋极化计算弛豫 Ni(111) 五层板顶部两层并重复表面能分析。
- 分析：用 py4vasp 测量弛豫前后层间距与键角变化，比较 (111) 与 (100) 的表面能与弛豫行为（密排面弛豫更小、表面能更低）。

**本部分涉及的 VASP 功能与标签汇总**：slab 建模（真空层、垂直方向单 k 点）、Selective Dynamics 分层弛豫、`IBRION = 1` RMM-DIIS、`EDIFFG` 默认值（10xEDIFF）、表面能公式与三步弛豫能方法、局域量计算（`RWIGS`、`LORBIT`、PROCAR、投影算符法）、`ICHARG = 11` 非自洽表面能带、slab 布里渊区与高对称点、带字符绘图、MP 展宽原理、py4vasp 几何测量。

---

## 6.2 Part 2：CO 吸附在 Ni(111) 表面（偶极修正、吸附能、功函数、PDOS 与振动）

> 来源：https://vasp.at/tutorials/latest/surface/part2/ ｜ 练习文件：`surface-part2.zip`

**示例 5：CO/Ni(111) 弛豫（100% 覆盖度与偶极修正）**

- 教学目标：弛豫表面上的分子（100% 覆盖度）；理解并应用偶极修正；评估 slab-吸附质体系的几何变化。
- 背景：吸附问题构型空间大、势能面平坦；CO 吸附是催化反应性测试案例。吸附质改变 slab 电荷分布 → 体系产生有限偶极矩 → 周期边界条件下真空中出现势的梯度（虚假电场）；偶极修正（dipole correction）校正势、能量与力，凡破坏反演对称、允许偶极矩的体系都推荐使用。
- 结构准备：用 ASE 构建——`fcc111('Ni',(1,1,5),a=3.53)` + `add_adsorbate(slab, CO, 1.8, 'ontop')`；ontop（顶位）是高覆盖度下最稳定吸附位；低覆盖度需超胞。

![Ni(111) 表面高对称吸附位示意](/imgs/2026-08-05/fig21.png)

- 输入要点：两步流程——先 INCAR.1 收敛电荷分布，再 INCAR.relax 做结构弛豫；偶极修正：`LDIPOL = .TRUE.`、`IDIPOL = 3`（沿 z）、`DIPOL = 0.5 0.5 z_COM`（反应坐标中心，可用 ASE 的质心计算）；`ENCUT = 400`（由 O 的 ENMAX 决定；对 CO 这类短键分子应测试硬势 C_h/O_h 的影响）；Ni(111) 有磁性但 CO 无磁矩，此处做非磁计算演示流程；长单胞+破缺反演对称易发生 charge sloshing，收敛困难时调节密度混合参数（`AMIX` 等）。
- 收敛性教学：不加偶极修正时总能对真空厚度 L 呈 $1/L^2$ 依赖，阻碍真空收敛，实际模拟的是"slab 层状块体"。
- 结果分析：OUTCAR 中搜 `dipolmoment` 与 `dipol+quadrupol energy correction`（金属 slab+弱极性吸附质时很小；带电吸附质如 K、O 自由基时大）；py4vasp 测得 CO–Ni 距离 1.76 Å、顶两层向外弛豫 3.2%、吸附后 C–O 键长增加 ~1%（1.142→1.154 Å）；比较几何量不要求各计算 ENCUT 相同（各自收敛即可），但比较总能必须相同。

**示例 6：吸附能与功函数**

- 教学目标：计算分子吸附能；计算并分析功函数。
- 吸附能：E_ads = E(CO/slab) - E(slab) - E(CO)，三个量需相同计算参数；教程强调半局域泛函低估表面形成能但高估吸附能。
- 功函数 $\varphi$ = 真空势 - 费米能：INCAR 设 `LVHAR = .TRUE.`（写出 Hartree 势，LOCPOT）；py4vasp `mycalc.workfunction` 给出沿表面法线的真空势平台与费米能（本例 $\varphi \approx 5.18$ eV，生产级计算 5.10 eV）。为什么用 LVHAR 而非 LVTOT：交换关联势在真空中数值上难评估、物理上应趋于零，真空能级不含 xc 势更准确（xc 效应已完全包含在自洽场与费米能中）；slab 两侧真空势略有差异是因为 slab 不对称（一侧弛豫、一侧未弛豫）。

**示例 7：CO/Ni(111) 分波 DOS（Blyholder 模型的验证）**

- 教学目标：计算表面上分子的分波 DOS；绘制并分析平面平均势与功函数变化。
- 理论（Blyholder 模型）：CO 的 $5\sigma$（HOMO，C 的 $sp_z$ 杂化）与 $2\pi^*$（LUMO）与金属表面作用形成成键/反键杂化轨道：$5\sigma$–金属反键态在费米能之上 → 电荷捐赠（donation）；$2\pi^*$–金属成键态被占据 → 反馈（back-donation）。
- 输入要点：从示例 5 的 CONTCAR/CHGCAR 重启（`ICHARG = 1`），加密 k 网格（13x13x1），`ISMEAR = -5`（四面体法，DOS 与精确总能推荐）、`NEDOS = 5001`、`EMIN/EMAX = ±10`、`LORBIT = 11`、`NBANDS = 48`、`LMAXMIX = 6`、偶极修正全套。
- 分析：LORBIT=11 的简单投影无法直接给出分子轨道，但按 LCAO/分子轨道理论组合原子投影即可：py4vasp 读取 C(pz)+C(s)（$5\sigma$ 特征）、O(pz)、Ni d 的投影 DOS——可见 $5\sigma$ 分裂为成键/反键两部分，费米能上方的小峰表明吸附后 $5\sigma$ 未全占（电荷转移到表面）；$1\pi$ 态以 O 为主、$2\pi^*$ 以 C 为主。教程演示 Hann 窗卷积平滑表面 DOS 的技巧，并警告平滑可能引入赝特征或抹掉真实峰，须个案决定。

**示例 8：CO/Ni(111) 的振动频率**

- 教学目标：用有限差分法计算吸附分子的振动频率。
- 输入要点：`IBRION = 5/6` + `NFREE` + 小 `POTIM`；POSCAR 选择性动力学只放开吸附质原子（C、O）的位移，slab 冻结，从而只对分子模式求 Hessian，大幅降低成本；从 OUTCAR 动力学矩阵部分读取频率，与气相 CO 伸缩频率对比（吸附后红移，与键长增大一致），并讨论表面模式的截断处理。

**本部分涉及的 VASP 功能与标签汇总**：ASE 表面/吸附质建模、偶极修正（`LDIPOL`/`IDIPOL`/`DIPOL`）、真空厚度收敛与 $1/L^2$ 效应、吸附能计算、功函数（`LVHAR`/LOCPOT、真空势平台）、`LMAXMIX = 6`（含 f 电子/磁性体系）、charge sloshing 与混合参数、Blyholder 模型的 PDOS 验证、DOS 平滑的注意事项、slab 上吸附质振动的有限差分法（选择性自由度冻结基底）。

---

## 6.3 表面科学（三）：STM 扫描隧道显微镜模拟（恒定高度 / 恒定电流模式）

> 来源：<https://vasp.at/tutorials/latest/surface__part3.html>
> 配套输入文件：[surface_part3.zip](https://vasp.at/tutorials/latest/surface__part3.zip)

本页通过两个例子演示如何用 VASP 计算**部分电荷（partial charge）**并结合 py4vasp 的 Tersoff–Hamann 模型模拟 STM（扫描隧道显微镜）图像：例 9 为石墨 (0001) 表面的**恒定高度模式**，例 10 为石墨烯的**恒定电流模式**。

### 例 9：恒定高度模式下的石墨 (0001) STM 图像

**教学目标**
- 掌握用 `LPARD`/`NBMOD`/`EINT` 计算能量窗口内的部分电荷；
- 理解 STM 模拟的物理基础（Bardeen 隧穿理论与 Tersoff–Hamann 近似）；
- 用 py4vasp 的 `to_stm` 在恒定高度模式下生成 STM 图像，并解释偏压极性、三角格子形貌的成因。

**物理背景**
- 在 Bardeen 隧穿理论框架下，隧穿电流正比于针尖处样品的局域态密度。Tersoff–Hamann 近似把针尖视为 s 波（球形对称）尖端，使 STM 图像正比于费米能附近的局域态密度，从而可以用 DFT 的部分电荷直接模拟。
- 若需要更高精度（如考虑针尖轨道取向），可采用 Chen 的导数规则。本教程的实现位于 py4vasp 中。

**第一步：静态自洽计算**
- 对石墨 (0001) 表面板模型做一次常规静态 SCF 计算，得到收敛的 WAVECAR/CHGCAR。
- 教程特别提醒 `PREC` 设置：使用 `PREC=Single` 时电荷密度的精细网格与粗网格一致（`NGX=NGXF`），可避免真空区电荷出现插值伪影；这是表面/真空体系的一个实用细节（教程也指出可用 py4vasp 的高斯平滑作为替代方案）。

**第二步：部分电荷后处理**
- 新建 `partial_charges/` 子目录，把 WAVECAR 复制进去，配合只读波函数的后处理 INCAR（教程给出的关键参数组合，示例文件内容此处略）：
  - `LPARD = .TRUE.`：开启部分电荷计算；
  - `LPARDH5 = .TRUE.`：把部分电荷写入 vaspout.h5，便于 py4vasp 读取；
  - `NBMOD = -3`：把 `EINT` 解释为相对费米能级的能量窗口；
  - `EINT = -0.25 0.0`：能量窗口为费米能以下 0.25 eV，相当于约 -0.25 V 的样品偏压，只积分占据态；
  - `LSEPK = .FALSE.`、`LSEPB = .FALSE.`：对所有 k 点和能带求和而不分开输出；py4vasp 做 STM 时要求二者均为 `.FALSE.`。
- 教程说明：对该体系实际只有第 18、19 条能带在 k 点 33 处落入能量窗口并贡献部分电荷。

**第三步：py4vasp 生成 STM 图像**
- 调用 `calc.partial_density.to_stm(supercell=7, tip_height=3)`：`supercell=7` 表示把表面原胞放大 7x7 以便观察周期形貌，`tip_height=3` 设定针尖距表面的高度（Å）。
- 结果解读：
  - **负偏压**对应占据态、**正偏压**对应未占据态的隧穿成像；
  - 石墨表面呈现**三角格子**形貌（由于 A/B 层堆垛，每层只有一个亚晶格位靠近表面被成像）。

![石墨 0001 表面结构](https://vasp.at/tutorials/latest/_images/Graphite_ambient_STM.jpg)

### 例 10：恒定电流模式下的石墨烯 STM 图像

**教学目标**
- 在同样的部分电荷框架下，切换到更接近实验的**恒定电流模式**；
- 理解两种成像模式的差别与各自适用场景。

**计算设置**
- 第一步静态 SCF 与第二步部分电荷计算与例 9 完全相同（同样的 `LPARD`/`NBMOD`/`EINT`/`LSEPK`/`LSEPB` 组合），示例文件内容此处略。

**py4vasp 生成恒流模式图像**
- 调用 `to_stm(selection="constant_current")`；默认隧穿电流为 **1 nA**。
- 算法：从真空中部开始向下扫描针尖高度，对平滑化后的部分电荷积分，直到局域隧穿电流等于设定电流值，记录该高度即得恒流形貌图。
- 与恒定高度模式相比，恒流模式更贴近真实 STM 实验（反馈回路保持电流恒定），而恒定高度模式计算更简单、适合平整表面快速预览。

### 功能与标签汇总
- 部分电荷：`LPARD`、`LPARDH5`、`NBMOD=-3`、`EINT`（相对费米能）、`LSEPK`/`LSEPB`
- 真空区精度细节：`PREC=Single` 与 NGX=NGXF、py4vasp 高斯平滑
- STM 理论：Bardeen 隧穿理论、Tersoff–Hamann 近似（s 波针尖）、Chen 导数规则
- py4vasp：`calc.partial_density.to_stm(supercell=..., tip_height=...)`、`selection="constant_current"`（默认 1 nA）
- 物理结论：偏压极性决定占据/未占据态成像；石墨呈三角格子形貌

---

## 7. 过渡态（Transition States）

化学反应在反应物与产物两个极小值之间进行，势能面上连接二者的能量极大值即过渡态。学习用静态与动态方法建模反应过渡态。

## 7.1 过渡态（一）：基础过渡态方法（NEB / IDM / 肘形图 / IRC）

> 来源：<https://vasp.at/tutorials/latest/transition_states/part1/>
> 配套输入文件：[transition_states-part1.zip](https://vasp.at/tutorials/latest/transition_states-part1.zip)
> 注：本教程使用的预训练机器学习力场（MLFF）可从教程页面给出的 Google Drive 链接下载，放入 `transition_states/mlff` 目录。

本页包含四个例子：例 1 用 NEB 求反应路径，例 2 用改进二聚体方法（IDM）优化过渡态，例 3 构建 H2/Cu(111) 解离吸附的肘形图（elbow plot），例 4 用内禀反应坐标（IRC）追踪 SN2 反应。

### 例 1：用 NEB 计算反应路径（Si 自扩散）

**教学目标**
- 准备并运行 NEB（nudged elastic band，弹性带）计算；
- 优化图像结构与反应路径；
- 可视化图像结构与能量路径。

**任务与方法原理**
- 任务：在 15 原子 Si 原胞超胞中，用 4 个 image 的 NEB 计算 Si 原子向空位扩散的反应路径。
- NEB 方法（J. Chem. Phys. 113, 9978 (2000)）：在反应物与产物结构之间插值一系列点（image，通常 4–20 个）来模拟最低能量路径（MEP）。image 之间以及 image 与初末态之间用弹簧相互作用连接，形成连续的"弹性带"。
- 带的受力分解为垂直于带的真实力 $F_{i,\perp}$ 和平行于带的弹簧力 $F_{i,\mathrm{spring}}$（即"nudging"处理），按合力 $F_i$ 弛豫即得反应路径。
- image 越多越精确但越昂贵、越难收敛；本例只用了 4 个。

![NEB 方法示意图](/imgs/2026-08-05/fig22.png)

![Si 空位超胞：蓝色 Si 原子向红色空位扩散](/imgs/2026-08-05/fig23.png)

**输入文件要点**
- 目录 `$TUTORIALS/transition_states/e01_NEB`；提供初态 `00/POSCAR`、4 个 image `01–04/POSCAR`、末态 `05/POSCAR`（均已预弛豫），并附 ASE 建模与插值脚本（内容略）。
- `INCAR.neb`（NEB 计算）与 `INCAR.sp`（单点）标签相同，关键标签（示例文件内容略）：
  - 电子最小化：`ENCUT`（平面波截断能）、`PREC`（精度）、`EDIFF`（电子自洽收敛判据；教程建议用负的力阈值，即单位 eV/Å 的力判据）、`ISMEAR`/`SIGMA`（展宽）；
  - 离子弛豫：`NSW`（最大离子步数，单点时取 0）、`EDIFFG`（离子弛豫收敛判据）、`IBRION=1`（准牛顿算法）、`POTIM`（离子步长；若出现 "BRIONS problems: POTIM should be increased" 可增大 POTIM 或换 IBRION 算法）、`ISIF=2`（只弛豫原子位置，晶胞固定）；
  - NEB 专用：`IMAGES=4`（image 数；MPI 进程数必须是 image 数的整数倍）、`SPRING=-5`（弹簧常数）。
- `KPOINTS` 取 $\Gamma$ 中心 k 网格；教程提示若要更高质量结果应加密 k 网格、增大超胞并把 `EDIFFG` 收紧到约 -0.01 eV/Å。

**计算流程与结果**
- NEB 必须从父目录 `e01_NEB` 启动，VASP 会在各子目录中并行计算每个 image；屏幕 stdout 来自 image 1，其余 image 的输出在各自目录。
- 初态与末态在 NEB 中固定不动，需分别做单点计算得到其能量，才能拼出完整能量剖面。
- 提升能量精度的建议：加密 k 网格（如 4x4x4），并建议最后用四面体方法做一次单点。
- 分析（py4vasp + 脚本）：把各 image 结构投影到二维截面，可清楚看到扩散的 Si 原子从初始位置直线移动到空位；能量剖面大致对称，峰值约 **0.7 eV**。

![各 NEB image 投影结构](/imgs/2026-08-05/fig24.png)

- 结论与经验：路径研究很复杂时可取中间 image 作为寻找过渡态的试探结构；image 少（如 $\leq 4$）时单次计算更易收敛，image 太多会显著增加开销。

### 例 2：改进二聚体方法（IDM）优化过渡态（H2/Cu(111) 解离）

**教学目标**
- 用振动计算求虚频模式；
- 用 IDM 优化 H2 在 Cu(111) 上解离的过渡态；
- 用虚频个数验证过渡态。

**方法原理**
- NEB 给路径，而过渡态（TS）本身——鞍点，只有一个虚频——需要用 IDM（J. Chem. Phys. 111, 7010 (1999)）等局部鞍点搜索方法。IDM 只用一阶力（而非 Hessian）：
  1. 从试探结构的最不稳定振动模式取初始方向 $u_\xi$（点 $q$）；
  2. 沿 $u_\xi$ 在 $q+\delta u_\xi$ 处取第二个点，两者构成"二聚体"；
  3. 二聚体绕 $q$ 在势能面上旋转角度 $\phi_1$；
  4. 再旋转 $u_\xi$ 一个 $\phi_{\min}$ 使负曲率 $\lambda$ 最小，定义搜索方向 $\bar N$，平移 $\epsilon$ 到 $q+\epsilon\bar N$；
  - 重复"旋转 + 平移"直到鞍点，即只有一个负振动模式的过渡态。

![IDM 四步弛豫示意图](/imgs/2026-08-05/fig25.png)

**输入文件要点**
- 目录 `e02_IDM`，子目录 `freq1/`（试探结构振动分析）、`idm/`（IDM 弛豫）、`freq2/`（过渡态验证）；`POSCAR.idm` 在原子坐标之后附加的行是**初速度/试探不稳定方向**。
- 关键标签（示例文件内容略）：
  - 频率计算：`IBRION=5`（有限差分法）、`NFREE`（位移点数）；
  - 机器学习力场加速：`ML_LMLFF=T`、`ML_MODE=run`（使用预训练 MLFF，需把 `ML_FFN` 文件软链接为 `ML_FF`）；
  - IDM：`IBRION=44`；`NELMIN`（最少电子步数）；
  - `NSW=20`（为教学速度只跑 20 步；真正收敛需数百到数千步）。
- 表面体系 KPOINTS 用二维 k 网格：垂直表面方向只取 $\Gamma$ 点，配合大真空避免周期性镜像相互作用；更高质量计算应收敛 k 网格、真空厚度与层数（通常 $\geq 4$ 层）。

**计算流程与结果**
1. **试探结构振动分析**（freq1）：得到振动频率，本例出现两个虚频；过渡态应当只有一个虚频，IDM 将沿较大的虚频模式弛豫、消除较小的那个。可用 `vibFreq.py` 生成 `hesseModes_shifted.molden`（用 JMol/molden 可视化），或用 py4vasp 提取频率并导出 molden 文件。

![H2 吸附在 Cu(111) 桥位，橙色箭头为 H-H 伸缩振动模式](/imgs/2026-08-05/fig26.png)

   该模式是 H–H 键伸缩——解离时正是此键断裂，因此沿此不稳定模式做 IDM。
2. **IDM 弛豫**（idm）：从 OUTCAR 提取最大虚频模式写入 POSCAR 的附加行（`POSCAR.idm` 已提供），复制为 POSCAR 后运行 `IBRION=44` 的弛豫。
3. **验证**（freq2）：再做一次振动分析，确认只剩**一个虚频**，即成功找到过渡态。
- 能量-步数曲线显示偶发的快速能量下降；教程展示了一个 ~600 步的更长 IDM 计算：前期缓慢衰减，约 600 步处快速下降后平台化——正是额外虚频被消除之处。

### 例 3：H2/Cu(111) 解离的肘形图（elbow plot）

**教学目标**
- 为肘形图准备结构；
- 可视化势能面（PES）的二维截面。

**方法原理**
- 双原子分子有 6 个自由度；固定表面后，取分子质心高度 $z$ 与键长 $r$ 两个自由度即可构造二维 PES 的等高线图，即"肘形图"（Phys. Rev. Lett. 73, 1400 (1994)；Phys. Rev. B 62, 8295）。

![肘形图定义：r 与 z 网格上逐点计算能量](/imgs/2026-08-05/fig27.png)

**输入文件要点**
- 目录 `e03_Elbow_plot`；新标签 `LCHARG`/`LWAVE` 设为 `.FALSE.`——由于要跑大量单点计算，不写出 CHGCAR/WAVECAR 以节省磁盘；`NSW=0`（单点）。
- KPOINTS 只用 $\Gamma$ 点（教学演示用，网格非常稀疏）。

**计算流程与结果**
- 用 ASE 脚本对 $z$、$r$ 网格（各 6+1=7 个点，共 $7^2=49$ 个结构）批量生成 POSCAR 并建目录，再用运行脚本批量提交。
- 因为只用 $\Gamma$ 点，使用 gamma-only 版本 `vasp_gam` 更高效（波函数为实数）。
- 用 py4vasp 提取各结构能量并画等高线图：

![肘形图（49 点，$\Gamma$ 点）](/imgs/2026-08-05/fig28.png)

![肘形图（121 点，3x3x1 k 网格）](/imgs/2026-08-05/fig29.png)

- 两图对比：121 点（11x11 网格、3x3x1 k 网格）的等高线更平滑，且因 k 网格更密能量整体更低（更负）。完全收敛还需更细网格。
- 图中两个极小值分别对应：桥位上的 H2 分子、以及解离后吸附在 fcc 与 hcp 位点的两个 H 原子。

### 例 4：内禀反应坐标（IRC）——SN2 取代反应

**教学目标**
- 沿 IRC 向前（产物）与向后（反应物）追踪；
- 可视化反应剖面与反应物/产物结构。

**任务**
- 在 12x12x12 Å$^3$ 超胞中，用机器学习力场沿 IRC 追踪氯甲烷 SN2 氯离子取代反应：从过渡态出发分别走向反应物与产物。

**输入文件要点（示例文件内容略）**
- IRC 计算用 `IBRION=40`；
- `IRC_DIRECTION`：初位移方向，`1` 平行于不稳定方向（朝向产物），`-1` 反平行（朝向反应物）；
- `IRC_STOP=20`：能量相对平稳的连续步数达到该值时终止；
- `IRC_VNORM0`：设定恒定速度矢量的模；
- `NELECT=22`：给体系多加一个电子（氯离子带负电荷）；
- MLFF 标签 `ML_LMLFF`/`ML_MODE` 同例 2；KPOINTS 只用 $\Gamma$ 点（气相反应足够）。

**计算流程与结果**
- 把预训练 MLFF（`MLFF_e04/ML_FFN`）软链接为 `ML_FF` 以节省磁盘；IRC 需要大量步数，MLFF 把计算从"数小时纯第一性原理"降到"几秒完成 ~450 步"。
- 需要两个目录：向后 `m/`（`IRC_DIRECTION=-1`）与向前 `p/`（`IRC_DIRECTION=1`），其余输入相同。
- 用 `ircShift.sh` 提取每步的能量与 IRC 坐标，画反应剖面：曲线近乎完美对称——氯离子从氯原子正对面进攻碳、使甲基翻转并弹出原来的 Cl$^-$；反应物与产物相同（自交换反应），故曲线对称。

![SN2 反应示意](/imgs/2026-08-05/fig30.png)

- 用 py4vasp 可分别查看反应物/产物结构；用 `pos2xyzAnim_rev.py`/`pos2xyzAnim_for.py` 每 10 步提取一次 IRC 坐标并拼成 `anim.xyz`，即可在 JMol 中播放整个 SN2 反应动画。

### 功能与标签汇总
- NEB：`IMAGES`、`SPRING`，父目录运行、`00`–`05` 目录结构，MPI 进程数须为 image 数整数倍；单点用四面体方法提精度
- IDM：`IBRION=44`，POSCAR 附加行给试探不稳定方向；过渡态判据为恰好一个虚频
- 振动分析：`IBRION=5`、`NFREE`；`vibFreq.py` → molden 可视化，py4vasp 提取频率
- 肘形图：`LCHARG`/`LWAVE=.FALSE.` 省磁盘、`NSW=0` 单点、`vasp_gam` 加速、ASE 批量建模
- IRC：`IBRION=40`、`IRC_DIRECTION=±1`、`IRC_STOP`、`IRC_VNORM0`；带电体系用 `NELECT`
- MLFF 加速：`ML_LMLFF`、`ML_MODE=run`、`ML_FFN` 软链接为 `ML_FF`

---

## 7.2 过渡态（二）：静态方法实战（沸石催化羰基化反应全流程）

> 来源：<https://vasp.at/tutorials/latest/transition_states/part2/>
> 配套输入文件：[transition_states-part2.zip](https://vasp.at/tutorials/latest/transition_states-part2.zip)
> 注：本页 MLFF 需另行下载（页面给出链接）放入 `transition_states/mlff`。

本页以真实催化体系为例，把第一部分学到的方法串成完整工作流：合成气在菱沸石（chabazite）沸石上催化合成乙醇，其中**决速步（RDS）是甲基化沸石的羰基化**（CO + CH3-Z → CH3CO$^+$ + Z$^-$）。五个例子依次为：例 5 试探过渡态的振动分析、例 6 IDM 优化过渡态、例 7 鞍点验证与过渡态自由能、例 8 IRC 反应路径、例 9 反应热力学与速率常数。

![菱沸石结构：红 O、蓝 Si、青 Al；沸石空腔提供巨大内表面](/imgs/2026-08-05/fig31.png)

**反应机理背景**：合成气制乙醇共六步元反应（CO+2H2→CH3OH；CH3OH+H-Z→H2O+CH3-Z；CO+CH3-Z→CH3CO$^+$+Z$^-$；CH3CO$^+$+Z$^-$→CH3CO-Z；CH3OH+CH3CO-Z→CH3COOCH3+H-Z；CH3COOCH3+2H2→CH3CH2OH+CH3OH），第 3 步羰基化为决速步（PCCP 13, 2603 (2011)）。本教程聚焦正向活化自由能 $\Delta A^{\ddagger}_{R\to P}$（$\ddagger$ 表示过渡态）。

### 例 5：试探过渡态的振动分析

**教学目标**
- 对试探过渡态做振动分析；
- 找到用于过渡态弛豫的试探不稳定方向。

**任务与设置**
- 试探过渡态的构建：平面 CH3 基团置于沸石 O 与 CO 的 C 中间，两个距离均取 2.1 Å。若猜测接近鞍点则应只有一个虚频；若有多个虚频，取最强虚频作为 IDM 的试探方向。
- 目录 `e05_Static_Vib_analysis`。INCAR 要点（示例文件内容略）：
  - `IBRION=5`（有限差分振动计算）、`POTIM`（位移量）、`NFREE`（位移数）；对 IBRION=5，`NSW` 只需为正，VASP 会将其重置为 3x原子数；
  - `ML_LMLFF`+`ML_MODE=run`：仅运行预训练 MLFF 不再训练，因此不需要电子结构标签；MLFF 使"几秒"完成振动分析（纯第一性原理需数小时）；
  - KPOINTS 只用 $\Gamma$ 点；POTCAR 含 Al/C/H/O/Si。
- 结果：126 个频率模式，前 5 个为虚频（输出以 `f/i` 标记）——说明结构还不是极小点。用 `vibFreq.py` 或 py4vasp 提取频率并生成 molden 可视化文件；下一步将把虚频从 5 个降到 1 个。

### 例 6：用 IDM 优化过渡态

**教学目标**
- 沿最"硬"（最大）虚频方向弛豫过渡态猜测结构；
- 比较收敛结构与初始猜测。

**输入要点**
- 目录 `e06_Static_TS_Opt`；`POSCAR.init` = 例 5 的 POSCAR + 从 `hesseModes.dat`（vibFreq.py 产生）提取的最强虚频力（写在坐标之后）。
- INCAR 要点（示例文件内容略）：`IBRION=44`（IDM）；`NSW` 此时为结构弛豫的离子步数；新标签 `EDIFFG`——负值表示以力（eV/Å）为收敛判据。

**计算与分析**
- IDM 回顾：沿二聚体轴 **N** 前后各取一点构成二聚体，绕中点旋转找 PES 上最小曲率 $\lambda$ 方向，再沿轴计算力并把中点反向平移向鞍点；每 IDM 步约 4 个离子步，迭代至收敛。
- 从 OUTCAR 提取最大梯度画图：前 100 步梯度快速下降，振荡收敛，200 步后基本不变。
- 再从 OUTCAR 提取曲率 $\lambda$ 画图（曲率每 4 步更新一次）：曲率先升后缓慢收敛并振荡，需 ~600 步才平台化——**曲率收敛才是限速因素**，这解释了为何力在 200 步内已收敛而计算仍跑到 ~600–800 步。
- 结构对比（py4vasp）：初始与收敛结构视觉差别很小，但 C–C 与 C–O 距离均从 ~2.1 Å 收缩到 ~2.0 Å，即过渡态中反应物相互靠拢。
- 效率对比：若用 Hessian 方法需每步做频率计算，代价巨大；IDM 每步只需几次单点计算。

### 例 7：鞍点验证与过渡态自由能

**教学目标**
- 验证弛豫结构为一阶鞍点（恰好一个虚频）；
- 得到 IRC 所需的不稳定试探方向；
- 计算过渡态在 440 K 的自由能。

**理论与输入**
- 一阶鞍点应只有一个虚频（扣除消失的平动模式）。虚频模式可用作 IRC 的不稳定试探方向。
- 振动热力学量：零点振动能
$$E_{ZPV} = \sum_i \frac{\hbar}{2}\nu_i$$
  以及振动焓 $H_{vib}(T)$、振动自由能 $A_{vib}$。此处自由能是 Helmholtz 自由能（正则 NVT 系综，$N$、$V$、$T$ 固定），而非 Gibbs 自由能。总 Helmholtz 自由能 $A_{tot}=E_{tot}+A_{vib}$ 是后续算速率常数所需。
- 目录 `e07_Static_Saddle_point`，输入与例 5 相同，POSCAR 为例 6 收敛结构（示例文件内容略）。

**计算与分析**
- 再做一次 IBRION=5 振动分析；脚本提取频率并构造 Hessian 与动力学矩阵、输出本征值/本征矢。
- 结果对比：IDM 前 100 cm$^{-1}$ 以上有 3 个虚频、另有 2 个更小虚频；IDM 后只剩 1 个 >100 cm$^{-1}$ 的虚频，其余为 <10 cm$^{-1}$ 的消失平动模式——**平动模式在计算振动热力学量前须手动剔除**。
- 虚频模式沿 O–C–C 轴方向（molden/JMol 可视化）。
- `vibpartition_eV.py`：剔除平动模式后计算 440 K（接近实验温度）的谐振振动贡献，输出 `thermo.dat`，其字段含义：

| 标签 | 含义 |
|---|---|
| `nDOF` | 自由度（=3x原子数） |
| `nVIB` | 去除虚频后的振动模式数 |
| nVIB 下一行 | 振动配分函数 $Q_{vib}$（科学计数法，无显式标签） |
| `T` | 温度（K） |
| `ZPV` | 零点振动能 |
| `Hvib - ZPV` | 振动能的热贡献（温度导致激发态布居） |
| `Hvib` | 振动焓 |
| `Svib` | 振动熵 |
| `Gvib` | 振动 Helmholtz 自由能 |

  也可用 ASE 的 `thermochemistry` 模块得到同样结果。
- $G_{vib}$ 加上 OUTCAR 总能即得过渡态 Helmholtz 自由能 $\Delta A^{\ddagger}_{R\to P}$。孤立自由能不能直接与实验对比，需构造自由能差（活化自由能）。

### 例 8：IRC 反应路径

**教学目标**
- 沿不稳定试探方向追踪 IRC；
- 找到并可视化反应物与产物态。

**方法与输入**
- 鞍点只有一个虚频；沿该模式正向传播到产物、反向传播到反应物，从而验证过渡态确实连接预期的反应物与产物。传播用阻尼速度 Verlet 算法（J. Chem. Phys. A 106, 165 (2002)）。
- 目录 `e08_Static_Reaction_pathway`；INCAR 相对例 6 改为 `IBRION=40`；新标签 `IRC_VNORM0`（速度矢量重标定的恒定值，提供阻尼）；`IRC_DIRECTION=1` 正向（产物）、`-1` 反向（反应物）（示例文件内容略）。

**计算与分析**
- 从例 7 复制 POSCAR（含不稳定模式试探方向）；在 `m/`（反应物）与 `p/`（产物）两个目录分别运行。
- `ircShift.sh` 输出 `ircShift.dat`（全部 IRC 点），画反应剖面：约 -3 Å 处为反应物（CO 与 CH3-Z）；CO 靠近、CH3 翻转，能量升至 0 Å 处极大值（过渡态）；自由乙酰阳离子形成后能量快速降至 ~4 Å 处的产物。
- 用 py4vasp 查看反应物/过渡态/产物结构；`pos2xyzAnim_rev.py`/`pos2xyzAnim_for.py` 每 10 步提取 IRC 坐标生成 `anim.xyz` 反应动画（JMol 播放）。
- IRC 终止条件：能量连续 `IRC_STOP`（默认 20）步单调上升；因此从计算末尾数第 20 个结构最接近势能极小。`pickMin_irc.sh` 从 XDATCAR 提取该结构为 `POSCAR.pick`。由于本例 `IRC_VNORM0=0.0005` eV/Å 很小，提取结构已非常接近极小点，无需再弛豫。

### 例 9：反应热力学（静态）与速率常数

**教学目标**
- 对 IRC 得到的反应物结构做振动分析、确认其为极小点；
- 比较谐振量子方法与半经典方法；
- 计算活化自由能与速率常数。

**理论**
- 总能由第一性原理算得，而实验测量的是自由能等热力学量；两者对比需计入零点振动能、振动热贡献与熵项。
- 活化自由能 $\Delta A^{\ddagger}_{R\to P}=A_{TS}-A_R$ 代入 Eyring 方程得速率常数：
$$k_{R\to P} = \frac{k_B T}{h} e^{-(\Delta A^{\ddagger}_{R\to P}/k_B T)}$$
  其中 $k_B$ 为 Boltzmann 常数、$h$ 为 Planck 常数、$T$ 为温度。
- Helmholtz 自由能与配分函数：$A = -k_B T \ln Q$；量子谐振振动配分函数
$$Q_{\mu,vib} = \prod_{i=a}^{N_{vib}} \frac{e^{-h\nu_i/(2k_BT)}}{1 - e^{-h\nu_i/(k_BT)}}$$
  半经典振动配分函数
$$Q^{SC}_{\mu,vib} = \prod_{i=a}^{N_{vib}} \frac{h\nu_i}{2k_B T}$$

**输入与计算**
- 目录 `e09_Static_Thermodynamics`；POSCAR 取例 8 的 `m/POSCAR.pick`；INCAR/KPOINTS/POTCAR 同例 5（示例文件内容略）。
- 振动分析确认：除消失的平动虚频外全为正频——确为能量极小点，即反应物。
- `vibpartition_eV.py` 在 440 K 计算谐振贡献（`thermo.dat` 格式同例 7）；也可用 ASE `thermochemistry`。
- 半经典对比：用 `vibpartitionClassic_eV.py` 由 $Q^{SC}_{\mu,vib}$ 计算半经典振动自由能。结果：半经典与量子 Helmholtz 自由能各自相差 ~1 eV，但两者的 $\Delta A^{\ddagger}$ 之差仅 ~0.02 eV——**量子效应很小，半经典近似成立**。这是使用 MD（半经典力、不含零点振动能）的重要验证。
- 用 Eyring 方程计算速率常数；振动能量差异使速率常数减小约 40%——小能量差经指数放大后对可测量影响显著。同理可算产物→过渡态方向的速率常数。

### 功能与标签汇总
- 工作流串联：猜测 TS → IBRION=5 振动分析 → IBRION=44 IDM（`EDIFFG` 力判据）→ 鞍点验证（唯一虚频）→ IBRION=40 IRC（`IRC_DIRECTION`、`IRC_VNORM0`、`IRC_STOP`）→ 振动热力学 → Eyring 速率常数
- 关键脚本：`vibFreq.py`（`hesseModes.dat`/molden）、`vibpartition_eV.py`/`vibpartitionClassic_eV.py`（`thermo.dat`）、`ircShift.sh`、`pickMin_irc.sh`、`pos2xyzAnim_*.py`
- 概念：一阶鞍点判据、消失平动模式的剔除、Helmholtz（NVT）自由能、谐振 vs 半经典配分函数、曲率收敛决定 IDM 步数

---

## 7.3 过渡态（三）：动态方法（反应坐标、Blue Moon 系综、动态热力学）

> 来源：<https://vasp.at/tutorials/latest/transition_states/part3/>
> 配套输入文件：[transition_states-part3.zip](https://vasp.at/tutorials/latest/transition_states-part3.zip)

本页继续同一沸石催化羰基化体系，改用**分子动力学（MD）系综**视角，三个例子依次为：例 10 用 MD 定义并检验反应坐标、例 11 用 Blue Moon 约束 MD 验证反应坐标、例 12 动态计算活化自由能与速率常数。

**为什么需要动态方法**：静态方法用单一驻点结构代表反应物/过渡态/产物，但真实 PES 上布满局部极小，多数体系需要**系综（结构集合）**来描述——例如酸性菱沸石上烯烃异构化势垒，系综平均可带来高达 ~20 kJ/mol 的差异（J. Catal. 373, 361 (2019)）；升温时熵甚至会改变反应机理本身（J. Chem. Phys. 131, 214508 (2009)），静态方法无法反映这些。

### 例 10：定义反应坐标

**教学目标**
- 用 MD 检验试探反应坐标；
- 确定试探反应物结构的概率密度。

**反应坐标的定义**
- 反应坐标 $\xi$ 应是原子位置 $\mathbf q$ 的连续函数，且无限缓慢地改变 $\xi$ 时过程可逆（Annu. Rev. Phys. Chem. 67, 669 (2016)）。
- 本例取 $\xi = r_2 - r_1$：$r_1$ 为 CO 的 C 与甲基 C 的距离（成键中），$r_2$ 为甲基 C 与沸石 O 的距离（断键中）。目标过渡态的反应坐标记为 $\xi^*$。

![反应坐标定义示意](/imgs/2026-08-05/fig32.png)

**输入要点**
- 目录 `e10_Dynamic_Coordinates`；POSCAR 为例 9 的反应物 CONTCAR。INCAR 要点（示例文件内容略）：
  - MD：`MDALGO=5`（CSVR 恒温器）、`CSVR_PERIOD`（恒温器时间尺度，以 MD 步数计）、`POTIM`（CSVR 时间步，fs）、`TEBEG`（初始/目标温度 440 K）；
  - MLFF：`ML_LMLFF=T`、`ML_MODE=run`；`ML_OUTPUT_MODE` 减少 MD 文件输出以提速；`ML_OUTBLOCK` 每 n 步输出一次；运行脚本自动加 `RANDOM_SEED`（随机数生成器定义初速度，保证可复现）。
  - KPOINTS 为 $\Gamma$ 点；`ICONST` 文件监测反应坐标：`R` 行定义键长（C(3)–O(25)、C(2)–C(3)），`S` 行为两个 R 之差；末位标志 `0` 表示约束、`7` 表示只监测（monitor）。
- 计算：运行脚本跑 50000 步（分 5 次 MD，共约 8–10 分钟；时间紧可只跑 1 次并用 refdata 补齐）。

**REPORT 文件解读**：开头是 MD 初始条件，随后每 10 步输出一次，字段含义：

| 标签 | 含义 |
|---|---|
| `MDALGO` | 恒温器类型（CSVR 即 5） |
| `CSVR_PERIOD` | CSVR 恒温器时间尺度（步数） |
| `SCALING` | 速度标度（本例未用） |
| `CNEXP` | ICONST `flag=D` 的设置（本例未用） |
| `RANDOM_SEED` | 定义初速度的随机数生成器 |
| `DOF` | 离子自由度 |
| `mc>` | 被监测的反应坐标 |
| `e_b>` | MD 能量：`E_tot`/`E_pot`/`E_kin`，及 Nosé 质量的势能 `EPS`、动能 `ES`（本例未用） |
| `tmprt>` | `T_sim` 初始/目标温度，`T_inst` 当前步瞬时温度 |

**分析流程**
- `monitorVal.sh` 从 report 文件提取 `mc>` 到 `mc.dat`；`normalizedHistogram.py` 画 100 个区间的概率密度直方图到 `mc.dat.hist`；`normalizedKDE.py` 用核密度估计（KDE，平滑宽度 0.05 Å）平滑直方图。
- 结果：概率密度最大处 $P(\xi)\approx 1.05$ Å$^{-1}$ 对应 $\xi\approx 2.4$ Å，即**反应物态**；从采样最好的区域取参考反应物结构 $\xi_{ref,R}$ 可使 $P(\xi)$ 的相对误差最小。
- 重要对比：静态方法给出的反应物 $\xi_{ref,R}=1.95$ Å，动态系综给出 2.40 Å，相差 0.45 Å——直观体现"单一结构 vs 结构分布"的差异。
- 用脚本对保存的 `POSCAR.*` 逐一计算 $\xi$，选出最接近概率峰的结构（参考数据为 POSCAR.2，$\xi\approx 2.32$ Å，$P\approx 0.87$ Å$^{-1}$）供下一步使用。

### 例 11：验证反应坐标（Blue Moon 约束 MD + slow-growth）

**教学目标**
- 用约束 MD 把反应物（R）转变为过渡态（TS）；
- 检验自由能梯度的连续性；
- （可选）检验 TS→R 的可逆性。

**输入要点**
- 目录 `e11_*`；POSCAR/POTCAR/KPOINTS/ICONST 同例 10；INCAR 新增两个标签（示例文件内容略）：
  - `INCREM`：slow-growth 模拟中控制转变速度（即 $\xi$ 每步增量）；取负值使 $\xi$ 从反应物走向过渡态；
  - `LBLUEOUT=T`：执行约束 MD（Blue Moon 系综），把自由能梯度写入 REPORT。
  - 约束 MD 使用 SHAKE 算法。

**计算与分析**
- 从例 10 选定的初始 POSCAR 出发（若出现 `Error: SHAKE algorithm did not converge!` 说明所选结构离反应物态太远，应重选更接近概率密度峰的结构）。
- Blue Moon 模拟分 6 段、每段 4000 步，沿 $\xi$ 扫描自由能梯度（检验其连续性）；结束后得到规则网格上 7 个 $\xi$ 的结构（`POSCAR.1`–`POSCAR.6` + 末态 CONTCAR）。
- `pos2xyzStatic.py` 把 7 个结构拼成 `anim.xyz` 动画（JMol 播放）：可见 CH3 与沸石 O（ZO）成键、CO 在胞内别处；随后 CO 与 CH3 靠近、ZO–CH3 断裂、CH3 翻转，最终形成乙酰阳离子 H3C–CO$^+$——机理符合预期。注意 7 幅图像是各自系综的代表结构。
- `fgrad.sh` 从 report 文件的 `cc>`（$\xi$）与 `b_m>`（自由能梯度）提取数据到 `grad.dat`；画图可见自由能梯度**没有明显不连续**——梯度连续性是对 $\xi$ 选择的第二重验证。
- （可选）可逆性检验：从最终 TS 结构出发，把 `INCREM` 符号从 `-1.1e-4` 改为 `+1.1e-4` 反向 slow-growth；正反向自由能梯度剖面定性相似，说明过程可逆。

### 例 12：反应热力学（动态）与速率常数

**教学目标**
- 用 Blue Moon 模拟得到 R→TS 沿反应坐标的自由能梯度；
- 积分得反应自由能；
- 在动态近似下计算速率常数，并与静态结果对比。

**理论公式**
- 活化自由能：
$$\Delta A^{\ddagger}_{R\to TS} = \Delta A_{\xi_{ref,R}\to\xi^*} - k_B T \ln\left[\frac{h}{k_B T}\frac{\langle|\dot\xi^*|\rangle}{2}P(\xi_{ref,R})\right]$$
  其中 $\Delta A_{\xi_{ref,R}\to\xi^*}$ 为把反应坐标从 $\xi_{ref,R}$ 移到 $\xi^*$ 所做的可逆功，$\langle|\dot\xi^*|\rangle$ 为过渡态构型的反应坐标速度系综平均，$P(\xi_{ref,R})$ 为反应物系综的概率密度。
- 可逆功由自由能梯度积分得到：
$$\Delta A_{\xi_{ref,R}\to\xi^*} = \int_{\xi_{ref,R}}^{\xi^*}\left(\frac{\partial A}{\partial \xi}\right)_{\xi^*} d\xi$$
- Blue Moon 系综的自由能梯度（约束系综平均）：
$$\left(\frac{\partial A}{\partial \xi}\right)_{\xi^*} = \frac{1}{\langle|Z|^{-1/2}\rangle_{\xi^*}} \langle |Z|^{-1/2}[\lambda + G k_B T]\rangle_{\xi^*}$$
  其中 $Z$ 为质量度规张量，$\lambda$ 为 SHAKE 约束的 Lagrange 乘子，$G$ 由 $|Z|$、$Z^{-1}$、$\nabla\xi$、$\nabla|Z|$ 导出。
- 速率常数仍用 Eyring 方程：$k_{R\to P} = \frac{k_BT}{h}e^{-\Delta A^{\ddagger}_{R\to P}/k_BT}$。

**输入与计算**
- 目录 `e12_Dynamic_Thermodynamics`；7 个 `POSCAR.n_init` 为例 11 生成的 7 个 $\xi$ 结构；INCAR 相对例 11 把 `NSW` 从 4000 提高到 10000（每点 10000 步约束 MD，7 点共 ~140 ps）；ICONST/KPOINTS/POTCAR 同例 10（示例文件内容略）。
- 完整计算约 25 分钟；可只跑第 1 个点，其余用 refdata 补齐。

**分析流程与结果**
- 从各段约束 MD 的 REPORT（第 4 列 $|z|^{-1/2}(\lambda+Gk_BT)$ 的系综平均除以第 3 列 $|z|^{-1/2}$ 的系综平均）得到各 $\xi$ 的自由能梯度；`fgradBM.sh` 生成 `grad.dat`（作图时把 $\xi$ 轴反转，使 R→TS 从左到右）。
- 梯度形状：$\xi$ 约 2.5–1.5 Å（反应物系综）区间大致平坦，随后降至 0.5 Å，再快速上升；TS 位于 $\xi\approx 0.0$ Å。
- `simpsonI_non_cum.py` 用 Simpson 法对 `grad.dat` 积分得 `dA.dat`；`splinesCubic.py` 做三次样条插值得 `dA_splines.dat`，即自由能 $A(\xi)$ 曲线（从 ~2.4 Å 反应物处上升至 TS）。
- TS 参考结构的选取：$\langle|Z|^{-1/2}\rangle$ 由约束 MD 系综平均完全决定；因 TS 对应 $\xi=0.043$ Å 并非积分点，`zetBM.sh` 提取各点的 $|z|^{-1/2}$ 作图——其在 TS 附近几乎不变，故用 POSCAR.6 很好地近似 TS。
- `aux_vel.sh` 从 TS 附近结构提取 $\langle|\dot\xi^*|\rangle$（参考值 $1.05285\times 10^{13}$ Å/s）。
- 汇总脚本读取 `dA_splines.dat`（$A_{TS}-A_R$）、`xi_dot.dat`（$\dot\xi^*$）、例 10 的 $\xi_{ref,R}$ 与 $P(\xi_{ref,R})$，代入公式得 $\Delta A^{\ddagger}_{R\to TS}\approx 1.14$ eV，比静态结果高 ~0.1 eV（远大于静态内量子/半经典之差的 0.02 eV）——原因是纳入了系综平均而非单一态。
- 速率常数对比：自由能仅差 ~0.1 eV，速率常数却相差约**一个数量级**——小能量差经指数放大。

### 功能与标签汇总
- MD：`MDALGO=5`（CSVR）、`CSVR_PERIOD`、`POTIM`、`TEBEG`、`RANDOM_SEED`；REPORT 输出字段（`mc>`、`cc>`、`b_m>`、`e_b>`、`tmprt>`）
- 约束/监测：`ICONST`（R=键长、S=组合坐标；0=约束、7=监测）
- Blue Moon / slow-growth：`LBLUEOUT=T`、`INCREM`（±号控制方向）、SHAKE 算法
- MLFF-MD 输出控制：`ML_OUTPUT_MODE`、`ML_OUTBLOCK`
- 分析脚本链：`monitorVal.sh` → `normalizedHistogram.py`/`normalizedKDE.py`（KDE）→ `fgrad.sh`/`fgradBM.sh` → `simpsonI_non_cum.py`（Simpson 积分）→ `splinesCubic.py`（样条）→ `zetBM.sh`/`aux_vel.sh` → Eyring 速率常数
- 概念：反应坐标的连续性与可逆性要求、概率密度定参考反应物态、系综平均对活化自由能/速率常数的影响

---

## 8. 杂化泛函（Hybrid Functionals）

杂化泛函是交换关联泛函的特殊类别，把一定比例的 Fock 交换混入密度泛函，灵感来自对无关联体系有精确解的 Hartree–Fock 方法。

## 8.1 杂化泛函（一）：可用泛函总览（PBE0 / B3LYP / HF / 屏蔽杂化 / HSE06 能带）

> 来源：<https://vasp.at/tutorials/latest/hybrids/part1/>
> 配套输入文件：[hybrids-part1.zip](https://vasp.at/tutorials/latest/hybrids-part1.zip)

本页四个例子：例 1 Si 带隙（PBE vs PBE0），例 2 Ar 带隙（PBE/B3LYP/HF），例 3 MgO 屏蔽杂化泛函带隙优化，例 4 CaS 的 HSE06 能带结构。

**背景**：LDA/GGA/meta-GGA 计算效率最高，但 LDA/GGA 无法正确描述带隙。杂化泛函把一定比例的 Fock 交换混入半局域泛函，以矫正 LDA/GGA 著名的带隙低估。非局域 Fock 势为

$$V_x^{\mathrm{HF}}(\mathbf r,\mathbf r') = -\frac{e^2}{2}\sum_{n\mathbf k} f_{n\mathbf k} \frac{\psi^*_{n\mathbf k}(\mathbf r')\psi_{n\mathbf k}(\mathbf r)}{|\mathbf r-\mathbf r'|}$$

它依赖 KS 轨道而非密度，因此杂化泛函严格说**不是**密度泛函。

### 例 1：Si 带隙（PBE 与 PBE0）

**教学目标**：解释什么是杂化泛函、设置 PBE0 计算、计算带隙。

**PBE0 的构成**（无经验拟合参数）：

$$E_{xc}^{\mathrm{PBE0}} = \frac{1}{4}E_x^{\mathrm{HF}} + \frac{3}{4}E_x^{\mathrm{PBE}} + E_c^{\mathrm{PBE}}$$

**输入要点**（目录 `e01_Si-gap`，示例文件内容略）
- PBE0 用 `LHFCALC` 开启 Fock 交换；只要 `AEXX`、`AGGAC`、`ALDAC` 保持默认，`GGA=PE` 即对应式中的系数（`AEXX`/`AGGAC`/`ALDAC` 即各分量的前置因子）。
- PBE 计算：`ISMEAR=-5`（四面体方法+Blöchl 修正）给出最平滑的 DOS；`LORBIT=11` 把投影 DOS 写入 DOSCAR 与 vaspout.h5。
- PBE0 计算：电子弛豫用阻尼算法（`ALGO=Damped`，直接优化、总能变分、比迭代法稳健，杂化泛函常需此稳健性）。由于 Blöchl 修正使能量非变分，收敛阶段改用高斯展宽 `ISMEAR=0`、`SIGMA=0.01`；收敛后再做单发 `ALGO=None`+`ISMEAR=-5`+`LORBIT=11` 拿 DOS。
- 收敛精度声明：带隙精度在 0.1 eV 内；要提高需对 ENCUT 与 k 点数做收敛测试。

**结果**
- 从 PBE 复制 WAVECAR 到 PBE0 目录再运行（stdout 每步多一行阻尼算法信息）。
- **不同交换关联泛函的总能不可直接比较**，只能比较带隙等物理量。
- PBE 带隙 0.62 eV，PBE0 带隙 1.84 eV：杂化泛函普遍给出更大带隙，部分矫正 GGA 低估。
- PBE0 的 DOS 不仅带隙变大，能带相对能量与带宽也不同于 PBE，这些细微变化会显著影响材料性质。

### 例 2：Ar 带隙（PBE / B3LYP / HF）

**教学目标**：使用 B3LYP 与 HF 方法、理解二者的局限。

**泛函构成**
- B3LYP 是为分子性质设计的半经验杂化泛函（参数拟合自原子化能、电子/质子亲和能、电离势等测试集）：

$$E_x^{\mathrm{B3LYP}} = 0.8E_x^{\mathrm{LDA}} + 0.2E_x^{\mathrm{HF}} + 0.72\Delta E_x^{\mathrm{B88}},\qquad E_c^{\mathrm{B3LYP}} = 0.19E_c^{\mathrm{VWN3}} + 0.81E_c^{\mathrm{LYP}}$$

- HF：全 Fock 交换（`AEXX=1`）、无关联（$E_c=0$）。
- **根本缺陷**：B3LYP 等非屏蔽泛函继承 HF 在均匀电子气极限的失败——Fock 势傅里叶变换在 $\mathbf q=0$ 处奇异，会完全压制巡游态；因此 HF/杂化泛函不推荐用于金属与小带隙半导体等巡游性强的体系（Paier et al. 2007）。

**输入要点**（目录 `e02_Ar-gap`，fcc Ar 晶体，示例文件内容略）
- HF：`LHFCALC` 开启、`AEXX=1`；此时默认 `GGA=PE`、`AGGAX=1-AEXX=0`，且 `AGGAC`/`ALDAC` 自动置零——只有 Fock 交换。
- B3LYP：用 `GGA` 标签切换半局域泛函并设置各交换/关联分量的比例系数。

**结果**：实验带隙 14.4 eV；PBE 8.5 eV（严重低估）、HF 17.6 eV、B3LYP 10.3 eV——三者都偏离实验 3 eV 以上，说明对大带隙绝缘体这些方法同样不精确。

### 例 3：MgO 屏蔽杂化泛函带隙优化

**教学目标**：理解屏蔽杂化泛函基础、HSE06 的基本假设、为复现带隙优化交换比例与屏蔽长度。

**理论**：屏蔽杂化泛函把 Fock 交换的 $1/r$ 分解为短程（SR）与长程（LR），引入屏蔽参数 $\mu$：

$$E_{xc}^{\mathrm{screened}} = \frac{1}{4}E_x^{\mathrm{HF,SR}}(\mu) + \frac{3}{4}E_x^{\mathrm{PBE,SR}}(\mu) + E_x^{\mathrm{PBE,LR}}(\mu) + E_c^{\mathrm{PBE}}$$

- SR 极限即 PBE0；截断 LR 部分避免了 Fock 交换的奇异（正是它压制金属态），并可通过 Hartree-Fock 算子降采样（downsampling）降低计算量。
- HSE06：$\mu=0.2$ Å$^{-1}$（经验优化），作用范围 $2/\mu\approx 10$ Å，约覆盖几个近邻；VASP 中用 `HFSCREEN` 控制。
- 教程强调：为复现实验带隙而拟合 AEXX 牺牲了第一性原理的预测性（更接近模型物理），文献中不应把拟合后的泛函笼统称为 HSE——HSE06 专指 `AEXX=0.25`、`HFSCREEN=0.2`。

**输入与计算**（目录 `e03_MgO-gap-optimization`，示例文件内容略）
- 相关标签：`HFSCREEN`、`HFRCUT`、`HFALPHA`。
- `loop_aexx.sh`：扫描 AEXX 值计算带隙，复现 MgO 实验带隙 7.7 eV 的最优 AEXX 在 **0.45–0.50** 之间。
- `loop_hfscreen.sh`：固定接近最优的 AEXX 再扫描 HFSCREEN——屏蔽长度越短，Fock 交换越多、带隙越大，结果对屏蔽参数敏感。

### 例 4：CaS 能带结构（PBE 与 HSE06）

**教学目标**：用杂化泛函计算能带结构。

**背景**：能带沿连接布里渊区高对称点的高对称线计算；VBM/CBM（HOMO/LUMO）若位于同一 k 点为直接带隙，否则为间接带隙。可用 SeeK-path 工具从 POSCAR 生成 k 路径。

![SeeK-path 生成的 CaS 第一布里渊区](/imgs/2026-08-05/fig33.png)

**输入要点**（目录 `e04_CaS-band`，fcc CaS，示例文件内容略）
- PBE 与 HSE06 的自洽计算都在均匀 k 网格上；PBE 能带用 line-mode KPOINTS 沿路径取点。
- **杂化泛函能带的关键技巧**：既需要高对称路径上的 k 点，又需要均匀 k 网格来评估 Fock 交换能。因此 `HSE06/band/KPOINTS` 包含 16 个带权重的不可约 k 点（张成第一布里渊区，权重按 IBZ 确定）+ 37 个零权重的路径 k 点（与 PBE 能带同路径）。

**计算要点**
- 画能带可用 py4vasp 并绘制能带特征（band character）。
- 比较 `HSE06/band/KPOINTS` 与 `IBZKPT`、PBE 的 band/OUTCAR：路径 k 点相同；杂化计算额外需要均匀网格，故路径点权重为零。
- `NELMIN` 的意义：从已有 WAVECAR 重启时总能已收敛，VASP 在 NELMIN 次迭代后停止更新轨道、直接对角化 KS 哈密顿量取本征值；附加 k 点的轨道若收敛不足会导致能带出现跳变，因此至少保证 4 次更新（`NELMIN`）。
- 导带仍可见轻微折线（kink），可用更大 NELMIN 或更多能带（`NBANDS`）改善。
- 结果：CaS 带隙为**间接带隙**（HOMO 在 $\Gamma$、LUMO 在 X），实验值 5.4 eV。

### 功能与标签汇总
- 杂化泛函开关与配比：`LHFCALC`、`AEXX`、`AGGAX`、`AGGAC`、`ALDAC`、`GGA`
- 屏蔽杂化：`HFSCREEN`（$\mu$）、`HFRCUT`、`HFALPHA`；downsampling 降本
- 收敛与算法：`ALGO=Damped`（变分要求）、`ISMEAR=-5` 与 Blöchl 修正不相容、`ALGO=None` 单发 DOS、`NELMIN`、`NBANDS`
- DOS/能带：`LORBIT=11`、line-mode KPOINTS、零权重路径 k 点技巧、SeeK-path、py4vasp 画 DOS/能带
- 结论性认识：杂化泛函非严格密度泛函；HF/非屏蔽杂化压制巡游态（金属慎用）；拟合 AEXX 牺牲预测性

---

## 9. 线性响应（Linear Response）

除基态性质外，VASP 还能计算体系对外界微扰的响应。微扰足够小时响应处于线性区，可应用线性响应理论。

## 9.1 响应函数（一）：静态与频率依赖的介电响应

> 来源：<https://vasp.at/tutorials/latest/response/part1/>
> 配套输入文件：[response-part1.zip](https://vasp.at/tutorials/latest/response-part1.zip)

本页三个例子：例 1 有限差分法计算 SiC 的静态响应，例 2 DFPT 计算 AlP 静态介电响应，例 3 NaCl 频率依赖介电函数。

### 例 1：有限差分法计算静态响应（SiC）

**教学目标**
- 说出并计算对小电场、应变、离子位移的静态响应；
- 掌握有限差分方法相关的 INCAR 标签与 PEAD 方法；
- 区分"离子冻结（ion-clamped）"与"离子弛豫（relaxed-ion）"版本的极化率、弹性模量与压电张量。

**理论框架**：总能 $E$ 对电场 $\mathcal E$、应变 $\eta$、离子位移 $u$ 做二阶展开：

$$E(\mathbf u,\mathcal E,\eta) = E^0 + \frac{\partial E}{\partial u_i}u_i + \frac{\partial E}{\partial \mathcal E_i}\mathcal E_i + \frac{\partial E}{\partial \eta_{ij}}\eta_{ij} + \frac{1}{2}\frac{\partial^2 E}{\partial u_i\partial u_j}u_iu_j + \cdots + \frac{\partial^2 E}{\partial \mathcal E_i\partial \eta_{jk}}\mathcal E_i\eta_{jk}$$

各系数的物理定义（重复指标求和）：
- 一阶：力 $F_i=-\Omega_0\,\partial E/\partial u_i$；极化 $P_i=-\partial E/\partial\mathcal E_i$（现代极化理论中为所有占据态 Berry 相位之和，$P_i = -\frac{fe}{(2\pi)^3}\sum_n\int_{BZ}d\mathbf k\,\langle u_{n\mathbf k}|\mathrm i\nabla_{\mathbf k}|u_{n\mathbf k}\rangle$）；应力 $\sigma_{ij}=\partial E/\partial\eta_{ij}$。力与应力由 Hellmann–Feynman 定理计算。
- 二阶：力常数矩阵（Hessian）$K_{ij}$、离子冻结介电极化率 $\chi_{ij}=\partial P_i/\partial\mathcal E_j$、零场离子冻结弹性模量 $C_{ijkm}$；
- 混合二阶：Born 有效电荷 $Z^*_{ij}=-(\Omega_0/e)\,\partial^2 E/\partial u_i\partial\mathcal E_j$、内应变张量 $\Lambda_{ijk}$、离子冻结压电张量 $e^{(0)}_{ijk}$。

**输入要点**（目录 `e01_SiC-static-finitedifferences`，示例文件内容略）
- 基态：GGA-PBE、高斯展宽（`ISMEAR`/`SIGMA`）；`ENCUT` 设得较高以避免 Pulay 应力；`EDIFF` 控制电子弛豫停止。
- `electric/`：在小静态电场下计算，`LCALCPOL=T` 计算电子极化（Berry 相位）；`LCALCEPS=T` 启用 PEAD（perturbation expression after discretization）方法，用有限差分算场极化 KS 轨道胞周期部分对 k 的导数：
$$\frac{\partial|u_{n\mathbf k_j}\rangle}{\partial k} = \frac{1}{2\Delta k}\sum_{m=1}^N\left[|u_{n\mathbf k_{j+1}}\rangle S^{-1}_{mn}(\mathbf k_j,\mathbf k_{j+1}) - |u_{n\mathbf k_{j+1}}\rangle S^{-1}_{mn}(\mathbf k_j,\mathbf k_{j-1})\right]$$
  其中 `EFIELD_PEAD` 为步长、`IPEAD` 为有限差分步数；设置时要考虑带隙大小与 k 点数（详见 LPEAD wiki）。
- `electric-ionic/`：在电场基础上加离子位移有限差分（`IBRION`、`POTIM` 步长、`NFREE` 步数）；`ISIF=3` 允许体积与晶胞形状变化（计算弹性模量所必需）。
- KPOINTS 为 $\Gamma$ 中心均匀网格；KPOINTS_OPT 为能带高对称点。

**计算与分析**
- 先算 DFT 基态并用 py4vasp 画能带（确认带隙大小，用于后续设置 EFIELD_PEAD）。
- 在 `electric/` 加小电场：若电场过强会警告可能发生 Zener 隧穿——可减少 k 点数或减小 `EFIELD_PEAD`（下一步即减小一个数量级）。
- `electric-ionic/` 从基态 WAVECAR 启动，可得式 (1.2)–(1.4) 全部静态响应量；py4vasp 可直接打印，OUTCAR 末尾也有同样信息。
- 实验对比需要**离子弛豫版本**的量，它们由离子冻结量加内禀修正得到：
$$\chi_{ij}|_{\eta} = \chi_{ij}|_{u,\eta} + \Omega_0^{-1}Z^*_{mi}(K^{-1})_{mn}Z^*_{nj}$$
$$C_{ijkm}|_{\mathcal E} = C_{ijkm}|_{u,\mathcal E} - \Omega_0^{-1}\Lambda_{nij}(K^{-1})_{no}\Lambda_{okm}$$
$$e^{(0)}_{ijk} = e^{(0)}_{ijk}|_{u} - \Omega_0^{-1}Z^*_{mi}(K^{-1})_{mn}\Lambda_{njk}$$
  VASP 不单独打印极化率，而是直接用它构造介电张量（下一节详述）。

### 例 2：DFPT 静态介电响应（AlP）

**教学目标**
- 把介电常数与极化率联系起来；
- 用 DFPT 计算静态介电性质；
- 解释长波近似与宏观/微观响应的区别；
- 说明 `LEPSILON` 与 `LCALCEPS` 的优缺点。

**理论**
- 介电常数把材料内电场与外加电场耦合：$\mathcal E_i(\mathbf r,\omega)=\int d\mathbf r'\,\epsilon^{-1}_{ij}(\mathbf r,\mathbf r',\omega)\mathcal E_{\mathrm{ext},j}(\mathbf r',\omega)$；动量空间中对 $\mathbf G,\mathbf G'$ 求和。
- **长波近似**：材料均匀或电场变化尺度远大于原子间距时，$\epsilon$ 只依赖 $\mathbf r-\mathbf r'$，宏观介电响应为 $\epsilon^\infty_{ij}(\hat{\mathbf q},\omega)\approx\lim_{\mathbf q\to 0}\epsilon_{00\,ij}(\mathbf q,\omega)$；若局域场效应不可忽略，则取逆矩阵的 $\mathbf q\to 0$ 极限：$1/\epsilon^\infty_{ij}=\lim_{\mathbf q\to 0}\epsilon^{-1}_{00\,ij}(\mathbf q,\omega)$。上标 $\infty$ 表示电子（高频）贡献，是极化激元（polariton）领域的常用记号。
- DFPT 通过线性 Sternheimer 方程求 KS 轨道对电场的线性响应 $|\xi_{n\mathbf k}\rangle$：
$$[H(\mathbf k)-\epsilon_{n\mathbf k}S(\mathbf k)]|\xi_{n\mathbf k}\rangle = -\Delta H_{\mathrm{SCF}}(\mathbf k)|u_{n\mathbf k}\rangle - \hat{\mathbf q}|\nabla_{\mathbf k}u_{n\mathbf k}\rangle$$
  其中 $\Delta H_{\mathrm{SCF}}$ 为局域场效应（微观 KS 哈密顿量变化）；`LRPA=F` 时令其为零。再用第二个 Sternheimer 方程求 $\nabla_{\mathbf k}|u_{n\mathbf k}\rangle$，最终静态宏观介电矩阵为
$$\hat{\mathbf q}\cdot\epsilon^\infty\cdot\hat{\mathbf q} = 1 - \frac{f\,8\pi e^2}{\Omega_0}\sum_{n\mathbf k}w_{\mathbf k}\langle\hat{\mathbf q}\nabla_{\mathbf k}u_{n\mathbf k}|\xi_{n\mathbf k}\rangle$$

**输入与对比要点**（目录 `e02_AlP-static-DFPT`，示例文件内容略）
- `LEPSILON=T` 开介电响应，`IBRION=8` 用微扰论算离子位移与应变（二阶导/Hessian/声子频率）。
- **LEPSILON vs LCALCEPS**：两者都对未占据态无求和；但 LEPSILON（DFPT）**适用于金属**，而 LCALCEPS（PEAD 有限差分）不适用；目前 LEPSILON 不能用于显含轨道的交换关联泛函（HF 类、杂化泛函）。

**计算**：先算 DFT 基态，再在 `electric-ionic/` 计算介电响应；用 py4vasp 打印静态介电响应。

### 例 3：频率依赖介电响应（NaCl）

**教学目标**
- 计算频率依赖介电函数；
- 用 Kramers–Kronig 关系联系介电函数实部与虚部；
- 把介电函数与光学性质联系起来；
- 计入离子对介电响应的贡献、计算声子频率并换算频率单位。

**理论**
- 完整频率依赖介电函数 = 电子贡献 + 离子贡献，两者分别在 DFT 基态之上计算后相加：$\varepsilon(\omega)=\varepsilon_{\mathrm{elec}}(\omega)+\varepsilon_{\mathrm{ion}}(\omega)$。
- 实部与虚部由 Kramers–Kronig 关系相联系；从介电函数可进一步得反射率、吸收、光导率等。
- 电子部分：`LOPTICS=T` 时用二阶微扰论计算 $\nabla_{\mathbf k}|u_{n\mathbf k}\rangle$，公式对**未占据态求和**：
$$\nabla_{\mathbf k}|u_{n\mathbf k}\rangle = \sum_{n'\neq n}\frac{|u_{n'\mathbf k}\rangle\langle u_{n'\mathbf k}|\frac{\partial[H(\mathbf k)-\epsilon_{n\mathbf k}S(\mathbf k)]}{\partial\mathbf k}|u_{n\mathbf k}\rangle}{\epsilon_{n\mathbf k}-\epsilon_{n'\mathbf k}}$$
  因此必须设置足够多的 `NBANDS` 并对 NBANDS 收敛极化率。

**输入要点**（目录 `e03_NaCl-optics`，示例文件内容略）
- 子目录：`electron/`（电子介电函数）、`ion-dfpt/`（DFPT 离子贡献）、`ion-finite-differences/`（有限差分离子贡献）。
- POTCAR 使用 **GW 推荐的赝势**：标准 DFT 中未占据带不贡献基态性质，而增大 NBANDS 时应选用对费米能级以上性质也经过检验的 GW 推荐赝势。

**计算与分析**
- 基态算完后复制 WAVECAR 到 `electron/` 计算电子介电函数；用 `get_dielectric_function.sh` 从 vasprun.xml 提取介电函数作图，或用 py4vasp 画。
- 离子贡献两种方法（DFPT 快、有限差分慢得多）结果对比。
- 离子位移给出振动方程 $m_{\mathrm{ion}}\omega^2\mathbf u = (\partial^2E/\partial u_i\partial u_j)\mathbf u$；从 OUTCAR 读声子频率。
- **偶极活性声子**：若某声子模式使涉及离子的电偶极矩在一个周期内变化，则它出现在频率依赖介电函数中（本例有 3 个偶极活性模式）。
- 另有 3 个虚频（f/i）模式——此处不代表动力学不稳定，而是有限平面波基组导致的数值效应（总能随全体离子质心在胞内位置略有变化）。
- 质量加权：动力学矩阵本征矢按质量加权，Na 与 Cl 位移矢量平行但长度不同，长度比恰为 $\sqrt{m_{\mathrm{Na}}/m_{\mathrm{Cl}}}$（质量取自 POTCAR 的 `POMASS`）。
- **频率单位换算**：VASP 声子频率单位 $\sqrt{\mathrm{eV/Å^2/a.m.u.}} = 9.822517\times 10^{13}$ s$^{-1}$，换算得
$$1\,\mathrm{THz} = 4.1357\,\mathrm{meV} = 33.356\,\mathrm{cm}^{-1}$$
$$1\,\mathrm{meV} = 0.242\,\mathrm{THz} = 8.066\,\mathrm{cm}^{-1}$$
$$1\,\mathrm{cm}^{-1} = 0.030\,\mathrm{THz} = 0.124\,\mathrm{meV}$$
- 离子介电函数的最大频率由声子频率自动决定（取 $1.2\times\omega_{max}$）。最后把 DFPT 离子贡献以 meV 为频率轴作图（脚本或 py4vasp）。

### 功能与标签汇总
- 有限差分响应：`LCALCPOL`、`LCALCEPS`（PEAD）、`EFIELD_PEAD`、`IPEAD`、`LPEAD`、`IBRION`+`POTIM`+`NFREE`、`ISIF=3`；Zener 隧穿警告的处理
- DFPT 响应：`LEPSILON=T`、`IBRION=8`、`LRPA`；Sternheimer 方程、局域场效应、长波近似、$\epsilon^\infty$ 记号
- 光学：`LOPTICS=T`（对未占据态求和→NBANDS 收敛）、GW 推荐赝势、`get_dielectric_function.sh`、py4vasp 介电函数绘图
- 离子贡献与声子：偶极活性模式、质量加权本征矢、`POMASS`、频率单位换算、最大频率 $1.2\omega_{max}$
- 概念：ion-clamped vs relaxed-ion 响应量及其 $K^{-1}$ 修正关系

---

## 10. GW 近似（GW Approximation）

GW 近似是多体微扰理论中的一种近似，自能取 Green 函数 G 与屏蔽库仑相互作用 W 之积。GW 通过准粒子能量给出体系的谱性质，是目前计算带隙最精确的多体方法之一。

## 10.1 GW 近似（一）：入门（$G_0W_0$、$GW_0$ 与 Wannier90 能带）

> 来源：<https://vasp.at/tutorials/latest/gw/part1/>
> 配套输入文件：[GW-part1.zip](https://vasp.at/tutorials/latest/GW-part1.zip)

本页含 GW 简介（第 0 节）、例 1 Si 的 $G_0W_0$ 带隙、例 2 用 Wannier90 得到 $GW_0$ 能带结构。

**GW 简介**：GW 近似是处理电子关联的高级方法。电子不再被视为独立粒子，而是在材料中运动时与其他电子相互作用、极化其环境——被自身诱导的极化云"着装（dressed）"的电子称为**准粒子**。准粒子态由准粒子方程求解：

$$(T+V_{ext}+V_h)\psi_{n\mathbf k}(\mathbf r) + \int d\mathbf r'\,\Sigma(\mathbf r,\mathbf r',E_{n\mathbf k})\psi_{n\mathbf k}(\mathbf r') = E_{n\mathbf k}\psi_{n\mathbf k}(\mathbf r)$$

与 KS 方程相比，交换关联效应由非局域、轨道与能量依赖的**自能** $\Sigma$ 引入。GW 中取 $\Sigma=GW$（Green 函数 $G$ 与屏蔽库仑作用 $W$ 之积）；$W$ 的介电屏蔽通常在随机相位近似（RPA）下计算——电子与环境的相互作用被严格限定为极化事件。原则上 $G$ 和 $W$ 都依赖准粒子本征态/本征能，故方程是自洽问题，但实践中几乎不做全自洽，常见两类做法：
- **$G_0W_0$（单发 GW）**：$G$、$W$ 均由 KS 本征态/本征能构造，准粒子方程只解一次；
- **$GW_0$（部分自洽，只对 Green 函数）**：首次单发后，用准粒子本征能（但仍用初始 KS 本征态）重建 $G$，再解新方程，重复若干次。

### 例 1：Si 的 $G_0W_0$ 带隙

**教学目标**：运行单发 GW（$G_0W_0$）、获得准粒子能量与重整化因子。

$G_0W_0$ 分三步：① DFT 基态；② 从自洽解（WAVECAR）重启，计算未占据态并写 WAVEDER；③ 以 ② 的 WAVECAR/WAVEDER 为起点做 $G_0W_0$。

**第 1 步：DFT 基态**（目录 `e01_Si-G0W0`，示例文件内容略）
- GGA-PBE；`ISMEAR=-5`（四面体+Blöchl 修正）得平滑 DOS；KPOINTS 为 6x6x6 $\Gamma$ 中心网格；KPOINTS_OPT 为高对称线（用于非自洽能带）。
- **POTCAR 要用 `_GW` 后缀版本**：在更宽能量范围内再现原子散射性质，计算大量未占据态时需要（可查 POTCAR 的 TITEL 字段）。
- 结果讨论：从 DOS/HOMO–LUMO 差与能带图得到的带隙可能不一致——Si 是间接带隙（HOMO 在 $\Gamma$，LUMO 在 $\Gamma$–X 线上略偏离 X 处），该点不在规则 k 网格上，只有 KPOINTS_OPT 的路径覆盖到。

**第 2 步：未占据 KS 轨道**（`unoccupied-states/`）
- 需要约 3–4 倍于默认数量的 KS 轨道（`NBANDS`）；`ALGO=Exact`+`NELM=1` 做一次精确对角化；`LOPTICS` 开启后 Bloch 态对 k 的导数写入 WAVEDER。
- 副产品：可用 py4vasp 画电子介电函数——注意这是独立粒子近似（IPA），不含局域场效应。

**第 3 步：$G_0W_0$**（`single-shot/`）
- `ALGO=EVGW0`（eigenvalue GW）+ `NELMGW=1`（单次电子更新）；`NBANDS` 必须与生成 WAVECAR/WAVEDER 时一致。
- OUTCAR 中 `QP-energies` 与重整化因子 `Z`：第 1 列为 KS 轨道带指标，第 3 列为准粒子能量 $E_{n\mathbf k}$，第 2/4/5/7 列分别为 KS 能量 $E^{(0)}_{n\mathbf k}$、自能对角元 $\langle\psi^{(0)}|\Sigma(\omega=E^{(0)})|\psi^{(0)}\rangle$、交换关联势能与重整化因子 $Z_{n\mathbf k}$。
- 结果：HOMO 在 $\Gamma$（k 点 1）、LUMO 在 X（k 点 13），$G_0W_0$ 带隙 $=6.2826-5.1235=1.1591$ eV。
- 对比 IPA 与 RPA 介电函数（`get_dielectric_function.sh` + 脚本或 py4vasp）：RPA 计入电子对环境的极化；VASP 在 GW 中默认用 RPA 计算屏蔽库仑作用，可观察 $\omega=0$ 处实部/虚部的变化。

### 例 2：Si 的 $GW_0$ 与 Wannier90 能带

**教学目标**：运行部分自洽 $GW_0$；用 Wannier90 近似 GW 能带结构。

**$GW_0$ 计算**（目录 `e02_Si-GW0-band`）
- INCAR 与 $G_0W_0$ 类似：`ALGO=EVGW0`，但 `NELMGW>1`（多次电子更新）；`NBANDS` 与上一步一致；`NBANDSGW=16` 增加计算准粒子态的数目，便于画能带（示例文件内容略）。
- 需要例 1 第 2 步的 WAVECAR/WAVEDER。
- OUTCAR 对每次本征能迭代都写 QP 能量：第一次迭代第 2 列是 KS 本征能，之后各次是上一轮的准粒子能量。
- 带隙对比：第一次迭代（$G_0W_0$）1.1591 eV；三次迭代后（$GW_0$）$6.2999-5.0779=1.2220$ eV。

**Wannier90 插值能带**（`band-structure/`）
- 局限：GW 中 VASP 只能算规则网格上的准粒子能量，不能沿任意高对称线取点；可用 Wannier90 近似。
- 本质是**近似能带**：对 Monkhorst-Pack 网格上的 $GW_0$ 准粒子能量（存于 WAVECAR）做 Fourier/Wannier 插值得到任意 k 点的能量。
- INCAR 要点：`ALGO=NONE`+`NELM=1` 让 VASP 读入后立即进入后处理；`LWANNIER90=T` 让 VASP 写出 Wannier90 输入文件；`NUM_WANN` 指定 Wannier 函数个数；`WANNIER90_WIN` 的内容原样复制到 `wannier90.win`；NBANDS 与 $GW_0$ 计算一致。
- 流程：VASP 后处理生成 `wannier90.*` 文件 → 调用 Wannier90 构造最大局域 Wannier 函数（检查 `wannier90.wout` 末尾确认成功）→ 脚本生成 `band.png`。
- 若存在 `e01_Si-G0W0/vaspout.h5`，可把 DFT 能带叠加到 GW 能带图上（`band_DFT_vs_GW.png`）。

### 功能与标签汇总
- 三步流程：DFT 基态 → `ALGO=Exact`+`NELM=1`+`LOPTICS`（多未占据带、写 WAVEDER）→ `ALGO=EVGW0`
- 关键标签：`NBANDS`（3–4 倍默认）、`NELMGW`（1=$G_0W_0$，>1=$GW_0$）、`NBANDSGW`、`ISMEAR=-5`
- 赝势：`_GW` 后缀 POTCAR
- 输出：OUTCAR 的 QP-energies/Z（列含义）、RPA 屏蔽 W、IPA vs RPA 介电函数对比
- Wannier90：`LWANNIER90=T`、`NUM_WANN`、`WANNIER90_WIN`、`ALGO=NONE`+`NELM=1` 后处理、最大局域 Wannier 函数插值 GW 能带

---

## 11. Bethe–Salpeter 方程（BSE）

VASP 可基于 Bethe–Salpeter 方程计算含激子效应的频率依赖介电函数，可在 DFT 基态、杂化泛函或 GW 近似之上进行。

## 11.1 BSE（一）：金刚石碳的光吸收（GW+BSE、TDDDH）

> 来源：<https://vasp.at/tutorials/latest/bse/part1/>
> 配套输入文件：[BSE-part1.zip](https://vasp.at/tutorials/latest/BSE-part1.zip)

本页含 BSE 形式理论简介（第 0 节）、例 1 金刚石 $G_0W_0$ 准备计算、例 2 BSE 光吸收、例 3 TDDDH 光吸收。

**BSE 形式理论**：Bethe–Salpeter 方程在介电函数中纳入激子效应。激发能是线性问题

$$\begin{pmatrix}\mathbf A & \mathbf B\\ \mathbf B^* & \mathbf A^*\end{pmatrix}\begin{pmatrix}\mathbf X_\lambda\\ \mathbf Y_\lambda\end{pmatrix} = \omega_\lambda\begin{pmatrix}-\mathbf 1 & \mathbf 0\\ \mathbf 0 & \mathbf 1\end{pmatrix}\begin{pmatrix}\mathbf X_\lambda\\ \mathbf Y_\lambda\end{pmatrix}$$

的本征值 $\omega_\lambda$。矩阵 $A$、$B$ 描述占据态 $v,v'$ 与未占据态 $c,c'$ 之间的跃迁：

$$A_{vc}^{v'c'} = (\varepsilon_v-\varepsilon_c)\delta_{vv'}\delta_{cc'} + \langle cv'|V|vc'\rangle - \langle cv'|W|c'v\rangle,\qquad B_{vc}^{v'c'} = \langle vv'|V|cc'\rangle - \langle vv'|W|c'c\rangle$$

相互作用含裸库仑 $V$ 与屏蔽库仑 $W$；能量与轨道通常取自 GW 计算。$A$ 的三项分别称对角项 $D$、交换项 $K^x$、直接项 $K^D$。$W$ 的屏蔽在 RPA 介电函数中计算，实践上取静态近似 $W(\omega=0)$ 为标准。**Tamm–Dancoff 近似（TDA）**忽略激发-退激发耦合（$B$、$B^*$），退化为 $AX_\lambda=\omega_\lambda X_\lambda$。宏观介电函数由 BSE 本征值/本征矢构造；光学极限 $q\to 0$ 下 $e^{i\mathbf q\mathbf r}\approx i\mathbf q\mathbf r$，矩阵元正比于 $\mathbf q\langle c\mathbf k|\nabla|v\mathbf k\rangle$。

**与 TDDFT 的联系**：Casida 形式的 TDDFT 给出同样的线性问题，但电子-空穴相互作用由交换关联核 $f_{xc}=\delta v_{xc}/\delta\rho$ 描述；含杂化泛函时 $\langle cv'|f_{xc}|c'v\rangle = \gamma\langle cv'|V|c'v\rangle + (1-\gamma)\langle cv'|f_x^{ALDA}|vc'\rangle + \langle cv'|f_c^{ALDA}|vc'\rangle$，$\gamma$ 由 Fock 交换比例决定。

**VASP 文件需求**：BSE/TDHF 需要前置基态计算的准粒子轨道 $\psi_{n\mathbf k}$、轨道导数 $\nabla_{\mathbf k}\psi_{n\mathbf k}$ 与能量（存于 WAVECAR 与 WAVEDER）；BSE 还需要 GW 产生的静态屏蔽势 $W_{\mathbf G,\mathbf G'}(\omega=0)$（写于 WFULLxxxx.tmp / Wxxxx.tmp）。

### 例 1：金刚石的准备性 $G_0W_0$ 计算

**教学目标**：设置单发 $G_0W_0$、找 QP 带隙、画 GW 计算中的 RPA 介电函数。

- 三步流程同 GW 教程：① PBE 基态（目录 `e01_C_GW/ground-state`）；② 精确对角化算未占据带（`unoccupied-states/`，含 60 个空轨道，足以使带隙收敛到 0.1 eV）；③ $G_0W_0$（50 个频率点，小于默认但满足 0.1 eV 收敛标准）。k 点较少但带隙收敛在 0.1 eV 内；小胞用精确对角化更快；POTCAR 用 `_GW` 版本（示例文件内容略）。
- 带隙概念：**直接带隙**为同一 k 点最高占据与最低未占据带之差；**基本带隙**为 VBM 与 CBM（可位于不同 k 点）之差。
- 金刚石实验基本带隙 5.85 eV（Phys. Rev. Materials 2, 073803 (2018)，已做零声子重整化 ZPR 修正——ZPR 会"打开"实验带隙，计算中不含，故比较时需修正实验值）；PBE 明显低估。
- 检查 OUTCAR 确认计算的带数正确、WAVEDER 存在；`LOPTICS` 开启时 VASP 在独立粒子近似下计算介电函数（py4vasp 可提取）。
- $G_0W_0$ 结果：直接带隙 7.40 eV、基本带隙 5.58 eV（接近实验 5.85 eV）；再画 RPA 介电函数。

### 例 2：BSE 光吸收

**教学目标**：计算含激子效应的光吸收谱。

- 具备 WAVECAR、WAVEDER、W*.tmp 后，在 TDA 下做 BSE：先用 QP 轨道与 RPA 屏蔽势构造 BSE 哈密顿量，解本征值/本征矢，最后按式 (9) 算 BSE 介电函数。
- INCAR 要点（目录 `e02_C_BSE`，示例文件内容略）：`ALGO=BSE`；对全部 4 个占据态与 4 个空态构造哈密顿量（最高跃迁能量约 37 eV）；`NBANDS` 必须与前置 GW 计算一致，并读取 WAVECAR/WAVEDER。
- 输出检查：vasprun.xml 含光学跃迁与振子强度、BSE 介电函数（AWK 脚本可提取 xx/yy/zz 平均）；OUTCAR 中 `Bands included in the BSE spin= 1 VB(min)= 1 VB(max)= 4 CB(min)= 5 CB(max)= 8` 给出参与能带，`BSE (scaLAPACK) attempting allocation of ... rank= 3456` 给出 BSE 矩阵秩。
- **排错提示**：若 BSE 介电函数全为零，检查 stdout 中是否有 `the WAVECAR file was read successfully`。
- 对比 BSE 与 RPA/IP 谱，观察电子-空穴相互作用对介电常数的影响。

### 例 3：TDDDH 光吸收（免 GW）

**教学目标**：确定模型介电函数参数；用含激子效应的 TDDFT 计算光吸收。

- 动机：GW 步很昂贵；用杂化泛函方法可达到相近精度。**介电依赖杂化泛函（DDH）**中 Fock 交换的屏蔽由模型函数描述：
$$\epsilon^{-1}(|\mathbf q+\mathbf G|) = 1 - (1-\epsilon_\infty^{-1})e^{-|\mathbf q+\mathbf G|^2/4\mu^2}$$
  同一函数也可屏蔽电子-空穴吸引库仑作用，即 **TDDDH**（Phys. Rev. Research 2, 032019 (2020)）。
- 两步流程（目录 `e03_C_TDHF`）：① DDH 基态（从 DFT WAVECAR 重启加速杂化收敛）；② DDH+精确对角化算未占据带。
- 参数拟合：把模型介电函数拟合到 GW 计算 vasprun.xml 中的 RPA 介电函数，得 $\varepsilon=5.7$、$\mu=1.79$（示例文件内容略）。
- 结果：DDH 直接带隙 7.5 eV、基本带隙 5.6 eV（与 GW 相当）；TDDDH 谱与 BSE 谱**高度吻合**——确认 BSE 哈密顿量秩与例 2 相同（3456）。
- 稀疏 k 网格下看不出 RPA 与含 e-h 作用谱的差异；激子效应的重要性参见 Phys. Rev. Lett. 107, 186401 (2011) 图 2。
- 物理结论：金刚石介电屏蔽强（$\varepsilon\sim 5.7$），激子结合弱——**Wannier–Mott 激子**；第二部分将研究弱屏蔽、强束缚 **Frenkel 激子**的 LiF。

### 功能与标签汇总
- 流程文件链：WAVECAR + WAVEDER（`LOPTICS` 产生轨道导数）+ W*.tmp（GW 的静态屏蔽势）→ `ALGO=BSE`
- 近似层级（INCaR 开关）：IP（`LHARTREE=.FALSE.`、`LADDER=.FALSE.`）、RPA（仅交换）、TDA（`LHARTREE=.TRUE.`+`LADDER=.TRUE.`，忽略 B/B*）、full BSE（`ANTIRES=2`）
- 参数：NBANDS 与 GW 一致、频率点数、空轨道数、OUTCAR 的 BSE 能带范围与矩阵秩
- TDDDH：`LMODELHF`/DDH、模型介电函数参数（$\epsilon_\infty$、$\mu$）拟合自 RPA 介电函数
- 输出：vasprun.xml（跃迁、振子强度、介电函数）、py4vasp 提取

---

## 11.2 BSE（二）：LiF 光吸收——逐层近似对比（IP / RPA / TDA / full BSE）

> 来源：<https://vasp.at/tutorials/latest/bse/part2/>
> 配套输入文件：[BSE-part2.zip](https://vasp.at/tutorials/latest/BSE-part2.zip)

与金刚石不同，LiF 的电子间相互作用屏蔽很弱，导致强束缚的 **Frenkel 激子**。本页用同一 $G_0W_0$ 起点，逐级开启 BSE 哈密顿量的各项，对比四个近似层级的光吸收谱：例 4 准备性 $G_0W_0$、例 5 独立粒子近似（BSE-IP）、例 6 BSE-RPA、例 7 BSE-TDA、例 8 超越 TDA 的完整 BSE。

### 例 4：LiF 的准备性 $G_0W_0$ 计算

- 三步流程与第一部分相同（目录 `e04_LiF_GW`：ground-state、unoccupied-states、GW；POTCAR 为 Li/F 的 `_GW` 版本，示例文件内容略）。
- 检查 OUTCAR：LiF 直接带隙约 **13.44 eV**。

### 例 5：独立粒子近似（BSE-IP）

- 设置：`LHARTREE=.FALSE.`、`LADDER=.FALSE.`——BSE 哈密顿量只剩对角项（$A=D$），粒子间相互作用被忽略。对屏蔽强、e-h 结合能小的体系，IP 近似通常已能合理描述光吸收（示例文件内容略）。
- 结果：与实验谱吻合差——吸收边（13.4 eV）处**没有强激子峰**，即电子与空穴未形成束缚态（束缚态本应降低跃迁能量）。
- 该谱等价于用相同轨道与能量做 `LOPTICS=.TRUE.` 的结果。

### 例 6：BSE-RPA

- 设置：`LADDER=.FALSE.`（不算梯子图）、`LHARTREE=.TRUE.`（只含交换作用），$A=D+K^x$（示例文件内容略）。
- 结果：与实验仍不吻合，且谱**略微蓝移**——来自排斥的交换作用。
- 注意：此处 RPA 只在 TDA 下含 bubble 图（忽略式 (1) 的 $B$、$B^*$），因此**不等价于 `ALGO=CHI`**（后者含激发-退激发 bubble 的耦合）。
- RPA 与 IP 有类似缺点（缺激子），但它计入与等离子激元激发的相互作用，常能很好描述 EELS 谱。

### 例 7：BSE-TDA

- 物理背景：LiF 是宽带隙、小介电常数（$\varepsilon_\infty=1.9$）的弱屏蔽体系，e-h 吸引强，激子强局域、结合能大（0.1–1 eV）；强 e-h 作用是碱卤化物与有机分子体系的典型特征。
- 设置：`LHARTREE=.TRUE.`+`LADDER=.TRUE.`（bubble 图与梯子图都计入），但忽略 $B$、$B^*$，哈密顿量 $A=D+K^x+K^D$（示例文件内容略）。
- 结果：与实验吻合**大幅改善**，第一个强激子峰被正确捕捉；部分特征强度仍与实验不符——第三部分将看到加密 k 网格可解决。

### 例 8：超越 TDA 的完整 BSE

- 求解完整 BSE，包含共振-反共振耦合 $H^{BSE}=\begin{pmatrix}A & B\\ B^* & A^*\end{pmatrix}$。
- `ANTIRES=2` 时纳入 $B$、$B^*$，原则上需对角化秩为 $2N_cN_vN_k$ 的矩阵；VASP 利用时间反演对称性把问题拆成两个秩 $N_cN_vN_k$ 的方程——**完整 BSE 只比 TDA 贵约一倍**（示例文件内容略）。
- 结果：用 py4vasp 画完整 BSE 介电函数并与 TDA 对比。

### 功能与标签汇总
- 近似层级开关矩阵：

| 近似 | LHARTREE | LADDER | ANTIRES | 哈密顿量 |
|---|---|---|---|---|
| IP | .FALSE. | .FALSE. | — | $A=D$ |
| RPA（TDA 内） | .TRUE. | .FALSE. | — | $A=D+K^x$ |
| TDA | .TRUE. | .TRUE. | — | $A=D+K^x+K^D$ |
| full BSE | .TRUE. | .TRUE. | 2 | 含 $B,B^*$ 耦合 |

- 认识：激子效应只由梯子图（直接项 $K^D$）带来；交换项导致蓝移；RPA-TDA $\neq$ ALGO=CHI；ANTIRES=2 借助时间反演对称性仅约 2 倍开销
- 物理：弱屏蔽→强束缚 Frenkel 激子（LiF、碱卤化物、有机分子）；强屏蔽→弱束缚 Wannier–Mott 激子（金刚石）

---

## 11.3 BSE（三）：高效布里渊区采样与激子分析

> 来源：<https://vasp.at/tutorials/latest/bse/part3/>
> 配套输入文件：[BSE-part3.zip](https://vasp.at/tutorials/latest/BSE-part3.zip)

### 例 9：高效布里渊区采样（平移网格）

**教学目标**：用平移网格（shifted grids）加速收敛的光吸收计算。

**动机与方法**：BSE/TDHF 需要对秩 $N_k\times N_v\times N_c$ 的矩阵对角化，k 点多时极其昂贵。本例用 8 次粗网格 $4\times4\times4$ 计算近似等效的 $16\times16\times16$ 稠密网格 TDDDH 谱，步骤：
1. 用 $4\times4\times4$ 网格的 DFT 计算确定不可约 k 点 $\mathbf k_n^{ir}$ 与权重 $W_n$（收集到 POINTS 文件）；
2. 对每个 $\mathbf k_n^{ir}$ 做一次平移后的 TDDDH 计算并提取介电函数；
3. 按权重 $W_n$ 加权求和得到最终介电函数。

**输入与计算**（目录含 INCAR.DFT/INCAR.DDH/INCAR.TDHF，示例文件内容略）
- 运行脚本一次跑完 8 个 TDDDH 计算（约 10 分钟）。
- **局限**：激子结合能随相邻 k 点间距 $\mathbf q=\mathbf k-\mathbf k'$ 线性收敛，多网格平均对谱是合理近似，但会**高估激子结合能**——结合能必须用显式布里渊区采样。
- 另一常用技巧：小的非对称平移 $\mathbf s=\{1/12,3/12,5/12\}$ 打破规则网格对称性，包含最多不等价点，改善谱的 k 收敛。
- 结果：LiF 有 1.15 eV 的大零声子重整化（ZPR），计算未含；为便于比较把 BSE 谱平移 1.15 eV，总体与实验吻合良好。

### 例 10：激子分析（BSEFATBAND）

**教学目标**：检查特定激子态的性质；用 fatband 可视化 BSE 本征矢。

**方法**：VASP 可把最低 `NBSEEIG` 个本征矢写入 `BSEFATBAND` 文件；fatband 图把耦合系数以圆圈大小画在对应能带与 k 点上。

**任务 10.1：DDH 基态**（目录 `e10_C_FATBAND/DFT` 与 `DDH`）
- 用长程分离的介电依赖杂化泛函（range-separated DDH）做基态：① DDH 基态；② DDH+精确对角化算未占据带（4 个空轨道，写 WAVEDER）（示例文件内容略）。
- 确认金刚石的直接/基本带隙。

**任务 10.2：TDDDH 并写激子本征矢**
- INCAR 设 `NBSEEIG=10` 写出前 10 个激子的本征矢（示例文件内容略）。
- `BSEFATBAND` 文件结构：`E_BSE`/`E_IP` 为 BSE 与 IP 跃迁能量，`Kx Ky Kz` 为 k 点坐标，`E_h`/`E_e` 为空穴/电子本征值，`X_BSE` 为本征矢分量，`W_k` 为 k 点权重，`NB_h NB_e` 为空穴/电子轨道编号。
- 找第一个**亮激子**：vasprun.xml 中第一个振子强度非零的跃迁是第 4 个态，故画 `4BSE` 本征态的 fatband（gnuplot 或脚本），横轴取 $\Gamma$–L 与 $\Gamma$–X 方向。
- 由该跃迁的 IP 能量与 BSE 本征值之差得第一亮激子结合能 **0.42 eV**（稀疏 k 网格下明显高估）。
- 结论：最大贡献来自 $\Gamma$ 点直接带隙附近的带顶，相邻 k 点贡献很小——该激子态在 **k 空间高度局域**。

### 功能与标签汇总
- 高效采样：不可约 k 点+权重（POINTS）→ 多个平移网格 TDDDH 计算 → 加权求和；非对称平移 $\{1/12,3/12,5/12\}$
- 局限：多网格平均会高估激子结合能（线性收敛于 k 间距）；ZPR 需手工平移谱对比
- 激子分析：`NBSEEIG`、`BSEFATBAND`（字段含义）、fatband 可视化、亮激子（振子强度非零）、k 空间局域性

---

## 12. X 射线吸收谱（XAS）

X 射线吸收谱涉及芯电子激发到导带，留下芯空穴与可形成束缚态（激子）的激发电子。VASP 有两种建模途径：基于 DFT 的超胞核空穴（SCH）方法，以及基于多体微扰理论、需求解 BSE 的方法。

## 12.1 XAS（一）：LiCl 的 X 射线吸收谱（超胞核空穴 / 全核空穴）

> 来源：<https://vasp.at/tutorials/latest/xas/part1/>
> 配套输入文件：[XAS-part1.zip](https://vasp.at/tutorials/latest/XAS-part1.zip)

本页两个例子：例 1 用超胞核空穴（SCH）方法计算 LiCl 中 Li 的 K 边 XAS，例 2 全核空穴（FCH）与激发电子核空穴（XCH）方法及与实验的对比规范。

**XAS 理论背景**：X 射线把电子从芯态激发到导带。X 射线波长远大于固体特征动量，故取长波极限（$\mathbf q=0$）介电函数虚部，它正比于吸收谱：

$$\epsilon^{(2)}_{\alpha\beta}(\omega,\mathbf q=0) = \frac{4\pi^2 e^2\hbar^2}{\Omega\omega^2 m_e^2}\sum_{c,v,\mathbf k}2w_{\mathbf k}\delta(\varepsilon_{c\mathbf k}-\varepsilon_{v\mathbf k}-\omega)M_\alpha^{v\to c}{M_\beta^{v\to c}}^*$$

PAW 形式下的动量矩阵元含全电子部分波 $|\phi_i\rangle$、赝部分波 $|\tilde\phi_i\rangle$ 与投影算符 $|\tilde p_i\rangle$。XAS 只考虑单一格点的核空穴，利用完备性关系可把矩阵元简化为单格点形式，并对带求和限制到导带 $c$：

$$M^{\mathrm{core}\to c\mathbf k}_\alpha = \sum_i\langle\tilde\psi_{c\mathbf k}|\tilde p_i\rangle\langle\phi_i|\nabla_\alpha|\phi_{\mathrm{core}}\rangle$$

### 例 1：超胞核空穴计算

**教学目标**：设置超胞核空穴计算；对超胞大小做收敛；理解激子效应在 XAS 中的重要性。

**输入要点**（目录 `e01_SCH`，示例文件内容略）
- 含核空穴的原子必须作为**单独的元素（species）**且只含 1 个原子：如 8 Li + 8 Cl 的 LiCl 算 Li K 边时变为 1+7+8（Li* Li Cl）三种元素，因此需要第三份 POTCAR（用 `cat` 合并）。
- 采用**终态近似**（`ICORELEVEL=2`，Phys. Rev. B 63, 205419 (2001)）；`CLNT` 选择做核空穴的元素（该元素只放 1 个原子，保证超胞中只有一个核空穴）；`CLN`/`CLL`/`CLZ` 指定激发芯态的主量子数、角量子数与电子数；`CH_NEDOS` 为介电函数能量网格点数；`CH_AMPLIFICATION` 用于标度谱（取超胞中 Li 原子数，补偿谱被超胞体积除导致的强度下降，使不同胞大小的谱可比）；`NBANDS` 需保证有足够未占据带可激发（教程取占据带数的两倍）。
- KPOINTS 只用单个 k 点（也需收敛）；提供 1x1x1 到 4x4x4 的 POSCAR。

**超胞收敛**：`run_conv.sh` 依次跑 2x2x2（16 原子）、3x3x3（54 原子）、4x4x4（128 原子）：峰的劈裂与相对强度随胞大小变化。**超胞的首要目的是最小化核空穴与其周期镜像的相互作用**；增大超胞同时改善 k 点表示，但一旦空穴间作用足够小，直接加密 k 点比继续加大超胞更便宜。

**核空穴效应可视化（激子效应）**
- 对比终态近似（`ICORELEVEL=2`）与**初态近似**（`ICORELEVEL=1`，价电子不随核空穴弛豫、无激子效应）：终态谱整体移向低能、态更局域（第一个峰尤其明显）——皆源于电子-空穴相互作用。
- 用 PARCHG（`INCAR.step1/step2_CH_PARCHG` 与 noCH 版本）画最低导带的部分电荷密度：有核空穴时电荷局域在核空穴原子 (0.5, 0.5, 0.5) 周围；无核空穴时第一未占据态完全离域。

### 例 2：全核空穴（FCH）计算与实验对比

**教学目标**：做 FCH 计算；掌握计算谱对实验的呈现方式。

**方法对比**
- 超胞核空穴法中，一个电子从芯态移出并加入导带——选 `ICORELEVEL=2`、`CLZ=1` 时自动完成，即**激发电子核空穴（XCH）**方法（Phys. Rev. Lett. 96, 215502 (2006)）。
- **全核空穴（FCH）**把激发电子放回均匀背景电荷，自洽中不含激发电子分布（J. Chem. Phys. 120, 8632 (2004)）；许多情形下 FCH 比 XCH 更接近实验（Unzog et al., Phys. Rev. B 106, 155133 (2022) 第 V 部分给出原因）。

**计算流程**（目录 `e02_FCH`，示例文件内容略）
- FCH 需要把激发电子当作均匀背景电荷：VASP 会自动用背景电荷补偿带电胞，因此只需用 `NELECT` 把电子数减一。先用 dry-run 模式查 OUTCAR 的 `NELECT = 513.0000 total number of electrons`，再写入 INCAR.FCH。
- `run_job.sh` 在 4x4x4 超胞、单 k 点上同时跑 FCH 与 XCH。

**与实验对比的规范**
- 芯电子能量总是与实验有系统性偏差，谱需要**平移**：常把实验第一峰移到零、计算谱平移使第一峰极大值对齐（选哪个峰是任意的）。
- 强度常以任意单位报告（谱仪展宽等参数未知）：把计算谱标度到与实验第一峰等高。即看 XAS 的"指纹"而非绝对位置与强度。
- 教程用较小展宽 0.1 eV，峰比实验窄得多；为便于对比，把 `CH_SIGMA` 改大后**复用已有 WAVECAR 重算介电函数**（不重跑整个 SCF），再作图对比。
- 展宽类型（Gaussian、Lorentzian、常数、能量依赖等）可自由选择——实验展宽来自谱仪等能量依赖因素，难以从纯理论定量再现。
- 平移+标度后计算与实验相似；相对峰位与强度要用收敛参数才能显著改善（收敛计算见 Unzog et al.）。

### 功能与标签汇总
- 核空穴方法：`ICORELEVEL`（1=初态近似无激子、2=终态近似）、`CLNT`/`CLN`/`CLL`/`CLZ`、XCH（CLZ=1 自动）vs FCH（`NELECT` 减一+背景电荷）
- 谱控制：`CH_NEDOS`、`CH_SIGMA`（展宽）、`CH_AMPLIFICATION`（强度标度）、NBANDS（足够空带）
- 超胞收敛：核空穴与周期镜像相互作用最小化优先于 k 点收敛
- 分析：PARCHG 部分电荷密度看激子局域化；谱的平移/标度规范

---

## 12.2 XAS（二）：BSE 方法计算 LiCl K 边 XAS

> 来源：<https://vasp.at/tutorials/latest/xas/part2/>
> 配套输入文件：[XAS-part2.zip](https://vasp.at/tutorials/latest/XAS-part2.zip)

本页用多体微扰理论（MBPT）模拟岩盐 LiCl 的 X 射线吸收：例 3 先做 $G_0W_0$ 得准粒子与 RPA 介电张量，例 4 用 BSE 计算含激子效应的 XAS、激子波函数并区分暗/亮激子。

（GW 背景回顾：准粒子方程、$\Sigma=GW$、RPA 屏蔽、单发 $G_0W_0$——与 GW 教程相同。）

### 例 3：LiCl 的 $G_0W_0$ 带隙

**教学目标**：运行单发 $G_0W_0$、得带隙、算 RPA 介电函数。

- 三步流程（目录 `e03_LiCl_G0W0`：dft → unoccupied-states → gw）：① PBE 基态（立方 LiCl；OUTCAR 的 `fundamental gap:` 给出 DFT 基本带隙 **6.4 eV**）；② `ALGO=Exact`+`NELM=1` 精确对角化算约 3–4 倍默认数量的未占据 KS 轨道（NBANDS），`LOPTICS` 写 WAVEDER；副产品为 IPA 介电函数（不含激子效应）；③ `ALGO=EVGW0`+`NELMGW=1`，NBANDS 与第 2 步一致（示例文件内容略）。
- $G_0W_0$ 结果：OUTCAR 的 QP-energies 给出基本准粒子带隙 **8.9 eV**。
- IPA vs RPA 介电函数对比（$\omega=0$ 实部）：实验介电常数 2.7，看哪个近似更接近。
- 强调：GW 需对空态数 `NBANDS` 与频率点数 `NOMEGA` 做收敛分析（本教程收敛度较低）。

### 例 4：BSE 计算 XAS

**教学目标**：模拟有/无激子效应的 XAS；计算并可视化激子电荷密度；区分暗激子与亮激子。

**BSE 形式（针对芯态）**：BSE 可写成双粒子极化率 $L$ 的 Dyson 方程 $L=L_0+L_0(v-W)L$；TDA 下 $AX_\lambda=\omega_\lambda X_\lambda$，双粒子哈密顿量

$$A_{vc}^{v'c'} = (\varepsilon_v-\varepsilon_c)\delta_{vv'}\delta_{cc'} + \langle cv'|V|vc'\rangle - \langle cv'|W|c'v\rangle$$

光学区 $v,v'$ 只含价带；**XAS 中 $v,v'$ 必须包含激发的芯态**。能量/轨道取自 GW 准粒子；$W$ 取静态近似。

**所需文件**：与 $G_0W_0$ 步一致的 POTCAR/KPOINTS/POSCAR，加上 $G_0W_0$ 产生的 WAVECAR（准粒子能量与轨道）、WFULLxxxx.tmp（RPA 介电张量）；WAVEDER 只在包含价轨道时需要。

**输入要点**（示例文件内容略）
- `ALGO=BSE`；与光学区 BSE 不同，芯态 BSE 只需指定导带数，价带数设为零：`NBANDSO=0`；芯态用与 $\Delta$SCF 相同的标签指定：`ICORELEVEL`、`CLNT`、`CLN`、`CLL`。
- 无激子对照（RPA）：关梯子图 `LADDER=.FALSE.`、保留裸库仑 `LHARTREE=.TRUE.`。
- OUTCAR 检查：BSE 哈密顿量秩 `rank= 3456` 应等于 $N_c\times N_v\times N_k$（芯态数x导带数x全布里渊区 k 点数）。

**结果与分析**
- vasprun.xml 的振子强度（偶极近似+选择定则计算）：第一个跃迁是**暗激子**（振子强度为零，违反选择定则、无光学活性）；更高能量的激子是**亮激子**；跃迁 2/3/4 简并且同能量。
- 排错：虚部介电函数全为零时检查 stdout 的 `the WAVECAR file was read successfully`。
- BSE 谱与实验（Shirley 2004）对比：主要特征与相对位置吻合良好；**绝对位置再现差**——PAW 中芯态冻结、直接取自 POTCAR，GW QP 移动只作用于价带与导带，芯态固定在 DFT 水平。
- 收敛提醒：6x6x6 网格未收敛（收敛结果见 Unzog et al. 2022、Olovsson et al. 2009）；XAS 激子强局域，介电函数需对导带数充分收敛。
- BSE vs RPA 介电函数对比：唯一差别是激子效应——清楚展示芯空穴-电子相互作用对 XAS 的重要性。第一个亮激子结合能即 RPA 与 BSE 本征值之差：$59.08-57.37=1.71$ eV。
- BSE 算法内的 RPA 实部与 $G_0W_0$ 步的 RPA 差别巨大：BSE 中只考虑芯电子，1s 芯态对屏蔽贡献极小（介电常数接近 1）；概念上还有 TDA 之别（BSE 内 RPA 用 TDA，GW 步未用，本例数值影响不大）。
- $\Delta$SCF 与 BSE 对比见 Unzog et al. (2022)：两者吻合，BSE 总体更接近实验。

**激子波函数可视化**
- 激子波函数是双坐标函数 $\psi_\lambda(\mathbf r_e,\mathbf r_h)=\sum_{vc}X^\lambda_{vc}\phi^*_v(\mathbf r_h)\phi_c(\mathbf r_e)$；3D 可视化需固定一个坐标——XAS 中芯空穴强局域，惯例是**固定空穴坐标**。
- INCAR：`NBSEEIG` 选要分析的本征矢数目；`BSEHOLE` 把空穴位置固定在激发原子坐标（示例文件内容略）。
- 结果：第一个亮激子三重简并，画波函数需把三个分量相加；其电荷密度很好地局域在两个最近邻 Li 原子之间、以 Li 原子为中心。py4vasp：`bse.exciton.density.to_ngl("1", isolevel=1, center=True)` 可选第一个暗激子。

### 功能与标签汇总
- 芯态 BSE：`ALGO=BSE`+`NBANDSO=0`+`ICORELEVEL`/`CLNT`/`CLN`/`CLL`；需要 WAVECAR + WFULLxxxx.tmp（+WAVEDER 若含价带）
- RPA 对照：`LADDER=.FALSE.`+`LHARTREE=.TRUE.`
- 激子分析：`NBSEEIG`、`BSEHOLE`（固定空穴坐标）、激子密度可视化（py4vasp `to_ngx`/`to_ngl`）、简并亮激子分量相加
- 概念：暗/亮激子与偶极选择定则、激子结合能=RPA-BSE 本征值差、芯态冻结导致绝对能量偏差、芯态对屏蔽贡献可忽略

---

## 13. 强关联体系（Strongly Correlated Systems）

强关联材料中电子–电子相互作用占主导，标准 DFT 的独立粒子近似不足以描述（通常含部分填充的局域轨道），可出现金属–绝缘体转变、磁性与非常规超导。教程覆盖 DFT+U、约束随机相位近似（cRPA）与动力学平均场理论（DMFT）等扩展。

## 13.1 强关联（一）：约束随机相位近似（cRPA）求 Hubbard 参数（NiO）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part1/>
> 配套输入文件：[strongly_correlated-part1.zip](https://vasp.at/tutorials/latest/strongly_correlated-part1.zip)

**背景**：强关联材料中电子-电子相互作用占主导，标准 DFT 的独立粒子近似无法充分描述（典型为局域的部分占据 $d$/$f$ 轨道元素），关联效应导致金属-绝缘体转变、磁性、非常规超导等。VASP 的扩展方法包括 cRPA（本部分）、DFT+U 与 DFT+DMFT（第二部分）。cRPA 用于计算模型哈密顿量的有效相互作用参数 $U$、$J$、$U'$（Phys. Rev. B 77, 085122 (2008)；Phys. Rev. B 112, 245102 (2025)）：核心思想是在 GW 的屏蔽库仑作用 $W$ 中**剔除目标态自身的屏蔽贡献**，得到的部分屏蔽库仑作用在张成目标空间的局域基（模型哈密顿量）下评估；目标空间通常低维（$\leq 5$ 态），便于用 DMFT 等高级理论处理。

cRPA 工作流五步：① DFT 基态（例 1）；② 未占据 KS 态（例 1）；③ DOS 投影与 fatband（例 1）；④ Wannier 投影（例 2）；⑤ cRPA 计算（例 3）。体系为原胞 NiO，12x12x12 k 网格。

### 例 1：能带结构、未占据态与 DOS 投影

**教学目标**：计算基态波函数；生成 cRPA 所需未占据态；投影 DOS（含 fatband）。

**步骤与要点**（目录 `e01_bands`，示例文件内容略）
- DFT 基态：GGA-PBE，标准标签（`SYSTEM`、`NBANDS`、`ISMEAR` 等）；POTCAR 用 `_GW` 版本（需大量未占据态）。结果：NiO 的 PBE 带隙很小，看似金属，但因强关联电子实为 **Mott 绝缘体**——直观展示标准 DFT 的失败。
- 未占据态（`unoccupied/`）：增大 `NBANDS`、`ALGO=Exact` 单次精确对角化；`LOPTICS=T` 计算频率依赖介电矩阵；`LFINITE_TEMPERATURE=T` 是计算全部光学矩阵元所需；WAVECAR/WAVEDER 作为 cRPA 起点。
- fatband（`fatbands/`）：`EMIN`/`EMAX` 定义 pDOS 能量范围、`NEDOS` 定网格点数。
- 分析：费米能附近约 -8 至 2 eV 有 8 条能带承载大部分 Ni $d$ 与 O $p$ 特征（带指标 2–9；PROCAR 含 spd 与格点投影的波函数特征）。
- **目标空间选择**：只取 Ni 的 5 个 $d$ 态会使基函数整体不够局域；纳入 O $p$ 态（离域到 O 上）允许混合、改善费米能附近 Wannier 轨道的描述。$s$ 态与更高 $p$ 态纠缠（~9 eV 处交叉）难以干净分离，且对费米能附近 5 条带贡献可忽略，故排除。py4vasp 可画 `dos.plot("d"/"p"/"Ni"/"O")` 验证：$d$ 主导态有 $p$ 贡献、~-5 eV 的 $p$ 主导态有 $d$ 贡献——说明必须含 $p$；~7 eV 处几乎无 $d$、0 eV 附近几乎无 $s$——说明可排除 $s$。

### 例 2：Wannier 投影

**教学目标**：用 Wannier90 生成 Wannier 函数；用 gnuplot 画其投影。

**理论**：Bloch 轨道 $\psi^\sigma_{n\mathbf k}(\mathbf r)=e^{i\mathbf{k\cdot r}}u^\sigma_{n\mathbf k}(\mathbf r)$；Wannier 函数由 Bloch 轨道线性组合得到

$$|w^\sigma_{i\mathbf R}\rangle = \frac{1}{N_k}\sum_{n\mathbf k}e^{-i\mathbf{k\cdot R}}T^{\sigma(\mathbf k)}_{in}|\psi^\sigma_{n\mathbf k}\rangle$$

其中 $T$ 为变换矩阵（此处用 $T$ 记号以避免与 DFT+$U$ 的 $U$ 混淆）。基组足够局域时可只在 $\mathbf R=0$ 原胞工作。模型哈密顿量（如 cRPA）只用一小部分 Bloch 函数——**目标态**通常围绕化学势/费米能且强局域于离子附近；只有目标态被 Wannier 基很好表示时模型才可解。

**输入要点**（目录 `e02_wannier_projection`，示例文件内容略）
- `LWRITE_WANPROJ`：是否写 WANPROJ 文件（含每个参与投影的 k 点与能带的 Wannier 投影矩阵 $T$）；
- `NUM_WANN`：Wannier 函数个数；`WANNIER90_WIN`：wannier90.win 内容（用 `exclude_bands` 只纳入 8 条带；另设 Wannier 函数数、投影带数与轨道类型）；`LWANNIER90_RUN`：让 VASP 运行 wannier90 库；
- `ALGO=NONE`：后处理步，不改动波函数；
- 用 gnuplot 文件画 Wannier fatband。

**结果**
- Wannier 能带着色显示轨道特征（Ni $d$ 红、O $p$ 蓝）：费米能附近 5 条 Wannier 带以 $d$ 特征为主、有重要 $p$ 贡献；~-7 eV 的 3 条低能带以 $p$ 为主、有显著 $d$ 贡献。
- 费米能附近的 5 个 $d$ 态即 cRPA 的**目标态**。wannier90.wout 显示前 5 个 Wannier 函数中心在 Ni 原子 (0,0,0)，后 3 个在 O 原子 (2.0855, 2.0855, 2.0855)。应选**展宽（spread）最小**的轨道：Ni 中心 $d$ 轨道 ~0.45，优于 O 中心 ~0.87（spread 定义见 Phys. Rev. B 56, 12847 (1997) 式 14）。
- （可选）仅用 $d$ 轨道：Wannier 态数 8→5，`exclude_bands=1-4` 排除 3 条 $p$ 带并去掉 `O:p` 投影；结果 $t_{2g}$ 与 $e_g$ 的 spread 分别为 0.95788 与 1.6914，比含 $p$ 时大 2–4 倍——说明生成 Wannier 基时纳入 $p$ 轨道非常重要。

### 例 3：cRPA 计算

**教学目标**：用粗/细 k 网格计算 cRPA Hubbard $U$、$U'$、$J$；计算离中心（off-center）项；把 Ohno 势拟合到 Hubbard $U$ 势并研究等离激元的频率依赖。

**理论**：在选定的 Wannier 轨道下计算 cRPA 相互作用矩阵（存于 WANPROJ/vaspout.h5）：

$$U^{\sigma\sigma'}_{ijkl} = \int d\mathbf r\int d\mathbf r'\,w_i^{*\sigma}(\mathbf r)w_j^{\sigma}(\mathbf r)U(\mathbf r,\mathbf r',\omega)w_k^{*\sigma'}(\mathbf r')w_l^{\sigma'}(\mathbf r')$$

常用 Hubbard–Kanamori 参数作为指标（单胞 $R=0$；对应 Vijkl/Uijkl 标签）：

$$\mathcal U^{\sigma\sigma'} = \frac{1}{N}\sum_{i\in\mathcal T}U^{\sigma\sigma'}_{iiii},\quad \mathcal U'^{\sigma\sigma'} = \frac{1}{N(N-1)}\sum_{i\neq j}U^{\sigma\sigma'}_{ijji},\quad \mathcal J^{\sigma\sigma'} = \frac{1}{N(N-1)}\sum_{i\neq j}U^{\sigma\sigma'}_{ijij}$$

扩展体系（$R>0$）对应 VRijkl/URijkl。实际 DFT+U 计算（Phys. Rev. B 44, 943 (1991)；Phys. Rev. B 57, 1505 (1998)）推荐球平均，在 DFT+DMFT 教程讨论。

**输入要点**（目录 `e03_CRPA`，示例文件内容略）
- `ALGO=CRPA`；`LOCALIZED_BASIS` 设置局域基——工作目录有 WANPROJ 时 VASP 直接用它（复用例 2 的 Wannier 投影）；
- `NTARGET_STATES` 控制在 cRPA 中**被排除（不屏蔽）**的 Wannier 态；
- `LFINITE_TEMPERATURE` 开启有限温度形式（要求 Fermi 展宽 `ISMEAR=-1`），零带隙或纠缠体系必需；
- `NBANDS` 与生成 WAVECAR 时一致；k 网格 4x4x4（粗）与 6x6x6（细，`fine_grid/`，必须与 DFT 步的 12x12x12 公度，8x8x8 不行）。

**结果与分析**
- OUTCAR 中找 `Hubbard U` 行；完整数据在 vaspout.h5。
- 粗→细 k 网格：裸相互作用几乎不变（已收敛）；$J$ 已收敛；$U$ 与 $U'$ 未收敛到 0.1 eV，还需更密 k 网格。DMFT 需要的是屏蔽后的 $U$、$U'$、$J$。
- **离中心项**：`LTWO_CENTER` 计算离中心库仑积分，vaspout.h5 的 `uijkl` 数据集多出 Wannier 中心 $R$ 维度：
$$U_{ijkl}(R) = \int d\mathbf r\int d\mathbf r'\,w_i(\mathbf r-\mathbf R)w_j(\mathbf r-\mathbf R)U(\mathbf r,\mathbf r',\omega=0)w_k(\mathbf r')w_l(\mathbf r')$$
  Hubbard–Kanamori $U$ 的空间衰减常用 Ohno 势描述：$U(R)/U(0)=\sqrt{\delta/(R+\delta)}$；py4vasp 可拟合（调节 `radius_max` 限制拟合点）。噪声主要来自 k 采样不足；教程演示剔除离群点（如 >10 Å 处高于 0.74 eV 的点）后重拟合，但**实际工作中应加密 k 网格与提高平面波截断能**（噪声源于过疏的 FFT 网格）。
- **等离激元频率依赖**（`plasmons/`）：cRPA 可用欧几里得时空算法评估大原胞（内部在虚时间轴计算，Minkowski 度规 $(-+++)$ 变为欧氏 $(++++)$，允许库仑核的稀疏表示）；小体系更耗内存/时间，大体系快得多（4 MPI 进程 + 4x4x4 网格约需 8 GB RAM；打印 `estimated memory requirement` 后有较长静默期属正常）。
- 用 AAA 拟合做解析延拓；改变 k 网格与频率点数检验峰位：第一个等离激元峰（5.5 eV，约对应 $d$–$p$ 中心间距）位置基本不变——是真实等离激元峰而非 Froissart 双峰；该峰之前有效库仑作用几乎恒定，说明对 NiO 的 Hubbard 类模型**静态近似是合适的**。高频（~50 eV 主等离激元峰，对应 GW 类近似的集体激发）难收敛、可靠性低；更高频处屏蔽势 $U$ 趋近裸势 $V$。
- `NOMEGA` 12→24：Minimax 等距方法已很好选择虚频点（Phys. Rev. B 101, 205145 (2020)），增加频率点对结果改善有限；反而加密 k 点能稳定高频结构、允许更小 `clean_up_tol`。关键结论：**解析延拓不是黑箱**——低频（模型哈密顿量物理最有趣处）峰可靠，应聚焦收敛；高频峰（~50 eV，接近 伽马射线能量）意义不大。
- （可选）仅 $d$ 轨道的 Wannier 化：裸 Hubbard U 从 26.3724 降到 22.8494（4x4x4），屏蔽 U 从 7.0806 降到 6.1218——大量相互作用被漏掉，再次说明 $p$ 轨道的必要性。

### 功能与标签汇总
- 工作流：DFT（`_GW` POTCAR）→ `ALGO=Exact`+`LOPTICS`+`LFINITE_TEMPERATURE`（未占据态/WAVEDER）→ fatband/PROCAR 选目标带 → `LWRITE_WANPROJ`+`NUM_WANN`+`WANNIER90_WIN`+`LWANNIER90_RUN`+`ALGO=NONE`（WANPROJ）→ `ALGO=CRPA`+`LOCALIZED_BASIS`+`NTARGET_STATES`
- 收敛：k 网格须与 DFT 公度；裸作用先收敛、$U$/$U'$ 后收敛；`NOMEGA` 与 AAA 解析延拓、Froissart 双峰识别
- 扩展输出：`LTWO_CENTER`（off-center $U(R)$）、Ohno 势拟合、vaspout.h5 的 Hubbard 数据、欧几里得时空算法（大胞）
- 物理：Mott 绝缘体、目标态/Wannier 基选择（spread 判据、$p$–$d$ 混合）、静态近似适用性

---

## 13.2 强关联（二）：DFT+U 与 DFT+DMFT（NiO 全流程）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part2/>
> 配套输入文件：[strongly_correlated-part2.zip](https://vasp.at/tutorials/latest/strongly_correlated-part2.zip)

本页先跑 VASP（DFT 与 DFT+U），再用 `solid_dmft`（TRIQS）以 VASP 的投影子（PLO）做 DFT+DMFT。五个例子：例 4 顺磁/反磁 DFT 与 DMFT 准备、例 5 反磁 DFT+U 与 DOS、例 6 单发 DFT+DMFT（Hartree 求解器）、例 7 电荷自洽（CSC）DFT+DMFT、例 8 反磁 vs 顺磁、Hartree vs QMC 对比。

### 例 4：DFT SCF 与 DMFT 准备

**教学目标**：设置 AFM 原胞 NiO 的 DFT 计算；控制 DMFT 所需 Ni-$d$ 投影子生成；了解 DMFT 相关输出文件；从 cRPA 提取 $U$、$J$。

**DFT+U 理论回顾**：半局域泛函的虚假自相互作用使电子人为离域、低估 Mott 绝缘体带隙；DFT+U 加 Hubbard 式在位库仑罚项。Dudarev 简化旋转不变形式（Phys. Rev. B 57, 1505 (1998)）：

$$E = E_{DFT} + \sum_I\left[\frac{U^I}{2}\sum_{m\sigma\neq m'\sigma'}n_m^{I\sigma}n_{m'}^{I\sigma'} - \frac{U^I}{2}n^I(n^I-1)\right]$$

第二项为**双计数修正**（减去 $E_{DFT}$ 已含的相互作用估计）。括号项对整数占据消失（原子极限不改变 $E_{DFT}$），但惩罚分数占据、驱动体系走向局域的整数解。DFT+U 便宜且常定性正确，但 $U$ 需谨慎选取（cRPA 或线性响应），且不能捕捉动态关联与 Kondo 物理——这正是 DFT+DMFT 的动机。

**输入要点**（目录 `e04_scf`，示例文件内容略）
- 收紧 `EDIFF`、增大 `NELMIN`/`NBANDS` 保证严格收敛（OSZICAR 中残差范数 `rms`）；并行：`NCORE`（每带 FFT 的核数）、`KPAR`（并行 k 点数）；`ISMEAR` 用四面体+Blöchl；`KSPACING` 按间距生成 k 网格；`LORBIT` 选原子球投影方法；`EMIN`/`EMAX`/`NEDOS` 控制 DOS；`LMAXMIX`（单中心 PAW 电荷密度通过混合器并写入 CHGCAR 的最大 $l$）；
- **`LOCPROJ`**：指定投影轨道的局域函数（Ni-$d$）；`LORBIT=14` 在 EMIN–EMAX 能量窗内自动优化 PAW 投影子通道选择；
- POTCAR 用 `Ni_pv` 显式处理与 $d$ 弱杂化的 $p$ 价态（强烈影响能带与晶格参数；理想是 `Ni_sv_GW`，为省时间未用）。
- 输出：vaspout.h5 与投影子输出（LOCPROJ/PROJCAR）。

**PLO 构造与转换**
- DFT+DMFT 需要关联子空间（此处 Ni $3d$）；本教程用**投影局域轨道（PLO）**。VASP 按 LDAU 相同方式把 KS 态投影到 PAW 局域轨道：
$$P^{\mathbf R}_{L,\nu}(\mathbf k) = \sum_i\langle\chi^{\mathbf R}_L|\phi_i\rangle\langle\tilde p_i|\tilde\Psi_{\nu\mathbf k}\rangle$$
  $|\chi^{\mathbf R}_L\rangle$ 为关联格点的局域基函数，$|\phi_i\rangle$ 全电子部分波，$|\tilde p_i\rangle$ PAW 投影算符；原始投影子正交归一化后连接 KS 基 $\nu$ 与局域轨道 $m$。
- `plo.cfg` 关键项：`EWINDOW`（找相关 KS 态的能量窗）、`BANDS`（投影的 VASP 带指标窗）、`LSHELL=2`（$d$ 壳层）、`IONS`（携带关联壳层的原子，AFM 为两个 Ni）、`CORR=True`（标记为 DMFT 处理的关联壳层）。
- TRIQS/dft_tools 的 `converter.py` 读 VASP 输出生成 `vasp.h5`（solid_dmft 输入）与投影子 DOS（`pdos_*.dat`）；log.converter 显示生成的杂质格点数与 5x5 占据矩阵 $n_m^{I\sigma}$。
- 合理性检查：VASP 投影 DOS 与 PLO 投影 DOS 吻合良好（两者正交归一化程序略有不同，~-5 eV 处高度有差异）。

![VASP 局域投影（蓝）与 PLO（青）投影 DOS 对比](/imgs/2026-08-05/fig34.png)

**从 cRPA 提取 U/J（UIJKL_cRPA）**
- 用 TRIQS 的 `fit_slater_fulld` 把 cRPA 的完整四指标库仑张量 $U_{ijkl}$ 拟合到 Slater 参数化相互作用：① 用 `U_matrix_slater` 由试验 Slater 参数（$U$,$J$ 或 $F^0,F^2,F^4$）构造旋转不变张量；② `reduce_4index_to_2index` 化为密度-密度 $U_{iijj}$ 与交换 $U_{ijij}$ 矩阵；③ `scipy.optimize.minimize` 最小二乘拟合：
$$\min_{U,J}\sum_{i,j}\left[(U^{cRPA}_{iijj}-U^{Slater}_{iijj})^2 + (U^{cRPA}_{ijij}-U^{Slater}_{ijij})^2\right]$$
- 得到的 ($U$,$J$) 一致用于后续 DFT+U（e05）与各 DMFT 配置（e06/e07/e08）。

### 例 5：AFM DFT+U 与 DOS

**任务**：用 cRPA 参数（$U=6.181$、$J=1.134$）计算 AFM NiO 的 DOS。

**输入要点**（目录 `e05_dftu`，示例文件内容略）
- `EFERMI=MIDGAP`（费米能取带隙中点）；`LASPH`（PAW 球内梯度修正的非球贡献，磁性计算应视为标准设置）；`ISPIN=2`（自旋极化）；`MAGMOM`（各原子初始磁矩）；
- `LDAU` 开启 DFT+U；`LDAUTYPE=1`（Liechtenstein 方案）；`LDAUL`（加在位相互作用的 $l$ 量子数）、`LDAUU`（在位库仑 $U$）、`LDAUJ`（在位交换 $J$）、`LDAUPRINT`（输出详略）。

**结果与物理**
- NiO 是原型 Mott 绝缘体，Ni$^{2+}$（$d^8$）：八面体晶场把 Ni $d$ 分裂为低能 $t_{2g}$ 三重态与高能 $e_g$ 双重态；8 个 $d$ 电子填满 $t_{2g}$（6 个），其余 2 个占据 $e_g$。
- 加 $U$ 并允许磁序后：交换劈裂驱动高自旋 $S=1$——多数自旋 $e_g$ 全占、少数自旋 $e_g$ 全空；DFT+U 把占据的多数自旋 $e_g$ 推向低能、空的少数自旋 $e_g$ 推向高能，打开带隙（~3 eV），DOS 在 ~2.5 eV 出现强 $e_g$ 峰；全满的 $t_{2g}$ 基本不受影响。
- 带隙具**电荷转移**特征（O $2p$ → Ni $e_g$）而非纯 Mott–Hubbard，但 DFT+U 恢复了接近实验 ~4 eV（Phys. Rev. Lett. 53, 2339 (1984)）的带隙。

![DFT（灰虚线）与 DFT+U（紫）DOS 对比，绿线标 VBM/CBM](/imgs/2026-08-05/fig35.png)

### 例 6：单发 DFT+DMFT（Hartree 求解器）

**教学目标**：用 solid_dmft 的静态 Hartree 杂质求解器做单发 DFT+DMFT；监控收敛诊断；后处理得关联谱函数并与 DFT+U 对比。

**理论**：DFT 收敛能量与密度，DMFT 收敛 Green 函数与自能。单粒子 Green 函数 $G(\mathbf k,\omega)$ 编码向相互作用体系加/减电子的概率幅，是 DMFT 的自然语言：谱函数 $A(\mathbf k,\omega)=-\frac{1}{\pi}\mathrm{Im}\,G(\mathbf k,\omega+i0^+)$ 直接对应 ARPES 强度；k 积分谱函数在无相互作用极限还原为 DOS。

**计算与监控**（目录 `e06_dmft_os`）
- 复制（或软链接）`e04_scf/vasp.h5` 后运行 solid_dmft（4 核约 6–10 分钟）。
- `dmft_config.toml` 分 `general`/`solver`/`dft`/`advanced` 节（toml 规范）：关键参数 `U`、求解器 `type`（此处 Hartree）、`n_iter_dmft`（DMFT 迭代数）、`jobname`、`dc_type`（双计数类型）；核对 $U$、$J$ 与 cRPA 提取值一致。
- 收敛诊断（`observables_imp0_up.dat`、`conv_imp0.dat`、h5 归档）：化学势 `mu` 收敛到定值；每自旋通道杂质占据；Weiss 场收敛 $\delta\mathcal G^0(i\omega_n)$；DMFT 自洽条件 $\delta G_{imp}=||G^{loc}-G^{imp}||$。杂质密度矩阵在求解器前后应一致。

![DMFT 收敛指标](/imgs/2026-08-05/fig36.png)

**后处理：谱函数**
- 晶格 Green 函数：$\hat G(k,\omega)=[\omega+\mu-\hat\epsilon(k)-(\hat\Sigma(\omega)^{imp}-\hat\Sigma^{dc})]^{-1}$（$\hat\epsilon(k)$ 为 KS 本征值，$\hat\Sigma^{dc}$ 双计数修正）；谱函数 $A(k,\omega)=-\frac{1}{\pi}\mathrm{Im}\,G(k,\omega)$；k 积分谱函数 $A(\omega)=-\frac{1}{\pi}\sum_k\mathrm{Im}\,G(k,\omega)$ 对应 PES。
- `calc_Aw_hartree.py`（MPI 并行）从 h5 归档读 $\Sigma^{imp}$、$\mu$、$\Sigma^{dc}$（`dc_energ`）与 `block_structure`，构造 SumK 对象并计算 k 积分谱函数；注意 Sigma 网格须覆盖能带能量范围（脚本会打印相关警告）。
- 结果：DFT+U 与 DFT+DMFT@Hartree 谱几乎相同——验证二者方法等价（静态极限下 DMFT@Hartree 还原为 DFT+U），带隙均 ~3 eV。

![DFT+U 与 DFT+DMFT@Hartree 谱函数对比](/imgs/2026-08-05/fig37.png)

### 例 7：电荷自洽（CSC）DFT+DMFT

**教学目标**：解释 CSC 条件及其重要性；以单发结果为起点配置 CSC 计算；评估 CSC 对谱函数与带隙的影响。

**CSC 条件**：除 DMFT 自洽 $|G^{imp}-G^{loc}|\to 0$ 外，DFT+DMFT 还需 DFT 密度与 DMFT 密度矩阵一致：由晶格 Green 函数得密度矩阵 $N_{\nu\nu'}(\mathbf k)=\frac{1}{\beta}\sum_{i\omega_n}G_{\nu\nu'}(\mathbf k,i\omega_n)$，与 DFT 密度 $\rho(\mathbf r)=\sum_{\mathbf k}\sum_{\nu\nu'}\langle\mathbf r|\Psi_{\nu\mathbf k}\rangle N_{\nu\nu'}(\mathbf k)\langle\Psi_{\nu'\mathbf k}|\mathbf r\rangle$ 须吻合——因此每个 DMFT 步后要更新电荷密度（`ICHARG=5`），DMFT 自能反馈进 DFT 循环。

**输入要点**（目录 `e07_dmft_csc`，示例文件内容略）
- INCAR：`NELMIN`/`NELM` 设得很大（防止 VASP 提前终止——由 TRIQS 控制终止）；`IMIX`（Broyden 2 + Pulay 混合）、`AMIX`（线性混合参数）、`BMIX`（混合截断波矢）；`LSYNCH5`（运行中实时同步 vaspout.h5，可用 py4vasp/solid_dmft 监控）；`LWAVE`/`LCHARG`。
- dmft_config.toml 新增：`load_sigma=true`+`path_to_sigma`（从 e06 加载自能做起点）；`n_iter_dmft_first=2`、`n_iter_dmft_per=2`（每次 CSC 更新跑 2 步 DMFT）、`n_iter_dmft=10`；另有完整 `dft` 节。
- 从收敛的 WAVECAR/CHGCAR 软链接启动；out 文件可见 VASP 与 solid_dmft 交替运行；`solid_dmft: Stopping VASP` 表示结束（约 10–15 分钟）。

**结果**：CSC 相对单发 DMFT 略微减小带隙（单发通常高估关联效应，CSC 使之降低）；带隙仍约 3 eV，略小于实验 4.3 eV。

![CSC vs 单发 DMFT 谱函数](/imgs/2026-08-05/fig38.png)

### 例 8：AFM vs 顺磁、Hartree vs QMC

**教学目标**：解读 Matsubara 与实频轴的 QMC 自能、识别强关联特征；区分 AFM 对称破缺带隙与顺磁动态 Mott 带隙；跨 DFT+U / 单发 Hartree DMFT / 顺磁 QMC 批判性比较谱函数。

**物理对比**
- Hartree 计算通过 AFM 长程序破缺自旋对称（两个 Ni 格点静态磁矩相反）；顺磁（PM）解平均上强制自旋对称（$\langle n_\uparrow\rangle=\langle n_\downarrow\rangle$），但保留静态平均场看不见的**动态局域磁矩涨落**——CT-HYB 通过 $\Sigma(\omega)$ 的完整频率依赖捕捉它：即使无长程序，强在位库仑 $U$ 压制双占据，费米能谱权重转移到上下 Hubbard 带，打开 **Mott 带隙**。这是 Néel 温度（$T_N\approx 520$ K）之上的关键物理——NiO 在 PM 相仍绝缘，不是靠磁序，而是关联驱动的局域化本身足以打开谱隙。

**QMC 求解器**（目录 `e08_dmft_os_qmc`，预计算结果在 `UcRPA-beta20-qmc1e+7/vasp.h5`；自己跑需数小时）
- CT-HYB（连续时间杂化展开 QMC，Rev. Mod. Phys. 83, 349 (2011)）：在虚时间 $\tau$ 采样杂质 Green 函数 $G^{imp}(\tau)$，按杂化函数展开杂质哈密顿量；`beta=20` 对应 $T\approx 580$ K；温度越高 $[0,\beta]$ 区间越短、越便宜；求解器 MPI 并行、C++ 实现。
- `[solver]` 精度参数：`n_warmup_cycles`（测量前热化，每个 MC walker/MPI 进程执行，高温简单问题 5000–10000 足够）；`n_cycles_tot`（测量总数，在所有 walker 间分配）；`length_cycle`（两次测量间的 MC 移动数；总移动数 = n_cycles_tot x length_cycle）。
- 若自行准备输入：先在 `scf_plo` 跑 VASP，再跑 converter，把 `vasp.h5` 复制到 e08 目录。

**自能分析**
- Matsubara 轴看 $-\mathrm{Im}\,\Sigma(i\omega_n)$：低频行为——费米液体中 $-\mathrm{Im}\,\Sigma(i\omega_n)\propto\omega_n\to 0$；大/发散值预示强散射或临近 Mott 行为。轨道分化：立方晶场分裂 $t_{2g}$（三重简并，指标 `up_0`）与 $e_g$（二重简并，`up_2`）；$e_g$ 与氧杂化弱、关联应更强。
- 实频自能：用最大熵方法（MaxEnt，`sigma_maxent.py` 后处理）把 $\Sigma(i\omega_n)$ 解析延拓到 $\Sigma(\omega)$，存于 h5 的 `Sigma_maxent_*`；可识别准粒子权重 $Z=(1-\partial\mathrm{Re}\,\Sigma/\partial\omega|_{\omega=0})^{-1}$ 与散射率 $\Gamma\propto-\mathrm{Im}\,\Sigma(\omega=0)$。注意解析延拓是病态反问题，高频细节需谨慎解读。

![t2g 与 eg 的实频自能 Re Σ 与 -Im Σ](/imgs/2026-08-05/fig39.png)

![实频自能用于实频网格的晶格 Green 函数](/imgs/2026-08-05/fig40.png)

**顺磁谱函数与总结**
- 观察要点：**Hubbard 带**（CT-QMC 的动态关联在准粒子流形上下产生非相干谱权重，静态 Hartree 没有）；**带隙性质**（Hartree/AFM 靠磁对称破缺开隙；PM QMC 无长程序而由强局域关联动态生成 Mott 隙）；**谱权重重分布**（QMC 谱更宽、更高束缚能处有额外非相干权重，与光发射谱一致）。

![DFT+U、DFT+DMFT@Hartree 与 PM CT-HYB 三层理论谱函数对比](/imgs/2026-08-05/fig41.png)

- 三者都打开带隙但物理原因不同：DFT+U 与 DMFT@Hartree 谱几乎相同（静态平均场、靠 AFM 对称破缺，~3 eV 锐隙+窄准粒子峰）——验证 DMFT@Hartree 实现在静态极限正确还原 DFT+U；PM CT-HYB 定性不同——无磁长程序仍开隙（顺磁 Mott 绝缘体的标志，计算温度 580 K > $T_N$ 520 K），谱特征更宽（$\mathrm{Im}\,\Sigma$ 给出有限准粒子寿命），更高束缚能（~6–8 eV）出现下 Hubbard 带——动态关联的直接指纹。

### 功能与标签汇总
- VASP 端：`LOCPROJ`（局域投影→LOCPROJ/PROJCAR）、`LORBIT=14`（能量窗内优化投影通道）、`LMAXMIX`、`ISPIN=2`/`MAGMOM`、`EFERMI=MIDGAP`、`LASPH`、DFT+U 全套（`LDAU`、`LDAUTYPE=1`、`LDAUL`、`LDAUU`、`LDAUJ`、`LDAUPRINT`）、CSC 相关（`ICHARG=5`、`IMIX`/`AMIX`/`BMIX`、`LSYNCH5`、大 `NELMIN`/`NELM`）、`Ni_pv` 赝势选择
- TRIQS/solid_dmft 端：plo.cfg（EWINDOW/BANDS/LSHELL/IONS/CORR）→ converter.py → vasp.h5；dmft_config.toml（general/solver/dft/advanced）；Hartree vs CT-HYB 求解器；CSC 参数（load_sigma、n_iter_dmft_first/per）；QMC 精度（n_warmup_cycles/n_cycles_tot/length_cycle）；MaxEnt 解析延拓
- 概念链：双计数修正、PLO 关联子空间、Slater 参数拟合 cRPA 张量、Weiss 场与 $||G^{loc}-G^{imp}||$ 收敛、电荷转移 vs Mott–Hubbard 带隙、对称破缺隙 vs 动态 Mott 隙、Hubbard 带

---

## 13.3 强关联（三）：PBE+U 起点的 $G_0W_0$ 与 BSE 光学性质（NiO）

> 来源：<https://vasp.at/tutorials/latest/strongly_correlated/part3/>
> 配套输入文件：[strongly_correlated-part3.zip](https://vasp.at/tutorials/latest/strongly_correlated-part3.zip)

本页用多体微扰理论（MBPT）矫正 AFM NiO 的带隙并模拟光学性质：例 9 以 PBE+U 波函数为起点做 $G_0W_0$，例 10 BSE 计算光吸收、EELS 与自旋翻转激发。

### 例 9：$G_0W_0$ 电子结构（PBE+U 起点）

**教学目标**：从 PBE+U 波函数出发运行 $G_0W_0$；计算未占据 KS 轨道与 IPA 光吸收；确定 NiO 的 $G_0W_0$ 带隙；计算 RPA 光吸收。

**为什么不用 PBE 起点**：$GW$ 自能的屏蔽由 RPA 极化率给出

$$\chi(1,2) = -i\int d(3,4)\,G(1,3)G(4,1)\Gamma(3,4;2)$$

（顶点函数 $\Gamma(3,4;2)=\delta(3-2)\delta(4-2)$）。以 PBE 为起点的 RPA 依赖**部分误差抵消**：PBE 低估带隙→高估屏蔽，缺顶点修正→低估屏蔽；对 $s$/$p$ 离域带主导的半导体这一抵消很好（Phys. Rev. Lett. 99, 246403 (2007)）。但对局域 $d$ 态主导的体系（如 AFM NiO，PBE 带隙低估更严重）误差抵消失效——RPA 强烈高估屏蔽（介电常数 ~15，实验 ~5.7，Phys. Rev. Materials 2, 073803 (2018)）；Mott 绝缘体（PBE 金属）更甚。修复途径：自洽解 Hedin 方程并加顶点修正（$QSG\tilde W$，昂贵），或**用改进的电子结构做起点**——本教程用 cRPA 确定 $U$ 的 PBE+U。

**三步计算**（目录 `e09_G0W0`，示例文件内容略）
1. **PBE+U 基态**（`dft/`）：`ISPIN=2`+`MAGMOM` 设 AFM；`LMAXMIX=4`（$d$ 态电荷密度须通过混合器写出，自旋极化+U 计算必需，也加速收敛）；`LDAU`+`LDAUTYPE`（Dudarev）+`LDAUL`/`LDAUU`/`LDAUJ`（$U=6.181$、$J=1.134$，来自第一部分的 cRPA）。k 网格 3x3x3（教学用，完整计算需收敛分析）。POTCAR 用 `_GW` 版本；精确结果可能还需把半芯 $s$/$p$ 态纳入价带（`_sv`/`_pv`）。
   - 结果：AFM NiO 实验基本带隙 4.3 eV；PBE 严重低估（~0.95 eV）；加 Hubbard 修正显著改善到 ~3.6 eV（仍低估）。
2. **未占据态**（`unoccupied/`）：`ALGO=Exact`（小胞算大量空带的优选方式）；`LOPTICS` 算频率依赖介电矩阵与轨道导数（WAVEDER）；增大 `NBANDS`。副产品：IPA 介电函数——$\omega=0$ 实部（介电常数）相对实验 5.7 被高估，这是带隙低估的后果；若带隙再现良好，IP 近似因缺电子-空穴相互作用反而低估介电常数。
3. **$G_0W_0$**（`g0w0/`）：NBANDS 与上一步一致；**注意 $G_0W_0$ 计算中不应再加 Hubbard 修正**；`NOMEGA`（极化率与自能的实频率点数）；`ALGO=EVGW0`+`NELM=1`；`PRECFOCK`/`ENCUTGW` 为教学降低（完整计算需收敛分析）。
   - 结果：OUTCAR `spin component 1` 后 `QP-energies` 列给出 $G_0W_0$ 带隙 **4.6 eV**（合理接近实验 4.3 eV；本计算未对 k 网格/频率数/带数完全收敛）。RPA 介电常数 5.76（略高于实验）；严格收敛后为 5.2 与带隙 4.25 eV——RPA 相对实验略低估主要因缺激子效应。结论：PBE+U 起点同时给出好的电子结构（基本带隙）与后续 BSE 所需的屏蔽作用 $W$。

### 例 10：NiO 的 BSE 光学性质

#### 10.1 光吸收

- INCAR 选 BSE 算法并指定纳入 BSE 核的价带/导带数（示例文件内容略）。
- **求解算法**：`IBSE=2`（精确对角化）最精确，直接得本征值与本征矢（算振子强度所需），但按 $N^3$ 标度（$N$ 为 BSE 矩阵大小）。
- 谱与实验合理吻合（位置与整体形状）；谱在 16 eV 截断由 `OMEGAMAX` 设定；教程附收敛良好的参考结果说明 3x3x3 网格已能正确再现主峰位置与谱形。
- **寿命展宽**：静态近似的 BSE 解为实值、寿命无限；实际激子因与其他激子（声子为主）耦合有有限寿命。方法忽略这些，用本征值的小复移动引入 Lorentz 展宽——由 `CSHIFT` 设定（此处 0.3 eV，通常调到与实验主要特征宽度匹配）。
- vasprun.xml 的 `<varray name="opticaltransitions">`：第一列激子态能量、第二列振子强度。前 12 个跃迁强度为零——违反选择定则的**暗激子**；第一个**亮激子**出现在 ~4 eV（吸收谱中可见）。

#### 10.2 EELS

- EELS 中入射电子携带有限动量，可激发有限动量转移 $q$ 的态（光子只能竖直激发）；$q$ 由 `KPOINT_BSE` 指定（只能算计算中包含的特定 $q$；对应 k 点指标见 OUTCAR `irreducible k-points` 之后）。
- EELS 信号来自 BSE 宏观介电函数：$\mathrm{EELS}(\mathbf q,\omega)=-\mathrm{Im}[\epsilon^{-1}(\mathbf q,\omega)]$。
- 光学区 TDA 很好，但**有限动量计算必须含共振-反共振耦合**（`ANTIRES=2`，Phys. Rev. B 91, 045209 (2015)）。
- 观察等离激元需宽能量范围（`OMEGAMAX=45`），BSE 矩阵大增，精确对角化太贵——改用迭代算法：`IBSE=1`（时间演化）或 `IBSE=3`（Lanczos，通常算介电函数最快，但当前实现不支持 full BSE），故用 `IBSE=1`；注意 IBSE=1 时 vasprun.xml 只写介电函数、没有 opticaltransitions 块。
- 实验常在小动量转移测谱，可能需要稠密 k 采样才能包含对应 $q$ 点；本例选第二个不可约 k 点（0.25, 0, 0），`KPOINT_BSE=2`。
- 结果：EELS 谱与实验合理一致；最显著特征是 20 eV 起损失函数陡升——对应纵向等离子体频率，即价电子的集体电荷振荡。

#### 10.3 自旋翻转激发

- 用 `LHARTREE=.FALSE.` 做 BSE 找自旋翻转解（目录 `e10_BSE/spin_flip`，示例文件内容略）。
- vasprun.xml 中第一个跃迁能量比吸收计算低 50 meV——该移动来自排斥裸库仑贡献的缺失（`LHARTREE=.FALSE.`）。
- 所得跃迁对应**自旋翻转激发**，违反光学偶极跃迁的自旋守恒，因此**无光学活性**；这些激发包含集体自旋激发即**磁振子（magnon）**，可在 BSE 框架内分析。

### 功能与标签汇总
- 起点策略：PBE+U（cRPA 的 $U$/$J$）→ $G_0W_0$（GW 步不加 LDAU）；`LMAXMIX=4` 对磁性+U 计算必需
- GW 控制：`ALGO=EVGW0`+`NELM=1`、`NOMEGA`、`PRECFOCK`、`ENCUTGW`、`_GW`（可加 `_sv`/`_pv`）POTCAR
- BSE 求解：`IBSE=2`（精确对角化，$N^3$）/`IBSE=1`（时间演化，支持 full BSE）/`IBSE=3`（Lanczos，最快但不支持 ANTIRES=2）；`ANTIRES=2`（有限动量必需）；`CSHIFT`（Lorentz 寿命展宽）；`OMEGAMAX`
- EELS：`KPOINT_BSE`（动量转移 $q$）、$-\mathrm{Im}\,\epsilon^{-1}$ 损失函数、20 eV 等离激元
- 自旋翻转：`LHARTREE=.FALSE.` 去裸库仑项、暗激发/磁振子、自旋守恒选择定则

---

## 14. 核磁共振（NMR）

核磁共振（NMR）谱学是研究分子、液体与固体化学环境的高灵敏技术。第一性原理模拟可计算化学屏蔽、四极耦合常数与超精细常数，帮助理解 NMR 谱。

## 14.1 NMR（一）：化学屏蔽（收敛、诱导电流、与实验对比、预测）

> 来源：<https://vasp.at/tutorials/latest/nmr/part1/>
> 配套输入文件：[nmr-part1.zip](https://vasp.at/tutorials/latest/nmr-part1.zip)

本页四个例子：例 1 金刚石中化学屏蔽的收敛性，例 2 诱导电流成像（苯的环电流），例 3 计算屏蔽与实验化学位移对比，例 4 用理论预测实验。

**NMR 物理背景**：非零自旋 $\mathbf I$ 的原子核有磁矩 $\boldsymbol\mu$，绕外磁场 $\mathbf B_{ext}$ 以 Larmor 频率 $\omega_L$ 进动（由 $\mathbf B_{ext}$ 强度与旋磁比 $\gamma$ 决定）。施加垂直于恒定外磁场的射频脉冲 $\omega_{rf}\approx\omega_L$ 时发生磁共振，磁矩从参考系翻转到横向系综，随后弛豫产生的信号即 NMR 测量对象。

![核磁矩绕外磁场进动示意](/imgs/2026-08-05/fig42.png)

同位素核虽 $\gamma$ 相同，但化学环境不同——核外电子壳层**屏蔽**原子核（如溴乙烷中 -CH2- 的 1H 比 -CH3 更靠近 Br），改变核处磁场与共振频率。相对标准参考（四甲基硅烷 TMS）的频率移动即**化学位移** $\delta$：

$$\delta_{ij}(\mathbf R) = \sigma^{ref}_{ij} - \sigma_{ij}(\mathbf R)$$

**化学屏蔽张量**由线性响应计算：

$$\sigma_{ij}(\mathbf R) = -\frac{\partial B^{ind}_i(\mathbf R)}{\partial B^{ext}_j}$$

### 例 1：金刚石中化学屏蔽的收敛

**教学目标**：计算化学屏蔽；对平面波截断能与 k 网格收敛；对比屏蔽与能量的收敛行为。

**输入要点**（目录 `e01_shielding`，2 原子金刚石原胞，示例文件内容略）
- `LCHIMAG`：线性响应计算化学屏蔽；`LNMR_SYM_RED`：丢弃与 k 空间导数线性响应计算方式不一致的对称操作；`NLSPLINE`：保证倒空间 PAW 投影算符可微；`KPAR`：并行 k 点数加速；KPOINTS 为 4x4x4。

**收敛流程**
- 截断能：`encut_shielding.sh` 把 ENCUT 从 400 eV 以 100 eV 步长扫到 900 eV，从 OUTCAR grep 13C 屏蔽。stdout 先显示 DFT 收敛 CHGCAR/WAVECAR，随后是沿磁场轴的线性响应：`BDIR` 为外加磁场方向、`IDIR` 为感应磁场方向（1:x、2:y、3:z）；`IQ` 是 `ICHIBARE` 定义的有限差分模板网格点（默认 ICHIBARE=1，IQ 从 -1 到 +1 共 3 点）。
- OUTCAR 屏蔽输出字段：

| 标签 | 含义 |
|---|---|
| `EXCLUDING G=0 CONTRIBUTION` | 不含宏观长程贡献的屏蔽（分子有用） |
| `INCLUDING G=0 CONTRIBUTION` | 含宏观长程贡献（尤其表面电流）的屏蔽（晶体有用） |
| `ion` | 按 POSCAR 排序的原子编号 |
| `BDIR` | 外加磁场方向（BDIR=1 即沿 x），每行再分 X/Y/Z，即打印每个离子的完整屏蔽张量 |
| `iso_shield` | 各向同性屏蔽 $\sigma_{iso}=\mathrm{Tr}(\sigma)/3=(\sigma_{11}+\sigma_{22}+\sigma_{33})/3$（括号内还给出 span 与 skew） |
| `span` | $\Omega=\sigma_{33}-\sigma_{22}$ |
| `skew` | $\kappa=3(\sigma_{iso}-\sigma_{22})/(\sigma_{33}-\sigma_{11})=3(\sigma_{iso}-\sigma_{22})/\Omega$ |

  （定义见 Solid State Nucl. Magn. Reson. 2, 285 (1993)；这些屏蔽已含芯贡献，是与实验比较的最终重要项。）
- **收敛结论**：总能量 500 eV 即收敛到 10 meV 内，但屏蔽收敛慢得多——500 eV 可到 1 ppm，0.1 ppm 需 $\geq 800$ eV。精度要求与元素有关：1H 化学位移范围 0–14 ppm、13C 0–250 ppm，而 195Pt 达 $\pm 7000$ ppm，0.1 ppm 不现实、应选 1 ppm。
- **不建议对总能收敛**（此处仅为对比）；应对能量差或目标性质收敛。参数应**顺序收敛**：先截断能，再用收敛的截断能测 k 网格（分子体系则测真空大小）。
- k 网格：`kpoints.sh` 把 KPOINTS 从 4x4x4 以步长 2 扫到 16x16x16。能量 6x6x6 收敛到 ~10 meV，而屏蔽 0.1 ppm 需 14x14x14。**屏蔽对 k 点的依赖远大于截断能**：400→900 eV 只差 ~1 ppm，4x4x4→16x16x16 差近 90 ppm（近两个数量级）——布里渊区的精细描述比基组完备性更影响化学屏蔽。

### 例 2：诱导电流成像（苯的环电流）

**教学目标**：理解诱导电流的三部分构成；可视化苯的环电流。

**理论**：PAW 框架下诱导电流为三项之和（Phys. Rev. B 76, 024401 (2007)）：

$$\mathbf j^{(1)}(\mathbf r') = \mathbf j^{(1)}_{bare}(\mathbf r') + \mathbf j^{(1)}_{\Delta d}(\mathbf r') + \mathbf j^{(1)}_{\Delta p}(\mathbf r')$$

赝化价波函数贡献 + 抗磁项 + 顺磁项（后两者修正核附近赝波函数与全电子波函数的偏差）。

**输入与计算**（目录 `e02_induced_current`，示例文件内容略）
- INCAR 同例 1 加 `WRT_NMRCUR`：把电流写出为 `NMRCURBX`（x 轴；y/z 轴另有文件）；
- 孤立分子：大胞 + $\Gamma$ 点 k 网格；真空大小也是需收敛的参数。
- `plot_nmrcur_slice.py` 四参数：电流文件（NMRCURBX/Y/Z）、网格取样间隔 n、yz 截面位置（x 轴直接坐标）、箭头缩放指数。

![苯分子诱导电流切片图](/imgs/2026-08-05/fig43.png)

- **物理**：苯平面上下有离域 $\pi$ 电子环（芳香性）；磁场下沿环流动形成**环电流**，在环外产生与外场同向、环内反向的强磁场——质子（H）因此强烈**去屏蔽**，这是芳香化合物 NMR 谱的特征。也可用 py4vasp 画 NMR 电流。

### 例 3：与实验化学位移对比

**教学目标**：计算含 O 系列的化学屏蔽；把实验化学位移对计算屏蔽作图并线性拟合。

**理论**：化学屏蔽不是可观测量；NMR 测的是进动频率 $\omega_L=-\gamma B_z$。现代谱仪磁场 1–20 T，场越强自旋态能量差越大、信噪比越高，但样品频率也随之变化；为跨谱仪比较，取参考化合物频率 $\omega_{ref}$ 计算相对化学位移：

$$\delta = \frac{\omega_{ref}-\omega_{sample}}{\omega_{sample}}\times 10^6$$

计算上由线性响应得屏蔽：诱导电流经 Biot–Savart 定律给出感应磁场，分别得 bare/抗磁/顺磁屏蔽，加上化学不变的芯贡献：

$$\sigma_{tot} = \sigma_{bare} + \sigma_{\Delta d} + \sigma_{\Delta p} + \sigma_{core}$$

![自旋能级在外磁场下的分裂](/imgs/2026-08-05/fig44.png)

**输入与计算**（目录 `e03_experiment`：BaSnO3、MgO、SrO、SrTiO3 四个子目录，INCAR 同例 1、6x6x6 $\Gamma$ 中心网格、experiment.dat，示例文件内容略）
- 建议只跑 MgO（其余耗时）；`O_shielding.sh` 可跑全部四体系并把各向同性屏蔽存入 `O_shielding.dat`。
- OUTCAR 中 `CSA tensor` 下表：前三行为屏蔽张量、末行各向同性屏蔽；固体取 `INCLUDING G=0 CONTRIBUTION`（该行第 4 个数）。
- 计算中没有标准参考物；单个参考计算可能因方法局限扭曲结果，故把多个实验位移对计算屏蔽做线性拟合（平均掉个别计算的误差，J. Chem. Phys. 146, 064115 (2017)）：
$$\delta^{exp}_{iso}[N] = \delta_{ref}[N] + m\,\sigma^{calc}_{iso}[N]$$
- 教程四化合物拟合结果：$\sigma_{ref}=210.94$ ppm、$m=-0.8639$、Pearson $r=-0.981$。若 $|r|$ 低于 ~0.95，检查两数据文件中化合物顺序是否一致。

### 例 4：用理论预测实验

**教学目标**：计算 BaZrO3 与 CaO 的屏蔽；用线性拟合预测其化学位移并与实验对比；与更大文献系列比较。

- 目录 `e04_comparison`（BaZrO3/CaO 子目录，INCAR 同例 1、6x6x6 网格，示例文件内容略）；建议只跑 CaO（BaZrO3 很耗时）。
- 预测结果表：

| | BaZrO3 | CaO |
|---|---|---|
| 计算 $\sigma$ / ppm | -161.697 | -148.444 |
| 预测 δ / ppm | 350.639 | 339.188 |
| 实验 δ / ppm | 376 | 294 |
| 差值 / ppm | 25.0 | -45.0 |

- 偏差来源：BaZrO3 因 5s/5p 半芯态未冻结（浅芯态在晶场中可极化）；CaO 因已知效应——空 Ca 3d 态在 DFT 中离价带顶太近（J. Chem. Phys. 146, 064115 (2017)）。
- 拟合质量对比表：

| | $\sigma_\mathrm{ref}$ / ppm | m | r |
|---|---|---|---|
| 四化合物拟合 | 210.94 | -0.8639 | -0.981 |
| 六化合物拟合 | 208.07 | -0.8592 | -0.9248 |
| Laskowski（标准） | 220.65 | -0.8724 | -0.995 |
| Laskowski（优化） | 216.67 | -0.8558 | -0.996 |

- 四化合物系列截距/斜率接近 Laskowski 优化系列（差 <1 ppm 与 0.02），是好的拟合；加入 BaZrO3/CaO 后偏离文献值、Pearson 系数下降。Laskowski O 系列含 10 个化合物："标准"组用 700 eV 与更硬的 `_h` POTCAR（更小芯半径），"优化"组用 900 eV 与特制改进 POTCAR。
- 理想情况 $m=-1$（计算与实验只差参考位移 $\sigma_{ref}$）；$m$ 偏离 -1 源于 DFT 对交换关联效应的捕捉局限。
- **赝势选择很重要**：多数情况普通 POTCAR 足够；少数情况（如 BaZrO3）需要更小芯半径的硬赝势。

### 功能与标签汇总
- NMR 核心：`LCHIMAG`（线性响应屏蔽）、`LNMR_SYM_RED`、`NLSPLINE`、`ICHIBARE`（有限差分模板）、`KPAR`
- 电流成像：`WRT_NMRCUR` → NMRCURBX/Y/Z；PAW 三项电流分解（bare/抗磁/顺磁）+ 芯贡献
- 输出解读：EXCLUDING/INCLUDING G=0、CSA tensor、iso_shield/span/skew
- 方法论：对目标性质（而非总能）顺序收敛（ENCUT→k 网格→真空）；屏蔽对 k 网格更敏感；多化合物线性拟合（$\sigma_\mathrm{ref}$、m、r）替代单一参考；PAW 赝势选择（半芯态、硬赝势）

---

## 14.2 NMR（二）：耦合常数与双中心修正

> 来源：<https://vasp.at/tutorials/latest/nmr/part2/>
> 配套输入文件：[NMR-part2.zip](https://vasp.at/tutorials/latest/NMR-part2.zip)

本页四个例子：例 5 超精细常数（$\mathrm{CH_3\cdot}$ 自由基），例 6 电场梯度（MAPbI3 的 EDIFF 依赖），例 7 四极耦合常数收敛，例 8 双中心屏蔽贡献（LiH）。

### 例 5：超精细常数

**教学目标**：做超精细耦合计算；与实验和文献对比；比较 GGA 与杂化泛函的影响。

**任务**：在 10x10x10 Å$^3$ 胞中用 PBE 与 HSE06 计算气相 CH3$^\bullet$ 自由基的超精细耦合张量。

**理论**：电子也有自旋；电子内禀磁场与核磁偶极矩相互作用使简并能级分裂——**超精细分裂**。哈密顿量 $\hat H_D=-\boldsymbol\mu_I\cdot\mathbf B$，磁场为电子轨道与自旋角动量贡献之和 $\mathbf B=\mathbf B^l_{el}+\mathbf B^s_{el}$；核自旋 $S^I$ 与电子自旋 $S^e$ 经超精细张量 $A^I$ 耦合（Phys. Rev. 88, 075202 (2013)）：

$$E = \sum_{ij}S^e_i A^I_{ij}S^I_j$$

**输入要点**（目录 `e05_hyperfine`，示例文件内容略）
- `LHYPERFINE=.TRUE.` 计算超精细耦合张量；`NGYROMAG` 给出各离子的核旋磁比（此处选 13C 与 1H）；有未配对电子故 `ISPIN=2`；
- INCAR.hse06：`LHFCALC=.TRUE.`+`GGA=PE`+`HFSCREEN=0.2` 设置 HSE06；`ALGO=Damped`（LHFCALC 推荐）；KPOINTS 为 $\Gamma$ 点。

**输出解读**（OUTCAR）
- 超精细参数分为各向同性的 **Fermi 接触**项与各向异性的**偶极**贡献。`Fermi contact (isotropic) hyperfine coupling parameter (MHz)` 之后按 POSCAR 离子列出，分量含义：

| 标签 | 含义 |
|---|---|
| `A_pw` | Fermi 接触项的平面波贡献 |
| `A_1PS` | 赝单中心贡献 |
| `A_1AE` | 全电子单中心贡献 |
| `A_1c` | 单中心芯贡献（Phys. Rev. 71, 115110 (2005)） |
| `A_tot` | Fermi 接触项总贡献（不含 A_1c） |

  随后是偶极贡献张量 $A_{ij}$，对角化后写出总超精细张量。$\mathrm{CH_3\cdot}$ 有三重对称，三个 1H 核等价（数值几乎相同）。
- **收敛与对比**：与化学屏蔽不同，超精细参数对截断能收敛很快（500 eV）；固体应用还需收敛 k 网格（如金刚石 NV- 缺陷）。
- 文献 $A_{tot}=182.5$ MHz；教程值与文献之差来自 `LASPH=.TRUE.`（PAW 球内梯度的非球贡献，文献未含；关掉即与文献一致，仅胞大小 10 vs 20 Å 有小差异）。
- 实验值（EPR 电子顺磁共振测得）107.4 MHz；不含芯贡献的 $A_{tot}$ 与实验差很多，加上芯贡献（$A_{tot}+A_{1c}$）文献值变 86.1 MHz（接近实验），教程值 ~70 MHz，仍差 20–30 MHz。
- **HSE06 显著改善**：文献 $A_{tot}+A_{1c}=101.9$ MHz，明显更接近实验 107.4 MHz。原因：杂化泛函含精确交换，有助于局域化孤电子或缺陷（J. Chem. Phys. 125, 224106 (2006)）。

### 例 6：电场梯度（EFG）

**教学目标**：计算 EFG；理解 EDIFF 对 EFG 的影响。

**任务**：计算钙钛矿太阳能电池材料 MAPbI3（甲基铵铅碘）中 I 与 N 的核电四极矩，考察收紧 SCF 收敛判据 EDIFF 的影响。

**理论**：EFG 是静电势 $V$ 在核位置的二阶导数 $V_{ij}=\partial^2 V/\partial x_i\partial x_j$。EFG 本身不可测，但核四极耦合常数 $C_q=eQV_{zz}/h$ 可测（NMR、EPR 或核四极共振 NQR）；$Q$ 为元素与同位素特定的四极矩。

**输入要点**（目录 `e06_efg_ediff`，示例文件内容略）
- INCAR 同例 1 加两个标签：`LEFG` 开启 EFG 计算；`QUAD_EFG` 按 POTCAR 原子类型定义核四极矩（单位毫靶恩 millbarn）；KPOINTS 仅 $\Gamma$ 点。

**计算与结果**
- `ediff_efg.sh` 把 EDIFF 从 1E-4 扫到 1E-8 eV；OUTCAR 中 `V_ii` 为沿 i 轴的 EFG 分量，再往下是四极耦合常数。
- MAPbI3 结构：I 分两类——轴向 I（与 Pb 同层、每层 2 个，与 CH3NH3 甲基铵 MA 层同层）与赤道 I（MA 层之间、每层 4 个），I 在 Pb 周围构成八面体。
- Cq 随 EDIFF 变化表（log EDIFF 4→8）：轴向/赤道 I 从 1E-4 到 1E-5 eV 变化很大，之后 <0.1 MHz 视为收敛；N 在 1E-4 已收敛（0.471 MHz 不变）。**EFG 常受 EDIFF 影响强烈**；耦合常数对技术参数的敏感度因体系而异，每组新计算都要做收敛测试。
- 另注：POTCAR 选择对四极耦合常数影响可能很大（GW 赝势、`_pv`/`_sv` 半芯电子），收敛前应先检查。

### 例 7：四极耦合常数的收敛

**教学目标**：对 k 点与截断能收敛 EFG；为不同元素确定合适的收敛限。

**任务**：对四方相 MAPbI3 计算 14N（丰度 99.9(6)%）与 127I（100%）的核四极矩，收敛 k 网格与截断能并尽量对比实验（NQR 可观测固体中这些核）。

**计算与结果**（目录 `e07_efg_convergence`，输入为例 6 的 INCAR + 2x2x2 $\Gamma$ 中心网格，示例文件内容略）
- 截断能：`encut_efg.sh` 把 ENCUT 从 200 eV 以 100 eV 步长扫到 600 eV（500/600 eV 结果已提供）。**注意**：实际计算绝不应使用低于 POTCAR 中 ENMAX 的 ENCUT，此处仅为演示收敛。
- 结果表显示：赤道 I 收敛快（300 eV 内到 ~0.1 MHz），轴向 I 类似；N 因数值小需要更严的限（0.01 MHz），400–500 eV 才完全收敛——推荐至少 400 eV、最好 500 eV。
- k 网格（400 eV 固定，实际应顺序收敛用 500 eV）：`kpoints_efg.sh` 扫 1x1x1、2x2x2（另加预制的 4x4x4、6x6x6）。$\Gamma$ 点远未收敛；N 在 2x2x2 收敛到 0.001 MHz 内，而 I 需要 6x6x6 才到 1 MHz 内——**不同元素收敛严格度不同，四极耦合常数的大小是关键因素**（N 0.001 MHz vs I 1 MHz）。
- **与实验对比的教训**：Franssen et al.（J. Phys. Chem. Lett. 8, 61 (2017)）NQR 测得 14N 的 Cq 从 -100 到 75 °C 由 0.11 降到 0.00 MHz，与计算值很不同；但论文用实验结构算得甲胺 Cq=0.45 MHz（接近教程值）。差异原因：甲胺在毫秒量级重新取向（NMR 时间尺度上很快），取四种构型平均后计算 Cq 随温度从 0.09 降到 0.00 MHz，与实验趋势一致——**收敛的计算只是与实验对比的一部分，还需考虑未预期的动力学**。

### 例 8：双中心屏蔽贡献

**教学目标**：计算 LiH 有/无双中心修正的化学屏蔽；对截断能收敛；比较双中心贡献对 Li 与 H 的影响。

**背景**：化学屏蔽通常被当作单中心性质（忽略相邻 PAW 球的增补电流贡献），一般成立；但**双中心贡献**有时很大（J. Chem. Phys. 139, 014109 (2013)）——尤其对周期表顶部元素（H、B、C、N、O、F）使用硬赝势（`*_h`）且键短时。

**输入与计算**（目录 `e08_two_center`，15x15x15 Å$^3$ 胞中的 LiH，示例文件内容略）
- INCAR 同例 1；`LLRAUG` 开启屏蔽张量的双中心贡献。
- `llraug.sh` 在 400 eV 分别跑 LLRAUG=.FALSE./.TRUE.（500/600 eV 结果后补），屏蔽写入 `llraug.dat`。
- 结果表（600 eV）：$\sigma$(Li) 88.8883→88.8918 ppm（差 0.0035 ppm，几乎可忽略）；$\sigma$(H) 24.6947→26.2643 ppm（差 ~1.1–1.6 ppm，大三个数量级）——凸显双中心贡献的重要性。
- 双中心修正只源于 PAW 球的使用；完整原子中心基组下屏蔽为 Li 88.12 ppm、H 26.34 ppm（Phys. Chem. Chem. Phys. 18, 21145 (2016)）——与含双中心修正的 PAW 结果吻合良好。

### 功能与标签汇总
- 超精细：`LHYPERFINE`、`NGYROMAG`（核旋磁比）、ISPIN=2；Fermi 接触/偶极分解（A_pw/A_1PS/A_1AE/A_1c/A_tot）；LASPH 的影响；杂化泛函局域化孤电子
- EFG/四极：`LEFG`、`QUAD_EFG`（毫靶恩）；EFG 对 EDIFF 极度敏感；Cq=eQVzz/h；POTCAR 选择（GW、_pv/_sv）
- 收敛认识：不同元素收敛限不同（按 Cq 量级定）；ENCUT 不应低于 ENMAX；动力学平均对实验对比的必要性
- 双中心：`LLRAUG`；顶部行元素+硬赝势+短键时重要；H 比 Li 大三个数量级

---

## 14.3 Part 3：芳香性（Aromaticity）

- 来源：<https://vasp.at/tutorials/latest/nmr/part3/>
- 教程输入文件下载：[NMR-part3.zip](https://vasp.at/tutorials/latest/NMR-part3.zip)

本部分把 NMR 屏蔽计算推广到**芳香性判据**：通过核独立化学位移（NICS, nucleus-independent chemical shift）与诱导电流的方向，从磁学角度区分芳香（苯）与反芳香（环丁二烯）分子。

---

### 示例 9：苯的 NICS 等高线图（NICS - benzene）

**教学目标**

- 计算核独立化学位移（NICS）；
- 绘制苯分子平面内的 NICS 等高线图；
- 用磁学判据判断芳香性。

**理论背景**

- 外加磁场 $\textbf{B}_{ext}$ 诱导出电流 $\textbf{j}^{(1)}$，该电流又产生诱导磁场 $\textbf{B}_{in}^{(1)}$。
- 芳香分子中共轭双键形成的环流（ring current）使环外质子受到强烈去屏蔽（这是芳香化合物的标志性 NMR 特征），而环内部则被强烈屏蔽。
- 空间中任意一点（不必在原子核上）的化学屏蔽称为 NICS，参见 [Chem. Rev. 105 (2005) 3842](https://doi.org/10.1021/cr030088+)。常用环中心单点 NICS 判断芳香性，但 NICS 等高线图能更直观地展示诱导磁场分布，参见 [J. Phys. Chem. A 117 (2013) 518](https://doi.org/10.1021/jp311536c)。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/NMR/e09_nics_benzene`）

- 在 $8\times8\times8$ Å$^3$ 盒子中计算苯分子平面内的 NICS 等高线图；
- INCAR 与示例 1 类似，额外设置：
  - `NUCIND = .TRUE.`：开启 NICS 计算；
  - `LPOSNICS = .TRUE.`：在 **POSNICS** 文件指定的位置计算 NICS；
- POSNICS 文件先给出 NICS 点数，然后逐行给出直接坐标（x、y、z）。本例取 xy 平面内、z 方向位于晶胞中央的一组网格点（教程提供脚本自动生成）；
- KPOINTS 使用 $\Gamma$ 点网格。

**结果分析**

- 从 OUTCAR 中提取各 NICS 点的屏蔽值，结合 POSNICS 中的坐标绘制等高线图：

![苯分子平面的 NICS 等高线图](/imgs/2026-08-05/fig45.png)

- 图中 C 核周围呈红色去屏蔽区，H 核周围呈蓝色屏蔽区；
- 最清晰的芳香性证据出现在环中心：环内呈现强屏蔽，这源于示例 2 中计算过的抗磁环流（diamagnetic current）。

**思考题**

1. 用哪个文件定义 NICS 的计算位置？（POSNICS）
2. NICS 结果在哪个输出文件中？（OUTCAR）

---

### 示例 10：反芳香电流（Antiaromatic current）

**教学目标**

- 计算并可视化环丁二烯的诱导电流；
- 利用电流方向判断环丁二烯的（反）芳香性。

**理论背景**

- 芳香性的一眼判据是数环上的 $\pi$ 电子数：芳香环有 $4n+2$ 个 $\pi$ 电子（苯有 6 个，$n=1$）；反芳香环有 $4n$ 个 $\pi$ 电子。
- $n=1$ 的反芳香体系即环丁二烯 C4H4（四元环）。芳香性使苯稳定为平面六元环；反芳香性使环丁二烯不稳定，从正方形畸变为矩形——$\pi$ 电子不再离域于整个环，而是定域在两个双键上。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/NMR/e10_current_c4h4`）

- 在 $8\times8\times8$ Å$^3$ 盒子中计算气相矩形环丁二烯的诱导电流；
- INCAR 标签与示例 2 完全相同（含 `WRT_NMRCUR` 等电流输出设置）；
- KPOINTS 使用 $\Gamma$ 点网格。

**结果分析**

- 用 py4vasp 查看分子结构，可与苯的六元环明确区分；
- 用 `plot_nmrcur_slice_c4h4.py`（与示例 2 的 `plot_nmrcur_slice.py` 相同，仅改原子位置）提取并绘制电流切片：

![环丁二烯电流密度切片](/imgs/2026-08-05/fig46.png)

- 黄色区域清楚显示双键对电流密度的影响（电流定域化）；
- 苯与环丁二烯的平面波贡献电流都看似逆时针绕环流动，但关键区别在于：环丁二烯实际产生**顺时针的顺磁（paratropic）电流**，苯产生**逆时针的抗磁（diatropic）电流**（[Chem. Commun. 21 (2001) 2220](https://doi.org/10.1039/B104847N)）。环流方向是芳香/反芳香的指示器：逆时针为芳香，顺时针为反芳香。
- 为什么图中看不到顺时针电流？因为绘图只画出了**平面波贡献**的电流（易于在精细 FFT 网格上表达）。PAW 方法把电流拆分为平面波部分与单中心（one-center, 1c）部分：抗磁 1c 贡献很小，但**顺磁 1c 贡献很大**，而 1c 贡献处于"局域"基组中、不在 FFT 网格上，因此未被画出。
- 若设置 `LNMRLEG = .TRUE.`，可在 OUTCAR 中查找 `plane wave contribution`、`one center paramagnetic contribution`、`one center diamagnetic contribution` 三项，直接比较各贡献的数值。

**思考题**

1. 电流结果在哪个输出文件中？
2. 电流在什么网格上计算？（精细 FFT 网格）
3. 芳香环的电流方向如何？属于什么电流？（逆时针，抗磁）
4. 反芳香环的电流方向如何？为什么图中看不到？

---

### 示例 11：环丁二烯的 NICS（NICS - cyclobutadiene）

**教学目标**

- 计算环丁二烯的 NICS 并绘制等高线图；
- 用磁学判据区分芳香性与反芳香性。

**理论背景**

- 芳香环中电流逆时针（抗磁），反芳香环中电流顺时针（顺磁）；由 [Biot–Savart 定律](https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law)可知两者诱导的磁场方向相反：

![芳香与反芳香环流产生的诱导磁场示意](/imgs/2026-08-05/fig47.png)

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/NMR/e11_nics_c4h4`）

- 在 $8\times8\times8$ Å$^3$ 盒子中计算环丁二烯分子平面的 NICS 等高线图；
- INCAR 与示例 9 基本相同，但**不使用 POSNICS/LPOSNICS**：只设 `NUCIND = .TRUE.`（等价于 `LNICSALL = .TRUE.`），NICS 将在**精细 FFT 网格的所有点**上计算，结果写入专门的 **NICS 输出文件**；
- KPOINTS 使用 $\Gamma$ 点网格。

**结果分析**

- 从 NICS 文件提取数值并绘制等高线图：

![环丁二烯分子平面的 NICS 等高线图](/imgs/2026-08-05/fig48.png)

- 与苯的 NICS 图并排对比：
  - H 核周围的蓝色屏蔽、C 核周围的红色去屏蔽两者相似；
  - 主要区别在化学键与环内：苯的 C=C 双键区域（去屏蔽晕之间的屏蔽区）因芳香稳定化而明显更屏蔽；
  - 环内差异最显著：**芳香的苯环内为屏蔽，反芳香的环丁二烯环内为去屏蔽**。
- 物理图像（Biot–Savart 定律）：环丁二烯的顺时针顺磁电流在环内诱导与外场同向的磁场 → 环内去屏蔽、环外屏蔽；苯则相反，逆时针抗磁电流在环内诱导与外场反向的磁场 → 环内屏蔽、环外去屏蔽。（注意：所绘电流仅为平面波贡献。）
- 两种 NICS 计算方式的取舍：
  - 精细 FFT 网格上计算 NICS 更快且可并行（`LNICSALL`）；
  - POSNICS 方式允许在感兴趣的局部区域（如化学键、氢键）布置比精细 FFT 更密的网格，适合研究局部化学环境对屏蔽的影响，例如聚合物中的氢键（[Macromolecules 49 (2016) 5548](https://doi.org/10.1021/acs.macromol.6b01051)）。

**思考题**

1. `LNICSALL = .TRUE.` 在哪些位置计算 NICS？（精细 FFT 网格全部点）
2. 为什么 `LNICSALL` 比 POSNICS 方式更快？（可并行）
3. POSNICS 自定义位置有什么优势？（可在关注区域使用更密网格）

---

### 功能与标签汇总（NMR Part 3）

| 功能 | 关键标签 / 文件 | 说明 |
|---|---|---|
| NICS 计算 | `NUCIND = .TRUE.` | 开启核独立化学位移计算 |
| 指定位置 NICS | `LPOSNICS = .TRUE.` + POSNICS | 在用户给定坐标（可更密网格）计算 NICS，结果在 OUTCAR |
| 全网格 NICS | `LNICSALL = .TRUE.`（有 NUCIND、无 POSNICS 时默认） | 在精细 FFT 网格所有点计算，写入 NICS 文件，更快可并行 |
| 磁化率/电流各贡献 | `LNMRLEG = .TRUE.` | OUTCAR 输出平面波、单中心顺磁/抗磁贡献 |
| 芳香性磁判据 | NICS 环中心值、环流方向 | 环内屏蔽+逆时针抗磁流=芳香；环内去屏蔽+顺时针顺磁流=反芳香 |

---

## 15. 声子（Phonons）

声子是扩展周期体系中原子核的集体激发。

## 15.1 Part 1：石墨烯声子（Graphene）

- 来源：<https://vasp.at/tutorials/latest/phonon/part1/>
- 教程输入文件下载：[phonon-part1.zip](https://vasp.at/tutorials/latest/phonon-part1.zip)

本部分以石墨烯为例，介绍声子计算的完整流程：结构弛豫 → 有限差分法计算力常数 → 声子色散关系，并讨论超胞尺寸、对称性与收敛性。

---

### 示例 1：静态晶格近似（Static-lattice approximation）

**教学目标**

- 列举静态晶格近似的局限；
- 对二维材料在固定体积下弛豫离子位置与晶胞形状；
- 设置 OUTCAR 的详细程度（verbosity）；
- 用 py4vasp 可视化晶体结构。

**理论背景**

- 标准第一性原理计算采用静态晶格模型：离子构成固定、刚性的周期结构，电子与离子自由度通过 Born–Oppenheimer 近似分离。典型流程：选取初始结构 → 量子力学处理电子 → 依据 Hellmann–Feynman 定理得到的力更新离子位置。
- 静态晶格模型的局限：无法描述热膨胀、热导、压电性、熔化、声波传播、X 射线衍射的温度依赖强度、中子散射的某些方面以及部分光学/介电性质（如反射率在远低于带隙处出现共振）；许多材料中声子对比热的贡献大于电子。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/phonon/e01_static-lattice`）

- 在固定体积下弛豫石墨烯的离子位置与晶胞形状；
- POSCAR：单层碳原子蜂窝晶格，z 方向留 8 Å 真空层抑制周期性镜像间相互作用；selective dynamics 固定 z 方向弛豫；
- INCAR：电子与离子自由度迭代更新（Kohn–Sham 基态 → Hellmann–Feynman 力与应力 → `IBRION` 决定的结构优化算法，本例为共轭梯度 → 重复直到满足 `EDIFFG` 收敛判据）；
  - `EDIFF` 必须足够小，否则电子密度收敛不足会导致力/应力不准，可能虚假收敛或收敛到错误解；
  - `ISIF = 4`：固定体积、允许改变晶胞形状与离子位置。虽然体积固定，晶格参数仍可变化，满足 `ENCUT` 的倒格矢数目可能改变 → 产生 **Pulay 应力**；对策是若干离子步后重启计算，或提高 `ENCUT`；
  - `NWRITE` 控制 OUTCAR 输出详细程度，本例借短计算机会观察额外输出；
- KPOINTS：$\Gamma$ 中心等距 k 网格；单层材料无需在表面法向细分 k 网格（法向 k 点只描述镜像层间相互作用）。

**结果分析**

- 用 py4vasp 可视化弛豫后的石墨烯，可绘制 3x3x1 超胞观察典型石墨烯层状图案；
- `EDIFFG` 默认为 10x`EDIFF`（即 $1\times10^{-7}$）；在 OUTCAR 中搜索 `TOTAL-FORCE (eV/Angst)` 可见每步力均为零，说明初始结构已是（亚）稳态；
- 后续示例沿用 POSCAR 中对称化的结构。

**思考题**

1. 静态晶格近似有哪些局限？
2. 垂直表面方向 k 网格应取多少细分？为什么？
3. `NWRITE` 控制什么？
4. 要变体积计算需改哪个标签？（`ISIF`）

---

### 示例 2：有限差分法计算力常数（Force constants using finite differences）

**教学目标**

- 为声子计算构建超胞；
- 给出力常数定义并用有限差分法计算；
- 判断何时设置 `LREAL = Auto`；
- 解释对称性在有限差分算法中的作用。

**理论背景**

- Born–Oppenheimer 能量面 $E(\{\mathbf{R}\})$ 依赖离子位置 $\mathbf{R}_I = \mathbf{R}^0_I + \mathbf{u}_I$。围绕平衡位置展开：

$$E(\{\mathbf{R}\}) = E(\{\mathbf{R}^0\}) - \sum_{I\alpha} F_{I\alpha}(\{\mathbf{R}^0\}) u_{I\alpha} + \sum_{I\alpha J\beta} \Phi_{I\alpha J\beta}(\{\mathbf{R}^0\}) u_{I\alpha} u_{J\beta} + \mathcal{O}(\mathbf{R}^3)$$

- 一阶系数为力 $F$，二阶系数为 Hessian / 力常数矩阵 $\Phi$（谐波近似截断到二阶）：

$$F_{I\alpha} = -\left. \frac{\partial E}{\partial R_{I\alpha}} \right|_{\mathbf{R}=\mathbf{R}^0}, \qquad \Phi_{I\alpha J\beta} = \left. \frac{\partial^2 E}{\partial R_{I\alpha} \partial R_{J\beta}} \right|_{\mathbf{R}=\mathbf{R}^0} = - \left. \frac{\partial F_{I\alpha}}{\partial R_{J\beta}} \right|_{\mathbf{R}=\mathbf{R}^0}$$

- 要捕捉波长大于原胞的声子模式，必须使用足够大的超胞；**所有声子计算都应对超胞尺寸收敛**。超胞放大 n 倍时，k 网格按同倍数缩小，保持采样密度不变。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/phonon/e02_fin-diff-force-constants`）

- 计算石墨烯的力常数；用 Python 脚本生成面内 3x3、4x4、5x5、6x6 超胞；
- INCAR：`IBRION = 6` 触发有限差分法计算力常数与声子模式；
- 讨论 `LREAL` 标签：投影算符（投影到局域的部分）涉及 aliasing 问题；本例应使用 `LREAL = Auto`（配合自动确定的投影网格）；
- KPOINTS：各超胞对应的等密度 k 网格。

**结果分析**

- py4vasp 可视化基态结构与有限差分所用的微扰结构：结构 0 为未微扰结构，其余结构中某一个原子在"抖动"；VASP 利用对称性，仅靠这些少量位移即可提取整个体系的力常数；
- 大超胞也只有 2 个自由度（而非 3x原子数）：石墨烯六方结构中面内晶格矢量等长、夹角 $60^\circ$，绕第三晶格矢量的 6 重旋转对称把 6 个自由度降到 4（每个位点仅一个面内自由度），反演对称再减 2（只需位移一个位点）；
- `IBRION = 5` 与 `IBRION = 6` 的区别：6 利用对称性，5 不利用；
- 力常数矩阵可经 py4vasp 读取。石墨烯等二维材料中，面内（x、y）与面外（z）力常数不耦合，相应对角块外元素为零；
- 把面外力常数、面内力常数（旋转到纵向 $\hat{\mathbf{r}}_L$ / 横向 $\hat{\mathbf{r}}_T$ / 法向的局域参考系 $\tilde{\Phi} = R^T \Phi R$）对原子间距作图：增大超胞尺寸可获得更多原子对的力常数；
- 由于远距离原子对的力常数幅度衰减，声子色散最终会随超胞尺寸收敛。两个著名例外：
  1. 带电离子位移形成偶极，其相互作用衰减缓慢（非解析修正，Part 2 处理）；
  2. 金属/半金属中核运动改变费米能级进而改变总能，形成长程的 [Kohn 反常](https://en.wikipedia.org/wiki/Kohn_anomaly)。
- 注意：出现少量负力常数不足以判定亚稳——某对离子位移可能降低能量，但其他正力常数可抵消之，整体仍稳定。

**思考题**

1. `IBRION = 5` 和 `6` 的区别？
2. 力常数的定义是什么？
3. `LREAL` 控制什么？
4. 超胞尺寸为何对精确声子色散重要？

---

### 示例 3：声子色散（Phonon dispersion）

**教学目标**

- 计算并绘制石墨烯声子色散。

**理论背景**

- 谐波近似下离子哈密顿量给出运动方程 $M_I \ddot{u}_{I\alpha} = -\Phi_{I\alpha J\beta} u_{J\beta}$。设平面波形式的解，求声子模式 $\varepsilon^{\mu}_{I\alpha}(\mathbf{q})$ 与频率 $\omega^\mu(\mathbf{q})$ 归结为动力学矩阵本征值问题：

$$\sum_{J\beta} \frac{1}{\sqrt{M_I M_J}} \Phi_{I\alpha J\beta} e^{i\mathbf{q}\cdot(\mathbf{R}_J-\mathbf{R}_I)} \varepsilon^{\mu}_{J\beta}(\mathbf{q}) = \omega^\mu(\mathbf{q})^2 \varepsilon^{\mu}_{I\alpha}(\mathbf{q})$$

- 超胞含 N 个离子时动力学矩阵维度为 3N；转换到原胞基（$\mathbf{R}_I = \mathbf{L}_i + \mathbf{r}_i$）后维度为 3n（n 为原胞原子数）；
- 与超胞**可公度**的 q 点（$\mathbf{q} = \frac{s_1}{n_1}\mathbf{b}_1 + \frac{s_2}{n_2}\mathbf{b}_2 + \frac{s_3}{n_3}\mathbf{b}_3$）上声子频率是精确计算的，其余 q 点由实空间力常数插值得到。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/phonon/e03_dispersion`）

- 计算石墨烯声子色散；
- 把示例 2 的 `vaspout.h5` 与 POSCAR 复制为 `vaspin.h5`（VASP 可读的力常数格式）；
- QPOINTS 文件定义色散路径上的 q 点序列；
- INCAR 中设置声子插值相关标签（示例要求在文件中逐个注释其含义）。

**结果分析**

- 用 py4vasp 绘制各超胞尺寸的声子色散并叠图比较：
  - 超胞过小时 $\Gamma$ 点附近出现鞍点（轻微虚频）；这类二维材料的二次型"软"声子模式收敛特别困难，需格外小心；
  - 负声子频率表示结构不稳定：平衡结构力为零，但沿某声子模式畸变会降低总能（亚稳）。虚频有时是计算不准（超胞不够、`ENCUT` 过低）或截断力常数插值的假象，有时有物理意义（沿该模式弛豫得到更低对称性、更低能量的结构）。本例 $\Gamma$ 附近的负模式通常指示计算问题；
  - 5x5 超胞的 $\Gamma$ 点光学支频率与其他超胞不同，因为 3x3/4x4/6x6 保持了等密度 k 采样而 5x5 只能使用更粗 k 网格；
- 倒空间布里渊区示意：对给定超胞，仅在蓝色采样点（$\Gamma$ 中心网格对应点）上频率精确，其余为插值：

![超胞尺寸对应的布里渊区采样点](/imgs/2026-08-05/fig49.png)

- 计算 M 点频率的最小超胞为 2x2，计算 K 点的最小超胞为 3x3；
- 用 [phononwebsite](https://henriquemiranda.github.io/phononwebsite) 交互式可视化声子模式：从 `vaspout.h5` 生成 `graphene.json`，在网页中可选择原胞重复次数、调节振荡幅度、显示原子运动方向矢量；点击色散曲线上任一点即可在中央面板播放对应声子模式的动画。

**思考题**

1. $\mathbf{q}=0$ 声子的波长在晶体中延伸多远？
2. 什么是声子色散？
3. 力常数与声子频率如何联系？

---

### 功能与标签汇总（Phonon Part 1）

| 功能 | 关键标签 / 文件 | 说明 |
|---|---|---|
| 固定体积形状弛豫 | `ISIF = 4` | 体积固定，允许晶胞形状与离子位置变化；注意 Pulay 应力 |
| OUTCAR 详细度 | `NWRITE` | 控制输出信息量 |
| 有限差分歧变 | `IBRION = 6` / `IBRION = 5` | 6 利用对称性减少位移数，5 不用对称性 |
| 投影局域化 | `LREAL = Auto` | 大超胞/有限差分推荐，配合自动投影网格 |
| 力常数读取 | `vaspout.h5` → `vaspin.h5` | h5 文件传递力常数 |
| 声子插值 | QPOINTS | 指定色散路径 q 点；可公度点精确、其余插值 |
| 模式可视化 | phononwebsite | 导出 json 后在线播放声子模式动画 |

---

## 15.2 Part 2：MgO 声子与 LO-TO 分裂

- 来源：<https://vasp.at/tutorials/latest/phonon/part2/>
- 教程输入文件下载：[phonon-part2.zip](https://vasp.at/tutorials/latest/phonon-part2.zip)

本部分以离子晶体 MgO 为例，在 Part 1 流程之上补充：晶格常数的 k 点收敛、声子态密度（DOS）、以及极性材料长程偶极–偶极相互作用的非解析处理（LO-TO 分裂）。

---

### 示例 4：晶格常数（Lattice parameter）

**教学目标**

- 执行晶胞体积弛豫；
- 监测晶格常数随 k 点网格密度的收敛。

**要点**

- 精确晶格常数需要足够大的基组（`ENCUT`）与足够的布里渊区采样；
- 推荐做法：观察目标量（本例为弛豫晶格常数）随 k 点密度的变化，选取能给出精确结果的最小密度，节省计算资源——因为后续超胞有限差分计算昂贵得多；
- 教程脚本对 N = 4、6、8、10、12 生成 NxNxN KPOINTS，逐个运行 VASP 并把 `vaspout.h5` 存为 `vaspout_$N.h5`；随后用 py4vasp 读取各文件的晶格常数并对 k 网格作图，判断收敛网格；
- 可直接从收敛计算的结果生成后续 POSCAR。

**思考题**：为何要对 k 点密度收敛晶格常数？

---

### 示例 5：有限差分法计算力常数

**教学目标**

- 为体材料声子计算构建超胞；
- 用有限差分法计算力常数。

**要点**

- 与二维材料不同，体材料沿三个笛卡尔方向扩大原胞；立方对称允许三方向共用同一缩放因子；
- 脚本生成 2x2x2、3x3x3、4x4x4 超胞，并生成对应于原胞 6x6x6 采样的等密度 KPOINTS；
- INCAR 设 `IBRION = 6`，对各超胞分别运行；
- 3x3x3 与 4x4x4 耗时较长，教程建议先用 2x2x2 继续后续章节（但要记住它对精确色散明显不足），大超胞在后台运行，最后可用更大力常数重跑第 6、8 节。

**思考题**：体材料与二维材料构建力常数超胞的主要区别是什么？

---

### 示例 6：声子色散与态密度（Phonon dispersion and DOS）

**教学目标**

- 计算 MgO 的声子色散与声子 DOS。

**要点**

- 利用上一步超胞力常数（写入 `vaspout.h5`）：把其复制为 `vaspin.h5` 后，INCAR 中设置 `LPHON_READ_FORCE_CONSTANTS`，VASP 会从 `vaspin.h5` 读入力常数，在 QPOINTS 指定的 q 点上通过傅里叶插值构建动力学矩阵、绘出色散后退出；
- **声子 DOS 的 q 点收敛**：分别在 10x10x10、20x20x20、30x30x30、40x40x40 q 网格上插值计算 DOS 并比较：
  - 10x10x10 与 20x20x20 之间仍有明显差别；30x30x30 与 40x40x40 差别可忽略 → DOS 已收敛；
  - 插值非常快：2x2x2 力常数在 40x40x40 q 点上仅需约 5 秒；
- **原子分解的声子 DOS**：Mg 原子在低频声子模式中权重最大，O 原子在高频部分权重最大。原因：O 质量小于 Mg，其参与振动的模式频率更高；频率大小还取决于力常数强度（由声子振动模式所形变化学键的强度决定）。

**思考题**

1. 与超胞不可公度的 q 点用什么插值得到声子模式？（力常数傅里叶插值）
2. 决定某声子模式频率大小的两个主要量是什么？（原子质量与力常数）

---

### 示例 7：长程偶极–偶极相互作用（LO-TO splitting）——计算介电张量与 Born 有效电荷

**教学目标**

- 计算 MgO 的静态介电张量与 Born 有效电荷。

**理论背景**

- 半导体/绝缘体中电子对离子的屏蔽不完全，产生**长程（LR）原子间力常数**。显式计算需要无穷大超胞；有限尺寸截断会在声子色散中引起 Gibbs 振荡；
- 解决办法（[Gonze and Lee 1997](https://doi.org/10.1103/PhysRevB.55.10355)）：把二阶力常数分解为短程与长程两部分：

$$\Phi_{I\alpha J\beta} = \Phi_{I\alpha J\beta}^\text{SR} + \Phi_{I\alpha J\beta}^\text{LR}$$

- 长程部分来自总能量中离子–离子贡献的长程部分的解析导数，用 Ewald 求和技术评估（实空间项对应短程、倒空间项对应长程）；长程力常数由 **Born 有效电荷张量** 与 **（离子钳制）静态介电张量** 完全描述。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/phonon/e07_mgo_dielectric_becs`）

- INCAR 设 `LEPSILON = .TRUE.`，计算钳制离子静态介电张量与 Born 有效电荷（配合 `PREC = Accurate`、`EDIFF = 1e-8`、`ISMEAR = 0`、`LREAL = .FALSE.` 等）；
- KPOINTS：$\Gamma$ 中心 6x6x6 网格。

**结果分析**

- 从 `vaspout.h5` 用 py4vasp 读取介电张量与 Born 有效电荷；
- 教程提供 Python 代码把 `vaspout.h5` 中的信息写成 INCAR 可用的 `PHON_DIELECTRIC` 与 `PHON_BORN_CHARGES` 标签，供下一节使用。

**思考题**

1. 力常数长程相互作用的起源是什么？（离子屏蔽不完全导致的偶极–偶极作用）
2. 描述长程力常数需要哪些量？（介电张量与 Born 有效电荷）

---

### 示例 8：含长程偶极–偶极处理的声子色散与 DOS

**教学目标**

- 利用静态介电张量与 Born 有效电荷，计算极性材料力常数的长程部分，得到含 LO-TO 分裂的声子色散。

**理论背景**

- 力常数的分解对应动力学矩阵的分解：

$$D_{I\alpha J\beta}(\mathbf{q}) = D^\text{SR}_{I\alpha J\beta}(\mathbf{q}) + D^\text{LR}_{I\alpha J\beta}(\mathbf{q})$$

- 长程部分（倒空间表达式，含 Born 有效电荷 $\mathbf{Z}^*$ 与高频介电张量 $\epsilon^\infty$，用 Ewald 参数 $\lambda$ 平滑收敛）：

$$D^\text{LR}_{i\alpha j\beta}(\mathbf{q}) = \frac{4\pi e^2}{\Omega_0} \sum_\mathbf{G} \sum_{l} \frac{\big[(\mathbf{G}+\mathbf{q})\cdot\mathbf{Z}^*_{i\alpha}\big]\big[(\mathbf{G}+\mathbf{q})\cdot\mathbf{Z}^*_{j\beta}\big]}{(\mathbf{G}+\mathbf{q})\cdot\epsilon^\infty\cdot(\mathbf{G}+\mathbf{q})} e^{i(\mathbf{q}+\mathbf{G})\cdot(\mathbf{l}+\mathbf{r}_i-\mathbf{r}_j)} e^{-\frac{(\mathbf{G}+\mathbf{q})\cdot\epsilon^\infty\cdot(\mathbf{G}+\mathbf{q})}{4\lambda^2}}$$

- 实用流程：
  1. 用有限尺寸超胞计算 $\Phi_{I\alpha J\beta}$；
  2. 减去解析长程部分得到 $\Phi^\text{SR}$；
  3. 用 $\Phi^\text{SR}$ 傅里叶插值得到 $D^\text{SR}(\mathbf{q})$；
  4. 在原胞中计算 $D^\text{LR}(\mathbf{q})$ 并相加。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/phonon/e08_*`）

- VASP 通过 `LPHON_POLAR = .TRUE.` 自动完成上述处理，并用 `PHON_DIELECTRIC` 提供介电张量、`PHON_BORN_CHARGES` 提供 Born 有效电荷；
- 把示例 5 的 `vaspout.h5` 复制为 `vaspin.h5` 后运行。

**结果分析**

- 对比示例 6 的色散图：加入长程处理后，极性材料在 $\Gamma$ 点附近出现 **LO-TO 分裂**（纵向光学支与横向光学支分开），且色散插值更加平滑；
- 同样在不同 q 网格上计算 DOS 并检查最密网格结果；
- 加分任务：用 3x3x3、4x4x4 超胞力常数重跑第 6、8 节，验证长程处理如何让插值更平滑。

**思考题**：加入长程力常数特殊处理后，声子色散最主要的变化是什么？（LO-TO 分裂 + 更平滑插值）

---

### 功能与标签汇总（Phonon Part 2）

| 功能 | 关键标签 / 文件 | 说明 |
|---|---|---|
| 体积弛豫 | `IBRION` + `ISIF` | 晶格常数对 k 网格收敛后使用 |
| 有限差分歧变 | `IBRION = 6` | 体材料三方向等比扩胞 |
| 读取力常数插值 | `LPHON_READ_FORCE_CONSTANTS` + `vaspin.h5` + QPOINTS | 计算色散/DOS 后退出 |
| 声子 DOS | QPOINTS q 网格 | 插值快，需对 q 网格收敛 |
| 介电张量 + Born 有效电荷 | `LEPSILON = .TRUE.` | 极性材料长程力常数所需 |
| LO-TO 分裂 | `LPHON_POLAR = .TRUE.` + `PHON_DIELECTRIC` + `PHON_BORN_CHARGES` | 自动分解短程/长程力常数并插值 |

---

## 16. 电子–声子相互作用（Electron-Phonon Interactions）

许多体系中电子与振动（声子）自由度可分别处理（电子远快于核运动）。该处理是近似的，可通过电子–声子耦合修正。电子–声子散射在诸多应用中占主导，如半导体迁移率、室温金属电导率。

## 16.1 Part 1：微扰理论计算带隙重正化（Bandgap renormalization）

- 来源：<https://vasp.at/tutorials/latest/electron-phonon/part1/>
- 教程输入文件下载：[electron-phonon-part1.zip](https://vasp.at/tutorials/latest/electron-phonon-part1.zip)

电子–声子相互作用计算至少分两步：第一步在**超胞**中用有限差分计算电子–声子势，第二步在**原胞**中计算电子–声子矩阵元 $g_{mn\nu}(\mathbf{k},\mathbf{q})$ 与带隙重正化。该流程把昂贵的超胞计算与原胞中的性质计算分离，兼顾效率与精度。

---

### 示例 1：金刚石的温度依赖带隙重正化

**教学目标**

- 从超胞计算电子–声子势；
- 计算声子诱导的带隙重正化；
- 比较零点与有限温度下的带隙。

**第一步：超胞计算电子–声子势**（输入位于 `$TUTORIALS/electron_phonon/e01_diamond/supercell`）

- POSCAR 为金刚石原胞的 2x2x2 重复（教程因算力限制用小超胞；实际应对目标性质做 2x2x2→3x3x3→4x4x4 的超胞收敛测试）；
- INCAR 两个关键标签：
  - `IBRION`：激活有限差分驱动器；
  - `ELPH_POT_GENERATE`：计算电子–声子势；
- VASP 先做超胞电子基态，然后逐原子沿三个笛卡尔方向位移、计算 Kohn–Sham 势的变化，结果写入 `phelel_params.hdf5`；KPOINTS 只有一个 k 点，可用 gamma 版 VASP；
- `ELPH_POT_GENERATE` 还会生成 `CONTCAR_ELPH`：与 CONTCAR 类似但包含**原胞**结构，应作为下一步的 POSCAR；不满意自动选的原胞可用 `ELPH_POT_LATTICE` 自定义；
- 电子–声子势在超胞精细 FFT 网格上计算，按 `ELPH_POT_FFT_MESH` 定义的原胞 FFT 网格输出；该网格默认必须与原胞电子–声子计算所用 FFT 网格完全一致。网格尺寸由 `ENCUT`、`PREC` 决定（也可直接设 `NGX/NGY/NGZ`）：两步保持相同的 `ENCUT` 与 `PREC`（推荐）则网格自动匹配。

**第二步：原胞计算带隙重正化**（输入位于 `unit cell` 目录）

- 把 `phelel_params.hdf5` 与 `CONTCAR_ELPH`（改名为 POSCAR）复制到原胞目录；
- INCAR 与超胞计算基本相同（`ENCUT`、`PREC` 完全一致），设 `ELPH_MODE = renorm`（元标签，自动设置合理默认值）：
  - `ELPH_RUN`：启动电子–声子计算；
  - `ELPH_SELFEN_GAPS`：自动计算带隙重正化；
  - `ELPH_SELFEN_FAN` 与 `ELPH_SELFEN_DW`：自能中包含 Fan–Migdal（FM）与 Debye–Waller（DW）两项贡献；
- `ELPH_SELFEN_TEMPS`：以 K 为单位列出温度列表，单次运行计算所有温度下的重正化；
- KPOINTS 只含少量 k 点用于自洽电子最小化；电子–声子物理量（矩阵元、电子能量、声子频率）在更密的 `KPOINTS_ELPH` 网格上计算。

**结果分析**

- 输出三处：标准输出（最简摘要，含直接/间接带隙的温度依赖重正化）、OUTCAR 末尾（详细）、`phelel_params.hdf5`（机器可读，供 py4vasp 等后处理）；
- 零温下由于量子涨落仍存在晶格振动（零点运动），产生**零点重正化（ZPR）**：即便基态计算也需要电子–声子物理才能得到正确带隙；
- 重正化符号为负：准粒子（重正化）带隙小于 Kohn–Sham 带隙，ZPR 通常使带隙收缩；
- 温度升高 → 可参与散射的声子态增多 → 电子态被额外散射过程重正化 → 带隙进一步收缩。用 py4vasp 可绘制直接带隙随温度的变化曲线。

**思考题**

1. 为什么零温也有重正化？（零点运动）
2. `CONTCAR_ELPH` 包含什么？（原胞结构）

---

### 示例 2：金刚石零点重正化的收敛研究

**教学目标**

- 对 q 点数 $N_q$、中间态数 $N_b$、展宽参数 $\delta$ 做收敛研究；
- 理解各收敛参数如何影响自能。

**理论背景**

- 只看 ZPR 即可：温度只通过费米/玻色分布函数进入自能，ZPR 收敛通常足以保证有限温度精度；
- $T=0$ 的 Fan–Migdal 自能：

$$\Sigma^{\text{FM}}_{n\mathbf{k}}(T=0) = \frac{1}{N_q} \sum_{q}^{N_q} \sum_{m}^{N_b} \sum_{\nu} \frac{|g_{mn\nu}(\mathbf{k},\mathbf{q})|^2}{\varepsilon_{n\mathbf{k}} - \varepsilon_{m\mathbf{k}+\mathbf{q}} \pm \omega_{\nu\mathbf{q}} + i\delta}$$

**三个收敛参数**

1. **q 点网格 $N_q$**（`KPOINTS_ELPH`）：在 4x4x4、6x6x6、8x8x8、10x10x10 网格上计算 ZPR。本例基本带隙在 10x10x10 即收敛，但直接带隙尚未收敛——实际中常需更高密度，每种材料、每个带隙都不同，务必复查；
2. **中间态数 $N_b$**（`ELPH_NBANDS_SUM`，可给列表、单次运行全部评估）：态求和收敛慢；本例 200 条带仍勉强收敛，且曲线非单调、难以判断；视体系与赝势可能需要数百至上千条中间带。`ELPH_MODE = renorm` 默认设 `ELPH_NBANDS = -2`（用满所有可用能带，数量=平面波系数数）。要进一步增多只能提高 `ENCUT`——但 FFT 网格会变，**必须用更高 ENCUT 重跑超胞计算**。因此推荐超胞计算一开始就用较高截断能，避免反复重跑昂贵的超胞步；
3. **展宽参数 $\delta$**（`ELPH_SELFEN_DELTA`，可给列表）：来自多体微扰理论的无穷小虚移，数值上充当处理（近）简并态的展宽。应尽量小（大值引入非物理贡献），太小则浮点噪声显现；默认 0.01 eV 通常是好的起点。本例中 1 meV–100 meV 范围内结果无显著变化。
   - 展宽应与 q 网格配套：更密 q 网格引入更多近能量态、自能分母更小。推荐 **q 网格 x δ 双收敛研究**（如在 q-points 目录的 INCAR 加 `ELPH_SELFEN_DELTA = 0.0001 0.001 0.01 0.1 1`）；
   - 对**极性材料**尤其重要：长程静电作用使纵向光学（LO）模式对应的电子–声子矩阵元在长波极限发散（[J. Chem. Phys. 143 (2015) 102813](https://doi.org/10.1063/1.4927081)），所需 q 网格可能非常密。

**思考题**

1. 什么是 q 点？
2. 电子–声子计算能带不够收敛怎么办？（提高 ENCUT 并重跑超胞）
3. 展宽参数为何影响 ZPR？

---

### 示例 3：极性材料 MgO 的带隙重正化

**教学目标**

- 计算 Born 有效电荷与介电张量；
- 计算电子–声子势；
- 计算声子诱导带隙重正化；
- 用外推法确定零点重正化。

**理论背景**

- 极性材料中长程静电作用使某些声子模式的电子–声子矩阵元发散（类似声子的 LO-TO 分裂），需要密得多的 q 网格与谨慎的收敛测试；与声子类似，存在把长程相互作用并入矩阵元的特殊修正方案；
- 长程（仅偶极）电子–声子矩阵元：

$$g_{mn\nu}^{\text{LR}}(\mathbf{k},\mathbf{q}) = i \frac{4\pi}{V_{\text{pc}}} \sum_{\kappa} \sqrt{\frac{\hbar}{2 m_\kappa \omega_{\nu\mathbf{q}}}} \sum_{\mathbf{G}} \frac{(\mathbf{q}+\mathbf{G})\cdot\mathbf{Z}^\star_\kappa\cdot\mathbf{e}_{\kappa,\nu\mathbf{q}}}{(\mathbf{q}+\mathbf{G})\cdot\boldsymbol{\epsilon}_\infty\cdot(\mathbf{q}+\mathbf{G})} \langle \psi_{m\mathbf{k}+\mathbf{q}} | e^{i(\mathbf{q}+\mathbf{G})\cdot(\mathbf{r}+\boldsymbol{\tau}_\kappa)} | \psi_{n\mathbf{k}} \rangle$$

- 声子频率/本征矢与电子 Bloch 轨道可由 VASP 在电子–声子计算中得出；Born 有效电荷与介电张量需事先计算并作为输入。
- "NAC"（non-analytic-term correction）术语指布里渊区中心动力学矩阵非解析行为的处理（[Phys. Rev. B 55 (1997) 10355](https://doi.org/10.1103/PhysRevB.55.10355)），电子–声子中有类似概念，故沿用该词。

**三步流程**

1. **介电性质**（`e03_MgO/nac`）：INCAR 设 `LEPSILON = .TRUE.`，用密度泛函微扰理论计算电场响应，得到 Born 有效电荷张量与离子钳制静态介电张量；结果可在 OUTCAR 末尾或用 py4vasp 从 `vaspout.h5` 读取（`.print()` 查看、`.read()` 用于脚本）；
2. **超胞电子–声子势**（`e03_MgO/supercell`）：与金刚石类似，但 INCAR 中必须填入上一步得到的 `PHON_DIELECTRIC` 与 `PHON_BORN_CHARGES`（教程提供脚本自动填充）；
3. **ZPR 与 q 点收敛**（`el-ph` 目录多个 q 网格子目录）：把 `phelel_params.hdf5` 与 `CONTCAR_ELPH` 复制到各子目录；因介电张量已在超胞步设置，此步 INCAR 与金刚石例子相同。MgO 只有直接带隙，输出中 "direct" 与 "fundamental" 相同。

**结果分析**

- 关注量是 ZPR，输出记为 `KS-QP gap (meV)`（准粒子带隙减 Kohn–Sham 带隙）；
- ZPR 随 q 网格明显不收敛，而且**即使继续加密 q 网格也不收敛**：LO 声子模式的电子–声子矩阵元在长波极限（$\mathbf{q}\to\mathbf{0}$）发散。好在发散可积，（非绝热）自能仍有限、ZPR 良定义；
- 由长程矩阵元公式：$\mathbf{q}\to\mathbf{0}$ 时分母行为 $g^{\text{LR}} \propto 1/|\mathbf{q}|$，自能积分 $\int d^3q\, 1/|\mathbf{q}|^2$ 在原点附近可积，因此可通过逐步加密 q 网格并**外推**到 $N_q \to \infty$ 得到 ZPR；
- 外推做法：把 ZPR 对**网格间距的倒数**作图，取线性区外推。注意：教程仅在线性区演示外推思想；实际需密得多（如 32x32x32 甚至更高）的 q 网格才能进入可靠的线性区。

**思考题**

1. 为什么要计算 Born 有效电荷与介电张量？（描述极性材料电子–声子相互作用的长程部分）
2. 为什么仅靠密 q 网格不足以收敛极性材料的 ZPR？（LO 模式矩阵元长波发散，需外推）

---

### 功能与标签汇总（Electron-phonon Part 1）

| 功能 | 关键标签 / 文件 | 说明 |
|---|---|---|
| 电子–声子势生成 | `IBRION` + `ELPH_POT_GENERATE` | 超胞有限差分，写 `phelel_params.hdf5` 与 `CONTCAR_ELPH` |
| 原胞选择 | `ELPH_POT_LATTICE` | 自定义 `CONTCAR_ELPH` 的原胞 |
| FFT 网格匹配 | `ELPH_POT_FFT_MESH`、`ENCUT`/`PREC`、`NGX/NGY/NGZ` | 两步须一致，否则重跑超胞 |
| 带隙重正化 | `ELPH_MODE = renorm`、`ELPH_RUN`、`ELPH_SELFEN_GAPS` | 元标签自动设默认 |
| 自能贡献 | `ELPH_SELFEN_FAN`、`ELPH_SELFEN_DW` | Fan–Migdal 与 Debye–Waller |
| 温度列表 | `ELPH_SELFEN_TEMPS` | 单次运行多温度 |
| 电子–声子 k/q 网格 | `KPOINTS_ELPH` | 密于自洽 KPOINTS |
| 中间态收敛 | `ELPH_NBANDS_SUM`、`ELPH_NBANDS = -2` | 列表单次评估；带数受平面波数限制 |
| 展宽收敛 | `ELPH_SELFEN_DELTA` | 与 q 网格配套做双收敛 |
| 极性修正 | `LEPSILON` + `PHON_DIELECTRIC` + `PHON_BORN_CHARGES` | Born 电荷/介电张量供长程偶极修正 |
| ZPR 外推 | ZPR 对网格间距倒数作图 | 极性材料必须外推到 $N_q\to\infty$ |

---

## 16.2 Part 2：统计采样法计算电子–声子相互作用

- 来源：<https://vasp.at/tutorials/latest/electron-phonon/part2/>
- 教程输入文件下载：[electron-phonon-part2.zip](https://vasp.at/tutorials/latest/electron-phonon-part2.zip)

与 Part 1 的微扰理论路线不同，本部分用**统计采样**方法：温度 $T$ 下的可观测量 $O(T)$ 是对不同坐标组 $x_T^{\text{MC},i}$ 采样后的统计平均。完整 Monte Carlo（MC）采样中每组位移由平衡位置加随机位移得到：

$$x_T^{\text{MC},i} = x_{\text{eq}} + \Delta\tau^{\text{MC},i}$$

位移由笛卡尔坐标、原子质量 $M_\kappa$、原子 $\kappa$ 上的声子本征模 $\nu$、本征频率 $\omega_\nu$ 与正态分布随机数决定（声子本征模与频率在谐波近似下计算）。MC 需要大量结构、计算昂贵。Zacharias 与 Giustino 提出 **one-shot（OS）方法**（ZG 构型）：只用一组位移（[Phys. Rev. B 94 (2016) 075125](https://doi.org/10.1103/PhysRevB.94.075125)）：

$$\Delta\tau^{\text{OS}} = \sqrt{\frac{1}{M_\kappa}} \sum_{\nu}^{3(N-1)} (-1)^{\nu-1} \varepsilon_{\kappa,\nu}\, \sigma_{\nu,T}$$

本征模按本征频率升序求和，位移幅度为

$$\sigma_{\nu,T} = \sqrt{(2n_{\nu,T}+1)\frac{\hbar}{2\omega_\nu}}, \qquad n_{\nu,T} = \left[\exp\left(\frac{\hbar\omega_\nu}{k_B T}\right) - 1\right]^{-1}$$

其中 $n_{\nu,T}$ 为 Bose–Einstein 占据数。one-shot 方法把 MC 求和约化为单次计算（[New J. Phys. 20 (2018) 123008](https://doi.org/10.1088/1367-2630/aaf53f)）。

---

### 示例 4：金刚石的 one-shot 计算

**教学目标**

- 用 one-shot 方法创建（含特殊构型的）结构；
- 解读电子–声子相互作用对带隙的影响。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/electron-phonon/e04_SIMPLE_ONE-SHOT_C`）

- 两套 INCAR：PBE 基态计算（`INCAR.dft`）与 one-shot 构型生成（`INCAR.os_disp`）；
- 关键标签：
  - `IBRION = 6`：获得动力学矩阵本征矢；
  - `PHON_LMC = .TRUE.`：启用 one-shot 与 Monte Carlo 计算；
  - `PHON_NSTRUCT > 0`：MC 采样的结构数；`PHON_NSTRUCT = 0` 切换为 one-shot；
  - `PHON_TLIST`：温度列表——一次计算可同时生成多个温度的 one-shot POSCAR（动力学矩阵只需算一次）；
  - `PHON_NWRITE`、`PHON_DOS`：绘制声子态密度；
- QPOINTS 仅声子 DOS 计算需要，one-shot/随机采样的电子–声子计算不需要。

**结果分析**

- 计算生成 0 K 的位移结构 `POSCAR.T=0.`（只有原子位置移动）；声子 DOS 峰对应本征模能量，正是 one-shot 位移所用；
- 分别对原胞与 one-shot 结构做标准 DFT：
  - 完整超胞：带隙附近因结构对称性本征能级强简并，带隙约 4.54 eV；
  - one-shot 结构（占据数按 Bose–Einstein 分布定义）：简并全部解除（对称性破缺），带隙约 3.57 eV，比完整结构低近 1 eV；
- 文献中常不用最高价带与最低导带之差作带隙，而是取原先简并的价带（band 30–32）与导带（band 33–38）的**平均 band energy** 之差，以便与实验更好比较。

**思考题**

1. one-shot 方法为何要求解动力学矩阵？
2. 占据数用什么分布计算？（Bose–Einstein）
3. 电子–声子相互作用对带隙附近本征态有什么影响？（解除简并、带隙收缩）
4. QPOINTS 文件为何使用？one-shot 计算是否必需？

---

### 示例 5：one-shot 方法的超胞收敛

**教学目标**

- 学会做简单的收敛测试；
- 体会超胞尺寸对 one-shot 与随机采样方法的重要性。

**要点**

- 声子 q 网格必须大到覆盖胞内所有声子；one-shot/MC 方法用**超胞尺寸**控制 q 网格大小。大超胞计算量增长快，因此做收敛测试、取最小够用尺寸；
- 对 2x2x2、3x3x3、4x4x4 超胞各做三步计算：完整超胞带隙 → 生成 one-shot POSCAR → 在位移结构上算重正化带隙；重正化 = 两者之差；
- **收敛设置的注意事项（重要）**：测试对象应是目标性质（带隙、力、声子模式）。求导会放大误差：力（一阶导）、动力学矩阵（二阶导）及其依赖量（声子模式）误差层层累积。one-shot 位移结构 `POSCAR.NxNxN_T0` 对初始误差敏感——松散设置下核数/机器差异可导致约 10 mÅ 的结构差异，进而使 ZPR 差约 0.05 eV。因此只收敛带隙不够，**还必须收敛声子模式**；改进方向：提高 `ENCUT`、`PREC`，把 `ALGO` 从 Fast 改为 Normal，收紧 `EDIFF`。小 k 网格+低截断会引入大 Pulay 应力（基组不完备→力不收敛→声子模式不收敛），可用 `grep "external pressure" OUTCAR.*` 检查；
- 用平均简并价带/导带能量之差计算带隙（文献标准做法，利于与实验收敛）；
- 可视化位移：红色为原超胞、青色为位移后超胞：

![one-shot 方法位移前后 2x2x2 超胞对比](/imgs/2026-08-05/fig50.png)

- 用 `grep 'fundamental gap' OUTCAR.xxx` 提取带隙；脚本输出存入 `gap_vs_supercell_size.dat` 并与参考数据对比作图；
- 带隙重正化随超胞尺寸收敛缓慢且不单调；用更严格设置后金刚石间接带隙可收敛（见 New J. Phys. 20 (2018) 123008 表 2、表 3）。

**思考题**

1. 为什么要做收敛测试？
2. 得到带隙重正化需要几次计算？（完整超胞 + one-shot 结构各一次，每个超胞尺寸共两次）

---

### 示例 6：金刚石带隙的温度依赖

**教学目标**

- 计算温度依赖的带隙重正化；
- 把温度依赖拟合到合适的经验曲线。

**要点**

- 以 4x4x4 超胞为起点，用 `PHON_TLIST` 在单次计算中生成多个温度下的 one-shot POSCAR（`POSCAR.T=xxx`），对每个结构做基态计算得带隙（0–700 K 范围，教程为省时注释掉了部分温度）；
- 半导体/绝缘体带隙温度依赖常用 **modified-Varshni 关系**（[Physica 34 (1967) 149](https://doi.org/10.1016/0031-8914(67)90062-6)）：

$$E_{\text{GAP}}(T) = E_{\text{GAP}}(0) - \frac{A T^{4}}{(B+T)^2}$$

其中 $A$、$B$ 为材料拟合参数。拟合结果显示计算带隙能被 Varshni 关系很好描述；
- **缺失的效应**：体积膨胀（热膨胀）贡献未包含，可通过准谐波计算（quasi-harmonic）获得。

**思考题**

1. 半导体/绝缘体带隙温度依赖用什么关系描述？（modified-Varshni）
2. 完整描述带隙温度依赖还需要什么？（体积膨胀/准谐波效应）

---

### 示例 7：Monte Carlo 采样 vs one-shot

**教学目标**

- 用 Monte Carlo 采样计算带隙；
- 与 one-shot 方法比较。

**要点**

- 用谐波分布的密度矩阵生成 30 个随机结构（[Phys. Status Solidi A 188 (2001) 1209](https://doi.org/10.1002/1521-396X(200112)188:4<1209::AID-PSSA1209>3.0.CO;2-2)），逐个做基态计算，MC 带隙取算术平均（[New J. Phys. 20 (2018) 123008](https://doi.org/10.1088/1367-2630/aaf53f)）；
- 输入位于 `$TUTORIALS/electron-phonon/e07_C_Monte_Carlo_sampling`，三套 INCAR：基态（`INCAR.dft`）、one-shot（`INCAR.elphon_oneshot`）、MC（`INCAR.elphon_mc`）；
  - `PHON_NSTRUCT` 设定 MC 结构数；与 one-shot 不同，MC 的温度由 `TEBEG` 选择；
- 流程：先生成 0 K one-shot 结构并计算（结果存为参考 `OUTCAR`），再生成 30 个编号结构 `POSCAR.T=0.X`，脚本批量计算后统计平均；
- 结果：20 个结构时 MC 已相对 30 结构结果收敛到第二位小数；one-shot 带隙与 MC 平均带隙非常接近，且超胞越大两者越趋于一致；
- **实践建议**：超胞足够大时两者结果相近，推荐**始终使用 one-shot**；MC 采样仅用于特殊场合，例如需要多个有限温度谐波采样结构作为分子动力学起点。

**思考题**

1. MC 采样得到含电子–声子作用带隙的计算步骤是什么？
2. one-shot 与 MC 结果何时互相收敛？（超胞足够大时）

---

### 功能与标签汇总（Electron-phonon Part 2）

| 功能 | 关键标签 / 文件 | 说明 |
|---|---|---|
| 动力学矩阵 | `IBRION = 6` | one-shot/MC 位移的基础 |
| 统计采样开关 | `PHON_LMC = .TRUE.` | 启用 one-shot 与 MC |
| 采样结构数 | `PHON_NSTRUCT` | 0 = one-shot；>0 = MC 结构数 |
| 温度列表 | `PHON_TLIST` | one-shot 多温度 POSCAR 一次生成 |
| MC 温度 | `TEBEG` | MC 采样用 TEBEG 指定温度 |
| 声子 DOS | `PHON_DOS`、`PHON_NWRITE`、QPOINTS | 仅 DOS 需要 QPOINTS |
| 带隙提取 | `grep 'fundamental gap'` | 简并带平均能量差更贴近实验 |
| 温度拟合 | modified-Varshni | 缺热膨胀项，需准谐波补充 |

---

## 16.3 Part 3：电子–声子矩阵元与 VASP+phelel 工作流

- 来源：<https://vasp.at/tutorials/latest/electron-phonon/part3/>
- 教程输入文件下载：[electron-phonon-part3.zip](https://vasp.at/tutorials/latest/electron-phonon-part3.zip)

---

### 示例 8：获取电子–声子矩阵元（ELPH_DRIVER = mels）

**教学目标**

- 计算电子–声子矩阵元；
- 沿 q 路径绘制矩阵元；
- 把结果图与电子、声子能带结构联系起来。

**理论背景**

- 电子与晶格振动的相互作用修正电子能量，由电子自能描述。二阶微扰理论把自能分为 Fan–Migdal（FM）与 Debye–Waller（DW）两项，它们都通过电子–声子矩阵元计算：

$$g_{mn\mathbf{k},\nu\mathbf{q}} = \left\langle \psi_{m\mathbf{k}+\mathbf{q}} \left| \partial_{\nu\mathbf{q}} V_{\text{KS}} \right| \psi_{n\mathbf{k}} \right\rangle$$

其中 $(\nu,\mathbf{q})$ 是声子模式指标与波矢，$m,n$ 是电子能带；
- 计算自能时矩阵元是即算即弃的（高效）；但若想单独分析或绘制矩阵元，需要显式输出。

**任务与输入要点**（示例文件内容略，输入位于 `$TUTORIALS/electron_phonon/e08_matrix_elements`）

- 用 4x4x4 k 网格，沿 $\Gamma$-X-U-K-$\Gamma$-L-W-X q 路径绘制金刚石（原胞）的平均电子–声子矩阵元；
- INCAR 电子–声子相关标签：
  - `ELPH_RUN`：启动电子–声子计算；
  - `ELPH_DRIVER = mels`：不计算自能，而是累积所有满足选择定则的矩阵元并写入 `vaspelph.h5` HDF5 文件；
  - `ELPH_SELFEN_IKPT`：指定计算自能的 k 点（本例仅 $\Gamma$）；
  - `ELPH_SELFEN_BAND_START` / `ELPH_SELFEN_BAND_STOP`：自能计算所用电子能带起止（本例 1 2 3 4）；
  - `ELPH_NBANDS`：电子–声子驱动器在密 k 网格上计算的能带数；
- 教程提供 `phelel_params.hdf5`（跳过超胞步）与绘图脚本 `elph_mels.py`。

**结果分析**

- 绘图脚本生成四联图：电子能带（1）、声子能带（2）、$\varepsilon_{m,\mathbf{k}+\mathbf{q}}$（$\mathbf{k}=\Gamma$）+ $\omega_{\nu,\mathbf{q}}$ + 矩阵元沿 q 路径（3）、按声子支分解的矩阵元（4）；
- 物理图像：$\Gamma$ 点电子态 $n,\mathbf{k}$ 与 $\mathbf{k}+\mathbf{q}$ 处电子态 $m$ 通过发射/吸收一支波矢 $\mathbf{q}$、支指标 $\nu$ 的声子发生耦合；
- 第三张图绘制的是 $m=n=4$ 所有简并态（含简并的 $n=2,3,4$）的平均矩阵元（平均方式见 `elph_mels.py` 或 [Phys. Rev. B 100, 174304 (2019)](https://doi.org/10.1103/PhysRevB.100.174304) 式 97）；
- 沿声子能带追踪 $\nu=6$：若按**带指标**追踪，声子支交叉处导数突变 → 矩阵元曲线不连续；若按**物理态**追踪，则在交叉处切换 $\nu$，曲线连续。不连续并非物理的。

**思考题**

1. 如何沿声子能带与矩阵元追踪一个声子态？
2. 按声子带指标追踪为何出现不连续？是否物理？

---

### 示例 9：金刚石电子能带（velph 入门）

**教学目标**

- 用 `velph` 生成 `velph.toml` 配置文件，准备能带与 DOS 计算目录；
- 认识 phelel/velph 自动化工作流。

**背景**

- `phelel` 是 VASP 之外的电子–声子工具，配套命令行工具 `velph` 可自动完成：
  - NAC 所需的 Born 有效电荷与介电张量计算设置；
  - 超胞生成与电子–声子势计算；
  - 原胞电子–声子计算准备；
  - 文件管理与各步设置一致性维护。
- 相比 VASP-only 流程，`phelel` 在原子位移选择上更灵活、构建电子–声子势时使用更多对称操作（可缩短运行时间）；
- 本例演示 C 金刚石的完整 VASP+phelel 工作流，起点只有原胞 POSCAR 与 POTCAR。

**工作流（输入位于 `$TUTORIALS/electron_phonon/e09_phelel_band`）**

![velph 工作流示意](/imgs/2026-08-05/fig51.png)

1. **初始化**：`velph init` 需要 POSCAR、项目目录与超胞尺寸（`--max-num-atoms`）；`--kspacing=0.6` 对应 INCAR 的 `KSPACING`（倒空间 k 点间距，自动生成 KPOINTS；0.6 给出与前面 MgO 例子相同的 2x2x2 超胞采样），生成 `velph.toml`；
   - `velph.toml` 含全部后续设置，大多与 INCAR 标签一一对应：`[vasp.*.incar]`、`[vasp.*.kpoints]` 小节，文件末尾还有 `[symmetry]`、`[unitcell.*]`、`[lattice]`；
   - `template.toml` 可提供 velph 默认之外的输入设置（如 `ENCUT`）；
2. **生成并运行能带计算**：`velph el_bands generate` 生成 `el_bands/bands` 与 `el_bands/dos` 两套文件（能带设置在 `[vasp.el_bands.incar]`），运行 VASP 后用 py4vasp 查看能带；
3. **DOS 与并排绘图**：运行 `dos` 目录计算，用 velph 把能带与 DOS 并排绘制：

![金刚石电子能带与态密度](/imgs/2026-08-05/fig52.png)

**思考题**

1. 使用 velph 时在哪里找 INCAR 设置？（`velph.toml` 的 `[vasp.*.incar]` 小节）
2. 这种结构化配置的优点是什么？

---

### 示例 10：极性材料 MgO 的介电性质与带隙重正化（velph 全流程）

**教学目标**

- 用 velph 计算极性材料的介电性质；
- 用集成工作流计算带隙重正化；
- 对比手动与自动化流程的结果。

**任务**

- 用 velph 生成 `velph.toml`（kspacing=0.6），计算 Born 有效电荷与介电张量，自动生成 2x2x2 超胞、计算电子–声子势，然后在原胞做电子–声子计算，得到 MgO 在 0–700 K 的带隙重正化。

**工作流（输入位于 `$TUTORIALS/electron_phonon/e10_phelel_MgO`）**

1. `velph init` 生成 `velph.toml`；
2. **NAC/介电性质**：查看 `vasp.nac` 小节（`vasp.nac.incar` 含收敛标签与 `LEPSILON = True`，`vasp.nac.kpoints` 定义 k 网格）。因教程初始化用了粗网格（4x4x4），需在 `velph.toml` 中手动改为 6x6x6 以与 Part 1 一致；`velph nac generate` 生成 `nac` 目录，运行 VASP 得到 Born 有效电荷与介电张量——与 Part 1 示例 3 手动结果**完全一致**，唯一区别是工作流；velph 随后自动知道去何处找这些张量并带入后续步骤；
3. **超胞电子–声子势**：`velph phelel init` 找出所需原子位移、生成 `phelel` 目录与 `phelel_disp.yml`（指明位移哪些原子、应用哪些对称性）；`velph phelel generate` 为每个位移（加未微扰结构）生成子目录与输入文件；每个子目录 INCAR 含新标签 `ELPH_PREPARE`：把局域势等必要信息写入 `vaspout.h5`；逐个子目录运行 VASP 后，回到根目录执行 `velph phelel differentiate`：自动调用 phelel 做有限差分、生成 `phelel_params.hdf5`。VASP+phelel 在位移选择上比 VASP-only 更灵活；
4. **原胞自能计算**：编辑 `velph.toml` 的 `[vasp.selfenergy.*]` 小节——`[vasp.selfenergy.kpoints_dense]` 改密网格为 8x8x8，`[vasp.selfenergy.incar]` 把 `ELPH_SELFEN_TEMPS` 设为 0–700 K、删除 `ELPH_NBANDS_SUM` 并设 `ELPH_NBANDS = -2`（用满所有能带）；`velph selfenergy generate` 生成 `selfenergy` 目录（核对 INCAR 与 KPOINTS_ELPH），运行 VASP；
5. **结果**：用 py4vasp 查看带隙重正化随温度的变化——0 K 点与手动流程一致，ZPR 随温度升高而减小（绝对值增大、带隙进一步收缩）。

**VASP-only vs VASP+phelel**：前者适合熟悉 VASP 计算的用户；后者适合熟悉 phonopy 风格工作流的用户，目前灵活性更高；选择取决于个人偏好。

**思考题**

1. VASP+phelel 工作流在哪些方面辅助计算？
2. 何时选 VASP-only、何时选 VASP+phelel？

---

### 功能与标签汇总（Electron-phonon Part 3）

| 功能 | 关键标签 / 工具 | 说明 |
|---|---|---|
| 矩阵元输出 | `ELPH_DRIVER = mels` | 矩阵元写入 `vaspelph.h5` |
| 自能 k 点/能带选择 | `ELPH_SELFEN_IKPT`、`ELPH_SELFEN_BAND_START/STOP` | 按需选择 |
| 密网格能带数 | `ELPH_NBANDS` | -2 = 全部平面波系数对应能带 |
| 工作流管理 | `velph init/generate/differentiate` + `velph.toml` | 单配置文件驱动全流程 |
| 超胞位移 | `phelel_disp.yml`、`ELPH_PREPARE` | 每位移一个子目录，局域势写入 `vaspout.h5` |
| NAC 集成 | `velph nac` | 自动衔接 Born 电荷/介电张量到后续步骤 |

---

## 16.4 Part 4：铁的电导率（Conductivity of iron）

- 来源：<https://vasp.at/tutorials/latest/electron-phonon/part4/>
- 教程输入文件下载：[electron-phonon-part4.zip](https://vasp.at/tutorials/latest/electron-phonon-part4.zip)

金属中电输运源于载流子（电子与空穴）在外电场下的运动。电导率 $\sigma$ 量化电流通过材料的难易程度，是电阻率 $\rho$ 的倒数 $\sigma = 1/\rho$。由线性化 Boltzmann 方程，输运分布函数为：

$$\mathcal{T}(\epsilon) = \frac{e^2}{N_\mathbf{k}\Omega} \sum_{n\mathbf{k}} \tau_{n\mathbf{k}}\, \mathbf{v}_{n\mathbf{k}} \otimes \mathbf{v}_{n\mathbf{k}}\, \delta(\epsilon_{n\mathbf{k}}-\epsilon)$$

其中 $\tau_{n\mathbf{k}}$ 是态 $(n,\mathbf{k})$ 的弛豫时间，$\mathbf{v}_{n\mathbf{k}}$ 为群速度。Onsager 输运系数（[J.M. Ziman, Electrons and Phonons (2001)](https://global.oup.com/academic/product/electrons-and-phonons-9780198507796)）：

$$\mathcal{L}_{ij} = \int d\epsilon\, \mathcal{T}(\epsilon)\, (\epsilon-\mu)^{i+j-2} \left(-\frac{\partial f^0}{\partial \epsilon}\right)$$

由此得电导率 $\sigma = \mathcal{L}_{11}$、Seebeck 系数 $S = \mathcal{L}_{11}/\mathcal{L}_{12}$、电子热导 $\kappa_e = \frac{1}{T}(\mathcal{L}_{22} - \mathcal{L}_{21}\mathcal{L}_{11}^{-1}\mathcal{L}_{12})$（晶格热导另行处理，此处省略）。

---

### 示例 11：常数弛豫时间近似（CRTA）

**教学目标**

- 在常数弛豫时间近似（CRTA）下计算简单金属的电导率并与实验对比；
- 测试电导率对 k 网格的收敛；
- 考虑铁的磁性做自旋极化计算，并与非磁结果比较。

**CRTA 基本计算**（输入位于 `$TUTORIALS/electron-phonon/e11_conductivity_iron_crta/crta`）

- CRTA 假设费米面上所有电子共享同一平均散射寿命 $\tau$：$\tau$ 成为费米面电子速度平均值的整体缩放因子，让 CRTA 能快速展示能带几何如何单独决定电导率。金属费米面电子本就可自由移动（无带隙），CRTA 是合理起点；
- INCAR 关键标签：
  - `ELPH_RUN`：启动电子–声子计算；
  - `ELPH_MODE = TRANSPORT`：设置输运计算（自动设定合理默认值）；
  - `ELPH_SELFEN_CARRIER_DEN`：载流子密度（空穴/电子）；
  - `ELPH_SCATTERING_APPROX`：弛豫时间近似方式，此处为 CRTA；
  - `TRANSPORT_RELAXATION_TIME`：输运弛豫时间。
- 单次 `vasp_std` 运行内完成五步：
  1. 电子最小化得基态；
  2. 在 `KPOINTS_ELPH`（或 `ELPH_KSPACING`）指定的另一套（通常更密的）k 网格上非自洽计算；
  3. 用位置算符与哈密顿量的对易子计算电子群速度；
  4. 在 Gauss–Legendre 网格上评估输运函数、计算 Onsager 系数；
  5. 计算输运系数。
- 输出关键量：`sigma S/m`（电导率）、`seebeck $\mu$V/K`（Seebeck 系数）、`kappa_e W/(m.K)`（电子热导）；
- 300 K 结果：$\sigma = 9.2\times10^6$ S/m；实验电阻率 97.1 nΩ m（CRC Handbook）对应 $1.03\times10^7$ S/m——第一次计算就不错！电子热导计算值 70.4 W/(m K) 对比实验 80.4 W/(m K) 也合理。但这种吻合是**巧合**：CRTA 恰好在 300 K 与实验一致。

**k 点收敛**（`kpoints/NxNxN` 目录）

- KPOINTS 用于自洽基态，KPOINTS_ELPH 用于非自洽与输运计算（也可用 `ELPH_KSPACING` 替代）；INC 新增 `ELPH_ISMEAR`：电子–声子计算前确定费米能级与化学势的展宽方法；
- 电导率对 k 采样极其敏感，因为：(1) 依赖费米面几何；(2) $\delta(\varepsilon_{n\mathbf{k}}-\varepsilon_F)$ 只挑选恰在费米能级的态；
- 在 10–60 网格上测试：收敛并不完全平滑（60x60x60 略有回升）。收敛慢的情形：复杂费米面拓扑、各向异性材料（bcc 铁不是）、轻有效质量（高曲率能带）；
- 生产计算中应监控电导率随 k 网格变化；收敛判据应参照收敛值而非实验值。教程取与收敛值偏差 5% 以内为足够 → bcc 铁需要 30x30x30 或更密的 KPOINTS_ELPH。

**磁性**（`kpoints_mag/NxNxN` 目录）

- 铁是磁性材料（在位磁矩约 2.2 $\mu_B$），忽略磁性是"因错误的原因得到正确结果"。重做共线自旋（铁磁）计算；
- 新增标签：`ISPIN`、`MAGMOM`（初始磁矩）、`KPAR`（k 点并行数）、`IBRION = 6`、`ELPH_POT_GENERATE`（从有限位移算电子–声子势）、`LPHON_DISPERSION`（沿 QPOINTS 路径算声子色散）；
- 结果：非磁结果似乎更贴近实验值，但对铁的物理准确描述必须包含磁性；复杂费米面金属需要极密 k 采样。

**思考题**

1. 输运计算用 `ELPH_MODE` 的哪个设置？（TRANSPORT）
2. 铁的电导率为何对 k 网格收敛慢？
3. 非磁解更接近实验，就能忽略铁的磁性吗？

---

### 示例 12：自能弛豫时间近似（SERTA）

**教学目标**

- 用 SERTA 计算简单金属电导率；
- 对比 CRTA 测试 k 网格收敛；
- 把散射寿命对 Kohn–Sham 能量作图，理解其对电导率的影响。

**理论**

- 从第一性原理计算弛豫时间：在超胞中计算 Kohn–Sham 势对离子位移的导数，以涵盖尽可能多的声子。输运计算改用 `ELPH_SCATTERING_APPROX = SERTA`（[Phys. Rev. B 97 (2018) 121201(R)](https://doi.org/10.1103/PhysRevB.97.121201)）：

$$\frac{1}{\tau^{\text{SERTA}}_{n\mathbf{k}}} = \sum_{n'\nu\mathbf{k}'} \frac{2\pi}{\hbar} w_{n\mathbf{k},n'\mathbf{k}'} |g^{\nu}_{n\mathbf{k},n'\mathbf{k}'}|^2 \left[ (n_{\nu\mathbf{q}}+1-f_{n'\mathbf{k}'})\, \delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{n'\mathbf{k}'}-\hbar\omega_{\nu\mathbf{q}}) + (n_{\nu\mathbf{q}}+f_{n'\mathbf{k}'})\, \delta(\varepsilon_{n\mathbf{k}}-\varepsilon_{n'\mathbf{k}'}+\hbar\omega_{\nu\mathbf{q}}) \right]$$

其中 $f$ 为 Fermi–Dirac 占据、$n_{\nu\mathbf{q}}$ 为 Bose–Einstein 占据；寿命进入输运分布函数 → Onsager 系数 → 电导率。

**流程**

1. **超胞电子–声子势**（`e12_conductivity_iron_serta/supercell`）：自旋极化；`ELPH_POT_GENERATE` + `LPHON_DISPERSION`（配 QPOINTS），生成 `phelel_params.hdf5`，可用 py4vasp 查看声子色散；
2. **声子限制输运**（`kpoints/NxNxN`）：与 CRTA 计算相同，但传入 `phelel_params.hdf5` 并设 `ELPH_SCATTERING_APPROX = SERTA`；教程只运行 10x10x10 与 15x15x15（更大网格提供预计算结果）。

**结果分析**

- 计入声子散射后电导率比实验值高 40%——"更多物理"有时只是"更多问题"；
- 把 300 K 散射寿命对 KS 能量作图（寿命由 Fan 自能虚部换算 $\tau = \hbar / (2|\text{Im}\,\Sigma^{\text{FM}}|)$）：
  - 紫线为 CRTA 假设：所有载流子同一寿命 $10^{-14}$ s；
  - 实际寿命分布很宽，高载流子能量处偏差尤大；多数寿命小于 $2\times10^{-14}$ s，但右上区域有不少载流子寿命高达约 $9\times10^{-14}$ s；
  - 寿命越长电导率越高。SERTA 包含这些高能、高迁移率载流子 → 电导率显著升高 → 与实验偏差更大；CRTA 忽略了它们，与实验吻合纯属巧合；
- k/q 网格由同一 KPOINTS_ELPH 决定：网格越密，$\varepsilon_{n\mathbf{k}}$ 与对应寿命越多，散射可能性与图像分辨率越高。

**思考题**

1. CRTA 与 SERTA 哪个用平均寿命？（CRTA）
2. 寿命范围展宽如何影响电导率？
3. 如何定义 k 与 q 网格？（KPOINTS_ELPH 同时决定）

---

### 示例 13：晶格常数、密度泛函与赝势的影响

**教学目标**

- 在实验晶格常数下计算铁的电导率；
- 理解不同密度泛函如何影响电子结构与输运性质；
- 比较不同铁赝势（standard 与 GW）的电导率。

**任务**：在实验晶格常数下生成电子–声子势，用 SERTA（10/20/30 q 网格）计算 bcc 铁电导率，同时变换赝势（标准 vs GW）与泛函（CA/LDA vs PBE）。

**晶格常数的影响**（从 PBE 计算值换到实验值）

1. **带宽变化**：晶格常数改变原子轨道重叠，影响导带带宽；
2. **费米面改变**：费米面大小与形状随体积变化，直接影响载流子数目与速度；
3. **态密度变化**：费米能级态密度 $N(\varepsilon_F)$ 随体积变化，影响电导率。

**赝势的影响**：不同赝势在以下方面不同——(1) 核–价电子划分（显式处理 vs 冻结在核内的电子数）；(2) 截断能需求（更"硬"的赝势需要更高平面波截断）；(3) 能带精度（对输运尤其重要的是费米能级附近）。

**结果分析**

- **赝势是主导因素**：GW 赝势结果更接近实验。GW 赝势更精确、专为考虑未占据态而制备，在高能区复现全电子散射性质；
- 泛函影响较小：CA（LDA）与 PBE 差别不大，两者都低估实验值，且收敛后更明显；无穷 k 极限下差别不大。LDA 对金属常表现良好，但用 LDA 弛豫结构会得到差的晶格常数，泛函选择仍须谨慎；
- **晶格常数影响巨大**：计算晶格常数（PBE）与实验晶格常数下的 SERTA 结果相差约 50%，足以把预测从低估变为高估。用实验晶格常数还是泛函自洽晶格常数，目前并无定论；
- 无穷 k 极限（外推）下：标准赝势低估电导率约一半，但数量级（$10^7$）正确；
- **结论**：示例 11 中 CRTA 与实验的接近是"因错误的原因"——常数弛豫时间、未收敛的 k 网格（收敛后反而变差）、不合适的标准赝势、非实验晶格常数，多重误差恰好抵消。与实验吻合很容易出于错误原因，必须通过仔细收敛测试判定每个因素的误差大小；
- 其他可考虑因素：超胞尺寸、不同磁相、缺陷影响、自洽 k 网格、电子–电子散射（post-DFT 方法）。

**思考题**

1. 如何获得电导率的"无穷 k 点极限"？（外推）
2. GW POTCAR 考虑了标准 POTCAR 缺失的哪些态？（未占据态/高能散射性质）

**Part 4 小结**：输运性质需要比基态计算密得多的 k 采样；CRTA 是电导率计算的好起点；晶格常数通过能带结构对电导率影响显著。

---

### 功能与标签汇总（Electron-phonon Part 4）

| 功能 | 关键标签 | 说明 |
|---|---|---|
| 输运计算 | `ELPH_MODE = TRANSPORT` | 自动设定输运默认值 |
| 散射近似 | `ELPH_SCATTERING_APPROX` | CRTA / SERTA |
| 常数弛豫时间 | `TRANSPORT_RELAXATION_TIME` | CRTA 的唯象寿命 |
| 载流子密度 | `ELPH_SELFEN_CARRIER_DEN` | 掺杂/额外载流子 |
| 输运展宽 | `ELPH_ISMEAR` | 费米能级/化学势确定 |
| 磁性输运 | `ISPIN` + `MAGMOM` | 铁磁 bcc 铁必须自旋极化 |
| 寿命分析 | Fan 自能虚部 | $\tau = \hbar/(2|\text{Im}\,\Sigma|)$ |

---

## 16.5 Part 5：半导体声子限制迁移率与 ZT 热电优值

- 来源：<https://vasp.at/tutorials/latest/electron-phonon/part5/>
- 教程输入文件下载：[electron-phonon-part5.zip](https://vasp.at/tutorials/latest/electron-phonon-part5.zip)

与金属不同，半导体在 $T\neq0$ 时本征无载流子，需靠掺杂（n 型加电子、p 型加空穴）提高电荷密度。载流子在材料中运动的难易程度用迁移率 $\mu_e$（电子）或 $\mu_h$（空穴）描述，是半导体器件性能的关键指标。把带求和限制在导带或价带，电导率分为电子与空穴贡献：

| 量 | 定义 | 载流子密度 |
|---|---|---|
| 电子迁移率 $\mu_e$ | $\mu_e = \sigma_{n\in\text{CB}} / n_e$ | $n_e = \frac{1}{\Omega N_\mathbf{k}} \sum_{\mathbf{k}n\in\text{CB}} f(\varepsilon,T,\mu)$ |
| 空穴迁移率 $\mu_h$ | $\mu_h = \sigma_{n\in\text{VB}} / n_h$ | $n_h = \frac{1}{\Omega N_\mathbf{k}} \sum_{\mathbf{k}n\in\text{VB}} [1-f(\varepsilon,T,\mu)]$ |

本部分聚焦用线性化 Boltzmann 输运方程（BTE）计算迁移率，重点是电子–声子相互作用的影响（[Rep. Prog. Phys. 83 036501 (2020)](https://doi.org/10.1088/1361-6633/ab6a43)）。

---

### 示例 14：CRTA 计算迁移率（硅）

**教学目标**

- 用 CRTA 计算硅的迁移率；
- 对能量积分网格与 k 点做数值收敛测试；
- 比较两种输运分布函数积分方案（线性网格 + Simpson vs Gauss–Legendre）。

**要点**

- CRTA 下迁移率只受能带几何影响（载流子轻重、谷各向异性、应变引起的能带重排），是"能带结构指纹"——描述无散射时晶体的行为；对半导体 CRTA 很粗糙（态密度行为尖锐）；
- 输入位于 `$TUTORIALS/electron-phonon/e14_mobility_si_crta/`，关键标签：
  - `ELPH_SCATTERING_APPROX = CRTA`：无需事先算超胞电子–声子势；
  - `ELPH_TRANSPORT_DRIVER = 2`：用 Gauss–Legendre 积分计算输运分布函数的能量积分（默认）；`= 1` 为线性网格；
  - `ELPH_SELFEN_TEMPS`：温度列表；
  - `ELPH_ISMEAR = -24`：电子–声子计算前确定费米能级/化学势的展宽方法，是四面体方法的更省内存/时间实现（`-15` 约慢 6–7 倍）；
  - `ELPH_SELFEN_CARRIER_DEN`：额外载流子数列表。**注意**：正值为 n 掺杂（加电子）、负值为 p 掺杂（加空穴），与文献习惯相反。

**能量网格收敛**

- 迁移率依赖载流子密度：高密度下是真实物理行为，低密度下应与密度无关；低密度区出现"鼓包"提示需对 `ELPH_TRANSPORT_NEDOS`（能量网格点数）收敛；
- CRTA 数值很差：实验电子/空穴迁移率约 1400 / 500 cm$^2$/(V s)（[Phys. Rev. B 97 121201(R) (2018)](https://doi.org/10.1103/PhysRevB.97.121201)），CRTA 只给出约 99 / 72——因为弛豫时间固定为 `1E-14`，调该参数可恢复实验值；

**线性网格 + Simpson 积分**（`linear_grid` 目录，51/501/5001/10001 点数）

- Onsager 系数 $L_{ij} = \int d\epsilon\, \mathcal{T}(\epsilon)\, (\epsilon-\mu)^{i+j-2}(-\partial f^0/\partial\epsilon)$ 用 Simpson 积分；
- 相关标签：`ELPH_TRANSPORT_DRIVER = 1`、`ELPH_FERMI_NEDOS` 与 `ELPH_TRANSPORT_NEDOS`（线性网格点数）、`ELPH_TRANSPORT_DFERMI_TOL`（Onsager 系数的能量窗口）；
- 蓝/青/红线重合 → 能量网格收敛。但**缺点**：能量范围由 `ELPH_TRANSPORT_DFERMI_TOL` 决定（排除化学势周围 Fermi–Dirac 导数积分的某一比例），只增加点数不够，还必须对 `ELPH_TRANSPORT_DFERMI_TOL` 双重收敛（如从 1e-6 收紧到 1e-8，低密度区鼓包消失）；
- CRTA 收敛的网格是 SERTA 的最小尝试起点。

**Gauss–Legendre 积分**（`gauss_legendre` 目录）

- `ELPH_TRANSPORT_DRIVER = 2`（默认）：按点数收敛较慢，但无需额外调能量范围参数（从两个参数减为一个，收敛更简单）；线性网格可加速但需要更多经验；
- 对 SERTA 尤其重要：寿命计算昂贵（每个 $\tau_{n\mathbf{k}}$ 要对 q 求和），能量范围用于挑选需要计算哪些 $\tau_{n\mathbf{k}}$；Gauss–Legendre 网格自动适应 Fermi 窗口宽度，数值高效，只需通过 `TRANSPORT_NEDOS` 定义点数（不必手动设 `ELPH_TRANSPORT_DFERMI_TOL`、`ELPH_TRANSPORT_EMIN/EMAX`）；
- 所用网格可在 OUTCAR 电子–声子段末尾（电子–声子累加器行）看到，标注为 `Linear grids` 或 `Gauss-Legendre grids`。

**k 点收敛**（q-mesh 12–48 目录）

- 从 12x12x12 起电子与空穴迁移率开始收敛，但完全收敛需要更密网格；
- 更方便的做法：把迁移率对**网格划分数倒数**作图并线性外推到无穷 k 采样；
- 真实体系中电离杂质（缺陷）也影响散射：掺杂同时增加载流子与缺陷，超过某点后收益递减、迁移率反被限制。

**思考题**

1. 用什么 INCAR 标签定义 Gauss–Legendre 能量网格？
2. 如何把迁移率外推到无穷 k 采样？

---

### 示例 15：SERTA 计算迁移率（硅）

**教学目标**

- 在超胞上计算电子–声子势；
- 在 SERTA 下对 k 与 q 网格联合收敛迁移率；
- 观察空穴与电子迁移率随 k 收敛的变化。

**流程**

1. **超胞电子–声子势**（`e15_mobility_si_serta/supercell`，2x2x2 Si 超胞，约 5 分钟）：教程提供 `phelel_params.hdf5`；
2. **k/q 联合收敛**（`kpoints` 目录，12–48 网格）：`ELPH_SCATTERING_APPROX = SERTA` 时从 `phelel_params.hdf5` 读取电子–声子势，并插值到连接 KPOINTS_ELPH 所有 k 点的 q 网格上——**改 KPOINTS_ELPH 同时调整 k 与 q 网格**，可系统研究迁移率对网格密度的收敛。教程只运行 12x12x12 与 18x18x18，其余用预计算结果。

**结果分析**

- 与 CRTA 曲线相比，数值问题消除：低密度下迁移率与载流子密度无关；
- k 网格加密，电子迁移率升高直至收敛；
- 电子迁移率远大于空穴，源于有效质量与平均弛豫时间的差异（[Phys. Rev. B 68, 125210 (2003)](https://doi.org/10.1103/PhysRevB.68.125210)）；空穴收敛所需动量低于电子；
- 动量为 q 的声子散射电子：q 网格越大，参与散射的声子模式越多，计算越昂贵——应在收敛前提下用尽可能小的 q 网格。

**思考题**

1. 用哪个 INCAR 标签定义载流子密度？（`ELPH_SELFEN_CARRIER_DEN`）
2. 电子–声子计算中 q 网格是哪种粒子的动量？（声子）

---

### 示例 16：热电材料与 ZT 优值（Bi2Te3）

**教学目标**

- 计算 Bi2Te3 有无自旋–轨道耦合（SOC）的能带结构；
- 在 CRTA 下计算各输运系数；
- 在给定温度下把 DOS、电导率、热电势、ZT 对化学势作图；
- 绘制 ZT 随载流子浓度与温度的热图。

**理论背景**

- 热电材料把热转化为电，效率由热电优值描述：

$$ZT = \frac{\sigma S^2 T}{\kappa_e + \kappa_l}$$

其中 $\sigma$ 为电导率、$S$ 为 Seebeck 系数、$\kappa_e$ 为电子热导、$\kappa_l$ 为晶格热导（可用 Müller–Plathe 方法或 phon3py 计算）；
- 典型热电材料 ZT 在 1–1.5；ZT > 2 的新材料对废热回收、无运动部件制冷器意义重大；
- ZT 把电导率、Seebeck、热导等输运系数（均可追溯到 Onsager 系数与输运分布函数）组合在一起；改进热电性能是经典优化问题（[PNAS 93, 7436 (1996)](https://www.pnas.org/doi/10.1073/pnas.93.15.7436)）。本例按 Sofo 等人的做法（[Phys. Rev. B 68, 125210 (2003)](https://doi.org/10.1103/PhysRevB.68.125210)）用 CRTA 计算 Bi2Te3。

**第一步：能带（含/不含 SOC）**（`e16_thermoelectric_bi2te3/bands` 与 `bands_soc`）

- 新标签：`LSORBIT`（开启自旋–轨道耦合）、`MAGMOM`（初始磁矩）；
- 结果：加入 SOC 后带隙改变，价带升至 0 的"pocket"从 $\Gamma$ 点移到 T 点附近 → 改变简并度与有效电荷，进而影响电导率。

**第二步：输运系数与 ZT 对化学势**（`crta` 目录）

- 新标签：`ELPH_SELFEN_MU_RANGE`：把化学势设为相对费米能级的能量移动范围（费米能级是 $T\to0$ 且无额外载流子时化学势的极限）；
- 300 K 下绘制 DOS、电导率、Seebeck、功率因子、ZT 对化学势的曲线：

![Bi2Te3 输运系数对化学势（300 K）](/imgs/2026-08-05/fig53.png)

- 各量均由 Onsager 系数联系，图形与 Sofo 等人 Fig. 6 非常相似（加密 k 网格、用实验晶格常数可进一步改善）；
- 电导率实际上是被温度平滑后的 DOS：载流子越多（离费米能越远）电导越大；
- Seebeck 与电导率相关：设 $\sigma = a(\varepsilon-\mu)+b$，则 $\sigma S^2 \propto a^2/b^2$；Seebeck 峰与功率因子峰不重合；ZT 有两个峰，分别在化学势正/负侧，对应空穴与电子载流子峰；
- 化学势与载流子密度直接（非线性）相关：

![化学势与载流子密度关系](/imgs/2026-08-05/fig54.png)

- 因此二者只能二选一设定：给 `ELPH_SELFEN_MU` 则计算达到该化学势移动所需载流子数；给 `ELPH_SELFEN_CARRIER_PER_CELL` 或 `ELPH_SELFEN_CARRIER_DEN` 则计算化学势移动。

**第三步：ZT 热图（载流子密度 x 温度）**

- 先按化学势移动计算不同温度下的 ZT 热图（热图中 $\kappa_l$ 取 1）：

![ZT 随化学势与温度热图](/imgs/2026-08-05/fig55.png)

- 理论计算常对化学势作图，但实验上更常直接控制载流子浓度（掺杂），因此对载流子密度作图更实用；
- 改用 `ELPH_SELFEN_CARRIER_DEN`：`ELPH_SELFEN_CARRIER_DEN_RANGE` 在对数网格上生成空穴（$-10^{22}$ 到 $-10^{17}$）与电子（$10^{17}$ 到 $10^{22}$）的载流子密度列表（也可用 `ELPH_SELFEN_CARRIER_DEN` 逐点定义）；
- 载流子密度版热图：

![ZT 随载流子密度与温度热图](/imgs/2026-08-05/fig56.png)

- 化学势与载流子密度关系高度非线性（尤其带边低密度区），对载流子密度作图便于与实验对比、寻找最优掺杂；两种表示包含相同物理信息，峰谷特征一一对应。

**思考题**

1. 用什么 INCAR 标签定义化学势范围？（`ELPH_SELFEN_MU_RANGE`）
2. 化学势–载流子密度关系如何影响 ZT 图的解读（尤其带边附近）？

---

### 功能与标签汇总（Electron-phonon Part 5）

| 功能 | 关键标签 | 说明 |
|---|---|---|
| 输运积分方案 | `ELPH_TRANSPORT_DRIVER = 1 / 2` | 线性网格+Simpson / Gauss–Legendre（默认） |
| 能量网格 | `TRANSPORT_NEDOS`、`ELPH_FERMI_NEDOS` | 网格点数 |
| 能量窗口 | `ELPH_TRANSPORT_DFERMI_TOL`、`ELPH_TRANSPORT_EMIN/EMAX` | 仅线性网格需要手动收敛 |
| 费米能级方法 | `ELPH_ISMEAR = -24` | 高效四面体法 |
| 载流子设定 | `ELPH_SELFEN_CARRIER_DEN`、`ELPH_SELFEN_CARRIER_DEN_RANGE`、`ELPH_SELFEN_CARRIER_PER_CELL` | 正 n 型 / 负 p 型（与文献相反） |
| 化学势设定 | `ELPH_SELFEN_MU`、`ELPH_SELFEN_MU_RANGE` | 与载流子密度二选一 |
| SOC 能带 | `LSORBIT` | Bi2Te3 等重元素材料必须 |
| 晶格热导 | Müller–Plathe / phon3py | ZT 分母需要 |

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

>本文档由 VASP 官方教程 （vasp.at/tutorials/latest） 整理而成，教程内容版权归 VASP 团队所有。*
<!--stackedit_data:
eyJoaXN0b3J5IjpbNzQ2MzA1Nzk4XX0=
-->