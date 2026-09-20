# 文献总结-高熵合金中的短程有序：理论表述及其在 Mo-Nb-Ta-V-W 体系中的应用-2017

## 一、文献信息

- **英文原文标题：** *Short-Range Order in High Entropy Alloys: Theoretical Formulation and Application to Mo-Nb-Ta-V-W System*
- **中文标题：** 高熵合金中的短程有序：理论表述及其在 Mo-Nb-Ta-V-W 体系中的应用
- **作者：** A. Fernández-Caballero、J. S. Wróbel、P. M. Mummery、D. Nguyen-Manh
- **期刊：** *Journal of Phase Equilibria and Diffusion*
- **年份、卷、期、页码：** 2017，38(4)，391–403
- **DOI：** [https://doi.org/10.1007/s11669-017-0582-3](https://doi.org/10.1007/s11669-017-0582-3)
- **投稿日期：** 2017-04-27
- **返修日期：** 论文原文未提供。
- **接收日期：** 2017-07-10
- **在线发表日期：** 2017-07-27
- **通讯作者：** D. Nguyen-Manh，duc.nguyen@ukaea.uk

### 数据与代码链接

- **arXiv 预印本：** [https://arxiv.org/abs/1705.01844](https://arxiv.org/abs/1705.01844)
- **论文专属数据仓库：** 论文未给出。
- **论文专属代码仓库：** 论文未给出。
- **所用软件：** 论文明确说明 DFT 计算使用 VASP，CE 拟合使用 ATAT；原文没有为本研究提供输入文件、拟合脚本、MC 脚本或可下载训练数据库的链接。

## 二、论文内容详细总结

### 2.1 研究目标

论文研究难熔体心立方高熵合金 Mo-Nb-Ta-V-W 的局域化学有序。作者关注的问题是：高温下近似随机的多组元固溶体在降温后会如何偏离理想随机分布，这些偏离如何用 Warren-Cowley 短程有序（SRO）参数描述，以及五元体系中的多体相互作用会如何改变从二元合金直觉得到的判断。

为此，作者把第一性原理数据库、团簇展开（CE）Hamiltonian、准正则 Monte Carlo（MC）模拟和解析 SRO 公式组合起来。论文的主要理论贡献是从五元 CE 的点关联函数与对关联函数出发，显式推导 10 个不同元素对的平均对概率及 Warren-Cowley SRO 参数；同时给出四元体系 6 个 SRO 参数的显式表达式。

### 2.2 理论与计算链条

作者先用 DFT 计算 428 个 bcc-like 结构的能量。这些结构来自不同的二元、三元和四元构型。DFT 使用 VASP 中的 PAW 方法，交换-相关采用 GGA-PBE；V、Nb、Ta 使用包含半芯态 $p$ 电子的 11 价电子 PAW 势，Mo 和 W 使用 12 价电子势。

随后使用 ATAT 的结构反演方法把 DFT 混合焓映射到 CE Hamiltonian。模型包含 5 个点 ECI、30 个对 ECI 和 40 个三体 ECI，DFT 与 CE 能量之间的交叉验证误差约为 $8\ \mathrm{meV/atom}$。论文指出，主导贡献来自第一和第二近邻对相互作用，第三近邻对和三体相互作用的量级明显更小；但多体项仍被保留，因此该模型超出了只包含第一近邻二元对参数的简单 Ising 模型。

MC 模拟在固定合金组成的准正则框架中进行，通过温度扫描获得点关联函数、对关联函数、混合焓及 SRO。较大体系的预熔构型由随机数生成；论文只说明较小模拟胞可以使用 SQS，并未给出本研究实际 MC 超胞的原子数、每个温度的 MC 步数、平衡步数或采样间隔。

### 2.3 SRO 参数的物理含义

论文采用：

$$
\alpha_{2,m}^{pq}=1-\frac{y_m^{pq}}{c_pc_q},
$$

其中，$c_p$ 和 $c_q$ 为元素平均浓度，$y_m^{pq}$ 为第 $m$ 近邻壳层中 $p-q$ 对的平均概率。由此：

- $\alpha_{2,m}^{pq}=0$ 表示随机分布；
- $\alpha_{2,m}^{pq}<0$ 表示异类 $p-q$ 原子对有序，即彼此近邻的概率高于随机值；
- $\alpha_{2,m}^{pq}>0$ 表示同类原子聚集或 $p-q$ 相互偏离，即具有偏析趋势。

对于 bcc 晶格，作者还把第一近邻和第二近邻合并成实验上更容易通过漫散射测得的平均参数：

$$
\alpha_{1+2}^{pq}=\frac{8\alpha_{2,1}^{pq}+6\alpha_{2,2}^{pq}}{14}.
$$

### 2.4 五元 Mo-Nb-Ta-V-W 的主要结果

在高温极限，所有 SRO 参数趋近 0，对应理想随机单相固溶体。低于约 750 K 后，混合焓曲线出现与有序-无序转变相联系的拐点，部分有序构型比等原子随机固溶体更稳定。

第一近邻中，Mo-Ta 的 SRO 最负，说明 Mo 周围出现 Ta 的概率显著高于随机值；V-W 也具有强烈负 SRO，Mo-Nb 次之。Mo-Ta 倾向与 bcc 二元 Mo-Ta 的负混合焓及 $A2\rightarrow B2$ 有序化一致。低温下 Nb-V、Mo-W 和 Ta-V 的 SRO 为正，表明这些元素对具有偏析倾向。

结果也显示多组元环境不能简单由二元相图或二元基态结构外推。例如，低温时五元体系中的 V-W 第一近邻 SRO 为负，而 Ta-W 为正，这与仅根据对应二元体系的 $B32$ 和 $B2_3$ 基态作出的直观判断并不一致。作者把这种复杂性归因于 CE 中多种对相互作用与三体相互作用的共同作用。

### 2.5 四元子体系的主要结果

在不含 Nb 的 Mo-Ta-V-W 中，Mo-Ta 与 V-W 仍是最有利的异类近邻对；Mo-W 和 Ta-V 的 SRO 在低温下为正。低温时 Ta-W 与 Mo-V 的 SRO 也转正。

在不含 Ta 的 Mo-Nb-V-W 中，Mo-Nb 取代 Mo-Ta 成为最显著的负 SRO，并与 V-W 的负 SRO 竞争。Nb-V 的 SRO 保持显著为正，与五元体系中的偏析趋势一致。

### 2.6 第二近邻效应

把第二近邻纳入平均 SRO 后，Mo-Ta 的负值绝对值明显减小，因为 Mo-Ta 的第二近邻贡献为正。作者据此认为，在 Mo-Ta 的 $B2$ 型局域环境中，第二近邻壳层由相同物种 Mo-Mo 或 Ta-Ta 构成，这与平均 SRO 的变化一致。

V-W 的行为不同：在大部分温度范围内，第二近邻 V-W SRO 仍为负，仅在低于约 100 K 的五元体系中转为正。作者认为，100 K 以上第二近邻异类 V-W 对的存在支持 $B32$ 型局域化学环境。在两个四元子体系中，$B32$-like 到 $B2$-like 的低温转变发生在 200 K 以下。

### 2.7 论文结论的适用边界

- 模型基于刚性 bcc 晶格，只处理化学占位偏离；未在 MC 中显式考虑热致离位和晶格振动。
- 论文未给出完整 MC 运行参数，因此无法仅凭本文精确复现采样长度和有限尺寸设置。
- 作者指出，经典分子动力学需要可靠的多组元原子间势，而大体系从头算分子动力学成本过高。
- 作者引用先前工作说明其相关 DFT 晶格畸变预测的原子坐标平均误差约为 1–2 pm，但该误差并不是本文当前刚性晶格 MC 的直接误差评估。
- 论文建议用中子衍射等实验对 Mo-Nb-Ta-V-W 的 SRO 进行更细致验证。

## 三、摘要完整翻译

**原文位置：PDF 第 391 页，标题和作者信息之后。**

在高熵合金（HEA）中，由于元素数目较多，从无序固溶体状态向偏析、析出和有序构型发展的局域化学涨落十分复杂。本工作建立了多组元合金体系的团簇展开（CE）Hamiltonian，以研究 HEA 化学有序随温度的变化；这种变化源于构型熵偏离理想固溶体。论文为五组元合金体系推导了 Warren-Cowley 短程有序（SRO）参数的解析表达式。作者把该理论表述用于研究 MoNbTaVW 及其四元子体系中 10 个不同 SRO 参数的演化；这些结果由 CE 与第一性原理形式体系相结合的 Monte Carlo 模拟获得。

预测表明，最强的化学 SRO 参数对应第一近邻 Mo-Ta 原子对，这与该二元体系 $B2$ 结构中较大的混合焓值一致。对于所考察 bcc HEA 中 Mo-Ta 对存在 $B2$ 相的预测，还得到第二近邻壳层对平均 SRO 的正贡献的支持。值得注意的是，与 Mo-Ta 对相比，V-W 对在第一和第二近邻壳层上的平均 SRO 参数也显著为负。HEA 中这一结果可以通过有序型 $B32$ 相的存在加以合理化和讨论；该相已被预测为等原子二元 bcc V-W 体系的基态结构。

## 四、引言完整翻译

**原文位置：第 1 节 Introduction，PDF 第 391–392 页。**

以多种主元素构成、并主要呈单一固溶体相的合金被称为高熵合金（HEA）。这类合金已经成为合金开发的新兴领域，并涉及许多基础概念，其中包括熵效应的起源，以及复杂浓集合金热力学分析与其微观组织性质之间的相互作用 [1-3]。这类新材料在 2004 年经 Cantor 等 [4] 和 Yeh 等 [5] 的工作首次受到广泛关注。其基本概念是：较高的混合构型熵应当使 fcc 或 bcc 等简单固溶体相相对于可能导致脆化的金属间相得到稳定。然而，近期实验研究对不同 HEA 中固溶体的熵稳定效应提出了质疑 [6-9]。特别是，在最初被认为是单一 fcc 相的 CrMnFeCoNi HEA 中，长期热处理后报道了两种不同的富 Cr 析出相，即 $M_{23}C_6$ 和 sigma 相 [7]。

尽管如此，HEA 不仅为合金设计提供了新颖而令人期待的方法，也因发现具有异常且有吸引力的物理、化学、磁学和力学性质的合金而受到关注 [3,9]。异常力学性质的一个实例是：bcc 结构难熔高熵合金 MoNbTaVW 的屈服强度高于其任一单质组元。此类 HEA 具有高熔点，并能在超高温下保持优异屈服强度；这通常归因于固溶强化机制和相关晶格应变 [10]。一般而言，HEA 的特征不仅是熵值较高，还包括因不同尺寸元素混合而产生的高原子尺度应力 [11,12]。

因此，为了从有序与无序角度理解 HEA 的局域原子结构，仍有许多基本问题需要解决。HEA 中的每种元素都倾向于占据能够使其位点能和键能最小的位置，而这些能量又取决于优选的局域化学环境、原子间相互作用以及不同组元物种的原子体积。为了通过实验确定元素分布，研究者已采用多种互补技术来确定 HEA 的局域结构，包括中子散射、高能同步辐射 X 射线衍射、原子探针层析，以及与 X 射线能量色散谱（XEDS）联用的透射电子显微镜 [9,13,14]。

近年来，HEA 因其辐照抗性优于常规单相 Fe-Cr-Ni 奥氏体不锈钢而受到广泛关注，这使它们成为高温裂变与聚变应用的潜在候选材料 [15,16]。级联事件后的辐照稳定性，可能归因于 HEA 中不同原子尺寸引起的高原子尺度应力，以及它们在辐照诱导热尖峰后以较高速率发生非晶化和再结晶的倾向 [17]。值得一提的是，在磁约束聚变堆设计中，W 基合金被认为是面向等离子体材料的优选方案 [18]。除高熔点之外，其他受到关注的原因还包括高热导率、低活化、低氚滞留和低溅射产额，这些性质与核聚变电站结构材料的辐照损伤有关 [19,20]。

本工作采用了此前为研究多组元合金相稳定性而开发的预测模型 [11,21]，把第一性原理计算与团簇展开（CE）方法结合起来，量化 HEA 中不同原子物种的短程有序。CE Hamiltonian 可以一致地描述多组元合金体系中固溶体和金属间相自由能的焓贡献与熵贡献。特别是，通过利用 ECI 的 Monte Carlo 模拟进行热力学积分，结果表明构型熵贡献通常强烈依赖温度，因此，HEA 研究中以唯象方式采用的理想随机固溶体熵表达式只在高温极限下有效 [21]。非随机构型会产生相分离或化学短程有序（SRO）倾向，这两种趋势都会使构型熵低于理想估计值。

本文安排如下。第 2 节介绍 CE 形式体系，并针对五元合金重点建立多组元体系的 Hamiltonian。第 3 节用平均对关联函数显式推导五组元合金的短程有序公式。第 4 节把这些 SRO 表达式应用于特定 HEA MoNbTaVW 及其相应的 5 个四元子体系。第 5 节讨论第二近邻对 SRO 性质的贡献。第 6 节总结本文结果。

## 五、方法部分完整翻译

论文没有单独设置 “Methods” 章节。与方法对应的是第 2 节 **Cluster Expansion Hamiltonian for Multicomponent Alloys** 和第 3 节 **Short Range Order Parameters for Multicomponent Alloy**，位于 PDF 第 392–397 页。以下按原文顺序完整翻译。

### 2. 多组元合金的团簇展开 Hamiltonian

多组元合金的相稳定性可以通过密度泛函理论（DFT）技术与基于 CE 形式体系的晶格统计模拟相结合来研究 [11,21]。DFT 计算使用 Vienna *Ab-initio* Simulation Package（VASP）中实现的投影缀加波（PAW）方法 [22-27]。作者使用包含半芯态 $p$ 电子贡献的 PAW 势：V、Nb 和 Ta 把 11 个电子作为价电子，Mo 和 W 把 12 个电子作为价电子。交换与关联相互作用使用广义梯度近似 GGA-PBE [28] 处理。

可以用 DFT 计算的五组元合金混合焓定义为：

$$
\begin{aligned}
\Delta H_{\mathrm{mix}}(\vec{\sigma})={}&
E_{\mathrm{tot}}^{\mathrm{lat}}
\left(A_{c_A}B_{c_B}C_{c_C}D_{c_D}E_{c_E},\vec{\sigma}\right)
-c_AE_{\mathrm{tot}}^{\mathrm{lat}}(A)
-c_BE_{\mathrm{tot}}^{\mathrm{lat}}(B)\\
&-c_CE_{\mathrm{tot}}^{\mathrm{lat}}(C)
-c_DE_{\mathrm{tot}}^{\mathrm{lat}}(D)
-c_EE_{\mathrm{tot}}^{\mathrm{lat}}(E).
\end{aligned}\tag{1}
$$

其中，$c_A$、$c_B$、$c_C$、$c_D$ 和 $c_E$ 分别是合金组元 A、B、C、D 和 E 的平均浓度，$E_{\mathrm{tot}}^{\mathrm{lattice}}$ 是所考虑结构的每原子总能。这里，向量 $\vec{\sigma}$ 定义给定晶格，例如体心立方（bcc）晶格，上的合金构型。

在团簇展开形式体系中 [29-32]，式 (1) 的构型混合焓可以用不同团簇相互作用能表示：

$$
\Delta H_{\mathrm{mix}}(\vec{\sigma})=
\sum_{\omega}J_{\omega}m_{\omega}
\left\langle\Gamma_{\omega'}(\vec{\sigma})\right\rangle_{\omega}.
\tag{2}
$$

式 (2) 对底层晶格中经对称操作后互不等价的所有团簇 $\omega$ 求和。$m_{\omega}$ 表示多重度，即与 $\omega$ 对称等价的团簇数目；$J_{\omega}$ 是团簇 $\omega$ 对应的有效团簇相互作用。$\left\langle\Gamma_{\omega'}(\vec{\sigma})\right\rangle_{\omega}$ 是团簇函数：它由特定团簇 $\omega$ 上占位变量的点函数 $\gamma_i(\sigma_p)$ 的乘积定义，并对与 $\omega$ 等价的原子团簇 $\omega'$ 取平均。对应合金构型 $\vec{\sigma}$ 的团簇关联函数一般表达式为：

$$
\left\langle
\Gamma_{|\omega|,m}^{(ijk\cdots)}(\vec{\sigma})
\right\rangle
=\sum_{pqr\cdots}
\gamma_i(\sigma_p)\gamma_j(\sigma_q)\gamma_k(\sigma_r)\cdots
y_m^{pqr\cdots}.
\tag{3}
$$

这里，$m$ 是与团簇壳层中各晶格点原子标记构型相对应的整数。$y_m^{pqr\cdots}$ 表示在由 $m$ 指定的壳层原子构型中找到原子物种 $pqr\cdots$ 的温度依赖概率。对于一般的 $n$ 组元体系，正交点函数 $\gamma$ 按参考文献 [31] 定义为：

$$
\gamma_{j,n}(\sigma_i)=
\begin{cases}
1, & j=0,\\
-\cos\left(2\pi\left\lceil\dfrac{j}{2}\right\rceil\dfrac{\sigma_i}{n}\right),
& j>0\ \text{且为奇数},\\
-\sin\left(2\pi\left\lceil\dfrac{j}{2}\right\rceil\dfrac{\sigma_i}{n}\right),
& j>0\ \text{且为偶数}.
\end{cases}
\tag{4}
$$

由式 (3) 和式 (4)，点团簇和对团簇的团簇关联函数分别为：

$$
\left\langle\Gamma_{1,m}^{(i)}(\vec{\sigma})\right\rangle
=\sum_p\gamma_i(\sigma_p)y_m^p
=\sum_p\gamma_i(\sigma_p)c_p,
\tag{5}
$$

$$
\left\langle\Gamma_{2,m}^{(ij)}(\vec{\sigma})\right\rangle
=\sum_{pq}\gamma_i(\sigma_p)\gamma_j(\sigma_q)y_m^{pq}.
\tag{6}
$$

在式 (5) 中，平均单点团簇函数可以用合金浓度 $c_p$ 表示。对于式 (6) 的平均成对团簇函数，指标 $m$ 表示第 $m$ 近邻原子对；该原子对以概率 $y_m^{pq}$ 出现，其中 $p$ 型原子位于第一个位点，$q$ 型原子位于第 $m$ 个原子对的第二个位点。需要注意，只有当式 (2) 的 Hamiltonian 包含对相互作用时，式 (6) 所述对概率 $y_m^{pq}$ 才相互独立并与 CE 形式体系一致。表 1 随后表明，与两体相互作用参数相比，三体相互作用参数的量值小得多。因此，尽管 CE Hamiltonian 包含三体贡献，在本研究中，成对概率仍在确定多组元合金准独立构型方面起主导作用。三体有效团簇相互作用对 Mo-Nb-Ta-V-W 合金短程有序分析的影响将在第 4 节和第 5 节讨论。

利用式 (5)，五组元合金的单点关联函数可以用体系内所有元素的平均原子浓度表示为：

$$
\left\langle\Gamma_{1,1}^{(0)}\right\rangle=1,
\tag{7a}
$$

$$
\left\langle\Gamma_{1,1}^{(1)}\right\rangle=
\frac{1}{4}\left[\phi_-\left(c_B+c_E\right)+
\phi_+\left(c_C+c_D\right)-4c_A\right],
\tag{7b}
$$

$$
\left\langle\Gamma_{1,1}^{(2)}\right\rangle=
\sqrt{\chi_-}(c_D-c_C)+\sqrt{\chi_+}(c_E-c_B),
\tag{7c}
$$

$$
\left\langle\Gamma_{1,1}^{(3)}\right\rangle=
\frac{1}{4}\left[\phi_-\left(c_C+c_D\right)+
\phi_+\left(c_B+c_E\right)-4c_A\right],
\tag{7d}
$$

$$
\left\langle\Gamma_{1,1}^{(4)}\right\rangle=
\sqrt{\chi_-}(c_E-c_B)+\sqrt{\chi_+}(c_C-c_D).
\tag{7e}
$$

这里使用记号 $\phi_{\pm}=1\pm\sqrt{5}$，$\chi_{\pm}=5/8\pm\sqrt{5}/8$。由式 (6)，成对团簇函数可以用对概率 $y_m^{pq}$ 表示如下。

<details>
<summary>展开查看式 (8a)–(8j)：五元体系的 10 个成对团簇函数</summary>

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{11}\right\rangle=\frac{1}{16}\big[&
2\phi_-\phi_+(y_m^{BC}+y_m^{BD}+y_m^{CE}+y_m^{DE})
-8\phi_-(y_m^{AB}+y_m^{AE})\\
&+\phi_-^2(y_m^{BB}+2y_m^{BE}+y_m^{EE})
-8\phi_+(y_m^{AC}+y_m^{AD})\\
&+\phi_+^2(y_m^{CC}+2y_m^{CD}+y_m^{DD})+16y_m^{AA}\big].
\end{aligned}\tag{8a}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{12}\right\rangle=\frac{1}{4}\big[&
\phi_-\sqrt{\chi_+}(y_m^{EE}-y_m^{BB})
+\sqrt{\chi_-}\phi_+(y_m^{DD}-y_m^{CC})\\
&+4\sqrt{\chi_-}(y_m^{AC}-y_m^{AD})
+\sqrt{\chi_-}\phi_-(-y_m^{BC}+y_m^{BD}-y_m^{CE}+y_m^{DE})\\
&+4\sqrt{\chi_+}(y_m^{AB}-y_m^{AE})
+\sqrt{\chi_+}\phi_+(-y_m^{BC}-y_m^{BD}+y_m^{CE}+y_m^{DE})\big].
\end{aligned}\tag{8b}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{13}\right\rangle=\frac{1}{16}\big[&
\phi_-\phi_+(y_m^{BB}+2y_m^{BE}+y_m^{CC}+2y_m^{CD}+y_m^{DD}+y_m^{EE})\\
&-4\phi_-(y_m^{AB}+y_m^{AC}+y_m^{AD}+y_m^{AE})
+\phi_-^2(y_m^{BC}+y_m^{BD}+y_m^{CE}+y_m^{DE})\\
&-4\phi_+(y_m^{AB}+y_m^{AC}+y_m^{AD}+y_m^{AE})
+\phi_+^2(y_m^{BC}+y_m^{BD}+y_m^{CE}+y_m^{DE})+16y_m^{AA}\big].
\end{aligned}\tag{8c}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{14}\right\rangle=\frac{1}{4}\big[&
\phi_-\sqrt{\chi_+}(y_m^{BC}-y_m^{DE})
+\sqrt{\chi_-}\phi_+(y_m^{DE}-y_m^{BC})\\
&+(\phi_-\sqrt{\chi_+}+\sqrt{\chi_-}\phi_+)(y_m^{CE}-y_m^{BD})
+4\sqrt{\chi_-}(y_m^{AB}-y_m^{AE})\\
&+\sqrt{\chi_-}\phi_-(y_m^{EE}-y_m^{BB})
-4\sqrt{\chi_+}(y_m^{AC}-y_m^{AD})
+\sqrt{\chi_+}\phi_+(y_m^{CC}-y_m^{DD})\big].
\end{aligned}\tag{8d}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{22}\right\rangle={}&
2\sqrt{\chi_-\chi_+}(y_m^{BC}-y_m^{BD}-y_m^{CE}+y_m^{DE})\\
&+\chi_-(y_m^{CC}-2y_m^{CD}+y_m^{DD})
+\chi_+(y_m^{BB}-2y_m^{BE}+y_m^{EE}).
\end{aligned}\tag{8e}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{23}\right\rangle=\frac{1}{4}\big[&
(\phi_-\sqrt{\chi_+}+\sqrt{\chi_-}\phi_+)(y_m^{DE}-y_m^{BC})\\
&+\phi_-\sqrt{\chi_+}(y_m^{CE}-y_m^{BD})
+\sqrt{\chi_-}\phi_+(y_m^{BD}-y_m^{CE})\\
&+4\sqrt{\chi_-}(y_m^{AC}-y_m^{AD})
+\sqrt{\chi_-}\phi_-(y_m^{DD}-y_m^{CC})\\
&+4\sqrt{\chi_+}(y_m^{AB}-y_m^{AE})
+\sqrt{\chi_+}\phi_+(y_m^{EE}-y_m^{BB})\big].
\end{aligned}\tag{8f}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{24}\right\rangle={}&
\sqrt{\chi_-\chi_+}(y_m^{BB}-2y_m^{BE}-y_m^{CC}+2y_m^{CD}-y_m^{DD}+y_m^{EE})\\
&+\chi_-(y_m^{BC}-y_m^{BD}-y_m^{CE}+y_m^{DE})
+\chi_+(-y_m^{BC}+y_m^{BD}+y_m^{CE}-y_m^{DE}).
\end{aligned}\tag{8g}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{33}\right\rangle=\frac{1}{16}\big[&
2\phi_-\phi_+(y_m^{BC}+y_m^{BD}+y_m^{CE}+y_m^{DE})
-8\phi_-(y_m^{AC}+y_m^{AD})\\
&+\phi_-^2(y_m^{CC}+2y_m^{CD}+y_m^{DD})
-8\phi_+(y_m^{AB}+y_m^{AE})\\
&+\phi_+^2(y_m^{BB}+2y_m^{BE}+y_m^{EE})+16y_m^{AA}\big].
\end{aligned}\tag{8h}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{34}\right\rangle=\frac{1}{4}\big[&
\sqrt{\chi_-}\phi_+(y_m^{EE}-y_m^{BB})
+\phi_-\sqrt{\chi_+}(y_m^{CC}-y_m^{DD})\\
&+4\sqrt{\chi_-}(y_m^{AB}-y_m^{AE})
+\sqrt{\chi_-}\phi_-(-y_m^{BC}-y_m^{BD}+y_m^{CE}+y_m^{DE})\\
&-4\sqrt{\chi_+}(y_m^{AC}-y_m^{AD})
+\sqrt{\chi_+}\phi_+(y_m^{BC}-y_m^{BD}+y_m^{CE}-y_m^{DE})\big].
\end{aligned}\tag{8i}
$$

$$
\begin{aligned}
\left\langle\Gamma_{2,m}^{44}\right\rangle={}&
\sqrt{\chi_-\chi_+}(-2y_m^{BC}+2y_m^{BD}+2y_m^{CE}-2y_m^{DE})\\
&+\chi_-(y_m^{BB}-2y_m^{BE}+y_m^{EE})
+\chi_+(y_m^{CC}-2y_m^{CD}+y_m^{DD}).
\end{aligned}\tag{8j}
$$

</details>

把式 (2) 按照式 (7) 和式 (8) 给出的平均点团簇函数与平均对团簇函数重写后，五元合金的构型混合焓可以表示为浓度 $c_p$ 和平均对概率 $y_m^{pq}$ 的函数：

$$
\begin{aligned}
\Delta H_{\mathrm{mix}}(\vec{\sigma})={}&
\sum_{i=0}^{4}J_{1,1}^{i}\left\langle\Gamma_{1,1}^{i}(\vec{\sigma})\right\rangle\\
&+\sum_{m\ \mathrm{pairs}}\bigg[
m_{2,m}^{11}J_{2,m}^{11}\left\langle\Gamma_{2,m}^{11}\right\rangle+
m_{2,m}^{12}J_{2,m}^{12}\left\langle\Gamma_{2,m}^{12}\right\rangle+
m_{2,m}^{13}J_{2,m}^{13}\left\langle\Gamma_{2,m}^{13}\right\rangle\\
&\quad+m_{2,m}^{14}J_{2,m}^{14}\left\langle\Gamma_{2,m}^{14}\right\rangle+
m_{2,m}^{22}J_{2,m}^{22}\left\langle\Gamma_{2,m}^{22}\right\rangle+
m_{2,m}^{23}J_{2,m}^{23}\left\langle\Gamma_{2,m}^{23}\right\rangle\\
&\quad+m_{2,m}^{24}J_{2,m}^{24}\left\langle\Gamma_{2,m}^{24}\right\rangle+
m_{2,m}^{33}J_{2,m}^{33}\left\langle\Gamma_{2,m}^{33}\right\rangle+
m_{2,m}^{34}J_{2,m}^{34}\left\langle\Gamma_{2,m}^{34}\right\rangle\\
&\quad+m_{2,m}^{44}J_{2,m}^{44}\left\langle\Gamma_{2,m}^{44}\right\rangle
\bigg]+\sum_{\mathrm{triplets}}\cdots.
\end{aligned}\tag{9}
$$

五组元 MoNbTaVW 体系的 ECI 参数 $J_{\omega}^{(s)}$，通过结构反演方法（SIM）[33,34] 把 428 个 bcc-like 结构的 DFT 能量映射到式 (9) 的 CE Hamiltonian 而得到；这些结构来自不同的二元、三元和四元体系。拟合过程使用 ATAT 软件包 [29] 完成，DFT 与 CE 能量之间的交叉验证误差约为 $8\ \mathrm{meV/atom}$。参考文献 [11] 最初报告的 5 个点 ECI、30 个对 ECI 和 40 个三体 ECI 的数值列于表 1。

![表 1](/imgs/2026-09-20/table1.png)

**表 1。** bcc 五元 Mo-Nb-Ta-V-W 合金的团簇大小 $|\omega|$、修饰 $(s)$、各点坐标、多重度 $m_{|\omega|,n}^{(s)}$ 和 ECI 参数 $J_{|\omega|,n}^{(s)}$，ECI 单位为 meV。

结果表明，成对能量的主导贡献来自第一和第二 bcc 近邻相互作用，而第三近邻对相互作用和三体相互作用则明显较小。值得注意的是，式 (9) 的多组元合金 CE 能量是 Ising-like Hamiltonian 的推广；后者只考虑不同物种之间的最近邻相互作用，例如四组元 HEA MoNbTaW 的情形 [35]。CE Hamiltonian 可用于执行准正则 Monte Carlo 模拟，以研究合金从无序到有序相变过程中的形成自由能。

自由能的评估使用热力学积分算法计算合金构型熵，具体方法此前已有详细说明 [21]。更关键的是，描述局域原子尺度占位相对于平均随机构型偏离程度的化学短程有序参数，可以由 CE 混合自由能求得，并与可用实验数据比较。

### 3. 多组元合金的短程有序参数

一般而言，多组元合金可以在短程和中程尺度上都表现出结构有序，后者的范围达到 2 nm 或更长。理论上，当前 CE Hamiltonian 方法与 Monte Carlo 模拟相结合，可以通过式 (3) 的多体团簇概率函数 $y_m^{pqr\cdots}$ 研究 HEA 的中程有序。然而，中程有序难以在实验中测量，也难以作无歧义解释。本文聚焦于化学短程有序（SRO）分析；该分析此前已成功用于三元 Fe-Cr-Ni 合金，并与可用实验数据进行了比较 [21]。

化学 SRO 程度会同时影响复杂合金的构型熵和混合焓。在原子物种随机占据位点的理想固溶体中，不存在化学 SRO。通常使用 Warren-Cowley 短程有序参数或对关联参数描述 SRO：

$$
\alpha_{2,m}^{pq}=1-\frac{y_m^{pq}}{c_pc_q}.
\tag{10}
$$

这里，$c_p$ 表示平均浓度，$y_m^{pq}$ 表示第 $m$ 近邻壳层的平均对概率，它们分别由式 (5) 和式 (6) 定义。$P_m^{pq}=y_m^{pq}/c_p$ 是条件概率，表示在原子 $p$ 周围的第 $m$ 配位壳层中找到原子 $q$ 的概率。当 $\alpha_{2,m}^{pq}=0$ 时，合金为随机合金，即原子对构型 $m$ 中的元素 $p$ 和 $q$ 以 $c_pc_q$ 的概率出现在合金中。当 $\alpha_{2,m}^{pq}>0$ 时，$p-p$ 与 $q-q$ 原子对具有聚集或偏析倾向；当 $\alpha_{2,m}^{pq}<0$ 时，异类 $p-q$ 原子对具有有序倾向。

短程有序参数可以用平均点关联函数和平均对关联函数表示。式 (10) 中 $y_{2,m}^{pq}$ 的表达式可以通过对式 (8) 和式 (7) 求逆得到。五组元合金体系有 10 个不同的对概率函数，其显式公式如下。

<details>
<summary>展开查看式 (11a)–(11j)：五元体系的 10 个对概率函数</summary>

令 $\kappa=\sqrt{5}$。为保持公式可读性，下式中的所有 $\Gamma$ 均表示原文定义的平均团簇关联函数。

$$
\begin{aligned}
y_{2,m}^{AB}=\frac{1}{50}\big[&
-(\kappa+3)\langle\Gamma_{1,1}^{1}\rangle
+\sqrt{2}\sqrt{\kappa+5}\{2(\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{1,1}^{2}\rangle\}\\
&+\sqrt{10-2\kappa}\{2(\langle\Gamma_{2,m}^{14}\rangle+\langle\Gamma_{2,m}^{34}\rangle)-\langle\Gamma_{1,1}^{4}\rangle\}
+(\kappa-3)\langle\Gamma_{1,1}^{3}\rangle\\
&+2(\kappa-1)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle
-2(\kappa+1)\langle\Gamma_{2,m}^{33}\rangle+2\big].
\end{aligned}\tag{11a}
$$

$$
\begin{aligned}
y_{2,m}^{AC}=\frac{1}{50}\big[&
(\kappa-3)\langle\Gamma_{1,1}^{1}\rangle
+\sqrt{10-2\kappa}\{2(\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{1,1}^{2}\rangle\}\\
&+\sqrt{2}\sqrt{\kappa+5}\{\langle\Gamma_{1,1}^{4}\rangle-2(\langle\Gamma_{2,m}^{14}\rangle+\langle\Gamma_{2,m}^{34}\rangle)\}
-(\kappa+3)\langle\Gamma_{1,1}^{3}\rangle\\
&-2(\kappa+1)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle
+2(\kappa-1)\langle\Gamma_{2,m}^{33}\rangle+2\big].
\end{aligned}\tag{11b}
$$

$$
\begin{aligned}
y_{2,m}^{AD}=\frac{1}{50}\big[&
(\kappa-3)\langle\Gamma_{1,1}^{1}\rangle
+\sqrt{10-2\kappa}\{\langle\Gamma_{1,1}^{2}\rangle-2(\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{23}\rangle)\}\\
&+\sqrt{2}\sqrt{\kappa+5}\{2(\langle\Gamma_{2,m}^{14}\rangle+\langle\Gamma_{2,m}^{34}\rangle)-\langle\Gamma_{1,1}^{4}\rangle\}
-(\kappa+3)\langle\Gamma_{1,1}^{3}\rangle\\
&-2(\kappa+1)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle
+2(\kappa-1)\langle\Gamma_{2,m}^{33}\rangle+2\big].
\end{aligned}\tag{11c}
$$

$$
\begin{aligned}
y_{2,m}^{AE}=\frac{1}{50}\big[&
-(\kappa+3)\langle\Gamma_{1,1}^{1}\rangle
+\sqrt{2}\sqrt{\kappa+5}\{\langle\Gamma_{1,1}^{2}\rangle-2(\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{23}\rangle)\}\\
&+\sqrt{10-2\kappa}\{\langle\Gamma_{1,1}^{4}\rangle-2(\langle\Gamma_{2,m}^{14}\rangle+\langle\Gamma_{2,m}^{34}\rangle)\}
+(\kappa-3)\langle\Gamma_{1,1}^{3}\rangle\\
&+2(\kappa-1)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle
-2(\kappa+1)\langle\Gamma_{2,m}^{33}\rangle+2\big].
\end{aligned}\tag{11d}
$$

$$
\begin{aligned}
y_{2,m}^{BC}=\frac{1}{25}\big[&
\langle\Gamma_{1,1}^{1}\rangle-\sqrt{2\kappa+5}(\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{2,m}^{14}\rangle)
+\sqrt{5-2\kappa}(\langle\Gamma_{1,1}^{4}\rangle-\langle\Gamma_{2,m}^{23}\rangle)\\
&+\langle\Gamma_{1,1}^{3}\rangle-\langle\Gamma_{2,m}^{11}\rangle
-\sqrt{10-2\kappa}\langle\Gamma_{2,m}^{12}\rangle+3\langle\Gamma_{2,m}^{13}\rangle\\
&+\kappa\langle\Gamma_{2,m}^{22}\rangle-\kappa\langle\Gamma_{2,m}^{24}\rangle
-\langle\Gamma_{2,m}^{33}\rangle+\sqrt{2}\sqrt{\kappa+5}\langle\Gamma_{2,m}^{34}\rangle
-\kappa\langle\Gamma_{2,m}^{44}\rangle+1\big].
\end{aligned}\tag{11e}
$$

$$
\begin{aligned}
y_{2,m}^{BD}=\frac{1}{100}\big[&
4\langle\Gamma_{1,1}^{1}\rangle-4\sqrt{5-2\kappa}\langle\Gamma_{1,1}^{2}\rangle
+4\langle\Gamma_{1,1}^{3}\rangle-4\sqrt{2\kappa+5}\langle\Gamma_{1,1}^{4}\rangle\\
&-4\langle\Gamma_{2,m}^{11}\rangle
+\sqrt{2}\sqrt{\kappa+5}\{-4\langle\Gamma_{2,m}^{12}\rangle+(\kappa-3)\langle\Gamma_{2,m}^{14}\rangle+(\kappa+1)\langle\Gamma_{2,m}^{23}\rangle\}\\
&+4\{3\langle\Gamma_{2,m}^{13}\rangle-\kappa\langle\Gamma_{2,m}^{22}\rangle+\kappa\langle\Gamma_{2,m}^{24}\rangle
-\langle\Gamma_{2,m}^{33}\rangle-\sqrt{10-2\kappa}\langle\Gamma_{2,m}^{34}\rangle+\kappa\langle\Gamma_{2,m}^{44}\rangle\}+4\big].
\end{aligned}\tag{11f}
$$

$$
\begin{aligned}
y_{2,m}^{BE}=\frac{1}{50}\big[&
-2(\kappa-1)\langle\Gamma_{1,1}^{1}\rangle+2(\kappa+1)\langle\Gamma_{1,1}^{3}\rangle
-(\kappa-3)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle\\
&-(\kappa+5)\langle\Gamma_{2,m}^{22}\rangle-4\kappa\langle\Gamma_{2,m}^{24}\rangle
+(\kappa+3)\langle\Gamma_{2,m}^{33}\rangle+(\kappa-5)\langle\Gamma_{2,m}^{44}\rangle+2\big].
\end{aligned}\tag{11g}
$$

$$
\begin{aligned}
y_{2,m}^{CD}=\frac{1}{50}\big[&
2(\kappa+1)\langle\Gamma_{1,1}^{1}\rangle-2(\kappa-1)\langle\Gamma_{1,1}^{3}\rangle
+(\kappa+3)\langle\Gamma_{2,m}^{11}\rangle-4\langle\Gamma_{2,m}^{13}\rangle\\
&+(\kappa-5)\langle\Gamma_{2,m}^{22}\rangle+4\kappa\langle\Gamma_{2,m}^{24}\rangle
-(\kappa-3)\langle\Gamma_{2,m}^{33}\rangle-(\kappa+5)\langle\Gamma_{2,m}^{44}\rangle+2\big].
\end{aligned}\tag{11h}
$$

$$
\begin{aligned}
y_{2,m}^{CE}=\frac{1}{100}\big[&
4\langle\Gamma_{1,1}^{1}\rangle+4\sqrt{5-2\kappa}\langle\Gamma_{1,1}^{2}\rangle
+4\langle\Gamma_{1,1}^{3}\rangle+4\sqrt{2\kappa+5}\langle\Gamma_{1,1}^{4}\rangle\\
&-4\langle\Gamma_{2,m}^{11}\rangle
-\sqrt{2}\sqrt{\kappa+5}\{-4\langle\Gamma_{2,m}^{12}\rangle+(\kappa-3)\langle\Gamma_{2,m}^{14}\rangle+(\kappa+1)\langle\Gamma_{2,m}^{23}\rangle\}\\
&+4\{3\langle\Gamma_{2,m}^{13}\rangle-\kappa\langle\Gamma_{2,m}^{22}\rangle+\kappa\langle\Gamma_{2,m}^{24}\rangle
-\langle\Gamma_{2,m}^{33}\rangle+\sqrt{10-2\kappa}\langle\Gamma_{2,m}^{34}\rangle+\kappa\langle\Gamma_{2,m}^{44}\rangle\}+4\big].
\end{aligned}\tag{11i}
$$

$$
\begin{aligned}
y_{2,m}^{DE}=\frac{1}{25}\big[&
\langle\Gamma_{1,1}^{1}\rangle+\sqrt{2\kappa+5}(\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{2,m}^{14}\rangle)
+\sqrt{5-2\kappa}(\langle\Gamma_{2,m}^{23}\rangle-\langle\Gamma_{1,1}^{4}\rangle)\\
&+\langle\Gamma_{1,1}^{3}\rangle-\langle\Gamma_{2,m}^{11}\rangle
+\sqrt{10-2\kappa}\langle\Gamma_{2,m}^{12}\rangle+3\langle\Gamma_{2,m}^{13}\rangle\\
&+\kappa\langle\Gamma_{2,m}^{22}\rangle-\kappa\langle\Gamma_{2,m}^{24}\rangle
-\langle\Gamma_{2,m}^{33}\rangle-\sqrt{2}\sqrt{\kappa+5}\langle\Gamma_{2,m}^{34}\rangle
-\kappa\langle\Gamma_{2,m}^{44}\rangle+1\big].
\end{aligned}\tag{11j}
$$

</details>

把式 (11) 代入式 (10)，就可以使用由带 ECI 的 Monte Carlo 模拟生成的点关联函数和对关联函数，计算全部 10 个不同的化学 SRO 参数 $\alpha_{2,m}^{pq}$。每个团簇关联函数的数值由 Hamiltonian 的团簇展开与随温度变化的 Monte Carlo 模拟相结合得到。等原子 Mo-Nb-Ta-V-W HEA 的 Warren-Cowley SRO 参数 $\alpha_{2,m}^{AB}$、$\alpha_{2,m}^{AC}$、$\ldots$、$\alpha_{2,m}^{ED}$ 随合金温度的变化将在下一节分析。

## 六、结果与讨论部分完整翻译

论文没有单独使用 “Results” 或 “Discussion” 标题。对应内容为第 4 节 **Application High Entropy Alloy Mo-Nb-Ta-V-W** 和第 5 节 **Second Nearest-Neighbor SRO Effects**，位于 PDF 第 397–401 页。以下完整翻译其正文、公式与图题。

### 4. 在高熵合金 Mo-Nb-Ta-V-W 中的应用

使用为五元 Mo-Nb-Ta-V-W 体系计算的有效团簇相互作用（ECI）[11]，通过准正则 Monte Carlo 模拟研究能量上有利的原子构型如何随温度和合金组成变化。图 1 给出等原子五元 Mo-Nb-Ta-V-W 体系以及两个四元子体系 Mo-Ta-V-W 和 Mo-Nb-V-W 的混合焓从 3000 K 开始随温度的演化。对于大型体系，这些预熔构型由随机数生成；对于较小的模拟胞，可以使用特殊准随机结构（SQS）[37]。

![图 1](/imgs/2026-09-20/fig9.png)

**图 1。** 3 种等原子 HEA 的混合焓随温度的变化：五元 Mo-Nb-Ta-V-W，以及四元 Mo-Ta-V-W 和 Mo-Nb-V-W。

混合焓温度曲线上的拐点表明，体系从固溶体相向不同非随机构型发生有序-无序相变。由图 1 可见，在低于 750 K 时，部分有序相在热力学上比等原子随机固溶体构型更稳定。对于两个子体系合金，含有 Mo-Ta 二元组合的 HEA，即 Mo-Ta-V-W，混合焓不仅显著低于不含 Ta 的 Mo-Nb-V-W，也低于五元 Mo-Nb-Ta-V-W 体系。

该合金体系在 $T=400\ \mathrm{K}$ 时生成的模拟结构如图 2 所示。从所给构型可以清楚看到 Nb，绿色，在原子胞边缘发生相分离或聚集，而 V，黄色，在胞中心聚集。这意味着 Nb 与 V 原子之间不存在化学吸引相互作用。该结果与作者对 bcc Nb-V 合金的 DFT 计算一致；这些计算表明，该二元体系在整个组成范围内的混合焓均为正。相反，在立方胞左侧面可以看到 Mo-Ta，红色和蓝色，的化学有序，这与 bcc 二元 Mo-Ta 体系强烈负的混合焓完全一致。随后对化学 SRO 参数的详细分析将证实上述观察。

![图 2](/imgs/2026-09-20/fig10.png)

**图 2。** 当前 MC 模拟在 400 K 获得的等原子 Mo-Nb-V-Ta-W HEA 原子构型。原图为彩色在线图。

对于五组元合金，可以应用上一节提出的 SRO 形式体系，研究高熵 Mo-Nb-Ta-V-W 体系的有序-无序趋势如何随组成和温度变化。等原子组成下，第一近邻（1NN）壳层中 10 个不同 SRO 参数随温度的演化见图 3。

![图 3](/imgs/2026-09-20/fig11.png)

**图 3。** 等原子 HEA Mo-Nb-Ta-V-W 中 1NN 短程有序参数随温度的演化。

在高温极限，图 3 清楚显示所有 SRO 参数都趋于 0；这对应理想随机单相固溶体，其混合构型熵为：

$$
\Delta S_{\mathrm{mix}}=-R\sum_p c_p\ln(c_p),
$$

其中 $R$ 为气体常数。该表达式通常用于 HEA 的定义 [3]。在低于 750 K 的温区，图 3 表明 Mo-Ta 对的 SRO 参数 $\alpha_{2,1}^{Mo-Ta}$ 变成最负的参数。这说明，在所考察等原子 bcc Mo-Nb-Ta-V-W 合金中，Mo 原子第一壳层周围出现 Ta 原子的概率很高。Mo-Ta 对具有强烈化学 SRO 的预测，与此前对等原子四元 Mo-Nb-Ta-W 体系中第一近邻 Mo-Ta 键相互作用的 DFT 研究 [35] 一致，也与本文 DFT 数据库中 bcc 二元 Mo-Ta 体系的负混合焓一致。

MC 模拟预测的下一个强负 SRO 参数对应 V-W 对，其后是 Mo-Nb 对。需要注意，Mo-Ta、V-W 和 Mo-Nb 的负 SRO 参数与元素周期表第 V 族和第 VI 族之间 bcc 二元体系的能量相稳定性趋势一致 [38]。图 3 还表明，Mo-Nb 对的 SRO 参数在极低温下转为正，说明在五组元合金中，Nb 和 Ta 原子竞争占据 Mo 原子周围的第一近邻壳层。在低温区，$\alpha_{2,1}^{NbV}$、$\alpha_{2,1}^{MoW}$ 和 $\alpha_{2,1}^{TaV}$ 为正，表明同属第 V 族或第 VI 族元素之间的这些原子对存在强烈偏析趋势。

还应指出，多组元 bcc HEA 中第一近邻 SRO 参数的行为远比相关二元体系的预测复杂。例如，根据图 3，在极低温下，V-W 具有负 $\alpha_{2,1}^{VW}$，表示有利的 1NN 化学键合；Ta-W 则具有正 $\alpha_{2,1}^{TaW}$，表示不利的 1NN 化学键合。这似乎与 DFT 对等原子 V-W 和 Ta-W 二元体系基态结构的预测不同；前者为 $B32$，后者为 $B2_3$ [39]。已知 $B32$ 结构中的 1NN 环境对化学相互作用并不像 $B2$ 或 $B2_3$ 结构那样完全有利。

还需再次强调，本研究给出的原子构型来自一组 ECI；这些 ECI 不仅包含对相互作用，也包含式 (9) 中的三体有效团簇展开贡献。这意味着，五组元 Mo-Nb-V-Ta-W 的混合焓，以及由此得到的构型熵和混合自由能，显然超出了由 10 个组元二元体系第一近邻相互作用参数构成的简单 Ising 模型。图 3 显示，$\alpha_{2,1}^{TaW}$、$\alpha_{2,1}^{NbW}$ 和 $\alpha_{2,1}^{MoV}$ 为正；然而，对应二元体系的混合焓为负，只是其绝对值小于占主导地位的 Mo-Ta 二元合金。10 个 SRO 参数随温度变化的复杂耦合行为，可以用稳定多组元 HEA 的团簇相互作用多体效应来解释。

为了交叉检验五组元体系 SRO 表达式的一致性，尤其是五元 Mo-Nb-Ta-V-W 中第 V 族过渡金属元素 Nb 与 Ta 的竞争，作者计算了等原子四元子体系 Mo-Ta-V-W 和 Mo-Nb-V-W 的 SRO 参数，分别见图 4 和图 5。四元体系 A-B-C-D 的 SRO 参数推导方式与五组元体系相似：对关联函数公式求逆以获得平均对概率，见式 (10)。四组元体系中 6 个参数 $\alpha_{2,m}^{AB}$、$\alpha_{2,m}^{AC}$、$\ldots$、$\alpha_{2,m}^{CD}$ 的显式结果如下。

<details>
<summary>展开查看式 (12a)–(12f)：四元体系的 6 个 SRO 参数</summary>

$$
\alpha_{2,m}^{AB}=1-
\frac{-2(\langle\Gamma_{1,1}^{1}\rangle+\langle\Gamma_{1,1}^{2}\rangle-2\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{13}\rangle-\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{2,m}^{33}\rangle+1}
{(-2\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{3}\rangle+1)(-2\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{1,1}^{3}\rangle+1)}.
\tag{12a}
$$

$$
\alpha_{2,m}^{AC}=1-
\frac{-2\langle\Gamma_{1,1}^{3}\rangle-4\langle\Gamma_{2,m}^{11}\rangle+\langle\Gamma_{2,m}^{33}\rangle+1}
{(1-\langle\Gamma_{1,1}^{3}\rangle)^2-4(\langle\Gamma_{1,1}^{1}\rangle)^2}.
\tag{12b}
$$

$$
\alpha_{2,m}^{AD}=1-
\frac{-2(\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{2}\rangle+2\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{13}\rangle+\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{2,m}^{33}\rangle+1}
{(-2\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{3}\rangle+1)(2\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{1,1}^{3}\rangle+1)}.
\tag{12c}
$$

$$
\alpha_{2,m}^{BC}=1-
\frac{2(\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{2}\rangle-2\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{13}\rangle+\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{2,m}^{33}\rangle+1}
{(2\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{3}\rangle+1)(-2\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{1,1}^{3}\rangle+1)}.
\tag{12d}
$$

$$
\alpha_{2,m}^{BD}=1-
\frac{2\langle\Gamma_{1,1}^{3}\rangle-4\langle\Gamma_{2,m}^{22}\rangle+\langle\Gamma_{2,m}^{33}\rangle+1}
{(\langle\Gamma_{1,1}^{3}\rangle+1)^2-4(\langle\Gamma_{1,1}^{2}\rangle)^2}.
\tag{12e}
$$

$$
\alpha_{2,m}^{CD}=1-
\frac{2(\langle\Gamma_{1,1}^{1}\rangle+\langle\Gamma_{1,1}^{2}\rangle+2\langle\Gamma_{2,m}^{12}\rangle+\langle\Gamma_{2,m}^{13}\rangle-\langle\Gamma_{2,m}^{23}\rangle)-\langle\Gamma_{2,m}^{33}\rangle+1}
{(2\langle\Gamma_{1,1}^{1}\rangle-\langle\Gamma_{1,1}^{3}\rangle+1)(2\langle\Gamma_{1,1}^{2}\rangle+\langle\Gamma_{1,1}^{3}\rangle+1)}.
\tag{12f}
$$

</details>

利用式 (12)，等原子 Mo-Ta-V-W 子体系中 6 个 SRO 参数的演化见图 4。结果表明，该四元子体系中最有利的化学键合发生在 Mo-Ta 和 V-W 对之间，与等原子 Mo-Nb-Ta-V-W 中相应 SRO 参数的结果一致。四元 Mo-V-Ta-W 合金在低温区 $T<750\ \mathrm{K}$ 时，Mo 与 W 以及 Ta 与 V 之间的相分离倾向，也由相应参数 $\alpha_{2,1}^{MoW}$ 和 $\alpha_{2,1}^{TaV}$ 的正值体现。低温下，$\alpha_{2,1}^{TaW}$ 和 $\alpha_{2,1}^{MoV}$ 也转为正，这与前面对五元 Mo-Nb-V-Ta-W 合金的讨论一致。重要的是，与 V-W 和 Ta-W 二元体系各自 $B32$ 与 $B2_3$ 基态结构的行为相比，1NN SRO 参数 $\alpha_{2,1}^{VW}$ 和 $\alpha_{2,1}^{TaW}$ 所表现出的相反特征，在四元 Mo-Ta-V-W 中仍然存在。

![图 4](/imgs/2026-09-20/fig12.png)

**图 4。** 不含 Nb 的四元子体系中 1NN 短程有序参数随温度的变化。

在不含 Ta 时，四元 Mo-Nb-V-W 子体系的 6 个 SRO 参数见图 5。图中最值得注意的结果是，Mo 与 Nb 之间的化学 SRO $\alpha_{2,1}^{MoNb}$ 现在占据主导的负值，并与 $\alpha_{2,1}^{VW}$ 强烈竞争，取代了 Mo-Ta 对的作用。该结果说明 Mo-Nb 对优选有序；这一行为不同于五元 Mo-Nb-V-Ta-W 和四元 Mo-Ta-V-W 合金，因为在后两种体系中，Mo-Nb 对在低于 200 K 时转为正。另一方面，V-W 对的 SRO 仍为负，与此前考察的五元和四元 HEA 情形相同。最后，Nb-V 对的 SRO 参数 $\alpha_{2,1}^{NbV}$ 在四元 Mo-Nb-V-W 体系中保持显著为正，与它在初始五组元合金中的行为一致。

![图 5](/imgs/2026-09-20/fig13.png)

**图 5。** 不含 Ta 的四元子体系中 1NN 短程有序参数随温度的变化。

### 5. 第二近邻 SRO 效应

图 3、图 4 和图 5 表示所考察 3 种 HEA 第一近邻壳层中的 SRO 参数。对于 bcc 合金，第一和第二近邻距离彼此非常接近，因此与合并配位壳层 $z_1+z_2$ 对应的平均 SRO 参数定义为：

$$
\alpha_{1+2}^{pq}=
\frac{z_1\alpha_{2,1}^{pq}+z_2\alpha_{2,2}^{pq}}{z_1+z_2},
\tag{13}
$$

其中 $z_1=8$，$z_2=6$。需要注意，平均 SRO 参数 $\alpha_{1+2}^{pq}$ 可以由中子漫散射测量以较高置信度在实验上获得 [40]。图 6、图 7 和图 8 分别给出等原子 HEA Mo-Nb-Ta-V-W、Mo-Ta-V-W 和 Mo-Nb-V-W 的全部平均 SRO 参数随温度的变化。

![图 6](/imgs/2026-09-20/fig14.png)

**图 6。** 五元体系平均短程有序参数随温度的变化。

把图 6 中同时包括 1NN 和 2NN 的平均 SRO 参数 $\alpha_{1+2}^{pq}$ 与图 3 中只包括第一壳层的结果比较，可以看到五元 HEA 中 Mo-Ta 对的 SRO 参数发生显著变化。Mo-Ta 对的第二近邻 SRO 贡献转为正，所得平均 SRO 参数 $\alpha_{1+2}^{MoTa}$ 在低温区的绝对值相对于 1NN SRO 参数 $\alpha_{2,1}^{MoTa}$ 减小。这证实，在等原子 Mo-Nb-Ta-V-W HEA 中，从短程有序角度看，$B2$ 相的存在是非常有利的；在该相中，2NN 壳层只包含相同化学物种的原子，即 Mo-Mo 或 Ta-Ta 对。不含 Nb 的四元合金中也一致预测到 Mo-Ta 平均 SRO 参数的类似减小，见图 7。

![图 7](/imgs/2026-09-20/fig15.png)

**图 7。** 不含 Nb 的四元子体系平均短程有序参数随温度的变化。

比较图 6 与图 3 还得到另一个有趣结果：考虑 2NN 贡献后，V-W 对的平均 SRO 参数变得更加显著，并与 Mo-W 的相应数值可比。分别分析 1NN 与 2NN 贡献可知，V-W 对的 2NN SRO 参数 $\alpha_{2,2}^{VW}$ 在很宽温度范围内为负，仅在低于 100 K 的极低温下转为正。因此，与 Mo-Ta 二元体系不同，所考察五元 HEA 中 V-W bcc 二元体系形成 $B2$-like 相，只在低温相变时有利。

在高于 100 K 时，负 SRO 参数表明 2NN 中存在异类 V-W 原子对。这一结果似乎支持有序 $B32$ 相的局域化学环境；DFT 已预测它是等原子二元体系的基态结构 [39]。对不含 Nb 和不含 Ta 的两个四元子体系中 V-W 对平均 SRO 参数进行类似分析，分别见图 7 和图 8，结果表明从有序型 $B32$ 构型向 $B2$-like 结构的相变发生在低于 200 K 的温度。

![图 8](/imgs/2026-09-20/fig16.png)

**图 8。** 不含 Ta 的四元子体系平均短程有序参数随温度的变化。

最后，在纳入 2NN 效应后，Ta-V 对的平均 SRO $\alpha_{1+2}^{TaV}$，见图 7，和 Nb-V 对的平均 SRO $\alpha_{1+2}^{NbV}$，见图 8，在四元 HEA Mo-Ta-V-W 和 Mo-Nb-V-W 的极低温区分别明确为正。这两个 SRO 参数在五元 Mo-Nb-Ta-V-W 体系中仍保持显著为正，见图 6，显示 Ta 与 V 之间以及 Nb 与 V 之间的偏析趋势。该预测趋势与图 2 所示原子构型一致。

## 七、结论完整翻译

**原文位置：第 6 节 Conclusion，PDF 第 401–402 页。**

本文系统建立了用于五组元和四组元合金体系的 Warren-Cowley 短程有序参数理论表述。其解析表达式用关联函数写出；这些关联函数可以由多体团簇展开 Hamiltonian 与使用 ECI 的 Monte Carlo 模拟相结合来计算。作为该形式体系的应用，作者一致分析了 bcc 晶格中等原子五元 Mo-Nb-Ta-V-W 及其两个四元子体系 Mo-Nb-V-W 和 Mo-Ta-V-W 的 SRO 温度依赖。

结果表明，在低于 750 K 的温度下，体系从理想固溶体的高温极限向不同化学物种的非随机分布发生有序-无序相变；其中，构型熵的过剩项来自 SRO 参数的非零值。一般而言，当 HEA 结构优选某一特定原子对成为第一近邻时，就存在 SRO。

本研究预测，五元 Mo-Nb-Ta-V-W 和四元 Mo-Ta-V-W HEA 中 Mo-Ta 对具有强负 SRO 参数。该现象起源于负混合焓，并由此对应作者为 CE Hamiltonian 建立的 DFT 数据库中 bcc 二元 Mo-Ta 体系的有效第一近邻相互作用。该结果也与此前关于无序 $A2$ 相向有序 $B2$ 相转变的从头算分析 [35] 以及四元 Mo-Nb-Ta-W 体系的电子结构计算 [41] 一致。

对于不含 Ta 原子的四元 Mo-Nb-V-W HEA，Mo-Nb 对的化学优选在第一近邻中占据主导，因此其 SRO 在低温下显著为负。对于所考察的全部 3 种 HEA，W-V 对的化学 SRO 参数也被预测为显著负值，并且在 bcc 晶格第一近邻中始终受到偏好。需要指出，3 个 SRO 参数 $\alpha_{2,1}^{MoTa}$、$\alpha_{2,1}^{MoNb}$ 和 $\alpha_{2,1}^{WV}$ 表征元素周期表第 V 族与第 VI 族过渡金属之间的化学键合；它们的负值可以通过相应二元体系的混合焓分析来理解。

等原子五元和四元体系中其他原子对的解析 Warren-Cowley SRO 参数在低温下为正，表明第一近邻之间存在不同程度的相分离。$\alpha_{2,1}^{MoW}$、$\alpha_{2,1}^{NbTa}$、$\alpha_{2,1}^{NbV}$ 和 $\alpha_{2,1}^{TaV}$ 的正值，可以用以下事实解释：对于同属第 V 族或第 VI 族的过渡金属元素，bcc 二元体系中的化学键合在热力学上不利。然而，HEA 中 Nb-W、Ta-W 和 Mo-V 原子对之间的相分离并不直接对应其二元体系的混合焓值。多组元合金中 SRO 参数的复杂趋势，只能由超越本研究最近邻成对近似的团簇展开 Hamiltonian 模型解释。

在 bcc 合金 SRO 参数研究中进一步加入 2NN 壳层，证实五元和四元体系中的 Mo-Ta 二元组合形成有利的 $B2$ 相。对 V-W 原子对的研究表明，其平均 SRO 参数变得与 Mo-Ta 对同样显著地为负。对于所考察 HEA 中的 V-W，只有在极低温下形成 1NN 和 2NN SRO 参数符号相反的 $B2$-like 构型：五元体系中低于 100 K，四元体系中低于 200 K。在更高温度下，V-W 对的 1NN 和 2NN SRO 参数都为负，表明存在与此前对相应 bcc 二元体系预测的基态 $B32$ 相相似的局域环境。1NN 与 2NN 对的平均 SRO 参数还明确证实，bcc 过渡金属系列第 V 族元素 Nb 与 V，或 Ta 与 V，之间存在偏析趋势。

最后需要指出，本文考察的 SRO 现象用化学占位描述局域尺度相对于平均状态的偏离。在当前基于刚性晶格的 Monte Carlo 模拟中，没有适当地考虑升温时离位展开效应对原子构型的影响。混合 Monte Carlo/分子动力学模拟可以通过偏径分布函数研究 SRO，从而帮助克服这一问题 [35]。然而，对于大原子胞，从头算分子动力学技术存在计算限制；对于经典分子动力学，多组元合金体系仍很难获得准确的原子间势 [42]。

尽管本文没有考虑由原子尺寸效应导致的离位偏移，广义的第一性原理 SRO 理论可以把它纳入考虑 [11,43,44]。特别是，作者此前的研究 [11] 已表明，由 Mo-Nb-Ta-V-W 合金中不同原子物种相互作用引起的晶格畸变，其预测原子坐标与完全弛豫 DFT 计算相比，平均误差可以达到 1–2 pm；这一结论覆盖所考察的全部组成和温度。作者认为，非常有必要把本文 Mo-Nb-Ta-V-W 的 SRO 理论数据与可由实验探测的局域化学有序详细分析进行比较，例如与参考文献 [45] 首次进行的中子衍射实验比较。

## 八、结果与讨论章节子标题对照

- **4 Application High Entropy Alloy Mo-Nb-Ta-V-W：** 4. 在高熵合金 Mo-Nb-Ta-V-W 中的应用
- **5 Second Nearest-Neighbor SRO Effects：** 5. 第二近邻 SRO 效应

## 九、计算细节汇总

| 项目 | 论文给出的具体信息 | 复现时需要注意的缺失信息 |
|---|---|---|
| 第一性原理程序 | VASP | 论文未给出 VASP 版本 |
| 电子结构方法 | PAW 方法；交换关联泛函采用 GGA-PBE | 未给出平面波截断能、$k$ 点网格及电子收敛阈值 |
| 价电子处理 | V、Nb、Ta 各取 11 个价电子；Mo、W 各取 12 个价电子 | 未列出所用 PAW 数据集的具体文件名或版本 |
| 结构数据库 | 428 个基于 bcc 晶格的结构，覆盖二元、三元和四元构型 | 未逐项公开 428 个结构及其 DFT 总能 |
| 团簇展开拟合 | 使用 ATAT；以结构反演法（SIM）拟合 5 个点团簇、30 个成对团簇和 40 个三体团簇的 ECI | 未提供 ATAT 输入文件、拟合脚本及完整训练集 |
| 拟合误差 | 交叉验证误差约为 $8\ \mathrm{meV/atom}$ | 未给出各结构的逐点残差 |
| Monte Carlo 方法 | 准正则 Monte Carlo；从 $3000\ \mathrm{K}$ 的随机占位构型开始降温 | 未给出超胞尺寸、每个温度的平衡步数、采样步数、温度步长及随机种子 |
| 初始构型 | 大超胞采用预熔化随机构型；较小超胞可由 SQS 作为初始构型 | 文中未明确报告本研究每个体系实际采用的原子数和 SQS 规格 |
| 构型示例 | 图 2 展示 $400\ \mathrm{K}$ 下等原子五元体系的 Monte Carlo 原子构型 | 未公开该构型文件 |
| SRO 输出 | 五元体系计算 10 种不同原子对；四元体系各计算 6 种不同原子对 | 原始温度-SRO 数值表未公开 |
| Warren-Cowley 符号 | $\alpha^{pq}<0$ 表示异类 $p$-$q$ 近邻富集，即化学有序；$\alpha^{pq}>0$ 表示相同元素环境或相分离倾向 | 需要结合近邻壳层解释，不能仅用符号判定具体有序结构 |
| bcc 近邻数 | 第一近邻配位数 $z_1=8$，第二近邻配位数 $z_2=6$ | 平均参数 $\alpha_{1+2}^{pq}$ 会掩盖 1NN 与 2NN 的符号差异 |
| 主要温标 | 约低于 $750\ \mathrm{K}$ 出现明显有序-无序变化；V-W 的 2NN 相关性在五元体系低于约 $100\ \mathrm{K}$、四元体系低于约 $200\ \mathrm{K}$ 时改变特征 | 这些温度来自 Monte Carlo 曲线及论文文字描述，论文没有给出统计误差条或有限尺寸分析 |

计算流程可以严格概括为：先用 DFT 为固定 bcc 晶格上的多组元占位构型建立能量数据库，再用多体团簇展开把构型能表示为 ECI 与关联函数的线性组合；随后以该 CE Hamiltonian 驱动准正则 Monte Carlo 采样，计算点关联函数和成对关联函数的热平均；最后用论文推导的五元与四元解析式，把这些关联函数转换为各元素对的 Warren-Cowley SRO 参数。论文未公开足够的输入文件和采样参数，因此可重建其理论路线和公式，但仅凭正文不能逐数值复现全部曲线。

## 十、逐图与表解读

### 表 1：五元团簇展开的 ECI

表 1 列出用于五元体系的点团簇、成对团簇和三体团簇 ECI。该表是 Monte Carlo Hamiltonian 的数值核心：点项控制各组元的单点能量基准，成对项刻画不同距离和装饰下的有效二体相互作用，三体项补充二体模型无法表达的多元协同效应。作者强调，最终采用 5 个点团簇、30 个成对团簇和 40 个三体团簇，并以约 $8\ \mathrm{meV/atom}$ 的交叉验证误差衡量预测能力。表中 ECI 不能直接逐项等同于某一具体元素对的键能，因为五元点函数基底中的元素占位信息通过多个装饰后的关联函数共同编码。

### 图 1：不同组元数下的三体团簇装饰

图 1 用最近邻等边三角形说明，组元数从二元增加到五元后，可区分的三体团簇装饰数量迅速增长。该图对应式 (1) 至式 (3) 中的装饰索引：二元体系只有一种非平凡点函数，五元体系则需要 4 种非平凡点函数及其组合。图的作用是解释为何高组元 CE 的参数空间和关联函数转换关系显著复杂于二元合金。

### 图 2：$400\ \mathrm{K}$ 下五元合金的 Monte Carlo 构型

图 2 是等原子 Mo-Nb-Ta-V-W 在 $400\ \mathrm{K}$ 下的代表性原子占位快照。作者据此指出体系不是完全随机固溶体：某些异类原子对优先相邻，另一些元素之间出现回避。该图提供的是局域化学分布的直观证据；具体强弱仍由图 3 和图 6 中的 SRO 参数定量给出，不能仅凭视觉聚集判断长程有序或宏观分相。

### 图 3：五元体系的 1NN SRO

图 3 同时比较五元合金 10 种原子对的第一近邻 SRO。低温下 Mo-Ta 和 V-W 显著为负，说明这两类异元素第一近邻受到偏好；Mo-W、Nb-Ta、Nb-V 和 Ta-V 等参数为正，说明相应元素之间存在回避或偏析倾向。所有曲线在高温侧趋向零，对应随机固溶体极限。约 $750\ \mathrm{K}$ 以下的系统性偏离说明局域化学有序开始明显增强，但论文没有把该温度作为经过临界标度确定的严格相变温度。

### 图 4：不含 Nb 的 Mo-Ta-V-W 的 1NN SRO

图 4 表明移除 Nb 后，Mo-Ta 与 V-W 仍是最强的负 SRO 对，因此它们的化学偏好对五元体系的主要有序趋势具有稳健性。Mo-W 与 Ta-V 在低温为正，Mo-V 和 Ta-W 也在低温转正。该结果说明四元环境中的竞争相互作用能够改变某些原子对的局部偏好，而 Mo-Ta 和 V-W 的第一近邻有序倾向仍保持主导。

### 图 5：不含 Ta 的 Mo-Nb-V-W 的 1NN SRO

图 5 显示移除 Ta 后，Mo-Nb 取代 Mo-Ta 成为最显著的负 SRO 对，并与 V-W 的负 SRO 竞争。Nb-V 保持明显正值，反映这两种第 V 族元素之间的偏析倾向。该图直接说明，多组元合金中的某一元素对不能只用孤立二元体系解释：改变其余组元后，化学环境和相互作用竞争会重新排序主要的 SRO 通道。

### 图 6：五元体系的 1NN+2NN 平均 SRO

图 6 将第一和第二近邻按 $z_1=8$、$z_2=6$ 加权平均。与图 3 相比，Mo-Ta 的负值绝对量减小，原因是其 2NN 贡献为正；这种 1NN 异类富集、2NN 同类富集的组合与局域 $B2$ 环境一致。V-W 的平均负 SRO 则更加突出，说明在较宽温区内第二近邻 V-W 异类配对同样受到偏好。Ta-V 和 Nb-V 的正平均 SRO 进一步支持相应元素间的偏析趋势。

### 图 7：Mo-Ta-V-W 的 1NN+2NN 平均 SRO

图 7 对不含 Nb 的四元合金给出相同的平均方式。Mo-Ta 平均 SRO 的绝对值相较其 1NN 值减小，再次支持 Mo-Ta 局域 $B2$ 型有序。V-W 在大部分温区保持显著负平均 SRO，而只有低于约 $200\ \mathrm{K}$ 时，1NN 与 2NN 的组合才转向论文所称的 $B2$-like 特征。Ta-V 的正平均 SRO 表明其偏析倾向在加入 2NN 后仍存在。

### 图 8：Mo-Nb-V-W 的 1NN+2NN 平均 SRO

图 8 显示不含 Ta 时，Mo-Nb 与 V-W 仍是主要负 SRO 通道。与五元体系相似，V-W 在大部分温区呈现与 $B32$ 局域环境相符的 1NN 和 2NN 异类配对特征，而低于约 $200\ \mathrm{K}$ 后才出现 $B2$-like 变化。Nb-V 的正平均 SRO 清楚表明二者回避，并与作者对第 V 族元素间偏析的解释一致。

## 十一、严格性说明

- 论文只明确给出投稿、接收和在线发表日期，没有列出返修日期；本文档将其标为“未提供”，没有推定日期。
- 论文没有提供专门的数据仓库、代码仓库、ATAT 输入文件、VASP 输入文件或 Monte Carlo 原始输出；arXiv 链接属于论文预印本，不等同于数据或代码开源。
- Monte Carlo 超胞尺寸、平衡与采样步数、温度步长、随机种子等关键复现参数未在正文中给出；本文档没有补全这些参数。
- 本文档保留论文的公式编号、物理量定义和文献编号，并以 [n] 形式表示引用。参考文献表本身不属于用户指定的 5 个正文部分，因此未逐条翻译。
- 原文在元素顺序的书写上同时出现 Mo-Nb-Ta-V-W 与 Mo-Nb-V-Ta-W；两者表示同一个五元元素集合。本文档在标题和书目信息中沿用正式标题的顺序，在对应原句的翻译中保留原文语境。
- 图题、表题和逐图解释均基于论文正文及对应图表；曲线的物理解释限于作者在正文中明确给出的结论，没有从图像自行读取或补造数值。
