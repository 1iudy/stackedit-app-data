# vasp实现AIMD相关（ab initial MD & MLFF）

内容主要参考vasp wiki有关分子动力学计算的相关内容：[分子动力学计算](https://vasp.at/wiki/Molecular-dynamics_calculations)
## **AIMD的基本原理**
VASP 实现 AIMD 的核心思路是：**在每个分子动力学时间步中，通过自洽求解 Kohn-Sham 方程获得电子基态，再由 Hellmann-Feynman 力驱动离子运动**。适用于静态基态计算不足以描述的情形——**热无序、扩散、结构涨落、相稳定性、熔化、以及体系在真实条件下如何弛豫平衡**。具体流程主要是：
-  给定离子位置 $\{\mathbf{R}_I\}$，求解 Kohn-Sham 方程得到电子基态密度 $\rho(\mathbf{r})$。
-   通过 Hellmann-Feynman 定理（加上 Pulay 力修正，因为 VASP 使用平面波+PAW 基组）计算每个离子上的力 $\mathbf{F}_I = -\nabla_{\mathbf{R}_I} E_{\text{KS}}$。
-   用牛顿运动方程（通常 Verlet 算法）更新离子位置和速度。
-   重复以上步骤，得到离子轨迹。 

计算分子动力学轨迹的过程中，除了直接通过DFT对原子间力进行计算，还可以通过**机器学习力场(MLFF)** 实现，通过对已计算的结构**插值**的方法进行训练，vasp使用MLFF可以进行即时训练，也可以采用外部的机器学习力场，可以通过**ML_MODE参数**进行选择。
即时机器学习力场生成方案如图所示
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-28/36eRy9Nt8SPrHbDZ.png)
## 机器学习力场MLFF的基本原理
MLFF的数据集包含布拉维晶格、原子坐标、DFT计算的总能量、力和应力张量，通过包含径向和角向分布信息的描述符识别每个原子周围的**局域构型**，从而使用模型学习对应的局域构型附近的力场。
### 训练结构判据
	
在输入文件中提前定义**力不确定性阈值**，当结构的力的不确定性超过该阈值时，则将该结构归入**训练集**，计算并提取其DFT信息。初始阈值由**ML_CTIFOR**（单位 eV/Å）给出，其后续如何自适应由**ML_ICRITERIA** 控制：`0` 不更新；`1` 用历史贝叶斯不确定性的均值更新；`2` 用滑动均值。推荐自动更新：**ML_ICRITERIA = 1**，此时阈值根据前几步中的**平均贝叶斯不确定性**进行更新，具体更新公式如下：
$$\texttt{ML\_CTIFOR} = \langle \text{stored Bayesian uncertainties} \rangle \times (1 + \texttt{ML\_CX})$$
力场不会每次遇到上述新结构时都进行重新训练，当遇到与已有结构差异不是太大的结构时，通常会进行**攒批(blocking)** 操作，收集待学习的数据，并在后续步骤中同时对所有结构进行学习，可以通过标签**ML_MCONF_NEW**设置学习块大小。一般根据力的贝叶斯不确定性大小进行更新操作：
-   **不确定性 > `ML_CDOUB` × `ML_CTIFOR`**：立刻采样该结构的局域参考构型，**当场做一次 DFT 并立即重建力场**。
-   **`ML_CTIFOR` < 不确定性 ≤ `ML_CDOUB` × `ML_CTIFOR`**：也做 DFT，但结构先列为"候选"，**攒够 `ML_MCONF_NEW` 个候选后**才一起并入训练集、统一更新力场。为避免采到太相似的结构，相邻候选之间还隔 `ML_NMDINT` 步。
-   **特例**：尚无任何力场时，第一个结构的所有原子都直接采样，建初始力场。
### 局部能量
VASP 的总能量可以写成居于原子能量之和：
$$ U = \sum_{i=1}^{N_a} U_i, \qquad U_i = F[\rho_i(\mathbf{r})] $$
局部能量$U_i$是局域原子密度的泛函，原子 $i$的能量只由它周围截断半径 $R_{\text{cut}}$ 内的邻域原子密度 $\rho_i(\mathbf{r})$ 决定。
$$\rho_i(\mathbf{r}) = \sum_{j=1}^{N_a} f_{\text{cut}}(r_{ij})\, g(\mathbf{r}-\mathbf{r}_{ij}), \quad r_{ij}=|\mathbf{r}_j-\mathbf{r}_i|$$
其中，$f_{\text{cut}}(r_{ij})$ 是截断函数，$r_{ij}>R_{\text{cut}}$ 时严格为 0，也就是仅考虑局域内原子的作用。$g(\mathbf{r})$是$\delta$ 函数，与SOAP相似，将$g(\mathbf{r})$经过高斯展宽后，得到：
$$g(\mathbf{r}) = \frac{1}{\sigma_{\text{atom}}\sqrt{2\pi}} \exp\!\left(-\frac{|\mathbf{r}|^2}{2\sigma_{\text{atom}}^2}\right)$$
### 描述符
由于能量$U_i$具有旋转不变性，即数值不跟随空间旋转发生改变，但是$\rho_i(\mathbf{r})$不具有旋转不变性，因此不能将$\rho_i(\mathbf{r})$直接作为模型描述符使用，因此需要构建具有旋转不变性的描述符。主要构建了径向描述符和角向描述符
**径向（两体）描述符**，对方向积分消去取向，仅保留距离r处有多少近邻原子，只含两体信息：
$$\rho_i^{(2)}(r) = \frac{1}{4\pi}\int \rho_i(r\hat{\mathbf{r}})\, d\hat{\mathbf{r}} $$
**角向（三体）描述符**，用 $\delta(\hat{\mathbf{r}}\cdot\hat{\mathbf{s}}-\cos\theta)$把两个方向的夹角锁成 $\theta$：

$$ \rho_i^{(3)}(r,s,\theta) = \iint d\hat{\mathbf{r}}\, d\hat{\mathbf{s}}\; \delta(\hat{\mathbf{r}}\cdot\hat{\mathbf{s}}-\cos\theta)\sum_{j}\sum_{k\ne j}\rho_{ik}(r\hat{\mathbf{r}})\,\rho_{ij}(s\hat{\mathbf{s}}) $$
### 基展开
为减少计算消耗，将$\rho_i(\mathbf{r})$投到一组基上——径向基 $\chi_{nl}(r)$（归一化球贝塞尔函数）乘球谐 $Y_{lm}(\hat{\mathbf{r}})$:
$$ \rho_i(\mathbf{r}) = \sum_{l=1}^{L_{\max}}\sum_{m=-l}^{l}\sum_{n=1}^{N_R^l} c_{nlm}^{i}\, \chi_{nl}(r)\, Y_{lm}(\hat{\mathbf{r}}) $$

系数 $c_{nlm}^{i} = \langle \chi_{nl} Y_{lm} \mid \rho_i \rangle$就是 $\rho_i$ 在这组基上的"坐标"，由 $L_{\max}$（`ML_LMAX2`）和径向截断 $N_R^l$ 截断。同理每个两体分布也有系数 $\rho_{ij}(\mathbf{r}) = \sum c_{nlm}^{ij}\chi_{nl}Y_{lm}$，且有线性性 $c_{nlm}^{i} = \sum_j c_{nlm}^{ij}$。 将展开代入$\rho_i^{(2)}(r)$，可得：
$$\rho_i^{(2)}(r) = \frac{1}{4\pi}\sum_{l,m,n} c_{nlm}^{i}\chi_{nl}(r)\int Y_{lm}\,d\hat{\mathbf{r}} = \frac{1}{\sqrt{4\pi}}\sum_n c_{n00}^{i}\chi_{n0}(r)$$
在此基础上加入角度信息，从而获得三体描述符
$$\rho_i^{(3)} = 2\pi \sum_l P_l(\cos\theta)\sum_{n,\nu}\chi_{nl}(r)\chi_{\nu l}(s)\sum_m \sum_{j,\,k\ne j} c_{nlm}^{ik}\, c_{\nu lm}^{ij*}$$
在许多情况下乘以角度滤波函数 $\eta$ ([ML_IAFILT2](https://vasp.at/wiki/ML_IAFILT2 "ML IAFILT2")），可以显著减少必要的基集大小，同时不损失计算准确性$\eta_{l,a}=1/(1+a[l(l+1)]^2)$（`ML_AFILT2`/`ML_IAFILT2`）
把对元素种类的二次求和改成"一种元素对全部元素求和"，使描述符数对元素数从二次降到近似线性：
$$ p_{n\nu l}^{iJ} = \sqrt{\frac{8\pi^2}{2l+1}}\sum_m c_{nlm}^{iJ}\sum_{J'} c_{\nu lm}^{iJ'} $$
### 核回归
对于需要学习的构型数据集$(\rho_{i_B},\, U_{i_B}^{\text{ref}})$，学习泛函关系$U_i = F[\rho_i(\mathbf{r})]$ ，使用核回归方法：将 $F$ 写成核的线性组合，即 $F$ 处在数据集构型所张成的空间之中：
$$U_i^{\alpha} = \sum_{i_B=1}^{N_B} w_{i_B}\, K(\mathbf{X}_i^{\alpha}, \mathbf{X}_{i_B})$$
其中核是衡量训练集中构型之间的相似度，用二体描述符和三体描述符表示
$$K(\hat{\mathbf{X}}_i, \hat{\mathbf{X}}_{i_B}) = \left[\, \beta\, \hat{\mathbf{X}}_i^{(2)}\!\cdot\hat{\mathbf{X}}_{i_B}^{(2)} + (1-\beta)\, \hat{\mathbf{X}}_i^{(3)}\!\cdot\hat{\mathbf{X}}_{i_B}^{(3)} \,\right]^{\zeta}$$
-   **点积** $\hat{\mathbf{X}}\cdot\hat{\mathbf{X}}_B$ 是旋转不变量之积的和，故核旋转不变——能量旋转不变由此保证。
-   **\(\beta\)**（`ML_W1`）权衡两体与三体相似度；正因为第五节把三体做成了"纯角向"，调节 $\beta$ 用来调节两体和三体相似度的权重。
-   **$\zeta$**（`ML_NHYP`）是核的幂次：$\zeta$ 越大，核越会在匹配度越高是产生相应；同时把括号展开后会出现描述符的高次乘积，对应**更高阶的多体相互作用**——所以 $\zeta$ 同时控制锐度与多体阶。
-   **归一化** $\hat{\mathbf{X}}=\mathbf{X}/\|\mathbf{X}_c\|$，其中 $\mathbf{X}_c=[\sqrt{\beta}\mathbf{X}^{(2)};\sqrt{1-\beta}\mathbf{X}^{(3)}]$。它除掉"邻居多少/密度高低"带来的尺度，让核只比较环境的**几何形状**，否则配位数不同的相似结构会被误判为不像。

### 贝叶斯线性回归
因为 $U=\sum_i U_i$ 且每个 $U_i$ 对 $w$ 线性，对它求导得到的力与应力**对 $w$ 仍然线性**：

$$ \mathbf{F}_I = -\frac{\partial U}{\partial \mathbf{R}_I} = -\sum_i \sum_{i_B} w_{i_B}\, \frac{\partial K(\mathbf{X}_i,\mathbf{X}_{i_B})}{\partial \mathbf{R}_I} $$

应力张量同理，是核对元胞坐标的导数。于是能量、$3N_a$ 个力、6 个应力分量可以用**同一方程组表示**

$$ \mathbf{Y} = \mathbf{\Phi}\,\mathbf{w} $$

设计矩阵 $\mathbf{\Phi}$ 按结构分块：每块首行是核（能量行），随后 $3N_a$ 行是核对原子坐标的导数（力行），末 6 行是核对元胞坐标的导数（应力行）。这就把"局部能量"无缝接到了上一节的贝叶斯线性回归：解出后验均值 $\bar{\mathbf{w}}$ 给预测，后验协方差 $\mathbf{\Sigma}$ 传播成预测的**不确定**——也就是 on-the-fly 每步用来判断"信力场还是回去算 DFT"的那个量。换句话说，**局部能量分解 + 核的线性结构，是 on-the-fly 能"每步廉价给出不确定"的根本原因**。

## AIMD & MLFF 设置
### POSCAR
需要一个包含足够大晶格的[POSCAR](https://vasp.at/wiki/POSCAR "POSCAR")。实际上，MD 通常需要相当数量的离子，这样轨迹才能采样到有意义的局域环境分布。如果晶胞太小，统计效果会很差，而且原子、缺陷或局域畸变可能会与其周期性镜像发生过于强烈的相互作用。

在进行分子动力学模拟之前可以考虑进行[结构优化](https://vasp.at/wiki/Structure_optimization "Structure optimization")。如果初始结构仍受应变或残留力较大，MD运行的初始步骤将用于消除异常应力（弛豫），而非采样所需的物理运动。这往往导致轨迹不稳定、温度控制变差，并且在 MLFF 训练工作流中产生低质量的训练数据。或者，你也可以用之前 MD 运行中的 [CONTCAR](https://vasp.at/wiki/CONTCAR "CONTCAR") 文件作为起点继续轨迹，因为除了结构本身，它已经包含了离子速度。
### INCAR
在INCAR文件中需要设置启用MD，此外，还需要定义轨迹步数，时间步长，温度设置和单元约束

**IBRION**
-  `IBRION = 0` **设置启用分子动力学模拟**

**NSW**

-  `NSW != 0` **设置MD运动的步数**。通常单点能计算中设置`NSW = 0`，但对于分子动力学必须设置最大离子步，对于每个离子步，最多可以执行  [NELM](https://vasp.at/wiki/NELM "NELM") 个电子步。如果在之前满足 [EDIFF](https://vasp.at/wiki/EDIFF "EDIFF") 设定的收敛标准，则数量会减少。力和应力根据 [ISIF](https://vasp.at/wiki/ISIF "ISIF") 对每个离子阶的设置计算。

**POTIM**
- `POTIM = [real]` **设置分子动力学的时间步长，单位为fs**。对于含有氢键或刚键的系统通常需要比重且弱结合的系统更短的时间步长。

**TEBEG、TEEND**
- `TEBEG = [real]` **设置起始温度，单位为K**。如果[POSCAR](https://vasp.at/wiki/POSCAR "POSCAR")文件中没有提供初始速度，则速度将根据初始温度TEBEG下的麦克斯韦-玻尔兹曼分布随机设定。速度仅用于分子动力学（[IBRION](https://vasp.at/wiki/IBRION "IBRION")=0）。
- `TEEND = [real]` **设置模拟终止温度，单位为K**。默认条件下 `TEEND = TEBEG`。

**MDALGO**
用于指定计算过程中的恒温器选择
- `MDALGO = 0 (Default)`
- `MDALGO = 1` **Andersen 恒温器**  可用于NVT系综和NVE系综。当用于NVE系综时，需要为[ANDERSEN_PROB](https://vasp.at/wiki/ANDERSEN_PROB "ANDERSEN PROB")设定合适的值。通常在加热到某个目标温度后进行。
- `MDALGO = 2` **Nosé-Hoover 恒温器** 仅适用于NVT系综，并且需要为SMASS设置合适的值（Nosé-Hoover 恒温器需要SMASS大于等于0）
- `MDALGO = 3` **Langevin恒温器** 适用于NVT系综、NpT系综和NpH系综。NVT系综需要为POSCAR文件中所有元素设置适当的摩擦系数（[LANGEVIN_GAMMA](https://vasp.at/wiki/LANGEVIN_GAMMA "LANGEVIN GAMMA")值），以启用 Langevin 恒温器。并且需要使用固定晶胞形状和体积（ISIF小于等于2）。NpT系综为实现晶格动力学，为晶格自由度（[LANGEVIN_GAMMA_L](https://vasp.at/wiki/LANGEVIN_GAMMA_L "LANGEVIN GAMMA L")）设置并指定一组独立的摩擦系数，以及为晶格自由度（[PMASS）](https://vasp.at/wiki/PMASS "PMASS")设定假质量。
- `MDALGO = 4` **Nosé-Hoover 链式恒温器** 相当于多个Nosé-Hoover 恒温器，使用[NHC_NCHAINS](https://vasp.at/wiki/NHC_NCHAINS "NHC NCHAINS") 设置恒温器的数量，并选择合适的恒温器参数[设置NHC_PERIOD](https://vasp.at/wiki/NHC_PERIOD "NHC PERIOD")。
-  `MDALGO = 5` **CSVR恒温器** 用于NVT系综，需要设置周期[CSVR_PERIOD](https://vasp.at/wiki/CSVR_PERIOD "CSVR PERIOD")。
-  `MDALGO = 13`**多个 Andersen 恒温器** 最多可将三个用户定义的原子子系统与独立的[Anderson 恒温器](https://vasp.at/wiki/Andersen_thermostat "Andersen thermostat")耦合

**ISIF**
用于设置应力张量计算和自由度
| ISIF | 力 | 应力张量 | 离子位置 | 晶胞形状 | 晶胞体积 | 典型用途 |
|---|---|---|---|---|---|---|
| 0 | ✔ | ✘ | ✔ | ✘ | ✘ | MD 默认；只动离子，不算应力 |
| 1 | ✔ | 仅迹 | ✔ | ✘ | ✘ | 只关心总压强 |
| 2 | ✔ | ✔ | ✔ | ✘ | ✘ | 固定晶胞的离子弛豫（最常用） |
| 3 | ✔ | ✔ | ✔ | ✔ | ✔ | 全自由度弛豫（离子+形状+体积） |
| 4 | ✔ | ✔ | ✔ | ✔ | ✘ | 固定体积下弛豫离子和形状 |
| 5 | ✔ | ✔ | ✘ | ✔ | ✘ | 只弛豫形状（离子、体积冻结） |
| 6 | ✔ | ✔ | ✘ | ✔ | ✔ | 只弛豫晶胞（离子冻结） |
| 7 | ✔ | ✔ | ✘ | ✘ | ✔ | 只弛豫体积（形状、离子冻结） |
| 8 | ✔ | ✔ | ✔ | ✘ | ✔ | 离子+体积，固定形状（VASP ≥ 6.4.1） |

**MDALGO和ISIF共同决定系综和恒温器**
| 系综 | Andersen | Nosé-Hoover | Langevin | Nosé-Hoover 链 | CSVR | 多 Andersen |
|---|---|---|---|---|---|---|
| 微正则 (NVE) | MDALGO=1, ANDERSEN_PROB=0.0 | | | | | |
| 正则 (NVT) | MDALGO=1<br>ISIF=2 | MDALGO=2<br>ISIF=2 | MDALGO=3<br>ISIF=2 | MDALGO=4<br>ISIF=2 | MDALGO=5<br>ISIF=2 | MDALGO=13<br>ISIF=2 |
| 等温等压 (NpT) | 不可用 | 不可用 | MDALGO=3<br>ISIF=3 | 不可用 | 不可用 | 不可用 |
| 等焓等压 (NpH) | | | MDALGO=3, ISIF=3, LANGEVIN_GAMMA=LANGEVIN_GAMMA_L=0.0 | | | |
根据你想采样的物理条件选择合奏。由于VASP通过结合[MDALGO](https://vasp.at/wiki/MDALGO "MDALGO")和[ISIF](https://vasp.at/wiki/ISIF "ISIF")来确定集成，你选择[的恒温](https://vasp.at/wiki/Thermostats "Thermostats")器直接影响你可用的细胞自由度。例如，虽然[朗之万恒温器](https://vasp.at/wiki/Langevin_thermostat "Langevin thermostat")（）灵活支持[NVT](https://vasp.at/wiki/NVT_ensemble "NVT ensemble")和[NpT](https://vasp.at/wiki/NpT_ensemble "NpT ensemble")模拟，但其他算法不支持单元的独立配置。`[MDALGO](https://vasp.at/wiki/MDALGO "MDALGO") = 3`

-   **正则系综NVT**利用它在固定粒子数（N）、固定体积（V）和恒定温度（T）下运行模拟。这里可以使用多种恒温器，包括[Andersen thermostat、Nosé-Hoover thermostat、Langevin thermostat、"Nosé-Hoover chain thermostat"、CSVR thermostat、以及多个安达森恒温器（）。保持该系综与[ISIF](https://vasp.at/wiki/ISIF "ISIF") < 3保持固定; 是常见的选择，因为它还报告了完整的应力张量。
-   **微正则系综NVE**仅在平衡后使用。该系综非常有用，因为原子仅由MLFF或DFT力传播。所以速度不会添加人工恒温器数据。如果对自相关函数感兴趣，这个集合可能会有帮助。例如，速度[自相关函数](https://vasp.at/wiki/Sampling_phonon_spectra_from_molecular-dynamics_simulations "Sampling phonon spectra from molecular-dynamics simulations")可能值得关注，因为可以从中获得[声子DOS](https://vasp.at/wiki/Computing_the_phonon_dispersion_and_DOS "Computing the phonon dispersion and DOS")。它被视为一种特殊情况，选择了一个恒温器但实际上关闭了。最简单的方法是用 和 。另一种选择是 ，该选项禁用 [Nosé–Hoover 恒温器](https://vasp.at/wiki/Nos%C3%A9-Hoover_thermostat "Nosé-Hoover thermostat")，并产生 [NVE](https://vasp.at/wiki/NVE_ensemble "NVE ensemble") 动态学。保持与[ISIF](https://vasp.at/wiki/ISIF "ISIF") 3的<保持联系。请注意，恒温器的选择将决定NVE仿真所采用的传播方案。
-   **等温等压系综NpT**当压力和体积波动是问题的一部分时，可以使用此方法。例如，当相变研究时，仿真箱会发生变化，这种情况就适用。在 VASP 中，这针对[朗之文动力学](https://vasp.at/wiki/Langevin_thermostat "Langevin thermostat")实现，即与 一起。[朗之文恒温](https://vasp.at/wiki/Langevin_thermostat "Langevin thermostat")器需要以下附加标签：离子用[LANGEVIN_GAMMA](https://vasp.at/wiki/LANGEVIN_GAMMA "LANGEVIN GAMMA")，晶格用[LANGEVIN_GAMMA_L](https://vasp.at/wiki/LANGEVIN_GAMMA_L "LANGEVIN GAMMA L")。[PMASS](https://vasp.at/wiki/PMASS "PMASS") 控制虚构晶格质量。`[MDALGO](https://vasp.at/wiki/MDALGO "MDALGO") = 3``[ISIF](https://vasp.at/wiki/ISIF "ISIF") = 3`
-   **[等焓–等压系（NpH）：](https://vasp.at/wiki/NpH_ensemble "NpH ensemble")**当你希望保持恒压且不触发恒温器时，可以使用此方法。该集合在研究结晶过程时颇具趣味性。结晶将势能转化为动能。一个（[NpT](https://vasp.at/wiki/NpT_ensemble "NpT ensemble")）恒温器人为消耗动能以保持温度平稳，破坏调节真实成核率的自然加热。同样，必须使用带有 和 的朗之万路径，但带有 和 ，这样可以同时关闭离子和晶格恒温器。`[MDALGO](https://vasp.at/wiki/MDALGO "MDALGO") = 3``[ISIF](https://vasp.at/wiki/ISIF "ISIF") = 3``[LANGEVIN_GAMMA](https://vasp.at/wiki/LANGEVIN_GAMMA "LANGEVIN GAMMA") = 0``[LANGEVIN_GAMMA_L](https://vasp.at/wiki/LANGEVIN_GAMMA_L "LANGEVIN GAMMA L") = 0`

对于大多数工作流程，先从[NVT](https://vasp.at/wiki/NVT_ensemble "NVT ensemble")开始，然后在[NVE](https://vasp.at/wiki/NVE_ensemble "NVE ensemble")中验证时间步和强制质量。只有在细胞波动时才使用[NpT](https://vasp.at/wiki/NpT_ensemble "NpT ensemble")。除非施加晶格约束，[NpT](https://vasp.at/wiki/NpT_ensemble "NpT ensemble")可能导致液体或长程有序有限系统的不可逆胞体变形。

**ISYM**
- `ISYM=0` **分子动力学设置ISYM=0**。VASP不使用对称性，但会假设$\Psi_{\mathbf{k}} = \Psi^*_{-\mathbf{k}}$，并相应减少布里渊区的采样。

**ML_LMFF & ML_MODE & ML_TYPE**
- **Ab initio MD**：所有力都从每一步的DFT计算出来。默认配置，不需要进行参数设置。
- **原生 on-the-fly VASP MLFF**：开启即时学习MLFF，
设置：`ML_LMFF =  0` `ML_MODE = train`
- **原生VASP MLFF，仅预测**：训练并重新调整力场后，仅用MLFF预测进行MD。该模式不生成新的从头数据，因此应仅在力场适用性确认后使用。设置：`ML_MODE = run`
- **外部预训练力场，以GRACE为例**：使用 GRACE 等外部模型进行MD预测.设置：`ML_MODE = run` `ML_LMLFF = .TRUE.` `ML_TYPE = grace`

**ML_IWEIGHT、ML_WTOTEN、ML_WITFOR、ML_WTSIF**
- 用于设置训练数据的归一化和权重，分别对应能量、力和应力

**ML_ISTART**


## MLFF文件
### ML_AB
该文件作为机器学习力场方法中的输入 ML_AB 和输出 ML_ABN 。它包含了之前计算的从头开始数据收集：布拉维矩阵、原子位置、能量、力和应力张量。
多个ML_AB文件可以手动合并，相关内容可在[ML_AB - VASP维基](https://vasp.at/wiki/ML_AB)查看
当任务中存在 ML_AB 文件时，将基于现有的结构的数据集生成力场，然后对于 POSCAR 结构进行 MD 模拟，使用 `ML_MODE = SELECT` 。不存在ML_AB文件时，将从零开始进行力场构建

### ML_FF
该文件包含机器学习的力场，该力场用于仅预测模式 `ML_MODE = run`

### ML_LOGFILE
该文件为日志文件，包含MLFF相关设置和运行结果。后处理相关内容基本以此文件为主。
#### 1. 内存使用情况
这通常是ML_LOGFILE的第一部分，包含基于启动时读取的VASP文件对内存需求的**估算**。
```
* MEMORY INFORMATION ***********************************************************************************************************************

Estimated memory consumption for ML force field generation (MB):

Persistent allocations for force field        :    516.9
|
|-- CMAT for basis                            :     20.3
|-- FMAT for basis                            :    458.5
|-- DESC for basis                            :      2.6
|-- DESC product matrix                       :      2.3

Persistent allocations for ab initio data     :      8.1
|
|-- Ab initio data                            :      7.8
|-- Ab initio data (new)                      :      0.3

Temporary allocations for sparsification      :    460.9
|
|-- SVD matrices                              :    460.7

Other temporary allocations                   :     15.5
|
|-- Descriptors                               :      4.7
|-- Regression                                :      6.5
|-- Prediction                                :      4.2

Total memory consumption                      :   1001.4

********************************************************************************************************************************************

```
#### 2. 机器学习设置
包含本次模拟中对于机器学习力场生成的 INCAR 设置内容。格式如下：
```
* MACHINE LEARNING SETTINGS ****************************************************************************************************************

This section lists the available machine-learning related settings with a short description, their
selected values and the INCAR tags. The column between the value and the INCAR tag may contain a
"state indicator" highlighting the origin of the value. Here is a list of possible indicators:

 *     : (empty) Tag was not provided in the INCAR file, a default value was chosen automatically.
 * (I) : Value was provided in the INCAR file.
 * (i) : Value was provided in the INCAR file, deprecated tag.
 * (!) : A value found in the INCAR file was overwritten by the contents of the ML_FF file.
 * (?) : The value for this tag was never set (please report this to the VASP developers).

Tag values with associated units are given here in Angstrom/eV, if not specified otherwise.

Please refer to the VASP online manual for a detailed description of available INCAR tags.


General settings
--------------------------------------------------------------------------------------------------------------------------------------------
Machine learning operation mode in strings (supertag)                                                 :         REFIT (I) ML_MODE       
Machine learning operation mode                                                                       :             4     ML_ISTART     
Precontraction of weights on Kernel for fast execution (ML_ISTART=2 only), but no error estimation    :             T     ML_LFAST      
Controls the verbosity of the output at each MD step when machine learning is used                    :             1     ML_OUTPUT_MODE
Sets the output frequency at various places for ML_ISTART=2                                           :             1     ML_OUTBLOCK

 
Descriptor settings
--------------------------------------------------------------------------------------------------------------------------------------------
Radial descriptors:
-------------------
Cutoff radius of radial descriptors                                                                   :   5.00000E+00     ML_RCUT1
Gaussian width for broadening the atomic distribution for radial descriptors                          :   5.00000E-01     ML_SION1
Number of radial basis functions for atomic distribution for radial descriptors                       :             8     ML_MRB1

Angular descriptors:
--------------------
Cutoff radius of angular descriptors                                                                  :   5.00000E+00     ML_RCUT2
Gaussian width for broadening the atomic distribution for angular descriptors                         :   5.00000E-01     ML_SION2
Number of radial basis functions for atomic distribution for angular descriptors                      :             8     ML_MRB2
Maximum angular momentum quantum number of spherical harmonics used to expand atomic distributions    :             4     ML_LMAX2
...
```
#### 3. 已有的第一性原理数据信息
总结了 ML_AB 文件中的 ab initial 数据信息
```
* AVAILABLE AB INITIO DATA *****************************************************************************************************************

Number of stored (maximum) ab initio structures:       114 (     1500)
 * System   1 :       114 , name: "Si cubic diamond 2x2x2 super cell"
 * System   2 :         0 , name: "Si cubic diamond 2x2x2 super cell"
Maximum number of atoms per element:
 * Element Si :        64

********************************************************************************************************************************************
```
#### 4. 主循环信息
根据 ML_MODE 的机器学习模式，它包含了VASP运行中所有时间步（或其他迭代方案）收集的数据。分为两部分：主循环头部的描述块解释可用数据，并以行列形式呈现其排列。然后，主循环主体包含了之前定义布局中的实际数据（主要是原始数字）。
**循环开头**
由多个块组成，每个块引入后面在循环正体中出现的行列。
```
* MAIN LOOP ********************************************************************************************************************************

# STATUS ###############################################################
# STATUS This line describes the overall status of each step.
# STATUS 
# STATUS nstep ..... MD time step or input structure counter
# STATUS state ..... One-word description of step action
# STATUS             - "accurate"  (1) : Errors are low, force field is used
# STATUS             - "threshold" (2) : Errors exceeded threshold, structure is sampled from ab initio
# STATUS             - "learning"  (3) : Stored configurations are used for training force field
# STATUS             - "critical"  (4) : Errors are high, ab initio sampling and learning is enforced
# STATUS             - "predict"   (5) : Force field is used in prediction mode only, no error checking
# STATUS is ........ Integer representation of above one-word description (integer in parenthesis)
# STATUS doabin .... Perform ab initio calculation (T/F)
# STATUS iff ....... Force field available (T/F, False after startup hints to possible convergence problems)
# STATUS nsample ... Number of steps since last reference structure collection (sample = T)
# STATUS ngenff .... Number of steps since last force field generation (genff = T)
# STATUS ###############################################################
# STATUS            nstep     state is doabin    iff   nsample    ngenff
# STATUS                2         3  4      5      6         7         8
# STATUS ###############################################################

# STDAB ####################################################################
# STDAB This line contains the standard deviation of the collected ab initio reference data.
# STDAB
# STDAB nstep ........ MD time step or input structure counter
# STDAB std_energy ... Standard deviation in energy (eV atom^-1)
# STDAB std_force .... Standard deviation in forces (eV Angst^-1)
# STDAB std_stress ... Standard deviation in stress (kB)
# STDAB ####################################################################
# STDAB             nstep       std_energy        std_force       std_stress
# STDAB                 2                3                4                5
# STDAB ####################################################################

...
```
**循环主体**
在头部之后，主环主体呈现了VASP运行中收集的信息时间序列。属于同一时间步的数据块被虚线围栏。第一行描述了结构的关键词，在下面的例子中，这行表示“学习”，即机器学习力场被重新训练。
```
...
--------------------------------------------------------------------------------
STATUS                 82 learning   3      T      T         0        72
LCONF                  82 Si      1222      1228
SPRSC                  82       129       129 Si      1228      1224
REGR                   82    1    1   1.27238822E+00   5.73175466E-02   7.83203623E-12 
REGR                   82    1    2   1.28510216E+00   5.73084508E-02   7.75332075E-12 
REGRF                  82    1    3   1.29486873E+00   5.73015362E-02   7.69391276E-12    2.23430718E+16   5.75166077E+09
STDAB                  82   1.28851006E-01   1.02791005E+00   1.07081172E+01
ERR                    82   1.21269596E-02   2.35740491E-01   4.40365370E+00
CFE                    82   2.71935242E-01   2.20681769E-01   7.30391193E-01
LASTE                  82   1.63070075E-02   2.66475855E-01   7.17595981E+00
BEE                    82   4.72039040E-05   1.03291046E-01   3.02999592E-02   9.56824349E-02   6.23077315E-01   4.66683801E-01
THRHIST                82    1   8.45535075E-02
THRHIST                82    2   8.99995395E-02
THRHIST                82    3   9.42765991E-02
THRHIST                82    4   9.37027237E-02
THRHIST                82    5   9.78682111E-02
THRHIST                82    6   1.02991465E-01
THRHIST                82    7   1.04972577E-01
THRHIST                82    8   1.02574658E-01
THRHIST                82    9   9.68150073E-02
THRHIST                82   10   8.90700596E-02
THRUPD                 82   9.54674570E-02   9.56824349E-02   6.60216623E-02   1.06906899E-02
BEEF                   82   4.58511233E-05   9.95065359E-02   2.94732909E-02   9.56824349E-02   6.03276708E-01   4.51396163E-01
--------------------------------------------------------------------------------
...
```
此外，关键词 accurate 代表预测的结构。
#### 5. 时间信息
最后一节提供了不同机器学习程序部分的时间间序（不考虑从头代码部分）。系统时钟（墙时）和CPU时间（进程所有线程的总和）各有独立列。
```
* TIMING INFORMATION ***********************************************************************************************************************

Program part                                         system clock (sec)       cpu time (sec)
---------------------------------------------------|--------------------|-------------------
Setup (file I/O, parameters,...)                   |              0.242 |              0.240
Descriptor and design matrix                       |             10.540 |             10.536
Sparsification of configurations                   |              9.183 |              9.177
Regression                                         |             14.778 |             14.770
Prediction                                         |             32.461 |             32.450
---------------------------------------------------|--------------------|-------------------
TOTAL                                              |             67.204 |             67.173

********************************************************************************************************************************************
```
## MLFF 后处理
通过使用 grep 处理 ML_LOGFILE 文件对机器学习结果进行分析。例如，**预测误差分析**：
```
grep ERR ML_LOGFILE
```
这将将主循环头部和正文的内容合并，得到以下结果：
```
# ERR ######################################################################
# ERR This line contains the RMSEs of the predictions with respect to ab initio results for the training data.
# ERR 
# ERR nstep ......... MD time step or input structure counter
# ERR rmse_energy ... RMSE of energies (eV atom^-1)
# ERR rmse_force .... RMSE of forces (eV Angst^-1)
# ERR rmse_stress ... RMSE of stress (kB)
# ERR ######################################################################
# ERR               nstep      rmse_energy       rmse_force      rmse_stress
# ERR                   2                3                4                5
# ERR ######################################################################
ERR                     2   8.77652825E-05   1.00592308E-02   2.68800480E-02
ERR                     3   3.01865279E-05   1.06283576E-02   5.81209819E-02
ERR                     4   1.52820686E-04   1.31384993E-02   1.10439716E-01
ERR                     5   1.62739008E-04   1.74252575E-02   1.40488725E-01
ERR                     6   2.97462508E-04   2.32615279E-02   1.79092561E-01
ERR                     7   2.10891509E-04   2.79123925E-02   1.94566420E-01
ERR                     8   3.26150852E-04   3.15081244E-02   1.76637577E-01
ERR                     9   7.03479132E-04   3.42249550E-02   1.66830771E-01
ERR                    10   2.41808229E-04   3.54422133E-02   1.80246157E-01
ERR                    11   2.46299647E-04   3.70102675E-02   2.01262013E-01
ERR                    12   3.57654922E-04   3.93143970E-02   2.20533745E-01
ERR                    14   1.95974374E-04   4.31813231E-02   2.44026531E-01
ERR                    15   4.94080997E-04   4.73774930E-02   2.74308998E-01
ERR                    16   9.62150633E-04   5.07005683E-02   3.17482301E-01
ERR                    18   1.31336233E-03   5.39222716E-02   3.25526268E-01
ERR                    21   1.07020831E-03   5.67663475E-02   3.04995023E-01
ERR                    24   9.88977484E-04   6.37987961E-02   3.83686143E-01
ERR                    26   9.63361971E-04   6.81972633E-02   4.92021943E-01
ERR                    29   1.81730719E-03   7.47758864E-02   6.38563225E-01
```
可以将得到的结果输出到文件中：
```
grep ERR ML_LOGFILE > err.dat
```
`grep ERR ML_LOGFILE|grep -v "#"|awk '{print $2, $4}' > ERR.dat`

`grep BEEF ML_LOGFILE|grep -v "#"|awk '{print $2, $4}' > BEEF.dat`

`grep BEEF ML_LOGFILE|grep -v "#"|awk '{print $2, $6}' > CTIFOR.dat`
除了均方误差分析，还有结构数量分析、能量分析等。
**结构数量分析**： `grep LCONF` 分析每个学习步骤的构型数量。
**贝叶斯误差分析**：`grep BEEF` 能量、力和应力的估计贝叶斯误差。
**贝叶斯误差阈值**：`grep THRUPD` 贝叶斯误差阈值参数 ML_CTIFOR 的更新记录。`grep THRHIST`  贝叶斯误差阈值参数 ML_CTIFOR 的历史记录。

## 参数优化
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE3MzY3ODMxNDQsLTE4MTA5NjI0NzgsMT
QzMTQzMDU0NiwtMTQ0NTUyMzgyMCwtMTAwNDE3NTM0Miw0OTU2
NDUyNDgsLTE2Nzg5Mzk1MzMsMTA3MDIyMjM3MiwxMjU0MjcwND
c4LDc5MzM5MTMzNCwtMzQwOTU0OTUyLDE0NDEyNDUzMzEsOTIw
NjAxMzg0LDE4NTk0MDU2NjksLTIxMjEwODk3MDMsMTgzMTM4MT
U4OCw3MTkyMDM1OTIsODMzNjU4Mzg1LC0xNTQwNTAwNjQ0LC0z
Mzk2NzA0MjBdfQ==
-->