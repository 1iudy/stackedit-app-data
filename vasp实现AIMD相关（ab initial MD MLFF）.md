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

## MLFF文件
### ML_AB
该文件作为机器学习力场方法中的输入 ML_AB 和输出 ML_ABN 。它包含了之前计算的从头开始数据收集：布拉维矩阵、原子位置、能量、力和应力张量。
多个ML_AB文件可以手动合并，相关内容可在[ML_AB - VASP维基](https://vasp.at/wiki/ML_AB)查看

### ML_FF
该文件包含机器学习的力场，该力场用于仅预测模式 `ML_MODE = run`
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTg1OTQwNTY2OSwtMjEyMTA4OTcwMywxOD
MxMzgxNTg4LDcxOTIwMzU5Miw4MzM2NTgzODUsLTE1NDA1MDA2
NDQsLTMzOTY3MDQyMCwtODQxNjg0MDU2LC0xMzg0NzI0MDIyLC
0xNjE2ODcxMzUwLC0xNTY4OTkxMTgsNTQ3MDY0NDczLDE3MjYy
NzE2NjIsNTgzNDY5MzMyLDIwNDk3MDYyNCwtOTI4MDc5ODI1LD
E5MjgxMDM3MjUsLTM4OTAyNjE5NCwtNDI3NTQ2MzAxLC0xNTAz
NDIyMzc0XX0=
-->