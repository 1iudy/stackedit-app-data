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

应力张量同理，是核对元胞坐标的导数。于是能量、$3N_a$ 个力、6 个应力分量可以**塞进同一个线性方程组**

\[ \mathbf{Y} = \mathbf{\Phi}\,\mathbf{w} \]

设计矩阵 \(\mathbf{\Phi}\) 按结构分块：每块首行是核（能量行），随后 \(3N_a\) 行是核对原子坐标的导数（力行），末 6 行是核对元胞坐标的导数（应力行）。这就把"局部能量"无缝接到了上一节的贝叶斯线性回归：解出后验均值 \(\bar{\mathbf{w}}\) 给预测，后验协方差 \(\mathbf{\Sigma}\) 传播成预测的**不确定**——也就是 on-the-fly 每步用来判断"信力场还是回去算 DFT"的那个量。换句话说，**局部能量分解 + 核的线性结构，是 on-the-fly 能"每步廉价给出不确定"的根本原因**。
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTkxMzcxNTk3LC00Mjc1NDYzMDEsLTE1MD
M0MjIzNzQsMjc4NjQ5MDQyLDEzMzAxNDM4NTQsLTIzMzcyNDgw
MiwxMDA2MzA4MjQyLC0xOTI4OTgxMzQ0LC0xMTg0MDg2MzM0LC
0xODUwMDY4NTk3LDI1NzM0OTAxNywtMTcxNTcwMTY5MiwtNjg0
NDc0MDI2LDE3MTM1NDA0NjEsMTI3NTM5Njc4NiwtMjEyMjgwMj
AyOCwtMTE2NzM3OTU0Niw5Njc1MzA4ODEsMzIxOTc2NTkyXX0=

-->