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
#### 训练结构判据
	
在输入文件中提前定义**力不确定性阈值**，当结构的力的不确定性超过该阈值时，则将该结构归入**训练集**，计算并提取其DFT信息。初始阈值由**ML_CTIFOR**（单位 eV/Å）给出，其后续如何自适应由**ML_ICRITERIA** 控制：`0` 不更新；`1` 用历史贝叶斯不确定性的均值更新；`2` 用滑动均值。推荐自动更新：**ML_ICRITERIA = 1**。
力场不会每次遇到上述新结构时都进行重新训练，一般根据，
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE3MDE1MTI1NjQsLTY4NDQ3NDAyNiwxNz
EzNTQwNDYxLDEyNzUzOTY3ODYsLTIxMjI4MDIwMjgsLTExNjcz
Nzk1NDYsOTY3NTMwODgxLDMyMTk3NjU5Ml19
-->