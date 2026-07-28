# vasp实现AIMD相关（ab initial MD & MLFF）

内容主要参考vasp wiki有关分子动力学计算的相关内容：[分子动力学计算](https://vasp.at/wiki/Molecular-dynamics_calculations)
## **AIMD的基本原理**
VASP 实现 AIMD 的核心思路是：**在每个分子动力学时间步中，通过自洽求解 Kohn-Sham 方程获得电子基态，再由 Hellmann-Feynman 力驱动离子运动**。适用于静态基态计算不足以描述的情形——**热无序、扩散、结构涨落、相稳定性、熔化、以及体系在真实条件下如何弛豫平衡**。具体流程主要是：
-  给定离子位置 $\{\mathbf{R}_I\}$，求解 Kohn-Sham 方程得到电子基态密度 $\rho(\mathbf{r})$。
-   通过 Hellmann-Feynman 定理（加上 Pulay 力修正，因为 VASP 使用平面波+PAW 基组）计算每个离子上的力 $\mathbf{F}_I = -\nabla_{\mathbf{R}_I} E_{\text{KS}}$。
-   用牛顿运动方程（通常 Verlet 算法）更新离子位置和速度。
-   重复以上步骤，得到离子轨迹。 


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEyOTQyODIxMTksOTY3NTMwODgxLDMyMT
k3NjU5Ml19
-->