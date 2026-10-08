# on-the-fly MLFF 熔点计算
关于2019年 [On-the-fly machine learning force field generation: Application to melting points | Phys. Rev. B](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.014105) 原文内容中关于熔点计算的讨论
熔点计算方法：固液界面平衡法
## MLFF训练
原文使用108原子固态Al、108原子液态Al和144原子的界面Al模型进行MLFF的训练
#### 1. 固体Al
固态Al初始结构先通过在0K下通过 DFT 获得最优结构，随后使用 NPT 系综在300K、0.1MPa下弛豫开启即时机器学习力场弛豫20ps，最终得到平衡的晶格常数。随后用这个初始结构进行300K和1000K下的力场训练。![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-09-16/qf6yHUABRwR8z1VB.png)

#### 2. 液体Al
液态结构是将构建的Al晶胞在2000K下弛豫20ps达到平衡，这个过程中也使用MLFF进行加速，随后进行100ps的训练，液体结构使用NPzT系综。
#### 3.界面Al
界面初始结构通过固定一半计算胞原子，在高温下熔化另一半原子构建，对原子施加简谐偏置势。使用液态和固态学习生成的力场进行加速，同样学习100ps。
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-09-18/BTbIQBjs9oAirL8B.png)
最终训练出336个构型的力场模型：
-   液态：149个构型
-   固态：83个
-   界面弛豫：18个（把构建界面结构过程中也学习了）
-   界面训练：86个
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-09-18/VIhWyb1P5CKe1URv.png)

相比原文献来说结构偏多。

## 钉扎界面法计算熔点

原文用集体密度区分固相和液相：

$$

Q=|\rho_{\mathbf q}|,\qquad

\rho_{\mathbf q}=\frac{1}{\sqrt N}\sum_{j=1}^{N}e^{-i\mathbf q\cdot\mathbf r_j},\qquad

\mathbf q=f_1\mathbf b_1+f_2\mathbf b_2+f_3\mathbf b_3.

$$

这里 $N$ 为原子数，$\mathbf r_j$ 为第 $j$ 个原子的位置，$\mathbf b_i$ 为倒格矢；Al 采用 $(f_1,f_2,f_3)=(8,0,0)$。**归一化是 $N^{-1/2}$，不是 $1/N$。** 对 512 原子结构及第一分量为 $x_j$ 的 Direct 坐标，有 $Q=|\sum_j e^{-2\pi i8x_j}|/\sqrt{512}$。

钉扎势写成 $U'=U+\kappa(Q-a)^2/2$。补充材料 Table S6 给 Al 的参数是 $\kappa=10\ \mathrm{eV}/Q^2$、$a=7.0$。$a$ 是运行前设定的目标序参量，原文取在纯固、纯液序参量之间，不是由熔点拟合得到。本项目 890 K 采用 $a=7.0$；920、940、960 K 采用 $a=7.96935194$，后者等于初始 512 原子 POSCAR 的 $Q$，并非文献值。正式对比最好统一 $a$，以减小不同固液比例的有限尺寸效应。

### 2. 纯相参考模拟

在候选熔点附近，分别对纯固相和纯液相用同一 MLFF 进行约 100 ps 的 NPT 模拟，统计平衡段的 $\langle Q\rangle_s$、$\langle Q\rangle_l$ 及每原子焓 $h_s$、$h_l$。这些是原文计算化学势差、熔化熵和 Newton 温度更新所需的数据。比较纯相与界面序参量时，须核对倒格矢、胞尺寸和归一化约定，不能把不同原子数体系的原始 $Q$ 不经处理直接相减。

如果只用多个温度的平均钉扎力零点求熔点，纯相焓不是必需输入；但这与原文同时计算 $\Delta\mu$ 和 $\Delta s$ 的 Newton 步骤不同。

### 3. 512 原子界面生产计算

拼接 256 原子固相和 256 原子液相，确保横向晶格参数 $L_x,L_y$ 匹配；界面平面在 $xy$，周期性边界下沿 z 有两个固液界面。生产阶段必须释放此前用于制备液相的原子固定约束。使用训练完成的同一 MLFF，在 NPzT 下只允许 z 轴长度变化。原文的 Al 界面模拟约 200 ps，时间步长 3 fs，压强为环境压强。

本项目 VASP 6.4.3 实际运行的主要 INCAR 设置如下；将 `T` 替换为各条轨迹的设定温度。`PSTRESS=0.001` 的单位为 kbar，对应 0.1 MPa。

```ini

IBRION = 0

NSW = 66667

POTIM = 3.0

MDALGO = 3

ISIF = 3

TEBEG = T

TEEND = T

LANGEVIN_GAMMA = 10

LANGEVIN_GAMMA_L = 10

PMASS = 1000

PSTRESS = 0.001

LATTICE_CONSTRAINTS = .FALSE. .FALSE. .TRUE.

ML_LMLFF = .TRUE.

ML_MODE = run

OFIELD_KAPPA = 10.0

OFIELD_K = 8 0 0

OFIELD_A = 7.0

NBLOCK = 10

```

此处应使用 `OFIELD_KAPPA`、`OFIELD_K`、`OFIELD_A`，**不是** `SPRING_K=10`、`SPRING_R0=7`：`SPRING_*` 是作用于 `ICONST` 所定义坐标的另一套偏置标签。本项目钉扎任务的空 `ICONST` 是预期现象。VASP 当前在线示例使用 $Q_6$，与原文的集体密度 $Q$ 不同；本机提交后应从 `OUTCAR` 确认 `External order field`、逐步 `Q=` 和 `order field energy` 均已输出。

### 4. 分析轨迹和求熔点

先检查 `XDATCAR` 的分层局域有序度，确认平衡段一直有固相、液相及两个界面。一种直接诊断是沿 z 分箱，逐帧计算每箱的 $C_k=|\sum_{j\in k}e^{-2\pi i8x_j}|/N_k$；有序区的 $C_k$ 明显高于无序区，但阈值应参考同温度纯相轨迹，不能用单帧噪声定界面位置。还需检查 $Q(t)$ 无长期漂移、温度与 z 向应力已平衡。不能只凭最终 `CONTCAR` 判断，也不能把逐步 $Q$ 当作独立样本。

`OUTCAR` 同行的 `Q=` 是瞬时值，`<Q>=` 是从开始到当前步的累计平均。若舍弃初期弛豫，必须对瞬时 `Q=` 重新平均。本项目各条钉扎轨迹实际输出 66660 步（199.98 ps）；当前分析舍弃前 16665 步（约 50 ps），对后 49995 步（约 150 ps）计算

$$ y(T)=\langle Q\rangle'_T-a.$$

令 $\Delta Q=Q_s-Q_l>0$，$\Delta\mu=\mu_s-\mu_l$，$N_a$ 为原文公式中的体系原子数。原文的关系为

$$\Delta\mu=-\kappa\,y(T)\,\frac{\Delta Q}{N_a}.
$$

因此 $y>0$ 表示固相更稳定，$y<0$ 表示液相更稳定，$y=0$ 对应熔点。钉扎轨迹在两侧温度均保留界面是正常的；应依据平均偏置力的符号，而非最后一帧的固相比例来定熔点。

原文以 $\Delta h=h_s-h_l$ 和 $\Delta s=(\Delta h-\Delta\mu)/T$，由 $\mathrm d\Delta\mu/\mathrm dT=-\Delta s$ 做 Newton 更新：$T_{\rm new}=T+\Delta\mu/\Delta s$，直至温度变化小于约 10 K。也可以在零点附近取多个温度，线性拟合 $y(T)$ 并解 $y(T_m)=0$；这种零点拟合不需要先计算 $\Delta Q$ 和 $\Delta h$，但必须检验线性区间和统计收敛。

### 5. 当前 Al 结果及下一步

当前四条钉扎轨迹采用同一初始 512 原子 POSCAR 和同一 `ML_FF`。后 150 ps 的结果为：

| 设定温度 (K) | $a$ | $\langle Q\rangle'$ | $y=\langle Q\rangle'-a$ |
| -- | -- | -- | -- |
| 890 | 7.00000000 | 7.00161698 | +0.00161698 |
| 920 | 7.96935194 | 7.96304245 | -0.00630949 |
| 940 | 7.96935194 | 7.95507996 | -0.01427198 |
| 960 | 7.96935194 | 7.94970412 | -0.01964782 |

当前可运行 `python3 Al_interface/analyze_pinning.py` 复算各温度的 $Q$；该脚本按 3 ps 分块误差加权，得到零点 896.16 K，但 3 ps 分块标准误可能偏小，最终误差还应比较 6、12、24 ps 分块及不同平衡段截取方案。四点普通最小二乘拟合得 $y(T)\simeq-3.1038\times10^{-4}(T-896.4\ \mathrm K)$，因此 **$T_m^{\rm MLFF}\approx896$ K（初步值）**。约 12 ps 时间块重采样给出的统计 95% 区间约为 886--904 K，不包含力场、有限尺寸、压强和温控的系统误差。另三条无偏置共存轨迹在 850、870、890 K 均于约 30 ps 内结晶，与熔点高于 890 K 附近的趋势相符，但不能单独给出精确熔点。

现阶段尚不能把 896 K 作为最终值：890 K 的正偏移接近统计噪声；四点采用不同的 $a$；`OSZICAR` 的长时平均动温比 `TEBEG` 普遍低约 13 K，应通过时间步长或温控对照检查，不能简单把熔点减去 13 K。建议固定同一 MLFF 及文献的 $a=7.0$，补做零点附近的 880、900、920 K 轨迹，并在关键温度进行独立重复。每条轨迹检查两相共存、$Q$ 平稳、分块误差及 z 向应力后，再拟合零点。

原文 Table I 中，**Al 的 PBE 力场熔点为 $871\pm8$ K**，经 DFT 热力学微扰修正后为 **$837\pm9$ K**；实验值约为 933 K。本项目的 896 K 是当前自训练 MLFF 的初步结果，不能视为原文力场的复现值或已经修正的第一性原理熔点。

## 参考

- [Jinnouchi、Karsai 和 Kresse, Phys. Rev. B 100, 014105 (2019)](https://doi.org/10.1103/PhysRevB.100.014105)，正文 Sec. III C、Table I；补充材料 Sec. S2、S4、Table S6。

- [VASP: Interface pinning calculations](https://vasp.at/wiki/Interface_pinning_calculations)，NPzT 和 $Q_6$ 示例；本项目使用的不是 $Q_6$。

- [VASP: SPRING_R0](https://vasp.at/wiki/SPRING_R0)，说明 `SPRING_*` 作用于 `ICONST` 中的偏置坐标。


<!--stackedit_data:
eyJoaXN0b3J5IjpbNjEyODk1ODY3LDEyODc1NzUxMjQsLTIwOT
IzOTQwNzQsLTEzMTU5MTU1NTksODUxOTcwNjgwLDc0ODg1OTIy
NSwtMTg5Mzk1OTc4MCwtNzQwNzkxNTY0LDEwODg1NzgxNzQsOD
E0ODc4ODgwLDczMzQwNzgwMCwyMDk3NTg3NjUsLTEzMDk1ODgz
MTEsOTAwNTc5ODA1LDEyMzY2OTQ2NzUsMjA0MDI5NzYyMl19
-->