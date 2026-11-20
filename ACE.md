## NbTaTiHf ACE 势函数
论文：[Ductility mechanisms in complex concentrated refractory alloys from atomistic fracture simulations | npj Computational Materials](https://www.nature.com/articles/s41524-025-01810-3)
代码仓库：[https://github.com/wenqing2344/NbTaTiHf_ACE_potential.git](https://github.com/wenqing2344/NbTaTiHf_ACE_potential.git)

公开数据集包含 **5,844 个训练构型、649 个测试构型**，合计 **6,493 个构型、263,810 个原子**，与论文 Methods 一致。训练集中，**二元固溶体占 60.8% 的构型，四元合金占 67.5% 的原子**。前者主要是 16 原子小胞，后者包含更多大超胞，并集中了高温 MD、自由表面、层错和弹性相关数据。

###  1. 数据集构型组成
#### 1.1 单质构型
构型类型 | Hf | Nb | Ta | Ti | 总计
|--|--|--|--|--|--|
BCC | 1 | 1 | 1 | 1 | 4
HCP (relaxed)  | 1 | 1 | 1 | 1 | 4
$\omega$ 相 (relaxed/unrelaxed) |  1/1 | 1/1 | 1/1 | 1/1 | 4/4
测试构型 |0|6|0|0|6
BCC 均匀应变| 120|120|120|120|480
BCC 自由表面| 6|6|6|6|24
**合计**|**130**|**136**|**130**|**130**|**526**

**均匀应变**：对于每个元素单质 **2 原子BCC晶胞**进行 xx、yy、zz、xy、xz、yz 六个方向、 形变量 -0.1-0.1 ，间隔 0.01 的计算，即每个方向 20 个计算构型。
**自由表面**：每种元素分别进行 $(100)$、$(110)$、$(111)$、$(112)$、$(113)$、$(322)$ 界面构型计算，对应的晶胞分别含 **32、36、24、24、88、88** 个原子。

#### 2. 二元体系
二元体系对各个元素通过将含 2 个原子的常规 BCC 晶胞扩成 $2\times2\times2$超胞，得到 16 个周期格点，进行**固定 16 原子 BCC 超胞的占位穷举，再按晶体对称性去重**。对于全部的对称不等价构型，分别进行**弛豫/非弛豫**计算。（该步骤构型可以通过icet等方式生成）
原文 Hf–Ta 和 Hf–Ti 的等比构型似乎存在问题。
#### 3. 等比三元体系
对 $3\times3\times3$ 常规 BCC 超胞先将 54 个格点分成三个各含 18 个位点的集合，再把三种元素分别分配给这三个集合，遍历全部六种分配，随后分别进行弛豫/非弛豫计算。`HfNbTa`、`HfNbTi`、`HfTaTi`、`NbTaTi` 各 6 个初始构型，进行弛豫/非弛豫计算。（）






<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIxMzA1NDAwMTYsODQzODk2Nzg5LDE5MT
E3OTk4NDIsMTI4MzAwNDA1NiwtNjgyMzgzODUwLDIwMzI4NTMz
NDIsLTExMjA0NzcxODYsLTEzNjY3ODg4OSwxMDc0Mzg4MDUsMj
U1MjEwMDYwLC0yMDU0MTEwNzAzLDE3NDU5OTEzNTUsLTEyOTUz
MTM1NzcsMjA0MDI5NzYyMl19
-->