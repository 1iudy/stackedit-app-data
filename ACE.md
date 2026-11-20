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







<!--stackedit_data:
eyJoaXN0b3J5IjpbLTExMjA0NzcxODYsLTEzNjY3ODg4OSwxMD
c0Mzg4MDUsMjU1MjEwMDYwLC0yMDU0MTEwNzAzLDE3NDU5OTEz
NTUsLTEyOTUzMTM1NzcsMjA0MDI5NzYyMl19
-->