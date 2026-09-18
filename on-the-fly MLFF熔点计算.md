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
最终训练出336个构型的力场模型
-   液态：149个构型
-   固态：83个
-   界面弛豫：18个
-   界面训练：86个


<!--stackedit_data:
eyJoaXN0b3J5IjpbMzE3MjYyODE0LDc0ODg1OTIyNSwtMTg5Mz
k1OTc4MCwtNzQwNzkxNTY0LDEwODg1NzgxNzQsODE0ODc4ODgw
LDczMzQwNzgwMCwyMDk3NTg3NjUsLTEzMDk1ODgzMTEsOTAwNT
c5ODA1LDEyMzY2OTQ2NzUsMjA0MDI5NzYyMl19
-->