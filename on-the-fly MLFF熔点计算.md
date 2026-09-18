# on-the-fly MLFF 熔点计算
关于2019年 [On-the-fly machine learning force field generation: Application to melting points | Phys. Rev. B](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.014105) 原文内容中关于熔点计算的讨论
熔点计算方法：固液界面平衡法
## MLFF训练
原文使用108原子固态Al、108原子液态Al和144原子的界面Al模型进行MLFF的训练
#### 1. 固体Al
固态Al初始结构先通过在0K下通过 DFT 获得最优结构，随后使用 NPT 系综在300K、0.1MPa下弛豫开启即时机器学习力场弛豫20ps，最终得到平衡的晶格常数。随后用这个初始结构进行300K和1000K下的力场训练。![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-09-16/qf6yHUABRwR8z1VB.png)

#### 2. 液体Al
液态ji
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTM3MDY4MCwxMDg4NTc4MTc0LDgxNDg3OD
g4MCw3MzM0MDc4MDAsMjA5NzU4NzY1LC0xMzA5NTg4MzExLDkw
MDU3OTgwNSwxMjM2Njk0Njc1LDIwNDAyOTc2MjJdfQ==
-->