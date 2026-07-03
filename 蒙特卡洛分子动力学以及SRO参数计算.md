# NbTaHfZr分子动力学和蒙特卡洛模拟
参考文章：[Assessing the Effects of Chemical Composition and Short-Range Ordering on the Tensile Deformation Behavior of Nanocrystalline High-Entropy Alloys Nb–Ta–Hf–Zr: A Combined Study on Molecular Dynamics and Monte Carlo Simulation](https://link.springer.com/article/10.1007/s11661-025-08003-z#Sec2)
## 1.NbTaHfZr单晶和纳米多晶模型构建
根据文献内容，Nb和Ta的含量始终相等，成分变化为5，共有45种不同的合金成分，但文章仅主要分析了Nb25Ta25Hf5Zr45、Nb25Ta25Hf45Zr5、Nb45Ta45Hf5Zr5、NbTaHfZr四种结构。
文献中单晶超胞约有14000个原子，因此对于单晶结构使用元胞**各方向扩胞20**；纳米多晶结构包含102000个原子，使用元胞**各方向扩胞47**。查询文献和数据库，TaNbHf高熵合金晶格常数约为3.36Å。
以Nb25Ta25Hf5Zr45为例，单晶结构构建命令如下：
```
atomsk --create bcc 3.36 Nb -duplicate 20 20 20 Nb_supercell.xsf
```
超胞结构共16000个原子，
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE4MjIxMTk1NTgsMTk3NzU4NDA4LDQ2NT
AxOTk4NCwxMjMzNTg3NjAyLDc4MTI4OTU3MiwxNjIxNzEwMjA4
LC0xMjgyNzg4NTY0LDEzNjEwMjc5NjMsMTM2MTAyNzk2MywtMj
U3ODE0OTMyLDE5MjM0MDA0NDAsLTE2Nzc0NzUwNjEsLTIwODg3
NDY2MTJdfQ==
-->