# NbTaHfZr分子动力学和蒙特卡洛模拟
参考文章：[Assessing the Effects of Chemical Composition and Short-Range Ordering on the Tensile Deformation Behavior of Nanocrystalline High-Entropy Alloys Nb–Ta–Hf–Zr: A Combined Study on Molecular Dynamics and Monte Carlo Simulation](https://link.springer.com/article/10.1007/s11661-025-08003-z#Sec2)
## 1.NbTaHfZr单晶和纳米多晶模型构建
根据文献内容，Nb和Ta的含量始终相等，成分变化为5，共有45种不同的合金成分，但文章仅主要分析了Nb25Ta25Hf5Zr45、Nb25Ta25Hf45Zr5、Nb45Ta45Hf5Zr5、NbTaHfZr四种结构。
文献中单晶超胞约有14000个原子，因此对于单晶结构使用元胞**各方向扩胞12**；纳米多晶结构包含102000个原子，使用元胞**各方向扩胞47**。查询文献和数据库，TaNbHf高熵合金晶格常数约为3.36Å，

<!--stackedit_data:
eyJoaXN0b3J5IjpbNjY4MDEzMTkzLDE5Nzc1ODQwOCw0NjUwMT
k5ODQsMTIzMzU4NzYwMiw3ODEyODk1NzIsMTYyMTcxMDIwOCwt
MTI4Mjc4ODU2NCwxMzYxMDI3OTYzLDEzNjEwMjc5NjMsLTI1Nz
gxNDkzMiwxOTIzNDAwNDQwLC0xNjc3NDc1MDYxLC0yMDg4NzQ2
NjEyXX0=
-->