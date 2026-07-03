# NbTaHfZr分子动力学和蒙特卡洛模拟
参考文章：[Assessing the Effects of Chemical Composition and Short-Range Ordering on the Tensile Deformation Behavior of Nanocrystalline High-Entropy Alloys Nb–Ta–Hf–Zr: A Combined Study on Molecular Dynamics and Monte Carlo Simulation](https://link.springer.com/article/10.1007/s11661-025-08003-z#Sec2)
## 1.NbTaHfZr单晶和纳米多晶模型构建
根据文献内容，Nb和Ta的含量始终相等，成分变化为5，共有45种不同的合金成分，但文章仅主要分析了Nb25Ta25Hf5Zr45、Nb25Ta25Hf45Zr5、Nb45Ta45Hf5Zr5、NbTaHfZr四种结构。
文献中单晶超胞约有14000个原子，因此对于单晶结构使用元胞**各方向扩胞19**；纳米多晶结构包含102000个原子，使用元胞**各方向扩胞47**。查询文献和数据库，TaNbHf高熵合金晶格常数约为3.36Å。
以Nb25Ta25Hf5Zr45为例，单晶结构构建命令如下：
```
atomsk --create bcc 3.36 Nb -duplicate 19 19 19 Nb_supercell.xsf
```
超胞结构共13718个原子，
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTU2Mjk0MDYxOSwxOTc3NTg0MDgsNDY1MD
E5OTg0LDEyMzM1ODc2MDIsNzgxMjg5NTcyLDE2MjE3MTAyMDgs
LTEyODI3ODg1NjQsMTM2MTAyNzk2MywxMzYxMDI3OTYzLC0yNT
c4MTQ5MzIsMTkyMzQwMDQ0MCwtMTY3NzQ3NTA2MSwtMjA4ODc0
NjYxMl19
-->