# NbTaHfZr分子动力学和蒙特卡洛模拟
参考文章：[Assessing the Effects of Chemical Composition and Short-Range Ordering on the Tensile Deformation Behavior of Nanocrystalline High-Entropy Alloys Nb–Ta–Hf–Zr: A Combined Study on Molecular Dynamics and Monte Carlo Simulation](https://link.springer.com/article/10.1007/s11661-025-08003-z#Sec2)
## 1. NbTaHfZr单晶和纳米多晶模型构建
根据文献内容，Nb和Ta的含量始终相等，成分变化为5，共有45种不同的单晶结构合金成分，但文章仅主要分析了Nb25Ta25Hf5Zr45、Nb25Ta25Hf45Zr5、Nb45Ta45Hf5Zr5、NbTaHfZr四种纳米多晶结构。
文献中单晶超胞约有14000个原子，因此对于单晶结构使用元胞**各方向扩胞20**；纳米多晶结构包含102000个原子，使用元胞**各方向扩胞47**。查询文献和数据库，TaNbHf高熵合金晶格常数约为3.36Å。
以Nb25Ta25Hf5Zr45为例，单晶结构构建命令如下：
```
atomsk --create bcc 3.36 Nb -duplicate 20 20 20 Nb_supercell.xsf
```
超胞结构共16000个原子，使用Ta、Hf、Zr分别替代原有Nb：
```
atomsk Nb_supercell.xsf -select random 4000 Nb -sub Nb Ta -select random 800 Nb -sub Nb Hf -select random 7200 Nb -sub Nb Zr Nb25Ta25Hf5Zr45.xsf
```
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-03/R3jzC4oFjNKwAFq8.png)
纳米多晶结构构建参考atomsk官方教程：
Atomsk 利用Voronoi镶嵌构造多晶体，首先在模拟盒内的给定位置引入节点，sui'ha节点与其邻居（红线）相连。使用周期性边界条件。

（c） 找到红线的法线（蓝线）。这些蓝线定义了未来晶粒的轮廓，即晶界。

（d） 原子“种子”（例如晶胞）放置在节点位置，按照给定的晶体取向。

（e） 种子在空间的三个方向上展开。去除粒子外的原子。

（f） 在所有种子都被扩展并在各自晶粒内切割后，得到最终的多晶体。
文献构建的结构尺寸为13×13×13nm³，假设晶粒数量为10，写入polycrystal.txt文件：
```
box 130 130 130
random 10
```
创建元胞种子文件：
```
atomsk --create bcc 3.36 Nb Nb_seed.xsf
```
使用polycrystal命令生成Nb多晶结构：
```
 atomsk --polycrystal Nb_seed.xsf polycrystal.txt Nb_polycrystal.cfg -wrap
 ```
 用Ta、Hf、Zr分别替代原有Nb：
 ```
  atomsk Nb_polycrystal.cfg -select random  25% Nb -sub Nb Ta -select random 6.6666% Nb -sub Nb Hf -select random 64.2857% Nb -sub Nb Zr Nb25Ta25Hf5Zr45.cfg
``` 
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-03/CXoCszlmptwuLmJz.png)![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-03/JV3TMRSabhOBWHt1.png)
 
 共有115841原子
 ## **2. MC/MD**
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDE4OTE3MDUzLDEwNzYzMTcxNTQsLTk2OD
E2MDc5Niw4NTU4ODk2ODcsLTQ0NTUyMzEyNSwtMjAwNDE2OTQz
NCwtMTc3ODY4MDcwMywyMzQ4MTYyMzQsMTE4NjM5NjI5NCwxNz
YwMDEwMzI4LDIxMzk4MDMwMjMsLTUyNjQyNjMzMSwxOTc3NTg0
MDgsNDY1MDE5OTg0LDEyMzM1ODc2MDIsNzgxMjg5NTcyLDE2Mj
E3MTAyMDgsLTEyODI3ODg1NjQsMTM2MTAyNzk2MywxMzYxMDI3
OTYzXX0=
-->