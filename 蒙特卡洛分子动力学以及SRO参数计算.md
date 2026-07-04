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
Atomsk 利用Voronoi镶嵌构造多晶体，首先在模拟盒内的给定位置引入节点，随后在周期性边界条件下将相邻节点相连，绘制出连线的法线作为未来多晶粒的晶界。将原子“种子”（例如晶胞）按照给定的晶体取向放置在节点位置。种子在空间的三个方向上展开并去除扩展至晶粒外的原子，当各个晶粒内的种子扩展满对应的晶粒并去除晶粒外的晶粒后得到最终的多晶体。
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
**单晶**
按照文献内容，在MC模拟进行之前先在1000K下50ps的NPT弛豫，in文件如下：
```
units metal
dimension 3
boundary p p p
timestep 0.001
  
atom_style atomic
read_data relaxed_Nb25Ta25Hf5Zr45.lmp
  
pair_style meam
pair_coeff * * ../../../library.meam Hf Nb Ta Ti Zr ../../../HfNbTaTiZr.meam Zr Nb Ta Hf
  
#目标温度1000K下在NPT系综下进行零压弛豫
fix 1 all npt temp 1000.0 1000.0 0.1 iso 0.0 0.0 1.0
thermo 1000
thermo_style custom step pe press etotal temp vol lx ly lz
  
run 50000
  
unfix 1
  
write_data NPT_Nb25Ta25Hf5Zr45.lmp
```

<!--stackedit_data:
eyJoaXN0b3J5IjpbNzk2MDQ1NzIxLDE4ODIzNDY3MjQsMTkzMD
gzOTY3OSwxMDc2MzE3MTU0LC05NjgxNjA3OTYsODU1ODg5Njg3
LC00NDU1MjMxMjUsLTIwMDQxNjk0MzQsLTE3Nzg2ODA3MDMsMj
M0ODE2MjM0LDExODYzOTYyOTQsMTc2MDAxMDMyOCwyMTM5ODAz
MDIzLC01MjY0MjYzMzEsMTk3NzU4NDA4LDQ2NTAxOTk4NCwxMj
MzNTg3NjAyLDc4MTI4OTU3MiwxNjIxNzEwMjA4LC0xMjgyNzg4
NTY0XX0=
-->