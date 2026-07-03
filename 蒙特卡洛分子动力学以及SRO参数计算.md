# NbTaHfZr分子动力学和蒙特卡洛模拟
参考文章：[Assessing the Effects of Chemical Composition and Short-Range Ordering on the Tensile Deformation Behavior of Nanocrystalline High-Entropy Alloys Nb–Ta–Hf–Zr: A Combined Study on Molecular Dynamics and Monte Carlo Simulation](https://link.springer.com/article/10.1007/s11661-025-08003-z#Sec2)
## 1.NbTaHfZr单晶和纳米多晶模型构建
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
文献构建的结构尺寸为13×13×13nm³，假设晶粒数量为10，写入polycrystal.txt文件：
```
box 130 130 130
random 10
```

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTk1ODc5NDE5MywyMzQ4MTYyMzQsMTE4Nj
M5NjI5NCwxNzYwMDEwMzI4LDIxMzk4MDMwMjMsLTUyNjQyNjMz
MSwxOTc3NTg0MDgsNDY1MDE5OTg0LDEyMzM1ODc2MDIsNzgxMj
g5NTcyLDE2MjE3MTAyMDgsLTEyODI3ODg1NjQsMTM2MTAyNzk2
MywxMzYxMDI3OTYzLC0yNTc4MTQ5MzIsMTkyMzQwMDQ0MCwtMT
Y3NzQ3NTA2MSwtMjA4ODc0NjYxMl19
-->