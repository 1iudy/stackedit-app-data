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
~~随后进行MC/MD模拟，一个周期内MC先运行12000步，每次交换一个原子，随后进行5000步MD弛豫，共进行80个周期，in文件如下：~~
```
units metal
dimension 3
boundary p p p
timestep 0.001
  
atom_style atomic
read_data NPT_Nb25Ta25Hf5Zr45.lmp
  
pair_style meam
pair_coeff * * ../../../library.meam Hf Nb Ta Ti Zr ../../../HfNbTaTiZr.meam Zr Nb Ta Hf
  
variable swap_fra equal 1.0/16000 #交换原子占比
  
#MC周期循环
variable period_count loop 80
label period_start
  
#MC设置
##fix 2 all sgcmc 1 ${swap_fra} 1000.0 0.0 0.0 0.0 #已通过修改lammps源代码的方式实现了多种原子进行atom/swap
fix 2 all atom/swap 1 1 89536 1000.0 types 1 2 3 4
fix 3 all nvt temp 1000.0 1000.0 0.1
  
thermo 2000
thermo_style custom step pe press temp vol lx ly lz
  
dump 2 all custom 4000 dump_MC_MD_Nb25Ta25Hf5Zr45.lammpstrj id type x y z
  
run 12000
  
unfix 2
undump 2
  
dump 3 all custom 2500 dump_MC_MD_Nb25Ta25Hf5Zr45.lammpstrj id type x y z
  
run 5000
unfix 3
undump 3
  
next period_count
jump in_MC_MD.lmp period_start
  
write_data MC_MD_Nb25Ta25Hf5Zr45.lmp
```
~~但fix sgcmc方法在不使用EAM势并添加atomic/energy yes的关键词下只能串行运行，只使用单核计算，速度较慢。~~（搞错了，原文中的Metropolis接受准则应该是通过fix atom/swap命令实现的）关于MC的Metropolis接受准则相关内容查看[【MC模拟入门】MC模拟原理：随机漫步与Metropolis算法 - CSDN文库](https://wenku.csdn.net/column/4mgwvxeokf)
fix atom/swap命令可以内置实现Metropolis接受准则，语法如下
```
fix 2 all atom/swap 1 1 12156 1000.0 types 1 2
```
分别定义MC间隔、单次交换原子对的数量、随机种子、模拟温度。~~对于高熵合金这种2种以上元素参与交换的情况，需要设置semi-grand yes使用半巨正则系综，通过types设置参与交换原子的种类，当使用semi-grand yes时需要设置各个元素的化学势，使用mu将各元素相对化学势设置为0，保证原子交换没有偏好性。~~（**semi grand半正则系综实际上是选定原子并对其进行替换，会在不改变原子总数的前提下改变原子种类的比例，跟本文的原子对之间进行交换并不相符**）
最终MC的实现方式是通过多个fix atom/swap实现多种原子对的交换，执行一步MC实际上对于每一种元素对均进行了一次交换，不知道是不是符合文献的描述“***Swapping of one atom type with an atom of another type is realized using the widely accepted MC procedure that follows the Metropolis acceptance criterion***”
输入文件相关部分如下：
```
fix 12 all atom/swap 1 1 12156 1000.0 types 1 2
fix 13 all atom/swap 1 1 23765 1000.0 types 1 3
fix 14 all atom/swap 1 1 38948 1000.0 types 1 4
fix 23 all atom/swap 1 1 66534 1000.0 types 2 3
fix 24 all atom/swap 1 1 24516 1000.0 types 2 4
fix 34 all atom/swap 1 1 99564 1000.0 types 3 4
fix 3 all nvt temp 1000.0 1000.0 0.1
  
thermo 2000
thermo_style custom step pe press temp vol lx ly lz
  
run 12000
  
unfix 12
unfix 13
unfix 14
unfix 23
unfix 24
unfix 34
```
## **3. Warren Cowley Parameter计算**
进行MC交换后的结构已经具有明显的局域有序现象：
![MC/MD完成后的Nb25Ta25Hf5Zr45结构](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-13/hlFk29bak3z2aHwa.png)
目前计算warren cowley的思路基本都是在第一邻近截断半径的范围内，统计每一个中心原子的近邻原子列表，并将统计结果代入计算公式得到最终的WCP:
$$\alpha_{IJ}^n = 1 - \frac{P_n^{J|I}}{c_J}$$
但是目前似乎并没有统一的Warren Cowley参数计算工具以供使用，相关文献对于Warren Cowley参数计算的实现也各不相同，目前使用了mdapy软件包、Warren Cowley Parameters软件包和AI辅助编程进行了计算，计算结果存在较多问题。
### 3.1 mdapy实现
mdapy是用于分析分子动力学（MD）模拟生成的原子轨迹的python软件包，可以直接处理lammps的dump和data文件。mdapy于今年进行了全面重写，发布了1.0.0及后续版本（目前已经更新至1.0.8a1）将原有依赖Tichi的部分通过c++实现，但是新版本暂时还没有开发wcp的计算模块，因此仍然使用[旧版本0.11.5](https://github.com/mushroomfire/mdapy/tree/mdapy_old)。
mdapy对于wcp计算结果以矩阵形式给出，比较方便在文献中展示，但是数值范围一般只能在-2~1之间，进行精细调控展示稍微麻烦。
Nb25Ta25Hf5Zr45计算结果如下：![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-13/q5ljNoYL4BDACvzZ.png)

### 3.2 Warren Cowley Parameter程序包实现
[Warren Cowley Parameter](https://github.com/killiansheriff/WarrenCowleyParameters)软件包可以作为ovito的modifier，也可以通过ovito的python接口进行实现。通过python接口自动计算，结果如下：
**1NN WC参数：**

| \ | Zr |	Nb |	Ta | Hf |	
|--|--|--|--|--|
| **Zr** | -0.36 |	0.32 |	0.41 | -0.26 |
| **Nb** |0.32|0.46| -1.06 |0.05|
| **Ta** | 0.39 |	-1.06 |	0.23 |0.45  |
| **Hf** | -0.27 |0.07|0.48|-0.15 |

**2NN WC参数：**
| \ | Zr |	Nb |	Ta | Hf |	
|--|--|--|--|--|
| **Zr** | -0.25 |	0.08 |	0.39 | -0.20 |
| **Nb** |0.07|-0.70| 0.49 |0.42|
| **Ta** | 0.41 |	0.49 |-1.19 |0.03  |
| **Hf** | -0.27 |0.39|0.00|-0.38 |

这种形式下的结果存在不对称的问题
> 此处计算的是未经过最小化弛豫的结构，经过最小化弛豫后使用这个程序包计算的结果与RDF计算结果基本相同
### 3.3 通过径向分布函数（RDF）实现
参考了[Effects of Chemical Short-Range Order and Temperature on Basic Structure Parameters and Stacking Fault Energies in Multi-Principal Element Alloys](https://www.mdpi.com/2673-3951/5/1/19)中的实现方式，先通过compute rdf命令计算各个原子对的径向分布函数，然后根据积分求配位数：
$$z_{AB} = \int_0^{r_{\text{cut}}} \rho \cdot g_{AB}(r) \cdot 4\pi r^2 \, dr$$
进而求出WC参数：
$$\alpha_{AB} = 1 - \frac{\displaystyle \int_0^{r_{\text{cut}}} \rho \cdot g_{AB}(r) \cdot 4\pi r^2 \, dr}{c_B \cdot \displaystyle \int_0^{r_{\text{cut}}} \rho \cdot g_A^{\text{total}}(r) \cdot 4\pi r^2 \, dr}$$
将得到的RDF结果通过python脚本处理得到Warren Cowley参数如下：
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-14/2EJZxr7pz8vqXKis.png)

### 小结
三种方式计算出来的WC参数虽然数值存在差异，但是元素间的偏聚/有序分布趋势的描述还是比较一致的，**但是最大的问题在于当前复现的warren cowley parameter跟原文献不符**，不清楚是MC随机性的问题还是其他什么问题
![原文献中WC参数的计算结果](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11661-025-08003-z/MediaObjects/11661_2025_8003_Fig6_HTML.png)
（原文结果只标明了等原子比和Hf45Zr5，其余两种没有标明）
**与原文献参数对比（以纳米多晶Nb25Ta25Hf45Zr5为例）**
|原子对| 计算值 |	文献值|
|--|--|--|
| Zr-Zr | -0.96 |-0.94|
|Zr-Nb|0.46|0.27|
|Zr-Ta|0.44|0.23|
Zr-Hf|-0.30|-0.21|
Nb-Nb|-0.30|-0.36|
Nb-Ta|0.06|0.04
Nb-Hf|0.09|0.14
Ta-Ta|-0.61|-0.41
Ta-Hf|0.25|0.17
Hf-Hf|-0.11|-0.17

**就WC参数来看，似乎对于16000原子的单晶结构，wc参数的计算偏差较大，而原子数较多的多晶结构相对来说偏差小一些**

另外，混合蒙特卡洛分子动力学前后的多晶结构对比如下：
MC/MD前：
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-22/oDzB6mDZJMprRQ3R.png)
MC/MD后：
![输入图片说明](https://raw.githubusercontent.com/1iudy/Learning_markdown_files/images/imgs/2026-07-22/uVj1b5x9VtfKmvKB.png)
晶界范围变大的现象在NPT弛豫过程后就已经发生，猜测是NbTaHfZr在1000K下的正常变化？
生成了45种不同的合金成分
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0MTI4MDM1MjQsMTkxOTE3NzE1NSwtOD
I2NjA3NzAyLC0yNTc0MzQzMjAsLTE2MDgzNTE5MjAsLTEyNTU3
NzkxNTIsMTQzMjkxMjU2MCwyMDA5OTUyNzQyLDE0MDEyNjQ3NT
AsLTE0NzkyMjg4MzIsMTM4NDA0MzcwMSwtMTk3MDI1MDE4Niwx
MjY5NTEzMDk2LC0xMDk1OTcwMjIwLC0xMTMxMzA2MTU1LC0xMT
YzMzMwMDc5LC05MzQ2ODc1MDgsLTExMDQ0NDcxOTIsMTkzNzEx
MjEzNywtMjk2ODg5MDQzXX0=
-->