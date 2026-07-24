#  短程有序对单晶CoCrFeMnNi高熵合金空隙增长的影响
参考文献：[Effects of short-range order on void growth in single-crystalline CoCrFeMnNi high-entropy alloys](https://www.sciencedirect.com/science/article/pii/S1359645426006142#b45)
>！！原文对于WC参数的计算疑似使用了mdapy程序包但是未加引用
## 1.CoCrFeMnNi模型构建
使用lammps构建，相对于atomsk构建，可以构建三坐标轴方向等长的模拟盒子（理论上讲atomsk也可以实现，不确定atomsk的-duplicate命令扩胞是否可以使用浮点数）
```
units metal
dimension 3
boundary p p p
atom_style atomic
  
lattice fcc 3.60 orient x 1 1 1 orient y 1 -1 0 orient z 1 1 -2 #确定晶胞构建方向
 
variable N equal 50.0	#总体的扩胞数量，各个方向上除以模长保证扩胞后三方向距离相等
variable Nx equal ${N}/sqrt(3.0) 
variable Ny equal ${N}/sqrt(2.0)
variable Nz equal 1.5*${N}/sqrt(6.0)
  
region box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice
create_box 5 box	#创建模拟盒子并用5种原子填充
mass 1 55.845 # Fe
mass 2 58.933 # Co
mass 3 51.996 # Cr
mass 4 58.693 # Ni
mass 5 54.938 # Mn
create_atoms 1 box
  
#对Fe原子随机替换
set type 1 type/ratio 2 0.2 66531
set type 1 type/ratio 3 0.25 66531
set type 1 type/ratio 4 0.33333 66531
set type 1 type/ratio 5 0.5 66531
  
write_data CoCrFeMnNi_111.lmp
```
目前使用修改过的lammps代码实现了fix atom/swap交换两种以上的元素类型，大幅提升了MC速度。


**7.24**
目前已经完成了**1264 / 2000**次MC周期，
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTg4MDQzNzI3OCwxMzEwMTU3MDc3LC0xOD
U4MzM5OTE3LC0xOTM0MTIwMjYyLC02ODc5NjIyMzcsLTk5Njkw
MzQ4OCwxMzkzODA3MTY2LC0xMDgyMDY5MjQzLDIxMTUwMzIzMj
MsLTk1MDk0ODQyMCwtMjk5NDM2NjIxLC04MzE2NDE3NjUsLTE2
Nzk2NzkyODFdfQ==
-->