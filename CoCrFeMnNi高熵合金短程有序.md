#  短程有序对单晶CoCrFeMnNi高熵合金空隙增长的影响
参考文献：[Effects of short-range order on void growth in single-crystalline CoCrFeMnNi high-entropy alloys](https://www.sciencedirect.com/science/article/pii/S1359645426006142#b45)
## 1.CoCrFeMnNi模型构建
使用lammps构建，相对于atomsk构建，可以构建三坐标轴方向等长的模拟盒子（理论上讲atomsk也可以实现，不确定atomsk的-duplicate命令扩胞是否可以使用浮点数）
```
units metal
dimension 3
boundary p p p
atom_style atomic
  
lattice fcc 3.60 orient x 1 1 1 orient y 1 -1 0 orient z 1 1 -2 #确定晶胞构建方向
 
variable N equal 50.0	#总体的扩胞数量
variable Nx equal ${N}/sqrt(3.0) 
variable Ny equal ${N}/sqrt(2.0)
variable Nz equal 1.5*${N}/sqrt(6.0)
  
region box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice
create_box 5 box
mass 1 55.845 # Fe
mass 2 58.933 # Co
mass 3 51.996 # Cr
mass 4 58.693 # Ni
mass 5 54.938 # Mn
create_atoms 1 box
  
set type 1 type/ratio 2 0.2 66531
set type 1 type/ratio 3 0.25 66531
set type 1 type/ratio 4 0.33333 66531
set type 1 type/ratio 5 0.5 66531
  
write_data CoCrFeMnNi_111.lmp
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTIwMTA4MzcxNiwtOTUwOTQ4NDIwLC0yOT
k0MzY2MjEsLTgzMTY0MTc2NSwtMTY3OTY3OTI4MV19
-->