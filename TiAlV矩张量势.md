# TiAlV矩张量势函数验证
参考文献：[Moment tensor potential and its application in the Ti-Al-V multicomponent system](https://doi.org/10.1103/PhysRevMaterials.9.063601)  
主要进行mtp势函数验证的部分
## 1. 基本性质
### 1.1 EOS曲线
通过改变晶格常数间隔计算EOS曲线，在lammps中实现方式为：通过循环生成不同晶格常数下的计算盒子，计算原子能量后将盒子清空然后进入下一个循环
```
variable a_start equal 3.80
variable a_end equal 4.30
variable da equal 0.05
variable nsteps equal round((${a_end}-${a_start})/${da})
  
print "# V/atom E/atom" file eos.dat screen no	
  
variable n loop 0 ${nsteps}	
label loop
  
units metal
dimension 3
boundary p p p
atom_style atomic
  
variable a equal ${a_start}+${n}*${da}
  
# 直接设为绝对尺寸，历史影响归零
lattice fcc ${a}
region box block 0 4 0 4 0 4 units lattice
create_box 1 box
create_atoms 1 box
  
pair_style mlip mlip.ini
pair_coeff * *
mass 1 26.98
  
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
  
run 0
  
variable natom equal count(all)
variable V equal (lx*ly*lz)
variable E equal pe
variable V_atom equal (lx*ly*lz)/count(all)
variable E_atom equal pe/${natom}
  
print "${natom} ${V} ${E}"
  
print "${V_atom} ${E_atom}" append eos.dat screen no
  
clear
  
next n
jump SELF loop
```

>注意：这个势函数存在较大的问题，问题如下：[Lammps+MLIP：Ti/Al/V机器学习势的安装与使用](https://zhuanlan.zhihu.com/p/1923794430491620846)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDAyNDM3MDIxLDE5ODY2MzMyNjMsLTM1MD
MyNzg3NywtMjQ5ODgyODIzLC03ODg0ODI0OTUsMjA0MDI5NzYy
Ml19
-->