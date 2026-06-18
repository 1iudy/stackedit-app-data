# **PAD位错模型构建和Peierls应力计算**

## 1.PAD位错模型构建
主要参考了atomsk的官方教程文档：[Atomsk - Tutorial - Edge dislocation in Al](https://atomsk.univ-lille.fr/tutorial_Al_edge.php)

1. **构建BCC结构Nb晶胞**
分别选择BCC晶格滑移方向<111>为x轴，滑移面法向<11-2>为y轴，<-110>方向为z轴，Nb晶格常数为3.294$\mathring{A}$，使用atomsk命令生成晶胞如下：

```  
atomsk --create bcc 3.294 Nb orient [111] [11-2] [-110] Nb_unitcell.xsf

```

2. **构建PAD模型下半部分**
后续计算Peierls应力时需要模型足够大，参考atomsk教程中，构建80×40的下半部分模型，沿x方向变形，刚度为0，命令如下：
```
atomsk Nb_unitcell.xsf -duplicate 80 40 1 -deform X 0.00625 0.0 bottom.xsf
```

3. **构建PAD模型上半部分**
与下半部分的构建方法相同，但构建上半部分时需要在x方向上增加一个版原子面：
```
atomsk Nb_unitcell.xsf -duplicate 81 40 1 -deform X 0.0061728 0.0 top.xsf
```
4. **拼接模型上下部分**
```
 atomsk --merge Y 2 bottom.xsf top.xsf Nb_pad.xsf
```
使用ovito查看得到的结构
![输入图片说明](/imgs/2026-06-17/rhDbPlzCB56Mf4tU.png)
（似乎x方向建的不够长）
按理来说需要进行结构弛豫才能得到合理的位错模型，但是计算periels应力的过程中会再进行弛豫，此处就省略掉了

## 2. Peierls应力计算
参考内容：[On the significance of model design in atomistic calculations of the Peierls stress in Nb](https://www.sciencedirect.com/science/article/pii/S0927025620306418#b0135)
文章中提供了三种加载方式，分别为位移加载、剪切加载、应力加载，并提供了相应的输入文件模板：[wrj2018/CMS_2020](https://github.com/wrj2018/CMS_2020/tree/master)
1. **位移加载方式**
通过固定上下各为6$\mathring{A}$大小边界层并向上边界施加位移实现，根据文献中的in文件修改如下：
```
variable sname index ms
log ${sname}.log

variable ustrain equal -2e-5
variable energyConv equal 160.21917 # conversion factor
variable Etol equal 1.0e-12
 
# ------------------------ INITIALIZATION ----------------------------
units metal
dimension 3
boundary p f p
atom_style atomic
read_data ../pad/Nb_pad.lmp
 
neighbor 2.5 bin
neigh_modify delay 10 every 2 check yes one 30000 page 300000
  
# Potentials
pair_style eam/alloy
pair_coeff * * ../Ackland_Nb.eam.fs Nb
  
# -------------------------RELAXATION SETTINGS ---------------------------------
variable tmp0 equal "ylo+6"
variable ylo0 equal ${tmp0}
variable tmp1 equal "yhi-6"
variable yhi0 equal ${tmp1}
  
region upper block INF INF ${yhi0} INF INF INF units box
region lower block INF INF INF ${ylo0} INF INF units box
  
group upper region upper
group lower region lower
group boundary union upper lower
group mobile subtract all boundary
  
# define the force to apply
variable tmp2 equal lx
variable LXX equal ${tmp2}
variable tmp3 equal ly
variable LYY equal ${tmp3}
variable tmp4 equal lz
variable LZZ equal ${tmp4}
  
variable udisp equal ${ustrain}*${LYY}
 

variable c loop 2501
label loopc
  
variable tstrain equal ${ustrain}*($c-1)
  
compute tfx upper reduce sum fx
  
thermo 1
thermo_style custom step c_tfx
  
run 0
  
variable tfxx equal c_tfx
  
 
variable sigmaxx equal -${tfxx}*${energyConv}/(${LXX}*${LZZ}) #####GPa
  

thermo 20
thermo_style custom step cpu etotal pe press pxx pyy pzz pxy pxz pyz vol density
  
displace_atoms upper move ${udisp} 0.0 0.0 units box
  
fix 1 upper setforce 0.0 0.0 NULL
fix 2 lower setforce 0.0 0.0 NULL
 
min_style cg
minimize ${Etol} ${Etol} 100000 100000
  
min_style fire
minimize ${Etol} ${Etol} 100000 100000
  
print "${tstrain} ${sigmaxx}" append strain-stress.txt
  
uncompute tfx
unfix 1
unfix 2
next c
jump Peierls_in_2.lmp loopc
  
print "All done"
```
根据结果绘制应力应变曲线：
![位移控制派纳力曲线](/imgs/2026-06-18/xX1UKRh3HBvhAgyz.png)
Peierls应力大约在507MPa左右收敛
2. **剪切加载方式**
类似地，通过对PAD模型施加剪切应变实现，输入文件如下：
```
# ------- initialization -------------------------------------------------
units metal
dimension 3
boundary p s p
atom_style atomic
neighbor 2 bin
  
# ------- create basic geometry -------------------------------------------------
  
read_data ../pad/Nb_pad.lmp
change_box all triclinic
  
 
# ------- EAM potentials -------------------------------------------------
  
pair_style eam/alloy
pair_coeff * * ../Ackland_Nb.eam.fs Nb
  
# ------- timestep & log -------------------------------------------------
  
thermo_style custom step cpu temp pxx pyy pzz pxy pxz pyz xy xz yz pe
thermo 100
  
# ------- fixed region -------------------------------------------------
variable tmp0 equal "ylo+6"
variable ylo0 equal ${tmp0}
variable tmp1 equal "yhi-6"
variable yhi0 equal ${tmp1}
  
region upper block INF INF ${yhi0} INF INF INF units box
region lower block INF INF INF ${ylo0} INF INF units box
  
group upper region upper
group lower region lower
group boundary union upper lower
group mobile subtract all boundary
  
# ------- energy minimization -------------------------------------------------
  
variable Etol equal 1.0e-12
  
min_style cg
fix 1 all box/relax x 0 z 0 nreset 1
minimize ${Etol} ${Etol} 100000 100000
unfix 1
  
min_style fire
minimize ${Etol} ${Etol} 100000 100000
  
# ------- MS Load ---------
  
fix freeze boundary setforce 0.0 NULL NULL
  
variable LY equal ly
  
variable Eel1 equal 0.00
variable Eel2 equal -0.02
  
variable theta equal PI/2 # in units of radian
  
variable Epzy equal yz/${LY}
variable Epxy equal xy/${LY}
  
variable Lzy1 equal (${Epzy}+${Eel1})*${LY}*cos(${theta})
variable Lzy2 equal (${Epzy}+${Eel2})*${LY}*cos(${theta})
  
variable Lxy1 equal (${Epxy}+${Eel1})*${LY}*sin(${theta})
variable Lxy2 equal (${Epxy}+${Eel2})*${LY}*sin(${theta})
  
variable N equal 2001
  
label loop
variable a loop ${N}
  
variable zyTilt equal ${Lzy1}+(${a}-1)/(${N}-1)*(${Lzy2}-${Lzy1})
  
variable xyTilt equal ${Lxy1}+(${a}-1)/(${N}-1)*(${Lxy2}-${Lxy1})
  
change_box all yz final ${zyTilt} remap units box
  
change_box all xy final ${xyTilt} remap units box
  
min_style cg
minimize ${Etol} ${Etol} 100000 100000
  
min_style fire
minimize ${Etol} ${Etol} 100000 100000
  
variable yTilt equal ${zyTilt}*cos(${theta})+${xyTilt}*sin(${theta})
variable strain equal ${yTilt}/${LY}
variable PY equal (pyz*cos(${theta})+pxy*sin(${theta}))*0.1 #####MPa
  
print "${yTilt} ${strain} ${PY}" append edge_pad.txt
print "${theta} ${yTilt} ${xyTilt}" append angle_log.txt
  
 
run 0
  
next a
jump Peierls_in.lmp loop
```
由于theta是90°固定值，yTilt其实跟xyTilt是相等的，可能源代码给的是更加通用的形式。通过分析结果文件得到的应力曲线：
![剪切控制派纳力曲线](/imgs/2026-06-18/oc8uwsEih1yRXGBB.png)
相对于位移加载控制方式下的曲线更加平滑
3. 
<!--stackedit_data:
eyJoaXN0b3J5IjpbMzE5OTQxODM4LDE0Nzk1NjY0MDcsMzEzMD
I1ODkyLC0xOTA4MTU2NTU1LDE3NjY0MzQzODMsLTQyNDg1ODUy
NywxNDMyNzM3MzMsLTY5NjkyNjQwNCwxNTM0NzUwMjIwLDExMD
M1OTkyNDMsMTY2Nzg3MzU4LDEzOTA2MDE5NDUsLTE3NTg3NzE0
MzNdfQ==
-->