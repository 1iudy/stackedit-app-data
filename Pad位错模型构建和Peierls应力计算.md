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


<!--stackedit_data:
eyJoaXN0b3J5IjpbMzEzMDI1ODkyLC0xOTA4MTU2NTU1LDE3Nj
Y0MzQzODMsLTQyNDg1ODUyNywxNDMyNzM3MzMsLTY5NjkyNjQw
NCwxNTM0NzUwMjIwLDExMDM1OTkyNDMsMTY2Nzg3MzU4LDEzOT
A2MDE5NDUsLTE3NTg3NzE0MzNdfQ==
-->