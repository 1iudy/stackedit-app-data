# **PAD位错模型构建和Peierls应力计算**

## 1.PAD位错模型构建
主要参考了atomsk的官方教程文档：[Atomsk - Tutorial - Edge dislocation in Al](https://atomsk.univ-lille.fr/tutorial_Al_edge.php)

1. **构建BCC结构Nb晶胞**
分别选择BCC晶格滑移方向<111>为x轴，滑移面法向<11-2>为y轴，<-110>方向为z轴，Nb晶格常数为3.294$\mathring{A}$，使用atomsk命令生成晶胞如下：

```  
atomsk --create bcc 3.294 Nb orient [111] [11-2] [-110] Nb_unitcell.xsf

```

2. **构建PAD模型下半部分**
后续计算Peierls应力时需要模型足够大，


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE2NzAzNTUzNjIsMTUzNDc1MDIyMCwxMT
AzNTk5MjQzLDE2Njc4NzM1OCwxMzkwNjAxOTQ1LC0xNzU4Nzcx
NDMzXX0=
-->