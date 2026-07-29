
# 液态 Al–Zn 合金冷却过程中原子间距的反常膨胀（AIMD和MLFF）
参考文献[Anomalous expansion of interatomic distance in liquid Al–Zn alloy during cooling | The Journal of Chemical Physics | AIP Publishing](https://pubs.aip.org/aip/jcp/article/163/5/054307/3357370/Anomalous-expansion-of-interatomic-distance-in)
## AIMD计算细节
使用atomsk分别构建Al–75Zn、Al–23Zn和Al–38Zn的原子结构，首先是200个原子的AIMD结构，使用atomsk构建FCC Al超胞并进行扩胞，根据质量百分数换算成原子百分数并使用随机替代的方式构建200原子的Al–75Zn超胞结构
```
 atomsk --create fcc 4.05 Al -duplicate 5 5 2 -sub Al Zn 44.7% POSCAR  
```
原文献进行了1600 K 熔化后逐步冷却至 700 K（步长 100 K）的AIMD模拟，每个温度分为两步：先在上个温度下进行2500步冷却至当前温度，随后在该温度下进行12000 MD 步弛豫，最后 4000 步统计


<!--stackedit_data:
eyJoaXN0b3J5IjpbMTA1OTgxOTczNiwtMTE5MDUzNjY0MiwxMD
Y0ODAwMjcyLC0zNzEwMTAyNzFdfQ==
-->