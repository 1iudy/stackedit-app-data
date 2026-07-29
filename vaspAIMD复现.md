
# 液态 Al–Zn 合金冷却过程中原子间距的反常膨胀（AIMD和MLFF）
参考文献[Anomalous expansion of interatomic distance in liquid Al–Zn alloy during cooling | The Journal of Chemical Physics | AIP Publishing](https://pubs.aip.org/aip/jcp/article/163/5/054307/3357370/Anomalous-expansion-of-interatomic-distance-in)
## AIMD计算细节
使用atomsk分别构建Al–75Zn、Al–23Zn和Al–38Zn的原子结构，首先是200个原子的AIMD结构，使用atomsk构建FCC Al超胞并进行扩胞，根据质量百分数换算成原子百分数并使用随机替代的方式构建200原子的Al–75Zn超胞结构
```
 atomsk --create fcc 4.05 Al -duplicate 5 5 2 -sub Al Zn 44.7% POSCAR  
```
原文献进行了AIMD部分分为两步，


<!--stackedit_data:
eyJoaXN0b3J5IjpbMjA5NjkxMzgxMywtMTE5MDUzNjY0MiwxMD
Y0ODAwMjcyLC0zNzEwMTAyNzFdfQ==
-->