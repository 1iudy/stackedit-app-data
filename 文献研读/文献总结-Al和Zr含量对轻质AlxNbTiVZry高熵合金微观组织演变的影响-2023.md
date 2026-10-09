# 文献总结-Al 和 Zr 含量对轻质 AlxNbTiVZry 高熵合金微观组织演变的影响-2023

## 阅读来源与翻译说明

本报告以用户 Zotero 文献库中的原论文 PDF 为主，并与同条目下保存的 MDPI 网页快照核对。Zotero 条目键：`BUP7HZ3R`；PDF 附件键：`9NDRQNZ4`；网页快照附件键：`V94XJ6HC`。PDF 共 12 页，正文及作者声明位于第 1–11 页，参考文献位于第 11–12 页。

原论文采用 CC BY 4.0 许可。本报告提供中文翻译和整理，图像沿用原论文；译文不是作者认可或审定的版本。原文出处：[完整 DOI 链接](https://doi.org/10.3390/ma16247581)；[许可说明](https://creativecommons.org/licenses/by/4.0/)。

下文分为详细总结、完整章节翻译、计算细节汇总、原文核对备注和参考文献。完整翻译保留原文的数值、判断强度和引用编号，不将原文中的不一致内容擅自统一。核对备注与译文分开。PDF 的文本提取中存在由图像文字层造成的重复段落，已对照可见页面及网页快照去重；没有删去实际正文内容。图题、表题只给中文，图中的原始标记和坐标保留。

## 一、文献信息

| 项目 | 原文信息 |
|---|---|
| 英文原文标题 | The Effects of the Al and Zr Contents on the Microstructure Evolution of Light-Weight AlxNbTiVZry High Entropy Alloy |
| 中文标题 | Al 和 Zr 含量对轻质 AlxNbTiVZry 高熵合金微观组织演变的影响 |
| 完整 DOI | [https://doi.org/10.3390/ma16247581](https://doi.org/10.3390/ma16247581) |
| 作者，按原文顺序 | Hongwei Yan；Rui Liu；Shenglong Li；Yong’an Zhang；Wei Xiao；Boyu Xue；Baiqing Xiong；Xiwu Li；Zhihui Li |
| 通讯作者 | Yong’an Zhang；Baiqing Xiong |
| 期刊 | Materials |
| 年份、卷、期、文章号 | 2023，16(24)，7581 |
| 投稿日期，Received | 2023 年 10 月 26 日 |
| 返修日期，Revised | 2023 年 11 月 28 日 |
| 接收日期，Accepted | 2023 年 12 月 2 日 |
| 在线发表日期，Published | 2023 年 12 月 9 日 |
| 学术编辑 | Jana Bidulská |
| 数据公开声明 | 原文写明：数据包含在本文内。 |
| 独立数据集或原始数据链接 | 原文未提供。 |
| 代码链接 | 原文未提供。 |
| 外部模型文件、输入文件或计算脚本链接 | 原文未提供。 |

## 二、详细总结

### 2.1 研究问题、对象与范围

论文研究 Al 和 Zr 含量如何影响 AlNbTiVZr 系列轻质耐火高熵合金的铸态微观组织、相组成、密度以及室温压缩性能。作者希望通过成分比较筛选合适的 Al、Zr 配比。本文实验对象是经真空电弧熔炼制备的铸态样品，未报告随后进行固溶、时效或长时间退火的实验，也未报告拉伸或高温压缩实验。

表 1 和各图实际采用的五个样品为：$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}$、$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$、$\mathrm{AlNbTiVZr}$、$\mathrm{Al}_{0.5}\mathrm{NbTiVZr}$、$\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$。摘要和引言所列第四个样品与这里不同，见后文核对备注。

### 2.2 制备与表征证据

作者采用纯度超过 99.9% 的块状或棒状纯金属，去除表面氧化层并超声处理 5 min 后称量。样品在 Ar 保护下电弧熔炼，以预留 Ti 锭吸收残余氧，并在铜坩埚中反复熔炼、翻转 7 次。

密度由阿基米德法测定，每种样品测量 3 次。XRD 用于分析相结构，作者报告使用 Rietveld 方法和 Jade 6.5 处理数据。SEM/BSE 和 EDS 分别用于观察组织与元素分布；EBSD 用于获取取向图和相含量；TEM/SAED 用于进一步确认 BCC 和 C14 相。室温压缩试样为直径 6 mm、高度 9 mm 的圆柱，每种样品重复试验 3 次，加载速度为 $1\ \mathrm{mm\,min^{-1}}$。

### 2.3 密度结果

表 1 中的实测密度依次为 5.337、5.291、5.599、5.870 和 $5.826\ \mathrm{g\,cm^{-3}}$，均小于 $6\ \mathrm{g\,cm^{-3}}$。最低值对应 $\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$。作者将 Al 添加降低密度与其在五种主元素中质量较轻联系起来。摘要和结论所写的密度范围上限是 $5.826\ \mathrm{g\,cm^{-3}}$，与表 1 的最大值不同。

### 2.4 XRD、BSE 与元素分布结果

第 3.2 节将主要相描述为有序 BCC 相和 C14-Laves 相，并把部分未标注的低强度衍射峰与 $\mathrm{Zr}_5\mathrm{Al}_3$ 析出联系起来。作者观察到提高 Zr 含量时出现额外 C14 衍射峰，低 Al、低 Zr 样品的 BCC 主峰最强。

BSE 图显示亮相常位于暗相边界附近，亮相中可见孔隙等铸造缺陷。低 Al 样品的亮相呈细长透镜状；高 Al 样品出现细小块状相。$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$ 样品中可见作者用红圈标出的裂纹。

作者依据 EDS 分布讨论 Al、Zr 的共同富集以及 V 的分布，并将部分 V 缺失、Al/Zr 富集区域与 $\mathrm{Zr}_5\mathrm{Al}_3$ 联系起来。作者还认为 Al 蒸发影响实际成分、局部浓度梯度及铸造孔隙。表 3 的二元混合焓被用于解释元素分配：Al–Zr 为 $-44\ \mathrm{kJ\,mol^{-1}}$，Nb–Zr 为 $+4\ \mathrm{kJ\,mol^{-1}}$，V–Zr 为 $-4\ \mathrm{kJ\,mol^{-1}}$。这些是原文采用的解释，不等于本文直接测量了凝固反应过程或扩散动力学。

### 2.5 EBSD 相含量与 TEM 相鉴定

| 样品简称 | BCC 含量，% | Laves 含量，% |
|---|---:|---:|
| $\mathrm{Al}_{1.5}$-Zr | 70.8 | 29.2 |
| $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ | 93.9 | 6.1 |
| Al-Zr | 76.9 | 23.1 |
| $\mathrm{Al}_{0.5}$-Zr | 98.4 | 1.6 |
| $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ | 99.3 | 0.7 |

数据来源：原文表 4。作者说明，由于部分晶内球状组织可能是同为六方结构的 $\mathrm{Zr}_5\mathrm{Al}_3$，实际 Laves 相含量可能略低于 EBSD 统计值。因而不能将表 4 的数值不加限定地视为经独立确认的纯 C14 体积分数。

TEM 重点比较低 Al 的样品 4 与样品 5，SAED 被用于将第二相鉴定为 ZrAlV 型 C14-Laves 相，采用的结构卡片为 PDF #05-0312。作者报告降低 Zr 含量后，C14 相平均宽度由约 1 μm 减小至 0.5 μm。

### 2.6 室温压缩结果及作者解释

作者认为 $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$ 在这五个样品中表现出最佳压缩强度与变形能力配合：峰值工程应力为 1783 MPa，工程应变为 28.8%，对应最低的 EBSD Laves 相含量 0.7%。这里的 1783 MPa 是峰值工程应力／压缩强度，不能改称为屈服强度或拉伸强度。

$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$ 在达到屈服阶段之前开裂，作者将其与已有铸造裂纹联系起来。作者还讨论了 Laves 相边界附近位错聚集、应力梯度和沿晶开裂，用以解释相分布对压缩性能的影响。本文未给出独立的位错动力学计算或原位变形跟踪实验来进一步量化这一解释。

### 2.7 C14 结构模型与键级数据

作者使用 CrystalMaker，依据 PDF #05-0312 重建图 7 的 ZrAlV 三维结构。按原文模型，Zr 占据 `4f` 位置，Al 和 V 随机分布于 `2a`、`6h` 位置。作者以 Zr、Al、V 的原子半径和与相关二元结构的相似性解释相稳定性，并得出 Zr 对该 Laves 相形成具有更主导作用的观点。

表 5 报告六种原子对的平均键级。ZrAlV 中 V–V 的值为 0.602，而等原子比 AlNbTiVZr 中为 0.331；Zr–Al、Zr–V 在 ZrAlV 中的值则较低。作者据此讨论相形成和稳定性。本文未交代这些键级的计算软件、计算模型或电子结构参数，不能补写为某种已确定的 DFT 计算流程。

### 2.8 结论的适用范围

原文结论限于所研究的五种铸态样品：Al 和 Zr 都影响组织与相含量；作者认为 Zr 对 ZrAlV 型 C14 相具有更主导作用；低 Al、低 Zr 样品取得最佳室温压缩性能。本文没有据此验证其他配比、其他加工状态或高温服役条件下也具有相同最优行为。

## 三、指定五个部分的完整中文翻译

以下翻译覆盖摘要、第 1 节引言、第 2 节实验过程、第 3 节结果与讨论（包括 3.1–3.4 全部子节）和第 4 节结论。这五部分在原文均有明确对应位置，无须用相似段落替代。原文中的“可能”“似乎”“表明”等判断强度予以保留；译者核对说明置于独立章节。

### 摘要

原文位置：PDF 第 1 页，Abstract。

为研究 Al 和 Zr 元素含量对 AlNbTiVZr 系列轻质耐火高熵合金（HEAs）微观组织演变的综合影响，本文研究了 5 个样品。不同成分的样品分别命名为 $\mathrm{Al}_{1.5}\mathrm{NbTiVZr}$、$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$、$\mathrm{AlNbTiVZr}$、$\mathrm{AlNbTiVZr}_{0.5}$ 和 $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$。结果表明，所研究 HEA 样品的实际密度范围为 $5.291$ 至 $5.826\ \mathrm{g\,cm^{-3}}$。这些 HEAs 的微观组织包含具有 BCC 结构的固溶体相和 Laves 相。通过 TEM 观察，进一步将 Laves 相鉴定为 ZrAlV 金属间化合物。AlNbTiVZr 系列 HEAs 的微观组织同时受到 Al 和 Zr 元素含量的影响；由于 Zr 原子占据 ZrAlV Laves 相（C14 结构）的核心位置，Zr 元素表现出更主导的作用。因此，铸态 $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$ 样品表现出最佳室温压缩性能，其压缩强度（$\sigma_p$）为 1783 MPa，工程应变为 28.8%；这是因为它具有最低的 ZrAlV 金属间化合物面积分数（0.7%），该面积分数由 EBSD 技术表征得到。

关键词：AlNbTiVZr 系列 HEAs；元素含量变化；微观组织演变；ZrAlV Laves 相。

### 1. 引言

原文位置：PDF 第 1–2 页，第 1 节。

Yeh 于 2004 年率先将高熵合金（HEAs）定义为至少含有 5 种主元素、每种元素原子百分比介于 5% 与 35% 之间的合金 [1,2]。HEAs 不同于传统合金，后者通常含有 1 种或 2 种作为所谓溶剂的元素。近期研究表明，由 4 种主元素组成的合金也可以归类为 HEAs [3-5]。鸡尾酒效应会产生不可预测的协同混合结果，是 HEAs 的四大核心效应之一。这为设计复杂合金提供了创新途径；与传统合金相比，这些合金具有高机械强度、耐高温、耐腐蚀和耐辐照等优异性能 [6-9]。尤其是，由耐火元素组成的 HEAs 常被认为是高温服役条件下静态结构部件的潜在结构材料；由于迟滞扩散效应，它们可能替代在高温服役条件下起稳定结构作用的传统 Ni 基高温合金。

一般而言，耐火元素是指熔点超过 1650 °C 的金属元素，包括 W、Mo、Nb、V、Ta 和 Zr。含耐火元素的 HEAs 因其固有的良好高温力学性能而引起研究者的兴趣。Senkov [10] 开发的 VNbMoTaW HEAs 即使在 1600 °C 下仍保持 477 MPa 的屈服强度，他是耐火 HEAs 研究的先驱。

NbTiVZr 系列 HEAs 也由 Senkov 于 2013 年首次开发。这一系列合金在 600 °C 下表现出 834 MPa 的屈服强度，而密度为 $6.52\ \mathrm{g\,cm^{-3}}$，使其比 VNbMoTaW HEAs 轻近 50% [11]。Liao [12] 结合虚晶近似、CALPHAD 建模和准谐 Debye–Grüneisen 模型，深入研究了 NbTiVZr HEAs 的晶体结构，并得出结论：与 FCC 和 HCP 结构相比，BCC 结构是最稳定的晶体结构。元素成分是影响 HEAs 微观组织和性能的关键因素。将真空电弧熔炼（VAM）法制备的铸态 $\mathrm{NbTiVZr}_{0.5}$ 合金与 NbTiVZr 和 $\mathrm{NbTiVZr}_2$ 合金比较可见，增加 Zr 元素含量会使单一 BCC 结构固溶体相之外出现额外的 Laves 相，这一点得到 XRD 试验的确认 [13]。为进一步降低耐火 HEAs 的密度，利用 Al 的低密度将其引入 NbTiVZr 系列 HEAs；$\mathrm{Al}_x\mathrm{NbTiVZr}$（$x=0.5,1,1.5$）合金的实际密度范围为 $5.5$ 至 $6.04\ \mathrm{g\,cm^{-3}}$ [14]。然而，由于 Al 与其他耐火金属元素之间高度负的混合焓所带来的成键能力，这些 AlxNbTiVZr HEAs 的相组成变得更加复杂。而且，与具有单一 BCC 结构的 NbTiVZr HEAs 相比，Laves 相的存在明显削弱了其塑性。压缩试验表明，$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}$ 合金在达到屈服阶段之前发生断裂。Qiu 等 [15] 利用密度泛函理论，通过定义形成热，揭示了 Al–Zr 键有助于 AlNbVTiZr 合金有序构型的结构稳定性。因此，一些研究关注 Zr 元素含量对 AlNbTiVZr 系列 HEAs 结构和力学性能的影响。降低 Zr 含量后，$\mathrm{AlNbTiVZr}_{0.5}$ 合金在室温下表现出 1430 MPa 的优异屈服强度，且只发现表面裂纹，即使圆柱试样被压缩至原高度的一半也是如此。在 BCC 晶粒内部观察到的部分独立 $\mathrm{Zr}_2\mathrm{Al}$ 型 Laves 相颗粒，可能有助于其良好的强度–延性平衡 [16]。然而，进一步降低 Zr 含量时，由于 Al 和 Zr 元素在晶界处偏聚以及 $\mathrm{Zr}_5\mathrm{Al}_3$ 相的形成，$\mathrm{AlNbTiVZr}_{0.25}$ HEA 的压缩应变下降至 6% [17]。

现有研究表明，AlNbTiVZr 系列 HEAs 的优异性能和理想微观组织特征与其成分密切相关。本研究采用真空电弧熔炼工艺制备了 5 个具有不同 Al 和 Zr 含量的样品。这些样品分别命名为 $\mathrm{Al}_{1.5}\mathrm{NbTiVZr}$、$\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$、$\mathrm{AlNbTiVZr}$、$\mathrm{AlNbTiVZr}_{0.5}$ 和 $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$，下标表示原子比。为筛选 AlNbTiVZr 系列 HEAs 的最佳 Al 和 Zr 成分，研究了 Al 和 Zr 元素含量对这 5 个铸态轻质 HEAs 的微观组织演变及相关力学性能的影响。

### 2. 实验过程

原文位置：PDF 第 2–3 页，第 2 节。

本研究用于制备 HEA 的原材料为纯度超过 99.9% 的块状或棒状纯金属。去除表面氧化层并随后进行 5 min 超声振荡后，使用精度为 0.0001 g 的电子天平对这些纯金属进行精确称量，以确保每种主元素的质量误差为 ±0.001 g。各轻质 HEA 样品的名义成分（at.%）见表 1。铸态 HEA 样品采用真空电弧熔炼炉（DHL400，Sky Technology Development Co., Ltd., Chinese Academy of Science，中国沈阳）制备，通入纯 Ar 作为保护气体，使炉内压力保持在 $8\times10^4\ \mathrm{Pa}$。先熔化预留的纯 Ti 锭以吸收残余氧，随后在铜坩埚中对样品反复熔炼并机械翻转 7 次，以确保各主元素均匀分布。

**表 1. 轻质 AlNbTiVZr HEAs 的详细元素组成与密度。**

| 样品 | 合金（简称） | 名义密度，$\mathrm{g\,cm^{-3}}$ | 实际密度，$\mathrm{g\,cm^{-3}}$ |
|---|---|---:|---:|
| 1 | $\mathrm{Al}_{1.5}\mathrm{NbTiVZr}$（$\mathrm{Al}_{1.5}$-Zr） | 5.482 | 5.337 |
| 2 | $\mathrm{Al}_{1.5}\mathrm{NbTiVZr}_{0.5}$（$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$） | 5.345 | 5.291 |
| 3 | $\mathrm{AlNbTiVZr}$（Al-Zr） | 5.739 | 5.599 |
| 4 | $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}$（$\mathrm{Al}_{0.5}$-Zr） | 6.049 | 5.870 |
| 5 | $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}_{0.5}$（$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$） | 5.975 | 5.826 |

首先，从铸态 HEA 铸锭中心区域采用电火花线切割加工边长为 10 mm 的立方体样品。然后，使用 800# 至 3000# 的 SiC 砂纸研磨这些立方体样品，并用金刚石研磨膏抛光 10 min。随后，按照 GB/T 1423-1996，采用阿基米德法对每种 HEA 样品的实际密度测量 3 次。样品采用 RIGAKU 衍射仪进行 X 射线衍射（XRD，Ultima IV，Rigaku Corporation，日本东京）分析，使用 Cu Kα 辐射，以确定相组成和结构。原始 XRD 数据采用 Rietveld 方法处理，随后依据已发表研究和相关 PDF 数据库卡片，在 Jade 6.5 软件中进行校准。使用配备电子背散射衍射（EBSD）和能量色散谱（EDS）探测器的扫描电子显微镜（SEM，JEOL JSM 7001F，JEOL Ltd.，日本东京），观察这些轻质 HEA 样品的微观组织、识别相含量并确定元素分布。此外，切取直径为 3 mm 的圆片样品，并机械减薄至 60 μm。这些样品在 -15 °C、30–35 V 条件下，使用 9% $\mathrm{HClO}_4$–甲醇溶液（体积分数）进行双喷电解抛光，为透射电子显微镜（TEM，JEM-2010，JEOL Ltd.，日本东京）观察做准备。EBSD 试验还在 TEM 样品透光区域附近的无应力区进行。压缩试验是评价材料力学性能的一种常用方法。因此，按照 ASTM D695，使用万能试验机（MTS 858，MTS Systems，美国明尼苏达州），以 $1\ \mathrm{mm\,min^{-1}}$ 的速度，对直径为 6 mm、高度为 9 mm 的圆柱形 HEA 样品进行 3 次压缩试验。

### 3. 结果与讨论

原文位置：PDF 第 3–10 页，第 3 节。

#### 3.1 AlNbTiVZr HEAs 的密度

原文位置：PDF 第 3 页。

密度是材料的一项基本物理属性，表征其占据空间的程度，并与物质的原子质量和晶格排列密切相关。考虑到这些轻质 HEAs 是无序固溶体，各样品的名义密度（ρ）采用式 (1) 计算。

$$
\rho=\frac{\sum c_i A_i}{\sum \frac{c_i A}{\rho_i}}\tag{1}
$$

其中，$c_i$、$A_i$ 和 $\rho_i$ 分别为这些 HEA 样品中第 $i$ 种主元素的原子百分比、原子量和密度。如表 1 所示，这些样品的实际密度均低于 $6\ \mathrm{g\,cm^{-3}}$，并与名义密度一致。$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ 样品的密度最低，为 $5.291\ \mathrm{g\,cm^{-3}}$。加入 Al 元素可能降低 AlNbTiVZr 系列 HEAs 的实际密度，因为 Al 是五种主元素中最轻的元素。

#### 3.2 AlNbTiVZr HEAs 的结构

原文位置：PDF 第 3–6 页。

在五种主元素中，Nb 和 V 单质在整个固态温度范围内均为单一 BCC 结构，而 Zr 和 Ti 单质在高温区间转变为 BCC 相。铸态轻质 AlNbTiVZr HEAs 的 XRD 结果如图 1 所示。这些样品的微观组织主要由有序 BCC 相和 C14-Laves 相组成，与先前研究一致 [14]。同时，在 BCC 相衍射峰之间出现一些未标注的低强度衍射峰，表明存在 $\mathrm{Zr}_5\mathrm{Al}_3$ 析出。与 $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ 和 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品相比，随着 Zr 元素含量增加，$\mathrm{Al}_{1.5}$-Zr 和 $\mathrm{Al}_{0.5}$-Zr 样品显示出额外的 C14-Laves 相衍射峰。就 Al 元素变化而言，$\mathrm{Al}_{1.5}$-Zr、Al-Zr 和 $\mathrm{Al}_{0.5}$-Zr 样品的 BCC 相主衍射峰强度持续增加。在 XRD 实验中，较高的衍射峰强度表明较好的结晶性能和更有序的晶格排列。Al 和 Zr 元素含量最低的 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品表现出最强的 BCC 结构衍射峰。有理由认为，这些 AlNbTiVZr 合金的 C14-Laves 相是富 Al–Zr 相。

![图 1](/imgs/2026-10-09/fig1.png)

**图 1. AlNbTiVZr 轻质 HEA 样品的 XRD 结果。**

图 2 展示了 AlNbTiVZr 样品的背散射电子（BSE）模式微观组织图像。总体而言，值得注意的是，亮相倾向于在暗相边界附近形成；如低倍图像所示，亮相内部存在孔隙等铸造缺陷。在高倍图像中，亮相表现出不均匀的组成：在低 Al 元素含量的 $\mathrm{Al}_{0.5}$-Zr 和 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品中，它呈细长透镜状结构；而在 Al 元素含量增加的 $\mathrm{Al}_{1.5}$-Zr 样品中，它转变为细小块状结构。如图 1 所示，$\mathrm{Al}_{1.5}$-Zr 样品具有大量亮相和胞状暗相。在高倍图像中，可以观察到细小块状相随机分布于两个暗相之间的区域。随着 Zr 元素含量降低，亮相含量低于 $\mathrm{Al}_{1.5}$-Zr 样品。细小块状相的平均尺寸从约 1 μm 增加至 2 μm。然而，沿块状亮相可见红圈标出的部分裂纹，表明 $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ 样品的延性欠佳。对于等原子比 Al-Zr 样品，随着 Al 元素含量降低，细小块状亮相连接在一起。图 2d,f 的高倍图像显示，暗相和亮相主要并排分布。由于 VAM 装置的水冷铜坩埚具有较快冷却速率，$\mathrm{Al}_{0.5}$-Zr 和 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品中的亮相沿传热方向分布；这是因为，在固液界面形核的相最终形成细长透镜状形貌，存在充足的成分过冷液相区域以散热，而在凝固过程中没有胞状共晶空间或组分限制。然而，对于 Zr 元素含量较低的 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品，这些细长透镜状相趋于共格。采用 Image J 软件（1.53g）测量，与 Al-$\mathrm{Zr}_{0.5}$ 样品相比，Al-Zr 和 $\mathrm{Al}_{0.5}$-Zr 样品的亮相面积分数分别由 35.5% 增加至 52.7% 和 61.4%。总体而言，随着 Zr 元素含量增加，暗相倾向于被亮相分隔，晶界逐渐呈现带胞状结构的球形，表明 Al 和 Zr 两种元素均影响 AlNbTiVZr HEAs 的相结构。

![图 2](/imgs/2026-10-09/fig2.png)

**图 2. AlNbTiVZr HEAs 的 BSE 图像：(a) 样品 1，$\mathrm{Al}_{1.5}$-Zr；(b) 样品 2，$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$；(c) 样品 3，Al-Zr；(d) 样品 4，$\mathrm{Al}_{0.5}$-Zr；(e) 样品 5，$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$。**

图 3 和表 2 提供了各 HEA 样品的元素分布及相应实际化学元素含量的信息。如图 3 所示，Al 和 Zr 原子之间的吸引相互作用在 AlNbTiVZr HEAs 的五种主元素中具有最高的亲和性，表明这种特定键中形成固溶体的倾向更大。V 元素的分布似乎均匀，说明 C14-Laves 相与 ZrAlV 金属间化合物相关。此外，在元素分布图（图 3b）中可以观察到一些 V 元素缺失、由 Al 和 Zr 元素填充的区域。因此，先前 XRD 测试中未标注的衍射峰可能归因于 $\mathrm{Zr}_5\mathrm{Al}_3$ 析出。在 VAM 过程中，由于 Al 与其他耐火元素之间的熔点差异相当大，会发生 Al 元素蒸发，导致这些样品的实际 Al 元素含量低于名义成分。值得注意的是，$\mathrm{Al}_{0.5}$-Zr 和 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品的实际成分与名义成分之间的 Al 元素差值，似乎小于其他样品。蒸发在凝固过程中诱导亮相出现局部浓度梯度，使低 Al 元素含量样品中的细长透镜状相沿这一富 Al 相的浓度梯度形成。然而，当提供过量稳定化 Al 原子时，细长透镜状结构在高 Al 元素含量样品中消失。另一方面，Al 元素蒸发也可能在亮色 C14-Laves 相内部引起孔隙等铸造缺陷。从能量角度出发，参考表 3 中 AlNbTiVZr 系列 HEAs 的二元混合焓（$\Delta H_{\mathrm{mix}}$），由于 Nb–Zr 键的 $\Delta H_{\mathrm{mix}}$ 为正值，Zr 元素表现出与最先进入凝固阶段的 Nb 元素强烈分离的倾向。同时，Zr 与另外两种元素（V、Al）混合时具有负焓。正的混合焓意味着化学键抑制有序构型形成，而负值有利于有序构型形成。因此，在 VAM 过程快速冷却引起的非平衡凝固期间，由于扩散时间不足，包晶反应很难完全进行，并发生 Zr 与 Nb 元素偏析以及 Zr、Al、V 元素富集，这有利于保留 Ti 元素的高温 BCC 相。

![图 3](/imgs/2026-10-09/fig3.png)

**图 3. AlNbTiVZr HEAs 的化学元素分布图像：(a) 样品 1，$\mathrm{Al}_{1.5}$-Zr；(b) 样品 2，$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$；(c) 样品 3，Al-Zr；(d) 样品 4，$\mathrm{Al}_{0.5}$-Zr；(e) 样品 5，$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$。**

**表 2. AlNbTiVZr HEA 样品的实际化学元素组成（at.%）。**

| 合金 | Al | Nb | Ti | V | Zr |
|---|---:|---:|---:|---:|---:|
| $\mathrm{Al}_{1.5}$-Zr | 22 | 25 | 15 | 17 | 21 |
| $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ | 25 | 25 | 17 | 19 | 14 |
| Al-Zr | 16 | 27 | 16 | 18 | 23 |
| $\mathrm{Al}_{0.5}$-Zr | 10 | 30 | 17 | 19 | 25 |
| $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ | 11 | 31 | 30 | 22 | 17 |

**表 3. AlNbTiVZr 系列 HEAs 的二元混合焓（$\Delta H_{\mathrm{mix}}$）[18]。**

| $\Delta H_{\mathrm{mix}}$，$\mathrm{kJ\,mol^{-1}}$ | Al | Nb | Ti | V | Zr |
|---|---:|---:|---:|---:|---:|
| Al | - | -18 | -30 | -16 | -44 |
| Nb | - | - | 2 | -1 | 4 |
| Ti | - | - | - | -2 | 0 |
| V | - | - | - | - | -4 |

#### 3.3 AlNbTiVZr HEAs 的压缩性能

原文位置：PDF 第 7 页。

图 4 展示了全部 5 个样品的工程应变–应力曲线，每条曲线在所进行的 3 次压缩试验中均可重复。可以得出，这些 HEA 样品在压缩试验中不能维持较长的屈服阶段；$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品表现出最均衡的强度–延性性能，峰值工程应力（$\sigma_p$）为 1783 MPa，工程应变为 28.8%。尤其是，由于图 2b 中红圈标出的铸造裂纹，$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ 样品在达到屈服阶段之前开裂。随着 Al 元素含量降低，$\sigma_p$ 降低，而工程应变倾向于增加。对于 Zr 元素含量变化，在 Al 元素原子比为 0.5 的较低 Al 含量下，$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品倾向于具有更高强度和更好延性。然而，对于 Al 元素原子比为 1.5 的较高 Al 含量样品，未观察到明显的力学性能变化。这说明 Al 和 Zr 两种元素均影响 AlNbTiVZr HEAs 的力学性能。

![图 4](/imgs/2026-10-09/fig4.png)

**图 4. AlNbTiVZr HEAs 的工程应力–应变曲线。**

#### 3.4 AlNbTiVZr HEAs 的相组成

原文位置：PDF 第 7–10 页。

依据 EBSD 试验，通过 Aztectcrystal 软件获得的反极图（IPF）和相组成图如图 5 所示。可以观察到，BCC 相是这些 AlNbTiVZr HEAs 的主要结构，而 Laves 相广泛分布于晶界。当 Al 元素含量降低时，将 $\mathrm{Al}_{1.5}$-Zr 样品与 Al-Zr 和 $\mathrm{Al}_{0.5}$-Zr 样品比较，Laves 相含量由 29.2% 显著降低至 23.1% 和 1.6%。此外，随着 Al 元素含量降低，Laves 相的宽度减小。图 5d,e 展示了一些具有多种晶体取向的相和细小黑色空区，随后被视为 Laves 相。对于 Al 和 Zr 元素含量最低的 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品，BCC 相含量高达 99.3%。然而，由于 ZrAlV 金属间化合物和 $\mathrm{Zr}_5\mathrm{Al}_3$ 析出物均具有六方结构，该晶内球状组织似乎是 $\mathrm{Zr}_5\mathrm{Al}_3$ 析出物。因此，实际 Laves 相含量略低于表 4 中 EBSD 结果所示的含量。

![图 5](/imgs/2026-10-09/fig5.png)

**图 5. AlNbTiVZr HEAs 的 IPF 和相组成图像：(a) $\mathrm{Al}_{1.5}$-Zr；(b) $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$；(c) Al-Zr；(d) $\mathrm{Al}_{0.5}$-Zr；(e) $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$。**

**表 4. 轻质 HEAs 的相含量。**

| 相含量，% | $\mathrm{Al}_{1.5}$-Zr | $\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ | Al-Zr | $\mathrm{Al}_{0.5}$-Zr | $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ |
|---|---:|---:|---:|---:|---:|
| BCC | 70.8 | 93.9 | 76.9 | 98.4 | 99.3 |
| Laves | 29.2 | 6.1 | 23.1 | 1.6 | 0.7 |

由于 $\mathrm{Al}_{0.5}$-Zr 和 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品表现出较低的 Laves 相含量，因此采用 TEM 技术进一步研究这些样品，以阐明 Zr 元素含量对 AlNbTiVZr HEAs 微观组织演变的影响。如图 6 所示，相应的 SAED 图样确定这些合金包含 BCC 结构相，其标准电子衍射图样沿 $[1\bar{1}\bar{1}]$ 晶带轴；并包含与 PDF #05-0312 卡片一致的 ZrAlV（C14-Laves 结构）金属间化合物。根据 Hume–Rothery 规则，合金只有在其主元素具有相似的晶体结构、原子尺寸和电化学特征时，才能形成简单固溶体。然而，Al 与其他主元素不同，因为 Al 具有 FCC 结构，这不利于形成简单固溶体。Zr、Al 和 V 元素的原子半径分别为 160、143 和 134 pm。此外，Zr 与（Al、V）的原子半径比约为 1.12 至 1.19，与 $AB_2$ Laves 相（C14 结构）的理论最佳形核原子半径比（$R_A/R_B$）1.225 接近。Zr 与（Al、V）原子的半径差异改善了晶格堆积并增强了它们之间的相互作用，最终促进 ZrAlV Laves 相的稳定性。

![图 6](/imgs/2026-10-09/fig6.png)

**图 6. AlNbTiVZr HEAs 的 TEM 图像及相应 SAED 图样：(a) 样品 4，$\mathrm{Al}_{0.5}$-Zr；(b) 样品 5，$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$。**

这种具有较高硬度和熔点的金属间化合物通常表现出较差的变形能力 [19,20]。在压缩过程中，六方结构 ZrAlV 相的滑移系容易被激活，随后在变形初期位错倾向于聚集到 Laves 相边界。随着压缩过程进行，位错逐渐在晶界附近堆积，造成晶内局部应力梯度。这种应力梯度可能导致沿 ZrAlV 金属间化合物发生沿晶开裂，最终降低这些样品的压缩性能。如图 6 所示，BCC 和 C14-Laves 结构表现出近乎平行的分布。随着 Zr 元素含量降低，C14-Laves 结构的平均宽度由 1 μm 减小至 0.5 μm，验证了此前 EBSD 试验的发现。BCC 相表现为不连续块状分布，被 ZrAlV 金属间化合物分隔。因此，位错需要更大的动力才能滑移，并在尖锐的 Laves 相边界处聚集，从而使 $\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品具有更长的弹性变形阶段。此外，BCC 相之间较小的间距使塑性恶化，令其在达到屈服阶段时易于断裂，这一点由压缩试验期间应力迅速降低得到证明。总体而言，在特定 Al 元素含量下，降低 Zr 含量会提高 AlNbTiVZr 系列 HEAs 的变形能力。

本研究的结果表明，Al 和 Zr 两种元素均显著影响 AlNbTiVZr 轻质 HEAs 的金属间相含量。而且，结果暗示 Zr 元素对 C14-Laves 结构（ZrAlV 相）的稳定具有更主导的影响。如图 7 所示，依据 PDF #05-0312 卡片，采用 CrystalMaker 重建了 XRD 结果中观察到的 C14-Laves 相的三维结构。Zr 原子占据中心 `4f` Wyckoff 位置，而 Al 和 V 原子随机分布在 `2a` 和 `6h` Wyckoff 位置。考虑到 ZrAlV 金属间化合物的原子结构与 $\mathrm{ZrAl}_2$ 和 $\mathrm{ZrV}_2$ 相相似，这一三元相的热力学行为可与这些二元相相比。一般而言，当两种主元素的原子数接近化学计量比时，Laves 相形成的可能性增加。根据 $\mathrm{ZrAl}_2$ 和 $\mathrm{ZrV}_2$ 二元相图 [21,22]，显然这些二元金属间化合物的凝固温度范围随 Zr 元素减少而变窄。因此，大量 Zr 基金属间化合物不能充分粗化，形成细小的 ZrAlV 颗粒；其中，二元合金原先稳定的 `2a` 和 `6h` Wyckoff 位置由 Al 和 V 原子替代。综上，与 Al 元素相比，中心 Zr 原子在 ZrAlV 金属间化合物形成中发挥更主导的作用。

![图 7](/imgs/2026-10-09/fig7.png)

**图 7. ZrAlV 金属间化合物的三维原子结构。**

此外，键级是决定原子间键强度的重要因素。键级概念由分子轨道理论引入，等于成键轨道电子数与反键轨道电子数之差的一半。根据结构化学中的化学键理论，较大的正键级值表明相邻原子之间的键更强、更稳定，更难断裂，从而提高构型的稳定性。相反，负键级值表明相邻原子之间的键不稳定，抑制结构构型的稳定性。ZrAlV Laves 相和等原子比 AlNbTiVZr HEA 中 Zr、Al、V 原子之间的键级列于表 5。Zr、Al、V 原子之间的所有键级均为正值。V–V 键的键级由 0.331 增加至 0.602，提供促进相形成并维持 ZrAlV 金属间化合物结构稳定性的结合动力。此外，与 AlNbTiVZr HEA 固溶体相比，Laves 相中 Zr–Al 和 Zr–V 键的键级降低，表明 Zr 元素对稳定 ZrAlV Laves 相具有主导影响。

**表 5. ZrAlV Laves 相和 AlNbTiVZr 合金的平均键级。**

| 平均键级 | Zr–Zr | Al–Al | V–V | Zr–V | Zr–Al | Al–V |
|---|---:|---:|---:|---:|---:|---:|
| ZrAlV Laves 相 | 0.227 | 0.285 | 0.602 | 0.228 | 0.217 | 0.360 |
| AlNbTiVZr 合金 | 0.354 | 0.236 | 0.331 | 0.273 | 0.243 | 0.255 |

### 4. 结论

原文位置：PDF 第 11 页，第 4 节。

本研究考察了 Al 和 Zr 元素含量对 AlNbTiVZr 系列 HEAs 微观组织演变及相关力学性能的综合影响。得出以下结论。

1. 所研究 AlNbTiVZr 系列 HEAs 的实际密度均低于 $6\ \mathrm{g\,cm^{-3}}$，范围为 $5.291$ 至 $5.826\ \mathrm{g\,cm^{-3}}$。这些合金可归类为轻质 HEAs。

2. 这些 AlNbTiVZr HEAs 主要表现出 BCC 结构和 C14-Laves 相，后者进一步被确认为 ZrAlV 金属间化合物。Al 和 Zr 元素含量均影响微观组织演变；由于 Zr 原子占据 ZrAlV Laves 相的中心位置，Zr 元素对相组成具有更主导的影响。将 $\mathrm{Al}_{1.5}$-Zr 样品与 Al 元素含量较低的 Al-Zr 和 $\mathrm{Al}_{0.5}$-Zr 样品比较，Laves 相含量由 29.2% 明显降低至 23.1% 和 1.6%。随着 Al 元素含量降低，Laves 相宽度减小。

3. 脆性 Laves 相容易聚集并分隔 BCC 结构，显著影响 AlNbTiVZr HEAs 的压缩力学性能。值得注意的是，由于具有最低的 ZrAlV Laves 相含量（0.7%），$\mathrm{Al}_{0.5}$-$\mathrm{Zr}_{0.5}$ 样品表现出最佳压缩力学性能，其压缩强度（$\sigma_p$）为 1783 MPa，工程应变为 28.8%。相比之下，$\mathrm{Al}_{1.5}$-$\mathrm{Zr}_{0.5}$ 样品因较高 Al 元素含量诱发的凝固裂纹而未能达到屈服阶段。

## 四、全文计算与数据处理细节汇总

本节汇总原文中出现的计算和定量处理，不补写原文没有交代的方法。

### 4.1 名义密度计算

位置：第 3.1 节，PDF 第 3 页，式 (1) 与表 1。

原文采用各主元素的原子百分比 $c_i$、原子量 $A_i$ 和密度 $\rho_i$ 计算样品的名义密度。原文式 (1) 已在译文中照录。原文没有列出逐项代入的元素原子量、元素密度、输入组成表或计算代码。表 1 给出的是最终名义密度和实际密度。实际密度采用阿基米德法测量，重复 3 次；正文未列出测量分散程度或误差条。

式 (1) 分母印为 $c_i A/\rho_i$，其中 $A$ 未带下标；随后的变量定义只解释了 $A_i$。此处未替作者补下标，详见核对备注。

### 4.2 XRD 数据处理

位置：第 2 节，PDF 第 3 页；结果见第 3.2 节及图 1。

原文报告：使用 Cu Ka 辐射获取 XRD 数据，对原始数据应用 Rietveld 方法，依据已有研究与 PDF 数据库卡片在 Jade 6.5 中校准。原文未报告扫描角度范围的文字参数、步长、计数时间、仪器展宽模型、精修拟合指标或完整精修文件。图 1 中的坐标范围是图示信息，不在本报告中补写为已明确给出的实验设置。本文也未提供 B2 长程有序参数的计算式或数值。

### 4.3 亮相面积分数测量

位置：第 3.2 节，PDF 第 4–5 页。

作者使用 Image J 1.53g 测量 BSE 图像中的亮相面积分数，正文报告 35.5%、52.7%、61.4%。相关比较中的样品名称原文存在不清楚之处，已按原文翻译并在备注中列出。原文未给出图像阈值、视场数、总统计面积或统计误差。BSE 亮相面积分数与表 4 的 EBSD Laves 相含量来自不同的表征与处理，不在本报告中合并为同一数据。

### 4.4 EBSD 相含量统计

位置：第 3.4 节，PDF 第 7–8 页，图 5 与表 4。

作者使用 Aztectcrystal 软件生成 IPF 和相组成图，并报告两类相的含量。最低 Laves 数值为 0.7%，最高为 29.2%。正文同时指出部分被识别为 Laves 的晶内球状组织可能是 $\mathrm{Zr}_5\mathrm{Al}_3$，实际 Laves 含量略低于 EBSD 结果。原文未给出步长、索引质量阈值、相分类判据、误索引修正流程或误差评估。

### 4.5 压缩试验中的定量结果

位置：第 2 节、第 3.3 节，PDF 第 3、7 页，图 4。

原文给出圆柱试样尺寸、加载速度、重复次数和试验标准。图 4 为工程应力–应变曲线，作者报告最佳样品的 $\sigma_p=1783\ \mathrm{MPa}$ 和工程应变 28.8%。原文未写出工程应力或工程应变的换算公式，也未报告逐样品力学性能数值表。本报告未从曲线另行读数、拟合或推算未给出的屈服强度。

### 4.6 原子半径比与结构重建

位置：第 3.4 节，PDF 第 8–10 页，图 6、7。

原文给出 Zr、Al、V 原子半径分别为 160、143、134 pm，Zr 与（Al、V）的半径比约 1.12–1.19，并与理论最佳比例 $R_A/R_B=1.225$ 比较。作者据此讨论尺寸差异与结构稳定性。本文没有自行增加原子半径来源、理论半径比推导或稳定性判据。

图 7 使用 CrystalMaker，根据 PDF #05-0312 卡片重建。原文列出 Zr 的 `4f` 占位以及 Al/V 在 `2a`、`6h` 上随机分布的模型描述；未报告结构文件、完整分数坐标、晶格常数、精修占位率或结构优化参数。

### 4.7 二元混合焓

位置：第 3.2 节、表 3，PDF 第 5–7 页。

表 3 的二元混合焓引自参考文献 [18]。作者用 Nb–Zr 的正值和 Zr–Al、Zr–V 的负值讨论元素分离与有序构型形成倾向。本文未给出五元合金整体混合焓的求和公式、整体计算值、自由能计算结果或 CALPHAD 相图计算。

### 4.8 键级的定义、结果与缺失计算细节

位置：第 3.4 节末尾，PDF 第 10 页，表 5。

原文文字定义为：“键级等于成键轨道电子数与反键轨道电子数之差的一半。”原文没有为这个定义单独列编号公式。表 5 给出两种结构中 Zr–Zr、Al–Al、V–V、Zr–V、Zr–Al、Al–V 六类平均键级，全部为正。作者重点比较 V–V 从 0.331 到 0.602 的增加，以及 Zr–V 从 0.273 到 0.228、Zr–Al 从 0.243 到 0.217 的降低。

原文没有说明：键级采用哪一种具体定义或布居分析方法；数据是本研究计算还是取自其他来源；使用何种软件；模型尺寸与成分；交换关联泛函；赝势或基组；能量截断；$k$ 点网格；结构弛豫与电子收敛阈值；平均键级如何按原子对计数。因此不能补写这些参数，也不能据此声称已复现表 5。

引言提及的虚晶近似、CALPHAD、准谐 Debye–Grüneisen 模型和 DFT 属于作者对参考文献 [12]、[15] 的介绍，不应转写为本论文实际采用的计算流程。

## 五、原文核对备注

以下均为对本文文字、图表和数值的核对，未用推测代替作者更正。译文保留原文；这些备注不属于原论文结论。

1. **样品 4 名称不一致。** 摘要及引言最后一段列出 $\mathrm{AlNbTiVZr}_{0.5}$；表 1、图 2、图 3、图 6 及后文讨论则使用 $\mathrm{Al}_{0.5}\mathrm{NbTiVZr}$／$\mathrm{Al}_{0.5}$-Zr。未将两者静默统一。

2. **密度范围上限不一致。** 摘要与结论为 $5.291$ 至 $5.826\ \mathrm{g\,cm^{-3}}$，但表 1 的样品 4 为 $5.870\ \mathrm{g\,cm^{-3}}$。表 1 的五个实际值最高是 5.870。这里仅指出数值不一致，不选择性修改摘要或表格。

3. **表 2 部分原子百分比不归一。** 按表中原数值做算术核对，五行之和依次为 100、100、100、101、111。原文未解释后两行偏离 100 的原因，未提供可用于纠正的原始数据；本报告保持原数值。不能把它们自动归一化后称为作者报告的数据。

4. **密度公式的下标不一致。** PDF 第 3 页式 (1) 分子为 $\sum c_i A_i$，分母中为 $\sum c_i A/\rho_i$，而变量说明只定义了 $A_i$。译文保持可见公式，不擅自补全。

5. **正文图号／分图号存在不对应。** 第 3.2 节 BSE 讨论中出现“如图 1 所示”，以及“图 2d,f”；可见图 2 仅有 (a)–(e)。译文保留这些指向，并在此提示，未推断作者原本想引用哪个分图。

6. **亮相面积分数比较中的样品简称不明确。** 第 3.2 节写有 Al-$\mathrm{Zr}_{0.5}$，但表 1 中的五个简称未包含这一项。译文保留该名称，不据此重新分配 35.5%、52.7%、61.4% 的样品对应关系。

7. **相有序性描述的限定。** 第 3.1 节称这些 HEAs 为无序固溶体，第 3.2 节又称主要组织含有序 BCC 相；本文未提供 B2 长程有序参数或明确的超点阵峰定量分析。本报告如实翻译，不把 BCC 主峰增强额外改写为已定量验证的 B2 有序程度增加。

8. **结构模型与直接占位测量的限定。** `4f`、`2a`、`6h` 的描述来自 CrystalMaker 按数据库卡片重建的模型；原文没有报告独立的元素占位精修结果。将“中心位置”按原文翻译，不另行将它解释为某一几何体心位置。

9. **键级结果的可复现性限定。** 表 5 的六类键级完整保留，但计算设置和数据来源没有明确交代。本报告不补写软件、泛函或布居分析方法。

10. **术语与归因的保留。** 引言把 $\mathrm{Zr}_2\mathrm{Al}$ 型颗粒称为 Laves 相；第 3.2 节使用“共格”表述；第 3.4 节包含由尺寸、二元相图和键级解释主导作用的语句。本报告保留这些原文表述，不增加结构判定、界面关系或凝固动力学的验证结论。

## 六、其他原文声明的翻译

原文位置：PDF 第 11 页。

**作者贡献：** H.Y.：数据整理、研究实施、验证、初稿撰写、项目管理。R.L.：审阅与编辑、数据整理。S.L.：写作、数据整理。Y.Z.：概念构思、审阅与编辑、方法学、指导。W.X.：资源、写作。B.X.（Boyu Xue）：方法学、正式分析。B.X.（Baiqing Xiong）：概念构思、审阅与编辑、指导。X.L.：方法学、资源。Z.L.：方法学、经费获取。所有作者均已阅读并同意稿件的发表版本。

**经费支持：** 本研究获得中国 GRIMAT Engineering Institute Co., Ltd. 科技创新基金项目的资助。

**机构审查委员会声明：** 不适用。

**知情同意声明：** 不适用。

**数据可获得性声明：** 数据包含在本文内。

**利益冲突：** 所有作者均受雇于 China GRINM Group Co., Ltd. 和 GRIMAT Engineering Institute Co., Ltd.。作者声明不存在利益冲突。

## 七、原文引用编号索引

以下按照原论文参考文献的著录列出，用于对应译文中的 [n]。没有补入原文未给出的 DOI，也没有将参考文献中的结论追加到本论文总结。作者姓名、年份、卷页及原文拼写按原文保留；此索引不表示已对所有被引论文全文进行独立核查。

[1] Cantor, B.; Chang, I.T.H.; Knight, P. Microstructural development in equiatomic multicomponent alloys. Mat. Sci. Eng. A-Struct. 2004, 376, 213–218.

[2] Yeh, J.; Chen, S.; Lin, S.; Gan, J.; Chin, T.; Shun, T.; Tsau, C.; Chang, S. Nanostructured high-entropy alloys with multiple principal elements: Novel alloy design concepts and outcomes. Adv. Eng. Mater. 2004, 6, 299–303.

[3] Senkov, O.N.; Miracle, D.B.; Chaput, K.J.; Couzinie, J.P. Development and exploration of refractory high entropy alloys-A review. J. Mater. Res. 2018, 33, 3092–3128.

[4] Tian, Y.; Zhou, W.; Tan, Q.; Wu, M.; Qiao, S.; Zhu, G.; Dong, A.; Shu, D.; Sun, B. A review of refractory high-entropy alloys. Trans. Nonferrous Met. Soc. 2022, 32, 3487–3515.

[5] Wang, Z.; Chen, S.; Yang, S.; Luo, Q.; Jin, Y.; Xie, W.; Zhang, L.; Li, Q. Light-weight refractory high-entropy alloys: A comprehensive review. J. Mater. Sci. Technol. 2023, 151, 41–65.

[6] Youssef, K.; Zddach, A.; Niu, C.; Irving, D.; Koch, C. A novel low-density, high-hardness, high-entropy alloy with close-packed single-phase nanocrystalline structures. Mater. Res. Lett. 2014, 3, 95–99.

[7] Tong, Y.; Chen, D.; Han, B.; Wang, J.; Feng, R.; Yang, T.; Zhao, C.; Zhao, Y.; Guo, W.; Shimizu, Y. Outstanding tensile properties of a precipitation-strengthened FeCoNiCrTi0.2 high-entropy alloy at room and cryogenic temperatures. Acta Mater. 2018, 165, 228–240.

[8] Shi, Y.; Yang, B.; Liaw, P.K. Corrosion-resistant high-entropy alloys: A review. Metals 2017, 7, 344.

[9] Pickering, E.J.; Carruthers, A.W.; Barron, P.J.; Middleburge, S.C.; Armstrong, D.E.J.; Gandy, A.S. High-entropy alloys for advanced nuclear applications. Entropy 2021, 23, 144.

[10] Senkov, O.N.; Wilks, G.B.; Scott, J.M.; Miracle, D.B. Mechanical properties of Nb25Mo25Ta25W25 and V20Nb20Mo20Ta20W20 refractory high entropy alloys. Intermetallics 2011, 19, 698–706.

[11] Senkov, O.N.; Senkova, S.V.; Woodward, C.; Miracle, D.B. Low-density, refractory multi-principal element alloys of the Cr-Nb-Ti-V-Zr system: Microstructure and phase analysis. Acta Mater. 2013, 61, 1545–1557.

[12] Liao, M.; Liu, Y.; Min, L.; Lai, Z.; Han, T.; Yang, D.; Zhu, J. Alloying effect on phase stability, elastic and thermodynamic properties of Nb-Ti-V-Zr high entropy alloy. Intermetallics 2018, 101, 152–164.

[13] King, D.J.M.; Cheung, S.T.Y.; Humphry-Baker, S.A.; Parkin, C.; Couet, A.; Cortie, M.B.; Lumpkin, G.R.; Middleburge, S.C.; Knowles, A.J. High temperature, low neutron cross-section high-entropy alloys in the Nb-Ti-V-Zr system. Acta Mater. 2019, 116, 435–446.

[14] Stepanov, N.D.; Yurchenko, N.Y.; Shaysultanov, D.G.; Salishchev, G.A.; Tikhonovsky, M.A. Effect of Al on structure and mechanical properties of AlxNbTiVZr (x = 0, 0.5, 1, 1.5) high entropy alloys. Mater. Sci. Technol. 2015, 31, 1184–1193.

[15] Qiu, S.; Chen, S.; Naihua, N.; Zhou, J.; Hu, Q.; Sun, Z. Structural stability and mechanical properties of B2 ordered refractory AlNbTiVZr high entropy alloys. J. Alloys Compd. 2021, 886, 161289.

[16] Stepanov, N.D.; Yurchrnko, N.Y.; Slkolovsky, V.S.; Tikhonovsky, M.A.; Salishchev, G.A. An AlNbTiVZr0.5 high-entropy alloy combining high specific strength and good ductility. Mater. Lett. 2015, 161, 136–139.

[17] Yurchenko, N.; Panina, E.; Belyakov, A.; Salishchev, G.; Zherebstov, S.; Stepanov, N. On the yield stress anomaly in a B2-ordered refractory AlNbTiVZr0.25 high-entropy alloy. Mater. Lett. 2021, 311, 131584.

[18] Takeuchi, A.; Inoue, A. Classification of bulk metallic glasses by atomic size difference, heat of mixing and period of constituent elements and its application to characterization of the main alloying element. Mater. Trans. 2005, 46, 2817–2829.

[19] Zhang, X.; Zhang, B.; Liu, S.; Zhang, X.; Liu, R. Microstructure and mechanical properties of novel Zr-Al-V alloys processed by hot rolling. Intermetallics 2020, 116, 106639.

[20] Junior, J.M.M.; Chaia, N.; Cotton, J.D.; Coelho, G.C.; Nunes, C.A. A multi-principal element alloy combining high specific strength and good ductility. Mater. Lett. 2022, 325, 132905.

[21] Wang, Q.; Chen, F.; Wu, J.; Qiang, J.; Dong, C.; Zhang, Y.; Xu, F.; Sun, L. Hydrogen absorption by Laves phase related BCC solid solution. Rare Matels 2006, 25, 252–255.

[22] Hu, Z.; Xu, Y.; Che, Y.; Schützendübe, P.; Wan, J.; Huang, Y.; Liu, Y.; Wang, Z. Anomalous formation of micrometer-thick amorphous oxide surficial layers during high-temperature oxidation of ZrAl2. J. Mater. Sci. Technol. 2019, 35, 1479–1484.

## 八、具体章节翻译覆盖说明

用户本次明确指定的摘要、引言、方法、结果、结论已完整覆盖；结果与讨论的 4 个子节标题均已翻译，原文 7 幅图和 5 张表均随相应段落输出。未另外指定其他章节，因此没有将参考文献中的内容或外部论文正文替代为“指定章节”翻译。
