# 相应方法：UTAG portraits、KRONOS、Eva、UTAG domains

这些方法用来对照「有名字的线性解释」和「只有向量、没有名字的 embedding」。

## 1. UTAG native domain portraits

面板主 UTAG 解释器。**不再**把 domain 比例送进 Lasso/SHAP。

每个数据集拟 10 个 domain，给每个 domain 一张画像：多数细胞类型 + centroid 上最高的 marker。再看与 motif 标签关联最强的 domain，是否命中预期细胞/marker，且空间任务不能只是纯肿瘤丰度。

表：`tables/utag_native_portraits_all_tasks.csv`，`tables/utag_cellular_selected_task_verdicts.csv`

选定任务 raw：

| 数据集 | 任务 | 最强 domain | 多数类型 | 画像 AUC | raw |
|---|---|---|---|---:|---|
| Jackson | tumor_high | Tumor (CK lo HR lo)__06 | 肿瘤 | 0.748 | **pass** |
| Jackson | 四个空间 motif | 几乎都是肿瘤 domain | 肿瘤 | 0.52–0.56 | fail |
| TNBC | cd8_high | CD79a^+Plasma__07 | 浆细胞，marker 含 CD8/CD3 | 0.846 | **pass** |
| HNC | cd8_clustering / immune_exclusion | Tumor__08 | 肿瘤 | 0.64 | fail（miss=肿瘤） |
| METABRIC | interface_immune | Myofibroblasts__01 | 肌成纤维 | 0.647 | pass（raw） |
| METABRIC | tumor_stroma_mixing | HR+ CK7-__09 | 肿瘤 | 0.532 | fail |

选定 control **2/2**（Jackson tumor_high，TNBC cd8_high）。空间 raw **1/8**（仅 METABRIC interface_immune）。HNC/METABRIC 的匹配 control 未同时过，协议空间 **0/4**。

结论：UTAG domain 能抓住「这块组织更像肿瘤 / 更像淋巴」这种粗分区，抓不住残差化之后的接触和聚类。空间任务的 top domain 经常就是肿瘤，和 abundance-only 失败模式相同。

## 2. KRONOS linear probe

冻结 KRONOS region embedding，Ridge（默认 `C=1, L2`）做 patient CV。没有特征名，只报 AUC。

Control 阈值 0.90，空间 0.60。

表：`tables/kronos_embedding_probe.csv`

选定任务：

| 任务 | AUC | raw |
|---|---:|---|
| Jackson tumor_high | 0.839 | fail control（&lt;0.90） |
| TNBC cd8_high | 0.924 | **pass control** |
| Jackson 四个空间 | 0.54–0.64 | 接触/混合各过一条，APC 和巨噬龛不过 |
| HNC clustering / exclusion | 0.613 / 0.706 | pass |
| METABRIC mixing / interface | 0.597 / 0.613 | mixing 不过，interface 过 |

空间 raw **5/8**。因 Jackson control 失败，且多数数据集两个 control 不能同时 ≥0.90，协议空间 **0/8**。

KRONOS 能区分部分空间标签，但不能指出是 mixing 还是 density，也不能过严格 control 门。

## 3. Eva linear probe

与 KRONOS 相同协议，换 Eva embedding。

表：`tables/eva_embedding_probe.csv`

选定任务空间 AUC：

| 任务 | AUC |
|---:|---:|
| apc_t_contact | 0.715 |
| cd8_tumor_contact | 0.804 |
| macrophage_tumor_niche | 0.666 |
| t_tumor_mixing | 0.812 |
| cd8_clustering | 0.616 |
| immune_exclusion | 0.721 |
| interface_immune | 0.754 |
| tumor_stroma_mixing | 0.888 |

空间 raw **8/8**（全部 ≥0.60）。Control：Jackson tumor_high 0.939 过，TNBC cd8_high 0.873 不过。协议空间仍是 **0/8**。

在「有没有空间信号」上，Eva 是 embedding 里最强的；在「解释对不对」上，它和 KRONOS 一样交白卷。

## 4. UTAG domains + Lasso/SHAP（仅 HNC 旧跑）

`results/pseudo_label_explanations/utag_domains/` 把 fold 内 domain 比例当表格特征，再用 tabular/UTAG 规则回收。

HNC 上：

- tumor_high：AUC 0.925，回收 1.00，control 过
- cd8_high：AUC 0.815，回收 0.13，control 不过
- 四个空间：AUC 0.49–0.73，回收低，判 abundance_only_or_fail

这再次说明：domain 比例更像组成向量，不是接触描述子。面板因此改成 native portraits，并且 **禁止** 把「UTAG + Lasso」当成一个解释器报分。

## 对照小结

| 问题 | 该看谁 |
|---|---|
| 选定丰度 control 过了没有？ | UTAG portraits 是 2/2；Lasso/SHAP 是 1/2 |
| 同数据集 control 过关时，空间机制能不能被点名？ | Lasso / SHAP 命名特征（Jackson / HNC） |
| 平滑 marker 是否携带丰度？ | UTAG-Linear embedding |
| 空间分区画像是否像预期细胞？ | UTAG portraits |
| 黑盒 embedding 有没有空间 AUC？ | Eva（强于 KRONOS） |
| 图模型节点是否落在对的细胞上？ | GNN Explainer（当前否） |
