# 临床终点上的 Global Lasso / SHAP

伪标签验证回答「解释器会不会指出已知生成特征」。临床任务没有这种 ground truth，只能报 **跨 fold 稳定、方向一致** 的特征。下面这些表不要和 motif 回收分数混用。

## Stability Lasso（组成）

`scripts/run_stability_lasso.py`：patient CV、多种子、patient bootstrap、L1 logistic。当前 **不做 survival**（没有 L1 Cox）。

输入可以是 composition、expression、二者拼接、UTAG、CytoCommunity TCN 等。

Top-5（按 `|coefficient_mean|`）见 `tables/clinical_stability_lasso_top5.csv`。

Jackson `response` 上反复出现的方向（组成 Lasso）：

- 正：B cell、Tumor (Apoptotic)、部分 SMA/Vimentin 间质
- 负：Stroma (Fibronectin hi)、Stroma (Large elongated)、Tumor (CK lo HR hi p53+)

多数特征的 bootstrap CI 仍跨 0。报告时应同时给 `seed_selection_frequency` 和 `bootstrap_ci_direction`，不要只报点估计。

## 组成 + 表达的线性 SHAP（Cox）

生存任务用 `--model cox`，SHAP 在 log-risk 空间，精确线性公式与分类相同。

Top-5（`direction_consistency == 1`）见 `tables/clinical_shap_comp_expr_cox_top5.csv`。

Jackson DFS 上 |SHAP| 最高且方向一致的包括：`expression__Vimentin`、`expression__Ecadherin`、`expression__c-erbB-2 - Her2`、`composition__Tumor (CK7+)`。这是全区域均值，不是空间邻域。

## CytoCommunity TCN 的 SHAP

TCN 比例是区域级组成，ID **不能跨数据集、甚至不能跨 region 直接比较**。生物标签来自该次运行里每个 TCN 的细胞类型富集，见：

`tables/clinical_shap_cytocommunity_tcn_biological.csv`

（源文件 `results/shap_cytocommunity_selected_tcn_group_features_ordered_biological.csv`）

筛选逻辑：精选 dataset×task×scheme 上，`direction_consistency == 1`，再按 `mean_abs_shap` 每组取 Top-1。

例子：

| 数据集 | 任务 | 方案 | 特征 |
|---|---|---|---|
| NSCLC-Aung | immunotherapy_response | Yale→YaleExt / Yale→UQ | TCN5: DC/B-cell/M1-enriched tumor niche |
| HNC-Wu | primary_outcome | CV 与 UPMC→DFCI | TCN4: dendritic/NK–proliferative tumor niche |
| Jackson | response | CV | TCN3: proliferative/p53+ tumor–macrophage niche |
| TNBC-Wang | pCR_all | CV | TCN5: neutrophil/B-cell/TCF1+CD4 T niche |

这些是 **临床关联的邻域组成**，不是伪标签验证过的接触 motif。TCN 编号换一次无监督聚类就会变。

## UTAG 临床 Lasso

`results/utag_all_datasets_lasso_feature_summary.csv` 里的列名是 `utag_0000` 这种匿名维。没有 marker/domain 名时，临床 UTAG-Linear 只有预测稳定性，没有可叙述的生物学。解释应回到 message-passing 的 mean/std 列名，或 portraits。

## 临床任务上怎么分工（可以各解释一部分）

伪标签没有给出「合格的全能解释器」，但给出了每家能答什么、不能答什么。临床没有 ground truth，更要按问题拆开用，不要拼成同一个空间故事。

| 临床问题 | 用谁 | 可以写进结果的句子 | 不要写成 |
|---|---|---|---|
| 哪些细胞类型 / marker 均值与终点相关？ | 组成或表达的 Stability Lasso；Cox 线性 SHAP | 「线性模型稳定选中 / 归因到这些命名特征」 | 「已验证的空间接触」 |
| 平滑后的哪些 marker 与终点相关？ | UTAG-Linear（`utag_mean__` / `utag_std__`） | 「message-passing 后的 marker 统计与终点相关」 | 「UTAG domain / 空间分区」 |
| 哪个粗分区 / 细胞邻域组成与终点相关？ | UTAG portraits；CytoCommunity TCN SHAP（带生物标签） | 「某 domain / TCN 画像与终点相关」 | 「CD8 聚类、免疫排斥、细胞接触」 |
| embedding 里有没有结局信号？ | Eva / KRONOS linear probe | 「冻结 embedding 的 CV AUC」 | 任何细胞类型或空间机制 |
| 图模型觉得哪些细胞重要？ | GNN Explainer | 目前不建议作为临床主解释 | 节点归因 = 生物学证据 |

推荐顺序：

1. 先做 **组成 ± 表达的 Lasso / 线性 SHAP**。这是临床解释的主结果：有列名、可报稳定性和方向。分类用 Lasso，生存用 Cox SHAP。
2. 若要补「分区」而不是「全区域均值」，再加 **UTAG portraits 或 TCN**，并写明 TCN ID 不可跨数据集比较。
3. 若要补「哪些平滑 marker」，用 **UTAG-Linear 的 mean/std**，不要和 portraits 捆成一个 UTAG 分数。
4. Eva/KRONOS 只回答有没有信号。GNN Explainer 在伪标签上定位失败，临床不要当主解释。

三条硬限制：

- Lasso/SHAP 在 TNBC `cd8_high` 上会避开组成、改抓 mixing。临床表里如果 Top 特征是 mixing，只能说「线性模型用了邻域统计」，不能说「恢复了空间 motif」。
- 它们在 `immune_exclusion` 上会走肿瘤面积捷径。临床密度/面积特征同样可能只是丰度代理。
- 各方法的 Top 特征对不上时，并列报告，不要互相「验证」。

## 和伪标签验证如何一起写

1. **可信度（伪标签）：** 没有方法同时过选定丰度 control 和空间回收。Lasso/SHAP 能在部分数据集上写出 mixing 列名，但选定 control 是 1/2。UTAG-Linear 能找回 marker。UTAG portraits 选定 control 2/2，空间几乎失败。GNN Explainer 当前不能。
2. **临床发现（无 ground truth）：** 只在上表允许的句子里报告稳定关联，并写明这是关联不是因果、也不是已验证的空间接触。
