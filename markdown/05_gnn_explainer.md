# GNN Explainer（SPACE-GM）

## 方法

GNN Explainer 解释的是 **SPACE-GM 细胞图分类器**，不是表格线性模型。

流程：

1. 每个选定 motif 训练 SPACE-GM（脚本默认只用 **seed 0 / fold 0** 一个 checkpoint）
2. `scripts/export_space_gm_node_importance.py --method gnn-explainer` 导出 `region_id, cell_type, importance`
3. `verify_pseudo_label_explanations.py --mode gnn-explainer`：每个 region 取 importance 最高的 10% 细胞，看预期细胞集合的比例是否相对背景富集 ≥ 1.25，且高于背景
4. 任务通过：≥50% 有标签的 region 过富集

预期细胞集合见 `GNN_EXPECTED_SETS`，例如：

- `tumor_high` → tumor
- `cd8_high` → cd8 或 t_cell
- `cd8_clustering` / `immune_exclusion` → cd8
- `apc_t_contact` → apc + t_cell
- `tumor_stroma_mixing` → tumor + stroma

导出 CSV 已齐（Jackson 299 MB，HNC 255 MB，TNBC 74 MB，METABRIC 50 MB）。本次把验证补跑完，写入 `results/pseudo_label_explanations_panel/gnn_explainer/`。

## 结果

表：`tables/gnn_explainer_task_summary.csv`

| 数据集 | 任务 | 类型 | 评分数 | 均值富集 | region 通过率 | 阈值 0.50 |
|---|---|---|---:|---:|---:|---|
| Jackson | tumor_high | control | 376 | 1.28 | 0.319 | fail |
| Jackson | t_tumor_mixing | spatial | 187 | 0.92 | 0.037 | fail |
| Jackson | cd8_tumor_contact | spatial | 188 | 1.10 | 0.202 | fail |
| Jackson | macrophage_tumor_niche | spatial | 194 | 0.92 | 0.062 | fail |
| Jackson | apc_t_contact | spatial | 180 | 1.12 | 0.300 | fail |
| TNBC | cd8_high | control | 766 | 1.40 | 0.477 | fail |
| HNC | cd8_clustering | spatial | 172 | 1.34 | 0.459 | fail |
| HNC | immune_exclusion | spatial | 175 | 1.27 | 0.463 | fail |
| METABRIC | interface_immune | spatial | 273 | 1.04 | 0.319 | fail |
| METABRIC | tumor_stroma_mixing | spatial | 307 | 0.99 | 0.010 | fail |

选定 control **0/2**，空间 raw **0/8**。协议分因此也是 0。

CSV 里没有同数据集的另一条 control 时，验证会 skip（例如 Jackson 未导出 `cd8_high`，HNC/METABRIC 未导出任何 control）。这是导出任务表造成的，不是验证脚本丢了文件。

## 怎么理解这些数字

均值富集 &gt; 1 只说明 **平均而言** top 节点里预期细胞略多，不表示多数 region 都过关。通过规则是逐 region 的：`enrichment ≥ 1.25` 且 `top_hit > background`。大量 region 略高于 1.0 就会把均值拉起来，通过率仍低于 0.5。

最接近阈值的是：

- TNBC `cd8_high`：0.477
- HNC `cd8_clustering`：0.459
- HNC `immune_exclusion`：0.463

空间接触类任务（T–肿瘤混合、巨噬细胞龛、肿瘤–间质混合）接近或低于背景，说明 importance 没有落在接触双方上。

## 限制（报告时必须写）

1. **单 checkpoint。** 不是 15 fold 平均。Lasso/SHAP 的 15 fold 稳定性这里没有。
2. **训练标签是 SPACE-GM 的 noisy v2 设定**，与表格面板的 `pseudo_labels/*_v2.csv` 不完全同一条流水线。
3. 没有 faithfulness（没有按节点 mask 做 MORF/LERF）。只有定位富集。
4. 细胞类型名要映射到 motif catalog 的 cell set；命名不一致会低估富集。

可选后续：`ALL_FOLDS=1` 平均 15 个 checkpoint；或对同一批图跑 IG / occlusion（`run_space_gm_gnn_explainer.sh` 里已预留）。

## 结论

以当前协议，**GNN Explainer 不能声称恢复了 motif 生成细胞。** 它在 CD8 相关任务上有弱定位倾向（通过率 ~0.46–0.48），在混合/接触任务上接近无效。和 Lasso/SHAP 的命名特征回收不在同一成熟度。
