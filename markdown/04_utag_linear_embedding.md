# UTAG-Linear embedding

## 方法是什么

UTAG 原论文是：细胞 × marker 做空间 message passing，再聚类成 domain。  
**UTAG-Linear embedding** 停在第一步之后：

1. 读缓存的 cell × marker 平滑矩阵；
2. 每个 region 对每个 marker 算 `utag_mean__*` 和 `utag_std__*`；
3. 把这条固定向量交给 Lasso / 线性 SHAP。

这和 KRONOS/Eva 的 linear probe 同类：都是 embedding → 线性分类。差别是 UTAG 列仍带着 marker 名，理论上可以做回收。

导出：`scripts/export_cached_model_features.py --source utag-message-passing`  
解释：`verify_pseudo_label_explanations.py --mode precomputed --rules utag`

当前完整 fold 结果只在 **HNC-Wu2022**（`results/pseudo_label_explanations/utag_message_passing/`）。四数据集面板跑的是 UTAG native portraits，不是这条 linear embedding。

## 一次必须纠正的打分错误

原始 `utag_message_passing/summary.csv` 里，control AUC 已经有 0.93，但 `frac_lasso_passed = 0`。原因是那次运行用了默认 `--rules tabular`。表格规则在找 `composition::tumor`，不会认 `utag_mean__PanCK`。

用 `UTAG_RULES` 对同一份 `feature_ranks.csv` 重打分后（top-5）：

表：`tables/utag_linear_hnc_rescored_summary.csv`

| 任务 | 类型 | Lasso AUC | 回收率（Lasso / SHAP） | Top-1 | 重打分结论 |
|---|---|---:|---|---|---|
| tumor_high | control | 0.933 | 0.93 / 1.00 | utag_mean__CD45（top-5 含 PanCK） | **通过** |
| cd8_high | control | 0.928 | 1.00 / 1.00 | utag_std__CD8 | **通过** |
| cd8_clustering | spatial | 0.578 | 0.00 / 0.00 | utag_std__Vimentin / utag_mean__HLA-DR | 失败 |
| immune_exclusion | spatial | 0.815 | 0.00 / 0.00 | utag_std__Vimentin | 失败 |
| tumor_stroma_mixing | spatial | 0.832 | 0.00 / 0.00 | utag_std__Vimentin | 失败 |
| interface_immune | spatial | 0.816 | 0.00 / 0.00 | utag_std__Vimentin | 失败 |

Faithfulness 在 control 和多数空间任务上接近 1.0：模型确实在用这些 marker 统计。失败的是 **空间回收规则**。

## 为什么空间一定失败

`benchmark/motifs/recovery.py` 里，UTAG 空间规则的 hit 几乎全是：

```text
utag_domain::<cell type>__kk
```

miss（不许当 Top-1 的丰度泄漏）则是：

```text
utag_mean__PanCK / CD8 / CD68 / ...
```

Linear embedding **根本没有 domain 列**。它能贡献的只有 mean/std。于是：

- control：规则允许 `utag_mean__PanCK` / `utag_mean__CD8` → 能过
- spatial：规则要 domain，mean/std 还被当成泄漏 → 不能过

这不是拟合坏了，是 **表示和回收规则不匹配**。要把 UTAG 当空间解释器，必须走 domain portraits 或 fold 内 KMeans domain，而不是 message-passing 统计。

## 和 UTAG portraits 的分工

| | UTAG-Linear embedding | UTAG native portraits |
|---|---|---|
| 向量 | marker mean/std | domain 比例 + 画像（多数细胞类型 + centroid marker） |
| 能回答 | 哪些平滑 marker 与标签相关 | 哪个空间 domain 与标签相关 |
| HNC control | 通过 | tumor_high 失败（画像偏肿瘤但关联不够） |
| 空间 motif | 0/4 | 面板上 raw 1/8，协议 0/4 |

二者不要合成「UTAG 解释器」一个分数。文档和 `summarize_interpreter_panel.py` 也明确：不做「UTAG domains + Lasso/SHAP」这种捆绑。

## 结论

UTAG-Linear embedding 是 **可命名的 marker embedding + 线性模型**。在 HNC 上它能稳定找回肿瘤/CD8 相关 marker，也能用 faithfulness 证明这些列驱动预测。它 **不能** 验证空间 motif，除非改回收规则（例如允许 `utag_std__Vimentin` 代表混合），或改用 domain 表示。

四数据集上的 UTAG 主结果仍然是 native portraits，见 [06_corresponding_methods.md](06_corresponding_methods.md)。
