# 表格说明

全部由 `ExplainingConclusions/_compile_tables.py` 从
`TME_modeling_benchmark/results/` 汇总。数值以这些 CSV 为准，markdown 中的分数已四舍五入。

## 伪标签面板（10 个选定 motif）

| 文件 | 内容 |
|---|---|
| `selected_motif_catalog.csv` | 10 个 motif 的定义、样本量、组成/表达/密度 AUC |
| `global_method_comparison.csv` | 方法总分：control、空间 raw、空间 protocol |
| `selected_tasks_all_methods.csv` | 每个选定任务 × 方法的 pass/fail |
| `lasso_shap_fold_summary_all_tasks.csv` | 含匹配 control 在内的全部 tabular 任务 |
| `lasso_selected_task_verdicts.csv` | Lasso 选定任务 |
| `shap_selected_task_verdicts.csv` | SHAP 选定任务 |
| `lasso_shap_dominant_top1.csv` | 每个选定任务最常见的 Top-1 特征 |
| `lasso_shap_top1_features.csv` | 全部 Top-1 特征及 fold 频率 |
| `utag_linear_hnc_rescored_summary.csv` | UTAG message-passing + 正确 UTAG_RULES |
| `utag_linear_hnc_original_tabular_rules.csv` | 原始错误规则下的汇总（勿作主结论） |
| `utag_linear_hnc_rescored_folds.csv` | 重打分的 fold 明细 |
| `utag_native_portraits_all_tasks.csv` | 四数据集 UTAG 画像 |
| `utag_cellular_selected_task_verdicts.csv` | UTAG portraits 选定任务 |
| `kronos_embedding_probe.csv` | KRONOS Ridge AUC |
| `eva_embedding_probe.csv` | Eva Ridge AUC |
| `gnn_explainer_task_summary.csv` | 节点富集、region 通过率 |
| `gnn_explainer_selected_task_verdicts.csv` | GNN 选定任务门控 |
| `tabular_family_method_comparison.csv` | 四族单独 Lasso/SHAP 的方法总分 |
| `tabular_family_selected_task_verdicts.csv` | 四族选定任务 raw / protocol |
| `tabular_family_selected_tasks_wide.csv` | 选定任务 × 族 × 解释器 |
| `tabular_family_fold_summary.csv` | 四族全部任务（含匹配 control） |
| `tabular_family_dominant_top1.csv` | 四族每个选定任务最常见 Top-1 |
| `tabular_family_auc_compare.csv` | 单族 AUC / 回收率对照联合表 |
| `tabular_pair_method_comparison.csv` | 六对两两 Lasso/SHAP 方法总分 |
| `tabular_pair_selected_task_verdicts.csv` | 两两选定任务 raw / protocol |
| `tabular_pair_selected_tasks_wide.csv` | 选定任务 × 两两组合 × 解释器 |
| `tabular_pair_fold_summary.csv` | 两两全部任务（含匹配 control） |
| `tabular_pair_dominant_top1.csv` | 两两每个选定任务最常见 Top-1 |
| `tabular_pair_auc_compare.csv` | 两两 AUC / 回收率对照联合表 |
| `tabular_recovery_regex.csv` | 每条 motif 的 hit / miss_as_top 正则（一行一条） |
| `tabular_recovery_regex_by_task.csv` | 按任务汇总的回收正则与通过条件 |

## 临床终点（无伪标签 ground truth）

| 文件 | 内容 |
|---|---|
| `clinical_stability_lasso_top5.csv` | 组成 Stability Lasso 每任务 Top-5 |
| `clinical_shap_comp_expr_cox_top5.csv` | 组成+表达 Cox SHAP，方向一致的 Top-5 |
| `clinical_shap_cytocommunity_tcn_biological.csv` | 精选任务上带生物标签的 TCN Top-1 |
