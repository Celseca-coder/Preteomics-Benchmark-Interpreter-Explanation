# Global Lasso

## 方法

Lasso 在这里不是事后解释器，而是 **带 L1 的 Logistic Regression**：在每个 patient fold 的训练集上拟合，在 held-out region 上预测，并用 `|β|` 排序做特征回收。

- 脚本：`scripts/verify_pseudo_label_explanations.py --mode tabular`
- 模型：`LinearClassifier(C=0.5, l1_ratio=1.0)`
- 特征：composition + density + mixing + type-specific point pattern
- 通过：该 fold 的 top-5 命中预期模式；空间任务还要求 Top-1 不是丰度泄漏
- 任务级通过：≥50% fold 通过
- top-k 细指标（hit@1 / P@3 / MRR / `|H|` / recall@k）见 [12_tabular_topk_precision.md](12_tabular_topk_precision.md)，不改 pass/fail

临床终点上的 Stability Lasso（真实 OS/response 等）见 [07_clinical_global.md](07_clinical_global.md)。本文只讨论 motif 伪标签验证。

## 选定任务结果

表：`tables/lasso_selected_task_verdicts.csv`，`tables/lasso_shap_fold_summary_all_tasks.csv`

| 数据集 | 任务 | 类型 | 均值 AUC | fold 回收率 | raw | 协议 |
|---|---|---|---:|---:|---|---|
| Jackson | tumor_high | control | 0.973 | 1.00 | pass | pass_control |
| Jackson | t_tumor_mixing | spatial | 0.883 | 1.00 | pass | pass_spatial |
| Jackson | cd8_tumor_contact | spatial | 0.881 | 1.00 | pass | pass_spatial |
| Jackson | macrophage_tumor_niche | spatial | 0.790 | 1.00 | pass | pass_spatial |
| Jackson | apc_t_contact | spatial | 0.684 | 1.00 | pass | pass_spatial |
| HNC | cd8_clustering | spatial | 0.948 | 1.00 | pass | pass_spatial |
| HNC | immune_exclusion | spatial | 0.979 | 0.00 | fail | fail_spatial |
| METABRIC | tumor_stroma_mixing | spatial | 0.874 | 1.00 | pass | fail_control_block |
| METABRIC | interface_immune | spatial | 0.585 | 0.00 | fail | fail_control_block |
| TNBC | cd8_high | control | 0.975 | 0.00 | fail | fail_control |

总分：选定 control **1/2**（方法级未过关），空间 raw **6/8**，空间协议 **5/6**。  
5/6 不是因为选定 control 过了，而是 Jackson / HNC **本数据集**的 `tumor_high`+`cd8_high` 都过了（Jackson 的 `cd8_high` 是匹配 control，不是面板选定的那条）。METABRIC 两条空间因本数据集 `cd8_high` 失败被挡住。

## 它实际选中了什么

表：`tables/lasso_shap_dominant_top1.csv`

| 任务 | 最常见 Top-1（15 fold 中的比例） |
|---|---|
| tumor_high | `composition::Stromal cells`（0.73；top-5 仍命中肿瘤组成/面积，故 control 通过） |
| t_tumor_mixing | `mixing::neighbor_fraction__T cell__to__T cell`（0.40） |
| cd8_tumor_contact | `mixing::neighbor_fraction__T cell__to__T cell`（0.47） |
| macrophage_tumor_niche | `mixing::neighbor_fraction__Macrophage__to__Stromal cells`（1.00） |
| apc_t_contact | `mixing::neighbor_fraction__Macrophage__to__T cell`（1.00） |
| cd8_clustering | `mixing::neighbor_fraction__CD8 T cell__to__CD8 T cell`（1.00） |
| immune_exclusion | `density::tumor_area_ratio`（0.73）← 丰度泄漏 |
| tumor_stroma_mixing | `mixing::shannon_entropy_normalized`（1.00） |
| interface_immune | `mixing::neighbor_fraction__HER2+__to__HRlow CKlow`（0.47），且 AUC 接近随机 |
| cd8_high（TNBC） | mixing（Treg→耗竭 CD8 等），不是 `composition::CD8^+T` |

## 怎么读这些失败

**TNBC `cd8_high`：** 标签是 CD8 丰度，但 TNBC 把 CD8 拆成多个功能亚型。Lasso 更愿意用「谁挨着 CD8」而不是总 CD8 比例。预测很准，对解释协议来说算点错名。

**HNC `immune_exclusion`：** 标签在残差化 CD8% 和肿瘤% 之后，看 CD8 是否在肿瘤多边形外。线性模型用肿瘤面积比就能把标签分开，于是解释停在丰度，没有走到 CD8 空间特征。

**METABRIC `interface_immune`：** 组成和表达本来就接近随机（目录里 composition AUC 0.58）。Lasso 同样弱（0.585），回收为 0。这是信号不足，不是解释器单独坏了。

**METABRIC `tumor_stroma_mixing`：** 解释是对的（熵/混合），但被数据集级 `cd8_high` 门控挡住。若只问「Lasso 会不会指出混合」，答案是会。

## 结论

不能说 Lasso 已经通过解释协议。选定 control 是 1/2，`control_ok=False`。

能说的更窄：在同数据集丰度 control 回收成功的数据集上（Jackson、HNC），它会把若干空间 motif 点成 mixing / Ripley，而不是细胞比例。这和「面板指定的 CD8 丰度 control 已通过」不是一回事。TNBC `cd8_high` 失败说明：特征库里一旦有 mixing，Lasso 不必、也没有回到 `composition::CD8`。

四族拆开后这句话可以直接对照：只给 composition 或 density 时，选定 control 是 2/2；只给 mixing 时，空间 raw 6/8、control 0/2。见 [08_tabular_families.md](08_tabular_families.md)。
