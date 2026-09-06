# Global 解释协议

## Global 指什么

Global 解释器在 **region 级** 工作：

1. 每个 region 得到一条固定长度特征，或每个细胞内的节点分数；
2. 下游是线性分类器、线性 SHAP，或对节点分数做富集；
3. 不把 region 切成 MIL window。

这和 Gated Attention MIL 不同。MIL 解释的是 bag 里的 window，不是这张 region×feature 表。

## 命名表格特征（Lasso / SHAP 共用）

`scripts/verify_pseudo_label_explanations.py --mode tabular`，特征源：

- `composition`：细胞类型比例
- `density`：组织/肿瘤密度、肿瘤面积比
- `mixing`：邻域分数、局部混合、熵
- `point-pattern`：按细胞类型分层的 Ripley K/L（`--by-type`）

Lasso 和 SHAP **共用特征表，但独立门控**：

- Lasso：`LinearClassifier(C=0.5, l1_ratio=1.0)`，按 `|β|` 排序
- SHAP：先拟合 Ridge（`C=1.0, l1_ratio=0.0`），再算 `φ_j = β_j (z_j)` 的 mean |φ| 排序
- Faithfulness：打乱 SHAP top-k 特征后的 AUC 下降，应大于打乱随机 k 个特征

默认 3 个 seed × 5 个 patient fold = 15 fold。回收看 top-5。空间任务若 Top-1 是已知丰度泄漏（例如 `composition::CD8`），即使后面有空间特征也判失败。

四族也曾 **单独** 各跑一次（各自的 `--output-dir`，避免读到联合 `tabular/` 缓存）。那是把每一族当作自己的解释器，不能用联合表的 Top-1 反推。结果见 [08_tabular_families.md](08_tabular_families.md)。expression 仍不在 `--mode tabular` 的 featurizer 里。

## 回收规则在问什么

Control（`tumor_high` / `cd8_high`）要求 top-5 里出现对应组成或密度列。

空间任务要求 top-5 里出现 mixing / type-specific K·L / tumor-compartment 密度，并且 Top-1 不能是单纯丰度。

这就是 HNC `immune_exclusion` 失败的原因：模型极强（AUC 0.98），但解释指向 `tumor_area_ratio`，没有指向 CD8 的空间特征。

## 两套 control

面板选定的 control 只有两条：Jackson `tumor_high`、TNBC `cd8_high`。方法级 `control_ok` 要求这两条都过。Lasso/SHAP 是 1/2。

空间任务另外用 **同数据集匹配 control**：每个数据集都附带自己的 `tumor_high` 和 `cd8_high`。`summarize_interpreter_panel.py` 检查这两条是否都过；不过则该数据集所有空间任务记 `fail_control_block`。

这会放大 control 失败：

- METABRIC 的 `cd8_high` 回收失败 → `tumor_stroma_mixing` 即使 Lasso 15/15 fold 通过，协议上仍被挡住
- Embedding 的 control 阈值是 AUC 0.90，绝大多数数据集过不了，于是 Eva 的 8 个空间 AUC 全部不算协议分

读结论时请同时看 **raw** 和 **protocol**。

## 标签

面板使用 `results/pseudo_labels/<dataset>_v2.csv`，`--label-version v2`。  
GNN 的 SPACE-GM 训练脚本写的是 noisy v2 标签；验证时仍用同一套 v2 伪标签协议，以便和 Lasso/SHAP 对齐。
