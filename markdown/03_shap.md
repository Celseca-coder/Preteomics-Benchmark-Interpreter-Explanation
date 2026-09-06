# Global 线性 SHAP

## 方法

这里的 SHAP 不是 TreeSHAP，也不是 DeepSHAP。它是 **Ridge 线性模型的精确 SHAP**：

\[
\varphi_j = \beta_j\, z_j
\]

\(z_j\) 是 held-out region 上该特征的标准化值。按 mean \(|\varphi_j|\) 排序，与 Lasso 的 `|β|` 排序独立门控。

Faithfulness 也挂在 Ridge 上：打乱 SHAP top-5 后的 AUC 下降，应大于打乱随机 5 个特征。

脚本与 Lasso 相同（`--mode tabular`），只是 `shap_passed` 单独计数。同一跑里可以换 `--feature-sources`：四族联合、单族、两两组合各有目录；SHAP 与 Lasso 共用特征表。

## 选定任务结果（四族联合）

特征：`composition + density + mixing + point-pattern`。  
目录：`results/pseudo_label_explanations_panel/tabular/`。  
表：`tables/shap_selected_task_verdicts.csv`。

与联合 Lasso **任务级结论完全一致**（同一 10 个选定 motif 的 pass/fail）。

| 数据集 | 任务 | Ridge AUC | SHAP 回收率 | Faithfulness | raw | 协议 |
|---|---|---:|---:|---:|---|---|
| Jackson | tumor_high | 0.972 | 1.00 | 1.00 | pass | pass_control |
| Jackson | t_tumor_mixing | 0.880 | 1.00 | 1.00 | pass | pass_spatial |
| Jackson | cd8_tumor_contact | 0.876 | 1.00 | 1.00 | pass | pass_spatial |
| Jackson | macrophage_tumor_niche | 0.791 | 1.00 | 1.00 | pass | pass_spatial |
| Jackson | apc_t_contact | 0.681 | 1.00 | 0.93 | pass | pass_spatial |
| HNC | cd8_clustering | 0.951 | 1.00 | 1.00 | pass | pass_spatial |
| HNC | immune_exclusion | 0.980 | 0.00 | 0.87 | fail | fail_spatial |
| METABRIC | tumor_stroma_mixing | 0.872 | 1.00 | 1.00 | pass | fail_control_block |
| METABRIC | interface_immune | 0.587 | 0.00 | 0.53 | fail | fail_control_block |
| TNBC | cd8_high | 0.977 | 0.20 | 0.93 | fail | fail_control |

总分：选定 control **1/2**，空间 raw **6/8**，空间 protocol **5/6**。

Faithfulness 几乎总是过：打乱 top 特征比打乱随机特征掉分更多。所以 SHAP 的失败来自 **回收规则**（点错生物学名字），不是「分数与模型无关」。

## Top-1 与联合 Lasso 的差别

表：`tables/lasso_shap_dominant_top1.csv`

多数任务上 SHAP 和 Lasso 指向同一类 mixing 特征，但具体列不完全相同：

| 任务 | Lasso Top-1 | SHAP Top-1 |
|---|---|---|
| tumor_high | composition::Stromal cells | density::tumor_area_ratio |
| apc_t_contact | Macrophage→T cell（1.00） | 同左（0.80） |
| macrophage_tumor_niche | Macrophage→Stromal（1.00） | 同左（0.80） |
| cd8_clustering | CD8→CD8（1.00） | 同左（1.00） |
| immune_exclusion | tumor_area_ratio | tumor_area_ratio |
| tumor_stroma_mixing | shannon_entropy_normalized | 同左 |
| cd8_high（TNBC） | Treg→耗竭 CD8 | M2 Mac→CD8^+T |

SHAP 把 `tumor_high` 更多地归到肿瘤面积/密度，Lasso 更多地罚间质比例。两者都能在 top-5 里碰到肿瘤相关列，所以 control 都过。

`apc_t_contact` 的 1/15 fold faithfulness 失败（打乱 top 后 AUC 几乎不降），是唯一反复出现的 faithfulness 缺口，与该任务本身较弱（AUC 0.68）一致。

## 单族 SHAP

与单族 Lasso 同一次跑（`--feature-sources` 只开一族）。  
表：`tables/tabular_family_method_comparison.csv`（`explainer=shap`），`tables/tabular_family_dominant_top1.csv`。

| 族 | 选定 control | 空间 raw | 空间 protocol | SHAP 常见 Top-1 |
|---|---|---|---|---|
| composition | 2/2 | 0/8 | 0/6 | TNBC `cd8_high` → `composition::CD8^+T`；空间任务全是组成列 |
| density | 2/2 | 2/8 | 2/8 | control → `tumor_area_ratio` / CD8 `tissue_density`；exclusion Top-1 仍是面积比 |
| mixing | 0/2 | 6/8 | 0/8 | 空间 → 邻域分数；TNBC `cd8_high` → 熵，过不了组成规则 |
| point-pattern | 0/2 | 6/8 | 0/8 | clustering → `CD8 T cell_L_r200`；control 点不到组成名 |

单族上 SHAP 与对应 Lasso **任务级完全一致**。更细的生物学读法见 [08_tabular_families.md](08_tabular_families.md)（文中默认写 Lasso，SHAP 同行）。

## 两两组合 SHAP

表：`tables/tabular_pair_method_comparison.csv`（`explainer=shap`），`tables/tabular_pair_dominant_top1.csv`，`tables/tabular_pair_selected_tasks_wide.csv`。

| 组合 | 选定 control | 空间 raw | 空间 protocol | 相对 Lasso |
|---|---|---|---|---|
| composition + density | 2/2 | 2/8 | 2/8 | 一致 |
| composition + mixing | **0/2** | 6/8 | **1/2** | **分叉**（见下） |
| density + mixing | 1/2 | 6/8 | 0/8 | 一致 |
| composition + point-pattern | 2/2 | 5/8 | 5/6 | 一致 |
| density + point-pattern | 2/2 | 6/8 | 6/8 | 一致 |
| mixing + point-pattern | 0/2 | 6/8 | 0/8 | 一致 |

**唯一分叉：`composition + mixing`。** Jackson `tumor_high` 的 SHAP 回收率 0.47（不过 0.5），Lasso 是 0.67。SHAP Top-1 仍常是 `composition::Stromal cells`，但 top-5 命中肿瘤列的 fold 不够一半。后果：

- 选定 control：SHAP **0/2**（Lasso 1/2）
- Jackson 匹配 control 因这条失败 → 四条 Jackson 空间 raw 虽过，协议全是 `fail_control_block`
- 空间 protocol 分母只剩 HNC 两条，过了 clustering → **1/2**（Lasso 是 5/6）

其余五对与 Lasso 同判。协议数字最好的仍是 **density + point-pattern（2/2 · 6/8）**，但 exclusion 的 SHAP Top-1 同样是 `tumor_area_ratio`（规则漏洞）。不含 mixing 且 control 齐、空间有协议分的是 **composition + point-pattern（2/2 · 5/6）**。明细见 [09_tabular_pairs.md](09_tabular_pairs.md)。

## 和临床 SHAP 的关系

伪标签 SHAP 验证的是：**已知生成过程时，线性归因会不会指出对的命名特征。**

临床 SHAP（composition/expression Cox，或 CytoCommunity TCN）没有这种 ground truth，只能报稳定的高 |SHAP| 特征。那部分见 [07_clinical_global.md](07_clinical_global.md)。不要把临床 Top-1 TCN 直接当成这里已经验证过的空间解释。

## 结论

- **四族联合**：选定 control 1/2，不能写成已验证；任务级与联合 Lasso 相同。
- **单族**：与单族 Lasso 完全同判；丰度族过 control、空间族过 raw 被 control 门挡住。
- **两两**：五对与 Lasso 一致；仅 `composition + mixing` 因 Jackson `tumor_high` 让 SHAP 更严（0/2 control、1/2 空间 protocol）。

线性 SHAP 的价值仍是：

1. 确认 Ridge（无稀疏）和 Lasso（稀疏）在绝大多数设定下指向同一套机制；
2. 用 perturbation drop 确认这些特征确实驱动预测；
3. 在 `composition + mixing` 上暴露出与 Lasso 不同的 control 边界（Ridge 归因更散，tumor 回收擦边不过）。

若只选一个 Global 表格解释器做主结果，联合表上 Lasso 和 SHAP 可以并列；报两两组合时要把 `composition + mixing` 的分叉写清楚，不要写成「SHAP 全程等于 Lasso」。
