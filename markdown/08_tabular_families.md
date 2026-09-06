# 单族表格：composition / density / mixing / point-pattern

联合 Lasso/SHAP（`tabular/`）把四族拼成一张表再做 L1 / Ridge。那次结果不能拆回「每一族自己当解释器会怎样」。这里是四次独立 `--feature-sources` 运行，目录：

- `results/pseudo_label_explanations_panel/tabular_composition/`
- `results/pseudo_label_explanations_panel/tabular_density/`
- `results/pseudo_label_explanations_panel/tabular_mixing/`
- `results/pseudo_label_explanations_panel/tabular_point_pattern/`

脚本、阈值、15 fold、`--rules tabular` 与联合跑相同。每个目录有自己的 `tabular_features.csv`，不会读到联合缓存。

表：`tables/tabular_family_method_comparison.csv`，`tables/tabular_family_selected_task_verdicts.csv`，`tables/tabular_family_dominant_top1.csv`，`tables/tabular_family_auc_compare.csv`。

Lasso 与 SHAP 在这四族上的任务级 pass/fail **完全一致**。下面分数用 Lasso；SHAP 只在 Top-1 列名不同时另写。

## 方法总分

和联合表一样：选定 control 是 Jackson `tumor_high` + TNBC `cd8_high`。空间 protocol 仍看同数据集匹配 control。

| 族 | 选定 control | 空间 raw | 空间 protocol | 一句话 |
|---|---|---|---|---|
| 联合四族（对照） | 1/2 | 6/8 | 5/6 | mixing 抢走 TNBC CD8 名字；空间靠 mixing |
| composition | **2/2** | 0/8 | 0/6 | control 第一次齐；空间全是组成，回收为 0 |
| density | **2/2** | 2/8 | 2/8 | 四个数据集匹配 control 都过；两条空间通过要打折 |
| mixing | 0/2 | 6/8 | 0/8 | 能点名邻域；丰度 control 全灭，空间被挡住 |
| point-pattern | 0/2 | 6/8 | 0/8 | 能点名 K/L；control 全灭；exclusion 擦边过 |

没有一族同时做到：选定 control 过关 **并且** 空间 protocol 像联合表那样拿到 5/6。

## 选定任务：回收

| 数据集 | 任务 | 类型 | composition | density | mixing | point-pattern | 联合 |
|---|---|---|---|---|---|---|---|
| Jackson | tumor_high | control | pass | pass | fail | fail | pass |
| TNBC | cd8_high | control | **pass** | **pass** | fail | fail | **fail** |
| Jackson | t_tumor_mixing | spatial | fail | fail | pass* | pass* | pass |
| Jackson | cd8_tumor_contact | spatial | fail | pass | pass* | pass* | pass |
| Jackson | macrophage_tumor_niche | spatial | fail | fail | pass* | pass* | pass |
| Jackson | apc_t_contact | spatial | fail | fail | pass* | pass* | pass |
| HNC | cd8_clustering | spatial | fail | fail | pass* | pass* | pass |
| HNC | immune_exclusion | spatial | fail | pass | fail* | pass* | fail |
| METABRIC | tumor_stroma_mixing | spatial | fail† | fail | pass* | fail* | pass† |
| METABRIC | interface_immune | spatial | fail† | fail | fail* | fail* | fail† |

`*`：该族匹配 control 失败，协议记 `fail_control_block`。  
`†`：METABRIC 本数据集 `tumor_high` 或 `cd8_high` 失败被挡住。density 是唯一四个数据集匹配 control 都过的族，所以它的空间分数没有 `fail_control_block`。

## 各族实际选中了什么

### composition

| 任务 | Lasso Top-1 |
|---|---|
| Jackson tumor_high | `composition::Stromal cells`（1.00；top-5 仍有 Tumor 亚型，故 control 过） |
| TNBC cd8_high | `composition::CD8^+PD1^+T_{Ex}`（0.73；SHAP 更常是 `composition::CD8^+T`） |
| Jackson 四条空间 | Macrophage / Tumor 亚型比例 |
| HNC 两条空间 | `composition::CD8 T cell` |
| METABRIC 两条空间 | `HR- CK7+` / Myoepithelial |

空间任务的回收规则要 mixing 或 type-specific K/L，并且组成列常被标成丰度泄漏。只给比例时，模型只能点组成，于是空间 raw **0/8**。

TNBC `cd8_high` 在联合表里失败（Top-1 是 mixing），在这里 15/15 fold 都回到 CD8 组成。联合 L1 不是「不会解释丰度」，而是 **特征库里有 mixing 时不必回到 `composition::CD8`**。

Jackson / HNC 的 `tumor_high` 能过，是因为类型名里有 `Tumor`。TNBC 和 METABRIC 的 `tumor_high` 回收是 0：Top-1 是 `composition::Fibroblasts`，亚型叫 HR+ / HER2+ 等，对不上 `composition::.*tumor`。AUC 仍然 ≈ 0.996。这是命名规则和细胞类型词表的错位，不是模型没把肿瘤和间质分开。

### density

| 任务 | Lasso Top-1 |
|---|---|
| Jackson tumor_high | `density::tumor_area_ratio`（1.00） |
| TNBC cd8_high | `density::tissue_density::CD8^+PD1^+T_{Ex}`（0.67） |
| cd8_tumor_contact | `density::tumor_density::T cell`（1.00） |
| t_tumor_mixing | 同上 `tumor_density::T cell`（1.00），但回收失败 |
| immune_exclusion | `density::tumor_area_ratio`（1.00） |
| cd8_clustering | `density::tissue_density::Tumor (Proliferating)`（0.73） |

两条空间「通过」都要按规则读，不能写成 density 已经提取到空间机制。

**`cd8_tumor_contact`：** 规则显式允许 `tumor_density::(cd8|t cells?)`。同一列 `tumor_density::T cell` 在 `t_tumor_mixing` 上 **不过**，因为那条规则的 hit 只有 mixing / Ripley。这是规则不对称，不是接触比混合更「像密度」。

**`immune_exclusion`：** 15/15 fold 的 Top-1 仍是 `tumor_area_ratio`，和联合表失败时一样。它能过，是因为 rank 2 稳定出现 `tissue_density::CD8 T cell`，而 `miss_as_top` 没有把 `tumor_area_ratio` 算进丰度泄漏。联合表里 mixing 把 CD8 密度挤出 top-5，所以同样的 Top-1 在联合跑里失败。

四个数据集的匹配 control 都过：`tumor_area_ratio` / `tissue_density::Tumor` 对得上 `tumor_high`，`tissue_density` 里的 T/CD8 对得上 `cd8_high`。所以 density 的空间 2/8 没有被 `fail_control_block` 挡掉。

### mixing

| 任务 | Lasso Top-1 |
|---|---|
| Jackson 四条空间 | 与联合表同类：T→T、Mac→Stroma、Mac→T |
| cd8_clustering | `CD8 T cell → CD8 T cell`（1.00） |
| tumor_stroma_mixing | `shannon_entropy_normalized`（0.67） |
| immune_exclusion | `CD8 T cell → Tumor`（0.53），但规则要的是 `tumor_density::cd8` / CD8 K/L，故 raw 失败 |
| TNBC cd8_high | `shannon_entropy_normalized`（0.47） |
| Jackson tumor_high | `Stromal → Stromal`（0.73） |

空间 raw **6/8**（四条 Jackson + clustering + METABRIC mixing）。`interface_immune` 和 `immune_exclusion` 不过。协议 **0/8**：没有组成/密度列，两条选定 control 和全部匹配 control 都是 0。

空间 AUC 已经接近联合表（clustering 0.940 vs 0.948；四条 Jackson 接触/混合只低 0.02–0.03）。联合表在空间上多出来的，主要是 exclusion 上的 `tumor_area_ratio`（联合 0.979，mixing 单独 0.929）。

### point-pattern

| 任务 | Lasso Top-1 |
|---|---|
| cd8_clustering | `point-pattern::CD8 T cell_L_r200`（1.00） |
| immune_exclusion | 分散：B cell / Tumor L，只有部分 fold 的 top-5 有 CD8 K/L（lasso 0.53，shap 0.73） |
| macrophage_tumor_niche | 常是 Tumor 的 L（规则 `tumor.*[_ ][kl]_r` 能过），AUC 只有 0.57 |
| TNBC cd8_high | `CD8^+PD1^+T_{Ex}_L_r10`（0.47），对不上组成规则 |

空间 raw 也是 **6/8**，但成员和 mixing 不同：多了擦边的 `immune_exclusion`，少了 `tumor_stroma_mixing`。`macrophage_tumor_niche` 回收 15/15，但 faithfulness 只有 0.47、AUC ≈ 0.57，是规则碰巧命中肿瘤 Ripley，不是稳定的巨噬细胞空间解释。control 全灭，协议 0/8。

## AUC：单族相对联合表

选定任务的 Lasso 均值 AUC（`tables/tabular_family_auc_compare.csv`）：

| 任务 | 联合 | composition | density | mixing | point-pattern |
|---|---:|---:|---:|---:|---:|
| Jackson tumor_high | 0.973 | 0.999 | 0.994 | 0.955 | 0.837 |
| TNBC cd8_high | 0.975 | 0.995 | 0.986 | 0.972 | 0.924 |
| t_tumor_mixing | 0.883 | 0.690 | 0.793 | 0.852 | 0.717 |
| cd8_tumor_contact | 0.881 | 0.676 | 0.792 | 0.855 | 0.692 |
| macrophage_tumor_niche | 0.790 | 0.529 | 0.667 | 0.802 | 0.566 |
| apc_t_contact | 0.684 | 0.572 | 0.627 | 0.678 | 0.592 |
| cd8_clustering | 0.948 | 0.604 | 0.681 | 0.940 | 0.914 |
| immune_exclusion | 0.979 | 0.638 | 0.950 | 0.929 | 0.858 |
| tumor_stroma_mixing | 0.874 | 0.594 | 0.821 | 0.820 | 0.747 |
| interface_immune | 0.585 | 0.580 | 0.631 | 0.517 | 0.623 |

组成几乎解释不了空间标签（除了本来就接近随机的 `interface_immune`）。mixing 单独就能接近联合空间 AUC。exclusion 的高 AUC 主要来自密度/面积，不是邻域分数。

## 怎么读这次拆开

1. **联合表的 TNBC `cd8_high` 失败，不是 CD8 丰度不可解释。** 只给 composition 或 density 时，两条选定 control 都是 2/2。失败发生在「同一张表里还有 mixing」的时候。
2. **联合表的空间 5/6 来自 mixing 列。** composition 空间 0/8。point-pattern 能过若干回收规则，但 AUC 和 faithfulness 明显弱于 mixing，而且 control 过不了，协议分为 0。
3. **density 的 2 条空间通过不能当成空间解释已经成立。** 一条是规则允许 `tumor_density::T cell`，同一特征在 mixing 任务上失败；一条是 Top-1 仍为 `tumor_area_ratio`，只是 rank 2 的 CD8 组织密度被算进 hit。
4. **没有单族 Lasso/SHAP 通过完整协议。** 丰度族过 control、不过空间；空间族过 raw、被 control 门挡住。联合表把两半拼在一起：control 变成 1/2，空间在匹配 control 过关的数据集上变成 5/6。

这和 [02_lasso.md](02_lasso.md) 的窄结论一致，只是把「特征库里一旦有 mixing，Lasso 不必回到 composition::CD8」变成了直接对照。

两两组合见 [09_tabular_pairs.md](09_tabular_pairs.md)：联合四族 ≈ composition + mixing；composition + point-pattern 是首个不含 mixing 却同时拿到 control 2/2 与空间协议分的组合。
