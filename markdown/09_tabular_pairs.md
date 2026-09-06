# 两两组合表格：六对 Lasso / SHAP

单族跑证明：丰度族过 control、不过空间；空间族过 raw、被 control 挡住。联合四族是 1/2 control + 5/6 空间。这里把四族两两拼成六对，各自独立 `--feature-sources`，目录：

- `tabular_composition_density/`
- `tabular_composition_mixing/`
- `tabular_density_mixing/`
- `tabular_composition_point_pattern/`
- `tabular_density_point_pattern/`
- `tabular_mixing_point_pattern/`

脚本：`scripts/run_tabular_family_combos.sh`（`MODE=pairs`）。阈值、15 fold、`--rules tabular` 与单族 / 联合表相同。

表：`tables/tabular_pair_method_comparison.csv`，`tables/tabular_pair_selected_task_verdicts.csv`，`tables/tabular_pair_dominant_top1.csv`，`tables/tabular_pair_auc_compare.csv`。

下面分数默认用 **Lasso**。`composition + mixing` 上 SHAP 与 Lasso 第一次在选定 control 上分叉，文末单独写。

## 方法总分

| 组合 | 选定 control | 空间 raw | 空间 protocol | 一句话 |
|---|---|---|---|---|
| 联合四族（对照） | 1/2 | 6/8 | 5/6 | mixing 抢走 TNBC CD8 |
| composition+density | **2/2** | 2/8 | 2/8 | 几乎等于 density 单族 |
| composition+mixing | 1/2 | 6/8 | 5/6 | **复现联合表** |
| density+mixing | 1/2 | 6/8 | **0/8** | tumor 过、CD8 不过 → 全数据集挡住 |
| composition+point-pattern | **2/2** | 5/8 | **5/6** | 首个「control 齐 + 空间有协议分」且不含 mixing |
| density+point-pattern | **2/2** | 6/8 | **6/8** | 协议分最高；exclusion 仍是面积比漏洞 |
| mixing+point-pattern | 0/2 | 6/8 | 0/8 | 等于 mixing 单族：空间有、control 无 |

## 选定任务：Lasso raw

| 数据集 | 任务 | comp+dens | comp+mix | dens+mix | comp+pp | dens+pp | mix+pp |
|---|---|---|---|---|---|---|---|
| Jackson | tumor_high | pass | pass | pass | pass | pass | fail |
| TNBC | cd8_high | **pass** | **fail** | **fail** | **pass** | **pass** | fail |
| Jackson | t_tumor_mixing | fail | pass | pass* | pass | pass | pass* |
| Jackson | cd8_tumor_contact | pass | pass | pass* | pass | pass | pass* |
| Jackson | macrophage_tumor_niche | fail | pass | pass* | pass | pass | pass* |
| Jackson | apc_t_contact | fail | pass | pass* | pass | pass | pass* |
| HNC | cd8_clustering | fail | pass | pass* | pass | pass | pass* |
| HNC | immune_exclusion | pass | fail | fail* | fail | pass | fail* |
| METABRIC | tumor_stroma_mixing | fail | pass† | pass* | fail† | fail | pass* |
| METABRIC | interface_immune | fail | fail† | fail* | fail† | fail | fail* |

`*`：匹配 control 失败 → protocol `fail_control_block`。  
`†`：METABRIC 被挡住（本数据集 control 未齐）。

## 各对在干什么

### composition + density → 还是丰度解释器

选定 control **2/2**。空间 raw / protocol 都是 **2/8**，和 density 单族同一对任务：`cd8_tumor_contact`、`immune_exclusion`。Top-1 几乎全是 density；TNBC `cd8_high` 回到 `composition::CD8^+PD1^+T_{Ex}`。加组成没有打开新的空间通道。

### composition + mixing → 联合表的最小充分集

选定 control **1/2**，空间 raw **6/8**，protocol **5/6**——和联合四族 **同一套数字**。TNBC `cd8_high` Top-1 仍是 mixing（熵 / 邻域），不是组成。Jackson 四条空间和 HNC clustering 的 Top-1 都是 mixing。

结论：联合表里 density 和 point-pattern **不是** 造成 1/2 + 5/6 的必要零件；**composition + mixing 已经够**。四族联合只是把这两族的行为包在更大的特征库里。

注意：这条上 **Lasso 与 SHAP 在 Jackson `tumor_high` 分叉**（Lasso 0.67 过、SHAP 0.47 不过），所以 SHAP 的选定 control 是 0/2，空间 protocol 也变成 1/2（Jackson 匹配 control 因 SHAP 的 tumor_high 失败而被挡住）。这是此前单族和联合表没有出现的差异。

### density + mixing → tumor 过、CD8 灭、空间全挡

Jackson / HNC / TNBC / METABRIC 的 `tumor_high` 都回收通过（Top-1 = `tumor_area_ratio`），但四个数据集的 `cd8_high` 都失败（Top-1 进了 mixing）。于是 **每个数据集的匹配 control 都不齐**，空间 raw 虽有 6/8，protocol **0/8**。

这比 composition+mixing 更糟：composition 至少能救 Jackson / HNC 的匹配 `cd8_high`；density 在和 mixing 同桌时 **救不了 CD8 名字**。

### composition + point-pattern → 首个「双过」且不含 mixing

选定 control **2/2**（TNBC 回到 CD8 组成）。空间 raw **5/8**，protocol **5/6**（METABRIC 仍挡）。Jackson 四条空间和 HNC clustering 的 Top-1 主要是 Ripley（如 `T cell_L_r10`、`CD8 T cell_L_r200`）。`immune_exclusion` 失败：Top-1 掉回 `composition::CD8 T cell`（丰度泄漏）。

这是第一次在 **没有 mixing** 的情况下同时做到：选定 control 齐，且空间 protocol 不是 0。AUC 仍明显低于 composition+mixing（Jackson 接触/混合约 0.74 vs 0.86）。

### density + point-pattern → 协议分最高，但 exclusion 仍是漏洞

选定 control **2/2**，空间 raw / protocol **6/8**。四个数据集匹配 control 都过，所以没有 `fail_control_block`。METABRIC 两条空间仍 raw 失败。

`immune_exclusion` 再次「通过」：15/15 fold Top-1 仍是 `density::tumor_area_ratio`，hit 来自 rank 2 的 `tissue_density::CD8 T cell`——和 density 单族同一规则漏洞。`t_tumor_mixing` 的 Top-1 是 `tumor_density::T cell`（规则对 mixing 任务不允许这列单族通过，但这里和 point-pattern 同跑时 raw 过了，多半是 top-5 里进了 Ripley）。读这条组合时，**不要把 6/8 写成干净的空间机制提取**。

### mixing + point-pattern → 纯空间特征库

选定 control **0/2**，空间 raw **6/8**，protocol **0/8**。Top-1 几乎全是 mixing；point-pattern 很少抢到第一名。行为接近 mixing 单族。

## AUC 对照（选定任务，Lasso 均值）

| 任务 | 联合 | c+d | c+m | d+m | c+pp | d+pp | m+pp |
|---|---:|---:|---:|---:|---:|---:|---:|
| Jackson tumor_high | 0.973 | 0.999 | 0.968 | 0.967 | 0.989 | 0.984 | 0.958 |
| TNBC cd8_high | 0.975 | 0.991 | 0.974 | 0.974 | 0.967 | 0.958 | 0.974 |
| t_tumor_mixing | 0.883 | 0.778 | 0.858 | 0.865 | 0.752 | 0.799 | 0.868 |
| cd8_tumor_contact | 0.881 | 0.774 | 0.860 | 0.868 | 0.737 | 0.799 | 0.869 |
| macrophage_tumor_niche | 0.790 | 0.652 | 0.808 | 0.808 | 0.568 | 0.634 | 0.779 |
| apc_t_contact | 0.684 | 0.656 | 0.691 | 0.716 | 0.610 | 0.632 | 0.648 |
| cd8_clustering | 0.948 | 0.689 | 0.939 | 0.958 | 0.906 | 0.912 | 0.937 |
| immune_exclusion | 0.979 | 0.953 | 0.929 | 0.978 | 0.890 | 0.925 | 0.962 |
| tumor_stroma_mixing | 0.874 | 0.820 | 0.828 | 0.833 | 0.757 | 0.762 | 0.864 |
| interface_immune | 0.585 | 0.667 | 0.524 | 0.513 | 0.627 | 0.630 | 0.586 |

含 mixing 的组合在空间 AUC 上最接近联合表。含 point-pattern、不含 mixing 的组合 AUC 更低，但 control 名字更干净。

## 怎么读这六对

1. **联合四族 ≈ composition + mixing。** 同样的 1/2 control 和 5/6 空间；TNBC CD8 被 mixing 抢走。density / Ripley 不是联合表协议分的必要来源。
2. **mixing 一旦进表，CD8 丰度 control 就危险。** composition+mixing、density+mixing、mixing+point-pattern、联合表都是如此。只有「丰度族 ± point-pattern、没有 mixing」时，选定 control 才能稳定 2/2。
3. **composition + point-pattern 是最干净的「双过」候选。** control 2/2，空间 protocol 5/6，Top-1 真是 K/L，不是 mixing。代价是 AUC 和联合表有差距，且 METABRIC 仍挡。
4. **density + point-pattern 的 6/8 协议分不能照单全收。** exclusion 的 Top-1 仍是肿瘤面积比；部分 Jackson 通过也掺了 `tumor_density::`。
5. **没有一对既（a）选定 control 2/2、（b）空间 Top-1 稳定是真正的邻域机制、（c）exclusion 也不靠面积比漏洞。** composition+mixing 有（b）没（a）；composition+point-pattern 有（a）和较弱的（b）；density+point-pattern 数字最好但（b）（c）有水分。

单族结论见 [08_tabular_families.md](08_tabular_families.md)。
