# 各 motif「通过」的回收正则（tabular Lasso / SHAP）

源码：`TME_modeling_benchmark/benchmark/motifs/recovery.py` 里的 `RULES`（`--rules tabular`）。  
匹配前特征名会做规范化：去首尾空白、压缩空白、**转小写**，再跑正则。

完整一行一条正则：`tables/tabular_recovery_regex.csv`  
按任务汇总：`tables/tabular_recovery_regex_by_task.csv`

## 怎样算一个 fold「通过」

| 类型 | 条件 |
|---|---|
| control | top-5 里 **至少一条** 特征名命中任一 `hit` 正则 |
| spatial | 同上，并且 **Top-1 不能** 命中任一 `miss_as_top` 正则（丰度泄漏） |

任务级：15 fold 中 ≥50% 的 `lasso_passed` / `shap_passed` 为真 → raw 过。

下面 `|` 只是表格里分隔多条正则，源码里是独立的 `tuple` 项（或关系）。

## Control

### `motif_tumor_high`

| 角色 | 正则 |
|---|---|
| hit | `composition::.*tumor(?!_stroma)` |
| hit | `tissue_density::.*tumor(?!_stroma)` |
| hit | `tumor_area_ratio` |
| hit | `(^\|::)tumor$` |
| miss_as_top | （无） |

能过的例子：`composition::Tumor`、`density::tissue_density::Tumor`、`density::tumor_area_ratio`。  
注意：特征带 `composition::` / `density::` 前缀时，中间的 `tissue_density::` / `tumor_area_ratio` 子串仍可被命中。

### `motif_cd8_high`

| 角色 | 正则 |
|---|---|
| hit | `composition::.*(cd8\|t cells?)` |
| hit | `tissue_density::.*(cd8\|t cells?)` |
| hit | `(^\|::)(cd8 t cells?\|t cells?)$` |
| miss_as_top | （无） |

能过的例子：`composition::CD8 T cell`、`composition::CD8^+T`、`density::tissue_density::T cell`。  
过不了的典型：只有 `mixing::...`、只有无关组成亚型且名字里没有 cd8/t cell。

## Spatial

### `motif_cd8_clustering`

| 角色 | 正则 |
|---|---|
| hit | `cd8.*[_ ][kl]_r` |
| hit | `cd8.*ripley` |
| hit | `cd8.*cluster` |
| hit | `mixing::neighbor_fraction__.*(cd8\|t cell).*__to__.*(cd8\|t cell)` |
| hit | `mixing::same_type_neighbor_fraction` |
| hit | `mixing::local_mixing` |
| miss_as_top | `composition::.*cd8` |
| miss_as_top | `tissue_density::.*cd8` |
| miss_as_top | `frac__cd8` |
| miss_as_top | `(^\|::)cd8 t cell$` |

能过：`mixing::neighbor_fraction__CD8 T cell__to__CD8 T cell`、`point-pattern::CD8 T cell_L_r200`。  
Top-1 若是 `composition::CD8 T cell` → 即使后面有 mixing 也 **fail**。

### `motif_tumor_stroma_mixing`

| 角色 | 正则 |
|---|---|
| hit | `stroma.*[_ ][kl]_r` |
| hit | `tumor.*[_ ][kl]_r` |
| hit | `mixing` |
| miss_as_top | `composition::.*(tumor\|stroma)$` |
| miss_as_top | `tissue_density::.*(tumor\|stroma)$` |
| miss_as_top | `(^\|::)(tumor\|stroma)$` |

`hit` 里的 `mixing` 很宽：任意带 `mixing` 子串的列（含熵、邻域分数）都算命中。

### `motif_interface_immune`

| 角色 | 正则 |
|---|---|
| hit | `tumor_density::(cd8\|cd4\|b cell\|macrophage\|immune)` |
| hit | `interface` |
| hit | `(cd8\|cd4\|b cell).*[_ ][kl]_r` |
| miss_as_top | `composition::.*(cd8\|cd4\|b cell\|macrophage\|immune)` |
| miss_as_top | `tissue_density::.*(cd8\|cd4\|b cell\|macrophage)` |

注意：一般 `mixing` **不算** 这条的 hit（与其它接触/混合任务不同）。

### `motif_immune_exclusion`

| 角色 | 正则 |
|---|---|
| hit | `tumor_density::cd8` |
| hit | `tissue_density::cd8` |
| hit | `cd8.*[_ ][kl]_r` |
| miss_as_top | `composition::.*cd8` |
| miss_as_top | `frac__cd8` |

**没有** 把 `tumor_area_ratio` 放进 `miss_as_top`。因此 Top-1 是面积比、但 top-5 里有 `tissue_density::CD8...` 时仍可 `passed=True`——这就是 dens 族 / dens+pp 上 exclusion「虚高通过」的规则来源。`mixing` 也不在 hit 里。

### `motif_t_tumor_mixing`

| 角色 | 正则 |
|---|---|
| hit | `(t cells?\|cd8).*[_ ][kl]_r` |
| hit | `tumor.*[_ ][kl]_r` |
| hit | `mixing` |
| miss_as_top | `composition::.*(tumor\|t cells?\|cd8)` |
| miss_as_top | `tissue_density::.*(tumor\|t cells?\|cd8)` |

`tumor_density::T cell` **不是** hit（除非名字里另有 mixing / Ripley）；单族 density 上这条常 raw 失败，同列在 `cd8_tumor_contact` 上却可以过。

### `motif_cd8_tumor_contact`

| 角色 | 正则 |
|---|---|
| hit | `(cd8\|t cells?).*[_ ][kl]_r` |
| hit | `tumor_density::(cd8\|t cells?)` |
| hit | `mixing` |
| miss_as_top | `composition::.*(tumor\|cd8\|t cells?)` |
| miss_as_top | `tissue_density::.*(tumor\|cd8\|t cells?)` |

显式允许 `density::tumor_density::T cell` 这类列。

### `motif_macrophage_tumor_niche`

| 角色 | 正则 |
|---|---|
| hit | `macrophage.*[_ ][kl]_r` |
| hit | `tumor.*[_ ][kl]_r` |
| hit | `mixing` |
| miss_as_top | `composition::.*(tumor\|macrophage)` |
| miss_as_top | `tissue_density::.*(tumor\|macrophage)` |

### `motif_apc_t_contact`

| 角色 | 正则 |
|---|---|
| hit | `(apc\|dendritic\|macrophage).*[_ ][kl]_r` |
| hit | `(t cells?\|cd8).*[_ ][kl]_r` |
| hit | `mixing` |
| miss_as_top | `composition::.*(apc\|dendritic\|macrophage\|t cells?\|cd8)` |
| miss_as_top | `tissue_density::.*(apc\|dendritic\|macrophage\|t cells?\|cd8)` |

## 和 UTAG / GNN 规则的关系

- 上表只适用于 **tabular** 命名特征（Lasso / SHAP）。  
- UTAG message-passing 用另一套 `UTAG_RULES`（要 `utag_domain::` / marker 名）。  
- GNN Explainer 用的是细胞集合富集，不是这些正则。  
- UTAG portraits 用画像 hit_sets / markers，也不是这些正则。

## 读结果时怎么对上号

`fold_recovery.csv` 里：

- `lasso_hit` / `shap_hit`：top-5 中命中 `hit` 的特征名，`;` 拼接  
- `lasso_top1`：是否踩中 `miss_as_top` 看 `shap_miss_as_top1` 等列（Lasso 的 miss 体现在 `lasso_passed=False`）  
- `lasso_passed` / `shap_passed`：上表逻辑的布尔结果

同一套正则上的 hit@1 / P@k / MRR / `|H|` / recall@k 见 [12_tabular_topk_precision.md](12_tabular_topk_precision.md)。宽 `mixing` 正则会让 Jackson / METABRIC stroma 的 `|H|` 涨到整张 mixing 词表。
