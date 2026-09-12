# Globals：top-1 / top-3 precision

离线重算，**没有重跑模型**。源数据是 **clean** 标签上 10 个 `tabular_*/feature_ranks.csv`（4 单族 + 6 对），每 fold 存了 Lasso / SHAP 的 top-20。正则与 [10_recovery_regex.md](10_recovery_regex.md) 相同。脚本：`ExplainingConclusions/_compile_topk_precision.py`。

选定任务仍是面板那 10 条（Jackson 5、HNC 2、METABRIC 2、TNBC 1）。每个 (任务 × 特征集 × 解释器) 平均 15 个 fold。

表：`tables/tabular_topk_precision_by_task.csv`、`tables/tabular_topk_precision_method_summary.csv`、`tables/tabular_topk_hit_candidates.csv`。

## 指标

对每个 fold 的有序列表算完再按任务平均。

| 指标 | 定义 |
|---|---|
| `hit@1` | rank-1 是否命中任一 `rule.hit` |
| `precision@3` | top-3 里命中条数 / 3 |
| `precision@5` | top-5 里命中条数 / 5 |
| `MRR` | 若第 *r* 名首次命中则为 1/*r*，从未命中则为 0 |
| `|H|` | 该特征集全部列里，有多少列能匹配 `rule.hit` |
| `recall@k` | hits@k / min(*k*, \|H\|) |

`|H|=0` 时 precision 和 recall 都记 0：特征集里根本没有「合格名字」，不是解释器没排对。

## 为什么必须报 \|H\|

若干空间任务的 hit 正则里有一条裸的 `mixing`。mixing 族里几乎每一列都叫 `mixing::...`，Jackson 上 `|H|=407`，METABRIC stroma 上 `|H|=491`。这时 precision@3=1.00 只说明 top-3 都带了 `mixing::` 前缀，**不能**说明点到了 CD8–肿瘤接触。

反例：HNC `cd8_clustering` 的 mixing hit 更严（`neighbor_fraction__*cd8*__to__*cd8*`、`same_type_neighbor_fraction`、`local_mixing`），`|H|=9`。SHAP precision@3 = **0.356**，和 recall@3 相同——top-3 大约只捞到 1 个合格列。若只看「mixing 空间 precision@3=0.54」会被那 407 列的任务抬上去。

composition 在 8/10 个选定任务上 `|H|=0`（空间规则不要组成列）。这时 hit@1=0 是规则与词表不交，不是 Lasso 坏了。

## 方法总分（10 个选定任务平均，SHAP）

Lasso 与 SHAP 几乎重合，这里只列 SHAP。完整两套见 `tabular_topk_precision_method_summary.csv`。

| 特征集 | hit@1 | P@3 | P@5 | MRR | R@3 | 读法 |
|---|---:|---:|---:|---:|---:|---|
| composition | 0.10 | 0.13 | 0.12 | 0.14 | 0.13 | 只在 control 上有合格列 |
| density | 0.32 | 0.21 | 0.17 | 0.37 | 0.29 | tumor_area_ratio / 组织密度能救一部分 |
| mixing | 0.60 | 0.54 | 0.53 | 0.60 | 0.54 | 被宽 `mixing` 正则抬高 |
| point-pattern | 0.39 | 0.36 | 0.35 | 0.47 | 0.36 | 空间任务上更「诚实」 |
| composition + density | 0.17 | 0.19 | 0.19 | 0.27 | 0.27 | 联合后 control 仍在，空间仍缺名字 |
| composition + mixing | 0.61 | 0.54 | 0.54 | 0.65 | 0.54 | 近似单族 mixing；TNBC `cd8_high` 被 mixing 挤掉 |
| density + mixing | **0.69** | 0.57 | 0.54 | **0.70** | 0.57 | 这批里 hit@1 / MRR 最高 |
| composition + point-pattern | 0.46 | 0.46 | 0.44 | 0.57 | 0.46 | control + Ripley |
| density + point-pattern | 0.53 | 0.43 | 0.40 | 0.65 | 0.43 | 选定 control 两边都能点到名 |
| mixing + point-pattern | 0.60 | 0.57 | 0.56 | 0.60 | 0.57 | 宽 mixing + clustering 的 Ripley |

## 分任务（SHAP，四个有代表性的特征集）

数字是 15 fold 平均。`—` 表示 `|H|=0`。

### mixing

| 任务 | \|H\| | hit@1 | P@3 | R@3 | MRR |
|---|---:|---:|---:|---:|---:|
| Jackson 四条接触 / mixing | 407 | 1.00 | 1.00 | 1.00 | 1.00 |
| METABRIC tumor_stroma_mixing | 491 | 1.00 | 1.00 | 1.00 | 1.00 |
| HNC cd8_clustering | **9** | 1.00 | **0.36** | 0.36 | 1.00 |
| Jackson tumor_high | 0 | — | — | — | — |
| TNBC cd8_high | 0 | — | — | — | — |
| HNC immune_exclusion | 0 | — | — | — | — |
| METABRIC interface_immune | 0 | — | — | — | — |

clustering 的 Top-1 稳定是 CD8→CD8 邻域分数（MRR=1），但 top-3 里另外两条经常是别的邻域对，所以 P@3 掉到 0.36。这才是「排对了名字」的分辨率。

### composition

| 任务 | \|H\| | hit@1 | P@3 | MRR |
|---|---:|---:|---:|---:|
| TNBC cd8_high | 4 | 1.00 | 1.00 | 1.00 |
| Jackson tumor_high | 14 | 0.00 | 0.33 | 0.38 |
| 其余 8 条空间 | 0 | — | — | — |

TNBC control 的 CD8 组成列在 Top-1。Jackson `tumor_high` 合格列在表里（14 个），但 Top-1 经常不是它们（hit@1=0，MRR=0.38，大约排在第 3 名附近）。

### point-pattern

| 任务 | \|H\| | hit@1 | P@3 | MRR |
|---|---:|---:|---:|---:|
| HNC cd8_clustering | 10 | 1.00 | 1.00 | 1.00 |
| Jackson t_tumor_mixing | 150 | 0.87 | 0.78 | 0.93 |
| Jackson macrophage_tumor_niche | 150 | 0.80 | 0.84 | 0.90 |
| Jackson apc_t_contact | 20 | 0.60 | 0.44 | 0.77 |
| HNC immune_exclusion | 10 | 0.40 | 0.22 | 0.52 |
| Jackson cd8_tumor_contact | 10 | 0.20 | 0.24 | 0.47 |
| METABRIC interface_immune | 20 | 0.07 | 0.04 | 0.14 |
| 两个 control、stroma_mixing | 0 | — | — | — |

`tumor.*[_ ][kl]_r` 会匹配大量肿瘤亚型的 Ripley 列，所以 macrophage / T–tumor 的 `|H|` 有 150，P@3 偏乐观。clustering 的 `|H|=10` 且 P@3=1.00，是这张表里最干净的一条。

### density + point-pattern

协议总分里这是少数选定 control 2/2 且空间也能过的组合。top-k 也对得上：Jackson `tumor_high` hit@1=1.00，TNBC `cd8_high` hit@1=0.93，HNC clustering hit@1=1.00、P@3=1.00。interface / stroma 仍然接近 0（stroma 的 `|H|=0`，interface 的 P@3=0.04）。

### composition + mixing（对照「联合表把 control 名字挤掉」）

TNBC `cd8_high`：单族 composition 的 hit@1=1.00，加上 mixing 后掉到 **0.13**（MRR=0.28）。Jackson `tumor_high` 的 hit@1 仍是 0。空间四条接触则回到宽 `mixing` 的 1.00。这和协议里「mixing 进表，TNBC control 就不稳」是同一件事，只是现在能量化到 hit@1。

## HNC `cd8_clustering` 跨特征集（SHAP）

这是唯一一条 hit 正则既严、又在多族里都有名字的空间任务。

| 特征集 | \|H\| | hit@1 | P@3 | R@3 |
|---|---:|---:|---:|---:|
| composition / density / 二者联合 | 0 | 0 | 0 | 0 |
| mixing、composition+mixing | 9 | 1.00 | 0.36 | 0.36 |
| density+mixing | 9 | 1.00 | 0.44 | 0.44 |
| point-pattern、density+point-pattern | 10 | 1.00 | 1.00 | 1.00 |
| composition+point-pattern | 10 | 1.00 | 0.91 | 0.91 |
| mixing+point-pattern | 19 | 1.00 | 0.73 | 0.73 |

结论：clustering 的 Top-1 只要词表里有空间列就能稳住；要把 top-3 也装满合格列，Ripley 比邻域分数矩阵干净。

## 和原协议 pass 的关系

原协议只问 top-5 里有没有命中、以及 Top-1 是不是丰度泄漏。top-k precision 更细：

- **hit@1 / MRR** 回答「第一个名字对不对、排第几」。
- **P@3 / P@5** 回答「前几名里掺了多少无关列」，必须和 `|H|`、`R@k` 一起看。
- 空间任务上 composition 的 `miss_as_top1` 约 0.36：Top-1 经常是组成泄漏。mixing / point-pattern 上这条是 0。

不要单独用 10 任务平均的 P@3 给特征集排序。先丢掉 `|H|=0` 的任务，或改报 `R@k`，宽 `mixing` 正则那几条再单独注明。
