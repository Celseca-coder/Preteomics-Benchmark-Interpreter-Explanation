# Global 解释方法总结

本目录整理 **Global（区域级）** 解释器，对应 `TME_modeling_benchmark` 中与 MIL 相对的一条线。Global 把每个 region 压成固定长度向量或节点归因，再解释；MIL 则解释局部 window bag。

覆盖四种核心方法，以及它们的对照方法：

| 方法 | 角色 | 输入 | 解释输出 |
|---|---|---|---|
| **Lasso** | 带 L1 的线性分类器，顺带特征选择 | 命名表格特征（组成/密度/混合/点过程） | 系数非零的特征及方向 |
| **SHAP** | 对 Ridge 线性模型做精确线性归因 | 同一套命名表格特征 | 每个特征对 held-out 预测的贡献 |
| **UTAG-Linear embedding** | UTAG message-passing 的区域均值/标准差 + 线性模型 | 平滑后的 marker 统计 | 哪些 `utag_mean__` / `utag_std__` 被选中 |
| **GNN Explainer** | SPACE-GM 图模型的节点归因 | 细胞图 | 每个细胞的 importance |

对照方法（相应方法）：

- **UTAG native portraits**：用 domain 画像解释，不再经过 Lasso/SHAP
- **KRONOS / Eva linear probe**：冻结 region embedding + Ridge，只能测信号，不能点名特征
- **UTAG domains + Lasso/SHAP**：HNC 上的 domain 比例线性模型

结果表在 `tables/`，分方法说明在 `markdown/`。

---

## 核心结论

1. **按选定 control 门控，Lasso / SHAP 并没有过关。** 面板指定的两条 control 是 Jackson `tumor_high`（过）和 TNBC `cd8_high`（不过），所以是 **1/2，方法级 `control_ok=False`**。不能据此写成「已经验证过的空间解释器」。后面的空间 5/6 用的是另一套门：同数据集匹配的 `tumor_high`+`cd8_high`。Jackson 上非选定的 `cd8_high` 是过的，TNBC 这条选定 control 没过。
2. **它们唯一多出来的能力是：在同数据集丰度 control 过关的数据集上，能写出 mixing / Ripley 这种特征名。** Jackson 四条空间和 HNC `cd8_clustering` 的 Top-1 是邻域分数，不是细胞比例。这是相对 UTAG / embedding / GNN 的差别，不是「选定 control 已通过」。UTAG portraits 反而是选定 control **2/2 都过**，但空间几乎全失败。
3. **丰度 control 和空间 motif 必须分开看。** 模型 AUC 高不等于解释对。HNC `immune_exclusion` 的 Lasso/SHAP AUC ≈ 0.98，但 Top-1 是 `density::tumor_area_ratio`。TNBC `cd8_high` AUC ≈ 0.98，但 Top-1 是 mixing 而不是 `composition::CD8`。
4. **UTAG-Linear embedding 适合丰度，不适合空间。** 用正确的 UTAG 回收规则重打分后：HNC `tumor_high` / `cd8_high` 回收通过，四个空间 motif 全部失败。原因是空间规则要的是 `utag_domain::`，而 linear embedding 只有 marker 的 mean/std。
5. **冻结 embedding（KRONOS、Eva）有预测信号，但过不了解释协议。** Eva 在 8/8 空间任务上 AUC ≥ 0.60，但 control 阈值是 0.90，同一数据集两个丰度 control 很少同时过，空间分数被全部挡住。Embedding 也没有可回收的特征名。
6. **GNN Explainer 目前不能当作已验证的 Global 解释器。** 导出用的是每个任务 seed0/fold0 一个 checkpoint。两个选定 control 的 region 通过率都 &lt; 0.50（TNBC `cd8_high` 最接近，0.477）。空间任务最高也只有 HNC 的 ~0.46。节点 importance 没有稳定富集到 motif 对应细胞类型。
7. **四族拆开跑之后，联合表的 1/2 control 和 5/6 空间可以对上，但没有单族通过完整协议。** composition / density 单独时选定 control 是 **2/2**（TNBC `cd8_high` 回到 CD8 组成或组织密度）；空间 raw 分别是 0/8 和 2/8（density 的两条通过要按规则打折）。mixing / point-pattern 空间 raw 都是 6/8，但 control 0/2，协议 0/8。联合 L1 有 mixing 时不会回到 `composition::CD8`。明细见 [markdown/08_tabular_families.md](markdown/08_tabular_families.md)。
8. **两两组合里，联合四族 ≈ composition + mixing（同为 1/2 + 5/6）。** 真正同时拿到选定 control 2/2 和空间协议分的是 **composition + point-pattern（5/6）** 和 **density + point-pattern（6/8）**；后者 exclusion 仍是 `tumor_area_ratio` 规则漏洞。mixing 一进表，TNBC `cd8_high` 名字就不稳。明细见 [markdown/09_tabular_pairs.md](markdown/09_tabular_pairs.md)。
9. **协议 pass 只问 top-5 有没有命中；top-k precision 更细，但必须带 `|H|`。** 10 任务平均上 density+mixing 的 SHAP hit@1 / MRR 最高（0.69 / 0.70），但 mixing 相关组合被裸 `mixing` 正则抬高（Jackson `|H|=407`）。诚实对照是 HNC `cd8_clustering`：mixing `|H|=9`、P@3=0.36；point-pattern `|H|=10`、P@3=1.00。composition 只在 control 上有合格列。明细见 [markdown/12_tabular_topk_precision.md](markdown/12_tabular_topk_precision.md)。

---

## 10 个选定 motif

每个 motif 只在其来源数据集上评分。同一数据集会附带 `tumor_high` / `cd8_high` 作为匹配 control。

| motif | 数据集 | 类型 | 定义（残差化后） |
|---|---|---|---|
| tumor_high | BC-Jackson2020 | control | 肿瘤比例 ≥ 发现队列中位数 |
| cd8_high | TNBC-Wang2023 | control | CD8（或 T 细胞代理）比例 ≥ 中位数 |
| t_tumor_mixing | BC-Jackson2020 | spatial | T 细胞与肿瘤邻域接触偏高 |
| cd8_tumor_contact | BC-Jackson2020 | spatial | CD8/T 与肿瘤接触偏高 |
| macrophage_tumor_niche | BC-Jackson2020 | spatial | 巨噬细胞–肿瘤邻域偏高 |
| apc_t_contact | BC-Jackson2020 | spatial | APC–T 接触偏高 |
| cd8_clustering | HNC-Wu2022 | spatial | CD8 自邻域偏高 |
| immune_exclusion | HNC-Wu2022 | spatial | CD8 落在肿瘤多边形外偏多 |
| tumor_stroma_mixing | BC-METABRIC | spatial | 肿瘤–间质邻域偏高 |
| interface_immune | BC-METABRIC | spatial | 免疫细胞集中在肿瘤边界 50 µm 内 |

完整定义见 `tables/selected_motif_catalog.csv`。

---

## 评分协议（读表时必须分清两列）

**原始通过（raw）**：该任务自己的回收/AUC/富集是否过阈值。

两套 control 不要混：

1. **选定 control（方法级，2 条）**：Jackson `tumor_high`、TNBC `cd8_high`。这决定表里的 `protocol_control` 和 `control_ok`。Lasso/SHAP 是 1/2，`control_ok=False`。
2. **同数据集匹配 control（空间门控）**：每个数据集自己的 `tumor_high`+`cd8_high`（包括非选定的那条）。这决定某条空间任务会不会被标 `fail_control_block`。

**协议通过（protocol）**：空间任务用的是第 2 套。任一匹配 control 失败，该数据集全部空间任务记为 `fail_control_block`，即使空间回收本身是对的。因此 Lasso 可以在选定 control 只有 1/2 的情况下，仍在 Jackson/HNC 上拿到空间 protocol 分。这不能回写成「选定 control 已通过」。

阈值：

| 方法 | control | spatial |
|---|---|---|
| Lasso / SHAP | ≥50% fold 回收到预期命名特征 | 同上，且 Top-1 不能是已知丰度泄漏 |
| Embedding probe | 均值 AUC ≥ 0.90 | 均值 AUC ≥ 0.60 |
| UTAG portraits | domain 画像命中预期细胞/marker | 命中空间画像且不是纯肿瘤丰度 |
| GNN Explainer | ≥50% region 的 top 10% 节点对预期细胞类型富集 ≥ 1.25 | 同左 |

因此表中会出现：Eva 空间 raw 8/8、protocol 0/8。这不是计算错误。

各 `*_method_comparison.csv` 另有两列**固定分母**总分，被匹配 control 挡住的空间任务计失败、不移出分母：

- `raw_overall = (控制通过 + 空间 raw 通过) / 10`
- `protocol_overall = (控制通过 + 空间协议通过) / 10`

选定面板永远是 2 个 control + 8 个空间。control 比空间容易过，2/10 可以全是丰度任务；`raw_overall` 也会被「空间 raw 过、匹配 control 没过」抬高（Eva 9/10 vs protocol 1/10）。**只用于方法排序，分项必须并排看。** 恒有 `raw_overall ≥ protocol_overall`。

---

## 方法总分

### Global 主方法（联合四族表格 = Lasso / SHAP）

| 方法 | 选定 control | 空间 raw | 空间 protocol | overall raw | overall protocol | 一句话 |
|---|---|---|---|---|---|---|
| Lasso（四族联合） | 1/2 | 6/8 | 5/6 | 7/10 | 6/10 | 选定 control 未齐；匹配 control 过关的数据集上能点名 mixing |
| SHAP（四族联合） | 1/2 | 6/8 | 5/6 | 7/10 | 6/10 | 与联合 Lasso 同一套 pass/fail |
| UTAG-Linear（HNC，重打分） | 2/2（仅 HNC） | 0/4（仅 HNC） | — | — | — | 能找回 marker，找不回空间 domain（只在 HNC 四任务上打分，不进 10 任务 overall） |
| UTAG portraits | 2/2 | 1/8 | 0/4 | 3/10 | 2/10 | control 过，空间几乎全是肿瘤 domain |
| KRONOS linear | 1/2 | 5/8 | 0/8 | 6/10 | 1/10 | 有弱空间信号，解释协议不过 |
| Eva linear | 1/2 | 8/8 | 0/8 | **9/10** | 1/10 | 空间 AUC 最好的 embedding，仍无特征名；overall raw 被空间 raw 抬高 |
| GNN Explainer | 0/2 | 0/8 | 0/8 | 0/10 | 0/10 | 节点归因未过 50% region 阈值 |

### 单族命名特征（同一次 `--mode tabular`，Lasso / SHAP 分列）

| 特征族 | Lasso control | Lasso 空间 raw | Lasso 空间 protocol | Lasso overall | SHAP control | SHAP 空间 raw | SHAP 空间 protocol | SHAP overall |
|---|---|---|---|---|---|---|---|---|
| composition | 2/2 | 0/8 | 0/6 | 2/10 · 2/10 | 2/2 | 0/8 | 0/6 | 2/10 · 2/10 |
| density | 2/2 | 2/8 | 2/8 | 4/10 · 4/10 | 2/2 | 2/8 | 2/8 | 4/10 · 4/10 |
| mixing | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 |
| point-pattern | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 |

单族上 Lasso 与 SHAP **任务级完全一致**。`overall` 列是 `raw_overall · protocol_overall`，分母固定 10。

### 两两组合（同上，Lasso / SHAP 分列）

| 组合 | Lasso control | Lasso 空间 raw | Lasso 空间 protocol | Lasso overall | SHAP control | SHAP 空间 raw | SHAP 空间 protocol | SHAP overall |
|---|---|---|---|---|---|---|---|---|
| composition + density | 2/2 | 2/8 | 2/8 | 4/10 · 4/10 | 2/2 | 2/8 | 2/8 | 4/10 · 4/10 |
| composition + mixing | 1/2 | 6/8 | 5/6 | 7/10 · 6/10 | **0/2** | 6/8 | **1/2** | 6/10 · 1/10 |
| density + mixing | 1/2 | 6/8 | 0/8 | 7/10 · 1/10 | 1/2 | 6/8 | 0/8 | 7/10 · 1/10 |
| composition + point-pattern | 2/2 | 5/8 | 5/6 | 7/10 · 7/10 | 2/2 | 5/8 | 5/6 | 7/10 · 7/10 |
| density + point-pattern | 2/2 | 6/8 | 6/8 | **8/10 · 8/10** | 2/2 | 6/8 | 6/8 | **8/10 · 8/10** |
| mixing + point-pattern | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 | 0/2 | 6/8 | 0/8 | 6/10 · 0/10 |

`overall` 列同样是 `raw_overall · protocol_overall`。按 protocol overall 排序，最高是 density + point-pattern（8/10），其次 composition + point-pattern（7/10）。唯一分叉：`composition + mixing` 上 Jackson `tumor_high`（Lasso 回收 0.67 过、SHAP 0.47 不过）→ SHAP 选定 control 变 0/2，Jackson 空间被匹配 control 挡住，协议分母变成 2，protocol overall 从 6/10 掉到 1/10。其余五对与对应 Lasso 一致。

单族 / 两两不是把联合表拆列，而是各自重跑。明细：`tables/global_method_comparison.csv`，`tables/tabular_family_method_comparison.csv`，`tables/tabular_pair_method_comparison.csv`。说明见 [08](markdown/08_tabular_families.md)、[09](markdown/09_tabular_pairs.md)。

### MIL noisy 特征组合（Attention；window 定位协议）

对应 `logs/mil_interpreter_panel_noisy_all_combos.nohup.log` 六组已跑完。协议是定位+忠实性，不是列名正则。明细：[11_mil_noisy_combos.md](markdown/11_mil_noisy_combos.md)。

| 特征组合 | 选定 control | 空间 raw | 空间 protocol | overall raw | overall protocol |
|---|---|---|---|---|---|
| composition | 2/2 | 2/8 | 2/6 | 4/10 | 4/10 |
| mixing | 2/2 | 0/8 | 0/4 | 2/10 | 2/10 |
| celltype_density | 2/2 | 2/8 | 0/4 | 4/10 | 2/10 |
| composition + mixing | 2/2 | 1/8 | 0/4 | 3/10 | 2/10 |
| composition + celltype_density | 1/2 | 2/8 | 0/4 | 3/10 | 1/10 |
| mixing + celltype_density | 2/2 | 1/8 | 0/4 | 3/10 | 2/10 |

Attention 最强仍是 composition（仅 HNC 两条空间 protocol 过）。同设定下 IG 在 composition 上可达 5/6，说明 attention 权重常对不齐几何证据。

---

## 文件索引

| 文件 | 内容 |
|---|---|
| [markdown/01_protocol.md](markdown/01_protocol.md) | 协议、特征、门控规则 |
| [markdown/02_lasso.md](markdown/02_lasso.md) | Lasso |
| [markdown/03_shap.md](markdown/03_shap.md) | 线性 SHAP |
| [markdown/04_utag_linear_embedding.md](markdown/04_utag_linear_embedding.md) | UTAG message-passing + 线性模型 |
| [markdown/05_gnn_explainer.md](markdown/05_gnn_explainer.md) | SPACE-GM GNN Explainer |
| [markdown/06_corresponding_methods.md](markdown/06_corresponding_methods.md) | UTAG portraits、KRONOS、Eva、UTAG domains |
| [markdown/07_clinical_global.md](markdown/07_clinical_global.md) | 真实临床终点上的 Global Lasso / SHAP |
| [markdown/08_tabular_families.md](markdown/08_tabular_families.md) | 四族拆开：composition / density / mixing / point-pattern |
| [markdown/09_tabular_pairs.md](markdown/09_tabular_pairs.md) | 六对两两组合 |
| [markdown/10_recovery_regex.md](markdown/10_recovery_regex.md) | 各 motif 回收通过的 hit / miss 正则 |
| [markdown/11_mil_noisy_combos.md](markdown/11_mil_noisy_combos.md) | MIL noisy 六组特征组合（Attention / IG 等） |
| [markdown/12_tabular_topk_precision.md](markdown/12_tabular_topk_precision.md) | Globals top-1 / top-3 precision、MRR、`|H|`、recall@k（clean tabular，不重跑） |
| `tables/` | 上述结论对应的 CSV |

源数据主要来自：

`TME_modeling_benchmark/results/pseudo_label_explanations_panel/`
