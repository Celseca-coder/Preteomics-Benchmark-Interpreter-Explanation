# ExplainingConclusions TODO

解释器协议的后续工作清单。不含排期。标签口径已定为 **noisy**（`--labels .../PseudoNoisyDataset/per_dataset/{dataset}_v2_noisy.csv`）。

代码仓库：`TME_modeling_benchmark`，当前分支 `interpreter/MIL`。结论目录：本目录。

---

## 0. 标签口径：后续一律 noisy

**上周 `mil_ig_*` 三组已经是 noisy**，log 叫 `mil_ig_features_*.log`，命令里带了 `--labels .../{dataset}_v2_noisy.csv`。那批结果不用因为「clean」再重跑。

容易混的是目录名，不是标签本身：

| 目录 / log | 实际标签 |
|---|---|
| `mil_ig_mixing` 等、`logs/mil_ig_features_*.log` | **noisy**（上周 IG 特征归因） |
| `mil_noisy_*`、`markdown/11_mil_noisy_combos.md` | **clean**（名字骗了人：脚本没传 `--labels`，走了默认 `{dataset}_v2.csv`） |
| 现有 `tabular_*` / embedding | **clean** |
| TODO 里曾写的 `mil_ig_clean_*` | 不要用；是误把口径定成 clean 时写的重跑命令 |

和旧 `mil_noisy_*` 比，noisy 上控制任务 bag AUC 从 ~0.96 掉到 ~0.86、过 0.90 门的 fold 从 0.94 掉到 0.13，这是噪声标签的真实效应，不是跑错了。特征级在剔除被门挡住的 fold 后仍可读（`abs` recovered 0.61–0.63）。协议 summary 里全是 `fail_control_block`，因为 AUC 门仍按 clean 的 0.90 设。noisy 线上要么下调控制 AUC 门，要么主报 raw 特征恢复、不让 control 把空间任务整表抹掉。

要做：

1. 把 working tree 里 5 个未提交文件先 commit。
2. **之后所有新跑（tabular 两两、MIL、IG）都显式传 noisy 标签**，不要省略 `--labels`。
3. 纠正文档：`mil_noisy_*` 实际是 clean；`11_mil_noisy_combos.md` 标题改口径。上周 `mil_ig_*` 保持 noisy 命名。
4. `.gitignore`：现有规则只匹配 `mil_noisy*`，`mil_ig_*` 不在内。

若还要再跑 IG（例如加上 3.2 的连续证据列），命令必须带 noisy 标签：

```bash
cd /autofs/nas8/tywang/tjzou/TME_modeling_benchmark
nohup bash -c '
PY=/autofs/nas8/tywang/tjzou/Miniconda3/envs/p3/bin/python
export CUDA_VISIBLE_DEVICES=0
for combo in "composition mixing" "mixing celltype_density"; do
  tag=$(echo "$combo" | tr " " "_")
  echo "======== START $tag $(date) ========"
  "$PY" -u scripts/verify_pseudo_label_explanations.py \
    --panel --mode mil --mil-explainers ig \
    --feature-groups $combo --mil-feature-aggs abs signed top \
    --labels "/autofs/nas8/tywang/tjzou/PseudoNoisyDataset/per_dataset/{dataset}_v2_noisy.csv" \
    --label-version v2 \
    --data-root /autofs/bal14/zqwu/CellularTables/TME_benchmark_data \
    --device cuda --seeds 0 1 2 \
    --output-dir "results/pseudo_label_explanations_panel/mil_ig_${tag}" \
    > "logs/mil_ig_features_${tag}.log" 2>&1
  echo "======== DONE $tag $(date) ========"
done
' > logs/mil_ig_features_all_combos.nohup.log 2>&1 &
```

注意：上周这两组已经跑完 noisy。只有改了定位 metric 或覆盖旧目录时才需要再发。

---

## 1. control / spatial 合并汇总（不分控制任务和空间任务）

现状：`ExplainingConclusions/_compile_tables.py` 的 `score_method()` 只吐分离分数（`protocol_control`、`raw_spatial`、`protocol_spatial`）。`protocol_spatial` 的分母会随被拦任务浮动。

要做：在 summary dict 里加两列，**固定分母**，被 control 拦掉的空间任务计为失败而不是移出分母：

- `raw_overall = (控制通过 + 空间 raw 通过) / (控制数 + 空间数)`
- `protocol_overall = (控制通过 + 空间协议通过) / (控制数 + 空间数)`

分项必须并排保留。面板里 control 只有 `tumor_high` / `cd8_high` 两类但覆盖 4 个数据集，合并分会被 control 抬高，只用于排序。文档写清这一点。

验收：`global_method_comparison.csv`、`tabular_family_method_comparison.csv`、`tabular_pair_method_comparison.csv`、`mil_noisy_combo_method_comparison.csv`（以及后续 MIL IG 表）都出现新列，且恒有 `raw_overall ≥ protocol_overall`。

---

## 2. Globals：top-1 / top-3 precision

完全不用重跑。10 个 `tabular_*/feature_ranks.csv` 已存每 fold 的 top-20（lasso 和 shap）。

对每个 fold 算完再按 task 平均：

| 指标 | 定义 |
|---|---|
| `hit@1` | rank-1 是否命中 `rule.hit`（top-1 metric） |
| `precision@3` | top-3 命中数 / 3 |
| `precision@5` | top-5 命中数 / 5（已有 `n_hit_in_top`） |
| `MRR` | 1 / `best_hit_rank`（区分「第 2 名命中」和「第 5 名险过」） |

必须同时输出每个 `(task, 特征集)` 的候选命中数 `|H|`，以及归一化 `recall@k = hits@k / min(k, |H|)`。否则 precision@k 被「该特征集里到底有几个合格列」卡死：例如 `motif_cd8_clustering` 在 by-type mixing 里可能只匹配到 1 列 `neighbor_fraction__*cd8*__to__*cd8*`，precision@3 上限就是 0.33，跨 task 横比没有意义。

**已写完（clean tabular，不重跑）**：`markdown/12_tabular_topk_precision.md`；表在 `tables/tabular_topk_precision_*.csv`、`tabular_topk_hit_candidates.csv`。后续若补 `02`/`03`/`08`/`09` 正文里的分列表，以这批 CSV 为准。MIL IG 文档改用 `13_mil_ig_features.md`，不要占 12。

---

## 3. MIL

### 3.1 特征重要性（= 原第 4 项，进行中）

代码已完成：`ig_attributions` 保留 `(window, feature)` 轴；三种聚合（`abs` / `signed` / `top`）；`MIL_FEATURE_RULES`；`summarize_mil_features`。验证方法和全局重要特征协议相同：top-k 命中 `hit` 正则，空间任务 Top-1 不能是纯局部丰度。

上周 noisy 三组已落盘。剩下的是汇总成 `markdown/13_mil_ig_features.md`（verdict 列在 noisy 上会被 control AUC 门挡住，主看 raw `frac_recovered` / `frac_passed`）。若 3.2 (c)(e) 加列后再跑，必须带 `--labels` noisy。`12` 已给 Globals top-k precision。

主推 `abs`。`signed` / `top` 作敏感性分析放附录。已有观察：剔除被 AUC 门挡住的 fold 后，`signed` 恢复率大约只有 `abs` 的一半（0.32 vs 0.62），说明同一 motif 在不同窗里的证据方向可以相反。

### 3.2 更强的空间定位 ground truth

现状：二值 `pos`/`neg` 掩码 + top-10% 富集 ≥ 1.25 单点阈值。标签是 motif 已知的，可以做更细的 metric。

按性价比：

**(a) 修口径不一致（先做，但会改 pass 判据，需要重跑）。** 标签侧用置换零模型校正过（`detect.py` 的 `_permute`；`motif_source_assignment.csv` 写着 "After residualizing on CD8%..."）。窗口证据侧是裸阈值 `mix >= 0.15`，没有丰度校正。在窗级做同样校正（对局部丰度回归取残差，或按窗内细胞数做置换），pos 窗才对应「超出丰度预期的聚集」。修完后 `faithful_but_wrong_place` 那批判决可能翻盘。

**(b) 尺度不匹配。** motif 邻域 50 µm，window 固定 100 µm、step 50 µm。至少做 window size 敏感性（50 / 100 / 200 µm）。

**(c) 用上已算出却被扔掉的连续证据。** `evidence_for_bag` 返回连续 `score`，`localization_metrics` 只用了二值 `pos`/`neg`。可加 Spearman ρ(解释器分数, 证据分数)、以证据分数为 gain 的 nDCG@k、加权富集。只改一个函数，不重提特征。**纯加列、不动 pass 判据。**

**(d) 判别性对照。** 用 A motif 的解释器分数去富集 B motif 的证据窗。如果一样富集，说明只是找到了「细胞多的窗」。目前协议唯一缺的阴性对照。

**(e) 单点换曲线。** 同时报 top-5% / 10% / 20%。**纯加列。**

(c) 和 (e) 若要进下一轮 IG，必须带 noisy `--labels`。 (a)(b)(d) 改判据或要额外实验，单独立项。

---

## 4. 生成并验证特征重要性解释器（进行中）

与 3.1 是同一件事。全局协议：ranked 特征名 → `rank_recovery` 对 `hit` / `miss_as_top`。MIL 侧 IG 已接到同一套规则（窗口特征用 `__` 分隔）。第 0 项重跑完成后，把 `ig_feature_summary.csv` 按 `abs` 主推写进结论。

---

## 5. Global 方法两两组合跑满伪标签数据集（进行中）

不是临床终点 pairwise（`feature/MethodCombination` 上的 `run_pairwise_feature_combinations.py` 跑的是 `ds.task_ids` / pCR，没有伪标签接线，且该分支删了 motif 基础设施）。也不移植那个脚本。

做法：扩 `verify_pseudo_label_explanations.py`。tabular 模式已经同时输出 `auc_lasso` / `auc_ridge` 和解释协议，fold / 标签口径天然一致。

1. `_featurizer()` 加 `expression` 分支（`MeanExpressionFeaturizer`），凑齐 6 组：composition / density / expression / spatial-distance / point-pattern / mixing。`spatial-distance` 代码已支持，从未被用过。
2. 扩 `scripts/run_tabular_family_combos.sh`：`FAMILIES` 到 6 个，pairs 到 15。已有 4 单族 + 6 对；还差单族 `expression`、`spatial-distance`，以及 9 对。point-pattern 相关组合排最后（现有脚本已如此）。
3. 每个组合约 8–11 分钟 CPU，不抢 GPU。新跑一律加 `--labels .../{dataset}_v2_noisy.csv`。已有 10 个 `tabular_*` 是 clean，不要和 noisy 新结果混在同一目录。
4. `compile_tabular_pairs` 的 `PAIR_DIRS` 从 6 扩到 15；更新 `markdown/09_tabular_pairs.md`。

---

## 已砍 / 不在本清单

- **再跑一遍 clean IG**：口径已改回 noisy；上周 `mil_ig_*` 就是 noisy。
- **`mixing` 单组 IG 特征级重跑**：词表与 control hit 正则无交集。
- **从 `feature/MethodCombination` 移植 pairwise 脚本**：validation API 分叉，且不接伪标签。

## 文档产出对照

| 项 | 主要产出 |
|---|---|
| 0 | 后续命令一律 noisy；文档纠正 `mil_noisy_*` 其实是 clean |
| 1 | 各 `*_method_comparison.csv` 增加 `raw_overall` / `protocol_overall`；README 方法总分表 |
| 2 | `markdown/12_tabular_topk_precision.md` + `tables/tabular_topk_*`（clean；`02`/`08`/`09` 已加交叉引用） |
| 3.1 / 4 | `markdown/13_mil_ig_features.md` |
| 3.2 (c)(e) | MIL fold / summary 新列；定位文档补连续证据与多点 k |
| 5 | 15 对 tabular 结果；更新 `09_tabular_pairs.md` |
