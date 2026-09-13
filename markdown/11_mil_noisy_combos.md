# MIL（noisy v2）特征组合：Attention + 实例解释器

对应 `logs/mil_interpreter_panel_noisy_all_combos.nohup.log` 串行跑完的六组：

| 顺序 | 特征组 | 输出目录 | 起止（日志） |
|---|---|---|---|
| 1 | composition | `mil_noisy_composition/` | 9/4 18:44 → 9/5 02:07 |
| 2 | mixing | `mil_noisy_mixing/` | 9/5 02:07 → 11:04 |
| 3 | celltype_density | `mil_noisy_celltype_density/` | 9/5 11:04 → 17:52 |
| 4 | composition + mixing | `mil_noisy_composition_mixing/` | 9/5 17:52 → 9/6 03:04 |
| 5 | composition + celltype_density | `mil_noisy_composition_celltype_density/` | 9/6 03:04 → 10:00 |
| 6 | mixing + celltype_density | `mil_noisy_mixing_celltype_density/` | 9/6 10:00 → 20:20 |

模型：Gated Attention MIL（window bag）。标签：noisy v2 伪标签。  
实例解释器：`attention`（基线）、`single`、`one_removed`、`ig`、`random`。

表：`tables/mil_noisy_combo_method_comparison.csv`，`tables/mil_noisy_combo_attention_wide.csv`，`tables/mil_noisy_combo_selected_metrics.csv`，`tables/mil_noisy_combo_selected_task_verdicts.csv`。

## 和 Global Lasso/SHAP 协议的差别

MIL **不**用特征名正则。每个 fold 要过三关：

1. **袋级 AUC**：control ≥ 0.90，spatial ≥ 0.60  
2. **定位（localization）**：window 分数在几何证据上的 top-10% 富集 ≥ 1.25  
3. **忠实性（faithfulness）**：MORF/LERF，打乱高分 window 应比随机更伤 AUC  

任务级：≥50% fold 过。control 任务要求 AUC+定位+忠实性都过（`pass_control`）。  
空间 raw：定位与忠实性都 ≥50% fold。  
空间 protocol：仍用 **同数据集匹配** `tumor_high`+`cd8_high` 都 `pass_control` 才计分（与 Global 文档同一套门；`protocol_spatial` 分母可变）。`raw_overall` / `protocol_overall` 分母固定为 10，被挡空间计失败。

选定 control 仍是 Jackson `tumor_high`、TNBC `cd8_high`。

## Attention 总分（主结果）

| 特征组合 | 选定 control | 空间 raw | 空间 protocol | overall raw | overall protocol | 一句话 |
|---|---|---|---|---|---|---|
| composition | **2/2** | **2/8** | **2/6** | **4/10** | **4/10** | 仅 HNC `cd8_clustering` / `immune_exclusion`；METABRIC 挡 |
| mixing | 2/2 | 0/8 | 0/4 | 2/10 | 2/10 | HNC `tumor_high` 定位失败 → HNC 空间全挡；空间 raw 也无双过 |
| celltype_density | 2/2 | 2/8 | 0/4 | 4/10 | 2/10 | HNC 两条 raw 过，但 HNC `tumor_high` 忠实性不够 → 协议挡 |
| composition + mixing | 2/2 | 1/8 | 0/4 | 3/10 | 2/10 | 仅 exclusion raw；HNC/METABRIC 匹配 control 不齐 |
| composition + celltype_density | **1/2** | 2/8 | 0/4 | 3/10 | 1/10 | TNBC `cd8_high` 定位 0.47 不过；HNC raw 过仍被挡 |
| mixing + celltype_density | 2/2 | 1/8 | 0/4 | 3/10 | 2/10 | 同 composition+mixing 量级 |

**没有一组 Attention 在空间 protocol 上超过 composition 的 2/6。** Jackson 四条空间在 Attention 下几乎全是「有 AUC / 有时 faithful，但定位不过」或反过来。

## Attention：选定任务明细

| 数据集 | 任务 | type | comp | mix | dens | c+m | c+d | m+d |
|---|---|---|---|---|---|---|---|---|
| Jackson | tumor_high | C | pass | pass | pass | pass | pass | pass |
| TNBC | cd8_high | C | pass | pass | pass | pass | **fail** | pass |
| Jackson | t_tumor_mixing | S | fail | fail | fail | fail | fail | fail |
| Jackson | cd8_tumor_contact | S | fail | fail | fail | fail | fail | fail |
| Jackson | macrophage_tumor_niche | S | fail | fail | fail | fail | fail | fail |
| Jackson | apc_t_contact | S | fail | fail | fail | fail | fail | fail |
| HNC | cd8_clustering | S | **pass** | fail* | pass* | fail* | pass* | fail* |
| HNC | immune_exclusion | S | **pass** | fail* | pass* | pass* | pass* | pass* |
| METABRIC | tumor_stroma_mixing | S | fail† | fail† | fail† | fail† | fail† | fail† |
| METABRIC | interface_immune | S | fail† | fail† | fail† | fail† | fail† | fail† |

`*` = raw 可能过，但协议 `fail_control_block`（该数据集匹配 control 未齐）。  
`†` = METABRIC 匹配 control 未齐（常见是 `cd8_high` 定位失败）。

## Attention AUC（选定任务）

| 任务 | comp | mix | dens | c+m | c+d | m+d |
|---|---:|---:|---:|---:|---:|---:|
| Jackson tumor_high | 0.971 | 0.940 | 0.978 | 0.943 | 0.981 | 0.949 |
| TNBC cd8_high | 0.965 | 0.972 | 0.967 | 0.973 | 0.969 | 0.974 |
| t_tumor_mixing | 0.801 | 0.663 | 0.821 | 0.678 | 0.832 | 0.663 |
| cd8_tumor_contact | 0.813 | 0.655 | 0.823 | 0.663 | 0.857 | 0.661 |
| macrophage_tumor_niche | 0.699 | 0.588 | 0.600 | 0.581 | 0.670 | 0.589 |
| apc_t_contact | 0.602 | 0.654 | 0.631 | 0.643 | 0.608 | 0.652 |
| cd8_clustering | 0.895 | 0.845 | 0.881 | 0.848 | 0.888 | 0.853 |
| immune_exclusion | 0.883 | 0.919 | 0.790 | 0.921 | 0.878 | 0.925 |
| tumor_stroma_mixing | 0.766 | 0.706 | 0.747 | 0.707 | 0.738 | 0.711 |
| interface_immune | 0.573 | 0.615 | 0.594 | 0.628 | 0.581 | 0.628 |

含 mixing 时 Jackson 接触/混合 AUC 反而低于 composition/density 单组——window mixing 特征没有帮 Attention 在这些任务上拉开差距，却容易把丰度 control 的定位打穿（HNC `tumor_high`）。

## 其它解释器（相对 Attention）

同一套 bag / 同一模型，换实例归因：

| 组合 | Attention 空间 protocol | IG | one_removed | single | random |
|---|---|---|---|---|---|
| composition | 2/6 | **5/6** | 4/6 | 2/6 | 0/8 |
| mixing | 0/4 | 2/6 | 2/6 | — | 0/8 |
| celltype_density | 0/4 | 4/6 | 4/6 | — | 0/8 |
| composition + mixing | 0/4 | 2/6 | 2/6 | — | 0/8 |
| composition + dens | 0/4 | **5/6** | **5/6** | 2/6 | 0/8 |
| mixing + dens | 0/4 | 2/6 | 2/6 | — | 0/8 |

（上表 `protocol_spatial` 分母随匹配 control 变化；`raw_overall` / `protocol_overall` 分母固定 10。完整数字见 `mil_noisy_combo_method_comparison.csv`。）

要点：

1. **`random` 全灭** → 定位/忠实性不是随机噪声。  
2. **Integrated Gradients（`ig`）在 composition（± density）上明显强于 Attention 权重**：空间 protocol 可到 5/6。说明「模型有信号」和「attention 权重对准几何证据」不是一回事。  
3. **`one_removed` 接近 IG**，`single` 接近 Attention。  
4. 即使 IG 协议分好看，也 **不要** 和 Global Lasso 的「点名 mixing 列」混为一谈：MIL 过的是 window 定位，不是列名回收。

## 怎么读

1. **Attention 主结果：composition 单组相对最好**（选定 control 2/2，HNC 两条空间 protocol 过），但 Jackson 空间接触类全军覆没。  
2. **加 mixing 没有救出 Jackson 空间**，还常搞挂 HNC/METABRIC 的匹配 control 定位。  
3. **若目标是「window 是否点在对的几何区域」**，应并列报告 IG / one_removed，不能只报 Attention。  
4. 与 Global 线对比：Global composition+mixing 能**写出**邻域列名；MIL composition 只能在 HNC 上证明 **bag 预测 + 部分定位**。两条线回答的问题不同。

源码：`benchmark/motifs/mil_run.py`，`scripts/verify_pseudo_label_explanations.py --mode mil`。
