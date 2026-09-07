# top100_network GO-MF TF 分析流程

## 1. 输入与目标

- 输入会话：`top100_network.pre_go_tf.cys`
- 物种：`Mus musculus`
- 原始 GeneMANIA 网络：117 个节点、2,229 条综合证据边
- 候选基因：原会话中 `node type=query` 的 97 个基因
- 分析目标：重新核实严格意义上的 DNA 结合型转录因子，更新 Cytoscape 节点表和视觉样式，并构建候选基因诱导互作网络。

## 2. TF 数据库与判定规则

使用的数据库类型为 **Gene Ontology - Molecular Function**：

- 小鼠基因注释：`org.Mm.eg.db 3.17.0`
- GO 结构与术语：`GO.db 3.17.0`
- TF 根条目：`GO:0003700 DNA-binding transcription factor activity`
- 纳入规则：基因具有 `GO:0003700` 或其任一后代 Molecular Function 注释
- 排除规则：仅作为转录共因子、染色质调节因子或一般转录调节因子的基因，不自动归为严格 TF

该标准强调可复现性和 DNA 结合型 TF，可能漏掉尚未被 GO 完整注释的转录因子。

## 3. 分析流程

1. 从 `.cys` ZIP 会话中读取 GeneMANIA 节点表、边表、网络 XGMML 和视图 XGMML。
2. 提取 97 个 `query` 候选基因。
3. 使用 `org.Mm.eg.db` 查询每个基因的 GO Molecular Function 注释。
4. 使用 `GO.db` 获得 `GO:0003700` 的全部后代条目。
5. 将符合规则的基因标记为 `TF`，其他候选基因标记为 `DEG`；GeneMANIA 添加的 20 个结果节点标记为 `BACKGROUND` 并在最终样式中隐藏。
6. 构建 97 个候选基因的诱导子图，保留它们之间的 1,181 条证据边。
7. 计算每个候选基因在诱导网络中的度数 `degree_val`。
8. DEG 节点大小固定为 8；TF 节点大小定义为：

   `display_size = 42 + 0.9 x degree_val`

9. 应用力导向布局和参考图样式：TF 为金黄色圆形，DEG 为紫红色小方块，互作边为浅蓝色，仅显示 TF 标签。
10. 更新 `.cys` 节点表、视觉映射、视图坐标和会话缩略图，并输出 PDF/PNG。

## 4. 节点表新增字段

- `node_type`：`TF`、`DEG` 或 `BACKGROUND`
- `tf_database_type`
- `tf_database_version`
- `tf_reference_rule`
- `tf_go_terms`
- `tf_go_term_names`
- `tf_evidence_codes`
- `tf_classification_date`
- `degree_val`
- `display_size`
- `label_text`

## 5. 可复现脚本

### 5.1 TF 分类

```bash
Rscript scripts/classify_tf_go.R \
  top100_network.pre_go_tf.cys \
  output/data/top100_tf_go_classification.csv
```

### 5.2 修改 Cytoscape 会话并生成报告

```bash
PYTHON=/Users/andychen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3

$PYTHON scripts/build_go_tf_artifacts.py \
  top100_network.pre_go_tf.cys \
  output/data/top100_tf_go_classification.csv \
  --output-cys top100_network.cys \
  --backup-cys top100_network.pre_go_tf.cys \
  --preview-png output/images/top100_network_GO_TF_preview.png \
  --report-pdf output/pdf/top100_network_GO_TF_analysis.pdf
```

### 5.3 单独导出网络 PDF

```bash
$PYTHON scripts/export_network_pdf.py \
  output/images/top100_network_GO_TF_preview.png \
  output/pdf/top100_network_GO_TF_network.pdf
```

## 6. TF 核实结果

| TF | 度数 | 度中心性 | 节点大小 | GO 条目 | 证据代码 |
|---|---:|---:|---:|---|---|
| Sox2 | 36 | 0.375 | 74.4 | GO:0000981; GO:0001228; GO:0003700 | IBA; IDA; ISO |
| Pax9 | 34 | 0.354 | 72.6 | GO:0000981; GO:0001228; GO:0003700 | IBA; IDA; TAS |
| Trp63 | 29 | 0.302 | 68.1 | GO:0000981; GO:0001228; GO:0003700 | IBA; IDA; ISO |
| Tcf4 | 26 | 0.271 | 65.4 | GO:0000981; GO:0001228; GO:0003700 | IBA; IDA; ISO |
| Sox21 | 17 | 0.177 | 57.3 | GO:0000981; GO:0001228; GO:0003700 | IBA; IDA |
| Etv1 | 5 | 0.052 | 46.5 | GO:0000981; GO:0001228; GO:0003700 | IBA; ISO; TAS |

证据代码：`IDA` 为直接实验；`IBA` 为基于系统发育的生物学属性推断；`ISO` 为基于序列正交性的推断；`TAS` 为可追溯作者陈述。

此前标注中的 `Nupr1` 和 `Eya1` 未满足严格的 `GO:0003700` 判据。`Etv1`、`Pax9` 和 `Sox21` 被新增为 TF。

## 7. 互作网络结果

| 互作类型 | 边数 | 占比 |
|---|---:|---:|
| Co-expression | 696 | 58.9% |
| Predicted | 281 | 23.8% |
| Co-localization | 175 | 14.8% |
| Other | 23 | 1.9% |
| Physical Interactions | 5 | 0.4% |
| Genetic Interactions | 1 | 0.1% |

## 8. 结果解析

- `Sox2`、`Pax9`、`Trp63` 和 `Tcf4` 的网络度数最高，是候选基因诱导网络中的主要 TF 枢纽。
- `Sox21` 具有中等连接度，可能位于较局部的调控模块。
- `Etv1` 满足严格 TF 判据，但在该诱导网络中的连接度较低；这表示当前 GeneMANIA 证据中与候选基因的直接关联较少，不表示其生物学重要性较低。
- 58.9% 的边来自共表达，23.8% 来自预测证据。因此网络更适合描述功能关联和候选调控模块，不能直接解释为 TF 到靶基因的方向性调控。
- 节点大小只表示网络连接度，不表示表达倍数、差异显著性、TF 活性或调控效应强度。

## 9. 输出文件

- 修改后的 Cytoscape 会话：`top100_network.cys`
- 原始会话备份：`top100_network.pre_go_tf.cys`
- TF 分类表：`output/data/top100_tf_go_classification.csv`
- 单页网络 PDF：`output/pdf/top100_network_GO_TF_network.pdf`
- 三页分析报告：`output/pdf/top100_network_GO_TF_analysis.pdf`
- 高分辨率网络图：`output/images/top100_network_GO_TF_preview.png`

