# KL-TRP63 肺腺癌向鳞癌转分化（AST）单细胞分析流程

基于 R/Seurat 的肺类器官单细胞 RNA 测序（scRNA-seq）分析流程，研究 TRP63 阳性细胞在肺腺癌（LUAD）向肺鳞癌（LUSC）转分化（Adeno-to-Squamous Transition, AST）过程中的动态规律，并运用 Dynamic Network Biomarker（DNB）理论识别表型转换临界点（Tipping Point），评估 NRG1 与 trajectory_up signature 的临床价值。

## 核心概念

| 概念 | 描述 |
|------|------|
| `tdTomato_ratio` | TRP63 表达比例（0 = LUAD 表型，1 = LUSC 表型），由 Amplicon-Seq 荧光 barcode 计算而来 |
| `GFP` / `tdTomato` | TRP63 低 / 高表达细胞群 |
| `AST` | 腺鳞转化（Adeno-to-Squamous Transition）轨迹亚群 |
| `DNB` | 动态网络生物标志物，CI = (PCC_in / PCC_out) × SD_in |
| `Tipping Point` | 表型转换临界状态；其前后分别为 Pre-Tipping / Post-Tipping |

## 目录结构

```
repository/
├── README.md                        # 本文件
├── R/                               # 工具函数与配色
│   ├── custom_functions.R           # 通用绘图与辅助函数
│   └── custom_colors.R              # 统一颜色配置
├── scripts/                         # 分析脚本（按编号顺序运行）
│   ├── 01_KLTN_process.R            # 本实验 KLTN 数据处理
│   ├── 02_integrated_analysis.R     # 与其他研究单细胞数据整合
│   ├── 03_trajectory_subset_analysis.R  # AST 轨迹亚群分析（含 FigS4 的 top100 关联基因聚类/GO/网络）
│   ├── 04_DNB_analysis.R            # DNB 临界点分析
│   ├── 05_survival_analysis.R       # 临床价值分析（含 FigS5 的 EGFR/KRAS signature 箱线图）
│   └── cellrank.ipynb               # CellRank 轨迹潜能分析（Python, 配合脚本 03）
├── top100_associate_cytoscape/       # FigS4 B 的 top100 关联基因 GeneMANIA/Cytoscape 网络工作流
├── data/                            # 输入数据（见下文"输入数据"）
└── output/                          # 图表与 RDS 输出
    ├── figures/                     # 图表（按分析模块分子目录）
    ├── object/                      # RDS 对象
    └── tables/                      # 分析表格
```

## 分析流程与运行顺序

脚本假定**从仓库根目录**运行（`DIR_ROOT <- getwd()`），请勿在工作目录不一致时 source。

| 顺序 | 脚本 | 内容 | 关键产出 |
|------|------|------|----------|
| 1 | `01_KLTN_process.R` | 读取 10X h5 → QC → DoubletFinder 去双胞 → 聚类 → 免疫/上皮精细注释 → Amplicon-Seq 荧光整合（`tdTomato_ratio`）→ CytoTRACE2 → inferCNV | `se.rds`、各子集 RDS、CNV 打分 |
| 2 | `02_integrated_analysis.R` | KLTN + SCC + KP-LUAD 三数据集统一 QC → SCT + Harmony 整合 → 参考注释映射（symphony K-NN）→ AST 过渡态验证 → GMT 通路评分 | 整合 RDS、cell_annotation 图表 |
| 3 | `03_trajectory_subset_analysis.R` | AST 轨迹亚群：tdRatio 时序基因 → slingshot 轨迹 → tradeSeq 伪时间 → SCENIC 调控网络 → signature 评分 → top100 关联基因聚类/GO/网络(FigS4) | `se_trajectory_subset.rds`、SCENIC 结果、top100 热图/GO 图 |
| 3.5 | `cellrank.ipynb`（Python） | PseudotimeKernel + ConnectivityKernel 组合核 → GPCCA 宏状态/命运概率 → C1/C2 谱系驱动基因 | `C1/C2_drivers.csv`（供脚本 03 第 5 章使用） |
| 4 | `04_DNB_analysis.R` | DNB 滑动窗口 → Tipping Point 识别 → Pre/Tipping/Post 状态 DEG → 动态 TF-靶基因网络 | `dnb_state` 列、DNB_final_df.csv |
| 5 | `05_survival_analysis.R` | NRG1 与 trajectory_up signature 的 KM 生存分析（TCGA-LUAD/LUSC、KRYSTAL-1）→ EGFR/KRAS 通路相关性 → EGFR/KRAS signature 箱线图(FigS5) | 生存曲线、Cox 结果、signature 箱线图 |

### 脚本间数据流（generated_data 说明）

本项目两个核心中间数据对象（原 `src/generated_data.R` 记录的路径）：

- `output/object/seurat/se.rds` — 脚本 01 产出（完整 KLTN 注释对象，含荧光与 CNV 信息）→ 供脚本 02、03 使用
- `output/object/seurat/se_trajectory_subset.rds` — 脚本 03 产出（AST 轨迹亚群）→ 供脚本 04 使用，并会随后续分析（脚本 02 第 10 章、脚本 04）追加评分列

注意运行顺序：脚本 02 的部分章节依赖脚本 03 的产出（`association_results.csv`、`se_trajectory_subset.rds`），建议完整顺序 01 → 03 → 02 → 04 → 05，或参照各脚本头部"输入"说明。

另注意：`cellrank.ipynb`（Python）需在脚本 03 的 slingshot 章节之后运行（读取 `se_trajectory_subset.h5ad` 与 slingshot 伪时间 CSV），其产出的 `C1/C2_drivers.csv` 是脚本 03 第 5 章（TCGA 生存分析）的输入。

## 输入数据

**完整数据清单、md5 校验与大文件下载链接见 [data/README.md](data/README.md)。**

```
data/
├── README.md                        # 数据清单与下载说明（38个文件的用途/大小/校验）
├── MD5SUMS.txt                      # 全部数据文件的 md5 校验和
├── 10xgenomics/                     # 本实验 KLTN 原始 10X 矩阵
├── GFP_algin_score.txt              # GFP 扩增子比对结果（>100MB, 外部托管）
├── tdTomato_algin_score.txt         # tdTomato 扩增子比对结果（>100MB, 外部托管）
├── geneset.hs.gmt                   # 人源基因集（GMT 格式）
├── scc_tong_2026/                   # SCC 数据集
├── gene_orthology/                  # 人源/鼠源 signature 列表（脚本 05 第2/5章）
├── clinical_validation/             # KRAS/EGFR 临床验证数据 + 小鼠类器官 TPM（脚本 05 第3/5章）
└── precomputed/                     # 上游产出的中间对象（镜像 output/object/ 结构,
                                      #   复制到 output/object/ 后脚本可从中间结果继续）
```

数据获取要点：

- **小文件（<100MB，29 个）**：随仓库直接分发
- **大文件（>100MB，9 个，约 14.3GB）**：超过 GitHub 单文件上限，已列入 `.gitignore`，
  需从外部托管下载（Zenodo 链接见 data/README.md），下载后放回 data/ 对应路径
- **中间对象使用**：`cp -r data/precomputed/* output/object/` 后，各脚本（开关 FALSE）
  即可从中间结果继续，无需从头运行上游步骤

外部公共数据（未随仓库分发，官方下载链接见 data/README.md）：

- KPTracer：`PATH_KP_LUAD_RDS`（脚本 02 顶部声明）— KPTracer-Data 的 `adata_processed.nt.rds`
- SCENIC mm10 ranking 数据库（feather）：脚本 03 配置的 `SCENIC_ref` 路径
- SCENIC 辅助 `aux_createViewerApp.R` 与 KRYSTAL-1 `mmc2/mmc3.xlsx`：已包含于 `data/precomputed/`

## 逻辑开关

重跑耗时段落时将对应开关设为 `TRUE`（默认 `FALSE` 时读取已保存结果）：

| 脚本 | 开关 | 作用 |
|------|------|------|
| 01 | `RUN_DOUBLETFINDER` | 从 h5 完整重跑 QC + DoubletFinder |
| 02 | `RUN_SCC_PREPROCESS` | 从 SCC 原始 h5 完整预处理 |
| 03 | `RUN_TCGA_DOWNLOAD` | 重新下载 TCGA 数据 |
| 05 | `RUN_KRYSTAL_PREPROCESS` / `RUN_EXPORT_SIGNATURE_LONG_CSV` | KRYSTAL-1 预处理 / 导出 signature |

## 依赖包

核心：Seurat, harmony, tidyverse(dplyr/ggplot2/tidyr/stringr/purrr), patchwork, clustree, DoubletFinder, future
轨迹：slingshot, tradeSeq, CytoTRACE2
调控网络：SCENIC(GENIE3/AUCell/RcisTarget), igraph
通路：GSVA/ssGSEA, clusterProfiler
临床：survival, survminer, timeROC, rstatix
其他：pheatmap, ggalluvial, reshape2, SCP, symphony, presto, scico, viridis, ggsci, RColorBrewer, org.Mm.eg.db

> 版本提示：原环境 R 4.4；igraph 与 Seurat 存在版本冲突时，先加载 igraph（脚本 02 已处理）。部分包（DoubletFinder、SCP、CytoTRACE2）需从 GitHub 安装。

Python（`scripts/cellrank.ipynb`）：cellrank, scanpy, scvelo, anndata, pandas, numpy；`cr.models.GAMR` 需 rpy2 并调用 R 的 mgcv 包。

## 维护说明

- 整理过程与各脚本详细变更记录见 `../整理进度索引.md`
- 原始未整理脚本保留于 `../src/`
- 各脚本头部含"输入/输出/上游脚本"声明，修改时请同步更新
