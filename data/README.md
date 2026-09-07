# data/ 数据说明与下载清单

本目录包含复现全部分析所需的数据文件（不含图表 PDF）。
`MD5SUMS.txt` 为全部文件的 md5 校验和。

## 目录结构

```
data/
├── 10xgenomics/          # 本实验 KLTN 原始 10X 矩阵（脚本 01）
├── GFP_algin_score.txt   # GFP 扩增子比对结果（脚本 01, Amplicon-Seq）※大文件
├── tdTomato_algin_score.txt  # tdTomato 扩增子比对结果 ※大文件
├── geneset.hs.gmt        # 人源基因集（脚本 02/03）
├── scc_tong_2026/        # SCC 数据集（脚本 02）※scc_epi_subset_obj.rds 为大文件
├── gene_orthology/       # 人源/鼠源 signature 列表（脚本 05）
├── clinical_validation/  # KRAS/EGFR 临床验证数据与类器官 TPM（脚本 05）※kras_se_obj / egfr_epi_subset 为大文件
└── precomputed/          # 上游脚本产出的中间对象（镜像 output/object/ 结构）
```

## 使用方法

- `data/` 根下的直接输入文件：脚本按相对路径直接读取，无需移动。
- `data/precomputed/` 下的中间对象：**复制到仓库根的 `output/object/` 对应相对路径**后，
  各脚本（开关为 FALSE 时）即可从中间结果继续，无需从头运行上游步骤。例如：

  ```bash
  # 将全部中间对象一次性就位（在仓库根目录执行）
  cp -r data/precomputed/* output/object/
  ```

## 文件清单

标记说明：✅ = 可直接提交 GitHub（<100MB）；☁ = 超过 GitHub 100MB 单文件上限，
已列入 `.gitignore`，需外部托管（见下文"大文件托管"）。

### 1. 直接输入数据（脚本从 data/ 读取）

| 文件 | 大小 | 标记 | 用途（脚本） | 来源 |
|------|------|------|--------------|------|
| `10xgenomics/sample_filtered_feature_bc_matrix.h5` | 62M | ✅ | KLTN 原始矩阵（01） | 本实验 10X 测序 |
| `GFP_algin_score.txt` | 1.2G | ☁ | GFP 扩增子比对（01 第7章） | 本实验 Amplicon-Seq |
| `tdTomato_algin_score.txt` | 1.6G | ☁ | tdTomato 扩增子比对（01 第7章） | 本实验 Amplicon-Seq |
| `geneset.hs.gmt` | 15K | ✅ | GMT 基因集（02 第9/10章） | 项目整理 |
| `scc_tong_2026/scc_epi_subset_obj.rds` | 1.2G | ☁ | SCC 上皮子集（02 第1章） | 合作数据 |
| `scc_tong_2026/filtered_feature_bc_matrix.h5` | 82M | ✅ | SCC 原始矩阵（02 开关TRUE） | 合作数据 |
| `gene_orthology/human_all_signature_without_pattern_list.rds` | 2.8K | ✅ | signature 列表（05 第2章） | 03/04 产出转人源 |
| `clinical_validation/KRAS_induced& suppressed.csv` | 11K | ✅ | KRAS 基因集（05 第3章） | MSigDB 整理 |
| `clinical_validation/KOBAYASHI_EGFR_SIGNALING_24HR_UP.v2026.1.Hs.gmt` | 719B | ✅ | EGFR 基因集（05 第3章） | MSigDB C3 |
| `clinical_validation/KOBAYASHI_EGFR_SIGNALING_24HR_DN.v2026.1.Hs.gmt` | 1.6K | ✅ | EGFR 基因集（05 第3章） | MSigDB C3 |
| `clinical_validation/kras_se_obj_with_scores.rds` | 336M | ☁ | KRAS 单细胞对象（05 第3章） | GSE137912 处理版 |
| `clinical_validation/egfr_epi_subset.rds` | 1.3G | ☁ | EGFR 单细胞对象（05 第3.6章） | scell_lung_adenocarcinoma 处理版（原始数据 PRJNA591860） |
| `clinical_validation/EGFR_Seu_metadata.csv` | 12M | ✅ | EGFR 单细胞元数据（05 第3.6章） | scell_lung_adenocarcinoma（原始数据 PRJNA591860） |

### 2. 中间对象（data/precomputed/，镜像 output/object/ 相对路径）

| 文件 | 大小 | 标记 | 产出脚本 | 消费脚本 |
|------|------|------|----------|----------|
| `precomputed/seurat/se.rds` | 2.1G | ☁ | 01 | 01/02/03 |
| `precomputed/seurat/se_trajectory_subset.rds` | 1.4G | ☁ | 03 | 01/02/03/04 |
| `precomputed/seurat/se_trajectory_subset.h5ad` | 289M | ☁ | 03（ExportToSeurat/SCE） | cellrank.ipynb |
| `precomputed/seurat/se_harmony_integrated_epi_subset_final_3_samples.rds` | 4.4G | ☁ | 02 第5章 | 02 第8/9章 |
| `precomputed/seurat/se_harmony_integrated_3samples_epithelium.rds` | 2.6G | ☁ | 02 | 04（需含 dnb_state 列） |
| `precomputed/seurat/immune_annotation.csv` | 58K | ✅ | 01 第4章 | 01 第6章 |
| `precomputed/seurat/epithelium_annotation.csv` | 263K | ✅ | 01 第5章 | 01 第6章 |
| `precomputed/slingshot/sce_full.rds` | 75M | ✅ | 03 第3章 | 03 第4章 |
| `precomputed/slingshot/slingshot_pseudotime.csv` | 322K | ✅ | 03 第3章 | cellrank.ipynb |
| `precomputed/slingshot/slingshot_weights.csv` | 198K | ✅ | 03 第3章 | cellrank.ipynb |
| `precomputed/slingshot/patternRes_filter.csv` | 667K | ✅ | 03 第4章 | 03 第4章 |
| `precomputed/cellrank/C1_drivers.csv` | 2.7M | ✅ | cellrank.ipynb | 03 第5章 |
| `precomputed/cellrank/C2_drivers.csv` | 2.7M | ✅ | cellrank.ipynb | 03 第5章 |
| `precomputed/td_ratio/association_results.csv` | 201K | ✅ | 03 第6章 | 02 第8章 |
| `precomputed/survival/TCGA-LUAD_vst_counts_with_symbols.rds` | 68M | ✅ | 03 第5章 | 03/05 |
| `precomputed/survival/TCGA-LUAD_clinical_data.rds` | 194K | ✅ | 03 第5章 | 03/05 |
| `precomputed/survival/TCGA-LUSC_vst_counts_with_symbols.rds` | 67M | ✅ | 03 第5章 | 05 第2章 |
| `precomputed/survival/TCGA-LUSC_clinical_data.rds` | 179K | ✅ | 03 第5章 | 05 第2章 |
| `precomputed/survival/human_all_signature_without_pattern_list.rds` | 2.3K | ✅ | 05 | 05 第2/3章 |
| `precomputed/survival/KRYSTAL-1_trial/vst_counts.rds` | 4.1M | ✅ | 05（开关生成） | 05 第1/2/3章 |
| `precomputed/survival/KRYSTAL-1_trial/meta_dat.csv` | 12K | ✅ | 05（开关生成） | 05 第1/2/3章 |
| `precomputed/survival/KRYSTAL-1_trial/mmc2 (1).xlsx` | 26K | ✅ | — | 05 预处理开关 |
| `precomputed/survival/KRYSTAL-1_trial/mmc3 (1).xlsx` | 7.2M | ✅ | — | 05 预处理开关 |
| `precomputed/SCENIC_trajectory/int/scenicOptions.Rds` | 1.4K | ✅ | 03 第7章 | 03 第7章 |
| `precomputed/SCENIC_trajectory/int/3.5_AUCellThresholds.Rds` | 1.8M | ✅ | 03 第7章 | 03 第7章 |
| `precomputed/SCENIC_trajectory/int/2.6_regulons_asGeneSet.Rds` | 1.0M | ✅ | 03 第7章 | 04 第12/13章 |
| `precomputed/SCENIC_td/aux_createViewerApp.R` | 12K | ✅ | SCENIC 辅助脚本 | 03 source |

## 大文件托管（☁ 标记，共 10 个，约 15.6GB）

以下文件超过 GitHub 单文件 100MB 上限，已列入 `.gitignore`，不会随仓库提交。
建议整体打包上传至 **Zenodo**（单文件最高 50GB，免费 DOI），并将链接填入下表：

| 文件 | 大小 | 下载链接 |
|------|------|----------|
| `GFP_algin_score.txt` | 1.2G | 【待上传后填写】 |
| `tdTomato_algin_score.txt` | 1.6G | 【待上传后填写】 |
| `scc_tong_2026/scc_epi_subset_obj.rds` | 1.2G | 【待上传后填写】 |
| `clinical_validation/kras_se_obj_with_scores.rds` | 336M | 【待上传后填写】 |
| `clinical_validation/egfr_epi_subset.rds` | 1.3G | 【待上传后填写】 |
| `precomputed/seurat/se.rds` | 2.1G | 【待上传后填写】 |
| `precomputed/seurat/se_trajectory_subset.rds` | 1.4G | 【待上传后填写】 |
| `precomputed/seurat/se_trajectory_subset.h5ad` | 289M | 【待上传后填写】 |
| `precomputed/seurat/se_harmony_integrated_epi_subset_final_3_samples.rds` | 4.4G | 【待上传后填写】 |
| `precomputed/seurat/se_harmony_integrated_3samples_epithelium.rds` | 2.6G | 【待上传后填写】 |

打包命令（在仓库根目录执行，产物可直接上传 Zenodo）：

```bash
cd data
tar -czf KL_TRP63_large_files.tar.gz \
  GFP_algin_score.txt tdTomato_algin_score.txt \
  scc_tong_2026/scc_epi_subset_obj.rds \
  clinical_validation/kras_se_obj_with_scores.rds \
  clinical_validation/egfr_epi_subset.rds \
  precomputed/seurat/se.rds precomputed/seurat/se_trajectory_subset.rds \
  precomputed/seurat/se_trajectory_subset.h5ad \
  precomputed/seurat/se_harmony_integrated_epi_subset_final_3_samples.rds \
  precomputed/seurat/se_harmony_integrated_3samples_epithelium.rds
md5sum KL_TRP63_large_files.tar.gz > KL_TRP63_large_files.tar.gz.md5
```

## 外部公共数据（未打包，按官方链接下载）

| 数据 | 用途 | 官方下载 |
|------|------|----------|
| KPTracer `adata_processed.nt.rds` (616M) | 02 第1章 LUAD 参考 | https://github.com/morris-lab/KPTracer-Data (`expression/adata_processed.nt.rds`)，下载后修改脚本 02 顶部 `PATH_KP_LUAD_RDS` |
| EGFR 单细胞数据 (scell_lung_adenocarcinoma) | 05 第3.6章 EGFR 通路验证（FigS5 D-F） | https://github.com/czbiohub-sf/scell_lung_adenocarcinoma （Data_input 文件夹内 Google Drive 链接）；原始测序 PRJNA591860，转化脚本见 data/clinical_validation/EGFR_seu_convertion.R |
| SCENIC mm10 ranking 数据库（feather） | 03 第7章 SCENIC | https://resources.aertslab.org/cistarget/ （mm10 `10bp_up_10bp_down` 与 `500bp_up_100bp_down` 两套 clustering 版），置于脚本 03 配置的 `SCENIC_ref` 路径 |
| TCGA-LUAD/LUSC 原始数据 | 03/05 生存分析 | 无需下载——本目录已提供处理后的 vst+clinical RDS（`precomputed/survival/`）；如需原始数据见 https://portal.gdc.cancer.gov |
| KRYSTAL-1 mmc2/mmc3 | 05 预处理开关 | 无需下载——已包含（`precomputed/survival/KRYSTAL-1_trial/`） |

## 校验

下载后请核对完整性：

```bash
cd data && md5sum -c MD5SUMS.txt
```
