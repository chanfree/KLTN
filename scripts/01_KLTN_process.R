# ============================================================================
# 脚本名称: 01_KLTN_process.R
# 描述: KLTN样本scRNA-seq数据处理全流程：10X数据读取与QC、DoubletFinder去双胞、
#       聚类与多级细胞注释（免疫/上皮精细注释）、Amplicon-Seq荧光数据整合
#       （计算tdTomato_ratio）、CytoTRACE2干性分析与inferCNV拷贝数变异分析
# 输入: data/10xgenomics/sample_filtered_feature_bc_matrix.h5   10X表达矩阵
#       data/GFP_algin_score.txt                                 GFP扩增子比对结果
#       data/tdTomato_algin_score.txt                            tdTomato扩增子比对结果
#       output/object/seurat/se.rds                              RUN_DOUBLETFINDER=FALSE时读取(已含双胞分类结果)
#       output/object/seurat/se_trajectory_subset.rds            下游脚本(02/03)生成的轨迹子集对象(inferCNV章节依赖)
# 输出: output/object/seurat/se.rds                              完整Seurat对象(含注释与tdTomato_ratio)
#       output/object/seurat/immune_subset_se.rds                免疫子集对象
#       output/object/seurat/immune_annotation.csv               免疫细胞注释表
#       output/object/seurat/epithelium_subset_se.rds            上皮子集对象
#       output/object/seurat/epithelium_annotation.csv           上皮细胞注释表
#       output/object/PREPARE/amplicon_processed_data.csv        扩增子荧光计数汇总表
#       output/object/cytotrace/cytotrace2_result.rds            CytoTRACE2结果
#       output/object/CNV/expFile.txt                            inferCNV表达矩阵输入
#       output/object/CNV_trajectory_HMM_i6/                     inferCNV运行结果目录
#       output/object/td_ratio/KLTN_cnv_scores.csv               各细胞CNV打分表
#       output/figures/下各子目录的QC/聚类/注释/分析PDF图
# 上游脚本: 无（原始10X数据）
# ============================================================================

# ---- 0. 环境与路径配置 ------------------------------------------------------
# 本脚本假定从仓库根目录运行（DIR_ROOT = getwd()）
set.seed(1)

DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "object")
fig_res <- 300

# 脚本特定子目录
DIR_FIG_SEURAT <- file.path(DIR_RES, "figures", "seurat")      # 主流程图表
DIR_OBJ_SEURAT <- file.path(DIR_OBJ, "seurat")                  # 主流程RDS对象
DIR_FIG_PREPARE <- file.path(DIR_RES, "figures", "PREPARE")     # Amplicon-Seq QC图表
DIR_OBJ_PREPARE <- file.path(DIR_OBJ, "PREPARE")                # Amplicon-Seq中间数据
DIR_FIG_cytotrace <- file.path(DIR_RES, "figures", "cytotrace") # CytoTRACE2图表
DIR_OBJ_cytotrace <- file.path(DIR_RES, "object", "cytotrace")  # CytoTRACE2对象
DIR_OBJ_CNV_INPUT <- file.path(DIR_OBJ, "CNV")                  # inferCNV输入文件(expFile.txt)
DIR_OBJ_CNV <- file.path(DIR_OBJ, "CNV_trajectory_HMM_i6")      # inferCNV输出目录
DIR_OBJ_TD_RATIO <- file.path(DIR_OBJ, "td_ratio")              # td_ratio分析输出
DIR_FIG_TD_RATIO <- file.path(DIR_FIG, "td_ratio")              # td_ratio分析图表

dir.create(path = DIR_FIG_SEURAT, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_SEURAT, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_PREPARE, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_PREPARE, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_cytotrace, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_cytotrace, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_CNV_INPUT, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_CNV, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_TD_RATIO, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_TD_RATIO, recursive = TRUE, showWarnings = FALSE)

# 加载R包
# 注: DoubletFinder 经 remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 安装，
#     原脚本曾指定个人库路径 lib.loc，此处改用默认库路径加载
library(Seurat)
library(stringr)
library(clustree)
library(patchwork)
library(dplyr)
library(pheatmap)
library(harmony)
library(future)
library(tidyverse)
library(DoubletFinder)
# SCP包提供 CellDimPlot()（原脚本未显式library，此处补齐显式加载）
library(SCP)
# CytoTRACE2 经 devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 安装
library(CytoTRACE2)
# inferCNV及相关分析包
library(infercnv)
library(AnnoProbe)
library(phylogram)
library(gridExtra)
library(grid)
library(dendextend)
library(ggthemes)
library(miscTools)

# 并行设置
plan(multisession, workers = 8)
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)

# 加载工具函数与配色（plotProportion等自定义函数、cols_*配色方案）
source(file.path(DIR_ROOT, "R", "custom_functions.R"))
source(file.path(DIR_ROOT, "R", "custom_colors.R"))

# ---- 1. 读取10X数据、质控与DoubletFinder去双胞 ------------------------------
# 逻辑开关：DoubletFinder的pK选择与双胞去除流程已运行过，结果已存入se.rds；
# 从10X h5完整重跑时设为TRUE（耗时较长：paramSweep扫描 + doubletFinder分类）
RUN_DOUBLETFINDER <- FALSE

if (RUN_DOUBLETFINDER) {
  ## (1) 读取10X h5并构建Seurat对象，计算线粒体基因比例
  se <- Read10X_h5(file.path(DIR_DATA, "10xgenomics", "sample_filtered_feature_bc_matrix.h5"))
  se <- CreateSeuratObject(se)
  se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^mt-|^MT-|^Mt-")

  ## (2) QC可视化：nFeature/nCount/percent.mt小提琴图与散点图
  Idents(se) <- 'orig.ident'
  se$orig.ident <- 'KLTN'
  p1 <- VlnPlot(se, features = c("nFeature_RNA"), ncol = 1, layer = 'counts', raster = TRUE) +
    theme(axis.title.x = element_blank(), legend.position = "none")
  p2 <- VlnPlot(se, features = c("nCount_RNA"), ncol = 1, layer = 'counts', raster = TRUE) +
    theme(axis.title.x = element_blank(), legend.position = "none")
  p3 <- VlnPlot(se, features = c("percent.mt"), ncol = 1, layer = 'counts', raster = TRUE) +
    theme(axis.title.x = element_blank(), legend.position = "none")
  plot1 <- FeatureScatter(se, feature1 = "nCount_RNA", feature2 = "percent.mt",
                          raster = TRUE, pt.size = 2) +
    theme(legend.text = element_blank(), legend.title = element_blank(), legend.position = "none")
  plot2 <- FeatureScatter(se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          raster = TRUE, pt.size = 2) +
    theme(legend.text = element_blank(), legend.title = element_blank(), legend.position = "none")
  pdf(file = file.path(DIR_FIG_SEURAT, "QC_plot_2.pdf"), width = 4, height = 4)
  print(p1)
  print(p2)
  print(p3)
  print(plot1)
  print(plot2)
  dev.off()

  ## (3) 过滤低质量细胞（nFeature_RNA<6500用于预去除双胞）、标准化与降维
  se <- subset(se, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 15)
  se <- NormalizeData(se)
  se <- FindVariableFeatures(se)
  all.genes <- rownames(se)
  se <- ScaleData(se, features = all.genes)
  se <- RunPCA(se)
  # 交互查看ElbowPlot确定数据维度，此处选择10个PC
  ElbowPlot(se)
  se <- FindNeighbors(se, dims = 1:10)
  se <- RunUMAP(se, dims = 1:10)
  # 注意：此处不要加载future相关冲突包，可能遇到未知错误（原文备注）
  se <- RunTSNE(se, dims = 1:10)
  se <- SCTransform(se)

  ## (4) DoubletFinder pK识别：paramSweep扫描后按BCmetric最大值选择（人工选定pK=0.11）
  sweep.res.list <- paramSweep(se, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn$pK))
  BCmetric <- bcmvn$BCmetric
  pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
  pK_choose <- 0.11
  par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)
  pk_dat <- data.frame(pk = pK, BCmetric = BCmetric)
  pk_plot <- ggplot(data = pk_dat, aes(x = pK, y = BCmetric)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_vline(xintercept = pK_choose, linetype = "dashed", color = "red") +
    annotate("text", x = pK_choose, y = max(BCmetric),
             label = as.character(pK_choose), color = "red", vjust = 1, hjust = -0.5, size = 5) +
    labs(title = "The BCmvn distributions", x = "pK", y = "BCmetric") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  pdf(file = file.path(DIR_FIG_SEURAT, "doubletfinder_pk_determination.pdf"), width = 8, height = 6)
  print(pk_plot)
  dev.off()
  saveRDS(pk_dat, file = file.path(DIR_OBJ_SEURAT, "doubletfinder_pk_dataframe.rds"))

  ## (5) 同源双胞比例估计（假设7.5%双胞形成率）并运行doubletFinder标记双胞
  nExp_poi <- round(0.075 * nrow(se@meta.data))
  se <- doubletFinder(se, PCs = 1:10, pN = 0.25, pK = 0.11,
                      nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
  # 保存含双胞分类结果的对象，供后续快速恢复
  SaveSeuratRds(se, file = file.path(DIR_OBJ_SEURAT, "se.rds"))
} else {
  # 快捷方式：读取已保存的se.rds（含DoubletFinder分类结果，跳过h5处理与paramSweep）
  se <- readRDS(file.path(DIR_OBJ_SEURAT, "se.rds"))
}

# 去除被DoubletFinder标记的双胞细胞
se <- subset(se, subset = DF.classifications_0.25_0.11_844 != 'Doublet')

# ---- 2. 全细胞聚类与聚类分辨率选择 -----------------------------------------
# 多分辨率聚类（0~1）后用clustree评估聚类稳定性，最终选定res=0.2
seq <- seq(0, 1, by = 0.1)
for (res in seq) {
  se[[sprintf("res_%s_clusters", res)]] <- NULL
  se <- FindClusters(se, verbose = FALSE, resolution = res,
                     cluster.name = sprintf("res_%s_clusters", res))
}

# clustree与UMAP用于交互检查分辨率选择（不自动保存PDF）
p1 <- clustree(se, prefix = "res_", suffix = '_clusters') + coord_flip()
p2 <- DimPlot(se, group.by = "res_0.2_clusters", label = TRUE, reduction = "umap")

se$seurat_clusters <- se$res_0.2_clusters
pdf(file = file.path(DIR_FIG_SEURAT, "KLTN_all_se_uamp.pdf"), width = 6, height = 5)
DimPlot(se, group.by = 'seurat_clusters', label = TRUE, raster = TRUE, pt.size = 2)
dev.off()

# 各cluster标记基因（用于初步鉴定细胞类型）
Idents(se) <- se$seurat_clusters
se.markers <- FindAllMarkers(se, only.pos = TRUE)
se.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
dot_plot <- DotPlot(se, features = unique(top5$gene), group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "unannotated_clusters_heatmap.pdf"), width = 20, height = 20)
DoHeatmap(se, features = top5$gene)
dev.off()

# ---- 3. 一级注释（免疫 vs 上皮） --------------------------------------------
# 经典谱系标记基因：Ptprc(免疫)、Epcam(上皮)、Pecam1(内皮)
genes_to_check <- c("Ptprc",   # immune
                    "Epcam",   # epithelium
                    "Pecam1")  # endothelial
se <- SetIdent(se, value = "seurat_clusters")
p <- DotPlot(se, features = genes_to_check)
pdf(file = file.path(DIR_FIG_SEURAT, "first_annotation_markers_dotplot.pdf"), width = 8, height = 6)
p
dev.off()

# 基于res_0.2聚类结果的一级注释（cluster 5/6/7/12为免疫，其余为上皮）
se$first_annotation <- case_when(
  se@meta.data$seurat_clusters == 0 ~ "epithelium",
  se@meta.data$seurat_clusters == 1 ~ "epithelium",
  se@meta.data$seurat_clusters == 2 ~ "epithelium",
  se@meta.data$seurat_clusters == 3 ~ "epithelium",
  se@meta.data$seurat_clusters == 4 ~ "epithelium",
  se@meta.data$seurat_clusters == 5 ~ "immune",
  se@meta.data$seurat_clusters == 6 ~ "immune",
  se@meta.data$seurat_clusters == 7 ~ "immune",
  se@meta.data$seurat_clusters == 8 ~ "epithelium",
  se@meta.data$seurat_clusters == 9 ~ "epithelium",
  se@meta.data$seurat_clusters == 10 ~ "epithelium",
  se@meta.data$seurat_clusters == 11 ~ "epithelium",
  se@meta.data$seurat_clusters == 12 ~ "immune"
)

# ---- 4. 免疫细胞子集精细注释 ------------------------------------------------
imu_se <- subset(se, first_annotation == 'immune')
# 快捷方式（已注释）：直接读取已保存的免疫子集对象（本脚本下方会重新生成该文件）
# imu_se <- readRDS(file.path(DIR_OBJ_SEURAT, "immune_subset_se.rds"))

# 免疫子集重新降维聚类（ElbowPlot交互检查后选10个PC）
imu_se <- RunPCA(imu_se)
ElbowPlot(imu_se)
imu_se <- FindNeighbors(imu_se, dims = 1:10)
imu_se <- RunUMAP(imu_se, dims = 1:10)
imu_se <- RunTSNE(imu_se)

seq <- seq(0, 2, by = 0.2)
for (res in seq) {
  imu_se <- FindClusters(imu_se, verbose = FALSE, resolution = res,
                         cluster.name = sprintf("immune_res_%s_clusters", res))
}

# clustree确定分辨率：最终选定immune_res_1
p1 <- clustree(imu_se, prefix = "immune_res_", suffix = '_clusters') + coord_flip()
p2 <- DimPlot(imu_se, group.by = "immune_res_1_clusters", label = TRUE, reduction = "umap")
p <- p1 + p2 + plot_layout(widths = c(3, 2))
pdf(file = file.path(DIR_FIG_SEURAT, "immnue_clustree_resolution_determination.pdf"),
    width = 16, height = 8)
p
dev.off()
imu_se$seurat_clusters <- imu_se$immune_res_1_clusters

# 各cluster标记基因（曾结合在线注释工具ACT辅助鉴定，JobID: 20250404125150UJR6XDDHKQQ44M）
Idents(imu_se) <- imu_se$seurat_clusters
imu_se.markers <- FindAllMarkers(imu_se, only.pos = TRUE)

imu_se.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50
imu_se.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
dot_plot <- DotPlot(imu_se, features = unique(top5$gene), group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_unannotated_clusters_heatmap.pdf"),
    width = 20, height = 20)
DoHeatmap(imu_se, features = top5$gene)
dev.off()
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_unannotated_clusters_dotplot.pdf"),
    width = 15, height = 6)
print(dot_plot)
dev.off()
write.csv(top50, file = file.path(DIR_OBJ_SEURAT, "immune_subset_unannotated_cluster_top50_markers.csv"),
          row.names = FALSE)

# 免疫细胞谱系经典标记基因
imu_marker <- c(
  "Cd79a", "Cd19", "Ms4a1", 'Ebf1',   # B cell
  "Cd3e", "Cd8a", 'Cd4', "Lef1",      # Naive thymus-derived CD8-positive, alpha-beta T cell
  "Cd3e", 'Cd3d', "Nkg7",             # Mature NK T cell
  'Ncr1', "Gzma", 'Klre1', "Klrb1a",  # Natural Killer
  "Retnlg", "S100a8", "S100a9", 'Il1rn', 'Siglecf',  # Neutrophil
  'Siglech', "H2-DMb1", "Cd68", "Cxcl16",            # Dendritic cell
  "Cd14", "F13a1", "Ly6c2", 'Csf1r',                 # Monocyte
  "Cd3e", "Il17a", 'Trdc',            # Gamma-delta T
  'Tjp1', 'Fermt2',                    # endothelial
  "Marco", "Msr1", "Mrc1", 'Trem2'     # Alveolar Macrophage/TAM
)
imu_marker <- unique(imu_marker)
imu_marker <- factor(imu_marker, levels = imu_marker)
dot_plot <- DotPlot(imu_se, features = imu_marker, group.by = 'seurat_clusters') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")

# 免疫精细注释（按immune_res_1聚类逐一鉴定细胞类型）
imu_se$immnue_annotation <- case_when(
  imu_se@meta.data$seurat_clusters == 0 ~ "B cell",
  imu_se@meta.data$seurat_clusters == 1 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 2 ~ "NKT",
  imu_se@meta.data$seurat_clusters == 3 ~ "NK",
  imu_se@meta.data$seurat_clusters == 4 ~ "B cell",
  imu_se@meta.data$seurat_clusters == 5 ~ "Neutrophil",
  imu_se@meta.data$seurat_clusters == 6 ~ "Dendritic cell",
  imu_se@meta.data$seurat_clusters == 7 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 8 ~ "Neutrophil",
  imu_se@meta.data$seurat_clusters == 9 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 10 ~ "Monocyte",
  imu_se@meta.data$seurat_clusters == 11 ~ "Monocyte"
)
dot_plot <- DotPlot(imu_se, features = imu_marker, group.by = 'immnue_annotation') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_act_auto_annotated_clusters_dotplot.pdf"),
    width = 15, height = 6)
print(dot_plot)
dev.off()
p1 <- DimPlot(imu_se, reduction = 'tsne', group.by = 'immnue_annotation')
p2 <- DimPlot(imu_se, reduction = 'umap', group.by = 'immnue_annotation')
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_act_auto_annotated_clusters_dimplot.pdf"),
    width = 12, height = 5)
p1 + p2
dev.off()

# 粗粒度版本免疫注释（coarse）：标记基因列表去掉Gamma-delta T/内皮/巨噬细胞组
imu_marker_2 <- c(
  "Cd79a", "Cd19", "Ms4a1", 'Ebf1',   # B cell
  "Cd3e", "Lef1",                      # Naive thymus-derived CD8-positive, alpha-beta T cell
  "Cd3e", 'Cd3d', "Nkg7",             # Mature NK T cell
  'Ncr1', "Gzma", 'Klre1', "Klrb1a",  # Natural Killer
  "Retnlg", "S100a8", "S100a9", 'Il1rn', 'Siglecf',  # Neutrophil
  'Siglech', "H2-DMb1", "Cd68", "Cxcl16",            # Dendritic cell
  "Cd14", "F13a1", "Ly6c2", 'Csf1r'                   # Monocyte
)
imu_marker_2 <- unique(imu_marker_2)
imu_marker_2 <- factor(imu_marker_2, levels = imu_marker_2)
imu_se$immnue_annotation_coarse <- case_when(
  imu_se@meta.data$seurat_clusters == 0 ~ "B cell",
  imu_se@meta.data$seurat_clusters == 1 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 2 ~ "NKT",
  imu_se@meta.data$seurat_clusters == 3 ~ "NK",
  imu_se@meta.data$seurat_clusters == 4 ~ "B cell",
  imu_se@meta.data$seurat_clusters == 5 ~ "Neutrophil",
  imu_se@meta.data$seurat_clusters == 6 ~ "Dendritic cell",
  imu_se@meta.data$seurat_clusters == 7 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 8 ~ "Neutrophil",
  imu_se@meta.data$seurat_clusters == 9 ~ "T cell",
  imu_se@meta.data$seurat_clusters == 10 ~ "Monocyte",
  imu_se@meta.data$seurat_clusters == 11 ~ "Monocyte"
)
dot_plot <- DotPlot(imu_se, features = imu_marker_2, group.by = 'immnue_annotation_coarse') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_coarse_annotated_clusters_dotplot.pdf"),
    width = 15, height = 6)
print(dot_plot)
dev.off()
p1 <- DimPlot(imu_se, reduction = 'tsne', group.by = 'immnue_annotation_coarse')
p2 <- DimPlot(imu_se, reduction = 'umap', group.by = 'immnue_annotation_coarse')
pdf(file = file.path(DIR_FIG_SEURAT, "immune_subset_coarse_annotated_clusters_dimplot.pdf"),
    width = 12, height = 5)
p1 + p2
dev.off()

# 保存免疫子集对象与注释表
SaveSeuratRds(imu_se, file = file.path(DIR_OBJ_SEURAT, "immune_subset_se.rds"))
imu_meta <- data.frame(imu_se$immnue_annotation, imu_se$immnue_annotation_coarse)
write.csv(imu_meta, file = file.path(DIR_OBJ_SEURAT, "immune_annotation.csv"),
          row.names = TRUE)

# ---- 5. 上皮细胞子集精细注释 ------------------------------------------------
epi_se <- subset(se, first_annotation != 'immune')

# 上皮子集重新降维聚类（UMAP/TSNE使用1:15 PCs，与原文一致）
epi_se <- RunPCA(epi_se)
ElbowPlot(epi_se)
epi_se <- FindNeighbors(epi_se, dims = 1:10)
epi_se <- RunUMAP(epi_se, dims = 1:15)
epi_se <- RunTSNE(epi_se, dims = 1:15)

seq <- seq(0, 2, by = 0.2)
for (res in seq) {
  epi_se <- FindClusters(epi_se, verbose = FALSE, resolution = res,
                         cluster.name = sprintf("epithelium_res_%s_clusters", res))
}

# clustree确定分辨率：最终选定epithelium_res_0.2
p1 <- clustree(epi_se, prefix = "epithelium_res_", suffix = '_clusters') + coord_flip()
p2 <- DimPlot(epi_se, group.by = "epithelium_res_0.2_clusters", label = TRUE, reduction = "umap")
# fluorescence列由后续Amplicon-Seq章节（或既有se.rds）提供；
# 从10X h5全新重跑(RUN_DOUBLETFINDER=TRUE)时该列尚不存在，此时跳过按荧光分组的UMAP图
p3 <- NULL
if ("fluorescence" %in% colnames(epi_se@meta.data)) {
  p3 <- DimPlot(epi_se, reduction = "umap",
                group.by = "epithelium_res_0.2_clusters", split.by = "fluorescence")
}
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_clustree_resolution_determination.pdf"),
    width = 20, height = 9)
p1 + p2 + plot_layout(widths = c(3, 2))
dev.off()
if (!is.null(p3)) {
  pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_unannotated_fluorescence_group.pdf"),
      width = 16, height = 6)
  p3
  dev.off()
}
epi_se$seurat_clusters <- epi_se$epithelium_res_0.2_clusters

# 各cluster标记基因
Idents(epi_se) <- epi_se$seurat_clusters
epi_se.markers <- FindAllMarkers(epi_se, only.pos = TRUE)
epi_se.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50
epi_se.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
dot_plot <- DotPlot(epi_se, features = unique(top5$gene), group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_subset_unannotated_clusters_heatmap.pdf"),
    width = 20, height = 20)
DoHeatmap(epi_se, features = top5$gene)
dev.off()
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_subset_unannotated_clusters_dotplot.pdf"),
    width = 20, height = 6)
print(dot_plot)
dev.off()
write.csv(top50, file = file.path(DIR_OBJ_SEURAT, "epithelium_subset_unannotated_cluster_top50_markers.csv"),
          row.names = FALSE)

# 上皮各亚型经典标记基因
epi_marker <- list(
  'AT1' = c("Ager", 'Hopx', "Clic5"),                          # AT1 markers
  'AT2' = c("Sftpc", "Sftpa1", "Sftpd", "Etv5"),               # AT2 cells
  'Basal' = c("Krt15", "Krt5", "Krt17", "Trp63"),              # Basal
  'Ciliated' = c("Pifo", "Foxj1", "Ccdc78"),                   # Ciliated
  'Club' = c("Scgb3a2", "Scgb1a1", "Cyp2f2"),                  # Club
  'Goblet' = c('Muc2', "Tff2", "Muc5b", "Spdef"),              # Goblet
  'Endothelium' = c('Pecam1', 'Cdh5', 'Cldn5', ''),            # Endothelial
  'Fibroblast' = c('Col1a2', 'Col1a1', 'Vim', 'Pdgfra'),       # Fibroblast
  'Malignant' = c("Kng2", "Meg3", "Dlk1", "Gkn2", "Lcn2")      # Malignant
)
dot_plot <- DotPlot(epi_se, features = epi_marker, group.by = 'seurat_clusters') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_subset_unannotated_cannonical_markers_dotplot.pdf"),
    width = 15, height = 6)
print(dot_plot)
dev.off()
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_subset_UMAP.pdf"), width = 6, height = 6)
DimPlot(epi_se)
dev.off()

# 上皮精细注释（按epithelium_res_0.2聚类逐一鉴定亚型）
epi_se$epithelium_annotation <- case_when(
  epi_se@meta.data$seurat_clusters == 0 ~ "AT2",
  epi_se@meta.data$seurat_clusters == 1 ~ "AT2",
  epi_se@meta.data$seurat_clusters == 2 ~ "Club/Goblet",
  epi_se@meta.data$seurat_clusters == 3 ~ "Endothelium",
  epi_se@meta.data$seurat_clusters == 4 ~ "AT1",
  epi_se@meta.data$seurat_clusters == 5 ~ "AT2",
  epi_se@meta.data$seurat_clusters == 6 ~ "Club/Goblet",
  epi_se@meta.data$seurat_clusters == 7 ~ "Fibroblast",
  epi_se@meta.data$seurat_clusters == 8 ~ "Ciliated",
  epi_se@meta.data$seurat_clusters == 9 ~ "Endothelium"
)
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_res_0.2_annotated_revised_clusters_dotplot.pdf"),
    width = 12, height = 6)
DimPlot(epi_se, reduction = 'umap', group.by = c('epithelium_annotation', 'seurat_clusters'))
dev.off()
dot_plot <- DotPlot(epi_se, features = epi_marker, group.by = 'epithelium_annotation') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "epithelium_subset_annotated_cannonical_markers_dotplot.pdf"),
    width = 15, height = 6)
print(dot_plot)
dev.off()

# 保存上皮子集对象与注释表
epi_meta <- data.frame(epi_se@meta.data$epithelium_annotation)
rownames(epi_meta) <- colnames(epi_se)
write.csv(epi_meta, file = file.path(DIR_OBJ_SEURAT, "epithelium_annotation.csv"),
          row.names = TRUE)
SaveSeuratRds(epi_se, file = file.path(DIR_OBJ_SEURAT, 'epithelium_subset_se.rds'))

# ---- 6. 最终注释整合 --------------------------------------------------------
# 将免疫(coarse版)与上皮注释合并回全对象，生成final_annotation
immu_meta <- read.csv(file.path(DIR_OBJ_SEURAT, "immune_annotation.csv"),
                      header = TRUE, row.names = 1)
immu_meta$imu_se.immnue_annotation <- NULL
colnames(immu_meta) <- 'annotation'
se <- AddMetaData(se, immu_meta, col.name = 'immune_annotation')
epi_meta <- read.csv(file.path(DIR_OBJ_SEURAT, "epithelium_annotation.csv"),
                     row.names = 1)
se <- AddMetaData(se, metadata = epi_meta, col.name = 'epi_annotation')
se$final_annotation <- ifelse(se$first_annotation == 'immune',
                              se$immune_annotation, se$epi_annotation)
pdf(file = file.path(DIR_FIG_SEURAT, "final_annotation_UMAP.pdf"), width = 6, height = 6)
DimPlot(se, group.by = 'final_annotation')
dev.off()

# 一级注释修订版：非免疫cluster统一为"epithelium/endothelium"
se$first_annotation <- case_when(
  se@meta.data$seurat_clusters == 0 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 1 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 2 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 3 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 4 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 5 ~ "immune",
  se@meta.data$seurat_clusters == 6 ~ "immune",
  se@meta.data$seurat_clusters == 7 ~ "immune",
  se@meta.data$seurat_clusters == 8 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 9 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 10 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 11 ~ "epithelium/endothelium",
  se@meta.data$seurat_clusters == 12 ~ "immune"
)

p1 <- CellDimPlot(se, group.by = "first_annotation", palcolor = c('pink', 'lightblue'),
                  reduction = "UMAP", theme_use = "theme_blank")
pdf(file = file.path(DIR_FIG_SEURAT, "UMAP_first_annotation_res_0.2_clusters.pdf"),
    width = 6, height = 5)
p1
dev.off()

# 全部细胞类型经典标记基因（上皮+免疫）
markers <- list(
  'AT1' = c("Ager", 'Hopx', "Clic5"),                          # AT1 markers
  'AT2' = c("Sftpc", "Sftpa1", "Sftpd", "Etv5"),               # AT2 cells
  'Basal' = c("Krt15", "Krt5", "Krt17", "Trp63"),              # Basal
  'Ciliated' = c("Pifo", "Foxj1", "Ccdc78"),                   # Ciliated
  'Club' = c("Scgb3a2", "Scgb1a1", "Cyp2f2"),                  # Club
  'Goblet' = c('Muc2', "Tff2", "Muc5b", "Spdef"),              # Goblet
  'Malignant' = c("Kng2", "Meg3", "Dlk1", "Gkn2", "Lcn2"),     # Malignant
  'Endothelial' = c('Pecam1', 'Cdh5', 'Cldn5', ''),            # Endothelial
  'Fibroblast' = c('Col1a2', 'Col1a1', 'Vim', 'Pdgfra'),       # Fibroblast
  "B cells" = c("Cd79a", "Cd19", "Ms4a1", 'Ebf1'),
  "T cells" = c("Cd3e", "Lef1"),
  "NKT" = c("Cd3d", "Nkg7"),
  "NK" = c('Ncr1', "Gzma", 'Klre1', "Klrb1a"),
  "Neutrophil" = c("Retnlg", "S100a8", "S100a9", 'Il1rn', 'Siglecf'),
  "Dentritic cell" = c('Siglech', "H2-DMb1", "Cd68", "Cxcl16"),
  "Monocyte" = c("Cd14", "F13a1", "Ly6c2", 'Csf1r')
)
dot_plot <- DotPlot(se, features = markers, group.by = 'res_0.2_clusters') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "unannotaed_final_cannonical_markers_dotplot.pdf"),
    width = 20, height = 6)
print(dot_plot)
dev.off()

#### Fig 1F：发表版最终注释UMAP图
pdf(file = file.path(DIR_FIG_SEURAT, "final_annotation_UMAP_publish.pdf"), width = 6, height = 6)
CellDimPlot(se, group.by = "final_annotation",
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()

# 发表版标记基因列表（去掉Basal与Malignant）并固定final_annotation因子顺序
markers <- list(
  'AT1' = c("Ager", 'Hopx', "Clic5"),                          # AT1 markers
  'AT2' = c("Sftpc", "Sftpa1", "Sftpd", "Etv5"),               # AT2 cells
  'Ciliated' = c("Pifo", "Foxj1", "Ccdc78"),                   # Ciliated
  'Club' = c("Scgb3a2", "Scgb1a1", "Cyp2f2"),                  # Club
  'Goblet' = c('Muc2', "Tff2", "Muc5b", "Spdef"),              # Goblet
  'Endothelial' = c('Pecam1', 'Cdh5', 'Cldn5', ''),            # Endothelial
  'Fibroblast' = c('Col1a2', 'Col1a1', 'Vim', 'Pdgfra'),       # Fibroblast
  "B cells" = c("Cd79a", "Cd19", "Ms4a1", 'Ebf1'),
  "T cells" = c("Cd3e", "Lef1"),
  "NKT" = c("Cd3d", "Nkg7"),
  "NK" = c('Ncr1', "Gzma", 'Klre1', "Klrb1a"),
  "Neutrophil" = c("Retnlg", "S100a8", "S100a9", 'Il1rn', 'Siglecf'),
  "Dentritic cell" = c('Siglech', "H2-DMb1", "Cd68", "Cxcl16"),
  "Monocyte" = c("Cd14", "F13a1", "Ly6c2", 'Csf1r')
)
se$final_annotation <- factor(se$final_annotation,
                              levels = c('AT1', "AT2", "Ciliated", "Club/Goblet",
                                         "Endothelium", "Fibroblast",
                                         "B cell", "T cell", "NKT", "NK", "Neutrophil",
                                         "Monocyte", "Dendritic cell"))
dot_plot <- DotPlot(se, features = markers, group.by = 'final_annotation') +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_SEURAT, "final_cannonical_markers_dotplot.pdf"),
    width = 20, height = 6)
print(dot_plot)
dev.off()

# ---- 7. Amplicon-Seq荧光数据整合（计算tdTomato_ratio） ----------------------
# Amplicon-Seq给出每个细胞中tdTomato与GFP的UMI计数；
# 详细比对流程见项目原始scripts/Amplicon-seq目录
## (1) GFP比对结果QC与barcode提取
GFP <- read.table(file.path(DIR_DATA, "GFP_algin_score.txt"), header = TRUE, sep = '\t')
GFP$Similarity <- str_extract(GFP$Similarity, "(?<=\\().+(?=%\\))") %>% as.numeric()  # 提取括号内数值
GFP$Gaps <- str_extract(GFP$Gaps, ".*(?=/)") %>% as.numeric()
GFP$Score <- as.numeric(GFP$Score)
# 比对质量分布直方图（红线为卡阈值：Similarity>65、Gaps<50、Score>600）
pdf(file = file.path(DIR_FIG_PREPARE, 'hist_QC_GFP.pdf'), width = 10, height = 4)
par(mfrow = c(1, 3))
hist(GFP$Similarity, breaks = 100)
abline(v = 65, col = "red", lwd = 1, lty = 2)
hist(GFP$Gaps, breaks = 100)
abline(v = 50, col = "red", lwd = 1, lty = 2)
hist(GFP$Score, breaks = 100)
abline(v = 600, col = "red", lwd = 1, lty = 2)
dev.off()

# 按阈值过滤后提取barcode（前16bp）与UMI（第17-28bp）
index_GFP <- which(GFP$Gaps < 50 & GFP$Similarity > 65 & GFP$Score > 600)
barcode_umi_GFP <- GFP$bseq[index_GFP]
barcode_GFP <- substr(barcode_umi_GFP, 1, 16)
dat_GFP <- data.frame(barcode = barcode_GFP, umi = substr(barcode_umi_GFP, 17, 28))

## (2) tdTomato比对结果QC与barcode提取（阈值：Gaps<=100、Similarity>=50、Score>=440）
tdTomato <- read.table(file.path(DIR_DATA, "tdTomato_algin_score.txt"), header = TRUE, sep = '\t')
tdTomato$Similarity <- str_extract(tdTomato$Similarity, "(?<=\\().+(?=%\\))") %>% as.numeric()
tdTomato$Gaps <- str_extract(tdTomato$Gaps, ".*(?=/)") %>% as.numeric()
tdTomato$Score <- as.numeric(tdTomato$Score)
pdf(file = file.path(DIR_FIG_PREPARE, 'hist_QC_tdTomato.pdf'), width = 10, height = 4)
par(mfrow = c(1, 3))
hist(tdTomato$Similarity, breaks = 100)
abline(v = 50, col = "red", lwd = 1, lty = 2)
hist(tdTomato$Gaps, breaks = 100)
abline(v = 100, col = "red", lwd = 1, lty = 2)
hist(tdTomato$Score, breaks = 100)
abline(v = 440, col = "red", lwd = 1, lty = 2)
dev.off()

index_tdTomato <- which(tdTomato$Gaps <= 100 & tdTomato$Similarity >= 50 & tdTomato$Score >= 440)
barcode_umi_tdTomato <- tdTomato$bseq[index_tdTomato]
barcode_tdTomato <- substr(barcode_umi_tdTomato, 1, 16)
dat_tdTomato <- data.frame(barcode = barcode_tdTomato, umi = substr(barcode_umi_tdTomato, 17, 28))

## (3) 统计每个barcode的GFP/tdTomato UMI种类数，并与单细胞barcode取交集
data1 <- aggregate(dat_GFP$umi, list(dat_GFP$barcode), function(x) {length(table(x))})
colnames(data1) <- c("BC", "GFP")
data2 <- aggregate(dat_tdTomato$umi, list(dat_tdTomato$barcode), function(x) {length(table(x))})
colnames(data2) <- c("BC", "tdTomato")
data <- merge(data1, data2, by = "BC", all.x = TRUE, all.y = TRUE)

# 交换单细胞对象的barcode前缀（取前16bp比对）
ref <- colnames(se) |> as.data.frame()
colnames(ref) <- 'barcodes'
ref1 <- substr(ref$barcodes, 1, 16)
data <- data[which(data$BC %in% ref1), ]
rownames(data) <- paste0(data$BC, '-1')

## (4) 计算tdTomato_ratio（加0.1伪计数避免除零）
data$GFP <- ifelse(is.na(data$GFP), 0, data$GFP)
data$GFP <- data$GFP + 0.1
data$tdTomato <- ifelse(is.na(data$tdTomato), 0, data$tdTomato)
data$tdTomato <- data$tdTomato + 0.1
data$ratio <- data$tdTomato / data$GFP
data$ratio_scale <- scale(data$ratio)
data[is.na(data)] <- 0
data$sum <- data$GFP + data$tdTomato
data$tdTomato_ratio <- data$tdTomato / data$sum

# 保存扩增子处理结果
write.csv(data, file = file.path(DIR_OBJ_PREPARE, "amplicon_processed_data.csv"))

## (5) 将荧光计数与tdTomato_ratio加入Seurat对象并保存
se <- AddMetaData(se, metadata = data)
se$tdTomato_ratio <- ifelse(is.na(se$tdTomato_ratio), 0, se$tdTomato_ratio)
p <- FeaturePlot(se, features = 'tdTomato_ratio')
pdf(file.path(DIR_FIG_PREPARE, 'umap_td_ratio.pdf'), width = 4, height = 4)
p + scale_color_gradientn(
  colors = viridis::viridis(256),
  limits = c(0, 1),
  oob = scales::squish  # 超出范围的值压缩到边界
)
dev.off()

SaveSeuratRds(se, file = file.path(DIR_OBJ_SEURAT, "se.rds"))

# ---- 8. CytoTRACE2干性分析 --------------------------------------------------
# 读取已保存的完整对象（上一步刚保存），对全部细胞评估分化潜能
se_obj <- readRDS(file.path(DIR_OBJ_SEURAT, "se.rds"))

cytotrace2_result <- cytotrace2(se_obj, is_seurat = TRUE, slot_type = "counts",
                                species = 'mouse', ncores = 32)

annotation <- data.frame(phenotype = se_obj@meta.data$final_annotation) %>%
  set_rownames(., colnames(se_obj))
plots <- plotData(cytotrace2_result = cytotrace2_result,
                  annotation = annotation,
                  is_seurat = TRUE)

# 输出CytoTRACE2各可视化结果
pdf(file = file.path(DIR_FIG_cytotrace, "CytoTRACE2_UMAP.pdf"), width = 8, height = 6)
plots$CytoTRACE2_UMAP
dev.off()
pdf(file = file.path(DIR_FIG_cytotrace, "CytoTRACE2_Relative_UMAP.pdf"), width = 8, height = 6)
plots$CytoTRACE2_Relative_UMAP
dev.off()
pdf(file = file.path(DIR_FIG_cytotrace, "annotated_cell_type_fluorescence_CytoTRACE2_Phenotype_UMAP.pdf"),
    width = 8, height = 6)
plots$Phenotype_UMAP
dev.off()
pdf(file = file.path(DIR_FIG_cytotrace, "annotated_cell_type_fluorescence_CytoTRACE2_Boxplot_byPheno.pdf"),
    width = 8, height = 6)
plots$CytoTRACE2_Boxplot_byPheno
dev.off()
pdf(file = file.path(DIR_FIG_cytotrace, "CytoTRACE2_Potency_UMAP.pdf"), width = 8, height = 6)
plots$CytoTRACE2_Potency_UMAP
dev.off()
saveRDS(cytotrace2_result, file = file.path(DIR_OBJ_cytotrace, "cytotrace2_result.rds"))

# ---- 9. inferCNV拷贝数变异分析 ----------------------------------------------
# 以免疫细胞为正常参考组，对肿瘤细胞进行CNV推断
se_obj <- readRDS(file.path(DIR_OBJ_SEURAT, "se.rds"))
# se_trajectory_subset.rds由下游脚本(02/03)生成，本脚本无法重新生成，保留读取
se_traj <- readRDS(file.path(DIR_OBJ_SEURAT, "se_trajectory_subset.rds"))
# 将轨迹子集对象的最后4列metadata（GSVA模块评分等）合并入全对象
n_cols <- ncol(se_traj@meta.data)
meta_data <- se_traj@meta.data[, (n_cols - 3):n_cols]
se_obj <- AddMetaData(se_obj, metadata = meta_data)
# 区分免疫（参考组）与肿瘤（观察组）
se_obj$traj_immu <- ifelse(se_obj$first_annotation == 'immune', 'immune', 'tumor')
se_obj$traj_immu <- factor(se_obj$traj_immu, levels = c("tumor", "immune"))

# 提取counts矩阵，并用AnnoProbe注释基因的染色体位置
dat <- GetAssayData(se_obj, assay = 'RNA', layer = 'count')
groupinfo <- data.frame(v1 = colnames(dat), v2 = se_obj@meta.data$traj_immu)
geneInfor <- annoGene(rownames(dat), "SYMBOL", "mouse")
# 移除线粒体与性染色体基因
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
# 提取染色体编号并按染色体与起始位置排序
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
geneInfor <- geneInfor[with(geneInfor, order(chr_num, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
# 仅保留有染色体注释的基因，并按geneInfor顺序重排表达矩阵
dat <- dat[rownames(dat) %in% geneInfor[, 1], ]
dat <- dat[match(geneInfor[, 1], rownames(dat)), ]

# 写出inferCNV三个输入文件（表达矩阵、细胞分组、基因排序）
expFile <- file.path(DIR_OBJ_CNV_INPUT, 'expFile.txt')
colnames(dat) <- gsub("-", "_", colnames(dat))  # barcode中的"-"替换为"_"，与分组文件保持一致
write.table(as.matrix(dat), file = expFile, sep = '\t', quote = FALSE)
groupFiles <- file.path(DIR_OBJ_CNV, 'groupFiles.txt')
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(groupinfo, file = groupFiles, sep = '\t',
            quote = FALSE, col.names = FALSE, row.names = FALSE)
geneFile <- file.path(DIR_OBJ_CNV, 'geneFile.txt')
write.table(geneInfor, file = geneFile, sep = '\t',
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# 创建inferCNV对象（immune为正常参考组）
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c('immune')
)
# 使用analysis_mode="subclusters"时建议设置scipen，避免hclust生成阶段报错
options(scipen = 100)

# 运行inferCNV（HMM i6模型，subclusters模式；cluster_by_groups=FALSE以便读取进化树）
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # 10X数据选择0.1（smart-seq选择1）
                              out_dir = DIR_OBJ_CNV,
                              cluster_by_groups = FALSE,
                              plot_steps = FALSE,
                              scale_data = TRUE,
                              denoise = TRUE,
                              noise_filter = 0.12,
                              analysis_mode = 'subclusters',
                              num_threads = 20,
                              HMM = TRUE,
                              HMM_type = 'i6',
                              write_phylo = TRUE  # 生成newick进化树
)
# 中值滤波降噪并绘制CNV热图
infercnv_obj_medianfiltered <- infercnv::apply_median_filtering(infercnv_obj)
infercnv::plot_cnv(infercnv_obj_medianfiltered,
                   out_dir = DIR_OBJ_CNV,
                   output_filename = 'infercnv.median_filtered',
                   x.range = "auto",
                   x.center = 1,
                   title = "infercnv",
                   color_safe_pal = FALSE)

# 快捷方式（已注释）：读取已完成的inferCNV最终对象，跳过上文运行步骤
# infercnv_obj <- readRDS(file.path(DIR_OBJ_CNV, "run.final.infercnv_obj"))

# ---- 10. CNV打分与可视化 ----------------------------------------------------
# 基于参考组（immune）表达分布确定拷贝数增益/丢失的阈值
expr <- infercnv_obj@expr.data
reference_groups <- c("immune")
tmp_list <- lapply(reference_groups, function(group) {
  expr[, infercnv_obj@reference_grouped_cell_indices[[group]]]
})
tmp <- do.call(cbind, tmp_list)
# 计算参考组均值±2SD得到中性区间，进一步推算拷贝数区间边界
down <- mean(rowMeans(tmp)) - 2 * mean(apply(tmp, 1, sd))
up <- mean(rowMeans(tmp)) + 2 * mean(apply(tmp, 1, sd))
oneCopy <- up - down
a1 <- down - 2 * oneCopy
a2 <- down - 1 * oneCopy
a3 <- up + 1 * oneCopy
a4 <- up + 2 * oneCopy

# 按表达值区间对每个基因-细胞打分类别（A完全丢失~F多拷贝增益）
cnv_score_table <- infercnv_obj@expr.data
cnv_score_mat <- as.matrix(cnv_score_table)
cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A"          # 完全丢失，2分
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B"      # 丢失一个拷贝，1分
cnv_score_table[cnv_score_mat >= down & cnv_score_mat < up] <- "C"      # 中性，0分
cnv_score_table[cnv_score_mat >= up & cnv_score_mat <= a3] <- "D"       # 增加一个拷贝，1分
cnv_score_table[cnv_score_mat > a3 & cnv_score_mat <= a4] <- "E"        # 增加两个拷贝，2分
cnv_score_table[cnv_score_mat > a4] <- "F"                              # 增加两个以上拷贝，2分

# 类别替换为CNV分数并求和得到每个细胞的总CNV得分
cnv_score_table_pts <- cnv_score_mat
rm(cnv_score_mat)
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2

cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
score <- cell_scores_CNV
rownames(score) <- gsub('_', '-', rownames(score))  # barcode分隔符还原为"-"

# 将CNV总分写入metadata并保存
meta <- se_obj@meta.data
meta$totalCNV <- score[match(colnames(se_obj), rownames(score)), 1]
write.csv(meta, file = file.path(DIR_OBJ_TD_RATIO, "KLTN_cnv_scores.csv"))

# 按聚类分组的CNV得分箱线图
meta$seurat_clusters <- factor(meta$seurat_clusters)
p <- ggplot(meta, aes(x = seurat_clusters, y = totalCNV, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "Seurat Clusters", y = "Total CNV score") +
  theme_classic() +
  theme(legend.position = "right")
pdf(file = file.path(DIR_FIG_TD_RATIO, "KLTN_cnv_scores_boxplot.pdf"), width = 6, height = 6)
print(p)
dev.off()
