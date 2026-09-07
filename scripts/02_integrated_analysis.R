# ============================================================================
# 脚本名称: 02_integrated_analysis.R
# 描述: 将本项目KLTN(AST)轨迹数据集与SCC及KPTracer-LUAD单细胞数据整合
#       (SCT标准化→Harmony批次校正→聚类注释→K-NN参考注释映射)，
#       并验证AST细胞处于LUAD与SCC之间的过渡态(签名得分/统计检验/通路评分)
# 输入: output/object/seurat/se.rds                                  01_KLTN_process.R产出(本项目数据, 本章作为trajectory数据集)
#       data/scc_tong_2026/scc_epi_subset_obj.rds                    SCC数据集上皮子集(RUN_SCC_PREPROCESS=FALSE时读取)
#       data/scc_tong_2026/filtered_feature_bc_matrix.h5             SCC原始10X矩阵(RUN_SCC_PREPROCESS=TRUE时从头预处理)
#       PATH_KP_LUAD_RDS(外部KPTracer数据库adata_processed.nt.rds)   KPTracer-LUAD参考数据集
#       data/geneset.hs.gmt                                          人类基因集GMT文件(第9/10章module评分)
#       output/object/td_ratio/association_results.csv               03轨迹tradeSeq关联结果(第8章轨迹上调signature)
#       output/object/seurat/se_trajectory_subset.rds                03产出的轨迹子集对象(第10-12章输入; 第10章会覆盖保存)
#       output/object/seurat/se_harmony_integrated_epi_subset_final_3_samples.rds  本章第5-6章产出(重跑第8章时读取)
# 输出: output/object/seurat/se_harmony_integrated_all_final_3_samples.rds          全细胞整合对象
#       output/object/seurat/se_harmony_integrated_epi_subset_final_3_samples.rds   上皮子集整合对象(含symphony注释)
#       output/object/seurat/se_harmony_integrated_epi_subset_with_symphony.rds     含K-NN注释的上皮子集对象
#       output/object/seurat/query_symphony_annotated.rds                           Query数据集注释对象
#       output/object/se_harmony_integrated_epi_subset_final_3_samples.rds          第9章module评分版(注意: 位于output/object/而非seurat/子目录, 原脚本如此)
#       output/object/se_harmony_integrated_epi_subset_final_3_samples_analysis_meta.csv  整合对象metadata导出
#       output/object/seurat/transition_score_stats.csv                             三组LUAD/SCC签名得分统计
#       output/object/seurat/se_trajectory_subset.rds                               第10章覆盖保存(追加Module_*评分列)
#       output/object/gsva/deg_td_vs_gfp_filter.csv                                 tdTomato vs GFP差异基因
#       output/cell_annotation/*.pdf, *.csv                                         K-NN注释结果图与统计表
#       output/figures/integration/raw_mixed/*.pdf                                  QC/PCA elbow/一级注释点图
#       output/figures/integration/transition_validation/*.pdf                      过渡态验证全部图表
#       output/figures/integration/*canonical_signature*.pdf                        通路signature箱线图
#       output/figures/gsva/*.pdf                                                   轨迹子集td组图表
#       output/figures/signature_score/*.pdf                                        通路评分箱线图/Nrg1相关性图
# 上游脚本: 01_KLTN_process.R (se.rds等); 03_trajectory_subset_analysis.R (se_trajectory_subset.rds, association_results.csv)
# ============================================================================

# ---- 0. 环境与路径配置 --------------------------------------------------------
# 本脚本假定从仓库根目录运行（DIR_ROOT = getwd()）
set.seed(123) # 固定随机种子保证可重复性

DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "object")
fig_res <- 300

# 脚本特定子目录
DIR_FIG_RAW_MIXED <- file.path(DIR_RES, "figures", "integration", "raw_mixed")  # 第1-4章: 整合前QC与全细胞注释图
DIR_FIG_INTEGRATION <- file.path(DIR_RES, "figures", "integration")             # 第8-9章: 整合对象signature图
FIG_TRANS <- file.path(DIR_FIG_INTEGRATION, "transition_validation")            # 第8-9章: 过渡态验证图
DIR_FIG_CELL_ANN <- file.path(DIR_RES, "cell_annotation")                       # 第6-7章: K-NN注释图
DIR_OBJ_SEURAT <- file.path(DIR_OBJ, "seurat")                                  # Seurat对象输出
DIR_FIG_SIG <- file.path(DIR_FIG, "signature_score")                            # 第10-11章: 通路评分图
DIR_OBJ_SIG <- file.path(DIR_OBJ, "signature_score")                            # 第10-11章: 评分对象
DIR_FIG_GSVA <- file.path(DIR_FIG, "gsva")                                      # 第10/12章: 轨迹子集图
DIR_OBJ_GSVA <- file.path(DIR_OBJ, "gsva")                                      # 第10章: DEG表

dir.create(path = DIR_FIG_RAW_MIXED, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_INTEGRATION, recursive = TRUE, showWarnings = FALSE)
dir.create(path = FIG_TRANS, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_CELL_ANN, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_SEURAT, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_SIG, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_SIG, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_FIG_GSVA, recursive = TRUE, showWarnings = FALSE)
dir.create(path = DIR_OBJ_GSVA, recursive = TRUE, showWarnings = FALSE)

# 外部公共数据集: KPTracer LUAD单细胞数据(https://github.com/morris-lab/KPTracer-Data,
# expression/adata_processed.nt.rds), 复现时请下载至本地并修改此路径
PATH_KP_LUAD_RDS <- "/home/data/t230405/database/KPTracer-Data/expression/adata_processed.nt.rds"

# 逻辑开关: SCC数据集从原始10X h5完整预处理(QC→标准化→聚类→上皮注释→取cluster 2/4)
# 已运行过, 结果已保存为 data/scc_tong_2026/scc_epi_subset_obj.rds, 重跑时设为TRUE
RUN_SCC_PREPROCESS <- FALSE

# 加载R包
# 注: igraph需在Seurat之前加载(版本依赖, 原脚本指定个人库路径lib.loc, 此处改用默认库路径)
library(igraph)
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(future)
library(clustree)
library(patchwork)
library(presto)
library(SCP)        # CellDimPlot/CellDensityPlot/FeatureDimPlot
library(tibble)
library(FNN)        # K-NN搜索
library(orthogene)  # 人源基因集→小鼠同源基因转换
library(ggpubr)     # stat_compare_means/ggarrange/stat_pvalue_manual/geom_signif
library(tidyr)
library(stringr)
library(pheatmap)
library(scales)
library(ggrastr)
library(rstatix)

# 并行设置: 主流程(SCTransform/NormalizeData等)用multisession可能报错, 故用sequential
plan(sequential)
options(future.globals.maxSize = 120000 * 1024^2, future.seed = TRUE)

# 加载工具函数与配色
source(file.path(DIR_ROOT, "R", "custom_functions.R"))
source(file.path(DIR_ROOT, "R", "custom_colors.R"))

# ---- 1. 三数据集加载与统一QC --------------------------------------------------
cat("=== Step 1: Loading and QC 3 datasets ===\n")

# 对象1: 本项目KLTN轨迹分析子集(已QC, 01脚本产出)
se_trajectory <- readRDS(file.path(DIR_OBJ_SEURAT, "se.rds"))
cat("Loaded trajectory dataset:", ncol(se_trajectory), "cells\n")

# 对象2: SCC 10X数据(已QC的上皮子集; 或从原始h5重新预处理)
if (RUN_SCC_PREPROCESS) {
  sc_data <- Read10X_h5(file.path(DIR_DATA, "scc_tong_2026", "filtered_feature_bc_matrix.h5"))
  sc_obj <- CreateSeuratObject(counts = sc_data,
                                project = "scc_2026",
                                min.cells = 3,
                                min.features = 200)

  # 计算QC指标(小鼠基因命名规范)
  sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
  sc_obj[["percent.ribo"]] <- PercentageFeatureSet(sc_obj, pattern = "^Rp[sl]")
  sc_obj[["percent.hb"]] <- PercentageFeatureSet(sc_obj, pattern = "^Hb[^(p)]")

  # QC小提琴图
  pdf(file.path(DIR_FIG_RAW_MIXED, "scc_qc_violin.pdf"), width = 10, height = 6)
  VlnPlot(sc_obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
          ncol = 5,
          pt.size = 0.01) +
    theme(axis.title.x = element_blank())
  dev.off()

  # 过滤低质量细胞
  scc_obj <- subset(sc_obj,
                    subset = nFeature_RNA > 200 &
                      nFeature_RNA < 6500 &
                      percent.mt < 15 &
                      percent.hb < 1)
  cat("Filtered SCC dataset:", ncol(scc_obj), "cells retained\n")

  # 降维聚类
  scc_obj <- NormalizeData(scc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  scc_obj <- FindVariableFeatures(scc_obj, selection.method = "vst", nfeatures = 2000)
  all_genes <- rownames(scc_obj)
  scc_obj <- ScaleData(scc_obj, features = all_genes)
  scc_obj <- RunPCA(scc_obj, features = VariableFeatures(object = scc_obj), verbose = FALSE)
  ElbowPlot(scc_obj, ndims = 30)
  scc_obj <- FindNeighbors(scc_obj, dims = 1:20)
  scc_obj <- FindClusters(scc_obj, resolution = 0.2)
  scc_obj <- RunUMAP(scc_obj, dims = 1:20, seed.use = 123)
  DimPlot(scc_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE, label.size = 4)

  # 经典谱系标记定大类
  markers <- c("Ptprc",    # immune
               "Pdgfra",   # stromal
               "Epcam",    # epithelium
               "Pecam1"    # endothelial
  )
  dot_plot <- DotPlot(scc_obj, features = markers, group.by = "seurat_clusters") +
    scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                          midpoint = 0) +
    theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
    labs(title = "Gene Expression Dot Plot")
  print(dot_plot)

  # 谱系注释(基于谱系标记点图逐cluster判定)
  scc_obj$first_annotation <- case_when(
    scc_obj@meta.data$seurat_clusters == 0 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 1 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 2 ~ "epithelium",
    scc_obj@meta.data$seurat_clusters == 3 ~ "endothelium",
    scc_obj@meta.data$seurat_clusters == 4 ~ "epithelium",
    scc_obj@meta.data$seurat_clusters == 5 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 6 ~ "stromal",
    scc_obj@meta.data$seurat_clusters == 7 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 8 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 9 ~ "immune",
    scc_obj@meta.data$seurat_clusters == 10 ~ "epithelium",
    scc_obj@meta.data$seurat_clusters == 11 ~ "epithelium",
    scc_obj@meta.data$seurat_clusters == 12 ~ "epithelium"
  )
  table(scc_obj$first_annotation)

  # LUAD/LUSC标记基因确认肿瘤上皮亚群
  tumor_marker <- list(luad_features = c("Ttf1", "Napsa", "Muc1", "Krt7", "Sftpa1", "Sftpc"),
                       lusc_features = c("Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt14", "Krt6a", "Krt6b", "Krt17", "Sox2", "Sfn", "Slc2a1"))
  dot_plot <- DotPlot(scc_obj, features = tumor_marker, group.by = "seurat_clusters") +
    scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                          midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Gene Expression Dot Plot")
  print(dot_plot)

  # 仅cluster 2和4为SCC特异性上皮亚群
  scc_obj <- subset(scc_obj, seurat_clusters %in% c(2, 4))
  table(scc_obj$first_annotation)
  SaveSeuratRds(scc_obj, file.path(DIR_DATA, "scc_tong_2026", "scc_epi_subset_obj.rds"))
} else {
  # 直接读取已保存的SCC上皮子集对象
  scc_obj <- readRDS(file.path(DIR_DATA, "scc_tong_2026", "scc_epi_subset_obj.rds"))
}

# 对象3: KPTracer LUAD数据(已有Cluster.Name/leiden_sub注释)
luad_obj <- readRDS(PATH_KP_LUAD_RDS)
table(luad_obj$Cluster.Name, luad_obj$leiden_sub) # 查看亚群构成以确定保留亚群

markers <- c("Ptprc",    # immune
             "Pdgfra",   # stromal
             "Epcam",    # epithelium
             "Pecam1"    # endothelial
)
dot_plot <- DotPlot(luad_obj, features = markers, group.by = "Cluster.Name") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file.path(DIR_FIG_RAW_MIXED, "kp_luad_raw_first_annotation.pdf"), width = 4, height = 8)
print(dot_plot)
dev.off()

# LUAD/LUSC标记确认肿瘤亚型
tumor_marker <- list(luad_features = c("Ttf1", "Napsa", "Muc1", "Krt7", "Sftpa1", "Sftpc"),
                     lusc_features = c("Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt14", "Krt6a", "Krt6b", "Krt17", "Sox2", "Sfn", "Slc2a1"))
dot_plot <- DotPlot(luad_obj, features = tumor_marker, group.by = "Cluster.Name") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file.path(DIR_FIG_RAW_MIXED, "kp_luad_raw_luad_lusc_markers.pdf"), width = 8, height = 8)
print(dot_plot)
dev.off()

# 仅保留leiden_sub=4(肿瘤上皮亚群)
luad_obj <- subset(luad_obj, leiden_sub %in% c(4))

# 统一使用RNA assay并添加批次标识
DefaultAssay(se_trajectory) <- "RNA"
DefaultAssay(scc_obj) <- "RNA"
DefaultAssay(luad_obj) <- "RNA"
se_trajectory$batch <- "trajectory_subset"
scc_obj$batch <- "scc_2026"
luad_obj$batch <- "kp_luad"

# ---- 2. SCT标准化与整合特征选择 ----------------------------------------------
cat("\n=== Step 2: SCTransform normalization ===\n")

# 创建数据集列表
obj_list <- list(
  trajectory = se_trajectory,
  scc = scc_obj,
  luad = luad_obj
)

# 对每个数据集独立运行SCTransform(官方推荐整合前预处理方法)
obj_list <- lapply(X = obj_list, FUN = SCTransform,
                   vst.flavor = "v2",
                   vars.to.regress = c("nCount_RNA"),
                   verbose = FALSE)

cat("\n=== Step 3: Selecting integration features ===\n")

# 选择用于整合的高变基因(官方推荐3000个)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

# 准备SCT整合
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
obj_list[["trajectory"]]@meta.data$orig.ident <- "ast"
obj_list[["luad"]]@meta.data$orig.ident <- "KP_laud"

# ---- 3. 合并、PCA与Harmony批次整合 -------------------------------------------
cat("\n=== Step 4: Merging and PCA ===\n")

# 合并3个数据集
se_merged <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  merge.data = TRUE
)

# PCA降维
VariableFeatures(se_merged) <- features
se_merged <- RunPCA(se_merged, features = features, verbose = FALSE)

# PCA elbow图选择合适的维度
pdf(file.path(DIR_FIG_RAW_MIXED, "pca_elbow.pdf"), width = 8, height = 6)
ElbowPlot(se_merged, ndims = 30) + ggtitle("PCA Elbow Plot")
dev.off()

# 选择前20个维度用于后续分析(根据elbow图调整)
dims_use <- 1:20

cat("\n=== Step 5: Harmony integration ===\n")

# Harmony批次校正(校正batch层)
se_merged <- RunHarmony(
  object = se_merged,
  group.by.vars = c("batch"),
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

# ---- 4. 全细胞聚类与一级注释 --------------------------------------------------
cat("\n=== Step 6: Dimensional reduction and clustering ===\n")

# UMAP降维(使用Harmony校正后的PCA)
se_merged <- RunUMAP(se_merged, dims = dims_use, reduction = "harmony",
                     min.dist = 0.3, spread = 1, verbose = FALSE)

# 近邻搜索与聚类
se_merged <- FindNeighbors(se_merged, dims = dims_use, reduction = "harmony", verbose = FALSE)
se_merged <- FindClusters(se_merged, resolution = 0.2, verbose = FALSE)

# 选择分辨率0.2作为默认
se_merged$seurat_clusters <- se_merged$SCT_snn_res.0.2

# 保存全细胞整合对象
SaveSeuratRds(se_merged, file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_all_final_3_samples.rds"))

# 一级注释: 经典谱系标记点图
markers <- c("Ptprc",    # immune
             "Pdgfra",   # stromal
             "Epcam",    # epithelium
             "Pecam1"    # endothelial
)
dot_plot <- DotPlot(se_merged, features = markers, group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file.path(DIR_FIG_RAW_MIXED, "dotplot_merged_seu_first_annotation.pdf"), width = 4, height = 8)
print(dot_plot)
dev.off()

# 谱系注释(基于谱系标记点图逐cluster判定)
se_merged$first_annotation <- case_when(
  se_merged@meta.data$seurat_clusters == 0 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 1 ~ "immune",
  se_merged@meta.data$seurat_clusters == 2 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 3 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 4 ~ "stromal",
  se_merged@meta.data$seurat_clusters == 5 ~ "stromal",
  se_merged@meta.data$seurat_clusters == 6 ~ "immune",
  se_merged@meta.data$seurat_clusters == 7 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 8 ~ "stromal",
  se_merged@meta.data$seurat_clusters == 9 ~ "endothelium",
  se_merged@meta.data$seurat_clusters == 10 ~ "stromal",
  se_merged@meta.data$seurat_clusters == 11 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 12 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 13 ~ "immune",
  se_merged@meta.data$seurat_clusters == 14 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 15 ~ "immune",
  se_merged@meta.data$seurat_clusters == 16 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 17 ~ "immune",
  se_merged@meta.data$seurat_clusters == 18 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 19 ~ "immune",
  se_merged@meta.data$seurat_clusters == 20 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 21 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 22 ~ "epithelium",
  se_merged@meta.data$seurat_clusters == 23 ~ "epithelium"
)

# 取上皮细胞子集
se_merged_epi <- subset(se_merged, first_annotation == "epithelium")

# ---- 5. 上皮子集二次Harmony整合与聚类 -----------------------------------------
# 上皮子集内重新PCA+Harmony, 消除上皮亚群间残余批次效应
VariableFeatures(se_merged_epi) <- features
se_merged_epi <- RunPCA(se_merged_epi, features = features, verbose = FALSE)

# 选择前20个维度用于后续分析(根据elbow图调整)
dims_use <- 1:20

se_merged_epi <- RunHarmony(
  object = se_merged_epi,
  group.by.vars = c("batch"),
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

# UMAP降维与聚类
se_merged_epi <- RunUMAP(se_merged_epi, dims = dims_use, reduction = "harmony",
                         min.dist = 0.3, spread = 1, verbose = FALSE)
se_merged_epi <- FindNeighbors(se_merged_epi, dims = dims_use, reduction = "harmony", verbose = FALSE)
se_merged_epi <- FindClusters(se_merged_epi, resolution = 0.2, verbose = FALSE)
se_merged_epi$seurat_clusters <- se_merged_epi$SCT_snn_res.0.2

# 批次因子顺序与组织学类型标注: AST=本项目轨迹子集, SCC, LUAD
se_merged_epi$batch <- factor(se_merged_epi$batch, levels = c("kp_luad", "trajectory_subset", "scc_2026"))
se_merged_epi$histology <- case_when(
  se_merged_epi$batch == "trajectory_subset" ~ "AST",
  se_merged_epi$batch == "scc_2026" ~ "SCC",
  TRUE ~ "LUAD"
)

# 保存上皮子集整合对象(本章及后续注释结果的基础对象)
SaveSeuratRds(se_merged_epi, file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_final_3_samples.rds"))
# 原脚本此处立即readRDS重读同一文件(子脚本边界), 线性运行时无需重读:
# se_merged_epi <- readRDS(file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_final_3_samples.rds"))

# ---- 6. Harmony空间K-NN参考注释映射 ------------------------------------------
# 以KP-LUAD为Reference(其Cluster.Name为已注释标签), 其余数据集为Query,
# 在Harmony降维空间中用K-NN投票转移标签, 距离加权置信度过滤低置信预测
cat("=== Step 1: Loading data ===\n")

# 原脚本此处readRDS读取第5章保存的对象; 线性运行时直接沿用内存中的se_merged_epi:
# se_merged_epi <- readRDS(file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_final_3_samples.rds"))

cat("Total cells in integrated object:", ncol(se_merged_epi), "\n")
cat("Batch distribution:\n")
print(table(se_merged_epi$batch))

cat("\nReference dataset (kp_luad) annotations:\n")
ref_cells <- se_merged_epi$batch == "kp_luad"
print(table(se_merged_epi$Cluster.Name[ref_cells]))

# 分离Reference和Query数据集
ref_seurat <- subset(se_merged_epi, batch == "kp_luad")
query_seurat <- subset(se_merged_epi, batch != "kp_luad")

cat("\nReference cells:", ncol(ref_seurat), "\n")
cat("Query cells:", ncol(query_seurat), "\n")
cat("Query batch distribution:\n")
print(table(query_seurat$batch))

cat("\n=== Step 2: K-NN label transfer in Harmony space ===\n")

# 检查是否存在Harmony降维
if (!"harmony" %in% names(se_merged_epi@reductions)) {
  stop("Harmony reduction not found in integrated object!")
}

# 提取Reference和Query的Harmony坐标
ref_coords <- Embeddings(ref_seurat, reduction = "harmony")
query_coords <- Embeddings(query_seurat, reduction = "harmony")

cat("Reference coordinates dimensions:", dim(ref_coords), "\n")
cat("Query coordinates dimensions:", dim(query_coords), "\n")

# K-NN搜索
k <- 30
cat("Using k =", k, "nearest neighbors for label transfer\n")

knn_result <- FNN::knnx.index(ref_coords, query_coords, k = k)

# Reference的细胞类型标签
ref_labels <- ref_seurat$Cluster.Name

# 每个Query细胞的k个最近邻标签
knn_labels <- matrix(ref_labels[knn_result], nrow = nrow(query_coords), ncol = k)

# 每个Query细胞到最近邻的距离
knn_distances <- FNN::knnx.dist(ref_coords, query_coords, k = k)

cat("\n=== Step 3: Computing confidence scores ===\n")

# 距离加权置信度得分(距离越近权重越高)
predict_scores <- numeric(nrow(query_coords))
predict_labels <- character(nrow(query_coords))

for (i in 1:nrow(query_coords)) {
  labels_i <- knn_labels[i, ]
  dists_i <- knn_distances[i, ]

  weights <- 1 / (dists_i + 1e-10)
  label_weight_sum <- tapply(weights, labels_i, sum)
  max_weight <- max(label_weight_sum)
  total_weight <- sum(label_weight_sum)

  predict_labels[i] <- names(label_weight_sum)[which.max(label_weight_sum)]
  predict_scores[i] <- max_weight / total_weight
}

# 置信度阈值, 低于阈值标记为Unknown
confidence_threshold <- 0.6
predict_labels_final <- ifelse(predict_scores >= confidence_threshold, predict_labels, "Unknown")

# 统计结果
cat("\nPrediction statistics:\n")
cat("Total Query cells:", length(predict_labels_final), "\n")
cat("Cells assigned to known types:", sum(predict_labels_final != "Unknown"), "\n")
cat("Cells labeled as Unknown:", sum(predict_labels_final == "Unknown"), "\n")
cat("Unknown percentage:", round(sum(predict_labels_final == "Unknown") / length(predict_labels_final) * 100, 2), "%\n")

cat("\nLabel distribution:\n")
print(table(predict_labels_final))

cat("\nConfidence score statistics:\n")
print(summary(predict_scores))

cat("\n=== Step 4: Writing results back ===\n")

# 预测标签与置信度写回Query与整合对象
query_seurat$symphony_prediction <- predict_labels_final
query_seurat$symphony_confidence <- predict_scores

se_merged_epi$symphony_prediction <- NA
se_merged_epi$symphony_confidence <- NA

se_merged_epi$symphony_prediction[colnames(query_seurat)] <- predict_labels_final
se_merged_epi$symphony_confidence[colnames(query_seurat)] <- predict_scores

# 保存注释结果(最终整合上皮对象, 下游脚本依赖)
saveRDS(query_seurat, file.path(DIR_OBJ_SEURAT, "query_symphony_annotated.rds"))
saveRDS(se_merged_epi, file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_with_symphony.rds"))
SaveSeuratRds(se_merged_epi, file = file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_final_3_samples.rds"))

cat("Results saved successfully!\n")

# ---- 6.1 注释结果可视化 -------------------------------------------------------
cat("\n=== Step 5: Visualization ===\n")

theme_custom <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
)

# 6.1.1 Query细胞UMAP-预测标签
pdf(file.path(DIR_FIG_CELL_ANN, "query_umap_predicted_labels.pdf"), width = 12, height = 10)
p1 <- DimPlot(query_seurat, reduction = "umap", group.by = "symphony_prediction",
              label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Query Cells - Predicted Labels (with confidence filter)") +
  theme_custom
print(p1)
dev.off()

# 6.1.2 Query细胞UMAP-置信度得分
pdf(file.path(DIR_FIG_CELL_ANN, "query_umap_confidence_score.pdf"), width = 12, height = 10)
p2 <- FeaturePlot(query_seurat, features = "symphony_confidence", reduction = "umap",
                  cols = c("lightgray", "darkblue"), min.cutoff = 0) +
  ggtitle("Query Cells - Mapping Confidence Score") +
  theme_custom
print(p2)
dev.off()

# 6.1.3 Query细胞UMAP-按批次显示预测标签
pdf(file.path(DIR_FIG_CELL_ANN, "query_umap_predicted_labels_by_batch.pdf"), width = 16, height = 10)
p3 <- DimPlot(query_seurat, reduction = "umap", group.by = "symphony_prediction",
              split.by = "batch", label = TRUE, label.size = 3) +
  ggtitle("Query Cells - Predicted Labels by Batch") +
  theme_custom
print(p3)
dev.off()

# 6.1.4 Reference细胞UMAP-原始注释
pdf(file.path(DIR_FIG_CELL_ANN, "reference_umap_original_labels.pdf"), width = 12, height = 10)
p4 <- DimPlot(ref_seurat, reduction = "umap", group.by = "Cluster.Name",
              label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Reference Cells (kp_luad) - Original Labels") +
  theme_custom
print(p4)
dev.off()

# 6.1.5 整合后UMAP-所有细胞预测结果
pdf(file.path(DIR_FIG_CELL_ANN, "integrated_umap_all_predictions.pdf"), width = 14, height = 12)
p5 <- DimPlot(se_merged_epi, reduction = "umap", group.by = "symphony_prediction",
              label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Integrated Dataset - Predicted Labels") +
  theme_custom
print(p5)
dev.off()

# 6.1.6 整合后UMAP-Reference vs Query
pdf(file.path(DIR_FIG_CELL_ANN, "integrated_umap_reference_vs_query.pdf"), width = 12, height = 10)
se_merged_epi$dataset_type <- ifelse(se_merged_epi$batch == "kp_luad", "Reference", "Query")
p6 <- DimPlot(se_merged_epi, reduction = "umap", group.by = "dataset_type") +
  ggtitle("Integrated Dataset - Reference vs Query") +
  theme_custom
print(p6)
dev.off()

# 6.1.7 置信度得分直方图
pdf(file.path(DIR_FIG_CELL_ANN, "confidence_score_histogram.pdf"), width = 10, height = 6)
p7 <- ggplot(query_seurat@meta.data, aes(x = symphony_confidence)) +
  geom_histogram(bins = 50, fill = "#4A90D9", color = "black") +
  geom_vline(xintercept = confidence_threshold, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = confidence_threshold + 0.05, y = max(table(cut(query_seurat$symphony_confidence, breaks = 50))) * 0.9,
           label = paste0("Threshold = ", confidence_threshold), color = "red", fontface = "bold") +
  xlab("Mapping Confidence Score") +
  ylab("Number of Cells") +
  ggtitle("Distribution of Mapping Confidence Scores") +
  theme_bw() +
  theme_custom
print(p7)
dev.off()

# 6.1.8 各标签的置信度得分箱线图
pdf(file.path(DIR_FIG_CELL_ANN, "confidence_score_by_label_boxplot.pdf"), width = 14, height = 8)
p8 <- ggplot(query_seurat@meta.data, aes(x = symphony_prediction, y = symphony_confidence, fill = symphony_prediction)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  geom_hline(yintercept = confidence_threshold, color = "red", linetype = "dashed") +
  xlab("Predicted Cell Type") +
  ylab("Mapping Confidence Score") +
  ggtitle("Confidence Score Distribution by Predicted Label") +
  theme_bw() +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p8)
dev.off()

# 6.1.9 各批次预测标签分布堆叠图
pdf(file.path(DIR_FIG_CELL_ANN, "label_distribution_by_batch_stacked.pdf"), width = 10, height = 6)
label_counts <- query_seurat@meta.data %>%
  group_by(batch, symphony_prediction) %>%
  summarise(count = n()) %>%
  group_by(batch) %>%
  mutate(percentage = count / sum(count) * 100)

p9 <- ggplot(label_counts, aes(x = batch, y = percentage, fill = symphony_prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  xlab("Query Batch") +
  ylab("Percentage of Cells") +
  ggtitle("Predicted Label Distribution by Query Batch") +
  theme_bw() +
  theme_custom +
  scale_fill_discrete(name = "Predicted Cell Type")
print(p9)
dev.off()

# 6.1.10 Reference和Query标签一致性对比图
pdf(file.path(DIR_FIG_CELL_ANN, "label_consistency_comparison.pdf"), width = 12, height = 10)
ref_labels_table <- table(ref_seurat$Cluster.Name)
query_labels_table <- table(query_seurat$symphony_prediction)

all_labels <- union(names(ref_labels_table), names(query_labels_table))

ref_counts <- ref_labels_table[all_labels]
ref_counts[is.na(ref_counts)] <- 0
query_counts <- query_labels_table[all_labels]
query_counts[is.na(query_counts)] <- 0

consistency_df <- data.frame(
  Label = rep(all_labels, 2),
  Dataset = rep(c("Reference", "Query"), each = length(all_labels)),
  Count = c(ref_counts, query_counts),
  Percentage = c(ref_counts / sum(ref_counts) * 100, query_counts / sum(query_counts) * 100)
)

p10 <- ggplot(consistency_df, aes(x = Label, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Percentage, 1)),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  xlab("Cell Type") +
  ylab("Percentage") +
  ggtitle("Label Distribution Consistency") +
  theme_bw() +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p10)
dev.off()

# ---- 6.2 注释统计表导出 -------------------------------------------------------
cat("\n=== Step 6: Exporting statistics ===\n")

# 6.2.1 预测标签统计
prediction_stats <- query_seurat@meta.data %>%
  group_by(symphony_prediction) %>%
  summarise(
    count = n(),
    percentage = n() / nrow(query_seurat@meta.data) * 100,
    mean_confidence = mean(symphony_confidence, na.rm = TRUE),
    median_confidence = median(symphony_confidence, na.rm = TRUE),
    min_confidence = min(symphony_confidence, na.rm = TRUE),
    max_confidence = max(symphony_confidence, na.rm = TRUE)
  ) %>%
  arrange(desc(count))

write.csv(prediction_stats, file.path(DIR_FIG_CELL_ANN, "prediction_statistics.csv"), row.names = FALSE)

# 6.2.2 各批次预测统计
batch_prediction_stats <- query_seurat@meta.data %>%
  group_by(batch, symphony_prediction) %>%
  summarise(
    count = n(),
    percentage = n() / sum(query_seurat$batch == batch) * 100,
    mean_confidence = mean(symphony_confidence, na.rm = TRUE)
  ) %>%
  arrange(batch, desc(count))

write.csv(batch_prediction_stats, file.path(DIR_FIG_CELL_ANN, "batch_prediction_statistics.csv"), row.names = FALSE)

# 6.2.3 细胞水平预测结果
cell_predictions <- query_seurat@meta.data %>%
  dplyr::select(orig.ident, batch, seurat_clusters, symphony_prediction, symphony_confidence)

write.csv(cell_predictions, file.path(DIR_FIG_CELL_ANN, "cell_level_predictions.csv"), row.names = TRUE)

cat("\n=== Analysis completed successfully! ===\n")
cat("Results saved to:", DIR_FIG_CELL_ANN, "\n")

# ---- 7. 参考注释结果补充图(figJ) ----------------------------------------------
# AST(本项目)细胞按symphony预测标签着色, LUAD/SCC参考细胞置灰作为背景
# 注: 原脚本使用group列判断(其值与histology一致, 由第8章写入持久化对象);
#     此处直接使用histology以保证本章可独立于第8章运行
se_merged_epi$plot_group <- ifelse(
  se_merged_epi$histology == "AST",
  se_merged_epi$symphony_prediction,
  "LUAD/SCC"
)

# 所有类别
all_levels <- unique(se_merged_epi$plot_group)

# 背景类别置灰
background_groups <- c("LUAD/SCC", "Unknown")
foreground_groups <- setdiff(all_levels, background_groups)

# 彩色类别用默认hue配色
foreground_cols <- scales::hue_pal()(length(foreground_groups))
names(foreground_cols) <- foreground_groups

plot_cols <- c(
  "Unknown" = "#BDBDBD",
  "LUAD/SCC" = "lightgrey",
  foreground_cols
)

# 按实际类别重新匹配颜色
plot_cols <- plot_cols[all_levels]

# 绘图: order仅放彩色组, 使其覆盖灰色背景点
p <- DimPlot(
  se_merged_epi,
  group.by = "plot_group",
  cols = plot_cols,
  order = foreground_groups,
  shuffle = FALSE
)
p_raster_on_save <- rasterise(p, layers = "Point", dpi = 1200)
pdf(file.path(DIR_FIG_CELL_ANN, "reference_based_KLTN_cell_annotation.pdf"), width = 7, height = 5)
print(p_raster_on_save)
dev.off()

# ---- 8. 整合对象可视化与AST过渡态验证 -----------------------------------------
# 读取第5-6章保存的整合上皮对象(含symphony注释);
# 原脚本此章节为独立子脚本, 从此处readRDS载入
se_merged <- readRDS(file.path(DIR_OBJ_SEURAT, "se_harmony_integrated_epi_subset_final_3_samples.rds"))

# 生成分组列: AST(本项目KLTN模型)/LUAD/SCC
# 注: 原脚本此处条件为'KLTN'(中途改名残留), 对象实际histology取值为AST, 已统一为AST
se_merged$group <- ifelse(se_merged$histology == "AST", "AST", se_merged$histology)

# ---- 8.1 Trp63靶基因可视化 ----------------------------------------------------
# Trp63靶基因列表在LUAD/AST/SCC三组的表达
feature_list <- c("B3galt2", "Baiap2l2", "Foxj1", "Krt4", "Lrp4",
                  "Nrg1", "Pifo", "Plcb1", "Rassf9", "Rgs12", "Traf3ip1", "Trp63")

# 组别因子顺序按生物学进程排列
se_merged$group <- factor(se_merged$group, levels = c("LUAD", "AST", "SCC"))

p_dot <- DotPlot(se_merged, features = feature_list, group.by = "group") +
  # Magma连续调色板(深紫→亮黄)
  scale_color_viridis_c(option = "magma", name = "Average\nExpression") +
  theme_classic(base_size = 14) +
  theme(
    # X轴(基因名)倾斜45度防斜体规范
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  labs(size = "Percent\nExpressed")

pdf(file = file.path(FIG_TRANS, "trp63_target_genes_dotplot.pdf"), width = 6, height = 4)
print(p_dot)
dev.off()

# Trp63与Nrg1两基因点图
feature_list <- c("Trp63", "Nrg1")
p_dot <- DotPlot(se_merged, features = feature_list, group.by = "group") +
  # Nature经典蓝红连续色板
  scale_color_gradient(low = "#2b83ba", high = "#d7191c", name = "Average\nExpression") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  labs(size = "Percent\nExpressed")

pdf(file = file.path(FIG_TRANS, "trp63_nrg1_genes_dotplot.pdf"), width = 4, height = 5)
print(p_dot)
dev.off()

# Nrg1的UMAP表达图
p <- FeaturePlot(se_merged, features = "Nrg1", order = TRUE, pt.size = 0.8) +
  scale_color_gradient2(
    low = "#2b83ba",     # 低表达蓝色
    mid = "white",       # 零值白色
    high = "#d7191c",    # 高表达红色
    midpoint = 0.5,
    name = "Nrg1\nExpression"
  ) +
  theme_void() +
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5))
pdf(file = file.path(FIG_TRANS, "umap_nrg1_expression_integrated_obj.pdf"), width = 6, height = 5)
print(p)
dev.off()

# ---- 8.2 轨迹相关signature验证 ------------------------------------------------
# 读取03章tradeSeq的关联结果, 取top50上调轨迹基因计算module得分
trajectory_associated_genes <- read.csv(file.path(DIR_OBJ, "td_ratio", "association_results.csv"),
                                        header = TRUE, row.names = 1)
top50_upregulated_trajectory_associated_genes <- trajectory_associated_genes |>
  filter(pvalue < 0.05) |>
  top_n(n = 50, wt = meanLogFC) |>
  rownames()
se_merged <- AddModuleScore(se_merged, features = list(top50_upregulated_trajectory_associated_genes), name = "trajectory_up")

pdf(file = file.path(FIG_TRANS, "umap_module_trajectory_up_integrated_obj.pdf"), width = 6, height = 5)
FeaturePlot(se_merged, features = "trajectory_up1", pt.size = 0.6) +
  scale_color_gradient2(
    low = "#2b83ba",
    mid = "#e0e7f2",
    high = "#d7191c",
    midpoint = 0,
    name = "Signature\nScore"
  ) +
  ggtitle("Module_trajectory_up") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right"
  )
dev.off()

# ---- 8.3 MIF信号轴基因可视化 --------------------------------------------------
feature_list <- c("Mif", "Cd74", "Cd44", "Ackr3", "Trp63", "Nrg1")
p_list <- list()
for (i in 1:length(feature_list)) {
  gene <- feature_list[[i]]

  p <- FeaturePlot(se_merged,
                   features = gene,
                   order = TRUE,
                   pt.size = 0.8) +
    # 浅蓝→白→朱红配色
    scale_color_gradientn(
      colors = c("#e0e7f2", "#FFFFFF", "#d7191c"),
      name = sprintf("%s\nExpression", gene)
    ) +
    ggtitle(gene) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 14),
      legend.position = "right"
    )

  p_list[[i]] <- p
}
pdf(file = file.path(FIG_TRANS, "umap_gene_list_integrated_obj.pdf"), width = 6, height = 5)
print(p_list)
dev.off()

p_dot <- DotPlot(se_merged, features = feature_list, group.by = "group") +
  scale_color_viridis_c(option = "magma", name = "Average\nExpression") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  labs(size = "Percent\nExpressed")
pdf(file = file.path(FIG_TRANS, "gene_list_genes_dotplot.pdf"), width = 6, height = 5)
print(p_dot)
dev.off()

# ---- 8.4 通路signature得分箱线图 ----------------------------------------------
# 依赖Module_*评分列(由第9章AddModuleScore写入持久化对象; 从头完整重跑时需先运行第9章)
module_scores_name <- grep("Module_", colnames(se_merged@meta.data), value = TRUE)
if (length(module_scores_name) > 0) {
  # 列名去掉"Module_"前缀
  new_names <- str_remove(module_scores_name, "^Module_")
  colnames(se_merged@meta.data)[colnames(se_merged@meta.data) %in% module_scores_name] <- new_names
  module_scores_name <- new_names
  se_merged$group <- ifelse(se_merged$group == "AST", "AST", se_merged$group)

  # 8.4.1 全部signature箱线图(带两两Wilcoxon检验)
  plot_data <- se_merged@meta.data %>%
    dplyr::select(group, all_of(module_scores_name))

  # 宽数据转长数据(ggplot2必须长格式)
  plot_data_long <- pivot_longer(
    data = plot_data,
    cols = -group,
    names_to = "Signature",
    values_to = "Expression"
  )
  plot_data_long$group <- factor(plot_data_long$group,
                                 levels = c("AST", "LUAD", "SCC"))
  # 自动生成所有两两比较组合
  my_comparisons <- combn(levels(plot_data_long$group), 2, simplify = FALSE)

  p_pval <- ggplot(plot_data_long, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~Signature, scales = "free_y", ncol = 5) +
    stat_compare_means(
      comparisons = my_comparisons,  # 两两比较组合列表
      method = "wilcox.test",        # Wilcoxon秩和检验
      label = "p.signif",            # 显示显著性标记(*, **, ***, ****)
      hide.ns = TRUE,                # 隐藏不显著(ns)的比较连线和标记
      step.increase = 0.1, vjust = 0.8
    ) +
    scale_fill_manual(values = c("LUAD" = "#377eb8", "AST" = "#e41a1c", "SCC" = "#4daf4a")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8),
      legend.position = "none",
      plot.margin = margin(t = 50, r = 5, b = 5, l = 5)
    ) +
    labs(
      x = "Annotation",
      y = "Signature Expression"
    )

  pdf(file.path(DIR_FIG_INTEGRATION, "all_canonical_signature_group_boxplot_with_pval.pdf"),
      width = 10, height = 12)
  print(p_pval)
  dev.off()

  # 8.4.2 选定signature箱线图(带检验)
  selectced_module <- c("AT2-like", "AT1/AT2", "Alveolar",
                        "Glandular", "Basal", "Squamous",
                        "Biosynthetic mixed identity",
                        "High-plasticity cell state")
  plot_data <- se_merged@meta.data %>%
    dplyr::select(group, all_of(selectced_module))

  plot_data_long <- pivot_longer(
    data = plot_data,
    cols = -group,
    names_to = "Signature",
    values_to = "Expression"
  )
  plot_data_long$group <- factor(plot_data_long$group,
                                 levels = c("AST", "LUAD", "SCC"))
  my_comparisons <- combn(levels(plot_data_long$group), 2, simplify = FALSE)

  p_pval <- ggplot(plot_data_long, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~Signature, scales = "free_y", ncol = 4) +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      step.increase = 0.1, vjust = 0.8
    ) +
    scale_fill_manual(values = c("LUAD" = "#377eb8", "AST" = "#e41a1c", "SCC" = "#4daf4a")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 8, color = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    ) +
    labs(
      x = "Annotation",
      y = "Signature Expression"
    )

  pdf(file.path(DIR_FIG_INTEGRATION, "selected_canonical_signature_group_boxplot_with_pval.pdf"),
      width = 12, height = 4)
  print(p_pval)
  dev.off()

  # 8.4.3 选定signature箱线图(不带检验, 主图版)
  p <- ggplot(plot_data_long, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~Signature, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("LUAD" = "#377eb8", "SCC" = "#e41a1c", "AST" = "#4daf4a")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 8, color = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black")
    ) +
    labs(
      x = "Annotation",
      y = "Signature Expression"
    )

  pdf(file.path(DIR_FIG_INTEGRATION, "selected_canonical_signature_group_boxplot_without_pval.pdf"),
      width = 12, height = 4)
  print(p)
  dev.off()
} else {
  message("未检测到Module_*评分列(需先运行第9章或读入含module评分的RDS), 跳过signature箱线图")
}

# ---- 8.5 LUAD/LUSC标记基因点图 ------------------------------------------------
tumor_marker <- list(luad_features = c("Ttf1", "Napsa", "Muc1", "Krt7", "Sftpa1", "Sftpc"),
                     lusc_features = c("Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt14", "Krt6a", "Krt6b", "Krt17", "Sox2", "Sfn", "Slc2a1"))
dot_plot <- DotPlot(se_merged, features = tumor_marker,
                    group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(FIG_TRANS, "mergerd_seurat_KLTN_genes_dotplot.pdf"), width = 6, height = 4)
print(dot_plot)
dev.off()

# ---- 8.6 AST过渡态验证 --------------------------------------------------------
cat("\n=== Step 11: Verifying KLTN as transition state between LUAD and SCC ===\n")

# 重新分组: LUAD组、SCC组、AST(本项目)组, 过滤无分组样本
se_merged$group <- se_merged$histology
se_merged <- subset(se_merged, !is.na(group))

# 8.6.1 验证1: UMAP空间分布验证
cat("\n1. UMAP spatial distribution validation...\n")

# 按分组着色的UMAP
p <- CellDimPlot(
  se_merged, pt.size = 2.5,
  group.by = "group",
  palcolor = c(
    "LUAD" = "#377eb8",
    "AST" = "#e41a1c",
    "SCC" = "#4daf4a"
  ),
  theme_use = "theme_blank",
  raster = TRUE,
  raster.dpi = c(1800, 1500)   # 6 x 5 inch, 约300 ppi
) +
  ggtitle("UMAP distribution of LUAD, KLTN, SCC groups")
ggsave(
  filename = file.path(FIG_TRANS, "umap_group.pdf"),
  plot = p,
  width = 6,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf
)
se_merged$group <- factor(se_merged$group, levels = c("AST", "LUAD", "SCC"))

# 8.6.2 密度分布叠加图, 看三组分布重叠情况
umap_all <- as.data.frame(Embeddings(se_merged, "umap"))
colnames(umap_all) <- c("UMAP_1", "UMAP_2")
umap_all$group <- se_merged$group

p1 <- ggplot(umap_all, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  stat_density_2d(aes(color = group), alpha = 0.7, linewidth = 0.6) +
  scale_color_manual(values = c("SCC" = "#4daf4a", "AST" = "#e41a1c", "LUAD" = "#377eb8")) +
  # 限制显示区间, 不截断密度计算
  scale_x_continuous(limits = c(-15, 15)) +
  scale_y_continuous(limits = c(-15, 15)) +
  theme_classic() +
  ggtitle("Density overlay of three groups")

pdf(file.path(FIG_TRANS, "density_overlay.pdf"), width = 5, height = 4)
print(p1)
dev.off()
pdf(file.path(FIG_TRANS, "density_overlay_wider.pdf"), width = 9, height = 5)
print(p1)
dev.off()

# 8.6.3 验证2: 分子特征得分验证
cat("\n2. Molecular signature score validation...\n")

# LUAD和SCC特征基因集(来自项目背景和文献)
luad_signature <- c("Mlph", "Tmc5", "Sfta3-ps", "Muc5b", "Ttf1", "Napsa", "Muc1", "St6galnac1", "Plin2", "Sftpc")
scc_signature <- c("Krt14", "Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt6a", "Trim29", "Slc7a5", "Spats2", "Clca2")

pdf(file.path(FIG_TRANS, "KLTN_related_genes.pdf"), width = 8, height = 8)
VlnPlot(se_merged, features = luad_signature, group.by = "group")
VlnPlot(se_merged, features = scc_signature, group.by = "group")
dev.off()

# 过滤掉不存在的基因
luad_signature <- luad_signature[luad_signature %in% rownames(se_merged)]
scc_signature <- scc_signature[scc_signature %in% rownames(se_merged)]

# 每个细胞的平均表达得分(使用SCT标准化数据)
DefaultAssay(se_merged) <- "SCT"
se_merged$LUAD_score <- colMeans(GetAssayData(se_merged)[luad_signature, ])
se_merged$SCC_score <- colMeans(GetAssayData(se_merged)[scc_signature, ])

# 得分分布对比
pdf(file.path(FIG_TRANS, "signature_score_distribution.pdf"), width = 6, height = 8)
p1 <- VlnPlot(se_merged, features = "LUAD_score", group.by = "group",
              cols = c("AST" = "#4daf4a", "SCC" = "#e41a1c", "LUAD" = "#377eb8"),
              pt.size = 0) + NoLegend() + ggtitle("LUAD signature score")
p2 <- VlnPlot(se_merged, features = "SCC_score", group.by = "group",
              cols = c("AST" = "#4daf4a", "SCC" = "#e41a1c", "LUAD" = "#377eb8"),
              pt.size = 0) + NoLegend() + ggtitle("SCC signature score")
p1 + p2
dev.off()
pdf(file.path(FIG_TRANS, "LUAD_signature_score_distribution.pdf"), width = 9, height = 5)
p1
dev.off()
pdf(file.path(FIG_TRANS, "LUSS_signature_score_distribution.pdf"), width = 9, height = 5)
p2
dev.off()

# 得分散点图, 看三组分布关系
p <- ggplot(
  se_merged@meta.data,
  aes(x = LUAD_score, y = SCC_score, color = group)
) +
  geom_point_rast(
    alpha = 0.5,
    size = 0.8,
    raster.dpi = 300
  ) +
  geom_density_2d(alpha = 0.8) +
  scale_color_manual(
    values = c(
      "SCC" = "#4daf4a",
      "AST" = "#e41a1c",
      "LUAD" = "#377eb8"
    )
  ) +
  theme_classic() +
  labs(
    title = "LUAD vs SCC signature score distribution",
    x = "LUAD signature score (normalized)",
    y = "SCC signature score (normalized)"
  )

ggsave(
  filename = file.path(FIG_TRANS, "signature_score_scatter.pdf"),
  plot = p,
  width = 5,
  height = 4,
  units = "in",
  device = grDevices::cairo_pdf
)

pdf(file.path(FIG_TRANS, "signature_score_scatter_wider.pdf"), width = 9, height = 5)
ggplot(se_merged@meta.data, aes(x = LUAD_score, y = SCC_score, color = group)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_density_2d(alpha = 0.8) +
  scale_color_manual(values = c("AST" = "#4daf4a", "SCC" = "#e41a1c", "LUAD" = "#377eb8")) +
  theme_classic() +
  labs(title = "LUAD vs SCC signature score distribution",
       x = "LUAD signature score (normalized)",
       y = "SCC signature score (normalized)")
dev.off()

# 8.6.4 统计检验: 验证AST得分显著介于LUAD和SCC之间
score_stats <- se_merged@meta.data %>%
  group_by(group) %>%
  summarise(
    LUAD_mean = mean(LUAD_score),
    LUAD_sd = sd(LUAD_score),
    SCC_mean = mean(SCC_score),
    SCC_sd = sd(SCC_score)
  )
write.csv(score_stats, file.path(DIR_OBJ_SEURAT, "transition_score_stats.csv"), row.names = FALSE)

# t检验验证AST与LUAD/SCC均有显著差异且处于中间
luad_ast_ttest <- t.test(se_merged$LUAD_score[se_merged$group == "LUAD"],
                         se_merged$LUAD_score[se_merged$group == "AST"])
scc_ast_ttest <- t.test(se_merged$SCC_score[se_merged$group == "SCC"],
                        se_merged$SCC_score[se_merged$group == "AST"])
cat("\nLUAD vs KLTN t-test p-value:", luad_ast_ttest$p.value, "\n")
cat("SCC vs KLTN t-test p-value:", scc_ast_ttest$p.value, "\n")

# 8.6.5 关键调控因子表达趋势: Trp63、Sox2等
p1 <- FeaturePlot(se_merged, features = c("Trp63", "Sox2", "Krt5", "Sftpc", "Napsa", "Ttf1"),
                  reduction = "umap", ncol = 3)
p2 <- VlnPlot(se_merged, features = c("Trp63", "Sox2", "Krt5", "Sftpc", "Napsa", "Ttf1"), group.by = "group", ncol = 3)
pdf(file.path(FIG_TRANS, "key_gene_expression_umap.pdf"), width = 12, height = 8)
print(p1)
dev.off()
pdf(file.path(FIG_TRANS, "key_gene_expression_vln.pdf"), width = 12, height = 8)
print(p2)
dev.off()

# ---- 9. GMT基因集module评分与通路箱线图(整合对象) ------------------------------
# 读取GMT基因集→人源转小鼠同源→AddModuleScore→三组箱线图
set.seed(1) # 原脚本该子章节独立运行时的随机种子
plan(multisession, workers = 8)
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)

cat("Loading data...\n")
# 本章沿用第8章内存中的se_merged(整合上皮对象)
cat("Loaded integrated epi object - cells:", ncol(se_merged), ", features:", nrow(se_merged), "\n")

# 读取GMT基因集
gmt_file <- file.path(DIR_DATA, "geneset.hs.gmt")
gmt_data <- readLines(gmt_file)
cat("Loaded geneset.hs.gmt -", length(gmt_data), "gene sets\n")

cat("\nParsing GMT file and converting to mouse genes...\n")

# 解析GMT文件(名称\t描述\t基因...)
parse_gmt <- function(gmt_lines) {
  gene_sets <- list()
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    name <- parts[1]
    genes <- parts[-c(1, 2)]
    genes <- genes[genes != ""]
    gene_sets[[name]] <- genes
  }
  return(gene_sets)
}

human_gene_sets <- parse_gmt(gmt_data)

# 人源基因集转换为小鼠同源基因
mouse_gene_sets <- list()
for (name in names(human_gene_sets)) {
  human_genes <- human_gene_sets[[name]]
  if (length(human_genes) >= 3) {
    tryCatch({
      mouse_genes <- convert_orthologs(
        gene_df = human_genes,
        input_species = "human",
        output_species = "mouse",
        method = "homologene"
      )
      mouse_gene_sets[[name]] <- rownames(mouse_genes)
    }, error = function(e) {})
  }
}
mouse_gene_sets <- mouse_gene_sets[sapply(mouse_gene_sets, length) >= 3]
cat("Loaded", length(mouse_gene_sets), "mouse gene sets\n")

# 全部基因集计算module得分
key_pathways <- names(mouse_gene_sets)
cat("Analyzing pathways:", paste(key_pathways, collapse = ", "), "\n")

gene_list_list <- lapply(key_pathways, function(p) mouse_gene_sets[[p]])
names(gene_list_list) <- key_pathways

se_merged <- AddModuleScore(
  se_merged,
  features = gene_list_list,
  name = paste0("Module_", key_pathways, "_"),
  ctrl = 100  # 100个控制基因
)

# 列名去掉AddModuleScore自动添加的数字后缀
module_scores_name <- grep("Module_", colnames(se_merged@meta.data), value = TRUE)
colnames(se_merged@meta.data) <- sub("^(Module_.+)_\\d+$", "\\1", colnames(se_merged@meta.data))
module_scores_name <- grep("Module_", colnames(se_merged@meta.data), value = TRUE)

# 9.1 各通路UMAP图
plot_list <- list()
for (pathway in module_scores_name) {
  p <- FeaturePlot(se_merged, features = pathway, reduction = "umap") +
    ggtitle(pathway) +
    theme(plot.title = element_text(size = 10))
  plot_list[[pathway]] <- p
}

pdf(file.path(FIG_TRANS, "umap_all_signature.pdf"), width = 24, height = 4)
ggarrange(plotlist = plot_list, ncol = 6)
dev.off()

# 9.2 全部通路得分箱线图
plot_data <- se_merged@meta.data %>%
  dplyr::select(group, all_of(module_scores_name))

plot_data_long <- pivot_longer(
  data = plot_data,
  cols = -group,
  names_to = "Signature",
  values_to = "Expression"
)

p <- ggplot(plot_data_long, aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~Signature, scales = "free_y", ncol = 5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "DNB State",
    y = "Signature Expression"
  )

pdf(
  file.path(FIG_TRANS, "all_canonical_signature_group_boxplot.pdf"),
  width = 12,
  height = 10
)
print(p)
dev.off()

# 9.3 主图选定通路箱线图
select_pathway <- grep("mixed|plasticity|AT2-like|Basal|Squamous|Glandular", module_scores_name, value = TRUE)
plot_data <- se_merged@meta.data %>%
  dplyr::select(group, all_of(select_pathway))

plot_data_long <- pivot_longer(
  data = plot_data,
  cols = -group,
  names_to = "Signature",
  values_to = "Expression"
)

p <- ggplot(plot_data_long, aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~Signature, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "DNB State",
    y = "Signature Expression"
  )

pdf(
  file.path(FIG_TRANS, "selected_signature_group_boxplot.pdf"),
  width = 10,
  height = 4
)
print(p)
dev.off()

# 保存含module评分的整合对象(注意: 原脚本此处保存于output/object/根目录而非seurat/子目录)
saveRDS(se_merged, file.path(DIR_OBJ, "se_harmony_integrated_epi_subset_final_3_samples.rds"))

# 导出最终metadata
write.csv(se_merged@meta.data,
          file.path(DIR_OBJ, "se_harmony_integrated_epi_subset_final_3_samples_analysis_meta.csv"), row.names = TRUE)

# ---- 10. AST轨迹子集: td组DEG与EMT/通路module评分箱线图 ------------------------
# 在本项目轨迹子集(se_trajectory_subset)上分析tdTomato/GFP两组的
# 差异基因、EMT评分及GMT通路module得分
set.seed(1) # 原脚本该子章节独立运行时的随机种子
plan(multisession, workers = 8)
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)

cat("Loading data...\n")

# 读取轨迹子集对象(03章产出; 本章末尾将覆盖保存以追加Module_*评分列)
se_obj_traj <- readRDS(file.path(DIR_OBJ_SEURAT, "se_trajectory_subset.rds"))
se_obj_traj$seurat_clusters <- droplevels(se_obj_traj$seurat_clusters)

# 聚类着色UMAP示意图
cols_seurat <- c("0" = "#6EB542", "1" = "pink", "2" = "#e0e896", "3" = "lightblue")
pdf(file = file.path(DIR_FIG_GSVA, "umap_illustration_to_show.pdf"), width = 6, height = 5)
DimPlot(se_obj_traj, group.by = "seurat_clusters", cols = cols_seurat,
        reduction = "umap") + theme_void()
dev.off()
cat("Loaded se_trajectory_subset.rds - cells:", ncol(se_obj_traj), ", features:", nrow(se_obj_traj), "\n")

# 10.1 按tdTomato_ratio阈值0.5分组: tdTomato(TRP63高) vs GFP(TRP63低)
se_obj_traj$td_group <- ifelse(se_obj_traj$tdTomato_ratio >= 0.5, "tdTomato", "GFP")
Idents(se_obj_traj) <- "td_group"

# 两组差异基因
td_vs_gfp <- FindMarkers(se_obj_traj, ident.1 = "tdTomato", ident.2 = "GFP")
td_vs_gfp_filter <- td_vs_gfp |>
  filter(p_val_adj < 0.05)
write.csv(td_vs_gfp_filter, file = file.path(DIR_OBJ_GSVA, "deg_td_vs_gfp_filter.csv"))

pdf(file = file.path(DIR_FIG_GSVA, "umap_td_gfp.pdf"), width = 6, height = 5)
CellDimPlot(se_obj_traj, group.by = "td_group", palcolor = c("#6EB542", "pink"),
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()

# 10.2 全AST对象(含免疫细胞)荧光分组UMAP
se_all <- readRDS(file.path(DIR_OBJ_SEURAT, "se.rds"))
se_all$td_group <- ifelse(se_all$tdTomato_ratio >= 0.5, "tdTomato", "GFP")
se_all$td_group_plot <- ifelse(se_all$first_annotation == "immune", "immune", se_all$td_group)
pdf(file = file.path(DIR_FIG_GSVA, "seu_all_umap_td_gfp.pdf"), width = 6, height = 5)
CellDimPlot(se_all, group.by = "td_group_plot", palcolor = c("#6EB542", "pink", "lightgrey"),
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()

# 10.3 EMT相关module评分(上皮标记 vs EMT标记)
gene_list <- list(Epithelium = c("Cdh1", "Epcam", "Cldn3", "Cldn4", "Krt8"),
                  EMT = c("Vim", "Fn1", "Cdh2", "Acta2", "S100a4"))
key_pathways <- c("Epithelium", "EMT")
se_obj_traj <- AddModuleScore(
  se_obj_traj,
  features = gene_list,
  name = paste0("Module_", key_pathways, "_"),
  ctrl = 100  # 100个控制基因
)
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
colnames(se_obj_traj@meta.data) <- sub("^(Module_.+)_\\d+$", "\\1", colnames(se_obj_traj@meta.data))
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)

pdf(file = file.path(DIR_FIG_GSVA, "umap_emt_signature_sore.pdf"), width = 12, height = 5)
FeatureDimPlot(se_obj_traj, features = module_scores_name, bg_cutoff = -0.1,
               theme_use = "theme_blank", legend.position = "bottom")
dev.off()

# 10.4 GMT基因集module评分(沿用第9章定义的parse_gmt与同源转换)
gmt_file <- file.path(DIR_DATA, "geneset.hs.gmt")
gmt_data <- readLines(gmt_file)

human_gene_sets <- parse_gmt(gmt_data)

# 人源基因集转换为小鼠同源基因
mouse_gene_sets <- list()
for (name in names(human_gene_sets)) {
  human_genes <- human_gene_sets[[name]]
  if (length(human_genes) >= 3) {
    tryCatch({
      mouse_genes <- convert_orthologs(
        gene_df = human_genes,
        input_species = "human",
        output_species = "mouse",
        method = "homologene"
      )
      mouse_gene_sets[[name]] <- rownames(mouse_genes)
    }, error = function(e) {})
  }
}
mouse_gene_sets <- mouse_gene_sets[sapply(mouse_gene_sets, length) >= 3]
cat("Loaded", length(mouse_gene_sets), "mouse gene sets\n")

# 全部基因集计算module得分
key_pathways <- names(mouse_gene_sets)
cat("Analyzing pathways:", paste(key_pathways, collapse = ", "), "\n")
se_obj_traj <- AddModuleScore(
  se_obj_traj,
  features = mouse_gene_sets,
  name = paste0("Module_", key_pathways, "_"),
  ctrl = 100  # 100个控制基因
)

module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
colnames(se_obj_traj@meta.data) <- sub("^(Module_.+)_\\d+$", "\\1", colnames(se_obj_traj@meta.data))
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)

# 覆盖保存轨迹子集对象(追加Module_*评分列; 该对象同时为04章DNB分析的输入)
SaveSeuratRds(se_obj_traj, file.path(DIR_OBJ_SEURAT, "se_trajectory_subset.rds"))

# 10.5 通路module得分箱线图
# 列名去掉"Module_"前缀
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
new_names <- str_remove(module_scores_name, "^Module_")
colnames(se_obj_traj@meta.data)[colnames(se_obj_traj@meta.data) %in% module_scores_name] <- new_names
module_scores_name <- new_names

############################################################
## 通路module得分箱线图
## 1. 全部通路: 一个PDF分面展示
## 2. 选定通路: 一个多页PDF, 每个pathway一页
## 3. Wilcoxon检验, 横线连接比较组, 星号居中
##
## 显著性标记:
## *     P < 0.05
## **    P < 0.01
## ***   P < 0.001
## ****  P < 0.0001
## hide_ns = TRUE: 不显示ns及不显著比较的横线
############################################################

# 参数设置
selected_pathways <- c(
  "AT1/AT2",
  "AT2-like",
  "SCC",
  "Gastric-like",
  "High-plasticity cell state",
  "Basal"
)

# 是否隐藏不显著结果(TRUE: 不显示ns及对应横线)
hide_ns <- TRUE

# Wilcoxon检验的P值校正方法("none": 原始P值)
p_adjust_method <- "none"

# 全部通路每行的分面数
all_plot_ncol <- 5

# 输出目录
output_dir <- DIR_FIG_SIG

# 准备绘图数据(td组x全部通路得分)
plot_data <- se_obj_traj@meta.data %>%
  dplyr::select(
    td_group,
    all_of(module_scores_name)
  ) %>%
  dplyr::filter(
    !is.na(td_group)
  )

# td_group显示顺序
if (is.factor(plot_data$td_group)) {
  group_levels <- levels(
    droplevels(plot_data$td_group)
  )
} else {
  group_levels <- unique(
    as.character(plot_data$td_group)
  )
}

plot_data$td_group <- factor(
  as.character(plot_data$td_group),
  levels = group_levels
)

if (length(group_levels) < 2) {
  stop("td_group中的有效分组少于两个, 无法进行比较。")
}

# 宽格式转长格式
plot_data_long <- plot_data %>%
  tidyr::pivot_longer(
    cols = -td_group,
    names_to = "Signature",
    values_to = "Expression"
  ) %>%
  dplyr::filter(
    !is.na(Expression)
  )

# 保持pathway在元数据中的原始顺序
plot_data_long$Signature <- factor(
  plot_data_long$Signature,
  levels = module_scores_name
)

# 分组颜色
base_colors <- c(
  "#6EB542",
  "#F4A6B8",
  "#4E79A7",
  "#F28E2B",
  "#59A14F",
  "#E15759",
  "#B07AA1",
  "#76B7B2"
)

if (length(group_levels) <= length(base_colors)) {
  group_colors <- base_colors[
    seq_along(group_levels)
  ]
} else {
  group_colors <- grDevices::colorRampPalette(
    base_colors
  )(
    length(group_levels)
  )
}

group_colors <- stats::setNames(
  group_colors,
  group_levels
)

# 全部两两比较组合
comparisons <- combn(
  group_levels,
  2,
  simplify = FALSE
)

# 统计检验函数: 分组Wilcoxon检验+显著性标记+横线坐标
calculate_pathway_statistics <- function(data) {
  stat_result <- data %>%
    dplyr::group_by(Signature) %>%
    rstatix::wilcox_test(
      Expression ~ td_group,
      comparisons = comparisons,
      p.adjust.method = p_adjust_method,
      exact = FALSE
    )

  # 根据原始P值或校正后的P值生成星号
  if (p_adjust_method == "none") {
    stat_result <- stat_result %>%
      rstatix::add_significance(
        p.col = "p"
      )

    stat_result$plot_significance <- stat_result$p.signif
  } else {
    stat_result <- stat_result %>%
      rstatix::add_significance(
        p.col = "p.adj"
      )

    stat_result$plot_significance <- stat_result$p.adj.signif
  }

  # 自动计算横线和星号的纵坐标(星号位于group1和group2的中央)
  stat_result <- stat_result %>%
    rstatix::add_xy_position(
      x = "td_group",
      fun = "max",
      step.increase = 0.12,
      scales = "free_y"
    )

  stat_result
}

# 对全部pathway进行统计检验
stat_all <- calculate_pathway_statistics(
  plot_data_long
)

# 通用箱线图主题
boxplot_theme <- theme_classic(
  base_size = 12
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      size = 9,
      face = "bold"
    ),
    strip.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.margin = margin(
      t = 15,
      r = 10,
      b = 10,
      l = 10
    )
  )

# 10.5.1 全部通路: 一个PDF分面展示
p_all <- ggplot(
  plot_data_long,
  aes(
    x = td_group,
    y = Expression,
    fill = td_group
  )
) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.80,
    linewidth = 0.55,
    colour = "black",
    outlier.shape = NA,  # 不显示离群点
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  facet_wrap(
    ~Signature,
    scales = "free_y",
    ncol = all_plot_ncol
  ) +
  # 比较横线和居中的星号
  ggpubr::stat_pvalue_manual(
    stat_all,
    label = "plot_significance",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    hide.ns = hide_ns,        # 不显著结果不显示, 也不留空横线
    tip.length = 0.01,        # 横线两端短竖线长度
    bracket.size = 0.40,      # 比较横线宽度
    size = 3.5,               # 星号大小
    vjust = 0.2,              # 星号相对横线的垂直位置
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.25)
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = "td_group",
    y = "Signature score"
  ) +
  boxplot_theme

# 根据pathway数量自动调整PDF高度
all_plot_nrow <- ceiling(
  length(module_scores_name) / all_plot_ncol
)
all_pdf_height <- max(
  6,
  all_plot_nrow * 3
)

ggsave(
  filename = file.path(
    output_dir,
    "all_canonical_signature_group_boxplot.pdf"
  ),
  plot = p_all,
  width = 14,
  height = all_pdf_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

# 10.5.2 选定通路: 检查pathway是否存在
selected_pathways_exist <- selected_pathways[
  selected_pathways %in% module_scores_name
]

selected_pathways_missing <- setdiff(
  selected_pathways,
  module_scores_name
)

if (length(selected_pathways_missing) > 0) {
  warning(
    paste0(
      "以下selected pathways未在元数据中找到, 将被忽略: ",
      paste(
        selected_pathways_missing,
        collapse = ", "
      )
    )
  )
}

if (length(selected_pathways_exist) == 0) {
  stop("selected_pathways中没有任何pathway存在于元数据中。")
}

# 按给定顺序排列选定通路
selected_data_long <- plot_data_long %>%
  dplyr::filter(
    as.character(Signature) %in%
      selected_pathways_exist
  ) %>%
  dplyr::mutate(
    Signature = factor(
      as.character(Signature),
      levels = selected_pathways_exist
    )
  )

# 选定通路横坐标标签: GFP -> ZsGreen⁺, tdTomato -> tdTomato⁺
selected_x_labels <- function(x) {
  labels <- lapply(x, function(group_name) {
    if (group_name == "GFP") {
      # "+"以上标形式显示
      expression(ZsGreen^"+")[[1]]
    } else if (group_name == "tdTomato") {
      expression(tdTomato^"+")[[1]]
    } else {
      group_name
    }
  })

  as.expression(labels)
}

# 10.5.3 每个选定pathway单独生成图
selected_plot_list <- list()

for (pathway in selected_pathways_exist) {
  pathway_data <- selected_data_long %>%
    dplyr::filter(
      as.character(Signature) == pathway
    ) %>%
    droplevels()

  # 当前pathway单独计算统计检验及标注位置
  pathway_stat <- calculate_pathway_statistics(
    pathway_data
  )

  p_pathway <- ggplot(
    pathway_data,
    aes(
      x = td_group,
      y = Expression,
      fill = td_group
    )
  ) +
    geom_boxplot(
      width = 0.62,
      alpha = 0.85,
      linewidth = 0.65,
      colour = "black",
      outlier.shape = NA,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = group_colors,
      drop = FALSE
    ) +
    # 选定pathway的横坐标显示
    scale_x_discrete(
      labels = selected_x_labels
    ) +
    # 比较横线及居中星号
    ggpubr::stat_pvalue_manual(
      pathway_stat,
      label = "plot_significance",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      hide.ns = hide_ns,
      tip.length = 0.015,
      bracket.size = 0.50,
      size = 5,
      vjust = 0.2,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.25)
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      title = pathway,
      x = NULL,
      y = "Signature score"
    ) +
    boxplot_theme +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        size = 12
      ),
      axis.text.y = element_text(
        size = 11
      ),
      axis.title = element_text(
        size = 12,
        face = "bold"
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5,
        face = "bold"
      )
    )

  selected_plot_list[[pathway]] <- p_pathway
}

# 10.5.4 选定pathway保存为多页PDF(每个pathway一页)
selected_pdf_file <- file.path(
  output_dir,
  "selected_canonical_signature_group_boxplot.pdf"
)

pdf(
  file = selected_pdf_file,
  width = 4,
  height = 6,
  onefile = TRUE,
  useDingbats = FALSE
)

for (pathway in selected_pathways_exist) {
  print(
    selected_plot_list[[pathway]]
  )
}

dev.off()

message(
  "All pathways PDF已保存至: ",
  file.path(
    output_dir,
    "all_canonical_signature_group_boxplot.pdf"
  )
)
message(
  "Selected pathways多页PDF已保存至: ",
  selected_pdf_file
)
message(
  "Selected PDF共包含 ",
  length(selected_pathways_exist),
  " 页, 每个pathway一页。"
)

# ---- 11. Nrg1与module评分相关性分析 -------------------------------------------
# 计算Nrg1基因表达与各module得分的Spearman相关性, 并以热图/柱状图可视化

# Nrg1表达量
nrg1_expr <- GetAssayData(se_obj_traj, layer = "data")["Nrg1", ]

# 全部module得分
module_scores_df <- se_obj_traj@meta.data[, module_scores_name, drop = FALSE]

# 整体Spearman相关性
nrg1_corr_overall <- sapply(module_scores_df, function(score) {
  cor(nrg1_expr, score, method = "spearman", use = "pairwise.complete.obs")
})
nrg1_corr_df_overall <- data.frame(
  Signature = names(nrg1_corr_overall),
  Correlation = nrg1_corr_overall
)
rownames(nrg1_corr_df_overall) <- nrg1_corr_df_overall$Signature
cat("Overall correlation between Nrg1 and module scores calculated\n")

# 按dnb_state分别计算相关性(dnb_state列来自04章DNB分析)
dnb_states <- unique(se_obj_traj$dnb_state)
dnb_states <- dnb_states[dnb_states != "Other"]
nrg1_corr_by_state <- list()
for (state in dnb_states) {
  cells_state <- colnames(se_obj_traj)[se_obj_traj$dnb_state == state]
  if (length(cells_state) < 30) next
  corr_state <- sapply(module_scores_df[cells_state, , drop = FALSE], function(score) {
    cor(nrg1_expr[cells_state], score, method = "spearman", use = "pairwise.complete.obs")
  })
  nrg1_corr_by_state[[state]] <- corr_state
}
nrg1_corr_matrix <- do.call(cbind, nrg1_corr_by_state)
nrg1_corr_matrix <- as.matrix(nrg1_corr_matrix)
colnames(nrg1_corr_matrix) <- names(nrg1_corr_by_state)

# 11.1 相关性热图(按DNB状态分组)
pdf(file.path(DIR_FIG_SIG, "nrg1_module_correlation_heatmap.pdf"),
    width = 8, height = 10)
pheatmap(nrg1_corr_matrix,
         color = colorRampPalette(c("#0765D7", "white", "#EAB404"))(100),
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 8,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Nrg1 Expression vs Module Scores\n(Spearman Correlation by DNB State)",
         border_color = NA,
         breaks = seq(-1, 1, length.out = 101))
dev.off()

# 11.2 整体相关性柱状图(选定通路)
nrg1_corr_df_overall_sorted <- nrg1_corr_df_overall[order(nrg1_corr_df_overall$Correlation), ]
nrg1_corr_df_overall_sorted$Signature <- factor(nrg1_corr_df_overall_sorted$Signature,
                                                levels = nrg1_corr_df_overall_sorted$Signature)
selected_pathways <- c(
  "AT1/AT2",
  "AT2-like",
  "SCC",
  "Glandular",
  "Gastric-like",
  "High-plasticity cell state",
  "Basal",
  "Ciliated",
  "High cycling",
  "Hypoxia",
  "Biosynthetic mixed identity",
  "Stressed",
  "SCC",
  "EMT",
  "pEMT",
  "Stress",
  "High-plasticity cell state"
)
nrg1_corr_df_overall_sorted <- nrg1_corr_df_overall_sorted[nrg1_corr_df_overall_sorted$Signature %in% selected_pathways, ]
p_nrg1_corr_bar <- ggplot(nrg1_corr_df_overall_sorted, aes(x = Signature, y = Correlation)) +
  geom_bar(stat = "identity", fill = "#2E5A8C", alpha = 0.8) +
  labs(title = "Nrg1 Correlation with Module Scores",
       subtitle = "Overall Spearman Correlation",
       x = "", y = "Spearman Correlation Coefficient") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "black"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

pdf(file.path(DIR_FIG_SIG, "nrg1_module_correlation_barplot.pdf"),
    width = 12, height = 8)
print(p_nrg1_corr_bar)
dev.off()

cat("Nrg1 correlation analysis completed\n")

# ---- 12. AST轨迹子集: LUAD/LUSC标记与签名得分 ---------------------------------
# 12.1 tdTomato/GFP两组的LUAD与LUSC标记基因点图
tumor_marker <- list("LUAD" = c("Mlph", "Tmc5", "Sfta3-ps", "Muc5b", "Ttf1", "Napsa", "Ceacam5", "Muc1", "St6galnac1", "Plin2", "Sftpc"),
                     "LUSC" = c("Krt14", "Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt6a", "Trim29", "Slc7a5", "Spats2", "Clca2"))

dot_plot <- DotPlot(se_obj_traj, features = tumor_marker,
                    group.by = "td_group") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_GSVA, "td_group_marker_gene_dotplot.pdf"), width = 6, height = 4)
print(dot_plot)
dev.off()

# 12.2 选定标记基因点图
tumor_marker <- list("LUAD" = c("Sfta3-ps", "Napsa", "Muc1", "Sftpc"),
                     "LUSC" = c("Dsg3", "Trp63", "Krt5", "Krt6a", "Krt14", "Spats2"))

dot_plot <- DotPlot(se_obj_traj, features = tumor_marker,
                    group.by = "td_group") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple",
                        midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_GSVA, "selected_td_group_marker_gene_dotplot.pdf"), width = 6, height = 4)
print(dot_plot)
dev.off()

# 12.3 LUAD/SCC特征基因小提琴图
luad_signature <- c("Mlph", "Tmc5", "Sfta3-ps", "Muc5b", "Ttf1", "Napsa", "Muc1", "St6galnac1", "Plin2", "Sftpc")
scc_signature <- c("Krt14", "Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt6a", "Trim29", "Slc7a5", "Spats2", "Clca2")
pdf(file.path(DIR_FIG_GSVA, "td_group_AST_related_genes.pdf"), width = 8, height = 8)
VlnPlot(se_obj_traj, features = luad_signature, group.by = "td_group")
VlnPlot(se_obj_traj, features = scc_signature, group.by = "td_group")
dev.off()

# 过滤不存在的基因
luad_signature <- luad_signature[luad_signature %in% rownames(se_obj_traj)]
scc_signature <- scc_signature[scc_signature %in% rownames(se_obj_traj)]

# 12.4 每个细胞的平均表达得分(RNA assay)
DefaultAssay(se_obj_traj) <- "RNA"
se_obj_traj$LUAD_score <- colMeans(GetAssayData(se_obj_traj)[luad_signature, ])
se_obj_traj$SCC_score <- colMeans(GetAssayData(se_obj_traj)[scc_signature, ])
p1 <- VlnPlot(se_obj_traj, features = "LUAD_score", group.by = "td_group",
              cols = c("#6EB542", "pink"),
              pt.size = 0) + NoLegend() + ggtitle("LUAD signature score")
p2 <- VlnPlot(se_obj_traj, features = "SCC_score", group.by = "td_group",
              cols = c("#6EB542", "pink"),
              pt.size = 0) + NoLegend() + ggtitle("SCC signature score")

pdf(file.path(DIR_FIG_GSVA, "signature_score_distribution.pdf"), width = 8, height = 10)
print(p1 + p2)
dev.off()

# 12.5 带P值的得分分布图(Wilcoxon检验)
# 绘图时坐标轴标签转换为数学表达式(上标+)
label_map <- c(
  "GFP" = expression(ZsGreen^"+"),
  "tdTomato" = expression(tdTomato^"+")
)

# 两组比较
comp <- list(c("GFP", "tdTomato"))

p1 <- VlnPlot(se_obj_traj,
              features = "LUAD_score",
              group.by = "td_group",
              cols = c("#6EB542", "pink"),
              pt.size = 0) +
  NoLegend() +
  ggtitle("LUAD signature score") +
  geom_signif(
    comparisons = comp,
    test = wilcox.test,       # Wilcoxon秩和检验
    map_signif_level = TRUE,  # 自动转为 */**/***
    y_position = 2,           # 星号高度, 根据分值上下调整
    textsize = 4
  ) + ylim(c(min(se_obj_traj$LUAD_score), 2.5)) +
  scale_x_discrete(labels = label_map)

p2 <- VlnPlot(se_obj_traj,
              features = "SCC_score",
              group.by = "td_group",
              cols = c("#6EB542", "pink"),
              pt.size = 0) +
  NoLegend() +
  ggtitle("SCC signature score") +
  geom_signif(
    comparisons = comp,
    test = wilcox.test,
    map_signif_level = TRUE,
    y_position = 1,
    textsize = 4
  ) + ylim(c(min(se_obj_traj$SCC_score), 1.2)) +
  scale_x_discrete(labels = label_map)

pdf(file.path(DIR_FIG_GSVA, "signature_score_distribution_with_p.pdf"), width = 4, height = 6)
print(p1)
print(p2)
dev.off()

# 12.6 得分UMAP分布图
pdf(file = file.path(DIR_FIG_GSVA, "umap_signature_sore.pdf"), width = 6, height = 5)
FeatureDimPlot(se_obj_traj, features = "LUAD_score", theme_use = "theme_blank", legend.position = "bottom")
FeatureDimPlot(se_obj_traj, features = "SCC_score", theme_use = "theme_blank", legend.position = "bottom")
dev.off()
