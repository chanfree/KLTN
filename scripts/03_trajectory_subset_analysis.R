# ============================================================================
# 脚本名称: 03_trajectory_subset_analysis.R
# 描述: AST(腺鳞转化)轨迹亚群核心分析：tdTomato_ratio相关时序基因、slingshot/tradeSeq
#       伪时间轨迹、TCGA生存验证、SCENIC调控网络、LUAD/LUSC signature评分, 以及
#       top100关联基因聚类/GO富集/调控网络(FigS4)
# 输入: output/object/seurat/se.rds
#       output/object/seurat/se_trajectory_subset.rds
#           (来源说明: 为 se.rds 中 seurat_clusters 属于 0/1/2/3 的上皮细胞子集)
#       output/object/cellrank/C1_drivers.csv, C2_drivers.csv (cellrank.ipynb 产物)
#       output/object/slingshot/sce_full.rds (tradeSeq GAM 拟合结果)
#       output/object/survival/TCGA-LUAD_vst_counts_with_symbols.rds
#       output/object/survival/TCGA-LUAD_clinical_data.rds
#       output/object/SCENIC_trajectory/ (int/scenicOptions.Rds,
#       int/3.5_AUCellThresholds.Rds, output/scenic.loom)
#       data/geneset.hs.gmt
#       ~/database/SCENIC_ref (外部SCENIC参考数据库, 需自行准备)
#       注: FigS4 B 的 top100 关联基因 GeneMANIA/Cytoscape 网络复现脚本见
#           top100_associate_cytoscape/ (独立于本 R 脚本)
# 输出: output/figures/td_ratio/, output/figures/slingshot/,
#       output/figures/survival/, output/figures/SCENIC_trajectory/,
#       output/figures/signature_score/, output/figures/gsva/ 下图表
#       (td_ratio/ 下含 FigS4 相关 heatmap_genes_cluster_*_go.pdf 与
#        heatmap_genes_cluster_goBP_top10_barplot.pdf)
#       output/object/ 对应子目录下RDS/CSV对象
# 上游脚本: 01_KLTN_process.R (se.rds等)
# ============================================================================

# ---- 0. 环境与路径配置 ----

set.seed(1)

DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "object")
fig_res <- 300

# 各分析子目录(沿用原输出目录结构)
DIR_FIG_td_ratio <- file.path(DIR_FIG, "td_ratio")
DIR_OBJ_td_ratio <- file.path(DIR_OBJ, "td_ratio")
DIR_FIG_slingshot <- file.path(DIR_FIG, "slingshot")
DIR_OBJ_slingshot <- file.path(DIR_OBJ, "slingshot")
DIR_FIG_survival <- file.path(DIR_FIG, "survival")
DIR_OBJ_survival <- file.path(DIR_OBJ, "survival")
DIR_FIG_SCENIC <- file.path(DIR_FIG, "SCENIC_trajectory")
DIR_OBJ_SCENIC <- file.path(DIR_OBJ, "SCENIC_trajectory")
DIR_FIG_sig <- file.path(DIR_FIG, "signature_score")
DIR_OBJ_sig <- file.path(DIR_OBJ, "signature_score")
DIR_FIG_gsva <- file.path(DIR_FIG, "gsva")
DIR_OBJ_gsva <- file.path(DIR_OBJ, "gsva")

for (d in c(DIR_FIG_td_ratio, DIR_OBJ_td_ratio, DIR_FIG_slingshot, DIR_OBJ_slingshot,
            DIR_FIG_survival, DIR_OBJ_survival, DIR_FIG_SCENIC, DIR_OBJ_SCENIC,
            DIR_FIG_sig, DIR_OBJ_sig, DIR_FIG_gsva, DIR_OBJ_gsva)) {
  dir.create(path = d, recursive = TRUE, showWarnings = FALSE)
}

# 逻辑开关(重跑控制) ---------------------------------------------------------
# TCGA-LUAD 数据下载与预处理: 已运行过, 结果已保存至 output/object/survival/;
# 首次运行或需重新下载数据时设为 TRUE
RUN_TCGA_DOWNLOAD <- FALSE

# 加载必要的包
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(future)
library(tradeSeq)
library(clusterExperiment)
library(pheatmap)
library(tibble)
library(stringr)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(orthogene)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(shiny)
library(patchwork)
library(ggrastr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggpubr)
library(rstatix)
library(tidyr)
library(SCP)

plan(multisession, workers = 8)
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)
bpparam <- BiocParallel::bpparam()
bpparam$workers <- 10

source(file.path(DIR_ROOT, "R", "custom_functions.R"))
source(file.path(DIR_ROOT, "R", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

# ---- 1. 读取数据与 tdTomato_ratio 可视化 ----

# 1.1 从完整对象提取上皮细胞(去除免疫细胞等非上皮群体)
seu_object <- readRDS(file.path(DIR_OBJ, "seurat", "se.rds"))
seu_object <- subset(seu_object, first_annotation == "epithelium")

# 1.2 tdTomato_ratio 分布直方图(叠加密度曲线与中位数虚线)
meta_dat <- seu_object@meta.data
data_median <- median(meta_dat$tdTomato_ratio)
pdf(file.path(DIR_FIG_td_ratio, "td_ratio_histo_distribution.pdf"), width = 6, height = 4)
ggplot(meta_dat, aes(x = tdTomato_ratio)) +
  geom_histogram(
    aes(y = after_stat(density)),  # y轴使用密度以匹配密度曲线
    bins = 30,
    color = "black",
    fill = "lightblue"
  ) +
  geom_density(
    color = "red",
    linewidth = 1
  ) +
  geom_vline(
    aes(xintercept = data_median),  # 中位数处虚线
    color = "black",
    linewidth = 1.5,
    show.legend = TRUE,
    linetype = "dashed"
  ) +
  labs(
    title = "Histogram with Density Curve",
    x = "tdTomato_ratio",
    y = "Density"
  ) +
  theme_minimal()
dev.off()

# 1.3 tdTomato_ratio 在UMAP上的空间分布(0.0=LUAD, 1.0=LUSC)
x_limits <- range(seu_object@reductions$umap@cell.embeddings[, 1], na.rm = TRUE)
y_limits <- range(seu_object@reductions$umap@cell.embeddings[, 2], na.rm = TRUE)

rel_order_umap <- FeaturePlot(seu_object, "tdTomato_ratio") +
  scale_colour_gradientn(colours = c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF",
                                     "#FE9F6DFF", "#FCFDBFFF"),
                         na.value = "transparent", limits = c(0, 1),
                         breaks = seq(0, 1, by = 0.2),
                         labels = c("0.0 (LUAD)", "0.2", "0.4", "0.6", "0.8", "1.0 (LUSC)"),
                         name = "tdTomato UMI Ratio",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("AST level") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                  margin = margin(b = 20)),
        aspect.ratio = 1) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
pdf(file = file.path(DIR_FIG_td_ratio, "td_ratio_umap.pdf"), width = 6, height = 5)
print(rel_order_umap)
dev.off()

# ---- 2. 沿AST过程的差异表达基因识别(Spearman相关) ----

# 检测与tdTomato_ratio显著相关的时序基因: 逐基因Spearman相关 + BH多重校正 +
# 筛选显著基因(p_val_adj<0.05且|correlation|>0.3) + 火山图/热图可视化;
# 全基因相关结果输出 correlation_all_genes_epi.csv
detect_and_visualize_temporal_genes_R <- function(seurat_obj, time_key = "time",
                                                  p_val_adj_cutoff = 0.05,
                                                  correlation_cutoff = 0.3) {
  # 检查time_key是否存在于metadata
  if (!time_key %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Time key '", time_key, "' not found in Seurat object metadata.", sep = ""))
  }

  # 提取归一化表达矩阵与时间向量
  cat("--- Starting Temporal Gene Analysis ---\n")
  expression_matrix <- GetAssayData(seurat_obj, slot = "data")
  expression_matrix <- as.matrix(expression_matrix)
  time_vector <- seurat_obj@meta.data[[time_key]]

  # 逐基因计算与时间向量的Spearman相关
  cat("Step 1/4: Calculating Spearman correlation for each gene...\n")
  correlation_results <- apply(expression_matrix, 1, function(gene_expression) {
    if (sd(gene_expression) > 0) {
      test_result <- cor.test(gene_expression, time_vector, method = "spearman", exact = FALSE)
      return(c(correlation = test_result$estimate, p_value = test_result$p.value))
    } else {
      return(c(correlation = 0, p_value = 1))
    }
  })

  results_df <- as.data.frame(t(correlation_results))
  colnames(results_df) <- c("correlation", "p_value")

  # BH校正p值
  cat("Step 2/4: Correcting for multiple testing...\n")
  results_df$p_val_adj <- p.adjust(results_df$p_value, method = "BH")

  cat("Step 3/4: Identifying significant genes...\n")
  results_df <- results_df %>%
    rownames_to_column(var = "gene") %>%
    na.omit()
  write.csv(results_df, file.path(DIR_OBJ_td_ratio, "correlation_all_genes_epi.csv"))
  significant_genes <- results_df %>%
    filter(p_val_adj < p_val_adj_cutoff, abs(correlation) > correlation_cutoff) %>%
    arrange(desc(correlation))

  # 火山图: 上调(红)/下调(蓝)/不显著(灰)
  results_df <- results_df %>%
    mutate(significance = ifelse(gene %in% significant_genes$gene,
                                 ifelse(correlation > 0, "Upregulated", "Downregulated"),
                                 "Not Significant"))

  volcano_plot <- ggplot(results_df, aes(x = correlation, y = -log10(p_val_adj))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue",
                                  "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(p_val_adj_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-correlation_cutoff, correlation_cutoff), linetype = "dashed") +
    labs(
      title = "Volcano Plot of Temporal Gene Expression",
      x = "Spearman Correlation with Time",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())

  print(volcano_plot)

  # 热图: 显著基因按时序排序的表达(行缩放)
  cat("\nGenerating heatmap...\n")
  heatmap_data <- GetAssayData(seurat_obj, slot = "data")[significant_genes$gene[10], ] |>
    as.matrix()
  cell_order <- order(seurat_obj@meta.data[[time_key]])
  heatmap_data_ordered <- heatmap_data[, cell_order]
  heatmap_data_scaled <- t(scale(t(heatmap_data_ordered)))

  pheatmap(
    heatmap_data_scaled,
    cluster_cols = FALSE,  # 列按时序排列, 不聚类
    cluster_rows = TRUE,
    show_colnames = FALSE,
    show_rownames = nrow(significant_genes) <= 50,
    main = "Expression of Significant Genes Over Time",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  )

  cat("\n--- Analysis Complete ---\n")
  return(significant_genes)
}

significant_results_r <- detect_and_visualize_temporal_genes_R(
  seurat_obj = seu_object,
  time_key = "tdTomato_ratio"
)

# ---- 3. slingshot 轨迹推断 ----

# 3.1 读取AST轨迹亚群对象
seu_object <- readRDS(file.path(DIR_OBJ, "seurat", "se_trajectory_subset.rds"))

# 3.2 各细胞类型经典marker基因DotPlot(验证亚群注释)
markers <- list(
  "AT1" = c("Ager", "Hopx", "Clic5"),
  "AT2" = c("Sftpc", "Sftpa1", "Sftpd", "Etv5"),
  "Basal" = c("Krt15", "Krt5", "Krt17", "Trp63"),
  "Ciliated" = c("Pifo", "Foxj1", "Ccdc78"),
  "Club" = c("Scgb3a2", "Scgb1a1", "Cyp2f2"),
  "Goblet" = c("Muc2", "Tff2", "Muc5b", "Spdef"),
  "KACs" = c("Cldn4", "Nfkb2", "Krt8", "Cdkn1a"),
  "Proliferative" = c("Birc5", "Top2a", "Mki67"),
  "endothelial" = c("Pecam1", "Cdh5", "Cldn5", ""),
  "Fibroblast" = c("Col1a2", "Col1a1", "Vim", "Pdgfra")
)
dot_plot <- DotPlot(seu_object, features = markers, group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple", midpoint = 0) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_slingshot, "cannonical_markers_dotplot.pdf"), width = 20, height = 6)
print(dot_plot)
dev.off()

# 3.3 CytoTRACE2 相对分化程度在UMAP上的分布
x_limits <- range(seu_object@reductions$umap@cell.embeddings[, 1], na.rm = TRUE)
y_limits <- range(seu_object@reductions$umap@cell.embeddings[, 2], na.rm = TRUE)
rel_order_umap <- FeaturePlot(seu_object, "CytoTRACE2_Relative") +
  scale_colour_gradientn(colours = c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF",
                                     "#FE9F6DFF", "#FCFDBFFF"),
                         na.value = "transparent", limits = c(0, 1),
                         breaks = seq(0, 1, by = 0.2),
                         labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8",
                                    "1.0 (Less diff.)"),
                         name = "Relative\norder \n",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("CytoTRACE 2") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                  margin = margin(b = 20)),
        aspect.ratio = 1) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
pdf(file = file.path(DIR_FIG_slingshot, "cytorace_relative_umap.pdf"), width = 6, height = 5)
print(rel_order_umap)
dev.off()

# 3.4 提取降维数据与聚类标签并转换为SingleCellExperiment
reduction_method <- "umap"
if (!(reduction_method %in% names(seu_object@reductions))) {
  stop(paste("Reduction method", reduction_method, "not found in Seurat object."))
}
reduced_dims <- Embeddings(seu_object, reduction = reduction_method)

cluster_column_name <- "seurat_clusters"
if (!(cluster_column_name %in% colnames(seu_object@meta.data))) {
  stop(paste("Cluster column", cluster_column_name, "not found in Seurat object metadata."))
}
cluster_labels <- Idents(seu_object)
if (!is.factor(cluster_labels)) cluster_labels <- factor(cluster_labels)

sce <- as.SingleCellExperiment(seu_object, assay = "RNA")
reducedDim(sce, toupper(reduction_method)) <- reduced_dims
colData(sce)$Cluster <- cluster_labels
sce_reduction_name <- toupper(reduction_method)
reducedDim(sce, sce_reduction_name) <- reduced_dims
colData(sce)$Cluster <- cluster_labels

# 3.5 运行slingshot(指定起始cluster为"0", 即AT2/LUAD端)
start_cluster <- "0"
if (!start_cluster %in% levels(colData(sce)$Cluster)) {
  stop(paste("Start cluster '", start_cluster, "' not found in the cluster labels.",
             "Available clusters are:",
             paste(levels(colData(sce)$Cluster), collapse = ", "), sep = ""))
}

tryCatch({
  sce_slingshot <- slingshot(sce,
                             clusterLabels = "Cluster",
                             reducedDim = sce_reduction_name,
                             start.clus = start_cluster,
                             approx_points = 150,
                             allow.breaks = FALSE)
  message("Slingshot analysis completed successfully.")
}, error = function(e) {
  message("Error during Slingshot analysis: ", e$message)
  message("Common issues to check:")
  message("- Is the 'start.clus' correctly specified and present in your 'Cluster' labels?")
  message("- Are the 'clusterLabels' and 'reducedDim' arguments pointing to the correct data in the SCE object?")
  message("- Does your data form a connected structure suitable for trajectory inference?")
})

# 3.6 轨迹可视化: ggplot版(细胞按cluster着色+轨迹曲线)
num_clusters <- length(levels(colData(sce_slingshot)$Cluster))
cluster_colors <- scales::hue_pal()(num_clusters)
names(cluster_colors) <- levels(colData(sce_slingshot)$Cluster)

slingshot_df <- as.data.frame(reducedDim(sce_slingshot, sce_reduction_name))
colnames(slingshot_df) <- c("Dim1", "Dim2")
slingshot_df$Cluster <- colData(sce_slingshot)$Cluster

sling_curves_list <- slingCurves(sce_slingshot)
curve_plot_df_list <- list()
for (i in seq_along(sling_curves_list)) {
  curve_df <- as.data.frame(sling_curves_list[[i]]$s[sling_curves_list[[i]]$ord, 1:2])
  colnames(curve_df) <- c("Dim1", "Dim2")
  curve_df$Lineage <- paste0("Lineage", i)
  curve_plot_df_list[[i]] <- curve_df
}
curve_plot_df <- bind_rows(curve_plot_df_list)

gg_traj <- ggplot(slingshot_df, aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = cluster_colors) +
  geom_path(data = curve_plot_df, aes(x = Dim1, y = Dim2, group = Lineage),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  labs(x = paste0(sce_reduction_name, "_1"), y = paste0(sce_reduction_name, "_2"),
       title = paste("Slingshot Trajectory on Seurat", sce_reduction_name)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", aspect.ratio = 1)
pdf(file.path(DIR_FIG_slingshot, "umap_slingshot_trajectory.pdf"), width = 6, height = 5)
print(gg_traj)
dev.off()

# 3.7 将各lineage的slingshot伪时间回填到Seurat对象metadata
pseudotime_slingshot <- slingPseudotime(sce_slingshot)
for (i in 1:ncol(pseudotime_slingshot)) {
  seu_object[[paste0("slingPseudotime_L", i)]] <- NA
  cells_in_lineage <- rownames(pseudotime_slingshot)[!is.na(pseudotime_slingshot[, i])]
  seu_object[[paste0("slingPseudotime_L", i)]][cells_in_lineage, ] <-
    pseudotime_slingshot[cells_in_lineage, i]
}

# 伪时间在UMAP上的分布(Lineage 1 / Lineage 2)
pdf(file.path(DIR_FIG_slingshot, "umap_slingshot_Pseudotime_Lineage1.pdf"), width = 6, height = 5)
print(FeaturePlot(seu_object, features = "slingPseudotime_L1", reduction = reduction_method) +
        labs(title = "Slingshot Pseudotime (Lineage 1)")) +
  scale_color_viridis_c(option = "mako")
dev.off()
pdf(file.path(DIR_FIG_slingshot, "umap_slingshot_Pseudotime_Lineage2.pdf"), width = 6, height = 5)
print(FeaturePlot(seu_object, features = "slingPseudotime_L2", reduction = reduction_method) +
        labs(title = "Slingshot Pseudotime (Lineage 2)")) +
  scale_color_viridis_c(option = "mako")
dev.off()

# ---- 4. tradeSeq 差异表达分析(slingshot轨迹) ----
# 注: cellrank 分析(cellrank drivers)见 scripts/cellrank.ipynb(Python, 已纳入本仓库),
#     其产出 C1/C2_drivers.csv 为本脚本第 5 章(TCGA 生存分析)的输入

# 4.1 提取计数矩阵
counts <- as.matrix(seu_object@assays$RNA@layers$counts)
rownames(counts) <- rownames(seu_object)
colnames(counts) <- colnames(seu_object)

# 4.2 构建lineage与主曲线(起始cluster "0")
lineages <- getLineages(data = reduced_dims, clusterLabels = cluster_labels, start.clus = "0")
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8,
                    allow.breaks = FALSE, shrink = 0.99)
saveRDS(curves, file.path(DIR_OBJ_slingshot, "curves.rds"))

# GAM拟合(耗时步骤): 结果已保存至 sce_full.rds, 重跑时取消下方注释
# sce <- fitGAM(counts = counts, sds = curves, parallel = TRUE, BPPARAM = bpparam)
# saveRDS(sce, file.path(DIR_OBJ_slingshot, "sce_full.rds"))
sce <- readRDS(file.path(DIR_OBJ_slingshot, "sce_full.rds"))

# 差异基因可视化函数: 轨迹上表达 + smoother拟合 并排展示
plot_differential_expression <- function(feature_id) {
  p1 <- plotGeneCount(curves, counts, gene = feature_id[1], clusters = cluster_labels,
                      models = sce) +
    ggplot2::theme(legend.position = "top")
  p2 <- plotSmoothers(sce, as.matrix(counts), gene = feature_id[1])
  print(p1 + p2)
}

# 4.3 基因表达与伪时间的关联检验(各lineage终点对比)
assoRes <- associationTest(sce, lineages = TRUE, contrastType = "end")
write.csv(assoRes, file.path(DIR_OBJ_slingshot, "assoRes.csv"))

# C1/C2 lineage各自top10关联基因的轨迹表达图
top10_c1 <- assoRes |>
  filter(pvalue_2 < 0.05) |>
  top_n(10, wt = waldStat_2) |>
  row.names()
plots <- lapply(top10_c1, function(i) plotGeneCount(curves, counts, gene = i))
pdf(file.path(DIR_FIG_slingshot, "C1_top10_associated_gene.pdf"), width = 30, height = 6)
grid.arrange(grobs = plots, ncol = 5, nrow = 2)
dev.off()

top10_c2 <- assoRes |>
  filter(pvalue_1 < 0.05) |>
  top_n(10, wt = waldStat_1) |>
  row.names()
plots <- lapply(top10_c2, function(i) plotGeneCount(curves, counts, gene = i))
pdf(file.path(DIR_FIG_slingshot, "C2_top10_associated_gene.pdf"), width = 30, height = 6)
grid.arrange(grobs = plots, ncol = 5, nrow = 2)
dev.off()

# 4.4 起始端(祖细胞)marker基因识别
startRes <- startVsEndTest(sce, lineages = TRUE)
startRes_filter <- startRes |> filter(pvalue < 0.05)
write.csv(startRes_filter, file.path(DIR_OBJ_slingshot, "startRes_filter.csv"))
oStart <- order(startRes_filter$waldStat_lineage2, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[1:10]]
plots <- lapply(sigGeneStart, function(i) plotGeneCount(curves, counts, gene = i))
pdf(file.path(DIR_FIG_slingshot, "C1_top10_starting_gene_umap.pdf"), width = 30, height = 6)
grid.arrange(grobs = plots, ncol = 5, nrow = 2)
dev.off()
plots <- lapply(sigGeneStart, function(i) plotSmoothers(sce, counts, gene = i))
pdf(file.path(DIR_FIG_slingshot, "C1_top10_starting_gene_smoother.pdf"), width = 20, height = 6)
grid.arrange(grobs = plots, ncol = 5, nrow = 2)
dev.off()

# 比较lineage内特定伪时间点(0 vs 15)的差异表达
customRes <- startVsEndTest(sce, pseudotimeValues = c(0, 15))
customRes_filter <- customRes |> filter(pvalue < 0.05)

# 4.5 不同lineage间表达模式差异基因(pattern test)
patternRes <- patternTest(sce)
patternRes_filter <- patternRes |> filter(pvalue < 0.05)
oPat <- order(patternRes_filter$waldStat, decreasing = TRUE)
write.csv(patternRes_filter, file.path(DIR_OBJ_slingshot, "patternRes_filter.csv"))
# 回读CSV以将行名恢复为基因符号(供后续按基因名取图)
patternRes_filter <- read.csv(file.path(DIR_OBJ_slingshot, "patternRes_filter.csv"),
                              row.names = 1)

# 4.6 终末分化细胞类型marker(diffEnd test)
endRes <- diffEndTest(sce)
endRes_filter <- endRes |> filter(pvalue < 0.05)
write.csv(endRes_filter, file.path(DIR_OBJ_slingshot, "endRes_filter.csv"))
o <- order(endRes_filter$waldStat, decreasing = TRUE)
sigGene <- rownames(endRes_filter)[o[1]]

# 4.7 分化早期驱动基因(相邻knots间的earlyDE test)
for (i in 1:3) {
  earlyDERes <- earlyDETest(sce, knots = c(i, i + 1))
  oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
  plot_genes <- head(rownames(earlyDERes)[oEarly])
  pdf(file = file.path(DIR_FIG_slingshot,
                       sprintf("early_different_expressed_gene_knot_%s_%s.pdf", i, i + 1)),
      width = 8, height = 4)
  for (n in plot_genes) {
    plot_differential_expression(n)
  }
  dev.off()
  write.csv(earlyDERes,
            file.path(DIR_OBJ_slingshot,
                      sprintf("earlyDERes_knot_%s_%s.csv", i, i + 1)))
}

# 4.8 AST相关基因的轨迹表达图(LUAD端与LUSC端marker)
pdf(file = file.path(DIR_FIG_slingshot, "AST_related_gene.pdf"), width = 8, height = 4)
plot_differential_expression("Sox2")
plot_differential_expression("Trp63")
plot_differential_expression("Ttf1")
plot_differential_expression("Napsa")
plot_differential_expression("Ager")
plot_differential_expression("Sftpc")
dev.off()

# 4.9 KAC(角质化腺细胞)相关基因的轨迹表达图
kac_genes <- c("Krt8", "Cldn4", "Nfkb2", "Cdkn1a")
pdf(file = file.path(DIR_FIG_slingshot, "KAC_related_gene.pdf"), width = 8, height = 4)
for (i in kac_genes) {
  plot_differential_expression(i)
}
dev.off()

# ---- 5. TCGA 生存分析(验证轨迹相关基因 signature) ----

# 5.1 读取cellrank drivers(来自 cellrank.ipynb), 取各lineage top50驱动基因作为签名
c1_driver <- read.csv(file.path(DIR_OBJ, "cellrank", "C1_drivers.csv"), header = TRUE,
                      row.names = 1)
c1_driver %>%
  filter(C1_3_qval < 0.05) %>%
  slice_max(order_by = C1_3_corr, n = 50) %>%
  rownames() -> c1_top50
mouse_gene_signature <- c1_top50

# c2_top50: CellRank lineage 2 (C2) 谱系驱动基因, 按 C2_corr(与 C2 命运概率的相关性)
# 降序取 top50, 并过滤 C2_qval < 0.05; 用于 Fig S2U 的 TCGA-LUAD KM 生存分层
c2_driver <- read.csv(file.path(DIR_OBJ, "cellrank", "C2_drivers.csv"), header = TRUE,
                      row.names = 1)
c2_driver %>%
  filter(C2_qval < 0.05) %>%
  slice_max(order_by = C2_corr, n = 50) %>%
  rownames() -> c2_top50

tcga_project <- "TCGA-LUAD"

# 5.2 TCGA数据下载与预处理(由 RUN_TCGA_DOWNLOAD 控制; 已运行过, 结果已保存)
if (RUN_TCGA_DOWNLOAD) {
  query <- GDCquery(
    project = tcga_project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # GDCdownload在当前工作目录创建GDCdata文件夹, 需切换到保存目录
  setwd(DIR_OBJ_survival)
  GDCdownload(query, method = "api", files.per.chunk = 10)

  data <- GDCprepare(query)

  clinical_data <- as.data.frame(colData(data))
  rna_counts <- assay(data, "unstranded")

  # 构建生存分析变量: status(Dead=1)与time_to_event(优先days_to_death)
  clinical_data$vital_status_numeric <- ifelse(clinical_data$vital_status == "Dead", 1, 0)
  clinical_data$time_to_event <- ifelse(
    !is.na(clinical_data$days_to_death),
    clinical_data$days_to_death,
    clinical_data$days_to_last_follow_up
  )
  clinical_data <- clinical_data[!is.na(clinical_data$time_to_event) &
                                   clinical_data$time_to_event > 0, ]

  rna_counts_filtered <- rna_counts[, rownames(clinical_data)]

  # VST标准化
  vst_counts <- DESeq2::varianceStabilizingTransformation(round(rna_counts_filtered),
                                                          blind = TRUE)

  # Ensembl ID转换为基因符号并去除重复
  gene_mapping <- as.data.frame(rowData(data))
  vst_counts_with_symbols <- merge(vst_counts, gene_mapping[, c("gene_id", "gene_name")],
                                   by.x = 0, by.y = "gene_id")
  vst_counts_with_symbols <- vst_counts_with_symbols[!duplicated(vst_counts_with_symbols$gene_name), ]
  rownames(vst_counts_with_symbols) <- vst_counts_with_symbols$gene_name
  vst_counts_with_symbols <- vst_counts_with_symbols[, -c(1, ncol(vst_counts_with_symbols))]
  saveRDS(vst_counts_with_symbols,
          file.path(DIR_OBJ_survival, "LUSC_vst_counts_with_symbols.rds"))
  saveRDS(clinical_data, file.path(DIR_OBJ_survival, "LUSC_clinical_data.rds"))
  setwd(DIR_ROOT)
}

# 5.3 读取已预处理的TCGA表达与临床数据
vst_counts_with_symbols <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_vst_counts_with_symbols.rds", tcga_project)))
clinical_data <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_clinical_data.rds", tcga_project)))

# 5.4 KM生存曲线绘制函数: 小鼠签名基因同源转换 -> 患者均分signature得分 ->
#     中位分组 -> Kaplan-Meier曲线(log-rank检验)
draw_KM_plot <- function(mouse_gene_signature, vst_counts_with_symbols, clinical_data) {
  require(orthogene)
  gene_signature <- convert_orthologs(
    gene_df = mouse_gene_signature,
    input_species = "mouse",
    output_species = "human",
    method = "homologene"
  )
  gene_signature <- unique(rownames(gene_signature))
  signature_expr <- vst_counts_with_symbols[rownames(vst_counts_with_symbols) %in%
                                              gene_signature, ]

  missing_genes <- setdiff(gene_signature, rownames(signature_expr))
  if (length(missing_genes) > 0) {
    warning(paste("The following genes were not found in the dataset:",
                  paste(missing_genes, collapse = ", ")))
  }

  # 每个患者的signature得分 = 签名基因平均表达
  signature_expr_t <- as.data.frame(t(signature_expr))
  signature_expr_t$signature_score <- rowMeans(signature_expr_t, na.rm = TRUE)

  survival_data <- merge(clinical_data, signature_expr_t[, "signature_score", drop = FALSE],
                         by.x = 0, by.y = 0)
  colnames(survival_data)[1] <- "patient_barcode"

  # 按中位得分分为高/低风险组
  median_score <- median(survival_data$signature_score)
  survival_data$risk_group <- ifelse(survival_data$signature_score >= median_score,
                                     "High_Risk", "Low_Risk")
  survival_data$risk_group <- factor(survival_data$risk_group,
                                     levels = c("Low_Risk", "High_Risk"))

  # 时间单位转换为月(平均30.44天/月)
  survival_data$time_to_event_months <- survival_data$time_to_event / 30.44

  fit <- survfit(Surv(time = time_to_event_months, event = vital_status_numeric) ~ risk_group,
                 data = survival_data)

  km_plot <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    risk.table.col = "strata",
    legend.labs = c("Low Risk", "High Risk"),
    legend.title = "Risk Group",
    palette = c("#0073C2FF", "#E7B800FF"),
    title = paste("Kaplan-Meier Curve for", tcga_project, "by Gene Signature"),
    xlab = "Time (Months)",
    ggtheme = theme_classic()
  )

  return(km_plot)
}

c1_plot <- draw_KM_plot(c1_top50, vst_counts_with_symbols = vst_counts_with_symbols,
                        clinical_data = clinical_data)
c2_plot <- draw_KM_plot(c2_top50, vst_counts_with_symbols = vst_counts_with_symbols,
                        clinical_data = clinical_data)
pdf(file.path(DIR_FIG_survival, "cellrank_c2_top50_TCGA_LUAD_KM.pdf"), width = 8, height = 5)
print(c2_plot)
dev.off()

# ---- 6. tradeSeq GAM 分析(以tdTomato_ratio为伪时间) ----

# 6.1 读取轨迹亚群对象并转换为SingleCellExperiment
seurat_obj <- readRDS(file.path(DIR_OBJ, "seurat", "se_trajectory_subset.rds"))
sce <- as.SingleCellExperiment(seurat_obj)

counts_matrix <- as.matrix(counts(sce))

# 以tdTomato_ratio为单轨迹伪时间, 细胞权重均为1
cell_pseudotime <- data.frame(curve1 = sce$tdTomato_ratio)
rownames(cell_pseudotime) <- colnames(sce)
cell_weights <- matrix(1, nrow = ncol(sce), ncol = 1)
rownames(cell_weights) <- colnames(sce)

cat("\nData prepared for tradeSeq analysis.\n")

# 6.2 拟合GAM模型(nknots=6), 将每个基因表达建模为伪时间的光滑函数
cat("\nFitting GAMs for each gene. This may take a while on large datasets...\n")
sce_tradeseq <- fitGAM(
  counts = as.matrix(counts(sce)),
  pseudotime = cell_pseudotime,
  cellWeights = cell_weights,
  nknots = 6,
  verbose = FALSE, parallel = TRUE, BPPARAM = bpparam
)
saveRDS(sce_tradeseq, file.path(DIR_OBJ_td_ratio, "sce_subset_tradeseq.rds"))
cat("GAM fitting complete.\n")

# 如需跳过耗时的GAM拟合, 可注释上方fitGAM并启用下述readRDS(与保存对象等价):
# sce_tradeseq <- readRDS(file.path(DIR_OBJ_td_ratio, "sce_subset_tradeseq.rds"))

# 6.3 关联检验: 识别与tdTomato_ratio显著相关的基因
cat("\nRunning association test to find significant genes...\n")
association_results <- associationTest(sce_tradeseq)
association_results$gene <- rownames(association_results)
write.csv(association_results, file.path(DIR_OBJ_td_ratio, "association_results.csv"))
# 回读CSV以将行名恢复为基因符号(供后续按基因名取图/取子集)
association_results <- read.csv(file.path(DIR_OBJ_td_ratio, "association_results.csv"),
                                row.names = 1)

# 6.4 AST相关基因的smoother拟合曲线图
genes_plot <- c("Sox2", "Trp63", "Ttf1", "Napsa", "Ager", "Sftpc")
pdf(file = file.path(DIR_FIG_td_ratio, "AST_related_gene.pdf"), width = 4, height = 4)
for (i in genes_plot) {
  p <- plotSmoothers(sce_tradeseq, counts_matrix, gene = i) + ggtitle(i)
  print(p)
}
dev.off()

# 6.5 top100显著基因的动态表达热图(predictSmooth预测100个伪时间点)
association_results |>
  filter(pvalue < 0.05) |>
  top_n(100, waldStat) -> top_100
significant_genes <- rownames(top_100)

smoothed_expression_matrix <- predictSmooth(sce_tradeseq,
                                            gene = significant_genes, nPoints = 100,
                                            tidy = FALSE)
P <- pheatmap(smoothed_expression_matrix, scale = "row", color = col_scale_cellrank,
              cluster_cols = FALSE,  # 列保持伪时间顺序
              cluster_rows = TRUE,
              show_colnames = FALSE,
              main = "Dynamic Genes Along Time",
              legend_labels = "Predict Expression",
              border_color = NA)
pdf(file.path(DIR_FIG_td_ratio, "top100_gene_heatmap.pdf"), width = 8, height = 15)
print(P)
dev.off()

# 6.6 top100动态基因热图聚类与逐类GO富集(FigS4 C/D)
#    将top100动态基因热图的行聚类划分为2类(沿tdRatio表达上升/下降),
#    分别进行GO过表达富集并绘制top10 BP条目条形图; clusters 1/2
#    分别对应上调/下调动态基因
row_clusters <- cutree(P$tree_row, h = 15)
table(row_clusters)
cluster_1 <- names(row_clusters[row_clusters == 1])  # 沿tdRatio表达上升(上调)
cluster_2 <- names(row_clusters[row_clusters == 2])  # 沿tdRatio表达下降(下调)

go_top_list <- list()
for (i in 1:2) {
  mygene <- names(row_clusters[row_clusters == i])
  write.csv(mygene, file = file.path(DIR_OBJ_td_ratio,
                                     sprintf("heatmap_genes_cluster_%s_list.csv", i)))
  # 基因符号转ENTREZID(org.Mm.eg.db)
  gene.df <- bitr(mygene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ego <- enrichGO(
    gene          = unique(gene.df$ENTREZID),
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",                 # BP/CC/MF
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  pdf(file = file.path(DIR_FIG_td_ratio, sprintf("heatmap_genes_cluster_%s_go.pdf", i)),
      width = 8, height = 8)
  print(dotplot(ego))
  dev.off()
  write.csv(ego@result, file.path(DIR_OBJ_td_ratio,
                                  sprintf("heatmap_genes_cluster_%s_go.csv", i)))
  # 收集top10 BP条目用于条形图
  go_top_list[[i]] <- ego@result |>
    filter(ONTOLOGY == "BP") |>
    arrange(p.adjust) |>
    head(10) |>
    mutate(cluster   = i,
           negLog10p = -log10(p.adjust))
}

# 6.6.2 top10 GO BP条形图(cluster_1=上调, cluster_2=下调)
col_up   <- "#EEC355"   # 金黄
col_down <- "#4A86E8"   # 蓝

plot_go_bar <- function(df, fill_col, panel_title) {
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot(df, aes(x = negLog10p, y = Description)) +
    geom_col(fill = fill_col, width = 0.85) +
    labs(title = panel_title,
         x = expression(-log[10](Adjusted~p-value)),
         y = "GO Term") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = 12) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
          axis.text.y  = element_text(size = 10, colour = "black"),
          axis.text.x  = element_text(size = 10, colour = "black"),
          axis.title.y = element_text(size = 12),
          axis.line.y  = element_blank(),
          axis.ticks.y = element_blank())
}

plot_list <- list()
if (nrow(go_top_list[[2]]) > 0) {
  plot_list[[1]] <- plot_go_bar(go_top_list[[2]], col_down, "Down-regulated Pathways (Cluster 2)")
}
if (nrow(go_top_list[[1]]) > 0) {
  plot_list[[2]] <- plot_go_bar(go_top_list[[1]], col_up,   "Up-regulated Pathways (Cluster 1)")
}

if (length(plot_list) > 0) {
  pdf(file = file.path(DIR_FIG_td_ratio, "heatmap_genes_cluster_goBP_top10_barplot.pdf"),
      width = 18, height = 8)
  grid.arrange(grobs = plot_list, ncol = length(plot_list),
               top = grid::textGrob("GO Biological Process Enrichment",
                                    gp = grid::gpar(fontsize = 16, fontface = "bold")))
  dev.off()
}

# 6.6.3 top100关联基因功能互作网络(FigS4 B, GeneMANIA/Cytoscape)
#     top100关联基因在Cytoscape GeneMANIA app中构建功能互作网络, 并基于
#     GO:0003700(DNA-binding transcription factor activity)及其子条目严格核实TF;
#     复现脚本见 top100_associate_cytoscape/ (classify_tf_go.R、build_go_tf_artifacts.py 等)

# 6.7 按表达模式对top500显著基因聚类(clusterExpressionPatterns)
association_results |>
  filter(pvalue < 0.05) |>
  top_n(500, wt = waldStat) |>
  rownames() -> signifcant_genes

nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce_tradeseq, nPoints = nPointsClus,
                                     genes = signifcant_genes)
saveRDS(clusPat, file = file.path(DIR_OBJ_td_ratio, "traj_clusPat.rds"))
# 与上方保存对象等价的快捷读取(需要时启用):
# clusPat <- readRDS(file.path(DIR_OBJ_td_ratio, "traj_clusPat.rds"))

clusterLabels <- primaryCluster(clusPat$rsec)
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1]  # 去除未聚类基因
top4_cluster <- 1:4
time <- seq(0, 1, length = nPointsClus)

# 各基因聚类的平滑表达模式图(展示top4聚类)
pdf(file.path(DIR_FIG_td_ratio, "smooth_cluster.pdf"), width = 4, height = 4)
for (xx in top4_cluster) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = time,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx), x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(time, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3"), breaks = c("0", "1"))
  print(p)
}
dev.off()

# 6.7 各聚类的KEGG通路富集分析
for (i in top4_cluster) {
  mygene <- colnames(clusPat$rsec)[which(clusterLabels == i)]
  # 利用org.Mm.eg.db整合基因名、ENTREZID、ENSEMBL
  gene.df <- bitr(mygene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"),
                  OrgDb = org.Mm.eg.db)

  ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = "mmu",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                      minGSSize = 10, maxGSSize = 500, use_internal_data = FALSE)
  ekegg@result$Description <- str_split(ekegg@result$Description, pattern = " - Mus",
                                        simplify = TRUE)[, 1]
  pdf(file = file.path(DIR_FIG_td_ratio, sprintf("association_genes_cluster_%s_kegg.pdf", i)),
      width = 6, height = 6)
  print(dotplot(ekegg))
  dev.off()
  write.csv(ekegg@result,
            file.path(DIR_OBJ_td_ratio, sprintf("association_genes_cluster_%s_kegg.csv", i)))
}

# 6.9 轨迹上调基因(trajectory_up signature): cluster 2的基因沿AST表达上升
genes_up <- colnames(clusPat$rsec)[which(clusterLabels == 2)]
cat(genes_up, sep = ",")

# ---- 7. SCENIC 转录因子调控网络分析 ----
# 注: SCENIC流程耗时较长, 内部函数(int/与output/目录)依赖相对路径,
# 需切换工作目录至 output/object/SCENIC_trajectory/ 运行

# 7.1 SCENIC设置初始化
se_obj <- readRDS(file.path(DIR_OBJ, "seurat", "se_trajectory_subset.rds"))
exprMat <- GetAssayData(se_obj, layer = "counts") %>% as.matrix()
cellInfo <- se_obj$seurat_clusters %>% as.data.frame()

org <- "mgi"  # 小鼠基因标注数据库
dbDir <- "~/database/SCENIC_ref"  # SCENIC参考数据库(需自行准备)
myDatasetTitle <- "SCENIC"
data(defaultDbNames)
defaultDbNames$mgi[1] <- "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
defaultDbNames$mgi[2] <- "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
data(list = "motifAnnotations_mgi_v9", package = "RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

scenicOptions <- initializeScenic(org = org,
                                  dbDir = dbDir,
                                  dbs = defaultDbNames[["mgi"]],
                                  datasetTitle = myDatasetTitle,
                                  nCores = 20)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123

setwd(DIR_OBJ_SCENIC)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

# 7.2 SCENIC主流程: 共表达网络 -> GRN构建 -> 细胞评分 -> binarize -> loom导出
source(file.path(DIR_ROOT, "output", "object", "SCENIC_td", "aux_createViewerApp.R"))

scenicOptions@settings$nCores <- 20

# 基因过滤与共表达网络(GENIE3)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered_log <- log2(exprMat_filtered + 1)
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat + 1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, minGenes = 3)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log,
                                        skipHeatmap = TRUE, skipTsne = TRUE)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

# 交互式AUCell阈值确定(shiny app), 选定阈值用于regulon活性二值化
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds", 1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file = getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, skipHeatmaps = TRUE)
tsneAUC(scenicOptions, aucType = "AUC")

# 导出loom文件供后续可视化
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

# 7.3 将Seurat对象完整metadata保存为SCENIC cellInfo
scenicOptions <- readRDS(file.path(DIR_OBJ_SCENIC, "int", "scenicOptions.Rds"))
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 20
scenicOptions@settings$seed <- 123
scenicOptions@inputDatasetInfo$cellInfo <- se_obj@meta.data
cellInfo <- scenicOptions@inputDatasetInfo$cellInfo
saveRDS(cellInfo, file = file.path(DIR_OBJ_SCENIC, "int", "cellInfo.Rds"))

setwd(DIR_ROOT)

# ---- 7.4 SCENIC结果可视化: Regulon活性与基因表达组合图(Fig1) ----

# 从loom读取regulon、AUC与阈值
loomPath <- file.path(DIR_OBJ_SCENIC, "output", "scenic.loom")
loom <- open_loom(loomPath)
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

# 从int目录读取scenicOptions与aucell_regulonAUC(loadInt依赖相对路径)
setwd(DIR_OBJ_SCENIC)
scenicOptions <- readRDS(file.path(DIR_OBJ_SCENIC, "int", "scenicOptions.Rds"))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulons_list <- loadInt(scenicOptions, "regulons")
# SCENIC自适应AUC阈值(每个regulon单独KDE/拐点确定)
aucell_thresholds_raw <- readRDS(file.path(DIR_OBJ_SCENIC, "int", "3.5_AUCellThresholds.Rds"))
setwd(DIR_ROOT)

# 目标转录因子: LUAD端(Nkx2-1/Foxa2)与LUSC端(Sox2/Trp63)关键TF
target_tfs <- c("Nkx2-1", "Foxa2", "Sox2", "Trp63")

# Regulon完整名映射(保留_extended后缀)
regulon_names_full <- onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))
regulon_names_map <- setNames(regulon_names_full[target_tfs], target_tfs)

regulon_gene_counts <- sapply(regulons, length)

# AUC阈值映射: 优先SCENIC/AUCell自适应结果, 次选模糊匹配, 兜底loom向量
build_auc_threshold_map <- function(target_tfs, regulon_names_map, aucell_thr_list,
                                    loom_thr_vec) {
  result <- c()
  for (tf in target_tfs) {
    full_name <- regulon_names_map[tf]
    threshold <- NA_real_
    # 1) 精确匹配int/3.5列表名(与full_name完全一致, 通常为"TF (Ng)")
    if (full_name %in% names(aucell_thr_list)) {
      threshold <- as.numeric(aucell_thr_list[[full_name]]$aucThr$selected)
    }
    # 2) 模糊匹配: 用base TF名搜索(自动区分extended/非extended)
    if (is.na(threshold)) {
      rgx <- paste0("^", gsub("([.-])", "\\\\\\1", sub("_.+$", "", tf)))
      hits <- names(aucell_thr_list)[grepl(rgx, names(aucell_thr_list))]
      if (length(hits) > 0) {
        chosen <- if (grepl("_extended", full_name)) {
          ex <- hits[grepl("_extended", hits)]
          ifelse(length(ex) > 0, ex[1], hits[1])
        } else {
          ne <- hits[!grepl("_extended", hits)]
          ifelse(length(ne) > 0, ne[1], hits[1])
        }
        threshold <- as.numeric(aucell_thr_list[[chosen]]$aucThr$selected)
      }
    }
    # 3) 回退: loom向量(名格式为"TF_(Ng)")
    if (is.na(threshold)) {
      loom_name <- gsub(" ", "_", full_name)
      if (loom_name %in% names(loom_thr_vec)) {
        threshold <- as.numeric(loom_thr_vec[loom_name])
      }
    }
    # 4) 最终兜底: 0.01(避免NA)
    if (is.na(threshold) || !is.finite(threshold)) threshold <- 0.01
    result[full_name] <- threshold
    message(sprintf("  TF=%-10s regulon=%-30s thr=%.5f", tf, full_name, threshold))
  }
  return(result)
}

auc_thresholds <- build_auc_threshold_map(target_tfs, regulon_names_map,
                                          aucell_thresholds_raw, regulonsAucThresholds)

# UMAP坐标与表达矩阵(复用上文exprMat_log)
umap_coord <- se_obj@reductions$umap@cell.embeddings
umap_df <- as.data.frame(umap_coord)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
exprMat_log_df <- as.data.frame(t(exprMat_log))
auc_mat_full <- getAUC(aucell_regulonAUC)

# Row 1: Regulon二值化AUC图(阈值上下分为Active/Inactive)
plot_regulon_binary <- function(tf_name) {
  full_name <- regulon_names_map[tf_name]
  if (is.na(full_name) || is.null(full_name)) {
    candidates <- rownames(aucell_regulonAUC)[grepl(
      paste0("^", gsub("([.-])", "\\\\\\1", tf_name)), rownames(aucell_regulonAUC))]
    if (length(candidates) > 0) full_name <- candidates[1] else full_name <- tf_name
  }
  auc_vals <- auc_mat_full[full_name, ]
  threshold <- as.numeric(auc_thresholds[full_name])
  if (is.na(threshold) || !is.finite(threshold)) threshold <- 0.01

  # 匹配regulon基因数(用于图题)
  n_genes <- NA_integer_
  if (full_name %in% names(regulon_gene_counts)) {
    n_genes <- regulon_gene_counts[full_name]
  }
  base_name <- sub("\\s+\\(\\d+g\\)$", "", full_name)
  if (is.na(n_genes) && base_name %in% names(regulon_gene_counts)) {
    n_genes <- regulon_gene_counts[base_name]
  }
  if (is.na(n_genes)) {
    nn <- grep(paste0("^", gsub("([.-])", "\\\\\\1", base_name)),
               names(regulon_gene_counts), value = TRUE)
    if (length(nn) > 0) n_genes <- regulon_gene_counts[nn[1]]
  }
  if (is.na(n_genes)) n_genes <- "?"

  if (grepl("_extended", full_name)) {
    display_name <- paste0(tf_name, "_extended (", n_genes, "g)")
  } else {
    display_name <- paste0(tf_name, " (", n_genes, "g)")
  }
  # signif(2)两位有效数字, 匹配参考图风格
  subtitle_text <- paste0("Cells with AUC > ",
                          format(signif(threshold, 2), scientific = FALSE, trim = TRUE))

  plot_df <- umap_df
  plot_df$auc <- auc_vals[plot_df$cell]
  plot_df$above_threshold <- plot_df$auc > threshold
  plot_df <- plot_df[!is.na(plot_df$above_threshold), ]

  thr_label <- format(signif(threshold, 2), scientific = FALSE, trim = TRUE)
  plot_df$status <- factor(
    ifelse(plot_df$above_threshold, "Active", "Inactive"),
    levels = c("Inactive", "Active")
  )

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = status, size = status, alpha = status)) +
    scale_color_manual(
      name = "Regulon Activity",
      values = c("Inactive" = "lightgrey", "Active" = "#1a478e"),
      labels = c(
        paste0("AUC \u2264 ", thr_label),
        paste0("AUC > ", thr_label)
      )
    ) +
    scale_size_manual(values = c("Inactive" = 0.5, "Active" = 0.8), guide = "none") +
    scale_alpha_manual(values = c("Inactive" = 0.6, "Active" = 0.8), guide = "none") +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(title = display_name, subtitle = subtitle_text) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )
  return(p)
}

# Row 2: TF基因表达图(log2计数)
plot_gene_expression <- function(tf_name) {
  expr_vals <- exprMat_log_df[[tf_name]]
  if (is.null(expr_vals)) {
    alt_name <- gsub("-", ".", tf_name)
    expr_vals <- exprMat_log_df[[alt_name]]
  }
  if (is.null(expr_vals)) {
    message("  Warning: ", tf_name, " not found in expression matrix, using zeros.")
    expr_vals <- rep(0, nrow(umap_df))
  }

  plot_df <- umap_df
  plot_df$expression <- expr_vals

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_gradientn(colors = col_scale_reds, na.value = "lightgrey") +
    labs(title = paste0(tf_name, " expression"), color = "log2(Expr+1)") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.spacing = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )
  return(p)
}

# 生成8张子图(4 TF x 2行)
p_regulon_list <- lapply(target_tfs, plot_regulon_binary)
p_expr_list <- lapply(target_tfs, plot_gene_expression)

# 左侧行标签(Regulon / Expression)
row_label_regulon <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Regulon",
           angle = 90, size = 6, fontface = "bold") +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(0, 5, 0, 5)
  )

row_label_expr <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Expression",
           angle = 90, size = 6, fontface = "bold") +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(0, 5, 0, 5)
  )

# 栅格化8个主体panel(行标签保持矢量), 按行组合为2x5布局
p_regulon_raster <- lapply(
  p_regulon_list,
  function(p) ggrastr::rasterise(p, dpi = fig_res)
)
p_expr_raster <- lapply(
  p_expr_list,
  function(p) ggrastr::rasterise(p, dpi = fig_res)
)

plot_list <- list(
  row_label_regulon,
  p_regulon_raster[[1]],
  p_regulon_raster[[2]],
  p_regulon_raster[[3]],
  p_regulon_raster[[4]],

  row_label_expr,
  p_expr_raster[[1]],
  p_expr_raster[[2]],
  p_expr_raster[[3]],
  p_expr_raster[[4]]
)

combined_plot <- wrap_plots(
  plot_list,
  nrow = 2,
  ncol = 5,
  byrow = TRUE
) +
  plot_layout(
    widths = c(0.07, 1, 1, 1, 1),
    heights = c(1, 1)
  )

# 输出PDF(栅格panel+矢量文字)与PNG
ggsave(
  filename = file.path(DIR_FIG_SCENIC, "Fig1_Regulon_Expression_2x4.pdf"),
  plot = combined_plot,
  width = 24,
  height = 10,
  units = "in",
  dpi = fig_res,
  device = cairo_pdf
)
ggsave(
  filename = file.path(DIR_FIG_SCENIC, "Fig1_Regulon_Expression_2x4.png"),
  plot = combined_plot,
  width = 24,
  height = 10,
  units = "in",
  dpi = fig_res
)

# ---- 8. signature 评分与 td_group 亚群比较 ----

# 8.1 读取轨迹亚群对象, 按tdTomato_ratio中位(0.5)划分荧光分组
se_obj_traj <- readRDS(file.path(DIR_OBJ, "seurat", "se_trajectory_subset.rds"))
table(se_obj_traj$seurat_clusters)
se_obj_traj$seurat_clusters <- droplevels(se_obj_traj$seurat_clusters)
cols_seurat <- c("0" = "#6EB542", "1" = "pink", "2" = "#e0e896", "3" = "lightblue")
pdf(file = file.path(DIR_FIG_gsva, "umap_illustration_to_show.pdf"), width = 6, height = 5)
DimPlot(se_obj_traj, group.by = "seurat_clusters", cols = cols_seurat,
        reduction = "umap") + theme_void()
dev.off()
cat("Loaded se_trajectory_subset.rds - cells:", ncol(se_obj_traj),
    ", features:", nrow(se_obj_traj), "\n")

se_obj_traj$td_group <- ifelse(se_obj_traj$tdTomato_ratio >= 0.5, "tdTomato", "GFP")
Idents(se_obj_traj) <- "td_group"

# 8.2 tdTomato vs GFP 差异表达基因(Wilcoxon)
td_vs_gfp <- FindMarkers(se_obj_traj, ident.1 = "tdTomato", ident.2 = "GFP")
td_vs_gfp_filter <- td_vs_gfp |>
  filter(p_val_adj < 0.05)
write.csv(td_vs_gfp_filter, file.path(DIR_OBJ_gsva, "deg_td_vs_gfp_filter.csv"))

pdf(file = file.path(DIR_FIG_gsva, "umap_td_gfp.pdf"), width = 6, height = 5)
CellDimPlot(se_obj_traj, group.by = "td_group", palcolor = c("#6EB542", "pink"),
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()

# 8.3 全对象(含免疫细胞)荧光分组UMAP: immune灰色, GFP绿, tdTomato粉
se_all <- readRDS(file.path(DIR_OBJ, "seurat", "se.rds"))
se_all$td_group <- ifelse(se_all$tdTomato_ratio >= 0.5, "tdTomato", "GFP")
se_all$td_group_plot <- ifelse(se_all$first_annotation == "immune", "immune",
                               se_all$td_group)
pdf(file = file.path(DIR_FIG_gsva, "seu_all_umap_td_gfp.pdf"), width = 6, height = 5)
CellDimPlot(se_all, group.by = "td_group_plot",
            palcolor = c("#6EB542", "pink", "lightgrey"),
            reduction = "UMAP", theme_use = "theme_blank")
dev.off()

# 8.4 EMT/上皮 module score(AddModuleScore, ctrl=100)
gene_list <- list(Epithelium = c("Cdh1", "Epcam", "Cldn3", "Cldn4", "Krt8"),
                  EMT = c("Vim", "Fn1", "Cdh2", "Acta2", "S100a4"))
key_pathways <- c("Epithelium", "EMT")
se_obj_traj <- AddModuleScore(
  se_obj_traj,
  features = gene_list,
  name = paste0("Module_", key_pathways, "_"),
  ctrl = 100  # 使用100个控制基因
)
# 清理AddModuleScore自动添加的数字后缀
colnames(se_obj_traj@meta.data) <- sub("^(Module_.+)_\\d+$", "\\1",
                                       colnames(se_obj_traj@meta.data))
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
pdf(file = file.path(DIR_FIG_gsva, "umap_emt_signature_sore.pdf"), width = 12, height = 5)
FeatureDimPlot(se_obj_traj, features = module_scores_name, bg_cutoff = -0.1,
               theme_use = "theme_blank", legend.position = "bottom")
dev.off()

# 8.5 读取GMT基因集文件并转换为小鼠同源基因
gmt_file <- file.path(DIR_DATA, "geneset.hs.gmt")
gmt_data <- readLines(gmt_file)

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

# 人类基因集转换为小鼠同源基因(orthologene/homologene)
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

# 对全部基因集计算module score并保存回Seurat对象
key_pathways <- names(mouse_gene_sets)
cat("Analyzing pathways:", paste(key_pathways, collapse = ", "), "\n")
se_obj_traj <- AddModuleScore(
  se_obj_traj,
  features = mouse_gene_sets,
  name = paste0("Module_", key_pathways, "_"),
  ctrl = 100  # 使用100个控制基因
)
colnames(se_obj_traj@meta.data) <- sub("^(Module_.+)_\\d+$", "\\1",
                                       colnames(se_obj_traj@meta.data))
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
SaveSeuratRds(se_obj_traj, file.path(DIR_OBJ, "seurat", "se_trajectory_subset.rds"))

# 去除module score列名前缀"Module_"
module_scores_name <- grep("Module_", colnames(se_obj_traj@meta.data), value = TRUE)
new_names <- str_remove(module_scores_name, "^Module_")
colnames(se_obj_traj@meta.data)[colnames(se_obj_traj@meta.data) %in% module_scores_name] <- new_names
module_scores_name <- new_names

# 8.6 各通路module得分在td_group间的箱线图比较(Wilcoxon检验)
#     1. 全部通路: 一个PDF分面展示; 2. 选定通路: 多页PDF每通路一页

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
# Wilcoxon检验P值校正方法
p_adjust_method <- "none"
# 全部通路图的每行分面数
all_plot_ncol <- 5
# 输出目录
output_dir <- DIR_FIG_sig

# 构建长格式绘图数据
plot_data <- se_obj_traj@meta.data %>%
  dplyr::select(
    td_group,
    all_of(module_scores_name)
  ) %>%
  dplyr::filter(
    !is.na(td_group)
  )

if (is.factor(plot_data$td_group)) {
  group_levels <- levels(droplevels(plot_data$td_group))
} else {
  group_levels <- unique(as.character(plot_data$td_group))
}
plot_data$td_group <- factor(as.character(plot_data$td_group), levels = group_levels)

if (length(group_levels) < 2) {
  stop("td_group中的有效分组少于两个，无法进行比较。")
}

plot_data_long <- plot_data %>%
  tidyr::pivot_longer(
    cols = -td_group,
    names_to = "Signature",
    values_to = "Expression"
  ) %>%
  dplyr::filter(!is.na(Expression))

# 保持pathway在元数据中的原始顺序
plot_data_long$Signature <- factor(plot_data_long$Signature, levels = module_scores_name)

# 分组颜色
base_colors <- c("#6EB542", "#F4A6B8", "#4E79A7", "#F28E2B", "#59A14F",
                 "#E15759", "#B07AA1", "#76B7B2")
if (length(group_levels) <= length(base_colors)) {
  group_colors <- base_colors[seq_along(group_levels)]
} else {
  group_colors <- grDevices::colorRampPalette(base_colors)(length(group_levels))
}
group_colors <- stats::setNames(group_colors, group_levels)

# 生成所有两两比较组合
comparisons <- combn(group_levels, 2, simplify = FALSE)

# 统计检验函数: 分通路Wilcoxon + 显著性星号 + 横线坐标
calculate_pathway_statistics <- function(data) {
  stat_result <- data %>%
    dplyr::group_by(Signature) %>%
    rstatix::wilcox_test(
      Expression ~ td_group,
      comparisons = comparisons,
      p.adjust.method = p_adjust_method,
      exact = FALSE
    )

  if (p_adjust_method == "none") {
    stat_result <- stat_result %>%
      rstatix::add_significance(p.col = "p")
    stat_result$plot_significance <- stat_result$p.signif
  } else {
    stat_result <- stat_result %>%
      rstatix::add_significance(p.col = "p.adj")
    stat_result$plot_significance <- stat_result$p.adj.signif
  }

  stat_result <- stat_result %>%
    rstatix::add_xy_position(
      x = "td_group",
      fun = "max",
      step.increase = 0.12,
      scales = "free_y"
    )
  stat_result
}

stat_all <- calculate_pathway_statistics(plot_data_long)

# 通用箱线图主题
boxplot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

# 全部通路: 一个PDF分面展示(高度按通路数自适应)
p_all <- ggplot(plot_data_long, aes(x = td_group, y = Expression, fill = td_group)) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.80,
    linewidth = 0.55,
    colour = "black",
    outlier.shape = NA,  # 不显示离群点
    na.rm = TRUE
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  facet_wrap(~Signature, scales = "free_y", ncol = all_plot_ncol) +
  ggpubr::stat_pvalue_manual(
    stat_all,
    label = "plot_significance",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    hide.ns = hide_ns,
    tip.length = 0.01,
    bracket.size = 0.40,
    size = 3.5,
    vjust = 0.2,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian(clip = "off") +
  labs(x = "td_group", y = "Signature score") +
  boxplot_theme

all_plot_nrow <- ceiling(length(module_scores_name) / all_plot_ncol)
all_pdf_height <- max(6, all_plot_nrow * 3)
ggsave(
  filename = file.path(output_dir, "all_canonical_signature_group_boxplot.pdf"),
  plot = p_all,
  width = 14,
  height = all_pdf_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

# 选定通路: 检查通路是否存在于元数据
selected_pathways_exist <- selected_pathways[selected_pathways %in% module_scores_name]
selected_pathways_missing <- setdiff(selected_pathways, module_scores_name)
if (length(selected_pathways_missing) > 0) {
  warning(paste0("以下selected pathways未在元数据中找到，将被忽略：",
                 paste(selected_pathways_missing, collapse = ", ")))
}
if (length(selected_pathways_exist) == 0) {
  stop("selected_pathways中没有任何pathway存在于元数据中。")
}

selected_data_long <- plot_data_long %>%
  dplyr::filter(as.character(Signature) %in% selected_pathways_exist) %>%
  dplyr::mutate(Signature = factor(as.character(Signature), levels = selected_pathways_exist))

# 选定通路横坐标标签: GFP -> ZsGreen+, tdTomato -> tdTomato+("+"上标显示)
selected_x_labels <- function(x) {
  labels <- lapply(x, function(group_name) {
    if (group_name == "GFP") {
      expression(ZsGreen^"+")[[1]]
    } else if (group_name == "tdTomato") {
      expression(tdTomato^"+")[[1]]
    } else {
      group_name
    }
  })
  as.expression(labels)
}

# 为每个选定通路单独生成图
selected_plot_list <- list()
for (pathway in selected_pathways_exist) {
  pathway_data <- selected_data_long %>%
    dplyr::filter(as.character(Signature) == pathway) %>%
    droplevels()

  pathway_stat <- calculate_pathway_statistics(pathway_data)

  p_pathway <- ggplot(pathway_data, aes(x = td_group, y = Expression, fill = td_group)) +
    geom_boxplot(
      width = 0.62,
      alpha = 0.85,
      linewidth = 0.65,
      colour = "black",
      outlier.shape = NA,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(labels = selected_x_labels) +
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    coord_cartesian(clip = "off") +
    labs(title = pathway, x = NULL, y = "Signature score") +
    boxplot_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
    )

  selected_plot_list[[pathway]] <- p_pathway
}

# 选定通路保存为多页PDF(每个pathway单独一页)
selected_pdf_file <- file.path(output_dir, "selected_canonical_signature_group_boxplot.pdf")
pdf(file = selected_pdf_file, width = 4, height = 6, onefile = TRUE, useDingbats = FALSE)
for (pathway in selected_pathways_exist) {
  print(selected_plot_list[[pathway]])
}
dev.off()

# 显著性标记: * P<0.05, ** P<0.01, *** P<0.001, **** P<0.0001;
# hide_ns=TRUE时不显示ns及对应横线
message("All pathways PDF已保存至：",
        file.path(output_dir, "all_canonical_signature_group_boxplot.pdf"))
message("Selected pathways多页PDF已保存至：", selected_pdf_file)
message("Selected PDF共包含 ", length(selected_pathways_exist), " 页，每个pathway一页。")

# ---- 8.7 Nrg1表达与各module得分的相关性分析 ----

# 提取Nrg1表达量与所有module得分
nrg1_expr <- GetAssayData(se_obj_traj, layer = "data")["Nrg1", ]
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
print(head(nrg1_corr_df_overall))

# 按dnb_state分别计算相关性(每个状态至少30个细胞)
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

# 相关性热图(按DNB状态分列)
pdf(file.path(DIR_FIG_sig, "nrg1_module_correlation_heatmap.pdf"), width = 8, height = 10)
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

# 整体相关性柱状图(选定通路)
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
nrg1_corr_df_overall_sorted <- nrg1_corr_df_overall_sorted[
  nrg1_corr_df_overall_sorted$Signature %in% selected_pathways, ]
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

pdf(file.path(DIR_FIG_sig, "nrg1_module_correlation_barplot.pdf"), width = 12, height = 8)
print(p_nrg1_corr_bar)
dev.off()

cat("Nrg1 correlation analysis completed\n")

# ---- 8.8 LUAD/LUSC marker与signature score可视化 ----
#
# 【LUAD/LUSC signature 来源与评分方式】
#
# 来源:
#   下方 luad_signature(10基因)与 scc_signature(11基因, 即LUSC signature)为人工精选的
#   肺腺癌/肺鳞癌谱系标志基因面板, 整理自文献报道的经典marker(原始代码注释"来自项目
#   背景和文献"): LUAD端以 Napsa/Ttf1(即Nkx2-1)/Sftpc/Muc1 等腺型分化基因为主,
#   LUSC端以 Trp63/Krt5/Krt14/Dsg3 等鳞型分化基因为主。
#   面板中多数基因与项目备用代码库 TCGA LUAD vs LUSC 差异表达signature大panel
#   (原 scripts/13.2.7_GSVA_analysis.R 的 tcga_luad_genes/tcga_lusc_genes,
#   ~100基因/谱系, 经 orthogene homologene 人→鼠同源转换)重叠;
#   该大panel的 ssGSEA 通路评分与本节的平均表达评分相互独立。
#
# 评分方式(注意区别于 AddModuleScore/ssGSEA):
#   1) 先过滤掉Seurat对象中不存在的基因(luad_signature/scc_signature过滤行);
#   2) 使用 RNA assay 的 log-normalized 表达矩阵(GetAssayData 默认 data slot);
#   3) 每个细胞的 signature score = 面板内基因表达的平均值(colMeans),
#      即简单平均表达得分, 不做背景基因集校正与中心化;
#   4) td_group(GFP vs tdTomato)组间差异用 Wilcoxon rank-sum 检验(geom_signif)。

# LUAD/LUSC marker基因DotPlot(按td_group分组)
tumor_marker <- list("LUAD" = c("Mlph", "Tmc5", "Sfta3-ps", "Muc5b", "Ttf1", "Napsa",
                                "Ceacam5", "Muc1", "St6galnac1", "Plin2", "Sftpc"),
                     "LUSC" = c("Krt14", "Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5",
                                "Krt6a", "Trim29", "Slc7a5", "Spats2", "Clca2"))
dot_plot <- DotPlot(se_obj_traj, features = tumor_marker, group.by = "td_group") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_gsva, "td_group_marker_gene_dotplot.pdf"), width = 6, height = 4)
print(dot_plot)
dev.off()

# 精选marker基因DotPlot
tumor_marker <- list("LUAD" = c("Sfta3-ps", "Napsa", "Muc1", "Sftpc"),
                     "LUSC" = c("Dsg3", "Trp63", "Krt5", "Krt6a", "Krt14", "Spats2"))
dot_plot <- DotPlot(se_obj_traj, features = tumor_marker, group.by = "td_group") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Dot Plot")
pdf(file = file.path(DIR_FIG_gsva, "selected_td_group_marker_gene_dotplot.pdf"),
    width = 6, height = 4)
print(dot_plot)
dev.off()

# LUAD与SCC signature基因的小提琴图
luad_signature <- c("Mlph", "Tmc5", "Sfta3-ps", "Muc5b", "Ttf1", "Napsa", "Muc1",
                    "St6galnac1", "Plin2", "Sftpc")
scc_signature <- c("Krt14", "Dsg3", "Dsc3", "Calml3", "Trp63", "Krt5", "Krt6a",
                   "Trim29", "Slc7a5", "Spats2", "Clca2")
pdf(file.path(DIR_FIG_gsva, "td_group_AST_related_genes.pdf"), width = 8, height = 8)
VlnPlot(se_obj_traj, features = luad_signature, group.by = "td_group")
VlnPlot(se_obj_traj, features = scc_signature, group.by = "td_group")
dev.off()

# 过滤掉对象中不存在的基因
luad_signature <- luad_signature[luad_signature %in% rownames(se_obj_traj)]
scc_signature <- scc_signature[scc_signature %in% rownames(se_obj_traj)]

# 计算每个细胞的LUAD/SCC signature平均表达得分
DefaultAssay(se_obj_traj) <- "RNA"
se_obj_traj$LUAD_score <- colMeans(GetAssayData(se_obj_traj)[luad_signature, ])
se_obj_traj$SCC_score <- colMeans(GetAssayData(se_obj_traj)[scc_signature, ])
p1 <- VlnPlot(se_obj_traj, features = "LUAD_score", group.by = "td_group",
              cols = c("#6EB542", "pink"), pt.size = 0) +
  NoLegend() + ggtitle("LUAD signature score")
p2 <- VlnPlot(se_obj_traj, features = "SCC_score", group.by = "td_group",
              cols = c("#6EB542", "pink"), pt.size = 0) +
  NoLegend() + ggtitle("SCC signature score")
pdf(file.path(DIR_FIG_gsva, "signature_score_distribution.pdf"), width = 8, height = 10)
print(p1 + p2)
dev.off()

# 带Wilcoxon显著性标记的得分分布图(坐标轴标签用上标+)
se_obj_traj$td_group_plot <- ifelse(se_obj_traj$td_group == "GFP", "ZsGreen+", "tdTomato+")
label_map <- c(
  "GFP" = expression(ZsGreen^"+"),
  "tdTomato" = expression(tdTomato^"+")
)

comp <- list(c("GFP", "tdTomato"))
p1 <- VlnPlot(se_obj_traj, features = "LUAD_score", group.by = "td_group",
              cols = c("#6EB542", "pink"), pt.size = 0) +
  NoLegend() +
  ggtitle("LUAD signature score") +
  geom_signif(
    comparisons = comp,
    test = wilcox.test,       # Wilcoxon rank-sum test
    map_signif_level = TRUE,  # 自动转为 */**/***
    y_position = 2,
    textsize = 4
  ) + ylim(c(min(se_obj_traj$LUAD_score), 2.5)) +
  scale_x_discrete(labels = label_map)
p2 <- VlnPlot(se_obj_traj, features = "SCC_score", group.by = "td_group",
              cols = c("#6EB542", "pink"), pt.size = 0) +
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
pdf(file.path(DIR_FIG_gsva, "signature_score_distribution_with_p.pdf"), width = 4, height = 6)
print(p1)
print(p2)
dev.off()

# LUAD/SCC signature得分在UMAP上的分布
pdf(file = file.path(DIR_FIG_gsva, "umap_signature_sore.pdf"), width = 6, height = 5)
FeatureDimPlot(se_obj_traj, features = "LUAD_score", theme_use = "theme_blank",
               legend.position = "bottom")
FeatureDimPlot(se_obj_traj, features = "SCC_score", theme_use = "theme_blank",
               legend.position = "bottom")
dev.off()
