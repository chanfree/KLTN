# ============================================================================
# 脚本名称: 04_DNB_analysis.R
# 描述: 基于Dynamic Network Biomarker (DNB)理论识别KLTN数据LUAD→LUSC转分化
#       轨迹中的Tipping Point，进行Pre/Post Tipping状态差异表达与差异变异分析、
#       GO/KEGG富集、Trp63核心调控网络及动态TF-靶基因网络构建
# 输入: output/object/seurat/se_trajectory_subset.rds (AST轨迹上皮细胞子集)
#       output/object/seurat/se_harmony_integrated_3samples_epithelium.rds
#         (3样本harmony整合上皮对象, 需已含group/dnb_state/histology/seurat_clusters列)
#       output/object/SCENIC_trajectory/int/2.6_regulons_asGeneSet.Rds (SCENIC regulons)
#       output/object/seurat/merged_obj/merged_obj_deg_tp_vs_luad.csv (可选, 存在则跳过重算)
#       output/object/seurat/merged_obj/merged_obj_deg_tp_vs_scc.csv (可选, 存在则跳过重算)
#       R/custom_functions.R, R/custom_colors.R
# 输出: output/object/td_ratio/DNB_final_df.csv及dnb_deg/dnb_top_var_drivers系列csv
#       output/object/seurat/se_trajectory_subset.rds (追加dnb_state列后回存)
#       output/figures/td_ratio/DNB系列图
#       output/tipping_points/ (GO富集与TF网络结果)
#       output/object/seurat/merged_obj/ (整合对象DEG/marker结果)
#       output/object/TF_regulation/, output/figures/TF_regulation/
#       output/figures/deg_network_core/dynamic_core_TF_network_LUAD_TP_SCC.pdf
# 上游脚本: 03_trajectory_subset_analysis.R (se_trajectory_subset.rds)
# ============================================================================

# ---- 0. 环境与路径配置 ----
set.seed(1)

# 脚本假定从仓库根目录运行
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "object")
fig_res <- 300

# 脚本特定子目录
DIR_FIG_TD <- file.path(DIR_FIG, "td_ratio")            # tdRatio/DNB分析图表
DIR_OBJ_TD <- file.path(DIR_OBJ, "td_ratio")            # tdRatio/DNB分析表格
DIR_FIG_INT <- file.path(DIR_FIG, "integration")        # 整合分析图表
DIR_OBJ_SEU <- file.path(DIR_OBJ, "seurat")             # Seurat对象
DIR_OBJ_MERGED <- file.path(DIR_OBJ_SEU, "merged_obj")  # 整合对象DEG结果
DIR_DEG <- file.path(DIR_RES, "tipping_points")         # tipping point富集/网络结果
DIR_FIG_TF <- file.path(DIR_FIG, "TF_regulation")       # TF调控分析图表
DIR_OBJ_TF <- file.path(DIR_OBJ, "TF_regulation")       # TF调控分析结果
DIR_FIG_NET <- file.path(DIR_FIG, "deg_network_core")   # 核心TF-靶基因网络图

dir.create(path = DIR_FIG_TD, recursive = TRUE)
dir.create(path = DIR_OBJ_TD, recursive = TRUE)
dir.create(path = DIR_FIG_INT, recursive = TRUE)
dir.create(path = DIR_OBJ_SEU, recursive = TRUE)
dir.create(path = DIR_OBJ_MERGED, recursive = TRUE)
dir.create(path = DIR_DEG, recursive = TRUE)
dir.create(path = DIR_FIG_TF, recursive = TRUE)
dir.create(path = DIR_OBJ_TF, recursive = TRUE)
dir.create(path = DIR_FIG_NET, recursive = TRUE)

# 加载R包(原脚本分散/重复加载的包已合并去重; 个人库路径lib.loc已移除)
library(spatstat.utils)  # 原以lib.loc指定个人库加载, 用于解决Seurat依赖版本冲突
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(future)
library(future.apply)
library(zoo)  # 滚动窗口函数
library(tradeSeq)
library(clusterExperiment)
library(rvcheck)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scRNAtoolVis)
library(enrichplot)
library(DOSE)
library(SCP)  # CellDimPlot/plotProportion/SaveSeuratRds
library(harmony)
library(patchwork)
library(ggsci)
library(presto)
library(ggrepel)
library(pheatmap)
library(igraph)
library(purrr)
library(stringr)
library(tibble)
library(ggraph)

# 并行设置(保持原脚本设置)
plan(multisession, workers = 8)
bpparam <- BiocParallel::bpparam()
bpparam$workers <- 10
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)

# 工具函数与配色
source(file.path(DIR_ROOT, "R", "custom_functions.R"))
source(file.path(DIR_ROOT, "R", "custom_colors.R"))

# ---- 1. 数据读取与tdRatio基础可视化 ----
# 读取AST轨迹上皮细胞子集(上游脚本03输出); obj为后续DNB分析的工作副本
seurat_obj <- readRDS(file.path(DIR_OBJ_SEU, "se_trajectory_subset.rds"))
obj <- seurat_obj

# 1.1 tdRatio分布直方图+密度曲线(中位数虚线标注)
meta_dat <- seurat_obj@meta.data
data_median <- median(meta_dat$tdTomato_ratio)
pdf(file.path(DIR_FIG_TD, "td_ratio_histo_distribution.pdf"), width = 8, height = 4)
ggplot(meta_dat, aes(x = tdTomato_ratio)) +
  geom_histogram(
    aes(y = after_stat(density), fill = after_stat(x)),  # 按tdRatio本身的值填色
    bins = 30,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#6EB542",   # ZsGreen
    high = "#f6c3cb",  # tdTomato
    name = "tdRatio"
  ) +
  geom_density(
    color = "red",
    linewidth = 1
  ) +
  geom_vline(
    aes(xintercept = data_median, linetype = "Median"),
    color = "black",
    linewidth = 1.5
  ) +
  scale_linetype_manual(name = "", values = c("Median" = "dashed")) +
  labs(
    title = "Histogram with Density Curve",
    x = "tdTomato_ratio",
    y = "Density"
  ) +
  theme_classic()
dev.off()

# 1.2 tdRatio在UMAP上的空间分布(三种配色版本)
x_limits <- range(seurat_obj@reductions$umap@cell.embeddings[, 1], na.rm = TRUE)
y_limits <- range(seurat_obj@reductions$umap@cell.embeddings[, 2], na.rm = TRUE)

rel_order_umap <- FeaturePlot(seurat_obj, "tdTomato_ratio") +
  scale_colour_gradientn(
    colours = c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF"),
    na.value = "transparent",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = c("0.0 (LUAD)", "0.2", "0.4", "0.6", "0.8", "1.0 (LUSC)"),
    name = "tdTomato UMI Ratio",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("AST level") +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
pdf(file.path(DIR_FIG_TD, "td_ratio_umap.pdf"), width = 6, height = 5)
print(rel_order_umap)
dev.off()

rel_order_umap_2 <- FeaturePlot(seurat_obj, "tdTomato_ratio", order = TRUE) +
  scale_colour_gradientn(
    colours = c("#6EB542", "#DDECCB", "#F8F5F2", "#F6D8DE", "#F6C3CB"),
    values = scales::rescale(c(0, 0.30, 0.50, 0.75, 1)),
    limits = c(0, 1),
    na.value = "grey90",
    breaks = seq(0, 1, 0.2),
    labels = c("0.0 (ZsGreen)", "0.2", "0.4", "0.6", "0.8", "1.0 (tdTomato)"),
    name = "tdRatio",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )
pdf(file.path(DIR_FIG_TD, "td_ratio_umap_ver_2.pdf"), width = 6, height = 5)
print(rel_order_umap_2)
dev.off()

rel_order_umap_3 <- FeaturePlot(seurat_obj, "tdTomato_ratio", order = TRUE) +
  scale_colour_gradientn(colours = c("#3B6BA5", "#F6F6F6", "#E6B4B4", "#B64A5A"))
pdf(file.path(DIR_FIG_TD, "td_ratio_umap_ver_3.pdf"), width = 6, height = 5)
print(rel_order_umap_3)
dev.off()

# ---- 2. DNB滑动窗口分析识别Tipping Point ----
# DNB核心公式: CI = (PCC_in / PCC_out) * SD_in
# 沿tdRatio排序的细胞上滑动窗口计算CI, CI峰值窗口即Tipping Point
plan(multisession, workers = 20)

# 计算单个窗口(表达矩阵)的DNB得分:
# 取top高变异基因->计算绝对Pearson相关矩阵->层次聚类划分模块->返回CI最大模块
calculate_dnb_window <- function(expr_mat, n_top_genes = 1000, n_clusters = 10) {
  # 1. 预过滤: 使用top变异基因以降低噪声与计算量(DNB理论关注"波动"基因)
  gene_vars <- apply(expr_mat, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_top_genes, length(gene_vars))]
  sub_mat <- expr_mat[top_genes, ]

  # 2. 计算每个基因的标准差(SD)
  sd_vector <- apply(sub_mat, 1, sd)

  # 3. 计算相关矩阵(绝对Pearson相关, DNB方法要求)
  cor_mat <- abs(cor(t(sub_mat), method = "pearson"))

  # 4. 层次聚类划分模块(距离 = 1 - correlation, 确定性聚类)
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "average")
  modules <- cutree(hc, k = n_clusters)

  # 5. 遍历每个模块, 计算CI并返回最大者
  results <- lapply(1:n_clusters, function(cluster_id) {
    in_genes <- names(modules)[modules == cluster_id]
    out_genes <- names(modules)[modules != cluster_id]

    # 模块基因数<3则跳过
    if (length(in_genes) < 3) return(NULL)

    # A. PCC_in: 组内平均绝对相关(取上三角以排除自身相关1.0)
    pcc_in_vals <- cor_mat[in_genes, in_genes]
    pcc_in <- mean(pcc_in_vals[upper.tri(pcc_in_vals)])

    # B. PCC_out: 组内外平均绝对相关
    pcc_out <- mean(cor_mat[in_genes, out_genes])

    # C. SD_in: 组内平均标准差
    sd_in <- mean(sd_vector[in_genes])

    # D. CI分数(分母加小epsilon避免除零)
    ci_score <- (pcc_in / (pcc_out + 1e-6)) * sd_in

    return(list(CI = ci_score, Module_Genes = in_genes))
  })

  # 过滤空模块并取CI最大者
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) return(list(Max_CI = 0, Top_Genes = NULL))

  ci_values <- sapply(results, function(x) x$CI)
  best_idx <- which.max(ci_values)

  return(list(Max_CI = ci_values[best_idx], Top_Genes = results[[best_idx]]$Module_Genes))
}

# 使用标准化后的data slot(非counts), 按tdRatio升序排列细胞
ordered_cells <- names(sort(obj$tdTomato_ratio))
expr_data <- GetAssayData(obj, slot = "data")

# 定义滑动窗口参数
window_size <- 300  # 窗口大小(约占总细胞数5-10%)
step_size <- 60     # 滑动步长
n_cells <- length(ordered_cells)
starts <- seq(1, n_cells - window_size + 1, by = step_size)

# DNB主循环(耗时步骤: 每个窗口做1000基因的层次聚类)
dnb_results <- list()
print("Starting DNB Sliding Window Analysis...")
pb <- txtProgressBar(min = 0, max = length(starts), style = 3)
for (i in seq_along(starts)) {
  start_idx <- starts[i]
  end_idx <- start_idx + window_size - 1

  # 当前窗口细胞及平均tdRatio(作图x轴)
  window_cells <- ordered_cells[start_idx:end_idx]
  avg_td <- mean(obj$tdTomato_ratio[window_cells])

  # 提取窗口表达矩阵, 剔除窗口内零方差基因避免报错
  win_mat <- as.matrix(expr_data[, window_cells])
  win_mat <- win_mat[apply(win_mat, 1, var) > 0, ]

  # 对窗口top 1000高变异基因、8个模块计算DNB
  res <- calculate_dnb_window(win_mat, n_top_genes = 1000, n_clusters = 8)

  dnb_results[[i]] <- data.frame(
    Window_ID = i,
    Avg_tdRatio = avg_td,
    CI_Score = res$Max_CI,
    Top_Genes = paste(res$Top_Genes, collapse = ",")
  )

  setTxtProgressBar(pb, i)
}
close(pb)

# 汇总结果并保存
final_df <- do.call(rbind, dnb_results)
write.csv(final_df, file = file.path(DIR_OBJ_TD, "DNB_final_df.csv"))
# 快捷方式: 已运行过DNB滑动窗口分析时, 可注释上方循环并直接读取缓存
# final_df <- read.csv(file.path(DIR_OBJ_TD, "DNB_final_df.csv"), row.names = 1)

# ---- 3. Tipping Point识别与CI轨迹可视化 ----
# CI最大值窗口即为Tipping Point
peak_row <- final_df[which.max(final_df$CI_Score), ]
print(paste("Tipping Point found at tdRatio:", round(peak_row$Avg_tdRatio, 3)))

# 3.1 CI随tdRatio变化曲线(标注Tipping Point)
pdf(file.path(DIR_FIG_TD, "DNB_tipping_point_window_plot.pdf"), width = 6, height = 4)
ggplot(final_df, aes(x = Avg_tdRatio, y = CI_Score)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(color = "black", size = 1.5) +
  geom_vline(
    xintercept = peak_row$Avg_tdRatio,
    color = "#E64B35",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = peak_row$Avg_tdRatio,
    y = 5.8,
    label = paste0("Tipping point\n(CI = ", round(peak_row$CI_Score, 2), ")"),
    color = "#E64B35",
    size = 4.5,
    hjust = -0.15
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = sprintf("%.2f", seq(0, 1, 0.25))
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, 2)
  ) +
  labs(
    x = "tdRatio",
    y = "Composite Index (CI)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )
dev.off()

# 3.2 出版风格CI轨迹图(高亮Tipping Point并标注公式)
pdf(file.path(DIR_FIG_TD, "DNB_tipping_point_window_plot_2.pdf"), width = 18, height = 6)
ggplot(final_df, aes(x = Avg_tdRatio, y = CI_Score)) +
  geom_line(color = "#2C5F8A", linewidth = 1) +
  geom_point(
    color = "#1A3A66",
    fill = "white",
    shape = 21,
    size = 2.2,
    stroke = 1,
    alpha = 0.75
  ) +
  geom_point(
    data = peak_row,
    aes(x = Avg_tdRatio, y = CI_Score),
    color = "#D62728",
    fill = "#D62728",
    shape = 21,
    size = 4
  ) +
  geom_vline(
    xintercept = peak_row$Avg_tdRatio,
    color = "#D62728",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  annotate(
    "label",
    x = peak_row$Avg_tdRatio,
    y = peak_row$CI_Score,
    label = paste0(
      "DNB tipping point\n",
      "tdRatio = ", round(peak_row$Avg_tdRatio, 3),
      "\nCI = ", round(peak_row$CI_Score, 2)
    ),
    color = "#B22222",
    fill = "white",
    size = 4,
    fontface = "bold",
    label.size = 0.3,
    hjust = -0.05,
    vjust = -0.3
  ) +
  labs(
    title = "Identification of DNB Tipping Point",
    subtitle = "Composite Index (CI) = (PCC_in / PCC_out) × SD_in",
    x = "tdRatio",
    y = "Composite Index (CI)"
  ) +
  scale_y_continuous(
    limits = c(min(final_df$CI_Score) * 0.85, max(final_df$CI_Score) * 1.15)
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.margin = margin(15, 20, 15, 15)
  )
dev.off()

# 3.3 查看驱动Tipping Point的DNB模块基因
print("Genes driving the Tipping Point module:")
print(strsplit(peak_row$Top_Genes, ",")[[1]])

# ---- 4. 细胞状态标注: Pre_Tipping / Tipping_Point / Post_Tipping ----
# 依据CI峰值窗口在排序细胞上的位置, 定义临界点前后各2个窗口的细胞状态
peak_window_row <- final_df[which.max(final_df$CI_Score), ]
peak_window_id <- peak_window_row$Window_ID

tipping_start_idx <- starts[peak_window_id]
tipping_end_idx <- tipping_start_idx + window_size - 1

# Pre-Tipping: 峰值窗口前2个窗口(不越界至索引1以下)
pre_start_idx <- max(1, tipping_start_idx - (window_size * 2))
pre_end_idx <- max(1, tipping_start_idx - 1)

# Post-Tipping: 峰值窗口后2个窗口(不超出细胞总数)
post_start_idx <- tipping_end_idx + 1
post_end_idx <- min(length(ordered_cells), post_start_idx + (window_size * 2))

# 映射到细胞名并写入metadata
cells_tip <- ordered_cells[tipping_start_idx:tipping_end_idx]
cells_pre <- ordered_cells[pre_start_idx:pre_end_idx]
cells_post <- ordered_cells[post_start_idx:post_end_idx]

obj$dnb_state <- "Other"
obj$dnb_state[cells_pre] <- "Pre_Tipping"
obj$dnb_state[cells_tip] <- "Tipping_Point"
obj$dnb_state[cells_post] <- "Post_Tipping"

# 回存至上游对象路径(覆盖保存, 新增dnb_state列; 注意此时dnb_state为字符型,
# 因子转换在第5章进行, 第11章依赖此处保存的字符型列)
SaveSeuratRds(obj, file = file.path(DIR_OBJ_SEU, "se_trajectory_subset.rds"))
Idents(obj) <- "dnb_state"

# ---- 5. dnb_state在UMAP上的可视化 ----
# 5.1 SCP CellDimPlot版本
# 因子水平顺序决定图层渲染顺序, Tipping_Point置于最后确保高光细胞绘制在顶层
obj$dnb_state <- factor(obj$dnb_state,
                        levels = c("Other", "Pre_Tipping", "Post_Tipping", "Tipping_Point"))

# 统一的Progression States配色(色盲友好且具过渡感)
state_cols <- c(
  "Other" = "#FEF1F1",
  "Pre_Tipping" = "#7570B3",
  "Post_Tipping" = "#1F78B4",
  "Tipping_Point" = "#E6AB02"
)

p_dnb <- CellDimPlot(
  obj,
  group.by = "dnb_state",
  pt.size = 0.5,
  palcolor = state_cols,
  theme_use = "theme_blank"
) +
  guides(color = guide_legend(
    override.aes = list(size = 4),
    title = "DNB State"
  )) +
  labs(title = "DNB Tipping Point States on UMAP") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right"
  )
pdf(file.path(DIR_FIG_TD, "DNB_state_umap.pdf"), width = 6, height = 5)
print(p_dnb)
dev.off()

# 5.2 分层绘图版本(手动控制绘制顺序: Other底层 -> Tipping_Point顶层)
umap_df <- cbind(
  as.data.frame(obj@reductions$umap@cell.embeddings),
  obj@meta.data[, "dnb_state", drop = F]
)

df_other <- umap_df %>% filter(dnb_state == "Other")
df_pre <- umap_df %>% filter(dnb_state == "Pre_Tipping")
df_post <- umap_df %>% filter(dnb_state == "Post_Tipping")
df_tip <- umap_df %>% filter(dnb_state == "Tipping_Point")

state_cols <- c(
  "Other" = "grey90",
  "Pre_Tipping" = "#7570B3",
  "Post_Tipping" = "#1F78B4",
  "Tipping_Point" = "#E6AB02"
)

p <- ggplot() +
  geom_point(data = df_other, aes(x = umap_1, y = umap_2, color = dnb_state), size = 0.4, alpha = 0.7) +
  geom_point(data = df_pre, aes(x = umap_1, y = umap_2, color = dnb_state), size = 0.4, alpha = 0.8) +
  geom_point(data = df_post, aes(x = umap_1, y = umap_2, color = dnb_state), size = 0.4, alpha = 0.8) +
  geom_point(data = df_tip, aes(x = umap_1, y = umap_2, color = dnb_state), size = 0.6, alpha = 1) +
  scale_color_manual(values = state_cols) +
  labs(title = "DNB Tipping Point States on umap", x = "umap 1", y = "umap 2", color = "DNB State") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "black"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )
pdf(file.path(DIR_FIG_TD, "DNB_state_umap_publish.pdf"), width = 6, height = 5)
print(p)
dev.off()

# ---- 6. DNB状态间差异表达与火山图 ----
Idents(obj) <- "dnb_state"

# 6.1 Tipping_Point vs Pre_Tipping差异基因火山图
de_genes <- FindMarkers(
  obj,
  group.by = "dnb_state",
  ident.1 = "Tipping_Point",
  ident.2 = "Pre_Tipping",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

de_genes$significant_0.05 <- de_genes$p_val_adj < 0.05
de_genes$sign_up_down <- ifelse(
  test = (de_genes$significant_0.05 & de_genes$avg_log2FC > 0),
  yes = "sign_up",
  no = ifelse((de_genes$significant_0.05 & de_genes$avg_log2FC < 0),
              yes = "sign_down", no = NA)
)
cols_sign_up_down <- setNames(
  c(cols_fluorescence[["tdTomato"]], cols_fluorescence[["GFP"]], "grey40"),
  nm = c("sign_up", "sign_down", "not_sign")
)

# LUAD/LUSC标志基因单独标注
luad_features_filtered <- c("Napsa", "Sftpa1", "Sftpc")
lusc_features_filtered <- c("Sox2", "Krt6a", "Trp63")
de_genes$luad_filter <- ifelse(rownames(de_genes) %in% luad_features_filtered, rownames(de_genes), "")
de_genes$lusc_filter <- ifelse(rownames(de_genes) %in% lusc_features_filtered, rownames(de_genes), "")

# 显著基因标签(p_adj<1e-4且|log2FC|>1.5)
de_genes$gene_labels <- ifelse((de_genes$p_val_adj < 0.0001 & abs(de_genes$avg_log2FC) > 1.5),
                               rownames(de_genes), "")

de_genes$negLog10_p_val_adj_capped <- ifelse(-log10(de_genes$p_val_adj) > 10, 10, -log10(de_genes$p_val_adj))
minmax_xilm <- round(max(abs(de_genes$avg_log2FC)) + 1)

p_volcano_all <- ggplot(de_genes, aes(x = avg_log2FC, y = negLog10_p_val_adj_capped,
                                      fill = sign_up_down, color = sign_up_down)) +
  geom_hline(yintercept = c(0, 10), linetype = "dashed", color = "grey") +
  geom_point(size = .5, alpha = 0.8, shape = 20) +
  scale_color_manual(values = cols_sign_up_down) +
  xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "DEGs Tipping vs Pre_tipping",
       x = "log2(fold-change)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8)
  )

p_volcano_all_txt <- p_volcano_all +
  ggrepel::geom_text_repel(
    data = subset(de_genes, gene_labels != "" & !is.na(gene_labels)),
    mapping = aes(x = avg_log2FC, y = negLog10_p_val_adj_capped, label = gene_labels),
    color = "black",
    size = 2,
    box.padding = 0.8,
    max.overlaps = Inf
  )

p_volcano_all <- ggrastr::rasterize(p_volcano_all, layers = "Point", dpi = 600)
p_volcano_all_txt <- ggrastr::rasterize(p_volcano_all_txt, layers = "Point", dpi = 600)

# LUAD(绿)/LUSC(红)标志基因标注版本
p_volcano_all_txt_luad_lusc <- p_volcano_all +
  ggrepel::geom_text_repel(
    data = subset(de_genes, luad_filter != "" & !is.na(luad_filter)),
    mapping = aes(x = avg_log2FC, y = negLog10_p_val_adj_capped, label = luad_filter),
    color = "darkgreen",
    size = 3,
    box.padding = 0.8,
    max.overlaps = Inf
  ) +
  ggrepel::geom_text_repel(
    data = subset(de_genes, lusc_filter != "" & !is.na(lusc_filter)),
    mapping = aes(x = avg_log2FC, y = negLog10_p_val_adj_capped, label = lusc_filter),
    color = "darkred",
    size = 3,
    box.padding = 0.8,
    max.overlaps = Inf
  )
pdf(file.path(DIR_FIG_TD, "DNB_tipping_vs_pre_volcano.pdf"), width = 6, height = 6)
p_volcano_all_txt_luad_lusc
p_volcano_all
p_volcano_all_txt
dev.off()

# 6.2 Tipping_Point vs 其余所有细胞差异基因火山图
de_genes <- FindMarkers(
  obj,
  group.by = "dnb_state",
  ident.1 = "Tipping_Point",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

de_genes$significant_0.05 <- de_genes$p_val_adj < 0.05
de_genes$sign_up_down <- ifelse(
  test = (de_genes$significant_0.05 & de_genes$avg_log2FC > 0),
  yes = "sign_up",
  no = ifelse((de_genes$significant_0.05 & de_genes$avg_log2FC < 0),
              yes = "sign_down", no = NA)
)
cols_sign_up_down <- setNames(
  c(cols_fluorescence[["tdTomato"]], cols_fluorescence[["GFP"]], "grey40"),
  nm = c("sign_up", "sign_down", "not_sign")
)

# 显著基因标签(p_adj<1e-4且|log2FC|>1)
de_genes$gene_labels <- ifelse((de_genes$p_val_adj < 0.0001 & abs(de_genes$avg_log2FC) > 1),
                               rownames(de_genes), "")

de_genes$negLog10_p_val_adj_capped <- ifelse(-log10(de_genes$p_val_adj) > 10, 10, -log10(de_genes$p_val_adj))
minmax_xilm <- round(max(abs(de_genes$avg_log2FC)) + 1)

p_volcano_all <- ggplot(de_genes, aes(x = avg_log2FC, y = negLog10_p_val_adj_capped,
                                      fill = sign_up_down, color = sign_up_down)) +
  geom_hline(yintercept = c(0, 10), linetype = "dashed", color = "grey") +
  geom_point(size = .5, alpha = 0.8, shape = 20) +
  scale_color_manual(values = cols_sign_up_down) +
  xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "DEGs Tipping vs all others",
       x = "log2(fold-change)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8)
  )

p_volcano_all_txt <- p_volcano_all +
  ggrepel::geom_text_repel(
    data = subset(de_genes, gene_labels != "" & !is.na(gene_labels)),
    mapping = aes(x = avg_log2FC, y = negLog10_p_val_adj_capped, label = gene_labels),
    color = "black",
    size = 2,
    box.padding = 0.8,
    max.overlaps = Inf
  )

p_volcano_all <- ggrastr::rasterize(p_volcano_all, layers = "Point", dpi = 600)
p_volcano_all_txt <- ggrastr::rasterize(p_volcano_all_txt, layers = "Point", dpi = 600)

pdf(file.path(DIR_FIG_TD, "DNB_tipping_vs_all_other_volcano.pdf"), width = 6, height = 6)
p_volcano_all_txt
p_volcano_all
dev.off()

# ---- 7. Tipping vs Pre: 差异表达、差异变异与富集分析 ----
# 分析A(标准DE)找"效应基因"; 分析B(差异变异)找"调控基因"(临界状态下闪烁的TF等)
table(obj$dnb_state)
Idents(obj) <- "dnb_state"

# 7.1 分析A: Tipping_Point vs Pre_Tipping标准差异表达
de_tip_vs_pre <- FindMarkers(
  obj,
  ident.1 = "Tipping_Point",
  ident.2 = "Pre_Tipping",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
head(de_tip_vs_pre)
write.csv(de_tip_vs_pre, file = file.path(DIR_OBJ_TD, "dnb_deg_tip_vs_pre.csv"))

# 7.2 分析B: 差异变异分析(Differential Variability)
# 计算两状态各基因方差, 筛选Tipping期变异倍数>2且有一定表达的"闪烁"基因
mat_tip <- GetAssayData(obj, slot = "data")[, cells_tip]
mat_pre <- GetAssayData(obj, slot = "data")[, cells_pre]

var_tip <- apply(mat_tip, 1, var)
var_pre <- apply(mat_pre, 1, var)

# 方差比值(Tipping/Pre), 加小epsilon避免除零
var_ratio <- (var_tip + 1e-6) / (var_pre + 1e-6)

# 过滤高变异驱动基因: 变异倍数>2 且 平均表达>0.1(避免dropout噪声)
mean_expr <- apply(mat_tip, 1, mean)
variable_drivers <- names(which(var_ratio > 2 & mean_expr > 0.1))

top_var_drivers <- sort(var_ratio[variable_drivers], decreasing = TRUE)
write.csv(top_var_drivers, file = file.path(DIR_OBJ_TD, "dnb_top_var_drivers_tip_vs_pre.csv"))
# 快捷方式: 已计算过时可直接读取缓存(注意缓存为data.frame, 需用rownames取基因名)
# top_var_drivers <- read.csv(file.path(DIR_OBJ_TD, "dnb_top_var_drivers_tip_vs_pre.csv"), row.names = 1)
print("Top Differential Variability Genes:")
print(head(top_var_drivers, 10))

# 7.3 top差异变异基因的KEGG/GO富集
# (原脚本从CSV读回data.frame用rownames取基因名; 此处直接使用计算得到的命名向量, 改用names)
mygene <- names(head(top_var_drivers, 100))
gene.df <- bitr(mygene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)

ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = "mmu",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, use_internal_data = F)
i <- "Top_Differential_Variability_Genes"
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_kegg.pdf", i)), width = 8, height = 8)
print(dotplot(ekegg))
dev.off()
write.csv(ekegg@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_kegg.csv", i)))

ego <- enrichGO(
  gene = unique(gene.df$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_go.pdf", i)), width = 8, height = 8)
print(dotplot(ego))
dev.off()
write.csv(ego@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_go.csv", i)))

# 7.4 Tipping vs Pre上调DEG的KEGG/GO富集
mygene <- de_tip_vs_pre |>
  filter(p_val_adj < 0.05, avg_log2FC > 0) |>
  row.names()
gene.df <- bitr(mygene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)

ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = "mmu",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, use_internal_data = F)
i <- "positive_DEG"
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_kegg.pdf", i)), width = 8, height = 8)
print(dotplot(ekegg))
dev.off()
write.csv(ekegg@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_kegg.csv", i)))

ego <- enrichGO(
  gene = unique(gene.df$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_go_MF.pdf", i)), width = 8, height = 8)
print(dotplot(ego))
dev.off()
write.csv(ego@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_go_MF.csv", i)))

# ---- 8. Tipping vs Post: 差异表达、差异变异与富集分析 ----
table(obj$dnb_state)
Idents(obj) <- "dnb_state"

# 8.1 分析A: Tipping_Point vs Post_Tipping标准差异表达
de_tip_vs_post <- FindMarkers(
  obj,
  ident.1 = "Tipping_Point",
  ident.2 = "Post_Tipping",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
head(de_tip_vs_post)
write.csv(de_tip_vs_post, file = file.path(DIR_OBJ_TD, "dnb_deg_tip_vs_post.csv"))

# 8.2 分析B: 差异变异分析(Tipping vs Post)
mat_tip <- GetAssayData(obj, slot = "data")[, cells_tip]
mat_post <- GetAssayData(obj, slot = "data")[, cells_post]

var_tip <- apply(mat_tip, 1, var)
var_post <- apply(mat_post, 1, var)

var_ratio <- (var_tip + 1e-6) / (var_post + 1e-6)

mean_expr <- apply(mat_tip, 1, mean)
variable_drivers <- names(which(var_ratio > 2 & mean_expr > 0.1))

top_var_drivers <- sort(var_ratio[variable_drivers], decreasing = TRUE)
write.csv(top_var_drivers, file = file.path(DIR_OBJ_TD, "dnb_top_var_drivers_tip_vs_post.csv"))

print("Top Differential Variability Genes:")
print(head(top_var_drivers, 10))

# 8.3 top差异变异基因的KEGG/GO富集
mygene <- names(head(top_var_drivers, 100))
gene.df <- bitr(mygene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)

ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = "mmu",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, use_internal_data = F)
i <- "Top_Differential_Variability_Genes_tipping_vs_post"
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_kegg.pdf", i)), width = 8, height = 8)
print(dotplot(ekegg))
dev.off()
write.csv(ekegg@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_kegg.csv", i)))

ego <- enrichGO(
  gene = unique(gene.df$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_go.pdf", i)), width = 8, height = 8)
print(dotplot(ego))
dev.off()
write.csv(ego@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_go.csv", i)))

# 8.4 Tipping vs Post上调DEG的KEGG/GO富集
mygene <- de_tip_vs_post |>
  filter(p_val_adj < 0.05, avg_log2FC > 0) |>
  row.names()
gene.df <- bitr(mygene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)

ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = "mmu",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, use_internal_data = F)
i <- "positive_DEG_tipping_vs_post"
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_kegg.pdf", i)), width = 8, height = 8)
print(dotplot(ekegg))
dev.off()
write.csv(ekegg@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_kegg.csv", i)))

ego <- enrichGO(
  gene = unique(gene.df$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
pdf(file.path(DIR_FIG_TD, sprintf("tipping_points_DNB_%s_go.pdf", i)), width = 8, height = 8)
print(dotplot(ego))
dev.off()
write.csv(ego@result, file.path(DIR_OBJ_TD, sprintf("tipping_points_DNB_%s_go.csv", i)))

# ---- 9. 整合队列: Tipping_Point vs LUAD/SCC差异基因 ----
set.seed(123)
plan(sequential)  # 差异分析阶段改为串行
options(future.globals.maxSize = 120000 * 1024^2, future.seed = TRUE)

deg_tp_vs_LUAD_file <- file.path(DIR_OBJ_MERGED, "merged_obj_deg_tp_vs_luad.csv")
deg_tp_vs_SCC_file <- file.path(DIR_OBJ_MERGED, "merged_obj_deg_tp_vs_scc.csv")

# 始终加载harmony整合对象(第11章仍需使用; 对象需已含dnb_state列)
se_merged <- readRDS(file.path(DIR_OBJ_SEU, "se_harmony_integrated_3samples_epithelium.rds"))

# DEG文件已存在则跳过重算(缓存机制)
if (file.exists(deg_tp_vs_LUAD_file) && file.exists(deg_tp_vs_SCC_file)) {
  cat("DEG文件已存在，跳过重新计算...\n")
} else {
  # AST样本细胞按dnb_state细分, 其余样本按原组别
  se_merged$tipping_group <- ifelse(se_merged$group == "AST", se_merged$dnb_state, se_merged$group)
  table(se_merged$tipping_group)
  se_merged <- subset(se_merged, !is.na(group))
  cat("预处理SCT数据用于差异分析...\n")
  se_merged <- PrepSCTFindMarkers(se_merged)

  cat("\n[1/3] 寻找Cluster Tipping_Point vs LUAD的差异基因...\n")
  Idents(se_merged) <- "tipping_group"
  deg_tp_vs_LUAD <- FindMarkers(se_merged, ident.1 = "Tipping_Point", ident.2 = "LUAD", assay = "SCT",
                                min.pct = 0.1, logfc.threshold = 0.25)
  deg_tp_vs_LUAD$gene <- rownames(deg_tp_vs_LUAD)
  cat("Cluster Tipping_Point vs LUAD DEG数量:", nrow(deg_tp_vs_LUAD), "\n")
  cat("Cluster Tipping_Point vs LUAD 上调基因:", sum(deg_tp_vs_LUAD$avg_log2FC > 0), "\n")
  write.csv(deg_tp_vs_LUAD, file = deg_tp_vs_LUAD_file, row.names = FALSE)

  cat("\n[1/3] 寻找Cluster Tipping_Point vs SCC的差异基因...\n")
  Idents(se_merged) <- "tipping_group"
  deg_tp_vs_SCC <- FindMarkers(se_merged, ident.1 = "Tipping_Point", ident.2 = "SCC", assay = "SCT",
                               min.pct = 0.1, logfc.threshold = 0.25)
  deg_tp_vs_SCC$gene <- rownames(deg_tp_vs_SCC)
  cat("Cluster Tipping_Point vs SCC DEG数量:", nrow(deg_tp_vs_SCC), "\n")
  cat("Cluster Tipping_Point vs SCC 上调基因:", sum(deg_tp_vs_SCC$avg_log2FC > 0), "\n")
  write.csv(deg_tp_vs_SCC, file = deg_tp_vs_SCC_file, row.names = FALSE)
}

dir.create(file.path(DIR_DEG, "GO_enrichment"), recursive = TRUE)
dir.create(file.path(DIR_DEG, "network"), recursive = TRUE)

# ---- 10. GO富集、鳞状通路分析与Trp63调控网络 ----
deg_datasets <- list(
  "tp_vs_LUAD" = deg_tp_vs_LUAD_file,
  "tp_vs_SCC" = deg_tp_vs_SCC_file
)

# 鳞状分化相关通路关键词(用于从GO结果中抽取鳞状相关条目)
squamous_keywords <- c("keratin", "epidermal", "desmosome", "epithelial", "squamous",
                       "cell junction", "tight junction", "adherens junction",
                       "cadherin", "cell fate", "differentiation", "morphogenesis",
                       "epithelium development", "skin development", "cornification")

# 10.1 GO富集函数: 对基因列表做BP/CC/MF富集并输出dotplot/barplot/网络图
run_go_enrichment_diverse <- function(genes, dataset_name, direction, ont_list = c("BP", "CC", "MF")) {
  cat(sprintf("\n[%s] 进行%s - %s GO富集分析...\n", dataset_name, direction, paste(ont_list, collapse = ",")))

  if (length(genes) < 3) {
    cat("基因数不足，跳过\n")
    return(NULL)
  }

  # SYMBOL转ENTREZID
  gene_ids <- tryCatch({
    suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  }, error = function(e) {
    cat("ID转换失败:", e$message, "\n")
    return(NULL)
  })

  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    cat("警告: 没有成功映射的基因ID\n")
    return(NULL)
  }
  cat(sprintf("  成功映射 %d/%d 个基因\n", nrow(gene_ids), length(genes)))

  results <- list()
  for (ont in ont_list) {
    cat(sprintf("  - %s富集分析 (%s)...\n", ont, direction))

    ego <- tryCatch({
      enrichGO(
        gene = unique(gene_ids$ENTREZID),
        OrgDb = org.Mm.eg.db,
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) {
      cat("  富集分析错误:", e$message, "\n")
      return(NULL)
    })

    if (!is.null(ego) && nrow(ego@result) > 0) {
      results[[ont]] <- ego

      write.csv(
        ego@result,
        file.path(DIR_DEG, "GO_enrichment", sprintf("GO_%s_%s_%s.csv", ont, dataset_name, direction)),
        row.names = FALSE
      )

      # dotplot
      cat(sprintf("  - 生成可视化 (%s)...\n", ont))
      p1 <- dotplot(ego, showCategory = 15) +
        ggtitle(sprintf("%s - %s (%s)", dataset_name, direction, ont)) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 10, hjust = 0.5))
      ggsave(file.path(DIR_DEG, "GO_enrichment", sprintf("%s_%s_%s_dotplot.pdf", dataset_name, direction, ont)),
             plot = p1, width = 10, height = 8, limitsize = FALSE)

      # barplot
      p2 <- barplot(ego, showCategory = 15) +
        ggtitle(sprintf("%s - %s (%s)", dataset_name, direction, ont)) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 10, hjust = 0.5))
      ggsave(file.path(DIR_DEG, "GO_enrichment", sprintf("%s_%s_%s_barplot.pdf", dataset_name, direction, ont)),
             plot = p2, width = 10, height = 8, limitsize = FALSE)

      # BP结果绘制GO条目聚类树状网络图
      if (ont == "BP" && nrow(ego@result) >= 5) {
        tryCatch({
          n_cluster <- 5
          cns_colors <- c(
            "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
            "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
          )
          cluster_colors <- cns_colors[seq_len(n_cluster)]

          p3 <- enrichplot::pairwise_termsim(ego)
          if (!is.null(p3) && !inherits(p3, "try-error")) {
            p4 <- treeplot(
              p3,
              showCategory = 20,
              nCluster = n_cluster,
              cluster_method = "ward.D",
              group_color = cluster_colors,
              color = "p.adjust",
              hilight = TRUE,
              fontsize_tiplab = 3.5,
              fontsize_cladelab = 4,
              label_format = 35,
              tiplab_offset = 0.2,
              cladelab_offset = 1
            ) +
              ggtitle(sprintf("%s - %s (%s) GO Network", dataset_name, direction, ont)) +
              theme_minimal() +
              theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
          }
          ggsave(file.path(DIR_DEG, "GO_enrichment", sprintf("%s_%s_%s_network.pdf", dataset_name, direction, ont)),
                 plot = p4, width = 14, height = 10, limitsize = FALSE)
        }, error = function(e) {
          cat("  网络图生成失败:", e$message, "\n")
        })
      }
    }
  }
  return(results)
}

# 10.2 从GO结果中抽取鳞状分化相关通路
analyze_squamous_pathways <- function(go_result, dataset_name, direction) {
  if (is.null(go_result) || nrow(go_result@result) == 0) return(NULL)

  squamous <- go_result@result %>%
    filter(grepl(paste(squamous_keywords, collapse = "|"), Description, ignore.case = TRUE))

  if (nrow(squamous) > 0) {
    cat(sprintf("\n[%s] %s中鳞状分化相关通路 (%d个):\n", dataset_name, direction, nrow(squamous)))
    print(squamous[, c("ID", "Description", "p.adjust", "Count")])

    write.csv(
      squamous,
      file.path(DIR_DEG, "GO_enrichment", sprintf("squamous_%s_%s.csv", dataset_name, direction)),
      row.names = FALSE
    )

    if (nrow(squamous) >= 3) {
      show_n <- min(nrow(squamous), 30)
      p <- dotplot(go_result, showCategory = show_n, color = "p.adjust") +
        ggtitle(sprintf("Squamous Differentiation - %s (%s)", dataset_name, direction)) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 9), plot.title = element_text(size = 11, hjust = 0.5))
      ggsave(file.path(DIR_DEG, "GO_enrichment", sprintf("squamous_%s_%s_dotplot.pdf", dataset_name, direction)),
             plot = p, width = 10, height = max(6, min(show_n, 20)), limitsize = FALSE)
    }
  }
  return(squamous)
}

# 10.3 Trp63相关基因网络分析: 基于SCENIC regulons计算DEG与TF靶基因重叠
analyze_trp63_network <- function(up_genes, down_genes, dataset_name, top_n_tf = 5) {
  cat(sprintf("\n[%s] Trp63相关基因分析...\n", dataset_name))

  regulons_file <- file.path(DIR_OBJ, "SCENIC_trajectory", "int", "2.6_regulons_asGeneSet.Rds")
  if (!file.exists(regulons_file)) {
    cat("  警告: regulons文件不存在，跳过网络分析\n")
    return(NULL)
  }

  regulons <- readRDS(regulons_file)
  cat(sprintf("  Loaded %d regulons\n", length(regulons)))

  all_tf_names <- gsub("_extended", "", gsub("[0-9]+$", "", names(regulons)))
  all_tf_names <- unique(all_tf_names)
  cat(sprintf("  共 %d 个唯一转录因子\n", length(all_tf_names)))

  # 按重叠比例(overlap/targets)寻找Top TF
  find_top_tfs_by_overlap <- function(degs, regulons, direction, top_n = top_n_tf) {
    tf_overlap_list <- list()

    for (tf_name in names(regulons)) {
      tf_clean <- gsub("_extended", "", gsub("[0-9]+$", "", tf_name))
      tf_targets <- regulons[[tf_name]]
      common_genes <- intersect(tf_targets, degs)
      if (length(common_genes) > 0) {
        tf_overlap_list[[tf_clean]] <- list(
          tf = tf_clean,
          original_name = tf_name,
          n_overlap = length(common_genes),
          overlap_genes = common_genes,
          n_targets = length(tf_targets)
        )
      }
    }

    if (length(tf_overlap_list) == 0) {
      return(NULL)
    }

    overlap_df <- do.call(rbind, lapply(names(tf_overlap_list), function(tf) {
      data.frame(
        TF = tf,
        OriginalName = tf_overlap_list[[tf]]$original_name,
        N_Overlap = tf_overlap_list[[tf]]$n_overlap,
        N_Targets = tf_overlap_list[[tf]]$n_targets,
        Overlap_Ratio = tf_overlap_list[[tf]]$n_overlap / tf_overlap_list[[tf]]$n_targets,
        Overlap_Genes = paste(tf_overlap_list[[tf]]$overlap_genes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))

    overlap_df <- overlap_df %>% arrange(desc(Overlap_Ratio))

    if (nrow(overlap_df) > top_n) {
      overlap_df <- overlap_df[1:top_n, ]
    }

    cat(sprintf("  [%s] Top %d TFs by Overlap Ratio (overlap/targets):\n", direction, nrow(overlap_df)))
    for (i in 1:nrow(overlap_df)) {
      cat(sprintf("    %d. %s: %d/%d (%.2f%%) - %s\n",
                  i, overlap_df$TF[i], overlap_df$N_Overlap[i], overlap_df$N_Targets[i],
                  overlap_df$Overlap_Ratio[i] * 100,
                  substr(overlap_df$Overlap_Genes[i], 1, 40)))
    }

    return(list(tf_info = overlap_df, tf_overlap_list = tf_overlap_list))
  }

  cat("\n  分析上调基因中的TF富集...\n")
  up_tf_result <- find_top_tfs_by_overlap(up_genes, regulons, "upregulated")

  cat("\n  分析下调基因中的TF富集...\n")
  down_tf_result <- find_top_tfs_by_overlap(down_genes, regulons, "downregulated")

  trp63_name <- grep("Trp63", names(regulons), value = TRUE)[1]
  trp63_targets <- if (!is.na(trp63_name) && length(trp63_name) > 0) regulons[[trp63_name]] else character()

  # 构建TF->DEG网络边
  build_tf_network_edges <- function(tf_overlap_list, selected_tfs, degs) {
    edges <- data.frame(Source = character(), Target = character())

    for (tf in selected_tfs) {
      tf_entry <- tf_overlap_list[[tf]]
      if (!is.null(tf_entry)) {
        tf_original_name <- tf_entry$original_name
        if (tf_original_name %in% names(regulons)) {
          tf_targets <- regulons[[tf_original_name]]
          common <- intersect(tf_targets, degs)
          if (length(common) > 0) {
            temp_df <- data.frame(Source = rep(tf, length(common)), Target = common)
            edges <- rbind(edges, temp_df)
          }
        }
      }
    }
    return(edges)
  }

  if (!is.null(up_tf_result)) {
    top_up_tfs <- up_tf_result$tf_info$TF
    up_edges <- build_tf_network_edges(up_tf_result$tf_overlap_list, top_up_tfs, up_genes)

    if (nrow(up_edges) > 0) {
      write.csv(up_edges, file.path(DIR_DEG, "network", sprintf("tf_network_%s_up.csv", dataset_name)), row.names = FALSE)
      cat(sprintf("  上调TF网络边数: %d\n", nrow(up_edges)))
    }

    write.csv(up_tf_result$tf_info,
              file.path(DIR_DEG, "network", sprintf("top_tfs_%s_up.csv", dataset_name)), row.names = FALSE)
  }

  if (!is.null(down_tf_result)) {
    top_down_tfs <- down_tf_result$tf_info$TF
    down_edges <- build_tf_network_edges(down_tf_result$tf_overlap_list, top_down_tfs, down_genes)

    if (nrow(down_edges) > 0) {
      write.csv(down_edges, file.path(DIR_DEG, "network", sprintf("tf_network_%s_down.csv", dataset_name)), row.names = FALSE)
      cat(sprintf("  下调TF网络边数: %d\n", nrow(down_edges)))
    }

    write.csv(down_tf_result$tf_info,
              file.path(DIR_DEG, "network", sprintf("top_tfs_%s_down.csv", dataset_name)), row.names = FALSE)
  }

  trp63_overlap_df <- data.frame(
    Dataset = dataset_name,
    Trp63_Targets_Total = length(trp63_targets),
    Up_TFs = ifelse(!is.null(up_tf_result), nrow(up_tf_result$tf_info), 0),
    Down_TFs = ifelse(!is.null(down_tf_result), nrow(down_tf_result$tf_info), 0),
    Top_Up_TFs = ifelse(!is.null(up_tf_result), paste(up_tf_result$tf_info$TF, collapse = ";"), ""),
    Top_Down_TFs = ifelse(!is.null(down_tf_result), paste(down_tf_result$tf_info$TF, collapse = ";"), ""),
    stringsAsFactors = FALSE
  )
  write.csv(trp63_overlap_df,
            file.path(DIR_DEG, "network", sprintf("tf_summary_%s.csv", dataset_name)), row.names = FALSE)

  cat("\n  构建Trp63上下游调控网络...\n")

  # Trp63上下游调控网络: 上游(调控Trp63的TF)+下游(Trp63靶基因与DEG关系)
  trp63_regulatory_network <- function(regulons, trp63_name, up_genes, down_genes, all_degs) {
    trp63_clean <- gsub("_extended", "", gsub("[0-9]+$", "", trp63_name))

    trp63_targets <- regulons[[trp63_name]]
    cat(sprintf("    Trp63直接靶基因数: %d\n", length(trp63_targets)))

    # 上游: 靶基因包含Trp63的TF
    upstream_tfs <- list()
    for (tf_name in names(regulons)) {
      tf_clean <- gsub("_extended", "", gsub("[0-9]+$", "", tf_name))
      if (tf_clean == trp63_clean) next
      tf_targets <- regulons[[tf_name]]
      if ("Trp63" %in% tf_targets || "Trp63" %in% gsub("_extended", "", gsub("[0-9]+$", "", tf_targets))) {
        upstream_tfs[[tf_clean]] <- list(
          tf = tf_clean,
          original_name = tf_name,
          n_targets = length(tf_targets)
        )
      }
    }

    cat(sprintf("    上游调控Trp63的TF数: %d\n", length(upstream_tfs)))

    edges_all <- data.frame(Source = character(), Target = character(), Type = character(), Regulation = character())
    edges_up <- data.frame(Source = character(), Target = character(), Type = character(), Regulation = character())
    edges_down <- data.frame(Source = character(), Target = character(), Type = character(), Regulation = character())

    for (tf_name in names(upstream_tfs)) {
      tf_original <- upstream_tfs[[tf_name]]$original_name
      tf_targets <- regulons[[tf_original]]
      if ("Trp63" %in% tf_targets) {
        temp_df <- data.frame(Source = tf_name, Target = "Trp63", Type = "activates", Regulation = "upstream", stringsAsFactors = FALSE)
        edges_all <- rbind(edges_all, temp_df)
        edges_up <- rbind(edges_up, temp_df)
        edges_down <- rbind(edges_down, temp_df)
      }
    }

    # 下游: Trp63靶基因按上下调归类
    for (target in trp63_targets) {
      target_in_up <- target %in% up_genes
      target_in_down <- target %in% down_genes

      if (target_in_up) {
        temp_df <- data.frame(Source = "Trp63", Target = target, Type = "DEG_up", Regulation = "up", stringsAsFactors = FALSE)
        edges_all <- rbind(edges_all, temp_df)
        edges_up <- rbind(edges_up, temp_df)
      } else if (target_in_down) {
        temp_df <- data.frame(Source = "Trp63", Target = target, Type = "DEG_down", Regulation = "down", stringsAsFactors = FALSE)
        edges_all <- rbind(edges_all, temp_df)
        edges_down <- rbind(edges_down, temp_df)
      } else {
        temp_df <- data.frame(Source = "Trp63", Target = target, Type = "target", Regulation = "other", stringsAsFactors = FALSE)
        edges_all <- rbind(edges_all, temp_df)
      }
    }

    cat(sprintf("    网络总边数: %d (上游: %d, 下游上调: %d, 下游下调: %d, 下游其他: %d)\n",
                nrow(edges_all),
                sum(edges_all$Regulation == "upstream"),
                sum(edges_all$Regulation == "up"),
                sum(edges_all$Regulation == "down"),
                sum(edges_all$Regulation == "other")))

    return(list(
      edges_all = edges_all,
      edges_up = edges_up,
      edges_down = edges_down,
      upstream_tfs = upstream_tfs,
      trp63_targets = trp63_targets,
      trp63_clean = trp63_clean
    ))
  }

  trp63_net_result <- trp63_regulatory_network(regulons, trp63_name, up_genes, down_genes, c(up_genes, down_genes))

  if (nrow(trp63_net_result$edges_all) > 0) {
    write.csv(trp63_net_result$edges_all,
              file.path(DIR_DEG, "network", sprintf("trp63_regulatory_network_all_%s.csv", dataset_name)), row.names = FALSE)
    cat(sprintf("    已保存: trp63_regulatory_network_all_%s.csv\n", dataset_name))
  }

  if (nrow(trp63_net_result$edges_up) > 0) {
    write.csv(trp63_net_result$edges_up,
              file.path(DIR_DEG, "network", sprintf("trp63_regulatory_network_up_%s.csv", dataset_name)), row.names = FALSE)
    cat(sprintf("    已保存: trp63_regulatory_network_up_%s.csv\n", dataset_name))
  }

  if (nrow(trp63_net_result$edges_down) > 0) {
    write.csv(trp63_net_result$edges_down,
              file.path(DIR_DEG, "network", sprintf("trp63_regulatory_network_down_%s.csv", dataset_name)), row.names = FALSE)
    cat(sprintf("    已保存: trp63_regulatory_network_down_%s.csv\n", dataset_name))
  }

  upstream_tfs_df <- if (length(trp63_net_result$upstream_tfs) > 0) {
    do.call(rbind, lapply(names(trp63_net_result$upstream_tfs), function(tf) {
      data.frame(
        TF = tf,
        OriginalName = trp63_net_result$upstream_tfs[[tf]]$original_name,
        N_Targets = trp63_net_result$upstream_tfs[[tf]]$n_targets,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(TF = character(), OriginalName = character(), N_Targets = integer(), stringsAsFactors = FALSE)
  }
  write.csv(upstream_tfs_df,
            file.path(DIR_DEG, "network", sprintf("trp63_upstream_tfs_%s.csv", dataset_name)), row.names = FALSE)
  cat(sprintf("    已保存: trp63_upstream_tfs_%s.csv\n", dataset_name))

  return(list(
    up_tf_result = up_tf_result,
    down_tf_result = down_tf_result,
    trp63_targets = trp63_targets,
    trp63_network = trp63_net_result
  ))
}

# 10.4 批量运行: 逐数据集做GO富集+鳞状通路抽取+Trp63网络分析
all_results <- list()

for (dataset_name in names(deg_datasets)) {
  deg_file <- deg_datasets[[dataset_name]]

  if (!file.exists(deg_file)) {
    cat(sprintf("\n警告: 文件不存在 %s，跳过\n", deg_file))
    next
  }

  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("处理数据集: %s\n", dataset_name))
  cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))

  deg_data <- read.csv(deg_file)
  cat(sprintf("DEG数据加载完成，共 %d 个基因\n", nrow(deg_data)))

  # 按log2FC方向+p_adj<=0.05划分上下调基因
  deg_data <- deg_data %>% mutate(change = ifelse(avg_log2FC > 0, "up", "down"))
  up_genes <- deg_data %>% filter(change == "up", p_val_adj <= 0.05) %>% pull(gene)
  down_genes <- deg_data %>% filter(change == "down", p_val_adj <= 0.05) %>% pull(gene)

  cat(sprintf("  上调基因数: %d\n", length(up_genes)))
  cat(sprintf("  下调基因数: %d\n", length(down_genes)))

  go_results <- list()

  go_up <- run_go_enrichment_diverse(up_genes, dataset_name, "upregulated")
  if (!is.null(go_up)) {
    go_results$up <- go_up
    if (!is.null(go_up$BP)) {
      analyze_squamous_pathways(go_up$BP, dataset_name, "upregulated")
    }
  }

  go_down <- run_go_enrichment_diverse(down_genes, dataset_name, "downregulated")
  if (!is.null(go_down)) {
    go_results$down <- go_down
    if (!is.null(go_down$BP)) {
      analyze_squamous_pathways(go_down$BP, dataset_name, "downregulated")
    }
  }

  all_results[[dataset_name]] <- list(go_results = go_results, up_genes = up_genes, down_genes = down_genes)
  all_results[[dataset_name]]$trp63_result <- analyze_trp63_network(up_genes, down_genes, dataset_name)
}

# 富集分析总结输出
cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("富集分析总结\n")
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))

for (dataset_name in names(all_results)) {
  cat(sprintf("\n[%s]\n", dataset_name))
  res <- all_results[[dataset_name]]

  for (direction in c("up", "down")) {
    if (!is.null(res$go_results[[direction]]$BP) && nrow(res$go_results[[direction]]$BP@result) > 0) {
      cat(sprintf("  %s基因Top 5 GO BP:\n", direction))
      top5 <- res$go_results[[direction]]$BP@result$Description[1:5]
      cat(sprintf("    %s\n", paste(top5, collapse = "\n    ")))
    }
  }

  if (!is.null(res$trp63_result)) {
    cat("  Top TFs by DEG overlap:\n")
    if (!is.null(res$trp63_result$up_tf_result)) {
      cat(sprintf("    上调中: %s\n", paste(res$trp63_result$up_tf_result$tf_info$TF, collapse = ", ")))
    }
    if (!is.null(res$trp63_result$down_tf_result)) {
      cat(sprintf("    下调中: %s\n", paste(res$trp63_result$down_tf_result$tf_info$TF, collapse = ", ")))
    }
  }
}

cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("分析完成！\n")
cat(sprintf("输出文件位置: %s\n", DIR_DEG))
cat("  - GO_enrichment/: GO富集结果和可视化\n")
cat("  - network/: TF基因网络图\n")

# ---- 11. 整合对象marker分析、组织学比较与dnb状态DEG ----
# 11.1 各细胞类型标志基因DotPlot(按seurat_clusters分组)
markers <- list(
  "AT2" = c("Sftpc", "Sftpa1", "Sftpd", "Etv5"),
  "Basal" = c("Krt15", "Krt5", "Krt17", "Trp63"),
  "Ciliated" = c("Pifo", "Foxj1", "Ccdc78"),
  "Club" = c("Scgb3a2", "Scgb1a1", "Cyp2f2"),
  "Goblet" = c("Muc2", "Tff2", "Muc5b", "Spdef"),
  "Proliferative" = c("Birc5", "Top2a", "Mki67"),
  "Endothelium" = c("Pecam1", "Cdh5", "Cldn5", ""),
  "Fibroblast" = c("Col1a2", "Col1a1", "Vim", "Pdgfra"),
  "Malignant" = c("Kng2", "Meg3", "Dlk1", "Gkn2", "Lcn2")
)

dot_plot <- DotPlot(se_merged, features = markers, group.by = "seurat_clusters") +
  scale_color_gradient2(low = "lightgray", mid = "lightblue", high = "purple", midpoint = 0) +
  theme(
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)
  ) +
  labs(title = "Gene Expression Dot Plot")
print(dot_plot)

# 11.2 对RNA assay重新SCT标准化后做全聚类marker分析
DefaultAssay(se_merged) <- "RNA"
se_merged <- SCTransform(se_merged, verbose = FALSE, method = "glmGamPoi")

# 恢复原有聚类分组
Idents(se_merged) <- se_merged$seurat_clusters
se_merged <- PrepSCTFindMarkers(se_merged)
all_markers <- FindAllMarkers(se_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, file = file.path(DIR_OBJ_MERGED, "merged_obj_all_seurat_clusters_markers.csv"))

# 11.3 组织学(histology)间差异比较
deg_ast_vs_luad <- FindMarkers(se_merged, group.by = "histology", ident.1 = "AST", ident.2 = "LUAD")
write.csv(deg_ast_vs_luad, file = file.path(DIR_OBJ_MERGED, "merged_obj_deg_ast_vs_luad.csv"))
deg_asy_vs_scc <- FindMarkers(se_merged, group.by = "histology", ident.1 = "ASY", ident.2 = "SCC")
write.csv(deg_asy_vs_scc, file = file.path(DIR_OBJ_MERGED, "merged_obj_deg_asy_vs_scc.csv"))

# 11.4 将轨迹对象的dnb_state按"_1"后缀映射到整合对象并做DEG
# 重新读取已回存的轨迹对象: 其dnb_state为字符型(保存于因子转换之前),
# 字符型保证下方case_when中"dnb_state > 0"判断可用(所有非NA字符均>"0")
seurat_obj <- readRDS(file.path(DIR_OBJ_SEU, "se_trajectory_subset.rds"))
meta_df <- seurat_obj@meta.data[, "dnb_state"] |> as.data.frame()
rownames(meta_df) <- paste0(rownames(seurat_obj@meta.data), "_1")
se_merged <- AddMetaData(se_merged, metadata = meta_df, col.name = "dnb_state")

# 合并分组: 有dnb_state的用状态, 无则按histology归入LUAD/SCC
se_merged@meta.data <- se_merged@meta.data %>%
  mutate(
    dnb_merged = case_when(
      dnb_state > 0 ~ as.character(dnb_state),
      is.na(dnb_state) & histology == "LUAD" ~ "LUAD",
      is.na(dnb_state) & histology == "SCC" ~ "SCC",
      TRUE ~ "Other"
    )
  )
table(se_merged$dnb_merged)

deg_tp_vs_luad <- FindMarkers(se_merged, group.by = "dnb_merged", ident.1 = "Tipping_Point", ident.2 = "LUAD")
write.csv(deg_tp_vs_luad, file = file.path(DIR_OBJ_MERGED, "merged_obj_deg_tp_vs_luad.csv"))
deg_tp_vs_scc <- FindMarkers(se_merged, group.by = "dnb_merged", ident.1 = "Tipping_Point", ident.2 = "SCC")
write.csv(deg_tp_vs_scc, file = file.path(DIR_OBJ_MERGED, "merged_obj_deg_tp_vs_scc.csv"))

# 11.5 每个cluster的top marker DotPlot与组织学比例图
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  ungroup()

DotPlot(
  se_merged,
  features = unique(top5_markers$gene),
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plotProportion(se_merged, meta.2 = "histology")

# ---- 12. Nrg1表达趋势与核心TF靶基因调控分析 ----
set.seed(1)

se_traj <- readRDS(file.path(DIR_OBJ_SEU, "se_trajectory_subset.rds"))
regulons <- readRDS(file.path(DIR_OBJ, "SCENIC_trajectory", "int", "2.6_regulons_asGeneSet.Rds"))
cat(sprintf("Loaded %d regulons\n", length(regulons)))

# 12.1 Nrg1沿AST轨迹的表达趋势(loess拟合)
gene_name <- "Nrg1"
tf_expr <- FetchData(se_traj, vars = c("Nrg1", "tdTomato_ratio", "dnb_state"))
colnames(tf_expr)[1] <- "expression"

p <- ggplot(tf_expr, aes(x = tdTomato_ratio, y = expression)) +
  geom_smooth(method = "loess", span = 0.3, colour = "red") +
  theme_classic() +
  labs(title = sprintf("%s Gene Expression Along AST", gene_name),
       x = "tdTomato_ratio", y = "Expression")
pdf(file.path(DIR_FIG_TF, "Nrg1_change_trend.pdf"), width = 8, height = 6)
print(p)
dev.off()

# CNS期刊风格版本
output_filename <- paste0(gene_name, "_expression_trend.pdf")

tf_expr <- FetchData(
  object = se_traj,
  vars = c(gene_name, "tdTomato_ratio", "dnb_state")
)
colnames(tf_expr)[1] <- "expression"

p <- ggplot(data = tf_expr, aes(x = tdTomato_ratio, y = expression)) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    color = "#E63946",
    fill = "#F8C8C8",
    linewidth = 1.1,
    se = TRUE
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_blank(),
    axis.title.x = element_text(family = "sans", size = 12),
    axis.title.y = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans", size = 10),
    plot.title = element_blank()
  ) +
  labs(
    x = "tdTomato signal ratio",
    y = paste0(gene_name, " normalized expression")
  )

pdf(
  file = file.path(DIR_FIG_TF, output_filename),
  width = 8,
  height = 6,
  useDingbats = FALSE
)
print(p)
dev.off()

# 12.2 提取核心TF的SCENIC靶基因
target_tfs <- c("Nkx2-1", "Foxa2", "Elf5", "Nfe2l2", "Trp63", "Sox2", "Sox21", "Elf5", "Zeb1", "Ezh2", "Foxp2")

cat("\nExtracting TF target genes...\n")
tf_targets <- list()
for (tf in target_tfs) {
  tf_idx <- grep(tf, names(regulons), value = TRUE)
  if (length(tf_idx) > 0) {
    tf_targets[[tf]] <- regulons[[tf_idx[1]]]
    cat(sprintf("%s: %d targets\n", tf, length(tf_targets[[tf]])))
  } else {
    cat(sprintf("%s: not found in regulons\n", tf))
  }
}

# 12.3 dnb_state间DEG(用于与TF靶基因取交集)
cat("\nFinding DEGs between states...\n")
Idents(se_traj) <- "dnb_state"
states <- c("Pre_Tipping", "Tipping_Point", "Post_Tipping")
se_states <- subset(se_traj, idents = states)
Idents(se_states) <- "dnb_state"

all_degs <- list()

de_pre_tip <- FindMarkers(se_states, ident.1 = "Tipping_Point", ident.2 = "Pre_Tipping", test.use = "wilcox")
de_pre_tip$gene <- rownames(de_pre_tip)
de_pre_tip$comparison <- "Tipping_vs_Pre"
all_degs[["Tipping_vs_Pre"]] <- de_pre_tip

de_tip_post <- FindMarkers(se_states, ident.1 = "Post_Tipping", ident.2 = "Tipping_Point", test.use = "wilcox")
de_tip_post$gene <- rownames(de_tip_post)
de_tip_post$comparison <- "Post_vs_Tipping"
all_degs[["Post_vs_Tipping"]] <- de_tip_post

de_pre_post <- FindMarkers(se_states, ident.1 = "Post_Tipping", ident.2 = "Pre_Tipping", test.use = "wilcox")
de_pre_post$gene <- rownames(de_pre_post)
de_pre_post$comparison <- "Post_vs_Pre"
all_degs[["Post_vs_Pre"]] <- de_pre_post

# 12.4 TF靶基因与各比较DEG的交集统计
cat("\nAnalyzing TF target overlaps with DEGs...\n")
tf_deg_overlap <- list()

for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  overlap_summary <- data.frame()

  for (comp in names(all_degs)) {
    deg_df <- all_degs[[comp]]
    sig_up <- deg_df %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
    sig_down <- deg_df %>% filter(p_val_adj < 0.05, avg_log2FC < 0) %>% pull(gene)

    overlap_up <- intersect(tf_genes, sig_up)
    overlap_down <- intersect(tf_genes, sig_down)

    overlap_summary <- rbind(overlap_summary, data.frame(
      TF = tf,
      Comparison = comp,
      TF_Targets = length(tf_genes),
      DEG_Up = length(sig_up),
      DEG_Down = length(sig_down),
      Overlap_Up = length(overlap_up),
      Overlap_Down = length(overlap_down),
      Overlap_Up_Genes = paste(overlap_up, collapse = ","),
      Overlap_Down_Genes = paste(overlap_down, collapse = ",")
    ))
  }
  tf_deg_overlap[[tf]] <- overlap_summary
}

overlap_all <- do.call(rbind, tf_deg_overlap)
write.csv(overlap_all, file.path(DIR_OBJ_TF, "TF_DEG_overlap.csv"), row.names = FALSE)

# 12.5 核心TF在各dnb_state的平均表达热图
cat("\nPlotting TF expression heatmap...\n")
tf_expr <- AverageExpression(se_states, features = target_tfs, group.by = "dnb_state", assays = "RNA")
tf_expr_mat <- as.matrix(tf_expr$RNA)

pdf(file.path(DIR_FIG_TF, "TF_expression_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  tf_expr_mat,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "TF Expression across dnb_states"
)
dev.off()

# 12.6 各TF靶基因(取top50高变异)在各状态的平均表达热图
cat("\nPlotting TF target gene expression...\n")
expr_mat_global <- GetAssayData(se_states, assay = "RNA", layer = "data")
for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  tf_genes_in_se <- tf_genes[tf_genes %in% rownames(se_states)]

  if (length(tf_genes_in_se) >= 5) {
    gene_vars <- apply(expr_mat_global[tf_genes_in_se, ], 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(tf_genes_in_se))]

    avg_expr <- AverageExpression(se_states, features = top_genes, group.by = "dnb_state", assays = "RNA")
    expr_mat <- as.matrix(avg_expr$RNA)

    pdf(file.path(DIR_FIG_TF, sprintf("%s_targets_heatmap.pdf", tf)), width = 10, height = 12)
    pheatmap(
      expr_mat,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = sprintf("%s Target Gene Expression", tf),
      fontsize = 11,
      fontsize_row = 9,
      fontsize_col = 10,
      border_color = "gray90",
      treeheight_row = 25,
      legend_breaks = seq(-2, 2, 1),
      legend_labels = c("Low", "-1", "0", "1", "High"),
      angle_col = 45,
      show_rownames = TRUE
    )
    dev.off()
  }
}

# 12.7 TF靶基因的GO富集
cat("\nPerforming GO enrichment for TF targets...\n")
for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  tryCatch({
    gene.df <- bitr(tf_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

    if (nrow(gene.df) > 0) {
      for (ont in c("BP", "CC", "MF")) {
        ego <- enrichGO(
          gene = unique(gene.df$ENTREZID),
          OrgDb = org.Mm.eg.db,
          keyType = "ENTREZID",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )

        if (nrow(ego@result) > 0) {
          write.csv(ego@result, file.path(DIR_OBJ_TF, sprintf("%s_GO_%s.csv", tf, ont)))

          show_n <- ifelse(nrow(ego@result) >= 10, 10, nrow(ego@result))
          p <- dotplot(ego, showCategory = show_n) + ggtitle(sprintf("%s GO %s", tf, ont))
          pdf(file.path(DIR_FIG_TF, sprintf("%s_GO_%s_dotplot.pdf", tf, ont)), width = 10, height = 8)
          print(p)
          dev.off()
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("Error processing %s: %s\n", tf, e$message))
  })
}

# 12.8 TF及其靶基因沿tdRatio的表达趋势图
cat("\nPlotting TF and target expression trends...\n")
se_with_ratio <- subset(se_traj, !is.na(tdTomato_ratio))
expr_mat_ratio <- GetAssayData(se_with_ratio, assay = "RNA", layer = "data")

# 绘制单个TF及其靶基因沿tdRatio的loess趋势(分面)
plot_tf_trend <- function(tf, target_genes, se_obj) {
  tf_expr <- FetchData(se_obj, vars = c(tf, "tdTomato_ratio", "dnb_state"))
  tf_expr$gene <- tf
  tf_expr$type <- "TF"
  colnames(tf_expr)[1] <- "expression"

  target_expr_list <- list()
  for (gene in target_genes) {
    if (gene %in% rownames(se_obj)) {
      expr <- FetchData(se_obj, vars = c(gene, "tdTomato_ratio", "dnb_state"))
      expr$gene <- gene
      expr$type <- "Target"
      colnames(expr)[1] <- "expression"
      target_expr_list[[gene]] <- expr
    }
  }

  if (length(target_expr_list) > 0) {
    all_expr <- rbind(tf_expr, do.call(rbind, target_expr_list))

    p <- ggplot(all_expr, aes(x = tdTomato_ratio, y = expression, color = type)) +
      geom_smooth(method = "loess", span = 0.3) +
      facet_wrap(~gene, scales = "free_y") +
      theme_bw() +
      labs(title = sprintf("%s and Target Genes Expression", tf),
           x = "tdTomato_ratio", y = "Expression")

    return(p)
  }
  return(NULL)
}

# 各TF取top6高变异靶基因绘制趋势图
for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  tf_genes_in_se <- tf_genes[tf_genes %in% rownames(se_with_ratio)]
  if (length(tf_genes_in_se) >= 3) {
    gene_vars <- apply(expr_mat_ratio[tf_genes_in_se, ], 1, var)
    top_targets <- names(sort(gene_vars, decreasing = TRUE))[1:min(6, length(tf_genes_in_se))]

    p <- plot_tf_trend(tf, top_targets, se_with_ratio)
    if (!is.null(p)) {
      pdf(file.path(DIR_FIG_TF, sprintf("%s_trend_plot.pdf", tf)), width = 14, height = 10)
      print(p)
      dev.off()
    }
  }
}

# Trp63全部靶基因趋势图(注意: 原脚本此处tf沿用上一循环的最后一个TF名,
# 图内标题显示该TF名而非Trp63; 原脚本此处if块缺少右括号导致无法解析, 已补全)
tf_genes <- tf_targets[["Trp63"]]
tf_genes_in_se <- tf_genes[tf_genes %in% rownames(se_with_ratio)]
if (length(tf_genes_in_se) >= 3) {
  top_targets <- tf_genes_in_se

  p <- plot_tf_trend(tf, top_targets, se_with_ratio)
  if (!is.null(p)) {
    pdf(file.path(DIR_FIG_TF, "trp63_all_targets_trend_plot.pdf"), width = 14, height = 10)
    print(p)
    dev.off()
  }
}

# 12.9 核心TF网络汇总: 共调控基因(>=2个TF调控)统计与热图
cat("\nCreating core TF network summary...\n")
all_target_genes <- unique(unlist(tf_targets))
tf_target_matrix <- matrix(0, nrow = length(target_tfs), ncol = length(all_target_genes))
rownames(tf_target_matrix) <- target_tfs
colnames(tf_target_matrix) <- all_target_genes

for (tf in names(tf_targets)) {
  tf_target_matrix[tf, tf_targets[[tf]]] <- 1
}

tf_target_overlap <- colSums(tf_target_matrix)
co_regulated_genes <- names(tf_target_overlap[tf_target_overlap >= 2])
cat(sprintf("Genes co-regulated by >=2 TFs: %d\n", length(co_regulated_genes)))

core_network <- list(
  TF_targets = tf_targets,
  Co_regulated_genes = co_regulated_genes,
  TF_target_matrix = tf_target_matrix
)
saveRDS(core_network, file.path(DIR_OBJ_TF, "core_TF_network.rds"))

if (length(co_regulated_genes) >= 5) {
  co_reg_matrix <- tf_target_matrix[, co_regulated_genes]
  co_reg_matrix <- co_reg_matrix[rowSums(co_reg_matrix) > 0, ]

  pdf(file.path(DIR_FIG_TF, "co_regulated_genes_heatmap.pdf"), width = 14, height = 10)
  pheatmap(
    co_reg_matrix,
    color = colorRampPalette(c("white", "darkblue"))(2),
    breaks = c(-0.5, 0.5, 1.5),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Co-regulated Genes by Core TFs"
  )
  dev.off()
}

cat("\n=== Analysis complete! ===\n")
cat(sprintf("Results saved to %s and %s\n", DIR_FIG_TF, DIR_OBJ_TF))

# ---- 13. 动态核心TF-靶基因网络图 (LUAD -> TP -> SCC) ----
# 本章原为独立脚本: 先清空环境再重建路径(依赖包已在第0章加载, 不受rm影响)
rm(list = ls())
set.seed(123)

DIR_ROOT <- getwd()
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures", "deg_network_core")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

REGULONS_FILE <- file.path(DIR_RES, "object", "SCENIC_trajectory", "int", "2.6_regulons_asGeneSet.Rds")
DEG_LUAD_FILE <- file.path(DIR_RES, "object", "seurat", "merged_obj", "merged_obj_deg_tp_vs_luad.csv")
DEG_SCC_FILE <- file.path(DIR_RES, "object", "seurat", "merged_obj", "merged_obj_deg_tp_vs_scc.csv")
OUTPUT_PDF <- file.path(DIR_FIG, "dynamic_core_TF_network_LUAD_TP_SCC.pdf")

# 参数: 核心TF/显著性阈值/每TF展示靶基因数/regulon类别
CORE_TFS <- c("Trp63", "Sox2")
PADJ_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.25
TOP_TARGETS_PER_TF <- 30
INCLUDE_EXTENDED_REGULONS <- TRUE

# 出版风格色盲友好配色:
# 红/蓝编码单调轨迹, 橙/青编码TP特异状态
TRAJECTORY_COLORS <- c(
  Monotonic_up = "#bd5b6a",
  Monotonic_down = "#3568A8",
  TP_peak = "#ddaa96",
  TP_trough = "#4F9D9D",
  Single_contrast = "#dce3e9",
  Not_significant = "#F7F7F7"
)

TRAJECTORY_LABELS <- c(
  Monotonic_up = "Monotonic increase",
  Monotonic_down = "Monotonic decrease",
  TP_peak = "Transition-state peak",
  TP_trough = "Transition-state trough",
  Single_contrast = "Significant in one contrast",
  Not_significant = "Not significant"
)

# 图中仅展示以下轨迹模式的靶基因
PLOTTED_TRAJECTORIES <- c("Monotonic_up", "TP_peak", "Single_contrast")

# 边颜色区分两个TF模块(图例隐藏, 由hub节点标签标识TF身份)
TF_EDGE_COLORS <- c(
  Trp63 = "#8A4B6E",
  Sox2 = "#3F7188"
)

PLOT_BACKGROUND <- "#FCFCFD"
LABEL_COLOR <- "#202124"

# 13.1 工具函数: 文件与列检查
check_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  invisible(TRUE)
}

check_required_columns <- function(df, required, object_name) {
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

# 标准化基因列名(兼容gene/Gene/symbol/X等)
standardize_gene_column <- function(df) {
  candidate_columns <- c("gene", "Gene", "symbol", "Symbol", "feature", "features", "X")

  matched <- candidate_columns[candidate_columns %in% colnames(df)]

  if (length(matched) == 0) {
    stop(
      "No gene column found. Expected one of: ",
      paste(candidate_columns, collapse = ", ")
    )
  }

  df$gene <- as.character(df[[matched[1]]])
  df
}

# DEG数据清洗: 统一列名/去重/标记显著性
prepare_deg <- function(deg_df, contrast_name, padj_cutoff = 0.05, logfc_cutoff = 0.25) {
  deg_df <- standardize_gene_column(deg_df)

  check_required_columns(
    deg_df,
    c("gene", "avg_log2FC", "p_val_adj"),
    contrast_name
  )

  deg_df %>%
    transmute(
      gene = as.character(gene),
      logFC = as.numeric(avg_log2FC),
      padj = as.numeric(p_val_adj)
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(logFC),
      !is.na(padj)
    ) %>%
    group_by(gene) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      significant = padj < padj_cutoff & abs(logFC) >= logfc_cutoff
    )
}

# 13.2 regulon名称解析: 提取TF名/调控方向/直接或扩展类别
parse_regulon_tf <- function(regulon_name) {
  regulon_name %>%
    str_remove("\\s*\\([0-9]+g\\)\\s*$") %>%
    str_remove("\\s*\\([+-]\\)\\s*$") %>%
    str_remove("_extended$") %>%
    str_trim()
}

parse_regulon_sign <- function(regulon_name) {
  case_when(
    str_detect(regulon_name, "\\(\\+\\)") ~ "activating",
    str_detect(regulon_name, "\\(-\\)") ~ "repressive",
    str_detect(str_to_lower(regulon_name), "positive|activating|activation") ~ "activating",
    str_detect(str_to_lower(regulon_name), "negative|repressive|repression") ~ "repressive",
    TRUE ~ "unknown"
  )
}

parse_regulon_class <- function(regulon_name) {
  ifelse(
    str_detect(regulon_name, "_extended"),
    "extended",
    "direct"
  )
}

# regulon列表转长表(Source->Target)
regulons_to_long <- function(regulons) {
  if (is.null(names(regulons))) {
    stop("The regulon list has no names.")
  }

  map_dfr(names(regulons), function(regulon_name) {
    targets <- unique(as.character(unlist(regulons[[regulon_name]], use.names = FALSE)))
    targets <- targets[!is.na(targets) & targets != ""]

    tibble(
      Source = parse_regulon_tf(regulon_name),
      Target = targets,
      regulon_name = regulon_name,
      edge_sign = parse_regulon_sign(regulon_name),
      regulon_class = parse_regulon_class(regulon_name)
    )
  }) %>%
    distinct()
}

# 13.3 联合TP-vs-LUAD与TP-vs-SCC两个对比, 划分基因的转换轨迹模式
build_transition_table <- function(deg_luad_clean, deg_scc_clean) {
  full_join(
    deg_luad_clean %>%
      rename(
        logFC_tp_vs_luad = logFC,
        padj_tp_vs_luad = padj,
        sig_luad = significant
      ),
    deg_scc_clean %>%
      rename(
        logFC_tp_vs_scc = logFC,
        padj_tp_vs_scc = padj,
        sig_scc = significant
      ),
    by = "gene"
  ) %>%
    mutate(
      sig_luad = coalesce(sig_luad, FALSE),
      sig_scc = coalesce(sig_scc, FALSE),
      transition_pattern = case_when(
        sig_luad & sig_scc &
          logFC_tp_vs_luad > 0 & logFC_tp_vs_scc < 0 ~ "Monotonic_up",
        sig_luad & sig_scc &
          logFC_tp_vs_luad < 0 & logFC_tp_vs_scc > 0 ~ "Monotonic_down",
        sig_luad & sig_scc &
          logFC_tp_vs_luad > 0 & logFC_tp_vs_scc > 0 ~ "TP_peak",
        sig_luad & sig_scc &
          logFC_tp_vs_luad < 0 & logFC_tp_vs_scc < 0 ~ "TP_trough",
        sig_luad | sig_scc ~ "Single_contrast",
        TRUE ~ "Not_significant"
      ),
      dynamic_score = abs(coalesce(logFC_tp_vs_luad, 0)) + abs(coalesce(logFC_tp_vs_scc, 0)),
      transition_pattern = factor(transition_pattern, levels = names(TRAJECTORY_COLORS))
    )
}

# 13.4 构建核心TF-target边表(按dynamic_score取每TF top靶基因)
build_transition_network <- function(regulon_long, transition_table, core_tfs,
                                     top_targets_per_tf = 30, include_extended_regulons = TRUE) {
  edge_candidates <- regulon_long %>%
    filter(Source %in% core_tfs, Target != Source)

  if (!include_extended_regulons) {
    edge_candidates <- edge_candidates %>%
      filter(regulon_class == "direct")
  }

  # 仅保留有转换信号的靶基因(或核心TF自身)
  edge_candidates <- edge_candidates %>%
    left_join(transition_table, by = c("Target" = "gene")) %>%
    mutate(
      transition_pattern = coalesce(as.character(transition_pattern), "Not_significant"),
      dynamic_score = coalesce(dynamic_score, 0)
    ) %>%
    filter(
      Target %in% core_tfs |
        transition_pattern != "Not_significant"
    )

  if (nrow(edge_candidates) == 0) {
    return(tibble())
  }

  # 汇总多重regulon支持的边
  edge_summary <- edge_candidates %>%
    group_by(Source, Target) %>%
    summarise(
      edge_sign = case_when(
        any(edge_sign == "activating") & !any(edge_sign == "repressive") ~ "activating",
        any(edge_sign == "repressive") & !any(edge_sign == "activating") ~ "repressive",
        TRUE ~ "unknown"
      ),
      regulon_class = ifelse(any(regulon_class == "direct"), "direct", "extended"),
      regulon_evidence = n_distinct(regulon_name),
      transition_pattern = first(transition_pattern),
      logFC_tp_vs_luad = first(logFC_tp_vs_luad),
      logFC_tp_vs_scc = first(logFC_tp_vs_scc),
      dynamic_score = first(dynamic_score),
      .groups = "drop"
    )

  core_edges <- edge_summary %>%
    filter(Target %in% core_tfs)

  target_edges <- edge_summary %>%
    filter(
      !Target %in% core_tfs,
      as.character(transition_pattern) %in% PLOTTED_TRAJECTORIES
    ) %>%
    group_by(Source) %>%
    slice_max(
      order_by = dynamic_score,
      n = top_targets_per_tf,
      with_ties = FALSE
    ) %>%
    ungroup()

  bind_rows(core_edges, target_edges) %>%
    distinct(Source, Target, .keep_all = TRUE) %>%
    mutate(
      edge_sign = factor(edge_sign, levels = c("activating", "repressive", "unknown")),
      regulon_class = factor(regulon_class, levels = c("direct", "extended")),
      transition_pattern = factor(transition_pattern, levels = names(TRAJECTORY_COLORS))
    )
}

# 13.5 构建节点表(核心TF/靶基因+轨迹模式+动态分数)
build_node_table <- function(edges, transition_table, core_tfs) {
  all_nodes <- unique(c(core_tfs, edges$Source, edges$Target))

  tibble(name = all_nodes) %>%
    left_join(
      transition_table %>%
        select(gene, transition_pattern, logFC_tp_vs_luad, logFC_tp_vs_scc, dynamic_score),
      by = c("name" = "gene")
    ) %>%
    mutate(
      is_core_tf = name %in% core_tfs,
      node_type = ifelse(is_core_tf, "Core_TF", "Target"),
      transition_pattern = coalesce(as.character(transition_pattern), "Not_significant"),
      transition_pattern = factor(transition_pattern, levels = names(TRAJECTORY_COLORS)),
      dynamic_score = coalesce(dynamic_score, 0)
    )
}

# 添加网络度统计并为所有节点生成标签
add_graph_statistics <- function(graph_object) {
  V(graph_object)$network_degree <- degree(
    graph_object,
    mode = "all",
    loops = FALSE
  )
  V(graph_object)$plot_label <- V(graph_object)$name
  graph_object
}

# 13.6 手工布局: 弧线坐标生成(靶基因按外向弧排列, 半径交替避免标签重叠)
make_arc_coordinates <- function(node_names, center_x, center_y, radius, angles_deg,
                                 radius_pattern = c(1, 1.10, 0.94)) {
  n_nodes <- length(node_names)

  if (n_nodes == 0) {
    return(tibble(name = character(), x = double(), y = double()))
  }

  if (length(angles_deg) != n_nodes) {
    stop("The number of angles must equal the number of nodes.")
  }

  radii <- radius * rep(radius_pattern, length.out = n_nodes)
  angles_rad <- angles_deg * pi / 180

  tibble(
    name = node_names,
    x = center_x + radii * cos(angles_rad),
    y = center_y + radii * sin(angles_rad)
  )
}

# 手工双hub参考布局: Trp63居左, Sox2居右, 靶基因沿各自外向弧排列
build_reference_layout <- function(graph_object, core_tfs) {
  vertex_table <- igraph::as_data_frame(graph_object, what = "vertices")

  edge_table <- igraph::as_data_frame(graph_object, what = "edges") %>%
    rename(Source = from, Target = to)

  # 两个TF hub对称置于同一水平轴
  if (all(c("Trp63", "Sox2") %in% vertex_table$name)) {
    hub_positions <- tibble(
      name = c("Trp63", "Sox2"),
      x = c(-0.58, 0.58),
      y = c(0, 0)
    )
  } else {
    available_core_tfs <- intersect(core_tfs, vertex_table$name)
    hub_positions <- tibble(
      name = available_core_tfs,
      x = seq(-0.58, 0.58, length.out = max(length(available_core_tfs), 1))[seq_along(available_core_tfs)],
      y = 0
    )
  }

  missing_core_tfs <- setdiff(intersect(core_tfs, vertex_table$name), hub_positions$name)

  if (length(missing_core_tfs) > 0) {
    hub_positions <- bind_rows(
      hub_positions,
      tibble(
        name = missing_core_tfs,
        x = seq(-0.70, 0.70, length.out = length(missing_core_tfs)),
        y = 1.35
      )
    )
  }

  # 按source TF收集靶基因
  target_source_map <- edge_table %>%
    filter(!Target %in% core_tfs) %>%
    distinct(Source, Target)

  order_targets <- function(target_names) {
    vertex_table %>%
      filter(name %in% target_names) %>%
      arrange(desc(dynamic_score), name) %>%
      pull(name)
  }

  # 小模块用单环, 稠密模块用双交错环(保持外向弧形并增大标签间距)
  place_module_targets <- function(target_names, center_x, center_y, side = c("left", "right")) {
    side <- match.arg(side)
    n_targets <- length(target_names)

    if (n_targets == 0) {
      return(tibble(name = character(), x = double(), y = double()))
    }

    angle_limits <- if (side == "left") c(72, 288) else c(-108, 108)

    if (n_targets <= 16) {
      angles <- if (n_targets == 1) {
        mean(angle_limits)
      } else {
        seq(angle_limits[1], angle_limits[2], length.out = n_targets)
      }

      return(
        make_arc_coordinates(
          node_names = target_names,
          center_x = center_x,
          center_y = center_y,
          radius = 1.02,
          angles_deg = angles,
          radius_pattern = c(1.00, 1.08, 0.96, 1.04)
        )
      )
    }

    outer_names <- target_names[seq(1, n_targets, by = 2)]
    inner_names <- target_names[seq(2, n_targets, by = 2)]

    make_angles <- function(n, limits, inset = 0) {
      if (n == 0) return(numeric())
      if (n == 1) return(mean(limits))
      seq(limits[1] + inset, limits[2] - inset, length.out = n)
    }

    bind_rows(
      make_arc_coordinates(
        node_names = outer_names,
        center_x = center_x,
        center_y = center_y,
        radius = 1.23,
        angles_deg = make_angles(length(outer_names), angle_limits, inset = 0),
        radius_pattern = c(1.00, 1.04, 0.98)
      ),
      make_arc_coordinates(
        node_names = inner_names,
        center_x = center_x,
        center_y = center_y,
        radius = 0.88,
        angles_deg = make_angles(length(inner_names), angle_limits, inset = 8),
        radius_pattern = c(1.00, 1.05, 0.96)
      )
    )
  }

  target_positions <- tibble(name = character(), x = double(), y = double())

  if ("Trp63" %in% hub_positions$name) {
    trp_center <- hub_positions %>% filter(name == "Trp63")
    trp_targets <- target_source_map %>%
      filter(Source == "Trp63") %>%
      pull(Target) %>%
      unique() %>%
      order_targets()

    target_positions <- bind_rows(
      target_positions,
      place_module_targets(
        target_names = trp_targets,
        center_x = trp_center$x,
        center_y = trp_center$y,
        side = "left"
      )
    )
  }

  if ("Sox2" %in% hub_positions$name) {
    sox_center <- hub_positions %>% filter(name == "Sox2")
    sox_targets <- target_source_map %>%
      filter(Source == "Sox2") %>%
      pull(Target) %>%
      unique() %>%
      order_targets()

    target_positions <- bind_rows(
      target_positions,
      place_module_targets(
        target_names = sox_targets,
        center_x = sox_center$x,
        center_y = sox_center$y,
        side = "right"
      )
    )
  }

  coordinate_table <- bind_rows(hub_positions, target_positions) %>%
    distinct(name, .keep_all = TRUE)

  # 未排布节点沿外围椭圆兜底放置
  unplaced_nodes <- setdiff(vertex_table$name, coordinate_table$name)
  if (length(unplaced_nodes) > 0) {
    fallback_angles <- seq(0, 360, length.out = length(unplaced_nodes) + 1)[-(length(unplaced_nodes) + 1)]

    coordinate_table <- bind_rows(
      coordinate_table,
      tibble(
        name = unplaced_nodes,
        x = 2.30 * cos(fallback_angles * pi / 180),
        y = 1.55 * sin(fallback_angles * pi / 180)
      )
    )
  }

  coordinate_table <- tibble(name = vertex_table$name) %>%
    left_join(coordinate_table, by = "name")

  if (anyNA(coordinate_table$x) || anyNA(coordinate_table$y)) {
    stop("The manual network layout contains missing coordinates.")
  }

  layout_data <- create_layout(
    graph_object,
    layout = "manual",
    x = coordinate_table$x,
    y = coordinate_table$y
  )

  x_range <- range(layout_data$x, na.rm = TRUE)
  y_range <- range(layout_data$y, na.rm = TRUE)
  x_padding <- max(diff(x_range) * 0.17, 0.55)
  y_padding <- max(diff(y_range) * 0.24, 0.55)

  attr(layout_data, "plot_xlim") <- x_range + c(-x_padding, x_padding)
  attr(layout_data, "plot_ylim") <- y_range + c(-y_padding, y_padding)

  layout_data
}

# 13.7 绘制网络图:
# 节点填充色=靶基因表达轨迹; 节点形状=核心TF或靶基因; 箭头方向=TF->target;
# 边线型=direct/extended regulon; 边宽=regulon支持数(图例隐藏); 全部基因带标签;
# Trp63左hub/Sox2右hub; 边色区分两个模块; 同时输出高分辨率PNG
plot_transition_network <- function(edges, nodes, output_pdf, core_tfs) {
  if (nrow(edges) == 0) {
    stop("No TF-target edges remained after filtering.")
  }

  # 保留TF身份作为边属性用于模块着色
  edges <- edges %>%
    mutate(
      edge_source = Source,
      edge_source = factor(edge_source, levels = names(TF_EDGE_COLORS))
    )

  # 适度的节点尺寸范围传达动态强度而不掩盖拓扑
  target_scores <- nodes$dynamic_score[!nodes$is_core_tf]
  target_scores <- target_scores[is.finite(target_scores)]

  if (length(target_scores) == 0 || diff(range(target_scores)) == 0) {
    nodes <- nodes %>%
      mutate(plot_size = ifelse(is_core_tf, 25.0, 10.0))
  } else {
    nodes <- nodes %>%
      mutate(
        target_size = scales::rescale(
          log1p(dynamic_score),
          to = c(8.6, 12.4),
          from = range(log1p(target_scores), na.rm = TRUE)
        ),
        plot_size = ifelse(is_core_tf, 25.0, target_size)
      ) %>%
      select(-target_size)
  }

  # 靶基因标签显示规则:
  # Sox2靶基因仅Monotonic_up标注; Trp63独有靶基因全部标注
  sox2_target_genes <- edges %>%
    filter(Source == "Sox2", !Target %in% core_tfs) %>%
    pull(Target) %>%
    unique()

  nodes <- nodes %>%
    mutate(
      show_target_label = case_when(
        is_core_tf ~ TRUE,
        name %in% sox2_target_genes ~ as.character(transition_pattern) == "Monotonic_up",
        TRUE ~ TRUE
      )
    )

  graph_object <- graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = nodes
  ) %>%
    add_graph_statistics()

  layout_data <- build_reference_layout(
    graph_object = graph_object,
    core_tfs = core_tfs
  )

  plot_xlim <- attr(layout_data, "plot_xlim")
  plot_ylim <- attr(layout_data, "plot_ylim")
  y_span <- diff(plot_ylim)

  plot_object <- ggraph(layout_data) +
    # 极浅模块背景色增强视觉分组
    annotate(
      "rect",
      xmin = plot_xlim[1], xmax = 0,
      ymin = plot_ylim[1], ymax = plot_ylim[2],
      fill = scales::alpha(TF_EDGE_COLORS[["Trp63"]], 0.035),
      colour = NA
    ) +
    annotate(
      "rect",
      xmin = 0, xmax = plot_xlim[2],
      ymin = plot_ylim[1], ymax = plot_ylim[2],
      fill = scales::alpha(TF_EDGE_COLORS[["Sox2"]], 0.035),
      colour = NA
    ) +
    annotate(
      "text",
      x = -0.58,
      y = plot_ylim[2] - 0.05 * y_span,
      label = "Trp63 regulon",
      colour = TF_EDGE_COLORS[["Trp63"]],
      fontface = "bold",
      size = 3.8
    ) +
    annotate(
      "text",
      x = 0.58,
      y = plot_ylim[2] - 0.05 * y_span,
      label = "Sox2 regulon",
      colour = TF_EDGE_COLORS[["Sox2"]],
      fontface = "bold",
      size = 3.8
    ) +
    geom_edge_link(
      aes(
        edge_colour = edge_source,
        edge_width = regulon_evidence,
        edge_linetype = regulon_class,
        edge_alpha = regulon_class
      ),
      start_cap = ggraph::circle(3.6, "mm"),
      end_cap = ggraph::circle(6.6, "mm"),
      arrow = grid::arrow(
        length = grid::unit(3.5, "mm"),
        type = "closed"
      ),
      lineend = "round",
      show.legend = c(
        edge_linetype = TRUE,
        edge_width = FALSE,
        edge_colour = FALSE,
        edge_alpha = FALSE,
        fill = FALSE,
        shape = FALSE
      )
    ) +
    # 白色光晕保证节点在稠密边束上可读
    geom_node_point(
      aes(shape = node_type, size = plot_size + 1.4),
      fill = "white",
      colour = "white",
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_node_point(
      aes(
        fill = transition_pattern,
        colour = transition_pattern,
        shape = node_type,
        size = plot_size
      ),
      stroke = 0,
      alpha = 1,
      show.legend = c(fill = TRUE, shape = TRUE, size = FALSE)
    ) +
    # 基因名标签: Arial字体, 黑色, 无装饰
    geom_node_label(
      aes(filter = is_core_tf, label = plot_label),
      nudge_y = 0.23,
      size = 4.2,
      fontface = "plain",
      colour = "black",
      family = "Arial",
      fill = NA,
      label.size = 0,
      label.padding = grid::unit(0, "lines"),
      label.r = grid::unit(0, "lines"),
      show.legend = FALSE
    ) +
    geom_node_label(
      aes(filter = !is_core_tf & show_target_label, label = plot_label),
      repel = TRUE,
      size = 2.65,
      fontface = "plain",
      colour = "black",
      family = "Arial",
      fill = NA,
      label.size = 0,
      label.padding = grid::unit(0, "lines"),
      label.r = grid::unit(0, "lines"),
      box.padding = 0.22,
      point.padding = 0.10,
      min.segment.length = 0,
      segment.colour = NA,
      segment.size = 0,
      force = 0.75,
      max.iter = 50000,
      max.overlaps = Inf,
      seed = 123,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_edge_colour_manual(
      values = TF_EDGE_COLORS,
      na.value = "#A7A7A7",
      guide = "none"
    ) +
    scale_edge_alpha_manual(
      values = c(direct = 0.82, extended = 0.42),
      guide = "none"
    ) +
    scale_edge_linetype_manual(
      values = c(direct = "solid", extended = "22"),
      breaks = c("direct", "extended"),
      labels = c(direct = "Direct", extended = "Extended"),
      drop = FALSE,
      name = "Regulon class"
    ) +
    # 不设硬上限: evidence>2的边应保持可见
    scale_edge_width_continuous(
      range = c(0.45, 1.45),
      trans = "sqrt",
      guide = "none"
    ) +
    scale_fill_manual(
      values = TRAJECTORY_COLORS,
      breaks = PLOTTED_TRAJECTORIES,
      labels = TRAJECTORY_LABELS[PLOTTED_TRAJECTORIES],
      drop = FALSE,
      name = "Expression trajectory"
    ) +
    scale_colour_manual(
      values = TRAJECTORY_COLORS,
      guide = "none"
    ) +
    scale_size_identity() +
    scale_shape_manual(
      values = c(Core_TF = 21, Target = 22),
      breaks = c("Core_TF", "Target"),
      labels = c(Core_TF = "Core TF", Target = "Target gene"),
      drop = FALSE,
      name = "Node type"
    ) +
    guides(
      edge_linetype = guide_legend(
        order = 1,
        keywidth = grid::unit(12, "mm"),
        override.aes = list(
          edge_colour = "#777777",
          edge_width = 0.9,
          edge_alpha = 1
        )
      ),
      fill = guide_legend(
        order = 2,
        ncol = 1,
        override.aes = list(
          shape = 22,
          size = 5.0,
          colour = "transparent",
          stroke = 0,
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 3,
        override.aes = list(
          size = c(7.2, 5.0),
          fill = "white",
          colour = "transparent",
          stroke = 0,
          alpha = 1
        )
      )
    ) +
    labs(
      title = "Dynamic TF–target network across the LUAD–TP–SCC transition",
      subtitle = paste0("Node colour summarizes joint TP-vs-LUAD and TP-vs-SCC expression behaviour"),
      caption = paste0(
        "Arrows indicate TF → target regulation. Solid edges are direct SCENIC regulons; ",
        "dashed edges are extended regulons. Target size reflects the dynamic score."
      )
    ) +
    coord_fixed(
      ratio = 1,
      xlim = plot_xlim,
      ylim = plot_ylim,
      clip = "off",
      expand = FALSE
    ) +
    theme_void(base_size = 10, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = PLOT_BACKGROUND, colour = NA),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.background = element_rect(fill = scales::alpha("white", 0.94), colour = NA),
      legend.margin = margin(8, 9, 8, 9),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(5.5, "mm"),
      legend.title = element_text(
        size = 10,
        face = "bold",
        colour = LABEL_COLOR,
        margin = margin(b = 3)
      ),
      legend.text = element_text(size = 8.5, colour = "#3C4043"),
      legend.spacing.y = grid::unit(2.0, "mm"),
      plot.title = element_text(
        hjust = 0,
        size = 16,
        face = "bold",
        colour = "#17191C",
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        hjust = 0,
        size = 10,
        colour = "#555B63",
        margin = margin(b = 10)
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 8.3,
        lineheight = 1.15,
        colour = "#626870",
        margin = margin(t = 10)
      ),
      plot.margin = margin(14, 16, 12, 16)
    )

  ggsave(
    filename = output_pdf,
    plot = plot_object,
    device = grDevices::cairo_pdf,
    width = 15.5,
    height = 10,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  # 高分辨率栅格副本便于幻灯片展示与快速查看
  output_png <- sub("\\.pdf$", ".png", output_pdf, ignore.case = TRUE)
  if (!identical(output_png, output_pdf)) {
    ggsave(
      filename = output_png,
      plot = plot_object,
      width = 15.5,
      height = 10,
      units = "in",
      dpi = 450,
      bg = "white",
      limitsize = FALSE
    )
  }

  invisible(plot_object)
}

# 13.8 读取并准备数据
check_file(REGULONS_FILE)
check_file(DEG_LUAD_FILE)
check_file(DEG_SCC_FILE)

regulons <- readRDS(REGULONS_FILE)

deg_luad_raw <- read.csv(
  DEG_LUAD_FILE,
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

deg_scc_raw <- read.csv(
  DEG_SCC_FILE,
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

deg_luad_raw$gene <- rownames(deg_luad_raw)
deg_scc_raw$gene <- rownames(deg_scc_raw)

deg_luad <- prepare_deg(
  deg_luad_raw,
  contrast_name = "TP_vs_LUAD",
  padj_cutoff = PADJ_CUTOFF,
  logfc_cutoff = LOGFC_CUTOFF
)

deg_scc <- prepare_deg(
  deg_scc_raw,
  contrast_name = "TP_vs_SCC",
  padj_cutoff = PADJ_CUTOFF,
  logfc_cutoff = LOGFC_CUTOFF
)

transition_table <- build_transition_table(
  deg_luad_clean = deg_luad,
  deg_scc_clean = deg_scc
)

regulon_long <- regulons_to_long(regulons)

# 13.9 构建整合的动态网络
edges_transition <- build_transition_network(
  regulon_long = regulon_long,
  transition_table = transition_table,
  core_tfs = CORE_TFS,
  top_targets_per_tf = TOP_TARGETS_PER_TF,
  include_extended_regulons = INCLUDE_EXTENDED_REGULONS
)

if (nrow(edges_transition) == 0) {
  stop(
    "No TF-target edges remained after filtering. ",
    "Check regulon names, gene symbols, and DEG thresholds."
  )
}

nodes_transition <- build_node_table(
  edges = edges_transition,
  transition_table = transition_table,
  core_tfs = CORE_TFS
)

# 13.10 生成网络图(dynamic_core_TF_network_LUAD_TP_SCC.pdf及同名PNG)
plot_transition_network(
  edges = edges_transition,
  nodes = nodes_transition,
  output_pdf = OUTPUT_PDF,
  core_tfs = CORE_TFS
)

cat("Figure generated:\n", OUTPUT_PDF, "\n", sep = "")
