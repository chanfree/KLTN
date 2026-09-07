# ============================================================================
# 脚本名称: custom_functions.R
# 描述: 项目通用绘图与辅助函数（QC 可视化、细胞比例统计、通路热图、富集气泡图）
# 依赖包: reshape2, ggalluvial, dplyr, org.Mm.eg.db, pheatmap, Seurat, ggplot2
# 注意: 部分函数依赖 custom_colors.R 中的颜色对象，请先 source 该文件
# ============================================================================

require(reshape2)
require(ggalluvial)
require(dplyr)
require(org.Mm.eg.db)
require(pheatmap)
require(Seurat)

# ----------------------------------------------------------------------------
# plotQC: 细胞 QC 指标可视化
#   参数:
#     seurat.object - Seurat 对象
#     fill.by       - 用于分面/填充的 metadata 列名 (默认 "batch")
#     label.by      - 用于分组的 metadata 列名 (默认 "label")
#     color.group   - 命名颜色向量 (默认 cols_fluorescence)
#   返回: 组合图 (基因数直方图 + UMI 数直方图 / 小提琴图)
# ----------------------------------------------------------------------------
plotQC <- function(seurat.object,
                   fill.by = "batch",
                   label.by = "label",
                   color.group = cols_fluorescence) {
  p_dat <- seurat.object@meta.data

  # 每细胞检测基因数分布
  p1 <- ggplot(p_dat, aes_string(x = "nFeature_RNA")) +
    geom_histogram(aes(y = ..density..), colour = NA, fill = "black", bins = 75) +
    geom_density(data = p_dat,
                 aes_string(x = "nFeature_RNA", y = "..density..", fill = fill.by),
                 alpha = .5) +
    geom_vline(aes(xintercept = mean(nFeature_RNA)),
               color = "orange", linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = color.group) +
    ggtitle("Unique genes per cells") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  # 每细胞 UMI 总数分布 (log10 坐标)
  p2 <- ggplot(p_dat, aes_string(x = "nCount_RNA")) +
    geom_histogram(aes(y = ..density..), colour = NA, fill = "black", bins = 75) +
    geom_density(data = p_dat,
                 aes_string(x = "nCount_RNA", y = "..density..", fill = fill.by),
                 alpha = .5) +
    geom_vline(aes(xintercept = mean(nFeature_RNA)),
               color = "orange", linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = color.group) +
    scale_x_log10() +
    ggtitle("Total counts per cells") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  # 按 label.by 分组、fill.by 分面的 QC 小提琴图
  p3 <- VlnPlot(seurat.object,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0,
                group.by = label.by,
                split.by = fill.by,
                cols = color.group) & theme_custom

  p <- (p1 + p2) / p3
  return(p)
}

# ----------------------------------------------------------------------------
# plotProportion: 两组 metadata 之间的细胞比例统计图
#   功能: 计算每个 meta.1 分组中 meta.2 各类别所占比例,
#         绘制堆叠条形图与 alluvial 冲积图
#   参数:
#     seurat.object - Seurat 对象
#     meta.1        - 分组列名 (默认 "seurat_clusters")
#     meta.2        - 类别列名 (默认 "annotation")
#   返回: 图列表 list(堆叠条形图, 冲积图)
# ----------------------------------------------------------------------------
plotProportion <- function(seurat.object,
                           meta.1 = "seurat_clusters",
                           meta.2 = "annotation") {
  metadata <- seurat.object[[]]

  # 统计 meta.1 x meta.2 交叉频数与比例
  dat1 <- metadata |>
    group_by(!!sym(meta.1), !!sym(meta.2)) |>
    summarise(count = n()) |>
    ungroup()
  dat2 <- dat1 |>
    group_by(!!sym(meta.1)) |>
    summarise(total = sum(count))
  dat_proportions <- dat1 |>
    left_join(dat2, by = meta.1) |>
    mutate(proportion = count / total)

  # 长宽数据转换
  data_long <- dcast(dat_proportions, reformulate(c(meta.1), response = meta.2),
                     value.var = "proportion", fill = 0)
  data_melted <- melt(data_long, id.vars = meta.2)
  names(data_melted) <- c(meta.2, meta.1, "proportion")

  # 堆叠条形图
  p1 <- ggplot(data = data_melted,
               aes(!!sym(meta.1), proportion, fill = !!sym(meta.2))) +
    geom_bar(stat = "identity", position = "fill",
             color = "black", width = 0.6, linewidth = 0.4)

  # 冲积图 (alluvial)
  p2 <- ggplot(data = data_melted,
               aes(!!sym(meta.1), proportion, fill = !!sym(meta.2),
                   stratum = !!sym(meta.2), alluvium = !!sym(meta.2))) +
    geom_stratum(color = "black", width = 0.3, linewidth = 0.5) +
    geom_flow(alpha = 0.5) +
    theme_minimal()

  return(list(p1, p2))
}

# ----------------------------------------------------------------------------
# plot_pathway_DEG: 指定通路中差异基因的表达热图
#   功能: 从 KEGG 等注释数据库提取目标通路的基因, 与差异基因取交集,
#         绘制各聚群平均表达的列标准化热图
#   参数:
#     se_obj           - Seurat 对象
#     interetsed_pathway - 目标通路名称向量
#     database         - 注释数据库名 (默认 "KEGG")
#     DEG              - 差异基因数据框, 需含 "gene" 列
#   返回: pheatmap 对象
# ----------------------------------------------------------------------------
plot_pathway_DEG <- function(se_obj, interetsed_pathway, database = "KEGG", DEG) {
  # 获取物种注释并筛选目标通路基因
  hsko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = database)
  gene_dat <- hsko@annot[hsko@annot$Annot %in% interetsed_pathway, ]
  gene_dat <- gene_dat[gene_dat$GeneID %in% DEG$gene, ]
  gene_dat <- arrange(gene_dat, Annot)

  # 各聚群平均表达矩阵
  avg.exp <- AverageExpression(se_obj)
  dat <- avg.exp[["RNA"]][gene_dat$GeneID, ]
  unqie_gene <- make.unique(gene_dat$GeneID, sep = " #")  # 重名基因加 # 后缀
  rownames(dat) <- unqie_gene
  dat <- t(dat)

  # 通路注释条
  col_annotation <- data.frame(Pathway = factor(gene_dat$Annot))
  rownames(col_annotation) <- colnames(dat)

  p <- pheatmap(
    dat,
    scale = "column",
    color = col_scale_bupu,           # 使用 custom_colors.R 中的配色
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = col_annotation,
    annotation_legend = TRUE,
    border_color = "grey60",
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = 10
  )
  return(p)
}

# ----------------------------------------------------------------------------
# erich2plot: 富集结果气泡图 (取 qvalue 最小的前 20 条)
#   参数:
#     data4plot - 富集结果数据框, 需含 qvalue / RichFactor / Description / Count 列
#   返回: ggplot 对象
# ----------------------------------------------------------------------------
erich2plot <- function(data4plot) {
  require(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE)[1:20], ]

  p <- ggplot(data4plot, aes(RichFactor, Description))
  p <- p + geom_point()

  pbubble <- p + geom_point(aes(size = Count, color = -1 * log10(qvalue)))

  pr <- pbubble +
    scale_colour_gradient(low = "#90EE90", high = "red") +
    labs(color = expression(-log[10](qvalue)), size = "observed.gene.count",
         x = "Richfactor", y = "term.description", title = "Enrichment Process") +
    theme_bw()
  pr
}
