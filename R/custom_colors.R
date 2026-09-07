# ============================================================================
# 脚本名称: custom_colors.R
# 描述: 项目统一使用的颜色配置（荧光分群、tdRatio、聚类、伪时间等配色方案）
# 依赖包: RColorBrewer, scico, viridis, ggsci, scales
# ============================================================================

require(RColorBrewer)
require(scico)
require(viridis)
require(ggsci)  # https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

# ---- 荧光分群颜色: GFP / tdTomato / none / immune ----
cols_fluorescence <- setNames(
  c("#5DA6C6", "#97DCAE", "#eb9898", "grey90"),
  nm = c("none", "GFP", "tdTomato", "immune")
)

# ---- tdTomato_ratio (TRP63 表达比例) 渐变色 ----
cols_td_ratio <- c("#97DCAE", "#F4D8A8", "#EDC962", "#e4ae69", "#eb9898")

# ---- CytoTRACE / CellRank 渐变色 ----
cols_cytotrace <- c("#000004FF", "#3B0F70FF", "#8C2981FF",
                    "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")
col_scale_cellrank <- viridis::viridis(256)

# ---- 常用连续型配色 ----
col_scale_reds    <- RColorBrewer::brewer.pal(9, "Reds")
col_scale_greens  <- RColorBrewer::brewer.pal(9, "Greens")
col_scale_blues   <- RColorBrewer::brewer.pal(9, "Blues")
col_scale_purples <- RColorBrewer::brewer.pal(9, "Purples")
col_scale_bupu    <- RColorBrewer::brewer.pal(9, "BuPu")
col_scale_purd    <- RColorBrewer::brewer.pal(9, "PuRd")
col_scale_mako    <- viridis::mako(22, direction = -1)

# ---- 伪时间渐变色 (深蓝 -> 浅蓝 -> 黄 -> 橙 -> 红) ----
my_pseudotime_colors <- colorRampPalette(
  c("darkblue", "lightblue", "yellow", "orange", "red"))(100)

# ---- 分类变量配色 ----
cols_cat_custom <- c("#5DA6C6", "#3873B9", "#215189",
                     "#E391BA", "#8D71A6",
                     "#F4D8A8", "#EDC962", "#E7B922", "#DB8712", "#CA6623", "#AB4A3D",
                     "#B1E4C2", "#78D6AE", "#49C1AD")

cols_cnv <- c("#5DA6C6", "#3873B9",
              "#E391BA", "#8D71A6",
              "#F4D8A8", "#EDC962",
              "#B1E4C2", "#78D6AE")

# ---- 聚类配色 ----
cols_cluster <- c(col_scale_reds[2:8],
                  col_scale_greens[3:7],
                  col_scale_blues[2:7],
                  col_scale_mako)

cols_seurat <- setNames(scales::hue_pal()(4), nm = c(0:3))  # Seurat 默认配色

# ---- 细胞类型注释配色 ----
cols_cluster_annotated <- setNames(
  c(col_scale_reds[2:7], col_scale_greens[3:7], col_scale_blues[2:8]),
  nm = c("AT1", "AT2", "Club", "AT2-like Club", "Goblet-like Club", "Cliated",
         "KAC", "Malignant", "Proliferating",
         "Endothelium", "Fibroblast",
         "B cell", "T cell", "NKT", "NK", "Granulocytes", "Monocyte", "cDCs")
)
