# ============================================================================
# 脚本名称: 05_survival_analysis.R
# 描述: 临床价值验证——NRG1及本研究全部signature(含ast_Traj_Up)在KRYSTAL-1试验与
#       TCGA-LUAD/LUSC队列中的Kaplan-Meier生存分析, 以及EGFR/KRAS通路活性与
#       signature评分的相关性分析(KRAS单细胞 + KRYSTAL-1 bulk)
# 输入: output/object/survival/KRYSTAL-1_trial/vst_counts.rds
#       output/object/survival/KRYSTAL-1_trial/meta_dat.csv
#       (源数据 mmc2 (1).xlsx / mmc3 (1).xlsx, 由 RUN_KRYSTAL_PREPROCESS 开关生成)
#       output/object/survival/TCGA-LUAD_vst_counts_with_symbols.rds
#       output/object/survival/TCGA-LUAD_clinical_data.rds
#       output/object/survival/TCGA-LUSC_vst_counts_with_symbols.rds
#       output/object/survival/TCGA-LUSC_clinical_data.rds
#       (TCGA数据由03脚本 RUN_TCGA_DOWNLOAD 下载并预处理)
#       data/gene_orthology/human_all_signature_without_pattern_list.rds
#       output/object/survival/human_all_signature_without_pattern_list.rds
#           (来源说明: 13 个鼠源 signature 转人源同源基因后的命名字符向量列表,
#            ast_Traj_Assoc/ast_Traj_Up/ast_td_vs_gfp/ast_DNB_DEG_Pre/Post/
#            ast_DNB_Var_Pre/Post/ast_DNB_Tip_Marker/ast_DNB_PEAK 及 4 个 integrated_deg_* 签名)
#       data/clinical_validation/KRAS_induced& suppressed.csv
#       data/clinical_validation/KOBAYASHI_EGFR_SIGNALING_24HR_UP.v2026.1.Hs.gmt
#       data/clinical_validation/KOBAYASHI_EGFR_SIGNALING_24HR_DN.v2026.1.Hs.gmt
#       data/clinical_validation/kras_se_obj_with_scores.rds
#       (KRAS单细胞对象原位于 clinical_analysis/results/KRAS_clinical_analysis/,
#        由外部脚本生成, 需复制到 data/clinical_validation/)
#       data/clinical_validation/egfr_epi_subset.rds
#       data/clinical_validation/EGFR_Seu_metadata.csv
#       (EGFR单细胞数据: 处理版来自 GitHub scell_lung_adenocarcinoma 的 Data_input,
#        原始测序数据 PRJNA591860, 转化脚本见 data/clinical_validation/EGFR_seu_convertion.R;
#        需复制到 data/clinical_validation/, 用于FigS5 D-F)
#       data/gene_orthology/gene_signatures_list_mouse.rds
#       data/clinical_validation/merged_tpm_organoid_processed.csv
#       (小鼠类器官TPM表达矩阵, 供FigS5 G的signature分组箱线图分析;
#        TPM合并处理脚本见 scr/process_tpm_organoid.R)
# 输出: output/figures/survival/ (KRYSTAL-1 NRG1 KM曲线)
#       output/figures/survival/human_signature/ (signature重叠UpSet图, 各队列整合KM图)
#       output/figures/survival/score_analysis/ (相关性热图/散点图/KM曲线,
#       signature表达箱线图, 见 FigS5)
#       output/figures/survival/score_analysis/mouse_organoid/ (FigS5 G类器官热图/箱线图/小提琴图)
#       output/tables/score_analysis/ (基因集信息与相关性/生存统计表)
#       output/tables/score_analysis/mouse_organoid/ (FigS5 G类器官得分/匹配率/ANOVA/两两Wilcoxon表)
#       output/object/survival/all_signature_list_human.rds
#       output/object/survival/score_analysis/KRAS_sc_analysis_results.rds
#       output/object/survival/score_analysis/EGFR_sc_analysis_results.rds
#       output/object/survival/score_analysis/KRYSTAL1_analysis_results.rds
#       output/object/survival/score_analysis/mouse_organoid_results.rds
# 上游脚本: 03_trajectory_subset_analysis.R / 04_DNB_analysis.R (signature相关)
# ============================================================================

# ---- 0. 环境与路径配置 ----

set.seed(1)

DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "output")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "object")
fig_res <- 300

# 生存分析各子目录(沿用原输出目录结构)
DIR_FIG_survival <- file.path(DIR_FIG, "survival")
DIR_FIG_survival_sig <- file.path(DIR_FIG, "survival", "human_signature")
DIR_OBJ_survival <- file.path(DIR_OBJ, "survival")
# EGFR/KRAS通路评分分析子目录(原位于 clinical_analysis/results/score_analysis/,
# 统一移至 output/ 下)
DIR_FIG_score <- file.path(DIR_FIG, "survival", "score_analysis")
DIR_TABLE_score <- file.path(DIR_RES, "tables", "score_analysis")
DIR_OBJ_score <- file.path(DIR_OBJ_survival, "score_analysis")
# 临床验证输入数据目录
DIR_DATA_CLINVAL <- file.path(DIR_DATA, "clinical_validation")

for (d in c(DIR_FIG_survival, DIR_FIG_survival_sig, DIR_OBJ_survival,
            DIR_FIG_score, DIR_TABLE_score, DIR_OBJ_score)) {
  dir.create(path = d, recursive = TRUE, showWarnings = FALSE)
}

# 逻辑开关(重跑控制) ---------------------------------------------------------
# KRYSTAL-1原始xlsx预处理(mmc2/mmc3 → vst_counts.rds + meta_dat.csv):
# 已运行过, 结果已保存, 重跑时设为TRUE
# (需DESeq2包, 且 mmc2 (1).xlsx / mmc3 (1).xlsx 位于
#  output/object/survival/KRYSTAL-1_trial/ 目录)
RUN_KRYSTAL_PREPROCESS <- FALSE
# signature长表CSV导出(一行一个基因, 供外部工具使用):
# 已运行过, 结果已保存, 重跑时设为TRUE
RUN_EXPORT_SIGNATURE_LONG_CSV <- FALSE

# 加载必要的包
# (原脚本另加载 TCGAbiolinks/SummarizedExperiment/viridis/readr,
#  本脚本未用到, 已移除; 原个人R库路径 .libPaths("/refdir/Rlib") 设置已删除)
# (pheatmap/RColorBrewer 为第5章FigS5 G小鼠类器官热图所需,
#  patchwork 为第4章FigS5 A-F箱线图排版所需)
library(readxl)
library(survival)
library(survminer)
library(orthogene)
library(UpSetR)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(future)
library(Matrix)
library(reshape2)
library(corrplot)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

plan(multisession, workers = 8)
options(future.globals.maxSize = 60000 * 1024^2, future.seed = TRUE)

# ---- 1. KRYSTAL-1试验: NRG1单基因KM生存曲线 ----

# KRYSTAL-1试验(KRAS G12C抑制剂Adagrasib)批量RNA-seq队列
# 数据来源: Tong, X. et al. Adeno-to-squamous transition drives resistance to KRAS
# inhibition in LKB1 mutant lung cancer. Cancer Cell 42, 413-428.e7 (2024).
tcga_project <- "KRYSTAL-1_trial"

# 1.1 预处理开关: 从论文附表xlsx生成VST标准化表达矩阵与临床信息表
if (RUN_KRYSTAL_PREPROCESS) {
  meta_dat <- read_excel(
    file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "mmc2 (1).xlsx"),
    sheet = "C-KRT6A and SCC signatures",
    skip = 2
  )
  raw_count <- read_excel(
    file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "mmc3 (1).xlsx"),
    skip = 1
  )
  raw_count <- tibble::column_to_rownames(raw_count, var = "Gene")

  # 借助KRT6A的log2CPM样本顺序, 将样本名对齐至临床表
  cpm <- 1e6 * sweep(raw_count, 2, colSums(raw_count), FUN = "/")
  log2cpm <- log2(cpm + 1)
  temp_dat <- t(log2cpm["KRT6A", ])
  meta_dat$sample <- rownames(temp_dat)

  raw_count <- as.matrix(raw_count)
  vst_counts <- DESeq2::varianceStabilizingTransformation(round(raw_count), blind = TRUE)
  # 删失编码转换: 原表cnsr 1=删失 → 转为0=删失(1=事件), 与Surv(event=)约定一致
  meta_dat$pfsind_cnsr <- ifelse(meta_dat$pfsind_cnsr == 1, 0, 1)
  meta_dat$os_cnsr <- ifelse(meta_dat$os_cnsr == 1, 0, 1)

  saveRDS(vst_counts, file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "vst_counts.rds"))
  write.csv(meta_dat, file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "meta_dat.csv"))
}

vst_counts <- readRDS(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "vst_counts.rds"))
meta_dat <- read.csv(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "meta_dat.csv"),
                     header = TRUE, row.names = 1)

# 1.2 构建生存分析数据: signature评分 + 按中位数分高低风险组
generate_survival_data <- function(mouse_gene_signature, vst_counts, clinical_data) {
  if (length(mouse_gene_signature) > 1) {
    # 多基因signature: 先将小鼠基因转换为人类同源基因
    require(orthogene)
    gene_signature <- convert_orthologs(
      gene_df = mouse_gene_signature,
      input_species = "mouse",
      output_species = "human",
      method = "homologene"
    )
    gene_signature <- unique(rownames(gene_signature))
    signature_expr <- vst_counts[rownames(vst_counts) %in% gene_signature, ] |> as.data.frame()
    # 检查未能匹配到的signature基因并给出警告
    missing_genes <- setdiff(gene_signature, rownames(signature_expr))
    if (length(missing_genes) > 0) {
      warning(paste("The following genes were not found in the dataset:",
                    paste(missing_genes, collapse = ", ")))
    }

    # 计算signature评分(每个样本取signature基因平均表达), 转置使基因为列、样本为行
    signature_expr_t <- as.data.frame(t(signature_expr))
    signature_expr_t$signature_score <- rowMeans(signature_expr_t, na.rm = TRUE)
  } else {
    # 单基因: 直接以该基因表达为评分
    gene_signature <- mouse_gene_signature
    signature_expr <- vst_counts[rownames(vst_counts) %in% gene_signature, ] |> as.data.frame()
    colnames(signature_expr) <- gene_signature
    signature_expr_t <- signature_expr
    signature_expr_t$signature_score <- signature_expr_t[, 1]
  }

  # 将signature评分合并至临床数据
  survival_data <- merge(clinical_data, signature_expr_t[, "signature_score", drop = FALSE],
                         by.x = 0, by.y = 0)
  colnames(survival_data)[1] <- "patient_barcode"

  # 按signature评分中位数将患者分为高/低风险组
  median_score <- median(survival_data$signature_score)
  survival_data$risk_group <- ifelse(survival_data$signature_score >= median_score,
                                     "High_Risk", "Low_Risk")
  survival_data$risk_group <- factor(survival_data$risk_group, levels = c("Low_Risk", "High_Risk"))
  return(survival_data)
}

# 1.3 NRG1单基因: OS(天)的KM曲线
gene_list <- "NRG1"
for (gene in gene_list) {
  survival_data <- generate_survival_data(gene, vst_counts = vst_counts, clinical_data = meta_dat)

  # 按风险组拟合生存模型并绘制KM曲线(含log-rank p值)
  fit <- survfit(Surv(time = survival_data$os, event = survival_data$os_cnsr) ~ risk_group,
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
    title = paste("Kaplan-Meier Curve for", tcga_project, "by Gene", gene),
    xlab = "Time (Days)", ylab = "OS",
    ggtheme = theme_classic()
  )
  pdf(file.path(DIR_FIG_survival, sprintf("KRYSTAL-1_trial_%s_mutant_km_os.pdf", gene)),
      width = 10, height = 6)
  print(km_plot)
  dev.off()
}

# 1.4 NRG1单基因: OS(月)的KM曲线, 图例标注各组人数
for (gene in gene_list) {
  survival_data <- generate_survival_data(
    gene,
    vst_counts = vst_counts,
    clinical_data = meta_dat
  )

  # OS: 天 → 月
  survival_data$os_month <- survival_data$os / 30.44

  fit <- survfit(
    Surv(time = os_month, event = os_cnsr) ~ risk_group,
    data = survival_data
  )
  # 统计分组人数并拼入图例: 组名(人数)
  group_cnt <- table(survival_data$risk_group)
  low_n <- group_cnt[["Low_Risk"]]
  high_n <- group_cnt[["High_Risk"]]
  legend_label_new <- c(paste0("Low Risk (n=", low_n, ")"),
                        paste0("High Risk (n=", high_n, ")"))
  km_plot <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    legend.labs = legend_label_new,
    legend.title = "Risk Group",
    palette = c("#0073C2FF", "#E7B800FF"),
    title = paste("Kaplan-Meier Curve for", tcga_project, "by Gene", gene),
    xlab = "Time (Months)",
    ylab = "OS",
    ggtheme = theme_classic()
  )
  pdf(file.path(DIR_FIG_survival,
                sprintf("KRYSTAL-1_trial_%s_mutant_km_os_month.pdf", gene)),
      width = 6, height = 5)
  print(km_plot)
  dev.off()
}

# ---- 2. 本研究全部signature在TCGA-LUAD/LUSC/KRYSTAL-1队列的KM验证 ----

# 2.1 读取人类同源signature列表并绘制重叠UpSet图
# (signature由小鼠单细胞分析得到, 人类同源转换在另一台计算机经babelgene完成,
#  转换结果存放于 data/gene_orthology/)
gene_signatures_list <- readRDS(file.path(DIR_DATA, "gene_orthology",
                                          "human_all_signature_without_pattern_list.rds"))
# 备份一份到输出目录
saveRDS(gene_signatures_list, file.path(DIR_OBJ_survival, "all_signature_list_human.rds"))
# 注: 原脚本此处再次readRDS output/object/survival/下的同名副本
# (与data/gene_orthology版本同名, 内容一致), 属冗余读取, 已注释;
# 如需完全复现原行为可取消下一行注释
# gene_signatures_list <- readRDS(file.path(DIR_OBJ_survival,
#                                           "human_all_signature_without_pattern_list.rds"))

# 可选: 导出signature长表CSV(一行一个基因)
if (RUN_EXPORT_SIGNATURE_LONG_CSV) {
  df_long <- do.call(rbind, lapply(names(gene_signatures_list), function(sig) {
    data.frame(
      signature = sig,
      gene = gene_signatures_list[[sig]]
    )
  }))
  write.csv(df_long,
            file.path(DIR_DATA, "gene_orthology", "human_all_signature_without_pattern_list.csv"),
            row.names = FALSE)
}

# 各signature间基因重叠UpSet图
p <- upset(fromList(gene_signatures_list),
           nsets = 23,
           order.by = "freq",
           main.bar.color = "black",
           sets.bar.color = "darkblue",
           text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1),
           mb.ratio = c(0.6, 0.4)
)
pdf(file.path(DIR_FIG_survival_sig, "overlap_all_sigantures_human.pdf"), width = 20, height = 10)
print(p)
dev.off()

# 2.2 TCGA-LUAD: 逐个signature计算评分、分组并绘制KM曲线, 汇总为一份PDF
tcga_project <- "TCGA-LUAD"
vst_counts_with_symbols <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_vst_counts_with_symbols.rds", tcga_project)))
clinical_data <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_clinical_data.rds", tcga_project)))

p_list <- lapply(names(gene_signatures_list), function(signature_name) {
  gene_signature <- unique(gene_signatures_list[[signature_name]])
  # 提取signature基因表达并计算每样本平均评分
  signature_expr <- vst_counts_with_symbols[rownames(vst_counts_with_symbols) %in% gene_signature, ]
  signature_expr_t <- as.data.frame(t(signature_expr))
  signature_expr_t$signature_score <- rowMeans(signature_expr_t, na.rm = TRUE)

  # 合并临床数据并按中位数分高低风险组
  survival_data <- merge(clinical_data, signature_expr_t[, "signature_score", drop = FALSE],
                         by.x = 0, by.y = 0)
  colnames(survival_data)[1] <- "patient_barcode"
  median_score <- median(survival_data$signature_score)
  survival_data$risk_group <- ifelse(survival_data$signature_score >= median_score,
                                     "High_Risk", "Low_Risk")
  survival_data$risk_group <- factor(survival_data$risk_group, levels = c("Low_Risk", "High_Risk"))

  # 图例标注各组人数
  group_cnt <- table(survival_data$risk_group)
  low_n <- group_cnt[["Low_Risk"]]
  high_n <- group_cnt[["High_Risk"]]
  legend_label_new <- c(paste0("Low Risk (n=", low_n, ")"),
                        paste0("High Risk (n=", high_n, ")"))

  fit <- survfit(Surv(time = survival_data$time_to_event,
                      event = survival_data$vital_status_numeric) ~ risk_group,
                 data = survival_data)
  km_plot <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    risk.table.col = "strata",
    legend.labs = legend_label_new,
    legend.title = "Risk Group",
    palette = c("#0073C2FF", "#E7B800FF"),
    title = paste("Kaplan-Meier Curve for", tcga_project,
                  "by Gene Signature", signature_name),
    xlab = "Time (Days)",
    ggtheme = theme_classic()
  )
  return(km_plot)
})
names(p_list) <- names(gene_signatures_list)
pdf(file.path(DIR_FIG_survival_sig, "KM_plot_integrated_all_signature_vs_TCGA-LUAD.pdf"),
    width = 8, height = 6)
print(p_list)
dev.off()

# 2.3 TCGA-LUSC: 同上流程
tcga_project <- "TCGA-LUSC"
vst_counts_with_symbols <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_vst_counts_with_symbols.rds", tcga_project)))
clinical_data <- readRDS(file.path(
  DIR_OBJ_survival, sprintf("%s_clinical_data.rds", tcga_project)))

p_list <- lapply(names(gene_signatures_list), function(signature_name) {
  gene_signature <- unique(gene_signatures_list[[signature_name]])
  signature_expr <- vst_counts_with_symbols[rownames(vst_counts_with_symbols) %in% gene_signature, ]

  # 检查未能匹配到的signature基因并给出警告
  missing_genes <- setdiff(gene_signature, rownames(signature_expr))
  if (length(missing_genes) > 0) {
    warning(paste("The following genes were not found in the dataset:",
                  paste(missing_genes, collapse = ", ")))
  }

  # 计算signature评分、合并临床数据并按中位数分组
  signature_expr_t <- as.data.frame(t(signature_expr))
  signature_expr_t$signature_score <- rowMeans(signature_expr_t, na.rm = TRUE)
  survival_data <- merge(clinical_data, signature_expr_t[, "signature_score", drop = FALSE],
                         by.x = 0, by.y = 0)
  colnames(survival_data)[1] <- "patient_barcode"
  median_score <- median(survival_data$signature_score)
  survival_data$risk_group <- ifelse(survival_data$signature_score >= median_score,
                                     "High_Risk", "Low_Risk")
  survival_data$risk_group <- factor(survival_data$risk_group, levels = c("Low_Risk", "High_Risk"))

  # 图例标注各组人数
  group_cnt <- table(survival_data$risk_group)
  low_n <- group_cnt[["Low_Risk"]]
  high_n <- group_cnt[["High_Risk"]]
  legend_label_new <- c(paste0("Low Risk (n=", low_n, ")"),
                        paste0("High Risk (n=", high_n, ")"))

  fit <- survfit(Surv(time = survival_data$time_to_event,
                      event = survival_data$vital_status_numeric) ~ risk_group,
                 data = survival_data)
  km_plot <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    risk.table.col = "strata",
    legend.labs = legend_label_new,
    legend.title = "Risk Group",
    palette = c("#0073C2FF", "#E7B800FF"),
    title = paste("Kaplan-Meier Curve for", tcga_project,
                  "by Gene Signature", signature_name),
    xlab = "Time (Days)",
    ggtheme = theme_classic()
  )
  return(km_plot)
})
names(p_list) <- names(gene_signatures_list)
pdf(file.path(DIR_FIG_survival_sig, "KM_plot_integrated_all_signature_vs_TCGA-LUSC.pdf"),
    width = 8, height = 6)
print(p_list)
dev.off()

# 2.4 KRYSTAL-1试验: 同上流程, 时间轴为OS(月)
tcga_project <- "KRYSTAL-1_trial"
vst_counts_with_symbols <- readRDS(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "vst_counts.rds"))
clinical_data <- read.csv(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "meta_dat.csv"),
                          header = TRUE, row.names = 1)

p_list <- lapply(names(gene_signatures_list), function(signature_name) {
  gene_signature <- unique(gene_signatures_list[[signature_name]])
  signature_expr <- vst_counts_with_symbols[rownames(vst_counts_with_symbols) %in% gene_signature, ]

  # 检查未能匹配到的signature基因并给出警告
  missing_genes <- setdiff(gene_signature, rownames(signature_expr))
  if (length(missing_genes) > 0) {
    warning(paste("The following genes were not found in the dataset:",
                  paste(missing_genes, collapse = ", ")))
  }

  # 计算signature评分、合并临床数据并按中位数分组
  signature_expr_t <- as.data.frame(t(signature_expr))
  signature_expr_t$signature_score <- rowMeans(signature_expr_t, na.rm = TRUE)
  survival_data <- merge(clinical_data, signature_expr_t[, "signature_score", drop = FALSE],
                         by.x = 0, by.y = 0)
  colnames(survival_data)[1] <- "patient_barcode"
  median_score <- median(survival_data$signature_score)
  survival_data$risk_group <- ifelse(survival_data$signature_score >= median_score,
                                     "High_Risk", "Low_Risk")
  survival_data$risk_group <- factor(survival_data$risk_group, levels = c("Low_Risk", "High_Risk"))

  # 图例标注各组人数
  group_cnt <- table(survival_data$risk_group)
  low_n <- group_cnt[["Low_Risk"]]
  high_n <- group_cnt[["High_Risk"]]
  legend_label_new <- c(paste0("Low Risk (n=", low_n, ")"),
                        paste0("High Risk (n=", high_n, ")"))

  # OS: 天 → 月
  survival_data$os_month <- survival_data$os / 30.44
  fit <- survfit(
    Surv(time = os_month, event = os_cnsr) ~ risk_group,
    data = survival_data
  )
  km_plot <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    risk.table.col = "strata",
    legend.labs = legend_label_new,
    legend.title = "Risk Group",
    palette = c("#0073C2FF", "#E7B800FF"),
    title = paste("Kaplan-Meier Curve for", tcga_project,
                  "by Gene Signature", signature_name),
    xlab = "Time (Months)",
    ggtheme = theme_classic()
  )
  return(km_plot)
})
names(p_list) <- names(gene_signatures_list)
pdf(file.path(DIR_FIG_survival_sig, "KM_plot_integrated_all_signature_vs_KRYSTAL-1.pdf"),
    width = 8, height = 6)
print(p_list)
dev.off()

# 单独保存ast_Traj_Up signature的KM曲线(原脚本依赖循环最后一次的km_plot,
# 即signature列表最后一个元素, ast_Traj_Up需位于列表末尾方可得到其曲线)
pdf(file.path(DIR_FIG_survival_sig, "KM_plot_ast_Traj_Up_signature_vs_KRYSTAL-1_os_month.pdf"),
    width = 6, height = 5)
print(km_plot)
dev.off()

# ---- 3. EGFR/KRAS通路活性与signature评分的相关性分析 ----

# 探究不同基因集与EGFR和KRAS信号通路活性的相关性:
#   - KRAS基因集: KRAS_induced& suppressed.csv
#   - EGFR基因集(GMT): KOBAYASHI_EGFR_SIGNALING_24HR_UP/DN.v2026.1.Hs.gmt
#     来源: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/
#           KOBAYASHI_EGFR_SIGNALING_24HR_UP.html (DN同理)
#   - KRAS基因集来源: https://www.nature.com/articles/s41586-019-1884-x#Sec44
#   - 单细胞数据集: KrasG12Ci_scRNA_GSE137912
#   - Bulk RNA数据集: KRYSTAL-1_trial

# 3.1 读取和准备基因集
cat("=== 步骤1: 读取和准备基因集 ===\n")

# 读取KRAS诱导/抑制基因集
cat("读取KRAS基因集...\n")
kras_genes_df <- read.csv(file.path(DIR_DATA_CLINVAL, "KRAS_induced& suppressed.csv"),
                          stringsAsFactors = FALSE)
kras_induced_genes <- unique(na.omit(kras_genes_df$induced))
kras_suppressed_genes <- unique(na.omit(kras_genes_df$suppressed))
kras_induced_genes <- kras_induced_genes[kras_induced_genes != ""]
kras_suppressed_genes <- kras_suppressed_genes[kras_suppressed_genes != ""]
cat(sprintf("  KRAS诱导基因: %d个\n", length(kras_induced_genes)))
cat(sprintf("  KRAS抑制基因: %d个\n", length(kras_suppressed_genes)))

# 读取EGFR基因集(GMT格式); 注意原脚本将UP文件读入EGFR_suppressed_genes、
# DN文件读入EGFR_induced_genes, 按原逻辑保留
parse_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  gene_set <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    name <- parts[1]
    genes <- parts[-(1:2)]
    genes <- genes[genes != ""]
    gene_set[[name]] <- genes
  }
  return(gene_set)
}

cat("读取EGFR基因集...\n")
EGFR_suppressed_genes <- parse_gmt(file.path(
  DIR_DATA_CLINVAL, "KOBAYASHI_EGFR_SIGNALING_24HR_UP.v2026.1.Hs.gmt"))
EGFR_induced_genes <- parse_gmt(file.path(
  DIR_DATA_CLINVAL, "KOBAYASHI_EGFR_SIGNALING_24HR_DN.v2026.1.Hs.gmt"))
EGFR_suppressed_genes <- EGFR_suppressed_genes[[1]]
EGFR_induced_genes <- EGFR_induced_genes[[1]]
cat(sprintf("  EGFR下调基因: %d个\n", length(EGFR_suppressed_genes)))
cat(sprintf("  EGFR上调基因: %d个\n", length(EGFR_induced_genes)))

# 读取单细胞分析signature列表
cat("读取单细胞分析基因集...\n")
signature_list <- readRDS(file.path(DIR_OBJ_survival,
                                    "human_all_signature_without_pattern_list.rds"))
cat(sprintf("  单细胞基因集数量: %d个\n", length(signature_list)))

# 整合通路基因集与单细胞signature
cat("整合所有基因集...\n")
all_gene_sets <- c(
  list(
    KRAS_induced = kras_induced_genes,
    KRAS_suppressed = kras_suppressed_genes,
    EGFR_induced = EGFR_induced_genes,
    EGFR_suppressed = EGFR_suppressed_genes
  ),
  signature_list
)
cat(sprintf("  总基因集数量: %d个\n", length(all_gene_sets)))

# 保存基因集信息表
gene_set_info <- data.frame(
  gene_set = names(all_gene_sets),
  n_genes = sapply(all_gene_sets, length),
  category = c(rep("Pathway", 4), rep("Single_cell", length(signature_list))),
  stringsAsFactors = FALSE
)
write.csv(gene_set_info, file.path(DIR_TABLE_score, "gene_set_info.csv"), row.names = FALSE)

# 3.2 辅助函数: 基因集评分与相关性检验

# 按基因集平均表达计算评分(基因名大写匹配, 支持稀疏矩阵)
calculate_gene_set_scores <- function(expr_mat, gene_sets, min_genes = 3) {
  expr_genes <- rownames(expr_mat)
  expr_genes_upper <- toupper(expr_genes)
  score_df <- data.frame(row.names = colnames(expr_mat))

  for (gs_name in names(gene_sets)) {
    genes_upper <- unique(toupper(na.omit(gene_sets[[gs_name]])))
    gene_idx <- unique(match(genes_upper, expr_genes_upper))
    gene_idx <- gene_idx[!is.na(gene_idx)]

    if (length(gene_idx) >= min_genes) {
      # Matrix::colMeans 可直接处理稀疏矩阵, 避免as.matrix造成内存膨胀
      scores <- as.numeric(Matrix::colMeans(expr_mat[gene_idx, , drop = FALSE]))
      score_df[[gs_name]] <- scores
    } else {
      cat(sprintf(
        "  警告: %s 匹配基因不足 (%d/%d)，跳过\n",
        gs_name, length(gene_idx), length(genes_upper)
      ))
    }
  }

  score_df
}

# 提取Seurat对象表达矩阵(兼容Assay5层结构)
get_expr_matrix <- function(se_obj, assay = "RNA", layer = "data") {
  if (inherits(se_obj[[assay]], "Assay5")) {
    return(LayerData(se_obj, assay = assay, layer = layer))
  }
  GetAssayData(se_obj, assay = assay, layer = layer)
}

# 安全z-score: 标准差为0或不可计算时返回全0并告警
safe_zscore <- function(x, variable_name = "score") {
  x <- as.numeric(x)
  x_sd <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(x_sd) || x_sd == 0) {
    warning(sprintf("%s 的标准差为0或不可计算，z-score暂设为0。", variable_name))
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / x_sd
}

# 基因名大写匹配至表达矩阵行号
match_gene_indices <- function(genes, expr_genes) {
  idx <- unique(match(unique(toupper(na.omit(genes))), toupper(expr_genes)))
  idx[!is.na(idx)]
}

# 按平均表达过滤基因集; 低表达过滤后基因不足时回退到全部匹配基因
filter_gene_set_by_expression <- function(
    expr_mat,
    genes,
    low_expr_cutoff = 0.1,
    min_genes = 3
) {
  expr_genes <- rownames(expr_mat)
  idx <- match_gene_indices(genes, expr_genes)

  if (length(idx) == 0) {
    return(character(0))
  }

  gene_avg_expr <- Matrix::rowMeans(expr_mat)
  keep_idx <- idx[is.finite(gene_avg_expr[idx]) & gene_avg_expr[idx] >= low_expr_cutoff]

  # 低表达过滤后不足时, 退回到所有匹配基因, 但仍保证组成分量和差值评分使用同一套基因
  if (length(keep_idx) < min_genes && length(idx) >= min_genes) {
    warning(sprintf(
      "低表达过滤后仅剩%d个基因，少于min_genes=%d；改用全部%d个匹配基因。",
      length(keep_idx), min_genes, length(idx)
    ))
    keep_idx <- idx
  }

  expr_genes[keep_idx]
}

# 基因集平均表达评分(带名称向量)
score_gene_set <- function(expr_mat, genes, min_genes = 3) {
  idx <- match_gene_indices(genes, rownames(expr_mat))

  if (length(idx) < min_genes) {
    return(NULL)
  }

  score <- as.numeric(Matrix::colMeans(expr_mat[idx, , drop = FALSE]))
  names(score) <- colnames(expr_mat)
  score
}

# 双向通路评分: induced与suppressed分量均使用同一套过滤后基因,
# 差值评分 = zscore(induced) - zscore(suppressed)
compute_bidirectional_pathway_score <- function(
    meta_data,
    expr_mat,
    induced_genes,
    suppressed_genes,
    prefix,
    low_expr_cutoff = 0.1,
    min_genes = 3
) {
  induced_col <- paste0(prefix, "_induced")
  suppressed_col <- paste0(prefix, "_suppressed")
  raw_col <- paste0(prefix, "_score_raw")
  score_col <- paste0(prefix, "_score")

  induced_filtered <- filter_gene_set_by_expression(
    expr_mat, induced_genes, low_expr_cutoff, min_genes
  )
  suppressed_filtered <- filter_gene_set_by_expression(
    expr_mat, suppressed_genes, low_expr_cutoff, min_genes
  )

  induced_score <- score_gene_set(expr_mat, induced_filtered, min_genes)
  suppressed_score <- score_gene_set(expr_mat, suppressed_filtered, min_genes)

  if (is.null(induced_score) || is.null(suppressed_score)) {
    stop(sprintf(
      "%s通路评分失败：induced或suppressed匹配基因不足%d个。",
      prefix, min_genes
    ))
  }

  # 保留此前使用全部匹配基因计算的原始分量, 便于追溯
  if (induced_col %in% colnames(meta_data)) {
    meta_data[[paste0(induced_col, "_unfiltered")]] <- meta_data[[induced_col]]
  }
  if (suppressed_col %in% colnames(meta_data)) {
    meta_data[[paste0(suppressed_col, "_unfiltered")]] <- meta_data[[suppressed_col]]
  }

  induced_score <- induced_score[rownames(meta_data)]
  suppressed_score <- suppressed_score[rownames(meta_data)]

  # 热图中的两个组成分量与差值评分使用完全相同的过滤后基因
  meta_data[[induced_col]] <- as.numeric(induced_score)
  meta_data[[suppressed_col]] <- as.numeric(suppressed_score)

  # 同时保存原始差值和标准化差值
  meta_data[[raw_col]] <- meta_data[[induced_col]] - meta_data[[suppressed_col]]
  meta_data[[score_col]] <-
    safe_zscore(meta_data[[induced_col]], induced_col) -
    safe_zscore(meta_data[[suppressed_col]], suppressed_col)

  gene_info <- data.frame(
    pathway = prefix,
    component = c("induced", "suppressed"),
    n_input = c(length(unique(induced_genes)), length(unique(suppressed_genes))),
    n_matched_filtered = c(length(induced_filtered), length(suppressed_filtered)),
    genes = c(
      paste(induced_filtered, collapse = ";"),
      paste(suppressed_filtered, collapse = ";")
    ),
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "  %s_induced: %d -> %d；%s_suppressed: %d -> %d\n",
    prefix, length(unique(induced_genes)), length(induced_filtered),
    prefix, length(unique(suppressed_genes)), length(suppressed_filtered)
  ))

  list(
    metadata = meta_data,
    induced_genes = induced_filtered,
    suppressed_genes = suppressed_filtered,
    gene_info = gene_info
  )
}

# 对评分数据框追加双向通路差值评分
add_bidirectional_score_df <- function(score_df, prefix) {
  induced_col <- paste0(prefix, "_induced")
  suppressed_col <- paste0(prefix, "_suppressed")

  if (!all(c(induced_col, suppressed_col) %in% colnames(score_df))) {
    warning(sprintf("%s评分缺少induced或suppressed分量，未构建差值评分。", prefix))
    return(score_df)
  }

  score_df[[paste0(prefix, "_score_raw")]] <-
    score_df[[induced_col]] - score_df[[suppressed_col]]
  score_df[[paste0(prefix, "_score")]] <-
    safe_zscore(score_df[[induced_col]], induced_col) -
    safe_zscore(score_df[[suppressed_col]], suppressed_col)

  score_df
}

# 安全Spearman相关检验(样本量/变异不足时返回NA)
safe_cor_test <- function(x, y, method = "spearman") {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.frame(n = length(x), correlation = NA_real_, p_value = NA_real_))
  }

  test <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  data.frame(
    n = length(x),
    correlation = unname(test$estimate),
    p_value = test$p.value
  )
}

# 偏Spearman相关(秩残差法校正协变量)
partial_spearman_test <- function(data, x_var, y_var, covariates) {
  vars <- unique(c(x_var, y_var, covariates))
  vars <- vars[vars %in% colnames(data)]
  covariates <- setdiff(vars, c(x_var, y_var))

  df <- data[, c(x_var, y_var, covariates), drop = FALSE]
  df <- df[complete.cases(df), , drop = FALSE]

  if (nrow(df) < 10 || length(unique(df[[x_var]])) < 2 || length(unique(df[[y_var]])) < 2) {
    return(data.frame(
      n = nrow(df), correlation = NA_real_, p_value = NA_real_,
      adjusted_for = paste(covariates, collapse = ";")
    ))
  }

  x_rank <- rank(df[[x_var]], ties.method = "average")
  y_rank <- rank(df[[y_var]], ties.method = "average")

  if (length(covariates) == 0) {
    result <- safe_cor_test(x_rank, y_rank, method = "pearson")
  } else {
    cov_df <- df[, covariates, drop = FALSE]
    cov_df[] <- lapply(cov_df, function(v) {
      if (is.numeric(v)) rank(v, ties.method = "average") else factor(v)
    })

    design <- model.matrix(~ ., data = cov_df)
    x_resid <- qr.resid(qr(design), x_rank)
    y_resid <- qr.resid(qr(design), y_rank)
    result <- safe_cor_test(x_resid, y_resid, method = "pearson")
  }

  result$adjusted_for <- paste(covariates, collapse = ";")
  result
}

# 单细胞评分诊断: QC相关矩阵、全局/偏相关、分组内相关性并输出表格
run_sc_score_diagnostics <- function(
    meta_data,
    signature_var,
    pathway_prefix,
    dataset_label,
    table_dir,
    group_vars = character(0),
    qc_vars = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    min_group_n = 20
) {
  pathway_vars <- c(
    paste0(pathway_prefix, "_induced"),
    paste0(pathway_prefix, "_suppressed"),
    paste0(pathway_prefix, "_score_raw"),
    paste0(pathway_prefix, "_score")
  )
  pathway_vars <- pathway_vars[pathway_vars %in% colnames(meta_data)]

  if (!signature_var %in% colnames(meta_data) || length(pathway_vars) == 0) {
    warning(sprintf("%s中缺少%s或%s评分，跳过诊断。", dataset_label, signature_var, pathway_prefix))
    return(invisible(NULL))
  }

  qc_present <- qc_vars[qc_vars %in% colnames(meta_data)]
  matrix_vars <- unique(c(signature_var, pathway_vars, qc_present))
  cor_matrix <- cor(
    meta_data[, matrix_vars, drop = FALSE],
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  write.csv(
    cor_matrix,
    file.path(table_dir, paste0(dataset_label, "_", signature_var, "_",
                                pathway_prefix, "_QC_correlation_matrix.csv"))
  )

  global_results <- do.call(rbind, lapply(pathway_vars, function(pathway_var) {
    result <- safe_cor_test(meta_data[[signature_var]], meta_data[[pathway_var]], "spearman")
    result$signature <- signature_var
    result$pathway_variable <- pathway_var
    result
  }))
  global_results$p_adj <- p.adjust(global_results$p_value, method = "BH")
  write.csv(
    global_results,
    file.path(table_dir, paste0(dataset_label, "_", signature_var, "_",
                                pathway_prefix, "_global_correlations.csv")),
    row.names = FALSE
  )

  partial_results <- do.call(rbind, lapply(pathway_vars, function(pathway_var) {
    result <- partial_spearman_test(
      meta_data, signature_var, pathway_var, covariates = qc_present
    )
    result$signature <- signature_var
    result$pathway_variable <- pathway_var
    result
  }))
  partial_results$p_adj <- p.adjust(partial_results$p_value, method = "BH")
  write.csv(
    partial_results,
    file.path(table_dir, paste0(dataset_label, "_", signature_var, "_",
                                pathway_prefix, "_partial_correlations.csv")),
    row.names = FALSE
  )

  group_vars <- group_vars[group_vars %in% colnames(meta_data)]
  grouped_results <- list()

  for (group_var in group_vars) {
    group_values <- as.character(meta_data[[group_var]])
    valid_groups <- unique(group_values[!is.na(group_values)])

    for (group_value in valid_groups) {
      idx <- which(!is.na(group_values) & group_values == group_value)
      if (length(idx) < min_group_n) next

      for (pathway_var in pathway_vars) {
        result <- safe_cor_test(
          meta_data[[signature_var]][idx],
          meta_data[[pathway_var]][idx],
          "spearman"
        )
        result$group_variable <- group_var
        result$group <- group_value
        result$signature <- signature_var
        result$pathway_variable <- pathway_var
        grouped_results[[length(grouped_results) + 1]] <- result
      }
    }
  }

  if (length(grouped_results) > 0) {
    grouped_results <- dplyr::bind_rows(grouped_results)
    grouped_results$p_adj <- p.adjust(grouped_results$p_value, method = "BH")
    write.csv(
      grouped_results,
      file.path(table_dir, paste0(dataset_label, "_", signature_var, "_",
                                  pathway_prefix, "_grouped_correlations.csv")),
      row.names = FALSE
    )
  }

  invisible(list(
    global = global_results,
    partial = partial_results,
    grouped = grouped_results,
    correlation_matrix = cor_matrix
  ))
}

# 两基因集重叠统计(交集大小、比例、Jaccard)
gene_set_overlap_stats <- function(set_a, set_b, name_a, name_b) {
  set_a <- unique(toupper(na.omit(set_a)))
  set_b <- unique(toupper(na.omit(set_b)))
  overlap <- intersect(set_a, set_b)
  union_set <- union(set_a, set_b)

  data.frame(
    gene_set_a = name_a,
    gene_set_b = name_b,
    n_a = length(set_a),
    n_b = length(set_b),
    overlap_n = length(overlap),
    overlap_fraction_a = ifelse(length(set_a) > 0, length(overlap) / length(set_a), NA_real_),
    overlap_fraction_b = ifelse(length(set_b) > 0, length(overlap) / length(set_b), NA_real_),
    jaccard = ifelse(length(union_set) > 0, length(overlap) / length(union_set), NA_real_),
    overlap_genes = paste(overlap, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

# 3.3 检查ast_Traj_Up与KRAS/EGFR组成基因的直接重叠
if ("ast_Traj_Up" %in% names(signature_list)) {
  overlap_diagnostics <- dplyr::bind_rows(
    gene_set_overlap_stats(signature_list[["ast_Traj_Up"]], kras_induced_genes,
                           "ast_Traj_Up", "KRAS_induced"),
    gene_set_overlap_stats(signature_list[["ast_Traj_Up"]], kras_suppressed_genes,
                           "ast_Traj_Up", "KRAS_suppressed"),
    gene_set_overlap_stats(kras_induced_genes, kras_suppressed_genes,
                           "KRAS_induced", "KRAS_suppressed"),
    gene_set_overlap_stats(signature_list[["ast_Traj_Up"]], EGFR_induced_genes,
                           "ast_Traj_Up", "EGFR_induced"),
    gene_set_overlap_stats(signature_list[["ast_Traj_Up"]], EGFR_suppressed_genes,
                           "ast_Traj_Up", "EGFR_suppressed")
  )
  write.csv(
    overlap_diagnostics,
    file.path(DIR_TABLE_score, "ast_Traj_Up_pathway_gene_overlap.csv"),
    row.names = FALSE
  )
}

# ---- 3.4 KRAS单细胞数据集分析 (KrasG12Ci_scRNA_GSE137912) ----

cat("\n=== 步骤3: KRAS单细胞数据集分析 ===\n")

# 读取已有的KRAS Seurat对象(原位于clinical_analysis/results/KRAS_clinical_analysis/,
# 由外部脚本生成, 需复制到data/clinical_validation/; 不存在则跳过本节)
kras_rds_path <- file.path(DIR_DATA_CLINVAL, "kras_se_obj_with_scores.rds")
if (file.exists(kras_rds_path)) {
  cat("读取已有的KRAS Seurat对象...\n")
  se_kras <- readRDS(kras_rds_path)
  cat(sprintf("  细胞数: %d, 基因数: %d\n", ncol(se_kras), nrow(se_kras)))
} else {
  cat("KRAS Seurat对象不存在，跳过KRAS单细胞分析\n")
  se_kras <- NULL
}

if (!is.null(se_kras)) {
  # 计算全部基因集评分并加入metadata
  cat("计算KRAS单细胞基因集评分...\n")
  expr_mat_kras <- get_expr_matrix(se_kras)
  scores_kras_sc <- calculate_gene_set_scores(expr_mat_kras, all_gene_sets)
  cat(sprintf("  成功计算 %d 个基因集评分\n", ncol(scores_kras_sc)))

  for (col in colnames(scores_kras_sc)) {
    se_kras@meta.data[[col]] <- scores_kras_sc[[col]]
  }

  # 用同一套过滤基因计算组成分量, 并以标准化差值构建方向性通路评分
  cat("KRAS单细胞: 计算一致且标准化的KRAS/EGFR通路评分...\n")
  low_expr_cutoff <- 0.1

  score_result_kras_kras <- compute_bidirectional_pathway_score(
    meta_data = se_kras@meta.data,
    expr_mat = expr_mat_kras,
    induced_genes = kras_induced_genes,
    suppressed_genes = kras_suppressed_genes,
    prefix = "KRAS",
    low_expr_cutoff = low_expr_cutoff
  )
  se_kras@meta.data <- score_result_kras_kras$metadata
  kras_induced_filtered_kras <- score_result_kras_kras$induced_genes
  kras_suppressed_filtered_kras <- score_result_kras_kras$suppressed_genes

  score_result_kras_egfr <- compute_bidirectional_pathway_score(
    meta_data = se_kras@meta.data,
    expr_mat = expr_mat_kras,
    induced_genes = EGFR_induced_genes,
    suppressed_genes = EGFR_suppressed_genes,
    prefix = "EGFR",
    low_expr_cutoff = low_expr_cutoff
  )
  se_kras@meta.data <- score_result_kras_egfr$metadata
  EGFR_induced_filtered_kras <- score_result_kras_egfr$induced_genes
  EGFR_suppressed_filtered_kras <- score_result_kras_egfr$suppressed_genes

  write.csv(
    dplyr::bind_rows(
      score_result_kras_kras$gene_info,
      score_result_kras_egfr$gene_info
    ),
    file.path(DIR_TABLE_score, "KRAS_sc_filtered_pathway_gene_info.csv"),
    row.names = FALSE
  )

  # 各基因集评分间Spearman相关矩阵热图
  cat("KRAS单细胞基因集相关性分析...\n")
  score_cols_kras <- c("KRAS_induced", "KRAS_suppressed", "KRAS_score_raw", "KRAS_score",
                       "EGFR_suppressed", "EGFR_induced", "EGFR_score_raw", "EGFR_score",
                       colnames(scores_kras_sc)[!colnames(scores_kras_sc) %in%
                                                   c("KRAS_induced", "KRAS_suppressed",
                                                     "EGFR_suppressed", "EGFR_induced")])
  score_cols_kras <- score_cols_kras[score_cols_kras %in% colnames(se_kras@meta.data)]

  cor_mat_kras_sc <- cor(se_kras@meta.data[, score_cols_kras],
                         use = "pairwise.complete.obs", method = "spearman")
  pdf(file.path(DIR_FIG_score, "KRAS_sc_correlation_heatmap.pdf"), width = 14, height = 12)
  corrplot(cor_mat_kras_sc, method = "color", type = "full",
           order = "hclust", tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("blue", "white", "red"))(100),
           title = "KRAS Single Cell - Gene Set Correlation",
           mar = c(0, 0, 2, 0))
  dev.off()
  write.csv(cor_mat_kras_sc, file.path(DIR_TABLE_score, "KRAS_sc_correlation_matrix.csv"))

  # 额外诊断: 共同表达/QC混杂、QC校正后的偏相关以及分组内相关性
  run_sc_score_diagnostics(
    meta_data = se_kras@meta.data,
    signature_var = "ast_Traj_Up",
    pathway_prefix = "KRAS",
    dataset_label = "KRAS_SC",
    table_dir = DIR_TABLE_score,
    group_vars = c("orig.ident", "Condition", "Hours", "cell_line", "time_point")
  )
  run_sc_score_diagnostics(
    meta_data = se_kras@meta.data,
    signature_var = "ast_Traj_Up",
    pathway_prefix = "EGFR",
    dataset_label = "KRAS_SC",
    table_dir = DIR_TABLE_score,
    group_vars = c("orig.ident", "Condition", "Hours", "cell_line", "time_point")
  )

  # KRAS/EGFR通路评分与各基因集的相关性(热图)
  pathway_scores_kras <- c("KRAS_induced", "KRAS_suppressed", "KRAS_score_raw", "KRAS_score",
                           "EGFR_induced", "EGFR_suppressed", "EGFR_score_raw", "EGFR_score")
  pathway_scores_kras <- pathway_scores_kras[pathway_scores_kras %in% colnames(se_kras@meta.data)]
  sig_scores_kras <- score_cols_kras[!score_cols_kras %in%
                                       c("KRAS_induced", "KRAS_suppressed", "KRAS_score_raw",
                                         "KRAS_score", "EGFR_suppressed", "EGFR_induced",
                                         "EGFR_score_raw", "EGFR_score")]

  if (length(pathway_scores_kras) > 0 && length(sig_scores_kras) > 0) {
    cor_with_pathway_kras_sc <- data.frame(gene_set = sig_scores_kras)
    for (ps in pathway_scores_kras) {
      cor_vals <- sapply(sig_scores_kras, function(gs) {
        cor(se_kras@meta.data[[ps]], se_kras@meta.data[[gs]],
            use = "pairwise.complete.obs", method = "spearman")
      })
      cor_with_pathway_kras_sc[[ps]] <- cor_vals
    }
    write.csv(cor_with_pathway_kras_sc,
              file.path(DIR_TABLE_score, "KRAS_sc_correlation_with_pathway.csv"),
              row.names = FALSE)

    cor_long_kras_sc <- reshape2::melt(cor_with_pathway_kras_sc, id.vars = "gene_set",
                                       variable.name = "pathway_score", value.name = "correlation")

    p_cor_kras_sc <- ggplot(cor_long_kras_sc, aes(x = pathway_score, y = gene_set, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-1, 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)) +
      labs(x = "Pathway Score", y = "Gene Set", fill = "Spearman Correlation",
           title = "KRAS Single Cell - Correlation with Pathway Scores")
    ggsave(file.path(DIR_FIG_score, "KRAS_sc_correlation_with_pathway.pdf"),
           p_cor_kras_sc, width = 8, height = 10, dpi = 300)
  }

  # 按细胞系和时间点展示评分动态变化
  if ("cell_line" %in% colnames(se_kras@meta.data) &&
      "time_point" %in% colnames(se_kras@meta.data)) {
    time_data <- se_kras@meta.data %>%
      select(all_of(c("cell_line", "time_point", score_cols_kras))) %>%
      pivot_longer(cols = -c(cell_line, time_point),
                   names_to = "gene_set", values_to = "score")

    p_time_kras <- ggplot(time_data %>%
                            filter(gene_set %in% c("KRAS_score", "EGFR_score",
                                                   head(sig_scores_kras, 6))),
                          aes(x = time_point, y = score, color = cell_line, group = cell_line)) +
      stat_summary(fun = mean, geom = "point") +
      stat_summary(fun = mean, geom = "line") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      facet_wrap(~ gene_set, scales = "free_y", ncol = 4) +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 8)) +
      labs(x = "Time Point", y = "Score", color = "Cell Line",
           title = "KRAS Single Cell - Scores Over Time")
    ggsave(file.path(DIR_FIG_score, "KRAS_sc_scores_over_time.pdf"),
           p_time_kras, width = 16, height = 12, dpi = 300)
  }

  # 基因集评分与通路评分的散点图(分面 + Spearman相关标注)
  cat("绘制KRAS单细胞基因集与通路评分散点图...\n")
  if (length(pathway_scores_kras) > 0 && length(sig_scores_kras) > 0) {
    scatter_data <- se_kras@meta.data %>%
      select(all_of(c(pathway_scores_kras, sig_scores_kras))) %>%
      pivot_longer(cols = all_of(sig_scores_kras),
                   names_to = "gene_set", values_to = "gene_set_score")

    for (ps in pathway_scores_kras) {
      plot_data <- scatter_data %>%
        select(all_of(c(ps, "gene_set", "gene_set_score"))) %>%
        rename(pathway_score = all_of(ps)) %>%
        filter(!is.na(pathway_score), !is.na(gene_set_score))

      p_scatter <- ggplot(plot_data, aes(x = gene_set_score, y = pathway_score)) +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 0.8) +
        stat_cor(method = "spearman", size = 3,
                 label.x.npc = "left", label.y.npc = "top") +
        facet_wrap(~ gene_set, scales = "free", ncol = 4) +
        theme_bw() +
        theme(axis.text = element_text(size = 8),
              strip.text = element_text(size = 8),
              panel.grid = element_line(linewidth = 0.2)) +
        labs(x = "Gene Set Score", y = paste(ps, "Score"),
             title = paste0("KRAS Single Cell - ", ps, " vs Gene Sets"))

      ggsave(file.path(DIR_FIG_score, paste0("KRAS_sc_scatter_", ps, ".pdf")),
             p_scatter, width = 16, height = 12, dpi = 300)
    }
    cat(sprintf("  成功生成 %d 个散点图\n", length(pathway_scores_kras)))
  }

  # 保存KRAS单细胞评分结果
  saveRDS(list(
    scores = se_kras@meta.data[, score_cols_kras, drop = FALSE],
    unfiltered_scores = scores_kras_sc,
    correlation_matrix = cor_mat_kras_sc,
    metadata = se_kras@meta.data
  ), file.path(DIR_OBJ_score, "KRAS_sc_analysis_results.rds"))
}

# ---- 3.5 KRYSTAL-1 Bulk RNA数据集分析 ----

cat("\n=== 步骤4: KRYSTAL-1 Bulk RNA数据集分析 ===\n")

# 读取表达与临床数据并匹配共同样本
cat("读取KRYSTAL-1数据...\n")
vst_counts_krystal <- readRDS(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "vst_counts.rds"))
clinical_krystal <- read.csv(file.path(DIR_OBJ_survival, "KRYSTAL-1_trial", "meta_dat.csv"),
                             header = TRUE, row.names = 1)
cat(sprintf("  样本数: %d, 基因数: %d\n",
            ncol(vst_counts_krystal), nrow(vst_counts_krystal)))
cat(sprintf("  临床数据样本数: %d\n", nrow(clinical_krystal)))

common_samples <- intersect(colnames(vst_counts_krystal), rownames(clinical_krystal))
cat(sprintf("  匹配样本数: %d\n", length(common_samples)))
vst_counts_krystal <- vst_counts_krystal[, common_samples]
clinical_krystal <- clinical_krystal[common_samples, ]

# 计算全部基因集评分及KRAS/EGFR双向通路差值评分
cat("计算KRYSTAL-1基因集评分...\n")
scores_krystal <- calculate_gene_set_scores(vst_counts_krystal, all_gene_sets)
cat(sprintf("  成功计算 %d 个基因集评分\n", ncol(scores_krystal)))

scores_krystal <- add_bidirectional_score_df(scores_krystal, "KRAS")
scores_krystal <- add_bidirectional_score_df(scores_krystal, "EGFR")

# 合并临床数据
scores_krystal$sample <- rownames(scores_krystal)
scores_krystal <- cbind(scores_krystal, clinical_krystal[rownames(scores_krystal), ])

# 各基因集评分间Spearman相关矩阵热图
cat("KRYSTAL-1基因集相关性分析...\n")
score_cols_krystal <- c("KRAS_induced", "KRAS_suppressed", "KRAS_score_raw", "KRAS_score",
                        "EGFR_suppressed", "EGFR_induced", "EGFR_score_raw", "EGFR_score",
                        colnames(scores_krystal)[!colnames(scores_krystal) %in%
                                                    c("KRAS_induced", "KRAS_suppressed",
                                                      "KRAS_score", "EGFR_suppressed",
                                                      "EGFR_induced", "EGFR_score", "sample",
                                                      colnames(clinical_krystal))])
score_cols_krystal <- score_cols_krystal[score_cols_krystal %in% colnames(scores_krystal)]

cor_mat_krystal <- cor(scores_krystal[, score_cols_krystal],
                       use = "pairwise.complete.obs", method = "spearman")
pdf(file.path(DIR_FIG_score, "KRYSTAL1_correlation_heatmap.pdf"), width = 14, height = 12)
corrplot(cor_mat_krystal, method = "color", type = "full",
         order = "hclust", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(100),
         title = "KRYSTAL-1 - Gene Set Correlation",
         mar = c(0, 0, 2, 0))
dev.off()
write.csv(cor_mat_krystal, file.path(DIR_TABLE_score, "KRYSTAL1_correlation_matrix.csv"))

# KRAS/EGFR通路评分与各基因集的相关性(热图)
pathway_scores_krystal <- c("KRAS_score_raw", "KRAS_score", "EGFR_score_raw", "EGFR_score")
pathway_scores_krystal <- pathway_scores_krystal[pathway_scores_krystal %in% colnames(scores_krystal)]
sig_scores_krystal <- score_cols_krystal[!score_cols_krystal %in%
                                           c("KRAS_induced", "KRAS_suppressed",
                                             "KRAS_score_raw", "KRAS_score",
                                             "EGFR_suppressed", "EGFR_induced",
                                             "EGFR_score_raw", "EGFR_score")]

if (length(pathway_scores_krystal) > 0 && length(sig_scores_krystal) > 0) {
  cor_with_pathway_krystal <- data.frame(gene_set = sig_scores_krystal)
  for (ps in pathway_scores_krystal) {
    cor_vals <- sapply(sig_scores_krystal, function(gs) {
      cor(scores_krystal[[ps]], scores_krystal[[gs]],
          use = "pairwise.complete.obs", method = "spearman")
    })
    cor_with_pathway_krystal[[ps]] <- cor_vals
  }
  write.csv(cor_with_pathway_krystal,
            file.path(DIR_TABLE_score, "KRYSTAL1_correlation_with_pathway.csv"),
            row.names = FALSE)

  cor_long_krystal <- reshape2::melt(cor_with_pathway_krystal, id.vars = "gene_set",
                                     variable.name = "pathway_score", value.name = "correlation")

  p_cor_krystal <- ggplot(cor_long_krystal,
                          aes(x = pathway_score, y = gene_set, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(x = "Pathway Score", y = "Gene Set", fill = "Spearman Correlation",
         title = "KRYSTAL-1 - Correlation with Pathway Scores")
  ggsave(file.path(DIR_FIG_score, "KRYSTAL1_correlation_with_pathway.pdf"),
         p_cor_krystal, width = 8, height = 10, dpi = 300)
}

# OS生存分析: 逐基因集按中位数分高低组, Cox回归计算HR并绘制KM曲线
cat("KRYSTAL-1生存分析 (OS)...\n")
if (all(c("os", "os_cnsr") %in% colnames(scores_krystal))) {
  survival_results_os <- data.frame()

  for (gs in score_cols_krystal) {
    gs_data <- scores_krystal %>%
      select(all_of(c(gs, "os", "os_cnsr"))) %>%
      filter(!is.na(os), !is.na(os_cnsr)) %>%
      mutate(high_low = ifelse(.[[gs]] >= median(.[[gs]], na.rm = TRUE), "High", "Low"))

    if (length(unique(gs_data$high_low)) >= 2 && nrow(gs_data) >= 10) {
      tryCatch({
        fit <- survfit(Surv(os, os_cnsr) ~ high_low, data = gs_data)
        cox_fit <- coxph(Surv(os, os_cnsr) ~ high_low, data = gs_data)
        cox_summary <- summary(cox_fit)

        survival_results_os <- rbind(survival_results_os, data.frame(
          gene_set = gs,
          HR = cox_summary$coefficients[1, "exp(coef)"],
          HR_low = cox_summary$conf.int[1, "lower .95"],
          HR_high = cox_summary$conf.int[1, "upper .95"],
          p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
          n_high = sum(gs_data$high_low == "High"),
          n_low = sum(gs_data$high_low == "Low")
        ))

        # 通路评分及前6个signature绘制KM曲线
        if (gs %in% c("KRAS_score", "EGFR_score") || gs %in% head(sig_scores_krystal, 6)) {
          p_km <- ggsurvplot(
            fit,
            data = gs_data,
            conf.int = TRUE,
            legend.labs = c("Low Score", "High Score"),
            legend = "right",
            main = paste("OS Kaplan-Meier -", gs),
            xlab = "OS (days)",
            ylab = "Survival Probability",
            palette = c("#E41A1C", "#377EB8")
          )
          ggsave(file.path(DIR_FIG_score, paste0("KRYSTAL1_OS_KM_", gs, ".pdf")),
                 p_km$plot, width = 10, height = 8, dpi = 300)
        }
      }, error = function(e) {
        cat(sprintf("  OS分析失败: %s - %s\n", gs, e$message))
      })
    }
  }

  if (nrow(survival_results_os) > 0) {
    survival_results_os$p_adj <- p.adjust(survival_results_os$p_value, method = "BH")
    write.csv(survival_results_os,
              file.path(DIR_TABLE_score, "KRYSTAL1_OS_survival_results.csv"),
              row.names = FALSE)
    cat(sprintf("  OS生存分析完成: %d 个基因集\n", nrow(survival_results_os)))
  }
}

# PFS生存分析: 同OS流程
cat("KRYSTAL-1生存分析 (PFS)...\n")
if (all(c("pfsind", "pfsind_cnsr") %in% colnames(scores_krystal))) {
  survival_results_pfs <- data.frame()

  for (gs in score_cols_krystal) {
    gs_data <- scores_krystal %>%
      select(all_of(c(gs, "pfsind", "pfsind_cnsr"))) %>%
      filter(!is.na(pfsind), !is.na(pfsind_cnsr)) %>%
      mutate(high_low = ifelse(.[[gs]] >= median(.[[gs]], na.rm = TRUE), "High", "Low"))

    if (length(unique(gs_data$high_low)) >= 2 && nrow(gs_data) >= 10) {
      tryCatch({
        fit <- survfit(Surv(pfsind, pfsind_cnsr) ~ high_low, data = gs_data)
        cox_fit <- coxph(Surv(pfsind, pfsind_cnsr) ~ high_low, data = gs_data)
        cox_summary <- summary(cox_fit)

        survival_results_pfs <- rbind(survival_results_pfs, data.frame(
          gene_set = gs,
          HR = cox_summary$coefficients[1, "exp(coef)"],
          HR_low = cox_summary$conf.int[1, "lower .95"],
          HR_high = cox_summary$conf.int[1, "upper .95"],
          p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
          n_high = sum(gs_data$high_low == "High"),
          n_low = sum(gs_data$high_low == "Low")
        ))

        if (gs %in% c("KRAS_score", "EGFR_score") || gs %in% head(sig_scores_krystal, 6)) {
          p_km <- ggsurvplot(
            fit,
            data = gs_data,
            conf.int = TRUE,
            legend.labs = c("Low Score", "High Score"),
            legend = "right",
            main = paste("PFS Kaplan-Meier -", gs),
            xlab = "PFS (days)",
            ylab = "Survival Probability",
            palette = c("#E41A1C", "#377EB8")
          )
          ggsave(file.path(DIR_FIG_score, paste0("KRYSTAL1_PFS_KM_", gs, ".pdf")),
                 p_km$plot, width = 10, height = 8, dpi = 300)
        }
      }, error = function(e) {
        cat(sprintf("  PFS分析失败: %s - %s\n", gs, e$message))
      })
    }
  }

  if (nrow(survival_results_pfs) > 0) {
    survival_results_pfs$p_adj <- p.adjust(survival_results_pfs$p_value, method = "BH")
    write.csv(survival_results_pfs,
              file.path(DIR_TABLE_score, "KRYSTAL1_PFS_survival_results.csv"),
              row.names = FALSE)
    cat(sprintf("  PFS生存分析完成: %d 个基因集\n", nrow(survival_results_pfs)))
  }
}

# 保存KRYSTAL-1分析结果
saveRDS(list(
  scores = scores_krystal,
  correlation_matrix = cor_mat_krystal,
  clinical_data = clinical_krystal
), file.path(DIR_OBJ_score, "KRYSTAL1_analysis_results.rds"))

# ---- 3.6 EGFR单细胞数据集分析 ----
#     (供FigS5 D-F的EGFR-induced/EGFR-suppressed/ast_Traj_Up得分比较;
#      EGFR单细胞数据来自奥希替尼治疗前后配对的EGFR突变肺腺癌样本)
#     EGFR数据文件 egfr_epi_subset.rds / EGFR_Seu_metadata.csv 需置于
#     data/clinical_validation/ (与kras_se_obj_with_scores.rds同为外部脚本生成)

cat("\n=== 步骤3.6: EGFR单细胞数据集分析 ===\n")

egfr_rds_path <- file.path(DIR_DATA_CLINVAL, "egfr_epi_subset.rds")
egfr_meta_path <- file.path(DIR_DATA_CLINVAL, "EGFR_Seu_metadata.csv")
if (file.exists(egfr_rds_path) && file.exists(egfr_meta_path)) {
  cat("读取EGFR单细胞数据...\n")
  se_egfr <- readRDS(egfr_rds_path)
  se_egfr <- UpdateSeuratObject(se_egfr)
  metadata_egfr <- read.csv(egfr_meta_path, row.names = 1)
  se_egfr@meta.data <- cbind(se_egfr@meta.data, metadata_egfr[colnames(se_egfr), , drop = FALSE])
  cat(sprintf("  细胞数: %d, 基因数: %d\n", ncol(se_egfr), nrow(se_egfr)))

  cat("计算EGFR单细胞基因集评分...\n")
  expr_mat_egfr <- get_expr_matrix(se_egfr)
  scores_egfr_sc <- calculate_gene_set_scores(expr_mat_egfr, all_gene_sets)
  cat(sprintf("  成功计算 %d 个基因集评分\n", ncol(scores_egfr_sc)))

  for (col in colnames(scores_egfr_sc)) {
    se_egfr@meta.data[[col]] <- scores_egfr_sc[[col]]
  }

  # 用同一套过滤基因计算KRAS/EGFR双向通路评分
  cat("EGFR单细胞: 计算一致且标准化的KRAS/EGFR通路评分...\n")
  low_expr_cutoff <- 0.1

  score_result_egfr_kras <- compute_bidirectional_pathway_score(
    meta_data = se_egfr@meta.data, expr_mat = expr_mat_egfr,
    induced_genes = kras_induced_genes, suppressed_genes = kras_suppressed_genes,
    prefix = "KRAS", low_expr_cutoff = low_expr_cutoff
  )
  se_egfr@meta.data <- score_result_egfr_kras$metadata
  kras_induced_filtered <- score_result_egfr_kras$induced_genes
  kras_suppressed_filtered <- score_result_egfr_kras$suppressed_genes

  score_result_egfr_egfr <- compute_bidirectional_pathway_score(
    meta_data = se_egfr@meta.data, expr_mat = expr_mat_egfr,
    induced_genes = EGFR_induced_genes, suppressed_genes = EGFR_suppressed_genes,
    prefix = "EGFR", low_expr_cutoff = low_expr_cutoff
  )
  se_egfr@meta.data <- score_result_egfr_egfr$metadata
  EGFR_induced_filtered <- score_result_egfr_egfr$induced_genes
  EGFR_suppressed_filtered <- score_result_egfr_egfr$suppressed_genes

  write.csv(
    dplyr::bind_rows(score_result_egfr_kras$gene_info, score_result_egfr_egfr$gene_info),
    file.path(DIR_TABLE_score, "EGFR_sc_filtered_pathway_gene_info.csv"), row.names = FALSE
  )

  # 各基因集评分间Spearman相关矩阵热图
  cat("EGFR单细胞基因集相关性分析...\n")
  score_cols_egfr <- c("KRAS_induced", "KRAS_suppressed", "KRAS_score_raw", "KRAS_score",
                       "EGFR_suppressed", "EGFR_induced", "EGFR_score_raw", "EGFR_score",
                       colnames(scores_egfr_sc)[!colnames(scores_egfr_sc) %in%
                                                   c("KRAS_induced", "KRAS_suppressed",
                                                     "EGFR_suppressed", "EGFR_induced")])
  score_cols_egfr <- score_cols_egfr[score_cols_egfr %in% colnames(se_egfr@meta.data)]

  cor_mat_egfr_sc <- cor(se_egfr@meta.data[, score_cols_egfr],
                         use = "pairwise.complete.obs", method = "spearman")
  pdf(file.path(DIR_FIG_score, "EGFR_sc_correlation_heatmap.pdf"), width = 14, height = 12)
  corrplot(cor_mat_egfr_sc, method = "color", type = "full", order = "hclust",
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("blue", "white", "red"))(100),
           title = "EGFR Single Cell - Gene Set Correlation", mar = c(0, 0, 2, 0))
  dev.off()
  write.csv(cor_mat_egfr_sc, file.path(DIR_TABLE_score, "EGFR_sc_correlation_matrix.csv"))

  # 额外诊断: 共同表达/QC混杂、QC校正后的偏相关以及分组内相关性
  run_sc_score_diagnostics(meta_data = se_egfr@meta.data, signature_var = "ast_Traj_Up",
                           pathway_prefix = "KRAS", dataset_label = "EGFR_SC",
                           table_dir = DIR_TABLE_score,
                           group_vars = c("orig.ident", "analysis", "treatment_status", "biopsy_time_status"))
  run_sc_score_diagnostics(meta_data = se_egfr@meta.data, signature_var = "ast_Traj_Up",
                           pathway_prefix = "EGFR", dataset_label = "EGFR_SC",
                           table_dir = DIR_TABLE_score,
                           group_vars = c("orig.ident", "analysis", "treatment_status", "biopsy_time_status"))

  # 保存EGFR单细胞评分结果
  saveRDS(list(scores = se_egfr@meta.data[, score_cols_egfr, drop = FALSE],
               unfiltered_scores = scores_egfr_sc,
               correlation_matrix = cor_mat_egfr_sc,
               metadata = se_egfr@meta.data),
          file.path(DIR_OBJ_score, "EGFR_sc_analysis_results.rds"))
} else {
  cat("EGFR单细胞数据文件不存在，跳过EGFR单细胞分析\n")
  se_egfr <- NULL
}

# ---- 4. signature 表达箱线图(FigS5) ----
#     FigS5 A-C: KRAS单细胞中 KRAS_induced / KRAS_suppressed / ast_Traj_Up 得分
#                在KRASG12C抑制剂处理 0/4/24/72 h 的比较;
#     FigS5 D-F: EGFR单细胞中 EGFR_induced / EGFR_suppressed / ast_Traj_Up 得分
#                在治疗前(naive) vs 奥希替尼治疗后(post-osimertinib)的比较;
#     组间两两比较采用Wilcoxon秩和检验; 显著性: * P<0.05, ** P<0.01, *** P<0.001, **** P<0.0001

# ---- 4.1 通用绘图辅助函数 ----
default_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
                    "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC")

make_group_colors <- function(group_levels) {
  n_groups <- length(group_levels)
  if (n_groups <= length(default_colors)) {
    colors <- default_colors[seq_len(n_groups)]
  } else {
    colors <- grDevices::colorRampPalette(default_colors)(n_groups)
  }
  stats::setNames(colors, group_levels)
}

save_pdf <- function(plot_object, filename, width = 16, height = 8) {
  ggsave(filename = file.path(DIR_FIG_score, filename), plot = plot_object,
         width = width, height = height, units = "in", device = "pdf", limitsize = FALSE)
}

# 通用分组箱线图: 不显示单细胞散点与离群点, 分组着色, 两两Wilcoxon检验并标注显著性
plot_pairwise_boxplot <- function(
    seurat_obj,
    group_var,
    features,
    group_order = NULL,
    comparisons = NULL,
    ncol = 4,
    test_method = "wilcox.test",
    hide_ns = TRUE,
    x_label = NULL,
    y_labels = NULL,
    x_tick_labels = NULL,
    title_labels = NULL,
    y_limit_method = c("whisker", "quantile", "max"),
    y_quantile = c(0.01, 0.99),
    bottom_padding_ratio = 0.05,
    bracket_step_ratio = 0.08,
    top_padding_ratio = 0.06
) {
  y_limit_method <- match.arg(y_limit_method)
  meta_data <- seurat_obj@meta.data

  if (!group_var %in% colnames(meta_data)) {
    stop(paste0("Seurat 元数据中不存在分组变量：", group_var))
  }
  missing_features <- setdiff(features, colnames(meta_data))
  if (length(missing_features) > 0) {
    stop(paste0("Seurat 元数据中不存在以下变量：", paste(missing_features, collapse = ", ")))
  }

  group_data <- meta_data[[group_var]]

  if (is.null(group_order)) {
    if (is.factor(group_data)) {
      group_order <- levels(droplevels(group_data))
    } else {
      group_order <- unique(as.character(group_data))
      group_order <- group_order[!is.na(group_order)]
    }
  } else {
    group_order <- as.character(group_order)
  }

  observed_groups <- unique(as.character(group_data))
  observed_groups <- observed_groups[!is.na(observed_groups)]
  group_order <- group_order[group_order %in% observed_groups]

  if (length(group_order) < 2) {
    stop(paste0(group_var, " 中有效分组少于两个，无法进行组间比较。"))
  }

  group_colors <- make_group_colors(group_order)

  get_custom_label <- function(label_vector, variable, default_value) {
    if (is.null(label_vector)) return(default_value)
    if (length(label_vector) == 1 && is.null(names(label_vector))) return(as.character(label_vector))
    if (!is.null(names(label_vector)) && variable %in% names(label_vector)) return(as.character(label_vector[[variable]]))
    default_value
  }

  plot_list <- lapply(features, function(feature) {
    plot_df <- FetchData(seurat_obj, vars = c(group_var, feature))
    colnames(plot_df) <- c(".group", ".value")
    plot_df <- plot_df[complete.cases(plot_df), , drop = FALSE]
    if (nrow(plot_df) == 0) {
      warning(paste0(feature, " 没有可用于绘图的数据。"))
      return(NULL)
    }
    plot_df$.group <- as.character(plot_df$.group)
    feature_groups <- group_order[group_order %in% unique(plot_df$.group)]
    if (length(feature_groups) == 0) return(NULL)
    plot_df <- plot_df[plot_df$.group %in% feature_groups, , drop = FALSE]
    plot_df$.group <- factor(plot_df$.group, levels = feature_groups)

    if (is.null(comparisons)) {
      current_comparisons <- combn(feature_groups, 2, simplify = FALSE)
    } else {
      current_comparisons <- Filter(
        function(x) length(x) == 2 && all(x %in% feature_groups),
        lapply(comparisons, as.character)
      )
    }

    values <- plot_df$.value
    if (y_limit_method == "whisker") {
      group_values <- split(values, plot_df$.group)
      whisker_limits <- lapply(group_values, function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) return(c(NA_real_, NA_real_))
        box_stats <- boxplot.stats(x, coef = 1.5)$stats
        c(lower = box_stats[1], upper = box_stats[5])
      })
      whisker_limits <- unlist(whisker_limits)
      y_data_min <- min(whisker_limits[grepl("lower", names(whisker_limits))], na.rm = TRUE)
      y_data_max <- max(whisker_limits[grepl("upper", names(whisker_limits))], na.rm = TRUE)
    } else if (y_limit_method == "quantile") {
      y_quantiles <- quantile(values, probs = y_quantile, na.rm = TRUE, names = FALSE)
      y_data_min <- y_quantiles[1]
      y_data_max <- y_quantiles[2]
    } else {
      y_data_min <- min(values, na.rm = TRUE)
      y_data_max <- max(values, na.rm = TRUE)
    }

    y_range <- y_data_max - y_data_min
    if (!is.finite(y_range) || y_range <= 0) {
      scale_reference <- max(abs(c(y_data_min, y_data_max)), na.rm = TRUE)
      if (!is.finite(scale_reference) || scale_reference == 0) scale_reference <- 1
      y_range <- scale_reference * 0.2
    }

    n_comparisons <- length(current_comparisons)
    if (n_comparisons > 0) {
      bracket_step <- y_range * bracket_step_ratio
      label_y_positions <- y_data_max + bracket_step * seq_len(n_comparisons)
      y_upper <- max(label_y_positions) + y_range * top_padding_ratio
    } else {
      label_y_positions <- NULL
      y_upper <- y_data_max + y_range * top_padding_ratio
    }
    y_lower <- y_data_min - y_range * bottom_padding_ratio

    x_title <- if (is.null(x_label)) group_var else x_label
    y_title <- get_custom_label(label_vector = y_labels, variable = feature, default_value = feature)

    p <- ggplot(plot_df, aes(x = .group, y = .value, fill = .group)) +
      geom_boxplot(width = 0.65, linewidth = 0.65, alpha = 0.90,
                   outlier.shape = NA, na.rm = TRUE, show.legend = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(title = NULL, x = x_title, y = y_title) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      theme_classic(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.margin = margin(t = 15, r = 8, b = 8, l = 8))

    if (!is.null(x_tick_labels)) {
      p <- p + scale_x_discrete(labels = x_tick_labels, drop = FALSE)
    }

    if (n_comparisons > 0) {
      p <- p + stat_compare_means(
        comparisons = current_comparisons, method = test_method,
        label = "p.signif", hide.ns = hide_ns, label.y = label_y_positions,
        step.increase = 0, tip.length = 0.01, bracket.size = 0.4, size = 4.5,
        symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                           symbols = c("****", "***", "**", "*", "ns"))
      )
    }

    p <- p + coord_cartesian(ylim = c(y_lower, y_upper), clip = "off")
    p
  })

  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) stop("没有可用于绘图的数据。")
  wrap_plots(plot_list, ncol = ncol)
}

# ---- 4.2 KRAS单细胞: 按Hours比较(FigS5 A-C) ----
if (!is.null(se_kras)) {
  features_kras <- c("KRAS_induced", "KRAS_suppressed", "ast_Traj_Up")
  hours_order <- c("0", "4", "24", "72")
  se_kras$Hours <- factor(as.character(se_kras$Hours), levels = hours_order)

  p_kras_hours <- plot_pairwise_boxplot(
    seurat_obj = se_kras,
    group_var = "Hours",
    features = features_kras,
    group_order = hours_order,
    ncol = 3,
    hide_ns = FALSE,
    y_limit_method = "max"
  )
  print(p_kras_hours)
  save_pdf(p_kras_hours, "KRAS_sc_signature_hours_boxplot.pdf", width = 12, height = 8)
}

# ---- 4.3 EGFR单细胞: 按treatment_status比较(FigS5 D-F) ----
if (!is.null(se_egfr)) {
  se_filter <- se_egfr[, !is.na(se_egfr$treatment_status)]
  target_order <- c("na•ve", "post-osimertinib")
  se_filter$treatment_status <- factor(se_filter$treatment_status, levels = target_order)

  features_egfr <- c("EGFR_induced", "EGFR_suppressed", "ast_Traj_Up")
  p_egfr_treatment <- plot_pairwise_boxplot(
    seurat_obj = se_filter,
    group_var = "treatment_status",
    x_tick_labels = c("post-osimertinib" = "Post-Osimertinib", "na•ve" = "Naive"),
    features = features_egfr,
    y_limit_method = "max",
    ncol = 3,
    hide_ns = TRUE
  )
  print(p_egfr_treatment)
  save_pdf(p_egfr_treatment, "EGFR_sc_signature_treatment_status_boxplot.pdf",
           width = 12, height = 8)
}

# ---- 5. 小鼠类器官signature表达分析 (FigS5 G) ----
# 来源脚本: clinical_analysis/customized_signature/scripts/01_mouse_organoid_expression.R
# 内容: 鼠源signature(ast_Traj_Up/ast_td_vs_gfp/Nrg1/Trp63)在小鼠肺类器官
#       A10/A8(LUAD样)/P1/A7(SCC样)样本中的表达差异;
#       signature得分 = 基因集成员 log2(TPM+1) 均值;
#       组间差异采用ANOVA(BH校正); 另输出两两Wilcoxon秩和检验表,
#       对应FigS5 G箱线图(发布版)的两两比较显著性标注
# 输入: data/gene_orthology/gene_signatures_list_mouse.rds
#       data/clinical_validation/merged_tpm_organoid_processed.csv
# 输出: output/figures/survival/score_analysis/mouse_organoid/
#       output/tables/score_analysis/mouse_organoid/
#       output/object/survival/score_analysis/mouse_organoid_results.rds

# 输出子目录
DIR_FIG_mouse  <- file.path(DIR_FIG_score, "mouse_organoid")
DIR_TABLE_mouse <- file.path(DIR_TABLE_score, "mouse_organoid")
dir.create(DIR_FIG_mouse, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLE_mouse, recursive = TRUE, showWarnings = FALSE)

# 分组配色与顺序(与原脚本一致)
GROUP_COLORS_MOUSE <- c(A10 = "#4DAF4A", P1 = "#984EA3",
                        A8 = "#377EB8", A7 = "#E41A1C")
GROUP_LEVELS_MOUSE <- c("A10", "A8", "P1", "A7")

# ---- 5.1 读取数据 ----
signature_list <- readRDS(file.path(DIR_DATA, "gene_orthology",
                                    "gene_signatures_list_mouse.rds"))
cat("基因集数量:", length(signature_list), "\n")
print(names(signature_list))

expr_data <- read.csv(file.path(DIR_DATA_CLINVAL, "merged_tpm_organoid_processed.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_data) <- expr_data[[1]]
expr_data[[1]] <- NULL
cat(sprintf("表达数据维度: %d 基因 x %d 样本\n", nrow(expr_data), ncol(expr_data)))

# ---- 5.2 样本分组 (按样本名前缀) ----
sample_info <- data.frame(
  sample = colnames(expr_data),
  group = case_when(
    grepl("^A7", colnames(expr_data)) ~ "A7",
    grepl("^A8", colnames(expr_data)) ~ "A8",
    grepl("^A10", colnames(expr_data)) ~ "A10",
    grepl("^P1", colnames(expr_data)) ~ "P1",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)
sample_info$group <- factor(sample_info$group, levels = GROUP_LEVELS_MOUSE)
print(sample_info)

# ---- 5.3 计算基因集表达得分 (log2(TPM+1) 均值) ----
expr_log2 <- log2(as.matrix(expr_data) + 1)

signature_scores <- matrix(NA_real_, nrow = length(signature_list), ncol = ncol(expr_data),
                           dimnames = list(names(signature_list), colnames(expr_data)))
overlap_info <- list()
for (sig_name in names(signature_list)) {
  genes <- signature_list[[sig_name]]
  genes_present <- intersect(genes, rownames(expr_log2))
  overlap_info[[sig_name]] <- data.frame(signature = sig_name,
                                         n_input = length(genes),
                                         n_found = length(genes_present),
                                         stringsAsFactors = FALSE)
  if (length(genes_present) > 0) {
    signature_scores[sig_name, ] <- colMeans(expr_log2[genes_present, , drop = FALSE])
  }
  cat(sprintf("基因集 %s: %d/%d 基因匹配\n", sig_name, length(genes_present), length(genes)))
}
write.csv(do.call(rbind, overlap_info), file.path(DIR_TABLE_mouse, "mouse_signature_overlap.csv"),
          row.names = FALSE)

signature_scores_df <- as.data.frame(signature_scores)
signature_scores_df$signature <- rownames(signature_scores_df)

signature_scores_long <- signature_scores_df %>%
  pivot_longer(cols = -signature, names_to = "sample", values_to = "score") %>%
  left_join(sample_info, by = "sample")

# ---- 5.4 可视化 ----
annotation_col <- data.frame(Group = sample_info$group, row.names = sample_info$sample)
ann_colors <- list(Group = GROUP_COLORS_MOUSE)

# 5.4.1 基因集表达热图 (按行标准化)
pdf(file.path(DIR_FIG_mouse, "mouse_signature_heatmap.pdf"), width = 8, height = 5)
pheatmap(signature_scores,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Mouse Organoid - Signature Expression Across Groups")
dev.off()

# 5.4.2 箱线图
p_box <- ggplot(signature_scores_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = GROUP_COLORS_MOUSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Group", y = "Signature Score (log2 TPM)", title = "Mouse Organoid - Signature Expression by Group")
ggsave(file.path(DIR_FIG_mouse, "mouse_signature_boxplot.pdf"), p_box,
       width = 12, height = 4, dpi = 300)

# 5.4.3 小提琴图
p_vio <- ggplot(signature_scores_long, aes(x = group, y = score, fill = group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = GROUP_COLORS_MOUSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Group", y = "Signature Score (log2 TPM)", title = "Mouse Organoid - Signature Expression Distribution")
ggsave(file.path(DIR_FIG_mouse, "mouse_signature_violin.pdf"), p_vio,
       width = 12, height = 4, dpi = 300)

# 5.4.4 组平均表达热图
mean_by_group <- signature_scores_long %>%
  group_by(signature, group) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  column_to_rownames("signature")

pdf(file.path(DIR_FIG_mouse, "mouse_signature_heatmap_by_group.pdf"), width = 6, height = 5)
pheatmap(as.matrix(mean_by_group),
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_cols = FALSE,
         main = "Mouse Organoid - Average Signature Expression by Group")
dev.off()

# ---- 5.5 统计分析 (ANOVA, 组间差异) ----
stat_results <- signature_scores_long %>%
  group_by(signature) %>%
  summarise(
    p_value = summary(aov(score ~ group, data = pick(group, score)))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = ifelse(p_adj < 0.05, "Significant (p_adj<0.05)", "Not significant"))
write.csv(stat_results, file.path(DIR_TABLE_mouse, "mouse_signature_anova_results.csv"),
          row.names = FALSE)
print(stat_results)

# 5.5b 两两Wilcoxon秩和检验(补充): 发布版FigS5 G箱线图的两两比较显著性标注
#      即以两两Wilcoxon检验(*号标注)为依据, 与原脚本的ANOVA结果并存
group_levels_v <- levels(sample_info$group)
pairwise_comps <- combn(group_levels_v, 2, simplify = FALSE)
pairwise_wilcox <- lapply(names(signature_list), function(sig_name) {
  sub_df <- signature_scores_long %>%
    filter(signature == sig_name, !is.na(score), !is.na(group))
  do.call(rbind, lapply(pairwise_comps, function(comp) {
    x <- sub_df$score[sub_df$group == comp[1]]
    y <- sub_df$score[sub_df$group == comp[2]]
    wt <- tryCatch(wilcox.test(x, y), error = function(e) NULL)
    data.frame(signature = sig_name,
               group1 = comp[1], group2 = comp[2],
               p_value = if (is.null(wt)) NA_real_ else wt$p.value,
               stringsAsFactors = FALSE)
  }))
})
pairwise_wilcox_df <- do.call(rbind, pairwise_wilcox)
pairwise_wilcox_df$p_adj <- p.adjust(pairwise_wilcox_df$p_value, method = "BH")
write.csv(pairwise_wilcox_df, file.path(DIR_TABLE_mouse, "mouse_signature_pairwise_wilcoxon.csv"),
          row.names = FALSE)
print(pairwise_wilcox_df)

# ---- 5.6 保存结果 ----
write.csv(signature_scores, file.path(DIR_TABLE_mouse, "mouse_signature_scores_matrix.csv"))
write.csv(signature_scores_long, file.path(DIR_TABLE_mouse, "mouse_signature_scores_detailed.csv"),
          row.names = FALSE)
saveRDS(list(scores = signature_scores,
             long = signature_scores_long,
             anova = stat_results,
             labels = sample_info),
        file.path(DIR_OBJ_score, "mouse_organoid_results.rds"))

cat("\n[5] 类器官 (mouse) 表达分析完成\n")
