#!/usr/bin/env Rscript

# MACS3峰值注释结果筛选脚本 - 提取启动子（Promoter）区域的峰值
# 作者：Claude
# 创建日期：2024-06-24

# 加载必要的包
suppressPackageStartupMessages({
  library(tidyverse)
  library(ChIPseeker)
})

# 设置工作目录
base_dir <- "D:/2023_Git/cut&tag/20230308结果/CutTag_local_analysis"
input_dir <- paste0(base_dir, "/peak_annotation_macs3/tables")
output_dir <- paste0(base_dir, "/peak_annotation_macs3/promoter_genes")

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 定义样本名称
samples <- c("H1_vs_I1", "H1_vs_I2",
            "H2_vs_I1", "H2_vs_I2",
            "H3_vs_I1", "H3_vs_I2")

# ===== 从注释文件中提取Promoter区域的峰值 =====
cat("开始筛选Promoter区域的峰值和基因...\n")

# 读取预先保存的峰值注释结果
peakAnno_list <- readRDS(paste0(input_dir, "/peak_annotation_list.rds"))
cat("已读取注释结果数据\n")

# 函数：从注释结果中提取启动子区域的峰值
extract_promoter_peaks <- function(anno_result) {
  # 获取完整注释数据
  anno_df <- as.data.frame(anno_result)

  # 筛选启动子区域的峰值
  # Promoter区域包括：Promoter (<=1kb), Promoter (1-2kb), Promoter (2-3kb),
  # Promoter (3-4kb), Promoter (4-5kb)
  promoter_df <- anno_df %>%
    filter(grepl("Promoter", annotation)) %>%
    # 添加一列提取Promoter的具体范围
    mutate(promoter_region = sub(".*Promoter \\(([^)]+)\\).*", "\\1", annotation))

  return(promoter_df)
}

# 为每个样本提取Promoter区域的峰值
promoter_list <- lapply(peakAnno_list, extract_promoter_peaks)

# 保存每个样本的Promoter峰值数据
for(sample in names(promoter_list)) {
  promoter_df <- promoter_list[[sample]]
  write.csv(promoter_df,
            paste0(output_dir, "/", sample, "_promoter_peaks.csv"),
            row.names = FALSE)

  cat(sample, ": 共", nrow(promoter_df), "个峰值位于启动子区域\n")
}

# ===== 提取唯一的基因信息 =====
cat("\n提取与启动子区域峰值相关的基因...\n")

# 提取与Promoter区域峰值相关的唯一基因
extract_promoter_genes <- function(promoter_df) {
  # 提取基因信息：基因ID、基因名称和描述
  genes_df <- promoter_df %>%
    select(geneId, SYMBOL, GENENAME) %>%
    distinct() %>%
    arrange(SYMBOL) %>%
    na.omit()

  return(genes_df)
}

# 为每个样本提取唯一的基因信息
genes_list <- lapply(promoter_list, extract_promoter_genes)

# 保存每个样本的基因信息
for(sample in names(genes_list)) {
  genes_df <- genes_list[[sample]]
  write.csv(genes_df,
            paste0(output_dir, "/", sample, "_promoter_genes.csv"),
            row.names = FALSE)

  cat(sample, ": 共", nrow(genes_df), "个唯一基因的启动子区域有峰值\n")
}

# ===== 提取样本间共有的基因 =====
cat("\n分析样本之间的共有基因...\n")

# 提取所有样本共有的基因
if(length(genes_list) > 1) {
  # 获取每个样本的基因ID列表
  gene_ids_list <- lapply(genes_list, function(df) df$geneId)

  # 找到所有样本共有的基因ID
  common_genes <- Reduce(intersect, gene_ids_list) %>% na.omit() %>% unique()

  # 如果有共有基因，则提取详细信息
  if(length(common_genes) > 0) {
    # 使用第一个样本的基因信息作为参考
    first_sample <- names(genes_list)[1]
    common_genes_df <- genes_list[[first_sample]] %>%
      filter(geneId %in% common_genes)

    # 保存共有基因信息
    write.csv(common_genes_df,
              paste0(output_dir, "/common_promoter_genes.csv"),
              row.names = FALSE)

    cat("所有样本共有", length(common_genes), "个基因的启动子区域有峰值\n")
  } else {
    cat("没有找到所有样本共有的基因\n")
  }
}

# ===== 生成Venn图显示样本间基因重叠关系 =====
cat("\n生成Venn图显示基因重叠关系...\n")

# 安装并加载VennDiagram包
if(!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
suppressPackageStartupMessages(library(VennDiagram))
# 设置venn图颜色
venn_colors <- RColorBrewer::brewer.pal(6, "Set2")

# 为不同组合的样本生成Venn图
# H1组vs H2组vs H3组
H_groups <- list(
  H1 = na.omit(unique(intersect(genes_list[["H1_vs_I1"]]$SYMBOL, genes_list[["H1_vs_I2"]]$SYMBOL))),
  H2 = na.omit(unique(intersect(genes_list[["H2_vs_I1"]]$SYMBOL, genes_list[["H2_vs_I2"]]$SYMBOL))),
  H3 = na.omit(unique(intersect(genes_list[["H3_vs_I1"]]$SYMBOL, genes_list[["H3_vs_I2"]]$SYMBOL)))
)

# 绘制H组间的Venn图
pdf(paste0(output_dir, "/H_groups_promoter_genes_venn.pdf"), width = 10, height = 8)
venn.plot <- venn.diagram(
  x = H_groups,
  filename = NULL,
  fill = venn_colors[1:3],
  alpha = 0.5,
  main = "H1 vs H2 vs H3 Groups - Promoter Genes",
  main.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

# 同一组不同IgG对照的比较
# H1组
H1_samples <- list(
  H1_vs_I1 = genes_list[["H1_vs_I1"]]$SYMBOL,
  H1_vs_I2 = genes_list[["H1_vs_I2"]]$SYMBOL
)

# 绘制H1组对照间的Venn图
pdf(paste0(output_dir, "/H1_samples_promoter_genes_venn.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = H1_samples,
  filename = NULL,
  fill = venn_colors[1:2],
  alpha = 0.5,
  main = "H1 vs I1 and H1 vs I2 - Promoter Genes",
  main.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

# H2组
H2_samples <- list(
  H2_vs_I1 = genes_list[["H2_vs_I1"]]$SYMBOL,
  H2_vs_I2 = genes_list[["H2_vs_I2"]]$SYMBOL
)

# 绘制H2组对照间的Venn图
pdf(paste0(output_dir, "/H2_samples_promoter_genes_venn.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = H2_samples,
  filename = NULL,
  fill = venn_colors[3:4],
  alpha = 0.5,
  main = "H2 vs I1 and H2 vs I2 - Promoter Genes",
  main.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

# H3组
H3_samples <- list(
  H3_vs_I1 = genes_list[["H3_vs_I1"]]$SYMBOL,
  H3_vs_I2 = genes_list[["H3_vs_I2"]]$SYMBOL
)

# 绘制H3组对照间的Venn图
pdf(paste0(output_dir, "/H3_samples_promoter_genes_venn.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = H3_samples,
  filename = NULL,
  fill = venn_colors[5:6],
  alpha = 0.5,
  main = "H3 vs I1 and H3 vs I2 - Promoter Genes",
  main.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

# ===== 生成汇总报告 =====
cat("\n生成汇总报告...\n")

# 创建汇总报告
sink(paste0(output_dir, "/promoter_peaks_summary.txt"))
cat("===== MACS3峰值启动子区域筛选分析报告 =====\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("每个样本的启动子区域峰值统计:\n")
cat("--------------------------------\n")
for(sample in names(promoter_list)) {
  # 计算不同启动子区域类型的峰值数量
  region_counts <- table(promoter_list[[sample]]$promoter_region)
  total_peaks <- nrow(promoter_list[[sample]])
  total_genes <- nrow(genes_list[[sample]])

  cat(sample, ":\n")
  cat("  - 总峰值数:", total_peaks, "\n")
  cat("  - 相关基因数:", total_genes, "\n")
  cat("  - 启动子区域分布:\n")
  for(region in names(region_counts)) {
    cat("      -", region, ":", region_counts[region],
        sprintf("(%.1f%%)", region_counts[region]/total_peaks*100), "\n")
  }
  cat("\n")
}

# 添加样本间比较信息
if(length(genes_list) > 1) {
  cat("样本间共有基因分析:\n")
  cat("--------------------------------\n")

  # 计算每组的共有基因（不同IgG对照的交集）
  H1_genes <- na.omit(unique(intersect(genes_list[["H1_vs_I1"]]$SYMBOL, genes_list[["H1_vs_I2"]]$SYMBOL)))
  H2_genes <- na.omit(unique(intersect(genes_list[["H2_vs_I1"]]$SYMBOL, genes_list[["H2_vs_I2"]]$SYMBOL)))
  H3_genes <- na.omit(unique(intersect(genes_list[["H3_vs_I1"]]$SYMBOL, genes_list[["H3_vs_I2"]]$SYMBOL)))

  # 输出各组基因数量
  cat("H1组共有基因数:", length(H1_genes), "\n")
  cat("H2组共有基因数:", length(H2_genes), "\n")
  cat("H3组共有基因数:", length(H3_genes), "\n")

  # 计算不同组合之间的重叠基因数量
  H1_H2_overlap <- length(intersect(H1_genes, H2_genes))
  H1_H3_overlap <- length(intersect(H1_genes, H3_genes))
  H2_H3_overlap <- length(intersect(H2_genes, H3_genes))
  H1_H2_H3_overlap <- length(Reduce(intersect, list(H1_genes, H2_genes, H3_genes)))

  cat("\n组间重叠基因数量:\n")
  cat("H1 ∩ H2:", H1_H2_overlap, "\n")
  cat("H1 ∩ H3:", H1_H3_overlap, "\n")
  cat("H2 ∩ H3:", H2_H3_overlap, "\n")
  cat("H1 ∩ H2 ∩ H3:", H1_H2_H3_overlap, "\n")
}
sink()

cat("\n分析完成，结果已保存到:", output_dir, "\n")
