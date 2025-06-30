#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2025-06-25
# 功能：读取RNA-seq下调基因和Cut&Tag靶基因，计算交集并绘制韦恩图

# 1. 加载包 ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(eulerr)
  library(grid)
})

setwd("D:/2023_Git/cut&tag/20230308结果")
source("./Integrative_analysis/functions/config_and_themes.R")

# 2. 定义路径和参数 ------------------------------------------------------------
rnaseq_file <- "./RNA-seq_analysis/results/tables/DE_results.csv"
cuttag_file <- "./CutTag_local_analysis/peak_annotation_macs3/proximal_promoter_genes/all_groups_common_proximal_genes.txt"

output_dir <- "./Integrative_analysis/results"
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# 3. 读取并处理基因列表 --------------------------------------------------------
cat("INFO: 读取并处理基因列表...\n")

down_in_cre_genes <- read.csv(rnaseq_file, row.names = 1) %>%
  filter(regulation == "Down_in_Cre") %>%
  rownames()

cat("INFO: Cre下调的差异基因:", length(down_in_cre_genes), "个\n")

if (!file.exists(cuttag_file)) {
  stop("ERROR: 找不到Cut&Tag靶基因文件，请检查路径:", cuttag_file)
}
cuttag_genes <- read.table(cuttag_file, header = FALSE, stringsAsFactors = FALSE)$V1
cat("INFO: Cut&Tag靶基因:", length(cuttag_genes), "个\n")

# 4. 计算交集并保存 ------------------------------------------------------------
cat("INFO: 计算基因交集...\n")
intersect_genes <- intersect(down_in_cre_genes, cuttag_genes)
cat("INFO: 交集基因:", length(intersect_genes), "个\n")

write.table(intersect_genes,
            file.path(output_dir, "intersect_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("INFO: 交集基因已保存至:", file.path(output_dir, "intersect_genes.txt"), "\n")

# 5. 绘制韦恩图 ---------------------------------------------------------------
cat("INFO: 绘制韦恩图...\n")
venn_list <- list(
  'Down in Cre' = down_in_cre_genes,
  'IKZF2 Targets' = cuttag_genes
)
fit <- euler(venn_list)

formatted_labels <- c(expression(atop("Down in", italic(Ikzf2^{"fl/fl"}*";Lck-cre"))), 'IKZF2 Targets')

pdf(file.path(plots_dir, "Venn_diagram_down_cre_vs_h123.pdf"), width = 8 / 2.54, height = 4 / 2.54)
plot(fit,
     fills = list(fill = c("#EFC000FF", "#0073C2FF"), alpha = 0.7),
     quantities = list(fontsize = 8),
     legend = list(labels = formatted_labels, fontsize = 8, font = 1))
dev.off()
cat("INFO: 韦恩图已保存至:", file.path(plots_dir, "Venn_diagram_down_cre_vs_h123.pdf"), "\n")

cat("INFO: 脚本 01_intersect_genes.R 执行完毕。\n")
