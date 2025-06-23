#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2025-04-17
# 分析内容：取Cre下调差异基因和H1_H2_H3_overlap基因的交集，进行功能富集分析

# 加载必要的包-------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(DOSE)
  library(RColorBrewer)
  library(ggplot2)
  library(eulerr)
  library(grid)
  library(ggbreak)
  library(Cairo)
})

# 设置工作目录
# setwd("RNA-seq_analysis")  # 如需要，请取消注释

# 设置输出目录
output_dir <- "enrichment_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 设置全局字体大小为8pt
my_theme <- theme_classic() +
  theme(
    text = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "lines")
  )

# 颜色设置
cols <- brewer.pal(8, "Dark2")

cat("开始读取RNA-seq差异基因数据...\n")

# 读取edgeR差异分析结果------------------
deg_data <- read.csv("./edgeR_paired_Cre.act_vs_Flox.act_results.csv", row.names = 1)

# 提取Cre下调的差异基因(Down_in_Cre)
down_in_cre_genes <- rownames(deg_data[deg_data$regulation == "Down_in_Cre", ])
cat("Cre下调的差异基因:", length(down_in_cre_genes), "个\n")

cat("开始读取H1_H2_H3共有基因数据...\n")

# 读取H1_H2_H3_overlap基因------------------
# 注意：这里需要根据实际文件路径和格式进行调整
h123_file <- "../local_analysis/peak_annotation_macs3/proximal_promoter_genes/all_groups_common_proximal_genes.txt"

if(file.exists(h123_file)) {
  h123_genes <- read.table(h123_file, header = FALSE, stringsAsFactors = FALSE)$V1
  cat("H1_H2_H3共有基因:", length(h123_genes), "个\n")
} else {
  stop("找不到H1_H2_H3共有基因文件，请检查路径")
}

# 找出交集基因-----------------------
intersect_genes <- intersect(down_in_cre_genes, h123_genes)
cat("交集基因:", length(intersect_genes), "个\n")

# 绘制韦恩图------------------
cat("开始绘制韦恩图...\n")
venn_list <- list(
  'Down in Cre' = down_in_cre_genes,
  'IKZF2 Targets' = h123_genes
)

# 使用eulerr创建面积成比例的韦恩图
fit <- euler(venn_list)

# 绘制并保存韦恩图
pdf(file.path(output_dir, "Venn_diagram_down_cre_vs_h123.pdf"), width = 4, height = 3)
plot(fit,
     fills = list(fill = c("#EFC000FF", "#0073C2FF"), alpha = 0.7),
     labels = list(font = 1, fontsize = 8),
     quantities = list(fontsize = 8)
)
# 手动添加标题以精确控制字体大小
grid.text("Overlap of Down-regulated Genes and IKZF2 Targets",
          y = unit(.95, "npc"), # 调整标题位置
          gp = gpar(fontsize = 8))
dev.off()
cat("韦恩图已保存至:", file.path(output_dir, "Venn_diagram_down_cre_vs_h123.pdf"), "\n")

# 将交集基因保存到文件------------------
write.table(intersect_genes,
            file.path(output_dir, "Down_in_Cre_and_H123_overlap_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 将基因符号转换为Entrez ID，用于后续富集分析------------------
cat("基因ID转换...\n")
gene_ids <- bitr(intersect_genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Mm.eg.db)
cat("成功映射的基因:", nrow(gene_ids), "个\n")

# 如果没有基因匹配，则终止分析
if(nrow(gene_ids) == 0) {
  stop("没有基因成功映射到Entrez ID，请检查基因名称格式")
}

# 保存基因ID对应关系
write.csv(gene_ids, file.path(output_dir, "gene_id_mapping.csv"))

# 开始富集分析------------------
cat("开始GO富集分析...\n")

# GO富集分析 - 合并BP、MF、CC------------------
go_all <- enrichGO(gene = gene_ids$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",       # 分析所有GO类别
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)

# 保存GO富集结果
if(nrow(go_all@result) > 0) {
  write.csv(go_all@result, file.path(output_dir, "GO_enrichment.csv"))
  cat("GO富集结果:", nrow(go_all@result), "个条目\n")

  # 按照GO子类别划分结果
  go_bp_results <- go_all@result[go_all@result$ONTOLOGY == "BP", ]
  go_mf_results <- go_all@result[go_all@result$ONTOLOGY == "MF", ]
  go_cc_results <- go_all@result[go_all@result$ONTOLOGY == "CC", ]

  cat("其中 BP 类别:", nrow(go_bp_results), "个条目\n")
  cat("其中 MF 类别:", nrow(go_mf_results), "个条目\n")
  cat("其中 CC 类别:", nrow(go_cc_results), "个条目\n")
} else {
  cat("GO富集结果为空\n")
}

# 保存GO富集对象为RData (移到if循环外)
save(go_all, file = file.path(output_dir, "GO_enrichment.RData"))
write.csv(go_all@result, file = file.path(output_dir, "GO_enrichment.csv"))
cat("GO富集对象已保存为RData文件\n")

# KEGG富集分析------------------
cat("开始KEGG富集分析...\n")
kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                 organism = 'mmu',
                 keyType = "ncbi-geneid",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kegg = setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# 保存KEGG富集结果
if(nrow(kegg@result) > 0) {
  write.csv(kegg@result, file.path(output_dir, "KEGG_enrichment.csv"))
  cat("KEGG富集结果:", nrow(kegg@result), "个条目\n")
} else {
  cat("KEGG富集结果为空\n")
}

# 保存KEGG富集对象为RData (移到if循环外)
save(kegg, file = file.path(output_dir, "KEGG_enrichment.RData"))
write.csv(kegg@result, file = file.path(output_dir, "KEGG_enrichment.csv"))
cat("KEGG富集对象已保存为RData文件\n")

# 可视化富集分析结果------------------
cat("开始绘制富集分析结果...\n")

# ===== GO分析可视化 =====

# 对所有GO结果按p值排序
go_all@result <- go_all@result[order(go_all@result$p.adjust), ]

# 找出包含"carbon"的词条
carbon_terms <- grep("carbon", go_all@result$Description, ignore.case = TRUE)
carbon_terms = carbon_terms[c(1, 3)]
# 选择前5个条目和包含carbon的词条
top10_indices <- 1:min(2, nrow(go_all@result))
selected_indices <- unique(c(top10_indices, carbon_terms))
selected_indices <- selected_indices[selected_indices <= nrow(go_all@result)]

# 提取选中的GO条目
show_terms <- go_all@result$Description[selected_indices]
go_selected <- go_all@result[selected_indices, ]

# 准备GO条形图数据：转换p.adjust为-log10格式
go_barplot_data <- go_selected %>%
  arrange(p.adjust) %>%
  mutate(log_padj = -log10(p.adjust),
         Description = factor(Description, levels = rev(Description)))

# GO条形图 - 使用ggplot2直接绘制########################
p1 <- ggplot(go_barplot_data, aes(x = log_padj, y = Description)) +
  geom_bar(stat = "identity", width = 0.6, fill = "steelblue") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "GO Enrichment Analysis",
       x = "-log10(p.adjust)",
       y = NULL) +
  my_theme +
  theme(axis.text.y = element_text(size = 8, hjust = 1))

# 检查是否需要断轴
max_log_padj_go <- max(go_barplot_data$log_padj, na.rm = TRUE)
if (is.finite(max_log_padj_go) && max_log_padj_go > 5) {
  p1 <- p1 + scale_x_cut(breaks = c(5), which = 2, scales = 0.4)
}

# 强制关闭所有现有图形设备，确保从一个干净的状态开始
while(!is.null(dev.list())) dev.off()
# 使用cairo_pdf设备，这是解决罕见绘图问题的最终方法
cairo_pdf(file.path(output_dir, "GO_barplot.pdf"), width = 5.5, height = length(show_terms) * 0.35)
print(p1)
dev.off()

# 准备GO点图数据
go_dotplot_data <- go_selected %>%
  arrange(p.adjust) %>%
  mutate(log_padj = -log10(p.adjust),
         GeneRatio = parse_ratio(GeneRatio),
         Description = factor(Description, levels = rev(Description)))

# GO点图 - 使用ggplot2直接绘制########################
p2 <- ggplot(go_dotplot_data, aes(x = log_padj, y = Description, size = Count, color = log_padj)) +
  geom_point() +
  scale_color_gradientn(colours = colorRampPalette(brewer.pal(9, "YlOrRd"))(30)) +
  labs(title = "GO Enrichment Analysis",
       x = "-log10(p.adjust)",
       y = NULL,
       size = "Gene Count") +
  my_theme +
  theme(axis.text.y = element_text(size = 8, hjust = 1))

# 使用ggsave保存图片，不再需要手动print
ggsave(file.path(output_dir, "GO_dotplot.pdf"), p2, width = 6, height = min(8, max(4, length(show_terms) * 0.3)))

# GO富集网络图 - 使用更新的参数语法并确保不超出边界########################
# 先检查showCategory参数是否超出可用的GO条目数量
if(length(show_terms) > 0) {
  # 确保用于网络图的条目数不会太多，避免下标出界错误
  max_terms_for_network <- min(10, length(show_terms))
  network_terms <- show_terms[1:max_terms_for_network]

  # 使用新的参数语法，避免过时警告
  p3 <- emapplot(pairwise_termsim(go_all),
               showCategory = network_terms,
               cex.params = list(category_label = 0.6),  # 新的参数语法
               color = "p.adjust") +
    labs(title = "GO Term Similarity Network") +
    my_theme

  # 使用ggsave保存图片，不再需要手动print
  ggsave(file.path(output_dir, "GO_emapplot.pdf"), p3, width = 7, height = 6)
} else {
  cat("警告: 没有足够的GO条目创建网络图\n")
}

# ===== KEGG分析可视化 =====

# 处理KEGG通路名称，移除"- Mus musculus (house mouse)"部分
kegg@result$clean_name <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg@result$Description)

# 找出包含"carbon"的KEGG通路
carbon_kegg <- grep("carbon", kegg@result$clean_name, ignore.case = TRUE)
carbon_kegg = carbon_kegg[1]
# 选择前10个条目和包含carbon的词条
top10_kegg <- 1:min(2, nrow(kegg@result))
selected_kegg <- unique(c(top10_kegg, carbon_kegg))
selected_kegg <- selected_kegg[selected_kegg <= nrow(kegg@result)]

# 提取选中的KEGG条目
show_kegg_terms <- kegg@result$clean_name[selected_kegg]
kegg_selected <- kegg@result[selected_kegg, ]

# 准备KEGG条形图数据
kegg_barplot_data <- kegg_selected %>%
  mutate(Description = clean_name) %>%
  arrange(p.adjust) %>%
  mutate(log_padj = -log10(p.adjust),
         Description = factor(Description, levels = rev(Description)))

# KEGG条形图 - 使用ggplot2直接绘制 #########################
p6 <- ggplot(kegg_barplot_data, aes(x = log_padj, y = Description)) +
  geom_bar(stat = "identity", width = 0.6, fill = "steelblue") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "KEGG Pathway Enrichment",
       x = "-log10(p.adjust)",
       y = NULL) +
  my_theme +
  theme(axis.text.y = element_text(size = 8, hjust = 1))

# 检查是否需要断轴
max_log_padj_kegg <- max(kegg_barplot_data$log_padj, na.rm = TRUE)
if (is.finite(max_log_padj_kegg) && max_log_padj_kegg > 5) {
  p6 <- p6 + scale_x_cut(breaks = c(3), which = 2, scales = 0.4)
}

# 强制关闭所有现有图形设备，确保从一个干净的状态开始
while(!is.null(dev.list())) dev.off()
# 使用cairo_pdf设备，这是解决罕见绘图问题的最终方法
cairo_pdf(file.path(output_dir, "KEGG_barplot.pdf"), width = 5.5, height = length(show_terms) * 0.35)
print(p6)
dev.off()

# 准备KEGG点图数据
kegg_dotplot_data <- kegg_selected %>%
  mutate(Description = clean_name) %>%
  arrange(p.adjust) %>%
  mutate(log_padj = -log10(p.adjust),
         GeneRatio = parse_ratio(GeneRatio),
         Description = factor(Description, levels = rev(Description)))

# KEGG点图 - 使用ggplot2直接绘制#########################
p8 <- ggplot(kegg_dotplot_data, aes(x = log_padj, y = Description, size = Count, color = log_padj)) +
  geom_point() +
  scale_color_gradientn(colours = colorRampPalette(brewer.pal(9, "YlOrRd"))(30)) +
  labs(title = "KEGG Pathway Enrichment",
       x = "-log10(p.adjust)",
       y = NULL,
       size = "Gene Count") +
  my_theme +
  theme(axis.text.y = element_text(size = 8, hjust = 1))

# 使用ggsave保存图片，不再需要手动print
ggsave(file.path(output_dir, "KEGG_dotplot.pdf"), p8, width = 6, height = min(8, max(4, length(show_kegg_terms) * 0.3)))

# KEGG网络图 - 同样修正参数和错误处理#########################
if(length(show_kegg_terms) > 0) {
  # 确保用于网络图的条目数不会太多
  max_kegg_terms_for_network <- min(10, length(show_kegg_terms))
  network_kegg_terms <- show_kegg_terms[1:max_kegg_terms_for_network]

  # 使用新的参数语法
  p9 <- emapplot(pairwise_termsim(kegg),
               showCategory = network_kegg_terms,
               cex.params = list(category_label = 0.6),  # 新的参数语法
               color = "p.adjust") +
    labs(title = "KEGG Pathway Similarity Network") +
    my_theme

  # 使用ggsave保存图片，不再需要手动print
  ggsave(file.path(output_dir, "KEGG_emapplot.pdf"), p9, width = 7, height = 6)
} else {
  cat("警告: 没有足够的KEGG通路创建网络图\n")
}

# 添加辅助函数用于处理GeneRatio
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(i) as.numeric(i[1])/as.numeric(i[2]))
}

# 创建富集结果汇总报告
cat("创建富集分析汇总报告...\n")

sink(file.path(output_dir, "enrichment_summary.txt"))
cat("===== 差异基因与启动子区域峰值基因交集的功能富集分析报告 =====\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("基因集统计:\n")
cat("--------------------------------\n")
cat("Cre下调的差异基因数量:", length(down_in_cre_genes), "\n")
cat("H1_H2_H3共有的近端启动子基因数量:", length(h123_genes), "\n")
cat("交集基因数量:", length(intersect_genes), "\n")
cat("成功映射到Entrez ID的基因数量:", nrow(gene_ids), "\n\n")

cat("GO富集分析结果统计:\n")
cat("--------------------------------\n")
cat("GO富集条目总数:", nrow(go_all@result), "\n")
cat("其中 BP 类别:", nrow(go_bp_results), "个条目\n")
cat("其中 MF 类别:", nrow(go_mf_results), "个条目\n")
cat("其中 CC 类别:", nrow(go_cc_results), "个条目\n\n")

if(nrow(go_all@result) > 0) {
  # 排序并获取前15个GO条目
  sorted_go_result <- go_all@result[order(go_all@result$p.adjust), ]
  top_go_all <- sorted_go_result[1:min(15, nrow(sorted_go_result)), ]

  cat("Top 15 GO条目 (所有类别):\n")
  for(i in 1:nrow(top_go_all)) {
    cat(i, ". [", top_go_all$ONTOLOGY[i], "] ", top_go_all$Description[i],
        " (p-adjust=", format(top_go_all$p.adjust[i], digits=3), ")\n", sep="")
  }
}
cat("\n")

cat("KEGG通路富集分析结果统计:\n")
cat("--------------------------------\n")
if(!is.null(kegg)) {
  cat("KEGG通路富集条目数:", nrow(kegg@result), "\n")
  if(nrow(kegg@result) > 0) {
    top_kegg <- kegg@result[1:min(10, nrow(kegg@result)), ]
    cat("Top 10 KEGG通路:\n")
    for(i in 1:nrow(top_kegg)) {
      cat(i, ". ", top_kegg$Description[i], " (p-adjust=", format(top_kegg$p.adjust[i], digits=3), ")\n", sep="")
    }
  }
} else {
  cat("KEGG通路富集结果为空\n")
}
cat("\n")

cat("生物学意义解读:\n")
cat("--------------------------------\n")
cat("这些交集基因既是IKZF2条件敲除中下调表达的基因，同时又在H1、H2和H3样本中的近端启动子区域存在峰值，\n")
cat("表明这些基因可能是转录因子直接调控的靶基因，其表达受到IKZF2条件敲除的正向调控。\n")
cat("这种调控模式暗示，IKZF2条件敲除影响了这些基因的转录活性，导致其表达水平下降。\n\n")
cat("上述富集的GO条目和KEGG通路揭示了这些基因可能参与的生物学过程和功能，\n")
cat("为理解IKZF2条件敲除的分子机制和潜在的下游效应提供了重要线索。\n")
sink()

cat("\n富集分析完成，结果已保存至:", output_dir, "\n")
