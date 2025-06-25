#!/usr/bin/env Rscript

# 加载必要的包
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(RColorBrewer)
})

# 设置工作目录和文件路径
base_dir = "D:/2023_Git/cut&tag/20230308结果/CutTag_local_analysis"
samples = c("H1_vs_I1", "H1_vs_I2", 
           "H2_vs_I1", "H2_vs_I2",
           "H3_vs_I1", "H3_vs_I2")  # SEACR结果文件名

# 创建输出目录
output_dir = paste0(base_dir, "/peak_annotation_seacr")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "/plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "/tables"), showWarnings = FALSE, recursive = TRUE)

# 加载TxDb
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

# ===== 第一部分：数据生成 =====
cat("===== 开始数据生成 =====\n")

# 1. 读取所有样本的峰文件
cat("读取峰文件...\n")
peak_files = paste0(base_dir, "/individual_peaks/", samples, ".stringent.bed")
names(peak_files) = samples
peaks_list = lapply(peak_files, readPeakFile)

cat("读取完成，各样本峰数量:\n")
for(sample in samples) {
  cat(sample, ":", length(peaks_list[[sample]]), "\n")
}

# 2. 获取启动子区域
cat("\n获取启动子区域...\n")
promoter = getPromoters(TxDb = txdb, upstream = 5000, downstream = 5000)
cat("启动子区域获取完成\n")

# 3. 计算所有样本峰在启动子区域的分布
cat("计算峰在启动子区域的分布...\n")
tagMatrix_list = lapply(peaks_list, getTagMatrix, windows = promoter)
cat("峰分布计算完成\n")

# 4. 注释所有样本的峰
cat("注释峰...\n")
peakAnno_list = lapply(peaks_list, annotatePeak, TxDb = txdb,
                       tssRegion = c(-5000, 5000),
                       annoDb = "org.Mm.eg.db")
cat("峰注释完成\n")

# ===== 第二部分：保存结果 =====
cat("\n===== 开始保存结果 =====\n")

# 1. 保存RDS文件
cat("保存RDS文件...\n")
# 保存峰注释结果
saveRDS(peakAnno_list, 
        paste0(output_dir, "/tables/peak_annotation_list.rds"))
# 保存TSS分布矩阵
saveRDS(tagMatrix_list,
        paste0(output_dir, "/tables/tag_matrix_list.rds"))
cat("RDS文件已保存\n")

# 2. 保存每个样本的峰信息
cat("保存峰信息...\n")
for(sample in samples) {
  peak_stats = as.data.frame(peaks_list[[sample]])
  write.csv(peak_stats,
            paste0(output_dir, "/tables/", sample, "_peak_info.csv"),
            row.names = FALSE)
}
cat("峰信息已保存\n")

# 3. 保存每个样本的注释结果
cat("保存注释结果...\n")
for(sample in samples) {
  anno_df = as.data.frame(peakAnno_list[[sample]])
  write.csv(anno_df,
            paste0(output_dir, "/tables/", sample, "_peak_annotation.csv"),
            row.names = FALSE)
}
cat("注释结果已保存\n")

# 3. 生成分析报告
cat("生成分析报告...\n")
sink(paste0(output_dir, "/peak_annotation_report.txt"))
cat("===== CUT&Tag SEACR峰注释分析报告 =====\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
for(sample in samples) {
  cat("样本:", sample, "\n")
  cat("峰数量:", length(peaks_list[[sample]]), "\n")
  cat("注释统计:\n")
  print(peakAnno_list[[sample]]@annoStat)
  cat("\n----------------------------------------\n\n")
}
cat("===== 分析完成 =====\n")
sink()
cat("分析报告已生成\n")

# ===== 第三部分：可视化 =====
cat("\n===== 开始生成可视化结果 =====\n")

# 1. 生成所有样本的TSS富集图
cat("生成TSS富集图...\n")
pdf_file = paste0(output_dir, "/plots/TSS_enrichment.pdf")
pdf(pdf_file, width = 8, height = 6)

# 生成颜色方案
colors = c(brewer.pal(3, "Reds")[3],    # H1_vs_I1
           brewer.pal(3, "Oranges")[3],  # H1_vs_I2
           brewer.pal(3, "Greens")[3],   # H2_vs_I1
           brewer.pal(3, "Blues")[3],    # H2_vs_I2
           brewer.pal(3, "Purples")[3],  # H3_vs_I1
           brewer.pal(3, "Greys")[3])    # H3_vs_I2

p = plotAvgProf(tagMatrix_list, xlim = c(-5000, 5000),
                xlab = "Genomic Region (5'->3')",
                ylab = "Read Count Frequency")

# 修改颜色和主题
p + scale_color_manual(values = colors) +
   theme_bw() +
   theme(legend.title = element_blank(),
         legend.position = "right",
         axis.text = element_text(size = 10),
         axis.title = element_text(size = 12))
dev.off()
cat("TSS富集图已保存到:", pdf_file, "\n")

# 2. 生成所有样本的注释饼图
cat("生成注释饼图...\n")
for(sample in samples) {
  pdf_file = paste0(output_dir, "/plots/", sample, "_annotation_pie.pdf")
  pdf(pdf_file, width = 8, height = 8)
  plotAnnoPie(peakAnno_list[[sample]])
  dev.off()
}
cat("饼图已保存\n")

cat("\n分析完成！结果保存在:", output_dir, "\n")
