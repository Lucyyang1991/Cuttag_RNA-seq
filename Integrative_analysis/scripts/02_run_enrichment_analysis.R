#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2025-06-25
# 功能：对交集基因进行GO和KEGG功能富集分析

# 1. 加载包 ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# 2. 定义路径 ----------------------------------------------------------------
setwd("D:/2023_Git/cut&tag/20230308结果")
input_file <- "./Integrative_analysis/results/intersect_genes.txt"
output_dir <- "./Integrative_analysis/results/enrichment"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 3. 读取交集基因 --------------------------------------------------------------
cat("INFO: 读取交集基因列表...\n")
if (!file.exists(input_file)) {
  stop("ERROR: 找不到交集基因文件，请先运行 01_intersect_genes.R 脚本。路径:", input_file)
}
intersect_genes <- read.table(input_file, header = FALSE)$V1
cat("INFO: 交集基因:", length(intersect_genes), "个\n")

# 4. 基因ID转换 ----------------------------------------------------------------
cat("INFO: 基因ID转换 (SYMBOL -> ENTREZID)...\n")
gene_ids <- bitr(intersect_genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Mm.eg.db,
                drop = TRUE)

if (nrow(gene_ids) == 0) {
  stop("ERROR: 没有基因成功映射到Entrez ID，请检查基因名称格式")
}
cat("INFO: 成功映射的基因:", nrow(gene_ids), "个\n")


# 5. GO富集分析 ----------------------------------------------------------------
cat("INFO: 开始GO富集分析...\n")
enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
) -> go_all

if (!is.null(go_all) && nrow(go_all@result) > 0) {
  write.csv(go_all@result, file.path(output_dir, "GO_enrichment.csv"), row.names = FALSE)
  save(go_all, file = file.path(output_dir, "GO_enrichment.RData"))
  cat("INFO: GO富集结果已保存。发现", nrow(go_all@result), "个条目\n")
} else {
  cat("INFO: GO富集结果为空\n")
}

# 6. KEGG富集分析 ---------------------------------------------------------------
cat("INFO: 开始KEGG富集分析...\n")
enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = 'mmu',
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
) -> kegg

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  kegg_readable <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  write.csv(kegg_readable@result, file.path(output_dir, "KEGG_enrichment.csv"), row.names = FALSE)
  save(kegg_readable, file = file.path(output_dir, "KEGG_enrichment.RData"))
  cat("INFO: KEGG富集结果已保存。发现", nrow(kegg@result), "个条目\n")
} else {
  cat("INFO: KEGG富集结果为空\n")
}

cat("INFO: 脚本 02_run_enrichment_analysis.R 执行完毕。\n") 
