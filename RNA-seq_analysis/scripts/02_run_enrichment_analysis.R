#!/usr/bin/env Rscript
# 脚本目的:
# 1. 读取差异分析结果。
# 2. 对差异基因进行ORA富集分析(GO & KEGG)。
# 3. 对所有基因按logFC排序后进行GSEA富集分析(GO & KEGG)。
# 4. 保存所有富集分析结果。

# 1. 加载包--------------------------------
suppressMessages({
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# 2. 定义输入输出路径--------------------------------
de_results_path <- "results/tables/DE_results.csv"

# 3. 读取差异分析结果--------------------------------
if (!file.exists(de_results_path)) {
  stop("差异分析结果文件不存在，请先运行 '01_run_DE_analysis.R'")
}
deg <- read.csv(de_results_path, row.names = 1)

# 4. ORA 富集分析 (GO & KEGG)--------------------------------
# 准备基因列表
up_genes   <- deg %>% filter(regulation == "Up_in_Cre") %>% rownames()
down_genes <- deg %>% filter(regulation == "Down_in_Cre") %>% rownames()
all_genes  <- rownames(deg)

# ID转换 (SYMBOL -> ENTREZID)
suppressMessages({
  up_entrez   <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  all_entrez  <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})

gene_list_ora <- list(Up_in_Cre = up_entrez, Down_in_Cre = down_entrez)
gene_list_ora <- gene_list_ora[sapply(gene_list_ora, function(x) length(x) > 0)]

if (length(gene_list_ora) > 0) {
  # GO ORA
  compare_go <- compareCluster(geneCluster = gene_list_ora, fun = "enrichGO", universe = all_entrez, OrgDb = org.Mm.eg.db, readable = TRUE, ont = "ALL")
  if (!is.null(compare_go) && nrow(as.data.frame(compare_go)) > 0) {
    write.csv(as.data.frame(compare_go), "results/tables/ORA_GO_results.csv", row.names = FALSE)
    save(compare_go, file = "results/RData/ORA_GO_results.RData")
    cat("ORA GO 分析完成。\n")
  }

  # KEGG ORA
  compare_kegg <- compareCluster(geneCluster = gene_list_ora, fun = "enrichKEGG", universe = all_entrez, organism = 'mmu')
  if (!is.null(compare_kegg) && nrow(as.data.frame(compare_kegg)) > 0) {
    compare_kegg <- setReadable(compare_kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    write.csv(as.data.frame(compare_kegg), "results/tables/ORA_KEGG_results.csv", row.names = FALSE)
    save(compare_kegg, file = "results/RData/ORA_KEGG_results.RData")
    cat("ORA KEGG 分析完成。\n")
  }
}

# 5. GSEA 富集分析 (GO & KEGG)--------------------------------
gene_list_gsea <- deg %>%
  rownames_to_column("SYMBOL") %>%
  dplyr::select(SYMBOL, logFC) %>%
  filter(!is.na(logFC)) %>%
  inner_join(
    suppressMessages(bitr(.$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)),
    by = "SYMBOL"
  ) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  pull(logFC, name = ENTREZID) %>%
  sort(decreasing = TRUE)
print(head(gene_list_gsea))

if(length(gene_list_gsea) > 0){
    # GSEA GO
    gsea_go <- gseGO(geneList = gene_list_gsea, OrgDb = org.Mm.eg.db, ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = T)
    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      gsea_go <- setReadable(gsea_go, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(as.data.frame(gsea_go), "results/tables/GSEA_GO_results.csv", row.names = FALSE)
      save(gsea_go, file = "results/RData/GSEA_GO_results.RData")
      cat("GSEA GO 分析完成。\n")
    }

    # GSEA KEGG
    gsea_kegg <- gseKEGG(geneList = gene_list_gsea, organism = 'mmu', pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = T)
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(as.data.frame(gsea_kegg), "results/tables/GSEA_KEGG_results.csv", row.names = FALSE)
      save(gsea_kegg, file = "results/RData/GSEA_KEGG_results.RData")
      cat("GSEA KEGG 分析完成。\n")
    }
}

cat("富集分析脚本执行完毕！\n")
