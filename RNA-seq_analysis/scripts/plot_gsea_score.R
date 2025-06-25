#!/usr/bin/env Rscript
# 脚本目的: 绘制GSEA score-plot (running score)

# 1. 加载包和函数--------------------------------
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(enrichplot)
  library(RColorBrewer)
  library(cowplot)
  library(ggtext)
  library(clusterProfiler)
})
source("functions/config_and_themes.R")

# 2. 定义路径--------------------------------
gsea_go_path <- "results/RData/GSEA_GO_results.RData"
plot_path <- "results/plots/GSEA_score_plot.pdf"

# 3. 绘制GSEA-plot--------------------------------
if (file.exists(gsea_go_path)) {
  load(gsea_go_path)

  # 选择与气泡图相同的GO terms，以保持一致性
  gsea_go_terms_to_plot <- c(
    "negative regulation of leukocyte mediated cytotoxicity", "acute inflammatory response",
    "negative regulation of interleukin-6 production", "negative regulation of natural killer cell mediated cytotoxicity",
    "negative regulation of T cell activation", "positive regulation of apoptotic process",
    "transferase activity, transferring one-carbon groups", "tetrahydrofolate metabolic process",
    "nuclear chromosome segregation", "DNA-templated DNA replication"
  )

  plot_data <- gsea_go %>% filter(Description %in% gsea_go_terms_to_plot) %>%
    arrange(NES, .by_group = T)

  if (nrow(as.data.frame(plot_data)) > 0) {
    p_gsea_score <- gseaplot2(
      plot_data,
      geneSetID = grep("DNA|tetrahydrofolate", plot_data@result$Description),
      base_size = 8,
      color = brewer.pal(3, "Dark2")[1:2],
      pvalue_table = F,
      rel_heights = c(1.5, .2, .5)
    )

    ggsave(plot_path, p_gsea_score, units = "cm", width = 10, height = 7)
    cat("GSEA score图已保存到:", plot_path, "\n")
  } else {
    cat("在富集结果中未找到指定的GSEA GO terms，无法绘图。\n")
  }
} else {
  cat("GSEA GO富集结果文件不存在，无法绘图。\n")
}
