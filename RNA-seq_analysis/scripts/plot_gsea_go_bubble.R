#!/usr/bin/env Rscript
# 脚本目的: 绘制GSEA GO富集分析的气泡图

# 1. 加载包和函数--------------------------------
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(cowplot)
  library(ggtext)
  library(RColorBrewer)
})
source("functions/config_and_themes.R")

# 2. 定义路径--------------------------------
gsea_go_path <- "results/RData/GSEA_GO_results.RData"
plot_path <- "results/plots/GSEA_GO_bubble_plot.pdf"

# 3. 绘制GSEA GO富集气泡图--------------------------------
if (file.exists(gsea_go_path)) {
  load(gsea_go_path)
  gsea_go_terms_to_plot <- c(
    "negative regulation of leukocyte mediated cytotoxicity", "acute inflammatory response",
    "negative regulation of interleukin-6 production", "negative regulation of natural killer cell mediated cytotoxicity",
    "negative regulation of T cell activation", "positive regulation of apoptotic process",
    "transferase activity, transferring one-carbon groups", "tetrahydrofolate metabolic process",
    "nuclear chromosome segregation", "DNA-templated DNA replication"
  )

  plot_data <- gsea_go %>%
    filter(Description %in% gsea_go_terms_to_plot) %>%
    arrange(NES, .by_group = T)

  if (nrow(as.data.frame(plot_data)) > 0) {
    lim <- round(max(abs(plot_data@result$NES))) + 1
    p_gsea_go <- dotplot(plot_data, showCategory=20, split=".sign", x = "NES", label_format = 40) +
      facet_grid(.sign~., scales = "free", space = "free") +
      scale_fill_distiller(palette = "YlOrRd", limits = c(NA, 0.05), breaks = c(0.01, 0.03, 0.05)) +
      scale_x_continuous(limits = c(-lim, lim)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = "GSEA GO Enrichment<br><i>Ikzf2<sup>fl/fl</sup>;Lck-cre</i> vs <i>Ikzf2<sup>fl/fl</sup></i>", x = "Normalized Enrichment Score (NES)") +
      my_theme
    p_gsea_go

    ggsave(plot_path, p_gsea_go, units = "cm", width = 10, height = 7)
    cat("GSEA GO气泡图已保存到:", plot_path, "\n")
  } else {
    cat("在富集结果中未找到指定的GSEA GO terms，无法绘图。\n")
  }
} else {
  cat("GSEA GO富集结果文件不存在，无法绘图。\n")
}
