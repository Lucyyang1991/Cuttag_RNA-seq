#!/usr/bin/env Rscript
# 脚本目的: 绘制ORA GO富集分析的气泡图

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
ora_go_path <- "results/RData/ORA_GO_results.RData"
plot_path <- "results/plots/ORA_GO_bubble_plot.pdf"

# 3. 绘制ORA GO富集气泡图--------------------------------
if (file.exists(ora_go_path)) {
  load(ora_go_path)
  ora_go_terms_to_plot <- c(
    "response to decreased oxygen levels", "pyruvate metabolic process", "inflammatory response",
    "positive regulation of programmed cell death", "negative regulation of leukocyte cell-cell adhesion",
    "negative regulation of T cell activation", "negative regulation of natural killer cell mediated cytotoxicity",
    "negative regulation of T cell differentiation", "negative regulation of leukocyte proliferation",
    "DNA replication", "cell cycle phase transition", "positive regulation of cell cycle process",
    "tetrahydrofolate metabolic process", "folic acid-containing compound metabolic process",
    "one-carbon metabolic process", "ATP-dependent activity", "DNA helicase activity"
  )

  plot_data <- compare_go %>% filter(Description %in% ora_go_terms_to_plot)
  plot_data@compareClusterResult <- plot_data@compareClusterResult %>%
    mutate(GeneRatio = enrichplot:::parse_ratio(GeneRatio)) %>%
    group_by(Cluster) %>%
    arrange(GeneRatio, .by_group = T)

  if (nrow(as.data.frame(plot_data)) > 0) {
    p_go <- dotplot(plot_data, font.size = 8, showCategory = 20, label_format = 40) +
      scale_x_discrete(labels = c("Up", "Down")) +
      scale_fill_distiller(palette = "YlOrRd", limits = c(NA, 0.05), breaks = c(0.01, 0.03, 0.05)) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) +
      labs(title = "GO Enrichment Analysis (ORA)<br><i>Ikzf2<sup>fl/fl</sup>;Lck-cre</i> vs <i>Ikzf2<sup>fl/fl</sup></i>", x = "") +
      my_theme

    ggsave(plot_path, p_go, units = "cm", width = 8, height = 12)
    cat("ORA GO气泡图已保存到:", plot_path, "\n")
  } else {
    cat("在富集结果中未找到指定的GO terms，无法绘图。\n")
  }
} else {
  cat("ORA GO富集结果文件不存在，无法绘图。\n")
}
