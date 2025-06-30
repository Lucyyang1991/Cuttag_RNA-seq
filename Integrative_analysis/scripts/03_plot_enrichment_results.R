#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2025-06-25
# 功能：可视化GO和KEGG富集分析的结果

# 1. 加载包 ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(enrichplot)
  library(ggbreak)
  library(Cairo)
  library(cowplot)
  library(ggtext)
  library(ggplot2)
})

# 2. 加载配置和定义路径 --------------------------------------------------------
setwd("D:/2023_Git/cut&tag/20230308结果")
source("./Integrative_analysis/functions/config_and_themes.R")

enrich_dir <- "./Integrative_analysis/results/enrichment"
plot_dir <- "./Integrative_analysis/results/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

go_rdata <- file.path(enrich_dir, "GO_enrichment.RData")
kegg_rdata <- file.path(enrich_dir, "KEGG_enrichment.RData")

# 3. 定义绘图函数 --------------------------------------------------------------

#' @description 绘制富集分析条形图
plot_enrich_bar <- function(data, title, file_path) {
  p <- ggplot(data, aes(x = log_padj, y = Description)) +
    geom_bar(stat = "identity", width = 0.6, fill = "steelblue") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = title, x = "-log<sub>10</sub><i>p.adjust</i>", y = NULL) +
    scale_y_discrete(labels = scales::label_wrap(width = 40)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    my_theme +
    theme(
      axis.text.y = element_text(size = 8, hjust = 1)
    )

  if (max(data$log_padj, na.rm = TRUE) > 5) {
    p <- p + scale_x_cut(breaks = c(4), which = 2)
  }

  cairo_pdf(file_path, width = 8 / 2.54, height = max(2, nrow(data) * 0.4 / 2.54))
  print(p)
  dev.off()
}


# 4. GO 可视化 ----------------------------------------------------------------
cat("INFO: 开始处理GO结果可视化...\n")
if (file.exists(go_rdata)) {
  load(go_rdata) # 加载 go_all 对象

  if (!is.null(go_all) && nrow(go_all@result) > 0) {
    # 筛选和准备数据
    go_plot_data <- go_all@result %>%
      arrange(p.adjust) %>%
      slice(unique(c(1:min(2, nrow(.)), grep("carbon", .$Description, ignore.case = TRUE)[c(1, 3)]))) %>%
      mutate(
        log_padj = -log10(p.adjust),
        Description = factor(Description, levels = rev(Description))
      )
    # 绘图
    plot_enrich_bar(go_plot_data, "GO Enrichment Analysis", file.path(plot_dir, "GO_barplot.pdf"))
    cat("INFO: GO条形图已保存。\n")

  } else {
    cat("INFO: GO富集结果为空，跳过可视化。\n")
  }
} else {
  cat("WARNING: 未找到GO_enrichment.RData文件，跳过GO可视化。\n")
}

# 5. KEGG 可视化 ---------------------------------------------------------------
cat("INFO: 开始处理KEGG结果可视化...\n")
if (file.exists(kegg_rdata)) {
  load(kegg_rdata) # 加载 kegg_readable 对象

  if (!is.null(kegg_readable) && nrow(kegg_readable@result) > 0) {
    # 筛选和准备数据
    kegg_plot_data <- kegg_readable@result %>%
      mutate(Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description)) %>%
      arrange(p.adjust) %>%
      slice(unique(c(1:min(2, nrow(.)), grep("carbon", .$Description, ignore.case = TRUE)[1]))) %>%
      filter(!is.na(ID)) %>%
      mutate(
        log_padj = -log10(p.adjust),
        Description = factor(Description, levels = rev(Description))
      )

    # 绘图
    plot_enrich_bar(kegg_plot_data, "KEGG Pathway Enrichment", file.path(plot_dir, "KEGG_barplot.pdf"))
    cat("INFO: KEGG条形图已保存。\n")

  } else {
    cat("INFO: KEGG富集结果为空，跳过可视化。\n")
  }
} else {
  cat("WARNING: 未找到KEGG_enrichment.RData文件，跳过KEGG可视化。\n")
}

cat("INFO: 脚本 03_plot_enrichment_results.R 执行完毕。\n")
