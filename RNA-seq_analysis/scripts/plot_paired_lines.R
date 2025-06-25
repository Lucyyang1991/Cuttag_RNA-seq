#!/usr/bin/env Rscript
# 脚本目的: 为感兴趣的基因绘制配对线图

# 1. 加载包和函数--------------------------------
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RColorBrewer)
  library(cowplot)
  library(ggtext)
})
source("functions/config_and_themes.R")

# 2. 定义路径--------------------------------
de_results_path <- "results/tables/DE_results.csv"
de_objects_path <- "results/RData/DE_objects.RData"
plot_path <- "results/plots/paired_lines_interest_genes.pdf"

# 3. 读取数据--------------------------------
if (!file.exists(de_results_path)) stop("DE results not found. Run 01_run_DE_analysis.R first.")
if (!file.exists(de_objects_path)) stop("DE objects not found. Run 01_run_DE_analysis.R first.")
deg <- read.csv(de_results_path, row.names = 1)
load(de_objects_path)

# 4. 准备绘图数据--------------------------------
if (length(genes_to_label) > 0) {
  plot_data <- deg %>%
    tibble::rownames_to_column("Gene") %>%
    filter(Gene %in% genes_to_label) %>%
    dplyr::select(Gene, PValue, all_of(act_samples)) %>%
    pivot_longer(cols = all_of(act_samples), names_to = "Sample", values_to = "Expression") %>%
    left_join(
      plot_meta_data %>% rownames_to_column("Sample"),
      by = "Sample"
    )

  # 准备显著性标记
  sig_data <- plot_data %>%
    group_by(Gene) %>%
    summarise(
      y_pos = max(Expression, na.rm = TRUE) * 1.05,
      PValue = first(PValue),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        PValue < 0.001 ~ "***",
        PValue < 0.01  ~ "**",
        PValue < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    )

  # 5. 绘制配对线图--------------------------------
  if (nrow(plot_data) > 0) {
    gene_count <- length(unique(plot_data$Gene))
    ncol_value <- min(5, gene_count)
    p_paired <- ggplot(plot_data, aes(x = Group, y = Expression, group = Pair, color = Pair)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.5) +
      scale_color_brewer(palette = "Dark2") +
      scale_x_discrete(labels = c("Flox" = "*Ikzf2*<sup>fl/fl</sup>", "Cre" = "*Ikzf2*<sup>fl/fl</sup>;*Lck*-cre")) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
      facet_wrap(~ Gene, ncol = ncol_value) +
      labs(x = "", y = "Expression (log~2~CPM)") +
      geom_text(
        data = sig_data,
        aes(x = 1.5, y = y_pos, label = label),
        inherit.aes = FALSE, size = 3
      ) +
      my_theme +
      theme(
        strip.text = element_text(face = "italic", size = 8),
        legend.position = "bottom",
        axis.text.x = element_markdown(angle = 60, hjust = 1)
      )

    ggsave(plot_path, p_paired, width = 10, height = 6, units = "cm")
    cat("配对线图已保存到:", plot_path, "\n")
  } else {
    cat("警告: 没有足够的配对数据用于创建线图\n")
  }
} else {
  cat("警告: 未在配置文件中指定任何感兴趣的基因\n")
}
