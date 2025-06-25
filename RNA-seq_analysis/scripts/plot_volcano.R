#!/usr/bin/env Rscript
# 脚本目的: 绘制火山图

# 1. 加载包和函数--------------------------------
suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(RColorBrewer)
  library(cowplot)
  library(ggtext)
})
source("functions/config_and_themes.R")

# 2. 定义路径--------------------------------
de_results_path <- "results/tables/DE_results.csv"
plot_path <- "results/plots/volcano_plot.pdf"

# 3. 读取数据--------------------------------
if (!file.exists(de_results_path)) stop("DE results not found. Run 01_run_DE_analysis.R first.")
deg <- read.csv(de_results_path, row.names = 1)

# 4. 准备绘图数据--------------------------------
plot_data <- deg %>%
  tibble::rownames_to_column("Gene") %>%
  mutate(
    Label = if_else(Gene %in% genes_to_label, Gene, ""),
    ToLabel = Gene %in% genes_to_label,
    regulation = factor(regulation, levels = c("Up_in_Cre", "Down_in_Cre", "Not_Sig"))
  )
print(head(plot_data))

# 计算差异基因数量
up_in_cre <- sum(plot_data$regulation == "Up_in_Cre")
down_in_cre <- sum(plot_data$regulation == "Down_in_Cre")

# 5. 绘制火山图--------------------------------
p_volcano <- ggplot(plot_data, aes(x = logFC, y = -log10(PValue), color = regulation)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(
    name = "Regulation",
    values = c("Up_in_Cre" = up_color, "Down_in_Cre" = down_color, "Not_Sig" = notsig_color),
    labels = c(paste0("Up (", up_in_cre, ")"), paste0("Down (", down_in_cre, ")"), "Not Sig")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Paired Comparison<br><i>Ikzf2<sup>fl/fl</sup>;Lck-cre</i> vs <i>Ikzf2<sup>fl/fl</sup></i>",
    x = "log~2~Fold Change",
    y = "-log~10~(*P*-value)"
  ) +
  my_theme +
  geom_text_repel(
    data = subset(plot_data, ToLabel == TRUE),
    aes(label = Label),
    color = "black", box.padding = 0.5,
    min.segment.length = 0, size = 2.5, fontface = "italic"
  )
p_volcano
ggsave(plot_path, p_volcano, units = "cm", width = 7, height = 6)

cat("火山图已保存到:", plot_path, "\n")
