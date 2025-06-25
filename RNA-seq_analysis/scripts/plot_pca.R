#!/usr/bin/env Rscript
# 脚本目的: 绘制PCA图

# 1. 加载包和函数--------------------------------
suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(RColorBrewer)
  library(cowplot)
  library(edgeR) # for cpm()
  library(ggtext)
})
source("functions/config_and_themes.R")

# 2. 定义路径--------------------------------
de_objects_path <- "results/RData/DE_objects.RData"
plot_path <- "results/plots/PCA_plot.pdf"

# 3. 读取数据--------------------------------
if (!file.exists(de_objects_path)) stop("DE objects not found. Run 01_run_DE_analysis.R first.")
load(de_objects_path)

# 4. 准备绘图数据--------------------------------
logcpm_all <- cpm(dge_filtered, log = TRUE)
pca <- prcomp(t(logcpm_all), scale = TRUE)

plot_data <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(
    Group = plot_meta_data$Group,
    Pair = plot_meta_data$Pair,
    Label = plot_meta_data$sample_rename
  )

# 5. 绘制PCA图--------------------------------
p_pca <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group, shape = Pair, label = Label)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, box.padding = 0.5, min.segment.length = 0) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "PCA of Paired Samples",
    x = paste0("PC1: ", round(summary(pca)$importance[2,1]*100, 1), "%"),
    y = paste0("PC2: ", round(summary(pca)$importance[2,2]*100, 1), "%")
  ) +
  my_theme
p_pca
ggsave(plot_path, p_pca, width = 10, height = 8, units = "cm")

cat("PCA图已保存到:", plot_path, "\n")
