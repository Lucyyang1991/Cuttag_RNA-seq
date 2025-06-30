# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
})

# 1. 全局绘图主题
my_theme <- theme_cowplot(font_size = 8) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = NA),
    plot.title = element_markdown(hjust = 0.5, size = 8),
    axis.title = element_markdown()
)

# 2. 颜色方案
up_color <- brewer.pal(9, "Set1")[1]    # 红色
down_color <- brewer.pal(9, "Set1")[2]  # 蓝色
notsig_color <- "grey70"

# 3. 标记基因列表
genes_to_label <- c("Shmt1", "Mthfd1", "Ikzf2")
