#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2025-06-25
# 描述：存放整合分析中共享的ggplot2主题和颜色配置

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  library(ggtext)
})

# 设置全局字体大小为8pt
my_theme <- theme_cowplot(font_size = 8) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = NA),
    plot.title = element_markdown(hjust = 0.5, size = 8),
    axis.title = element_markdown()
  )


# 颜色设置
cols <- brewer.pal(8, "Dark2")
