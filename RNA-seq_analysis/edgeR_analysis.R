#!/usr/bin/env Rscript
# 作者：Claude AI Assistant
# 创建日期：2024-06-24
# 分析内容：使用edgeR对Cre.act和Flox.act组样本进行配对比较差异分析

# 加载所需的包
library(edgeR)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggrepel)  # 用于在图上添加文本标签
library(RColorBrewer) # 用于配色
library(tibble)

# 设置工作目录
# setwd("RNA-seq_analysis")  # 如需要，请取消注释

# 定义自定义主题，设置全局字体大小为8pt
my_theme <- theme_bw() +
  theme(
    text = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "lines"),
    strip.text = element_text(size = 8)
  )

# 读取数据
counts <- read.csv("readcount.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("coldata.csv", row.names = 1)

# 选择.act样本
act_samples <- grep("\\.act", colnames(counts), value = TRUE)
counts_act <- counts[, act_samples]

# 创建分组和配对信息
group_info <- coldata[act_samples, "type"]
names(group_info) <- act_samples

# 创建配对信息，从sample_rename列获取对应的配对编号(act1, act2, act3)
pair_info <- sub(".*\\.act", "", coldata[act_samples, "sample_rename"])
names(pair_info) <- act_samples

# 获取样本重命名信息，用于图表标记
sample_rename <- coldata[act_samples, "sample_rename"]
names(sample_rename) <- act_samples

# 打印分组和配对信息以便验证
cat("Sample grouping information:\n")
print(data.frame(Sample = names(group_info), Group = group_info, Pair = pair_info, Rename = sample_rename))

# 创建DGEList对象
dge <- DGEList(counts = counts_act, group = group_info)

# 过滤低表达基因（每组中至少有3个样本CPM>1）
keep <- filterByExpr(dge, group = group_info)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes before filtering: ", nrow(dge), "\n")
cat("Genes after filtering: ", nrow(dge_filtered), "\n")

# 标准化
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

# 修改：调整设计矩阵，使logFC大于0表示在Cre中上调
# 方法1：修改对比方向
# 设置参考水平为Flox（即将Flox作为基线）
group_info <- factor(group_info, levels = c("Flox", "Cre"))

# 设计矩阵 - 使用配对设计
# 使用~pair+group模型，考虑配对效应
design <- model.matrix(~pair_info+group_info)
colnames(design) <- c("Intercept", "Pair2", "Pair3", "GroupCre")

# 估计离散度
dge_filtered <- estimateDisp(dge_filtered, design)

# 拟合模型
fit <- glmQLFit(dge_filtered, design)

# 差异表达分析 - 比较GroupCre（即Cre对Flox的差异）
# 在这种设计下，正的logFC表示在Cre中上调
qlf <- glmQLFTest(fit, coef = "GroupCre")
results <- topTags(qlf, n = Inf)
deg <- as.data.frame(results)

# 添加上调/下调标签
# 现在logFC>0表示在Cre中上调，logFC<0表示在Cre中下调
deg$regulation <- ifelse(deg$logFC > 0 & deg$PValue < 0.05, "Up_in_Cre",
                        ifelse(deg$logFC < 0 & deg$PValue < 0.05, "Down_in_Cre", "Not_Sig"))

# 提取标准化的CPM值，添加到差异基因结果中
logcpm <- cpm(dge_filtered, log = TRUE)
sample_expr <- logcpm[rownames(deg), ]
deg <- cbind(deg, sample_expr)

# 保存结果
write.csv(deg, "edgeR_paired_Cre.act_vs_Flox.act_results.csv")

# 输出显著差异基因数量
up_in_cre <- sum(deg$regulation == "Up_in_Cre")
down_in_cre <- sum(deg$regulation == "Down_in_Cre")
cat("Up-regulated genes in Cre: ", up_in_cre, "\n")
cat("Down-regulated genes in Cre: ", down_in_cre, "\n")
cat("Total significant DEGs (PValue < 0.05): ", up_in_cre + down_in_cre, "\n")

# 定义RColorBrewer颜色方案
up_color <- brewer.pal(9, "Set1")[1]    # 红色
down_color <- brewer.pal(9, "Set1")[2]  # 蓝色
notsig_color <- "grey70"
color_palette <- brewer.pal(8, "Dark2")

# 可视化

# 要标记的基因列表
genes_to_label <- c(
    "Mthfr", "Shmt1", "Shmt2", "Mthfd1", "Mthfd2", "Mthfd1l", "Mthfd2l",
    "Dhfr", "Mtr", "Mtrr", "Tyms", "Aldh1l1", "Aldh1l2", "Gart", "Amt",
    "Gcsh", "Gldc", "Phgdh", "Psat1", "Psph", "Adsl", "Atic", "Mat2a",
    "Ahcy", "Ikzf2", "Brd4", "Gzmb", "Prf1", "Mki67", "Bcl2"
)

# 1. 火山图 (使用ggplot2)(重点)
# 创建火山图数据
volcano_data <- deg
volcano_data$Gene <- rownames(volcano_data)
volcano_data$Label <- ifelse(volcano_data$Gene %in% genes_to_label, volcano_data$Gene, "")
volcano_data$ToLabel <- volcano_data$Gene %in% genes_to_label

# 创建火山图
p2 <- ggplot(volcano_data, aes(x = logFC, y = -log10(PValue), color = regulation)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Up_in_Cre" = up_color, "Down_in_Cre" = down_color, "Not_Sig" = notsig_color)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Volcano Plot: Paired Comparison Cre.act vs Flox.act",
       x = "Expression Difference (logFC, Cre vs Flox)",
       y = "-log10(PValue)") +
  my_theme

# 添加基因标签 (使用斜体)
p2 <- p2 + geom_text_repel(
  data = subset(volcano_data, ToLabel == TRUE),
  aes(label = Label),
  fontface = "italic", # 使基因名为斜体
  size = 2.5,          # 更小的字体，约8pt
  box.padding = 0.35,
  point.padding = 0.2,
  force = 2,
  max.overlaps = 30,
  segment.color = "grey50",
  segment.size = 0.2
)
p2
ggsave("Volcano_plot_paired_Cre.act_vs_Flox.act.pdf", p2, width = 5, height = 4)

# 2. 热图（只展示感兴趣的基因）
# 从感兴趣的基因中过滤出存在于我们数据集中的基因
interest_genes_in_data <- genes_to_label[genes_to_label %in% rownames(logcpm)]

if(length(interest_genes_in_data) > 0) {
  # 提取感兴趣基因的表达值
  interest_expr <- logcpm[interest_genes_in_data, ]

  # 对数据进行缩放，便于可视化
  scaled_expr <- t(scale(t(interest_expr)))

  # 准备样本注释
  sample_anno <- data.frame(
    Group = group_info,
    Pair = pair_info,
    row.names = names(group_info)
  )

  # 自定义颜色
  anno_colors <- list(
    Group = c(Cre = brewer.pal(8, "Set2")[1], Flox = brewer.pal(8, "Set2")[2]),
    Pair = c("1" = brewer.pal(8, "Dark2")[1], "2" = brewer.pal(8, "Dark2")[2], "3" = brewer.pal(8, "Dark2")[3])
  )

  # 绘制热图 - 使用sample_rename作为列名
  # 创建带有重命名的表达矩阵
  scaled_expr_renamed <- scaled_expr
  colnames(scaled_expr_renamed) <- sample_rename[colnames(scaled_expr)]

  # 更新sample_anno的行名以匹配新的列名
  sample_anno_renamed <- sample_anno
  rownames(sample_anno_renamed) <- sample_rename[rownames(sample_anno)]

  # 绘制热图
  pheatmap(scaled_expr_renamed,
           main = "Expression of Genes of Interest (PValue < 0.05)",
           cluster_rows = TRUE,
           cluster_cols = FALSE,    # 不对样品进行聚类
           show_rownames = TRUE,
           fontsize = 8,         # 字体大小8pt
           fontsize_row = 8,     # 行名字体大小
           border_color = "white",
           cellwidth = 8,
           cellheight = 8,
           annotation_col = sample_anno_renamed,
           annotation_colors = anno_colors,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           filename = "Heatmap_genes_of_interest.pdf",
           width = 5,
           height = 6)

  # 创建ggplot版本的热图数据
  heatmap_data <- as.data.frame(scaled_expr) %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

  # 添加重命名的样本信息
  heatmap_data$SampleName <- sample_rename[heatmap_data$Sample]
  heatmap_data$Group <- group_info[heatmap_data$Sample]
  heatmap_data$Pair <- pair_info[heatmap_data$Sample]

  # 根据层次聚类确定基因顺序
  hc_genes <- hclust(dist(scaled_expr))
  gene_order <- rownames(scaled_expr)[hc_genes$order]
  heatmap_data$Gene <- factor(heatmap_data$Gene, levels = gene_order)

  # 保持样本原始顺序，不进行聚类
  # 根据Group排序
  group_order <- names(sort(table(group_info), decreasing = TRUE))
  sample_by_group <- lapply(group_order, function(g) names(group_info[group_info == g]))
  ordered_samples <- unlist(sample_by_group)
  renamed_sample_order <- sample_rename[ordered_samples]
  heatmap_data$SampleName <- factor(heatmap_data$SampleName, levels = renamed_sample_order)

  # 使用ggplot2创建热图
  p3 <- ggplot(heatmap_data, aes(x = SampleName, y = Gene, fill = Expression)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
    facet_grid(. ~ Group, scales = "free", space = "free") +
    labs(title = "Expression of Genes of Interest (PValue < 0.05)",
         x = "", y = "") +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8, face = "italic"), # 基因名斜体
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      legend.key.size = unit(0.8, "lines"),
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(size = 8),
      panel.grid = element_blank()
    )
  p3
  ggsave("Heatmap_genes_of_interest_ggplot.pdf", p3, width = 3, height = 5)
}

# 3. PCA分析 (使用sample_rename作为标签)
logcpm_all <- cpm(dge_filtered, log = TRUE)
pca <- prcomp(t(logcpm_all), scale = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])
pca_data$group <- group_info
pca_data$pair <- pair_info
pca_data$sample <- rownames(pca_data)
pca_data$sample_rename <- sample_rename[pca_data$sample]

# 使用RColorBrewer配色
p4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = pair, label = sample_rename)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, box.padding = 0.3, point.padding = 0.2, max.overlaps = 20) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PCA: Paired Comparison Cre.act vs Flox.act",
       x = paste0("PC1: ", round(summary(pca)$importance[2,1]*100, 1), "%"),
       y = paste0("PC2: ", round(summary(pca)$importance[2,2]*100, 1), "%")) +
  my_theme
p4
ggsave("PCA_plot_paired.pdf", p4, width = 4, height = 3)

# 4. 感兴趣基因的箱线图 (使用sample_rename)
# 从指定的感兴趣基因中选择存在于数据集中的
interest_genes <- genes_to_label[genes_to_label %in% rownames(logcpm)]

if(length(interest_genes) > 0) {
  # 准备箱线图数据
  boxplot_data <- data.frame()
  
  for(gene in interest_genes) {
    gene_expr <- logcpm[gene, ]
    temp_data <- data.frame(
      Gene = gene,
      Expression = gene_expr,
      Sample = names(gene_expr),
      SampleName = sample_rename[names(gene_expr)],
      Group = group_info[names(gene_expr)],
      Pair = pair_info[names(gene_expr)]
    )
    boxplot_data <- rbind(boxplot_data, temp_data)
  }
  print(head(boxplot_data))
  
  # 将Gene转换为因子并设置为斜体
  boxplot_data$Gene <- factor(boxplot_data$Gene)
  
  # 确保Gene列有值
  if(length(unique(boxplot_data$Gene)) == 0) {
    cat("警告: 找不到任何匹配的基因数据\n")
  } else {
    # 使用所有感兴趣的基因而不是限制为前12个
    boxplot_data_subset <- boxplot_data
    
    # 检查数据
    cat("基因数量:", length(unique(boxplot_data_subset$Gene)), "\n")
    cat("样本数量:", length(unique(boxplot_data_subset$Sample)), "\n")
    cat("组别数量:", length(unique(boxplot_data_subset$Group)), "\n")
    
    # 确保Gene列有效
    if(length(unique(boxplot_data_subset$Gene)) > 0) {
      # 创建箱线图，根据基因数量调整列数
      gene_count <- length(unique(boxplot_data_subset$Gene))
      # 调整列数，确保每行不超过5个基因
      ncol_value <- min(5, gene_count)
      
      # 添加显著性信息
      sig_data <- data.frame()
      for(gene in unique(boxplot_data_subset$Gene)) {
        # 获取差异分析结果中的P值
        if(gene %in% rownames(deg)) {
          p_value <- deg[gene, "PValue"]
          sig_level <- ifelse(p_value < 0.001, "***", 
                             ifelse(p_value < 0.01, "**", 
                                   ifelse(p_value < 0.05, "*", "ns")))
          
          # 计算显著性标记的位置
          gene_data <- boxplot_data_subset[boxplot_data_subset$Gene == gene, ]
          y_max <- max(gene_data$Expression, na.rm = TRUE)
          y_pos <- y_max + 0.1 * (max(gene_data$Expression, na.rm = TRUE) - min(gene_data$Expression, na.rm = TRUE))
          
          temp_sig <- data.frame(
            Gene = gene,
            y = y_pos,
            label = sig_level,
            p_value = p_value
          )
          sig_data <- rbind(sig_data, temp_sig)
        }
      }
      
      # 创建箱线图
      p6 <- ggplot(boxplot_data_subset, aes(x = Group, y = Expression, fill = Group)) +
        geom_boxplot(width = 0.5, outlier.shape = NA) +
        geom_point(aes(color = Pair), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), size = 1) +
        scale_fill_brewer(palette = "Pastel1") +
        scale_color_brewer(palette = "Dark2") +
        facet_wrap(~ Gene, scales = "free_y", ncol = ncol_value) +
        labs(title = "Expression of Genes of Interest",
             x = "", y = "Expression (logCPM)") +
        theme_bw(base_size = 8) +
        theme(
          strip.text = element_text(face = "italic", size = 8), # 基因名斜体
          legend.position = "bottom",
          legend.box = "vertical",
          legend.key.size = unit(0.8, "lines"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      
      # 添加显著性标记
      if(nrow(sig_data) > 0) {
        p6 <- p6 + geom_text(
          data = sig_data,
          aes(x = 1.5, y = y, label = label),
          size = 2.5,
          inherit.aes = FALSE
        )
      }
      
      print(p6)
      
      # 自动调整保存的图像大小，基于基因数量
      width_value <- min(12, max(6, ncol_value * 1.5))
      height_value <- min(15, max(8, ceiling(gene_count/ncol_value) * 1.5))
      
      ggsave("Boxplot_interest_genes.pdf", p6, width = width_value, height = height_value)
      
      # 创建配对线图 - 使用sample_rename作为标签
      # 重整数据为配对格式
      paired_data <- data.frame()
      
      for(gene in unique(boxplot_data_subset$Gene)) {
        for(pair_id in unique(boxplot_data_subset$Pair)) {
          gene_subset <- boxplot_data_subset[boxplot_data_subset$Gene == gene & boxplot_data_subset$Pair == pair_id, ]
          
          if(nrow(gene_subset) == 2) { # 确保有配对
            paired_data <- rbind(paired_data, gene_subset)
          }
        }
      }
      
      # 检查配对数据是否有效
      if(nrow(paired_data) > 0) {
        # 点和线显示配对样本
        p7 <- ggplot(paired_data, aes(x = Group, y = Expression, group = interaction(Gene, Pair), color = Pair)) +
          geom_point(size = 1.5) +
          geom_line(linewidth = 0.5) +
          scale_color_brewer(palette = "Dark2") +
          facet_wrap(~ Gene, scales = "free_y", ncol = ncol_value) +
          labs(title = "Paired Expression of Genes of Interest",
               x = "", y = "Expression (logCPM)") +
          theme_bw(base_size = 8) +
          theme(
            strip.text = element_text(face = "italic", size = 8), # 基因名斜体
            legend.position = "bottom",
            legend.key.size = unit(0.8, "lines"),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
        
        # 添加显著性标记到配对图
        if(nrow(sig_data) > 0) {
          p7 <- p7 + geom_text(
            data = sig_data,
            aes(x = 1.5, y = y, label = label),
            size = 3,
            inherit.aes = FALSE
          )
        }
        
        print(p7)
        ggsave("Paired_lines_interest_genes.pdf", p7, width = width_value, height = height_value)
      } else {
        cat("警告: 无法创建配对线图，没有足够的配对数据\n")
      }
    } else {
      cat("警告: 没有有效的Gene列数据进行分面\n")
    }
  }
} else {
  cat("警告: 指定的基因都不在表达数据中\n")
}

# 5. 额外添加：探索我们感兴趣的基因在差异基因中的位置
# 创建函数来检查基因是否在差异基因中
check_genes_interest <- function(gene_list, deg_df) {
  result <- data.frame(
    Gene = gene_list,
    In_DEGs = gene_list %in% rownames(deg_df),
    logFC = NA,
    PValue = NA,
    Regulation = NA,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(result)) {
    if (result$In_DEGs[i]) {
      gene_idx <- which(rownames(deg_df) == result$Gene[i])
      result$logFC[i] <- deg_df$logFC[gene_idx]
      result$PValue[i] <- deg_df$PValue[gene_idx]
      result$Regulation[i] <- deg_df$regulation[gene_idx]
    }
  }

  result <- result[order(result$PValue), ]
  return(result)
}

# 检查我们感兴趣的基因
interest_genes_result <- check_genes_interest(genes_to_label, deg)
write.csv(interest_genes_result, "One_carbon_metabolism_genes_results.csv")

# 输出这些基因的情况
cat("\nResults for genes of interest:\n")
print(interest_genes_result)

cat("Paired comparison analysis completed! Results saved to current directory.\n")
cat("Note: Significant genes were filtered using only PValue < 0.05 criterion.\n")
cat("logFC direction has been adjusted: now positive logFC means up-regulated in Cre mice.\n")

