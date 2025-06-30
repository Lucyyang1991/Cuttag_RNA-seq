#!/usr/bin/env Rscript
# 作者：Claude AI Assistant & Your Name
# 创建日期：2025-06-25
# 功能：使用ComplexHeatmap绘制特定基因集的热图

# 1. 加载必要的包 -------------------------------------------------------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(circlize)
    library(gridtext)
})

# 2. 定义路径和参数 ---------------------------------------------------------------
setwd("D:/2023_Git/cut&tag/20230308结果")
kegg_data_path <- "./Integrative_analysis/results/enrichment/KEGG_enrichment.RData"
rnaseq_data_path <- "./RNA-seq_analysis/results/RData/DE_objects.RData"
output_dir <- "./Integrative_analysis/results/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pathway_name <- "One carbon pool by folate"
font_size <- 8

# 3. 加载数据 ----------------------------------------------------------------------
cat("INFO: 加载KEGG富集结果...\n")
if (!file.exists(kegg_data_path)) {
    stop("错误: 找不到文件: ", kegg_data_path)
}
load(kegg_data_path, verbose = FALSE) # 加载 kegg_readable 对象

cat("INFO: 加载RNA-seq结果...\n")
if (!file.exists(rnaseq_data_path)) {
    stop("错误: 找不到文件: ", rnaseq_data_path)
}
load(rnaseq_data_path, verbose = FALSE) # 加载 deg, plot_meta_data, act_samples 对象

# 4. 提取并准备基因列表 -----------------------------------------------------------------
# 从KEGG结果中提取基因
one_carbon_genes <- kegg_readable %>%
    as.data.frame() %>%
    filter(grepl(pathway_name, Description, ignore.case = TRUE)) %>%
    pull(geneID) %>%
    str_split("/") %>%
    unlist() %>%
    unique()

if (length(one_carbon_genes) == 0) {
    warning("警告: 未在KEGG结果中找到 '", pathway_name, "' 通路。")
} else {
    cat("INFO: 成功提取 '",
        pathway_name,
        "' 通路的基因:",
        length(one_carbon_genes),
        "个\n")
}

# 定义增殖相关基因
proliferation_markers <- c("Mki67",
                           "Ccnb1",
                           "Cdk1",
                           "Plk1",
                           "Aurkb",
                           "Cdc20",
                           "Bub1b",
                           "Birc5",
                           "Mad2l1",
                           "Cdt1")
cat("INFO: 定义了经典增殖标记基因:", length(proliferation_markers), "个\n")

# 5. 准备用于热图的表达矩阵 -------------------------------------------------------
# 按Group（Cre在前）和Pair排序样本元数据
ordered_meta_data <- plot_meta_data %>%
    arrange(factor(Group, levels = c("Cre", "Flox")), Pair)

# 提取排序后的原始样本名，这对应deg矩阵的列名
ordered_samples <- rownames(ordered_meta_data)

# 确定所有需要展示的基因，且这些基因存在于表达矩阵中
all_available_genes <- intersect(c(one_carbon_genes, proliferation_markers), rownames(deg))

# 从deg数据框中，按排好的样本顺序提取表达量列，并筛选需要的基因
heatmap_matrix <- deg[all_available_genes, ordered_samples]

# 对数据进行Z-score标准化（按行）
scaled_matrix <- heatmap_matrix %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    t() %>%
    na.omit() # na.omit处理由scale产生的NaN（当行中所有值相同时）

# 按照基因名（行名）字母顺序排序
scaled_matrix <- scaled_matrix[sort(rownames(scaled_matrix)), ]

# 将列名更新为配对信息
colnames(scaled_matrix) <- paste0("Pair", ordered_meta_data$Pair)

cat("INFO: 表达矩阵准备和标准化完成。\n")

# 6. 创建热图注释 -----------------------------------------------------------------
# 创建列注释
col_annotation <- HeatmapAnnotation(
    Group = ordered_meta_data$Group,
    Pair = ordered_meta_data$Pair,
    col = list(
        Group = c("Flox" = "#1B9E77", "Cre" = "#D95F02"),
        Pair = c(
            "1" = "#7570B3",
            "2" = "#E7298A",
            "3" = "#66A61E"
        )
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = font_size),
    simple_anno_size = unit(0.3, "cm"),
    annotation_legend_param = list(
        Group = list(
            title = "Group",
            labels = c(
                expression(italic(Ikzf2)^{fl/fl}*";Lck-cre"), expression(italic(Ikzf2)^{fl/fl})
            ),
            title_gp = gpar(fontsize = font_size),
            labels_gp = gpar(fontsize = font_size)
        ),
        Pair = list(
            title = "Pair",
            title_gp = gpar(fontsize = font_size),
            labels_gp = gpar(fontsize = font_size)
        )
    )
)

# 创建行分割（基因集来源）
row_split <- ifelse(
    rownames(scaled_matrix) %in% one_carbon_genes,
    "One carbon pool\nby folate",
    "Proliferation"
)

# 7. 绘制热图 ----------------------------------------------------------------------
cat("INFO: 开始绘制热图...\n")

# 定义颜色映射
spectral_pal <- brewer.pal(11, "Spectral")
col_fun <- colorRamp2(seq(-2, 2, length.out = 11), rev(spectral_pal))

# 绘制热图
ht <- Heatmap(
    scaled_matrix,
    name = "Z-score",
    col = col_fun,

    # 列设置
    top_annotation = col_annotation,
    cluster_columns = FALSE,
    column_split = ordered_meta_data$Group,
    column_title = NULL,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = font_size),

    # 行设置
    row_split = row_split,
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = font_size),
    row_names_gp = gpar(fontsize = font_size, fontface = "italic"),
    cluster_rows = FALSE,
    show_row_dend = TRUE,

    # 图例设置
    heatmap_legend_param = list(
        title = "Z-score",
        title_gp = gpar(fontsize = font_size),
        labels_gp = gpar(fontsize = font_size),
        legend_height = unit(1.5, "cm")
    )
)

# 8. 保存热图 ----------------------------------------------------------------------
output_pdf_path <- file.path(output_dir, "heatmap_one_carbon_and_proliferation.pdf")
pdf(output_pdf_path,
    width = 8 / 2.54,
    height = 8 / 2.54)
draw(ht,
     merge_legend = TRUE,
     heatmap_legend_side = "right")
dev.off()

cat("INFO: 热图已保存到:", output_pdf_path, "\n")
cat("INFO: 脚本 04_plot_heatmap.R 执行完毕。\n")
