#!/usr/bin/env Rscript
# 脚本目的:
# 1. 使用edgeR对配对样本进行差异表达分析。
# 2. 保存差异表达结果表格。
# 3. 保存后续分析和绘图所需的核心对象。

# 1. 加载包--------------------------------
suppressMessages({
  library(edgeR)
  library(dplyr)
  library(tibble)
})

# 2. 定义输入输出路径--------------------------------
counts_path <- "readcount.csv"
coldata_path <- "coldata.csv"
de_results_path <- "results/tables/DE_results.csv"
de_objects_path <- "results/RData/DE_objects.RData"

# 3. 创建结果目录--------------------------------
dir.create(dirname(de_results_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(de_objects_path), recursive = TRUE, showWarnings = FALSE)

# 4. 读取并准备数据--------------------------------
counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_path, row.names = 1)

# 筛选.act样本并准备元数据
act_samples <- grep("\\.act", colnames(counts), value = TRUE)
counts_act <- counts[, act_samples]

coldata_act <- coldata[act_samples, ] %>%
  mutate(
    Pair = sub(".*\\.act", "", sample_rename),
    Group = factor(type, levels = c("Flox", "Cre")) # 设置"Flox"为对照组
  )

group_info <- setNames(coldata_act$Group, rownames(coldata_act))
pair_info <- setNames(coldata_act$Pair, rownames(coldata_act))

# 5. DGEList对象创建、过滤和标准化--------------------------------
dge <- DGEList(counts = counts_act, group = group_info)
keep <- filterByExpr(dge, group = group_info)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
cat("基因过滤与标准化完成。\n")

# 6. 设计矩阵与差异分析--------------------------------
design <- model.matrix(~pair_info + group_info)
colnames(design) <- c("Intercept", "Pair2", "Pair3", "GroupCre")

dge_filtered <- estimateDisp(dge_filtered, design)
fit <- glmQLFit(dge_filtered, design)
qlf <- glmQLFTest(fit, coef = "GroupCre")
results <- topTags(qlf, n = Inf)

# 7. 整理并保存差异分析结果--------------------------------
logcpm <- cpm(dge_filtered, log = TRUE)
deg <- as.data.frame(results) %>%
  mutate(
    regulation = case_when(
      logFC > 0 & PValue < 0.05 ~ "Up_in_Cre",
      logFC < 0 & PValue < 0.05 ~ "Down_in_Cre",
      TRUE ~ "Not_Sig"
    )
  ) %>%
  merge(logcpm, by = "row.names") %>%
  column_to_rownames("Row.names")
print(head(deg))
write.csv(deg, de_results_path)
cat("差异分析结果已保存到:", de_results_path, "\n")

# 8. 保存核心对象以供下游分析--------------------------------
plot_meta_data <- coldata_act %>%
    select(Group = "type", Pair, sample_rename)

save(
  deg,
  dge_filtered,
  plot_meta_data,
  act_samples,
  file = de_objects_path
)
cat("核心对象已保存到:", de_objects_path, "\n")

cat("差异分析脚本执行完毕！\n")
