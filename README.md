# Cut&Tag测序数据分析流程

## 1. 项目概述

本项目是一个整合了Cut&Tag测序和RNA-seq数据的分析流程，旨在鉴定转录因子结合位点，并结合基因表达变化分析其调控功能。项目包含Cut&Tag数据的完整处理流程（从质控、比对到峰值检测和注释）、RNA-seq的差异表达分析，并最终通过`Integrative_analysis`模块对两者进行整合分析。

## 2. 目录结构

项目的建议目录结构如下：

```
project_root/
├── CutTag_Linux/             # 原始数据及处理脚本 (服务器端运行)
│   ├── scripts/              # 各步骤处理脚本
│   ├── raw_data/             # 原始FASTQ文件存放目录 (示例)
│   ├── genome/               # 参考基因组及注释文件
│   └── ...                   # 其他分析结果目录
├── CutTag_local_analysis/    # CUT&Tag下游本地分析 (本地环境运行)
│   └── *.R                   # R分析脚本 (如峰值注释、筛选)
├── RNA-seq_analysis/         # RNA-seq 数据分析 (本地环境运行)
│   ├── scripts/              # RNA-seq分析脚本
│   ├── results/              # RNA-seq分析结果
│   └── ...
├── Integrative_analysis/     # 整合分析 (本地环境运行)
│   ├── scripts/              # 整合分析脚本
│   └── results/              # 整合分析结果
└── reference_repos/          # 参考分析流程 (学习和参考用途)
```

详细的目录结构规范请参考项目级规则（如 `.cursor/rules`）。

## 3. 分析流程

本项目的分析流程主要分为服务器端分析和本地分析两部分。

### 3.1 服务器端分析流程 (位于 `CutTag_Linux` 目录)

服务器端分析流程主要处理原始测序数据，包括以下步骤：

1.  **数据准备**:
    *   `scripts/config/00_config.sh`: 配置全局变量、样本信息和路径。
    *   `scripts/00_subsample.sh` (可选): 对原始数据进行抽样，用于流程测试。
2.  **数据预处理**:
    *   `scripts/01_fastqc.sh` (可选): 使用 FastQC 对原始数据进行质量评估。
    *   `scripts/01_fastqc_stats.sh`: 汇总FastQC结果并生成统计报告。
    *   `scripts/02_trimming.sh`: 使用 Cutadapt 去除测序接头和低质量序列。
    *   `scripts/02_trimming_stats.sh`: 汇总去接头结果并生成统计报告。
3.  **序列比对**:
    *   `scripts/03_alignment.sh`: 使用 Bowtie2 将处理后的序列比对到参考基因组，并进行BAM文件处理（排序、标记重复、过滤）。
    *   `scripts/03_alignment_stats.sh`: 生成详细的比对质量报告。
4.  **峰值检测 (Peak Calling)**:
    *   `scripts/04_bedgraph_generation.sh`: 生成标准化的 BedGraph 文件，作为Peak Calling的输入。
    *   `scripts/04_peak_calling_macs3.sh`: 使用 MACS3 进行峰值检测（推荐）。
    *   `scripts/04_peak_calling_seacr.sh`: (可选) 使用 SEACR 进行峰值检测。
    *   `scripts/04_idr_analysis.sh`: (可选) 使用 IDR (Irreproducible Discovery Rate) 分析生物学重复样本间Peak的一致性，筛选高置信度的Peak。
    *   `scripts/04_peak_calling_stats.sh`: 生成详细的Peak Calling统计报告。
5.  **下游可视化准备**:
    *   `scripts/05_bam_to_bigwig.sh`: 将比对产生的BAM文件转换为BigWig格式，用于在IGV等基因组浏览器中可视化。
6.  **高级分析 (可选)**:
    *   `scripts/06_generate_heatmap.sh`: 使用 deepTools 生成TSS区域的信号富集热图。
    *   `scripts/07_sample_correlation.sh`: 使用 deepTools 分析样本间的相关性。
7.  **IGV风格可视化**:
    *   `scripts/08_plot_igv_style.sh`: 使用 pyGenomeTracks 生成IGV风格的基因组浏览器图，展示特定基因区域的信号富集模式。

**运行方式**:
服务器端的脚本设计为通过SGE (Sun Grid Engine) 作业调度系统提交和运行。具体提交方式请参考项目级规则。

### 3.2 本地分析流程 (位于 `CutTag_local_analysis`, `RNA-seq_analysis`, `Integrative_analysis` 目录)

本地分析流程主要对服务器端产生的Cut&Tag峰值文件等结果，以及RNA-seq数据进行下游的注释、差异表达分析、功能富集和整合分析。

1.  **`CutTag_local_analysis`**: Cut&Tag峰值注释与筛选
    *   主要使用R包 `ChIPseeker` 对服务器产生的peak文件进行注释，分析其在基因组上的分布特征。
    *   `05_peak_annotation_macs3.R`: 对MACS3产生的peak进行注释。
    *   `05_peak_annotation_idr.R`: 对IDR分析产生的peak进行注释。
    *   `06_filter_promoter_peaks.R`: 筛选启动子区域的peak。
    *   `07_filter_proximal_promoter_genes.R`: 筛选近端启动子区域的峰值及其关联基因。

2.  **`RNA-seq_analysis`**: 差异表达与功能富集分析
    *   `scripts/01_run_DE_analysis.R`: 使用 `edgeR` 等工具进行差异表达分析。
    *   `scripts/02_run_enrichment_analysis.R`: 对差异表达基因进行GO和KEGG等功能富集分析。
    *   包含多种绘图脚本 (`plot_volcano.R`, `plot_pca.R` 等) 对分析结果进行可视化。

3.  **`Integrative_analysis`**: 整合Cut&Tag和RNA-seq数据
    *   `scripts/01_intersect_genes.R`: 整合Cut&Tag的峰值数据（特别是与基因关联的峰）和RNA-seq的差异表达基因，以鉴定潜在的直接调控靶基因。
    *   `scripts/02_run_enrichment_analysis.R`: 对整合分析得到的基因列表进行功能富集分析。
    *   `scripts/03_plot_enrichment_results.R`: 对富集分析结果进行可视化。
    *   `scripts/04_plot_heatmap.R`: 使用 `ComplexHeatmap` 对特定的基因集（如KEGG通路基因和自定义标记基因）绘制表达热图，进行可视化分析。

**运行方式**:
本地分析脚本（主要是R脚本）在本地计算机环境中运行。

## 4. 环境配置

### 4.1 服务器端环境

*   **Conda 环境**: `cuttag`
*   **激活方式**: 在执行服务器端脚本前，使用 `source activate cuttag` 命令激活环境。
*   **主要工具**: FastQC, MultiQC, Cutadapt, Bowtie2, SAMtools, Picard, MACS3, SEACR, BEDtools, deepTools, GNU Parallel。
    确保这些工具已在 `cuttag` 环境中正确安装。脚本中通常包含依赖检查步骤。

*   **IGV可视化工具**: 在 `cuttag` 环境中安装 `pyGenomeTracks`
*   **安装命令**: `conda install -c bioconda -c conda-forge pygenometracks` (在cuttag环境中执行)
*   **用途**: 生成IGV风格的基因组浏览器图（通过 `08_plot_igv_style.sh` 脚本）

### 4.2 本地环境

*   **R 环境**: R (建议版本 4.0 或更高)
*   **主要R包**: `tidyverse` (包含 `ggplot2`, `dplyr`, `tidyr` 等), `ChIPseeker`, `clusterProfiler`, `edgeR`, `RColorBrewer`。
    请确保这些R包已在本地R环境中安装。

## 5. 文件命名与结果组织

请严格遵守项目级规则中关于文件命名和结果文件组织的规范。这有助于保持项目的一致性和可维护性。

## 6. 参考流程

本项目在设计和开发过程中参考了以下优秀的公开流程和教程：

*   **诺唯赞流程**: `reference_repos/诺唯赞流程/cuttag_cutrun分析流程v2414/` (主要参考)
*   **CebolaLab流程**: `reference_repos/CebolaLab_CUTandTAG_ Analysis pipeline for CUT&TAG data.html`
*   **nf-core/cutandrun**: `reference_repos/nf-core_cutandrun/`
*   **Henikoff Lab CUT&Tag Tutorial**: `reference_repos/Henikoff_CUTTag_tutorial/`

## 7. 注意事项

*   **环境分离**: `CutTag_Linux` 目录下的脚本仅能在服务器环境运行，`CutTag_local_analysis` 和 `Integrative_analysis` 目录下的脚本仅能在本地环境运行。
*   **配置文件**: 核心配置（如文件路径、样本列表）集中在 `CutTag_Linux/scripts/config/00_config.sh` 中管理。
*   **SGE提交**: 服务器端计算密集型任务通过SGE提交，遵循特定的脚本格式和依赖管理，详见项目级规则。
*   **抽样测试**: 在进行完整分析前，建议使用 `00_subsample.sh` 对少量数据进行抽样测试，确保流程稳定。测试输出与正式分析输出严格分离。

## 8. 作者与贡献

请在各脚本文件头部根据项目级规则中的作者署名规范，明确标注作者和贡献者信息。

---

*本文档将根据项目进展持续更新。请以实际的项目结构和脚本为准。* 