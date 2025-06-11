# Cut&Tag测序数据分析流程

## 1. 项目概述

本项目是一个整合了Cut&Tag测序和RNA-seq数据的分析流程，旨在鉴定转录因子结合位点，并结合基因表达变化分析其调控功能。项目包含Cut&Tag数据的完整处理流程（从质控、比对到峰值检测和注释）以及RNA-seq的差异表达分析，并最终对两者进行整合分析。

## 2. 目录结构

项目的建议目录结构如下：

```
project_root/
├── LBFC20230270/             # 原始数据及处理脚本 (服务器端运行)
│   ├── config/               # 配置文件 (如 00_config.sh)
│   ├── scripts/              # 各步骤处理脚本 (如 01_fastqc.sh, 02_trimming.sh 等)
│   ├── raw_data/             # 原始FASTQ文件存放目录
│   ├── genome/               # 参考基因组及注释文件
│   ├── analysis_results/     # 正式分析结果输出目录
│   │   ├── 01_fastqc/
│   │   ├── 02_trimmed/
│   │   ├── 03_alignment/
│   │   ├── 04_peaks/
│   │   ├── 05_visualization/
│   │   └── logs/
│   └── test_analysis_results/ # 测试分析结果输出目录 (结构同analysis_results)
├── local_analysis/           # CUT&Tag下游本地分析 (本地环境运行)
│   ├── peak_annotation/      # 峰值注释结果
│   ├── functional_analysis/  # (CUT&Tag峰值相关基因的)功能富集分析结果
│   ├── visualization/        # 可视化结果图表
│   └── *_*.R                 # R分析脚本
├── RNA-seq_analysis/         # RNA-seq 数据分析 (本地环境运行)
│   ├── analysis_results/     # RNA-seq 分析结果 (如差异表达基因列表)
│   └── enrichment_results/   # (可能包含RNA-seq或整合分析的)富集结果
└── reference_repos/          # 参考分析流程 (学习和参考用途)
    ├── CebolaLab_CUTandTAG_ Analysis pipeline for CUT&TAG data_files/
    ├── Henikoff_CUTTag_tutorial/
    ├── nf-core_cutandrun/
    └── 诺唯赞流程/
```

详细的目录结构规范请参考 `.cursor\rules\cuttag_workflow.mdc` 中的第3节和第9节。

## 3. 分析流程

本项目的分析流程主要分为服务器端分析和本地分析两部分。

### 3.1 服务器端分析流程 (位于 `LBFC20230270` 目录)

服务器端分析流程主要处理原始测序数据，包括以下步骤：

1.  **数据准备**:
    *   `00_config.sh`: 配置全局变量、样本信息和路径。
    *   `00_subsample.sh` (可选): 对原始数据进行抽样，用于流程测试。
2.  **数据预处理**:
    *   `01_fastqc.sh` (可选): 使用 FastQC 对原始数据进行质量评估。
    *   `01_fastqc_stats.sh`: 汇总FastQC结果并生成统计报告。
    *   `02_trimming.sh`: 使用 Cutadapt 去除测序接头和低质量序列。
    *   `02_trimming_stats.sh`: 汇总去接头结果并生成统计报告。
3.  **序列比对**:
    *   `03_alignment.sh`: 使用 Bowtie2 将处理后的序列比对到参考基因组，并进行BAM文件处理（排序、标记重复、过滤）。
    *   `03_alignment_stats.sh`: 生成详细的比对质量报告。
4.  **峰值检测 (Peak Calling)**:
    *   `04_bedgraph_generation.sh`: 生成标准化的 BedGraph 文件，作为Peak Calling的输入。
    *   `04_peak_calling_macs3.sh`: 使用 MACS3 进行峰值检测（推荐）。
    *   `04_peak_calling_seacr.sh`: (可选) 使用 SEACR 进行峰值检测。
    *   `04_peak_calling_stats.sh`: 生成详细的Peak Calling统计报告。
5.  **下游可视化准备**:
    *   `05_bam_to_bigwig.sh`: 将比对产生的BAM文件转换为BigWig格式，用于在IGV等基因组浏览器中可视化。
6.  **高级分析 (可选)**:
    *   `06_generate_heatmap.sh`: 使用 deepTools 生成TSS区域的信号富集热图。
    *   `07_sample_correlation.sh`: 使用 deepTools 分析样本间的相关性。
7.  **IGV风格可视化**:
    *   `08_plot_igv_style.sh`: 使用 pyGenomeTracks 生成IGV风格的基因组浏览器图，展示特定基因区域的信号富集模式。

**运行方式**:
服务器端的脚本设计为通过SGE (Sun Grid Engine) 作业调度系统提交和运行。具体提交方式请参考 `.cursor\rules\cuttag_workflow.mdc` 中的第6节"并行处理规范"和第7节中的 `qsub.sh` 脚本说明。

### 3.2 本地分析流程 (位于 `local_analysis` 和 `RNA-seq_analysis` 目录)

本地分析流程主要对服务器端产生的Cut&Tag峰值文件等结果，以及RNA-seq数据进行下游的注释、差异表达分析、功能富集和整合分析。

1.  **Cut&Tag峰值注释 (主要在 `local_analysis` 目录)**:
    *   使用R包 `ChIPseeker` 对MACS3或SEACR产生的peak文件进行注释，分析其在基因组上的分布特征。
    *   `07_filter_proximal_promoter_genes.R` (示例): 筛选近端启动子区域的峰值及其关联基因。
2.  **RNA-seq差异表达分析 (主要在 `RNA-seq_analysis` 目录)**:
    *   使用 `edgeR` (例如通过 `edgeR_analysis.R` 脚本) 或其他适用工具进行配对或分组比较，分析差异表达基因。
    *   生成火山图、热图等可视化结果。
    *   根据项目规则，可能特别关注某些基因集（如一碳代谢相关基因）的表达变化。
3.  **整合分析 (结果可能存放于 `local_analysis` 或 `RNA-seq_analysis/enrichment_results`)**: 
    *   整合Cut&Tag的峰值数据（特别是与基因关联的峰，如启动子区域）和RNA-seq的差异表达基因。
    *   例如，筛选近端启动子区域存在峰值且其关联基因差异表达的基因列表，以鉴定潜在的直接调控靶基因。
4.  **功能富集分析 (结果可能存放于 `local_analysis/functional_analysis` 或 `RNA-seq_analysis/enrichment_results`)**:
    *   `do_enrichment_analysis.R` (示例): 使用R包 `clusterProfiler` 对Cut&Tag关联基因、差异表达基因或整合分析得到的基因列表进行GO和KEGG等功能富集分析。
    *   分析关键生物学通路。
5.  **结果可视化 (分布于各分析子目录)**:
    *   使用 `ggplot2` 等R包对上述各步骤的分析结果进行可视化。

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

请严格遵守 `.cursor\rules\cuttag_workflow.mdc` 中关于文件命名（第2节）和结果文件组织（第9节）的规范。这有助于保持项目的一致性和可维护性。

## 6. 参考流程

本项目在设计和开发过程中参考了以下优秀的公开流程和教程：

*   **诺唯赞流程**: `reference_repos/诺唯赞流程/cuttag_cutrun分析流程v2414/` (主要参考)
*   **CebolaLab流程**: `reference_repos/CebolaLab_CUTandTAG_ Analysis pipeline for CUT&TAG data.html`
*   **nf-core/cutandrun**: `reference_repos/nf-core_cutandrun/`
*   **Henikoff Lab CUT&Tag Tutorial**: `reference_repos/Henikoff_CUTTag_tutorial/`

## 7. 注意事项

*   **环境分离**: `LBFC20230270` 目录下的脚本仅能在服务器环境运行，`local_analysis` 目录下的脚本仅能在本地环境运行。
*   **配置文件**: 核心配置（如文件路径、样本列表）集中在 `LBFC20230270/config/00_config.sh` 中管理。
*   **SGE提交**: 服务器端计算密集型任务通过SGE提交，遵循特定的脚本格式和依赖管理，详见 `.cursor\rules\cuttag_workflow.mdc`。
*   **抽样测试**: 在进行完整分析前，建议使用 `00_subsample.sh` 对少量数据进行抽样测试，确保流程稳定。测试输出与正式分析输出严格分离。

## 8. 作者与贡献

请在各脚本文件头部根据 `.cursor\rules\cuttag_workflow.mdc` 第4节的作者署名规范，明确标注作者和贡献者信息。

---

*本文档根据 `.cursor\rules\cuttag_workflow.mdc` 生成。请以 `.cursor\rules\cuttag_workflow.mdc` 中的详细规范为准。* 