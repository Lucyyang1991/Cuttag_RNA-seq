#!/bin/bash

# CUT&TAG数据分析流程 - 主配置文件
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")

#===============================
# 工作目录配置
#===============================
WORK_DIR="/data/home/lucyyang/Project/IKZF2_CUT_TAG"

# 创建主要目录
mkdir -p ${WORK_DIR}/{scripts,raw_data/{test_samples,production},genome,analysis_results,test_analysis_results}

# 创建标准分析结果目录
mkdir -p ${WORK_DIR}/analysis_results/01_fastqc/raw \
         ${WORK_DIR}/analysis_results/02_trimmed \
         ${WORK_DIR}/analysis_results/03_alignment/metrics \
         ${WORK_DIR}/analysis_results/04_peaks/individual_peaks \
         ${WORK_DIR}/analysis_results/05_visualization/heatmaps \
         ${WORK_DIR}/analysis_results/05_visualization/correlation \
         ${WORK_DIR}/analysis_results/logs/fastqc \
         ${WORK_DIR}/analysis_results/logs/trimming \
         ${WORK_DIR}/analysis_results/logs/alignment \
         ${WORK_DIR}/analysis_results/logs/peak_calling

# 创建测试分析结果目录
mkdir -p ${WORK_DIR}/test_analysis_results/01_fastqc/raw \
         ${WORK_DIR}/test_analysis_results/02_trimmed \
         ${WORK_DIR}/test_analysis_results/03_alignment/metrics \
         ${WORK_DIR}/test_analysis_results/04_peaks/individual_peaks \
         ${WORK_DIR}/test_analysis_results/05_visualization/heatmaps \
         ${WORK_DIR}/test_analysis_results/05_visualization/correlation \
         ${WORK_DIR}/test_analysis_results/logs/fastqc \
         ${WORK_DIR}/test_analysis_results/logs/trimming \
         ${WORK_DIR}/test_analysis_results/logs/alignment \
         ${WORK_DIR}/test_analysis_results/logs/peak_calling

# 设置目录路径
SCRIPTS_DIR="${WORK_DIR}/scripts"
CONFIG_DIR="${SCRIPTS_DIR}/config"
RAW_DATA_DIR="${WORK_DIR}/raw_data/production"
TEST_DATA_DIR="${WORK_DIR}/raw_data/test_samples"
RESULT_DIR="${WORK_DIR}/analysis_results"
TEST_RESULT_DIR="${WORK_DIR}/test_analysis_results"

# 设置结果子目录（正式）
FASTQC_DIR="${RESULT_DIR}/01_fastqc"
TRIMMED_DIR="${RESULT_DIR}/02_trimmed"
ALIGNMENT_DIR="${RESULT_DIR}/03_alignment"
PEAKS_DIR="${RESULT_DIR}/04_peaks"
VISUALIZATION_DIR="${RESULT_DIR}/05_visualization"

# 设置结果子目录（测试）
TEST_FASTQC_DIR="${TEST_RESULT_DIR}/01_fastqc"
TEST_TRIMMED_DIR="${TEST_RESULT_DIR}/02_trimmed"
TEST_ALIGNMENT_DIR="${TEST_RESULT_DIR}/03_alignment"
TEST_PEAKS_DIR="${TEST_RESULT_DIR}/04_peaks"
TEST_VISUALIZATION_DIR="${TEST_RESULT_DIR}/05_visualization"

# 设置日志目录（正式）
LOG_DIR="${RESULT_DIR}/logs"
FASTQC_LOG="${LOG_DIR}/fastqc"
TRIMMING_LOG="${LOG_DIR}/trimming"
ALIGNMENT_LOG="${LOG_DIR}/alignment"
PEAKCALLING_LOG="${LOG_DIR}/peak_calling"
VISUALIZATION_LOG="${LOG_DIR}/visualization"

# 设置日志目录（测试）
TEST_LOG_DIR="${TEST_RESULT_DIR}/logs"
TEST_FASTQC_LOG="${TEST_LOG_DIR}/fastqc"
TEST_TRIMMING_LOG="${TEST_LOG_DIR}/trimming"
TEST_ALIGNMENT_LOG="${TEST_LOG_DIR}/alignment"
TEST_PEAKCALLING_LOG="${TEST_LOG_DIR}/peak_calling"
TEST_VISUALIZATION_LOG="${TEST_LOG_DIR}/visualization"
#===============================
# 参考基因组配置
#===============================
GENOME_DIR="${WORK_DIR}/genome"
BOWTIE2_INDEX="${GENOME_DIR}/bowtie2_index/mm10"
GENOME_FA="${GENOME_DIR}/mm10.fa"
GTF_FILE="${GENOME_DIR}/mm10.refGene.gtf"
GENOME_CHROM_SIZES="${GENOME_DIR}/mm10.chrom.sizes"  # 添加染色体大小文件

#===============================
# 实验参数配置
#===============================
# 实验类型设置
EXPERIMENT_TYPE="CUT&Tag"  # 可选: CUT&Tag 或 CUT&RUN

# 接头序列配置
if [ "${EXPERIMENT_TYPE}" == "CUT&Tag" ]; then
    ADAPTER1="CTGTCTCTTATAC"
    ADAPTER2="CTGTCTCTTATAC"
else
    ADAPTER1="AGATCGGAAGAGCA"
    ADAPTER2="AGATCGGAAGAGCA"
fi
POLYG_ADAPTER="GGGGGGGGGGGGX"

# 质控参数
MIN_LENGTH=18
QUALITY_CUTOFF=30
MAX_N=0.05
ERROR_RATE=0.2
ADAPTER_OVERLAP=3

#===============================
# SEACR峰值检测参数
#===============================
# 根据 Henikoff 实验室的 CUTTag_tutorial 和 nf-core/cutandrun 的设置
SEACR_NORM="norm"           # 由于没有使用spike-in标准化，使用norm进行IgG对照标准化
SEACR_STRINGENCY="stringent"  # 使用严格模式，减少假阳性
SEACR_THRESHOLD="0.01"     # 当没有对照时使用的阈值

#===============================
# MACS3峰值检测参数
#===============================
# 根据多个CUT&Tag分析流程的综合建议设置
MACS3_QVALUE="0.05"        # 使用较宽松的q值阈值捕获更多可能的结合位点
MACS3_SHIFT="-75"          # 偏移修正值，调整Tn5转座酶插入偏好
MACS3_EXTSIZE="150"        # 延伸大小，对应于Cut&Tag片段的平均大小
MACS3_GENOME="mm"          # 基因组大小参数，mm为小鼠，hs为人类
MACS3_MODE="narrow"        # 窄峰模式，适用于转录因子，宽峰模式适用于组蛋白修饰
MACS3_BROAD_CUTOFF="0.1"   # 宽峰模式的峰值阈值，仅在MACS3_MODE="broad"时使用

#===============================
# BigWig生成参数
#===============================
# deepTools bamCoverage参数
BIGWIG_BIN_SIZE=50
BIGWIG_NORMALIZE="RPKM"
# BIGWIG_SMOOTH_LENGTH=50
# BIGWIG_EXTEND_READS=150

#===============================
# Heatmap参数
#===============================
HEATMAP_BIN_SIZE=50
HEATMAP_UPSTREAM=5000
HEATMAP_DOWNSTREAM=5000
HEATMAP_BODY_LENGTH=5000
# HEATMAP_COLOR="blue,yellow,red"
# PROFILE_COLORS="blue red green purple orange"

# 如果直接执行此脚本，输出配置信息
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # 输出配置信息
    echo "===== CUT&TAG 分析基础配置信息 ====="
    echo "工作目录: ${WORK_DIR}"
    echo "正式样本目录: ${FASTQ_DIR}"
    echo "测试样本目录: ${TEST_FASTQ_DIR}"
    echo "实验类型: ${EXPERIMENT_TYPE}"
    echo "参考基因组: ${GENOME_FA}"
    echo "==============================="

    # 输出预处理配置信息
    echo "===== 预处理配置信息 ====="
    echo "使用的接头序列:"
    echo "  Adapter1: ${ADAPTER1}"
    echo "  Adapter2: ${ADAPTER2}"
    echo "  PolyG: ${POLYG_ADAPTER}"
    echo "质控参数:"
    echo "  最小长度: ${MIN_LENGTH}"
    echo "  质量阈值: ${QUALITY_CUTOFF}"
    echo "  最大N比例: ${MAX_N}"
    echo "  错误率: ${ERROR_RATE}"
    echo "  接头重叠: ${ADAPTER_OVERLAP}"
    echo "========================="
fi