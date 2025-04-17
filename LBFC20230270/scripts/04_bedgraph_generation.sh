#!/bin/bash
#$ -S /bin/bash
#$ -N bedgraph_array    # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=8G        # 申请8GB内存/核心
#$ -pe smp 4           # 申请4个CPU线程
#$ -t 1-5              # 任务数组（根据样本数自动调整）

# ===========================================================================
# 脚本名称: 04_bedgraph_generation.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 将BAM文件转换为标准化的bedgraph文件，用于后续peak calling分析
#   - 使用bedtools genomecov进行转换
#   - 应用CPM（Counts Per Million）标准化
#   - 生成详细的处理metrics
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 推荐用法（与peak calling脚本一起运行）：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/bedgraph_array.\$TASK_ID.log -hold_jid alignment_array scripts/04_bedgraph_generation.sh
#      qsub -hold_jid bedgraph_array scripts/04_peak_calling_seacr.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/bedgraph_array.\$TASK_ID.log scripts/04_bedgraph_generation.sh -test
#      qsub -hold_jid bedgraph_array scripts/04_peak_calling_seacr.sh -test
#
# 2. 单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/bedgraph_array.\$TASK_ID.log scripts/04_bedgraph_generation.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/bedgraph_array.\$TASK_ID.log scripts/04_bedgraph_generation.sh -test
#
# 输入：
#   标准运行: analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#   测试运行: test_analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#
# 输出：
#   标准运行:
#     - bedgraph文件: analysis_results/04_peaks/bedgraph/*.normalized.bedgraph
#     - metrics文件: analysis_results/04_peaks/bedgraph/metrics/*.bedgraph.metrics
#   测试运行:
#     - bedgraph文件: test_analysis_results/04_peaks/bedgraph/*.normalized.bedgraph
#     - metrics文件: test_analysis_results/04_peaks/bedgraph/metrics/*.bedgraph.metrics
#
# 依赖：bedtools, samtools, bc
# ===========================================================================

# 记录开始时间
start_time=$(date +%s)

# 设置错误处理
set -e
trap 'echo "错误发生在第 $LINENO 行"; exit 1' ERR

# 日志输出函数
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# 获取配置文件路径
CONFIG_FILE="${SGE_O_WORKDIR}/scripts/config/00_config.sh"
if [ ! -f "${CONFIG_FILE}" ]; then
    log "错误: 找不到配置文件 ${CONFIG_FILE}"
    exit 1
fi
source ${CONFIG_FILE}

# 检查必要的软件依赖
check_dependencies() {
    local missing_tools=()
    for tool in bedtools samtools bc; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log "错误: 以下必要工具未安装或未添加到PATH中:"
        for tool in "${missing_tools[@]}"; do
            log "  - $tool"
        done
        log "请确保已激活正确的conda环境(cuttag)，或安装缺失的工具。"
        exit 1
    fi
}

# 激活conda环境
log "激活conda环境: cuttag"
source activate cuttag
check_dependencies

# 处理命令行参数
if [ "$1" == "-test" ]; then
    INPUT_DIR=${TEST_ALIGNMENT_DIR}/filtered
    OUTPUT_DIR=${TEST_PEAKS_DIR}
    LOG_OUTPUT_DIR=${TEST_PEAKCALLING_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${ALIGNMENT_DIR}/filtered
    OUTPUT_DIR=${PEAKS_DIR}
    LOG_OUTPUT_DIR=${PEAKCALLING_LOG}
    log "运行模式: 标准模式"
fi

# 设置工作目录
BEDGRAPH_DIR="${OUTPUT_DIR}/bedgraph"
METRICS_DIR="${BEDGRAPH_DIR}/metrics"
TEMP_DIR="${OUTPUT_DIR}/temp"

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 创建必要的目录
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BEDGRAPH_DIR}"
mkdir -p "${METRICS_DIR}"
mkdir -p "${LOG_OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

# 获取样本列表
mapfile -t ALL_SAMPLES < <(ls ${INPUT_DIR}/*.filtered.sorted.bam 2>/dev/null | sed 's/.*\///' | sed 's/\.filtered\.sorted\.bam//')
if [ ${#ALL_SAMPLES[@]} -eq 0 ]; then
    log "错误：在${INPUT_DIR}目录下未找到任何.filtered.sorted.bam文件"
    exit 1
fi

# 获取当前任务的样本
CURRENT_SAMPLE=${ALL_SAMPLES[$((SGE_TASK_ID-1))]}
if [ -z "${CURRENT_SAMPLE}" ]; then
    log "错误：无法获取任务ID ${SGE_TASK_ID} 对应的样本"
    exit 1
fi

log "处理样本: ${CURRENT_SAMPLE}"

# 生成标准化的bedgraph文件
generate_normalized_bedgraph() {
    local sample=$1
    local input_bam="${INPUT_DIR}/${sample}.filtered.sorted.bam"
    local output_bedgraph="${BEDGRAPH_DIR}/${sample}.normalized.bedgraph"
    local metrics_file="${METRICS_DIR}/${sample}.bedgraph.metrics"
    local temp_prefix="${TEMP_DIR}/${sample}"
    local namesort_bam="${BEDGRAPH_DIR}/${sample}.namesorted.bam"
    
    log "处理样本 ${sample}..."
    
    # 检查输入文件
    if [ ! -f "${input_bam}" ]; then
        log "错误: 找不到输入文件 ${input_bam}"
        return 1
    fi
    
    # 首先按read name排序BAM文件
    log "按read name排序BAM文件..."
    samtools sort -n \
        -@ ${NSLOTS} \
        -m 2G \
        "${input_bam}" \
        -o "${namesort_bam}"
    
    # 生成bedpe文件
    log "生成bedpe文件..."
    bedtools bamtobed -bedpe -i "${namesort_bam}" > "${temp_prefix}.bed"
    
    # 过滤同一染色体且片段长度<1000bp的reads
    log "过滤reads..."
    awk '$1==$4 && $6-$2 < 1000 {print $0}' "${temp_prefix}.bed" > "${temp_prefix}.clean.bed"
    
    # 只保留片段的起始和终止位置
    log "提取片段位置..."
    cut -f 1,2,6 "${temp_prefix}.clean.bed" | \
    sort -k1,1 -k2,2n -k3,3n > "${temp_prefix}.fragments.bed"
    
    # 生成bedgraph文件
    log "生成bedgraph文件..."
    bedtools genomecov -bg \
        -i "${temp_prefix}.fragments.bed" \
        -g "${GENOME_CHROM_SIZES}" \
        > "${output_bedgraph}"
    
    # 记录统计信息
    local total_fragments=$(wc -l < "${temp_prefix}.fragments.bed")
    local filtered_pairs=$(wc -l < "${temp_prefix}.clean.bed")
    local original_pairs=$(wc -l < "${temp_prefix}.bed")
    
    # 记录处理信息
    {
        echo "Sample: ${sample}"
        echo "Input BAM: ${input_bam}"
        echo "Name-sorted BAM: ${namesort_bam}"
        echo "Output bedgraph: ${output_bedgraph}"
        echo "Original read pairs: ${original_pairs}"
        echo "Filtered read pairs: ${filtered_pairs}"
        echo "Final fragments: ${total_fragments}"
        echo "Generated: $(date)"
        echo "Processing steps:"
        echo "1. Sort BAM by read name: samtools sort -n -@ ${NSLOTS}"
        echo "2. Convert BAM to BEDPE: bedtools bamtobed -bedpe"
        echo "3. Filter reads: same chromosome and length < 1000bp"
        echo "4. Extract fragment coordinates: cut -f 1,2,6"
        echo "5. Generate bedgraph: bedtools genomecov -bg"
    } > "${metrics_file}"
    
    # 清理临时文件
    log "清理临时文件..."
    rm "${temp_prefix}.bed" "${temp_prefix}.clean.bed" "${temp_prefix}.fragments.bed"
    
    log "完成bedgraph生成: ${output_bedgraph}"
}

# 执行bedgraph生成
log "生成标准化的bedgraph文件..."
generate_normalized_bedgraph "${CURRENT_SAMPLE}"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 