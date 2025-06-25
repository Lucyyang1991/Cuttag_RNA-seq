#!/bin/bash
#$ -S /bin/bash
#$ -N alignment_array   # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=16G       # 申请16GB内存/核心
#$ -pe smp 8           # 申请8个CPU线程
#$ -t 1-5              # 任务数组（5个样本）

# ===========================================================================
# 脚本名称: 03_alignment.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 对Cut&Tag测序数据进行比对分析，包括：
#   - 使用Bowtie2进行基因组比对
#   - 对BAM文件进行排序
#   - 标记重复序列
#   - 过滤低质量和未比对的reads
#   - 生成比对统计信息
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 单独运行比对：
#    标准运行: 
#      qsub -o analysis_results/logs/alignment/alignment_array.\$TASK_ID.log scripts/03_alignment.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/alignment/alignment_array.\$TASK_ID.log scripts/03_alignment.sh -test
#
# 2. 运行比对并自动生成统计报告（推荐）：
#    标准运行: 
#      # 第一步：运行比对分析
#      qsub -o analysis_results/logs/alignment/alignment_array.\$TASK_ID.log scripts/03_alignment.sh
#      # 第二步：生成统计报告
#      qsub -o analysis_results/logs/alignment/alignment_stats.log -hold_jid alignment_array scripts/03_alignment_stats.sh
#
#    测试运行: 
#      # 第一步：运行比对分析
#      qsub -o test_analysis_results/logs/alignment/alignment_array.\$TASK_ID.log scripts/03_alignment.sh -test
#      # 第二步：生成统计报告
#      qsub -o test_analysis_results/logs/alignment/alignment_stats.log -hold_jid alignment_array scripts/03_alignment_stats.sh -test
#
# 输入：
#   标准运行: analysis_results/02_trimming/*_R{1,2}_trimmed.fastq.gz
#   测试运行: test_analysis_results/02_trimming/*_R{1,2}_trimmed.fastq.gz
# 输出：
#   标准运行: 
#     - analysis_results/03_alignment/raw/*.sorted.bam
#     - analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#     - analysis_results/03_alignment/metrics/*.{bowtie2.log,markdup.metrics,flagstat,idxstats}
#   测试运行: 
#     - test_analysis_results/03_alignment/raw/*.sorted.bam
#     - test_analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#     - test_analysis_results/03_alignment/metrics/*.{bowtie2.log,markdup.metrics,flagstat,idxstats}
# 依赖：bowtie2, samtools, picard
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
    for tool in bowtie2 samtools picard; do
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
    INPUT_DIR=${TEST_TRIMMED_DIR}
    OUTPUT_DIR=${TEST_ALIGNMENT_DIR}
    LOG_OUTPUT_DIR=${TEST_ALIGNMENT_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${TRIMMED_DIR}
    OUTPUT_DIR=${ALIGNMENT_DIR}
    LOG_OUTPUT_DIR=${ALIGNMENT_LOG}
    log "运行模式: 标准模式"
fi

# 设置输出子目录
BAM_RAW_DIR="${OUTPUT_DIR}/raw"
BAM_FILTERED_DIR="${OUTPUT_DIR}/filtered"
METRICS_DIR="${OUTPUT_DIR}/metrics"

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 检查输入目录是否存在
if [ ! -d "${INPUT_DIR}" ]; then
    log "错误: 输入目录不存在: ${INPUT_DIR}"
    log "请确保已经运行了02_trimming.sh脚本"
    exit 1
fi

# 创建必要的目录
for dir in "${OUTPUT_DIR}" "${BAM_RAW_DIR}" "${BAM_FILTERED_DIR}" "${METRICS_DIR}" "${LOG_OUTPUT_DIR}"; do
    mkdir -p "${dir}"
done

# 检查bowtie2索引文件
if [ ! -f "${BOWTIE2_INDEX}.1.bt2" ]; then
    log "错误: Bowtie2索引文件不存在: ${BOWTIE2_INDEX}"
    exit 1
fi

# 获取样本文件
log "获取样本文件列表..."
# 获取所有*R1_trimmed.fastq.gz文件
R1_FILES=($(ls ${INPUT_DIR}/*_R1_trimmed.fastq.gz 2>/dev/null || true))
if [ ${#R1_FILES[@]} -eq 0 ]; then
    log "错误: ${INPUT_DIR}目录中没有找到*_R1_trimmed.fastq.gz文件"
    exit 1
fi

# 获取当前任务的样本文件（R1和R2）
CURRENT_R1=${R1_FILES[$((SGE_TASK_ID-1))]}
CURRENT_R2=${CURRENT_R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}
CURRENT_SAMPLE=$(basename "${CURRENT_R1}" _R1_trimmed.fastq.gz)

# 检查文件是否存在
if [ ! -f "${CURRENT_R1}" ]; then
    log "错误: 找不到R1文件 ${CURRENT_R1}"
    exit 1
fi
if [ ! -f "${CURRENT_R2}" ]; then
    log "错误: 找不到R2文件 ${CURRENT_R2}"
    exit 1
fi

log "开始处理双端文件:"
log "R1: ${CURRENT_R1}"
log "R2: ${CURRENT_R2}"
log "总样本数: ${#R1_FILES[@]}, 当前处理第 ${SGE_TASK_ID} 对"

# 1. Bowtie2比对
log "运行bowtie2比对..."
bowtie2 \
    -p ${NSLOTS} \
    --very-sensitive-local \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 -X 700 \
    -x ${BOWTIE2_INDEX} \
    -1 ${CURRENT_R1} \
    -2 ${CURRENT_R2} \
    2> "${METRICS_DIR}/${CURRENT_SAMPLE}.bowtie2.log" \
    | samtools view -bS - > "${BAM_RAW_DIR}/${CURRENT_SAMPLE}.bam"

# 2. 排序BAM文件
log "对BAM文件进行排序..."
samtools sort \
    -@ ${NSLOTS} \
    -o "${BAM_RAW_DIR}/${CURRENT_SAMPLE}.sorted.bam" \
    "${BAM_RAW_DIR}/${CURRENT_SAMPLE}.bam"

# 3. 标记重复序列
log "标记重复序列..."
picard MarkDuplicates \
    INPUT="${BAM_RAW_DIR}/${CURRENT_SAMPLE}.sorted.bam" \
    OUTPUT="${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.marked.bam" \
    METRICS_FILE="${METRICS_DIR}/${CURRENT_SAMPLE}.markdup.metrics" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

# 4. 过滤低质量和未比对的reads
log "过滤低质量reads..."
samtools view \
    -@ ${NSLOTS} \
    -F 1804 \
    -q 30 \
    -b "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.marked.bam" \
    > "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.bam"

# 5. 对过滤后的BAM文件进行排序和索引
log "对过滤后的BAM文件排序和索引..."
samtools sort \
    -@ ${NSLOTS} \
    -o "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.sorted.bam" \
    "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.bam"

samtools index "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.sorted.bam"

# 6. 生成比对统计信息
log "生成比对统计信息..."
samtools flagstat \
    "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.sorted.bam" \
    > "${METRICS_DIR}/${CURRENT_SAMPLE}.flagstat"

# 7. 计算测序深度
log "计算测序深度..."
samtools depth \
    "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.sorted.bam" \
    | awk '{sum+=$3} END {print "Average depth: "sum/NR}' \
    >> "${METRICS_DIR}/${CURRENT_SAMPLE}.flagstat"

# 8. 线粒体比对统计
log "统计线粒体比对情况..."
samtools idxstats "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.sorted.bam" \
    > "${METRICS_DIR}/${CURRENT_SAMPLE}.idxstats"

# 清理中间文件
log "清理中间文件..."
rm -f "${BAM_RAW_DIR}/${CURRENT_SAMPLE}.bam"
rm -f "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.filtered.bam"
rm -f "${BAM_FILTERED_DIR}/${CURRENT_SAMPLE}.marked.bam"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0