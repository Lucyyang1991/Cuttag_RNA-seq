#!/bin/bash
#$ -S /bin/bash
#$ -N bigwig_array      # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=8G        # 申请8GB内存/核心
#$ -pe smp 4           # 申请4个CPU线程
#$ -t 1-5              # 任务数组（根据样本数自动调整）

# ===========================================================================
# 脚本名称: 05_bam_to_bigwig.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 将BAM文件转换为BigWig格式用于IGV可视化
#   - 使用deeptools的bamCoverage进行转换
#   - 支持多种标准化方法（RPKM/CPM等）
#   - 生成IGV track列表和统计报告
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 标准运行：
#    qsub -o analysis_results/logs/visualization/bigwig_array.\$TASK_ID.log -hold_jid alignment_array scripts/05_bam_to_bigwig.sh
#
# 2. 测试运行：
#    qsub -o test_analysis_results/logs/visualization/bigwig_array.\$TASK_ID.log scripts/05_bam_to_bigwig.sh -test
#
# 输入：
#   标准运行: analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#   测试运行: test_analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#
# 输出：
#   标准运行:
#     - BigWig文件: analysis_results/visualization/bigwig/*.bw
#     - IGV track列表: analysis_results/visualization/bigwig/igv_tracks.txt
#     - 统计报告: analysis_results/logs/visualization/bigwig_statistics_*.txt
#   测试运行:
#     - BigWig文件: test_analysis_results/visualization/bigwig/*.bw
#     - IGV track列表: test_analysis_results/visualization/bigwig/igv_tracks.txt
#     - 统计报告: test_analysis_results/logs/visualization/bigwig_statistics_*.txt
#
# 依赖：deeptools, samtools, parallel
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
    for tool in bamCoverage samtools parallel; do
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
    OUTPUT_DIR=${TEST_VISUALIZATION_DIR}
    LOG_OUTPUT_DIR=${TEST_VISUALIZATION_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${ALIGNMENT_DIR}/filtered
    OUTPUT_DIR=${VISUALIZATION_DIR}
    LOG_OUTPUT_DIR=${VISUALIZATION_LOG}
    log "运行模式: 标准模式"
fi

# 设置工作目录
BIGWIG_DIR="${OUTPUT_DIR}/bigwig"

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 创建必要的目录
mkdir -p "${BIGWIG_DIR}"
mkdir -p "${LOG_OUTPUT_DIR}"

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

# 生成BigWig文件
generate_bigwig() {
    local sample=$1
    local input_bam="${INPUT_DIR}/${sample}.filtered.sorted.bam"
    local output_bw="${BIGWIG_DIR}/${sample}.bw"
    
    # 检查输入文件
    if [ ! -f "${input_bam}" ]; then
        log "错误: 找不到输入文件 ${input_bam}"
        return 1
    fi
    
    # 检查BAM索引
    if [ ! -f "${input_bam}.bai" ]; then
        log "为BAM文件创建索引..."
        samtools index ${input_bam}
    fi
    
    # 生成BigWig文件
    log "生成BigWig文件..."
    bamCoverage --bam ${input_bam} \
        --outFileName ${output_bw} \
        --outFileFormat bigwig \
        --binSize ${BIGWIG_BIN_SIZE} \
        --normalizeUsing ${BIGWIG_NORMALIZE} \
        # --smoothLength ${BIGWIG_SMOOTH_LENGTH} \
        # --extendReads ${BIGWIG_EXTEND_READS} \
        # --skipNonCoveredRegions \
        # --numberOfProcessors 4
        
    log "已生成BigWig文件: ${output_bw}"
}

# 执行BigWig生成
log "生成BigWig文件..."
generate_bigwig "${CURRENT_SAMPLE}"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 