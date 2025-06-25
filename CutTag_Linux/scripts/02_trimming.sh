#!/bin/bash
#$ -S /bin/bash
#$ -N trimming_array    # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=8G        # 申请8GB内存/核心
#$ -pe smp 4           # 申请4个CPU线程
#$ -t 1-5              # 任务数组（5个样本）

# ===========================================================================
# 脚本名称: 02_trimming.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 对Cut&Tag测序数据进行接头和低质量序列去除，包括：
#   - 去除测序接头序列
#   - 去除低质量碱基
#   - 去除polyG序列
#   - 去除过短的reads
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 单独运行Trimming：
#    标准运行: 
#      qsub -o analysis_results/logs/trimming/trimming_array.\$TASK_ID.log scripts/02_trimming.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/trimming/trimming_array.\$TASK_ID.log scripts/02_trimming.sh -test
#
# 2. 运行Trimming并自动生成统计报告（推荐）：
#    标准运行: 
#      # 第一步：运行Trimming分析
#      qsub -o analysis_results/logs/trimming/trimming_array.\$TASK_ID.log scripts/02_trimming.sh
#      # 第二步：生成统计报告
#      qsub -o analysis_results/logs/trimming/trimming_stats.log -hold_jid trimming_array scripts/02_trimming_stats.sh
#
#    测试运行: 
#      # 第一步：运行Trimming分析
#      qsub -o test_analysis_results/logs/trimming/trimming_array.\$TASK_ID.log scripts/02_trimming.sh -test
#      # 第二步：生成统计报告
#      qsub -o test_analysis_results/logs/trimming/trimming_stats.log -hold_jid trimming_array scripts/02_trimming_stats.sh -test
#
# 输入：
#   标准运行: analysis_results/01_fastqc/raw/*.fastq.gz
#   测试运行: test_analysis_results/01_fastqc/raw/*.fastq.gz
# 输出：
#   标准运行: 
#     - analysis_results/02_trimming/*_R{1,2}_trimmed.fastq.gz
#     - analysis_results/logs/trimming/*_cutadapt.log
#   测试运行: 
#     - test_analysis_results/02_trimming/*_R{1,2}_trimmed.fastq.gz
#     - test_analysis_results/logs/trimming/*_cutadapt.log
# 依赖：cutadapt
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
    for tool in cutadapt; do
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
    INPUT_DIR=${TEST_DATA_DIR}
    OUTPUT_DIR=${TEST_TRIMMED_DIR}
    LOG_OUTPUT_DIR=${TEST_TRIMMING_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${RAW_DATA_DIR}
    OUTPUT_DIR=${TRIMMED_DIR}
    LOG_OUTPUT_DIR=${TRIMMING_LOG}
    log "运行模式: 标准模式"
fi

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 检查输入目录是否存在
if [ ! -d "${INPUT_DIR}" ]; then
    log "错误: 输入目录不存在: ${INPUT_DIR}"
    exit 1
fi

# 创建输出目录
log "创建输出目录"
mkdir -p ${OUTPUT_DIR} ${LOG_OUTPUT_DIR}

# 获取样本文件
log "获取样本文件列表..."
# 获取所有*R1.fastq.gz文件
R1_FILES=($(ls ${INPUT_DIR}/*.R1.fastq.gz 2>/dev/null || true))
if [ ${#R1_FILES[@]} -eq 0 ]; then
    log "错误: ${INPUT_DIR}目录中没有找到*R1.fastq.gz文件"
    exit 1
fi

# 获取当前任务的样本文件（R1和R2）
CURRENT_R1=${R1_FILES[$((SGE_TASK_ID-1))]}
CURRENT_R2=${CURRENT_R1/R1.fastq.gz/R2.fastq.gz}
CURRENT_SAMPLE=$(basename "${CURRENT_R1}" .R1.fastq.gz)

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

# 运行cutadapt
log "运行Cutadapt分析..."
cutadapt \
    -a ${ADAPTER1} -A ${ADAPTER2} \
    -a ${POLYG_ADAPTER} -A ${POLYG_ADAPTER} \
    -m ${MIN_LENGTH} \
    -q ${QUALITY_CUTOFF},${QUALITY_CUTOFF} \
    --max-n ${MAX_N} \
    -e ${ERROR_RATE} \
    -O ${ADAPTER_OVERLAP} \
    -n 2 \
    --trim-n \
    --cores=${NSLOTS} \
    -o ${OUTPUT_DIR}/${CURRENT_SAMPLE}_R1_trimmed.fastq.gz \
    -p ${OUTPUT_DIR}/${CURRENT_SAMPLE}_R2_trimmed.fastq.gz \
    ${CURRENT_R1} \
    ${CURRENT_R2} \
    2>&1 | tee ${LOG_OUTPUT_DIR}/${CURRENT_SAMPLE}_cutadapt.log

# 检查输出文件
if [ ! -s "${OUTPUT_DIR}/${CURRENT_SAMPLE}_R1_trimmed.fastq.gz" ] || \
   [ ! -s "${OUTPUT_DIR}/${CURRENT_SAMPLE}_R2_trimmed.fastq.gz" ]; then
    log "错误: 输出文件不存在或为空"
    exit 1
fi

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0