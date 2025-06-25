#!/bin/bash
#$ -S /bin/bash
#$ -N macs3_array       # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=16G        # 申请16GB内存/核心，MACS3需要更多内存
#$ -pe smp 4            # 申请4个CPU线程
#$ -t 1-3               # 任务数组（根据处理组样本数自动调整）
#$ -hold_jid alignment_array  # 等待比对完成

# ===========================================================================
# 脚本名称: 04_peak_calling_macs3.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-05-11
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 使用MACS3对Cut&Tag数据进行peak calling分析
#   - 使用IgG作为对照样本
#   - 支持多个IgG对照的并行处理
#   - 针对Cut&Tag数据优化的参数设置
#   - 自动处理样本配对关系
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行
#
# 1. 推荐用法（在比对脚本后运行）：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/macs3_array.\$TASK_ID.log -hold_jid alignment_array scripts/04_peak_calling_macs3.sh
#      qsub -hold_jid macs3_array scripts/04_peak_calling_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/macs3_array.\$TASK_ID.log -hold_jid alignment_array scripts/04_peak_calling_macs3.sh -test
#      qsub -hold_jid macs3_array scripts/04_peak_calling_stats.sh -test
#
# 2. 单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/macs3_array.\$TASK_ID.log scripts/04_peak_calling_macs3.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/macs3_array.\$TASK_ID.log scripts/04_peak_calling_macs3.sh -test
#
# 输入：
#   标准运行:
#     - BAM文件: analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#     - 样本信息: scripts/config/sample_info.csv
#   测试运行:
#     - BAM文件: test_analysis_results/03_alignment/filtered/*.filtered.sorted.bam
#     - 样本信息: scripts/config/sample_info.csv
#
# 输出：
#   标准运行:
#     - MACS3 peaks: analysis_results/04_peaks/individual_peaks/*.macs3_peaks.narrowPeak
#   测试运行:
#     - MACS3 peaks: test_analysis_results/04_peaks/individual_peaks/*.macs3_peaks.narrowPeak
#
# 依赖：macs3, bedtools, awk
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
    for tool in macs3 bedtools awk; do
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
    OUTPUT_DIR=${TEST_PEAKS_DIR}
    INPUT_DIR=${TEST_ALIGNMENT_DIR}
    LOG_OUTPUT_DIR=${TEST_PEAKCALLING_LOG}
    log "运行模式: 测试模式"
else
    OUTPUT_DIR=${PEAKS_DIR}
    INPUT_DIR=${ALIGNMENT_DIR}
    LOG_OUTPUT_DIR=${PEAKCALLING_LOG}
    log "运行模式: 标准模式"
fi

# 设置工作目录
BAM_DIR="${INPUT_DIR}/filtered"
PEAK_OUTPUT_DIR="${OUTPUT_DIR}/individual_peaks"
TEMP_DIR="${OUTPUT_DIR}/temp"
SAMPLE_INFO="${SGE_O_WORKDIR}/scripts/config/sample_info.csv"

log "输入目录: ${BAM_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "Peak输出目录: ${PEAK_OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 检查样本信息文件
if [ ! -f "${SAMPLE_INFO}" ]; then
    log "错误: 找不到样本信息文件 ${SAMPLE_INFO}"
    exit 1
fi

# 创建必要的目录
mkdir -p "${PEAK_OUTPUT_DIR}" "${LOG_OUTPUT_DIR}" "${TEMP_DIR}"

# 获取处理组样本列表
mapfile -t TREATMENT_SAMPLES < <(awk -F',' 'NR>1 && $2=="treatment" {print $1}' "${SAMPLE_INFO}" | grep -v '^$')
if [ ${#TREATMENT_SAMPLES[@]} -eq 0 ]; then
    log "错误：在样本信息文件中未找到任何处理组样本"
    exit 1
fi

# 获取IgG对照样本列表
mapfile -t IGG_SAMPLES < <(awk -F',' 'NR>1 && $3=="IgG" {print $1}' "${SAMPLE_INFO}" | grep -v '^$')
if [ ${#IGG_SAMPLES[@]} -eq 0 ]; then
    log "错误：在样本信息文件中未找到任何IgG对照样本"
    exit 1
fi

# 获取当前任务的处理组样本
CURRENT_SAMPLE=${TREATMENT_SAMPLES[$((SGE_TASK_ID-1))]}
if [ -z "${CURRENT_SAMPLE}" ]; then
    log "错误：无法获取任务ID ${SGE_TASK_ID} 对应的样本"
    exit 1
fi

log "处理样本: ${CURRENT_SAMPLE}"

# 确定基因组大小参数
# 对于mm10使用适当的基因组大小值
GENOME_SIZE="${MACS3_GENOME}"  # 从配置文件读取

# 运行MACS3 peak calling
run_macs3() {
    local treatment=$1
    local control=$2
    local output_prefix="${PEAK_OUTPUT_DIR}/${treatment}_vs_${control}"
    local treatment_bam="${BAM_DIR}/${treatment}.filtered.sorted.bam"
    local control_bam="${BAM_DIR}/${control}.filtered.sorted.bam"
    
    # 检查输入文件
    if [ ! -f "${treatment_bam}" ]; then
        log "错误: 找不到处理组BAM文件 ${treatment_bam}"
        return 1
    fi
    if [ ! -f "${control_bam}" ]; then
        log "错误: 找不到对照组BAM文件 ${control_bam}"
        return 1
    fi
    
    log "运行MACS3: ${treatment} vs ${control}"
    log "  处理组: ${treatment_bam}"
    log "  对照组: ${control_bam}"
    log "  输出前缀: ${output_prefix}"
    
    # 运行MACS3 - Cut&Tag优化参数
    # 根据配置文件设置参数
    
    # 根据模式设置宽峰/窄峰参数
    local peak_mode_params=""
    if [ "${MACS3_MODE}" == "broad" ]; then
        peak_mode_params="--broad --broad-cutoff ${MACS3_BROAD_CUTOFF}"
    else
        peak_mode_params="--call-summits"
    fi
    
    macs3 callpeak \
        -t "${treatment_bam}" \
        -c "${control_bam}" \
        -f BAMPE \
        -g ${GENOME_SIZE} \
        -n "${treatment}_vs_${control}.macs3" \
        --outdir "${PEAK_OUTPUT_DIR}" \
        --keep-dup all \
        --nomodel \
        --shift ${MACS3_SHIFT} \
        --extsize ${MACS3_EXTSIZE} \
        -q ${MACS3_QVALUE} \
        ${peak_mode_params} \
        2> "${LOG_OUTPUT_DIR}/${treatment}_vs_${control}.macs3.log"
    
    log "MACS3分析完成: ${treatment} vs ${control}"
    
    # 复制并重命名结果文件，使其具有一致的命名格式
    if [ "${MACS3_MODE}" == "broad" ]; then
        cp "${PEAK_OUTPUT_DIR}/${treatment}_vs_${control}.macs3_peaks.broadPeak" "${PEAK_OUTPUT_DIR}/${treatment}_vs_${control}.macs3_peaks.bed"
    else
        cp "${PEAK_OUTPUT_DIR}/${treatment}_vs_${control}.macs3_peaks.narrowPeak" "${PEAK_OUTPUT_DIR}/${treatment}_vs_${control}.macs3_peaks.bed"
    fi
}

# 对当前处理组样本进行peak calling，使用每个IgG对照
for igg_sample in "${IGG_SAMPLES[@]}"; do
    log "使用IgG对照 ${igg_sample} 进行peak calling..."
    run_macs3 "${CURRENT_SAMPLE}" "${igg_sample}"
done

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 