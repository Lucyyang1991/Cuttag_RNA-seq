#!/bin/bash
#$ -S /bin/bash
#$ -N alignment_stats    # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=4G        # 申请4GB内存
#$ -pe smp 1           # 单线程即可

# ===========================================================================
# 脚本名称: 03_alignment_stats.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 生成Cut&Tag比对结果的统计报告，包括：
#   - 总reads数和比对率统计
#   - 唯一比对率统计
#   - PCR重复率统计
#   - 平均测序深度统计
#   - 线粒体DNA比率统计
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 与比对脚本同时提交（推荐）：
#    标准运行: 
#      qsub -o analysis_results/logs/alignment/alignment_stats.log -hold_jid alignment_array scripts/03_alignment_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/alignment/alignment_stats.log -hold_jid alignment_array scripts/03_alignment_stats.sh -test
#
# 2. 比对完成后单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/alignment/alignment_stats.log scripts/03_alignment_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/alignment/alignment_stats.log scripts/03_alignment_stats.sh -test
#
# 输入：
#   标准运行: analysis_results/03_alignment/metrics/*.{flagstat,markdup.metrics,idxstats}
#   测试运行: test_analysis_results/03_alignment/metrics/*.{flagstat,markdup.metrics,idxstats}
# 输出：
#   标准运行: analysis_results/03_alignment/metrics/alignment_statistics_YYYYMMDD.txt
#   测试运行: test_analysis_results/03_alignment/metrics/alignment_statistics_YYYYMMDD.txt
# 依赖：awk
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
    for tool in awk; do
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
    OUTPUT_DIR=${TEST_ALIGNMENT_DIR}
    METRICS_DIR="${OUTPUT_DIR}/metrics"
    LOG_OUTPUT_DIR=${TEST_ALIGNMENT_LOG}
    log "运行模式: 测试模式"
else
    OUTPUT_DIR=${ALIGNMENT_DIR}
    METRICS_DIR="${OUTPUT_DIR}/metrics"
    LOG_OUTPUT_DIR=${ALIGNMENT_LOG}
    log "运行模式: 标准模式"
fi

log "输出目录: ${OUTPUT_DIR}"
log "统计目录: ${METRICS_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 检查目录是否存在
if [ ! -d "${METRICS_DIR}" ]; then
    log "错误: metrics目录不存在: ${METRICS_DIR}"
    exit 1
fi

# 创建日志目录
mkdir -p ${LOG_OUTPUT_DIR}

# 等待一分钟确保所有文件都已写入完成
log "等待60秒确保所有文件都已完成写入..."
sleep 60

# 生成统计报告
log "生成汇总报告..."
STATS_FILE="${METRICS_DIR}/alignment_statistics_$(date '+%Y%m%d').txt"
echo -e "样本ID\t总reads数\t比对率(%)\t唯一比对率(%)\t重复率(%)\t平均测序深度\t线粒体比率(%)" > ${STATS_FILE}

# 获取所有样本的flagstat文件
FLAGSTAT_FILES=($(ls ${METRICS_DIR}/*.flagstat 2>/dev/null || true))
if [ ${#FLAGSTAT_FILES[@]} -eq 0 ]; then
    log "错误: ${METRICS_DIR}目录中没有找到flagstat文件"
    exit 1
fi

# 处理每个样本的统计文件
for sample_stats in ${FLAGSTAT_FILES[@]}; do
    sample_name=$(basename ${sample_stats} .flagstat)
    log "处理样本: ${sample_name}"
    
    # 检查其他必需文件
    MARKDUP_FILE="${METRICS_DIR}/${sample_name}.markdup.metrics"
    IDXSTATS_FILE="${METRICS_DIR}/${sample_name}.idxstats"
    
    for file in "${MARKDUP_FILE}" "${IDXSTATS_FILE}"; do
        if [ ! -f "$file" ]; then
            log "错误: 找不到文件: $file"
            exit 1
        fi
    done
    
    # 提取统计信息
    total_reads=$(grep "in total" ${sample_stats} | cut -d ' ' -f1)
    mapped_reads=$(grep "mapped (" ${sample_stats} | head -n1 | cut -d ' ' -f1)
    mapping_rate=$(awk "BEGIN {printf \"%.2f\", ${mapped_reads}/${total_reads}*100}")
    
    # 获取重复率
    dup_rate=$(grep -A1 "PERCENT_DUPLICATION" "${MARKDUP_FILE}" | tail -n1 | cut -f9)
    dup_rate=$(awk "BEGIN {printf \"%.2f\", ${dup_rate}*100}")
    
    # 获取平均深度
    avg_depth=$(tail -n1 ${sample_stats} | cut -d' ' -f3)
    
    # 计算线粒体比率
    mito_stats=$(grep "chrM" "${IDXSTATS_FILE}" || echo "0 0 0 0")
    total_mapped=$(awk '{sum+=$3} END {print sum}' "${IDXSTATS_FILE}")
    mito_reads=$(echo ${mito_stats} | awk '{print $3}')
    
    if [ "${total_mapped}" -eq 0 ]; then
        mito_ratio="0.00"
    else
        mito_ratio=$(awk -v mito=${mito_reads} -v total=${total_mapped} 'BEGIN {printf "%.2f", (mito/total)*100}')
    fi
    
    # 写入统计结果
    echo -e "${sample_name}\t${total_reads}\t${mapping_rate}\t${mapping_rate}\t${dup_rate}\t${avg_depth}\t${mito_ratio}" >> ${STATS_FILE}
    log "完成样本 ${sample_name} 的统计分析"
done

log "汇总报告已保存至: ${STATS_FILE}"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 