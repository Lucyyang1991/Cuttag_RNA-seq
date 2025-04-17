#!/bin/bash
#$ -S /bin/bash
#$ -N trimming_stats    # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=4G        # 申请4GB内存
#$ -pe smp 1           # 单线程即可

# ===========================================================================
# 脚本名称: 02_trimming_stats.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 生成Trimming结果的统计报告，包括：
#   - 每个样本的原始reads数量
#   - 过滤后的reads数量
#   - reads保留率统计
#   - 接头去除和低质量序列过滤情况
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 与Trimming脚本同时提交（推荐）：
#    标准运行: 
#      qsub -o analysis_results/logs/trimming/trimming_stats.log -hold_jid trimming_array scripts/02_trimming_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/trimming/trimming_stats.log -hold_jid trimming_array scripts/02_trimming_stats.sh -test
#
# 2. Trimming完成后单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/trimming/trimming_stats.log scripts/02_trimming_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/trimming/trimming_stats.log scripts/02_trimming_stats.sh -test
#
# 输入：
#   标准运行: analysis_results/logs/trimming/*_cutadapt.log
#   测试运行: test_analysis_results/logs/trimming/*_cutadapt.log
# 输出：
#   标准运行: analysis_results/02_trimming/trimming_statistics_YYYYMMDD.txt
#   测试运行: test_analysis_results/02_trimming/trimming_statistics_YYYYMMDD.txt
# 依赖：Python3
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
    for tool in python3; do
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
    OUTPUT_DIR=${TEST_TRIMMED_DIR}
    LOG_DIR=${TEST_TRIMMING_LOG}
    log "运行模式: 测试模式"
else
    OUTPUT_DIR=${TRIMMED_DIR}
    LOG_DIR=${TRIMMING_LOG}
    log "运行模式: 标准模式"
fi

log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_DIR}"

# 检查目录是否存在
if [ ! -d "${OUTPUT_DIR}" ]; then
    log "错误: 输出目录不存在: ${OUTPUT_DIR}"
    exit 1
fi

if [ ! -d "${LOG_DIR}" ]; then
    log "错误: 日志目录不存在: ${LOG_DIR}"
    exit 1
fi

# 等待一分钟确保所有文件都已写入完成
log "等待60秒确保所有文件都已完成写入..."
sleep 60

# 检查是否有cutadapt日志文件
LOG_FILES=($(ls ${LOG_DIR}/*_cutadapt.log 2>/dev/null || true))
if [ ${#LOG_FILES[@]} -eq 0 ]; then
    log "错误: ${LOG_DIR}目录中没有找到cutadapt日志文件"
    exit 1
fi

# 生成统计报告
log "生成数据预处理统计报告..."
STATS_FILE="${OUTPUT_DIR}/trimming_statistics_$(date '+%Y%m%d').txt"
echo -e "样本ID\t原始reads数\t过滤后reads数\t保留率" > ${STATS_FILE}

# 处理每个样本的日志文件
for log_file in ${LOG_FILES[@]}; do
    sample_name=$(basename ${log_file} _cutadapt.log)
    
    # 检查输出文件是否存在
    R1_FILE="${OUTPUT_DIR}/${sample_name}_R1_trimmed.fastq.gz"
    R2_FILE="${OUTPUT_DIR}/${sample_name}_R2_trimmed.fastq.gz"
    
    if [ ! -f "${R1_FILE}" ] || [ ! -f "${R2_FILE}" ]; then
        log "警告: 样本 ${sample_name} 的输出文件不存在"
        echo -e "${sample_name}\tNA\tNA\tNA" >> ${STATS_FILE}
        continue
    fi
    
    # 提取统计信息并去除逗号
    total_pairs=$(grep "Total read pairs processed:" ${log_file} | awk '{print $5}' | tr -d ',')
    pairs_written=$(grep "Pairs written" ${log_file} | awk '{print $5}' | tr -d ',')
    
    if [ ! -z "${total_pairs}" ] && [ ! -z "${pairs_written}" ]; then
        # 使用去除逗号后的数字计算保留率
        retention=$(awk "BEGIN {printf \"%.2f\", ${pairs_written}/${total_pairs}*100}")
        # 为了输出格式美观，给数字添加千位分隔符
        total_pairs_formatted=$(printf "%'d" ${total_pairs})
        pairs_written_formatted=$(printf "%'d" ${pairs_written})
        echo -e "${sample_name}\t${total_pairs_formatted}\t${pairs_written_formatted}\t${retention}%" >> ${STATS_FILE}
        log "处理完成: ${sample_name} (保留率: ${retention}%)"
    else
        log "警告: 无法从日志文件中提取统计信息: ${sample_name}"
        echo -e "${sample_name}\tNA\tNA\tNA" >> ${STATS_FILE}
    fi
done

log "统计报告已保存至: ${STATS_FILE}"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 