#!/bin/bash
#$ -S /bin/bash
#$ -N seacr_stats       # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=4G        # 申请4GB内存
#$ -pe smp 1           # 单线程即可

# ===========================================================================
# 脚本名称: 04_peak_calling_stats.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 生成SEACR peak calling的统计报告，包括：
#   - 每个样本的peak数量统计
#   - Peak长度分布分析（最小值、最大值、平均值、中位数）
#   - Bedgraph生成统计信息汇总
#   - 样本间peak重叠分析
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 与peak calling脚本同时提交（推荐）：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/seacr_stats.log -hold_jid seacr_array scripts/04_peak_calling_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/seacr_stats.log -hold_jid seacr_array scripts/04_peak_calling_stats.sh -test
#
# 2. peak calling完成后单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/peak_calling/seacr_stats.log scripts/04_peak_calling_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/peak_calling/seacr_stats.log scripts/04_peak_calling_stats.sh -test
#
# 输入：
#   标准运行:
#     - SEACR peaks: analysis_results/04_peaks/individual_peaks/*.peaks.stringent.bed
#     - bedgraph metrics: analysis_results/04_peaks/bedgraph/metrics/*.bedgraph.metrics
#   测试运行:
#     - SEACR peaks: test_analysis_results/04_peaks/individual_peaks/*.peaks.stringent.bed
#     - bedgraph metrics: test_analysis_results/04_peaks/bedgraph/metrics/*.bedgraph.metrics
#
# 输出：
#   标准运行:
#     - 统计报告: analysis_results/04_peaks/statistics/peak_statistics_YYYYMMDD.txt
#   测试运行:
#     - 统计报告: test_analysis_results/04_peaks/statistics/peak_statistics_YYYYMMDD.txt
#
# 依赖：bedtools, awk
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
    for tool in bedtools awk; do
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
    LOG_OUTPUT_DIR=${TEST_PEAKCALLING_LOG}
    log "运行模式: 测试模式"
else
    OUTPUT_DIR=${PEAKS_DIR}
    LOG_OUTPUT_DIR=${PEAKCALLING_LOG}
    log "运行模式: 标准模式"
fi

# 设置工作目录
PEAK_INPUT_DIR="${OUTPUT_DIR}/individual_peaks"
METRICS_DIR="${OUTPUT_DIR}/bedgraph/metrics"
STATS_DIR="${OUTPUT_DIR}/statistics"

log "输出目录: ${OUTPUT_DIR}"
log "统计目录: ${STATS_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 创建必要的目录
mkdir -p "${STATS_DIR}" "${LOG_OUTPUT_DIR}"

# 等待一分钟确保所有文件都已写入完成
log "等待60秒确保所有文件都已完成写入..."
sleep 60

# 生成统计报告
log "生成统计报告..."
STATS_FILE="${STATS_DIR}/peak_statistics_$(date '+%Y%m%d').txt"

# 创建统计报告头部
cat > ${STATS_FILE} << EOF
SEACR Peak Calling 统计报告
==========================
生成时间: $(date)
运行模式: $([ "$1" == "-test" ] && echo "测试模式" || echo "正式模式")

1. Peak数量统计
==============
EOF

# 统计每个样本的peak数量
log "统计peak数量..."
for peak_file in ${PEAK_INPUT_DIR}/*.stringent.bed; do
    if [ -f "$peak_file" ]; then
        sample_name=$(basename "$peak_file" .stringent.bed)
        peak_count=$(wc -l < "$peak_file")
        echo "${sample_name}: ${peak_count} peaks" >> ${STATS_FILE}
    fi
done

# 添加peak长度分布统计
echo -e "\n2. Peak长度分布统计\n====================" >> ${STATS_FILE}

for peak_file in ${PEAK_INPUT_DIR}/*.stringent.bed; do
    if [ -f "$peak_file" ]; then
        sample_name=$(basename "$peak_file" .stringent.bed)
        echo -e "\n${sample_name}:" >> ${STATS_FILE}
        
        # 计算peak长度统计
        awk -v OFS="\t" '{len=$3-$2; print len}' "$peak_file" | \
        sort -n | \
        awk '
        BEGIN {min=999999999; max=0; sum=0; count=0}
        {
            sum+=$1; count++
            if($1<min) min=$1
            if($1>max) max=$1
            all[count]=$1
        }
        END {
            mean=sum/count
            if(count%2) median=all[int(count/2)+1]
            else median=(all[int(count/2)]+all[int(count/2)+1])/2
            printf "  最小长度: %d bp\n", min
            printf "  最大长度: %d bp\n", max
            printf "  平均长度: %.2f bp\n", mean
            printf "  中位长度: %d bp\n", median
            printf "  总peak数: %d\n", count
        }' >> ${STATS_FILE}
    fi
done

# 添加bedgraph生成统计
echo -e "\n3. Bedgraph生成统计\n=====================" >> ${STATS_FILE}

for metrics_file in ${METRICS_DIR}/*.bedgraph.metrics; do
    if [ -f "$metrics_file" ]; then
        sample_name=$(basename "$metrics_file" .bedgraph.metrics)
        echo -e "\n${sample_name}:" >> ${STATS_FILE}
        cat "$metrics_file" | sed 's/^/  /' >> ${STATS_FILE}
    fi
done

# 添加样本间比较
echo -e "\n4. 样本间Peak重叠分析\n======================" >> ${STATS_FILE}

# 获取所有样本的peak文件
peak_files=(${PEAK_INPUT_DIR}/*.stringent.bed)
n_samples=${#peak_files[@]}

# 如果样本数大于1，进行两两比较
if [ $n_samples -gt 1 ]; then
    echo -e "\n样本间peak重叠数量:" >> ${STATS_FILE}
    for ((i=0; i<n_samples-1; i++)); do
        for ((j=i+1; j<n_samples; j++)); do
            sample1=$(basename "${peak_files[i]}" .stringent.bed)
            sample2=$(basename "${peak_files[j]}" .stringent.bed)
            
            # 计算重叠peaks
            overlap_count=$(bedtools intersect -a "${peak_files[i]}" -b "${peak_files[j]}" -u | wc -l)
            echo "  ${sample1} vs ${sample2}: ${overlap_count} peaks" >> ${STATS_FILE}
        done
    done
fi

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "统计报告已生成: ${STATS_FILE}"
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 