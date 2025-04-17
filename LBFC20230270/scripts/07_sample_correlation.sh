#!/bin/bash
#$ -S /bin/bash
#$ -N correlation_array   # 作业名称
#$ -cwd                  # 使用当前目录
#$ -j y                  # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=16G        # 申请16GB内存/核心
#$ -pe smp 4            # 申请4个CPU线程

# ===========================================================================
# 脚本名称: 07_sample_correlation.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 使用deeptools分析样本间相关性
#   - 使用multiBamSummary计算样本间相关性矩阵
#   - 生成相关性热图和散点图
#   - 提供详细的统计报告
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 标准运行：
#    qsub -o analysis_results/logs/visualization/correlation.\$JOB_ID.log -hold_jid alignment_array scripts/07_sample_correlation.sh [参数]
#
# 2. 测试运行：
#    qsub -o test_analysis_results/logs/visualization/correlation.\$JOB_ID.log scripts/07_sample_correlation.sh -test [参数]
#
# 参数说明：
#   -test          使用测试数据
#   -c METHOD      指定相关性计算方法（可选，默认：pearson）：
#                  - pearson：皮尔逊相关系数
#                  - spearman：斯皮尔曼相关系数
#
# 示例：
#   1. 使用Pearson相关系数（测试模式）：
#      qsub -o test_analysis_results/logs/visualization/correlation.\$JOB_ID.log scripts/07_sample_correlation.sh -test
#
#   2. 使用Spearman相关系数（标准模式）：
#      qsub -o analysis_results/logs/visualization/correlation.\$JOB_ID.log scripts/07_sample_correlation.sh -c spearman
#
# 输入：
#   标准运行:
#     - BAM文件: analysis_results/alignment/filtered/*.filtered.sorted.bam
#   测试运行:
#     - BAM文件: test_analysis_results/alignment/filtered/*.filtered.sorted.bam
#
# 输出：
#   标准运行:
#     - 相关性矩阵: analysis_results/visualization/correlation/*.npz
#     - 相关性热图: analysis_results/visualization/correlation/*.pdf
#     - 统计报告: analysis_results/logs/visualization/correlation_stats_*.txt
#   测试运行:
#     - 相关性矩阵: test_analysis_results/visualization/correlation/*.npz
#     - 相关性热图: test_analysis_results/visualization/correlation/*.pdf
#     - 统计报告: test_analysis_results/logs/visualization/correlation_stats_*.txt
#
# 依赖：deeptools (multiBamSummary, plotCorrelation), samtools
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
    for tool in multiBamSummary plotCorrelation samtools; do
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

# 解析命令行参数
TEST_MODE=false
CORRELATION_METHOD="pearson"

while [[ $# -gt 0 ]]; do
    case $1 in
        -test)
            TEST_MODE=true
            shift
            ;;
        -c)
            CORRELATION_METHOD="$2"
            shift 2
            ;;
        *)
            log "错误：未知参数 $1"
            exit 1
            ;;
    esac
done

# 设置目录
if [ "$TEST_MODE" = true ]; then
    INPUT_DIR=${TEST_ALIGNMENT_DIR}/filtered
    OUTPUT_DIR=${TEST_VISUALIZATION_DIR}/correlation
    LOG_OUTPUT_DIR=${TEST_VISUALIZATION_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${ALIGNMENT_DIR}/filtered
    OUTPUT_DIR=${VISUALIZATION_DIR}/correlation
    LOG_OUTPUT_DIR=${VISUALIZATION_LOG}
    log "运行模式: 标准模式"
fi

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"
log "相关性计算方法: ${CORRELATION_METHOD}"

# 创建必要的目录
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_OUTPUT_DIR}"

# 获取BAM文件列表
mapfile -t BAM_FILES < <(ls ${INPUT_DIR}/*.filtered.sorted.bam 2>/dev/null)
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    log "错误：在${INPUT_DIR}目录下未找到任何过滤后的BAM文件"
    exit 1
fi

# 检查BAM文件索引
log "检查BAM文件索引..."
for bam_file in "${BAM_FILES[@]}"; do
    if [ ! -f "${bam_file}.bai" ]; then
        log "为文件创建索引: ${bam_file}"
        samtools index ${bam_file}
    fi
done

# 构建样本标签列表
SAMPLE_LABELS=""
for bam_file in "${BAM_FILES[@]}"; do
    sample_name=$(basename "${bam_file}" .filtered.sorted.bam)
    SAMPLE_LABELS="${SAMPLE_LABELS} ${sample_name}"
done

# 设置输出文件名
BASE_NAME="correlation_$(date '+%Y%m%d')"
MATRIX_FILE="${OUTPUT_DIR}/${BASE_NAME}_matrix.npz"
HEATMAP_FILE="${OUTPUT_DIR}/${BASE_NAME}_heatmap.pdf"
SCATTERPLOT_FILE="${OUTPUT_DIR}/${BASE_NAME}_scatterplot.pdf"
STATS_FILE="${LOG_OUTPUT_DIR}/correlation_stats_$(date '+%Y%m%d').txt"

# 生成相关性矩阵
log "生成相关性矩阵..."
multiBamSummary bins \
    --bamfiles ${BAM_FILES[@]} \
    --labels ${SAMPLE_LABELS} \
    --minMappingQuality 30 \
    --binSize 10000 \
    --distanceBetweenBins 0 \
    --numberOfProcessors 4 \
    -o ${MATRIX_FILE}

# 生成相关性热图
log "生成相关性热图..."
plotCorrelation \
    --corData ${MATRIX_FILE} \
    --corMethod ${CORRELATION_METHOD} \
    --whatToPlot heatmap \
    --colorMap viridis \
    --plotNumbers \
    --plotTitle "样本间${CORRELATION_METHOD}相关性热图" \
    --plotHeight 10 \
    --plotWidth 10 \
    --plotFile ${HEATMAP_FILE}

# 生成散点图
log "生成相关性散点图..."
plotCorrelation \
    --corData ${MATRIX_FILE} \
    --corMethod ${CORRELATION_METHOD} \
    --whatToPlot scatterplot \
    --plotTitle "样本间${CORRELATION_METHOD}相关性散点图" \
    --plotHeight 10 \
    --plotWidth 10 \
    --plotFile ${SCATTERPLOT_FILE}

# 生成统计报告
log "生成统计报告..."
{
    echo "样本相关性分析统计报告"
    echo "生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    echo "分析参数:"
    echo "  相关性方法: ${CORRELATION_METHOD}"
    echo "  Bin大小: 10000bp"
    echo "  最小比对质量: 30"
    echo ""
    echo "样本信息:"
    echo "  总样本数: ${#BAM_FILES[@]}"
    echo "  样本列表:"
    for bam_file in "${BAM_FILES[@]}"; do
        echo "    - $(basename ${bam_file} .filtered.sorted.bam)"
    done
    echo ""
    echo "输出文件:"
    echo "  相关性矩阵: ${MATRIX_FILE}"
    echo "  相关性热图: ${HEATMAP_FILE}"
    echo "  相关性散点图: ${SCATTERPLOT_FILE}"
} > "${STATS_FILE}"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 