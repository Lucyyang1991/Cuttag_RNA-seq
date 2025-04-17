#!/bin/bash
#$ -S /bin/bash
#$ -N heatmap_array    # 作业名称
#$ -cwd                # 使用当前目录
#$ -j y                # 合并标准输出和错误输出
#$ -q all.q            # 队列指定（自动分配节点）
#$ -l h_vmem=16G      # 申请16GB内存/核心
#$ -pe smp 4          # 申请4个CPU线程

# ===========================================================================
# 脚本名称: 06_generate_heatmap.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 使用deeptools生成TSS区域和其他感兴趣区域的信号热图
#   - 支持多种区域类型：TSS、基因体、增强子
#   - 生成热图和信号曲线图
#   - 提供详细的统计报告
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 标准运行：
#    qsub -o analysis_results/logs/visualization/heatmap.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r TSS
#    qsub -o analysis_results/logs/visualization/heatmap.gene_body.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r gene_body
#    qsub -o analysis_results/logs/visualization/heatmap.DN.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r CD4-CD8-
#    qsub -o analysis_results/logs/visualization/heatmap.Thymus.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r Thymus
#
# 2. 测试运行：
#    qsub -o test_analysis_results/logs/visualization/heatmap.\$JOB_ID.log scripts/06_generate_heatmap.sh -test [参数]
#
# 参数说明：
#   -test          使用测试数据
#   -r REGION      指定区域类型（必需）：
#                  - TSS：转录起始位点区域
#                  - gene_body：基因体区域
#                  - CD4-CD8-：CD4-CD8-增强子区域
#                  - Thymus：胸腺增强子区域
#
# 示例：
#   1. 分析TSS区域（测试模式）：
#      qsub -o test_analysis_results/logs/visualization/heatmap.\$JOB_ID.log scripts/06_generate_heatmap.sh -test -r TSS
#
#   2. 分析基因体区域（标准模式）：
#      qsub -o analysis_results/logs/visualization/heatmap.\$JOB_ID.log scripts/06_generate_heatmap.sh -r gene_body
#
# 区域类型说明：
#   1. TSS区域：
#      - 以转录起始位点为中心，上下游各5000bp
#      - 参考点为TSS
#      - 适用于分析转录因子在转录起始位点附近的富集情况
#
#   2. 基因体：
#      - 从转录起始位点到转录终止位点
#      - 将基因体缩放为统一长度（5000bp）
#      - 适用于分析转录因子在整个基因区域的分布
#
#   3. CD4-CD8-增强子：
#      - 使用CD4-CD8-特异的增强子区域
#      - 以增强子中心为参考点，上下游各5000bp
#      - 文件：genome/modified/CD4-CD8-_modified.bed
#
#   4. Thymus增强子：
#      - 使用胸腺特异的增强子区域
#      - 以增强子中心为参考点，上下游各5000bp
#      - 文件：genome/modified/Thymus_modified.bed
#
# 输入：
#   标准运行:
#     - BigWig文件: analysis_results/visualization/bigwig/*.bw
#     - 基因注释文件: genome/mm10.refGene.gtf
#   测试运行:
#     - BigWig文件: test_analysis_results/visualization/bigwig/*.bw
#     - 基因注释文件: genome/mm10.refGene.gtf
#
# 输出：
#   标准运行:
#     - 热图矩阵: analysis_results/visualization/heatmaps/*.gz
#     - 热图: analysis_results/visualization/heatmaps/*.pdf
#     - 信号曲线图: analysis_results/visualization/heatmaps/*.pdf
#   测试运行:
#     - 热图矩阵: test_analysis_results/visualization/heatmaps/*.gz
#     - 热图: test_analysis_results/visualization/heatmaps/*.pdf
#     - 信号曲线图: test_analysis_results/visualization/heatmaps/*.pdf
#
# 依赖：deeptools, parallel
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
    for tool in computeMatrix plotHeatmap plotProfile parallel; do
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
REGION_TYPE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -test)
            TEST_MODE=true
            shift
            ;;
        -r)
            REGION_TYPE="$2"
            shift 2
            ;;
        *)
            log "错误：未知参数 $1"
            exit 1
            ;;
    esac
done

# 检查是否提供了区域类型
if [ -z "${REGION_TYPE}" ]; then
    log "错误：必须使用 -r 参数指定区域类型（TSS/gene_body/CD4-CD8-/Thymus）"
    exit 1
fi

# 设置目录
if [ "$TEST_MODE" = true ]; then
    INPUT_DIR=${TEST_VISUALIZATION_DIR}/bigwig
    OUTPUT_DIR=${TEST_VISUALIZATION_DIR}/heatmaps
    LOG_OUTPUT_DIR=${TEST_VISUALIZATION_LOG}
    log "运行模式: 测试模式"
else
    INPUT_DIR=${VISUALIZATION_DIR}/bigwig
    OUTPUT_DIR=${VISUALIZATION_DIR}/heatmaps
    LOG_OUTPUT_DIR=${VISUALIZATION_LOG}
    log "运行模式: 标准模式"
fi

log "输入目录: ${INPUT_DIR}"
log "输出目录: ${OUTPUT_DIR}"
log "日志目录: ${LOG_OUTPUT_DIR}"

# 创建必要的目录
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_OUTPUT_DIR}"

# 根据区域类型设置参数
log "处理区域类型: ${REGION_TYPE}"
case ${REGION_TYPE} in
    TSS)
        REFERENCE_POINT="TSS"
        REGION_LABEL="Transcription Start Site (TSS)"
        MATRIX_TYPE="reference-point"
        COMMAND="computeMatrix reference-point"
        GTF_FILE="${GENOME_DIR}/mm10.refGene.gtf"
        ;;
    gene_body) # 基因体区域，时间最长
        REFERENCE_POINT=""
        REGION_LABEL="Gene Body"
        MATRIX_TYPE="scale-regions"
        COMMAND="computeMatrix scale-regions"
        GTF_FILE="${GENOME_DIR}/mm10.refGene.gtf"
        ;;
    CD4-CD8-)
        ENHANCER_FILE="${GENOME_DIR}/modified/CD4-CD8-_modified.bed"
        REFERENCE_POINT="center"
        REGION_LABEL="CD4-CD8- Enhancers"
        MATRIX_TYPE="reference-point"
        COMMAND="computeMatrix reference-point"
        GTF_FILE="${ENHANCER_FILE}"
        ;;
    Thymus)
        ENHANCER_FILE="${GENOME_DIR}/modified/Thymus_modified.bed"
        REFERENCE_POINT="center"
        REGION_LABEL="Thymus Enhancers"
        MATRIX_TYPE="reference-point"
        COMMAND="computeMatrix reference-point"
        GTF_FILE="${ENHANCER_FILE}"
        ;;
    *)
        log "错误：未知的区域类型 ${REGION_TYPE}"
        log "可用的区域类型: TSS, gene_body, CD4-CD8-, Thymus"
        exit 1
        ;;
esac

# 检查输入文件
if [ ! -f "${GTF_FILE}" ]; then
    log "错误: 找不到参考文件 ${GTF_FILE}"
    exit 1
fi

# 获取BigWig文件列表
mapfile -t BIGWIG_FILES < <(ls ${INPUT_DIR}/*.bw 2>/dev/null)
if [ ${#BIGWIG_FILES[@]} -eq 0 ]; then
    log "错误：在${INPUT_DIR}目录下未找到任何.bw文件"
    exit 1
fi

# 构建参数字符串
BIGWIG_LIST=""
SAMPLE_LABELS=""
for bw_file in "${BIGWIG_FILES[@]}"; do
    BIGWIG_LIST="${BIGWIG_LIST} ${bw_file}"
    sample_name=$(basename "${bw_file}" .bw)
    SAMPLE_LABELS="${SAMPLE_LABELS} ${sample_name}"
done

# 设置输出文件名
BASE_NAME="${REGION_TYPE}_$(date '+%Y%m%d')"
MATRIX_FILE="${OUTPUT_DIR}/${BASE_NAME}_matrix.gz"
HEATMAP_FILE="${OUTPUT_DIR}/${BASE_NAME}_heatmap.pdf"
PROFILE_FILE="${OUTPUT_DIR}/${BASE_NAME}_profile.pdf"

# 生成热图矩阵
log "生成热图矩阵..."
if [ "${MATRIX_TYPE}" == "reference-point" ]; then
    ${COMMAND} \
        --referencePoint ${REFERENCE_POINT} \
        -S ${BIGWIG_LIST} \
        -R ${GTF_FILE} \
        -b ${HEATMAP_UPSTREAM} -a ${HEATMAP_DOWNSTREAM} \
        --skipZeros \
        --sortRegions descend \
        --binSize ${HEATMAP_BIN_SIZE} \
        --missingDataAsZero \
        --numberOfProcessors 4 \
        -o ${MATRIX_FILE}
else
    ${COMMAND} \
        -S ${BIGWIG_LIST} \
        -R ${GTF_FILE} \
        -b ${HEATMAP_UPSTREAM} -a ${HEATMAP_DOWNSTREAM} \
        --skipZeros \
        --sortRegions descend \
        --binSize ${HEATMAP_BIN_SIZE} \
        --regionBodyLength ${HEATMAP_BODY_LENGTH} \
        --missingDataAsZero \
        --numberOfProcessors 4 \
        -o ${MATRIX_FILE}
fi

# 生成热图
log "生成热图..."
plotHeatmap \
    -m ${MATRIX_FILE} \
    --colorMap RdYlBu_r \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel "${REGION_LABEL}" \
    --zMin 0 \
    --yAxisLabel "Signal" \
    --heatmapHeight 15 \
    --heatmapWidth 8 \
    --dpi 300 \
    --plotTitle "CUT&Tag Signal at ${REGION_LABEL}" \
    -o ${HEATMAP_FILE}

# 生成信号曲线图
log "生成信号曲线图..."
plotProfile \
    -m ${MATRIX_FILE} \
    --perGroup \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel "${REGION_LABEL}" \
    --plotHeight 10 \
    --plotWidth 10 \
    --dpi 300 \
    --plotTitle "CUT&Tag Signal Profile at ${REGION_LABEL}" \
    -o ${PROFILE_FILE}

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0