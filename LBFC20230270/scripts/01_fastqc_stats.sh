#!/bin/bash
#$ -S /bin/bash
#$ -N fastqc_stats      # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=4G        # 申请4GB内存
#$ -pe smp 1           # 单线程即可

# ===========================================================================
# 脚本名称: 01_fastqc_stats.sh
# 作者：Lucy Yang & Claude
# 创建日期：2024-01-20
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 生成FastQC质控结果的汇总报告，包括测序数据质量统计和可视化
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 与FastQC脚本同时提交（推荐）：
#    标准运行: 
#      qsub -o analysis_results/logs/fastqc/fastqc_stats.log -hold_jid fastqc_array scripts/01_fastqc_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/fastqc/fastqc_stats.log -hold_jid fastqc_array scripts/01_fastqc_stats.sh -test
#
# 2. FastQC完成后单独运行：
#    标准运行: 
#      qsub -o analysis_results/logs/fastqc/fastqc_stats.log scripts/01_fastqc_stats.sh
#    测试运行: 
#      qsub -o test_analysis_results/logs/fastqc/fastqc_stats.log scripts/01_fastqc_stats.sh -test
#
# 输入：
#   标准运行: analysis_results/01_fastqc/raw/*.zip
#   测试运行: test_analysis_results/01_fastqc/raw/*.zip
# 输出：
#   标准运行: analysis_results/01_fastqc/raw/multiqc_report.html
#   测试运行: test_analysis_results/01_fastqc/raw/multiqc_report.html
# 依赖：multiqc
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
    for tool in multiqc; do
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
    FASTQC_DIR="${TEST_FASTQC_DIR}/raw"
    LOG_DIR="${TEST_FASTQC_LOG}"
    log "运行模式: 测试模式"
else
    FASTQC_DIR="${FASTQC_DIR}/raw"
    LOG_DIR="${FASTQC_LOG}"
    log "运行模式: 标准模式"
fi

log "FastQC结果目录: ${FASTQC_DIR}"
log "日志输出目录: ${LOG_DIR}"

# 创建日志目录
mkdir -p ${LOG_DIR}

# 检查FastQC输出目录是否存在
if [ ! -d "${FASTQC_DIR}" ]; then
    log "错误: FastQC输出目录不存在: ${FASTQC_DIR}"
    exit 1
fi

# 检查是否有FastQC结果文件
FASTQC_FILES=($(ls ${FASTQC_DIR}/*.zip 2>/dev/null || true))
if [ ${#FASTQC_FILES[@]} -eq 0 ]; then
    log "错误: ${FASTQC_DIR}目录中没有找到FastQC结果文件(*.zip)"
    exit 1
fi

# 等待一分钟确保所有文件都已写入完成
log "等待60秒确保所有文件都已完成写入..."
sleep 60

# 运行MultiQC生成汇总报告
log "运行MultiQC生成汇总报告..."
multiqc \
    -f \
    -o ${FASTQC_DIR} \
    -n "multiqc_report.html" \
    --title "FastQC Quality Report" \
    ${FASTQC_DIR}

# 检查MultiQC报告是否生成
if [ ! -f "${FASTQC_DIR}/multiqc_report.html" ]; then
    log "错误: MultiQC报告生成失败"
    exit 1
fi

log "MultiQC报告已生成: ${FASTQC_DIR}/multiqc_report.html"

# 计算运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"

exit 0 