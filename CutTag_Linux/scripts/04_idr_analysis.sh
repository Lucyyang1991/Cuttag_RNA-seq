#!/bin/bash
#$ -S /bin/bash
#$ -N idr_macs3
#$ -cwd
#$ -j y
#$ -o logs/idr_macs3.o$JOB_ID
#$ -l h_vmem=16G
#$ -pe smp 4

# 作者：Lucy Yang & Claude
# 创建日期：2024-06-22
# 最后修改：2024-06-22

# 脚本功能：使用IDR分析MACS3 peaks的重复性
# 输入：MACS3产生的narrowPeak文件
# 输出：IDR分析结果和可视化
# 依赖：idr工具

# 使用方式：
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 标准运行：
#    qsub -o analysis_results/logs/idr/idr_macs3.log scripts/04_idr_analysis.sh
#
# 注意：IDR工具安装在单独的conda环境中(idr)，不在cuttag环境中。
# 请确保miniconda已安装并且idr环境已创建。
#
# 输入：
#   - MACS3 peaks: analysis_results/04_peaks/individual_peaks/*.macs3_peaks.narrowPeak
#   - 样本信息: scripts/config/sample_info.csv
#
# 输出：
#   - IDR分析结果: analysis_results/04_peaks/idr/*.idr.txt
#   - IDR通过阈值的峰: analysis_results/04_peaks/idr/*_idr_passed.txt
#   - IDR报告: analysis_results/04_peaks/idr/idr_summary.txt
#
# 依赖：idr, awk, bc

# ===== 日志函数 =====
log() {
    local message=$1
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

log_error() {
    local message=$1
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo -e "[ERROR] [${timestamp}] ${message}"
}

log_success() {
    local message=$1
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo -e "[SUCCESS] [${timestamp}] ${message}"
}

# 获取配置文件路径
CONFIG_FILE="${SGE_O_WORKDIR}/scripts/config/00_config.sh"
if [ ! -f "${CONFIG_FILE}" ]; then
    log_error "错误: 找不到配置文件 ${CONFIG_FILE}"
    exit 1
fi
source ${CONFIG_FILE}

# 创建日志目录
mkdir -p "${WORK_DIR}/analysis_results/04_peaks/idr"
mkdir -p "${WORK_DIR}/analysis_results/logs/idr"

# 日志文件
LOG_FILE="${WORK_DIR}/analysis_results/logs/idr/idr_analysis_$(date +%Y%m%d).log"

# 重定向输出到日志文件
exec > >(tee -a ${LOG_FILE}) 2>&1

# 设置输入/输出目录
PEAK_INPUT_DIR="${WORK_DIR}/analysis_results/04_peaks/individual_peaks"
IDR_OUTPUT_DIR="${WORK_DIR}/analysis_results/04_peaks/idr"

log "===== IDR分析开始 ====="
log "脚本路径: $0"
log "当前工作目录: ${SGE_O_WORKDIR}"
log "输入目录: ${PEAK_INPUT_DIR}"
log "输出目录: ${IDR_OUTPUT_DIR}"

# ===== 激活conda环境 =====
activate_conda_env() {
    log "激活conda环境..."
    
    # 检查conda是否可用
    if ! command -v conda &> /dev/null; then
        log_error "未找到conda命令。请确保miniconda或anaconda已安装并添加到PATH中。"
        exit 1
    fi
    
    # 激活idr环境
    log "激活idr环境..."
    
    # 使用CONDA_PREFIX检测是否已在conda环境中
    if [[ -n $CONDA_PREFIX && $(basename "$CONDA_PREFIX") == "idr" ]]; then
        log "已在idr环境中，无需重新激活"
    else
        # 确保conda命令可用于当前shell
        CONDA_BASE=$(conda info --base)
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        
        # 激活idr环境
        conda activate idr
        
        if [ $? -ne 0 ]; then
            log_error "激活idr环境失败。请确保idr环境已创建。"
            log "提示: 可以使用以下命令创建环境: conda create -n idr -c bioconda idr"
            exit 1
        fi
    fi
    
    log_success "conda环境激活成功: $(conda info --envs | grep '*' | awk '{print $1}')"
}

# ===== 检查IDR工具是否安装 =====
check_idr() {
    log "检查IDR工具..."
    if ! command -v idr &> /dev/null; then
        log_error "未找到IDR工具。请确保已在idr环境中正确安装IDR工具。"
        log "提示: 可以使用以下命令安装: conda install -c bioconda idr"
        exit 1
    fi
    
    idr_version=$(idr --version 2>&1 | head -n 1)
    log "IDR版本: ${idr_version}"
    log_success "IDR工具已安装"
}

# ===== 排序所有narrowPeak文件 =====
sort_all_peaks() {
    log "开始对所有峰文件进行排序..."
    
    # 创建排序输出目录
    local sorted_dir="${IDR_OUTPUT_DIR}/sorted_peaks"
    mkdir -p "${sorted_dir}"
    
    # 定义所有需要分析的样本
    local all_samples=("H1_vs_I1" "H1_vs_I2" "H2_vs_I1" "H2_vs_I2" "H3_vs_I1" "H3_vs_I2")
    
    # 对每个样本文件进行排序
    for sample in "${all_samples[@]}"; do
        local input_file="${PEAK_INPUT_DIR}/${sample}.macs3_peaks.narrowPeak"
        local output_file="${sorted_dir}/${sample}.sorted.narrowPeak"
        
        # 检查输入文件是否存在
        if [ ! -f "${input_file}" ]; then
            log_error "排序失败: 输入文件不存在 ${input_file}"
            continue
        fi
        
        log "对文件进行排序: ${input_file}"
        log "排序标准: 第8列 (-log10(p-value))"
        
        # 使用-log10(p-value)，第8列进行排序，降序排列（-rn参数）
        sort -k8,8nr "${input_file}" > "${output_file}"
        
        if [ $? -eq 0 ]; then
            log_success "文件排序完成: ${output_file}"
        else
            log_error "文件排序失败: ${input_file}"
            exit 1  # 排序失败是严重错误，终止脚本执行
        fi
    done
    
    log_success "所有峰文件排序完成"
    return 0
}

# ===== 运行IDR分析 =====
run_idr() {
    local sample1=$1
    local sample2=$2
    local output_prefix=$3
    
    # 使用排序后的峰文件
    local sorted_dir="${IDR_OUTPUT_DIR}/sorted_peaks"
    local sorted_peak1="${sorted_dir}/${sample1}.sorted.narrowPeak"
    local sorted_peak2="${sorted_dir}/${sample2}.sorted.narrowPeak"
    local output="${IDR_OUTPUT_DIR}/${output_prefix}_idr.txt"
    local plot="${IDR_OUTPUT_DIR}/${output_prefix}_idr.png"
    
    # 检查排序后的输入文件
    if [ ! -f "$sorted_peak1" ]; then
        log_error "找不到排序后的文件: $sorted_peak1"
        return 1
    fi
    
    if [ ! -f "$sorted_peak2" ]; then
        log_error "找不到排序后的文件: $sorted_peak2"
        return 1
    fi
    
    log "比较样本 ${sample1} 和 ${sample2}..."
    log "排序后输入文件1: ${sorted_peak1}"
    log "排序后输入文件2: ${sorted_peak2}"
    log "输出文件: ${output}"
    
    # 运行IDR分析，使用排序后的文件
    idr --samples ${sorted_peak1} ${sorted_peak2} \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file ${output} \
        --plot \
        --log-output-file ${IDR_OUTPUT_DIR}/${output_prefix}_idr_log.txt
    
    if [ $? -eq 0 ]; then
        log_success "IDR比较完成: ${output_prefix}"
        
        # IDR阈值0.05对应的分数是540
        # 第5列是缩放后的IDR值: min(int(log2(-125*IDR)), 1000)
        # IDR=0.05 -> score=int(-125*log2(0.05))=540
        # 因此提取第5列值大于等于540的行
        local idr_threshold=540
        local peaks_count=$(awk -v threshold=$idr_threshold '$5 >= threshold' ${output} | wc -l)
        log "通过IDR阈值(0.05，对应分数>=${idr_threshold})的峰数量: ${peaks_count}"
        
        # 保存通过IDR阈值的峰
        awk -v threshold=$idr_threshold '$5 >= threshold' ${output} > ${IDR_OUTPUT_DIR}/${output_prefix}_idr_passed.txt
    else
        log_error "IDR比较失败: ${output_prefix}"
    fi
}

# ===== 生成汇总报告 =====
generate_summary() {
    log "生成IDR分析汇总报告..."
    
    local summary_file="${IDR_OUTPUT_DIR}/idr_summary.txt"
    
    echo "===== MACS3 IDR 分析汇总报告 =====" > ${summary_file}
    echo "生成时间: $(date "+%Y-%m-%d %H:%M:%S")" >> ${summary_file}
    echo "" >> ${summary_file}
    echo "IDR阈值: 0.05 (对应分数值: 540)" >> ${summary_file}
    echo "注: IDR输出文件第5列是缩放后的IDR值，计算公式为min(int(log2(-125*IDR)), 1000)" >> ${summary_file}
    echo "    IDR=0.05对应分数=540，值越高表示峰的可重复性越好" >> ${summary_file}
    echo "" >> ${summary_file}
    echo "| 比较组合 | 总峰数 | 通过IDR阈值的峰数 | 百分比 |" >> ${summary_file}
    echo "|----------|--------|------------------|--------|" >> ${summary_file}
    
    for file in ${IDR_OUTPUT_DIR}/*_idr.txt; do
        local comparison=$(basename ${file} _idr.txt)
        local total=$(wc -l < ${file})
        local passed=$(wc -l < ${IDR_OUTPUT_DIR}/${comparison}_idr_passed.txt)
        local percent=$(echo "scale=2; ${passed}*100/${total}" | bc)
        
        echo "| ${comparison} | ${total} | ${passed} | ${percent}% |" >> ${summary_file}
    done
    
    log_success "汇总报告已生成: ${summary_file}"
}

# ===== 清理环境 =====
cleanup() {
    log "清理环境..."
    
    # 如果在conda环境中，退出环境
    if [[ -n $CONDA_PREFIX ]]; then
        conda deactivate
        log "已退出conda环境"
    fi
    
    log_success "清理完成"
}

# ===== 主函数 =====
main() {
    log "开始IDR分析..."
    
    # 激活conda环境
    activate_conda_env
    
    # 检查IDR是否安装
    check_idr
    
    # 对所有峰文件进行排序
    sort_all_peaks
    
    # 比较不同条件下的重复
    # H1组样本比较
    run_idr "H1_vs_I1" "H1_vs_I2" "H1_replicates"
    
    # H2组样本比较
    run_idr "H2_vs_I1" "H2_vs_I2" "H2_replicates"
    
    # H3组样本比较
    run_idr "H3_vs_I1" "H3_vs_I2" "H3_replicates"
    
    # 跨组比较
    # H1 vs H2
    run_idr "H1_vs_I1" "H2_vs_I1" "H1_H2_vs_I1"
    run_idr "H1_vs_I2" "H2_vs_I2" "H1_H2_vs_I2"
    
    # H1 vs H3
    run_idr "H1_vs_I1" "H3_vs_I1" "H1_H3_vs_I1"
    run_idr "H1_vs_I2" "H3_vs_I2" "H1_H3_vs_I2"
    
    # H2 vs H3
    run_idr "H2_vs_I1" "H3_vs_I1" "H2_H3_vs_I1"
    run_idr "H2_vs_I2" "H3_vs_I2" "H2_H3_vs_I2"
    
    # 生成汇总报告
    generate_summary
    
    # 清理环境
    cleanup
    
    log_success "IDR分析完成"
}

# 执行主函数
main

exit 0 