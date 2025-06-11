#!/bin/bash
#$ -S /bin/bash
#$ -N igv_plot          # 作业名称
#$ -cwd                 # 使用当前目录
#$ -j y                 # 合并标准输出和错误输出
#$ -q all.q             # 队列指定（自动分配节点）
#$ -l h_vmem=4G         # 申请4GB内存/核心
#$ -pe smp 2            # 申请2个CPU线程

# ===========================================================================
# 脚本名称: 08_plot_igv_style.sh
# 作者：Lucy Yang & Claude
# 创建日期：2025-6-10
# 最后修改：$(date +"%Y-%m-%d")
#
# 功能描述: 使用pyGenomeTracks生成IGV风格的基因组浏览器图
#   - 展示特定基因区域的IKZF2蛋白富集模式
#   - 同时显示BigWig信号轨道和peak注释
#   - 生成PNG和PDF格式的高质量图片
# 
# 使用方法: 
# 注意：所有命令都应在项目根目录下运行（~/Project/IKZF2_CUT_TAG/）
#
# 1. 标准运行（生成默认Mthfd1基因图）：
#    qsub -o analysis_results/logs/visualization/igv_plot.log -hold_jid bigwig_array scripts/08_plot_igv_style.sh
#
# 2. 自定义基因和区域：
#    qsub -o analysis_results/logs/visualization/igv_plot.log -hold_jid bigwig_array scripts/08_plot_igv_style.sh -gene Mthfd1 -chr chr6 -start 90976000 -end 90978000
#
# 3. 测试运行：
#    qsub -o test_analysis_results/logs/visualization/igv_plot.log scripts/08_plot_igv_style.sh -test
#
# 输入：
#   标准运行: 
#     - BigWig文件: analysis_results/05_visualization/bigwig/*.bw
#     - Peak文件: analysis_results/04_peaks/*.narrowPeak
#     - GTF文件: genome/mm10.refGene.gtf
#   测试运行:
#     - BigWig文件: test_analysis_results/05_visualization/bigwig/*.bw
#     - Peak文件: test_analysis_results/04_peaks/*.narrowPeak
#
# 输出：
#   标准运行:
#     - PNG图片: analysis_results/05_visualization/igv_plots/{gene}_{chr}_{start}_{end}.png
#     - PDF图片: analysis_results/05_visualization/igv_plots/{gene}_{chr}_{start}_{end}.pdf
#     - 配置文件: analysis_results/05_visualization/igv_plots/{gene}_track_config.ini
#   测试运行:
#     - PNG图片: test_analysis_results/05_visualization/igv_plots/{gene}_{chr}_{start}_{end}.png
#     - PDF图片: test_analysis_results/05_visualization/igv_plots/{gene}_{chr}_{start}_{end}.pdf
#
# 依赖：pyGenomeTracks (需要在cuttag环境中安装)
# ===========================================================================

# 记录开始时间
start_time=$(date +%s)

# 设置错误处理
set -e
trap 'log_error "错误发生在第 $LINENO 行"; exit 1' ERR

# 日志输出函数
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# 错误日志输出函数
log_error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
}

# 获取配置文件路径
CONFIG_FILE="${SGE_O_WORKDIR}/scripts/config/00_config.sh"
if [ -z "${SGE_O_WORKDIR}" ]; then
    log_error "错误: SGE_O_WORKDIR环境变量未设置。请通过qsub提交此脚本。"
    log_error "如果需要直接运行进行测试, 请确保在项目根目录下执行: bash scripts/08_plot_igv_style.sh"
    # 为直接运行提供备用路径
    CONFIG_FILE="./scripts/config/00_config.sh"
fi

if [ ! -f "${CONFIG_FILE}" ]; then
    log_error "错误: 找不到配置文件 ${CONFIG_FILE}"
    exit 1
fi
source ${CONFIG_FILE}

# 检查必要的软件依赖
check_dependencies() {
    local missing_tools=()
    # 检查pyGenomeTracks和python3
    for tool in pyGenomeTracks python3; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "错误: 以下必要工具未安装或未添加到PATH中:"
        for tool in "${missing_tools[@]}"; do
            log_error "  - $tool"
        done
        exit 1
    fi
}

# 激活conda环境
log "激活conda环境: cuttag"
source activate cuttag
check_dependencies

# 设置默认参数
GENE_NAME="Mthfd1"
CHROMOSOME=""
START_POS=""
END_POS=""
TEST_MODE=false

# 处理命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -gene)
            GENE_NAME="$2"
            shift 2
            ;;
        -chr)
            CHROMOSOME="$2"
            shift 2
            ;;
        -start)
            START_POS="$2"
            shift 2
            ;;
        -end)
            END_POS="$2"
            shift 2
            ;;
        -test)
            TEST_MODE=true
            shift
            ;;
        *)
            log "未知参数: $1"
            exit 1
            ;;
    esac
done

# 根据运行模式设置路径
if [ "$TEST_MODE" = true ]; then
    BIGWIG_DIR="${TEST_VISUALIZATION_DIR}/bigwig"
    PEAKS_DIR="${TEST_PEAKS_DIR}"
    OUTPUT_DIR="${TEST_VISUALIZATION_DIR}/igv_plots"
    LOG_OUTPUT_DIR="${TEST_VISUALIZATION_LOG}"
    log "运行模式: 测试模式"
else
    BIGWIG_DIR="${VISUALIZATION_DIR}/bigwig"
    PEAKS_DIR="${PEAKS_DIR}"
    OUTPUT_DIR="${VISUALIZATION_DIR}/igv_plots"
    LOG_OUTPUT_DIR="${VISUALIZATION_LOG}"
    log "运行模式: 标准模式"
fi

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_OUTPUT_DIR}"

log "开始为基因 ${GENE_NAME} 生成IGV风格可视化图..."

# 检查并创建目录，并提取基因信息
setup_and_validate() {
    log "检查输出和日志目录是否存在..."
    mkdir -p "${OUTPUT_DIR}" "${LOG_OUTPUT_DIR}"

    log "正在提取基因 ${GENE_NAME} 的坐标和链信息..."
    # 使用一个健壮的awk命令进行精确匹配
    # 1. IGNORECASE=1: 忽略大小写
    # 2. pattern: 构建一个正则表达式，要求gene_name "..." 后必须是分号或行尾，以确保精确匹配
    # 3. $3 == "transcript": 根据refGene GTF格式，使用'transcript'作为基因级别的条目
    local line_info=$(awk -v gene="$GENE_NAME" '
        BEGIN {
            FS="\t"; OFS=" "; IGNORECASE=1;
            pattern = "gene_name \\\"" gene "\\\"(;|$)";
        }
        $3 == "transcript" && $9 ~ pattern {
            print $1, $4, $5, $7;
        }' ${GTF_FILE} | head -n 1)
    
    if [ -z "${line_info}" ]; then
        log_error "在GTF文件中未找到基因: ${GENE_NAME}"
        exit 1
    fi
    
    # 将基因的真实坐标和链信息读入全局变量
    read -r CHROMOSOME START_POS END_POS STRAND <<< "${line_info}"
    
    log "找到基因坐标: ${CHROMOSOME}:${START_POS}-${END_POS}, 链: ${STRAND}"

    # 根据基因方向，计算仅包含上游10kb的非对称可视化区域
    local vis_start=""
    local vis_end=""
    if [ "${STRAND}" == "+" ]; then
        # 对于正链基因，TSS是START_POS，上游是坐标减小的方向
        vis_start=$((START_POS - 10000))
        vis_end=${END_POS}
        log "基因在正链上，设定可视化区域为 TSS上游10kb 到 基因末端。"
    else
        # 对于负链基因，TSS是END_POS，上游是坐标增加的方向
        vis_start=${START_POS}
        vis_end=$((END_POS + 10000))
        log "基因在负链上，设定可视化区域为 基因开端 到 TSS上游10kb。"
    fi
    
    [ ${vis_start} -lt 0 ] && vis_start=0 # 防止坐标为负

    REGION="${CHROMOSOME}:${vis_start}-${vis_end}"
    log "设定最终可视化区域: ${REGION}"
}

# 提取指定区域内的所有GTF条目
extract_region_gtf() {
    local source_gtf=$1
    local region=$2
    local output_gtf=$3

    log "正在为区域 ${region} 提取GTF条目以供核查..."

    # 解析区域字符串
    local chrom=$(echo "${region}" | cut -d: -f1)
    local range=$(echo "${region}" | cut -d: -f2)
    local start=$(echo "${range}" | cut -d- -f1)
    local end=$(echo "${range}" | cut -d- -f2)

    # 使用 awk 过滤 GTF 文件，提取指定区域内的所有条目
    # $1 == chrom: 染色体匹配
    # $4 < end: 特征的起始位置小于区域的结束位置
    # $5 > start: 特征的结束位置大于区域的起始位置
    # 这确保了任何与该区域有重叠的特征都会被提取
    awk -v chrom="${chrom}" -v start="${start}" -v end="${end}" '
        BEGIN { FS="\t" }
        $1 == chrom && $4 < end && $5 > start {
            print $0
        }
    ' "${source_gtf}" > "${output_gtf}"

    if [ -s "${output_gtf}" ]; then
        log "成功提取并生成区域GTF文件: ${output_gtf}"
    else
        log_warning "未能在区域 ${region} 中提取任何GTF条目。"
    fi
}

# 提取指定区域内所有peak
extract_region_peaks() {
    local source_bed=$1
    local region=$2
    local output_bed=$3

    log "正在为区域 ${region} 提取可见的peak以供核查..."

    if [ ! -f "${source_bed}" ]; then
        log_warning "源peak文件不存在: ${source_bed}。跳过可见peak提取。"
        return
    fi

    local chrom=$(echo "${region}" | cut -d: -f1)
    local range=$(echo "${region}" | cut -d: -f2)
    local start=$(echo "${range}" | cut -d- -f1)
    local end=$(echo "${range}" | cut -d- -f2)

    awk -v chrom="${chrom}" -v start="${start}" -v end="${end}" '
        BEGIN { FS="\t" }
        $1 == chrom && $2 < end && $3 > start {
            print $0
        }
    ' "${source_bed}" > "${output_bed}"

    if [ -s "${output_bed}" ]; then
        log "成功提取并生成可见peak的BED文件: ${output_bed}"
    else
        log_warning "未能在区域 ${region} 中提取任何可见的peak。"
    fi
}

# 新增函数：创建TSS位置的BED文件
create_tss_bed_file() {
    local chrom=$1
    local start=$2
    local end=$3
    local strand=$4
    local output_bed=$5

    local tss_pos=""
    if [ "$strand" == "+" ]; then
        tss_pos=$start
    else
        tss_pos=$end
    fi

    # 创建一个1bp的BED条目代表TSS
    # BED格式是0-based, half-open
    local tss_start=$((tss_pos - 1))
    
    echo -e "${chrom}\t${tss_start}\t${tss_pos}\t${GENE_NAME}_TSS" > "${output_bed}"
    if [ $? -eq 0 ]; then
        log "成功创建TSS BED文件: ${output_bed}"
    else
        log_error "创建TSS BED文件失败"
        exit 1
    fi
}

# 新增函数：创建高亮区域的BED文件
create_highlight_bed() {
    local peaks_bed=$1
    local tss_bed=$2
    local output_bed=$3

    log "正在查找与TSS重叠的peak以创建高亮区域..."

    if ! command -v bedtools &> /dev/null; then
        log_error "bedtools 未安装或不在您的PATH中。无法创建高亮区域。"
        touch "${output_bed}" # 创建空文件以防后续步骤失败
        return 1
    fi
    
    if [ ! -s "${peaks_bed}" ] || [ ! -s "${tss_bed}" ]; then
        log_warning "Peak文件或TSS文件不存在或为空，无法创建高亮区域。"
        touch "${output_bed}"
        return 1
    fi

    # 使用bedtools intersect -wa, -a参数在前，确保我们得到的是a文件中完整的 peak 区域
    bedtools intersect -wa -a "${peaks_bed}" -b "${tss_bed}" > "${output_bed}"

    if [ -s "${output_bed}" ]; then
        log "成功创建高亮区域BED文件: ${output_bed}"
        return 0
    else
        log_warning "未找到与TSS重叠的peak，不创建高亮区域。"
        return 1
    fi
}

# 合并、排序并融合所有的narrowPeak文件
merge_peak_files() {
    local output_merged_bed=$1
    log "开始合并所有的.narrowPeak文件..."

    if [ ! -d "${PEAKS_DIR}" ]; then
        log_warning "Peak目录不存在: ${PEAKS_DIR}。跳过peak轨道。"
        return 1
    fi

    local peak_files=($(find "${PEAKS_DIR}" -name "*.narrowPeak" 2>/dev/null))
    if [ ${#peak_files[@]} -eq 0 ]; then
        log_warning "在 ${PEAKS_DIR} 中未找到.narrowPeak文件。跳过peak轨道。"
        return 1
    fi

    log "找到 ${#peak_files[@]} 个peak文件，正在合并..."
    # 合并所有文件，只取前3列（chr, start, end），排序，然后融合
    cat "${peak_files[@]}" | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' | bedtools sort -i - | bedtools merge -i - > "${output_merged_bed}"

    if [ -s "${output_merged_bed}" ]; then
        log "成功合并并生成融合后的peak文件: ${output_merged_bed}"
        return 0
    else
        log_warning "未能成功生成合并的peak文件。"
        return 1
    fi
}

# 获取所有bigwig文件在指定区域的统一y轴最大值
get_global_max_value() {
    log "正在使用 get_max_signal_value.py 计算所有BigWig文件在区域内的统一Y轴最大值..."

    local helper_script_path=""
    # 根据是否在SGE环境中确定辅助脚本的路径
    if [ -n "${SGE_O_WORKDIR}" ]; then
        helper_script_path="${SGE_O_WORKDIR}/scripts/get_max_signal_value.py"
    else
        # 兼容直接运行: 假设从项目根目录执行 `bash scripts/08...`
        helper_script_path="./scripts/get_max_signal_value.py"
    fi
    
    if [ ! -f "${helper_script_path}" ]; then
        log_error "错误: 未找到辅助脚本 ${helper_script_path}"
        exit 1
    fi

    # 调用Python脚本计算最大值
    local max_val=$(python3 "${helper_script_path}" \
        --bwfiles "${BIGWIG_FILES[@]}" \
        --region "${CHROMOSOME}:${START_POS}-${END_POS}")
    
    if [[ -z "$max_val" || ! "$max_val" =~ ^[0-9]+$ ]]; then
        log_error "错误: get_max_signal_value.py未能返回一个有效的数值。得到的值: '${max_val}'"
        exit 1
    fi

    GLOBAL_MAX_VALUE=${max_val}
    log "所有轨道的统一Y轴最大值 (由Python脚本计算并向上取整): ${GLOBAL_MAX_VALUE}"
}

# 查找BigWig文件
find_bigwig_files() {
    # 检查BigWig文件
    if [ ! -d "${BIGWIG_DIR}" ] || [ -z "$(ls -A ${BIGWIG_DIR}/*.bw 2>/dev/null)" ]; then
        log "错误: BigWig目录不存在或为空: ${BIGWIG_DIR}"
        log "请先运行05_bam_to_bigwig.sh生成BigWig文件"
        exit 1
    fi
    
    # 检查GTF文件
    if [ ! -f "${GTF_FILE}" ]; then
        log "错误: GTF文件不存在: ${GTF_FILE}"
        exit 1
    fi
    
    log "输入文件检查完成"

    # 获取所有BigWig文件列表
    mapfile -t BIGWIG_FILES < <(find ${BIGWIG_DIR} -name "*.bw" | sort)
    if [ ${#BIGWIG_FILES[@]} -eq 0 ]; then
        log_error "错误：未找到任何BigWig文件"
        exit 1
    fi
    log "找到 ${#BIGWIG_FILES[@]} 个BigWig文件"
}

# 生成配置文件
generate_config() {
    local config_file=$1
    local peaks_file=$2
    local highlight_file=$3
    local highlight_status=$4
    local tss_file=$5

    log "==> 进入 generate_config 函数..."
    log "    检查输出文件路径: ${config_file}"
    
    # 检查目录是否存在且可写
    local config_dir=$(dirname "${config_file}")
    if [ ! -d "${config_dir}" ]; then
        log_error "    错误: 配置文件的目录不存在: ${config_dir}"
        return 1
    fi
    if [ ! -w "${config_dir}" ]; then
        log_error "    错误: 配置文件的目录不可写: ${config_dir}"
        return 1
    fi
    log "    目录存在且可写。"
    
    # 检查GTF文件是否存在
    log "    检查GTF文件: ${GTF_FILE}"
    if [ ! -f "${GTF_FILE}" ]; then
        log_error "    错误: GTF文件未找到: ${GTF_FILE}"
        return 1
    fi
    log "    GTF文件存在。"

    log "    执行第一个cat命令，创建配置文件..."
    cat > "${config_file}" << EOF
EOF
    if [ $? -ne 0 ]; then
        log_error "    错误: 第一个cat命令执行失败。"
        return 1
    fi
    log "    第一个cat命令执行成功。"

    # 如果高亮文件成功创建，则添加vhighlight轨道
    if [ "${highlight_status}" -eq 0 ]; then
        log "    添加 vhighlight 轨道..."
        cat >> "${config_file}" << EOF
[vhighlight]
file = ${highlight_file}
type = vhighlight
color = #FFFACD
alpha = 0.4
zorder = -100

EOF
    fi

    log "    添加 BigWig 轨道..."
    # 添加每个BigWig文件的轨道
    for bw_file in "${BIGWIG_FILES[@]}"; do
        local sample_name=$(basename "$bw_file" | sed 's/\.bw$//' | sed 's/\.bigwig$//')
        
        local display_name=""
        local rep_num
        
        # 从样本名中提取数字部分
        rep_num=$(echo "${sample_name}" | grep -o '[0-9]*$')
        
        # 根据样本组（H或I）和数字后缀创建新的显示名称
        if [[ "${sample_name}" == H* ]]; then
            display_name="IKZF2_rep${rep_num}"
        elif [[ "${sample_name}" == I* ]]; then
            display_name="IgG_rep${rep_num}"
        else
            display_name=${sample_name} # 如果不匹配，则使用原始名称
        fi

        local color
        # 根据样本名决定颜色
        if [[ ${sample_name} == H* ]]; then
            color="#d62728"  # 红色
        elif [[ ${sample_name} == I* ]]; then
            color="#808080"  # 灰色
        else
            color="#1f77b4"  # 默认蓝色
        fi

        cat >> "${config_file}" << EOF

[${sample_name}]
file = ${bw_file}
title = ${display_name}
color = ${color}
alpha = 0.8
min_value = 0
max_value = ${GLOBAL_MAX_VALUE}
show_data_range = true
fontsize = 8
type = fill
height = 2
EOF
    done

    # 如果peaks_file存在且不为空，则添加peaks轨道
    if [ -n "$peaks_file" ] && [ -s "$peaks_file" ]; then
        log "    添加 All Peaks 轨道..."
        cat >> "${config_file}" << EOF

[Peaks]
file = ${peaks_file}
title = All Peaks
color = black
display = collapsed
fontsize = 8
height = 0.5

EOF
    fi

    # 最后添加基因轨道和x轴
    log "    添加基因、TSS标记和X轴..."
    cat >> "${config_file}" << EOF

[genes]
file = ${GTF_FILE}
file_type = gtf
height = 1
color = navy
color_arrow = red
arrow_interval = 15
style = UCSC
labels = true
merge_transcripts = true
fontsize = 8
fontstyle = italic
color_backbone = navy
labels_in_margin = true
display = collapsed
gene_rows = 6

[vlines]
file = ${tss_file}
type = vlines
color = #FFA500
line_style = dashed

[x-axis]
where = bottom
fontsize = 8
EOF
    if [ $? -ne 0 ]; then
        log_error "    错误: 添加x-axis失败。"
        return 1
    fi

    log "==> 成功完成 generate_config 函数。"
    return 0
}

# 生成可视化图片
generate_plots() {
    local config_file=$1
    local output_prefix=$2
    local height=$3
    
    # 强制matplotlib使用非交互式后端 (Agg)，这是在SGE等无头环境运行的关键
    export MPLBACKEND=Agg
    unset DISPLAY

    # 记录将要执行的命令
    local cmd="pyGenomeTracks --tracks ${config_file} --region ${REGION} --outFileName ${output_prefix}.png --dpi 300 --trackLabelFraction 0.2 --width 12 --height ${height} --fontSize 8"
    log "准备执行PNG生成命令:"
    log "${cmd}"

    # 生成PNG格式图片
    log "生成PNG格式图片..."
    if pyGenomeTracks --tracks ${config_file} \
                   --region ${REGION} \
                   --outFileName ${output_prefix}.png \
                   --dpi 300 \
                   --trackLabelFraction 0.2 \
                   --width 12 \
                   --height ${height} \
                   --fontSize 8; then
        log "PNG图片生成成功: ${output_prefix}.png"
    else
        log_error "错误：PNG图片生成失败，退出码: $?"
        exit 1
    fi
    
    # 记录将要执行的命令
    cmd="pyGenomeTracks --tracks ${config_file} --region ${REGION} --outFileName ${output_prefix}.pdf --dpi 300 --trackLabelFraction 0.2 --width 12 --height ${height} --fontSize 8"
    log "准备执行PDF生成命令:"
    log "${cmd}"

    # 生成PDF格式图片
    log "生成PDF格式图片..."
    if pyGenomeTracks --tracks ${config_file} \
                   --region ${REGION} \
                   --outFileName ${output_prefix}.pdf \
                   --dpi 300 \
                   --trackLabelFraction 0.2 \
                   --width 12 \
                   --height ${height} \
                   --fontSize 8; then
        log "PDF图片生成成功: ${output_prefix}.pdf"
    else
        log_error "错误：PDF图片生成失败，退出码: $?"
        exit 1
    fi
    
    log "已生成图片: ${output_prefix}.png 和 ${output_prefix}.pdf"
}

# 主函数
main() {
    setup_and_validate
    
    log "步骤 1: 创建TSS位置的BED文件和区域GTF文件..."
    local tss_bed_file="${OUTPUT_DIR}/${GENE_NAME}_TSS.bed"
    create_tss_bed_file "${CHROMOSOME}" "${START_POS}" "${END_POS}" "${STRAND}" "${tss_bed_file}"
    local region_gtf_file="${OUTPUT_DIR}/${GENE_NAME}_region_for_verification.gtf"
    extract_region_gtf "${GTF_FILE}" "${REGION}" "${region_gtf_file}"

    log "步骤 2: 合并所有peak文件..."
    local merged_peaks_file="${OUTPUT_DIR}/all_merged_peaks.bed"
    merge_peak_files "${merged_peaks_file}"
    local merge_success=$?

    log "步骤 3: 提取可视化区域内的Peaks用于核查..."
    local visible_peaks_file="${OUTPUT_DIR}/${GENE_NAME}_visible_peaks.bed"
    if [ ${merge_success} -eq 0 ]; then
        extract_region_peaks "${merged_peaks_file}" "${REGION}" "${visible_peaks_file}"
    fi

    log "步骤 4: 创建与TSS重叠的Peak高亮区域..."
    local highlight_bed_file="${OUTPUT_DIR}/${GENE_NAME}_highlight_peaks.bed"
    create_highlight_bed "${visible_peaks_file}" "${tss_bed_file}" "${highlight_bed_file}"
    local highlight_success=$?

    log "步骤 5: 查找 BigWig 文件..."
    find_bigwig_files
    
    log "步骤 6: 计算所有 BigWig 文件在目标区域的最大信号值..."
    get_global_max_value

    local region_underscore="${REGION//:/_}"
    region_underscore="${region_underscore//-/_}"
    local output_prefix="${OUTPUT_DIR}/${GENE_NAME}_${region_underscore}"
    local config_file="${OUTPUT_DIR}/${GENE_NAME}_config.ini"

    log "步骤 7: 生成配置文件..."
    if [ ${merge_success} -eq 0 ]; then
         generate_config "${config_file}" "${merged_peaks_file}" "${highlight_bed_file}" "${highlight_success}" "${tss_bed_file}"
    else
         generate_config "${config_file}" "" "" 1 "${tss_bed_file}"
    fi

    log "步骤 8: 开始绘制 IGV 风格图..."
    
    # 生成PNG和PDF两种格式，并为height提供默认值
    generate_plots "${config_file}" "${output_prefix}" "10" # 增加了高度以容纳新轨道

    log_final_summary "${output_prefix}" "${region_gtf_file}" "${visible_peaks_file}" "${highlight_bed_file}"
}

# 生成统计报告和最终日志
log_final_summary() {
    local output_prefix=$1
    local region_gtf_path=$2
    local visible_peaks_path=$3
    local highlight_bed_path=$4
    local duration=$((SECONDS - start_time))
    
    local txt_report_path="${LOG_OUTPUT_DIR}/igv_plot_${GENE_NAME}_$(date +%Y%m%d_%H%M%S).txt"
    
    cat > ${txt_report_path} << EOF
IGV风格可视化报告
==================
生成时间: $(date)
基因名称: ${GENE_NAME}
基因组区域: ${CHROMOSOME}:${START_POS}-${END_POS}
运行模式: $([ "$TEST_MODE" = true ] && echo "测试模式" || echo "标准模式")

输入文件:
- GTF文件: ${GTF_FILE}
- BigWig目录: ${BIGWIG_DIR}
- Peak目录: ${PEAKS_DIR}

输出文件:
EOF

    # 列出生成的文件
    for file in ${OUTPUT_DIR}/${GENE_NAME}_${CHROMOSOME}_${START_POS}_${END_POS}*; do
        if [ -f "$file" ]; then
            echo "- $(basename $file): $(ls -lh $file | awk '{print $5}')" >> ${txt_report_path}
        fi
    done
    
    echo "==============================================================================" >> "${txt_report_path}"
    log "                              任务完成                                  "
    log "=============================================================================="
    log "可视化图生成完成！"
    log "输出文件位置："
    log "  - PNG图片: ${output_prefix}.png"
    log "  - PDF图片: ${output_prefix}.pdf"
    log "  - 区域GTF (用于校验): ${region_gtf_path}"
    if [ -s "${visible_peaks_path}" ]; then
        log "  - 可见Peaks (用于校验): ${visible_peaks_path}"
    fi
    if [ -s "${highlight_bed_path}" ]; then
        log "  - 高亮Peaks (用于校验): ${highlight_bed_path}"
    fi
    log "任务完成，运行时间: $(printf '%dh:%dm:%ds\n' $((duration/3600)) $((duration%3600/60)) $((duration%60)))"
}

# 执行主函数
main "$@"

# 生成统计报告
generate_report() {
    local report_file="${LOG_OUTPUT_DIR}/igv_plot_${GENE_NAME}_$(date +%Y%m%d_%H%M%S).txt"
    
    cat > ${report_file} << EOF
IGV风格可视化报告
==================
生成时间: $(date)
基因名称: ${GENE_NAME}
基因组区域: ${CHROMOSOME}:${START_POS}-${END_POS}
运行模式: $([ "$TEST_MODE" = true ] && echo "测试模式" || echo "标准模式")

输入文件:
- GTF文件: ${GTF_FILE}
- BigWig目录: ${BIGWIG_DIR}
- Peak目录: ${PEAKS_DIR}

输出文件:
EOF

    # 列出生成的文件
    for file in ${OUTPUT_DIR}/${GENE_NAME}_${CHROMOSOME}_${START_POS}_${END_POS}*; do
        if [ -f "$file" ]; then
            echo "- $(basename $file): $(ls -lh $file | awk '{print $5}')" >> ${report_file}
        fi
    done
    
    log "统计报告已生成: ${report_file}"
}

generate_report

exit 0 