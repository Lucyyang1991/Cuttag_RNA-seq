#!/bin/bash

# 脚本功能：从原始fastq文件中随机抽取1万条reads用于测试
# 输入：merge.file目录下的fastq.gz文件
# 输出：test_samples目录下的抽样后的fastq.gz文件
# 依赖：seqtk工具
# 环境：服务器端conda环境cuttag
# 执行位置：服务器端

# 使用方式：
# sh 00_subsample.sh [参数]
# 
# 示例：
#   sh 00_subsample.sh                  # 从merge.file目录抽取样本到test_samples目录
#   sh 00_subsample.sh 20000            # 指定抽取20000条reads (默认10000)

# 激活conda环境
source activate cuttag

# 设置错误时退出
set -e

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${SCRIPT_DIR}

# 导出配置变量
export MERGE_DIR="${SCRIPT_DIR}/merge.file"
export TEST_DIR="${SCRIPT_DIR}/test_samples"
export SAMPLE_SIZE=10000
export LOG_DIR="${SCRIPT_DIR}/logs/subsample"

# 创建必要的目录
mkdir -p ${TEST_DIR}
mkdir -p ${LOG_DIR}

# 开始日志记录
MASTER_LOG="${LOG_DIR}/master_$(date '+%Y%m%d_%H%M%S').log"
echo "[$(date)] 开始执行抽样脚本" | tee ${MASTER_LOG}
echo "[$(date)] 使用的conda环境: $(conda info --envs | grep '*' | awk '{print $1}')" | tee -a ${MASTER_LOG}
echo "[$(date)] seqtk版本: $(seqtk 2>&1 | grep Version || echo '无法获取版本信息')" | tee -a ${MASTER_LOG}

# 获取样本名列表（不包括.R1.fastq.gz和.R2.fastq.gz后缀）
SAMPLES=($(ls ${MERGE_DIR}/*.R1.fastq.gz | sed 's/.*\///' | sed 's/\.R1\.fastq\.gz//'))
echo "[$(date)] 检测到以下样本: ${SAMPLES[@]}" | tee -a ${MASTER_LOG}

# 定义抽样函数
subsample_fastq() {
    local sample=$1
    local timestamp=$(date "+%Y%m%d_%H%M%S")
    local log_file="${LOG_DIR}/${sample}_${timestamp}.log"
    
    echo "[$(date)] 开始处理样本: ${sample}" | tee -a ${log_file}
    
    # 检查输入文件是否存在
    if [ ! -f "${MERGE_DIR}/${sample}.R1.fastq.gz" ] || [ ! -f "${MERGE_DIR}/${sample}.R2.fastq.gz" ]; then
        echo "[$(date)] 错误: 输入文件不存在" | tee -a ${log_file}
        return 1
    fi
    
    # 处理R1
    echo "[$(date)] 抽样处理 ${sample}.R1.fastq.gz" | tee -a ${log_file}
    seqtk sample -s100 ${MERGE_DIR}/${sample}.R1.fastq.gz ${SAMPLE_SIZE} | gzip > ${TEST_DIR}/${sample}.R1.fastq.gz 2>> ${log_file}
    
    # 处理R2
    echo "[$(date)] 抽样处理 ${sample}.R2.fastq.gz" | tee -a ${log_file}
    seqtk sample -s100 ${MERGE_DIR}/${sample}.R2.fastq.gz ${SAMPLE_SIZE} | gzip > ${TEST_DIR}/${sample}.R2.fastq.gz 2>> ${log_file}
    
    # 验证文件大小
    if [ -s "${TEST_DIR}/${sample}.R1.fastq.gz" ] && [ -s "${TEST_DIR}/${sample}.R2.fastq.gz" ]; then
        echo "[$(date)] ${sample} 抽样完成" | tee -a ${log_file}
    else
        echo "[$(date)] 错误: ${sample} 抽样失败" | tee -a ${log_file}
        return 1
    fi
    
    # 统计reads数量
    echo "[$(date)] 验证抽样数量:" | tee -a ${log_file}
    echo "[$(date)] 原始R1 reads数: $(zcat ${MERGE_DIR}/${sample}.R1.fastq.gz | wc -l | awk '{print $1/4}')" | tee -a ${log_file}
    echo "[$(date)] 抽样R1 reads数: $(zcat ${TEST_DIR}/${sample}.R1.fastq.gz | wc -l | awk '{print $1/4}')" | tee -a ${log_file}
    echo "[$(date)] 原始R2 reads数: $(zcat ${MERGE_DIR}/${sample}.R2.fastq.gz | wc -l | awk '{print $1/4}')" | tee -a ${log_file}
    echo "[$(date)] 抽样R2 reads数: $(zcat ${TEST_DIR}/${sample}.R2.fastq.gz | wc -l | awk '{print $1/4}')" | tee -a ${log_file}
}

# 导出函数供parallel使用
export -f subsample_fastq
export MERGE_DIR
export TEST_DIR
export SAMPLE_SIZE
export LOG_DIR

# 并行处理所有样本
echo "[$(date)] 开始并行处理样本..." | tee -a ${MASTER_LOG}
parallel --will-cite -j 5 subsample_fastq ::: "${SAMPLES[@]}"

# 检查是否所有样本都处理成功
EXPECTED_COUNT=${#SAMPLES[@]}
ACTUAL_COUNT=$(ls ${TEST_DIR}/*.R1.fastq.gz 2>/dev/null | wc -l)
if [ "$ACTUAL_COUNT" -eq "$EXPECTED_COUNT" ]; then
    echo "[$(date)] 所有样本处理成功: $ACTUAL_COUNT/$EXPECTED_COUNT" | tee -a ${MASTER_LOG}
else
    echo "[$(date)] 警告: 部分样本处理失败: $ACTUAL_COUNT/$EXPECTED_COUNT" | tee -a ${MASTER_LOG}
fi

# 生成抽样报告
echo "[$(date)] 生成抽样报告..." | tee -a ${MASTER_LOG}
REPORT_FILE="${LOG_DIR}/subsample_report_$(date '+%Y%m%d').txt"
echo "Cut&Tag数据抽样统计报告" > ${REPORT_FILE}
echo "=========================" >> ${REPORT_FILE}
echo "生成时间: $(date)" >> ${REPORT_FILE}
echo "服务器: $(hostname)" >> ${REPORT_FILE}
echo "执行用户: $(whoami)" >> ${REPORT_FILE}
echo "Conda环境: $(conda info --envs | grep '*' | awk '{print $1}')" >> ${REPORT_FILE}
echo "=========================" >> ${REPORT_FILE}
echo "抽样数量: ${SAMPLE_SIZE} reads" >> ${REPORT_FILE}
echo "处理样本数: ${#SAMPLES[@]}" >> ${REPORT_FILE}
echo "=========================" >> ${REPORT_FILE}
echo "原始文件信息:" >> ${REPORT_FILE}
for sample in "${SAMPLES[@]}"; do
    echo "$sample:" >> ${REPORT_FILE}
    echo "  R1大小: $(du -h ${MERGE_DIR}/${sample}.R1.fastq.gz | cut -f1)" >> ${REPORT_FILE}
    echo "  R2大小: $(du -h ${MERGE_DIR}/${sample}.R2.fastq.gz | cut -f1)" >> ${REPORT_FILE}
done
echo "=========================" >> ${REPORT_FILE}
echo "抽样文件信息:" >> ${REPORT_FILE}
for sample in "${SAMPLES[@]}"; do
    if [ -f "${TEST_DIR}/${sample}.R1.fastq.gz" ]; then
        echo "$sample:" >> ${REPORT_FILE}
        echo "  R1大小: $(du -h ${TEST_DIR}/${sample}.R1.fastq.gz | cut -f1)" >> ${REPORT_FILE}
        echo "  R2大小: $(du -h ${TEST_DIR}/${sample}.R2.fastq.gz | cut -f1)" >> ${REPORT_FILE}
    else
        echo "$sample: 处理失败" >> ${REPORT_FILE}
    fi
done

echo "[$(date)] 抽样完成！请查看 ${REPORT_FILE} 获取详细信息。" | tee -a ${MASTER_LOG} 