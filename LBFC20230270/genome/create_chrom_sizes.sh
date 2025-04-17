#!/bin/bash

# 生成mm10.chrom.sizes文件
# 作者：lucyyang
# 合作者：Claude AI Assistant
# 创建日期：$(date +"%Y-%m-%d")
# 最后修改：$(date +"%Y-%m-%d")

# 激活conda环境
source activate cuttag

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 设置输入输出文件
GENOME_FA="${SCRIPT_DIR}/mm10.fa"
OUTPUT_FILE="${SCRIPT_DIR}/mm10.chrom.sizes"

# 检查输入文件是否存在
if [ ! -f "${GENOME_FA}" ]; then
    echo "错误：找不到参考基因组文件：${GENOME_FA}"
    conda deactivate
    exit 1
fi

echo "开始生成染色体大小文件..."

# 使用samtools创建索引并生成染色体大小文件
echo "1. 创建基因组索引..."
samtools faidx ${GENOME_FA}

echo "2. 提取染色体名称和大小信息..."
cut -f1,2 ${GENOME_FA}.fai > ${OUTPUT_FILE}

# 验证生成的文件
if [ -f "${OUTPUT_FILE}" ]; then
    echo "3. 验证生成的文件..."
    echo "染色体数量：$(wc -l < ${OUTPUT_FILE})"
    echo "文件内容预览："
    head -n 5 ${OUTPUT_FILE}
    echo "..."
    echo "成功生成文件：${OUTPUT_FILE}"
else
    echo "错误：文件生成失败"
    conda deactivate
    exit 1
fi

conda deactivate 