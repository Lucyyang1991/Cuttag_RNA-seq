#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
处理增强子文件，添加标识符和来源信息
数据来源：EnhancerAtlas 2.0 (http://www.enhanceratlas.org/downloadv2.php)
作者：Claude AI Assistant
创建日期：2024-04-07
"""

import os

def process_enhancer_file(input_file, output_file, prefix):
    """处理增强子文件，添加标识符"""
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # 添加文件头注释
        fout.write('# 数据来源: EnhancerAtlas 2.0 (http://www.enhanceratlas.org/downloadv2.php)\n')
        fout.write('# 格式: chr\tstart\tend\tname\tscore\n')
        
        # 处理每一行
        for i, line in enumerate(fin, 1):
            if line.startswith('#'):
                continue
            
            # 解析行
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
                
            chr_name = parts[0]
            start = parts[1]
            end = parts[2]
            score = parts[3]
            
            # 创建新的标识符
            enhancer_id = f"{prefix}_{str(i).zfill(6)}"
            
            # 写入新行
            fout.write(f"{chr_name}\t{start}\t{end}\t{enhancer_id}\t{score}\n")

def main():
    """主函数"""
    # 设置输入输出路径
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # CD4-CD8-增强子
    cd4cd8_input = os.path.join(base_dir, "CD4-CD8-.bed")
    cd4cd8_output = os.path.join(base_dir, "modified", "CD4-CD8-_modified.bed")
    process_enhancer_file(cd4cd8_input, cd4cd8_output, "CD4CD8N_Enh")
    
    # Thymus增强子
    thymus_input = os.path.join(base_dir, "Thymus.bed")
    thymus_output = os.path.join(base_dir, "modified", "Thymus_modified.bed")
    process_enhancer_file(thymus_input, thymus_output, "Thymus_Enh")

if __name__ == "__main__":
    main() 