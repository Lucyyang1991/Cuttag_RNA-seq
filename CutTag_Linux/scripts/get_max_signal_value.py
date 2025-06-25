#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
计算一个或多个bigWig文件中给定区域的最大信号值。
"""

import argparse
import pyBigWig
import numpy as np
import sys

def get_max_value(bigwig_files, region):
    """
    遍历所有bigWig文件，找到指定区域中的最大信号值。

    Args:
        bigwig_files (list): bigWig文件路径列表。
        region (str): 格式为 'chr:start-end' 的基因组区域。

    Returns:
        float: 所有文件在该区域内的最大信号值。
    """
    try:
        chrom, coords = region.split(':')
        start, end = map(int, coords.split('-'))
    except ValueError:
        print(f"错误: 区域格式不正确 '{region}'. 请使用 'chr:start-end' 格式。", file=sys.stderr)
        sys.exit(1)

    overall_max = 0.0

    for bw_file in bigwig_files:
        try:
            with pyBigWig.open(bw_file) as bw:
                # 确保染色体存在于文件中
                if chrom not in bw.chroms():
                    # print(f"警告: 在文件 {bw_file} 中找不到染色体 {chrom}，跳过。", file=sys.stderr)
                    continue
                
                # 获取区域内的最大值，如果区域无效或无信号，.stats返回None
                region_max = bw.stats(chrom, start, end, type='max')
                
                # region_max 是一个列表，我们取第一个元素
                if region_max and region_max[0] is not None:
                    if region_max[0] > overall_max:
                        overall_max = region_max[0]

        except Exception as e:
            print(f"错误: 处理文件 {bw_file} 时出错: {e}", file=sys.stderr)
            continue
            
    return overall_max

def main():
    parser = argparse.ArgumentParser(description="从多个bigWig文件中计算给定区域的最大信号值。")
    parser.add_argument('--bwfiles', required=True, nargs='+', help="一个或多个bigWig文件的路径，用空格分隔。")
    parser.add_argument('--region', required=True, type=str, help="格式为 'chr:start-end' 的基因组区域。")
    
    args = parser.parse_args()
    
    max_value = get_max_value(args.bwfiles, args.region)
    
    # 将最终的最大值向上取整后打印到标准输出
    print(int(np.ceil(max_value)))

if __name__ == '__main__':
    main() 