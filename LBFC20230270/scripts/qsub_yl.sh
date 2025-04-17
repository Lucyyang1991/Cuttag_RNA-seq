# 使用fastqc进行质量控制
qsub -o analysis_results/logs/trimming/trimming_array.\$TASK_ID.log scripts/02_trimming.sh
qsub -o analysis_results/logs/trimming/trimming_stats.log -hold_jid trimming_array scripts/02_trimming_stats.sh
# 使用bowtie2进行alignment
qsub -o analysis_results/logs/alignment/alignment_array.\$TASK_ID.log -hold_jid trimming_array scripts/03_alignment.sh
qsub -o analysis_results/logs/alignment/alignment_stats.log -hold_jid alignment_array scripts/03_alignment_stats.sh

# 使用seacr进行peak calling
qsub -o analysis_results/logs/peak_calling/bedgraph_array.\$TASK_ID.log -hold_jid alignment_array scripts/04_bedgraph_generation.sh
qsub -o analysis_results/logs/peak_calling/seacr_array.\$TASK_ID.log -hold_jid bedgraph_array scripts/04_peak_calling_seacr.sh
qsub -o analysis_results/logs/peak_calling/seacr_stats.log -hold_jid seacr_array scripts/04_peak_calling_stats.sh

# 使用macs3进行peak calling
qsub -o analysis_results/logs/peak_calling/macs3_array.\$TASK_ID.log -hold_jid alignment_array scripts/04_peak_calling_macs3.sh
# 使用idr进行peak calling重复性分析
qsub -o analysis_results/logs/idr/idr_macs3.log scripts/04_idr_analysis.sh

# 使用bam_to_bigwig进行bigwig文件生成
qsub -o analysis_results/logs/visualization/bigwig_array.\$TASK_ID.log -hold_jid alignment_array scripts/05_bam_to_bigwig.sh
# 使用generate_heatmap进行热图生成
qsub -o analysis_results/logs/visualization/heatmap.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r TSS
qsub -o analysis_results/logs/visualization/heatmap.gene_body.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r gene_body
qsub -o analysis_results/logs/visualization/heatmap.DN.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r CD4-CD8-
qsub -o analysis_results/logs/visualization/heatmap.Thymus.\$JOB_ID.log -hold_jid bigwig_array scripts/06_generate_heatmap.sh -r Thymus
# 使用sample_correlation进行样本相关性分析
qsub -o analysis_results/logs/visualization/correlation.\$JOB_ID.log -hold_jid alignment_array scripts/07_sample_correlation.sh

#删除空的文件夹
find ~/Project/IKZF2_CUT_TAG -type d -empty -delete