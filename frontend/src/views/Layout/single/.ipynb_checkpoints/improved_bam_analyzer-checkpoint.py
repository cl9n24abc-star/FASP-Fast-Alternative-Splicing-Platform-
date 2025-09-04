#!/usr/bin/env python3
"""
BAM分析器测试调用脚本
用于测试改进后的BAM分析器功能
"""

import sys
import os
import tempfile
from pathlib import Path

# 添加当前目录到Python路径，以便导入模块
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

def create_test_bam_list(bam_files, filename):
    """创建测试用的BAM文件列表"""
    with open(filename, 'w') as f:
        for bam_file in bam_files:
            f.write(f"{bam_file}\n")
    print(f"Created BAM list file: {filename}")

def test_bam_analyzer():
    """测试BAM分析器"""
    print("=" * 60)
    print("BAM分析器测试")
    print("=" * 60)
    
    # 你的BAM文件路径
    sample1_bams = ["/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam"]
    sample2_bams = ["/mnt/disk1/rna_seq_analysis/step0/SRR12125134_Aligned.sortedByCoord.out.bam"]
    
    # 创建临时的BAM列表文件
    sample1_list = "sample1_test.txt"
    sample2_list = "sample2_test.txt"
    
    create_test_bam_list(sample1_bams, sample1_list)
    create_test_bam_list(sample2_bams, sample2_list)
    
    # 输出文件
    output_file = "bam_analysis_test_result.json"
    
    try:
        print("\n1. 验证BAM文件是否存在...")
        for bam_file in sample1_bams + sample2_bams:
            if os.path.exists(bam_file):
                print(f"  ✓ {bam_file}")
            else:
                print(f"  ✗ {bam_file} - 文件不存在")
                return
        
        print("\n2. 创建BAM分析器实例...")
        # 使用较小的区域大小进行快速测试 (50kb)
        analyzer = BAMAnalyzer(
            sample1_bams=sample1_bams,
            sample2_bams=sample2_bams, 
            output_file=output_file,
            region_size=50000  # 50kb per chromosome
        )
        
        print("\n3. 运行分析...")
        result = analyzer.run_analysis()
        
        print("\n4. 分析结果摘要:")
        print(f"  - 输出文件: {output_file}")
        print(f"  - 分析时间: {result['analysis_info']['analysis_time_seconds']} 秒")
        print(f"  - 分析染色体数: {result['analysis_info']['total_chromosomes_analyzed']}")
        print(f"  - 每染色体区域大小: {result['analysis_info']['region_size_per_chromosome']} bp")
        
        print("\n5. 样本统计:")
        print(f"  Group 1:")
        print(f"    - 文件大小: {result['group1']['fileSize']}")
        print(f"    - 映射率: {result['group1']['mappingRate']}%")
        print(f"    - 平均MAPQ: {result['group1']['avgMapQ']}")
        print(f"    - 平均插入大小: {result['group1']['avgInsertSize']} bp")
        
        print(f"  Group 2:")
        print(f"    - 文件大小: {result['group2']['fileSize']}")
        print(f"    - 映射率: {result['group2']['mappingRate']}%")
        print(f"    - 平均MAPQ: {result['group2']['avgMapQ']}")
        print(f"    - 平均插入大小: {result['group2']['avgInsertSize']} bp")
        
        print("\n6. 图表数据验证:")
        chart_data = result['chart_data']
        
        # 饼图数据
        pie_data = chart_data['pie_chart']
        print(f"  饼图数据: {len(pie_data)} 个类别")
        for item in pie_data:
            print(f"    - {item['name']}: {item['value']}%")
        
        # MAPQ直方图
        mapq_data = chart_data['mapq_histogram']
        print(f"  MAPQ直方图: {len(mapq_data)} 个区间")
        
        # 热力图数据
        heatmap_data = chart_data['chromosome_heatmap']
        print(f"  染色体热力图: {len(heatmap_data)} 个数据点")
        
        # 统计每组的覆盖度
        group1_coverage = [d[2] for d in heatmap_data if d[1] == 0]
        group2_coverage = [d[2] for d in heatmap_data if d[1] == 1]
        
        if group1_coverage:
            print(f"    Group 1 覆盖度: 平均={sum(group1_coverage)/len(group1_coverage):.1f}%, 范围={min(group1_coverage):.1f}-{max(group1_coverage):.1f}%")
        if group2_coverage:
            print(f"    Group 2 覆盖度: 平均={sum(group2_coverage)/len(group2_coverage):.1f}%, 范围={min(group2_coverage):.1f}-{max(group2_coverage):.1f}%")
        
        # 火山图数据
        volcano_data = chart_data['volcano_plot']
        print(f"  火山图: {len(volcano_data)} 个基因点")
        
        print(f"\n✓ 分析完成! 结果已保存到 {output_file}")
        print(f"  你可以将此文件放到前端public目录供Vue组件使用")
        
    except Exception as e:
        print(f"\n✗ 分析失败: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        # 清理临时文件
        for temp_file in [sample1_list, sample2_list]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"清理临时文件: {temp_file}")

def quick_test():
    """快速测试 - 只分析少数几个染色体"""
    print("=" * 60)
    print("快速测试模式 - 仅分析前5个染色体")
    print("=" * 60)
    
    # 首先检查染色体命名
    print("检查染色体命名格式...")
    test_bam = "/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam"
    try:
        import subprocess
        result = subprocess.run(['samtools', 'view', '-H', test_bam], capture_output=True, text=True)
        if result.returncode == 0:
            lines = [line for line in result.stdout.split('\n') if line.startswith('@SQ')]
            print(f"找到 {len(lines)} 个染色体序列")
            print("前几个染色体名称:")
            for line in lines[:5]:
                if '\tSN:' in line:
                    chr_name = line.split('SN:')[1].split('\t')[0]
                    print(f"  - {chr_name}")
                    
            # 判断使用哪种命名
            first_chr = lines[0].split('SN:')[1].split('\t')[0] if lines else ""
            use_chr_prefix = first_chr.startswith('chr')
            print(f"染色体命名格式: {'chr前缀' if use_chr_prefix else '无chr前缀'}")
        else:
            print("无法读取BAM头信息")
            use_chr_prefix = True  # 默认使用chr前缀
    except Exception as e:
        print(f"检查染色体命名失败: {e}")
        use_chr_prefix = True
    
    # 修改BAMAnalyzer类进行快速测试
    sample1_bams = ["/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam"]
    sample2_bams = ["/mnt/disk1/rna_seq_analysis/step0/SRR12125134_Aligned.sortedByCoord.out.bam"]
    
    class QuickBAMAnalyzer(BAMAnalyzer):
        def analyze_basic_stats(self, bam_files, group_name):
            """修改基础统计分析，使用正确的染色体命名"""
            import pysam
            import numpy as np
            
            print(f"分析 {group_name} 基础统计...")
            
            total_reads = 0
            mapped_reads = 0
            duplicate_reads = 0
            mapq_scores = []
            insert_sizes = []
            
            # 根据检测结果使用正确的染色体命名
            if use_chr_prefix:
                test_regions = ['chr1:10000000-10100000', 'chr2:10000000-10100000']
            else:
                test_regions = ['1:10000000-10100000', '2:10000000-10100000']
            
            for bam_file in bam_files:
                try:
                    with pysam.AlignmentFile(bam_file, "rb") as bam:
                        for region in test_regions:
                            try:
                                for read in bam.fetch(region=region):
                                    total_reads += 1
                                    
                                    if not read.is_unmapped:
                                        mapped_reads += 1
                                        mapq_scores.append(read.mapping_quality)
                                        
                                        if read.is_proper_pair and read.template_length > 0:
                                            insert_sizes.append(abs(read.template_length))
                                    
                                    if read.is_duplicate:
                                        duplicate_reads += 1
                                    
                                    if total_reads > 10000:  # 减少到1万reads
                                        break
                                
                                if total_reads > 10000:
                                    break
                                    
                            except Exception as region_error:
                                print(f"    跳过区域 {region}: {region_error}")
                                continue
                        
                        if total_reads > 10000:
                            break
                            
                except Exception as e:
                    print(f"    处理文件错误 {bam_file}: {e}")
                    continue
            
            # 计算统计信息
            mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0
            avg_mapq = np.mean(mapq_scores) if mapq_scores else 0
            avg_insert_size = int(np.mean(insert_sizes)) if insert_sizes else 0
            
            print(f"    读取reads: {total_reads}, 映射率: {mapping_rate:.1f}%")
            
            return {
                'total_reads': total_reads,
                'mapped_reads': mapped_reads,
                'duplicate_reads': duplicate_reads,
                'mapping_rate': round(mapping_rate, 1),
                'avg_mapq': round(avg_mapq, 1),
                'avg_insert_size': avg_insert_size,
                'mapq_scores': mapq_scores,
                'insert_sizes': insert_sizes
            }
            
        def analyze_all_chromosomes_coverage(self, bam_files, group_index):
            """只分析前5个染色体，使用正确的命名格式"""
            import subprocess
            
            print(f"Quick test: analyzing first 5 chromosomes for group {group_index + 1}")
            
            if use_chr_prefix:
                chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
            else:
                chromosomes = ['1', '2', '3', '4', '5']
                
            coverage_data = []
            
            for chr_idx, chr_name in enumerate(chromosomes):
                region = f"{chr_name}:10000000-10200000"  # 增加到200kb区域
                
                try:
                    cmd = ['samtools', 'coverage', '-r', region] + bam_files
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                    
                    if result.returncode == 0:
                        lines = result.stdout.strip().split('\n')
                        if len(lines) >= 2:
                            data = lines[1].split('\t')
                            if len(data) >= 6:
                                coverage_percentage = float(data[5])
                            else:
                                coverage_percentage = 0
                        else:
                            coverage_percentage = 0
                    else:
                        print(f"    警告 {region}: {result.stderr.strip()}")
                        coverage_percentage = 0
                        
                except Exception as e:
                    print(f"    Error in {region}: {e}")
                    coverage_percentage = 0
                
                coverage_data.append([chr_idx, group_index, round(coverage_percentage, 1)])
                print(f"    {chr_name}: {coverage_percentage:.1f}%")
            
            return coverage_data
    
    try:
        analyzer = QuickBAMAnalyzer(
            sample1_bams=sample1_bams,
            sample2_bams=sample2_bams,
            output_file="quick_test_result.json",
            region_size=500000
        )
        
        result = analyzer.run_analysis()
        
        print(f"快速测试完成!")
        print(f"分析了 5 个染色体，每个 50kb 区域")
        print(f"结果保存到: quick_test_result.json")
        
    except Exception as e:
        print(f"快速测试失败: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    print("BAM分析器测试脚本")
    print("1 - 完整测试 (24个染色体)")
    print("2 - 快速测试 (5个染色体)")
    print("0 - 退出")
    
    while True:
        choice = input("\n请选择测试模式 (1/2/0): ").strip()
        
        if choice == "1":
            test_bam_analyzer()
            break
        elif choice == "2":
            quick_test()
            break
        elif choice == "0":
            print("退出测试")
            break
        else:
            print("无效选择，请输入 1, 2 或 0")