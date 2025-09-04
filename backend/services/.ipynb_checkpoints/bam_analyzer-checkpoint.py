#!/usr/bin/env python3
"""
BAM File Analyzer - Improved Version
改进版BAM文件分析器 - 修复饼图逻辑并分析所有染色体小区域
"""

import pysam
import numpy as np
import json
import os
import sys
from collections import defaultdict, Counter
import subprocess
import argparse
from pathlib import Path
import time
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BAMAnalyzer:
    def __init__(self, sample1_bams, sample2_bams, output_file, region_size=3000000, gtf_file=None):
        """
        初始化BAM分析器
        
        Args:
            sample1_bams: Sample 1的BAM文件路径列表
            sample2_bams: Sample 2的BAM文件路径列表
            output_file: 输出JSON文件路径
            region_size: 每个染色体分析的区域大小 (默认100kb)
        """
        self.sample1_bams = sample1_bams
        self.sample2_bams = sample2_bams
        self.output_file = output_file
        self.region_size = region_size
        self.gtf_file = gtf_file
        # 验证文件存在性
        self._validate_files()
        
        logger.info(f"Initialized BAM analyzer")
        logger.info(f"Sample 1: {len(sample1_bams)} BAM files")
        logger.info(f"Sample 2: {len(sample2_bams)} BAM files")
        logger.info(f"Region size per chromosome: {region_size} bp")
    def parse_gtf_genes(self):
        if not self.gtf_file or not os.path.exists(self.gtf_file):
            return []
        
        genes = []
        try:
            with open(self.gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and parts[2] == 'gene':
                        chrom = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        
                        # 提取基因名
                        attributes = parts[8]
                        gene_name = 'Unknown'
                        for attr in attributes.split(';'):
                            if 'gene_name' in attr:
                                gene_name = attr.split('"')[1]
                                break
                        
                        genes.append({
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'name': gene_name,
                            'length': end - start
                        })
        except Exception as e:
            logger.warning(f"Failed to parse GTF: {e}")
            return []
        
        return genes[:500]  # 限制分析基因数量

    def count_reads_in_region(self, bam_files, chrom, start, end):
        """统计指定区域的读段数"""
        total_reads = 0
        
        for bam_file in bam_files:
            try:
                with pysam.AlignmentFile(bam_file, "rb") as bam:
                    # 使用samtools count快速统计
                    reads = bam.count(chrom, start, end)
                    total_reads += reads
            except Exception as e:
                logger.warning(f"Failed to count reads in {bam_file}: {e}")
                continue
        
        return total_reads
    
    def _validate_files(self):
        """验证BAM文件是否存在且有索引"""
        all_bams = self.sample1_bams + self.sample2_bams
        for bam_file in all_bams:
            if not os.path.exists(bam_file):
                raise FileNotFoundError(f"BAM file not found: {bam_file}")
            
            # 检查索引文件
            index_file = bam_file + ".bai"
            if not os.path.exists(index_file):
                logger.warning(f"Index file missing for {bam_file}, creating...")
                pysam.index(bam_file)
    
    def get_file_sizes(self, bam_files):
        """计算文件总大小"""
        total_size = sum(os.path.getsize(f) for f in bam_files)
        return f"{total_size / (1024**3):.1f} GB"
    
    def analyze_basic_stats(self, bam_files, group_name):
        """分析基础统计信息 - 修复MAPQ采样"""
        logger.info(f"Analyzing basic stats for {group_name}...")
        
        total_reads = 0
        mapped_reads = 0
        duplicate_reads = 0
        mapq_scores = []
        insert_sizes = []
        
        for bam_file in bam_files:
            try:
                # 使用samtools flagstat获取全局统计（保持不变）
                result = subprocess.run(['samtools', 'flagstat', bam_file], 
                                      capture_output=True, text=True, timeout=120)
                if result.returncode == 0:
                    lines = result.stdout.strip().split('\n')
                    for line in lines:
                        if 'in total' in line:
                            total_reads += int(line.split()[0])
                        elif 'mapped (' in line and 'primary' not in line:
                            mapped_reads += int(line.split()[0])
                        elif 'duplicates' in line:
                            duplicate_reads += int(line.split()[0])
                
                # 修复采样方法 - 只分析主要比对
                with pysam.AlignmentFile(bam_file, "rb") as bam:
                    sample_count = 0
                    mapq_debug = {}
                    
                    for read in bam.fetch():
                        if sample_count >= 100000:  # 减少到10万样本
                            break
                        
                        # 只处理主要比对的reads
                        if (not read.is_unmapped and 
                            not read.is_secondary and 
                            not read.is_supplementary):  # 添加这个检查
                            
                            mapq = read.mapping_quality
                            mapq_debug[mapq] = mapq_debug.get(mapq, 0) + 1
                            mapq_scores.append(mapq)
                            
                            # 获取插入片段大小
                            if read.is_proper_pair and read.template_length > 0:
                                insert_sizes.append(abs(read.template_length))
                        
                        sample_count += 1
                    
                    # 输出调试信息（只显示主要比对的MAPQ）
                    primary_count = sum(mapq_debug.values())
                    logger.info(f"  Primary alignments in sample: {primary_count}/{sample_count}")
                    logger.info(f"  MAPQ distribution (primary only):")
                    for mapq in sorted(mapq_debug.keys()):
                        count = mapq_debug[mapq]
                        percent = count/primary_count*100 if primary_count > 0 else 0
                        logger.info(f"    MAPQ {mapq}: {count} reads ({percent:.1f}%)")
                            
            except Exception as e:
                logger.error(f"Error processing {bam_file}: {e}")
                continue
        
        # 计算统计信息（保持不变）
        mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0
        avg_mapq = np.mean(mapq_scores) if mapq_scores else 0
        avg_insert_size = int(np.mean(insert_sizes)) if insert_sizes else 0
        
        logger.info(f"  {group_name} MAPQ samples: {len(mapq_scores)}, avg MAPQ: {avg_mapq:.1f}")
        
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
    def generate_pie_chart_data(self, group1_stats, group2_stats):
        """生成修正的映射统计饼图数据"""
        # 合并两组统计数据
        total_reads = group1_stats['total_reads'] + group2_stats['total_reads']
        mapped_reads = group1_stats['mapped_reads'] + group2_stats['mapped_reads']
        duplicate_reads = group1_stats['duplicate_reads'] + group2_stats['duplicate_reads']
        
        # 计算百分比
        mapped_percent = (mapped_reads / total_reads * 100) if total_reads > 0 else 0
        unmapped_percent = ((total_reads - mapped_reads) / total_reads * 100) if total_reads > 0 else 0
        duplicate_percent = (duplicate_reads / total_reads * 100) if total_reads > 0 else 0
        
        logger.info(f"Pie chart stats: Total={total_reads}, Mapped={mapped_percent:.1f}%, Unmapped={unmapped_percent:.1f}%, Duplicates={duplicate_percent:.1f}%")
        
        return [
            {
                "value": round(mapped_percent, 1),
                "name": "Successfully Mapped",
                "itemStyle": {"color": "#5470c6"}
            },
            {
                "value": round(unmapped_percent, 1), 
                "name": "Unmapped",
                "itemStyle": {"color": "#ee6666"}
            },
            {
                "value": round(duplicate_percent, 1),
                "name": "Duplicates", 
                "itemStyle": {"color": "#fac858"}
            }
        ]
    
    def generate_mapq_histogram(self, mapq_scores):
        """生成MAPQ直方图数据"""
        ranges = [(0, 10), (11, 20), (21, 30), (31, 40), (41, 50), (51, 60)]
        histogram = []
        
        for min_val, max_val in ranges:
            count = sum(1 for score in mapq_scores if min_val <= score <= max_val)
            histogram.append({
                'range': f'{min_val}-{max_val}',
                'count': count
            })
        
        return histogram
    
    def generate_insert_size_distribution(self, group1_sizes, group2_sizes):
        """生成插入片段大小分布"""
        x_values = list(range(100, 601, 10))  # 100-600bp，间隔10bp
        
        def calculate_density(sizes, x_vals):
            if not sizes:
                return [0] * len(x_vals)
            
            # 使用核密度估计
            hist, bin_edges = np.histogram(sizes, bins=50, range=(100, 600), density=True)
            # 插值到目标x值
            density = np.interp(x_vals, bin_edges[:-1], hist)
            return [max(0, d * 1000) for d in density]  # 放大以便显示
        
        return {
            'x_values': x_values,
            'group1': calculate_density(group1_sizes, x_values),
            'group2': calculate_density(group2_sizes, x_values)
        }
    
    def analyze_all_chromosomes_coverage(self, bam_files, group_index):
        """分析所有染色体的小区域覆盖度"""
        logger.info(f"Analyzing all chromosomes coverage for group {group_index + 1}...")
        
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        coverage_data = []
        
        for chr_idx, chrom in enumerate(chromosomes):
            chr_name = chrom
            # 分析每个染色体的一个小区域
            start_pos = 10000000  # 从10M开始避开端粒
            end_pos = start_pos + self.region_size
            region = f"{chr_name}:{start_pos}-{end_pos}"
            
            try:
                # 使用samtools coverage分析指定区域
                cmd = ['samtools', 'coverage', '-r', region] + bam_files
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                
                if result.returncode == 0:
                    lines = result.stdout.strip().split('\n')
                    if len(lines) >= 2:
                        data = lines[1].split('\t')
                        if len(data) >= 6:
                            # samtools coverage输出格式第6列是覆盖度百分比
                            coverage_percentage = float(data[5])
                        else:
                            coverage_percentage = 0
                    else:
                        coverage_percentage = 0
                else:
                    logger.warning(f"Failed to analyze {region}: {result.stderr}")
                    coverage_percentage = 0
                    
            except subprocess.TimeoutExpired:
                logger.warning(f"Timeout analyzing {region}")
                coverage_percentage = 0
            except Exception as e:
                logger.error(f"Error analyzing {region}: {e}")
                coverage_percentage = 0
            
            # 添加到热力图数据 [染色体索引, 组别索引, 覆盖度]
            coverage_data.append([chr_idx, group_index, round(coverage_percentage, 1)])
            
            # 打印进度
            logger.info(f"  {chr_name}: {coverage_percentage:.1f}% coverage")
        
        return coverage_data
    
    def generate_gene_density_analysis(self, group1_stats, group2_stats):
        """生成真实的基因读段密度分析数据"""
        
        # 如果没有GTF文件，回退到模拟数据
        if not self.gtf_file:
            logger.warning("No GTF file provided, using simulated data")
            return self._generate_simulated_density_data()
        
        genes = self.parse_gtf_genes()
        if not genes:
            logger.warning("No genes parsed from GTF, using simulated data")
            return self._generate_simulated_density_data()
        
        density_data = []
        logger.info(f"Analyzing read density for {len(genes)} genes...")
        
        for i, gene in enumerate(genes):
            try:
                # 统计两组样本的读段数
                reads1 = self.count_reads_in_region(self.sample1_bams, gene['chrom'], gene['start'], gene['end'])
                reads2 = self.count_reads_in_region(self.sample2_bams, gene['chrom'], gene['start'], gene['end'])
                
                # 避免除零错误
                if reads1 == 0 and reads2 == 0:
                    continue
                if reads1 == 0:
                    reads1 = 1
                if reads2 == 0:
                    reads2 = 1
                
                # 计算密度比值和总读段数
                density_ratio = reads2 / reads1
                log2_ratio = np.log2(density_ratio)
                total_reads = reads1 + reads2
                log_total_reads = np.log10(total_reads) if total_reads > 0 else 0
                
                # 分类 - 匹配前端的分类名称
                category = 'similar'
                if abs(log2_ratio) > np.log2(1.5):  # 1.5倍差异阈值
                    if log2_ratio > 0:
                        category = 'sample2-enriched'
                    else:
                        category = 'sample1-enriched'
                elif total_reads < 100:
                    category = 'low-coverage'
                
                density_data.append([log2_ratio, log_total_reads, category, gene['name']])
                
                # 打印进度
                if (i + 1) % 50 == 0:
                    logger.info(f"  Processed {i + 1}/{len(genes)} genes")
                    
            except Exception as e:
                logger.warning(f"Error analyzing gene {gene['name']}: {e}")
                continue
        
        logger.info(f"Generated density data for {len(density_data)} genes")
        return density_data
    
    def _generate_simulated_density_data(self):
        """生成模拟的基因读段密度数据（用于没有GTF文件时）"""
        np.random.seed(42)
        
        gene_names = [
            'BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS', 'PIK3CA', 'PTEN', 'RB1', 'APC', 'VHL',
            'BRAF', 'PIK3R1', 'ARID1A', 'CTNNB1', 'SMAD4', 'FBXW7', 'NRAS', 'PPP2R1A'
        ]
        
        density_data = []
        for i in range(200):
            # 模拟读段数据
            reads1 = max(1, int(np.random.exponential(500)))
            reads2 = max(1, int(np.random.exponential(500)))
            
            density_ratio = reads2 / reads1
            log2_ratio = np.log2(density_ratio)
            total_reads = reads1 + reads2
            log_total_reads = np.log10(total_reads)
            
            # 分类
            category = 'similar'
            if abs(log2_ratio) > np.log2(1.5):
                if log2_ratio > 0:
                    category = 'sample2-enriched'
                else:
                    category = 'sample1-enriched'
            elif total_reads < 200:
                category = 'low-coverage'
            
            gene_name = gene_names[i] if i < len(gene_names) else f'Gene_{i+1}'
            density_data.append([log2_ratio, log_total_reads, category, gene_name])
        
        return density_data
    def run_analysis(self):
        """运行完整分析"""
        start_time = time.time()
        logger.info("Starting BAM analysis...")
        
        # 分析两个样本组的基础统计
        group1_stats = self.analyze_basic_stats(self.sample1_bams, "Group 1")
        group2_stats = self.analyze_basic_stats(self.sample2_bams, "Group 2")
        
        # 生成图表数据
        logger.info("Generating chart data...")
        
        # 1. 修正的饼图数据
        pie_data = self.generate_pie_chart_data(group1_stats, group2_stats)
        
        # 2. MAPQ直方图
        mapq_histogram = self.generate_mapq_histogram(
            group1_stats['mapq_scores'] + group2_stats['mapq_scores']
        )
        
        # 3. 插入片段大小分布
        insert_size_dist = self.generate_insert_size_distribution(
            group1_stats['insert_sizes'], 
            group2_stats['insert_sizes']
        )
        
        # 4. 染色体覆盖度热图 - 分析所有染色体
        logger.info("Analyzing chromosome coverage for all chromosomes...")
        coverage_data = []
        coverage_data.extend(self.analyze_all_chromosomes_coverage(self.sample1_bams, 0))
        coverage_data.extend(self.analyze_all_chromosomes_coverage(self.sample2_bams, 1))
        
        # 5. 火山图数据
        volcano_data = self.generate_gene_density_analysis(group1_stats, group2_stats)
        
        # 构建最终结果
        result = {
            "group1": {
                "fileSize": self.get_file_sizes(self.sample1_bams),
                "mappingRate": group1_stats['mapping_rate'],
                "avgMapQ": group1_stats['avg_mapq'],
                "avgInsertSize": group1_stats['avg_insert_size']
            },
            "group2": {
                "fileSize": self.get_file_sizes(self.sample2_bams),
                "mappingRate": group2_stats['mapping_rate'],
                "avgMapQ": group2_stats['avg_mapq'],
                "avgInsertSize": group2_stats['avg_insert_size']
            },
            "chart_data": {
                "pie_chart": pie_data,
                "mapq_histogram": mapq_histogram,
                "insert_size_distribution": insert_size_dist,
                "chromosome_heatmap": coverage_data,
                "volcano_plot": volcano_data
            },
            "analysis_info": {
                "region_size_per_chromosome": self.region_size,
                "total_chromosomes_analyzed": 24,
                "analysis_date": time.strftime("%Y-%m-%d %H:%M:%S"),
                "total_samples": len(self.sample1_bams) + len(self.sample2_bams),
                "analysis_time_seconds": round(time.time() - start_time, 2)
            }
        }
        
        # 保存结果
        logger.info(f"Saving results to {self.output_file}")
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        with open(self.output_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        
        # 打印统计摘要
        logger.info(f"Analysis completed in {time.time() - start_time:.2f} seconds")
        logger.info(f"Results saved to: {self.output_file}")
        
        # 打印覆盖度统计
        group1_coverage = [d[2] for d in coverage_data if d[1] == 0]
        group2_coverage = [d[2] for d in coverage_data if d[1] == 1]
        
        if group1_coverage and group2_coverage:
            logger.info(f"Group 1 coverage: avg={np.mean(group1_coverage):.1f}%, range={min(group1_coverage):.1f}-{max(group1_coverage):.1f}%")
            logger.info(f"Group 2 coverage: avg={np.mean(group2_coverage):.1f}%, range={min(group2_coverage):.1f}-{max(group2_coverage):.1f}%")
        
        return result

def read_bam_list(file_path):
    """读取BAM列表文件"""
    bam_files = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                if os.path.exists(line):
                    bam_files.append(line)
                else:
                    logger.warning(f"BAM file not found: {line}")
    return bam_files

def main():
    parser = argparse.ArgumentParser(description='Analyze BAM files for visualization (improved version)')
    parser.add_argument('--sample1', required=True, help='Sample 1 BAM list file')
    parser.add_argument('--sample2', required=True, help='Sample 2 BAM list file')
    parser.add_argument('--output', default='../frontend/public/bam_analysis_result.json', help='Output JSON file')
    parser.add_argument('--region-size', type=int, default=100000, help='Region size per chromosome (default: 100kb)')
    
    args = parser.parse_args()
    
    try:
        # 读取BAM文件列表
        sample1_bams = read_bam_list(args.sample1)
        sample2_bams = read_bam_list(args.sample2)
        
        if not sample1_bams:
            raise ValueError(f"No valid BAM files found in {args.sample1}")
        if not sample2_bams:
            raise ValueError(f"No valid BAM files found in {args.sample2}")
        
        # 创建分析器并运行
        analyzer = BAMAnalyzer(sample1_bams, sample2_bams, args.output, args.region_size)
        result = analyzer.run_analysis()
        
        print(f"Analysis completed successfully!")
        print(f"Results saved to: {args.output}")
        print(f"Region size per chromosome: {args.region_size} bp")
        print(f"Sample 1 mapping rate: {result['group1']['mappingRate']}%")
        print(f"Sample 2 mapping rate: {result['group2']['mappingRate']}%")
        print(f"Analysis time: {result['analysis_info']['analysis_time_seconds']} seconds")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()