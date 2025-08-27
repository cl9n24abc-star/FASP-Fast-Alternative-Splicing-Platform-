#!/usr/bin/env python3
"""
Direct test of BAM analyzer script
直接测试BAM分析脚本
"""

import subprocess
import os
import json
import sys
from pathlib import Path

def test_bam_analyzer():
    """直接测试BAM分析器"""
    
    # 配置参数
    sample1_list = "/home/cl9n24abc/step/sample1_bams.txt"
    sample2_list = "/home/cl9n24abc/step/sample2_bams.txt"
    output_file = "../frontend/public/bam_analysis_result.json"  # 修改为相对路径
    chromosome = "1"  # 改为不带chr前缀的格式
    
    print("=" * 60)
    print("Direct BAM Analysis Test")
    print("=" * 60)
    print(f"Sample 1 list: {sample1_list}")
    print(f"Sample 2 list: {sample2_list}")
    print(f"Output file: {output_file}")
    print(f"Chromosome: {chromosome}")
    print()
    
    # 检查输入文件
    if not os.path.exists(sample1_list):
        print(f"ERROR: Sample 1 BAM list not found: {sample1_list}")
        return False
    
    if not os.path.exists(sample2_list):
        print(f"ERROR: Sample 2 BAM list not found: {sample2_list}")
        return False
    
    print("Input files found, starting analysis...")
    print()
    
    # 查找分析脚本
    script_locations = [
        "./bam_analyzer.py",
        "./services/bam_analyzer.py",  # 添加services目录
        "../bam_analyzer.py", 
        "./scripts/bam_analyzer.py",
        "../scripts/bam_analyzer.py",
        "../services/bam_analyzer.py"
    ]
    
    analyzer_script = None
    for location in script_locations:
        if os.path.exists(location):
            analyzer_script = location
            break
    
    if not analyzer_script:
        print("ERROR: bam_analyzer.py script not found in expected locations:")
        for loc in script_locations:
            print(f"  - {loc}")
        return False
    
    print(f"Using analyzer script: {analyzer_script}")
    print()
    
    # 构建命令 - 使用简单的位置参数
    cmd = [
        sys.executable,  # 使用当前Python解释器
        analyzer_script,
        sample1_list,
        sample2_list
    ]
    
    print("Running command:")
    print(" ".join(cmd))
    print()
    print("-" * 40)
    
    try:
        # 运行分析
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800  # 30分钟超时
        )
        
        print("STDOUT:")
        print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        print("-" * 40)
        print(f"Return code: {result.returncode}")
        
        if result.returncode == 0:
            print("SUCCESS: Analysis completed!")
            
            # 检查多个可能的输出文件位置
            possible_outputs = [
                output_file,  # ./bam_analysis_result.json
                "/home/cl9n24abc/step/bam_analysis_result.json",  # 从日志看到的实际位置
                os.path.join(os.path.dirname(sample1_list), "bam_analysis_result.json"),  # 输入文件目录
                "bam_analysis_result.json"  # 当前目录
            ]
            
            found_output = None
            for possible_output in possible_outputs:
                if os.path.exists(possible_output):
                    found_output = possible_output
                    break
            
            if found_output:
                file_size = os.path.getsize(output_file)
                print(f"Output file generated: {output_file}")
                print(f"File size: {file_size / 1024:.1f} KB")
                
                # 显示结果摘要
                try:
                    with open(output_file, 'r') as f:
                        data = json.load(f)
                    
                    print()
                    print("=" * 40)
                    print("ANALYSIS SUMMARY")
                    print("=" * 40)
                    
                    # 基本信息
                    analysis_info = data.get('analysis_info', {})
                    print(f"Chromosome analyzed: {analysis_info.get('chromosome_analyzed', 'N/A')}")
                    print(f"Analysis date: {analysis_info.get('analysis_date', 'N/A')}")
                    print(f"Analysis time: {analysis_info.get('analysis_time_seconds', 'N/A')}s")
                    print(f"Total samples: {analysis_info.get('total_samples', 'N/A')}")
                    print()
                    
                    # 组别统计
                    group1 = data.get('group1', {})
                    group2 = data.get('group2', {})
                    
                    print("Group Statistics:")
                    print(f"  Sample 1 - File size: {group1.get('fileSize', 'N/A')}")
                    print(f"  Sample 1 - Mapping rate: {group1.get('mappingRate', 'N/A')}%")
                    print(f"  Sample 1 - Avg MAPQ: {group1.get('avgMapQ', 'N/A')}")
                    print(f"  Sample 1 - Avg insert size: {group1.get('avgInsertSize', 'N/A')} bp")
                    print()
                    print(f"  Sample 2 - File size: {group2.get('fileSize', 'N/A')}")
                    print(f"  Sample 2 - Mapping rate: {group2.get('mappingRate', 'N/A')}%")
                    print(f"  Sample 2 - Avg MAPQ: {group2.get('avgMapQ', 'N/A')}")
                    print(f"  Sample 2 - Avg insert size: {group2.get('avgInsertSize', 'N/A')} bp")
                    print()
                    
                    # 图表数据统计
                    chart_data = data.get('chart_data', {})
                    if chart_data:
                        print("Chart Data Generated:")
                        if chart_data.get('pie_chart'):
                            print(f"  - Pie chart: {len(chart_data['pie_chart'])} categories")
                        if chart_data.get('mapq_histogram'):
                            print(f"  - MAPQ histogram: {len(chart_data['mapq_histogram'])} bins")
                        if chart_data.get('volcano_plot'):
                            print(f"  - Volcano plot: {len(chart_data['volcano_plot'])} data points")
                        if chart_data.get('chromosome_heatmap'):
                            print(f"  - Chromosome heatmap: {len(chart_data['chromosome_heatmap'])} data points")
                    
                    print()
                    print("Analysis completed successfully!")
                    print("You can now copy this JSON file to your frontend public directory.")
                    
                except json.JSONDecodeError:
                    print("WARNING: Output file exists but contains invalid JSON")
                    return False
                except Exception as e:
                    print(f"WARNING: Error reading output file: {e}")
                    return False
                    
            else:
                print("ERROR: Analysis completed but no output file generated")
                return False
                
        else:
            print("ERROR: Analysis failed")
            print(f"Error output: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("ERROR: Analysis timed out (>30 minutes)")
        return False
    except FileNotFoundError:
        print("ERROR: Python interpreter or script not found")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}")
        return False
    
    return True

if __name__ == "__main__":
    success = test_bam_analyzer()
    if success:
        print()
        print("=" * 60)
        print("TEST PASSED")
        print("=" * 60)
        sys.exit(0)
    else:
        print()
        print("=" * 60)
        print("TEST FAILED")
        print("=" * 60)
        sys.exit(1)