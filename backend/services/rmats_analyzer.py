#!/usr/bin/env python3
"""
rMATS快速分析器 - 10分钟内生成Vue组件需要的JSON数据
使用方法: python rmats_analyzer.py /path/to/rmats/output
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
import glob

def find_rmats_files(input_path):
    """查找rMATS结果文件"""
    event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
    files = {}
    
    # 查找文件的多种模式
    patterns = [
        f"{input_path}/*.MATS.JC.txt",
        f"{input_path}/**/*.MATS.JC.txt",
        f"{input_path}/*.MATS.JCEC.txt"
    ]
    
    for pattern in patterns:
        for file_path in glob.glob(pattern, recursive=True):
            filename = os.path.basename(file_path)
            for event_type in event_types:
                if filename.startswith(event_type + "."):
                    files[event_type] = file_path
                    break
    
    print(f"找到rMATS文件: {list(files.keys())}")
    return files

def quick_analyze_rmats(input_path, output_path=None):
    """快速分析rMATS结果"""
    
    # 设置输出路径
    if output_path is None:
        output_path = input_path
    
    # 查找文件
    rmats_files = find_rmats_files(input_path)
    if not rmats_files:
        raise ValueError(f"未找到rMATS结果文件在: {input_path}")
    
    # 初始化结果
    result = {
        "rmatsData": {
            "totalEvents": 0,
            "significantEvents": 0,
            "maxPsiIncrease": 0.0,
            "maxPsiDecrease": 0.0,
            "group1Samples": 3,
            "group2Samples": 3,
            "eventCounts": {}
        },
        "chartData": {
            "eventDistribution": [],
            "significanceStats": {
                "categories": [],
                "significantData": [],
                "nonSignificantData": []
            },
            "volcanoPlot": []
        },
        "metadata": {
            "analysisDate": pd.Timestamp.now().isoformat(),
            "inputPath": input_path,
            "fdrThreshold": 0.05
        }
    }
    
    all_psi_values = []
    event_colors = {
        'SE': '#5470c6',
        'A5SS': '#91cc75', 
        'A3SS': '#fac858',
        'MXE': '#ee6666',
        'RI': '#73c0de'
    }
    
    event_names = {
        'SE': 'SE (Skipped Exon)',
        'A5SS': "A5SS (Alt 5' Splice Site)",
        'A3SS': "A3SS (Alt 3' Splice Site)", 
        'MXE': 'MXE (Mutually Exclusive Exon)',
        'RI': 'RI (Retained Intron)'
    }
    
    # 处理每个事件类型文件
    for event_type, file_path in rmats_files.items():
        print(f"处理 {event_type} 文件: {file_path}")
        
        try:
            # 读取文件
            df = pd.read_csv(file_path, sep='\t')
            print(f"  加载了 {len(df)} 个 {event_type} 事件")
            
            # 基础统计
            total_events = len(df)
            significant_events = 0
            
            if 'FDR' in df.columns:
                significant_events = len(df[df['FDR'] < 0.05])
            elif 'PValue' in df.columns:
                significant_events = len(df[df['PValue'] < 0.05])
            
            # 更新统计
            result["rmatsData"]["totalEvents"] += total_events
            result["rmatsData"]["significantEvents"] += significant_events
            result["rmatsData"]["eventCounts"][event_type] = total_events
            
            # 添加到图表数据
            result["chartData"]["significanceStats"]["categories"].append(event_type)
            result["chartData"]["significanceStats"]["significantData"].append(significant_events)
            result["chartData"]["significanceStats"]["nonSignificantData"].append(total_events - significant_events)
            
            # 饼图数据
            if total_events > 0:
                result["chartData"]["eventDistribution"].append({
                    "name": event_names.get(event_type, event_type),
                    "value": total_events,
                    "itemStyle": {"color": event_colors.get(event_type, '#d0d0d0')}
                })
            
            # PSI值统计
            if 'IncLevelDifference' in df.columns:
                psi_values = df['IncLevelDifference'].dropna()
                all_psi_values.extend(psi_values.tolist())
            
            # 样本数量推断(只做一次)
            if result["rmatsData"]["group1Samples"] == 3:  # 还是默认值
                try:
                    if 'IncLevel1' in df.columns and len(df) > 0:
                        inc1 = str(df['IncLevel1'].iloc[0])
                        inc2 = str(df['IncLevel2'].iloc[0])
                        if ',' in inc1 and ',' in inc2:
                            result["rmatsData"]["group1Samples"] = len(inc1.split(','))
                            result["rmatsData"]["group2Samples"] = len(inc2.split(','))
                except:
                    pass
            
            # 火山图数据 (每个类型取前50个)
            volcano_sample = df.head(50)
            for _, row in volcano_sample.iterrows():
                try:
                    if 'FDR' in row and 'IncLevelDifference' in row:
                        fdr = float(row['FDR'])
                        delta_psi = float(row['IncLevelDifference'])
                        
                        if pd.notna(fdr) and pd.notna(delta_psi) and fdr > 0:
                            gene_name = str(row.get('geneSymbol', f'{event_type}_event'))
                            
                            category = 'significant' if fdr < 0.05 else 'non-significant'
                            
                            result["chartData"]["volcanoPlot"].append({
                                "eventName": gene_name,
                                "deltaPsi": round(delta_psi, 4),
                                "fdr": fdr,
                                "negLog10Fdr": round(-np.log10(fdr), 4),
                                "category": category,
                                "eventType": event_type
                            })
                except:
                    continue
                    
        except Exception as e:
            print(f"  警告: 处理 {event_type} 文件时出错: {e}")
            continue
    
    # 计算PSI范围
    if all_psi_values:
        result["rmatsData"]["maxPsiIncrease"] = round(max(all_psi_values), 3)
        result["rmatsData"]["maxPsiDecrease"] = round(min(all_psi_values), 3)
    
    output_file = "../../frontend/public/rmats_analysis_result.json"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ 分析完成!")
    print(f"📊 总事件数: {result['rmatsData']['totalEvents']}")
    print(f"⭐ 显著事件: {result['rmatsData']['significantEvents']}")
    print(f"🧬 事件类型: {list(result['rmatsData']['eventCounts'].keys())}")
    print(f"📁 结果文件: {output_file}")
    
    return output_file

if __name__ == "__main__":
    # 命令行使用
    if len(sys.argv) < 2:
        print("使用方法: python rmats_analyzer.py <rMATS输出目录> [输出目录]")
        print("示例: python rmats_analyzer.py /data/rmats_output /data/analysis_results")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else input_path
    
    try:
        result_file = quick_analyze_rmats(input_path, output_path)
        print(f"\n🎉 成功生成JSON文件: {result_file}")
        print("现在可以将此文件用于Vue组件!")
        
    except Exception as e:
        print(f"❌ 分析失败: {e}")
        sys.exit(1)