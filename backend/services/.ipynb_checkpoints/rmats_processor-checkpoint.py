import os
import json
from typing import Dict, Any, Optional
from config import Config

class RMATSProcessor:
    """rMATS数据处理服务类"""
    
    @classmethod
    def process_single_analysis(cls, gtf_file: str, rmats_dir: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """处理单分析模式"""
        try:
            print("🔄 开始单分析模式处理...")
            
            # 如果没有指定输出目录，使用前端public目录
            if output_dir is None:
                output_dir = Config.get_frontend_public_path()
                Config.ensure_frontend_public_exists()
            
            print(f"📁 输出目录: {output_dir}")
            
            # 检查GTF文件是否提供
            if not gtf_file or gtf_file.strip() == '':
                print("⚠️ 未提供GTF文件，执行简化处理模式")
                return cls._process_rmats_without_gtf(rmats_dir, output_dir)
            else:
                print(f"📄 使用GTF文件: {gtf_file}")
                return cls._process_rmats_with_gtf(gtf_file, rmats_dir, output_dir)
                
        except Exception as e:
            import traceback
            print(f"❌ 单分析处理出错: {e}")
            print(traceback.format_exc())
            return {
                'success': False,
                'message': f'单分析处理失败: {str(e)}',
                'error': str(e)
            }
    
    @classmethod
    def _process_rmats_with_gtf(cls, gtf_file: str, rmats_dir: str, output_dir: str) -> Dict[str, Any]:
        """使用GTF文件的完整处理模式"""
        try:
            print("🧬 执行完整处理模式（包含基因注释）")
            
            # 导入并调用genestructure处理函数
            from .genestructure import process_rmats_data
            
            result = process_rmats_data(
                gtf_file=gtf_file,
                rmats_dir=rmats_dir,
                output_dir=output_dir,
                split_by_type=True
            )
            
            # 检查处理结果
            if result:
                event_count = len(result)
                print(f"✅ 完整处理成功! 处理了 {event_count} 个可变剪切事件")
                
                # 读取生成的摘要文件
                summary_data = cls._read_summary_file(output_dir)
                
                return {
                    'success': True,
                    'message': '完整分析模式处理成功，包含基因注释信息',
                    'analysis_type': 'single_with_gtf',
                    'event_count': event_count,
                    'summary': summary_data,
                    'output_directory': output_dir,
                    'public_files_ready': True,
                    'gtf_processed': True
                }
            else:
                return {
                    'success': False,
                    'message': '完整处理失败：未返回处理结果',
                    'error': 'No processing result returned'
                }
                
        except ImportError as e:
            return {
                'success': False,
                'message': '导入genestructure模块失败',
                'error': f'Import error: {str(e)}'
            }
    
    @classmethod
    def _process_rmats_without_gtf(cls, rmats_dir: str, output_dir: str) -> Dict[str, Any]:
        """不使用GTF文件的简化处理模式"""
        try:
            print("📊 执行简化处理模式（仅rMATS数据统计）")
            
            # 直接处理rMATS文件，生成基础统计信息
            stats = cls._generate_rmats_basic_stats(rmats_dir)
            
            if not stats['success']:
                return stats
            
            # 保存基础统计结果到输出目录
            output_file = os.path.join(output_dir, 'rmats_basic_stats.json')
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(stats['data'], f, indent=2, ensure_ascii=False)
            
            print(f"✅ 简化处理完成! 处理了 {stats['data']['total_events']} 个可变剪切事件")
            
            return {
                'success': True,
                'message': '简化分析模式处理成功，生成基础rMATS统计信息',
                'analysis_type': 'single_without_gtf',
                'event_count': stats['data']['total_events'],
                'summary': stats['data']['event_summary'],
                'output_directory': output_dir,
                'output_file': output_file,
                'gtf_processed': False,
                'note': '未使用GTF文件，仅包含rMATS基础统计信息'
            }
            
        except Exception as e:
            import traceback
            print(f"❌ 简化处理出错: {e}")
            print(traceback.format_exc())
            return {
                'success': False,
                'message': f'简化处理失败: {str(e)}',
                'error': str(e)
            }
    
    @classmethod
    def _generate_rmats_basic_stats(cls, rmats_dir: str) -> Dict[str, Any]:
        """生成rMATS基础统计信息（不依赖GTF）"""
        try:
            import pandas as pd
            
            event_files = {
                'SE': f"{rmats_dir}/SE.MATS.JC.txt",
                'MXE': f"{rmats_dir}/MXE.MATS.JC.txt", 
                'A5SS': f"{rmats_dir}/A5SS.MATS.JC.txt",
                'A3SS': f"{rmats_dir}/A3SS.MATS.JC.txt",
                'RI': f"{rmats_dir}/RI.MATS.JC.txt"
            }
            
            total_events = 0
            significant_events = 0
            event_summary = {}
            
            for event_type, file_path in event_files.items():
                try:
                    if os.path.exists(file_path):
                        df = pd.read_csv(file_path, sep='\t')
                        event_count = len(df)
                        
                        # 计算显著事件数（P值 < 0.05 且 FDR < 0.05）
                        significant_count = len(df[(df['PValue'] < 0.05) & (df['FDR'] < 0.05)])
                        
                        event_summary[event_type] = {
                            'total': event_count,
                            'significant': significant_count,
                            'file_exists': True
                        }
                        
                        total_events += event_count
                        significant_events += significant_count
                        
                        print(f"  {event_type}: {event_count} 事件 ({significant_count} 显著)")
                    else:
                        event_summary[event_type] = {
                            'total': 0,
                            'significant': 0,
                            'file_exists': False
                        }
                        print(f"  {event_type}: 文件不存在")
                        
                except Exception as e:
                    print(f"处理 {event_type} 文件时出错: {e}")
                    event_summary[event_type] = {
                        'total': 0,
                        'significant': 0,
                        'file_exists': False,
                        'error': str(e)
                    }
            
            return {
                'success': True,
                'data': {
                    'total_events': total_events,
                    'significant_events': significant_events,
                    'significance_rate': round(significant_events / total_events * 100, 2) if total_events > 0 else 0,
                    'event_summary': event_summary,
                    'processing_mode': 'basic_stats_only',
                    'gtf_used': False
                }
            }
            
        except Exception as e:
            return {
                'success': False,
                'message': f'生成基础统计失败: {str(e)}',
                'error': str(e)
            }
    
    @classmethod
    def process_comparative_analysis(cls, gtf_file: str, original_rmats_dir: str, 
                                   compare_rmats_dir: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """处理对比分析模式"""
        try:
            print("🔄 开始对比分析模式处理...")
            
            # 如果没有指定输出目录，使用前端public目录
            if output_dir is None:
                output_dir = Config.get_frontend_public_path()
                Config.ensure_frontend_public_exists()
            
            print(f"📁 输出目录: {output_dir}")
            print(f"📂 原始rMATS目录: {original_rmats_dir}")
            print(f"📂 对比rMATS目录: {compare_rmats_dir}")
            
            # TODO: 实现真正的对比分析逻辑
            # 目前先处理原始数据作为示例
            from .genestructure import process_rmats_data
            
            result = process_rmats_data(
                gtf_file=gtf_file,
                rmats_dir=original_rmats_dir,  # 先处理原始数据
                output_dir=output_dir,
                split_by_type=True
            )
            
            # 检查处理结果
            if result:
                event_count = len(result)
                print(f"✅ 对比分析处理成功! 处理了 {event_count} 个可变剪切事件")
                
                # 读取生成的摘要文件
                summary_data = cls._read_summary_file(output_dir)
                
                return {
                    'success': True,
                    'message': '对比分析模式处理成功，文件已保存到public目录',
                    'analysis_type': 'compare',
                    'event_count': event_count,
                    'summary': summary_data,
                    'output_directory': output_dir,
                    'public_files_ready': True,
                    'note': 'Currently processing original data only. Full comparative analysis to be implemented.'
                }
            else:
                return {
                    'success': False,
                    'message': '对比分析处理失败：未返回处理结果',
                    'error': 'No processing result returned'
                }
                
        except ImportError as e:
            return {
                'success': False,
                'message': '导入genestructure模块失败',
                'error': f'Import error: {str(e)}'
            }
        except Exception as e:
            import traceback
            print(f"❌ 对比分析处理出错: {e}")
            print(traceback.format_exc())
            return {
                'success': False,
                'message': f'对比分析处理失败: {str(e)}',
                'error': str(e)
            }
    
    @classmethod
    def check_rmats_files(cls, rmats_dir: str) -> Dict[str, Any]:
        """检查rMATS目录中的文件"""
        try:
            if not os.path.exists(rmats_dir) or not os.path.isdir(rmats_dir):
                return {
                    'valid': False,
                    'error': f'rMATS directory not found: {rmats_dir}'
                }
            
            file_status = {}
            for filename in Config.REQUIRED_RMATS_FILES:
                file_path = os.path.join(rmats_dir, filename)
                file_status[filename] = {
                    'exists': os.path.exists(file_path),
                    'size': os.path.getsize(file_path) if os.path.exists(file_path) else 0
                }
            
            all_exist = all(status['exists'] for status in file_status.values())
            
            return {
                'valid': all_exist,
                'files': file_status,
                'missing_files': [name for name, status in file_status.items() if not status['exists']]
            }
            
        except Exception as e:
            return {
                'valid': False,
                'error': f'Error checking rMATS files: {str(e)}'
            }
    
    @classmethod
    def get_processing_status(cls, output_dir: str) -> Dict[str, Any]:
        """获取处理状态"""
        try:
            # 检查输出目录是否存在
            if not os.path.exists(output_dir):
                return {
                    'status': 'not_started',
                    'message': 'Output directory does not exist'
                }
            
            # 检查是否有生成的文件
            expected_files = [f'frontend_data_{event_type}.json' for event_type in ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']]
            summary_file = 'frontend_data_summary.json'
            basic_stats_file = 'rmats_basic_stats.json'
            
            existing_files = []
            for filename in expected_files + [summary_file, basic_stats_file]:
                file_path = os.path.join(output_dir, filename)
                if os.path.exists(file_path):
                    existing_files.append(filename)
            
            if len(existing_files) == 0:
                status = 'not_started'
            elif summary_file in existing_files and len(existing_files) >= 5:
                status = 'completed'
            elif basic_stats_file in existing_files:
                status = 'completed_basic'
            else:
                status = 'processing'
            
            return {
                'status': status,
                'existing_files': existing_files,
                'expected_files': expected_files + [summary_file],
                'completion_rate': len(existing_files) / len(expected_files + [summary_file])
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'error': str(e)
            }
    
    # ==================== 私有辅助方法 ====================
    
    @classmethod
    def _read_summary_file(cls, output_dir: str) -> Optional[Dict[str, Any]]:
        """读取处理摘要文件"""
        try:
            summary_file = os.path.join(output_dir, 'frontend_data_summary.json')
            if os.path.exists(summary_file):
                with open(summary_file, 'r', encoding='utf-8') as f:
                    summary_data = json.load(f)
                print(f"📋 读取摘要文件成功，总事件数: {summary_data.get('total_events', 0)}")
                return summary_data
            else:
                print("⚠️ 摘要文件不存在")
                return None
        except Exception as e:
            print(f"⚠️ 读取摘要文件失败: {e}")
            return None
    
    @classmethod
    def _validate_processing_inputs(cls, gtf_file: str, rmats_dir: str) -> Dict[str, Any]:
        """验证处理输入参数"""
        errors = []
        
        # 检查GTF文件
        if gtf_file and not os.path.exists(gtf_file):
            errors.append(f"GTF file not found: {gtf_file}")
        
        # 检查rMATS目录
        if not os.path.exists(rmats_dir):
            errors.append(f"rMATS directory not found: {rmats_dir}")
        elif not os.path.isdir(rmats_dir):
            errors.append(f"rMATS path is not a directory: {rmats_dir}")
        else:
            # 检查必需文件
            missing_files = []
            for filename in Config.REQUIRED_RMATS_FILES:
                file_path = os.path.join(rmats_dir, filename)
                if not os.path.exists(file_path):
                    missing_files.append(filename)
            
            if missing_files:
                errors.append(f"Missing rMATS files: {', '.join(missing_files)}")
        
        return {
            'valid': len(errors) == 0,
            'errors': errors
        }