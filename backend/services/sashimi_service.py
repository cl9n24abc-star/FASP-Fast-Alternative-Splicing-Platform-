"""
SashimiPlot生成服务
处理rmats2sashimiplot相关的所有逻辑
"""

import subprocess
import os
import json
import logging
import tempfile
from typing import Dict, List, Optional, Tuple
from datetime import datetime

class SashimiPlotService:
    """SashimiPlot生成服务类"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def create_compatible_gff(self, original_gff: str, chromosome: str, start: int, end: int) -> str:
        """
        为指定染色体区域创建兼容的GFF文件
        
        Args:
            original_gff: 原始GFF文件路径
            chromosome: 染色体名称（数字格式，如 '21'）
            start: 起始位置
            end: 结束位置
            
        Returns:
            临时兼容GFF文件路径
        """
        try:
            # 创建临时文件
            temp_fd, temp_gff = tempfile.mkstemp(suffix='.gff3', prefix=f'chr{chromosome}_')
            
            self.logger.info(f"创建染色体{chromosome}的兼容GFF文件: {temp_gff}")
            
            # 扩大查找范围，确保包含足够的注释信息
            search_start = max(0, start - 10000)  # 向前扩展10kb
            search_end = end + 10000  # 向后扩展10kb
            
            with os.fdopen(temp_fd, 'w') as f_out:
                with open(original_gff, 'r') as f_in:
                    records_found = 0
                    
                    for line in f_in:
                        # 保留注释行
                        if line.startswith('#'):
                            f_out.write(line)
                            continue
                        
                        # 处理染色体记录
                        if line.startswith(f'chr{chromosome}\t'):
                            fields = line.strip().split('\t')
                            if len(fields) >= 5:
                                try:
                                    record_start = int(fields[3])
                                    record_end = int(fields[4])
                                    
                                    # 检查是否与目标区域重叠
                                    if (record_start <= search_end and record_end >= search_start):
                                        # 转换染色体命名：chr21 -> 21
                                        converted_line = line.replace(f'chr{chromosome}', chromosome, 1)
                                        f_out.write(converted_line)
                                        records_found += 1
                                        
                                except ValueError:
                                    # 如果坐标解析失败，仍然包含这行
                                    converted_line = line.replace(f'chr{chromosome}', chromosome, 1)
                                    f_out.write(converted_line)
                                    records_found += 1
            
            self.logger.info(f"兼容GFF文件创建成功，包含{records_found}条记录")
            
            if records_found == 0:
                self.logger.warning(f"警告：在染色体{chromosome}的区域{search_start}-{search_end}中未找到任何注释记录")
            
            return temp_gff
            
        except Exception as e:
            self.logger.error(f"创建兼容GFF文件失败: {e}")
            # 如果创建失败，返回原始文件
            return original_gff
    
    def validate_sashimi_params(self, params: Dict) -> Tuple[bool, str]:
        """
        验证SashimiPlot生成所需的参数（事件文件模式）
        
        Args:
            params: 前端传来的参数字典
            
        Returns:
            (is_valid, error_message): 验证结果和错误信息
        """
        # 检查必需参数
        required_fields = [
            'output_directory', 'event_type', 'events_file',
            'sample1_label', 'sample2_label'
        ]
        
        for field in required_fields:
            if not params.get(field):
                return False, f"缺少必需参数: {field}"
        
        # 验证事件类型
        valid_event_types = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
        if params['event_type'] not in valid_event_types:
            return False, f"无效的事件类型: {params['event_type']}，支持的类型: {', '.join(valid_event_types)}"
        
        # 验证事件文件存在性
        events_file = params['events_file']
        if not os.path.exists(events_file):
            return False, f"事件文件不存在: {events_file}"
        
        # 如果指定了BAM文件，检查其存在性
        if params.get("file_type") == "bam":
            bam_files = [params.get("sample1_bam"), params.get("sample2_bam")]
            for i, bam_file in enumerate(bam_files, 1):
                if bam_file and not os.path.exists(bam_file):
                    return False, f"样本{i} BAM文件不存在: {bam_file}"
        
        return True, ""
    
    def validate_coordinate_params(self, params: Dict) -> Tuple[bool, str]:
        """
        验证坐标模式SashimiPlot生成所需的参数
        
        Args:
            params: 前端传来的参数字典
            
        Returns:
            (is_valid, error_message): 验证结果和错误信息
        """
        required_fields = [
            'output_directory',
            'coordinate',
            'sample1_label', 'sample2_label',
            'sample1_bam', 'sample2_bam'
        ]
        
        for field in required_fields:
            if not params.get(field):
                return False, f"缺少必需参数: {field}"
        
        # 验证坐标格式 (chromosome:strand:start:end:gff_file)
        coordinate = params['coordinate']
        parts = coordinate.split(':')
        if len(parts) != 5:
            return False, f"坐标格式错误，应为 'chr:strand:start:end:gff_file'，收到: {coordinate}"
        
        chr_name, strand, start, end, gff_file = parts
        
        # 验证链信息
        if strand not in ['+', '-']:
            return False, f"无效的链方向: {strand}，应为 + 或 -"
        
        # 验证坐标
        try:
            start_pos = int(start)
            end_pos = int(end)
            if start_pos >= end_pos:
                return False, f"起始位置应小于结束位置: {start} >= {end}"
        except ValueError:
            return False, f"坐标必须为整数: start={start}, end={end}"
        
        # 验证GFF文件存在性
        if not os.path.exists(gff_file):
            return False, f"GFF文件不存在: {gff_file}"
        
        # 验证BAM文件存在性
        bam_files = [params.get("sample1_bam"), params.get("sample2_bam")]
        for i, bam_file in enumerate(bam_files, 1):
            if bam_file and not os.path.exists(bam_file):
                return False, f"样本{i} BAM文件不存在: {bam_file}"
        
        return True, ""
    
    def build_rmats2sashimi_command(self, params: Dict) -> List[str]:
        """
        构建rmats2sashimiplot命令行参数（事件文件模式）
        
        Args:
            params: 参数字典
            
        Returns:
            完整的命令参数列表
        """
        cmd = ["rmats2sashimiplot"]
        
        # 必需参数
        cmd.extend(["-o", params["output_directory"]])
        cmd.extend(["--event-type", params["event_type"]])
        cmd.extend(["-e", params["events_file"]])
        
        # 样本标签
        cmd.extend(["--l1", params["sample1_label"]])
        cmd.extend(["--l2", params["sample2_label"]])
        
        # BAM文件参数
        if params.get("file_type") == "bam":
            if params.get("sample1_bam"):
                cmd.extend(["--b1", params["sample1_bam"]])
            if params.get("sample2_bam"):
                cmd.extend(["--b2", params["sample2_bam"]])
        
        # 可选参数映射
        optional_params = {
            "--exon_s": params.get("exon_scale"),
            "--intron_s": params.get("intron_scale"),
            "--min-counts": params.get("min_counts"),
            "--color": params.get("color"),
            "--font-size": params.get("font_size"),
            "--fig-height": params.get("figure_height"),
            "--fig-width": params.get("figure_width"),
            "--group-info": params.get("group_info")
        }
        
        # 添加有效的可选参数
        for param, value in optional_params.items():
            if value is not None and str(value).lower() not in ["none", ""]:
                cmd.extend([param, str(value)])
        
        return cmd
    
    def build_coordinate_command(self, params: Dict) -> List[str]:
        """
        构建坐标模式的rmats2sashimiplot命令行参数
        
        Args:
            params: 参数字典
            
        Returns:
            完整的命令参数列表
        """
        cmd = ["rmats2sashimiplot"]
        
        # 必需参数
        cmd.extend(["-o", params["output_directory"]])
        cmd.extend(["-c", params["coordinate"]])  # 坐标参数
        
        # 样本标签
        cmd.extend(["--l1", params["sample1_label"]])
        cmd.extend(["--l2", params["sample2_label"]])
        
        # BAM文件参数
        cmd.extend(["--b1", params["sample1_bam"]])
        cmd.extend(["--b2", params["sample2_bam"]])
        
        # 修改：更宽松的参数设置
        optional_params = {
            "--exon_s": params.get("exon_scale", 1),
            "--intron_s": params.get("intron_scale", 1),
            "--min-counts": params.get("min_counts", 0),  # 改为0，最宽松
            "--color": params.get("color"),
            "--font-size": params.get("font_size", 10),
            "--fig-height": params.get("figure_height", 7),
            "--fig-width": params.get("figure_width", 8),
            "--group-info": params.get("group_info")
        }
        
        # 添加有效的可选参数
        for param, value in optional_params.items():
            if value is not None and str(value).lower() not in ["none", ""]:
                cmd.extend([param, str(value)])
        
        return cmd
    
    def generate_sashimi_plot(self, params: Dict) -> Dict:
        """
        生成SashimiPlot的主要方法（事件文件模式）
        
        Args:
            params: 前端传来的完整参数字典
            
        Returns:
            包含执行结果的字典
        """
        start_time = datetime.now()
        
        # 1. 参数验证
        self.logger.info("开始验证SashimiPlot参数...")
        is_valid, error_msg = self.validate_sashimi_params(params)
        if not is_valid:
            self.logger.error(f"参数验证失败: {error_msg}")
            return {
                "success": False,
                "message": "参数验证失败",
                "error": error_msg,
                "timestamp": start_time.isoformat()
            }
        
        # 2. 构建命令
        self.logger.info("构建rmats2sashimiplot命令...")
        cmd = self.build_rmats2sashimi_command(params)
        
        return self._execute_sashimi_command(cmd, params, start_time)
    
    def generate_sashimi_by_coordinate(self, params: Dict) -> Dict:
        """
        使用坐标模式生成SashimiPlot的主要方法
        
        Args:
            params: 前端传来的完整参数字典
            
        Returns:
            包含执行结果的字典
        """
        start_time = datetime.now()
        temp_gff_file = None
        
        try:
            # 设置默认输出目录
            if not params.get('output_directory'):
                params['output_directory'] = '../frontend/public/sashimi_output'
                self.logger.info(f"使用默认输出目录: {params['output_directory']}")
            
            # 解析坐标信息
            coordinate = params.get('coordinate', '')
            parts = coordinate.split(':')
            if len(parts) == 5:
                chr_name, strand, start_str, end_str, original_gff = parts
                
                # 去除可能的chr前缀，确保使用数字格式
                if chr_name.startswith('chr'):
                    chr_name = chr_name[3:]
                
                try:
                    start_pos = int(start_str)
                    end_pos = int(end_str)
                    
                    self.logger.info(f"解析坐标: 染色体{chr_name}, 位置{start_pos}-{end_pos}")
                    
                    # 创建兼容的GFF文件
                    if os.path.exists(original_gff):
                        temp_gff_file = self.create_compatible_gff(original_gff, chr_name, start_pos, end_pos)
                        
                        # 更新coordinate参数，使用新的GFF文件
                        new_coordinate = f"{chr_name}:{strand}:{start_pos}:{end_pos}:{temp_gff_file}"
                        params['coordinate'] = new_coordinate
                        
                        self.logger.info(f"已创建兼容GFF文件，更新坐标: {new_coordinate}")
                    else:
                        self.logger.warning(f"原始GFF文件不存在: {original_gff}")
                        
                except ValueError as e:
                    self.logger.error(f"坐标解析失败: {e}")
            
            # 1. 参数验证
            self.logger.info("开始验证坐标模式SashimiPlot参数...")
            is_valid, error_msg = self.validate_coordinate_params(params)
            if not is_valid:
                self.logger.error(f"参数验证失败: {error_msg}")
                return {
                    "success": False,
                    "message": "坐标模式参数验证失败",
                    "error": error_msg,
                    "timestamp": start_time.isoformat()
                }
            
            # 2. 构建命令
            self.logger.info("构建坐标模式rmats2sashimiplot命令...")
            cmd = self.build_coordinate_command(params)
            
            # 3. 执行命令
            result = self._execute_sashimi_command(cmd, params, start_time)
            
            return result
            
        except Exception as e:
            self.logger.error(f"生成SashimiPlot时发生错误: {e}")
            return {
                "success": False,
                "message": "生成SashimiPlot时发生错误",
                "error": str(e),
                "timestamp": start_time.isoformat()
            }
        finally:
            # 清理临时文件
            if temp_gff_file and temp_gff_file != params.get('coordinate', '').split(':')[-1]:
                try:
                    if os.path.exists(temp_gff_file):
                        os.unlink(temp_gff_file)
                        self.logger.info(f"已清理临时GFF文件: {temp_gff_file}")
                except Exception as e:
                    self.logger.warning(f"清理临时文件失败: {e}")
    
    def generate_coordinate_sashimi_plot(self, params: Dict) -> Dict:
        """
        生成SashimiPlot的坐标模式方法（别名方法，保持兼容性）
        
        Args:
            params: 前端传来的完整参数字典
            
        Returns:
            包含执行结果的字典
        """
        # 直接调用新的方法名
        return self.generate_sashimi_by_coordinate(params)
    
    def _execute_sashimi_command(self, cmd: List[str], params: Dict, start_time: datetime) -> Dict:
        """
        执行rmats2sashimiplot命令的通用方法
        
        Args:
            cmd: 完整的命令参数列表
            params: 参数字典
            start_time: 开始时间
            
        Returns:
            包含执行结果的字典
        """
        # 创建输出目录
        output_dir = params["output_directory"]
        try:
            abs_output_dir = os.path.abspath(output_dir)
            os.makedirs(abs_output_dir, exist_ok=True)
            self.logger.info(f"输出目录已准备: {abs_output_dir}")
        except Exception as e:
            return {
                "success": False,
                "message": "无法创建输出目录",
                "error": str(e)
            }
        
        # 执行rmats2sashimiplot命令
        try:
            cmd_str = ' '.join(cmd)
            self.logger.info(f"执行命令: {cmd_str}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=600  # 10分钟超时
            )
            
            # 查找生成的文件
            generated_files = self._find_generated_files(abs_output_dir)
            metadata = self.get_sashimi_metadata(abs_output_dir)
            
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            self.logger.info(f"SashimiPlot生成成功，耗时: {duration:.2f}秒")
            
            return {
                "success": True,
                "message": "SashimiPlot生成成功",
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": cmd_str,
                "output_directory": abs_output_dir,
                "generated_files": generated_files,
                "metadata": metadata,
                "duration_seconds": duration,
                "start_time": start_time.isoformat(),
                "end_time": end_time.isoformat()
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"rmats2sashimiplot执行失败: {e}")
            return {
                "success": False,
                "message": "rmats2sashimiplot执行失败",
                "error": e.stderr,
                "stdout": e.stdout,
                "command": ' '.join(cmd),
                "return_code": e.returncode
            }
        except subprocess.TimeoutExpired as e:
            self.logger.error("rmats2sashimiplot执行超时")
            return {
                "success": False,
                "message": "rmats2sashimiplot执行超时",
                "error": "命令执行超过10分钟",
                "command": ' '.join(cmd)
            }
        except FileNotFoundError:
            self.logger.error("rmats2sashimiplot命令未找到")
            return {
                "success": False,
                "message": "rmats2sashimiplot命令未找到",
                "error": "请确保rmats2sashimiplot已正确安装并在PATH中",
                "suggestion": "可以尝试: conda install -c bioconda rmats2sashimiplot"
            }
        except Exception as e:
            self.logger.error(f"生成SashimiPlot时发生未知错误: {e}")
            return {
                "success": False,
                "message": "生成SashimiPlot时发生未知错误",
                "error": str(e)
            }
    
    def _find_generated_files(self, output_dir: str) -> List[str]:
        """
        查找输出目录中生成的图片文件
        
        Args:
            output_dir: 输出目录路径
            
        Returns:
            生成的文件路径列表
        """
        generated_files = []
        if not os.path.exists(output_dir):
            return generated_files
        
        # SashimiPlot常见的输出文件格式
        image_extensions = ['.pdf', '.png', '.svg', '.eps']
        
        try:
            # 检查主目录和Sashimi_plot子目录
            search_dirs = [output_dir]
            sashimi_plot_dir = os.path.join(output_dir, 'Sashimi_plot')
            if os.path.exists(sashimi_plot_dir):
                search_dirs.append(sashimi_plot_dir)
            
            for search_dir in search_dirs:
                for file in os.listdir(search_dir):
                    file_path = os.path.join(search_dir, file)
                    if os.path.isfile(file_path):
                        # 检查文件扩展名
                        for ext in image_extensions:
                            if file.lower().endswith(ext):
                                generated_files.append(file_path)
                                break
        except Exception as e:
            self.logger.error(f"扫描输出文件时出错: {e}")
        
        return generated_files
    
    def get_sashimi_metadata(self, output_dir: str) -> Dict:
        """
        获取SashimiPlot生成结果的元数据
        
        Args:
            output_dir: 输出目录路径
            
        Returns:
            包含元数据的字典
        """
        metadata = {
            "plot_count": 0,
            "files": [],
            "total_size_bytes": 0,
            "total_size_mb": 0,
            "file_types": {},
            "directories": []
        }
        
        if not os.path.exists(output_dir):
            return metadata
        
        try:
            # 检查主目录和Sashimi_plot子目录
            search_dirs = [output_dir]
            sashimi_plot_dir = os.path.join(output_dir, 'Sashimi_plot')
            if os.path.exists(sashimi_plot_dir):
                search_dirs.append(sashimi_plot_dir)
                metadata["directories"].append("Sashimi_plot")
            
            for search_dir in search_dirs:
                for file in os.listdir(search_dir):
                    file_path = os.path.join(search_dir, file)
                    
                    if os.path.isfile(file_path):
                        # 检查是否是图片文件
                        if any(file.lower().endswith(ext) for ext in ['.pdf', '.png', '.svg', '.eps']):
                            file_size = os.path.getsize(file_path)
                            file_type = file.split('.')[-1].lower()
                            
                            # 文件信息
                            file_info = {
                                "name": file,
                                "path": file_path,
                                "directory": os.path.basename(search_dir),
                                "size_bytes": file_size,
                                "size_mb": round(file_size / (1024 * 1024), 2),
                                "type": file_type,
                                "modified_time": datetime.fromtimestamp(os.path.getmtime(file_path)).isoformat()
                            }
                            
                            metadata["files"].append(file_info)
                            metadata["total_size_bytes"] += file_size
                            metadata["plot_count"] += 1
                            
                            # 统计文件类型
                            if file_type in metadata["file_types"]:
                                metadata["file_types"][file_type] += 1
                            else:
                                metadata["file_types"][file_type] = 1
            
            # 转换总大小为MB
            metadata["total_size_mb"] = round(metadata["total_size_bytes"] / (1024 * 1024), 2)
            
        except Exception as e:
            self.logger.error(f"获取元数据时出错: {e}")
            metadata["error"] = str(e)
        
        return metadata
    
    def check_rmats2sashimi_availability(self) -> Dict:
        """
        检查rmats2sashimiplot工具是否可用
        
        Returns:
            检查结果字典
        """
        try:
            result = subprocess.run(
                ["rmats2sashimiplot", "--help"],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            return {
                "available": True,
                "version_info": result.stdout[:200] if result.stdout else "工具可用",
                "message": "rmats2sashimiplot工具检查通过"
            }
            
        except FileNotFoundError:
            return {
                "available": False,
                "message": "rmats2sashimiplot工具未找到",
                "suggestion": "请安装: conda install -c bioconda rmats2sashimiplot"
            }
        except Exception as e:
            return {
                "available": False,
                "message": f"检查rmats2sashimiplot时出错: {str(e)}"
            }
    
    def clean_output_directory(self, output_dir: str, keep_recent: int = 5) -> Dict:
        """
        清理输出目录中的旧文件（可选功能）
        
        Args:
            output_dir: 要清理的目录
            keep_recent: 保留最近的文件数量
            
        Returns:
            清理结果
        """
        try:
            if not os.path.exists(output_dir):
                return {"success": True, "message": "目录不存在，无需清理"}
            
            # 获取所有图片文件并按修改时间排序
            files = []
            search_dirs = [output_dir]
            sashimi_plot_dir = os.path.join(output_dir, 'Sashimi_plot')
            if os.path.exists(sashimi_plot_dir):
                search_dirs.append(sashimi_plot_dir)
            
            for search_dir in search_dirs:
                for file in os.listdir(search_dir):
                    file_path = os.path.join(search_dir, file)
                    if os.path.isfile(file_path) and any(file.lower().endswith(ext) for ext in ['.pdf', '.png', '.svg']):
                        files.append((file_path, os.path.getmtime(file_path)))
            
            # 按修改时间排序（新的在前）
            files.sort(key=lambda x: x[1], reverse=True)
            
            # 删除多余的文件
            deleted_count = 0
            if len(files) > keep_recent:
                for file_path, _ in files[keep_recent:]:
                    try:
                        os.remove(file_path)
                        deleted_count += 1
                    except Exception as e:
                        self.logger.warning(f"删除文件失败 {file_path}: {e}")
            
            return {
                "success": True,
                "message": f"清理完成，删除了 {deleted_count} 个文件",
                "deleted_count": deleted_count,
                "remaining_count": len(files) - deleted_count
            }
            
        except Exception as e:
            return {
                "success": False,
                "message": f"清理目录时出错: {str(e)}"
            }