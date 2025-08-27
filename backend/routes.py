from flask import request, jsonify
import json
import os
import subprocess
import tempfile
from pathlib import Path
from config import Config
import sys 
from services.file_cleanup import cleanup_json_only, cleanup_public_files

def register_routes(app):
    """注册所有API路由到Flask应用"""
    
    # ==================== 路径验证相关路由 ====================
    
    @staticmethod
    def validate_bam_list_file(path):
        """验证BAM列表文件"""
        try:
            import os
            
            # 检查文件是否存在
            if not os.path.exists(path):
                return {
                    'valid': False,
                    'error': f'BAM list file does not exist: {path}',
                    'exists': False
                }
            
            # 检查是否是文件
            if not os.path.isfile(path):
                return {
                    'valid': False,
                    'error': f'Path is not a file: {path}',
                    'exists': True
                }
            
            # 检查文件扩展名
            if not path.lower().endswith('.txt'):
                return {
                    'valid': False,
                    'error': f'BAM list file must be a .txt file: {path}',
                    'exists': True
                }
            
            # 读取并验证文件内容
            bam_files = []
            missing_bams = []
            
            with open(path, 'r', encoding='utf-8') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    line = line.strip()
                    
                    # 跳过空行和注释
                    if not line or line.startswith('#'):
                        continue
                    
                    # 检查BAM文件是否存在
                    if os.path.exists(line):
                        bam_files.append(line)
                    else:
                        missing_bams.append(f"Line {line_num}: {line}")
            
            # 计算文件大小
            file_size = os.path.getsize(path)
            size_mb = round(file_size / (1024 * 1024), 2)
            
            # 构建返回结果
            if not bam_files:
                return {
                    'valid': False,
                    'error': 'No valid BAM files found in the list',
                    'exists': True,
                    'size_mb': size_mb,
                    'bam_count': 0,
                    'missing_bams': missing_bams
                }
            
            # 如果有缺失的BAM文件，给出警告但仍然标记为有效
            result = {
                'valid': True,
                'message': f'BAM list file contains {len(bam_files)} valid BAM files',
                'exists': True,
                'size_mb': size_mb,
                'bam_count': len(bam_files),
                'total_lines': line_num,
                'sample_paths': bam_files[:3]  # 显示前3个路径作为样本
            }
            
            if missing_bams:
                result['warnings'] = f'{len(missing_bams)} BAM files not found'
                result['missing_bams'] = missing_bams[:5]  # 显示前5个缺失的文件
            
            return result
            
        except PermissionError:
            return {
                'valid': False,
                'error': f'Permission denied accessing file: {path}',
                'exists': True
            }
        except UnicodeDecodeError:
            return {
                'valid': False,
                'error': f'File encoding error - file may be binary or corrupted: {path}',
                'exists': True
            }
        except Exception as e:
            return {
                'valid': False,
                'error': f'Error validating BAM list file: {str(e)}',
                'exists': os.path.exists(path) if path else False
            }
    @app.route('/api/cleanup/generated-files', methods=['POST'])
    def cleanup_files():
        print(f"Config.get_frontend_public_path(): {Config.get_frontend_public_path()}")
        print(f"Path exists: {os.path.exists(Config.get_frontend_public_path())}")
        try:
            # 方法1: 使用现有的便捷函数
            public_path = Config.get_frontend_public_path()
            result = cleanup_public_files(public_path)  # 使用已导入的函数
            
            return jsonify(result)
            
        except Exception as e:
            import traceback
            print(f"Cleanup error: {e}")
            print(traceback.format_exc())
            return jsonify({
                "success": False,
                "message": f"Cleanup failed: {str(e)}",
                "errors": [str(e)]
            }), 500
    @app.route('/api/validate-path', methods=['POST', 'OPTIONS'])
    def validate_path():
        """路径验证API"""
        
        # 处理CORS预检请求
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        try:
            data = request.json
            path = data.get('path', '').strip()
            path_type = data.get('type', '')
            
            if not path:
                return jsonify({
                    'valid': False,
                    'error': 'Path cannot be empty'
                }), 400
            
            print(f"🔍 验证路径: {path} (类型: {path_type})")
            
            # 导入验证服务并调用
            from services.path_validator import PathValidator
            
            if path_type == 'rmats_dir':
                result = PathValidator.validate_rmats_directory(path)
            elif path_type == 'output_dir':
                result = PathValidator.validate_output_directory(path)
            elif path_type == 'gtf_file':
                result = PathValidator.validate_gtf_file(path)
            elif path_type == 'bam_file':
                result = PathValidator.validate_bam_file(path)
            else:
                return jsonify({
                    'valid': False,
                    'error': f'Unknown path type: {path_type}'
                }), 400
            
            print(f"✅ 验证结果: {result}")
            return jsonify(result), 200
            
        except Exception as e:
            print(f"❌ 验证出错: {e}")
            return jsonify({
                'valid': False,
                'error': f'Validation error: {str(e)}'
            }), 500
    
    @app.route('/api/validate-bam-lists', methods=['POST', 'OPTIONS'])
    def validate_bam_lists():
        """同时验证两个BAM列表文件"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        try:
            data = request.json
            sample1_path = data.get('sample1_path', '').strip()
            sample2_path = data.get('sample2_path', '').strip()
            
            print(f"🧬 验证双BAM列表文件:")
            print(f"  Sample 1: {sample1_path}")
            print(f"  Sample 2: {sample2_path}")
            
            # 验证两个路径都不为空
            if not sample1_path or not sample2_path:
                return jsonify({
                    'success': False,
                    'error': 'Both Sample 1 and Sample 2 BAM list paths are required',
                    'sample1_result': {
                        'valid': False,
                        'error': 'Sample 1 path is required' if not sample1_path else 'Valid path provided',
                        'exists': False
                    },
                    'sample2_result': {
                        'valid': False,
                        'error': 'Sample 2 path is required' if not sample2_path else 'Valid path provided',
                        'exists': False
                    }
                }), 400
            
            # 分别验证两个BAM列表文件
            sample1_result = validate_bam_list_file(sample1_path)
            sample2_result = validate_bam_list_file(sample2_path)
            
            # 构建响应
            response_data = {
                'success': sample1_result.get('valid', False) and sample2_result.get('valid', False),
                'message': 'BAM lists validation completed',
                'sample1_result': sample1_result,
                'sample2_result': sample2_result,
                'summary': {
                    'sample1_valid': sample1_result.get('valid', False),
                    'sample2_valid': sample2_result.get('valid', False),
                    'sample1_bam_count': sample1_result.get('bam_count', 0),
                    'sample2_bam_count': sample2_result.get('bam_count', 0),
                    'total_bam_files': sample1_result.get('bam_count', 0) + sample2_result.get('bam_count', 0)
                }
            }
            
            if response_data['success']:
                print("✅ 双BAM列表验证成功")
                response_data['message'] = f"Both BAM lists validated successfully. Total BAM files: {response_data['summary']['total_bam_files']}"
            else:
                print("❌ 双BAM列表验证失败")
                errors = []
                if not sample1_result.get('valid', False):
                    errors.append(f"Sample 1: {sample1_result.get('error', 'Unknown error')}")
                if not sample2_result.get('valid', False):
                    errors.append(f"Sample 2: {sample2_result.get('error', 'Unknown error')}")
                response_data['message'] = '; '.join(errors)
            
            return jsonify(response_data), 200
            
        except Exception as e:
            print(f"❌ 双BAM列表验证出错: {e}")
            return jsonify({
                'success': False,
                'error': str(e),
                'message': f'Validation failed: {str(e)}',
                'sample1_result': {
                    'valid': False,
                    'error': str(e),
                    'exists': False
                },
                'sample2_result': {
                    'valid': False,
                    'error': str(e),
                    'exists': False
                }
            }), 500

    # ==================== BAM分析相关路由 ====================
        
    @app.route('/api/analyze-bam-files', methods=['POST', 'OPTIONS'])
    def analyze_bam_files():
        """分析BAM文件并生成可视化数据"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        try:
            data = request.json
            sample1_bam_list = data.get('sample1_bam_list', '').strip()
            sample2_bam_list = data.get('sample2_bam_list', '').strip()
            chromosome = data.get('chromosome', 'chr1')
            
            print(f"🧬 收到BAM分析请求:")
            print(f"  Sample 1 BAM列表: {sample1_bam_list}")
            print(f"  Sample 2 BAM列表: {sample2_bam_list}")
            print(f"  分析染色体: {chromosome}")
            
            # 验证输入参数
            if not sample1_bam_list or not sample2_bam_list:
                return jsonify({
                    'success': False,
                    'message': 'Both sample1_bam_list and sample2_bam_list are required',
                    'error': 'Missing BAM list files'
                }), 400
            
            # 验证文件存在
            if not os.path.exists(sample1_bam_list):
                return jsonify({
                    'success': False,
                    'message': f'Sample 1 BAM list file not found: {sample1_bam_list}',
                    'error': 'File not found'
                }), 400
            
            if not os.path.exists(sample2_bam_list):
                return jsonify({
                    'success': False,
                    'message': f'Sample 2 BAM list file not found: {sample2_bam_list}',
                    'error': 'File not found'
                }), 400
            
            # 确定输出文件路径
            output_file = os.path.join(Config.get_frontend_public_path(), 'bam_analysis_result.json')
            
            # 找到BAM分析脚本 - 修复路径查找逻辑
            current_dir = os.path.dirname(os.path.abspath(__file__))
            
            # 扩展搜索路径，包含services目录
            possible_paths = [
                os.path.join(current_dir, 'services', 'bam_analyzer.py'),  # 最可能的位置
                os.path.join(current_dir, '..', 'services', 'bam_analyzer.py'),
                os.path.join(current_dir, 'bam_analyzer.py'),
                os.path.join(current_dir, '..', 'bam_analyzer.py'),
                os.path.join(current_dir, '..', 'scripts', 'bam_analyzer.py'),
                os.path.join(os.getcwd(), 'services', 'bam_analyzer.py'),
                os.path.join(os.getcwd(), 'bam_analyzer.py'),
                os.path.join(os.getcwd(), 'scripts', 'bam_analyzer.py')
            ]
            
            analyzer_script = None
            for path in possible_paths:
                if os.path.exists(path):
                    analyzer_script = path
                    break
            
            if not analyzer_script:
                print(f"❌ BAM分析脚本未找到，搜索路径:")
                for path in possible_paths:
                    print(f"  - {path} {'✓' if os.path.exists(path) else '✗'}")
                
                return jsonify({
                    'success': False,
                    'message': 'BAM analyzer script not found',
                    'error': f'Script not found in any expected locations',
                    'searched_paths': possible_paths
                }), 500
            
            print(f"📊 使用分析脚本: {analyzer_script}")
            print(f"📁 输出文件: {output_file}")
            
            # 确保使用正确的Python环境
            python_executable = sys.executable or 'python3'
            
            # 构建命令 - 使用旧格式（位置参数）
            cmd = [
                python_executable, analyzer_script,
                sample1_bam_list,
                sample2_bam_list
            ]
            
            print(f"🔧 执行命令: {' '.join(cmd)}")
            
            # 执行BAM分析
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=1800,  # 30分钟超时
                    cwd=os.path.dirname(analyzer_script)
                )
                
                print(f"🔄 命令执行完成，返回码: {result.returncode}")
                print(f"📤 标准输出: {result.stdout}")
                if result.stderr:
                    print(f"📥 错误输出: {result.stderr}")
                
                if result.returncode == 0:
                    print("✅ BAM分析成功完成")
                    
                    # 验证输出文件是否生成
                    if os.path.exists(output_file):
                        # 读取并返回结果摘要
                        try:
                            with open(output_file, 'r', encoding='utf-8') as f:
                                analysis_result = json.load(f)
                            
                            return jsonify({
                                'success': True,
                                'message': 'BAM analysis completed successfully',
                                'output_file': output_file,
                                'analysis_info': analysis_result.get('analysis_info', {}),
                                'summary': {
                                    'sample1_mapping_rate': analysis_result.get('group1', {}).get('mappingRate', 0),
                                    'sample2_mapping_rate': analysis_result.get('group2', {}).get('mappingRate', 0),
                                    'chromosome_analyzed': analysis_result.get('analysis_info', {}).get('chromosome_analyzed', chromosome)
                                }
                            }), 200
                        except json.JSONDecodeError as e:
                            return jsonify({
                                'success': False,
                                'message': 'Analysis completed but output file is corrupted',
                                'error': f'JSON decode error: {str(e)}'
                            }), 500
                    else:
                        return jsonify({
                            'success': False,
                            'message': 'Analysis completed but output file not found',
                            'error': f'Output file not generated: {output_file}',
                            'stdout': result.stdout,
                            'stderr': result.stderr
                        }), 500
                        
                else:
                    print(f"❌ BAM分析失败 (返回码: {result.returncode})")
                    
                    return jsonify({
                        'success': False,
                        'message': 'BAM analysis failed',
                        'error': result.stderr or result.stdout or 'Unknown error',
                        'return_code': result.returncode,
                        'stdout': result.stdout,
                        'stderr': result.stderr
                    }), 500
                    
            except subprocess.TimeoutExpired:
                return jsonify({
                    'success': False,
                    'message': 'BAM analysis timed out',
                    'error': 'Analysis took longer than 30 minutes'
                }), 500
            except FileNotFoundError as e:
                return jsonify({
                    'success': False,
                    'message': 'Python executable or script not found',
                    'error': f'FileNotFoundError: {str(e)}',
                    'python_executable': python_executable,
                    'script_path': analyzer_script
                }), 500
                
        except Exception as e:
            import traceback
            print(f"❌ BAM分析处理过程中发生错误: {e}")
            print(traceback.format_exc())
            return jsonify({
                'success': False,
                'message': f'BAM analysis processing error: {str(e)}',
                'error': str(e),
                'traceback': traceback.format_exc()
            }), 500
            
    @app.route('/api/create-directory', methods=['POST', 'OPTIONS'])
    def create_directory():
        """创建目录API"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        try:
            data = request.json
            path = data.get('path', '').strip()
            
            if not path:
                return jsonify({
                    'success': False,
                    'error': 'Path cannot be empty'
                }), 400
            
            print(f"📁 创建目录: {path}")
            
            # 导入服务并调用
            from services.path_validator import PathValidator
            result = PathValidator.create_directory(path)
            
            return jsonify(result), 200 if result['success'] else 500
            
        except Exception as e:
            return jsonify({
                'success': False,
                'error': str(e)
            }), 500
    
    # ==================== rMATS处理相关路由 ====================
    @app.route('/api/process-rmats', methods=['POST', 'OPTIONS'])
    def process_rmats():
        """处理rMATS数据并生成前端格式"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        try:
            data = request.json
            
            print("\n" + "="*50)
            print("收到rMATS处理请求")
            print("="*50)
            for key, value in data.items():
                print(f"{key}: {value}")
            print("-"*50)
            
            # 提取5个参数
            rmats_dir = data.get('rmats_dir', '').strip()
            output_dir = data.get('output_dir', '').strip()
            sample1_bam_list = data.get('sample1_bam_list', '').strip()
            sample2_bam_list = data.get('sample2_bam_list', '').strip()
            gtf_file = data.get('gtf_file', '').strip()  # 可选
            sashimi_config = data.get('sashimi_config')  # 可选
            
            # 新增：存储 Sashimi 配置到应用配置中
            if sashimi_config:
                app.config['SASHIMI_CONFIG'] = sashimi_config
                print(f"📊 Sashimi配置已保存: {sashimi_config}")
                print(f"   - BAM文件: {sashimi_config.get('bam_files', {})}")
                print(f"   - GFF文件: {sashimi_config.get('gff_file', 'N/A')}")
                print(f"   - 显示选项: {sashimi_config.get('display_options', {})}")
            else:
                # 清除之前的配置
                if 'SASHIMI_CONFIG' in app.config:
                    del app.config['SASHIMI_CONFIG']
                print("📊 未提供Sashimi配置或已清除之前的配置")
            
            # 验证必选参数（3个必选项）
            missing_params = []
            if not rmats_dir:
                missing_params.append('rmats_dir')
            if not output_dir:
                missing_params.append('output_dir')
            if not sample1_bam_list:
                missing_params.append('sample1_bam_list')
            if not sample2_bam_list:
                missing_params.append('sample2_bam_list')
            
            if missing_params:
                return jsonify({
                    'success': False,
                    'message': f'缺少必选参数: {", ".join(missing_params)}',
                    'error': f'Required parameters missing: {", ".join(missing_params)}'
                }), 400
            
            # 验证必选参数的路径是否存在
            if not os.path.exists(rmats_dir):
                return jsonify({
                    'success': False,
                    'message': f'rMATS目录不存在: {rmats_dir}',
                    'error': f'rMATS directory not found: {rmats_dir}'
                }), 400
            
            if not os.path.exists(sample1_bam_list):
                return jsonify({
                    'success': False,
                    'message': f'Sample 1 BAM列表文件不存在: {sample1_bam_list}',
                    'error': f'Sample 1 BAM list file not found: {sample1_bam_list}'
                }), 400
            
            if not os.path.exists(sample2_bam_list):
                return jsonify({
                    'success': False,
                    'message': f'Sample 2 BAM列表文件不存在: {sample2_bam_list}',
                    'error': f'Sample 2 BAM list file not found: {sample2_bam_list}'
                }), 400
            
            # 创建输出目录（如果不存在）
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                    print(f"创建输出目录: {output_dir}")
                except Exception as e:
                    return jsonify({
                        'success': False,
                        'message': f'无法创建输出目录: {output_dir}',
                        'error': f'Failed to create output directory: {str(e)}'
                    }), 500
            
            print("\n" + "="*30 + " 开始处理 " + "="*30)
            results = {
                'rmats_processing': None,
                'bam_analysis': None,
                'gtf_processing': None,
                'sashimi_processing': None,
                'processing_summary': {
                    'modules_executed': [],
                    'modules_skipped': [],
                    'total_processing_time': 0
                }
            }
            
            import time
            start_time = time.time()
            
            # 1. 必选：处理rMATS数据（使用单分析模式）
            print("步骤1: 处理rMATS数据...")
            try:
                if gtf_file:
                    print(f"使用GTF文件: {gtf_file}")
                    # 验证GTF文件
                    if not os.path.exists(gtf_file):
                        return jsonify({
                            'success': False,
                            'message': f'GTF文件不存在: {gtf_file}',
                            'error': f'GTF file not found: {gtf_file}'
                        }), 400
                    
                    if not os.path.isfile(gtf_file):
                        return jsonify({
                            'success': False,
                            'message': f'GTF路径不是文件: {gtf_file}',
                            'error': f'GTF path is not a file: {gtf_file}'
                        }), 400
                    
                    print("GTF文件验证通过，将进行完整处理")
                else:
                    print("未提供GTF文件，将跳过基因注释处理")
                    gtf_file = None
                
                from services.rmats_processor import RMATSProcessor
                rmats_result = RMATSProcessor.process_single_analysis(
                    gtf_file=gtf_file,
                    rmats_dir=rmats_dir,
                    output_dir=Config.get_frontend_public_path()
                )
                
                results['rmats_processing'] = rmats_result
                results['processing_summary']['modules_executed'].append('rMATS_processing')
                
                if rmats_result.get('success'):
                    print("rMATS数据处理完成")
                else:
                    print(f"rMATS数据处理失败: {rmats_result.get('message')}")
                    
            except Exception as e:
                import traceback
                print(f"rMATS处理出错: {e}")
                print(traceback.format_exc())
                results['rmats_processing'] = {
                    'success': False,
                    'message': f'rMATS处理失败: {str(e)}',
                    'error': str(e)
                }
            
            # 2. 必选：BAM文件分析
            print("\n步骤2: BAM文件分析...")
            try:
                # 直接调用BAM分析逻辑（复用现有代码）
                output_file = os.path.join(Config.get_frontend_public_path(), 'bam_analysis_result.json')
                
                # 找到BAM分析脚本
                current_dir = os.path.dirname(os.path.abspath(__file__))
                possible_paths = [
                    os.path.join(current_dir, 'services', 'bam_analyzer.py'),
                    os.path.join(current_dir, '..', 'services', 'bam_analyzer.py'),
                    os.path.join(current_dir, 'bam_analyzer.py'),
                    os.path.join(os.getcwd(), 'services', 'bam_analyzer.py'),
                ]
                
                analyzer_script = None
                for path in possible_paths:
                    if os.path.exists(path):
                        analyzer_script = path
                        break
                
                if analyzer_script:
                    print(f"使用BAM分析脚本: {analyzer_script}")
                    
                    python_executable = sys.executable or 'python3'
                    cmd = [python_executable, analyzer_script, sample1_bam_list, sample2_bam_list]
                    
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=1800,
                        cwd=os.path.dirname(analyzer_script)
                    )
                    
                    if result.returncode == 0 and os.path.exists(output_file):
                        with open(output_file, 'r', encoding='utf-8') as f:
                            analysis_result = json.load(f)
                        
                        results['bam_analysis'] = {
                            'success': True,
                            'message': 'BAM分析完成',
                            'output_file': output_file,
                            'summary': {
                                'sample1_mapping_rate': analysis_result.get('group1', {}).get('mappingRate', 0),
                                'sample2_mapping_rate': analysis_result.get('group2', {}).get('mappingRate', 0)
                            }
                        }
                        results['processing_summary']['modules_executed'].append('BAM_analysis')
                        print("BAM文件分析完成")
                    else:
                        results['bam_analysis'] = {
                            'success': False,
                            'message': 'BAM分析失败',
                            'error': result.stderr or result.stdout or 'Unknown error'
                        }
                        print(f"BAM分析失败: {result.stderr}")
                else:
                    results['bam_analysis'] = {
                        'success': False,
                        'message': 'BAM分析脚本未找到',
                        'error': 'BAM analyzer script not found'
                    }
                    print("BAM分析脚本未找到")
                    
            except Exception as e:
                print(f"BAM分析出错: {e}")
                results['bam_analysis'] = {
                    'success': False,
                    'message': f'BAM分析失败: {str(e)}',
                    'error': str(e)
                }
            
            # 3. 可选：GTF处理（如果在rMATS处理中已包含，这里可以记录状态）
            if gtf_file:
                results['gtf_processing'] = {
                    'success': True,
                    'message': 'GTF文件已在rMATS处理中使用',
                    'gtf_file': gtf_file
                }
                results['processing_summary']['modules_executed'].append('GTF_annotation')
            else:
                results['gtf_processing'] = {
                    'success': False,
                    'message': '未提供GTF文件，跳过基因注释',
                    'skipped': True
                }
                results['processing_summary']['modules_skipped'].append('GTF_annotation')
            
            # 4. Sashimi配置存储状态（不在这里生成，只是记录配置状态）
            if sashimi_config:
                results['sashimi_processing'] = {
                    'success': True,
                    'message': 'Sashimi配置已保存，可在分析页面中使用',
                    'config_stored': True,
                    'config_summary': {
                        'has_bam_files': bool(sashimi_config.get('bam_files')),
                        'has_gff_file': bool(sashimi_config.get('gff_file')),
                        'has_display_options': bool(sashimi_config.get('display_options'))
                    }
                }
                results['processing_summary']['modules_executed'].append('SashimiConfig_storage')
            else:
                results['sashimi_processing'] = {
                    'success': False,
                    'message': '未提供SashimiPlot配置，跳过配置存储',
                    'skipped': True
                }
                results['processing_summary']['modules_skipped'].append('SashimiConfig_storage')
            
            # 处理完成，生成最终结果
            end_time = time.time()
            results['processing_summary']['total_processing_time'] = round(end_time - start_time, 2)
            
            # 判断整体成功状态
            rmats_success = results['rmats_processing'] and results['rmats_processing'].get('success', False)
            bam_success = results['bam_analysis'] and results['bam_analysis'].get('success', False)
            
            overall_success = rmats_success and bam_success  # 必选项都成功才算成功
            
            print(f"\n" + "="*30 + " 处理完成 " + "="*30)
            print(f"总处理时间: {results['processing_summary']['total_processing_time']}秒")
            print(f"已执行模块: {', '.join(results['processing_summary']['modules_executed'])}")
            print(f"已跳过模块: {', '.join(results['processing_summary']['modules_skipped'])}")
            print(f"整体状态: {'成功' if overall_success else '部分失败'}")
            
            return jsonify({
                'success': overall_success,
                'message': f"处理完成 - 执行了 {len(results['processing_summary']['modules_executed'])} 个模块",
                'results': results,
                'sashimi_enabled': sashimi_config is not None,
                'output_directory': output_dir
            }), 200 if overall_success else 207  # 207 Multi-Status for partial success
            
        except Exception as e:
            import traceback
            print(f"处理请求时出错: {e}")
            print(traceback.format_exc())
            return jsonify({
                'success': False,
                'message': f'处理请求时出错: {str(e)}',
                'error': str(e)
            }), 500
    # ==================== SashimiPlot 相关路由 ====================
    
    # 初始化SashimiPlot服务（在路由外部，避免重复初始化）
    try:
        from services.sashimi_service import SashimiPlotService
        sashimi_service = SashimiPlotService()
        print("✅ SashimiPlot服务初始化成功")
    except ImportError as e:
        print(f"⚠️ SashimiPlot服务导入失败: {e}")
        sashimi_service = None
    except Exception as e:
        print(f"⚠️ SashimiPlot服务初始化失败: {e}")
        sashimi_service = None
    
    @app.route('/api/sashimi/generate', methods=['POST', 'OPTIONS'])
    def generate_sashimi_plot():
        """生成 SashimiPlot（事件文件模式）"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'success': False,
                'message': 'SashimiPlot服务未正确加载',
                'error': 'Service not initialized'
            }), 500
        
        try:
            data = request.json
            print("\n" + "🍣" * 20)
            print("📌 收到SashimiPlot生成请求（rMATS事件文件模式）")
            print("🍣" * 20)
            print(json.dumps(data, indent=2, ensure_ascii=False))
            print("-" * 50)
            
            # 使用rMATS模式处理
            result = sashimi_service.generate_sashimi_by_rmats(data)
            
            if result['success']:
                print("✅ rMATS模式SashimiPlot生成成功")
                return jsonify(result)
            else:
                print("❌ rMATS模式SashimiPlot生成失败")
                return jsonify(result), 400
            
        except Exception as e:
            import traceback
            print(f"❌ rMATS模式SashimiPlot处理过程中发生错误: {e}")
            print(traceback.format_exc())
            return jsonify({
                'success': False,
                'message': f'rMATS模式处理过程中发生错误: {str(e)}',
                'error': str(e)
            }), 500
    
    @app.route('/api/sashimi/generate-by-coordinate', methods=['POST', 'OPTIONS'])
    def generate_coordinate_sashimi_plot():
        """生成SashimiPlot（坐标模式）- 单个位点可视化"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'success': False,
                'message': 'SashimiPlot服务未正确加载',
                'error': 'Service not initialized'
            }), 500
        
        try:
            data = request.json
            print("\n" + "🎯" * 20)
            print("📌 收到SashimiPlot坐标模式生成请求")
            print("🎯" * 20)
            print(json.dumps(data, indent=2, ensure_ascii=False))
            print("-" * 50)
            
            # 获取预存的 Sashimi 配置
            stored_config = app.config.get('SASHIMI_CONFIG', {})
            
            if stored_config:
                print(f"🔧 使用预存的Sashimi配置: {stored_config}")
                
                # 合并 BAM 文件路径
                bam_files = stored_config.get('bam_files', {})
                if bam_files.get('sample1') and not data.get('sample1_bam'):
                    data['sample1_bam'] = bam_files['sample1']
                if bam_files.get('sample2') and not data.get('sample2_bam'):
                    data['sample2_bam'] = bam_files['sample2']
                    
                # 合并显示选项
                display_options = stored_config.get('display_options', {})
                for key, value in display_options.items():
                    if key not in data and value is not None:
                        data[key] = value
            
            # 使用坐标模式处理
            result = sashimi_service.generate_sashimi_by_coordinate(data)
            
            if result['success']:
                print("✅ 坐标模式SashimiPlot生成成功")
                return jsonify(result)
            else:
                print("❌ 坐标模式SashimiPlot生成失败")
                return jsonify(result), 400
                
        except Exception as e:
            import traceback
            print(f"❌ 坐标模式生成过程中发生错误: {e}")
            print(traceback.format_exc())
            return jsonify({
                'success': False,
                'message': '坐标模式生成过程中出错',
                'error': str(e)
            }), 500
    # 新增：前端期望的rmats模式API端点
    @app.route('/api/sashimi/generate-rmats', methods=['POST', 'OPTIONS'])
    def generate_rmats_sashimi_plot():
        """生成SashimiPlot（前端调用的rmats模式）- 与坐标模式相同"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'success': False,
                'message': 'SashimiPlot服务未正确加载',
                'error': 'Service not initialized'
            }), 500
        
        try:
            data = request.json
            print("\n" + "🔧" * 20)
            print("📌 收到前端rmats2sashimiplot模式生成请求")
            print("🔧" * 20)
            print(json.dumps(data, indent=2, ensure_ascii=False))
            print("-" * 50)
            
            # 实际使用坐标模式处理（因为前端发送的是坐标格式）
            result = sashimi_service.generate_sashimi_by_coordinate(data)
            
            if result['success']:
                print("✅ 前端rmats2sashimiplot模式生成成功")
                return jsonify(result)
            else:
                print("❌ 前端rmats2sashimiplot模式生成失败")
                return jsonify(result), 400
                
        except Exception as e:
            import traceback
            print(f"❌ 前端rmats2sashimiplot模式生成过程中发生错误: {e}")
            print(traceback.format_exc())
            return jsonify({
                'success': False,
                'message': '前端rmats2sashimiplot模式生成过程中出错',
                'error': str(e)
            }), 500
    
    
    # 新增：前端期望的rmats工具检查端点
    @app.route('/api/sashimi/check-rmats', methods=['GET', 'OPTIONS'])
    def check_rmats_tool_availability():
        """检查 rmats2sashimiplot 工具是否可用（前端调用版本）"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'available': False,
                'message': 'SashimiPlot服务未正确加载',
                'error': 'Service not initialized'
            }), 500
        
        try:
            print("🔍 检查rmats2sashimiplot工具可用性（前端调用）")
            
            result = sashimi_service.check_rmats2sashimi_availability()
            
            if result.get('available'):
                print("✅ rmats2sashimiplot工具可用（前端调用）")
            else:
                print("❌ rmats2sashimiplot工具不可用（前端调用）")
            
            return jsonify(result)
            
        except Exception as e:
            print(f"❌ 前端工具检查过程中发生错误: {e}")
            return jsonify({
                'available': False,
                'message': f'前端工具检查过程中发生错误: {str(e)}',
                'error': str(e)
            }), 500
    
    @app.route('/api/sashimi/metadata', methods=['POST', 'OPTIONS'])
    def get_sashimi_metadata():
        """获取指定目录的 SashimiPlot 元数据"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'success': False,
                'message': 'SashimiPlot服务未正确加载'
            }), 500
        
        try:
            data = request.json
            output_dir = data.get('output_dir')
            
            if not output_dir:
                return jsonify({
                    'success': False,
                    'message': '缺少 output_dir 参数'
                }), 400
            
            print(f"📊 获取SashimiPlot元数据: {output_dir}")
            
            metadata = sashimi_service.get_sashimi_metadata(output_dir)
            
            return jsonify({
                'success': True,
                'metadata': metadata
            })
            
        except Exception as e:
            print(f"❌ 获取元数据失败: {e}")
            return jsonify({
                'success': False,
                'message': f'获取元数据失败: {str(e)}',
                'error': str(e)
            }), 500
    
    @app.route('/api/sashimi/clean', methods=['POST', 'OPTIONS'])
    def clean_sashimi_output():
        """清理输出目录中的旧SashimiPlot文件"""
        
        if request.method == 'OPTIONS':
            return jsonify({'status': 'ok'}), 200
        
        if not sashimi_service:
            return jsonify({
                'success': False,
                'message': 'SashimiPlot服务未正确加载'
            }), 500
        
        try:
            data = request.json
            output_dir = data.get('output_dir')
            keep_recent = data.get('keep_recent', 5)
            
            if not output_dir:
                return jsonify({
                    'success': False,
                    'message': '缺少 output_dir 参数'
                }), 400
            
            print(f"🧹 清理SashimiPlot输出目录: {output_dir} (保留最近{keep_recent}个文件)")
            
            result = sashimi_service.clean_output_directory(output_dir, keep_recent)
            
            if result.get('success'):
                print(f"✅ 清理完成: {result.get('message')}")
            else:
                print(f"❌ 清理失败: {result.get('message')}")
            
            return jsonify(result)
            
        except Exception as e:
            print(f"❌ 清理失败: {e}")
            return jsonify({
                'success': False,
                'message': f'清理失败: {str(e)}',
                'error': str(e)
            }), 500
    
    # ==================== 辅助函数 ====================
    
    def _process_single_analysis(data):
        """处理单分析"""
        print("🔍 检测到: 单分析模式")
        
        gtf_file = data.get('gtf_file', '').strip()
        rmats_dir = data.get('rmats_dir', '').strip()
        
        print(f"GTF文件路径: {gtf_file if gtf_file else '(未提供，将跳过GTF相关处理)'}")
        print(f"rMATS目录: {rmats_dir}")
        
        # 验证必需参数 - 只有rmats_dir是必需的
        if not rmats_dir:
            return {
                'success': False,
                'message': 'rMATS目录路径不能为空', 
                'error': 'rmats_dir parameter is required and cannot be empty'
            }
        
        # 验证rMATS目录是否存在
        if not os.path.exists(rmats_dir):
            return {
                'success': False,
                'message': f'rMATS目录不存在: {rmats_dir}',
                'error': f'rMATS directory not found: {rmats_dir}'
            }
        
        if not os.path.isdir(rmats_dir):
            return {
                'success': False,
                'message': f'rMATS路径不是目录: {rmats_dir}',
                'error': f'rMATS path is not a directory: {rmats_dir}'
            }
        
        # 验证GTF文件（仅当提供时）
        if gtf_file:
            if not os.path.exists(gtf_file):
                return {
                    'success': False,
                    'message': f'GTF文件不存在: {gtf_file}',
                    'error': f'GTF file not found: {gtf_file}'
                }
            
            if not os.path.isfile(gtf_file):
                return {
                    'success': False,
                    'message': f'GTF路径不是文件: {gtf_file}',
                    'error': f'GTF path is not a file: {gtf_file}'
                }
            print("✅ GTF文件验证通过")
        else:
            gtf_file = None  # 明确设置为None以便后续处理
            print("⚠️ 未提供GTF文件，将跳过基因注释相关功能")
        
        print("📂 输出目录: 自动检测到前端public目录")
        
        # 导入rMATS处理服务
        try:
            from services.rmats_processor import RMATSProcessor
            result = RMATSProcessor.process_single_analysis(
                gtf_file=gtf_file,
                rmats_dir=rmats_dir,
                output_dir=Config.get_frontend_public_path()
            )
            
            print(f"✅ 单分析处理完成")
            return result
            
        except Exception as e:
            import traceback
            print(f"❌ 单分析处理出错: {e}")
            print(traceback.format_exc())
            return {
                'success': False,
                'message': f'rMATS处理过程中发生错误: {str(e)}',
                'error': str(e),
                'traceback': traceback.format_exc()
            }
    
    # ==================== 调试和测试路由 ====================
    
    @app.route('/api/test')
    def test_api():
        """测试API连通性"""
        return jsonify({
            'message': 'API is working!',
            'config': {
                'backend_dir': str(Config.BACKEND_DIR),
                'frontend_public_dir': Config.get_frontend_public_path(),
                'required_rmats_files': Config.REQUIRED_RMATS_FILES
            },
            'available_endpoints': [
                'GET /api/test - API测试',
                'POST /api/validate-path - 路径验证',
                'POST /api/validate-bam-lists - 双BAM列表验证',
                'POST /api/analyze-bam-files - BAM文件分析',
                'POST /api/create-directory - 创建目录',
                'POST /api/process-rmats - rMATS数据处理',
                'POST /api/sashimi/generate - 生成SashimiPlot（rMATS事件文件模式）',
                'POST /api/sashimi/generate-by-coordinate - 生成SashimiPlot（坐标模式）',
                'POST /api/sashimi/generate-rmats - 生成SashimiPlot（前端rmats模式）',
                'GET /api/sashimi/check - 检查SashimiPlot工具',
                'GET /api/sashimi/check-rmats - 检查SashimiPlot工具（前端版本）',
                'POST /api/sashimi/metadata - 获取SashimiPlot元数据',
                'POST /api/sashimi/clean - 清理SashimiPlot文件'
            ]
        })
    
    @app.route('/api/sashimi/test')
    def test_sashimi_api():
        """测试SashimiPlot API连通性"""
        try:
            if not sashimi_service:
                return jsonify({
                    'message': 'SashimiPlot服务未加载',
                    'service_status': 'not_loaded',
                    'error': 'Service initialization failed'
                }), 500
            
            # 检查工具可用性
            availability = sashimi_service.check_rmats2sashimi_availability()
            
            return jsonify({
                'message': 'SashimiPlot API is working!',
                'service_status': 'loaded',
                'tool_availability': availability,
                'endpoints': [
                    'POST /api/sashimi/generate - 生成SashimiPlot（rMATS事件文件模式）',
                    'POST /api/sashimi/generate-by-coordinate - 生成SashimiPlot（坐标模式）',
                    'POST /api/sashimi/generate-rmats - 生成SashimiPlot（前端rmats模式）',
                    'GET /api/sashimi/check - 检查工具可用性',
                    'GET /api/sashimi/check-rmats - 检查工具可用性（前端版本）',
                    'POST /api/sashimi/metadata - 获取元数据',
                    'POST /api/sashimi/clean - 清理文件',
                    'GET /api/sashimi/test - 测试API'
                ]
            })
        except Exception as e:
            return jsonify({
                'message': 'SashimiPlot API has issues!',
                'error': str(e),
                'service_status': 'error'
            }), 500
    
    @app.route('/api/analyze-rmats', methods=['POST'])
    def analyze_rmats_results():
        """分析rMATS结果接口"""
        try:
            data = request.get_json()
            input_path = data.get('input_path')
            output_path = data.get('output_path')  # 可选
            fdr_threshold = data.get('fdr_threshold', 0.05)
            
            if not input_path:
                return jsonify({"error": "input_path is required"}), 400
            
            # 执行分析并生成JSON文件
            result_file = rmats_analyzer.analyze_and_save(input_path, output_path, fdr_threshold)
            
            return jsonify({
                "status": "success", 
                "message": "rMATS analysis completed",
                "result_file": result_file
            })
            
        except Exception as e:
            return jsonify({"error": str(e)}), 500
    
    print("📝 所有API路由已注册 (包括BAM分析功能)")