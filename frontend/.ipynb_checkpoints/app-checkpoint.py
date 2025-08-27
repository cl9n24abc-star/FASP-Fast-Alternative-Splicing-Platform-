from flask import Flask, jsonify, request
import os
import sys
import subprocess
import json
from flask import Flask, jsonify, request
from flask_cors import CORS  # 需要安装: pip install flask-cors
from genestructure import process_rmats_data
 
app = Flask(__name__)

@app.route('/')
def index():
    """API服务状态页面"""
    return jsonify({
        'status': 'running',
        'message': 'Flask API服务正在运行',
        'available_endpoints': [
            'GET / - API状态',
            'GET /api/check-file - 文件检查',
            'POST /api/rmats2sashimiplot - rMATS绘图',
            'POST /api/noisyrinput - Noisyr路径接收',
            'POST /api/process-rmats - rMATS数据处理'  # 新增
        ],
        'working_directory': os.getcwd()
    })

@app.route('/api/check-file')
def check_file():
    """检查文件是否存在的API"""
    file_path = request.args.get('path')
    if file_path and os.path.exists(file_path):
        return jsonify({
            'exists': True,
            'size': os.path.getsize(file_path),
            'modified': os.path.getmtime(file_path)
        })
    else:
        return jsonify({'exists': False})

@app.route('/api/rmats2sashimiplot', methods=['POST', 'OPTIONS'])
def rmats2sashimiplot():
    """接收rMATS2sashimiplot参数并处理"""
    
    # 处理OPTIONS请求（CORS预检）
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'}), 200
    
    try:
        # 获取请求数据
        data = request.json
        
        # 打印全部发送过来的信息
        print("收到的完整请求数据:")
        print(json.dumps(data, indent=2, ensure_ascii=False))
        
        # 调用sashimi plot生成函数
        sashimi_result = run_rmats2sashimiplot_from_request(data)
        
        # 返回结果
        if sashimi_result['success']:
            return jsonify({
                'success': True,
                'message': 'Sashimi plot生成成功',
                'sashimi_result': sashimi_result,
                'received_data': data
            })
        else:
            return jsonify({
                'success': False,
                'message': 'Sashimi plot生成失败',
                'error': sashimi_result.get('error', '未知错误'),
                'sashimi_result': sashimi_result
            }), 500
        
    except Exception as e:
        print(f"错误: {e}")
        return jsonify({
            'success': False,
            'message': f'处理过程中发生错误: {str(e)}'
        }), 500

@app.route('/api/noisyrinput', methods=['POST', 'OPTIONS'])
def noisyrinput():
    """接收noisyr对比分析路径参数"""
    
    # 处理OPTIONS请求（CORS预检）
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'}), 200
    
    try:
        # 获取请求数据
        data = request.json
        
        # 打印全部发送过来的信息
        print("收到的完整请求数据:")
        print(json.dumps(data, indent=2, ensure_ascii=False))
        
        # 获取路径参数
        original_path = data.get('original_rmats_path')
        noisyr_path = data.get('noisyr_rmats_path')
        
        print(f"原始rMATS路径: {original_path}")
        print(f"Noisyr处理路径: {noisyr_path}")
        
        # 返回成功结果
        return jsonify({
            'success': True,
            'message': '路径接收成功',
            'job_id': 'noisyr_001',
            'received_data': data
        })
        
    except Exception as e:
        print(f"错误: {e}")
        return jsonify({
            'success': False,
            'message': f'处理过程中发生错误: {str(e)}'
        }), 500
        
# 在函数外面添加sashimi处理函数
def run_rmats2sashimiplot_from_request(request_data):
    """从请求数据运行rmats2sashimiplot"""
    import subprocess
    
    data = request_data
    
    # 构建基础命令
    cmd = ["rmats2sashimiplot"]
    
    # 必需参数
    cmd.extend(["-o", data["output_directory"]])
    cmd.extend(["--event-type", data["event_type"]])
    cmd.extend(["-e", data["events_file"]])
    
    # BAM文件参数
    if data["file_type"] == "bam":
        cmd.extend(["--b1", data["sample1_bam"]])
        cmd.extend(["--b2", data["sample2_bam"]])
    
    # 标签
    cmd.extend(["--l1", data["sample1_label"]])
    cmd.extend(["--l2", data["sample2_label"]])
    
    # 可选参数（只有非None值才添加）
    optional_params = {
        "--exon_s": data.get("exon_scale"),
        "--intron_s": data.get("intron_scale"),
        "--min-counts": data.get("min_counts"),
        "--color": data.get("color"),
        "--font-size": data.get("font_size"),
        "--fig-height": data.get("figure_height"),
        "--fig-width": data.get("figure_width"),
        "--group-info": data.get("group_info")
    }
    
    for param, value in optional_params.items():
        if value and value != "None":
            cmd.extend([param, str(value)])
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, 
                              capture_output=True, 
                              text=True, 
                              check=True)
        
        return {
            "success": True,
            "message": "Sashimi plot generated successfully",
            "stdout": result.stdout,
            "command": ' '.join(cmd)
        }
        
    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "message": "Failed to generate sashimi plot",
            "error": e.stderr,
            "command": ' '.join(cmd)
        }
    except FileNotFoundError:
        return {
            "success": False,
            "message": "rmats2sashimiplot command not found",
            "error": "Please ensure rmats2sashimiplot is installed and accessible"
        }
@app.route('/api/process-rmats', methods=['POST', 'OPTIONS'])
def process_rmats():
    """处理rMATS数据并生成前端格式，自动输出到public目录"""
    
    # 处理OPTIONS请求（CORS预检）
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'}), 200
    
    try:
        # 获取请求数据
        data = request.json
        
        # 打印收到的所有参数
        print("\n" + "="*50)
        print("📌 收到rMATS处理请求")
        print("="*50)
        for key, value in data.items():
            print(f"{key}: {value}")
        print("-"*50)
        
        # 获取分析类型
        analysis_type = data.get('analysis_type')
        
        # 判断分析模式
        if analysis_type == 'single':
            print("🔍 检测到: 单分析模式")
            
            # 获取参数（不需要output_dir，让genestructure自动检测）
            gtf_file = data.get('gtf_file')
            rmats_dir = data.get('rmats_dir')
            
            print(f"GTF文件路径: {gtf_file}")
            print(f"rMATS目录: {rmats_dir}")
            print("📂 输出目录: 自动检测Vue项目public目录")
            
            # 调用处理函数 - 不传output_dir，让它自动检测public目录
            print("🔄 开始调用process_rmats_data函数...")
            result = process_rmats_data(
                gtf_file=gtf_file,
                rmats_dir=rmats_dir,
                output_dir=None,  # 让genestructure.py自动检测public目录
                split_by_type=True
            )
            
            print(f"✅ 处理成功! 处理了 {len(result) if result else 0} 个可变剪切事件")
            
            # 检查生成的文件（genestructure会自动找到public目录）
            try:
                # 尝试找到生成的摘要文件来确认输出位置
                from genestructure import find_vue_public_dir
                detected_public_dir = find_vue_public_dir()
                summary_file = os.path.join(detected_public_dir, 'frontend_data_summary.json')
                
                summary_data = None
                if os.path.exists(summary_file):
                    with open(summary_file, 'r') as f:
                        summary_data = json.load(f)
                    print(f"📋 总事件数: {summary_data.get('total_events', 0)}")
                    print(f"📁 文件已保存到: {detected_public_dir}")
                
                return jsonify({
                    'success': True,
                    'message': '单分析模式处理成功，文件已自动保存到public目录',
                    'analysis_type': 'single',
                    'event_count': len(result) if result else 0,
                    'summary': summary_data,
                    'output_directory': detected_public_dir,  # 告诉前端实际的输出位置
                    'public_files_ready': True,  # 标记public文件已准备好
                    'received_data': data
                }), 200
                
            except Exception as file_check_error:
                print(f"⚠️ 无法检查生成的文件: {file_check_error}")
                return jsonify({
                    'success': True,
                    'message': '单分析模式处理成功',
                    'analysis_type': 'single',
                    'event_count': len(result) if result else 0,
                    'public_files_ready': True,
                    'received_data': data
                }), 200
            
        elif analysis_type == 'compare':
            print("🔍 检测到: 对比分析模式")
            
            # 获取对比分析的参数
            gtf_file = data.get('gtf_file')
            original_rmats_dir = data.get('rmats_dir')  # 原始数据目录
            compare_rmats_dir = data.get('compare_rmats_dir')  # 对比数据目录
            
            print(f"GTF文件路径: {gtf_file}")
            print(f"原始rMATS目录: {original_rmats_dir}")
            print(f"对比rMATS目录: {compare_rmats_dir}")
            print("📂 输出目录: 自动检测Vue项目public目录")
            
            # TODO: 实现对比分析逻辑
            print("🔄 开始处理对比分析...")
            result = process_rmats_data(
                gtf_file=gtf_file,
                rmats_dir=original_rmats_dir,  # 暂时先处理原始数据
                output_dir=None,  # 自动检测
                split_by_type=True
            )
            
            try:
                from genestructure import find_vue_public_dir
                detected_public_dir = find_vue_public_dir()
                
                return jsonify({
                    'success': True,
                    'message': '对比分析模式处理成功，文件已自动保存到public目录',
                    'analysis_type': 'compare',
                    'event_count': len(result) if result else 0,
                    'output_directory': detected_public_dir,
                    'public_files_ready': True,
                    'received_data': data
                }), 200
                
            except Exception as file_check_error:
                return jsonify({
                    'success': True,
                    'message': '对比分析模式处理成功',
                    'analysis_type': 'compare',
                    'event_count': len(result) if result else 0,
                    'public_files_ready': True,
                    'received_data': data
                }), 200
            
        else:
            print(f"❌ 未知的分析类型: {analysis_type}")
            return jsonify({
                'success': False,
                'message': f'未知的分析类型: {analysis_type}',
                'error': '分析类型必须是 single 或 compare'
            }), 400
        
    except Exception as e:
        import traceback
        print(f"❌ 处理请求时出错: {e}")
        print(traceback.format_exc())
        return jsonify({
            'success': False,
            'message': f'处理请求时出错: {str(e)}',
            'error': str(e)
        }), 500
if __name__ == '__main__':
    print("启动Flask API服务器...")
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)