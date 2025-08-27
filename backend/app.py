from flask import Flask, jsonify
from flask_cors import CORS
from config import Config
import sys
import os

def create_app():
    """Flask应用工厂函数"""
    app = Flask(__name__)
    
    # 配置CORS（允许跨域请求）
    CORS(app)
    
    # 添加当前目录到Python路径，确保能导入本地模块
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir not in sys.path:
        sys.path.insert(0, current_dir)
    
    # 注册所有API路由
    from routes import register_routes
    register_routes(app)
    
    # 基础路由
    @app.route('/')
    def index():
        """API首页 - 显示服务状态"""
        return jsonify({
            'status': 'running',
            'message': 'Flask API服务正在运行',
            'version': '1.0.0',
            'available_endpoints': [
                'GET / - 服务状态',
                'GET /health - 健康检查', 
                'POST /api/validate-path - 路径验证',
                'POST /api/create-directory - 创建目录',
                'POST /api/process-rmats - rMATS数据处理'
            ]
        })
    
    @app.route('/health')
    def health_check():
        """健康检查端点"""
        try:
            frontend_exists = Config.FRONTEND_PUBLIC_DIR.exists()
            return jsonify({
                'status': 'healthy',
                'backend_dir': str(Config.BACKEND_DIR),
                'frontend_public_dir': Config.get_frontend_public_path(),
                'frontend_public_exists': frontend_exists,
                'config_loaded': True
            })
        except Exception as e:
            return jsonify({
                'status': 'error',
                'error': str(e)
            }), 500
    
    return app


def main():
    """主函数 - 启动Flask应用"""
    # 打印配置信息
    print("🚀 启动Flask API服务器...")
    Config.print_config_info()
    
    # 确保前端public目录存在
    try:
        Config.ensure_frontend_public_exists()
        print("✅ 前端public目录检查完成")
    except Exception as e:
        print(f"⚠️ 前端public目录检查失败: {e}")
    
    # 创建并启动应用
    app = create_app()
    
    print(f"🌐 API服务地址: http://{Config.HOST}:{Config.PORT}")
    print("📋 可用端点:")
    print("   GET  /        - 服务状态")
    print("   GET  /health  - 健康检查")
    print("   POST /api/*   - API接口")
    
    # 启动服务器
    app.run(
        host=Config.HOST,
        port=Config.PORT,
        debug=Config.DEBUG,
        use_reloader=False  # 避免重复启动
    )

if __name__ == '__main__':
    main()