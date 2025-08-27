import os
from pathlib import Path

class Config:
    """后端配置类"""
    
    # 基础路径配置
    BACKEND_DIR = Path(__file__).parent                    # backend目录
    PROJECT_ROOT = BACKEND_DIR.parent                      # 项目根目录  
    FRONTEND_DIR = PROJECT_ROOT / "frontend"               # 前端目录
    FRONTEND_PUBLIC_DIR = FRONTEND_DIR / "public"          # 前端public目录
    
    # Flask服务器配置
    HOST = '0.0.0.0'
    PORT = 5000
    DEBUG = True
    
    # rMATS必需文件配置
    REQUIRED_RMATS_FILES = [
        'SE.MATS.JC.txt',
        'MXE.MATS.JC.txt', 
        'A5SS.MATS.JC.txt',
        'A3SS.MATS.JC.txt',
        'RI.MATS.JC.txt'
    ]
    
    # GTF文件配置
    VALID_GTF_EXTENSIONS = [
        '.gtf', '.gtf.gz', 
        '.gff', '.gff.gz', 
        '.gff3', '.gff3.gz'
    ]
    
    # 文件大小限制 (16MB)
    MAX_FILE_SIZE = 16 * 1024 * 1024
    
    @classmethod
    def get_frontend_public_path(cls):
        """获取前端public目录的绝对路径字符串"""
        return str(cls.FRONTEND_PUBLIC_DIR.absolute())
    
    @classmethod
    def ensure_frontend_public_exists(cls):
        """确保前端public目录存在"""
        cls.FRONTEND_PUBLIC_DIR.mkdir(parents=True, exist_ok=True)
        return cls.get_frontend_public_path()
    
    @classmethod
    def print_config_info(cls):
        """打印配置信息（用于调试）"""
        print(f"📁 后端目录: {cls.BACKEND_DIR}")
        print(f"📁 项目根目录: {cls.PROJECT_ROOT}")
        print(f"📁 前端目录: {cls.FRONTEND_DIR}")
        print(f"📁 前端public目录: {cls.get_frontend_public_path()}")
        print(f"🌐 服务器: http://{cls.HOST}:{cls.PORT}")
        print(f"🔧 调试模式: {cls.DEBUG}")


# 创建一个全局配置实例（可选）
config = Config()

# 如果直接运行这个文件，打印配置信息
if __name__ == '__main__':
    Config.print_config_info()