# backend/services/__init__.py
"""
Services模块初始化文件
用于将services目录变成Python包，方便其他模块导入
"""

# 这个文件可以是空的，或者导入一些常用的服务类

# 可选：导入主要的服务类，方便外部使用
# from .path_validator import PathValidator
# from .rmats_processor import RMATSProcessor

# 如果需要，可以在这里定义一些公共的服务函数或常量
__version__ = "1.0.0"
__author__ = "Your Name"

# 导出的服务列表（当使用 from services import * 时）
__all__ = [
    # 'PathValidator',
    # 'RMATSProcessor',
]