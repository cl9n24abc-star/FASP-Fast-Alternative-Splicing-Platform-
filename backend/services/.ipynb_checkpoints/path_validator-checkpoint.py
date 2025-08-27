import os
import gzip
from typing import Dict, List, Any, Tuple
from config import Config

class PathValidator:
    """路径验证服务类"""
    
    @classmethod
    def validate_rmats_directory(cls, rmats_path: str) -> Dict[str, Any]:
        """验证rMATS目录并检查必需的文件"""
        try:
            # 检查目录是否存在
            if not os.path.exists(rmats_path):
                return cls._create_error_result(
                    f'Directory does not exist: {rmats_path}',
                    exists=False
                )
            
            if not os.path.isdir(rmats_path):
                return cls._create_error_result(
                    f'Path is not a directory: {rmats_path}',
                    exists=True
                )
            
            # 检查rMATS必需文件
            found_files = []
            missing_files = []
            file_details = []
            
            for filename in Config.REQUIRED_RMATS_FILES:
                file_path = os.path.join(rmats_path, filename)
                if os.path.exists(file_path) and os.path.isfile(file_path):
                    found_files.append(filename)
                    file_size = os.path.getsize(file_path)
                    file_details.append({
                        'name': filename,
                        'size': file_size,
                        'size_mb': round(file_size / (1024 * 1024), 2),
                        'exists': True
                    })
                else:
                    missing_files.append(filename)
                    file_details.append({
                        'name': filename,
                        'size': 0,
                        'size_mb': 0,
                        'exists': False
                    })
            
            is_valid = len(missing_files) == 0
            
            return {
                'valid': is_valid,
                'exists': True,
                'files': found_files,
                'missing': missing_files,
                'file_details': file_details,
                'total_files': len(found_files),
                'required_files': len(Config.REQUIRED_RMATS_FILES),
                'message': f'Found {len(found_files)}/{len(Config.REQUIRED_RMATS_FILES)} required rMATS files',
                'error': None if is_valid else f'Missing {len(missing_files)} required files: {", ".join(missing_files)}'
            }
            
        except Exception as e:
            return cls._create_error_result(
                f'Error accessing directory: {str(e)}',
                exists=False
            )
    
    @classmethod
    def validate_output_directory(cls, output_path: str) -> Dict[str, Any]:
        """验证输出目录"""
        try:
            if os.path.exists(output_path):
                if not os.path.isdir(output_path):
                    return {
                        'valid': False,
                        'exists': True,
                        'writable': False,
                        'error': f'Path exists but is not a directory: {output_path}',
                        'can_create': False,
                        'message': 'Path exists but is not a directory'
                    }
                
                # 检查是否可写
                writable = cls._check_directory_writable(output_path)
                
                return {
                    'valid': writable,
                    'exists': True,
                    'writable': writable,
                    'can_create': False,
                    'error': None if writable else 'Directory exists but is not writable',
                    'message': 'Directory exists and is writable' if writable else 'Directory exists but is not writable'
                }
            else:
                # 目录不存在，检查是否可以创建
                parent_dir = os.path.dirname(output_path)
                parent_exists = os.path.exists(parent_dir)
                can_create = parent_exists and os.access(parent_dir, os.W_OK)
                
                return {
                    'valid': False,  # 目录不存在所以当前无效
                    'exists': False,
                    'writable': False,
                    'can_create': can_create,
                    'parent_exists': parent_exists,
                    'error': f'Directory does not exist: {output_path}',
                    'message': 'Directory can be created' if can_create else 'Cannot create directory (parent directory issues)'
                }
                
        except Exception as e:
            return {
                'valid': False,
                'exists': False,
                'writable': False,
                'can_create': False,
                'error': f'Error checking output directory: {str(e)}',
                'message': 'Error occurred while checking directory'
            }
    
    @classmethod
    def validate_gtf_file(cls, gtf_path: str) -> Dict[str, Any]:
        """验证GTF文件"""
        try:
            if not os.path.exists(gtf_path):
                return {
                    'valid': False,
                    'exists': False,
                    'error': f'GTF file does not exist: {gtf_path}',
                    'size': 0,
                    'size_mb': 0,
                    'is_compressed': False,
                    'message': 'File not found'
                }
            
            if not os.path.isfile(gtf_path):
                return {
                    'valid': False,
                    'exists': True,
                    'error': f'Path is not a file: {gtf_path}',
                    'size': 0,
                    'size_mb': 0,
                    'is_compressed': False,
                    'message': 'Path is not a file'
                }
            
            # 检查文件扩展名
            is_valid_extension = any(gtf_path.lower().endswith(ext) for ext in Config.VALID_GTF_EXTENSIONS)
            
            # 获取文件信息
            file_size = os.path.getsize(gtf_path)
            is_compressed = gtf_path.lower().endswith('.gz')
            
            # 验证文件格式
            format_valid, sample_lines = cls._validate_gtf_format(gtf_path, is_compressed)
            
            overall_valid = is_valid_extension and format_valid and file_size > 0
            
            return {
                'valid': overall_valid,
                'exists': True,
                'size': file_size,
                'size_mb': round(file_size / (1024 * 1024), 2),
                'is_compressed': is_compressed,
                'valid_extension': is_valid_extension,
                'format_valid': format_valid,
                'sample_lines': sample_lines,
                'error': None if overall_valid else cls._get_gtf_error_message(is_valid_extension, format_valid, file_size),
                'message': f'Valid GTF file ({file_size / (1024*1024):.1f} MB)' if overall_valid else 'GTF file validation failed'
            }
            
        except Exception as e:
            return {
                'valid': False,
                'exists': False,
                'error': f'Error validating GTF file: {str(e)}',
                'size': 0,
                'size_mb': 0,
                'is_compressed': False,
                'message': 'Error occurred during validation'
            }
    @classmethod
    def validate_bam_file(cls, bam_path: str) -> Dict[str, Any]:
        """验证单个BAM文件"""
        try:
            if not os.path.exists(bam_path):
                return {
                    'valid': False,
                    'exists': False,
                    'error': f'BAM file does not exist: {bam_path}',
                    'size': 0,
                    'size_mb': 0,
                    'message': 'File not found'
                }
            
            if not os.path.isfile(bam_path):
                return {
                    'valid': False,
                    'exists': True,
                    'error': f'Path is not a file: {bam_path}',
                    'size': 0,
                    'size_mb': 0,
                    'message': 'Path is not a file'
                }
            
            # 检查文件扩展名
            if not bam_path.lower().endswith('.bam'):
                return {
                    'valid': False,
                    'exists': True,
                    'error': f'File must have .bam extension: {bam_path}',
                    'size': os.path.getsize(bam_path),
                    'size_mb': round(os.path.getsize(bam_path) / (1024 * 1024), 2),
                    'message': 'Invalid file extension'
                }
            
            # 获取文件信息
            file_size = os.path.getsize(bam_path)
            
            return {
                'valid': True,
                'exists': True,
                'size': file_size,
                'size_mb': round(file_size / (1024 * 1024), 2),
                'error': None,
                'message': f'Valid BAM file ({file_size / (1024*1024):.1f} MB)'
            }
            
        except Exception as e:
            return {
                'valid': False,
                'exists': False,
                'error': f'Error validating BAM file: {str(e)}',
                'size': 0,
                'size_mb': 0,
                'message': 'Error occurred during validation'
            }
    @classmethod
    def create_directory(cls, output_path: str) -> Dict[str, Any]:
        """创建输出目录"""
        try:
            # 使用makedirs创建目录（包括父目录）
            os.makedirs(output_path, exist_ok=True)
            
            # 验证目录是否成功创建
            if os.path.exists(output_path) and os.path.isdir(output_path):
                return {
                    'success': True,
                    'message': f'Directory created successfully: {output_path}',
                    'path': output_path,
                    'error': None
                }
            else:
                return {
                    'success': False,
                    'message': f'Failed to create directory: {output_path}',
                    'path': output_path,
                    'error': 'Directory creation failed for unknown reason'
                }
                
        except PermissionError:
            return {
                'success': False,
                'message': f'Permission denied: cannot create directory {output_path}',
                'path': output_path,
                'error': 'Permission denied'
            }
        except Exception as e:
            return {
                'success': False,
                'message': f'Failed to create directory: {output_path}',
                'path': output_path,
                'error': str(e)
            }
    
    # ==================== 私有辅助方法 ====================
    
    @classmethod
    def _create_error_result(cls, error_message: str, **kwargs) -> Dict[str, Any]:
        """创建错误结果的辅助方法"""
        result = {
            'valid': False,
            'error': error_message,
            'files': [],
            'missing': [],
            'message': 'Validation failed'
        }
        result.update(kwargs)
        return result
    
    @classmethod
    def _check_directory_writable(cls, directory_path: str) -> bool:
        """检查目录是否可写"""
        try:
            # 尝试在目录中创建一个测试文件
            test_file = os.path.join(directory_path, '.write_test_temp')
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
            return True
        except:
            return False
    
    @classmethod
    def _validate_gtf_format(cls, gtf_path: str, is_compressed: bool) -> Tuple[bool, List[str]]:
        """验证GTF文件格式"""
        try:
            open_func = gzip.open if is_compressed else open
            sample_lines = []
            format_valid = False
            
            with open_func(gtf_path, 'rt', encoding='utf-8', errors='ignore') as f:
                line_count = 0
                data_line_count = 0
                
                for line in f:
                    line = line.strip()
                    line_count += 1
                    
                    # 收集样本行（包括注释）
                    if len(sample_lines) < 3:
                        display_line = line[:100] + '...' if len(line) > 100 else line
                        sample_lines.append(display_line)
                    
                    # 检查数据行格式（非注释行）
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 8:  # GTF标准格式至少8列
                            format_valid = True
                            data_line_count += 1
                            if data_line_count >= 3:  # 检查前3行数据行就够了
                                break
                    
                    # 避免读取过多行
                    if line_count >= 50:
                        break
            
            return format_valid, sample_lines
            
        except Exception as e:
            return False, [f"Error reading file: {str(e)}"]
    
    @classmethod
    def _get_gtf_error_message(cls, is_valid_extension: bool, format_valid: bool, file_size: int) -> str:
        """生成GTF文件错误消息"""
        errors = []
        
        if not is_valid_extension:
            errors.append("invalid file extension")
        if not format_valid:
            errors.append("invalid GTF format")
        if file_size == 0:
            errors.append("empty file")
        
        return f"GTF file validation failed: {', '.join(errors)}"