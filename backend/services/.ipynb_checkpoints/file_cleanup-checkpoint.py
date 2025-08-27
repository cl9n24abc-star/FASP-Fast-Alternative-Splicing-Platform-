import os
import shutil
import json
import logging
from pathlib import Path
from typing import List, Dict, Any

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class FileCleanupService:
    """文件清理服务类"""
    
    def __init__(self, base_public_path: str = "../../frontend/public"):
        """
        初始化文件清理服务
        
        Args:
            base_public_path: public文件夹的相对路径
        """
        self.base_public_path = Path(base_public_path).resolve()
        
    def cleanup_json_files(self) -> Dict[str, Any]:
        """
        清理public文件夹中的JSON文件
        
        Returns:
            Dict: 包含清理结果的字典
        """
        result = {
            "success": True,
            "deleted_files": [],
            "errors": [],
            "message": ""
        }
        
        try:
            if not self.base_public_path.exists():
                result["success"] = False
                result["message"] = f"Public folder not found: {self.base_public_path}"
                logger.warning(result["message"])
                return result
            
            # 查找并删除JSON文件
            json_files = list(self.base_public_path.glob("*.json"))
            
            for json_file in json_files:
                try:
                    json_file.unlink()
                    result["deleted_files"].append(str(json_file))
                    logger.info(f"Deleted JSON file: {json_file}")
                except Exception as e:
                    error_msg = f"Failed to delete {json_file}: {str(e)}"
                    result["errors"].append(error_msg)
                    logger.error(error_msg)
            
            if result["deleted_files"]:
                result["message"] = f"Successfully deleted {len(result['deleted_files'])} JSON file(s)"
            else:
                result["message"] = "No JSON files found to delete"
                
            # 如果有错误但也有成功删除的文件，标记为部分成功
            if result["errors"] and result["deleted_files"]:
                result["success"] = True
                result["message"] += f", but encountered {len(result['errors'])} error(s)"
            elif result["errors"] and not result["deleted_files"]:
                result["success"] = False
                result["message"] = "Failed to delete any files"
                
        except Exception as e:
            result["success"] = False
            result["message"] = f"Cleanup operation failed: {str(e)}"
            result["errors"].append(str(e))
            logger.error(f"Cleanup operation failed: {str(e)}")
        
        return result
    
    def cleanup_generated_folders(self, folders_to_clean: List[str] = None) -> Dict[str, Any]:
        """
        清理指定的生成文件夹
        
        Args:
            folders_to_clean: 要清理的文件夹列表，默认清理常见的生成文件夹
            
        Returns:
            Dict: 包含清理结果的字典
        """
        if folders_to_clean is None:
            folders_to_clean = [
                "sashimi_output",
                "generated_plots",
                "analysis_results",
                "temp_files"
            ]
        
        result = {
            "success": True,
            "deleted_folders": [],
            "deleted_files": [],
            "errors": [],
            "message": ""
        }
        
        try:
            if not self.base_public_path.exists():
                result["success"] = False
                result["message"] = f"Public folder not found: {self.base_public_path}"
                logger.warning(result["message"])
                return result
            
            for folder_name in folders_to_clean:
                folder_path = self.base_public_path / folder_name
                
                if folder_path.exists() and folder_path.is_dir():
                    try:
                        # 记录删除的文件
                        for file_path in folder_path.rglob("*"):
                            if file_path.is_file():
                                result["deleted_files"].append(str(file_path))
                        
                        # 删除整个文件夹
                        shutil.rmtree(folder_path)
                        result["deleted_folders"].append(str(folder_path))
                        logger.info(f"Deleted folder: {folder_path}")
                        
                    except Exception as e:
                        error_msg = f"Failed to delete folder {folder_path}: {str(e)}"
                        result["errors"].append(error_msg)
                        logger.error(error_msg)
                else:
                    logger.info(f"Folder {folder_path} does not exist, skipping")
            
            # 生成结果消息
            deleted_folders_count = len(result["deleted_folders"])
            deleted_files_count = len(result["deleted_files"])
            
            if deleted_folders_count > 0:
                result["message"] = f"Successfully deleted {deleted_folders_count} folder(s) and {deleted_files_count} file(s)"
            else:
                result["message"] = "No folders found to delete"
            
            # 处理错误情况
            if result["errors"]:
                if result["deleted_folders"]:
                    result["message"] += f", but encountered {len(result['errors'])} error(s)"
                else:
                    result["success"] = False
                    result["message"] = "Failed to delete any folders"
                    
        except Exception as e:
            result["success"] = False
            result["message"] = f"Folder cleanup operation failed: {str(e)}"
            result["errors"].append(str(e))
            logger.error(f"Folder cleanup operation failed: {str(e)}")
        
        return result
    
    def cleanup_all(self, folders_to_clean: List[str] = None) -> Dict[str, Any]:
        """
        执行完整的清理操作，包括JSON文件和生成的文件夹
        
        Args:
            folders_to_clean: 要清理的文件夹列表
            
        Returns:
            Dict: 包含清理结果的字典
        """
        logger.info("Starting complete cleanup operation")
        
        # 清理JSON文件
        json_result = self.cleanup_json_files()
        
        # 清理生成的文件夹
        folder_result = self.cleanup_generated_folders(folders_to_clean)
        
        # 合并结果
        combined_result = {
            "success": json_result["success"] and folder_result["success"],
            "json_cleanup": json_result,
            "folder_cleanup": folder_result,
            "total_deleted_files": len(json_result["deleted_files"]) + len(folder_result["deleted_files"]),
            "total_deleted_folders": len(folder_result["deleted_folders"]),
            "message": ""
        }
        
        if combined_result["success"]:
            combined_result["message"] = f"Complete cleanup successful: {combined_result['total_deleted_files']} files and {combined_result['total_deleted_folders']} folders removed"
        else:
            combined_result["message"] = "Cleanup completed with some errors"
        
        logger.info(f"Cleanup operation completed: {combined_result['message']}")
        return combined_result

# 便捷函数
def cleanup_public_files(public_path: str = "../../frontend/public", 
                        folders_to_clean: List[str] = None) -> Dict[str, Any]:
    """
    便捷函数：清理public文件夹中的文件
    
    Args:
        public_path: public文件夹路径
        folders_to_clean: 要清理的文件夹列表
        
    Returns:
        Dict: 清理结果
    """
    service = FileCleanupService(public_path)
    return service.cleanup_all(folders_to_clean)

def cleanup_json_only(public_path: str = "../../frontend/public") -> Dict[str, Any]:
    """
    便捷函数：只清理JSON文件
    
    Args:
        public_path: public文件夹路径
        
    Returns:
        Dict: 清理结果
    """
    service = FileCleanupService(public_path)
    return service.cleanup_json_files()

# 测试函数
if __name__ == "__main__":
    # 测试清理服务
    service = FileCleanupService()
    
    print("Testing JSON file cleanup...")
    result = service.cleanup_json_files()
    print(json.dumps(result, indent=2))
    
    print("\nTesting complete cleanup...")
    result = service.cleanup_all()
    print(json.dumps(result, indent=2))