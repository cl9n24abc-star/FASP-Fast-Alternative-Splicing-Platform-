#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 直接测试不带chr前缀的格式
def quick_test_direct():
    print("直接测试模式")
    
    bam_file = "/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam"
    
    # 测试几个染色体区域
    regions = ["1:10000000-15000000", "2:10000000-15000000", "3:10000000-15000000"]
    
    import subprocess
    
    for region in regions:
        print(f"测试区域: {region}")
        try:
            cmd = ['samtools', 'coverage', '-r', region, bam_file]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                print(result.stdout)
            else:
                print(f"错误: {result.stderr}")
                
        except Exception as e:
            print(f"异常: {e}")

if __name__ == "__main__":
    quick_test_direct()
