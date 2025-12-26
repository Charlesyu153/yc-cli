#!/usr/bin/env python3
"""
检查哪些样本缺少 CD66、CD68、CD86 的 marker 列，或这些列的值全为 NA。
"""

from pathlib import Path
import pandas as pd

ROOT = Path("/home/jacekyu/PCF")
CLEANDATA_ROOT = ROOT / "cleandata"

MARKERS_TO_CHECK = ["CD66", "CD68", "CD86"]

def check_sample(sample_dir: Path) -> dict:
    """
    检查单个样本的 marker 列情况。
    
    Returns
    -------
    dict
        {
            "sample_id": str,
            "markers": {
                "CD66": {"exists": bool, "all_na": bool, "has_data": bool},
                "CD68": {...},
                "CD86": {...}
            }
        }
    """
    sample_id = sample_dir.name
    result = {
        "sample_id": sample_id,
        "markers": {}
    }
    
    # 查找 objects.tsv 文件
    tsv_files = list(sample_dir.glob("*_objects.tsv"))
    if not tsv_files:
        result["error"] = "未找到 *_objects.tsv 文件"
        return result
    
    tsv_path = tsv_files[0]
    
    try:
        # 读取文件（只读前几行来检查列名）
        df = pd.read_csv(tsv_path, sep="\t", nrows=1000)
        
        for marker in MARKERS_TO_CHECK:
            marker_info = {
                "exists": False,
                "all_na": False,
                "has_data": False,
                "non_na_count": 0,
                "total_count": 0
            }
            
            if marker in df.columns:
                marker_info["exists"] = True
                marker_info["total_count"] = len(df)
                marker_info["non_na_count"] = df[marker].notna().sum()
                
                if marker_info["non_na_count"] == 0:
                    marker_info["all_na"] = True
                else:
                    marker_info["has_data"] = True
            else:
                marker_info["exists"] = False
            
            result["markers"][marker] = marker_info
            
    except Exception as e:
        result["error"] = str(e)
    
    return result


def main():
    """主函数：检查所有样本。"""
    # 发现所有样本目录
    sample_dirs = sorted([
        d for d in CLEANDATA_ROOT.iterdir() 
        if d.is_dir() and list(d.glob("*_objects.tsv"))
    ])
    
    if not sample_dirs:
        print(f"在 {CLEANDATA_ROOT} 下未找到任何样本目录")
        return
    
    print("=" * 80)
    print(f"检查 {len(sample_dirs)} 个样本的 CD66、CD68、CD86 marker 列情况")
    print("=" * 80)
    print()
    
    results = []
    for sample_dir in sample_dirs:
        result = check_sample(sample_dir)
        results.append(result)
    
    # 汇总统计
    print("详细结果：")
    print("-" * 80)
    
    for result in results:
        sample_id = result["sample_id"]
        if "error" in result:
            print(f"{sample_id:10s} - 错误: {result['error']}")
            continue
        
        print(f"\n样本: {sample_id}")
        for marker in MARKERS_TO_CHECK:
            info = result["markers"][marker]
            if not info["exists"]:
                status = "❌ 列不存在"
            elif info["all_na"]:
                status = "⚠️  列存在但全为 NA"
            else:
                status = f"✅ 有数据 (非NA: {info['non_na_count']}/{info['total_count']})"
            print(f"  {marker:6s}: {status}")
    
    # 汇总表格
    print("\n" + "=" * 80)
    print("汇总表格：")
    print("=" * 80)
    print(f"{'样本ID':<10s} {'CD66':<15s} {'CD68':<15s} {'CD86':<15s}")
    print("-" * 80)
    
    for result in results:
        if "error" in result:
            continue
        
        sample_id = result["sample_id"]
        row = [sample_id]
        
        for marker in MARKERS_TO_CHECK:
            info = result["markers"][marker]
            if not info["exists"]:
                status = "列不存在"
            elif info["all_na"]:
                status = "全为NA"
            else:
                status = f"有数据({info['non_na_count']})"
            row.append(status)
        
        print(f"{row[0]:<10s} {row[1]:<15s} {row[2]:<15s} {row[3]:<15s}")
    
    # 统计缺失情况
    print("\n" + "=" * 80)
    print("缺失统计：")
    print("=" * 80)
    
    for marker in MARKERS_TO_CHECK:
        missing_count = 0
        all_na_count = 0
        has_data_count = 0
        
        for result in results:
            if "error" in result:
                continue
            info = result["markers"][marker]
            if not info["exists"]:
                missing_count += 1
            elif info["all_na"]:
                all_na_count += 1
            else:
                has_data_count += 1
        
        print(f"{marker}:")
        print(f"  列不存在: {missing_count} 个样本")
        print(f"  全为 NA:  {all_na_count} 个样本")
        print(f"  有数据:   {has_data_count} 个样本")
        print()


if __name__ == "__main__":
    main()

