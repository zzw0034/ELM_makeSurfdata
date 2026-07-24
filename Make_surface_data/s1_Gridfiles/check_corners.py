#!/usr/bin/env python3
# %%
"""
检查 SCRIP 网格文件中格子角点 (corners) 的绕行方向。

背景 / 原理
-----------
SCRIP / ESMF 要求每个网格单元的四个角点按「逆时针」(counter-clockwise, CCW)
顺序排列（从球外看）。若顺序是顺时针 (CW)，重网格化权重可能出错或符号相反。

本脚本读取 SCRIP 文件中的 grid_corner_lat / grid_corner_lon，取出指定格子
的四个角点，用「鞋带公式」(Shoelace formula) 计算平面上的有符号面积：

    area = Σ (lon_n * lat_{n+1} - lon_{n+1} * lat_n)   (n 循环取模)

  - area > 0  → 逆时针 (CCW)，符合 SCRIP 标准
  - area < 0  → 顺时针 (CW)，需要翻转角点顺序

用法
----

cd Make_surface_data/s1_Gridfiles
python check_corners.py scr_out/SEUS_1_24deg_grids.nc


  python check_corners.py <scrip_file.nc>
  python check_corners.py scr_out/SEUS_1_24deg_grids.nc --cell 0

批量检查所有格子请用：
  s2_Mappingfiles/scripts/check_scrip_vertex_order.py <scrip_file.nc>
"""

import argparse
from netCDF4 import Dataset


def check_corners(file_path, cell_idx=0):
    """检查单个格子的角点绕向，打印结果。"""
    with Dataset(file_path) as f:
        print(f)
        latv = f.variables["grid_corner_lat"][:]
        lonv = f.variables["grid_corner_lon"][:]

    print(latv.shape)
    print(lonv.shape)

    if cell_idx < 0 or cell_idx >= latv.shape[0]:
        raise ValueError(f"cell index {cell_idx} out of range (0..{latv.shape[0] - 1})")

    lat = latv[cell_idx, :]
    lon = lonv[cell_idx, :]

    # 鞋带公式：正值 = 逆时针；负值 = 顺时针
    area = 0.0
    for n in range(4):
        nxt = (n + 1) % 4
        area += lon[n] * lat[nxt] - lon[nxt] * lat[n]

    print(f"\nFile: {file_path}")
    print(f"Cell index: {cell_idx}")
    print("Corners (lat, lon):")
    for n in range(4):
        print(f" corner {n}: {lat[n]:.4f}, {lon[n]:.4f}")

    print(f"Approx signed area: {area}")
    if area > 0:
        print("→ Counter-clockwise (SCRIP standard)")
    else:
        print("→ Clockwise (needs flipping)")

    return area


def main():
    parser = argparse.ArgumentParser(
        description="Check SCRIP grid corner winding order for one cell."
    )
    parser.add_argument("file", help="Path to SCRIP NetCDF file")
    parser.add_argument(
        "--cell",
        type=int,
        default=0,
        help="Cell index to check (default: 0)",
    )
    args = parser.parse_args()
    check_corners(args.file, args.cell)


if __name__ == "__main__":
    main()

# %%
# 交互式使用时，取消注释并修改路径：
# check_corners("scr_out/SEUS_1_24deg_grids.nc", cell_idx=0)
