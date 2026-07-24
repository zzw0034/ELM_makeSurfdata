#!/usr/bin/env python3
"""
生成 SCRIP 格式的规则经纬度网格文件，供 ESMF_RegridWeightGen / mksurfdata 使用。

背景
----
SCRIP 是气候模式重网格化常用的网格描述格式。本脚本在给定经纬度边界框 (bbox)
内，按固定分辨率（默认 1/24° ≈ 4.17 km）铺规则 lat/lon 格子，写出 NetCDF3
文件，包含格心坐标、四角点坐标、维度信息和掩膜。

主要输出变量（ESMF 读取 SCRIP 所需）：
  - grid_dims        : [nlon, nlat]
  - grid_center_lat/lon : 格心坐标
  - grid_corner_lat/lon : 每格 4 个角点，逆时针 (CCW) 排列
  - grid_imask       : 1=参与重网格，0=不参与

角点顺序（SCRIP/ESMF 要求 CCW）
  SW(0) → SE(1) → NE(2) → NW(3)，从球外看为逆时针。
  可用 check_corners.py 或 check_scrip_vertex_order.py 验证。
"""

import numpy as np
from netCDF4 import Dataset


def snap(val, dx, mode="floor"):
    """将边界值对齐到分辨率的整数倍，避免浮点累积误差导致格点偏移。"""
    k = val / dx
    return (np.floor(k) if mode == "floor" else np.ceil(k)) * dx


def make_scrip(
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    dlon=1 / 24,
    dlat=1 / 24,
    lon_360=False,
    outfile="SEUS_1_24deg_SCRIP.nc",
):
    # 将 bbox 四边 snap 到 dlon/dlat 的整数倍（min 向下、max 向上）
    lon_min = snap(lon_min, dlon, "floor")
    lon_max = snap(lon_max, dlon, "ceil")
    lat_min = snap(lat_min, dlat, "floor")
    lat_max = snap(lat_max, dlat, "ceil")

    # 1D 格心坐标：从 (边界 + 半格宽) 起，步长为分辨率
    lons = np.arange(lon_min + dlon / 2, lon_max, dlon)
    lats = np.arange(lat_min + dlat / 2, lat_max, dlat)

    nlon = lons.size
    nlat = lats.size
    grid_size = nlon * nlat  # 总格数
    grid_corners = 4
    grid_rank = 2  # 2D 逻辑矩形网格

    # 2D 格心场
    LON, LAT = np.meshgrid(lons, lats)

    # 四角点：CCW，顺序 SW → SE → NE → NW（ESMF/SCRIP 要求）
    lon_bnds = np.stack(
        [
            LON - dlon / 2,  # SW lon
            LON + dlon / 2,  # SE lon
            LON + dlon / 2,  # NE lon
            LON - dlon / 2,  # NW lon
        ],
        axis=-1,
    )

    lat_bnds = np.stack(
        [
            LAT - dlat / 2,  # SW lat
            LAT - dlat / 2,  # SE lat
            LAT + dlat / 2,  # NE lat
            LAT + dlat / 2,  # NW lat
        ],
        axis=-1,
    )

    # 展平为 SCRIP 一维 cell 列表 (grid_size,)
    clat = LAT.ravel()
    clon = LON.ravel()
    clat_b = lat_bnds.reshape(grid_size, grid_corners)
    clon_b = lon_bnds.reshape(grid_size, grid_corners)

    # 可选：经度从 [-180,180] 转为 [0,360)
    if lon_360:
        clon = np.mod(clon, 360.0)
        clon_b = np.mod(clon_b, 360.0)

    grid_dims = np.array([nlon, nlat], dtype="i4")

    # 掩膜：bbox 内全部有效（1）；若需裁掉海洋/边界可在此修改
    imask = np.ones(grid_size, dtype="i4")

    # 球面矩形格元面积（m²），当前未写入文件，仅作参考
    # area = R² × |Δλ| × |sin(φ_n) - sin(φ_s)|，λ/φ 为弧度
    R = 6371000.0
    lon_w = clon_b[:, 0]
    lon_e = clon_b[:, 2]
    lat_s = clat_b[:, 0]
    lat_n = clat_b[:, 2]
    dlam = np.deg2rad(lon_e - lon_w)
    dphi = np.abs(np.sin(np.deg2rad(lat_n)) - np.sin(np.deg2rad(lat_s)))
    area = (R**2) * np.abs(dlam) * dphi

    # 写入 SCRIP NetCDF（NETCDF3_CLASSIC，与 ESMF 兼容性最好）
    nc = Dataset(outfile, "w", format="NETCDF3_CLASSIC")
    try:
        nc.createDimension("grid_size", grid_size)
        nc.createDimension("grid_corners", grid_corners)
        nc.createDimension("grid_rank", grid_rank)

        v_dims = nc.createVariable("grid_dims", "i4", ("grid_rank",))
        v_imask = nc.createVariable("grid_imask", "i4", ("grid_size",))
        v_center_lat = nc.createVariable("grid_center_lat", "f8", ("grid_size",))
        v_center_lon = nc.createVariable("grid_center_lon", "f8", ("grid_size",))
        v_corner_lat = nc.createVariable(
            "grid_corner_lat", "f8", ("grid_size", "grid_corners")
        )
        v_corner_lon = nc.createVariable(
            "grid_corner_lon", "f8", ("grid_size", "grid_corners")
        )
        # v_area = nc.createVariable("grid_area", "f8", ("grid_size",))

        v_dims[:] = grid_dims
        v_imask[:] = imask
        v_center_lat[:] = clat
        v_center_lon[:] = clon
        v_corner_lat[:, :] = clat_b
        v_corner_lon[:, :] = clon_b
        # v_area[:] = area

        v_center_lat.units = "degrees"
        v_center_lon.units = "degrees"
        v_corner_lat.units = "degrees"
        v_corner_lon.units = "degrees"
        # v_area.units = "m2"
        nc.title = "SCRIP grid"
        nc.history = "Created by make_scrip_1_24deg.py"
        nc.conventions = "SCRIP"
    finally:
        nc.close()

    print(f"Wrote {outfile}")
    print(f"nlon x nlat = {nlon} x {nlat}  (grid_size={grid_size})")
    return outfile


if __name__ == "__main__":
    # 示例：美国东南部 (SEUS) 区域
    # 经度：degrees_east（西经为负）；纬度：degrees_north
    import os

    _dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(_dir, "scr_out"), exist_ok=True)
    outfile = os.path.join(_dir, "scr_out", "SCRIPgrid_SEUS_1_24deg.nc")

    make_scrip(
        lon_min=-95.0,
        lon_max=-74.0,
        lat_min=24.0,
        lat_max=37.5,
        dlon=1 / 24,
        dlat=1 / 24,
        lon_360=False,
        outfile=outfile,
    )

# 小范围测试网格示例（取消注释即可）：
# outfile = os.path.join(_dir, "scr_out", "small_SEUS_1_24deg_grids.nc")
# make_scrip(
#     lon_min=-85.0,
#     lon_max=-82.0,
#     lat_min=30.0,
#     lat_max=33,
#     dlon=1 / 24,
#     dlat=1 / 24,
#     lon_360=False,
#     outfile=outfile,
# )
