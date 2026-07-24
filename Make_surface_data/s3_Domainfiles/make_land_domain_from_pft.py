"""
从 PFT/植被 NetCDF 与 SCRIP 目标网格生成 ELM 陆面 domain（landfrac + mask）。

做什么
  - 读 SCRIP（grid_center_* / grid_corner_* / grid_dims），按角点算格元面积（球面矩形近似，单位 rad²）。
  - 从 PFT 文件估计「是否有自然植被陆地」：优先 PCT_NATVEG，否则 PCT_NAT_PFT（natpft 求和），
    可选 PFT_FRAC（pft 维求和）；得到 0/1 的 landfrac_pft。
  - 用 scipy.griddata(..., method='nearest') 把 PFT 格点上的 landfrac 映到 SCRIP 格心；SCRIP 超出
    PFT 经纬包络的格点置 landfrac=0、mask=0。
  - 写出 NetCDF3_64BIT：xc, yc, xv, yv, area, frac, mask（dims nj×ni），供 E3SM/ELM 域文件流程使用。

输入 / 输出（见下方 SCRIP / PFT / OUT）
  - SCRIP：与 mksurfdata、surfdata 使用的目标网格一致（如 small_SEUS_1_24deg_grids.nc）。
  - PFT：与 LUT/surfdata 同源网格为佳（常见为 NLCD→ELM 的 1/24° elmpft）；也支持带 LONGXY/LATIXY 的 LUT。


依赖
  - numpy, xarray, scipy, matplotlib（诊断图）；脚本内 os.chdir 指向本脚本目录，换环境请改路径。

运行
  - 在 IDE 中按单元运行,或:python make_land_domain_from_pft.py。

每次运行修改以下变量：
SCRIP = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s1_Gridfiles/scr_out/SCRIPgrid_SEUS_1_24deg.nc"
OUT = "domain.lnd.SEUS_1_24deg.nc"
"""

# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.interpolate import griddata

# %%
os.chdir("/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s3_Domainfiles")

# SCRIP：目标域 SCRIP（与 Step 1 网格一致）
SCRIP = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s1_Gridfiles/scr_out/SCRIPgrid_SEUS_1_24deg.nc"
# PFT：用于陆面掩膜的 ELM PFT/植被场（1D lat/lon 或 LUT 风格坐标）
PFT = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
OUT = "domain.lnd.SEUS_1_24deg.nc"
# %%
# ---- read grid ----
g = xr.open_dataset(SCRIP)


def _get_pft_lon_lat(ds):
    """PFT/surfdata 经纬度：LONGXY/LATIXY, lon2D/lat2D, 或 1D lat/lon (meshgrid)."""
    if "LONGXY" in ds.variables and "LATIXY" in ds.variables:
        return np.asarray(ds["LONGXY"].values).ravel(), np.asarray(
            ds["LATIXY"].values
        ).ravel()
    if "lon2D" in ds.variables and "lat2D" in ds.variables:
        return np.asarray(ds["lon2D"].values).ravel(), np.asarray(
            ds["lat2D"].values
        ).ravel()
    # Regular 1D lat/lon (e.g. elmpft_from_nlcd_frac_pred_*)
    if "lat" in ds.variables and "lon" in ds.variables:
        lat = np.asarray(ds["lat"].values)
        lon = np.asarray(ds["lon"].values)
        if lat.ndim == 1 and lon.ndim == 1:
            lon_2d, lat_2d = np.meshgrid(lon, lat)
            return lon_2d.ravel(), lat_2d.ravel()
    return None, None


def get_lon_lat_extent(ds, name):
    """从 dataset 得到 lon/lat 的 min/max，用于比较范围。"""
    if "LONGXY" in ds.variables and "LATIXY" in ds.variables:
        lon = ds["LONGXY"].values
        lat = ds["LATIXY"].values
    elif "lon2D" in ds.variables and "lat2D" in ds.variables:
        lon = ds["lon2D"].values
        lat = ds["lat2D"].values
    elif "grid_center_lon" in ds.variables and "grid_center_lat" in ds.variables:
        lon = ds["grid_center_lon"].values
        lat = ds["grid_center_lat"].values
    elif "grid_corner_lon" in ds.variables and "grid_corner_lat" in ds.variables:
        lon = ds["grid_corner_lon"].values.ravel()
        lat = ds["grid_corner_lat"].values.ravel()
    elif "lat" in ds.variables and "lon" in ds.variables:
        lat = np.asarray(ds["lat"].values)
        lon = np.asarray(ds["lon"].values)
        if lat.ndim == 1 and lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
    else:
        return None
    lon, lat = np.asarray(lon).ravel(), np.asarray(lat).ravel()
    return {
        "name": name,
        "lon_min": float(np.nanmin(lon)),
        "lon_max": float(np.nanmax(lon)),
        "lat_min": float(np.nanmin(lat)),
        "lat_max": float(np.nanmax(lat)),
    }


# ---- 比较两个 data 的空间范围（SCRIP 与 PFT 可不同） ----
pft = xr.open_dataset(PFT)
ext_pft = get_lon_lat_extent(pft, "PFT")
ext_scrip = get_lon_lat_extent(g, "SCRIP")
if ext_pft and ext_scrip:
    print("范围比较 (lon/lat 单位: degrees):")
    print(
        f"  PFT:   lon [{ext_pft['lon_min']:.4f}, {ext_pft['lon_max']:.4f}], lat [{ext_pft['lat_min']:.4f}, {ext_pft['lat_max']:.4f}]"
    )
    print(
        f"  SCRIP: lon [{ext_scrip['lon_min']:.4f}, {ext_scrip['lon_max']:.4f}], lat [{ext_scrip['lat_min']:.4f}, {ext_scrip['lat_max']:.4f}]"
    )
    scrip_inside = (
        ext_scrip["lon_min"] >= ext_pft["lon_min"]
        and ext_scrip["lon_max"] <= ext_pft["lon_max"]
        and ext_scrip["lat_min"] >= ext_pft["lat_min"]
        and ext_scrip["lat_max"] <= ext_pft["lat_max"]
    )
    print(f"  SCRIP 是否完全在 PFT 范围内: {scrip_inside}")
    if not scrip_inside:
        print("  → 超出 PFT 范围的 SCRIP 格点将设 landfrac=0, mask=0")
pft.close()

# %%
clon = g["grid_center_lon"].values
clat = g["grid_center_lat"].values
lonc = g["grid_corner_lon"].values
latc = g["grid_corner_lat"].values

# ---- compute area (radian^2) ----
lon_w = lonc[:, 0]
lon_e = lonc[:, 2]
lat_s = latc[:, 0]
lat_n = latc[:, 2]

dlam = np.deg2rad(lon_e - lon_w)
dphi = np.sin(np.deg2rad(lat_n)) - np.sin(np.deg2rad(lat_s))
area = np.abs(dlam * dphi)

# ---- PFT: 经纬度与陆地比例（允许与 SCRIP 不同范围/分辨率） ----
p = xr.open_dataset(PFT)
pft_lon, pft_lat = _get_pft_lon_lat(p)
if pft_lon is None:
    raise ValueError(
        "PFT file must have LONGXY/LATIXY, lon2D/lat2D, or 1D lat/lon for coordinates."
    )

# 植被/陆地比例：优先 PCT_NAT_PFT (dim natpft)，否则 PFT_FRAC (dim pft)
if "PCT_NATVEG" in p.variables:
    v = p["PCT_NATVEG"]
    if "time" in v.dims:
        v = v.isel(time=0)
    veg = v.values.ravel()
    # PCT_NATVEG 为百分比 0–100；>0.01% 视为有植被
    landfrac_pft = (veg > 1e-6).astype(np.float64)  # 0/1 landfrac
elif "PCT_NAT_PFT" in p.variables:
    v = p["PCT_NAT_PFT"]
    if "time" in v.dims:
        v = v.isel(time=0)
    veg = v.sum(dim="natpft").values.ravel()
    # PCT_NAT_PFT 为百分比 0–100；>0.01% 视为有植被
    landfrac_pft = (veg > 1e-2).astype(np.float64)  # 0/1 landfrac
elif "PFT_FRAC" in p.variables:
    v = p["PFT_FRAC"]
    if "time" in v.dims:
        v = v.isel(time=0)
    veg = v.sum(dim="pft").values.ravel()
    landfrac_pft = (veg > 1e-6).astype(np.float64)
else:
    raise KeyError(
        "PFT file must have 'PCT_NATVEG', 'PCT_NAT_PFT' (dim natpft), or 'PFT_FRAC' (dim pft). "
        f"Variables found: {list(p.variables)}"
    )
"""
#%%
import matplotlib.pyplot as plt

# Plot landfrac_pft as a map
plt.figure(figsize=(8, 6))
sc = plt.scatter(
    pft_lon,
    pft_lat,
    c=landfrac_pft,
    cmap="viridis",
    s=1,
    vmin=0.0,
    vmax=1.0,
    rasterized=True,
)
plt.title("landfrac_pft (source map)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(sc, label="landfrac_pft")
plt.tight_layout()
plt.show()
# %% test------------------------------------------------------------------------------
p["PCT_NAT_PFT"].isel(time=0).plot(x="lon", y="lat", col="natpft", col_wrap=3)
# %% test nlcd------------------------------------------------------------------------------
v = p["PCT_NAT_PFT"]
if "time" in v.dims:
    v = v.isel(time=0)
    # Exclude natpft==0 from the vegetation sum.

# v_use = v.where(v["natpft"] != 0, drop=True)
# veg = v_use.sum(dim="natpft").values.ravel()
# # PCT_NAT_PFT 为百分比 0–100；>0.01% 视为有植被
# landfrac_pft = (veg > 1e-2).astype(np.float64)
import matplotlib.pyplot as plt

# Plot landfrac_pft as a map
plt.figure(figsize=(8, 6))
sc = plt.scatter(
    pft_lon,
    pft_lat,
    c=landfrac_pft,
    cmap="viridis",
    s=1,
    vmin=0.0,
    vmax=1.0,
    rasterized=True,
)
plt.title("landfrac_pft (source map)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(sc, label="landfrac_pft")
plt.tight_layout()
plt.show()
"""
# %%
# 将 PFT 映射到 SCRIP 格点：griddata 插值或最近邻；SCRIP 超出 PFT 范围则 landfrac=0
pft_points = np.column_stack([pft_lon, pft_lat])
scrip_points = np.column_stack([clon, clat])
landfrac = griddata(
    pft_points,
    landfrac_pft,
    scrip_points,
    method="nearest",
    fill_value=0.0,
    rescale=False,
)


# 超出 PFT 范围或插值得 nan 的格点置 0
if ext_pft:
    outside_pft = (
        (clon < ext_pft["lon_min"])
        | (clon > ext_pft["lon_max"])
        | (clat < ext_pft["lat_min"])
        | (clat > ext_pft["lat_max"])
    )
else:
    outside_pft = np.zeros(clon.shape, dtype=bool)
landfrac = np.atleast_1d(np.asarray(landfrac, dtype=np.float64))
landfrac[outside_pft] = 0.0
np.nan_to_num(landfrac, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
landfrac = np.clip(landfrac, 0.0, 1.0)
mask = (landfrac > 0).astype(np.int32)

# ---- reshape to (nj, ni) ----
nlon, nlat = g["grid_dims"].values
nj, ni = nlat, nlon

# %% quick map check -----------------------------------------------------------------
mask_map = xr.DataArray(
    mask.reshape((nj, ni)),
    dims=("nj", "ni"),
    coords={
        "lon": (("nj", "ni"), clon.reshape((nj, ni))),
        "lat": (("nj", "ni"), clat.reshape((nj, ni))),
    },
    name="mask",
)
landfrac_map = xr.DataArray(
    landfrac.reshape((nj, ni)),
    dims=("nj", "ni"),
    coords={
        "lon": (("nj", "ni"), clon.reshape((nj, ni))),
        "lat": (("nj", "ni"), clat.reshape((nj, ni))),
    },
    name="landfrac",
)
mask_map.plot(x="lon", y="lat", cmap="gray_r")
landfrac_map.plot(x="lon", y="lat")
# %% quick map check -----------------------------------------------------------------


def reshape2(a):
    return a.reshape((nj, ni))


ds = xr.Dataset(
    data_vars=dict(
        xc=(("nj", "ni"), reshape2(clon)),
        yc=(("nj", "ni"), reshape2(clat)),
        xv=(("nj", "ni", "nv"), lonc.reshape(nj, ni, 4)),
        yv=(("nj", "ni", "nv"), latc.reshape(nj, ni, 4)),
        area=(("nj", "ni"), reshape2(area)),
        frac=(("nj", "ni"), reshape2(landfrac)),
        mask=(("nj", "ni"), reshape2(mask)),
    ),
    coords=dict(
        ni=np.arange(ni),
        nj=np.arange(nj),
        nv=np.arange(4),
    ),
)

ds["area"].attrs["units"] = "radian2"
ds["frac"].attrs["units"] = "unitless"

ds.to_netcdf(OUT, format="NETCDF3_64BIT")

p.close()
g.close()
print("Domain written:", OUT)

# %%
# Prepare data for plotting
mask_plot = ds["mask"].values
landfrac_plot = ds["frac"].values

xc = ds["xc"].values
yc = ds["yc"].values

fig, axs = plt.subplots(1, 3, figsize=(18, 6))

im0 = axs[0].pcolormesh(xc, yc, mask_plot, shading="auto", cmap="gray")
axs[0].set_title("mask")
plt.colorbar(im0, ax=axs[0])

im1 = axs[1].pcolormesh(xc, yc, landfrac_plot, shading="auto", cmap="viridis")
axs[1].set_title("landfrac")
plt.colorbar(im1, ax=axs[1])

# Plot source PFT land fraction on its native grid/points.
valid = np.isfinite(pft_lon) & np.isfinite(pft_lat) & np.isfinite(landfrac_pft)
sc = axs[2].scatter(
    pft_lon[valid],
    pft_lat[valid],
    c=landfrac_pft[valid],
    s=1,
    cmap="viridis",
    vmin=0.0,
    vmax=1.0,
    rasterized=True,
)
axs[2].set_title("landfrac_pft (source)")
plt.colorbar(sc, ax=axs[2])

for ax in axs:
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")

plt.tight_layout()
plt.savefig(OUT + ".png", dpi=300)
plt.close(fig)
print("Figure saved: ", OUT + ".png")
print(
    "Domain file path: \n/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s3_Domainfiles"
)
