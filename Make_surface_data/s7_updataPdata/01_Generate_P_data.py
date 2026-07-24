# Generate P data based on README.md in s7_updataPdata folder
# reference Yang, X. et al. (2013).
# 05/28/2026


# %% load library
import os

import matplotlib.pyplot as plt
import numpy as np
import rasterio
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap
from rasterio.transform import from_origin
from scipy.spatial import cKDTree

# %%
# Step 1: Calculate TPS
# TPS = 0.01 · D · ρ_P · C_P · (1 − PPDI) / (ε + 1)

#  0.01 ·D · ρp · CP = rockP_den(g P m⁻²)
# read rockP_den.nc unit (kg P/m²)
# To convert to the g P m⁻² used in Yang et al. (2013) Eq. 2:
# g P m⁻² = rockP_den (kg P m⁻²) × 1000
# read rockP_den.nc; .load() materializes into memory so we can close the file
with xr.open_dataset("/projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/P_rock_den.nc") as ds:
    rockP_den = ds["rockP_den"].load() * 1000  # g P m⁻²

# %%
# Resample rockP_den (g P/m², global 0.5°) onto the 1/24° SEUS surfdata grid.
#
# surfdata stores the grid as two 2D variables LATIXY(lsmlat, lsmlon) and
# LONGXY(lsmlat, lsmlon), not as 1D coords. The dim-only labels `lsmlat` /
# `lsmlon` are just virtual integer indices (0..323, 0..503) and have no
# physical meaning, so we cannot pass them to .interp directly — doing so
# would map them as "lat=0..323", fall outside the source range, and return
# all-NaN. We must use the real LATIXY / LONGXY values as the interpolation
# target. Passing them as 2D DataArrays with the same dims as SOIL_ORDER
# makes xarray do pointwise interpolation and return a result on (lsmlat,
# lsmlon).
surfdata_path = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/"
    "Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc"
)
with xr.open_dataset(surfdata_path) as ds:
    soil_order = ds["SOIL_ORDER"].load()
    target_lat = ds["LATIXY"].load()  # 2D, (lsmlat, lsmlon), -180..180 lon frame
    target_lon = ds["LONGXY"].load()  # 2D, (lsmlat, lsmlon)
# %%
# Linear interp is appropriate for a continuous field (rockP_den, g P/m²).
# Do NOT use this for SOIL_ORDER itself — categorical integers need "nearest".
#
# kwargs={"fill_value": None} → scipy.interpolate.RegularGridInterpolator
# extrapolates instead of returning NaN at out-of-range points. We want this
# because some SEUS coastal land cells (Florida Keys, Outer Banks, Mississippi
# Delta) sit right on the 0.5° land/ocean boundary; without extrapolation
# they would get NaN from rockP_den and lose all P. Linear extrapolation
# from the nearest valid cells is the least-bad fallback for those cells.
# (To suppress extrapolation instead, use kwargs={"fill_value": np.nan}.)
rockP_den_resampled = rockP_den.interp(
    lat=target_lat,
    lon=target_lon,
    method="linear",
    kwargs={"fill_value": None},
).assign_coords(LATIXY=target_lat, LONGXY=target_lon)

# %%
# Fill rockP_den_resampled where SOIL_ORDER says "land" but rockP is NaN.
#
# Why we need this: the 0.5° rockP_den.nc has NaN on every ocean cell. After
# resampling onto SEUS 1/24° many coastal land cells (Florida peninsula,
# Mississippi Delta, Outer Banks, Chesapeake) inherit NaN from those ocean
# neighbours even though SOIL_ORDER (derived from a 1 km USDA soil map) does
# carry a real value there. SOIL_ORDER is our authoritative land mask: any
# pixel with a real USDA soil order (codes 1..12) is a land cell that must
# have a rockP value. Codes 0 (Water), 13 (shiftingsand), 14 (rockland),
# 15 (iceglacier) and the -9999 fill are treated as non-soil / no-data.
# See SOIL_ORDER_NAMES below for the full code → name lookup.
#
# Strategy: nearest-neighbour fill in geographic space — for each land cell
# whose rockP is NaN, copy the value from the closest land cell that does
# have rockP. cKDTree on (lat, lon) handles this in one vectorised query.
SOIL_ORDER_NONLAND_CODES = {0, 13, 14, 15}  # Water, shiftingsand, rockland, iceglacier
SOIL_ORDER_FILL = -9999

so_arr = np.asarray(soil_order.values)
land_mask = (
    np.isfinite(so_arr)
    & ~np.isin(so_arr, list(SOIL_ORDER_NONLAND_CODES))
    & (so_arr != SOIL_ORDER_FILL)
)

rp_arr = np.asarray(rockP_den_resampled.values, dtype=np.float64)
need_fill = land_mask & ~np.isfinite(rp_arr)
have_data = land_mask & np.isfinite(rp_arr)

print(
    f"land cells           : {int(land_mask.sum()):>7d}\n"
    f"land cells with rockP: {int(have_data.sum()):>7d}\n"
    f"land cells to fill   : {int(need_fill.sum()):>7d}"
)

if need_fill.any() and have_data.any():
    lat_arr = np.asarray(target_lat.values)
    lon_arr = np.asarray(target_lon.values)
    tree = cKDTree(np.column_stack([lat_arr[have_data], lon_arr[have_data]]))
    _, nn = tree.query(np.column_stack([lat_arr[need_fill], lon_arr[need_fill]]))
    rp_arr[need_fill] = rp_arr[have_data][nn]
    rockP_den_resampled = rockP_den_resampled.copy(data=rp_arr)
    print(f"filled {int(need_fill.sum())} coastal land cells via nearest neighbour")
elif need_fill.any():
    print("WARNING: no valid rockP cells available to draw nearest neighbours from")

# %%
# Quick visual check: global rockP_den (source, 0.5°), the SEUS-resampled
# version (1/24°), and SOIL_ORDER on the same SEUS grid. Useful to confirm
# the longitude convention / extent are right and there are no all-NaN holes.
fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

# 1) Global source rockP_den, with SEUS bounding box overlay
seus_lon_min = float(target_lon.min())
seus_lon_max = float(target_lon.max())
seus_lat_min = float(target_lat.min())
seus_lat_max = float(target_lat.max())
m0 = axes[0].pcolormesh(
    rockP_den["lon"], rockP_den["lat"], rockP_den, shading="auto", cmap="viridis"
)
axes[0].plot(
    [seus_lon_min, seus_lon_max, seus_lon_max, seus_lon_min, seus_lon_min],
    [seus_lat_min, seus_lat_min, seus_lat_max, seus_lat_max, seus_lat_min],
    color="red",
    lw=1.2,
)
axes[0].set_title("rockP_den (global, 0.5°)\ng P m⁻²")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")
fig.colorbar(m0, ax=axes[0], shrink=0.85)

# 2) Resampled rockP_den on SEUS 1/24° grid
m1 = axes[1].pcolormesh(
    target_lon,
    target_lat,
    rockP_den_resampled,
    shading="auto",
    cmap="viridis",
    vmin=float(rockP_den_resampled.min(skipna=True)),
    vmax=float(rockP_den_resampled.max(skipna=True)),
)
axes[1].set_title("rockP_den_resampled (SEUS, 1/24°)\ng P m⁻²")
axes[1].set_xlabel("Longitude")
axes[1].set_ylabel("Latitude")
axes[1].set_aspect("equal", adjustable="box")
fig.colorbar(m1, ax=axes[1], shrink=0.85)

# 3) SOIL_ORDER on SEUS grid, discrete categorical colors
so_vals = np.rint(soil_order.values).astype(np.int64)
so_vals = so_vals[np.isfinite(soil_order.values)]
classes = np.unique(so_vals) if so_vals.size else np.array([0])
base = plt.get_cmap("tab20", max(classes.size, 1))
cmap = ListedColormap([base(i % base.N) for i in range(classes.size)])
boundaries = np.concatenate([classes - 0.5, [classes[-1] + 0.5]])
norm = BoundaryNorm(boundaries, cmap.N)
m2 = axes[2].pcolormesh(
    target_lon, target_lat, soil_order, shading="auto", cmap=cmap, norm=norm
)
axes[2].set_title("SOIL_ORDER (SEUS, 1/24°)\nUSDA class")
axes[2].set_xlabel("Longitude")
axes[2].set_ylabel("Latitude")
axes[2].set_aspect("equal", adjustable="box")
cb = fig.colorbar(m2, ax=axes[2], ticks=classes.tolist(), shrink=0.85)
cb.ax.set_yticklabels([str(int(c)) for c in classes])

plt.show()
os.makedirs("intermediate", exist_ok=True)
fig.savefig(
    "intermediate/check_rockP_den_vs_resampled_vs_soil_order.png",
    dpi=200,
    bbox_inches="tight",
)
plt.close(fig)


# %%
# Export rockP_den_resampled and soil_order as GeoTIFFs for QGIS inspection.
#
# surfdata's LATIXY / LONGXY are stored as 2D but for the SEUS file the grid
# is actually a regular rectangular lon/lat grid (LATIXY varies only with
# lsmlat, LONGXY only with lsmlon). That means we can collapse them to 1D and
# write a standard georeferenced GeoTIFF (EPSG:4326).
# GeoTIFF row 0 must be the *north* edge, so we flip the array vertically
# before writing and use a negative pixel-height in the affine transform.
def _write_seus_geotiff(
    da: xr.DataArray,
    latixy: xr.DataArray,
    longxy: xr.DataArray,
    out_path: str,
    dtype: str = "float32",
    nodata: float | int | None = None,
) -> None:
    lat1d = latixy.values[:, 0]  # increases south -> north
    lon1d = longxy.values[0, :]  # increases west  -> east

    if not np.allclose(latixy.values - lat1d[:, None], 0):
        raise ValueError("LATIXY is not constant along lsmlon; grid is not regular.")
    if not np.allclose(longxy.values - lon1d[None, :], 0):
        raise ValueError("LONGXY is not constant along lsmlat; grid is not regular.")

    dlat = float(np.mean(np.diff(lat1d)))
    dlon = float(np.mean(np.diff(lon1d)))
    west_edge = float(lon1d[0]) - dlon / 2.0
    north_edge = float(lat1d[-1]) + dlat / 2.0
    transform = from_origin(west_edge, north_edge, dlon, dlat)

    arr = np.asarray(da.values, dtype=dtype)
    arr = np.flipud(arr)  # so row 0 is the northernmost row

    profile = {
        "driver": "GTiff",
        "height": arr.shape[0],
        "width": arr.shape[1],
        "count": 1,
        "dtype": dtype,
        "crs": "EPSG:4326",
        "transform": transform,
        "compress": "deflate",
        "tiled": True,
    }
    if nodata is not None:
        profile["nodata"] = nodata
        if np.issubdtype(arr.dtype, np.floating):
            arr = np.where(np.isnan(arr), nodata, arr)

    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(arr, 1)
    print(f"wrote {out_path}  ({arr.shape}, dtype={dtype})")


os.makedirs("intermediate", exist_ok=True)

_write_seus_geotiff(
    rockP_den_resampled,
    target_lat,
    target_lon,
    "intermediate/rockP_den_resampled_SEUS_1_24deg.tif",
    dtype="float32",
    nodata=-9999.0,
)

_write_seus_geotiff(
    soil_order,
    target_lat,
    target_lon,
    "intermediate/SOIL_ORDER_SEUS_1_24deg.tif",
    dtype="int16",
    nodata=-9999,
)
# %%
# ---------------------------------------------------------
# TPS = 0.01 · D · ρ_P · C_P · (1 − PPDI) / (ε + 1)
# calculate  (1 − PPDI) and (ε + 1)
# ---------------------------------------------------------
# %%
# ELM soil order map to USDA soil order name
# Integer code -> USDA soil order name used in the SOIL_ORDER raster
# (see s6_updateSoildata/update_Soil_surfdat.py SOIL_ORDER source TIF
# legend). 0 = water; 13–15 are non-soil cover classes that should be
# treated as no-data for the PPDI / TPS calculation.
SOIL_ORDER_NAMES = {
    0: "Water",
    1: "Andisols",
    2: "Gelisols",
    3: "Histosols",
    4: "Entisols",
    5: "Inceptisols",
    6: "Aridsols",
    7: "Vertisols",
    8: "Mollisols",
    9: "Alfisols",
    10: "Spodosols",
    11: "Ultisols",
    12: "Oxisols",
    13: "shiftingsand",
    14: "rockland",
    15: "iceglacier",
}

# Inverse lookup (name -> integer code).
SOIL_ORDER_CODES = {name: code for code, name in SOIL_ORDER_NAMES.items()}

# make PPDI table based on Weathering stages
# Weathering stages based on soil order
# USDA soil orders grouped by weathering stage — Yang et al. (2013), Sect. 3.1
# (Gelisols, Histosols, Andisols are excluded as "special" soil orders.)
SOIL_ORDERS_BY_STAGE = {
    "Slight": ["Entisols", "Inceptisols"],
    "Intermediate": ["Aridsols", "Mollisols", "Vertisols", "Alfisols"],
    "Strong": ["Spodosols", "Ultisols", "Oxisols"],
}

# PPDI (%, mean) by weathering stage / soil order — Yang et al. (2013), Table 2
# Uniform nested format: every stage maps to a {soil_order: PPDI} dict, even
# when all values within a stage are equal. This keeps the lookup loop
# simple (always PPDI_PCT[stage][order]) and makes it trivial to refine a
# single order later without restructuring the dict.
PPDI_PCT = {
    "Slight": {
        "Entisols": 12,
        "Inceptisols": 12,
    },
    "Intermediate": {
        "Aridsols": 54,
        "Mollisols": 54,
        "Vertisols": 54,
        "Alfisols": 54,
    },
    "Strong": {
        "Spodosols": 65,
        "Ultisols": 70,
        "Oxisols": 90,
    },
}


# Table 3 — Strain values (ε, mean) in 50 cm soils, Yang et al. (2013).
# Slight weathering splits by parent material; Strong splits by soil order
# (Spodosols vs Ultisols & Oxisols). ± sd dropped.
#
# Parent-material assumption for "Slight" (Non-carbonate vs Carbonate):
# SOIL_ORDERS_BY_STAGE cannot distinguish carbonate from non-carbonate parent
# material — that information comes from a lithology map (Yang et al. 2013
# used Dürr et al. 2005 1 km lithology), which is not available here (see
# README, Step A 2026-05-19). For the SE US domain almost all slightly
# weathered cells (Entisols / Inceptisols) sit on coastal-plain siliciclastic
# sediments and Piedmont/Appalachian saprolite, i.e. non-carbonate parent
# material. Known carbonate exceptions (the Florida limestone platform and
# the Appalachian Valley & Ridge limestone valleys) are either small or
# dominated by Histosols, which are excluded from the PPDI calculation
# anyway. We therefore default every "Slight" cell to Non-carbonate
# (ε = 0.66) and revisit with a real lithology mask (e.g. GLiM) later if
# the TPS map shows artifacts in those regions.
DEFAULT_SLIGHT_PARENT_MATERIAL = "Non_carbonate"
# Uniform nested format like PPDI_PCT: every stage maps to a dict whose
# inner key is the *soil order* — except "Slight", whose inner key is
# parent-material lithology (Non_carbonate / Carbonate) because Yang et al.
# (2013) Table 3 only resolves Slight ε that way. The pixel-level lookup
# must therefore special-case "Slight".
STRAIN_EPSILON = {
    "Slight": {
        "Non_carbonate": 0.66,
        "Carbonate": -0.63,
    },
    "Intermediate": {
        "Aridsols": -0.06,
        "Mollisols": -0.06,
        "Vertisols": -0.06,
        "Alfisols": -0.06,
    },
    "Strong": {
        "Spodosols": -0.13,
        "Ultisols": -0.12,
        "Oxisols": -0.12,
    },
}

# ---------------------------------------------------------
# %%
# Build per-pixel PPDI / ε / weathering-stage maps by looking up each
# SOIL_ORDER integer code through the four reference dicts:
#   code → SOIL_ORDER_NAMES → SOIL_ORDERS_BY_STAGE → PPDI_PCT / STRAIN_EPSILON
#
# Strategy: build a length-16 numpy lookup table indexed by the SOIL_ORDER
# code, then evaluate the whole 324×504 map in one fancy-indexing call
# (vectorised in C, runs in ~5 ms instead of seconds).
#
# For every code that is NOT one of the 9 Yang et al. (2013) weathering-stage
# soil orders (i.e. Water, Andisols, Gelisols, Histosols, shiftingsand,
# rockland, iceglacier), we set PPDI=0 and ε=0 → Eq. 2 collapses to
# TPS = rockP_den, treating the parent-material rockP as the final P
# density without pedogenic depletion or strain correction.
#
# Slight weathering's ε depends on parent-material lithology
# (Non_carbonate vs Carbonate); see the comment by STRAIN_EPSILON for why
# we default to Non_carbonate for the whole SEUS domain (Option 3).

STAGE_NAMES = ("Slight", "Intermediate", "Strong", "Other")
STAGE_OTHER = 3  # int code used in stage_map for any non-stage soil

code_to_ppdi = np.zeros(16, dtype=np.float32)
code_to_eps = np.zeros(16, dtype=np.float32)
code_to_stage = np.full(16, STAGE_OTHER, dtype=np.int8)

for stage_int, stage_name in enumerate(("Slight", "Intermediate", "Strong")):
    for order_name in SOIL_ORDERS_BY_STAGE[stage_name]:
        code = SOIL_ORDER_CODES[order_name]
        code_to_stage[code] = stage_int
        code_to_ppdi[code] = PPDI_PCT[stage_name][order_name]
        # Slight ε keyed by parent material, all other stages keyed by order.
        if stage_name == "Slight":
            code_to_eps[code] = STRAIN_EPSILON["Slight"][DEFAULT_SLIGHT_PARENT_MATERIAL]
        else:
            code_to_eps[code] = STRAIN_EPSILON[stage_name][order_name]

so_int = np.rint(np.asarray(soil_order.values)).astype(np.int64)
so_int = np.clip(so_int, 0, 15)  # guard against rogue codes / fill values

ppdi_arr = code_to_ppdi[so_int]
eps_arr = code_to_eps[so_int]
stage_arr = code_to_stage[so_int]

ppdi_pct_map = soil_order.copy(data=ppdi_arr).astype("float32")
ppdi_pct_map.attrs = {
    "long_name": "Phosphorus Pedogenic Depletion Index",
    "units": "percent",
    "reference": "Yang et al. (2013) Table 2",
}

eps_map = soil_order.copy(data=eps_arr).astype("float32")
eps_map.attrs = {
    "long_name": "Volumetric soil strain in top 50 cm",
    "units": "1",
    "reference": "Yang et al. (2013) Table 3",
    "note": "Slight uses Non_carbonate value (SEUS Option 3)",
}

stage_map = soil_order.copy(data=stage_arr).astype("int8")
stage_map.attrs = {
    "long_name": "Weathering stage class",
    "flag_values": [0, 1, 2, 3],
    "flag_meanings": "Slight Intermediate Strong Other",
}

# Sanity print: stage counts + unique values per map
print("stage histogram:")
for i, name in enumerate(STAGE_NAMES):
    n = int((stage_arr == i).sum())
    print(f"  {name:<12s}: {n:>7d} cells")
print(f"  {'total':<12s}: {stage_arr.size:>7d} cells")
print(f"ppdi unique values: {np.unique(ppdi_arr)}")
print(f"eps  unique values: {np.unique(eps_arr)}")

# Save as a single small NetCDF holding all three variables (good for
# scripted reuse) and as GeoTIFFs (good for QGIS overlay).
# Directory convention:
#   outputs/      = final products consumed downstream (TPS, P_forms).
#   intermediate/ = transient artifacts for inspection / sanity checking
#                   (rockP_resampled, ppdi/eps/stage maps, TIFFs, PNGs).
os.makedirs("intermediate", exist_ok=True)
xr.Dataset({"PPDI_pct": ppdi_pct_map, "eps": eps_map, "stage": stage_map}).to_netcdf(
    "intermediate/ppdi_eps_stage_SEUS_1_24deg.nc"
)

os.makedirs("intermediate", exist_ok=True)
_write_seus_geotiff(
    ppdi_pct_map,
    target_lat,
    target_lon,
    "intermediate/PPDI_pct_SEUS_1_24deg.tif",
    dtype="float32",
)
_write_seus_geotiff(
    eps_map,
    target_lat,
    target_lon,
    "intermediate/EPS_SEUS_1_24deg.tif",
    dtype="float32",
)
_write_seus_geotiff(
    stage_map,
    target_lat,
    target_lon,
    "intermediate/STAGE_SEUS_1_24deg.tif",
    dtype="int16",
)

# %%
# Calculate TPS (total P in top 50 cm soil, g P m⁻²) via Yang et al. (2013) Eq. 2:
#
#   TPS = 0.01 · D · ρ_P · C_P · (1 − PPDI/100) / (ε + 1)
#
# rockP_den_resampled already encodes 0.01 · D · ρ_P · C_P in g P m⁻²
# (P_rock_den source kg P m⁻² × 1000; see top of script). Combined with the
# per-pixel PPDI / ε maps built above, every pixel reduces to:
#
#   TPS_pixel = rockP_den_resampled · (1 − PPDI/100) / (ε + 1)
#
# Numerics:
# - For "Other" pixels (PPDI=0, ε=0) the formula collapses to rockP_den,
#   i.e. raw parent-material P density with no pedogenic depletion or strain.
# - For non-land pixels rockP_den_resampled is NaN, so TPS is NaN by NaN
#   propagation — no need to mask explicitly.
# - We sanity-check that (ε + 1) is never ≤ 0 (would otherwise blow up or
#   flip sign). Yang Table 3 ε values are all > -1 so this should hold, but
#   we assert it to catch any future data corruption.
eps_plus_1 = eps_map + 1.0
if float(eps_plus_1.min()) <= 0.0:
    raise ValueError(
        f"(ε + 1) has non-positive values (min={float(eps_plus_1.min()):.3f}); "
        f"TPS would be ill-defined."
    )

tps_map = (rockP_den_resampled * (1.0 - ppdi_pct_map / 100.0) / eps_plus_1).astype(
    "float32"
)
tps_map.attrs = {
    "long_name": "Total phosphorus in top 50 cm soil (TPS)",
    "units": "g P m-2",
    "reference": "Yang et al. (2013) Eq. 2",
    "formula": "rockP_den * (1 - PPDI/100) / (eps + 1)  "
    "[rockP_den already encodes 0.01 * D * rho_P * C_P in g P m-2]",
}

# Sanity stats per weathering stage
print("TPS (g P m⁻²) summary by weathering stage:")
print(f"  {'stage':<12s}  {'N':>7s}  {'min':>8s}  {'mean':>8s}  {'max':>8s}")
stage_arr_np = np.asarray(stage_map.values)
tps_np = np.asarray(tps_map.values)
for i, name in enumerate(STAGE_NAMES):
    mask = (stage_arr_np == i) & np.isfinite(tps_np)
    if mask.any():
        v = tps_np[mask]
        print(
            f"  {name:<12s}  {mask.sum():>7d}  "
            f"{v.min():>8.2f}  {v.mean():>8.2f}  {v.max():>8.2f}"
        )
    else:
        print(f"  {name:<12s}  {0:>7d}  {'--':>8s}  {'--':>8s}  {'--':>8s}")
nan_frac = float(np.isnan(tps_np).mean())
print(f"  NaN fraction (non-land + missing rockP): {nan_frac:.3f}")

# %%
# Save TPS: NetCDF (single-variable, in outputs/ as a final product) +
# GeoTIFF (in intermediate/ for QGIS).
os.makedirs("outputs", exist_ok=True)
os.makedirs("intermediate", exist_ok=True)
tps_map.to_dataset(name="TPS").to_netcdf("outputs/TPS_SEUS_1_24deg.nc")
_write_seus_geotiff(
    tps_map,
    target_lat,
    target_lon,
    "intermediate/TPS_SEUS_1_24deg.tif",
    dtype="float32",
)

# %%
# Quick visual sanity plot: rockP_den_resampled vs TPS (same color scale
# would distort because rockP is the pre-Eq.2 input; instead each gets its
# own scale). The two panels should look similar in spatial pattern, with
# TPS visibly reduced wherever PPDI is large (Spodosols/Ultisols/Oxisols
# corridor in SE US).
fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

rp_vmin = float(rockP_den_resampled.min(skipna=True))
rp_vmax = float(rockP_den_resampled.max(skipna=True))
m0 = axes[0].pcolormesh(
    target_lon,
    target_lat,
    rockP_den_resampled,
    shading="auto",
    cmap="viridis",
    vmin=rp_vmin,
    vmax=rp_vmax,
)
axes[0].set_title("rockP_den_resampled (g P m⁻²)\n[input to Eq. 2]")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")
axes[0].set_aspect("equal", adjustable="box")
fig.colorbar(m0, ax=axes[0], shrink=0.85)

tps_vmin = float(np.nanmin(tps_np))
tps_vmax = float(np.nanmax(tps_np))
m1 = axes[1].pcolormesh(
    target_lon,
    target_lat,
    tps_map,
    shading="auto",
    cmap="viridis",
    vmin=tps_vmin,
    vmax=tps_vmax,
)
axes[1].set_title("TPS (g P m⁻²)\n[Eq. 2 output: rockP·(1−PPDI/100)/(ε+1)]")
axes[1].set_xlabel("Longitude")
axes[1].set_ylabel("Latitude")
axes[1].set_aspect("equal", adjustable="box")
fig.colorbar(m1, ax=axes[1], shrink=0.85)

plt.show()
os.makedirs("intermediate", exist_ok=True)
fig.savefig("intermediate/check_TPS_vs_rockP.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# %%

# Step 2: Calculate P form fractions
# %%
# Hedley fractionation P fractions (mean, %) by USDA soil order.
# - 9 weathering-stage orders: Yang et al. (2013) Table S2 (± sd dropped).
# - Histosols: not in Yang 2013 Table S2. Values estimated from Fig. 3 of
#   Yang & Post (2011), "Phosphorus transformations as a function of
#   pedogenesis: A synthesis of soil phosphorus data using Hedley
#   fractionation method", Biogeosciences 8, 2907–2916,
#   doi:10.5194/bg-8-2907-2011. Read off bar heights; sums to ~100%.

# ELM soilorder 2:Gelisol(Permafrost),
# 1:Andisolnot(volcanic ash soils) are not in Yang 2013 Table S2
# They are not in SE US region, so we do not need to consider them.

HEDLEY_P_FRACTIONS_PCT = {
    #               LABILE_P, SECONDARY_P, APATITE_P, OCCLUDED_P, organic_P
    "Entisols": {
        "LABILE_P": 11,
        "SECONDARY_P": 5,
        "APATITE_P": 47,
        "OCCLUDED_P": 22,
        "organic_P": 15,
    },
    "Inceptisols": {
        "LABILE_P": 12,
        "SECONDARY_P": 7,
        "APATITE_P": 17,
        "OCCLUDED_P": 23,
        "organic_P": 41,
    },
    "Aridsols": {
        "LABILE_P": 8,
        "SECONDARY_P": 6,
        "APATITE_P": 64,
        "OCCLUDED_P": 17,
        "organic_P": 5,
    },
    "Vertisols": {
        "LABILE_P": 6,
        "SECONDARY_P": 6,
        "APATITE_P": 29,
        "OCCLUDED_P": 47,
        "organic_P": 12,
    },
    "Mollisols": {
        "LABILE_P": 5,
        "SECONDARY_P": 4,
        "APATITE_P": 28,
        "OCCLUDED_P": 44,
        "organic_P": 19,
    },
    "Alfisols": {
        "LABILE_P": 7,
        "SECONDARY_P": 11,
        "APATITE_P": 19,
        "OCCLUDED_P": 38,
        "organic_P": 25,
    },
    "Spodosols": {
        "LABILE_P": 7,
        "SECONDARY_P": 12,
        "APATITE_P": 9,
        "OCCLUDED_P": 28,
        "organic_P": 44,
    },
    "Ultisols": {
        "LABILE_P": 7,
        "SECONDARY_P": 14,
        "APATITE_P": 3,
        "OCCLUDED_P": 50,
        "organic_P": 26,
    },
    "Oxisols": {
        "LABILE_P": 6,
        "SECONDARY_P": 14,
        "APATITE_P": 1,
        "OCCLUDED_P": 59,
        "organic_P": 20,
    },
    # Histosols — from Yang & Post (2011) Fig. 3 (bar-height estimates);
    # sum = 10.1 + 5.8 + 8.3 + 40.8 + 35.0 = 100.0
    "Histosols": {
        "LABILE_P": 10.1,
        "SECONDARY_P": 5.8,
        "APATITE_P": 8.3,
        "OCCLUDED_P": 40.8,
        "organic_P": 35.0,
    },
}

# %%
# Build a length-16 × 5 lookup table indexed by SOIL_ORDER code.
# Same pattern as code_to_ppdi / code_to_eps, but the per-code value is a
# vector of 5 fraction percentages instead of a scalar. Codes that don't
# appear in HEDLEY_P_FRACTIONS_PCT (Water=0, Andisols=1, Gelisols=2 and the
# 3 non-soil cover classes 13/14/15) keep the default 0 row, so all 5
# P-form pools come out 0 for them — correct behavior since none of those
# cells exist as real soils in SE US.
P_FORM_NAMES = ("LABILE_P", "SECONDARY_P", "APATITE_P", "OCCLUDED_P", "organic_P")

code_to_fracs_pct = np.zeros((16, len(P_FORM_NAMES)), dtype=np.float32)
for order_name, fracs in HEDLEY_P_FRACTIONS_PCT.items():
    code = SOIL_ORDER_CODES[order_name]
    code_to_fracs_pct[code, :] = [fracs[p] for p in P_FORM_NAMES]

# %%
# Vectorised per-pixel multiplication:
#   P_form_pixel[k] = TPS_pixel * fraction_pct(code, k) / 100
# fancy indexing -> fracs_arr.shape = (lsmlat, lsmlon, 5)
fracs_arr = code_to_fracs_pct[so_int]  # (324, 504, 5)
p_form_arr = (
    tps_np[..., None] * fracs_arr / 100.0  # broadcast (324,504,1)*(324,504,5)
).astype("float32")  # (324, 504, 5)

# Wrap each pool in its own xr.DataArray so we get clean per-variable
# NetCDF / GeoTIFF outputs and downstream ELM-style naming.
p_form_maps = {}
for k, pname in enumerate(P_FORM_NAMES):
    da = soil_order.copy(data=p_form_arr[..., k]).astype("float32")
    da.attrs = {
        "long_name": f"Soil P fraction: {pname.replace('_', ' ')}",
        "units": "g P m-2",
        "reference": (
            "Yang et al. (2013) Table S2 (9 weathering-stage orders); "
            "Yang & Post (2011) Fig. 3 for Histosols."
        ),
        "formula": f"TPS * HEDLEY_P_FRACTIONS_PCT[order]['{pname}'] / 100",
    }
    p_form_maps[pname] = da

# Sanity: sum of 5 fractions should equal TPS within rounding (only on
# pixels where the soil order has a defined fraction table).
sum_p = sum(da.values for da in p_form_maps.values())
fraction_table_codes = np.array([SOIL_ORDER_CODES[n] for n in HEDLEY_P_FRACTIONS_PCT])
has_fracs = np.isin(so_int, fraction_table_codes) & np.isfinite(tps_np)
if has_fracs.any():
    diff = (sum_p - tps_np)[has_fracs]
    rel = np.abs(diff) / np.abs(tps_np[has_fracs] + 1e-9)
    print(
        f"sum(5 P forms) vs TPS on {int(has_fracs.sum())} fractioned cells: "
        f"max |abs diff| = {float(np.abs(diff).max()):.3g} g P/m², "
        f"max rel error = {float(rel.max()):.3g}"
    )

print("\nP form fraction summary (g P m⁻²) on fractioned land cells:")
print(f"  {'pool':<14s}  {'N':>7s}  {'min':>10s}  {'mean':>10s}  {'max':>10s}")
for pname, da in p_form_maps.items():
    v = da.values[has_fracs]
    print(
        f"  {pname:<14s}  {has_fracs.sum():>7d}  "
        f"{v.min():>10.2f}  {v.mean():>10.2f}  {v.max():>10.2f}"
    )

# %%
# Save: one combined NetCDF (in outputs/, the final product consumed
# downstream by the surfdata write-back step) with all 5 P-form variables +
# TPS for traceability, and 5 separate GeoTIFFs (in intermediate/) for QGIS.
os.makedirs("outputs", exist_ok=True)
os.makedirs("intermediate", exist_ok=True)
xr.Dataset({"TPS": tps_map, **p_form_maps}).to_netcdf("outputs/P_forms_SEUS_1_24deg.nc")
for pname, da in p_form_maps.items():
    _write_seus_geotiff(
        da,
        target_lat,
        target_lon,
        f"intermediate/Pform_{pname}_SEUS_1_24deg.tif",
        dtype="float32",
    )

# %%
# Sanity plot: TPS + 5 P pools on a 2×3 grid, shared color scale for the
# 5 pools so their relative magnitudes are visible at a glance.
fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
axes = axes.ravel()

# Panel 0: TPS
m = axes[0].pcolormesh(
    target_lon,
    target_lat,
    tps_map,
    shading="auto",
    cmap="viridis",
    vmin=float(np.nanmin(tps_np)),
    vmax=float(np.nanmax(tps_np)),
)
axes[0].set_title("TPS (g P m⁻²)\n[total soil P]")
fig.colorbar(m, ax=axes[0], shrink=0.85)

# Panels 1-5: 5 P pools with shared color scale
pool_vals = np.concatenate(
    [da.values[np.isfinite(da.values)].ravel() for da in p_form_maps.values()]
)
pool_vmin = float(np.nanmin(pool_vals)) if pool_vals.size else 0.0
pool_vmax = float(np.nanmax(pool_vals)) if pool_vals.size else 1.0
for i, (pname, da) in enumerate(p_form_maps.items(), start=1):
    m = axes[i].pcolormesh(
        target_lon,
        target_lat,
        da,
        shading="auto",
        cmap="viridis",
        vmin=pool_vmin,
        vmax=pool_vmax,
    )
    axes[i].set_title(f"{pname.replace('_', ' ')} (g P m⁻²)")
    fig.colorbar(m, ax=axes[i], shrink=0.85)

for ax in axes:
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="box")

plt.show()
os.makedirs("intermediate", exist_ok=True)
fig.savefig("intermediate/check_P_forms.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# %%
