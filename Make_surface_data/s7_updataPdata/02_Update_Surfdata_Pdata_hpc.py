# HPC (Pathfinder) version of s7_updataPdata/02_Update_Surfdata_Pdata.py.
# Same logic as the mac script -- update LABILE_P, SECONDARY_P, APATITE_P,
# OCCLUDED_P in the s6-updated surfdata with the per-pixel values from
# 01_Generate_P_data.py -- only the path block differs.
#
# Target grid: SEUS 1/24°, (lsmlat=324, lsmlon=504). Write-back is
# pointwise on (lsmlat, lsmlon) after an explicit LATIXY / LONGXY
# alignment check.
# 06/02/2026

# %% load library
import datetime as _dt
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import rasterio
import xarray as xr
from rasterio.transform import from_origin
from scipy.spatial import cKDTree


def _write_seus_geotiff(
    arr: np.ndarray,
    latixy: np.ndarray,
    longxy: np.ndarray,
    out_path: Path | str,
    dtype: str = "float32",
    nodata: float | int | None = -9999.0,
) -> None:
    """
    Write a (lsmlat, lsmlon) array as a georeferenced GeoTIFF (EPSG:4326)
    so it can be loaded straight into QGIS.

    Mirrors _write_seus_geotiff in 01_Generate_P_data.py but takes plain
    numpy arrays (we're working with surfdata values cast via .values, not
    DataArrays). The SEUS 1/24° grid is regular in lon/lat -- we collapse
    LATIXY / LONGXY to 1D and assemble a standard north-up GTiff.
    """
    latixy = np.asarray(latixy)
    longxy = np.asarray(longxy)
    lat1d = latixy[:, 0]  # increases south -> north
    lon1d = longxy[0, :]  # increases west  -> east
    if not np.allclose(latixy - lat1d[:, None], 0):
        raise ValueError("LATIXY is not constant along lsmlon; grid is not regular.")
    if not np.allclose(longxy - lon1d[None, :], 0):
        raise ValueError("LONGXY is not constant along lsmlat; grid is not regular.")

    dlat = float(np.mean(np.diff(lat1d)))
    dlon = float(np.mean(np.diff(lon1d)))
    west_edge = float(lon1d[0]) - dlon / 2.0
    north_edge = float(lat1d[-1]) + dlat / 2.0
    transform = from_origin(west_edge, north_edge, dlon, dlat)

    arr2d = np.asarray(arr, dtype=dtype)
    arr2d = np.flipud(arr2d)  # GTiff row 0 must be the north edge

    profile = {
        "driver": "GTiff",
        "height": arr2d.shape[0],
        "width": arr2d.shape[1],
        "count": 1,
        "dtype": dtype,
        "crs": "EPSG:4326",
        "transform": transform,
        "compress": "deflate",
        "tiled": True,
    }
    if nodata is not None:
        profile["nodata"] = nodata
        if np.issubdtype(arr2d.dtype, np.floating):
            arr2d = np.where(np.isnan(arr2d), nodata, arr2d)

    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(arr2d, 1)
    print(f"wrote {out_path}  ({arr2d.shape}, dtype={dtype})")


# %%
# Paths (Make_surface_data layout).
#
# Anchored to this script's directory so the script doesn't care where
# SLURM (or you) launches it from. Actual layout:
#
#   /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/  (THIS_DIR.parent)
#       s6_updateSoildata/
#           Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc   (s6 output -> s7 input)
#       s7_updataPdata/                                              (THIS_DIR)
#           01_Generate_P_data.py
#           02_Update_Surfdata_Pdata_hpc.py                          (this script)
#           run_generate_P.slurm   run_update_P.slurm
#           outputs/
#               P_forms_SEUS_1_24deg.nc                              (s7 step 01 output -> input here)
#           surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c<today>.nc (surfdata output -> here)
#           landuse.timeseries_..._UpdatedP_c<today>.nc              (landuse output -> here, only if UPDATE_LANDUSE_TIMESERIES=True)
#           intermediate/check_surfdata_P_update.png                 (sanity plot -> here)
THIS_DIR = Path(__file__).resolve().parent

SURF_IN = (
    THIS_DIR.parent
    / "s6_updateSoildata"
    / "Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc"
)
P_SRC = THIS_DIR / "outputs" / "P_forms_SEUS_1_24deg.nc"

_today = _dt.date.today().strftime("%y%m%d")
SURF_OUT = THIS_DIR / f"surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c{_today}.nc"

if SURF_OUT.resolve() == SURF_IN.resolve():
    raise ValueError("SURF_OUT must differ from SURF_IN (would overwrite source).")
for p, name in [(SURF_IN, "SURF_IN"), (P_SRC, "P_SRC")]:
    if not p.is_file():
        raise FileNotFoundError(f"{name} not found: {p}")

print(f"SURF_IN : {SURF_IN}")
print(f"P_SRC   : {P_SRC}")
print(f"SURF_OUT: {SURF_OUT}")

P_VARS = ("LABILE_P", "SECONDARY_P", "APATITE_P", "OCCLUDED_P")

# %%
# Open target surfdata (read-only) to extract the grid and land mask.
# We use LATIXY / LONGXY (2D, the authoritative geographic coords) and
# PFTDATA_MASK (1 = land, 0 = ocean / non-land) -- this is the same
# pattern used by s6_updateSoildata/update_Soil_surfdat.py.
with xr.open_dataset(SURF_IN) as ds_surf:
    latixy = ds_surf["LATIXY"].load()
    longxy = ds_surf["LONGXY"].load()
    if "PFTDATA_MASK" in ds_surf.variables:
        pft_mask = np.asarray(ds_surf["PFTDATA_MASK"].values, dtype=np.int32) == 1
    else:
        # If PFTDATA_MASK is missing we fall back to "all land", matching the
        # s6 behaviour. In SEUS the mask is always present, so this branch
        # is defensive only.
        pft_mask = np.ones(latixy.shape, dtype=bool)
        print("WARNING: PFTDATA_MASK missing -- treating every cell as land.")

    surf_p_orig = {v: ds_surf[v].load() for v in P_VARS}

nlat, nlon = latixy.shape
print(f"surfdata grid: ({nlat}, {nlon})  land cells = {int(pft_mask.sum())}")

# %%
# Open the P source file and verify it is on the exact same (lsmlat, lsmlon)
# SEUS 1/24° grid as the target surfdata. The P_forms file was built from
# a sibling surfdata, which shares the SEUS template with the current
# target; but we cannot assume identical grids -- we verify.
#
# The P_forms NetCDF stores P pools as plain (lsmlat, lsmlon) DataArrays
# WITHOUT LATIXY / LONGXY coords (see 01_Generate_P_data.py: variables are
# built via soil_order.copy(data=...) which drops the 2D geo coords). So
# we compare shapes directly, and additionally compare LATIXY / LONGXY only
# if the P file happens to carry them.
with xr.open_dataset(P_SRC) as ds_p:
    missing = [v for v in P_VARS if v not in ds_p.variables]
    if missing:
        raise KeyError(f"P_forms file is missing variables: {missing}")

    p_src = {v: ds_p[v].load() for v in P_VARS}

    for v in P_VARS:
        if p_src[v].shape != (nlat, nlon):
            raise ValueError(
                f"{v}: source shape {p_src[v].shape} != target {(nlat, nlon)}. "
                f"P_forms must be regenerated against the new surfdata grid."
            )

    if "LATIXY" in ds_p.variables and "LONGXY" in ds_p.variables:
        if not np.allclose(ds_p["LATIXY"].values, latixy.values, atol=1e-9):
            raise ValueError("LATIXY mismatch between P_forms and target surfdata.")
        if not np.allclose(ds_p["LONGXY"].values, longxy.values, atol=1e-9):
            raise ValueError("LONGXY mismatch between P_forms and target surfdata.")
        print("LATIXY / LONGXY match exactly between P_forms and target surfdata.")
    else:
        print(
            "P_forms file has no LATIXY / LONGXY coords; alignment relies on "
            "identical (lsmlat, lsmlon) shape (verified above)."
        )

# %%
# Fill NaN on land cells via nearest neighbour.
#
# The P_forms file may have NaN on a few coastal cells where SOIL_ORDER is
# defined but rockP_den_resampled (and hence TPS) ended up NaN. ELM cannot
# carry NaN in a surfdata field, and PFTDATA_MASK==1 declares the cell as
# land -- so we must give it a value. We borrow from the closest land cell
# that does have data, using cKDTree on (lat, lon). Same strategy as the
# rockP fill in 01_Generate_P_data.py, applied independently to each pool.
lat_arr = np.asarray(latixy.values)
lon_arr = np.asarray(longxy.values)
filled = {}
for v in P_VARS:
    arr = np.asarray(p_src[v].values, dtype=np.float64)
    need_fill = pft_mask & ~np.isfinite(arr)
    have_data = pft_mask & np.isfinite(arr)
    if need_fill.any() and have_data.any():
        tree = cKDTree(np.column_stack([lat_arr[have_data], lon_arr[have_data]]))
        _, nn = tree.query(np.column_stack([lat_arr[need_fill], lon_arr[need_fill]]))
        arr[need_fill] = arr[have_data][nn]
        print(
            f"{v}: filled {int(need_fill.sum())} land NaN cells via nearest neighbour"
        )
    elif need_fill.any():
        raise ValueError(f"{v}: all land cells are NaN; cannot find any donor.")
    # Force any non-land cells that are still NaN back to 0 -- they will
    # be re-masked by PFTDATA_MASK on the write, but having a finite value
    # in the buffer keeps downstream tooling happy.
    arr[~np.isfinite(arr)] = 0.0
    filled[v] = arr

# %%
# Make a fresh copy of the target surfdata, then overwrite only the 4 P
# pools in place via netCDF4 (r+). Same write pattern as s6 -- this avoids
# the round-trip-through-xarray that would otherwise rewrite attributes,
# compression and fill values for every variable in the file.
SURF_OUT.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(SURF_IN, SURF_OUT)
print(f"copied {SURF_IN.name} -> {SURF_OUT.name}")

with nc.Dataset(SURF_OUT, "r+") as ds:
    for v in P_VARS:
        var = ds.variables[v]
        if var.shape != (nlat, nlon):
            raise ValueError(f"{v}: surfdata var shape {var.shape} != {(nlat, nlon)}")

        old = np.asarray(var[:], dtype=np.float64)
        new = filled[v].astype(np.float64)
        out = old.copy()
        out[pft_mask] = new[pft_mask]
        var[:, :] = out.astype(var.dtype)

        # Stamp provenance on the variable so it's clear where the values
        # came from. Keep existing long_name / units intact.
        var.setncattr(
            "history",
            f"Updated {_today} from {P_SRC.name} (land cells only, "
            f"PFTDATA_MASK==1); see s7_updataPdata/02_Update_Surfdata_Pdata_hpc.py.",
        )

    # File-level history append.
    prev_hist = getattr(ds, "history", "")
    ds.history = (
        f"{_dt.datetime.now().isoformat(timespec='seconds')}: "
        f"02_Update_Surfdata_Pdata_hpc.py wrote {', '.join(P_VARS)} "
        f"from {P_SRC.name}.\n"
    ) + (prev_hist if prev_hist else "")

print(f"wrote {SURF_OUT}")

# %%
# Sanity stats: before vs after on land cells.
print(
    f"\n{'pool':<12s}  {'orig mean':>10s}  {'new mean':>10s}  "
    f"{'orig max':>10s}  {'new max':>10s}"
)
with xr.open_dataset(SURF_OUT) as ds_chk:
    for v in P_VARS:
        old_vals = surf_p_orig[v].values[pft_mask]
        new_vals = ds_chk[v].values[pft_mask]
        print(
            f"{v:<12s}  {old_vals.mean():>10.3f}  {new_vals.mean():>10.3f}  "
            f"{old_vals.max():>10.3f}  {new_vals.max():>10.3f}"
        )

# %%
# Export old / new P pools as GeoTIFFs for QGIS inspection.
# Two tiffs per pool (8 total): *_old.tif from the original surfdata,
# *_new.tif from the just-written SURF_OUT. nodata = -9999 so NaN cells
# (ocean per PFTDATA_MASK, plus any leftover _FillValue cells from the
# original surfdata) show up as no-data in QGIS instead of as numeric 0.
tif_dir = THIS_DIR / "intermediate"
tif_dir.mkdir(parents=True, exist_ok=True)
with xr.open_dataset(SURF_OUT) as ds_chk:
    for v in P_VARS:
        old = surf_p_orig[v].values.astype(np.float32)
        new = ds_chk[v].values.astype(np.float32)
        _write_seus_geotiff(
            old,
            latixy.values,
            longxy.values,
            tif_dir / f"Pform_{v}_old.tif",
            dtype="float32",
            nodata=-9999.0,
        )
        _write_seus_geotiff(
            new,
            latixy.values,
            longxy.values,
            tif_dir / f"Pform_{v}_new.tif",
            dtype="float32",
            nodata=-9999.0,
        )

# %%
# Quick visual check: old vs new vs (new - old), one row per P pool.
# Shared color scale per row makes the spatial pattern shift obvious; the
# diff column uses a diverging map centered at zero.
# On a compute node there's no display, so we don't call plt.show() here
# (the SLURM script also exports MPLBACKEND=Agg as a belt-and-braces).
with xr.open_dataset(SURF_OUT) as ds_chk:
    fig, axes = plt.subplots(
        len(P_VARS), 3, figsize=(15, 3.6 * len(P_VARS)), constrained_layout=True
    )
    if len(P_VARS) == 1:
        axes = axes[None, :]
    for i, v in enumerate(P_VARS):
        old = surf_p_orig[v].values.astype(float)
        new = ds_chk[v].values.astype(float)
        diff = new - old
        # Mask non-land for cleaner display.
        old_m = np.where(pft_mask, old, np.nan)
        new_m = np.where(pft_mask, new, np.nan)
        diff_m = np.where(pft_mask, diff, np.nan)

        vmax = np.nanmax([np.nanmax(old_m), np.nanmax(new_m)])
        vmin = 0.0
        dmax = float(np.nanmax(np.abs(diff_m))) if np.isfinite(diff_m).any() else 1.0

        m0 = axes[i, 0].pcolormesh(
            longxy, latixy, old_m, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax
        )
        axes[i, 0].set_title(f"{v} (original)")
        fig.colorbar(m0, ax=axes[i, 0], shrink=0.85)

        m1 = axes[i, 1].pcolormesh(
            longxy, latixy, new_m, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax
        )
        axes[i, 1].set_title(f"{v} (new)")
        fig.colorbar(m1, ax=axes[i, 1], shrink=0.85)

        m2 = axes[i, 2].pcolormesh(
            longxy,
            latixy,
            diff_m,
            shading="auto",
            cmap="RdBu_r",
            vmin=-dmax,
            vmax=dmax,
        )
        axes[i, 2].set_title(f"{v} (new - original)")
        fig.colorbar(m2, ax=axes[i, 2], shrink=0.85)

        for ax in axes[i]:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_aspect("equal", adjustable="box")

    out_png = THIS_DIR / "intermediate" / "check_surfdata_P_update.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")

# %%
# Propagate the same 4 P pools into the landuse.timeseries file.
#
# landuse.timeseries_..._c260415.nc carries APATITE_P, LABILE_P, OCCLUDED_P,
# SECONDARY_P as time-invariant (lsmlat, lsmlon) fields on the exact same
# SEUS 1/24° grid as the surfdata. We reuse the already-NN-filled `filled`
# dict (and the same `pft_mask`) so the two outputs are numerically
# identical on every land cell -- no second nearest-neighbour pass, no risk
# of drift.
#
# PF layout: the landuse.timeseries file sits one directory above this
# script (in surfdata_results/, alongside postprocsSoil_onPF and
# postprocsP_onPF). The updated copy is written next to the SURF_OUT.
#
# DISABLED by default. Inspection of canonical CLM/ELM inputs on Pathfinder
# shows the P pools in landuse.timeseries are uniformly zero, e.g.:
#   /projects/hpcl-cli185/proj-shared/xyk/inputdata/lnd/clm2/surfdata_map/
#       landuse.timeseries_360x720cru_hist_simyr1850-2015_c180220.nc
#   /projects/hpcl-cli185/proj-shared/zdr/trendy_2024/landuse/
#       landuse.timeseries_360x720cru_hist_TRENDY_simyr1700-2024_c240826.nc
# i.e. ELM consumes P pools from surfdata, not from landuse.timeseries.
# We leave this block in place so that if that assumption ever changes,
# flipping UPDATE_LANDUSE_TIMESERIES to True is the only edit required.
UPDATE_LANDUSE_TIMESERIES = False

if UPDATE_LANDUSE_TIMESERIES:
    LU_IN = (
        THIS_DIR.parent
        / "landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260415.nc"
    )
    LU_OUT = (
        THIS_DIR
        / f"landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_UpdatedP_c{_today}.nc"
    )
    if LU_OUT.resolve() == LU_IN.resolve():
        raise ValueError("LU_OUT must differ from LU_IN (would overwrite source).")
    if not LU_IN.is_file():
        raise FileNotFoundError(f"LU_IN not found: {LU_IN}")
    print(f"\nLU_IN  : {LU_IN}")
    print(f"LU_OUT : {LU_OUT}")

    # Grid alignment check -- the surfdata and landuse.timeseries should both
    # be the SEUS 1/24° template, but we verify rather than trust the filename.
    with xr.open_dataset(LU_IN) as ds_lu:
        if not np.allclose(ds_lu["LATIXY"].values, latixy.values, atol=1e-9):
            raise ValueError(
                "LATIXY mismatch between landuse.timeseries and surfdata."
            )
        if not np.allclose(ds_lu["LONGXY"].values, longxy.values, atol=1e-9):
            raise ValueError(
                "LONGXY mismatch between landuse.timeseries and surfdata."
            )
        lu_p_orig = {v: ds_lu[v].load() for v in P_VARS}
    print("LATIXY / LONGXY match between landuse.timeseries and surfdata.")

    LU_OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LU_IN, LU_OUT)
    print(f"copied {LU_IN.name} -> {LU_OUT.name}")

    with nc.Dataset(LU_OUT, "r+") as ds:
        for v in P_VARS:
            var = ds.variables[v]
            if var.shape != (nlat, nlon):
                raise ValueError(
                    f"{v}: landuse.timeseries var shape {var.shape} != {(nlat, nlon)}"
                )
            old = np.asarray(var[:], dtype=np.float64)
            new = filled[v].astype(np.float64)
            out = old.copy()
            out[pft_mask] = new[pft_mask]
            var[:, :] = out.astype(var.dtype)
            var.setncattr(
                "history",
                f"Updated {_today} from {P_SRC.name} (land cells only, "
                f"PFTDATA_MASK==1); see s7_updataPdata/02_Update_Surfdata_Pdata_hpc.py.",
            )

        prev_hist = getattr(ds, "history", "")
        ds.history = (
            f"{_dt.datetime.now().isoformat(timespec='seconds')}: "
            f"02_Update_Surfdata_Pdata_hpc.py wrote {', '.join(P_VARS)} "
            f"from {P_SRC.name}.\n"
        ) + (prev_hist if prev_hist else "")

    print(f"wrote {LU_OUT}")

    # Sanity stats on landuse.timeseries P pools, on land cells only.
    print(
        f"\n[landuse.timeseries] {'pool':<12s}  {'orig mean':>10s}  "
        f"{'new mean':>10s}  {'orig max':>10s}  {'new max':>10s}"
    )
    with xr.open_dataset(LU_OUT) as ds_chk:
        for v in P_VARS:
            old_vals = lu_p_orig[v].values[pft_mask]
            new_vals = ds_chk[v].values[pft_mask]
            print(
                f"{v:<12s}  {old_vals.mean():>10.3f}  {new_vals.mean():>10.3f}  "
                f"{old_vals.max():>10.3f}  {new_vals.max():>10.3f}"
            )
else:
    print(
        "\nSkipping landuse.timeseries P-pool update "
        "(UPDATE_LANDUSE_TIMESERIES=False; ELM reads P pools from surfdata only)."
    )

# %%
