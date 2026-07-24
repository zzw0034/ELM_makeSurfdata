"""
update_Soil_surfdat.py

Update soil / topography fields in an ELM surfdata NetCDF by re-deriving
them from 1 km source datasets. Runs after ``s5_mksurfdata_map``, which
produces a surfdata file where several soil / topo fields are still at
~0.5 deg resolution (and look like big squares on the map); this script
replaces them with values regridded from 1 km sources on land cells only.

The input surfdata is COPIED to a new file and then modified; the original
file is left untouched. Default output is ``surfdata_Updatedsoil_<basename>.nc``
next to the input. Use ``--output`` to override (pass the same path as
``--surfdata`` for an in-place edit). Variables are overwritten on land
cells only (see "Land-cell guarantee" below).
#_______________________________________________#
Inputs (sources)
Global 1 km soil / topo variables in NetCDF under a directory conventionally named global_cf_float (ESSD-style land parameters; see script header for literature and HPC paths).
    location:
    -mac:/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/Soil_Properties/data/global_cf_float
    -Pathfinder: /projects/hpcl-cli185/proj-shared/ywo/Soil_Properties/data/global_cf_float
North America 1 km soil order as a GeoTIFF (Soil_Order_NA_1km.tif by default).
    location:
    -mac:/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/Soil_Properties/data/Soil_Order_NA_1km.tif
    -Pathfinder: /projects/hpcl-cli185/proj-shared/ywo/Soil_Properties/data/Soil_Order_NA_1km.tif
ELM surface data (surfdata_*.nc)
    location:
    -mac:/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260415.nc
    -Pathfinder: /projects/hpcl-cli185/proj-shared/zw5/elm-surface-dataset/SEUS_1_24deg/surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260415.nc

Inputs
------
1. Global 1 km soil / topo NetCDFs under a directory conventionally named
   ``global_cf_float``.
       mac:        /Users/zw5/ORNL_workplace/.../global_cf_float
       Pathfinder: /projects/hpcl-cli185/proj-shared/ywo/.../global_cf_float
2. North America 1 km USDA soil order as a GeoTIFF
   (``Soil_Order_NA_1km.tif`` by default).
3. The ELM surfdata target (output of s5_mksurfdata_map) whose horizontal
   coords are ``LATIXY`` and ``LONGXY`` on dims ``(lsmlat, lsmlon)``.

Variables updated (source name -> surfdata name)
------------------------------------------------
    ELEVATION  -> TOPO         bilinear   (continuous, 1 km -> target)
    SLOPE      -> SLOPE        bilinear   (continuous)
    ORGANIC    -> ORGANIC      bilinear   (continuous, 10 layers)
    PCT_CLAY   -> PCT_CLAY     bilinear   (continuous, 10 layers, %)
    PCT_SAND   -> PCT_SAND     bilinear   (continuous, 10 layers, %)
    SOIL_ORDER -> SOIL_ORDER   mode       (categorical USDA 1-15; created if
                                           missing; water=0 declared as nodata
                                           so coastal land cells survive)
    STDEV_ELEV -> STD_ELEV     total-std  (sub-grid topographic heterogeneity:
                                           sqrt(mean(STDEV_ELEV_1km^2)
                                           + var(ELEVATION_1km)) inside each
                                           target cell)
    SOIL_COLOR -> SOIL_COLOR   N x N mode smoothing IN-PLACE on the existing
                                           surfdata field (the 1 km
                                           SOIL_COLOR_*.nc is NOT used;
                                           default window=3, configurable via
                                           --soil-color-window). This step is
                                           OPTIONAL (opt-out): enabled by
                                           default and toggled via
                                           --smooth-soil-color /
                                           --no-smooth-soil-color. When
                                           disabled, SOIL_COLOR is left
                                           untouched and --soil-color-window
                                           is ignored.

Pipeline summary
----------------
Bilinear / mode regrids use ``rioxarray.reproject_match`` (GDAL warp); the
clip_box extent is taken from PFTDATA_MASK == 1 cells only, so memory is
bounded by the target's land bbox and not the full source raster.

STDEV_ELEV uses the law of total variance. For each target cell, the
between-1km component is ``var(ELEVATION_1km)`` and the within-1km
component is ``mean(STDEV_ELEV_1km**2)``. The final value is the square root
of their sum. Both 1 km sources are sliced by integer index before binning
(the global rasters are never loaded in full).

SOIL_COLOR smoothing uses a vectorised N x N majority filter
(``scipy.ndimage.convolve`` on one-hot encoded classes), restricted to
land cells; the result stays on integer USDA colour classes.

Land-cell guarantee
-------------------
Write-back is land-only: cells with ``PFTDATA_MASK != 1`` keep the existing
surfdata values; only ``PFTDATA_MASK == 1`` cells are overwritten. For
every variable, after the warp / bin / filter step, any land cell that
came back NaN is filled from its nearest land cell with a value using
``scipy.ndimage.distance_transform_edt`` -- so every land cell is
guaranteed a value. Newly created SOIL_ORDER variables still use
``fill_value`` on non-land cells. If ``PFTDATA_MASK`` is absent, all
cells are treated as land and a warning is emitted.

Outputs
-------
- The output surfdata NetCDF (a copy of the input with the variables above
  updated on land cells).
- Optional diagnostic PNGs under ``SOIL_CHECK_PLOT_DIR`` if ``--plot`` is set:
      check_<varname>.png              original 1 km vs new surfgrid panels
      check_SOIL_COLOR_smoothing.png   3-panel before / after / changed cells
                                       (only written when SOIL_COLOR smoothing
                                       runs, i.e. --smooth-soil-color; skipped
                                       when --no-smooth-soil-color)

CLI examples
------------
# Default (writes ``surfdata_Updatedsoil_<basename>.nc`` next to the input)
cd /Users/zw5/ORNL_workplace/ELM_makeSurfdata
uv run python Make_surface_data/s6_updateSoildata/update_Soil_surfdat.py \
  --surfdata /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260415.nc \
  --global-cf /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/Soil_Properties/data/global_cf_float \
  --soil-order-tif /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/Soil_Properties/data/Soil_Order_NA_1km.tif \
  --plot \
  --soil-color-window 9 \
  --output /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata/surfdata_UpdatedSoil_SEUS_1_24deg_simyr1850_c260602.nc


# Explicit output path
... --output /path/to/my_output.nc

# In-place edit (overwrites the input!)
... --surfdata /a/x.nc --output /a/x.nc

# Bigger SOIL_COLOR smoothing window
... --soil-color-window 5

# Skip SOIL_COLOR smoothing entirely (leave SOIL_COLOR untouched;
# --soil-color-window is ignored)
... --no-smooth-soil-color

References
----------
- ESSD (1 km global soil/topo): Wang et al., 2024 -- "Global 1 km land
  surface parameters for kilometer-scale Earth system modeling",
  https://essd.copernicus.org/articles/16/2007/2024/
- OpenLandMap USDA soil order:
  https://developers.google.com/earth-engine/datasets/catalog/OpenLandMap_SOL_SOL_GRTGROUP_USDA-SOILTAX_C_v01#bands
- Adapted from Yaoping's process_merge_surfdat.py
  https://github.com/ypwong22/soil_properties/blob/main/process_merge_surfdat.py
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import warnings
from netCDF4 import Dataset
import numpy as np
import rasterio as rio
import rioxarray  # noqa: F401  (registers the .rio accessor on xarray)
import xarray as xr
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from scipy.ndimage import convolve, distance_transform_edt
from scipy.stats import binned_statistic_2d

# -----------------------------------------------------------------------------
# Defaults (CLI, or env SURFDATA_PATH, SOIL_GLOBAL_CF_DIR, SOIL_ORDER_TIF)
# -----------------------------------------------------------------------------
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_MAKE_SURFACE_DATA = os.path.normpath(os.path.join(_THIS_DIR, ".."))

DEFAULT_SURFDATA = os.path.join(
    _MAKE_SURFACE_DATA, "surfdata_small_SEUS_1_24deg_simyr1850_c260410.nc"
)
DEFAULT_PATH_GLOBAL_CF = os.path.join(
    _THIS_DIR, "Soil_Properties", "data", "global_cf_float"
)
DEFAULT_PATH_SOIL_ORDER_TIF = os.path.join(
    _THIS_DIR, "Soil_Properties", "data", "Soil_Order_NA_1km.tif"
)
DEFAULT_CHECK_PLOT_DIR = os.path.join(_THIS_DIR, "intermediate", "merged_surfdat")

# Source variable names on disk -> surfdata variable names.
# Specialised paths (NOT in this loop):
#   - SOIL_ORDER  : NA 1 km TIF + mode warp           (interpolate_soil_order_tif)
#   - STDEV_ELEV  : total std from 1 km mean/std topo (coarsen_stdev_elev_from_1km_stats)
#   - SOIL_COLOR  : in-place NxN mode smoothing       (smooth_soil_color_in_place)
SOIL_NC_VARS = [
    "ELEVATION",
    "ORGANIC",
    "PCT_CLAY",
    "PCT_SAND",
    "SLOPE",
]

# Per-variable GDAL warp resampling for `interpolate_1km_nc_to_targets`.
# All entries are continuous fields downsampled with bilinear. The three
# categorical / statistical variables above have their own functions and do
# not appear here.
_RESAMPLING_BY_VAR: dict[str, Resampling] = {
    "ELEVATION": Resampling.bilinear,
    "ORGANIC": Resampling.bilinear,
    "PCT_CLAY": Resampling.bilinear,
    "PCT_SAND": Resampling.bilinear,
    "SLOPE": Resampling.bilinear,
}


def _path_global_cf() -> str:
    return os.environ.get("SOIL_GLOBAL_CF_DIR", DEFAULT_PATH_GLOBAL_CF)


def _path_soil_order_tif() -> str:
    return os.environ.get("SOIL_ORDER_TIF", DEFAULT_PATH_SOIL_ORDER_TIF)


def _check_plot_dir() -> str:
    return os.environ.get("SOIL_CHECK_PLOT_DIR", DEFAULT_CHECK_PLOT_DIR)


def _default_output_path(input_path: str) -> str:
    """
    Default output filename: insert ``Updatedsoil_`` after a leading
    ``surfdata_`` prefix, keeping the directory and rest of the basename.

    Examples
    --------
    surfdata_SEUS_..._c260415.nc -> surfdata_Updatedsoil_SEUS_..._c260415.nc
    foo.nc                       -> foo_Updatedsoil.nc  (no surfdata_ prefix)
    """
    d, base = os.path.split(input_path)
    if base.startswith("surfdata_"):
        new_base = "surfdata_Updatedsoil_" + base[len("surfdata_") :]
    else:
        stem, ext = os.path.splitext(base)
        new_base = f"{stem}_Updatedsoil{ext}"
    return os.path.join(d, new_base)


# -----------------------------------------------------------------------------
# Variable / filename mapping helpers
# -----------------------------------------------------------------------------
def surfdata_var_name(src_name: str) -> str:
    if src_name == "STDEV_ELEV":
        return "STD_ELEV"
    if src_name == "ELEVATION":
        return "TOPO"
    return src_name


def _nc_filename_for_var(varname: str) -> str:
    if varname in ("PCT_CLAY", "PCT_SAND", "ORGANIC"):
        return f"{varname}_10layer_1k_c230606.nc"
    return f"{varname}_1k_c230606.nc"


# -----------------------------------------------------------------------------
# Regrid / aggregate pipelines for the four "from 1 km source" variables.
# (SOIL_COLOR smoothing is in the surfdata-I/O section below, since it reads
#  and writes the same surfdata file rather than a 1 km source.)
#   _fill_missing_land_with_edt         -- shared land-cell NaN fill helper
#   interpolate_1km_nc_to_targets       -- bilinear regrid for continuous fields
#   interpolate_soil_order_tif          -- mode regrid for categorical SOIL_ORDER
#   _bin_edges_from_centres             -- helper for binned_statistic_2d
#   coarsen_stdev_elev_from_1km_stats   -- total std from 1 km mean/std fields
# -----------------------------------------------------------------------------
def _fill_missing_land_with_edt(
    plane: np.ndarray, land_2d: np.ndarray, varname: str
) -> np.ndarray:
    """
    Build a (nlat, nlon) result where:
      * non-land cells are NaN,
      * land cells with a finite warped value keep it,
      * land cells with NaN are filled from the nearest land cell that did
        receive a value (Euclidean distance transform on the 2-D mask).
    """
    nlat, nlon = land_2d.shape
    out = np.full((nlat, nlon), np.nan, dtype=np.float64)
    out[land_2d] = plane[land_2d]
    have = land_2d & np.isfinite(out)
    need = land_2d & ~have
    if need.any():
        if not have.any():
            raise RuntimeError(
                f"{varname}: no land cell received a valid value from the "
                "resampled raster; check land_mask vs raster coverage."
            )
        _, (rr, cc) = distance_transform_edt(~have, return_indices=True)
        out[need] = out[rr[need], cc[need]]
    return out


def interpolate_1km_nc_to_targets(
    varname: str,
    lat_tgt: np.ndarray,
    lon_tgt: np.ndarray,
    *,
    grid_shape: tuple[int, int],
    land_mask: np.ndarray,
    path_src: str | None = None,
    pad_deg: float = 0.25,
    resampling: Resampling | None = None,
) -> np.ndarray:
    """
    Regrid one continuous variable from the 1 km NetCDF stack onto the
    surfdata grid using rioxarray's GDAL-warp backend.

    Used for the bilinear variables only (ELEVATION, ORGANIC, PCT_CLAY,
    PCT_SAND, SLOPE). STDEV_ELEV and SOIL_COLOR have their own pipelines
    (``coarsen_stdev_elev_from_1km_stats`` and ``smooth_soil_color_in_place``)
    and are NOT handled here.

    The resampling method is chosen per-variable from :data:`_RESAMPLING_BY_VAR`
    (all entries currently bilinear); pass ``resampling=`` to override.

    The clip_box extent is taken from the land cells only (much tighter than
    the full target bbox when land covers a small fraction of the rectangle).
    Every land cell is guaranteed to receive a value: any land cell that came
    back NaN after warp is filled from its nearest *land* neighbour that has
    a value (Euclidean distance transform on the 2-D mask). Non-land cells
    are returned as NaN.

    Returns a 1-D float64 array of length ``nlat*nlon`` for 2-D source
    variables, or a 2-D array of shape ``(nlayer, nlat*nlon)`` for 3-D
    source variables. Output is in C-order matching ``LATIXY.ravel()``.
    """
    path_src = path_src or _path_global_cf()
    filename = _nc_filename_for_var(varname)
    fpath = os.path.join(path_src, filename)
    if not os.path.isfile(fpath):
        raise FileNotFoundError(f"Missing soil property file: {fpath}")

    if resampling is None:
        resampling = _RESAMPLING_BY_VAR.get(varname, Resampling.bilinear)

    nlat, nlon = int(grid_shape[0]), int(grid_shape[1])
    n_pts = nlat * nlon

    lat_tgt = np.asarray(lat_tgt, dtype=np.float64).ravel()
    lon_tgt = np.asarray(lon_tgt, dtype=np.float64).ravel()
    if lat_tgt.size != n_pts or lon_tgt.size != n_pts:
        raise ValueError(
            f"lat/lon length {lat_tgt.size}/{lon_tgt.size} does not match "
            f"grid_shape {grid_shape}"
        )

    lat2d = lat_tgt.reshape(nlat, nlon)
    lon2d = lon_tgt.reshape(nlat, nlon)
    lat1d = lat2d[:, 0]
    lon1d = lon2d[0, :]
    if not (np.allclose(lat2d, lat1d[:, None]) and np.allclose(lon2d, lon1d[None, :])):
        raise ValueError(
            "LATIXY/LONGXY is not a regular lat/lon grid; "
            "rioxarray.reproject_match requires a rectilinear target."
        )

    land_2d = np.asarray(land_mask, dtype=bool)
    if land_2d.shape != (nlat, nlon):
        raise ValueError(
            f"land_mask shape {land_2d.shape} != grid_shape {(nlat, nlon)}"
        )

    # ---- Determine land bbox before touching the source file ---------------
    if not land_2d.any():
        # No land in this surfdata; return all NaN without opening the source.
        with xr.open_dataset(fpath, mask_and_scale=False, decode_times=False) as _ds:
            band_dims = [d for d in _ds[varname].dims if d not in ("lat", "lon")]
        if not band_dims:
            return np.full(n_pts, np.nan, dtype=np.float64)
        with xr.open_dataset(fpath, decode_times=False) as _ds:
            nband = int(_ds.sizes[band_dims[0]])
        return np.full((nband, n_pts), np.nan, dtype=np.float64)

    lat_min = float(lat2d[land_2d].min()) - pad_deg
    lat_max = float(lat2d[land_2d].max()) + pad_deg
    lon_min = float(lon2d[land_2d].min()) - pad_deg
    lon_max = float(lon2d[land_2d].max()) + pad_deg

    # ---- IO-defensive: integer .isel() FIRST, then promote to float64 ------
    # Using clip_box on a 3D *_10layer NC (6 GB) does NOT always push the
    # slice down to the netCDF backend reliably, so the full variable can
    # get promoted to float64 (-> tens of GB) and the OS kills the process
    # (exit 137 / SIGKILL). Pinning the IO with .isel() on raw integer
    # indices is the same idiom used in `coarsen_stdev_elev_from_1km_stats`.
    print(f"[{varname}] reading 1km slab from {os.path.basename(fpath)}...", flush=True)
    ds = xr.open_dataset(fpath, mask_and_scale=True, decode_times=False)
    try:
        src_lat = np.asarray(ds["lat"].values, dtype=np.float64)
        src_lon = np.asarray(ds["lon"].values, dtype=np.float64)
        ilat = np.where((src_lat >= lat_min) & (src_lat <= lat_max))[0]
        ilon = np.where((src_lon >= lon_min) & (src_lon <= lon_max))[0]
        if ilat.size == 0 or ilon.size == 0:
            raise ValueError(
                f"{varname}: target land extent is outside the source raster "
                f"footprint (source lat [{src_lat.min():.3f},{src_lat.max():.3f}], "
                f"lon [{src_lon.min():.3f},{src_lon.max():.3f}])."
            )
        i0, i1 = int(ilat[0]), int(ilat[-1]) + 1
        j0, j1 = int(ilon[0]), int(ilon[-1]) + 1
        sub_lat = src_lat[i0:i1]
        sub_lon = src_lon[j0:j1]
        var = ds[varname].isel(lat=slice(i0, i1), lon=slice(j0, j1))
        print(
            f"[{varname}]   slab shape {tuple(var.shape)} "
            f"(lat[{i0}:{i1}] lon[{j0}:{j1}])",
            flush=True,
        )
        non_spatial = [d for d in var.dims if d not in ("lat", "lon")]
        if len(non_spatial) > 1:
            raise ValueError(
                f"{varname}: more than one non-spatial dim ({non_spatial}); "
                "rioxarray.reproject_match supports at most one band dim."
            )
        band_dim = non_spatial[0] if non_spatial else None
        # Read only the small slab and promote to float64 in-memory.
        sub_vals = np.asarray(var.values, dtype=np.float64)
        if band_dim is not None:
            band_axis = list(var.dims).index(band_dim)
            if band_axis != 0:
                sub_vals = np.moveaxis(sub_vals, band_axis, 0)
            nband = int(sub_vals.shape[0])
    finally:
        ds.close()

    # ---- Rebuild as a CRS-aware in-memory DataArray ------------------------
    if band_dim is None:
        src = xr.DataArray(
            sub_vals,
            coords={"y": sub_lat, "x": sub_lon},
            dims=("y", "x"),
        )
    else:
        src = xr.DataArray(
            sub_vals,
            coords={"band": np.arange(nband), "y": sub_lat, "x": sub_lon},
            dims=("band", "y", "x"),
        )
    src = src.rio.write_crs("EPSG:4326")
    src = src.rio.write_nodata(np.nan, inplace=False)

    # ---- Target grid as a CRS-aware DataArray ------------------------------
    target = xr.DataArray(
        np.zeros((nlat, nlon), dtype=np.float32),
        coords={"y": lat1d, "x": lon1d},
        dims=("y", "x"),
    ).rio.write_crs("EPSG:4326")

    # ---- Warp (source is already the small clipped slab) -------------------
    print(
        f"[{varname}]   warping slab {tuple(sub_vals.shape)} "
        f"-> target {(nlat, nlon)} via rioxarray ({resampling.name})...",
        flush=True,
    )
    warped = src.rio.reproject_match(target, resampling=resampling)

    arr = np.asarray(warped.values, dtype=np.float64)

    if arr.ndim == 2:
        return _fill_missing_land_with_edt(arr, land_2d, varname).ravel()
    if arr.ndim == 3:
        nlayer = arr.shape[0]
        out_layers = np.empty((nlayer, n_pts), dtype=np.float64)
        for k in range(nlayer):
            out_layers[k] = _fill_missing_land_with_edt(
                arr[k], land_2d, f"{varname}[{k}]"
            ).ravel()
        return out_layers
    raise ValueError(f"{varname}: unexpected warped ndim {arr.ndim}")


def interpolate_soil_order_tif(
    lat_tgt: np.ndarray,
    lon_tgt: np.ndarray,
    *,
    grid_shape: tuple[int, int],
    land_mask: np.ndarray,
    tif_path: str | None = None,
    pad_deg: float = 0.25,
) -> np.ndarray:
    """
    Regrid the categorical SOIL_ORDER GeoTIFF (USDA 1-15, 0 = water) onto the
    surfdata grid using majority (mode) resampling.

    The clip_box extent is taken from the land cells only (much tighter than the
    full target bbox when land covers a small fraction of the rectangle). Every
    cell where ``land_mask`` is True is guaranteed to receive a valid soil order
    (any land cell that came back NoData after warp is filled from its nearest
    *land* neighbour that has a value). Non-land cells are returned as NaN so
    that :func:`ensure_soil_order_variable` can stamp them with ``fill_value``.

    Returns a 1-D float64 array of length ``nlat * nlon`` (C-order, matching
    ``LATIXY.ravel()``).
    """
    tif_path = tif_path or _path_soil_order_tif()
    if not os.path.isfile(tif_path):
        raise FileNotFoundError(f"Missing soil order raster: {tif_path}")

    nlat, nlon = int(grid_shape[0]), int(grid_shape[1])
    n_pts = nlat * nlon

    lat_tgt = np.asarray(lat_tgt, dtype=np.float64).ravel()
    lon_tgt = np.asarray(lon_tgt, dtype=np.float64).ravel()
    if lat_tgt.size != n_pts or lon_tgt.size != n_pts:
        raise ValueError(
            f"lat/lon length {lat_tgt.size}/{lon_tgt.size} does not match "
            f"grid_shape {grid_shape}"
        )

    lat2d = lat_tgt.reshape(nlat, nlon)
    lon2d = lon_tgt.reshape(nlat, nlon)
    lat1d = lat2d[:, 0]
    lon1d = lon2d[0, :]
    if not (np.allclose(lat2d, lat1d[:, None]) and np.allclose(lon2d, lon1d[None, :])):
        raise ValueError(
            "LATIXY/LONGXY is not a regular lat/lon grid; "
            "rioxarray.reproject_match requires a rectilinear target."
        )

    land_2d = np.asarray(land_mask, dtype=bool)
    if land_2d.shape != (nlat, nlon):
        raise ValueError(
            f"land_mask shape {land_2d.shape} != grid_shape {(nlat, nlon)}"
        )
    if not land_2d.any():
        # No land in this surfdata, nothing to regrid.
        return np.full(n_pts, np.nan, dtype=np.float64)

    # ---- Target grid as a CRS-aware xarray.DataArray ------------------------
    target = xr.DataArray(
        np.zeros((nlat, nlon), dtype=np.uint8),
        coords={"y": lat1d, "x": lon1d},
        dims=("y", "x"),
    ).rio.write_crs("EPSG:4326")

    # ---- Source: load TIF, mark 0 as nodata so mode ignores water -----------
    src = (
        rioxarray.open_rasterio(tif_path, masked=False)
        .squeeze("band", drop=True)
        .rio.write_nodata(0, inplace=False)
    )
    if src.rio.crs is None:
        src = src.rio.write_crs("EPSG:4326")

    # ---- clip_box extent driven by land cells (tighter than full bbox) ------
    lat_land = lat2d[land_2d]
    lon_land = lon2d[land_2d]
    lat_min = float(lat_land.min()) - pad_deg
    lat_max = float(lat_land.max()) + pad_deg
    lon_min = float(lon_land.min()) - pad_deg
    lon_max = float(lon_land.max()) + pad_deg

    # Intersect with source bounds (replacement for the old
    # `Window.intersection(0,0,src.width,src.height)` guard).
    s_left, s_bottom, s_right, s_top = src.rio.bounds()
    minx = max(lon_min, s_left)
    maxx = min(lon_max, s_right)
    miny = max(lat_min, s_bottom)
    maxy = min(lat_max, s_top)
    if not (minx < maxx and miny < maxy):
        # Target land sits entirely outside the source raster footprint.
        return np.full(n_pts, np.nan, dtype=np.float64)

    if (minx, miny, maxx, maxy) != (s_left, s_bottom, s_right, s_top):
        src = src.rio.clip_box(minx=minx, miny=miny, maxx=maxx, maxy=maxy)

    # ---- Mode resampling onto the target grid -------------------------------
    warped = src.rio.reproject_match(target, resampling=Resampling.mode)
    arr = np.asarray(warped.values, dtype=np.float64)
    # Anything not a valid USDA order (0, fill, NaN, negative) -> NaN
    arr[~np.isfinite(arr) | (arr <= 0)] = np.nan

    # ---- Mask + guarantee land cells have a value ---------------------------
    # Land cells whose 1 km footprint contained only water/nodata receive
    # NaN here and are filled from the nearest land cell with a real value
    # by the shared EDT helper.
    return _fill_missing_land_with_edt(arr, land_2d, "SOIL_ORDER").ravel()


def _bin_edges_from_centres(centres: np.ndarray) -> np.ndarray:
    """
    Build N+1 monotonically increasing bin edges from N (sorted-ascending)
    cell-centre coordinates by placing edges midway between consecutive
    centres; the two outer edges extend by half of the adjacent spacing.
    """
    c = np.asarray(centres, dtype=np.float64)
    if c.size < 2:
        raise ValueError("Need at least 2 centres to derive bin edges.")
    edges = np.empty(c.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (c[1:] + c[:-1])
    edges[0] = c[0] - 0.5 * (c[1] - c[0])
    edges[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    return edges


def coarsen_stdev_elev_from_1km_stats(
    lat_tgt: np.ndarray,
    lon_tgt: np.ndarray,
    *,
    grid_shape: tuple[int, int],
    land_mask: np.ndarray,
    path_src: str | None = None,
    pad_deg: float = 0.25,
) -> np.ndarray:
    """
    Sub-grid topographic heterogeneity from the 1 km mean elevation and
    1 km within-cell elevation standard deviation fields.

    For each target 1/24-deg cell, combine:

      * between-1km variability: variance of 1 km ELEVATION pixels whose
        centres fall inside the target cell, and
      * within-1km variability: mean of STDEV_ELEV_1km squared.

    This is the law of total variance:

        sigma_total^2  =  E[sigma_within^2]  +  Var[mean_within]
                       =  mean(STDEV_ELEV_1km^2)  +  Var(ELEVATION_1km)

    Do not average source STDEV_ELEV directly: mean(std) != std(samples).

    The land-cell bbox is read by integer slices before binning so the
    global 1 km files are never loaded in full. Bins with no valid paired
    ELEVATION / STDEV_ELEV samples are NaN and the EDT helper fills any
    land cell that came out NaN from its nearest land neighbour with a value.

    Returns
    -------
    np.ndarray
        1-D float64 of length nlat*nlon, in C-order matching ``LATIXY.ravel()``.
    """
    path_src = path_src or _path_global_cf()
    elev_filename = _nc_filename_for_var("ELEVATION")
    stdev_filename = _nc_filename_for_var("STDEV_ELEV")
    elev_fpath = os.path.join(path_src, elev_filename)
    stdev_fpath = os.path.join(path_src, stdev_filename)
    if not os.path.isfile(elev_fpath):
        raise FileNotFoundError(
            f"Missing 1 km ELEVATION file: {elev_fpath} (needed for STDEV_ELEV)"
        )
    if not os.path.isfile(stdev_fpath):
        raise FileNotFoundError(
            f"Missing 1 km STDEV_ELEV file: {stdev_fpath} (needed for STDEV_ELEV)"
        )

    nlat, nlon = int(grid_shape[0]), int(grid_shape[1])
    n_pts = nlat * nlon

    lat_tgt = np.asarray(lat_tgt, dtype=np.float64).ravel()
    lon_tgt = np.asarray(lon_tgt, dtype=np.float64).ravel()
    if lat_tgt.size != n_pts or lon_tgt.size != n_pts:
        raise ValueError(
            f"lat/lon length {lat_tgt.size}/{lon_tgt.size} does not match "
            f"grid_shape {grid_shape}"
        )

    lat2d = lat_tgt.reshape(nlat, nlon)
    lon2d = lon_tgt.reshape(nlat, nlon)
    lat1d = lat2d[:, 0]
    lon1d = lon2d[0, :]
    if not (np.allclose(lat2d, lat1d[:, None]) and np.allclose(lon2d, lon1d[None, :])):
        raise ValueError(
            "LATIXY/LONGXY is not a regular lat/lon grid; "
            "total-variance aggregation requires a rectilinear target."
        )

    land_2d = np.asarray(land_mask, dtype=bool)
    if land_2d.shape != (nlat, nlon):
        raise ValueError(
            f"land_mask shape {land_2d.shape} != grid_shape {(nlat, nlon)}"
        )
    if not land_2d.any():
        return np.full(n_pts, np.nan, dtype=np.float64)

    # binned_statistic_2d needs monotonically increasing bin edges. If the
    # target lat is descending (north -> south), flip for binning and flip
    # the result back at the end so it matches LATIXY's order.
    lat_descending = lat1d[0] > lat1d[-1]
    lat_for_bins = lat1d[::-1] if lat_descending else lat1d
    lon_descending = lon1d[0] > lon1d[-1]
    lon_for_bins = lon1d[::-1] if lon_descending else lon1d
    lat_edges = _bin_edges_from_centres(lat_for_bins)
    lon_edges = _bin_edges_from_centres(lon_for_bins)

    # Source ELEVATION and STDEV_ELEV, clipped to the land bbox + pad.
    lat_min = float(lat2d[land_2d].min()) - pad_deg
    lat_max = float(lat2d[land_2d].max()) + pad_deg
    lon_min = float(lon2d[land_2d].min()) - pad_deg
    lon_max = float(lon2d[land_2d].max()) + pad_deg

    ds_elev = xr.open_dataset(elev_fpath, mask_and_scale=True, decode_times=False)
    ds_stdev = xr.open_dataset(stdev_fpath, mask_and_scale=True, decode_times=False)
    try:
        src_lat = np.asarray(ds_elev["lat"].values, dtype=np.float64)
        src_lon = np.asarray(ds_elev["lon"].values, dtype=np.float64)
        stdev_lat = np.asarray(ds_stdev["lat"].values, dtype=np.float64)
        stdev_lon = np.asarray(ds_stdev["lon"].values, dtype=np.float64)
        if not (
            src_lat.shape == stdev_lat.shape
            and src_lon.shape == stdev_lon.shape
            and np.allclose(src_lat, stdev_lat)
            and np.allclose(src_lon, stdev_lon)
        ):
            raise ValueError(
                "STDEV_ELEV: ELEVATION and STDEV_ELEV 1 km files do not "
                "share the same lat/lon grid; cannot combine variances."
            )
        ilat = np.where((src_lat >= lat_min) & (src_lat <= lat_max))[0]
        ilon = np.where((src_lon >= lon_min) & (src_lon <= lon_max))[0]
        if ilat.size == 0 or ilon.size == 0:
            raise ValueError(
                f"STDEV_ELEV: no 1 km ELEVATION overlap with target land bbox "
                f"(lat [{lat_min:.3f},{lat_max:.3f}], "
                f"lon [{lon_min:.3f},{lon_max:.3f}])"
            )
        i0, i1 = int(ilat[0]), int(ilat[-1]) + 1
        j0, j1 = int(ilon[0]), int(ilon[-1]) + 1
        sub_lat = src_lat[i0:i1]
        sub_lon = src_lon[j0:j1]
        elev = (
            ds_elev["ELEVATION"]
            .isel(lat=slice(i0, i1), lon=slice(j0, j1))
            .values.astype(np.float64)
        )
        stdev = (
            ds_stdev["STDEV_ELEV"]
            .isel(lat=slice(i0, i1), lon=slice(j0, j1))
            .values.astype(np.float64)
        )
    finally:
        ds_elev.close()
        ds_stdev.close()

    # Flatten source pixels and drop invalid paired samples before binning.
    lon_grid, lat_grid = np.meshgrid(sub_lon, sub_lat)
    samples_elev = elev.ravel()
    samples_stdev = stdev.ravel()
    finite = (
        np.isfinite(samples_elev) & np.isfinite(samples_stdev) & (samples_stdev >= 0)
    )
    samples_lat = lat_grid.ravel()[finite]
    samples_lon = lon_grid.ravel()[finite]
    samples_elev = samples_elev[finite]
    samples_stdev_sq = np.square(samples_stdev[finite])
    if samples_elev.size == 0:
        return np.full(n_pts, np.nan, dtype=np.float64)

    # Law of total variance per target cell:
    #   total_var = mean(within_1km_variance) + variance(1km_mean_elevation)
    elev_std = binned_statistic_2d(
        samples_lat,
        samples_lon,
        samples_elev,
        statistic="std",
        bins=[lat_edges, lon_edges],
    ).statistic.astype(np.float64)
    within_var = binned_statistic_2d(
        samples_lat,
        samples_lon,
        samples_stdev_sq,
        statistic="mean",
        bins=[lat_edges, lon_edges],
    ).statistic.astype(np.float64)

    stat = np.sqrt(np.maximum(within_var + np.square(elev_std), 0.0))

    # Undo the flips so axes match LATIXY/LONGXY ordering.
    if lat_descending:
        stat = stat[::-1, :]
    if lon_descending:
        stat = stat[:, ::-1]

    return _fill_missing_land_with_edt(stat, land_2d, "STDEV_ELEV").ravel()


# -----------------------------------------------------------------------------
# Surfdata I/O: grid shape, land mask, land-only write-back, and the
# in-place SOIL_COLOR mode-filter smoothing pipeline.
# -----------------------------------------------------------------------------
def _lsmlat_lon_shape(ds: Dataset) -> tuple[int, int]:
    """2D land grid size (lsmlat, lsmlon), same order as LATIXY/LONGXY."""
    latixy = ds.variables["LATIXY"]
    if len(latixy.shape) != 2:
        raise ValueError(
            f"Expected LATIXY to be 2D (lsmlat, lsmlon); got shape {latixy.shape}."
        )
    return int(latixy.shape[0]), int(latixy.shape[1])


def _land_mask_2d(ds: Dataset) -> np.ndarray:
    """
    True where surfdata says ELM should run land physics (update soil there).

    Uses PFTDATA_MASK == 1 when present (0 = ocean / no land in this file).
    If absent, all cells are treated as land (with a runtime warning).
    """
    nlat, nlon = _lsmlat_lon_shape(ds)
    if "PFTDATA_MASK" in ds.variables:
        m = np.asarray(ds.variables["PFTDATA_MASK"][:], dtype=np.int32)
        if m.shape != (nlat, nlon):
            raise ValueError(f"PFTDATA_MASK shape {m.shape} != LATIXY {(nlat, nlon)}")
        return m == 1

    warnings.warn(
        "PFTDATA_MASK not found; updating every (lsmlat, lsmlon) cell including ocean.",
        UserWarning,
        stacklevel=2,
    )
    return np.ones((nlat, nlon), dtype=bool)


def _write_surf_var(
    ds: Dataset,
    var_name: str,
    data: np.ndarray,
    land_mask: np.ndarray,
) -> None:
    """
    Write interpolated columns in the same C-order as LATIXY.ravel():
    last index (lsmlon) varies fastest -> reshape to (lsmlat, lsmlon) or
    (nlevsoi, lsmlat, lsmlon).
    """
    v = ds.variables[var_name]
    arr = np.asarray(data)
    out_shape = v.shape
    nlat, nlon = _lsmlat_lon_shape(ds)
    if land_mask.shape != (nlat, nlon):
        raise ValueError(f"land_mask shape {land_mask.shape} != {(nlat, nlon)}")
    if out_shape[-2:] != (nlat, nlon):
        raise ValueError(
            f"Variable {var_name!r} spatial tail {out_shape[-2:]} "
            f"does not match LATIXY {(nlat, nlon)}."
        )

    old = np.asarray(v[:])

    if len(out_shape) == 2:
        flat = arr.ravel()
        if flat.size != nlat * nlon:
            raise ValueError(
                f"{var_name}: expected {nlat * nlon} points, got {flat.size}"
            )
        new = flat.reshape(nlat, nlon).astype(v.dtype)
        out = np.array(old, dtype=new.dtype, copy=True)
        out[land_mask] = new[land_mask]
        v[:, :] = out
    elif len(out_shape) == 3:
        nz = out_shape[0]
        if arr.ndim == 2:
            if arr.shape != (nz, nlat * nlon):
                raise ValueError(
                    f"{var_name}: expected shape ({nz}, {nlat * nlon}), got {arr.shape}"
                )
            arr = arr.reshape(nz, nlat, nlon)
        elif arr.ndim == 3:
            if arr.shape != out_shape:
                raise ValueError(
                    f"{var_name}: expected shape {out_shape}, got {arr.shape}"
                )
        new = arr.astype(v.dtype)
        out = np.array(old, dtype=new.dtype, copy=True)
        out[:, land_mask] = new[:, land_mask]
        v[:, :, :] = out
    else:
        flat = arr.ravel()
        if flat.size != int(np.prod(out_shape)):
            raise ValueError(
                f"{var_name}: cannot map length {flat.size} to shape {out_shape}"
            )
        new = flat.reshape(out_shape).astype(v.dtype)
        out = np.array(old, dtype=new.dtype, copy=True)
        lm = land_mask
        for _ in range(len(out_shape) - 2):
            lm = np.expand_dims(lm, 0)
        lm = np.broadcast_to(lm, out_shape)
        out[lm] = new[lm]
        v[:] = out


def ensure_soil_order_variable(
    ds: Dataset,
    data: np.ndarray,
    land_mask: np.ndarray,
    fill_value: float = -9999.0,
) -> None:
    """Write SOIL_ORDER; only land_mask cells get new values (see _land_mask_2d)."""
    nlat, nlon = _lsmlat_lon_shape(ds)
    if land_mask.shape != (nlat, nlon):
        raise ValueError(f"land_mask shape {land_mask.shape} != {(nlat, nlon)}")
    data_flat = np.asarray(data, dtype=np.float64).ravel()
    if data_flat.size != nlat * nlon:
        raise ValueError(
            f"SOIL_ORDER: expected {nlat * nlon} points, got {data_flat.size}"
        )
    data_i = np.rint(data_flat).astype(np.int32)
    data_i[~np.isfinite(data_flat)] = int(fill_value)
    grid2 = data_i.reshape(nlat, nlon)

    if "SOIL_ORDER" in ds.variables:
        old = np.asarray(ds.variables["SOIL_ORDER"][:], dtype=np.int32)
        out = old.copy()
        out[land_mask] = grid2[land_mask]
        ds.variables["SOIL_ORDER"][:, :] = out
        return

    template = ds.variables.get("SOIL_COLOR")
    if template is None:
        template = ds.variables.get("TOPO") or ds.variables.get("STD_ELEV")
    if template is None:
        raise KeyError(
            "Cannot create SOIL_ORDER: need SOIL_COLOR, TOPO, or STD_ELEV "
            "as dimension template."
        )
    dims = template.dimensions
    v = ds.createVariable(
        "SOIL_ORDER",
        np.int32,
        dims,
        zlib=True,
        complevel=4,
        fill_value=int(fill_value),
    )
    v.long_name = "soil order"
    v.units = "unitless"
    v.missing_value = int(fill_value)
    out = np.full((nlat, nlon), int(fill_value), dtype=np.int32)
    out[land_mask] = grid2[land_mask]
    v[:, :] = out


def _mode_filter_int(enc: np.ndarray, window: int) -> np.ndarray:
    """
    Vectorised per-pixel mode over a window-by-window neighbourhood for a
    2-D integer raster of small label cardinality.

    Values ``<= 0`` in ``enc`` are treated as *invalid* and never voted
    for; cells whose entire window contains only invalid values come out
    as 0 (caller is expected to backfill those, e.g. via EDT).

    Implementation: one-hot encode each class, sum each class plane with
    a uniform ``window x window`` kernel via :func:`scipy.ndimage.convolve`,
    then take ``argmax`` across classes. O(N * K) where K is the number of
    distinct classes (~20 for USDA soil colour), far faster than the
    Python-callback ``generic_filter`` route at continental sizes.
    """
    if window < 1 or window % 2 == 0:
        raise ValueError(f"window must be a positive odd integer, got {window}")
    enc = np.asarray(enc, dtype=np.int32)
    nlat, nlon = enc.shape
    max_val = int(enc.max(initial=0))
    if max_val <= 0:
        return np.zeros_like(enc, dtype=np.int32)
    kernel = np.ones((window, window), dtype=np.int32)
    counts = np.zeros((max_val, nlat, nlon), dtype=np.int32)
    for k in range(max_val):
        cls = k + 1  # classes 1..max_val
        is_cls = (enc == cls).astype(np.int32)
        if is_cls.any():
            counts[k] = convolve(is_cls, kernel, mode="reflect")
    smoothed = counts.argmax(axis=0) + 1  # back to 1..max_val
    smoothed[counts.sum(axis=0) == 0] = 0  # window had no valid neighbours
    return smoothed.astype(np.int32)


def smooth_soil_color_in_place(
    ds: Dataset,
    land_mask: np.ndarray,
    *,
    window: int = 3,
    var_name: str = "SOIL_COLOR",
) -> None:
    """
    Smooth the existing categorical SOIL_COLOR variable in-place using a
    ``window x window`` majority (mode) filter, restricted to land cells.

    Motivation: SOIL_COLOR in surfdata is often upsampled from a much
    coarser source by nearest-neighbour, leaving visible stair-step
    blocks. A mode filter over a small window smooths these boundaries
    without introducing fractional codes that aren't valid USDA colour
    classes (which is what bilinear would do).

    Steps:
        1. Read SOIL_COLOR; values ``<= 0`` are treated as invalid.
        2. ``_mode_filter_int`` computes the per-cell window mode, ignoring
           invalid neighbours.
        3. Commit the smoothed value only on ``land_mask`` cells; ocean
           cells keep whatever they had.
        4. Any land cell still without a value (entire window was invalid)
           is filled from its nearest land neighbour with a value (EDT).
    """
    if var_name not in ds.variables:
        raise KeyError(f"{var_name} not found in surfdata; nothing to smooth.")

    nlat, nlon = _lsmlat_lon_shape(ds)
    land_2d = np.asarray(land_mask, dtype=bool)
    if land_2d.shape != (nlat, nlon):
        raise ValueError(f"land_mask shape {land_2d.shape} != {(nlat, nlon)}")

    v = ds.variables[var_name]
    raw = np.ma.filled(np.ma.asarray(v[:]), 0).astype(np.int32)
    enc = np.where(raw > 0, raw, 0).astype(np.int32)

    smoothed = _mode_filter_int(enc, window=window)

    out = raw.copy()
    out[land_2d] = smoothed[land_2d]

    # EDT backfill: land cells with no valid neighbours in their window.
    have = land_2d & (out > 0)
    need = land_2d & ~have
    if need.any():
        if not have.any():
            raise RuntimeError(
                f"{var_name}: no land cell has a valid value after mode "
                "filtering; check land_mask vs raw SOIL_COLOR coverage."
            )
        _, (rr, cc) = distance_transform_edt(~have, return_indices=True)
        out[need] = out[rr[need], cc[need]]

    v[:, :] = out.astype(v.dtype)


# -----------------------------------------------------------------------------
# Diagnostic plots (skip when matplotlib is not available; cartopy NOT used)
# -----------------------------------------------------------------------------
def check_plot(
    varname,
    data_new,
    lat_new,
    lon_new,
    *,
    do_plot: bool = True,
    path_global_cf: str | None = None,
    path_soil_order_tif: str | None = None,
    grid_shape: tuple[int, int] | None = None,
) -> None:
    """
    Verify interpolation by plotting original 1 km source vs new surfgrid.

    Plain lon/lat matplotlib axes -- no cartopy projection / coastlines so
    plots are fast even for the 10-layer ORGANIC / PCT_CLAY / PCT_SAND
    variables. Saved as ``check_<varname>.png`` at dpi=300 under
    ``SOIL_CHECK_PLOT_DIR``.
    """
    if not do_plot:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"Skipping check plot for {varname}: missing matplotlib ({e})")
        return

    print(f"Plotting verification for {varname}...", flush=True)
    os.makedirs(_check_plot_dir(), exist_ok=True)

    if data_new.ndim == 1:
        data_new = data_new.reshape(1, -1)

    nlayer = data_new.shape[0]

    vmin = float(np.nanmin(data_new))
    vmax = float(np.nanmax(data_new))

    fig, ax = plt.subplots(
        nlayer,
        2,
        figsize=(12, 4 * nlayer),
        squeeze=False,
    )

    for n in range(nlayer):
        val_new = data_new[n, :]
        use_raster = (
            grid_shape is not None
            and len(grid_shape) == 2
            and int(np.prod(grid_shape)) == lon_new.size
        )
        if use_raster:
            lon2d = np.asarray(lon_new).reshape(grid_shape)
            lat2d = np.asarray(lat_new).reshape(grid_shape)
            val2d = np.asarray(val_new).reshape(grid_shape)
            sc = ax[n, 1].pcolormesh(
                lon2d,
                lat2d,
                val2d,
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
                shading="auto",
            )
        else:
            sc = ax[n, 1].scatter(
                lon_new,
                lat_new,
                c=val_new,
                s=1,
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )
        ax[n, 1].set_title(f"New surfgrid (Layer {n + 1})")
        ax[n, 1].set_xlabel("Longitude")
        ax[n, 1].set_ylabel("Latitude")
        ax[n, 1].set_aspect("equal", adjustable="box")
        plt.colorbar(sc, ax=ax[n, 1], fraction=0.046, pad=0.04)

    try:
        path_src_1km = (
            path_global_cf if path_global_cf is not None else _path_global_cf()
        )

        if varname in SOIL_NC_VARS or varname == "STDEV_ELEV":
            fname = _nc_filename_for_var(varname)
            fpath = os.path.join(path_src_1km, fname)

            # IO-defensive: integer .isel() BEFORE any astype / mask, same as
            # `interpolate_1km_nc_to_targets` and `coarsen_stdev_elev_from_1km_stats`.
            # clip_box on 3D *_10layer NCs (6 GB) sometimes doesn't push down
            # to the backend, then the full var gets promoted to float64 and
            # the OS SIGKILLs us (exit 137).
            pad = 0.25
            lat_min = float(np.nanmin(lat_new)) - pad
            lat_max = float(np.nanmax(lat_new)) + pad
            lon_min = float(np.nanmin(lon_new)) - pad
            lon_max = float(np.nanmax(lon_new)) + pad

            ds_src = xr.open_dataset(fpath, mask_and_scale=True, decode_times=False)
            try:
                src_lat = np.asarray(ds_src["lat"].values, dtype=np.float64)
                src_lon = np.asarray(ds_src["lon"].values, dtype=np.float64)
                ilat = np.where((src_lat >= lat_min) & (src_lat <= lat_max))[0]
                ilon = np.where((src_lon >= lon_min) & (src_lon <= lon_max))[0]
                if ilat.size == 0 or ilon.size == 0:
                    raise ValueError(
                        f"{varname}: target bbox has no source overlap "
                        f"(src lat [{src_lat.min():.3f},{src_lat.max():.3f}], "
                        f"lon [{src_lon.min():.3f},{src_lon.max():.3f}])"
                    )
                i0, i1 = int(ilat[0]), int(ilat[-1]) + 1
                j0, j1 = int(ilon[0]), int(ilon[-1]) + 1
                y_coord = src_lat[i0:i1]
                x_coord = src_lon[j0:j1]
                var_da = ds_src[varname].isel(lat=slice(i0, i1), lon=slice(j0, j1))
                print(
                    f"  reading 1km slab for {varname} "
                    f"(shape={tuple(var_da.shape)}, lat[{i0}:{i1}] lon[{j0}:{j1}])...",
                    flush=True,
                )
                arr_src = np.asarray(var_da.values, dtype=np.float64)
                non_spatial = [d for d in var_da.dims if d not in ("lat", "lon")]
                if non_spatial:
                    axis = list(var_da.dims).index(non_spatial[0])
                    if axis != 0:
                        arr_src = np.moveaxis(arr_src, axis, 0)
            finally:
                ds_src.close()

            extent = [
                float(x_coord.min()),
                float(x_coord.max()),
                float(y_coord.min()),
                float(y_coord.max()),
            ]
            # Source axis direction is preserved by .isel(); pick the right
            # matplotlib origin so the image isn't flipped vertically.
            origin = (
                "upper" if (y_coord.size > 1 and y_coord[0] > y_coord[-1]) else "lower"
            )

            for n in range(nlayer):
                if arr_src.ndim == 3:
                    if n >= arr_src.shape[0]:
                        continue
                    data_orig = arr_src[n]
                else:
                    if n != 0:
                        # 2D source: only the first row gets the "original" panel;
                        # leave the rest blank rather than re-plotting the same field.
                        continue
                    data_orig = arr_src

                im = ax[n, 0].imshow(
                    data_orig,
                    extent=extent,
                    origin=origin,
                    cmap="viridis",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax[n, 0].set_title(f"Original 1km ({fname}) Layer {n + 1}")
                ax[n, 0].set_xlabel("Longitude")
                ax[n, 0].set_ylabel("Latitude")
                ax[n, 0].set_aspect("equal", adjustable="box")
                plt.colorbar(im, ax=ax[n, 0], fraction=0.046, pad=0.04)

        elif varname == "SOIL_ORDER":
            fpath = (
                path_soil_order_tif
                if path_soil_order_tif is not None
                else _path_soil_order_tif()
            )
            if os.path.exists(fpath):
                with rio.open(fpath) as src:
                    window = from_bounds(
                        lon_new.min(),
                        lat_new.min(),
                        lon_new.max(),
                        lat_new.max(),
                        src.transform,
                    )
                    data_orig = src.read(1, window=window, masked=True)
                    bounds = src.window_bounds(window)
                    extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
                    for n in range(nlayer):
                        im = ax[n, 0].imshow(
                            data_orig,
                            extent=extent,
                            origin="upper",
                            cmap="viridis",
                            vmin=vmin,
                            vmax=vmax,
                        )
                        ax[n, 0].set_title("Original 1km (SOIL_ORDER)")
                        ax[n, 0].set_xlabel("Longitude")
                        ax[n, 0].set_ylabel("Latitude")
                        ax[n, 0].set_aspect("equal", adjustable="box")
                        plt.colorbar(im, ax=ax[n, 0], fraction=0.046, pad=0.04)

    except Exception as e:
        print(f"Could not plot original data for {varname}: {e}")

    out_png = os.path.join(_check_plot_dir(), f"check_{varname}.png")
    print(f"  saving {out_png}", flush=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def check_plot_soil_color_smoothing(
    soil_color_before: np.ndarray,
    soil_color_after: np.ndarray,
    lat_new: np.ndarray,
    lon_new: np.ndarray,
    *,
    do_plot: bool = True,
    grid_shape: tuple[int, int] | None = None,
    varname: str = "SOIL_COLOR",
) -> None:
    """
    Three-panel diagnostic for the SOIL_COLOR in-place mode-filter smoothing:

        | before | after | cells that changed |

    The first two panels share a colour scale (1..max valid USDA class);
    the third panel highlights every cell whose value changed and reports
    the count in its title. Saved as ``check_<varname>_smoothing.png`` in
    ``SOIL_CHECK_PLOT_DIR`` at dpi=300. Plain lon/lat matplotlib axes --
    no cartopy projection. No-op when ``do_plot=False`` or matplotlib is
    not importable.
    """
    if not do_plot:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"Skipping smoothing plot for {varname}: missing matplotlib ({e})")
        return

    print(f"Plotting before/after for {varname} smoothing...")
    os.makedirs(_check_plot_dir(), exist_ok=True)

    before = np.asarray(soil_color_before, dtype=np.float64)
    after = np.asarray(soil_color_after, dtype=np.float64)
    if before.shape != after.shape:
        raise ValueError(
            f"{varname}: before/after shape mismatch {before.shape} vs {after.shape}"
        )

    # Mask invalid (0 / fill) for display; arithmetic comparison below uses raw values.
    before_disp = np.where(before > 0, before, np.nan)
    after_disp = np.where(after > 0, after, np.nan)
    finite_before = np.isfinite(before_disp)
    finite_after = np.isfinite(after_disp)
    if not (finite_before.any() or finite_after.any()):
        print(f"{varname}: no valid cells to plot.")
        return

    vmin = float(np.nanmin([np.nanmin(before_disp), np.nanmin(after_disp)]))
    vmax = float(np.nanmax([np.nanmax(before_disp), np.nanmax(after_disp)]))

    use_raster = (
        grid_shape is not None
        and len(grid_shape) == 2
        and int(np.prod(grid_shape)) == np.asarray(lon_new).size
    )
    if use_raster:
        lon2d = np.asarray(lon_new).reshape(grid_shape)
        lat2d = np.asarray(lat_new).reshape(grid_shape)
        before_2d = before_disp.reshape(grid_shape)
        after_2d = after_disp.reshape(grid_shape)
        # "changed" mask uses the raw arrays so 0==0 (ocean) is *not* highlighted.
        changed_2d = np.where(before != after, 1.0, np.nan).reshape(grid_shape)
    else:
        before_2d = before_disp
        after_2d = after_disp
        changed_2d = np.where(before != after, 1.0, np.nan)

    n_changed = int(np.sum(before != after))

    fig, ax = plt.subplots(1, 3, figsize=(18, 5))

    def _draw(ax_, data2d, title, *, cmap, vmin_=None, vmax_=None):
        if use_raster:
            sc = ax_.pcolormesh(
                lon2d,
                lat2d,
                data2d,
                cmap=cmap,
                vmin=vmin_,
                vmax=vmax_,
                shading="auto",
            )
        else:
            sc = ax_.scatter(
                lon_new,
                lat_new,
                c=np.asarray(data2d).ravel(),
                s=1,
                cmap=cmap,
                vmin=vmin_,
                vmax=vmax_,
            )
        ax_.set_title(title)
        ax_.set_xlabel("Longitude")
        ax_.set_ylabel("Latitude")
        ax_.set_aspect("equal", adjustable="box")
        plt.colorbar(sc, ax=ax_, fraction=0.046, pad=0.04)

    _draw(
        ax[0],
        before_2d,
        f"{varname} (before smoothing)",
        cmap="viridis",
        vmin_=vmin,
        vmax_=vmax,
    )
    _draw(
        ax[1],
        after_2d,
        f"{varname} (after smoothing)",
        cmap="viridis",
        vmin_=vmin,
        vmax_=vmax,
    )
    _draw(
        ax[2],
        changed_2d,
        f"Cells changed (n={n_changed})",
        cmap="Reds",
        vmin_=0.0,
        vmax_=1.0,
    )

    plt.savefig(
        os.path.join(_check_plot_dir(), f"check_{varname}_smoothing.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


# -----------------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------------
def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Regrid 1 km soil / topo NetCDFs and SOIL_ORDER GeoTIFF onto ELM surfdata "
            "(LATIXY/LONGXY), land-only write using PFTDATA_MASK."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--surfdata",
        type=str,
        default=None,
        metavar="PATH",
        help="Input surfdata NetCDF (else $SURFDATA_PATH or built-in default). "
        "This file is NOT modified; see --output.",
    )
    p.add_argument(
        "--output",
        type=str,
        default=None,
        metavar="PATH",
        help="Output surfdata NetCDF. Default: surfdata_Updatedsoil_<basename> "
        "in the same directory as the input. Pass the same path as --surfdata "
        "for in-place edits.",
    )
    p.add_argument(
        "--global-cf",
        type=str,
        default=None,
        dest="global_cf",
        metavar="DIR",
        help="Directory with 1 km *_1k_c230606.nc / *_10layer_1k_*.nc "
        "(else $SOIL_GLOBAL_CF_DIR or built-in default)",
    )
    p.add_argument(
        "--soil-order-tif",
        type=str,
        default=None,
        dest="soil_order_tif",
        metavar="PATH",
        help="Soil order GeoTIFF (else $SOIL_ORDER_TIF or built-in default)",
    )
    p.add_argument(
        "--soil-color-window",
        type=int,
        default=3,
        dest="soil_color_window",
        metavar="N",
        help="Side length (odd) of the SOIL_COLOR mode-filter smoothing window; "
        "smooths blocky upsampled SOIL_COLOR in-place on land cells.",
    )
    p.add_argument(
        "--smooth-soil-color",
        action=argparse.BooleanOptionalAction,
        default=True,
        dest="smooth_soil_color",
        help=(
            "Smooth the existing SOIL_COLOR field in-place using a mode filter "
            "on land cells. Use --no-smooth-soil-color to skip this step "
            "entirely (then --soil-color-window is ignored and SOIL_COLOR is "
            "left untouched)."
        ),
    )
    p.add_argument(
        "--plot",
        action="store_true",
        help="Write diagnostic PNGs (dpi=300, plain lon/lat axes) under "
        "SOIL_CHECK_PLOT_DIR (needs matplotlib only; cartopy is not used).",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    in_path = args.surfdata or os.environ.get("SURFDATA_PATH", DEFAULT_SURFDATA)
    out_path = args.output or _default_output_path(in_path)
    global_cf = args.global_cf or os.environ.get(
        "SOIL_GLOBAL_CF_DIR", DEFAULT_PATH_GLOBAL_CF
    )
    soil_order_tif = args.soil_order_tif or os.environ.get(
        "SOIL_ORDER_TIF", DEFAULT_PATH_SOIL_ORDER_TIF
    )
    soil_color_window = int(args.soil_color_window)
    smooth_soil_color = bool(args.smooth_soil_color)
    do_plot = bool(args.plot)

    if not os.path.isfile(in_path):
        raise FileNotFoundError(f"SURFDATA not found: {in_path}")

    in_abs = os.path.abspath(in_path)
    out_abs = os.path.abspath(out_path)
    if in_abs == out_abs:
        # Explicit in-place edit; warn so callers don't lose the original silently.
        print(f"WARNING: --output equals --surfdata; editing IN PLACE: {in_abs}")
    else:
        os.makedirs(os.path.dirname(out_abs) or ".", exist_ok=True)
        if os.path.isfile(out_abs):
            print(f"Overwriting existing output file: {out_abs}")
        print(f"Copying {in_abs}\n     -> {out_abs}")
        shutil.copyfile(in_abs, out_abs)

    nc_path = out_abs

    with Dataset(nc_path, mode="a", format="NETCDF4") as ds:
        grid_shape = tuple(int(v) for v in ds.variables["LATIXY"].shape)
        lat_tgt = np.asarray(ds.variables["LATIXY"][:]).ravel()
        lon_tgt = np.asarray(ds.variables["LONGXY"][:]).ravel()
        land_mask = _land_mask_2d(ds)

        soil_order = interpolate_soil_order_tif(
            lat_tgt,
            lon_tgt,
            grid_shape=grid_shape,
            land_mask=land_mask,
            tif_path=soil_order_tif,
        )

        for var in SOIL_NC_VARS:
            var_out = surfdata_var_name(var)
            data = interpolate_1km_nc_to_targets(
                var,
                lat_tgt,
                lon_tgt,
                grid_shape=grid_shape,
                land_mask=land_mask,
                path_src=global_cf,
            )
            check_plot(
                var,
                data,
                lat_tgt,
                lon_tgt,
                do_plot=do_plot,
                path_global_cf=global_cf,
                path_soil_order_tif=soil_order_tif,
                grid_shape=grid_shape,
            )
            _write_surf_var(ds, var_out, data, land_mask)
            ds.sync()

        check_plot(
            "SOIL_ORDER",
            soil_order.astype(float),
            lat_tgt,
            lon_tgt,
            do_plot=do_plot,
            path_global_cf=global_cf,
            path_soil_order_tif=soil_order_tif,
            grid_shape=grid_shape,
        )
        ensure_soil_order_variable(ds, soil_order, land_mask)
        ds.sync()

        # STDEV_ELEV: physically-correct sub-grid topographic heterogeneity
        # via the law of total variance:
        #   mean(STDEV_ELEV_1km**2) + var(ELEVATION_1km)
        # Do NOT directly average STDEV_ELEV (mean(std) != std(samples)).
        stdev_elev = coarsen_stdev_elev_from_1km_stats(
            lat_tgt,
            lon_tgt,
            grid_shape=grid_shape,
            land_mask=land_mask,
            path_src=global_cf,
        )
        check_plot(
            "STDEV_ELEV",
            stdev_elev,
            lat_tgt,
            lon_tgt,
            do_plot=do_plot,
            path_global_cf=global_cf,
            path_soil_order_tif=soil_order_tif,
            grid_shape=grid_shape,
        )
        _write_surf_var(ds, surfdata_var_name("STDEV_ELEV"), stdev_elev, land_mask)
        ds.sync()

        # SOIL_COLOR: smooth the existing (often coarse-upsampled, blocky)
        # field in-place with an N x N mode filter, restricted to land cells.
        # Snapshot before the in-place edit so we can plot before/after.
        if smooth_soil_color:
            soil_color_before = np.asarray(
                ds.variables["SOIL_COLOR"][:], dtype=np.float64
            ).copy()
            smooth_soil_color_in_place(ds, land_mask, window=soil_color_window)
            ds.sync()
            soil_color_after = np.asarray(
                ds.variables["SOIL_COLOR"][:], dtype=np.float64
            )
            check_plot_soil_color_smoothing(
                soil_color_before,
                soil_color_after,
                lat_tgt,
                lon_tgt,
                do_plot=do_plot,
                grid_shape=grid_shape,
            )
        else:
            print("Skipping SOIL_COLOR smoothing (--no-smooth-soil-color).")

    print("Done.")
    print(f"  Input  (unchanged): {in_abs}")
    print(f"  Output (updated):   {out_abs}")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
