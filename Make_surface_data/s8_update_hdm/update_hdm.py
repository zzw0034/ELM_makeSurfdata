"""
update_hdm.py

Bilinearly regrid the ``hdm`` (population density) field from a coarse elmforc
HDM NetCDF (e.g. 0.5 deg) onto the 1/24 deg SEUS grid defined by a surfdata
file, then fill any land cell that still has no valid value (NaN) from its
nearest land neighbour.

**EDT gap-fill** (Euclidean Distance Transform): for each land cell that is
still NaN after regridding, scipy finds the closest land cell that *does* have
a valid value and copies that value over. It is a nearest-neighbour fill in
2-D grid space (not bilinear).

The target grid, ``LATIXY`` / ``LONGXY``, ``EDGE*``, and ``LANDMASK`` all come
from the surfdata file (``PFTDATA_MASK == 1`` defines land). Ocean cells are
not gap-filled; only land cells are updated.

Pass ``--plot`` / ``--tif`` to write QA outputs for the first processed time
slice: a 3-panel PNG (original nearest-regrid | bilinear | abs diff) and
GeoTIFFs for QGIS.

Usage
-----
python update_hdm.py /path/to/elmforc_hdm_in.nc
python update_hdm.py /path/to/elmforc_hdm_in.nc --output /path/to/out.nc
python update_hdm.py /path/to/elmforc_hdm_in.nc --year 1850 --plot --tif
python update_hdm.py /path/to/elmforc_hdm_in.nc \\
    --surfdata /path/to/surfdata_SEUS_1_24deg_simyr1850.nc

Defaults
--------
--surfdata    s7_updataPdata/surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc
--output      elmforc_Updatedhdm_<basename> next to the input
--fill-gaps   on  (fill land NaNs after regrid; --no-fill-gaps to skip)
--plot        off (3-panel before/after/diff PNG for first time slice)
--tif         off (original + bilinear GeoTIFFs for QGIS)
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import xarray as xr
from netCDF4 import Dataset as NcDataset
from scipy.ndimage import distance_transform_edt

VAR_NAME = "hdm"
NODATA = -9999.0

# Target grid + land mask come from the s7 output surfdata under Make_surface_data.
DEFAULT_SURFDATA = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/"
    "s7_updataPdata/surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc"
)


def _default_output_path(input_path: str) -> str:
    """``elmforc_Updatedhdm_<basename>`` in the same directory as the input."""
    d, base = os.path.split(input_path)
    if base.startswith("elmforc."):
        new_base = "elmforc_Updatedhdm_" + base[len("elmforc.") :]
    else:
        stem, ext = os.path.splitext(base)
        new_base = f"{stem}_Updatedhdm{ext}"
    return os.path.join(d, new_base)


def _edges_from_grid(lat_1d: np.ndarray, lon_1d: np.ndarray) -> dict[str, float]:
    """Domain edge scalars from 1-D cell-center coordinates."""
    dlat = float(lat_1d[1] - lat_1d[0])
    dlon = float(lon_1d[1] - lon_1d[0])
    return {
        "EDGEW": float(lon_1d[0] - dlon / 2.0),
        "EDGEE": float(lon_1d[-1] + dlon / 2.0),
        "EDGES": float(lat_1d[0] - dlat / 2.0),
        "EDGEN": float(lat_1d[-1] + dlat / 2.0),
    }


def _load_target_grid_from_surfdata(
    surfdata_path: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """
    Read the 1/24 deg target grid and land mask from surfdata.

    Returns
    -------
    new_lat, new_lon, lat2d, lon2d, land_mask_bool, landmask_f32, edges
    """
    if not os.path.isfile(surfdata_path):
        raise FileNotFoundError(f"Surfdata not found: {surfdata_path}")

    sd = xr.open_dataset(surfdata_path, decode_times=False)
    for req in ("PFTDATA_MASK", "LATIXY", "LONGXY"):
        if req not in sd:
            raise KeyError(f"{req} not found in surfdata: {surfdata_path}")

    lat2d = np.asarray(sd["LATIXY"].values, dtype=np.float32)
    lon2d = np.asarray(sd["LONGXY"].values, dtype=np.float32)
    new_lat = lat2d[:, 0].astype(np.float32)
    new_lon = lon2d[0, :].astype(np.float32)

    pft_mask = np.asarray(sd["PFTDATA_MASK"].values, dtype=np.int32)
    land_mask = pft_mask == 1
    landmask_f32 = pft_mask.astype(np.float32)
    edges = _edges_from_grid(new_lat, new_lon)
    sd.close()

    print(
        f"Target grid from surfdata: lat={new_lat.size}, lon={new_lon.size}  "
        f"({edges['EDGEW']:.4f}..{edges['EDGEE']:.4f}, "
        f"{edges['EDGES']:.4f}..{edges['EDGEN']:.4f})"
    )
    print(f"Land cells: {int(land_mask.sum())} / {land_mask.size}")
    return new_lat, new_lon, lat2d, lon2d, land_mask, landmask_f32, edges


def fill_missing_hdm(
    field: np.ndarray,
    land_mask: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Fill land NaNs from the nearest land neighbour (EDT). Returns (field, n_filled)."""
    land_2d = np.asarray(land_mask, dtype=bool)
    raw = np.asarray(field, dtype=np.float64)
    have = land_2d & np.isfinite(raw)
    need = land_2d & ~have
    n_filled = int(need.sum())
    if n_filled == 0:
        return raw, 0
    if not have.any():
        raise RuntimeError(
            f"{VAR_NAME}: no land cell has a valid value to fill gaps from."
        )
    out = raw.copy()
    _, (rr, cc) = distance_transform_edt(~have, return_indices=True)
    out[need] = out[rr[need], cc[need]]
    return out, n_filled


def _regrid_hdm_bilinear(
    da: xr.DataArray,
    new_lat: np.ndarray,
    new_lon: np.ndarray,
) -> np.ndarray:
    """
    Bilinearly regrid one ``hdm`` time slice onto the target grid.

    ``fill_value=None`` extrapolates the thin edge band beyond the source cell
    centers so no NaN gaps appear at the domain boundary (same as
    ``regrid_hdm_1_24deg.py``).
    """
    out = da.interp(
        lat=new_lat,
        lon=new_lon,
        method="linear",
        kwargs={"fill_value": None},
    )
    return out.clip(min=0.0).astype(np.float32).values


def _regrid_hdm_nearest(
    da: xr.DataArray,
    new_lat: np.ndarray,
    new_lon: np.ndarray,
) -> np.ndarray:
    """Nearest-neighbour regrid for QA comparison (blocky coarse on fine grid)."""
    out = da.interp(lat=new_lat, lon=new_lon, method="nearest")
    return out.clip(min=0.0).astype(np.float32).values


def _write_check_plot(
    orig: np.ndarray,
    bilinear: np.ndarray,
    lat2d: np.ndarray,
    lon2d: np.ndarray,
    land_mask: np.ndarray,
    *,
    year_label: str,
    out_png: str,
) -> None:
    """3-panel PNG: original (nearest) | bilinear | abs diff on land."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    land = np.asarray(land_mask, dtype=bool)
    vmin = float(np.nanmin([np.nanmin(orig[land]), np.nanmin(bilinear[land])]))
    vmax = float(np.nanmax([np.nanmax(orig[land]), np.nanmax(bilinear[land])]))
    diff = np.where(land, np.abs(orig.astype(np.float64) - bilinear.astype(np.float64)), np.nan)
    mean_abs_diff = float(np.nanmean(diff))
    max_abs_diff = float(np.nanmax(diff))

    fig, axes = plt.subplots(1, 3, figsize=(21, 7))
    for ax, field, title in [
        (axes[0], orig, f"{VAR_NAME} original (nearest, {year_label})"),
        (axes[1], bilinear, f"{VAR_NAME} bilinear ({year_label})"),
    ]:
        im = ax.pcolormesh(
            lon2d, lat2d, field, cmap="viridis", vmin=vmin, vmax=vmax, shading="auto"
        )
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    im2 = axes[2].pcolormesh(lon2d, lat2d, diff, cmap="magma", shading="auto")
    axes[2].set_title(
        f"abs diff (land only)\nmean={mean_abs_diff:.4g}, max={max_abs_diff:.4g}"
    )
    axes[2].set_xlabel("Longitude")
    axes[2].set_ylabel("Latitude")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(
        f"[{VAR_NAME}] plot ({year_label}): mean abs diff={mean_abs_diff:.4g}, "
        f"max abs diff={max_abs_diff:.4g} -> {out_png}"
    )


def _write_seus_geotiff(
    arr: np.ndarray,
    lat2d: np.ndarray,
    lon2d: np.ndarray,
    out_path: str,
) -> None:
    """Write a (lat, lon) array as a north-up GeoTIFF (EPSG:4326)."""
    from rasterio.crs import CRS
    from rasterio.transform import from_origin
    import rasterio

    lat1d = np.asarray(lat2d)[:, 0]
    lon1d = np.asarray(lon2d)[0, :]
    dlat = float(np.mean(np.diff(lat1d)))
    dlon = float(np.mean(np.diff(lon1d)))
    west_edge = float(lon1d[0]) - dlon / 2.0
    north_edge = float(lat1d[-1]) + dlat / 2.0
    transform = from_origin(west_edge, north_edge, dlon, dlat)

    arr2d = np.asarray(arr, dtype=np.float32)
    arr2d = np.where(np.isnan(arr2d), NODATA, arr2d).astype(np.float32)
    arr2d = np.flipud(arr2d)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with rasterio.open(
        out_path,
        "w",
        driver="GTiff",
        height=arr2d.shape[0],
        width=arr2d.shape[1],
        count=1,
        dtype="float32",
        crs=CRS.from_epsg(4326),
        transform=transform,
        compress="deflate",
        tiled=True,
        nodata=NODATA,
    ) as dst:
        dst.write(arr2d, 1)


def _prepare_geotiff_array(arr: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """Mask ocean cells to NaN for cleaner QGIS overlays."""
    out = np.asarray(arr, dtype=np.float64)
    out[~np.asarray(land_mask, dtype=bool)] = np.nan
    return out


def _write_check_geotiffs(
    orig: np.ndarray,
    bilinear: np.ndarray,
    lat2d: np.ndarray,
    lon2d: np.ndarray,
    land_mask: np.ndarray,
    *,
    year_label: str,
    tif_dir: str,
) -> None:
    """Write original (nearest) and bilinear GeoTIFFs for QGIS."""
    stem = f"check_{VAR_NAME}_{year_label}"
    out_orig = os.path.join(tif_dir, f"{stem}_original.tif")
    out_bil = os.path.join(tif_dir, f"{stem}_bilinear.tif")
    _write_seus_geotiff(_prepare_geotiff_array(orig, land_mask), lat2d, lon2d, out_orig)
    _write_seus_geotiff(_prepare_geotiff_array(bilinear, land_mask), lat2d, lon2d, out_bil)
    print(f"[{VAR_NAME}] GeoTIFF original: {out_orig}")
    print(f"[{VAR_NAME}] GeoTIFF bilinear: {out_bil}")


def _build_output_slice_dataset(
    ds: xr.Dataset,
    ti: int,
    hdm_slab: np.ndarray,
    *,
    new_lat: np.ndarray,
    new_lon: np.ndarray,
    edges: dict[str, float],
    landmask: np.ndarray,
    lat2d: np.ndarray,
    lon2d: np.ndarray,
) -> xr.Dataset:
    """One-time-slice output dataset on the surfdata grid."""
    data_vars: dict[str, xr.DataArray] = {
        VAR_NAME: xr.DataArray(
            hdm_slab[np.newaxis, ...],
            dims=("time", "lat", "lon"),
            attrs=dict(ds[VAR_NAME].attrs),
        ),
        "year": ds["year"].isel(time=ti).expand_dims("time"),
        "time_bnds": ds["time_bnds"].isel(time=ti).expand_dims("time"),
    }
    for edge_name, edge_val in edges.items():
        dim = ds[edge_name].dims if edge_name in ds else ("scalar",)
        attrs = dict(ds[edge_name].attrs) if edge_name in ds else {}
        data_vars[edge_name] = xr.DataArray(
            np.array([edge_val], dtype="float32"),
            dims=dim,
            attrs=attrs,
        )
    landmask_attrs = dict(ds["LANDMASK"].attrs) if "LANDMASK" in ds else {}
    data_vars["LANDMASK"] = xr.DataArray(
        landmask,
        dims=("lat", "lon"),
        attrs=landmask_attrs,
    )
    longxy_attrs = dict(ds["LONGXY"].attrs) if "LONGXY" in ds else {}
    latixy_attrs = dict(ds["LATIXY"].attrs) if "LATIXY" in ds else {}
    data_vars["LONGXY"] = xr.DataArray(lon2d, dims=("lat", "lon"), attrs=longxy_attrs)
    data_vars["LATIXY"] = xr.DataArray(lat2d, dims=("lat", "lon"), attrs=latixy_attrs)

    lat_attrs = dict(ds["lat"].attrs) if "lat" in ds else {}
    lon_attrs = dict(ds["lon"].attrs) if "lon" in ds else {}
    return xr.Dataset(
        data_vars=data_vars,
        coords={
            "lat": ("lat", new_lat, lat_attrs),
            "lon": ("lon", new_lon, lon_attrs),
            "time": ds["time"].isel(time=ti).expand_dims("time"),
        },
        attrs=dict(ds.attrs),
    )


def _append_time_slice(out_path: str, slice_ds: xr.Dataset, time_index: int) -> None:
    """Append one time slice to an existing elmforc output file."""
    with NcDataset(out_path, "a") as nc:
        nc.variables[VAR_NAME][time_index, :, :] = np.ascontiguousarray(
            slice_ds[VAR_NAME].values[0]
        )
        nc.variables["year"][time_index] = int(slice_ds["year"].values[0])
        nc.variables["time"][time_index] = float(slice_ds["time"].values[0])
        nc.variables["time_bnds"][time_index, :] = np.ascontiguousarray(
            slice_ds["time_bnds"].values[0]
        )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            f"Bilinearly regrid {VAR_NAME} from a coarse elmforc HDM file onto "
            "the surfdata 1/24 deg grid; gap-fill land NaNs (EDT). "
            "Writes a new output file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "input",
        type=str,
        metavar="HDM_NC",
        help="Input elmforc HDM NetCDF (not modified).",
    )
    p.add_argument(
        "--output",
        type=str,
        default=None,
        metavar="PATH",
        help="Output NetCDF. Default: elmforc_Updatedhdm_<basename>.",
    )
    p.add_argument(
        "--surfdata",
        type=str,
        default=DEFAULT_SURFDATA,
        metavar="PATH",
        help="Surfdata NetCDF that defines the target grid and land mask "
        "(PFTDATA_MASK, LATIXY, LONGXY).",
    )
    p.add_argument(
        "--year",
        type=int,
        default=None,
        metavar="YEAR",
        help="Process only this calendar year (default: all time slices).",
    )
    p.add_argument(
        "--fill-gaps",
        action=argparse.BooleanOptionalAction,
        default=True,
        dest="fill_gaps",
        help=f"After regrid, fill land cells with no valid {VAR_NAME} (NaN).",
    )
    p.add_argument(
        "--plot",
        action="store_true",
        help="Write a 3-panel check PNG (nearest original | bilinear | abs diff) "
        "for the first processed time slice.",
    )
    p.add_argument(
        "--plot-dir",
        type=str,
        default=None,
        dest="plot_dir",
        metavar="DIR",
        help="[--plot] Directory for check PNG. Default: check_plots/ next to --output.",
    )
    p.add_argument(
        "--tif",
        action="store_true",
        help="Write check GeoTIFFs (original nearest + bilinear) for QGIS.",
    )
    p.add_argument(
        "--tif-dir",
        type=str,
        default=None,
        dest="tif_dir",
        metavar="DIR",
        help="[--tif] Directory for GeoTIFFs. Default: <plot-dir>/tif.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    in_path = args.input
    out_path = args.output or _default_output_path(in_path)
    fill_gaps = bool(args.fill_gaps)

    if not os.path.isfile(in_path):
        raise FileNotFoundError(f"Input HDM file not found: {in_path}")

    in_abs = os.path.abspath(in_path)
    out_abs = os.path.abspath(out_path)
    surfdata_abs = os.path.abspath(args.surfdata)
    if in_abs == out_abs:
        raise ValueError("--output must differ from the input path.")

    os.makedirs(os.path.dirname(out_abs) or ".", exist_ok=True)
    if os.path.isfile(out_abs):
        print(f"Overwriting existing output file: {out_abs}")

    new_lat, new_lon, lat2d, lon2d, land_mask, landmask, edges = (
        _load_target_grid_from_surfdata(surfdata_abs)
    )

    ds = xr.open_dataset(in_abs, decode_times=False)
    if VAR_NAME not in ds:
        raise KeyError(f"{VAR_NAME} not found in {in_abs}")

    if args.year is not None:
        tidx = np.where(ds["year"].values == args.year)[0]
        if tidx.size == 0:
            raise ValueError(f"year {args.year} not found in file.")
        time_indices = [int(tidx[0])]
        print(f"Processing year {args.year} (time index {time_indices[0]}).")
    else:
        time_indices = list(range(ds.sizes["time"]))
        print(f"Processing all {len(time_indices)} time slice(s).")

    encoding = {
        v: {"zlib": True, "complevel": 4}
        for v in [VAR_NAME, "LANDMASK", "LATIXY", "LONGXY"]
    }
    total_gap_filled = 0
    qa_done = False
    plot_dir = args.plot_dir or os.path.join(os.path.dirname(out_abs) or ".", "check_plots")
    tif_dir = args.tif_dir or os.path.join(plot_dir, "tif")

    for n_done, ti in enumerate(time_indices, start=1):
        hdm_da = ds[VAR_NAME].isel(time=ti)
        slab_bilinear = _regrid_hdm_bilinear(hdm_da, new_lat, new_lon)

        if fill_gaps:
            slab_bilinear, n_filled = fill_missing_hdm(slab_bilinear, land_mask)
            total_gap_filled += n_filled

        if (args.plot or args.tif) and not qa_done:
            year_label = str(int(ds["year"].values[ti]))
            slab_nearest = _regrid_hdm_nearest(hdm_da, new_lat, new_lon)
            if args.plot:
                out_png = os.path.join(plot_dir, f"check_{VAR_NAME}_{year_label}.png")
                _write_check_plot(
                    slab_nearest,
                    slab_bilinear,
                    lat2d,
                    lon2d,
                    land_mask,
                    year_label=year_label,
                    out_png=out_png,
                )
            if args.tif:
                _write_check_geotiffs(
                    slab_nearest,
                    slab_bilinear,
                    lat2d,
                    lon2d,
                    land_mask,
                    year_label=year_label,
                    tif_dir=tif_dir,
                )
            qa_done = True

        slice_ds = _build_output_slice_dataset(
            ds,
            ti,
            slab_bilinear,
            new_lat=new_lat,
            new_lon=new_lon,
            edges=edges,
            landmask=landmask,
            lat2d=lat2d,
            lon2d=lon2d,
        )
        if n_done == 1:
            hist = slice_ds.attrs.get("history", "")
            slice_ds.attrs["history"] = (
                f"Bilinearly regridded {VAR_NAME} onto surfdata 1/24 deg grid "
                f"(land-mask gap-fill via update_hdm.py); " + hist
            )
            slice_ds.to_netcdf(
                out_abs,
                engine="netcdf4",
                encoding=encoding,
                unlimited_dims=["time"],
            )
        else:
            _append_time_slice(out_abs, slice_ds, n_done - 1)

        if n_done % 25 == 0 or n_done == len(time_indices):
            print(f"  wrote time index {ti} ({n_done}/{len(time_indices)})")

    ds.close()

    if fill_gaps:
        print(f"Gap-filled {total_gap_filled} land cell(s) across all time slices.")
    else:
        print("Skipped gap fill (--no-fill-gaps).")
    print("Done.")
    print(f"  Input    (unchanged): {in_abs}")
    print(f"  Surfdata (grid/mask): {surfdata_abs}")
    print(f"  Output   (updated):   {out_abs}")


if __name__ == "__main__":
    main()
