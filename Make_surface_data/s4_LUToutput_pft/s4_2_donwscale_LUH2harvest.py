"""Downscale LUH2 wood harvest on 0.25 deg grid to 1/24 deg ELM PFT grid.

This script follows the workflow in
`/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/
AIprompt_s4_2_donwscale_LUH2harvest.md`
- Read LUH2 harvest variables on 0.25 deg grid.
- Use ELM PFT fractions to compute tree-based downscaling weights.
- Conserve harvested AREA inside every coarse LUH2 grid cell (area-weighted).
- Convert gridcell harvest fractions to LUT-style vegetated-unit fractions.
- Write yearly LUT-like NetCDF files (1850-2023).
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import xarray as xr


# -------------------------
# Paths and run parameters
# -------------------------
PFT_PATH = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
)
LUH2_PATH = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/luh/transitions.nc"
)
# Recent period (calendar year >= 2015): continue harvest with LUH2 v2f SSP2-4.5.
# Same 0.25 deg grid / lat-order / harvest variables as v2h; calendar 2015-2099.
# Scenarios are nearly identical before ~2030, so SSP2-4.5 stands in as near-reality.
LUH2_SSP_PATH = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/luh/ssp2rcp45_transitions.nc"
)
TEMPLATE_PATH = Path(
    "/projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/rawdata/LUT_LUH2_historical_04082019/LUT_LUH2_historical_2014_04082019.nc"
)
OUT_DIR = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest"
)

START_YEAR = 1850
END_YEAR = 2023
LUH2_YEAR0 = 850
LUH2_LAST_YEAR = 2014  # LUH2 v2h transitions file ends at 2014
LUH2_SSP_YEAR0 = 2015  # SSP2-4.5 v2f transitions starts at 2015
LUH2_SSP_LAST_YEAR = 2099
PFT_YEAR0 = 1850

# NetCDF missing-value sentinel (float32); valid data must not equal this.
FILL_VALUE = np.float32(-999.0)
TREE_PFT_IDXS = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int64)

# LUH2 -> LUT output naming
HARVEST_VAR_MAP = {
    "primf_harv": "HARVEST_VH1",
    "primn_harv": "HARVEST_VH2",
    "secmf_harv": "HARVEST_SH1",
    "secyf_harv": "HARVEST_SH2",
    "secnf_harv": "HARVEST_SH3",
}

LONG_NAMES = {
    "HARVEST_VH1": "harvest from primary forest",
    "HARVEST_VH2": "harvest from primary non-forest",
    "HARVEST_SH1": "harvest from secondary mature-forest",
    "HARVEST_SH2": "harvest from secondary young-forest",
    "HARVEST_SH3": "harvest from secondary non-forest",
    "GRAZING": "grazing of herbacous pfts",
}


def _safe_float32(arr: np.ndarray) -> np.ndarray:
    """Convert to float32 and replace NaN/inf with 0."""
    out = np.asarray(arr, dtype=np.float32)
    np.nan_to_num(out, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    return out


def _fraction_if_percent(arr: np.ndarray) -> np.ndarray:
    """Convert percent to fraction when values look like percent."""
    out = _safe_float32(arr)
    maxv = float(np.nanmax(out)) if out.size else 0.0
    if maxv > 1.5:
        out = out / 100.0
    return np.clip(out, 0.0, 1.0)


def _compute_grid_edges(
    lat_1d: np.ndarray, lon_1d: np.ndarray
) -> tuple[float, float, float, float]:
    """Return (EDGEN, EDGEE, EDGES, EDGEW)."""
    return (
        float(np.nanmax(lat_1d)),
        float(np.nanmax(lon_1d)),
        float(np.nanmin(lat_1d)),
        float(np.nanmin(lon_1d)),
    )


def _build_hr_to_coarse_index(
    lat_hr: np.ndarray,
    lon_hr: np.ndarray,
    lat_coarse: np.ndarray,
    lon_coarse: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map each high-res cell to a coarse cell id.

    Returns
    -------
    coarse_id_2d : (nlat_hr, nlon_hr) int64
        Flat coarse-cell id for each high-res cell.
    inside_2d : (nlat_hr, nlon_hr) bool
        True where high-res cell center falls in coarse grid extent.
    coarse_shape : (2,) int64
        [nlat_coarse, nlon_coarse]
    """
    nlat_c = lat_coarse.size
    nlon_c = lon_coarse.size

    dlat = float(np.median(np.diff(lat_coarse)))
    dlon = float(np.median(np.diff(lon_coarse)))

    lat0 = float(lat_coarse[0] - 0.5 * dlat)
    lon0 = float(lon_coarse[0] - 0.5 * dlon)

    lat_idx_1d = np.floor((lat_hr - lat0) / dlat).astype(np.int64)
    lon_idx_1d = np.floor((lon_hr - lon0) / dlon).astype(np.int64)

    lat_ok = (lat_idx_1d >= 0) & (lat_idx_1d < nlat_c)
    lon_ok = (lon_idx_1d >= 0) & (lon_idx_1d < nlon_c)

    lat_idx_2d = lat_idx_1d[:, None]
    lon_idx_2d = lon_idx_1d[None, :]
    inside_2d = lat_ok[:, None] & lon_ok[None, :]

    coarse_id_2d = lat_idx_2d * nlon_c + lon_idx_2d
    coarse_id_2d = coarse_id_2d.astype(np.int64)
    coarse_id_2d[~inside_2d] = 0

    coarse_shape = np.array([nlat_c, nlon_c], dtype=np.int64)
    return coarse_id_2d, inside_2d, coarse_shape


def _year_to_pft_index(year: int) -> int:
    idx = year - PFT_YEAR0
    if idx < 0:
        raise ValueError(f"PFT year index < 0 for year={year}.")
    return idx


def _year_to_luh2_index(year: int) -> int:
    y = min(year, LUH2_LAST_YEAR)
    idx = y - LUH2_YEAR0
    if idx < 0:
        raise ValueError(f"LUH2 year index < 0 for year={year}.")
    return idx


def _distribute_conservatively(
    coarse_2d: np.ndarray,
    weights_2d: np.ndarray,
    area_2d: np.ndarray,
    coarse_id_2d: np.ndarray,
    inside_2d: np.ndarray,
    ncoarse: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Distribute a coarse per-area fraction to high-res, conserving harvested AREA.

    The coarse value is a fraction of cell area (an intensity), so the quantity that
    must be conserved inside each coarse cell is the harvested AREA,
    ``sum(frac_i * area_i) = F_coarse * sum(area_i)``, not the sum of the fractions.
    With ``frac_i`` proportional to weight ``w_i`` this gives

        frac_i = F_coarse * w_i * sum(area) / sum(w * area)

    (The previous implementation used ``w_i / sum(w)``, i.e. it conserved the sum of
    fractions, which underestimated harvest by the coarse/fine cell-count ratio ~36x.)

    Returns
    -------
    hr_2d : (nlat_hr, nlon_hr)
        Area-conservatively distributed fraction on the high-res grid.
    coarse_area_sums : (ncoarse,)
        ``sum(hr_2d * area)`` within each coarse cell (placed harvested area, for
        diagnostics).
    """
    coarse_flat = _safe_float32(coarse_2d).reshape(-1)
    w_flat = _safe_float32(weights_2d).reshape(-1)
    a_flat = _safe_float32(area_2d).reshape(-1).astype(np.float64)
    id_flat = coarse_id_2d.reshape(-1)
    inside_flat = inside_2d.reshape(-1)

    w_flat = np.where(inside_flat, w_flat, 0.0)
    a_flat = np.where(inside_flat, a_flat, 0.0)
    sum_wa = np.bincount(id_flat, weights=w_flat * a_flat, minlength=ncoarse)
    sum_a = np.bincount(id_flat, weights=a_flat, minlength=ncoarse)

    hr_flat = np.zeros_like(w_flat, dtype=np.float32)
    denom = sum_wa[id_flat]
    valid = inside_flat & (denom > 0.0)
    hr_flat[valid] = (
        coarse_flat[id_flat[valid]]
        * sum_a[id_flat[valid]]
        * w_flat[valid]
        / denom[valid]
    ).astype(np.float32)

    coarse_area_sums = np.bincount(
        id_flat, weights=hr_flat.astype(np.float64) * a_flat, minlength=ncoarse
    )
    hr_2d = hr_flat.reshape(weights_2d.shape)
    return hr_2d, coarse_area_sums


def _build_output_dataset(
    year: int,
    lat_hr: np.ndarray,
    lon_hr: np.ndarray,
    pft_frac: np.ndarray,
    veg_frac_landmask: np.ndarray,
    harvest_vars: dict[str, np.ndarray],
    template_attrs: dict,
) -> xr.Dataset:
    """Create one LUT-style output dataset for one year."""
    lat2d, lon2d = np.meshgrid(lat_hr, lon_hr, indexing="ij")
    edgen, edgee, edges, edgew = _compute_grid_edges(lat_hr, lon_hr)

    ds = xr.Dataset(attrs=dict(template_attrs))
    ds.attrs["title"] = (
        f"Downscaled 0.25 deg LUH2 harvest on ELM 1/24 deg grid ({year})"
    )
    ds.attrs["creation_date"] = str(date.today())
    ds.attrs["source"] = "LUH2 transitions + ELM PFT fractions"

    ds["EDGEN"] = (
        (),
        edgen,
        {"long_name": "northern edge of surface grid", "units": "degrees north"},
    )
    ds["EDGEE"] = (
        (),
        edgee,
        {"long_name": "eastern edge of surface grid", "units": "degrees east"},
    )
    ds["EDGES"] = (
        (),
        edges,
        {"long_name": "southern edge of surface grid", "units": "degrees north"},
    )
    ds["EDGEW"] = (
        (),
        edgew,
        {"long_name": "western edge of surface grid", "units": "degrees east"},
    )

    ds["LAT"] = (
        ["lat"],
        lat_hr.astype(np.float32),
        {"long_name": "lat", "units": "degrees north"},
    )
    ds["LON"] = (
        ["lon"],
        lon_hr.astype(np.float32),
        {"long_name": "lon", "units": "degrees east"},
    )
    ds["LATIXY"] = (
        ["lat", "lon"],
        _safe_float32(lat2d),
        {"long_name": "latitude-2d", "units": "degrees north"},
    )
    ds["LONGXY"] = (
        ["lat", "lon"],
        _safe_float32(lon2d),
        {"long_name": "longitude-2d", "units": "degrees east"},
    )

    landmask = (veg_frac_landmask > 0.0).astype(np.float32)
    ds["LANDMASK"] = (
        ["lat", "lon"],
        landmask,
        {"long_name": "land mask", "units": "unitless", "_FillValue": FILL_VALUE},
    )

    ds["PCT_PFT"] = (
        ["pft", "lat", "lon"],
        _safe_float32(100.0 * pft_frac),
        {"long_name": "percent pft", "_FillValue": FILL_VALUE},
    )

    for name, arr in harvest_vars.items():
        ds[name] = (
            ["lat", "lon"],
            _safe_float32(arr),
            {"long_name": LONG_NAMES.get(name, name), "_FillValue": FILL_VALUE},
        )

    ds["GRAZING"] = (
        ["lat", "lon"],
        np.ones((lat_hr.size, lon_hr.size), dtype=np.float32),
        {"long_name": LONG_NAMES["GRAZING"], "_FillValue": FILL_VALUE},
    )

    ds = ds.assign_coords(
        lat=("lat", lat_hr.astype(np.float32)),
        lon=("lon", lon_hr.astype(np.float32)),
        pft=("pft", np.arange(pft_frac.shape[0], dtype=np.int32)),
    )
    ds["year"] = ((), np.int32(year))
    return ds


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Reading PFT: {PFT_PATH}")
    print(f"Reading LUH2 v2h    : {LUH2_PATH}")
    print(f"Reading LUH2 v2f SSP: {LUH2_SSP_PATH}")

    with xr.open_dataset(TEMPLATE_PATH, decode_times=False) as template_ds:
        template_attrs = dict(template_ds.attrs)

    with (
        xr.open_dataset(PFT_PATH, decode_times=False) as pft_ds,
        xr.open_dataset(LUH2_PATH, decode_times=False) as luh2_ds,
        xr.open_dataset(LUH2_SSP_PATH, decode_times=False) as ssp_ds,
    ):
        # SSP2-4.5 must share the v2h grid/lat-order so the coarse index and
        # lat sort built from v2h apply unchanged.
        if (ssp_ds["lat"].size, ssp_ds["lon"].size) != (
            luh2_ds["lat"].size,
            luh2_ds["lon"].size,
        ) or not np.allclose(ssp_ds["lat"].values, luh2_ds["lat"].values):
            raise RuntimeError("SSP2-4.5 grid does not match LUH2 v2h grid.")
        lat_hr = _safe_float32(pft_ds["lat"].values)
        lon_hr = _safe_float32(pft_ds["lon"].values)

        # Static 1850 vegetated fraction: used ONLY for the time-invariant LANDMASK.
        veg_frac_1850 = _fraction_if_percent(pft_ds["PCT_NATVEG"].isel(time=0).values)

        # High-res cell area (km^2) for area-conserving downscaling (depends on lat).
        _dlat_hr = float(np.median(np.diff(lat_hr)))
        _dlon_hr = float(np.median(np.diff(lon_hr)))
        _n = np.deg2rad(lat_hr + 0.5 * abs(_dlat_hr))
        _s = np.deg2rad(lat_hr - 0.5 * abs(_dlat_hr))
        _band = (6371.0**2) * np.deg2rad(abs(_dlon_hr)) * (np.sin(_n) - np.sin(_s))
        area_hr = np.repeat(_band[:, None], lon_hr.size, axis=1).astype(np.float64)

        # LUH2 latitude is descending in this dataset; sort to ascending for index mapping.
        lat_luh_raw = _safe_float32(luh2_ds["lat"].values)
        lon_luh = _safe_float32(luh2_ds["lon"].values)
        lat_sort_idx = np.argsort(lat_luh_raw)
        lat_luh = lat_luh_raw[lat_sort_idx]
        ncoarse = lat_luh.size * lon_luh.size

        coarse_id_2d, inside_2d, coarse_shape = _build_hr_to_coarse_index(
            lat_hr, lon_hr, lat_luh, lon_luh
        )
        if tuple(coarse_shape) != (lat_luh.size, lon_luh.size):
            raise RuntimeError("Unexpected coarse-shape mapping failure.")

        # Coarse (LUH2) cell area and the domain mask (coarse cells receiving at least
        # one high-res cell), for the area-conservation diagnostics in the log.
        _dlatc = float(np.median(np.diff(lat_luh)))
        _dlonc = float(np.median(np.diff(lon_luh)))
        _nc = np.deg2rad(lat_luh + 0.5 * abs(_dlatc))
        _sc = np.deg2rad(lat_luh - 0.5 * abs(_dlatc))
        _bandc = (6371.0**2) * np.deg2rad(abs(_dlonc)) * (np.sin(_nc) - np.sin(_sc))
        area_coarse = (
            np.repeat(_bandc[:, None], lon_luh.size, axis=1).astype(np.float64).reshape(-1)
        )
        domain_coarse = (
            np.bincount(
                coarse_id_2d.reshape(-1),
                weights=inside_2d.reshape(-1).astype(np.float64),
                minlength=ncoarse,
            )
            > 0
        )

        print("Beginning yearly downscaling...")
        for year in range(START_YEAR, END_YEAR + 1):
            pft_idx = _year_to_pft_index(year)
            # E3SM LUT convention (doc p3): the file labelled `year` carries the
            # harvest rate that occurred in the prior calendar year (cal = year-1).
            cal_year = year - 1
            # Source splice: v2h historical through cal 2014, then v2f SSP2-4.5.
            if cal_year <= LUH2_LAST_YEAR:
                harv_ds = luh2_ds
                harv_idx = cal_year - LUH2_YEAR0
                harv_src = f"v2h {cal_year}"
            else:
                capped = min(cal_year, LUH2_SSP_LAST_YEAR)
                harv_ds = ssp_ds
                harv_idx = capped - LUH2_SSP_YEAR0
                harv_src = f"ssp245 {capped}"

            # PFT fractions expected in %, convert to [0,1]
            pft_frac = _fraction_if_percent(
                pft_ds["PCT_NAT_PFT"].isel(time=pft_idx).values
            )
            # Annual vegetated fraction of the gridcell (previously frozen at 1850).
            veg_frac = _fraction_if_percent(
                pft_ds["PCT_NATVEG"].isel(time=pft_idx).values
            )
            tree_frac_veg = np.sum(
                pft_frac[TREE_PFT_IDXS, :, :], axis=0, dtype=np.float32
            )
            # Weight = forest area fraction of the whole gridcell; zero => no harvest.
            tree_frac_grid = _safe_float32(tree_frac_veg * veg_frac)
            weights_2d = np.where(tree_frac_grid > 0.0, tree_frac_grid, 0.0).astype(
                np.float32
            )

            # Coarse cells that can actually receive harvest this year (have forest).
            w_flat = np.where(inside_2d.reshape(-1), weights_2d.reshape(-1), 0.0)
            sum_w_coarse = np.bincount(
                coarse_id_2d.reshape(-1), weights=w_flat, minlength=ncoarse
            )
            alloc_coarse = domain_coarse & (sum_w_coarse > 0.0)

            # 1) Area-conservative downscaling of each harvest category (grid fraction).
            hr_grids: dict[str, np.ndarray] = {}
            luh_area_alloc = 0.0  # LUH2 harvested area in allocatable coarse cells
            hr_area_placed = 0.0  # area actually placed on the high-res grid
            for luh_name, out_name in HARVEST_VAR_MAP.items():
                coarse = _safe_float32(
                    harv_ds[luh_name].isel(time=harv_idx).values[lat_sort_idx, :]
                )  # reorder lat to ascending
                hr_grid, coarse_area_sums = _distribute_conservatively(
                    coarse_2d=coarse,
                    weights_2d=weights_2d,
                    area_2d=area_hr,
                    coarse_id_2d=coarse_id_2d,
                    inside_2d=inside_2d,
                    ncoarse=ncoarse,
                )
                hr_grids[out_name] = hr_grid
                cflat = coarse.reshape(-1).astype(np.float64)
                luh_area_alloc += float((cflat * area_coarse)[alloc_coarse].sum())
                hr_area_placed += float(coarse_area_sums[alloc_coarse].sum())

            # 2) Convert to LUT units (fraction of vegetated land unit) and cap the
            #    per-cell 5-category total to <= 1 (option A: clip, discard the excess),
            #    so ELM never sees > 100% annual tree removal.
            lut_grids = {
                out_name: np.divide(
                    hr,
                    veg_frac,
                    out=np.zeros_like(hr, dtype=np.float32),
                    where=veg_frac > 0.0,
                )
                for out_name, hr in hr_grids.items()
            }
            lut_sum = sum(lut_grids.values())
            over = lut_sum > 1.0
            scale = np.ones_like(lut_sum, dtype=np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                scale[over] = (1.0 / lut_sum[over]).astype(np.float32)
            clip_loss_area = float(
                (
                    np.where(over, (lut_sum - 1.0) * veg_frac, 0.0).reshape(-1)
                    * area_hr.reshape(-1)
                )[inside_2d.reshape(-1)].sum()
            )
            out_vars = {
                out_name: _safe_float32(np.clip(lut * scale, 0.0, 1.0))
                for out_name, lut in lut_grids.items()
            }

            # Build and write yearly LUT-style output
            out_ds = _build_output_dataset(
                year=year,
                lat_hr=lat_hr,
                lon_hr=lon_hr,
                pft_frac=pft_frac,
                veg_frac_landmask=veg_frac_1850,
                harvest_vars=out_vars,
                template_attrs=template_attrs,
            )

            out_path = (
                OUT_DIR
                / f"LUT_nlcd2elm_luh2_historical_{year}_{date.today().strftime('%m%d%Y')}.nc"
            )
            enc = {v: {"zlib": True, "complevel": 2} for v in out_ds.data_vars}
            out_ds.to_netcdf(out_path, encoding=enc)
            out_ds.close()

            # Validation summary: area-conservation ratio over allocatable coarse cells
            # (should be ~1 after the fix; was ~1/36 before) and the harvest fraction
            # discarded by the <=1 cap.
            hsum = sum(out_vars.values())
            cons = hr_area_placed / luh_area_alloc if luh_area_alloc > 0 else float("nan")
            clip_pct = (
                100.0 * clip_loss_area / hr_area_placed if hr_area_placed > 0 else 0.0
            )
            print(
                f"Year {year} ({harv_src}) wrote {out_path.name} | "
                f"sum_lut(max)={np.nanmax(hsum):.4g} | "
                f"area_placed/alloc={cons:.4f} | clip_loss={clip_pct:.3f}%"
            )

    print("All done.")


if __name__ == "__main__":
    main()
