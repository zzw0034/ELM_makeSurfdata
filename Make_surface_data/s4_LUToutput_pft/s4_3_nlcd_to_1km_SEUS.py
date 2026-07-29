#!/usr/bin/env python
"""s4_3 -- Annual NLCD -> 1 km ELM PFT over the SEUS domain, for one year.

Why this script exists
----------------------
`scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc` (step s4_1) is on the
1/24 deg (~4 km) grid, so it cannot produce a 1 km figure. The 1 km product is
built here straight from the Annual NLCD land cover that Dan reprojected to a
regular lat/lon grid:

    <NLCD_DIR>/nlcd_<year>_CONUS1.1800deg.nc      var `nlcd`, (band, y, x)
                                                  dy ~ 0.00068 deg, dx ~ 0.00061 deg

Pipeline (same idea as Dan's `regrid_nlcd.py`, one difference noted below):

  1. bin the ~60-70 m NLCD classes into 0.01 deg (~1 km) cells over the SEUS box,
     giving a percentage per NLCD class per cell;
  2. translate class percentages to ELM PFTs with the exact mapping used by
     `s4_1_convert_nlcd_to_PFT.py` (Dan's `nlcd_to_elmpft.py`, fc4 = 0.2);
  3. write a NetCDF and draw the figures.

Difference from `regrid_nlcd.py`: that script bins by an integer pixel count and
then does a nearest-neighbour `interp` onto the exact target grid ("_corr"
file), because target_res / pixel_size is not an integer. Here each source pixel
is instead assigned directly to the target cell that contains its centre, so the
binning lands on the exact 0.01 deg grid in one step with no nearest-neighbour
resampling and no truncated edge cells.

"1 km" means 0.01 deg here, matching the project's existing `0.01x0.01` grid
name (see the map_0.01x0.01_*_to_SEUS_1_24deg_* mapping files). 0.01 deg is
~1.11 km N-S and ~0.9 km E-W across the SEUS domain.

Default domain is the project SEUS grid box, i.e. the outer bounds of
`s1_Gridfiles/scr_out/SCRIPgrid_SEUS_1_24deg.nc` (504 x 324 cells at 1/24 deg):
lon [-95, -74], lat [24, 37.5].

Usage
-----
    python s4_3_nlcd_to_1km_SEUS.py --year 2023
    python s4_3_nlcd_to_1km_SEUS.py --year 2023 --plot-only   # re-draw figures
    python s4_3_nlcd_to_1km_SEUS.py --year 2015 --res 0.005   # finer target grid

Run it through `run_s4_3_nlcd_1km_SEUS.slurm`, not on a login node: it reads
~2.6 GB of NLCD per year.
"""

from __future__ import annotations

import argparse
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

NLCD_DIR = "/projects/hpcl-cli185/proj-shared/zdr/hires_data/landcover/NLCD"
NLCD_TEMPLATE = "nlcd_{year}_CONUS1.1800deg.nc"

OUT_DIR = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/"
    "s4_LUToutput_pft/scr_out/nlcd1km_SEUS"
)

# SEUS_1_24deg box (outer bounds of SCRIPgrid_SEUS_1_24deg.nc).
SEUS_BBOX = (-95.0, -74.0, 24.0, 37.5)  # lon_min, lon_max, lat_min, lat_max
TARGET_RES = 0.01  # deg; the project's "1 km" grid spacing

# NLCD classes kept. Anything else (0, 250, NaN, ...) counts as no-data.
NLCD_CODES = [11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95]

NLCD_LABELS = {
    11: "Open water",
    12: "Perennial ice/snow",
    21: "Developed, open space",
    22: "Developed, low intensity",
    23: "Developed, medium intensity",
    24: "Developed, high intensity",
    31: "Barren land",
    41: "Deciduous forest",
    42: "Evergreen forest",
    43: "Mixed forest",
    51: "Dwarf scrub",
    52: "Shrub/scrub",
    71: "Grassland/herbaceous",
    72: "Sedge/herbaceous",
    73: "Lichens",
    74: "Moss",
    81: "Pasture/hay",
    82: "Cultivated crops",
    90: "Woody wetlands",
    95: "Emergent herbaceous wetlands",
}

# Vegetated classes summed into PCT_NATVEG -- identical list to s4_1.
VEG_CODES = [21, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95]

FC4 = 0.2  # C4 fraction of the herbaceous classes, as in s4_1

ELM_PFT_NAMES = (
    "Bare_Ground",                          # 0
    "needleleaf_evergreen_temperate_tree",  # 1
    "needleleaf_evergreen_boreal_tree",     # 2
    "needleleaf_deciduous_boreal_tree",     # 3
    "broadleaf_evergreen_tropical_tree",    # 4
    "broadleaf_evergreen_temperate_tree",   # 5
    "broadleaf_deciduous_tropical_tree",    # 6
    "broadleaf_deciduous_temperate_tree",   # 7
    "broadleaf_deciduous_boreal_tree",      # 8
    "broadleaf_evergreen_shrub",            # 9
    "broadleaf_deciduous_temperate_shrub",  # 10
    "broadleaf_deciduous_boreal_shrub",     # 11
    "c3_arctic_grass",                      # 12
    "c3_non-arctic_grass",                  # 13
    "c4_grass",                             # 14
    "crop",                                 # 15
    "irrigated_crop",                       # 16
)
NPFT = len(ELM_PFT_NAMES)

# Aggregated NLCD classes for the dominant-land-cover figure (Dan's grouping in
# AI_nlcd/animate_landuse.py).
AGG_CLASSES = {
    "water": [11, 12],
    "urban": [21, 22, 23, 24],
    "barren": [31],
    "forest": [41, 42, 43],
    "shrub": [51, 52],
    "grass": [71, 72, 73, 74, 81],
    "crop": [82],
    "wetland": [90, 95],
}
AGG_COLORS = ["#1f6feb", "#c62828", "#8d6e63", "#1b5e20", "#ef8f00", "#9ccc65", "#fdd835", "#7e57c2"]


# --------------------------------------------------------------------------
# step 1 -- bin NLCD onto the 1 km target grid
# --------------------------------------------------------------------------


def build_target_grid(bbox, res):
    """Cell-centre coordinates of the regional target grid, lat ascending."""
    lon_min, lon_max, lat_min, lat_max = bbox
    nlon = int(round((lon_max - lon_min) / res))
    nlat = int(round((lat_max - lat_min) / res))
    lon = lon_min + (np.arange(nlon) + 0.5) * res
    lat = lat_min + (np.arange(nlat) + 0.5) * res
    return lon, lat


def bin_nlcd_to_grid(src_path, bbox, res, chunk_rows):
    """Percentage of each NLCD class within every target cell.

    Returns (lc_pct, valid_frac): lc_pct is (nclass, nlat, nlon) in percent of
    the *classified* source pixels in the cell; valid_frac is the fraction of
    source pixels in the cell that carried a real NLCD class.
    """
    lon_min, lon_max, lat_min, lat_max = bbox
    tlon, tlat = build_target_grid(bbox, res)
    nlon, nlat = tlon.size, tlat.size
    ncell = nlat * nlon
    nclass = len(NLCD_CODES)

    # value -> compact class index, -1 for no-data
    lut = np.full(256, -1, dtype=np.int16)
    for i, code in enumerate(NLCD_CODES):
        lut[code] = i

    ds = Dataset(src_path, "r")
    var = ds.variables["nlcd"]
    var.set_auto_maskandscale(False)
    ysrc = ds.variables["y"][:]  # descending
    xsrc = ds.variables["x"][:]  # ascending

    # Source pixel index range covering the box (y is descending -> reverse for
    # searchsorted, then map back).
    x0 = int(np.searchsorted(xsrc, lon_min, side="left"))
    x1 = int(np.searchsorted(xsrc, lon_max, side="right"))
    y_asc = ysrc[::-1]
    ya0 = int(np.searchsorted(y_asc, lat_min, side="left"))
    ya1 = int(np.searchsorted(y_asc, lat_max, side="right"))
    y0 = ysrc.size - ya1
    y1 = ysrc.size - ya0
    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"requested box {bbox} does not overlap {src_path}")
    print(f"  source window: y[{y0}:{y1}] ({y1 - y0} rows), x[{x0}:{x1}] ({x1 - x0} cols)")

    # Target column index of every source column (fixed for all chunks).
    ix = np.floor((xsrc[x0:x1] - lon_min) / res).astype(np.int64)
    np.clip(ix, 0, nlon - 1, out=ix)

    counts = np.zeros(nclass * ncell, dtype=np.int64)
    total = np.zeros(ncell, dtype=np.int64)

    t0 = time.time()
    for r0 in range(y0, y1, chunk_rows):
        r1 = min(r0 + chunk_rows, y1)
        block = np.asarray(var[0, r0:r1, x0:x1])

        iy = np.floor((ysrc[r0:r1] - lat_min) / res).astype(np.int64)
        np.clip(iy, 0, nlat - 1, out=iy)
        cell = iy[:, None] * nlon + ix[None, :]

        vals = np.nan_to_num(block, nan=0.0)
        np.clip(vals, 0, 255, out=vals)
        cls = lut[vals.astype(np.uint8)]

        ok = cls >= 0
        cell_ok = cell[ok]
        code = cls[ok].astype(np.int64) * ncell + cell_ok

        counts += np.bincount(code, minlength=nclass * ncell)
        total += np.bincount(cell.ravel(), minlength=ncell)

        done = (r1 - y0) / (y1 - y0)
        print(f"  rows {r0 - y0:>6d}-{r1 - y0:>6d}  {done * 100:5.1f}%  "
              f"{time.time() - t0:6.1f}s", flush=True)

    ds.close()

    counts = counts.reshape(nclass, nlat, nlon)
    total = total.reshape(nlat, nlon)
    classified = counts.sum(axis=0)

    lc_pct = np.zeros((nclass, nlat, nlon), dtype=np.float32)
    good = classified > 0
    lc_pct[:, good] = counts[:, good] / classified[good] * 100.0
    valid_frac = np.where(total > 0, classified / np.maximum(total, 1), 0.0).astype(np.float32)

    return tlon, tlat, lc_pct, valid_frac


# --------------------------------------------------------------------------
# step 2 -- NLCD class percentages -> ELM PFT
# --------------------------------------------------------------------------


def nlcd_pct_to_elmpft(lc_pct):
    """Exact mapping of s4_1_convert_nlcd_to_PFT.py, applied to one year."""
    idx = {code: i for i, code in enumerate(NLCD_CODES)}
    nlat, nlon = lc_pct.shape[1:]

    def lc(code):
        i = idx.get(code)
        return lc_pct[i] if i is not None else np.zeros((nlat, nlon), dtype=np.float32)

    pft = np.zeros((NPFT, nlat, nlon), dtype=np.float64)
    pft[0] = lc(31)                                         # bare ground
    pft[1] = lc(42) + (lc(43) + lc(90)) * 0.5               # NE temperate tree
    pft[7] = lc(41) + (lc(43) + lc(90)) * 0.5               # BD temperate tree
    pft[9] = (lc(52) + lc(51)) * 0.5                        # BE shrub
    pft[10] = (lc(52) + lc(51)) * 0.5                       # BD temperate shrub
    pft[13] = (1.0 - FC4) * (lc(71) + lc(72) + lc(81)) + lc(95)   # C3 grass
    pft[14] = FC4 * (lc(71) + lc(72) + lc(81))              # C4 grass
    pft[15] = lc(82)                                        # crop
    # developed open space treated as a 50/50 tree/grass mix
    pft[7] += lc(21) * 0.5
    pft[13] += lc(21) * 0.5

    sums = pft.sum(axis=0)
    mask = sums > 1e-6
    pft[:, mask] = pft[:, mask] / sums[mask] * 100.0
    pft[:, ~mask] = 0.0
    pft[0, ~mask] = 100.0  # only read where PCT_NATVEG > 0

    natveg = np.zeros((nlat, nlon), dtype=np.float64)
    for code in VEG_CODES:
        natveg += lc(code)

    urban = np.stack([lc(22), lc(23), lc(24)]).astype(np.float64)

    return pft.astype(np.float32), natveg.astype(np.float32), urban.astype(np.float32)


# --------------------------------------------------------------------------
# NetCDF I/O
# --------------------------------------------------------------------------


def write_netcdf(path, year, res, lon, lat, lc_pct, valid_frac, pft, natveg, urban):
    if os.path.exists(path):
        os.remove(path)
    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("lat", lat.size)
    nc.createDimension("lon", lon.size)
    nc.createDimension("natpft", NPFT)
    nc.createDimension("numurbl", 3)
    nc.createDimension("nlcdclass", len(NLCD_CODES))

    v = nc.createVariable("lat", "f8", ("lat",)); v.units = "degrees_north"; v[:] = lat
    v = nc.createVariable("lon", "f8", ("lon",)); v.units = "degrees_east"; v[:] = lon
    v = nc.createVariable("nlcd_class", "i4", ("nlcdclass",))
    v.long_name = "NLCD land cover class code"
    v[:] = np.array(NLCD_CODES, dtype=np.int32)

    def add(name, dims, data, units, long_name):
        var = nc.createVariable(name, "f4", dims, zlib=True, complevel=4)
        var.units = units
        var.long_name = long_name
        var[...] = data
        return var

    add("PCT_LC", ("nlcdclass", "lat", "lon"), lc_pct, "percent",
        "percentage of each NLCD class within the cell (of classified pixels)")
    add("PCT_NAT_PFT", ("natpft", "lat", "lon"), pft, "percent",
        "ELM natural PFT percentage, normalized to 100% of the vegetated fraction")
    add("PCT_NATVEG", ("lat", "lon"), natveg, "percent",
        "percentage of the cell that is naturally vegetated (sum of NLCD veg classes)")
    add("PCT_URBAN", ("numurbl", "lat", "lon"), urban, "percent",
        "NLCD developed low / medium / high intensity (classes 22, 23, 24)")
    add("LANDFRAC", ("lat", "lon"), valid_frac, "1",
        "fraction of source NLCD pixels in the cell carrying a valid class")

    nc.setncattr("title", f"NLCD {year} binned to {res} deg and mapped to ELM PFTs, SEUS domain")
    nc.setncattr("source", NLCD_TEMPLATE.format(year=year))
    nc.setncattr("pft_mapping", "identical to s4_1_convert_nlcd_to_PFT.py (Dan nlcd_to_elmpft.py), fc4=0.2")
    nc.setncattr("history", f"created by s4_3_nlcd_to_1km_SEUS.py on {time.strftime('%Y-%m-%d %H:%M:%S')}")
    nc.close()
    print(f"wrote {path}")


def read_netcdf(path):
    nc = Dataset(path, "r")
    out = tuple(
        np.asarray(nc.variables[name][:])
        for name in ("lon", "lat", "PCT_LC", "LANDFRAC", "PCT_NAT_PFT", "PCT_NATVEG", "PCT_URBAN")
    )
    nc.close()
    return out


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------


def _map_axes(ax, extent):
    """Plain lon/lat axes -- cartopy is not installed in make_surfdata_pf."""
    lon_min, lon_max, lat_min, lat_max = extent
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(0.5 * (lat_min + lat_max))))
    ax.tick_params(labelsize=7)


def fig_dominant_landcover(lc_pct, extent, year, out_png):
    idx = {code: i for i, code in enumerate(NLCD_CODES)}
    names = list(AGG_CLASSES)
    agg = np.stack([
        np.sum([lc_pct[idx[c]] for c in codes if c in idx], axis=0) for codes in AGG_CLASSES.values()
    ])
    dominant = np.argmax(agg, axis=0).astype(float)
    dominant[agg.sum(axis=0) <= 0] = np.nan

    cmap = ListedColormap(AGG_COLORS)
    cmap.set_bad("white")
    norm = BoundaryNorm(np.arange(-0.5, len(names)), cmap.N)

    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(dominant, origin="lower", extent=extent, cmap=cmap, norm=norm,
                   interpolation="nearest")
    _map_axes(ax, extent)
    ax.set_title(f"Dominant aggregated NLCD land cover, {year} — SEUS 0.01° (~1 km)")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    cb = fig.colorbar(im, ax=ax, ticks=range(len(names)), shrink=0.8)
    cb.ax.set_yticklabels(names)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def fig_pft_panels(pft, natveg, extent, year, out_png):
    """All 17 PFTs. Empty ones are kept so the panel index equals the PFT index.

    Cells with no vegetation are masked out. Without that, panel 0 renders the
    whole ocean at 100%: `nlcd_pct_to_elmpft` (like s4_1) parks a sentinel
    `PCT_NAT_PFT[0] = 100` wherever the PFT sum is zero, which is only
    meaningful because everything downstream gates on `PCT_NATVEG > 0`.
    """
    veg = natveg > 0
    fig, axes = plt.subplots(5, 4, figsize=(17, 17), constrained_layout=True)
    for k, ax in enumerate(axes.ravel()):
        if k >= NPFT:
            ax.axis("off")
            continue
        field = np.where(veg & (pft[k] > 0), pft[k], np.nan)
        im = ax.imshow(field, origin="lower", extent=extent, cmap="YlGn", vmin=0, vmax=100,
                       interpolation="nearest")
        _map_axes(ax, extent)
        empty = "  (empty)" if not np.isfinite(field).any() else ""
        ax.set_title(f"{k}: {ELM_PFT_NAMES[k]}{empty}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.75)
    fig.suptitle(
        f"PCT_NAT_PFT (% of natural vegetation), NLCD {year} — SEUS 0.01° (~1 km)\n"
        "panels 9 and 10 are identical by construction (NLCD 51/52 split 50/50); "
        "cells with PCT_NATVEG = 0 are masked",
        fontsize=13)
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def fig_natveg_urban(natveg, urban, landfrac, extent, year, out_png):
    panels = [
        (natveg, "PCT_NATVEG", "YlGn", 100),
        (urban.sum(axis=0), "PCT_URBAN total (22+23+24)", "OrRd", 100),
        (urban[0], "PCT_URBAN low intensity (22)", "OrRd", 100),
        (urban[1], "PCT_URBAN medium intensity (23)", "OrRd", 100),
        (urban[2], "PCT_URBAN high intensity (24)", "OrRd", 100),
        (landfrac * 100.0, "valid NLCD pixels in cell", "Blues", 100),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    for (field, title, cmap, vmax), ax in zip(panels, axes.ravel()):
        im = ax.imshow(np.where(field > 0, field, np.nan), origin="lower", extent=extent,
                       cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")
        _map_axes(ax, extent)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.8, label="%")
    fig.suptitle(f"Vegetated / urban fractions, NLCD {year} — SEUS 0.01° (~1 km)", fontsize=15)
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def fig_dominant_pft(pft, natveg, extent, year, out_png):
    """Dominant natural PFT.

    PFT 9 and PFT 10 are equal *by construction* -- the shrub classes (51, 52)
    are split 50/50 between them -- so a plain argmax would hand every shrub
    cell to PFT 9 on numpy's lower-index tie-break. Those cells get their own
    "shrub 9/10 tie" category instead. Same check is applied generically, so any
    other exact tie is reported rather than silently broken.
    """
    dominant = np.argmax(pft, axis=0)
    top = np.max(pft, axis=0)
    n_tied = (pft == top[None, :, :]).sum(axis=0)
    tie = (n_tied > 1) & (top > 0)

    cats = [(k, ELM_PFT_NAMES[k], c) for k, c in [
        (0, "#9e9e9e"), (1, "#1b5e20"), (7, "#8bc34a"), (13, "#fff59d"),
        (14, "#ffb74d"), (15, "#a1887f"),
    ]]
    codes = np.full(pft.shape[1:], -1, dtype=np.int16)
    labels, colors = [], []
    for slot, (k, name, color) in enumerate(cats):
        codes[(dominant == k) & ~tie] = slot
        labels.append(f"{k}: {name}")
        colors.append(color)
    codes[tie] = len(labels)
    labels.append("tie (shrub 9/10 split 50/50)")
    colors.append("#ef8f00")
    # anything not covered above
    other = (codes < 0) & (natveg > 0)
    codes[other] = len(labels)
    labels.append("other PFT")
    colors.append("#000000")

    codes_f = codes.astype(float)
    codes_f[natveg <= 0] = np.nan

    cmap = ListedColormap(colors)
    cmap.set_bad("white")
    norm = BoundaryNorm(np.arange(-0.5, len(labels)), cmap.N)

    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(codes_f, origin="lower", extent=extent, cmap=cmap, norm=norm,
                   interpolation="nearest")
    _map_axes(ax, extent)
    ax.set_title(f"Dominant ELM natural PFT, NLCD {year} — SEUS 0.01° (~1 km)")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    cb = fig.colorbar(im, ax=ax, ticks=range(len(labels)), shrink=0.85)
    cb.ax.set_yticklabels(labels, fontsize=8)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")
    tied_pct = 100.0 * tie.sum() / max((natveg > 0).sum(), 1)
    print(f"  exact ties: {tie.sum()} cells ({tied_pct:.1f}% of vegetated cells)")


# --------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--year", type=int, default=2023, help="NLCD year (1985-2023)")
    p.add_argument("--res", type=float, default=TARGET_RES,
                   help="target resolution in degrees (default 0.01 = the project's 1 km grid)")
    p.add_argument("--bbox", type=float, nargs=4, default=list(SEUS_BBOX),
                   metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
                   help="domain bounds; default is the SEUS_1_24deg box")
    p.add_argument("--nlcd-dir", default=NLCD_DIR)
    p.add_argument("--outdir", default=OUT_DIR)
    p.add_argument("--chunk-rows", type=int, default=2048,
                   help="source rows read per pass (file chunking is 512)")
    p.add_argument("--plot-only", action="store_true",
                   help="skip the binning and re-draw from the existing NetCDF")
    p.add_argument("--no-plot", action="store_true", help="write the NetCDF only")
    args = p.parse_args()

    bbox = tuple(args.bbox)
    os.makedirs(args.outdir, exist_ok=True)
    tag = f"{args.year}_SEUS_{args.res:g}deg"
    nc_path = os.path.join(args.outdir, f"nlcd_elmpft_{tag}.nc")

    print(f"year        : {args.year}")
    print(f"domain      : lon [{bbox[0]}, {bbox[1]}], lat [{bbox[2]}, {bbox[3]}]")
    print(f"resolution  : {args.res} deg")
    print(f"output      : {nc_path}")

    if args.plot_only:
        lon, lat, lc_pct, landfrac, pft, natveg, urban = read_netcdf(nc_path)
        print(f"re-plotting from {nc_path}")
    else:
        src = os.path.join(args.nlcd_dir, NLCD_TEMPLATE.format(year=args.year))
        if not os.path.exists(src):
            raise SystemExit(f"source not found: {src}")
        print(f"source      : {src}")
        t0 = time.time()
        lon, lat, lc_pct, landfrac = bin_nlcd_to_grid(src, bbox, args.res, args.chunk_rows)
        print(f"binning done in {time.time() - t0:.1f}s -> grid {lat.size} x {lon.size}")
        pft, natveg, urban = nlcd_pct_to_elmpft(lc_pct)
        write_netcdf(nc_path, args.year, args.res, lon, lat, lc_pct, landfrac,
                     pft, natveg, urban)

    if args.no_plot:
        return

    extent = (bbox[0], bbox[1], bbox[2], bbox[3])
    fig_dominant_landcover(lc_pct, extent, args.year,
                           os.path.join(args.outdir, f"fig1_dominant_landcover_{tag}.png"))
    fig_dominant_pft(pft, natveg, extent, args.year,
                     os.path.join(args.outdir, f"fig2_dominant_pft_{tag}.png"))
    fig_natveg_urban(natveg, urban, landfrac, extent, args.year,
                     os.path.join(args.outdir, f"fig3_natveg_urban_{tag}.png"))
    fig_pft_panels(pft, natveg, extent, args.year,
                   os.path.join(args.outdir, f"fig4_pct_nat_pft_panels_{tag}.png"))


if __name__ == "__main__":
    main()
