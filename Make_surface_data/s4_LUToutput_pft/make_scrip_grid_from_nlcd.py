#!/usr/bin/env python3
"""
Create a SCRIP-format grid file from the lat/lon of an NLCD (or LUT) netCDF.
Use this as the SOURCE grid when building mapping files: NLCD grid (North America 1/24 deg) -> target grid.

Input:  NLCD file with 1D lat, lon (e.g. elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc)
Output: SCRIP grid netCDF in scr_out/nlcd_NA_1_24deg_grids.nc
"""

# %%
import os
import numpy as np
from netCDF4 import Dataset

# %%
# NLCD file that defines the grid (same grid as your LUT output)
NLCD_FILE = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "scr_out")
OUT_FILE = os.path.join(OUT_DIR, "SCRPgrid_NLCD_PFTs_NA_1_24deg.nc")
# %%


def _derive_grid_imask(nc: Dataset, nlat: int, nlon: int) -> np.ndarray:
    """Return SCRIP grid_imask flattened as 1=active, 0=inactive.

    Priority:
      1. LANDMASK if present
      2. PCT_NATVEG(time=0) > 0
      3. all ones
    """
    if "LANDMASK" in nc.variables:
        landmask = np.asarray(nc.variables["LANDMASK"][:], dtype=np.float64)
        if landmask.shape != (nlat, nlon):
            raise ValueError(
                f"LANDMASK shape {landmask.shape} does not match (nlat, nlon)=({nlat}, {nlon})"
            )
        return (landmask > 0).astype("i4").ravel()

    if "PCT_NATVEG" in nc.variables:
        pct_natveg = nc.variables["PCT_NATVEG"]
        veg = np.asarray(
            pct_natveg[0, :, :] if getattr(pct_natveg, "ndim", 0) == 3 else pct_natveg[:],
            dtype=np.float64,
        )
        if veg.shape != (nlat, nlon):
            raise ValueError(
                f"PCT_NATVEG shape {veg.shape} does not match (nlat, nlon)=({nlat}, {nlon})"
            )
        return (veg > 0).astype("i4").ravel()

    return np.ones(nlat * nlon, dtype="i4")


def make_scrip_from_lat_lon(
    lat_1d: np.ndarray, lon_1d: np.ndarray, outfile: str, imask: np.ndarray | None = None
) -> str:
    """
    Build SCRIP grid from 1D lat, lon (cell centers).
    Corners are inferred from cell spacing (half-width to neighbors).
    """
    lat_1d = np.asarray(lat_1d, dtype=np.float64)
    lon_1d = np.asarray(lon_1d, dtype=np.float64)
    nlat = lat_1d.size
    nlon = lon_1d.size
    grid_size = nlat * nlon
    grid_corners = 4
    grid_rank = 2

    # Spacing: use absolute median step so boundary half-cells extend the right way
    # when lat is descending (north -> south), np.diff(lat) is negative — still OK here.
    dlat = np.median(np.diff(lat_1d)) if nlat > 1 else 1.0 / 24.0
    dlon = np.median(np.diff(lon_1d)) if nlon > 1 else 1.0 / 24.0
    dlat_h = 0.5 * (abs(dlat) if abs(dlat) > 1e-9 else 1.0 / 24.0)
    dlon_h = 0.5 * (abs(dlon) if abs(dlon) > 1e-9 else 1.0 / 24.0)

    # Cell edges (interfaces between centers); midpoints are correct for asc or desc lat/lon.
    lat_edges = np.zeros(nlat + 1)
    if nlat > 1:
        lat_edges[1:-1] = (lat_1d[:-1] + lat_1d[1:]) / 2
        lat_edges[0] = lat_1d[0] + dlat_h * np.sign(lat_1d[0] - lat_edges[1])
        lat_edges[-1] = lat_1d[-1] + dlat_h * np.sign(lat_1d[-1] - lat_edges[-2])
    else:
        lat_edges[0] = float(lat_1d[0]) - dlat_h
        lat_edges[1] = float(lat_1d[0]) + dlat_h

    lon_edges = np.zeros(nlon + 1)
    if nlon > 1:
        lon_edges[1:-1] = (lon_1d[:-1] + lon_1d[1:]) / 2
        lon_edges[0] = lon_1d[0] + dlon_h * np.sign(lon_1d[0] - lon_edges[1])
        lon_edges[-1] = lon_1d[-1] + dlon_h * np.sign(lon_1d[-1] - lon_edges[-2])
    else:
        lon_edges[0] = float(lon_1d[0]) - dlon_h
        lon_edges[1] = float(lon_1d[0]) + dlon_h

    # 2D centers (lat index i, lon index j) -> same order as meshgrid(lon, lat, indexing='ij') so lat is first
    LAT, LON = np.meshgrid(lat_1d, lon_1d, indexing="ij")
    clat = LAT.ravel()
    clon = LON.ravel()

    # Corners: SW -> SE -> NE -> NW (counterclockwise for ESMF)
    # Cell (i,j): lat_edges[i]..lat_edges[i+1], lon_edges[j]..lon_edges[j+1]
    grid_corner_lat = np.zeros((grid_size, grid_corners))
    grid_corner_lon = np.zeros((grid_size, grid_corners))
    for i in range(nlat):
        for j in range(nlon):
            gc = i * nlon + j
            # Geographic bounds: lat_edges[i] and lat_edges[i+1] need not be south<north
            # when the file stores lat descending; min/max fixes CCW corners for ESMF.
            lat_lo = min(lat_edges[i], lat_edges[i + 1])
            lat_hi = max(lat_edges[i], lat_edges[i + 1])
            lon_lo = min(lon_edges[j], lon_edges[j + 1])
            lon_hi = max(lon_edges[j], lon_edges[j + 1])
            # SW, SE, NE, NW (CCW, increasing lon along south edge)
            grid_corner_lat[gc, 0] = lat_lo
            grid_corner_lat[gc, 1] = lat_lo
            grid_corner_lat[gc, 2] = lat_hi
            grid_corner_lat[gc, 3] = lat_hi
            grid_corner_lon[gc, 0] = lon_lo
            grid_corner_lon[gc, 1] = lon_hi
            grid_corner_lon[gc, 2] = lon_hi
            grid_corner_lon[gc, 3] = lon_lo

    grid_dims = np.array([nlon, nlat], dtype="i4")
    if imask is None:
        imask = np.ones(grid_size, dtype="i4")
    else:
        imask = np.asarray(imask, dtype="i4").reshape(grid_size)

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
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

        v_dims[:] = grid_dims
        v_imask[:] = imask
        v_center_lat[:] = clat
        v_center_lon[:] = clon
        v_corner_lat[:, :] = grid_corner_lat
        v_corner_lon[:, :] = grid_corner_lon

        v_center_lat.units = "degrees"
        v_center_lon.units = "degrees"
        v_corner_lat.units = "degrees"
        v_corner_lon.units = "degrees"
        nc.title = "SCRIP grid from NLCD lat/lon (North America 1/24 deg)"
        nc.conventions = "SCRIP"
    finally:
        nc.close()

    print(f"Wrote {outfile}")
    print(f"  nlon={nlon}, nlat={nlat}, grid_size={grid_size}")
    print(
        f"  lat [{lat_1d.min():.4f}, {lat_1d.max():.4f}], lon [{lon_1d.min():.4f}, {lon_1d.max():.4f}]"
    )
    return outfile


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Create SCRIP grid from NLCD file lat/lon (for mapping: NLCD grid -> target)."
    )
    p.add_argument(
        "--nlcd", default=NLCD_FILE, help="Path to NLCD netCDF (with lat, lon)"
    )
    p.add_argument("-o", "--out", default=OUT_FILE, help="Output SCRIP grid path")
    args = p.parse_args()
    if not os.path.isfile(args.nlcd):
        raise FileNotFoundError(
            f"NLCD file not found: {args.nlcd}\nUse --nlcd /path/to/file.nc"
        )
    nc = Dataset(args.nlcd, "r")
    try:
        lat = np.asarray(nc.variables["lat"][:])
        lon = np.asarray(nc.variables["lon"][:])
        imask = _derive_grid_imask(nc, lat.size, lon.size)
    finally:
        nc.close()

    if lat.size > 1 and lat[0] > lat[1]:
        print(
            "Note: lat is descending (north to south); corner bounds use min/max per cell."
        )

    print(
        f"Using grid_imask with active={int(np.sum(imask > 0))}, inactive={int(np.sum(imask == 0))}"
    )
    make_scrip_from_lat_lon(lat, lon, args.out, imask=imask)
