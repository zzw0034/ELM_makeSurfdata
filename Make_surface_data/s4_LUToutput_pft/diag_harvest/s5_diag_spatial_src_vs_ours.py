"""Spatial-pattern comparison: source (s4_1 NLCD) vs our final ELM product.

Per-cell 'fraction of gridcell' for each veg type:
  source: comp * PCT_NATVEG_src/100                     (NA grid, SEUS box)
  ours  : comp * PCT_NATVEG_landunit/100 * LANDFRAC_PFT  (SEUS grid)
Source is nearest-neighbour sampled onto the SEUS grid (half-cell offset, same res)
so the two can be differenced cell-by-cell. Reports spatial correlation + maps.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
LU = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/surfdata_results/landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260723.nc"
SURF = "/projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map/surfdata_SEUS_1_24deg_simyr1850_c260723.nc"
DIAG = Path("/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest")
TYPES = {"tree": list(range(1, 9)), "grass": [13, 14], "crop": [15]}
YEARS = {1850: 0, 2023: 173}


def rd(ds, v, *idx):
    a = ds[v][idx] if idx else ds[v][:]
    return np.asarray(np.ma.filled(np.ma.masked_invalid(a), 0.0), dtype=np.float64)


def main() -> None:
    s = nc.Dataset(SRC); o = nc.Dataset(LU); sf = nc.Dataset(SURF)
    slat = rd(s, "lat"); slon = rd(s, "lon")
    olat = rd(sf, "LATIXY")[:, 0]; olon = rd(sf, "LONGXY")[0, :]
    # nearest-neighbour index maps: for each SEUS cell, nearest source cell
    li = np.abs(olat[:, None] - slat[None, :]).argmin(axis=1)
    oi = np.abs(olon[:, None] - slon[None, :]).argmin(axis=1)

    nv_lu = rd(sf, "PCT_NATVEG") / 100.0
    lfrac = rd(sf, "LANDFRAC_PFT")
    ours_mask = lfrac > 1e-6

    def src_frac(t, idxs):
        pft = rd(s, "PCT_NAT_PFT", t, slice(None), slice(None), slice(None))
        nv = rd(s, "PCT_NATVEG", t, slice(None), slice(None)) / 100.0
        f = pft[idxs].sum(0) / 100.0 * nv               # (NA lat,lon)
        return f[np.ix_(li, oi)]                        # sampled to SEUS grid

    def ours_frac(t, idxs):
        pft = rd(o, "PCT_NAT_PFT", t, slice(None), slice(None), slice(None))
        return pft[idxs].sum(0) / 100.0 * nv_lu * lfrac  # (SEUS lat,lon)

    print(f"{'year':6s}{'type':6s}{'spatial r':>10s}{'mean|diff|(frac)':>18s}{'bias(ours-src)':>16s}")
    fields = {}
    for y, t in YEARS.items():
        for name, idxs in TYPES.items():
            fs = src_frac(t, idxs); fo = ours_frac(t, idxs)
            m = ours_mask & (fs + fo > 1e-9)
            r = np.corrcoef(fs[m], fo[m])[0, 1] if m.sum() > 10 else np.nan
            mad = np.abs(fo[m] - fs[m]).mean()
            bias = (fo[m] - fs[m]).mean()
            print(f"{y:<6d}{name:6s}{r:10.3f}{mad:18.4f}{bias:+16.4f}")
            fields[(y, name)] = (fs, fo)
    s.close(); o.close(); sf.close()

    # maps for 2023 tree & crop
    ext = [olon.min(), olon.max(), olat.min(), olat.max()]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for row, name in enumerate(["tree", "crop"]):
        fs, fo = fields[(2023, name)]
        fs = np.where(ours_mask, fs, np.nan); fo = np.where(ours_mask, fo, np.nan)
        im = ax[row, 0].imshow(fs, origin="lower", extent=ext, vmin=0, vmax=1, cmap="YlGn", aspect="auto")
        ax[row, 0].set_title(f"source {name} frac, 2023"); plt.colorbar(im, ax=ax[row, 0])
        im = ax[row, 1].imshow(fo, origin="lower", extent=ext, vmin=0, vmax=1, cmap="YlGn", aspect="auto")
        ax[row, 1].set_title(f"ours {name} frac, 2023"); plt.colorbar(im, ax=ax[row, 1])
        d = fo - fs
        im = ax[row, 2].imshow(d, origin="lower", extent=ext, vmin=-0.3, vmax=0.3, cmap="RdBu_r", aspect="auto")
        ax[row, 2].set_title(f"ours - source {name}, 2023"); plt.colorbar(im, ax=ax[row, 2])
    fig.tight_layout()
    fig.savefig(DIAG / "spatial_src_vs_ours_2023.png", dpi=140, bbox_inches="tight")
    print("wrote spatial_src_vs_ours_2023.png")


if __name__ == "__main__":
    main()
