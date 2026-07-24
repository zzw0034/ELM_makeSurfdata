"""Compare the fixed s4_2 output (_07232026) against the old buggy files (_07112026).

Old files were moved to output_downscaled_luh2_harvest/prev_buggy_c0711/.
Only HARVEST_* changed; PCT_PFT / LANDMASK / GRAZING are identical (verified
elsewhere). This makes a 3-panel figure: domain-total harvest (old vs new),
the year-by-year ratio, and 2000 spatial maps (new field + new/old ratio map)
to show the spatial pattern is preserved and only rescaled.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest"
)
OLD_SUB = BASE / "prev_buggy_c0711"
DIAG_DIR = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest"
)
OLD_TAG, NEW_TAG = "07112026", "07232026"
HARV = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3"]
Y0, Y1 = 1850, 2023
R = 6371.0


def rd(v):
    return np.asarray(np.ma.filled(np.ma.masked_invalid(v[:]), 0.0), dtype=np.float64)


def hsum(path):
    d = nc.Dataset(path)
    tot = sum(rd(d[v]) for v in HARV)
    d.close()
    return tot


def main() -> None:
    years = np.arange(Y0, Y1 + 1)

    d0 = nc.Dataset(BASE / f"LUT_nlcd2elm_luh2_historical_2000_{NEW_TAG}.nc")
    lat = rd(d0["LAT"]); lon = rd(d0["LON"])
    d0.close()
    dlat = abs(np.median(np.diff(lat))); dlon = abs(np.median(np.diff(lon)))
    n = np.deg2rad(lat + 0.5 * dlat); s = np.deg2rad(lat - 0.5 * dlat)
    A = np.repeat(((R**2) * np.deg2rad(dlon) * (np.sin(n) - np.sin(s)))[:, None], lon.size, 1)

    old_tot = np.zeros(years.size)
    new_tot = np.zeros(years.size)
    for iy, y in enumerate(years):
        o = hsum(OLD_SUB / f"LUT_nlcd2elm_luh2_historical_{y}_{OLD_TAG}.nc")
        nw = hsum(BASE / f"LUT_nlcd2elm_luh2_historical_{y}_{NEW_TAG}.nc")
        old_tot[iy] = (o * A).sum()
        new_tot[iy] = (nw * A).sum()
        if y % 25 == 0 or y == Y1:
            print(f"{y}: old={old_tot[iy]:10.1f} new={new_tot[iy]:10.1f} ratio={new_tot[iy]/old_tot[iy] if old_tot[iy]>0 else float('nan'):.2f}")

    # 2000 spatial fields
    o2000 = hsum(OLD_SUB / f"LUT_nlcd2elm_luh2_historical_2000_{OLD_TAG}.nc")
    n2000 = hsum(BASE / f"LUT_nlcd2elm_luh2_historical_2000_{NEW_TAG}.nc")

    fig = plt.figure(figsize=(13, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.7, 1.3])

    ax = fig.add_subplot(gs[0, :])
    ax.plot(years, old_tot, color="#7f7f7f", lw=2, label="old (buggy, _07112026)")
    ax.plot(years, new_tot, color="#d62728", lw=2, label="new (fixed, _07232026)")
    ax.set_yscale("log"); ax.set_ylabel("$\\Sigma$ HARVEST $\\times$ cell area")
    ax.set_title("s4_2 harvest, old vs new (domain total, 5 categories summed)", fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

    ax = fig.add_subplot(gs[1, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(old_tot > 0, new_tot / old_tot, np.nan)
    ax.plot(years, ratio, color="#1f77b4", lw=2, label="new / old")
    ax.axhline(36, color="k", ls=":", lw=1.5, label="coarse/fine cell ratio = 36")
    ax.set_ylabel("ratio"); ax.set_xlabel("year")
    ax.set_ylim(0, 45)
    ax.set_title("Magnitude correction factor (area-conservation fix + veg_frac + 1-yr shift)", fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # SEUS-ish crop of the map for legibility
    ml = (lat >= 25) & (lat <= 50); mo = (lon >= -125) & (lon <= -65)
    ext = [lon[mo].min(), lon[mo].max(), lat[ml].min(), lat[ml].max()]

    def box(f):
        return f[np.ix_(ml, mo)]

    axn = fig.add_subplot(gs[2, 0])
    fld = box(n2000); fld = np.where(fld > 0, fld, np.nan)
    im = axn.imshow(fld, origin="lower", extent=ext, cmap="magma_r",
                    norm=matplotlib.colors.LogNorm(vmin=1e-5, vmax=np.nanmax(fld)), aspect="auto")
    axn.set_title("new HARVEST sum, 2000 (fraction of veg unit)", fontsize=10)
    plt.colorbar(im, ax=axn)

    axr = fig.add_subplot(gs[2, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        rmap = np.where(box(o2000) > 1e-9, box(n2000) / box(o2000), np.nan)
    im = axr.imshow(rmap, origin="lower", extent=ext, cmap="RdBu_r", vmin=30, vmax=42, aspect="auto")
    axr.set_title("new / old ratio, 2000  (~uniform => pattern preserved)", fontsize=10)
    plt.colorbar(im, ax=axr)

    fig.tight_layout()
    fig.savefig(DIAG_DIR / "old_vs_new_harvest.png", dpi=140, bbox_inches="tight")
    print("wrote old_vs_new_harvest.png")
    print(f"\n1850-2023 total: old={old_tot.sum():.1f}  new={new_tot.sum():.1f}  new/old={new_tot.sum()/old_tot.sum():.2f}")
    rr = ratio[np.isfinite(ratio)]
    print(f"per-year ratio: min={rr.min():.2f} median={np.median(rr):.2f} max={rr.max():.2f}")


if __name__ == "__main__":
    main()
