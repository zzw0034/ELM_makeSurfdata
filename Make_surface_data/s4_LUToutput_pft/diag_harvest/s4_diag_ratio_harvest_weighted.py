"""How much does the LUH2/NLCD forest disagreement actually matter for harvest?

The per-cell ratio LUH2_forest / NLCD_tree scales the local harvest intensity that
ELM ends up applying, 1:1. Unweighted per-cell statistics (r, mean|diff|) treat a
desert cell with 5% tree cover the same as a Southeast timber cell with 85%, but
harvest is overwhelmingly concentrated in the latter. So the decision-relevant
quantity is the distribution of that ratio *weighted by where the harvest is*.

Reports, per year: the harvest-weighted mean ratio, its weighted percentiles, and
the share of all harvest landing in cells where the implied intensity is off by
more than 2x in either direction (plus the share landing on zero NLCD tree cover,
which s4_2 discards outright).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft"
)
PFT_PATH = BASE / "scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
LUH_DIR = Path("/projects/hpcl-cli185/proj-shared/zw5/luh")
DIAG_DIR = BASE / "diag_harvest"

START_YEAR, END_YEAR = 1850, 2014
LUH2_YEAR0, PFT_YEAR0 = 850, 1850
R_EARTH = 6371.0
TREE = [1, 2, 3, 4, 5, 6, 7, 8]
HARV = ["primf_harv", "primn_harv", "secmf_harv", "secyf_harv", "secnf_harv"]


def read(var, *idx):
    arr = var[idx] if idx else var[:]
    return np.asarray(np.ma.filled(np.ma.masked_invalid(arr), 0.0), dtype=np.float64)


def cell_area(lat_1d, dlat, dlon, nlon):
    n = np.deg2rad(lat_1d + 0.5 * abs(dlat))
    s = np.deg2rad(lat_1d - 0.5 * abs(dlat))
    return np.repeat((R_EARTH**2 * np.deg2rad(abs(dlon)) * (np.sin(n) - np.sin(s)))[:, None], nlon, 1)


def wpct(values, weights, qs):
    """Weighted percentiles."""
    o = np.argsort(values)
    v, w = values[o], weights[o]
    c = np.cumsum(w)
    if c[-1] <= 0:
        return np.full(len(qs), np.nan)
    return np.interp(np.array(qs) / 100.0 * c[-1], c, v)


def main() -> None:
    pft = nc.Dataset(PFT_PATH, "r")
    st = nc.Dataset(LUH_DIR / "states.nc", "r")
    tr = nc.Dataset(LUH_DIR / "transitions.nc", "r")

    lat_hr = np.asarray(pft["lat"][:], float)
    lon_hr = np.asarray(pft["lon"][:], float)
    area_hr = cell_area(lat_hr, np.median(np.diff(lat_hr)), np.median(np.diff(lon_hr)), lon_hr.size)

    lat_raw = np.asarray(st["lat"][:], float)
    lon_c = np.asarray(st["lon"][:], float)
    srt = np.argsort(lat_raw)
    lat_c = lat_raw[srt]
    nlat_c, nlon_c = lat_c.size, lon_c.size
    area_c = cell_area(lat_c, np.median(np.diff(lat_c)), np.median(np.diff(lon_c)), nlon_c).reshape(-1)

    dlat_c, dlon_c = np.median(np.diff(lat_c)), np.median(np.diff(lon_c))
    li = np.floor((lat_hr - (lat_c[0] - 0.5 * dlat_c)) / dlat_c).astype(np.int64)
    oi = np.floor((lon_hr - (lon_c[0] - 0.5 * dlon_c)) / dlon_c).astype(np.int64)
    inside = ((li >= 0) & (li < nlat_c))[:, None] & ((oi >= 0) & (oi < nlon_c))[None, :]
    cid = (li[:, None] * nlon_c + oi[None, :]).astype(np.int64)
    cid[~inside] = 0
    id_flat, in_flat = cid.reshape(-1), inside.reshape(-1)
    a_hr_flat = area_hr.reshape(-1) * in_flat
    a_in = np.bincount(id_flat, weights=a_hr_flat, minlength=nlat_c * nlon_c)
    domain = a_in > 0

    years = np.arange(START_YEAR, END_YEAR + 1)
    out = {k: np.zeros(years.size) for k in
           ["wmean", "p10", "p25", "p50", "p75", "p90", "off2x", "zero_tree", "unweighted_mean"]}
    keep = {}

    for iy, year in enumerate(years):
        ti, pi = year - LUH2_YEAR0, year - PFT_YEAR0
        luh_forest = (read(st["primf"], ti) + read(st["secdf"], ti))[srt, :].reshape(-1)
        harv = sum(read(tr[v], ti) for v in HARV)[srt, :].reshape(-1)

        p = read(pft["PCT_NAT_PFT"], pi, TREE, slice(None), slice(None)) / 100.0
        veg = read(pft["PCT_NATVEG"], pi, slice(None), slice(None)) / 100.0
        tree_hr = (p.sum(axis=0) * veg).reshape(-1)
        num = np.bincount(id_flat, weights=tree_hr * a_hr_flat, minlength=nlat_c * nlon_c)
        nlcd_tree = np.zeros_like(num)
        np.divide(num, a_in, out=nlcd_tree, where=a_in > 0)

        # harvest weight = harvested area in the cell
        w = (harv * area_c) * domain
        has_tree = nlcd_tree > 1e-6
        out["zero_tree"][iy] = 100 * w[~has_tree].sum() / w.sum() if w.sum() > 0 else 0.0

        m = domain & has_tree & (w > 0) & (luh_forest > 1e-6)
        ratio = luh_forest[m] / nlcd_tree[m]
        ww = w[m]
        out["wmean"][iy] = np.average(ratio, weights=ww)
        out["p10"][iy], out["p25"][iy], out["p50"][iy], out["p75"][iy], out["p90"][iy] = wpct(
            ratio, ww, [10, 25, 50, 75, 90]
        )
        out["off2x"][iy] = 100 * ww[(ratio > 2) | (ratio < 0.5)].sum() / ww.sum()
        mu = domain & has_tree & (luh_forest > 1e-6)
        out["unweighted_mean"][iy] = np.median(luh_forest[mu] / nlcd_tree[mu])

        if year == 2014:
            keep = {"ratio": ratio, "w": ww}
        if year % 25 == 0 or year == END_YEAR:
            print(f"{year}: w-mean={out['wmean'][iy]:.3f} median={out['p50'][iy]:.3f} "
                  f"[p10={out['p10'][iy]:.2f}, p90={out['p90'][iy]:.2f}] "
                  f"off-by->2x={out['off2x'][iy]:5.2f}%  onto-zero-tree={out['zero_tree'][iy]:5.2f}%")

    pft.close(); st.close(); tr.close()

    fig, ax = plt.subplots(2, 1, figsize=(11, 9))
    a = ax[0]
    a.fill_between(years, out["p10"], out["p90"], color="#1f77b4", alpha=0.18, label="harvest-weighted p10-p90")
    a.fill_between(years, out["p25"], out["p75"], color="#1f77b4", alpha=0.32, label="p25-p75")
    a.plot(years, out["p50"], color="#1f77b4", lw=2, label="harvest-weighted median")
    a.plot(years, out["unweighted_mean"], color="#7f7f7f", lw=1.5, ls="--", label="unweighted median (all cells)")
    a.axhline(1.0, color="k", ls=":", lw=1.5)
    a.set_yscale("log"); a.set_ylabel("LUH2 forest / NLCD tree")
    a.set_title("Per-cell forest disagreement, weighted by where the harvest actually is\n"
                "(this ratio scales the local harvest intensity ELM applies, 1:1)", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3, which="both")

    a = ax[1]
    a.plot(years, out["off2x"], color="#d62728", lw=2, label="harvest landing where intensity is off by >2x")
    a.plot(years, out["zero_tree"], color="#ff7f0e", lw=2, label="harvest landing on zero NLCD tree (discarded)")
    a.set_xlabel("year"); a.set_ylabel("% of total harvest")
    a.set_title("Share of harvest affected", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(DIAG_DIR / "ratio_harvest_weighted.png", dpi=140, bbox_inches="tight")

    print("\n" + "=" * 72)
    r, w = keep["ratio"], keep["w"]
    print("2014, weighted by harvested area:")
    for q in (5, 10, 25, 50, 75, 90, 95):
        print(f"   p{q:<3d} = {wpct(r, w, [q])[0]:.3f}")
    print(f"   within 2x    : {100*w[(r>=0.5)&(r<=2)].sum()/w.sum():.2f}% of harvest")
    print(f"   within 1.5x  : {100*w[(r>=1/1.5)&(r<=1.5)].sum()/w.sum():.2f}% of harvest")
    print(f"   unweighted, same cells: median={np.median(r):.3f}  within 2x={100*np.mean((r>=0.5)&(r<=2)):.2f}% of cells")
    print("=" * 72)


if __name__ == "__main__":
    main()
