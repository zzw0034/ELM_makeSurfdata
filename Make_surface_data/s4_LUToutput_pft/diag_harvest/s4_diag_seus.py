"""Same harvest diagnostics as the full-domain runs, restricted to the SEUS pilot box.

Box is the one HARMONIZATION_SEUS_PILOT.md 13.0 uses: lon -95..-74, lat 24..37.5
(intersected with the 1/24 deg grid, which starts at 25N). SEUS is the region the
harmonization work actually targets, and it is the most intensively managed forest
region in the country, so the numbers here matter more than the CONUS-wide ones.

Reports, per year:
  - forest area, LUH2 (primf+secdf) vs NLCD tree, aggregated to 0.25 deg
  - harvest-weighted distribution of LUH2_forest / NLCD_tree (scales local intensity)
  - harvest that cannot be placed, split into no-NLCD-data vs genuinely treeless
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

SEUS = dict(lon0=-95.0, lon1=-74.0, lat0=24.0, lat1=37.5)
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


def wpct(v, w, qs):
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    return np.interp(np.array(qs) / 100.0 * c[-1], c, v) if c[-1] > 0 else np.full(len(qs), np.nan)


def main() -> None:
    pft = nc.Dataset(PFT_PATH, "r")
    st = nc.Dataset(LUH_DIR / "states.nc", "r")
    tr = nc.Dataset(LUH_DIR / "transitions.nc", "r")

    lat_hr = np.asarray(pft["lat"][:], float)
    lon_hr = np.asarray(pft["lon"][:], float)
    area_hr = cell_area(lat_hr, np.median(np.diff(lat_hr)), np.median(np.diff(lon_hr)), lon_hr.size)
    # restrict the hi-res grid to the SEUS box; everything downstream follows
    box = ((lat_hr >= SEUS["lat0"]) & (lat_hr <= SEUS["lat1"]))[:, None] & (
        (lon_hr >= SEUS["lon0"]) & (lon_hr <= SEUS["lon1"])
    )[None, :]
    print(f"SEUS box: lat {SEUS['lat0']}-{SEUS['lat1']}, lon {SEUS['lon0']}-{SEUS['lon1']}")
    print(f"hi-res cells in box: {int(box.sum())} of {box.size}")

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
    inside &= box
    cid = (li[:, None] * nlon_c + oi[None, :]).astype(np.int64)
    cid[~inside] = 0
    id_flat, in_flat = cid.reshape(-1), inside.reshape(-1)
    a_hr_flat = area_hr.reshape(-1) * in_flat
    a_in = np.bincount(id_flat, weights=a_hr_flat, minlength=nlat_c * nlon_c)
    domain = a_in > 0
    ndom = int(domain.sum())
    print(f"SEUS coarse cells: {ndom}")

    def agg(f):
        num = np.bincount(id_flat, weights=f * a_hr_flat, minlength=nlat_c * nlon_c)
        o = np.zeros_like(num)
        np.divide(num, a_in, out=o, where=a_in > 0)
        return o

    years = np.arange(START_YEAR, END_YEAR + 1)
    R = {k: np.zeros(years.size) for k in
         ["luh_for", "nlcd_tree", "p10", "p50", "p90", "off2x", "wmean",
          "n_harv", "n_nodata", "n_treeless", "h_nodata", "h_treeless", "harv"]}
    snap = {}

    for iy, year in enumerate(years):
        ti, pi = year - LUH2_YEAR0, year - PFT_YEAR0
        luh_for = (read(st["primf"], ti) + read(st["secdf"], ti))[srt, :].reshape(-1)
        harv = sum(read(tr[v], ti) for v in HARV)[srt, :].reshape(-1)
        hA = harv * area_c * domain

        p = read(pft["PCT_NAT_PFT"], pi, TREE, slice(None), slice(None)) / 100.0
        veg = read(pft["PCT_NATVEG"], pi, slice(None), slice(None)) / 100.0
        tree = agg((p.sum(axis=0) * veg).reshape(-1))
        natveg = agg(veg.reshape(-1))

        R["luh_for"][iy] = (luh_for * area_c)[domain].sum()
        R["nlcd_tree"][iy] = (tree * area_c)[domain].sum()
        R["harv"][iy] = hA[domain].sum()

        has_h = domain & (hA > 0)
        nodata = has_h & (natveg <= 1e-9)
        treeless = has_h & (natveg > 1e-9) & (tree <= 1e-6)
        R["n_harv"][iy], R["n_nodata"][iy], R["n_treeless"][iy] = has_h.sum(), nodata.sum(), treeless.sum()
        tot = hA[has_h].sum()
        R["h_nodata"][iy] = 100 * hA[nodata].sum() / tot if tot > 0 else 0
        R["h_treeless"][iy] = 100 * hA[treeless].sum() / tot if tot > 0 else 0

        m = has_h & (tree > 1e-6) & (luh_for > 1e-6)
        if m.sum() > 0:
            ratio, w = luh_for[m] / tree[m], hA[m]
            R["wmean"][iy] = np.average(ratio, weights=w)
            R["p10"][iy], R["p50"][iy], R["p90"][iy] = wpct(ratio, w, [10, 50, 90])
            R["off2x"][iy] = 100 * w[(ratio > 2) | (ratio < 0.5)].sum() / w.sum()
            if year == 2014:
                snap = {"ratio": ratio, "w": w}

        if year % 25 == 0 or year == END_YEAR:
            print(f"{year}: LUH2 for={R['luh_for'][iy]/1e3:7.1f}e3 NLCD tree={R['nlcd_tree'][iy]/1e3:7.1f}e3 km2 "
                  f"ratio={R['nlcd_tree'][iy]/R['luh_for'][iy]:.3f} | w-median={R['p50'][iy]:.3f} "
                  f"off>2x={R['off2x'][iy]:5.2f}% | nodata={int(R['n_nodata'][iy]):3d}c/{R['h_nodata'][iy]:5.2f}% "
                  f"treeless={int(R['n_treeless'][iy]):3d}c/{R['h_treeless'][iy]:.2f}%")

    pft.close(); st.close(); tr.close()

    fig, ax = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    a = ax[0]
    a.plot(years, R["luh_for"] / 1e3, color="#1f77b4", lw=2, label="LUH2 forest (primf+secdf)")
    a.plot(years, R["nlcd_tree"] / 1e3, color="#d62728", lw=2, label="NLCD tree cover")
    a.set_ylabel("area (10$^3$ km$^2$)")
    a.set_title(f"SEUS pilot box (lon {SEUS['lon0']}..{SEUS['lon1']}, lat {SEUS['lat0']}..{SEUS['lat1']}), "
                f"{ndom} coarse cells", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    a = ax[1]
    a.fill_between(years, R["p10"], R["p90"], color="#1f77b4", alpha=0.25, label="harvest-weighted p10-p90")
    a.plot(years, R["p50"], color="#1f77b4", lw=2, label="harvest-weighted median")
    a.plot(years, R["luh_for"] / R["nlcd_tree"], color="#2ca02c", lw=1.5, ls="--", label="ratio of SEUS totals")
    a.axhline(1.0, color="k", ls=":", lw=1.5)
    a.set_ylabel("LUH2 forest / NLCD tree")
    a.set_title("Local harvest-intensity scaling", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    a = ax[2]
    a.plot(years, R["off2x"], color="#d62728", lw=2, label="harvest where intensity off by >2x")
    a.plot(years, R["h_nodata"], color="#ff7f0e", lw=2, label="harvest onto no-NLCD-data cells")
    a.plot(years, R["h_treeless"], color="#9467bd", lw=2, label="harvest onto vegetated-but-treeless cells")
    a.set_xlabel("year"); a.set_ylabel("% of SEUS harvest")
    a.legend(fontsize=9); a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "seus_harvest_diag.png", dpi=140, bbox_inches="tight")

    i = -1
    print("\n" + "=" * 72)
    print(f"SEUS 2014: {int(R['n_harv'][i])} of {ndom} coarse cells carry harvest")
    print(f"  LUH2 forest = {R['luh_for'][i]/1e3:.1f}e3 km2   NLCD tree = {R['nlcd_tree'][i]/1e3:.1f}e3 km2"
          f"   ratio = {R['nlcd_tree'][i]/R['luh_for'][i]:.3f}")
    print(f"  harvest intensity: intended={100*R['harv'][i]/R['luh_for'][i]:.4f}%/yr  "
          f"applied={100*R['harv'][i]/R['nlcd_tree'][i]:.4f}%/yr")
    if snap:
        r, w = snap["ratio"], snap["w"]
        print("  harvest-weighted ratio percentiles:",
              "  ".join(f"p{q}={wpct(r,w,[q])[0]:.3f}" for q in (10, 25, 50, 75, 90)))
        print(f"  within 2x = {100*w[(r>=0.5)&(r<=2)].sum()/w.sum():.2f}% of harvest, "
              f"within 1.5x = {100*w[(r>=1/1.5)&(r<=1.5)].sum()/w.sum():.2f}%")
    print(f"  no-data cells = {int(R['n_nodata'][i])} ({R['h_nodata'][i]:.2f}% of harvest)")
    print(f"  treeless cells= {int(R['n_treeless'][i])} ({R['h_treeless'][i]:.2f}% of harvest)")
    print("=" * 72)


if __name__ == "__main__":
    main()
