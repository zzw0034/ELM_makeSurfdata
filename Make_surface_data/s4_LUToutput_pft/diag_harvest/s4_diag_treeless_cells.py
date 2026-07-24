"""Where does the harvest that s4_2 discards actually go?

s4_2 drops harvest in any coarse cell whose NLCD tree weight is zero. That loss
splits into two very different causes:

  no-data   : NLCD natveg == 0 for the whole cell. The 1/24 deg product covers
              CONUS only, but the domain runs to 50N / -125E, so southern Canada
              and northern Mexico sit inside it with no land cover at all. LUH2
              still supplies harvest there. This is a domain-extent artefact,
              not a data disagreement.
  treeless  : NLCD natveg > 0 but tree cover == 0. NLCD sees the cell, vegetated,
              with no trees, while LUH2 reports wood harvest. This is a genuine
              product disagreement and the only part worth acting on.

Counts cells and harvest shares for both, per year, plus a 2014 map.
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
    ndom = int(domain.sum())
    print(f"domain coarse cells: {ndom}")

    def agg(field_hr_flat):
        num = np.bincount(id_flat, weights=field_hr_flat * a_hr_flat, minlength=nlat_c * nlon_c)
        o = np.zeros_like(num)
        np.divide(num, a_in, out=o, where=a_in > 0)
        return o

    years = np.arange(START_YEAR, END_YEAR + 1)
    rec = {k: np.zeros(years.size) for k in
           ["n_nodata", "n_treeless", "n_harv", "h_nodata", "h_treeless", "h_total"]}
    snap = {}

    for iy, year in enumerate(years):
        ti, pi = year - LUH2_YEAR0, year - PFT_YEAR0
        harv = sum(read(tr[v], ti) for v in HARV)[srt, :].reshape(-1)
        hA = harv * area_c * domain

        p = read(pft["PCT_NAT_PFT"], pi, TREE, slice(None), slice(None)) / 100.0
        veg = read(pft["PCT_NATVEG"], pi, slice(None), slice(None)) / 100.0
        tree = agg((p.sum(axis=0) * veg).reshape(-1))
        natveg = agg(veg.reshape(-1))

        has_h = domain & (hA > 0)
        nodata = has_h & (natveg <= 1e-9)                      # NLCD sees nothing at all
        treeless = has_h & (natveg > 1e-9) & (tree <= 1e-6)    # vegetated but no trees

        rec["n_harv"][iy] = has_h.sum()
        rec["n_nodata"][iy] = nodata.sum()
        rec["n_treeless"][iy] = treeless.sum()
        rec["h_total"][iy] = hA[has_h].sum()
        rec["h_nodata"][iy] = hA[nodata].sum()
        rec["h_treeless"][iy] = hA[treeless].sum()

        if year == 2014:
            snap = {"nodata": nodata.copy(), "treeless": treeless.copy(),
                    "has_h": has_h.copy(), "hA": hA.copy()}
        if year % 25 == 0 or year == END_YEAR:
            print(f"{year}: cells w/ harvest={int(rec['n_harv'][iy]):5d} | "
                  f"no-data={int(rec['n_nodata'][iy]):4d} ({100*rec['h_nodata'][iy]/rec['h_total'][iy]:5.2f}% of harvest) | "
                  f"treeless={int(rec['n_treeless'][iy]):4d} ({100*rec['h_treeless'][iy]/rec['h_total'][iy]:5.2f}%)")

    pft.close(); st.close(); tr.close()

    fig, ax = plt.subplots(2, 1, figsize=(11, 9))
    a = ax[0]
    a.plot(years, 100 * rec["h_nodata"] / rec["h_total"], color="#ff7f0e", lw=2,
           label="outside NLCD coverage (Canada / Mexico) - domain artefact")
    a.plot(years, 100 * rec["h_treeless"] / rec["h_total"], color="#d62728", lw=2,
           label="NLCD vegetated but treeless - real disagreement")
    a.plot(years, 100 * (rec["h_nodata"] + rec["h_treeless"]) / rec["h_total"], color="k", lw=1.2, ls="--",
           label="total discarded by s4_2")
    a.set_ylabel("% of LUH2 harvest in domain")
    a.set_title("What s4_2 throws away, split by cause", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    a = ax[1]
    a.plot(years, rec["n_nodata"], color="#ff7f0e", lw=2, label="no-data cells")
    a.plot(years, rec["n_treeless"], color="#d62728", lw=2, label="treeless cells")
    a.plot(years, rec["n_harv"], color="#7f7f7f", lw=1.2, ls=":", label="all cells with harvest")
    a.set_xlabel("year"); a.set_ylabel("number of 0.25$\\degree$ cells")
    a.set_title(f"Cell counts (domain has {ndom} coarse cells)", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "treeless_cells.png", dpi=140, bbox_inches="tight")

    # 2014 map
    ml = (lat_c >= lat_hr.min()) & (lat_c <= lat_hr.max())
    mo = (lon_c >= lon_hr.min()) & (lon_c <= lon_hr.max())
    ext = [lon_c[mo].min(), lon_c[mo].max(), lat_c[ml].min(), lat_c[ml].max()]
    cat = np.zeros(nlat_c * nlon_c)
    cat[snap["has_h"]] = 1
    cat[snap["treeless"]] = 2
    cat[snap["nodata"]] = 3
    img = cat.reshape(nlat_c, nlon_c)[np.ix_(ml, mo)]
    fig, axm = plt.subplots(figsize=(11, 5.5))
    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#9ecae1", "#d62728", "#ff7f0e"])
    im = axm.imshow(img, origin="lower", extent=ext, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto")
    cb = plt.colorbar(im, ax=axm, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["no harvest", "harvest placed OK", "treeless (real)", "no NLCD data"])
    axm.set_title("2014: where LUH2 harvest cannot be placed", fontsize=11)
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "treeless_cells_map.png", dpi=140, bbox_inches="tight")

    i = -1
    print("\n" + "=" * 72)
    print(f"2014: {int(rec['n_harv'][i])} coarse cells carry harvest "
          f"({100*rec['n_harv'][i]/ndom:.1f}% of the {ndom}-cell domain)")
    print(f"  no-data  : {int(rec['n_nodata'][i]):4d} cells "
          f"({100*rec['n_nodata'][i]/rec['n_harv'][i]:.1f}% of them), "
          f"{100*rec['h_nodata'][i]/rec['h_total'][i]:.2f}% of harvest")
    print(f"  treeless : {int(rec['n_treeless'][i]):4d} cells "
          f"({100*rec['n_treeless'][i]/rec['n_harv'][i]:.1f}% of them), "
          f"{100*rec['h_treeless'][i]/rec['h_total'][i]:.2f}% of harvest")
    print("=" * 72)


if __name__ == "__main__":
    main()
