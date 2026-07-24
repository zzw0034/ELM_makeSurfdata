"""Compare the LUH2 forest/natural-land area against the NLCD-derived one, 1850-2014.

s4_2 takes wood harvest from LUH2 but applies it to NLCD-derived land cover. LUH2's
harvest is calibrated against LUH2's own forest, so if the two products disagree on
how much forest exists, the harvest rate ELM actually sees is biased by that ratio.
This quantifies the disagreement.

Everything is aggregated NLCD -> 0.25 deg (fine to coarse, area-weighted), never the
other way, and compared only over coarse cells that overlap the 1/24 deg domain.

Compared pairs
  forest    : LUH2 primf+secdf                vs  NLCD tree cover (PFT 1-8)
  natural   : LUH2 primf+secdf+primn+secdn    vs  NLCD natveg excluding cropland
              (the second pair is the harvest normalisation denominator the E3SM
               doc specifies: "normalized based on the sum of LUH primary and
               secondary land")

Caveat carried into the plots: LUH2 secdf is "potentially forested secondary land",
a land-use category, not observed canopy. It is an upper bound on actual tree cover,
so some of the gap is definitional rather than error.
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

START_YEAR, END_YEAR = 1850, 2014  # transitions.nc ends 2014
LUH2_YEAR0, PFT_YEAR0 = 850, 1850
R_EARTH = 6371.0

TREE = [1, 2, 3, 4, 5, 6, 7, 8]
CROP = 15
HARV = ["primf_harv", "primn_harv", "secmf_harv", "secyf_harv", "secnf_harv"]
LUH_CROP = ["c3ann", "c4ann", "c3per", "c4per", "c3nfx"]


def read(var, *idx) -> np.ndarray:
    arr = var[idx] if idx else var[:]
    return np.asarray(np.ma.filled(np.ma.masked_invalid(arr), 0.0), dtype=np.float64)


def cell_area(lat_1d, dlat, dlon, nlon):
    lat_n = np.deg2rad(lat_1d + 0.5 * abs(dlat))
    lat_s = np.deg2rad(lat_1d - 0.5 * abs(dlat))
    band = R_EARTH**2 * np.deg2rad(abs(dlon)) * (np.sin(lat_n) - np.sin(lat_s))
    return np.repeat(band[:, None], nlon, axis=1)


def build_index(lat_hr, lon_hr, lat_c, lon_c):
    dlat = float(np.median(np.diff(lat_c)))
    dlon = float(np.median(np.diff(lon_c)))
    lat0, lon0 = lat_c[0] - 0.5 * dlat, lon_c[0] - 0.5 * dlon
    li = np.floor((lat_hr - lat0) / dlat).astype(np.int64)
    oi = np.floor((lon_hr - lon0) / dlon).astype(np.int64)
    ok_l = (li >= 0) & (li < lat_c.size)
    ok_o = (oi >= 0) & (oi < lon_c.size)
    inside = ok_l[:, None] & ok_o[None, :]
    cid = (li[:, None] * lon_c.size + oi[None, :]).astype(np.int64)
    cid[~inside] = 0
    return cid, inside, lat_c.size * lon_c.size


def main() -> None:
    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    pft = nc.Dataset(PFT_PATH, "r")
    st = nc.Dataset(LUH_DIR / "states.nc", "r")
    tr = nc.Dataset(LUH_DIR / "transitions.nc", "r")

    lat_hr = np.asarray(pft["lat"][:], dtype=np.float64)
    lon_hr = np.asarray(pft["lon"][:], dtype=np.float64)
    area_hr = cell_area(
        lat_hr, np.median(np.diff(lat_hr)), np.median(np.diff(lon_hr)), lon_hr.size
    )

    lat_raw = np.asarray(st["lat"][:], dtype=np.float64)
    lon_c = np.asarray(st["lon"][:], dtype=np.float64)
    srt = np.argsort(lat_raw)
    lat_c = lat_raw[srt]
    nlat_c, nlon_c = lat_c.size, lon_c.size
    area_c = cell_area(lat_c, np.median(np.diff(lat_c)), np.median(np.diff(lon_c)), nlon_c)
    area_c_flat = area_c.reshape(-1)

    cid, inside, ncoarse = build_index(lat_hr, lon_hr, lat_c, lon_c)
    id_flat, in_flat = cid.reshape(-1), inside.reshape(-1)
    a_hr_flat = area_hr.reshape(-1) * in_flat
    # area of hi-res cells landing in each coarse cell (the aggregation denominator)
    a_in_coarse = np.bincount(id_flat, weights=a_hr_flat, minlength=ncoarse)
    domain = a_in_coarse > 0
    print(f"coarse cells in domain: {int(domain.sum())}")

    def to_coarse(field_hr_flat):
        """Area-weighted mean of a hi-res fraction field onto the coarse grid."""
        num = np.bincount(id_flat, weights=field_hr_flat * a_hr_flat, minlength=ncoarse)
        out = np.zeros(ncoarse)
        np.divide(num, a_in_coarse, out=out, where=a_in_coarse > 0)
        return out

    years = np.arange(START_YEAR, END_YEAR + 1)
    keys = ["luh_forest", "luh_nat", "nlcd_tree", "nlcd_nat", "luh_crop", "luh_past", "harv"]
    ts = {k: np.zeros(years.size) for k in keys}
    snap = {}

    for iy, year in enumerate(years):
        ti, pi = year - LUH2_YEAR0, year - PFT_YEAR0

        sub = read(st["primf"], ti) + read(st["secdf"], ti)
        luh_forest = sub[srt, :].reshape(-1)
        luh_nat = (
            sub + read(st["primn"], ti) + read(st["secdn"], ti)
        )[srt, :].reshape(-1)
        luh_crop = sum(read(st[v], ti) for v in LUH_CROP)[srt, :].reshape(-1)
        luh_past = (read(st["pastr"], ti) + read(st["range"], ti))[srt, :].reshape(-1)
        harv = sum(read(tr[v], ti) for v in HARV)[srt, :].reshape(-1)

        p = read(pft["PCT_NAT_PFT"], pi, TREE + [CROP], slice(None), slice(None)) / 100.0
        veg = read(pft["PCT_NATVEG"], pi, slice(None), slice(None)) / 100.0
        tree_hr = (p[:8].sum(axis=0) * veg).reshape(-1)
        # NLCD "natural" = vegetated fraction minus its cropland part
        nat_hr = ((1.0 - p[8]) * veg).reshape(-1)

        nlcd_tree = to_coarse(tree_hr)
        nlcd_nat = to_coarse(nat_hr)

        for k, f in [
            ("luh_forest", luh_forest), ("luh_nat", luh_nat), ("nlcd_tree", nlcd_tree),
            ("nlcd_nat", nlcd_nat), ("luh_crop", luh_crop), ("luh_past", luh_past),
            ("harv", harv),
        ]:
            ts[k][iy] = (f * area_c_flat)[domain].sum()

        if year in (1850, 1950, 2014):
            snap[year] = {
                "luh_forest": luh_forest.copy(), "nlcd_tree": nlcd_tree.copy(),
            }
        if year % 25 == 0 or year == END_YEAR:
            print(
                f"{year}: LUH2 forest={ts['luh_forest'][iy]/1e6:6.3f}e6  "
                f"NLCD tree={ts['nlcd_tree'][iy]/1e6:6.3f}e6 km2  "
                f"ratio={ts['nlcd_tree'][iy]/ts['luh_forest'][iy]:.3f}"
            )

    pft.close(); st.close(); tr.close()

    hdr = ["year"] + keys
    np.savetxt(
        DIAG_DIR / "forest_luh2_vs_nlcd_km2.csv",
        np.column_stack([years] + [ts[k] for k in keys]),
        delimiter=",", header=",".join(hdr), comments="", fmt="%.6g",
    )

    # ---------------- figure 1: time series ----------------
    fig, ax = plt.subplots(3, 1, figsize=(11, 12), sharex=True)

    a = ax[0]
    a.plot(years, ts["luh_forest"] / 1e6, color="#1f77b4", lw=2, label="LUH2 forest (primf+secdf)")
    a.plot(years, ts["nlcd_tree"] / 1e6, color="#d62728", lw=2, label="NLCD tree cover (PFT 1-8)")
    a.plot(years, ts["luh_nat"] / 1e6, color="#1f77b4", lw=1.4, ls="--", label="LUH2 natural (prim+sec, f+n)")
    a.plot(years, ts["nlcd_nat"] / 1e6, color="#d62728", lw=1.4, ls="--", label="NLCD natveg excl. cropland")
    a.plot(years, ts["luh_crop"] / 1e6, color="#7f7f7f", lw=1.2, ls=":", label="LUH2 cropland")
    a.set_ylabel("area (10$^6$ km$^2$)")
    a.set_title(
        "Forest / natural land over the 1/24$\\degree$ domain (25-50$\\degree$N, 125-65$\\degree$W)\n"
        "NLCD aggregated to 0.25$\\degree$; LUH2 secdf is a land-use category, an upper bound on canopy",
        fontsize=11,
    )
    a.legend(fontsize=9); a.grid(alpha=0.3)

    a = ax[1]
    a.plot(years, ts["nlcd_tree"] / ts["luh_forest"], color="#d62728", lw=2, label="NLCD tree / LUH2 forest")
    a.plot(years, ts["nlcd_nat"] / ts["luh_nat"], color="#1f77b4", lw=2, label="NLCD natveg / LUH2 natural")
    a.axhline(1.0, color="k", ls=":", lw=1.5)
    a.set_ylabel("ratio")
    a.set_title("Denominator disagreement: <1 means NLCD has less land than LUH2 assumes", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    a = ax[2]
    a.plot(years, 100 * ts["harv"] / ts["luh_forest"], color="#1f77b4", lw=2,
           label="LUH2 harvest / LUH2 forest  (intended)")
    a.plot(years, 100 * ts["harv"] / ts["nlcd_tree"], color="#d62728", lw=2,
           label="LUH2 harvest / NLCD tree  (what ELM effectively applies)")
    a.set_xlabel("year"); a.set_ylabel("% of forest harvested per year")
    a.set_title("Implied harvest intensity (after the 36x area bug is fixed)", fontsize=11)
    a.legend(fontsize=9); a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(DIAG_DIR / "forest_luh2_vs_nlcd_timeseries.png", dpi=140, bbox_inches="tight")

    # ---------------- figure 2: maps + scatter at 2014 ----------------
    ml = (lat_c >= lat_hr.min()) & (lat_c <= lat_hr.max())
    mo = (lon_c >= lon_hr.min()) & (lon_c <= lon_hr.max())
    ext = [lon_c[mo].min(), lon_c[mo].max(), lat_c[ml].min(), lat_c[ml].max()]

    def box(flat):
        return flat.reshape(nlat_c, nlon_c)[np.ix_(ml, mo)]

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    s = snap[2014]
    for j, (arr, ttl) in enumerate(
        [(s["luh_forest"], "LUH2 primf+secdf, 2014"), (s["nlcd_tree"], "NLCD tree cover, 2014")]
    ):
        im = ax[0, j].imshow(box(arr), origin="lower", extent=ext, vmin=0, vmax=1, cmap="YlGn", aspect="auto")
        ax[0, j].set_title(ttl, fontsize=11)
        plt.colorbar(im, ax=ax[0, j], label="fraction of cell")

    d = box(s["nlcd_tree"]) - box(s["luh_forest"])
    im = ax[1, 0].imshow(d, origin="lower", extent=ext, vmin=-0.6, vmax=0.6, cmap="RdBu_r", aspect="auto")
    ax[1, 0].set_title("NLCD tree $-$ LUH2 forest, 2014", fontsize=11)
    plt.colorbar(im, ax=ax[1, 0], label="fraction of cell")

    x, y = s["luh_forest"][domain], s["nlcd_tree"][domain]
    h = ax[1, 1].hist2d(x, y, bins=60, range=[[0, 1], [0, 1]], cmap="magma", norm=matplotlib.colors.LogNorm())
    ax[1, 1].plot([0, 1], [0, 1], "w--", lw=1.5)
    ax[1, 1].set_xlabel("LUH2 primf+secdf"); ax[1, 1].set_ylabel("NLCD tree cover")
    ax[1, 1].set_title(f"per-cell, 2014 (r = {np.corrcoef(x, y)[0,1]:.3f})", fontsize=11)
    plt.colorbar(h[3], ax=ax[1, 1], label="cells")

    fig.tight_layout()
    fig.savefig(DIAG_DIR / "forest_luh2_vs_nlcd_maps.png", dpi=140, bbox_inches="tight")

    # ---------------- summary ----------------
    print("\n" + "=" * 70)
    for y0 in (1850, 1950, 2014):
        i = y0 - START_YEAR
        print(
            f"{y0}: LUH2 forest={ts['luh_forest'][i]/1e6:.3f}e6  NLCD tree={ts['nlcd_tree'][i]/1e6:.3f}e6 km2 "
            f"| ratio={ts['nlcd_tree'][i]/ts['luh_forest'][i]:.3f} "
            f"| natural ratio={ts['nlcd_nat'][i]/ts['luh_nat'][i]:.3f}"
        )
    i = -1
    print(
        f"\n2014 harvest intensity: intended={100*ts['harv'][i]/ts['luh_forest'][i]:.4f}%/yr  "
        f"applied={100*ts['harv'][i]/ts['nlcd_tree'][i]:.4f}%/yr  "
        f"bias={ts['luh_forest'][i]/ts['nlcd_tree'][i]:.3f}x"
    )
    x, y = snap[2014]["luh_forest"][domain], snap[2014]["nlcd_tree"][domain]
    print(f"per-cell 2014: r={np.corrcoef(x,y)[0,1]:.3f}  "
          f"mean|diff|={np.abs(y-x).mean():.4f}  bias(NLCD-LUH2)={np.mean(y-x):+.4f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
