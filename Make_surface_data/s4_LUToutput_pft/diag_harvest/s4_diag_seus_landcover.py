"""Full land-cover comparison in the SEUS box: LUH2 states vs the NLCD-derived product.

Compares every major cover class, not just forest, on both spatial pattern and
interannual behaviour. NLCD is aggregated to 0.25 deg (area-weighted); LUH2 is
never downscaled.

Class pairing (LUH2 is land *use*, the NLCD product is land *cover*, so the pairing
is approximate by construction):

  tree       primf + secdf                      <->  PFT 1-8
  crop       c3ann+c4ann+c3per+c4per+c3nfx      <->  PFT 15,16
  open       primn + secdn + pastr + range      <->  PFT 9-14 (shrub + grass)
  (bare)     no LUH2 counterpart                     PFT 0
  (urban)    urban                                   PCT_URBAN, frozen at 1850

Note the columns do not have to close: LUH2 states sum to 1 over land, while the
NLCD classes sum to PCT_NATVEG, which excludes water, ice and developed land.

Outputs: area time series, year-to-year increments (temporal texture), per-year
spatial correlation, and 2014 maps.
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

GRP = {"tree": [1, 2, 3, 4, 5, 6, 7, 8], "open": [9, 10, 11, 12, 13, 14],
       "crop": [15, 16], "bare": [0]}
LUH = {"tree": ["primf", "secdf"], "open": ["primn", "secdn", "pastr", "range"],
       "crop": ["c3ann", "c4ann", "c3per", "c4per", "c3nfx"], "urban": ["urban"]}
CLASSES = ["tree", "open", "crop"]


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

    lat_hr = np.asarray(pft["lat"][:], float)
    lon_hr = np.asarray(pft["lon"][:], float)
    area_hr = cell_area(lat_hr, np.median(np.diff(lat_hr)), np.median(np.diff(lon_hr)), lon_hr.size)
    box = ((lat_hr >= SEUS["lat0"]) & (lat_hr <= SEUS["lat1"]))[:, None] & (
        (lon_hr >= SEUS["lon0"]) & (lon_hr <= SEUS["lon1"]))[None, :]

    lat_raw = np.asarray(st["lat"][:], float)
    lon_c = np.asarray(st["lon"][:], float)
    srt = np.argsort(lat_raw)
    lat_c = lat_raw[srt]
    nlat_c, nlon_c = lat_c.size, lon_c.size
    area_c = cell_area(lat_c, np.median(np.diff(lat_c)), np.median(np.diff(lon_c)), nlon_c).reshape(-1)

    dlat_c, dlon_c = np.median(np.diff(lat_c)), np.median(np.diff(lon_c))
    li = np.floor((lat_hr - (lat_c[0] - 0.5 * dlat_c)) / dlat_c).astype(np.int64)
    oi = np.floor((lon_hr - (lon_c[0] - 0.5 * dlon_c)) / dlon_c).astype(np.int64)
    inside = ((li >= 0) & (li < nlat_c))[:, None] & ((oi >= 0) & (oi < nlon_c))[None, :] & box
    cid = (li[:, None] * nlon_c + oi[None, :]).astype(np.int64)
    cid[~inside] = 0
    id_flat = cid.reshape(-1)
    a_hr_flat = area_hr.reshape(-1) * inside.reshape(-1)
    a_in = np.bincount(id_flat, weights=a_hr_flat, minlength=nlat_c * nlon_c)
    domain = a_in > 0
    print(f"SEUS coarse cells: {int(domain.sum())}")

    def agg(f):
        num = np.bincount(id_flat, weights=f * a_hr_flat, minlength=nlat_c * nlon_c)
        o = np.zeros_like(num)
        np.divide(num, a_in, out=o, where=a_in > 0)
        return o

    years = np.arange(START_YEAR, END_YEAR + 1)
    A = {f"{s}_{c}": np.zeros(years.size) for s in ("luh", "nlcd") for c in CLASSES + ["bare"]}
    A["luh_urban"] = np.zeros(years.size)
    A["nlcd_natveg"] = np.zeros(years.size)
    corr = {c: np.zeros(years.size) for c in CLASSES}
    mad = {c: np.zeros(years.size) for c in CLASSES}
    snap = {}

    for iy, year in enumerate(years):
        ti, pi = year - LUH2_YEAR0, year - PFT_YEAR0
        L, N = {}, {}
        for c in CLASSES:
            L[c] = sum(read(st[v], ti) for v in LUH[c])[srt, :].reshape(-1)
        A["luh_urban"][iy] = (read(st["urban"], ti)[srt, :].reshape(-1) * area_c)[domain].sum()

        p = read(pft["PCT_NAT_PFT"], pi, slice(None), slice(None), slice(None)) / 100.0
        veg = read(pft["PCT_NATVEG"], pi, slice(None), slice(None)) / 100.0
        for c in CLASSES + ["bare"]:
            N[c] = agg((p[GRP[c]].sum(axis=0) * veg).reshape(-1))
        A["nlcd_natveg"][iy] = (agg(veg.reshape(-1)) * area_c)[domain].sum()

        for c in CLASSES:
            A[f"luh_{c}"][iy] = (L[c] * area_c)[domain].sum()
            A[f"nlcd_{c}"][iy] = (N[c] * area_c)[domain].sum()
            x, y = L[c][domain], N[c][domain]
            corr[c][iy] = np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan
            mad[c][iy] = np.abs(x - y).mean()
        A["nlcd_bare"][iy] = (N["bare"] * area_c)[domain].sum()

        if year in (1850, 2014):
            snap[year] = {f"luh_{c}": L[c].copy() for c in CLASSES}
            snap[year].update({f"nlcd_{c}": N[c].copy() for c in CLASSES})

        if year % 25 == 0 or year == END_YEAR:
            print(f"{year}: " + "  ".join(
                f"{c}: L={A[f'luh_{c}'][iy]/1e3:6.1f} N={A[f'nlcd_{c}'][iy]/1e3:6.1f} r={corr[c][iy]:.2f}"
                for c in CLASSES))

    pft.close(); st.close()

    col = {"tree": "#2ca02c", "open": "#ff7f0e", "crop": "#8c564b"}

    # ---- figure A: areas, increments, spatial correlation ----
    fig, ax = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    a = ax[0]
    for c in CLASSES:
        a.plot(years, A[f"luh_{c}"] / 1e3, color=col[c], lw=2, label=f"LUH2 {c}")
        a.plot(years, A[f"nlcd_{c}"] / 1e3, color=col[c], lw=1.6, ls="--", label=f"NLCD {c}")
    a.set_ylabel("area (10$^3$ km$^2$)")
    a.set_title("SEUS land cover: LUH2 states (solid) vs NLCD-derived product (dashed)", fontsize=11)
    a.legend(fontsize=8, ncol=3); a.grid(alpha=0.3)

    a = ax[1]
    for c in CLASSES:
        a.plot(years[1:], np.diff(A[f"luh_{c}"]) / 1e3, color=col[c], lw=1.6, label=f"LUH2 {c}")
        a.plot(years[1:], np.diff(A[f"nlcd_{c}"]) / 1e3, color=col[c], lw=1.2, ls="--", label=f"NLCD {c}")
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("year-to-year change (10$^3$ km$^2$)")
    a.set_title("Temporal texture: how each product distributes change across years", fontsize=11)
    a.legend(fontsize=8, ncol=3); a.grid(alpha=0.3)

    a = ax[2]
    for c in CLASSES:
        a.plot(years, corr[c], color=col[c], lw=2, label=f"{c}  (r)")
        a.plot(years, mad[c], color=col[c], lw=1.2, ls=":", label=f"{c}  (mean|diff|)")
    a.set_xlabel("year"); a.set_ylabel("r  /  mean |diff| (fraction)")
    a.set_ylim(0, 1)
    a.set_title("Per-year spatial agreement across the 0.25$\\degree$ cells", fontsize=11)
    a.legend(fontsize=8, ncol=3); a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "seus_landcover_timeseries.png", dpi=140, bbox_inches="tight")

    # ---- figure B: 2014 maps ----
    ml = (lat_c >= max(lat_hr.min(), SEUS["lat0"])) & (lat_c <= min(lat_hr.max(), SEUS["lat1"]))
    mo = (lon_c >= SEUS["lon0"]) & (lon_c <= SEUS["lon1"])
    ext = [lon_c[mo].min(), lon_c[mo].max(), lat_c[ml].min(), lat_c[ml].max()]
    dom2 = domain.reshape(nlat_c, nlon_c)[np.ix_(ml, mo)]

    def box2(flat):
        return np.where(dom2, flat.reshape(nlat_c, nlon_c)[np.ix_(ml, mo)], np.nan)

    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    s = snap[2014]
    for i, c in enumerate(CLASSES):
        for j, (k, ttl) in enumerate([(f"luh_{c}", "LUH2"), (f"nlcd_{c}", "NLCD")]):
            im = ax[i, j].imshow(box2(s[k]), origin="lower", extent=ext, vmin=0, vmax=1,
                                 cmap="YlGn", aspect="auto")
            ax[i, j].set_title(f"{ttl} {c}, 2014", fontsize=10)
            plt.colorbar(im, ax=ax[i, j])
        d = box2(s[f"nlcd_{c}"]) - box2(s[f"luh_{c}"])
        im = ax[i, 2].imshow(d, origin="lower", extent=ext, vmin=-0.6, vmax=0.6,
                             cmap="RdBu_r", aspect="auto")
        r = corr[c][-1]
        ax[i, 2].set_title(f"NLCD $-$ LUH2 {c}  (r={r:.2f})", fontsize=10)
        plt.colorbar(im, ax=ax[i, 2])
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "seus_landcover_maps.png", dpi=140, bbox_inches="tight")

    hdr = ["year"] + sorted(A) + [f"r_{c}" for c in CLASSES]
    np.savetxt(DIAG_DIR / "seus_landcover_km2.csv",
               np.column_stack([years] + [A[k] for k in sorted(A)] + [corr[c] for c in CLASSES]),
               delimiter=",", header=",".join(hdr), comments="", fmt="%.6g")

    print("\n" + "=" * 78)
    for y0 in (1850, 1950, 2014):
        i = y0 - START_YEAR
        print(f"--- {y0} ---")
        for c in CLASSES:
            l, n = A[f"luh_{c}"][i], A[f"nlcd_{c}"][i]
            print(f"  {c:6s} LUH2={l/1e3:7.1f}  NLCD={n/1e3:7.1f} e3 km2  "
                  f"N/L={n/l if l>0 else np.nan:5.3f}  r={corr[c][i]:.3f}  mean|diff|={mad[c][i]:.4f}")
    print(f"\nNLCD bare (no LUH2 counterpart), 2014: {A['nlcd_bare'][-1]/1e3:.1f}e3 km2")
    print(f"LUH2 urban 2014: {A['luh_urban'][-1]/1e3:.1f}e3 km2   "
          f"(NLCD urban is frozen at 1850 and excluded from PCT_NATVEG)")
    print(f"NLCD natveg total 2014: {A['nlcd_natveg'][-1]/1e3:.1f}e3 km2")
    for c in CLASSES:
        d = np.abs(np.diff(A[f"luh_{c}"])).mean() / 1e3, np.abs(np.diff(A[f"nlcd_{c}"])).mean() / 1e3
        print(f"mean |year-to-year change| {c:6s}: LUH2={d[0]:.3f}  NLCD={d[1]:.3f} e3 km2/yr")
    print("=" * 78)


if __name__ == "__main__":
    main()
