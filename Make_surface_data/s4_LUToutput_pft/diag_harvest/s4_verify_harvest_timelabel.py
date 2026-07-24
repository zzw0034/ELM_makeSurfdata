"""Verify the harvest year-labelling convention against the reference LUT files.

E3SM doc p3 claims: the file labelled year L carries the harvest that occurred in
model year L-1 (i.e. LUH2 calendar year L-1). s4_2 currently writes LUH2 year L into
the file labelled L (no shift). Before changing s4_2 we confirm which is right by
comparing the *reference* LUT files (built by the real LUT code) against our LUH2
transitions.nc.

Method (grid-independent): the reference LUT harvest = LUH2 harvest / veg_frac, so the
two global-mean time series share temporal shape but differ in magnitude. Cross-
correlate them (raw and first-differenced) and report the lag k that maximises
correlation, where series are aligned as  LUT(label=Y) <-> LUH2(cal=Y+k).
  k = 0  -> label year == LUH2 data year  (current s4_2 is correct)
  k = -1 -> LUT(label Y) holds LUH2 year Y-1 (doc is correct; shift needed)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF_DIR = Path(
    "/projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/rawdata/LUT_LUH2_historical_04082019"
)
REF_GLOB = "LUT_LUH2_historical_{year}_04082019.nc"
LUH2_PATH = Path("/projects/hpcl-cli185/proj-shared/zw5/luh/transitions.nc")
DIAG_DIR = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest"
)

Y0, Y1 = 1850, 2014
LUH2_YEAR0 = 850
LUT_HARV = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3"]
LUH_HARV = ["primf_harv", "primn_harv", "secmf_harv", "secyf_harv", "secnf_harv"]


def read(var, *idx):
    arr = var[idx] if idx else var[:]
    return np.asarray(np.ma.filled(np.ma.masked_invalid(arr), 0.0), dtype=np.float64)


def area_weights(lat):
    w = np.cos(np.deg2rad(lat))
    return w / w.sum()


def main() -> None:
    years = np.arange(Y0, Y1 + 1)

    # --- reference LUT global-mean harvest per label year ---
    lut_series = np.full(years.size, np.nan)
    lat_w_ref = None
    for iy, y in enumerate(years):
        f = REF_DIR / REF_GLOB.format(year=y)
        if not f.exists():
            continue
        d = nc.Dataset(f, "r")
        latname = "lat" if "lat" in d.variables else ("LAT" if "LAT" in d.variables else None)
        lat = read(d[latname]) if latname else np.linspace(-89.75, 89.75, d.dimensions["lat"].size)
        if lat_w_ref is None:
            lat_w_ref = area_weights(lat)
        tot = sum(read(d[v]) for v in LUT_HARV)  # (lat, lon)
        lut_series[iy] = (tot.mean(axis=1) * lat_w_ref).sum()
        d.close()

    # --- LUH2 global-mean harvest per calendar year (Y0-1 .. Y1) ---
    lu = nc.Dataset(LUH2_PATH, "r")
    lat_luh = read(lu["lat"])
    lat_w_luh = area_weights(lat_luh)
    cal_years = np.arange(Y0 - 2, Y1 + 1)
    luh_by_cal = {}
    for y in cal_years:
        ti = y - LUH2_YEAR0
        tot = sum(read(lu[v], ti, slice(None), slice(None)) for v in LUH_HARV)
        luh_by_cal[y] = (tot.mean(axis=1) * lat_w_luh).sum()
    lu.close()

    def z(a):
        a = np.asarray(a, float)
        m = np.isfinite(a)
        return (a - a[m].mean()) / a[m].std(), m

    # align LUT(label=Y) with LUH2(cal=Y+k) for several lags
    print("lag k :  corr(raw)   corr(diff)   [LUT(label Y) vs LUH2(cal Y+k)]")
    results = {}
    for k in (-2, -1, 0, 1, 2):
        luh_aligned = np.array([luh_by_cal.get(y + k, np.nan) for y in years])
        za, ma = z(lut_series)
        zb, mb = z(luh_aligned)
        m = ma & mb
        c_raw = np.corrcoef(lut_series[m], luh_aligned[m])[0, 1]
        # differenced series (sharpens temporal features)
        da = np.diff(lut_series)
        db = np.diff(luh_aligned)
        md = np.isfinite(da) & np.isfinite(db)
        c_diff = np.corrcoef(da[md], db[md])[0, 1]
        results[k] = (c_raw, c_diff)
        print(f"  k={k:+d} :   {c_raw:.5f}     {c_diff:.5f}")

    best_raw = max(results, key=lambda k: results[k][0])
    best_diff = max(results, key=lambda k: results[k][1])
    print(f"\nbest lag by raw corr : k={best_raw:+d}")
    print(f"best lag by diff corr: k={best_diff:+d}")
    print("interpretation: k=0 -> current s4_2 correct; k=-1 -> LUT(label Y)=LUH2(Y-1), doc correct, shift needed")

    # --- figure: normalized overlay at the best lag vs k=0 ---
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    za, _ = z(lut_series)
    for axi, k in zip(ax, [0, best_diff]):
        luh_aligned = np.array([luh_by_cal.get(y + k, np.nan) for y in years])
        zb, _ = z(luh_aligned)
        axi.plot(years, za, color="#d62728", lw=2, label="reference LUT (label year), z-score")
        axi.plot(years, zb, color="#1f77b4", lw=1.6, ls="--",
                 label=f"LUH2 (cal year {'Y' if k==0 else f'Y{k:+d}'}), z-score")
        axi.set_title(f"lag k={k:+d}   raw r={results[k][0]:.4f}, diff r={results[k][1]:.4f}", fontsize=11)
        axi.legend(fontsize=9); axi.grid(alpha=0.3); axi.set_ylabel("z-score")
    ax[1].set_xlabel("year")
    fig.suptitle("Harvest year-label verification: reference LUT vs LUH2 transitions", fontsize=12)
    fig.tight_layout()
    fig.savefig(DIAG_DIR / "verify_harvest_timelabel.png", dpi=140, bbox_inches="tight")

    # --- a couple of explicit spot checks around distinctive years ---
    print("\nspot check (global-mean harvest):")
    print("  year Y | refLUT(label Y) | LUH2(Y) | LUH2(Y-1) | LUH2(Y+1)")
    for y in (1900, 1945, 1950, 2000, 2014):
        iy = y - Y0
        print(f"   {y}  |  {lut_series[iy]:.4e}  | {luh_by_cal.get(y):.4e} | "
              f"{luh_by_cal.get(y-1):.4e} | {luh_by_cal.get(y+1, np.nan):.4e}")


if __name__ == "__main__":
    main()
