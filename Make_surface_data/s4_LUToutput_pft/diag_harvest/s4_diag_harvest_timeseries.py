"""Compare downscaled harvest (s4_2 output) against the LUH2 source, 1850-2023.

For every year this computes the *area-weighted* harvested area (km2) inside the
1/24 deg domain, from three angles:

  luh2_total  : LUH2 harvest summed over all coarse cells that overlap the domain.
  luh2_alloc  : same, but only over coarse cells that have non-zero tree weight,
                i.e. the part s4_2 is even able to place somewhere.
  ours        : the s4_2 output, converted back to grid fraction (HARVEST * PCT_NATVEG)
                and integrated over true cell areas.

luh2_total / ours therefore separates into two independent effects:
  luh2_total -> luh2_alloc : harvest dropped in treeless coarse cells (by design).
  luh2_alloc -> ours       : the missing coarse/fine area ratio in
                             _distribute_conservatively (sum vs area-weighted mean).

Writes a CSV of the yearly numbers plus a 3-panel figure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -------------------------
# Paths (absolute: SLURM runs from a spool copy)
# -------------------------
BASE = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft"
)
PFT_PATH = BASE / "scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
LUH2_PATH = Path("/projects/hpcl-cli185/proj-shared/zw5/luh/transitions.nc")
LUH2_SSP_PATH = Path("/projects/hpcl-cli185/proj-shared/zw5/luh/ssp2rcp45_transitions.nc")
OUT_GLOB = "output_downscaled_luh2_harvest/LUT_nlcd2elm_luh2_historical_{year}_*.nc"
DIAG_DIR = BASE / "diag_harvest"

START_YEAR = 1850
END_YEAR = 2023
LUH2_YEAR0 = 850
LUH2_LAST_YEAR = 2014
LUH2_SSP_YEAR0 = 2015
LUH2_SSP_LAST_YEAR = 2099
PFT_YEAR0 = 1850

R_EARTH = 6371.0  # km
TREE_PFT_IDXS = np.arange(1, 9)  # same as s4_2 TREE_PFT_IDXS

HARVEST_VAR_MAP = {
    "primf_harv": "HARVEST_VH1",
    "primn_harv": "HARVEST_VH2",
    "secmf_harv": "HARVEST_SH1",
    "secyf_harv": "HARVEST_SH2",
    "secnf_harv": "HARVEST_SH3",
}


def read(var, *idx) -> np.ndarray:
    """Read a netCDF variable as float64 with masked/NaN values turned into 0.

    netCDF4 hands back a MaskedArray when _FillValue is set; np.asarray on that
    silently returns the raw fill (1e20 for LUH2), so fill it explicitly.
    """
    arr = var[idx] if idx else var[:]
    arr = np.ma.filled(np.ma.masked_invalid(arr), 0.0)
    return np.asarray(arr, dtype=np.float64)


def cell_area(lat_1d: np.ndarray, dlat: float, dlon: float, nlon: int) -> np.ndarray:
    """Exact spherical cell area (km2), shape (nlat, nlon)."""
    lat_n = np.deg2rad(lat_1d + 0.5 * abs(dlat))
    lat_s = np.deg2rad(lat_1d - 0.5 * abs(dlat))
    band = R_EARTH**2 * np.deg2rad(abs(dlon)) * (np.sin(lat_n) - np.sin(lat_s))
    return np.repeat(band[:, None], nlon, axis=1)


def build_hr_to_coarse_index(lat_hr, lon_hr, lat_coarse, lon_coarse):
    """Identical mapping to s4_2 _build_hr_to_coarse_index."""
    nlat_c, nlon_c = lat_coarse.size, lon_coarse.size
    dlat = float(np.median(np.diff(lat_coarse)))
    dlon = float(np.median(np.diff(lon_coarse)))
    lat0 = float(lat_coarse[0] - 0.5 * dlat)
    lon0 = float(lon_coarse[0] - 0.5 * dlon)

    lat_idx = np.floor((lat_hr - lat0) / dlat).astype(np.int64)
    lon_idx = np.floor((lon_hr - lon0) / dlon).astype(np.int64)
    lat_ok = (lat_idx >= 0) & (lat_idx < nlat_c)
    lon_ok = (lon_idx >= 0) & (lon_idx < nlon_c)

    inside_2d = lat_ok[:, None] & lon_ok[None, :]
    coarse_id_2d = (lat_idx[:, None] * nlon_c + lon_idx[None, :]).astype(np.int64)
    coarse_id_2d[~inside_2d] = 0
    return coarse_id_2d, inside_2d, nlat_c * nlon_c


def main() -> None:
    DIAG_DIR.mkdir(parents=True, exist_ok=True)

    pft = nc.Dataset(PFT_PATH, "r")
    luh = nc.Dataset(LUH2_PATH, "r")
    luh_ssp = nc.Dataset(LUH2_SSP_PATH, "r")  # v2f SSP2-4.5 for cal year >= 2015

    lat_hr = np.asarray(pft["lat"][:], dtype=np.float64)
    lon_hr = np.asarray(pft["lon"][:], dtype=np.float64)
    dlat_hr = float(np.median(np.diff(lat_hr)))
    dlon_hr = float(np.median(np.diff(lon_hr)))
    area_hr = cell_area(lat_hr, dlat_hr, dlon_hr, lon_hr.size)

    # s4_2 now uses the ANNUAL vegetated fraction (read per year inside the loop),
    # so "ours" is reconstructed with the same annual divisor the production run used.

    # LUH2 lat is descending; sort ascending as s4_2 does.
    lat_luh_raw = np.asarray(luh["lat"][:], dtype=np.float64)
    lon_luh = np.asarray(luh["lon"][:], dtype=np.float64)
    lat_sort = np.argsort(lat_luh_raw)
    lat_luh = lat_luh_raw[lat_sort]
    dlat_luh = float(np.median(np.diff(lat_luh)))
    dlon_luh = float(np.median(np.diff(lon_luh)))
    area_luh = cell_area(lat_luh, dlat_luh, dlon_luh, lon_luh.size)
    area_luh_flat = area_luh.reshape(-1)

    coarse_id_2d, inside_2d, ncoarse = build_hr_to_coarse_index(
        lat_hr, lon_hr, lat_luh, lon_luh
    )
    id_flat = coarse_id_2d.reshape(-1)
    inside_flat = inside_2d.reshape(-1)

    # Domain = coarse cells that actually receive at least one hi-res cell.
    # Using this instead of a lat/lon box removes any edge ambiguity.
    n_hr_per_coarse = np.bincount(
        id_flat, weights=inside_flat.astype(np.float64), minlength=ncoarse
    )
    domain = n_hr_per_coarse > 0
    print(f"hi-res grid : {lat_hr.size} x {lon_hr.size} @ {dlat_hr:.6f} deg")
    print(f"LUH2 grid   : {lat_luh.size} x {lon_luh.size} @ {dlat_luh:.6f} deg")
    print(f"coarse cells overlapping domain: {int(domain.sum())}")
    print(f"hi-res cells per coarse cell: {n_hr_per_coarse[domain].mean():.2f} (mean)")

    area_hr_flat = area_hr.reshape(-1) * inside_flat
    years = np.arange(START_YEAR, END_YEAR + 1)
    nvar = len(HARVEST_VAR_MAP)
    luh2_total = np.zeros((years.size, nvar))
    luh2_alloc = np.zeros((years.size, nvar))
    ours = np.zeros((years.size, nvar))

    for iy, year in enumerate(years):
        pft_idx = year - PFT_YEAR0
        # File labelled `year` carries harvest of the prior calendar year (s4_2 P3 fix),
        # spliced v2h (<=2014) -> v2f SSP2-4.5 (>=2015), matching s4_2.
        cal_year = year - 1
        if cal_year <= LUH2_LAST_YEAR:
            harv_ds, luh_idx = luh, cal_year - LUH2_YEAR0
        else:
            harv_ds = luh_ssp
            luh_idx = min(cal_year, LUH2_SSP_LAST_YEAR) - LUH2_SSP_YEAR0

        # Annual vegetated fraction (s4_2 no longer freezes it at 1850).
        veg_frac = read(pft["PCT_NATVEG"], pft_idx, slice(None), slice(None))
        if veg_frac.max() > 1.5:
            veg_frac = veg_frac / 100.0

        # Tree weight, rebuilt exactly as s4_2 does it.
        tree = read(
            pft["PCT_NAT_PFT"], pft_idx, TREE_PFT_IDXS, slice(None), slice(None)
        ).sum(axis=0)
        if tree.max() > 1.5:
            tree = tree / 100.0
        w_flat = np.clip(tree, 0.0, 1.0).reshape(-1) * veg_frac.reshape(-1) * inside_flat
        sum_w = np.bincount(id_flat, weights=w_flat, minlength=ncoarse)
        has_tree = sum_w > 0.0

        matches = sorted(BASE.glob(OUT_GLOB.format(year=year)))
        if len(matches) != 1:
            raise RuntimeError(f"expected 1 output file for {year}, got {len(matches)}")
        out = nc.Dataset(matches[0], "r")

        for iv, (luh_name, out_name) in enumerate(HARVEST_VAR_MAP.items()):
            coarse = read(harv_ds[luh_name], luh_idx, slice(None), slice(None))[
                lat_sort, :
            ].reshape(-1)
            harv_area = coarse * area_luh_flat
            luh2_total[iy, iv] = harv_area[domain].sum()
            luh2_alloc[iy, iv] = harv_area[domain & has_tree].sum()

            # LUT variable is a fraction of the vegetated unit -> back to grid fraction.
            lut = read(out[out_name])
            ours[iy, iv] = ((lut * veg_frac).reshape(-1) * area_hr_flat).sum()

        out.close()
        if year % 25 == 0 or year == END_YEAR:
            print(
                f"{year}: luh2_total={luh2_total[iy].sum():10.1f} "
                f"luh2_alloc={luh2_alloc[iy].sum():10.1f} "
                f"ours={ours[iy].sum():8.1f} km2"
            )

    pft.close()
    luh.close()
    luh_ssp.close()

    # ---- CSV ----
    names = list(HARVEST_VAR_MAP.values())
    header = ["year", "luh2_total", "luh2_alloc", "ours"]
    header += [f"luh2_{n}" for n in names] + [f"ours_{n}" for n in names]
    table = np.column_stack(
        [years, luh2_total.sum(1), luh2_alloc.sum(1), ours.sum(1), luh2_total, ours]
    )
    csv_path = DIAG_DIR / "harvest_timeseries_km2.csv"
    np.savetxt(
        csv_path, table, delimiter=",", header=",".join(header), comments="", fmt="%.6g"
    )
    print(f"wrote {csv_path}")

    # ---- figure ----
    tot_l, tot_a, tot_o = luh2_total.sum(1), luh2_alloc.sum(1), ours.sum(1)
    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)

    ax = axes[0]
    ax.plot(years, tot_l, color="#1f77b4", lw=2, label="LUH2 source (all coarse cells)")
    ax.plot(
        years, tot_a, color="#2ca02c", lw=1.6, ls="--",
        label="LUH2 restricted to cells with trees (allocatable)",
    )
    ax.plot(years, tot_o, color="#d62728", lw=2, label="s4_2 downscaled output")
    ax.set_yscale("log")
    ax.set_ylabel("harvested area (km$^2$ yr$^{-1}$)")
    ax.set_title(
        "Total wood harvest over the 1/24$\\degree$ domain (25-50$\\degree$N, 125-65$\\degree$W)\n"
        "harvest: v2h historical (<=2014) spliced to v2f SSP2-4.5 (>=2015)",
        fontsize=11,
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        r_tot = np.where(tot_o > 0, tot_l / tot_o, np.nan)
        r_all = np.where(tot_o > 0, tot_a / tot_o, np.nan)
    ax.plot(years, r_tot, color="#1f77b4", lw=2, label="LUH2 total / ours")
    ax.plot(years, r_all, color="#2ca02c", lw=1.6, ls="--", label="LUH2 allocatable / ours")
    ax.axhline(
        n_hr_per_coarse[domain].mean(),
        color="k", ls=":", lw=1.5,
        label=f"mean hi-res cells per coarse cell = {n_hr_per_coarse[domain].mean():.0f}",
    )
    ax.set_ylabel("ratio")
    ax.set_title("Deficit factor: green ~ the area-ratio bug alone, blue adds dropped harvest", fontsize=11)
    ax.legend(loc="best", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[2]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for iv, n in enumerate(names):
        if luh2_total[:, iv].max() <= 0 and ours[:, iv].max() <= 0:
            continue
        ax.plot(years, luh2_total[:, iv], color=colors[iv], lw=2, label=f"{n} LUH2")
        ax.plot(years, ours[:, iv], color=colors[iv], lw=1.6, ls="--", label=f"{n} ours")
    ax.set_yscale("log")
    ax.set_xlabel("year")
    ax.set_ylabel("harvested area (km$^2$ yr$^{-1}$)")
    zeros = [n for iv, n in enumerate(names) if luh2_total[:, iv].max() <= 0]
    ax.set_title(
        "By harvest type (solid = LUH2, dashed = ours)"
        + (f"  |  identically zero in LUH2 over this domain: {', '.join(zeros)}" if zeros else ""),
        fontsize=11,
    )
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig_path = DIAG_DIR / "harvest_timeseries_luh2_vs_downscaled.png"
    fig.savefig(fig_path, dpi=140, bbox_inches="tight")
    print(f"wrote {fig_path}")

    # ---- summary ----
    print("\n" + "=" * 64)
    print(f"{'':16s}{'1850-2023 total (km2)':>24s}")
    print(f"{'LUH2 source':16s}{tot_l.sum():24.1f}")
    print(f"{'LUH2 allocatable':16s}{tot_a.sum():24.1f}")
    print(f"{'s4_2 output':16s}{tot_o.sum():24.1f}")
    print(f"\ndeficit LUH2_total/ours   = {tot_l.sum()/tot_o.sum():.2f}")
    print(f"deficit LUH2_alloc/ours   = {tot_a.sum()/tot_o.sum():.2f}  <- pure area-ratio bug")
    print(f"dropped in treeless cells = {100*(1-tot_a.sum()/tot_l.sum()):.2f}% of LUH2")
    print("=" * 64)


if __name__ == "__main__":
    main()
