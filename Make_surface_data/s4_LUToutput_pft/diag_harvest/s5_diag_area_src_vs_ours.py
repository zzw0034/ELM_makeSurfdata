"""Compare vegetation-type AREA (km2) in the source (s4_1 NLCD) vs our final ELM product.

Source (s4_1, NA 1/24 grid), restricted to the SEUS box:
   area_type = Σ [ (Σ PCT_NAT_PFT over type)/100 * PCT_NATVEG/100 * cellarea ]
   (PCT_NATVEG here = NLCD vegetated fraction of the gridcell, time-varying)

Ours (final product on SEUS 1/24):
   area_type = Σ [ (Σ PCT_NAT_PFT over type)/100 * PCT_NATVEG_landunit/100 * LANDFRAC_PFT * cellarea ]
   PCT_NAT_PFT from landuse.timeseries (time-varying);
   PCT_NATVEG_landunit + LANDFRAC_PFT from surfdata (STATIC).

This isolates how much the pipeline (regrid + static landunit substitution + landfrac)
changes the actual areas, beyond the composition (which we already showed matches).
"""
from __future__ import annotations
import numpy as np
import netCDF4 as nc

SRC = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"
LU = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/surfdata_results/landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260723.nc"
SURF = "/projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map/surfdata_SEUS_1_24deg_simyr1850_c260723.nc"
SEUS = dict(lon0=-95.0, lon1=-74.0, lat0=24.0, lat1=37.5)
R = 6371.0
TREE = list(range(1, 9)); GRASS = [13, 14]; CROP = [15]
YEARS = {1850: 0, 2023: 173}


def sl(ds, v, *idx):
    a = ds[v][idx] if idx else ds[v][:]
    return np.asarray(np.ma.filled(np.ma.masked_invalid(a), 0.0), dtype=np.float64)


def area_km2(lat, lon):
    dlat = abs(np.median(np.diff(lat))); dlon = abs(np.median(np.diff(lon)))
    n = np.deg2rad(lat + 0.5 * dlat); s = np.deg2rad(lat - 0.5 * dlat)
    band = (R ** 2) * np.deg2rad(dlon) * (np.sin(n) - np.sin(s))
    return np.repeat(band[:, None], lon.size, axis=1)


def main() -> None:
    # ---- source (NA), restrict to SEUS ----
    s = nc.Dataset(SRC)
    slat = sl(s, "lat"); slon = sl(s, "lon")
    sA = area_km2(slat, slon)
    smask = ((slat >= SEUS["lat0"]) & (slat <= SEUS["lat1"]))[:, None] & (
        (slon >= SEUS["lon0"]) & (slon <= SEUS["lon1"]))[None, :]

    # ---- ours ----
    o = nc.Dataset(LU); sf = nc.Dataset(SURF)
    olat = sl(sf, "LATIXY")[:, 0]; olon = sl(sf, "LONGXY")[0, :]
    oA = area_km2(olat, olon)
    nv_lu = sl(sf, "PCT_NATVEG") / 100.0          # static landunit fraction
    lfrac = sl(sf, "LANDFRAC_PFT")                # static land fraction

    def type_area_src(t, idxs):
        pft = sl(s, "PCT_NAT_PFT", t, slice(None), slice(None), slice(None))  # (17,lat,lon)
        comp = pft[idxs].sum(0) / 100.0
        nv = sl(s, "PCT_NATVEG", t, slice(None), slice(None)) / 100.0
        return float((comp * nv * sA * smask).sum())

    def type_area_ours(t, idxs):
        pft = sl(o, "PCT_NAT_PFT", t, slice(None), slice(None), slice(None))  # (17,lat,lon)
        comp = pft[idxs].sum(0) / 100.0
        return float((comp * nv_lu * lfrac * oA).sum())

    print(f"{'':10s}{'type':6s}{'source km2':>14s}{'ours km2':>14s}{'ours/src':>10s}")
    for y, t in YEARS.items():
        for name, idxs in [("tree", TREE), ("grass", GRASS), ("crop", CROP)]:
            asrc = type_area_src(t, idxs); aours = type_area_ours(t, idxs)
            r = aours / asrc if asrc > 0 else float("nan")
            print(f"{y:<10d}{name:6s}{asrc:14.1f}{aours:14.1f}{r:10.3f}")
    # total vegetated (natveg) area, source vs ours
    print("\n--- total natveg area (all PFT) ---")
    for y, t in YEARS.items():
        allidx = list(range(17))
        asrc = type_area_src(t, allidx); aours = type_area_ours(t, allidx)
        print(f"  {y}: source={asrc:.1f}  ours={aours:.1f}  ours/src={aours/asrc:.3f}")
    s.close(); o.close(); sf.close()


if __name__ == "__main__":
    main()
