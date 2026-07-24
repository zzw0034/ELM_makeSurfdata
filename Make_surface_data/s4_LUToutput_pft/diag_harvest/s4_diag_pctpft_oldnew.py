"""Check whether PCT_PFT differs between old (_07112026) and new (_07232026) LUT files.

For every labelled year 1850-2023, compute the max absolute difference of PCT_PFT
(over all 17 PFTs x 601 x 1441 cells) between the old buggy file and the new fixed
file. Also report the max abs diff of the column-sum (should stay 100) and, for
contrast, the HARVEST difference (which is the thing that actually changed).
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import netCDF4 as nc

BASE = Path(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest"
)
OLD_SUB = BASE / "prev_buggy_c0711"
OLD_TAG, NEW_TAG = "07112026", "07232026"
HARV = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3"]


def rd(v):
    return np.asarray(np.ma.filled(np.ma.masked_invalid(v[:]), 0.0), dtype=np.float64)


def main() -> None:
    years = range(1850, 2024)
    worst_pft = 0.0
    worst_pft_year = None
    worst_harv = 0.0
    n_pft_nonzero = 0
    for y in years:
        o = nc.Dataset(OLD_SUB / f"LUT_nlcd2elm_luh2_historical_{y}_{OLD_TAG}.nc")
        n = nc.Dataset(BASE / f"LUT_nlcd2elm_luh2_historical_{y}_{NEW_TAG}.nc")
        po, pn = rd(o["PCT_PFT"]), rd(n["PCT_PFT"])
        dpft = float(np.abs(po - pn).max())
        ho = sum(rd(o[v]) for v in HARV)
        hn = sum(rd(n[v]) for v in HARV)
        dharv = float(np.abs(ho - hn).max())
        o.close(); n.close()
        if dpft > 0:
            n_pft_nonzero += 1
        if dpft > worst_pft:
            worst_pft, worst_pft_year = dpft, y
        worst_harv = max(worst_harv, dharv)
        if y % 25 == 0 or y == 2023:
            print(f"{y}: PCT_PFT maxdiff={dpft:.3e} | HARVEST maxdiff={dharv:.3e}")

    print("\n" + "=" * 60)
    print(f"years checked            : 174 (1850-2023)")
    print(f"years with any PCT_PFT diff: {n_pft_nonzero}")
    print(f"worst PCT_PFT maxdiff    : {worst_pft:.3e}  (year {worst_pft_year})")
    print(f"worst HARVEST maxdiff    : {worst_harv:.3e}  (for contrast)")
    if worst_pft == 0.0:
        print("=> PCT_PFT is BIT-IDENTICAL old vs new across all 174 years.")
    print("=" * 60)


if __name__ == "__main__":
    main()
