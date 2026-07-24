"""Verify the regenerated landuse.timeseries (c260723) vs the old one (c260712).

Expectation after the harvest-only fix:
  - PCT_NAT_PFT : identical (harvest fix does not touch PFT composition)
  - GRAZING     : identical (hard-coded ones upstream)
  - HARVEST_*   : changed (~36x area fix, + SSP2-4.5 for recent years)
Reports max abs diff of PFT/GRAZING and the domain-total harvest ratio new/old.
"""

from __future__ import annotations
import numpy as np
import netCDF4 as nc

D = "/projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map"
OLD = f"{D}/landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260712.nc"
NEW = f"{D}/landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260723.nc"
HARV = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3"]


def rd(v):
    return np.asarray(np.ma.filled(np.ma.masked_invalid(v[:]), 0.0), dtype=np.float64)


def main() -> None:
    o = nc.Dataset(OLD)
    n = nc.Dataset(NEW)
    print("old vars:", sorted(o.variables))
    print("time old/new:", o.dimensions["time"].size if "time" in o.dimensions else "?",
          n.dimensions["time"].size if "time" in n.dimensions else "?")

    # PCT_NAT_PFT identity
    dp = float(np.abs(rd(o["PCT_NAT_PFT"]) - rd(n["PCT_NAT_PFT"])).max())
    print(f"\nPCT_NAT_PFT max|diff| = {dp:.3e}   (expect 0)")
    if "GRAZING" in o.variables and "GRAZING" in n.variables:
        dg = float(np.abs(rd(o["GRAZING"]) - rd(n["GRAZING"])).max())
        print(f"GRAZING     max|diff| = {dg:.3e}   (expect 0)")

    # harvest: domain-total per year (sum over the 5 vars), ratio new/old
    print("\nyear-index :  old_sum      new_sum     new/old")
    ho = sum(rd(o[v]) for v in HARV)  # (time, lat, lon)
    hn = sum(rd(n[v]) for v in HARV)
    nt = ho.shape[0]
    to = ho.reshape(nt, -1).sum(1)
    tn = hn.reshape(nt, -1).sum(1)
    for it in (0, 25, 50, 100, 150, 165, nt - 1):
        r = tn[it] / to[it] if to[it] > 0 else float("nan")
        print(f"  t={it:3d}   : {to[it]:.4e}  {tn[it]:.4e}  {r:.2f}")
    tot_o, tot_n = to.sum(), tn.sum()
    print(f"\nall-time harvest sum: old={tot_o:.4e}  new={tot_n:.4e}  new/old={tot_n/tot_o:.2f}")
    o.close(); n.close()


if __name__ == "__main__":
    main()
