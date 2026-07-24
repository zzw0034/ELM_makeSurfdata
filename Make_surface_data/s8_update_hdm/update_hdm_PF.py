"""
update_hdm_PF.py  --  Pathfinder (HPC) launcher for update_hdm.py
=================================================================

Bilinearly regrid coarse elmforc ``hdm`` onto the SEUS 1/24 deg surfdata
grid (``PFTDATA_MASK`` land mask), with optional QA PNG / GeoTIFF outputs.

This is a thin wrapper around ``update_hdm.py`` (same directory) that fills
in the Make_surface_data default paths so it can be run with no arguments.

Layout (project hpcl-cli185)
----------------------------
    HDM source (read-only shared elmforc input, NOT under Make_surface_data):
        /projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/firedata/
            elmforc.Li_20181205_mod_hist_SSP2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c240906.nc

    /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/
        s7_updataPdata/
            surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc  (grid + PFTDATA_MASK)
        s8_update_hdm/                                               (WORKDIR)
            update_hdm_PF.py                     (this launcher)
            update_hdm.py                          (regrid logic)
            run_update_hdm.slurm
            elmforc.Li_hdm_1_24x1_24_bilinear_SEUS_simyr1850-2100.nc  (output)
            check_plots/check_hdm_<year>.png        (--plot)
            check_plots/tif/check_hdm_<year>_{original,bilinear}.tif  (--tif)

Conda env: ``make_surfdata_pf`` (xarray, netCDF4, scipy, matplotlib, rasterio).

Usage on Pathfinder
-------------------
    cd ${WORKDIR}
    sbatch run_update_hdm.slurm

Or interactively (defaults below):
    python update_hdm_PF.py
    python update_hdm_PF.py --year 1850 --plot --tif
    python update_hdm_PF.py --help
"""

from __future__ import annotations

import os
import sys

# update_hdm.py now lives in the same directory as this launcher.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from update_hdm import main  # noqa: E402

# ---------------------------------------------------------------------------
# Default paths (Make_surface_data layout, hpcl-cli185)
# ---------------------------------------------------------------------------
WORKDIR = _SCRIPT_DIR

HDM_INPUT = (
    "/projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/firedata/"
    "elmforc.Li_20181205_mod_hist_SSP2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c240906.nc"
)
SURFDATA = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/"
    "s7_updataPdata/surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc"
)
OUTPUT_NC = os.path.join(
    WORKDIR,
    "elmforc.Li_hdm_1_24x1_24_bilinear_SEUS_simyr1850-2100.nc",
)
PLOT_DIR = os.path.join(WORKDIR, "check_plots")


def _build_default_argv() -> list[str]:
    os.makedirs(os.path.dirname(OUTPUT_NC) or ".", exist_ok=True)
    return [
        sys.argv[0],
        HDM_INPUT,
        "--surfdata",
        SURFDATA,
        "--output",
        OUTPUT_NC,
        "--year",
        "1850",
        "--plot",
        "--plot-dir",
        PLOT_DIR,
        "--tif",
        "--tif-dir",
        os.path.join(PLOT_DIR, "tif"),
    ]


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv = _build_default_argv()
    main()
