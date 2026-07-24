"""
HPC wrapper for ``update_Soil_surfdat.py`` (project hpcl-cli185 on Pathfinder).

Thin shim around :func:`update_Soil_surfdat.main`: forwards every CLI flag
unchanged but substitutes shared Pathfinder paths for ``--global-cf`` /
``--soil-order-tif`` (and a local ``surfdata*.nc`` for ``--surfdata``) when
the user omits them. All real work -- SOIL_ORDER mode regrid, bilinear
regrids, STDEV_ELEV from the law of total variance, SOIL_COLOR mode-filter
smoothing, EDT land-cell fill, copy-to-new-file output, diagnostic plots --
is done by the main script. See ``update_Soil_surfdat.py --help`` and
``update_Soil_surfdat_explain.md`` for the full pipeline.

Usage on Pathfinder
-------------------
# Run from a directory that has your input surfdata*.nc next to it.
# The wrapper finds the surfdata automatically; output is written next to
# it as ``surfdata_Updatedsoil_<basename>.nc`` (main script's default).
python update_Soil_surfdat_hpc.py --plot

# Override any of the shared paths with the usual main-script flags.
python update_Soil_surfdat_hpc.py \\
    --surfdata        /path/to/surfdata.nc \\
    --output          /path/to/my_output.nc \\
    --global-cf       /other/global_cf_float \\
    --soil-order-tif  /other/Soil_Order_NA_1km.tif \\
    --soil-color-window 5 \\
    --plot
"""

from __future__ import annotations

import os
import sys

from update_Soil_surfdat import main as _main

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Pathfinder (project hpcl-cli185) defaults for the shared 1 km inputs.
# Override with --global-cf / --soil-order-tif on the command line.
DEFAULT_GLOBAL_CF_DIR = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/Soil_Properties/global_cf_float"
)
DEFAULT_SOIL_ORDER_TIF = (
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/Soil_Properties/Soil_Order_NA_1km.tif"
)


def _default_surfdata_in_script_dir() -> str | None:
    """
    Return the first ``surfdata*.nc`` next to this wrapper, or ``None``.

    Prefers the small SEUS example file when present so a no-arg invocation
    has a deterministic target; otherwise falls back to the alphabetically
    first ``surfdata*.nc`` in the script directory.
    """
    preferred = os.path.join(
        _THIS_DIR, "surfdata_small_SEUS_1_24deg_simyr1850.nc"
    )
    if os.path.isfile(preferred):
        return preferred
    candidates = sorted(
        os.path.join(_THIS_DIR, name)
        for name in os.listdir(_THIS_DIR)
        if name.startswith("surfdata") and name.endswith(".nc")
    )
    return candidates[0] if candidates else None


def _has_flag(argv: list[str], flag: str) -> bool:
    """True if ``--flag`` or ``--flag=value`` already appears in argv."""
    eq = f"{flag}="
    return any(a == flag or a.startswith(eq) for a in argv)


def _inject_defaults(argv: list[str]) -> list[str]:
    """Add the Pathfinder defaults for any flag the user did not pass."""
    out = list(argv)
    if not _has_flag(out, "--global-cf"):
        out += ["--global-cf", DEFAULT_GLOBAL_CF_DIR]
    if not _has_flag(out, "--soil-order-tif"):
        out += ["--soil-order-tif", DEFAULT_SOIL_ORDER_TIF]
    if not _has_flag(out, "--surfdata"):
        fallback = _default_surfdata_in_script_dir()
        if fallback is not None:
            out += ["--surfdata", fallback]
    return out


def main(argv: list[str] | None = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    forwarded = _inject_defaults(list(argv))
    print(f"[hpc wrapper] forwarding to update_Soil_surfdat.main: {forwarded}")
    _main(forwarded)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
