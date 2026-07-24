#!/usr/bin/env python3
"""
Generate mksrf_fdynuse list file for mksurfdata.
Format: one full path per line, padded to 195 chars, then space, then 4-digit year (year at column 197).

usage:

cd /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft

python3 make_fdynuse_list.py \
  --dir /gpfs/wolf2/cades/cli185/proj-shared/zw5/NLCD2ELM/LUToutput/output_downscaled_luh2_harvest \
  --pattern 'LUT_nlcd2elm_luh2_historical_{year}_04082026.nc' \
  --start 1850 \
  --end 2023 \
  -o LUT_nlcd2elm_luh2_historical_list.txt

"""

import os
import argparse

# Default: downscaled LUH2 harvest LUT files produced by s4_2_donwscale_LUH2harvest.py
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.join(SCRIPT_DIR, "output_downscaled_luh2_harvest")
DEFAULT_PATTERN = "LUT_nlcd2elm_luh2_historical_{year}_07232026.nc"
PATH_WIDTH = (
    195  # mksurfdata expects year at character 197 (1-based), so path + padding = 195
)


def main():
    p = argparse.ArgumentParser(
        description="Generate mksrf_fdynuse list file (path + padding, year at col 197)."
    )
    p.add_argument("--dir", default=DEFAULT_DIR, help="Directory containing LUT files")
    p.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help="Filename pattern with {year} placeholder",
    )
    p.add_argument("--start", type=int, default=1850, help="Start year")
    p.add_argument("--end", type=int, default=2023, help="End year (inclusive)")
    p.add_argument(
        "-o",
        "--out",
        default=None,
        help="Output list file (default: LUT_nlcd2elm_luh2_historical_list.txt in --dir)",
    )
    args = p.parse_args()

    out = args.out
    if out is None:
        out = os.path.join(args.dir, "LUT_nlcd2elm_luh2_historical_list.txt")

    lines = []
    for year in range(args.start, args.end + 1):
        fname = args.pattern.format(year=year)
        path = os.path.join(args.dir, fname)
        # Pad path to PATH_WIDTH so year starts at column 197 (1-based)
        line = path.ljust(PATH_WIDTH) + " " + str(year)
        lines.append(line)

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {out} ({len(lines)} years, {args.start}-{args.end})")


if __name__ == "__main__":
    main()
