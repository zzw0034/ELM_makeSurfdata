#!/bin/bash
# Submit all generated mapping jobs in a case directory.
# Usage: ./submitall.sh [case_dir] [hgrid_name]
#   case_dir    defaults to the current directory
#   hgrid_name  defaults to case_dir's basename (case dirs are named after the hgrid,
#               and the generated scripts are named <hgrid>.map_*.run)
# Typical use: cd into the case dir, then  bash ../submitall.sh

CASEDIR=${1:-$PWD}
HGRID=${2:-$(basename "$CASEDIR")}

cd "$CASEDIR" || exit 1
mkdir -p done

for f in ${HGRID}.map_*.run; do
    [ -f "$f" ] || { echo "No ${HGRID}.map_*.run files in $CASEDIR"; exit 1; }
    echo "Submitting $f ..."
    sbatch "$f"
    sleep 1   # optional small delay between submissions
done
echo "All ${HGRID} mapping jobs submitted!"
