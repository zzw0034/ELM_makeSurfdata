#!/bin/bash
#SBATCH -J s5_area
#SBATCH -A hpcl-cli185
#SBATCH -p hpcl-cli185
#SBATCH -q hpcl-cli185
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64g
#SBATCH -t 02:00:00
#SBATCH -o /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest/%x_%j.out
#SBATCH -e /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest/%x_%j.err

unset SLURM_EXPORT_ENV
set -euo pipefail

DIAG_DIR=/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/diag_harvest
cd "$DIAG_DIR"
PY=/projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf/bin/python

echo "Start : $(date)  Node: $(hostname)"
$PY "$DIAG_DIR/s5_diag_area_src_vs_ours.py"
echo "End   : $(date)"
