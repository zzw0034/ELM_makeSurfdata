#!/bin/bash -e
#SBATCH -p serial
#SBATCH -q normal
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32gb
#SBATCH -J updateSoilSurfdat
#SBATCH -o updateSoilSurfdat-%J.out
#SBATCH -e updateSoilSurfdat-%J.err

# Layout on Pathfinder:
#   .../SEUS_1_24deg/surfdata_results/
#       surfdata_SEUS_1_24deg_simyr1850_c260415.nc      (input, one level up)
#       postprocsSoil_onPF/
#           run_update_soil.slurm                       (this file)
#           update_Soil_surfdat.py
#           update_Soil_surfdat_hpc.py
#           environment.yml
#           surfdata_Updatedsoil_*.nc                   (output -> here)
#           intermediate/merged_surfdat/check_*.png     (plots -> here)
WORKDIR="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata"
INPUT_NC="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260712.nc"
OUTPUT_NC="Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc"

cd "${WORKDIR}"

# Activate the conda env created from environment.yml.
# Pathfinder doesn't set OLCF_MINIFORGE3_ROOT, so source the base conda
# directly (path matches `conda info --base` from the login node).
source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

# The HPC wrapper auto-injects Pathfinder defaults for --global-cf and
# --soil-order-tif (under /projects/hpcl-cli185/proj-shared/ywo/...).
python "${WORKDIR}/update_Soil_surfdat_hpc.py" \
  --surfdata "${INPUT_NC}" \
  --output   "${OUTPUT_NC}" \
  --smooth-soil-color \
  --soil-color-window 9 \
  --plot
