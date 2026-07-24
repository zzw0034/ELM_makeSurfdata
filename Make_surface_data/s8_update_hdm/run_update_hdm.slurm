#!/bin/bash -e
#SBATCH -p serial
#SBATCH -q normal
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32gb
#SBATCH -J updHdm
#SBATCH -o updHdm-%J.out
#SBATCH -e updHdm-%J.err

# Bilinearly regrid elmforc hdm (0.5 deg) onto the SEUS 1/24 deg surfdata grid.
# Land-mask gap-fill + optional QA PNG / GeoTIFF (first processed time slice).
#
# Layout (project hpcl-cli185):
#   HDM source (read-only shared elmforc input, NOT under Make_surface_data):
#     /projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/firedata/
#         elmforc.Li_20181205_mod_hist_SSP2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c240906.nc
#
#   .../ELM_makeSurfdata/Make_surface_data/
#       s7_updataPdata/
#           surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc  (grid + PFTDATA_MASK)
#       s8_update_hdm/                                               (WORKDIR)
#           update_hdm.py
#           update_hdm_PF.py
#           run_update_hdm.slurm                                     (this file)
#           elmforc.Li_hdm_1_24x1_24_bilinear_SEUS_simyr1850-2100.nc (output)
#           check_plots/check_hdm_1850.png
#           check_plots/tif/check_hdm_1850_{original,bilinear}.tif
#
# Submit:
#   cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s8_update_hdm
#   sbatch run_update_hdm.slurm

WORKDIR="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s8_update_hdm"
HDM_INPUT="/projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/firedata/elmforc.Li_20181205_mod_hist_SSP2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c240906.nc"
SURFDATA="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s7_updataPdata/surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c260712.nc"
OUTPUT_NC="${WORKDIR}/elmforc.Li_hdm_1_24x1_24_bilinear_SEUS_simyr1850-2100.nc"
PLOT_DIR="${WORKDIR}/check_plots"
TIF_DIR="${PLOT_DIR}/tif"

# Set YEAR=1850 for a single-year test + QA plots; leave empty for all 251 years.
YEAR=

mkdir -p "${WORKDIR}" "${PLOT_DIR}" "${TIF_DIR}"
cd "${WORKDIR}"

source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

ARGS=(
  "${HDM_INPUT}"
  --surfdata "${SURFDATA}"
  --output "${OUTPUT_NC}"
  --plot --plot-dir "${PLOT_DIR}"
  --tif --tif-dir "${TIF_DIR}"
)
if [[ -n "${YEAR}" ]]; then
  ARGS+=(--year "${YEAR}")
fi

python "${WORKDIR}/update_hdm.py" "${ARGS[@]}"
