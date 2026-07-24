#!/bin/bash -e
#SBATCH -p serial
#SBATCH -q normal
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32gb
#SBATCH -J generatePdata
#SBATCH -o generatePdata-%J.out
#SBATCH -e generatePdata-%J.err

# Generate the P-form NetCDF products consumed by 02_Update_Surfdata_Pdata_hpc.py.
#
# Inputs read by 01_Generate_P_data.py (absolute paths hard-coded in the script):
#   /projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/P_rock_den.nc
#       parent-material rock P density (kg P m-2, global 0.5 deg)
#   .../s6_updateSoildata/Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc
#       s6 output; supplies SOIL_ORDER + LATIXY/LONGXY (SEUS 1/24 deg grid)
#
# Outputs written under WORKDIR (relative paths in the script):
#   outputs/TPS_SEUS_1_24deg.nc            total soil P (TPS), final product
#   outputs/P_forms_SEUS_1_24deg.nc        TPS + 5 Hedley P-form pools, final product
#                                          (this is the file 02 reads)
#   intermediate/*.nc, *.tif, *.png        rockP_resampled / ppdi / eps / stage
#                                          maps + sanity plots for QGIS inspection
WORKDIR="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s7_updataPdata"

cd "${WORKDIR}"

# Reuse the env created for s6 (make_surfdata_pf); it already has
# xarray / scipy / netcdf4 / matplotlib-base / rasterio -- everything s7 needs.
source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

# Matplotlib has no display on a compute node -- force the Agg backend so
# the diagnostic PNGs still write cleanly and no GUI hooks are touched.
export MPLBACKEND=Agg

python "${WORKDIR}/01_Generate_P_data.py"
