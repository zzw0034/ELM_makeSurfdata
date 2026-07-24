#!/bin/bash -e
#SBATCH -p serial
#SBATCH -q normal
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32gb
#SBATCH -J updatePsurfdat
#SBATCH -o updatePsurfdat-%J.out
#SBATCH -e updatePsurfdat-%J.err

# Layout (project hpcl-cli185):
#   .../ELM_makeSurfdata/Make_surface_data/
#       s6_updateSoildata/
#           Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc   (s6 output, s7 input)
#       s7_updataPdata/                                              (WORKDIR)
#           run_update_P.slurm                                       (this file)
#           02_Update_Surfdata_Pdata_hpc.py
#           outputs/P_forms_SEUS_1_24deg.nc                          (from s7 step 01)
#           surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c<today>.nc (surfdata output -> here)
#           landuse.timeseries_..._UpdatedP_c<today>.nc              (landuse output -> here, only if UPDATE_LANDUSE_TIMESERIES=True)
#           intermediate/check_surfdata_P_update.png                 (sanity plot -> here)
WORKDIR="/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s7_updataPdata"

cd "${WORKDIR}"

# Reuse the env created for s6 (make_surfdata_pf); it already has
# xarray / scipy / netcdf4 / matplotlib-base -- everything s7 needs.
# See ../postprocsSoil_onPF/environment.yml for how it was built.
source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

# Matplotlib has no display on a compute node -- force the Agg backend so
# the diagnostic PNG still writes cleanly and no GUI hooks are touched.
export MPLBACKEND=Agg

python "${WORKDIR}/02_Update_Surfdata_Pdata_hpc.py"
