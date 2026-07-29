#!/bin/bash
#SBATCH -J s4_3_nlcd1km_seus
#SBATCH -A hpcl-cli185
#SBATCH -p hpcl-cli185
#SBATCH -q hpcl-cli185
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o scr_out/%x-%j.out
#SBATCH -e scr_out/%x-%j.err

# Bin Annual NLCD (~60-70 m lat/lon) to 0.01 deg (~1 km) over the SEUS domain,
# map to ELM PFTs, and draw the figures. One year per run.
#
#   sbatch --export=NONE run_s4_3_nlcd_1km_SEUS.slurm          # defaults to 2023
#   sbatch --export=NONE,YEAR=2015 run_s4_3_nlcd_1km_SEUS.slurm

set -euo pipefail

YEAR="${YEAR:-2023}"

WORKDIR=/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft
cd "$WORKDIR"
mkdir -p scr_out

# The script is not multithreaded; keep the numeric libs from oversubscribing.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

echo "Job              : $SLURM_JOB_NAME ($SLURM_JOB_ID)"
echo "Node             : $(hostname)"
echo "Python           : $(which python)"
echo "Year             : $YEAR"
echo "Start            : $(date)"
echo "-------------------------------------------------------------"

python -u s4_3_nlcd_to_1km_SEUS.py --year "$YEAR"

echo "-------------------------------------------------------------"
echo "End              : $(date)"
