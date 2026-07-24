#!/bin/bash
#SBATCH -J make_scrip_1_24
#SBATCH -A hpcl-cli185
#SBATCH -p hpcl-cli185
#SBATCH -q hpcl-cli185
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o scr_out/%x-%j.out
#SBATCH -e scr_out/%x-%j.err

set -euo pipefail

# --- Working directory ---
WORKDIR=/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s1_Gridfiles
cd "$WORKDIR"
mkdir -p scr_out

# --- Keep numeric libs single-threaded (avoids OpenBLAS oversubscription) ---
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- Activate conda environment ---
source /software/baseline/nsp/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate /projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf

echo "Job    : $SLURM_JOB_NAME ($SLURM_JOB_ID)"
echo "Node   : $(hostname)"
echo "Python : $(which python)"
echo "Start  : $(date)"
echo "-------------------------------------------------------------"

python -u make_scrip_1_24_deg.py

echo "-------------------------------------------------------------"
echo "End    : $(date)"
