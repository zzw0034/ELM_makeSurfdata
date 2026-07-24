#!/bin/sh
# Pathfinder version (2026-07). Baseline original kept in git / see
# software/esmf_pathfinder_setup.md for the full Baseline->Pathfinder migration notes.
#
# Generates ${hgrid_name}.map_XX.run batch scripts from the map_*.run templates
# in this directory. Templates are Baseline-style; the -pathfinder edits
# (enabled by default) rewrite them for Pathfinder:
#   - partition/account/QOS: -p hpcl-cli185 -A hpcl-cli185 -q hpcl-cli185
#   - explicit -c 1 and --mem=460g (Pathfinder Slurm rejects jobs without them)
#   - /gpfs/wolf2/cades/cli185/... -> /projects/hpcl-cli185/...
#   - LD_LIBRARY_PATH for the Baseline-built ESMF (RPATH points to dead paths)
#   - UCX_TLS=rc,ud,dc,self,tcp (required: default UCX shm transports SIGBUS
#     on Pathfinder compute nodes; verified 2026-07-09, jobs 394169/394172/394173)
#   - -d <dst grid> rewritten to ${scrip_filepath}/${scrip_filename}

YYMMDD=`date +"%y%m%d"`
#scrip_file_name=northamericax4v1pg2_scrip.nc
#hgrid_name=northamericax4v1pg2

scrip_filename=SCRIPgrid_SEUS_1_24deg.nc
# Destination grid is produced by s1_Gridfiles; read the freshly generated copy
# from scr_out (NOT the stale world-shared one).
scrip_filepath=/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s1_Gridfiles/scr_out
hgrid_name=SEUS_1_24deg
# MPI-enabled ESMF_RegridWeightGen (Baseline-built copy, verified on Pathfinder).
ESMF_APP=/projects/hpcl-cli185/proj-shared/zw5/software/esmf-8.8.1-openmpi-gcc12/bin/binO/Linux.gfortran.64.openmpi.default/ESMF_RegridWeightGen
ESMF_LIBDIR=/projects/hpcl-cli185/proj-shared/zw5/software/esmf-8.8.1-openmpi-gcc12/lib/libO/Linux.gfortran.64.openmpi.default
# Pathfinder edits are enabled by default.
pathfinder=1
verbose=0

##################################################
# The command line help
##################################################
display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "   -hgrid_name     <name>             The hgrid name (e.g. SEUS_1_24deg)"
    echo "   -scrip_filename <netcdf_filename>  The destination SCRIP filename (e.g. SCRIPgrid_SEUS_1_24deg.nc)"
    echo "   -scrip_filepath <path>             The path to the SCRIP file"
    echo "   -esmf_app       <path>             Optional. Embed ESMF_APP=path and invoke \"\${ESMF_APP}\"."
    echo "   -pathfinder                        Apply ORNL Pathfinder edits (enabled by default)."
    echo "   -no-pathfinder                     Disable Pathfinder edits (keep raw Baseline templates)."
    echo "   -v, --verbose                      Set verbosity option true"
    echo
    echo "Example: "
    echo "   ./create_SEUS_mapping_scripts.sh   \\"
    echo "   -hgrid_name SEUS_1_24deg           \\"
    echo "   -scrip_filename SCRIPgrid_SEUS_1_24deg.nc \\"
    echo "   -scrip_filepath /projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/mappingdata/grids"
    echo
    exit 1
}


##################################################
# Get command line arguments
##################################################
while [ $# -gt 0 ]
do
  case "$1" in
    -hgrid_name )    hgrid_name="$2"; shift ;;
    -scrip_filename) scrip_filename="$2"; shift ;;
    -scrip_filepath)   scrip_filepath="$2"; shift ;;
    -esmf_app)        ESMF_APP="$2"; shift ;;
    -pathfinder)      pathfinder=1;;
    -no-pathfinder)   pathfinder=0;;
    -baseline)        echo "Note: -baseline is deprecated; Pathfinder edits stay ON. Use -no-pathfinder to disable."; pathfinder=1;;
    -no-baseline)     pathfinder=0;;
    -v | --verbose)  verbose=1;;
    -*)
      echo "Unknown option: $1"
      display_help
      exit 0
      ;;
    -h | --help)
      display_help
      exit 0
      ;;
    *)  break;;	# terminate while loop
  esac
  shift
done

if [ $verbose -eq 1 ]
then
  echo "Verbosity: On"
  echo " "
fi

if [ -z $hgrid_name ]
then
  echo "hgrid_name is not specified"
  display_help
  exit 0
fi

if [ -z $scrip_filename ]
then
  echo "scrip_filename is not specified"
  display_help
  exit 0
fi

if [ ! -f "${scrip_filepath}/${scrip_filename}" ]
then
  echo "WARNING: destination SCRIP grid not found: ${scrip_filepath}/${scrip_filename}"
fi

rm -rf $hgrid_name
mkdir -p $hgrid_name

echo "Creating batch scripts:"
# Templates: map_01.run, map_02.run, ... in this directory.
# Outputs:    ${hgrid_name}.map_01.run (e.g. SEUS_1_24deg.map_01.run)
for filename in map_*.run; do
  if [ ! -f "$filename" ]; then
    echo "No template map_*.run files in $(pwd); nothing to do."
    break
  fi
  out="$hgrid_name/${hgrid_name}.${filename}"
  echo "  $out"
  cp "$filename" "$out"
  sed -i "s/YYMMDD/${YYMMDD}/g"    "$out"
  sed -i "s/HGRID_NAME/${hgrid_name}/g" "$out"
  sed -i 's/ESMF_RegridWeightGen \\/"${ESMF_APP}" \\/g' "$out"
  # Insert ESMF_APP assignment after netcdf-fortran/4.6.1-mpi-h5f in the output script
  sed -i "/module load netcdf-fortran\/4.6.1-mpi-h5f/a ESMF_APP=${ESMF_APP}" "$out"

  if [ $pathfinder -eq 1 ]
  then
    # --- Slurm header: Baseline -> Pathfinder ---
    sed -i 's/^#SBATCH -p .*/#SBATCH -p hpcl-cli185/'            "$out"
    sed -i 's/^#SBATCH -A CLI185/#SBATCH -A hpcl-cli185/'        "$out"
    sed -i '/^#SBATCH -A hpcl-cli185/a #SBATCH -q hpcl-cli185'   "$out"
    # Pathfinder Slurm rejects jobs without explicit -c and --mem
    sed -i '/^#SBATCH --ntasks=/a #SBATCH -c 1\n#SBATCH --mem=460g' "$out"

    # --- Filesystem paths: Baseline -> Pathfinder ---
    sed -i 's|/gpfs/wolf2/cades/cli185|/projects/hpcl-cli185|g'  "$out"

    # --- Destination grid: honor -scrip_filepath / -scrip_filename ---
    sed -i "s|^-d /.*|-d ${scrip_filepath}/${scrip_filename} \\\\|" "$out"

    # --- Runtime env: ESMF libs (stale RPATH) + UCX shm SIGBUS workaround ---
    sed -i "/^ESMF_APP=/a export LD_LIBRARY_PATH=${ESMF_LIBDIR}:\${LD_LIBRARY_PATH:-}\nexport UCX_TLS=rc,ud,dc,self,tcp" "$out"
  fi

done
