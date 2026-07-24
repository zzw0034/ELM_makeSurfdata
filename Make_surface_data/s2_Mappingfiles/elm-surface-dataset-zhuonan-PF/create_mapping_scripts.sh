#!/bin/sh
# Generic generator: creates ${hgrid_name}.map_XX.run batch scripts from the
# map_*.run templates in this directory. Templates are ORNL-Baseline-style.
#
# Pathfinder (2026-07): pass -pathfinder to rewrite the generated scripts for
# ORNL Pathfinder (partition/account/QOS, explicit -c/--mem, /gpfs->/projects
# paths, ESMF LD_LIBRARY_PATH, UCX_TLS SIGBUS fix, -d dst grid substitution).
# See software/esmf_pathfinder_setup.md for the background of every edit.
# For SEUS grids prefer the pre-configured wrappers:
#   create_SEUS_mapping_scripts.sh / create_smallSEUS_mapping_scripts.sh

YYMMDD=`date +"%y%m%d"`
#scrip_file_name=northamericax4v1pg2_scrip.nc
#hgrid_name=northamericax4v1pg2

scrip_filename=
scrip_filepath=
hgrid_name=
verbose=0
# Optional: embed MPI-enabled ESMF_RegridWeightGen path (see -esmf_app).
# -pathfinder sets the verified Pathfinder default if unset.
ESMF_APP=
ESMF_LIBDIR=/projects/hpcl-cli185/proj-shared/zw5/software/esmf-8.8.1-openmpi-gcc12/lib/libO/Linux.gfortran.64.openmpi.default
# If set (via -pathfinder), apply ORNL Pathfinder edits (off by default here;
# the SEUS/smallSEUS wrapper scripts enable them by default).
pathfinder=0

##################################################
# The command line help
##################################################
display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "   -hgrid_name     <name>             The hgrid name (e.g. northamericax4v1pg2)"
    echo "   -scrip_filename <netcdf_filename>  The destination SCRIP filename (e.g. northamericax4v1pg2_scrip.nc)"
    echo "   -scrip_filepath <path>             The path to SCRIP file (e.g. /global/cfs/cdirs/e3sm/inputdata on NERSC)"
    echo "   -esmf_app       <path>             Optional. Embed ESMF_APP=path and invoke \"\${ESMF_APP}\" instead of bare ESMF_RegridWeightGen."
    echo "   -pathfinder                        ORNL Pathfinder: SBATCH header (-p/-A/-q/-c/--mem), /gpfs->/projects paths,"
    echo "                                      LD_LIBRARY_PATH + UCX_TLS env, -d dst grid substitution; sets default -esmf_app if omitted."
    echo "   -baseline                          Deprecated alias for -pathfinder (Baseline is retired)."
    echo "   -v, --verbose                      Set verbosity option true"
    echo
    echo "Example: "
    echo "   ./create_mapping_scripts.sh                  \\"
    echo "   -hgrid_name northamericax4v1pg2              \\"
    echo "   -scrip_filename northamericax4v1pg2_scrip.nc \\"
    echo "   -scrip_filepath ~/data                       \\"
    echo "   -pathfinder"
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
    -baseline)        echo "Note: -baseline is deprecated (Baseline retired); applying Pathfinder edits."; pathfinder=1;;
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

if [ "$pathfinder" -eq 1 ] && [ -z "$ESMF_APP" ]; then
  ESMF_APP=/projects/hpcl-cli185/proj-shared/zw5/software/esmf-8.8.1-openmpi-gcc12/bin/binO/Linux.gfortran.64.openmpi.default/ESMF_RegridWeightGen
fi

if [ "$pathfinder" -eq 1 ] && [ -n "$scrip_filepath" ] && [ ! -f "${scrip_filepath}/${scrip_filename}" ]; then
  echo "WARNING: destination SCRIP grid not found: ${scrip_filepath}/${scrip_filename}"
fi

rm -rf $hgrid_name
mkdir -p $hgrid_name

echo "Creating batch scripts:"
# Templates: map_01.run, ...  Outputs: ${hgrid_name}/${hgrid_name}.map_01.run
for filename in map_*.run; do
  if [ ! -f "$filename" ]; then
    echo "No template map_*.run files in $(pwd); nothing to do."
    break
  fi
  out="${hgrid_name}/${hgrid_name}.${filename}"
  echo "  $out"
  cp "$filename" "$out"
  sed -i "s/YYMMDD/${YYMMDD}/g"     "$out"
  sed -i "s/HGRID_NAME/${hgrid_name}/g" "$out"

  if [ -n "$ESMF_APP" ]; then
    sed -i 's/^ESMF_RegridWeightGen \\/"${ESMF_APP}" \\/' "$out"
    # Insert ESMF_APP assignment after the last module load line
    sed -i "/module load netcdf-fortran\/4.6.1-mpi-h5f/a ESMF_APP=${ESMF_APP}" "$out"
  fi

  if [ "$pathfinder" -eq 1 ]
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
    if [ -n "$scrip_filepath" ]; then
      sed -i "s|^-d /.*|-d ${scrip_filepath}/${scrip_filename} \\\\|" "$out"
    fi

    # --- Runtime env: ESMF libs (stale RPATH) + UCX shm SIGBUS workaround ---
    sed -i "/^ESMF_APP=/a export LD_LIBRARY_PATH=${ESMF_LIBDIR}:\${LD_LIBRARY_PATH:-}\nexport UCX_TLS=rc,ud,dc,self,tcp" "$out"
  fi

done
