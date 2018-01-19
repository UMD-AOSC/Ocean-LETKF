#!/bin/bash
set -ex
date=/bin/date

#module load nco
#module load netcdf

# ncmax $var_nm $fl_nm : What is maximum of variable?
function ncmax {
  ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} tmp.max.nc | cut -f 3- -d ' ' ;
  rm -f tmp.max.nc
}
# ncmin $var_nm $fl_nm : What is minimum of variable?
function ncmin {
  ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} tmp.min.nc | cut -f 3- -d ' ' ;
  rm -f tmp.min.nc
}

source params.sh

# Link the observations to the subgrid working directories
SUBGRIDS=$GLOBAL_num_subgrids
ISLOT=1  # Use if there are multiple observation time bins

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME

SRCDIR=$GLOBAL_PERM/DATA/STATIC
DSTDIR0=$SCRATCH/$EXPNAME/WORK
mkdir -p $DSTDIR0

# Create a subdirectory in the working directory for each subgrid analysis
SG=$GLOBAL_subgrid_start
while [ $SG -lt $SUBGRIDS ]; do
  SG4=`printf "%04d" $SG`
  echo "Linking to subgrid number $SG out of $SUBGRIDS..."

  DSTDIR=$DSTDIR0/$SG4
  mkdir -p $DSTDIR
  cd $DSTDIR

  # Get sample restart file sugbrid information
  # Need to read: DOMAIN_position_first, DOMAIN_position_last 
  RESTART_SAMPLE=anal001.restart.nc
  cdlfile=anal001.header.cdl
  ncdump -h $RESTART_SAMPLE > $cdlfile

  gid=`grep 'DOMAIN_number ' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`
  ngrids=`grep 'DOMAIN_number_total' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`
  echo "domain $gid out of $ngrids"

  nlon=`grep 'DOMAIN_size_global' $cdlfile | grep -Eo '[0-9]*,' | grep -Eo '[0-9]*'`
  nlat=`grep 'DOMAIN_size_global' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`
  echo "global domain size : $nlon, $nlat"

  nlon_local=`grep 'DOMAIN_size_local' $cdlfile | grep -Eo '[0-9]*,' | grep -Eo '[0-9]*'`
  nlat_local=`grep 'DOMAIN_size_local' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`
  echo "local domain size : $nlon_local, $nlat_local"

  lon0=`grep 'DOMAIN_position_first' $cdlfile | grep -Eo '[0-9]*,' | grep -Eo '[0-9]*'`
  lat0=`grep 'DOMAIN_position_first' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`

  lon1=`grep 'DOMAIN_position_last' $cdlfile | grep -Eo '[0-9]*,' | grep -Eo '[0-9]*'`
  lat1=`grep 'DOMAIN_position_last' $cdlfile | grep -Eo '[0-9]*[[:blank:]];' | grep -Eo '[0-9]*'`

  echo "local lon indices : $lon0, $lon1"
  echo "local lat indices : $lat0, $lat1"

  # Subset the mesh_mask.nc file to these dimensions
  ncks -F -O -d x,$lon0,$lon1 -d y,$lat0,$lat1 $SRCDIR/mesh_mask.nc $DSTDIR/mesh_mask.nc

  var=mbathy  # for some reason 'tmask' doesn't work, maybe because it's in byte format? 
  mask_max=`ncmax $var mesh_mask.nc`
  mask_max="$(echo -e "${mask_max}" | tr -d '[:space:]')"
  mask_min=`ncmin $var mesh_mask.nc`
  mask_min="$(echo -e "${mask_min}" | tr -d '[:space:]')"

  echo "Min and max $var: $mask_min, $mask_max"

  if [ "$mask_max" -eq "0" ]; then
    echo "This grid contains no ocean points! Flagging..."
    touch no_ocean.skip
  fi

  SG=`expr $SG + 1`
done

