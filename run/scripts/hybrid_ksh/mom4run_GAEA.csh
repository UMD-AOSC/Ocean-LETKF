#!/bin/csh -f
#PBS -r y                                                              #This job is restartable
#PBS -S /bin/sh                                                        #Do not change this - it keeps your job from issuing a false alarm
#PBS -E                                                                #Do not change this - it gives your job more and more useful Moab environment variables
# -- Request 120 cores
#PBS -l size=120
# 
# -- Specify a maximum wallclock of 4 hours
#PBS -l walltime=0:30:00
#
# -- Specify under which account a job should run
#PBS -A cpo_orr
#
# -- Set the name of the job, or moab will default to STDIN
#PBS -N cpo_ens
#
# -- Set the queue: debug, batch, novel, bigmem
#PBS -q batch
# 
# -- Set the partition (for Gaea)
#PBS -l partition=c1:c2
#
# -- Set this as working directory
#PBS -d .
#
# -- Send output and error to same file
#PBS -j oe
#

#Sample call: ./mom4run.csh <full path to working directory> mk3p51_001 mom4p1_solo 4
echo "File = $0"
echo "Work Directory  = $1"
echo "Experiment Name = $2 (make sure to include ensemble member id subscript)"
echo "Model Type      = $3"
echo "Number of Procs = $4"
#echo "llfile          = $5 (mom4p1.ll)"

#Source the system initialization scripts (but not the user's so as to avoid interactive settings that can poison us)
source /etc/csh.cshrc

set platform      = ftn      # A unique identifier for your platform
set name          = $2 #00            # Name of the experiment you want to run
                                   # One of box1, box_channel1, bowl1, dome1, gyre1, iom1, mk3p51, symmetric_box1, torus1
set npes          = $4             # number of processors
                                   # Note: If you change npes you may need to change
                                   # the layout in the corresponding namelist
set type          = $3 #mom4p1_solo_prod            # type of the experiment
set workdir       = $1             # where the model is run and model output is produced
                                   # This is recommended to be a link to the $WORKDIR of the platform.
set expdir        = $workdir/$name 
set inputDataDir  = $expdir/INPUT   # This is path to the directory that contains the input data for this experiment.
                                    # You should have downloaded and untared this directory from MOM4p1 FTP site.
set diagtable     = $inputDataDir/diag_table  # path to diagnositics table
set datatable     = $inputDataDir/data_table  # path to the data override table.
set fieldtable    = $inputDataDir/field_table # path to the field table
set namelist      = $inputDataDir/input.nml   # path to namelist file

set root          = /lustre/f1/unswept/Steve.Penny/mom4p1 #$cwd:h         # The directory in which you checked out src
set executable    = $root/exec_$platform/$type/fms_$type.x      # executable created after compilation
set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set time_stamp    = $root/bin/time_stamp.csh          # path to cshell to generate the date
#STEVE: copy all of these to the experiment directory prior to running:
#       (this is to avoid referencing the $HOME directory)
#set executable    = $expdir/fms_$type.x      # executable created after compilation
#set mppnccombine  = $expdir/mppnccombine.$platform  # path to executable mppnccombine
#set time_stamp    = $expdir/time_stamp.csh          # path to cshell to generate the date

#===========================================================================
# The user need not change any of the following
#===========================================================================

#
# Users must ensure the correct environment file exists for their platform.
#
#source $root/bin/environs.$platform  # environment variables and loadable modules

# Check if the user has extracted the input data
  if ( ! -d $inputDataDir ) then
    echo "ERROR: the experiment directory '$inputDataDir' does not exist or does not contain input and preprocessing data directories!"
    exit 1
  endif

set echo

# setup directory structure
  if ( ! -d $expdir )         mkdir -p $expdir
  if ( ! -d $expdir/RESTART ) mkdir -p $expdir/RESTART

#
#Check the existance of essential input files
#
   if ( ! -e $inputDataDir/grid_spec.nc ) then
     echo "ERROR: required input file does not exist $inputDataDir/grid_spec.nc "
     exit 1
   endif
   if ( ! -e $inputDataDir/ocean_temp_salt.res.nc ) then
     echo "ERROR: required input file does not exist $inputDataDir/ocean_temp_salt.res.nc "
     exit 1
   endif

#Change to expdir

  cd $expdir
  echo "pwd:: `pwd`"
  echo "Should be: $expdir"

# Create INPUT directory. Make a link instead of copy
# 
if ( ! -d $expdir/INPUT   ) mkdir -p $expdir/INPUT

  if ( ! -e $namelist ) then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  endif
  if ( ! -e $datatable ) then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  endif
  if ( ! -e $diagtable ) then
    echo "ERROR: required input file does not exist $diagtable "
    exit 1
  endif
  if ( ! -e $fieldtable ) then
    echo "ERROR: required input file does not exist $fieldtable "
    exit 1
  endif

  cp $namelist   $expdir/input.nml
  cp $datatable  $expdir/data_table
  cp $diagtable  $expdir/diag_table
  cp $fieldtable $expdir/field_table 


#   --- run the model ---

  cd $expdir

  if($npes > 1) then
#    mpiexec_mpt -np $npes $executable > $expdir/fmt.out
     echo "Calling: aprun -n $npes $executable > $expdir/fmt.out"
     cp $executable $expdir/
     aprun -n $npes $expdir/fms_$type.x > $expdir/fmt.out
  else
     $executable:t > fms.out
  endif

#----------------------------------------------------------------------------------------------
# generate date for file names ---
    set begindate = `$time_stamp -bf digital`
    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`
    set enddate = `$time_stamp -ef digital`
    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
    if ( -f time_stamp.out ) rm -f time_stamp.out
#----------------------------------------------------------------------------------------------
# get a tar restart file
  cd RESTART
  cp $expdir/input.nml $expdir/RESTART 
  cp $expdir/*_table $expdir/RESTART
# combine netcdf files
  if ( $npes > 1 ) then
    set file_previous = ""
    set multires = (`ls *.nc.????`)
    foreach file ( $multires )
	if ( $file:r != $file_previous:r ) then
	    set input_files = ( `ls $file:r.????` )
              if ( $#input_files > 0 ) then
                 $mppnccombine $file:r $input_files
                 if ( $status != 0 ) then
                   echo "ERROR: in execution of mppnccombine on restarts"
                   exit 1
                 endif
                 rm $input_files
              endif
           else
              continue
           endif
           set file_previous = $file
       end
  endif

  cd $expdir
  mkdir -p history
  mkdir -p ascii
#----------------------------------------------------------------------------------------------
# rename ascii files with the date
  foreach out (`ls *.out`)
     mv $out ascii/$begindate.$out
  end

#----------------------------------------------------------------------------------------------
# combine netcdf files
  if ( $npes > 1 ) then
    set file_previous = ""
    set multires = (`ls *.nc.????`)
    foreach file ( $multires )
	if ( $file:r != $file_previous:r ) then
	    set input_files = ( `ls $file:r.????` )
              if ( $#input_files > 0 ) then
                 $mppnccombine $file:r $input_files
                 if ( $status != 0 ) then
                   echo "ERROR: in execution of mppnccombine on restarts"
                   exit 1
                 endif
                 rm $input_files
              endif
           else
              continue
           endif
           set file_previous = $file
       end
  endif

#----------------------------------------------------------------------------------------------
# rename nc files with the date
  foreach ncfile (`/bin/ls *.nc`)
     mv $ncfile history/$begindate.$ncfile
  end

  unset echo

#Get the ctl files
#cp $root/data/$name/*.ctl $workdir/$name/RESTART
#cd $workdir

echo end_of_run
echo "NOTE: Natural end-of-script."

exit 0
  
