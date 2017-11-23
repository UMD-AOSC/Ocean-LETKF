PROGRAM obsop_ecmwf_tprof
!===============================================================================
! PROGRAM: obsop_ecmwf_tprof
! 
! USES:
!  use common
!  use params_model
!  use vars_model
!  use common_oceanmodel
!  use params_obs
!  use vars_obs
!  use common_obs_oceanmodel
!
!
! DESCRIPTION: 
!   This program acts as the observation operator for temperature profiles. 
!   It inputs observations of temperature profiles (yo)
!   and a single member's model equivalent H(xf) used to form the 
!   innovations (yo-H(xf)) associated with that member.
!
! USAGE:
!   A separate instance is run independently for each member and timeslot
!   All observations are read in and converted to the letkf obs2 format w/ H(x) data 
!   for each member added in a new column.
! 
! !REVISION HISTORY:
!   10/27/2017 Steve Penny set up for reading ECMWF NEMOVAR feedback (fdbk) files
!   03/17/2016 Steve Penny separated obsop.f90 into separate files
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
!   04/03/2013 Takemasa Miyoshi created for SPEEDY atmospheric model.
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny, Takemasa Miyoshi $
!===============================================================================

  USE common
  USE params_model
  USE vars_model
  USE common_oceanmodel
  USE params_obs,                ONLY: nobs, id_t_obs, id_s_obs, id_u_obs, id_v_obs, id_eta_obs
  USE params_obs,                ONLY: DO_REMOVE_65N
  USE vars_obs
  USE common_obs_oceanmodel
  USE read_ecmwf_fdbk,           ONLY: argo_data, read_fdbk_nc
  ! For dynamic instantiation and use of namelists:
  USE input_nml_oceanmodel,      ONLY: read_input_namelist

  IMPLICIT NONE

  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: obsinfile2='obsin2.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)

  ! For storing observation data:
  TYPE(argo_data), ALLOCATABLE, DIMENSION(:) :: obs_data
  
  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  REAL(r_size), ALLOCATABLE :: obhr(:)
  INTEGER     , ALLOCATABLE :: oqc(:)
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size) :: dk,tg,qg
  REAL(r_size) :: ri,rj,rk
  INTEGER :: i,j,k,n

  ! For climatological salinity:
  REAL(r_size), ALLOCATABLE :: selm(:)
  REAL(r_size), ALLOCATABLE :: slon(:)
  REAL(r_size), ALLOCATABLE :: slat(:)
  REAL(r_size), ALLOCATABLE :: slev(:)
  REAL(r_size), ALLOCATABLE :: sdat(:)

  ! For observation error modification:
  REAL(r_size) :: obserr_scaling=1.0d0 !STEVE: use this to scale the input observations
  REAL(r_size) :: obs_randselect=1.0d0 !STEVE: use this to select random input observations
  REAL(r_size), DIMENSION(1) :: rand

  LOGICAL :: remap_obs_coords = .false.  ! Remap to mom4p1's -280 to 80 longitude range

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  LOGICAL :: debug_obsfilter = .false.
  LOGICAL :: debug_hdxf_0 = .false.   ! This error occured because there was not a model representation 
                                      ! of the observed value (i.e. SST obs with no SST model field).
                                      ! Solution was to populate a SST model field (v2d) with 
                                      ! surface temp data from the model (v3d(:,:,1)).
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
  LOGICAL :: print1st = .true.
  LOGICAL :: dodebug1 = .false.

  INTEGER :: cnt_obs_u=0, cnt_obs_v=0, cnt_obs_t=0, cnt_obs_s=0
  INTEGER :: cnt_obs_ssh=0, cnt_obs_sst=0, cnt_obs_sss=0, cnt_obs_eta=0
  INTEGER :: cnt_obs_x=0, cnt_obs_y=0, cnt_obs_z=0
  INTEGER :: cnt_thin=0
  INTEGER :: cnt_largedep=0 ! (NEMO)
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0, cnt_oerlt0=0, cnt_altlatrng=0
  INTEGER :: cnt_nl1=0, cnt_nl2=0, cnt_nl3=0

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  LOGICAL :: DO_REMOVE_BLACKSEA=.false.
  INTEGER :: cnt_blacksea=0
  REAL(r_size) :: dep_thresh = 100
  CHARACTER(slen) :: obs_name='POTM_OBS' !(NEMO) default. Change via command line argument
  CHARACTER(slen) :: hxb_name='POTM_Hx'  !(NEMO) default. Change via command line argument

  ! BEGIN
  
  !----------------------------------------------------------------------------
  ! Read in namelist parameters
  !----------------------------------------------------------------------------
  CALL read_input_namelist !STEVE: this has been moved to input_nml_{oceanmodel}.f90 since it needed to be slightly different for each model

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
  CALL set_common_oceanmodel !(gets the model grid information)
  CALL process_command_line  !(example options: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  !-----------------------------------------------------------------------------
  ! Read observations from ECMWF NEMOVAR feedback (fdbk) file
  !-----------------------------------------------------------------------------
  print *, "obsop_ecmwf_tprof.f90:: calling read_fdbk_nc for observations..."
  CALL read_fdbk_nc(obsinfile,id_t_obs,obs_data,nobs,trim(obs_name),trim(hxb_name),1) !STEVE: use ,2) if obs errors can be read in
  print *, "obsop_ecwmf_tprof.f90:: finished read_fdbk_nc."

  !-----------------------------------------------------------------------------
  ! Allocate memory for observation data
  !-----------------------------------------------------------------------------
  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  ALLOCATE( obhr(nobs) )

  !-----------------------------------------------------------------------------
  ! Assign observation and model equivalent data to data arrays in letkf format
  !-----------------------------------------------------------------------------
  print *, "obsop_tprof.f90:: starting nobs = ", nobs
  do i=1,nobs
    elem(i) = obs_data(i)%typ
    rlon(i) = obs_data(i)%x_grd(1)
    rlat(i) = obs_data(i)%x_grd(2)
    rlev(i) = obs_data(i)%x_grd(3)
    odat(i) = obs_data(i)%value
    oerr(i) = obs_data(i)%oerr
    ohx(i)  = obs_data(i)%hxb
    oqc(i)  = 0  ! default = don't use the ob - update below during qc checks
    obhr(i) = obs_data(i)%hour
  enddo
  if (print1st) then 
    print *, "elem(1),rlon(1),rlat(1),obhr(1) = ", elem(1),rlon(1),rlat(1),obhr(1)
    print *, "odat(1:40) = ", odat(1:40)
    print *, "oerr(1:40) = ", oerr(1:40)
  endif

  DEALLOCATE(obs_data)

  !-----------------------------------------------------------------------------
  ! Update the coordinate to match the model grid
  ! (Stolen from common_obs_mom4.f90::read_obs
  !-----------------------------------------------------------------------------
  if (remap_obs_coords) then
    CALL center_obs_coords(rlon,oerr,nobs)
  endif

  !-----------------------------------------------------------------------------
  ! Cycle through all observations and do basic quality control checks
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Cycling through observations..."
  oqc=0
  DO n=1,nobs
    !---------------------------------------------------------------------------
    ! Count bad obs errors associated with observations, and skip the ob
    !---------------------------------------------------------------------------
    if (oerr(n) <= 0) then
      !STEVE: this occurred for a few synthetic obs, perhaps due to Dave's code generating obs errors
      cnt_oerlt0 = cnt_oerlt0 + 1
      CYCLE
    endif

    !---------------------------------------------------------------------------
    ! Count bad departures, usually due to hxb being 'missing' where observations exist (NEMO) (NEMOVAR)
    !---------------------------------------------------------------------------
    if (abs(odat(n) - ohx(n)) > dep_thresh) then
      !STEVE: this occurred for a few synthetic obs, perhaps due to Dave's code generating obs errors
      cnt_largedep = cnt_largedep + 1
      CYCLE
    endif

    !---------------------------------------------------------------------------
    ! Filter in the tripolar region until localization is examined in the arctic !(ISSUE)
    !---------------------------------------------------------------------------
    if (DO_REMOVE_65N .and. rlat(n) > 65) then
      if (verbose) WRITE(6,'(A)') "Latitude above 65.0N, in tripolar region. Removing observation..."
      cnt_triout = cnt_triout + 1
      CYCLE
    endif
    
    !---------------------------------------------------------------------------
    ! Filter out observations in the black sea, which may cause problems
    !---------------------------------------------------------------------------
    if (DO_REMOVE_BLACKSEA) then
      if (40 < rlat(n) .and. rlat(n) < 50 .and. &
          25 < rlon(n) .and. rlon(n) < 45) then
        cnt_blacksea = cnt_blacksea + 1
        CYCLE
      endif
    endif

    !---------------------------------------------------------------------------
    ! (OPTIONAL) Scale the observation errors
    !---------------------------------------------------------------------------
    scale_obs : if (abs(obserr_scaling - 1.0d0) .gt. TINY(1.0d0) ) then
      oerr(n) = oerr(n) * obserr_scaling
    endif scale_obs

    !---------------------------------------------------------------------------
    ! (OPTIONAL) Thin the observations
    !---------------------------------------------------------------------------
    ! Select random subset of observations to speed up processing
    thin_obs : if (abs(obs_randselect - 1) > TINY(1.0d0)) then
      CALL com_rand(1,rand)
      if (rand(1) > obs_randselect) then
        cnt_thin = cnt_thin+1
        CYCLE
      endif
    endif thin_obs

    ! Debug output
    if (dodebug1 .AND. (abs(odat(n)-ohx(n)) > oerr(n)*5.0)) then
        print *, "======================================"
        print *, "n = ", n
        print *, "rlon(n),rlat(n),rlev(n) = ", rlon(n),rlat(n),rlev(n)
        print *, "ohx(n)  = ", ohx(n)
        print *, "odat(n) = ", odat(n)
        print *, "odep(n) = ", odat(n) - ohx(n)
        print *, "======================================"
    endif

    oqc(n) = 1
  enddo !1:nobs

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_tprof.f90:: observations removed for:"
  WRITE(6,*) "cnt_oerlt0 = ", cnt_oerlt0
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_zout = ", cnt_zout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_nl1 = ", cnt_nl1
  WRITE(6,*) "cnt_nl2 = ", cnt_nl2
  WRITE(6,*) "cnt_altlatrng = ", cnt_altlatrng
  WRITE(6,*) "cnt_thin = ", cnt_thin
  WRITE(6,*) "cnt_largedep = ", cnt_largedep !(NEMO)
  if (DO_REMOVE_BLACKSEA) WRITE(6,*) "cnt_blacksea = ", cnt_blacksea

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated model equivalent to file
  !-----------------------------------------------------------------------------
  CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)

  if (ALLOCATED(elem)) DEALLOCATE(elem)
  if (ALLOCATED(rlon)) DEALLOCATE(rlon)
  if (ALLOCATED(rlat)) DEALLOCATE(rlat)
  if (ALLOCATED(rlev)) DEALLOCATE(rlev)
  if (ALLOCATED(odat)) DEALLOCATE(odat)
  if (ALLOCATED(oerr)) DEALLOCATE(oerr)
  if (ALLOCATED(ohx )) DEALLOCATE(ohx)
  if (ALLOCATED(oqc )) DEALLOCATE(oqc)
  if (ALLOCATED(obhr)) DEALLOCATE(obhr)
  if (ALLOCATED(v3d )) DEALLOCATE(v3d)
  if (ALLOCATED(v2d )) DEALLOCATE(v2d)

CONTAINS

SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In obsop_tprof.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
!   case('-nlon')
!     CALL GET_COMMAND_ARGUMENT(i+1,arg2)
!     PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
!     read (arg2,*) nlon
!   case('-nlat')
!     CALL GET_COMMAND_ARGUMENT(i+1,arg2)
!     PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
!     read (arg2,*) nlat
!   case('-nlev')
!     CALL GET_COMMAND_ARGUMENT(i+1,arg2)
!     PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
!     read (arg2,*) nlev
    case('-obs_name')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obs_name = arg2
    case('-hxb_name')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      hxb_name = arg2
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-obsin2')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile2 = arg2
    case('-gues')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      guesfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case('-rmBlackSea')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_BLACKSEA
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug1
    case('-remap')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) remap_obs_coords
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM obsop_ecmwf_tprof
