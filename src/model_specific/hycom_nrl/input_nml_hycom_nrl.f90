MODULE input_nml_oceanmodel

USE params_model
USE params_letkf
USE params_obs
USE vars_model, ONLY: z_lvl

IMPLICIT NONE

PUBLIC :: read_input_namelist, read_ncoda_namelists

PRIVATE

INTEGER :: m,n,kko

  ! Namelist inputs:  
  NAMELIST /params_model_nml/ gridfile_lon,&  ! NCODA file containing model grid information
                              gridfile_lat,&  !
                              gridfile_lev,&  !
                              gridfile_kmt,&  !
                              glon, &      ! global number of longitude grid points (Can be specified via namelist)
                              glat, &      ! global number of latitude grid points
                              glev, &      ! global number of model levels
                              nlon, &      ! local number of longitude grid points (Can be specified via namelist or input file)
                              nlat, &      ! local number of latitude grid points
                              nlev, &      ! local number of model levels
                              istart, &    ! LAM start grid i index
                              iend, &      ! LAM end grid i index
                              jstart, &    ! LAM start grid j index
                              jend, &      ! LAM end grid j index
                              SSHclm_file  ! model ssh climatology for altimetry assimilation

  NAMELIST /params_obs_nml/   obs1nrec, &  ! number of records in obs.dat type file
                              obs2nrec     ! number of records in obs2.dat type file

  NAMELIST /params_letkf_nml/ nbv, &                 ! Number of ensemble members
                              nslots, &              ! Number of time slots for 4D assimilation
                              nbslot, &              ! Index of base timeslot (time at which to form analysis)
                              sigma_obs, &           ! Sigma-radius (half-width) for horizontal localization at the equator (m)
                              sigma_obs0, &          ! Sigma-radius for horizontal localization at the poles (m)
                              sigma_obsv, &          ! Sigma-radius for vertical localization (m)
                              sigma_obst, &          ! Sigma-radius for temporal localization (not activated)
                              gross_error, &         ! number of standard deviations for quality control (all outside removed)
                              DO_DRIFTERS, &         ! logical flag to do lagrangian drifters assimilation
                              DO_ALTIMETRY, &        ! logical flag to do altimetry data assimilation
                              DO_SLA, &              ! logical flag to use SLA for altimetry
                              DO_ADT, &              ! logical flag to use Absolute Dynamic Topography (ADT) for altimetry
                              DO_NO_VERT_LOC, &      ! logical flag to skip all vertical localization and project weights (default)
                              DO_MLD, &              ! logical flag to use 2-layer vertical localization: (1) SFC and mixed layer, (2) below mixed layer
                              DO_MLD_MAXSPRD, &      ! logical flag to use maximum spread instead of mixed layer depth do to steps for DO_MLD 
                              DO_REMOVE_65N, &       ! option to remove points above 65ÂN in letkf instead of in observation operators
                              DO_QC_MEANDEP, &       ! option to quality control observations based on mean departure
                              DO_QC_MAXDEP, &        ! option to quality control observation based on maximum departure across ensemble members
                              localization_method, & ! localization method to be used in letkf_local.f90
                              cov_infl_mul, &        ! multiplicative inflation factor (default=1.0, i.e. none)
                              sp_infl_add            ! additive inflation factor (default none)

  NAMELIST /gridnl/           m, &                   ! glon, global grid longitude dimension from NCODA
                              n, &                   ! glat, ""
                              kko                    ! glev, ""
  
  NAMELIST /oanl/             z_lvl                  ! z-levels from NCODA



CONTAINS

SUBROUTINE read_input_namelist
!===============================================================================
! Subroutine to read the parameter from an input namelist (input.nml)
!===============================================================================

  LOGICAL :: ex
  INTEGER :: fid=21

  INQUIRE(FILE="input.nml", EXIST=ex)
  if (ex) then
    OPEN(fid,file="input.nml", status='OLD') !, delim='APOSTROPHE')
    READ(fid,nml=params_model_nml)
    READ(fid,nml=params_obs_nml)
    READ(fid,nml=params_letkf_nml)
  endif

  if (DO_ADT .or. DO_SLA) DO_ALTIMETRY=.true.

  WRITE(6,*) "================================================================="
  WRITE(6,*) "Namelist inputs:"
  WRITE(6,*) "================================================================="
  WRITE(6,params_model_nml)
  WRITE(6,params_obs_nml)
  WRITE(6,params_letkf_nml)
  WRITE(6,*) "================================================================="

END SUBROUTINE read_input_namelist

SUBROUTINE read_ncoda_namelists
!===============================================================================
! Subroutine to read the NCODA grid parameters from input namelists
!===============================================================================

  LOGICAL :: ex
  INTEGER :: fid=21

  INQUIRE(FILE="gridnl", EXIST=ex)
  if (ex) then
    OPEN(fid,file="gridnl", status='OLD')
    READ(fid,nml=gridnl)
    CLOSE(fid)
  else
    WRITE(6,*) 'Please supply NCODA gridnl namelist file. EXITING...'
  endif

  INQUIRE(FILE="oanl", EXIST=ex)
  if (ex) then
    OPEN(fid,file="oanl", status='OLD')
    READ(fid,nml=oanl)
    CLOSE(fid)
  else
    WRITE(6,*) 'Please supply NCODA oanl namelist file. EXITING...'
  endif

  WRITE(6,*) "================================================================="
  WRITE(6,*) "NCODA Namelist inputs:"
  WRITE(6,*) "================================================================="
  WRITE(6,gridnl)
  WRITE(6,oanl)
  WRITE(6,*) "================================================================="

  !-----------------------------------------------------------------------------
  ! Get the grid dimensions and the reference density levels:
  !-----------------------------------------------------------------------------
  glon = m
  glat = n
  glev = kko

  WRITE(6,*) "m = ", m
  WRITE(6,*) "n = ", n
  WRITE(6,*) "kko = ", kko

END SUBROUTINE read_ncoda_namelists

END MODULE input_nml_oceanmodel
