MODULE input_nml_oceanmodel

USE params_model
USE params_letkf
USE params_obs
!USE input_nml_letkf  !STEV: consider generalizing the more common items

PUBLIC :: read_input_namelist

PRIVATE
  ! Namelist inputs:  
 
#ifdef DYNAMIC
  ! Grid dimensions are set in params_model.f90, but only used in a namelist if compiled for dynamic arrays:
  NAMELIST /grid_dimensions_nml/ & 
                              nlon, &      ! number of longitude grid points (Can be specified via namelist or in netcdf gridfile)
                              nlat, &      ! number of latitude grid points
                              nlev         ! number of model levels
#endif

  NAMELIST /params_model_nml/ gridfile, &  ! Model grid definition file (e.g. NEMO mesh_mask.nc)
                              tsbase,& !
                              grid_lon2d_name,& !
                              grid_lat2d_name,& !
                              grid_lev_name,& !
                              grid_temp_name,& !
                              grid_salt_name,& !
                              grid_u_name,& !
                              grid_v_name,& !
                              grid_lsmask_name,& ! land/sea mask
                              grid_depth_name,& !
                              diag_lon_name,& !
                              diag_lat_name,& !
                              diag_lev_name,& !
                              diag_temp_name,& !
                              diag_salt_name,& !
                              diag_u_name,& !
                              diag_v_name,& !
                              diag_ssh_name,& !
                              rsrt_lon_name,& !
                              rsrt_lat_name,& !
                              rsrt_lev_name,& !
                              rsrt_temp_name,& !
                              rsrt_salt_name,& !
                              rsrt_u_name,& !
                              rsrt_v_name,& !
                              rsrt_ssh_name,& !
                              SSHclm_file,&   ! model ssh climatology for altimetry assimilation
                              DO_SUBGRID,&    ! logical flag to identify whether using a subgrid of the full global grid (e.g. prior to MPPNCCOMBINE)
                              do_physlimit,&  ! specify .true. (default) to enable bounds checks on analysis output
                              max_t, min_t, &
                              max_s, min_s, &
                              max_u, min_u, &
                              max_v, min_v, &
                              max_eta, min_eta

  NAMELIST /params_obs_nml/   obs1nrec, &            ! number of records in obs.dat type file
                              obs2nrec               ! number of records in obs2.dat type file

  NAMELIST /params_letkf_nml/ nbv, &                 ! Number of ensemble members
                              nslots, &              ! Number of time slots for 4D assimilation
                              nbslot, &              ! Index of base timeslot (time at which to form analysis)
                              sigma_obs, &           ! Sigma-radius (half-width) for horizontal localization at the equator (m)
                              sigma_obs0, &          ! Sigma-radius for horizontal localization at the poles (m)
                              sigma_obsv, &          ! Sigma-radius for vertical localization (m)
                              sigma_obst, &          ! Sigma-radius for temporal localization (not activated)
                              gross_error, &         ! number of standard deviations for quality control (all outside removed)
                              DO_WRITE_ENS_MEAN_SPRD, &  ! logical flag to write the ensemble mean and spread for the forecast and analysis fields
                              DO_DRIFTERS, &         ! logical flag to do lagrangian drifters assimilation
                              DO_ALTIMETRY, &        ! logical flag to do altimetry data assimilation
                              DO_SLA, &              ! logical flag to use SLA for altimetry
                              DO_ADT, &              ! logical flag to use Absolute Dynamic Topography (ADT) for altimetry
                              DO_NO_VERT_LOC, &      ! logical flag to skip all vertical localization and project weights (default)
                              DO_MLD, &              ! logical flag to use 2-layer vertical localization: (1) SFC and mixed layer, (2) below mixed layer
                              DO_MLD_MAXSPRD, &      ! logical flag to use maximum spread instead of mixed layer depth do to steps for DO_MLD 
                              DO_REMOVE_65N, &       ! option to remove points above 65N in letkf instead of in observation operators
                              DO_QC_MEANDEP, &       ! option to quality control observations based on mean departure
                              DO_QC_MAXDEP, &        ! option to quality control observation based on maximum departure across ensemble members
                              DO_UPDATE_H, &         ! option to update model layer thicknesses based on assimilation of observations
                              localization_method, & ! localization method to be used in letkf_local.f90 (default=1, latitude-dependent)
                              cov_infl_mul, &        ! multiplicative inflation factor (default=1.0, i.e. none)
                              sp_infl_add            ! additive inflation factor (default none)



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
#ifdef DYNAMIC
    READ(fid,nml=grid_dimensions_nml)
#endif
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

END MODULE input_nml_oceanmodel
