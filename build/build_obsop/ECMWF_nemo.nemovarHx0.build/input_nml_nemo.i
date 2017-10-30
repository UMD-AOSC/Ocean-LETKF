# 1 "input_nml_nemo.f90"
MODULE input_nml_oceanmodel

USE params_model
USE params_letkf
USE params_obs

PUBLIC :: read_input_namelist

PRIVATE
  ! Namelist inputs:  
 








# 20
  NAMELIST /params_model_nml/ gridfile, &  ! MOM4 grid_spec.nc file 
                              gridfile1,&  !
                              gridfile2,&  !
                              gridfile3,&  !
                              tsbase,& !
                              uvbase,& !
                              grid_lon_name,& !
                              grid_lat_name,& !
                              grid_lev_name,& !
                              grid_temp_name,& !
                              grid_salt_name,& !
                              grid_u_name,& !
                              grid_v_name,& !
                              grid_h_name,& !
                              grid_lon2d_name,& !
                              grid_lat2d_name,& !
                              grid_wet_name,& !
                              grid_depth_name,& !
                              grid_height_name,& !
                              diag_lon_name,& !
                              diag_lat_name,& !
                              diag_lev_name,& !
                              diag_temp_name,& !
                              diag_salt_name,& !
                              diag_u_name,& !
                              diag_v_name,& !
                              diag_h_name,& !
                              diag_ssh_name,& !
                              diag_height_name,& !
                              rsrt_lon_name,& !
                              rsrt_lat_name,& !
                              rsrt_lev_name,& !
                              rsrt_temp_name,& !
                              rsrt_salt_name,& !
                              rsrt_u_name,& !
                              rsrt_v_name,& !
                              rsrt_h_name,& !
                              rsrt_ssh_name,& !
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
                              DO_UPDATE_H, &         ! option to update model layer thicknesses based on assimilation of observations
                              localization_method, & ! localization method to be used in letkf_local.f90
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



# 104
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
