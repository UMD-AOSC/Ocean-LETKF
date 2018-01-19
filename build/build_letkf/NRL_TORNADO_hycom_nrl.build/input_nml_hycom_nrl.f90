MODULE input_nml_oceanmodel

USE params_model
USE params_letkf
USE params_obs
USE vars_model, ONLY: z_lvl

IMPLICIT NONE

PUBLIC :: read_input_namelist, read_ncoda_namelists

PRIVATE

INTEGER :: m,n,kko
! Dummy variables for gridnl:
INTEGER :: nnest,nproj
REAL :: rlat
REAL, DIMENSION(1) :: delx, dely
! Dummy variables for oanl:
LOGICAL               :: adj_run, argo_bias, fcst_err, conflict, direct, dh_fcst, global, himem, ice_asm
LOGICAL               :: isop, isp_xtnd, mds_xtnd, mds_rej_sal, modas, pertob, pt_anl, ssh_asm, sst_asm
LOGICAL, DIMENSION(4) :: fcst
LOGICAL, DIMENSION(8) :: debug
INTEGER               :: n_cold_anl, n_it, n_pass, prf_hrs, ssh_hrs, upd_cyc
INTEGER, DIMENSION(2) :: corr_mdl
INTEGER, DIMENSION(4) :: fgat, redu
INTEGER, DIMENSION(7) :: deny
REAL                  :: bv_chk, dv_dz, lvl_nmo, rscl_cap, sal_std, ssh_del
REAL                  :: ssh_lim, ssh_std, tmp_std, topo_mn, topo_mx, vc_adj
REAL, DIMENSION(4)    :: cluster, offset, rscl, vscl
REAL, DIMENSION(10)   :: prf_slct
REAL, DIMENSION(172)  :: obs_err_scl
CHARACTER(4)          :: mask_opt, mld_src, prf_time, ssh_mean, ssh_time, vc_mdl, vc_src
CHARACTER(5)          :: topo_src
CHARACTER(12)         :: model
CHARACTER(4), DIMENSION(4) :: hc_mdl

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
                              SSHclm_file,&! model ssh climatology for altimetry assimilation
                              do_physlimit,& ! Apply physical bounds to the analysis output
                              max_t, &     ! Maximum temperature in analysis output
                              min_t, &     ! Minimum temperature in analysis output
                              max_s, &     ! Maximum salinity in analysis output
                              min_s, &     ! Minimum salinity in analysis output
                              max_uv, &    ! Maximum currents in analysis output
                              min_uv, &    ! Minimum currents in analysis output
                              max_ssh, &   ! Maximum sea surface height variation in analysis output
                              min_ssh, &   ! Minimum sea surface height variation in analysis output
                              max_h, &     ! Maximum layer height in analysis output
                              min_h, &     ! Minimum layer height in analysis output
                              DO_WRITE_TILE ! Logical to write tiled subgrid instead of updating a single global grid file
                              

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
                              DO_WRITE_ENS_MEAN_SPRD, &
                              localization_method, & ! localization method to be used in letkf_local.f90
                              cov_infl_mul, &        ! multiplicative inflation factor (default=1.0, i.e. none)
                              sp_infl_add            ! additive inflation factor (default none)

  NAMELIST /gridnl/           m, &                   ! glon, global grid longitude dimension from NCODA
                              n, &                   ! glat, ""
                              kko, &                 ! glev, ""
                              delx,&                 ! (NCODA) unused, placeholder only
                              dely,&                 ! ""
                              nnest,&                ! ""
                              nproj,&                ! ""
                              rlat                   ! ""
  
  NAMELIST /oanl/             z_lvl,&                ! z-levels from NCODA
                              adj_run,&              ! (NCODA) unused, placeholder only
                              argo_bias,&            ! ""
                              bv_chk,&               ! ""
                              fcst,&                 ! ""
                              fcst_err,&             ! ""
                              fgat,&                 ! ""
                              cluster,&              ! ""
                              conflict,&             ! ""
                              corr_mdl,&             ! ""
                              debug,&                ! ""
                              deny,&                 ! ""
                              direct,&               ! ""
                              dh_fcst,&              ! ""
                              dv_dz,&                ! ""
                              global,&               ! ""
                              hc_mdl,&               ! ""
                              himem,&                ! ""
                              ice_asm,&              ! ""
                              isop,&                 ! ""
                              isp_xtnd,&             ! ""
                              lvl_nmo,&              ! ""
                              mask_opt,&             ! ""
                              mds_xtnd,&             ! ""
                              mds_rej_sal,&          ! ""
                              modas,&                ! ""
                              model,&                ! ""
                              mld_src,&              ! ""
                              obs_err_scl,&          ! ""
                              offset,&               ! ""
                              n_cold_anl,&           ! ""
                              n_it,&                 ! ""
                              n_pass,&               ! ""
                              pertob,&               ! ""
                              prf_slct,&             ! ""
                              prf_time,&             ! ""
                              prf_hrs,&              ! ""
                              pt_anl,&               ! ""
                              redu,&                 ! ""
                              rscl,&                 ! ""
                              rscl_cap,&             ! ""
                              sal_std,&              ! ""
                              ssh_asm,&              ! ""
                              ssh_del,&              ! ""
                              ssh_hrs,&              ! ""
                              ssh_mean,&             ! ""
                              ssh_lim,&              ! ""
                              ssh_std,&              ! ""
                              ssh_time,&             ! ""
                              sst_asm,&              ! ""
                              tmp_std,&              ! ""
                              topo_mn,&              ! ""
                              topo_mx,&              ! ""
                              topo_src,&             ! ""
                              upd_cyc,&              ! ""
                              vc_adj,&               ! ""
                              vc_mdl,&               ! ""
                              vc_src,&               ! ""
                              vscl                   ! ""



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
