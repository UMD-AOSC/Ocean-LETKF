# 1 "params_model.f90"
MODULE params_model

USE common, ONLY: r_size, slen

IMPLICIT NONE

PUBLIC








  ! NEMO 1/4-degree, ECMWF ORAS5
# 17
  INTEGER,PARAMETER :: nlon=1442
  INTEGER,PARAMETER :: nlat=1021
  INTEGER,PARAMETER :: nlev=75


# 22
  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
  INTEGER,PARAMETER :: nv2d=1 ! ssh                  !(OCEAN)
  INTEGER,PARAMETER :: nv4d=0 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS

  ! Initialize via subroutine below:
  INTEGER,SAVE :: nij0           ! Number of gridpoints handled by myrank processor
  INTEGER,SAVE :: nlevall        ! Total number of variables and levels (3d + 2d)
  INTEGER,SAVE :: ngpv           ! Total number of gridpoints, including nij0*nlevall

  ! Identifiers for each analyzed variable in 3D
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                !(OCEAN)
  INTEGER,PARAMETER :: iv3d_h=5                !(OCEAN) (MOM6)

  ! Identifiers for each analyzed variable in 2D
  INTEGER,PARAMETER :: iv2d_ssh=1              !(OCEAN) ! time averaged thickness of top model grid cell (m)
  INTEGER,PARAMETER :: iv2d_sst=2              !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model
  INTEGER,PARAMETER :: iv2d_sss=3              !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models
  INTEGER,PARAMETER :: iv2d_eta=4              !(OCEAN) ! model-derived eta sea surface perturbation
  INTEGER,PARAMETER :: iv2d_mld=5              !(OCEAN) ! mixed layer depth

  ! Identifiers for each analyzed variable in 4D
  INTEGER,PARAMETER :: iv4d_x=1                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_y=2                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_z=3                !(OCEAN) (DRIFTERS)

  !
  ! Elements
  !
  CHARACTER(4) :: element(nv3d+nv2d+nv4d)
! element(iv3d_u) = 'U   '
! element(iv3d_v) = 'V   '
! element(iv3d_t) = 'T   '
! element(iv3d_s) = 'S   '               !(OCEAN)
! element(nv3d+iv2d_ssh) = 'SSH '        !(OCEAN)
! element(nv3d+iv2d_sst) = 'SST '        !(OCEAN)
! element(nv3d+iv2d_sss) = 'SSS '        !(OCEAN)
! if (DO_ALTIMETRY) then
!   element(nv3d+iv2d_eta) = 'eta '      !(OCEAN)
! endif
! if (DO_MLD) then
!   element(nv3d+iv2d_mld) = 'mld'       !(OCEAN)
! endif
! if (DO_DRIFTERS) then
!   element(nv3d+nv2d+iv4d_x) = 'X   '             !(OCEAN) (DRIFTERS)
!   element(nv3d+nv2d+iv4d_y) = 'Y   '             !(OCEAN) (DRIFTERS)
!   element(nv3d+nv2d+iv4d_z) = 'Z   '             !(OCEAN) (DRIFTERS)
! endif

  ! If using SLA assimilation calibrated to altimetry climatology:
  CHARACTER(13) :: SSHclm_file = 'not specified'

  !0001_nrt_20161228_000000_restart.nc
  CHARACTER(12) :: gridfile  = 'mesh_mask.nc'   !'MOM.res.nc'
  CHARACTER(12) :: gridfile1 = 'mesh_mask.nc'   !'MOM.res_1.nc'     (NEMO) n/a
  CHARACTER(14) :: gridfile2 = 'bathy_meter.nc' !'ocean_topog.nc'   (NEMO) n/a
  CHARACTER(14) :: gridfile3 = 'bathy_meter.nc' !'ocean_hgrid.nc'   (NEMO) n/a

  ! variable names in gridfile:
  CHARACTER(7) :: grid_lon_name = 'nav_lon'      ! (NEMO) n/a
  CHARACTER(7) :: grid_lat_name = 'nav_lat'      ! (NEMO) n/a
  CHARACTER(7) :: grid_lev_name = 'nav_lev'
  CHARACTER(2) :: grid_temp_name = 'tn'
  CHARACTER(2) :: grid_salt_name = 'sn'
  CHARACTER(2) :: grid_u_name = 'un'
  CHARACTER(2) :: grid_v_name = 'vn'
  CHARACTER(13):: grid_h_name = 'not specified'  ! (NEMO) n/a

  CHARACTER(7) :: grid_lon2d_name = 'nav_lon'
  CHARACTER(7) :: grid_lat2d_name = 'nav_lat'

  CHARACTER(5) :: grid_lsmask_name    = 'tmask'     ! land/sea mask (assuming >0 is ocean)
  CHARACTER(6) :: grid_depth_name  = 'mbathy'    ! (in mesh_mask file) !'Bathymetry'  ! In bathymetry file

  CHARACTER(3) :: grid_dx_name = 'e1t'
  CHARACTER(3) :: grid_dy_name = 'e2t'
  CHARACTER(5) :: grid_dz_name = 'e3t_0'

  ! Diagnostic file filenames 
  !STEVE: these aren't used for NEMOVAR, instead files are specified explicitly by obsop_xxx.f90 routine
  CHARACTER(13) :: diag_tsbase = 'not specified'   !(and u, and h)
  CHARACTER(13) :: diag_uvbase = 'not specified'   !(v and ave_ssh/sfc)
  CHARACTER(slen) :: diag_hbase
  ! variable names in diag file:
  CHARACTER(7) :: diag_lon_name = 'nav_lon'
  CHARACTER(7) :: diag_lat_name = 'nav_lat'
  CHARACTER(7) :: diag_lev_name = 'nav_lev'
  CHARACTER(2) :: diag_temp_name = 'tn'
  CHARACTER(2) :: diag_salt_name = 'sn'
  CHARACTER(2) :: diag_u_name = 'un'
  CHARACTER(2) :: diag_v_name = 'vn'
  CHARACTER(1) :: diag_h_name = 'h'
  CHARACTER(4) :: diag_ssh_name = 'sshn'
  CHARACTER(10):: diag_height_name = 'col_height'

  ! Restart filenames
  CHARACTER(slen) :: rsrt_tsbase = 'restart.nc' !'MOM.res.nc'   !(and u, and h)
  CHARACTER(slen) :: rsrt_uvbase = 'restart.nc' !(v and ave_ssh/sfc)
  CHARACTER(slen) :: rsrt_hbase  = 'restart.nc'
  ! Analysis variable names in restart file:
  CHARACTER(7) :: rsrt_lon_name = 'nav_lon'
  CHARACTER(7) :: rsrt_lat_name = 'nav_lat'
  CHARACTER(7) :: rsrt_lev_name = 'nav_lev'
  CHARACTER(2) :: rsrt_temp_name = 'tn'
  CHARACTER(2) :: rsrt_salt_name = 'tn'
  CHARACTER(2) :: rsrt_u_name = 'un'
  CHARACTER(2) :: rsrt_v_name = 'vn'
  CHARACTER(4) :: rsrt_ssh_name = 'sshn'
  CHARACTER(1) :: rsrt_h_name = 'h'              ! (NEMO) n/a

  !STEVE: unused:
  CHARACTER(slen) :: drbase   ! (DRIFTERS)   !(NEMO) n/a

  !STEVE: needed to read in ocean_hgrid.nc with supergrid format
  INTEGER :: nlon2d ! = 2*nlon  !STEVE: set below in initialize_params_model
  INTEGER :: nlat2d ! = 2*nlat  !STEVE: set below in initialize_params_model

  ! Bounds checking (for output by common_nemo.f90::write_restart)
  LOGICAL :: do_physlimit = .true.
  REAL(r_size) :: max_t = 40.0d0 ! ºC
  REAL(r_size) :: min_t = -4.0d0 ! ºC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu
  REAL(r_size) :: max_u = 99.0d0 ! m/s
  REAL(r_size) :: min_u =-99.0d0 ! m/s
  REAL(r_size) :: max_v = 99.0d0 ! m/s
  REAL(r_size) :: min_v =-99.0d0 ! m/s
! REAL(r_size) :: max_h = 90.0d3 ! =90000 m
! REAL(r_size) :: min_h =  0.0d0 ! m
  REAL(r_size) :: max_eta = 99.0d0 ! m
  REAL(r_size) :: min_eta =-99.0d0 ! m

  LOGICAL,SAVE :: params_model_initialized = .false.

CONTAINS


SUBROUTINE initialize_params_model
!===============================================================================
! Subroutine to initialize the parameters of the model
!===============================================================================

  IMPLICIT NONE

  if (params_model_initialized) then
    WRITE(6,*) "initialize_params_model:: already initialized, RETURNING..."
    RETURN
  endif

  nij0=nlon*nlat
  nlevall=nlev*nv3d+nv2d
  ngpv=nij0*nlevall

  !STEVE: needed to read in ocean_hgrid.nc with supergrid format for MOM6
  nlon2d = nlon !(NEMO) !2*nlon+1 !(MOM6)
  nlat2d = nlat !(NEMO) !2*nlat+1 !(MOM6)

  params_model_initialized = .true.

END SUBROUTINE initialize_params_model


END MODULE params_model
