MODULE params_model

USE common, ONLY: r_size, slen

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
!
  ! Now definable via namelist at runtime:

  ! HYCOM 1/4º global grid, 32 levels:
#ifdef DYNAMIC
  INTEGER :: nlon=563
  INTEGER :: nlat=351
  INTEGER :: nlev=32
#else
  INTEGER,PARAMETER :: nlon=563
  INTEGER,PARAMETER :: nlat=351
  INTEGER,PARAMETER :: nlev=32
#endif

  !INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
  INTEGER,PARAMETER :: nv3d=5 ! u,v,t,s,h             !(OCEAN)(MOM6)(HYCOM)
  INTEGER,PARAMETER :: nv4d=3 ! x,y,z                 !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
  !INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss          !(OCEAN)
  INTEGER,PARAMETER :: nv2d=3 ! ssh,ubt,vbt           !(OCEAN) !(ALTIMETRY)(HYCOM)

  INTEGER,PARAMETER :: ilev_sfc=1

  ! Initialize via subroutine below:
  INTEGER,SAVE :: nij0           ! Number of gridpoints handled by myrank processor
  INTEGER,SAVE :: nlevall        ! Total number of variables and levels (3d + 2d)
  INTEGER,SAVE :: ngpv           ! Total number of gridpoints, including nij0*nlevall

  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                !(OCEAN)
  INTEGER,PARAMETER :: iv3d_h=5                !(OCEAN)

  INTEGER,PARAMETER :: iv2d_ssh=1                    !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_ubt=2                    !(OCEAN) ! Barotropic zonal velocity (HYCOM)
  INTEGER,PARAMETER :: iv2d_vbt=3                    !(OCEAN) ! Barotropic meridional velocity (HYCOM)
  INTEGER,PARAMETER :: iv2d_sst=4                    !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model (MOM4p1)
  INTEGER,PARAMETER :: iv2d_sss=5                    !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models (MOM4p1)
  INTEGER,PARAMETER :: iv2d_eta=6              !(OCEAN) ! eta sea surface perturbation from mom4's ocean_barotropic.res.nc restart file
  INTEGER,PARAMETER :: iv2d_mld=7              !(OCEAN) ! mixed layer depth
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

  ! In place of mom4p1's grid_spec.nc data file, used to define the model grid
  CHARACTER(6)  :: gridfile  = '3zt.nc'
  CHARACTER(12) :: gridfile1 = 'blkdat.input'
  CHARACTER(6)  :: gridfile2 = 'bot.nc'
  CHARACTER(6)  :: gridfile3 = '3dh.nc'

  !For input/output of model binary files (converted form HYCOM ab-format):
  CHARACTER(8) :: base  = '.fsd.bin'! (HYCOM)(All variables assembled in one binary file)
  CHARACTER(2) :: base_a  = '.a'! (HYCOM)(All variables assembled in one
  CHARACTER(2) :: base_b  = '.b'! (HYCOM)(All variables assembled in one

  !model input and output file direct/sequential
  INTEGER, PARAMETER :: hycom_io_access = 1  ! 0 == direct, 1 == sequential

  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
! CHARACTER(32) :: ts_basefile = 'ocean_temp_salt.res.nc'
! CHARACTER(32) :: uv_basefile = 'ocean_velocity.res.nc'
! CHARACTER(32) :: sf_basefile = 'ocean_sbc.res.nc'
! CHARACTER(32) :: sh_basefile = 'ocean_barotropic.res.nc'
! CHARACTER(32) :: hs_basefile = 'ocean_TS.nc'
!
! ! For grid_spec.nc data file:
! CHARACTER(12) :: gridfile = 'grid_spec.nc'

  !For input/output of model binary files (converted form HYCOM ab-format):
! CHARACTER(slen) :: base  = 'fsd.bin'! (HYCOM)(All variables assembled in one binary file)
  !Otherwise, use netcdf files:
  CHARACTER(6) :: tbase = '3zt.nc' ! Temperature
  CHARACTER(6) :: sbase = '3zs.nc' ! Salinity
  CHARACTER(6) :: ubase = '3zu.nc' ! Zonal Current
  CHARACTER(6) :: vbase = '3zv.nc' ! Merional Current
  CHARACTER(6) :: hbase = '3dh.nc' ! Layer Thickness
  CHARACTER(6) :: ebase = 'fsd.nc' ! Sea Surface Height
  CHARACTER(6) :: drbase           ! (DRIFTERS)

  ! Bounds checking (for output by common_hycom.f90::write_restart)
  LOGICAL :: do_physlimit=.true.
  REAL(r_size) :: max_t = 50.0d0 ! degC
  REAL(r_size) :: min_t = 0.0d0 ! degC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu
  REAL(r_size) :: max_uv = 10.0 ! psu
  REAL(r_size) :: min_uv =  -10.0d0 ! psu
  REAL(r_size) :: max_h = 100000.0 ! psu
  REAL(r_size) :: min_h =  0.0d0 ! psu
  
  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: ncundef = 1.267650600228229E+030 !1.0e18

  LOGICAL,SAVE :: params_model_initialized = .false.

CONTAINS

SUBROUTINE initialize_params_model

  IMPLICIT NONE

  if (params_model_initialized) then
    WRITE(6,*) "initialize_params_model:: already initialized, RETURNING..."
    RETURN
  endif

  nij0=nlon*nlat
  nlevall=nlev*nv3d+nv2d
  ngpv=nij0*nlevall

  params_model_initialized = .true.

END SUBROUTINE initialize_params_model

END MODULE params_model
