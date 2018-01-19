MODULE params_model

USE common, ONLY: r_size, slen

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
!
  ! Now definable via namelist at runtime:

  ! HYCOM 1/12º global model grid dimensions:
  INTEGER :: glon=4500
  INTEGER :: glat=3298
  INTEGER :: glev=41

  ! HYCOM 1/12º local area model grid dimensions 
  ! (updated at runtime by grid specification file)
  INTEGER :: nlon
  INTEGER :: nlat
  INTEGER :: nlev
  ! LAM grid index:
  INTEGER :: istart,iend
  INTEGER :: jstart,jend

  !STEVE: going with z-coordinate in the vertical for now:
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s               !(OCEAN)(NCODA)
! INTEGER,PARAMETER :: nv3d=5 ! u,v,t,s,h             !(OCEAN)(MOM6)(HYCOM)
  INTEGER,PARAMETER :: nv4d=3 ! x,y,z                 !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
! INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss           !(OCEAN)
! INTEGER,PARAMETER :: nv2d=3 ! ssh,ubt,vbt           !(OCEAN) !(ALTIMETRY)(HYCOM)
  INTEGER,PARAMETER :: nv2d=1 ! ssh                   !(OCEAN) !(ALTIMETRY)(HYCOM)

  INTEGER,PARAMETER :: ilev_sfc=1

  ! Initialize via subroutine below:
  INTEGER,SAVE :: nij0           ! Number of gridpoints handled by myrank processor
  INTEGER,SAVE :: nlevall        ! Total number of variables and levels (3d + 2d)
  INTEGER,SAVE :: ngpv           ! Total number of gridpoints, including nij0*nlevall

  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                
  INTEGER,PARAMETER :: iv3d_h=5                !(NCODA)(n/a)

  INTEGER,PARAMETER :: iv2d_ssh=1              !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_ubt=2              !(OCEAN) ! Barotropic zonal velocity (HYCOM)
  INTEGER,PARAMETER :: iv2d_vbt=3              !(OCEAN) ! Barotropic meridional velocity (HYCOM)
  INTEGER,PARAMETER :: iv2d_sst=4              !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv2d_sss=5              !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv2d_eta=6              !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv2d_mld=7              !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv4d_x=1                !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv4d_y=2                !(NCODA)(n/a)
  INTEGER,PARAMETER :: iv4d_z=3                !(NCODA)(n/a)

  !
  ! Elements (labels for each model variable)
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

  ! In place of GFDL grid_spec.nc data file (e.g. as used for MOM6), used to define the model grid
  CHARACTER(8)  :: gridfile_lon  = 'glon.dat'
  CHARACTER(8)  :: gridfile_lat  = 'glat.dat'
  CHARACTER(8)  :: gridfile_lev  = 'glev.dat'
  CHARACTER(8)  :: gridfile_kmt  = 'kmt.dat'

  CHARACTER(14)   :: SSHclm_file = 'aEtaCds9399.nc'  ! Used for MOM4p1 GODAS (HYCOM n/a)
  CHARACTER(slen) :: rsrt_tbase = 't.dat'
  CHARACTER(slen) :: rsrt_sbase = 's.dat'
  CHARACTER(slen) :: rsrt_ubase = 'u.dat'
  CHARACTER(slen) :: rsrt_vbase = 'v.dat'
  CHARACTER(slen) :: rsrt_sshbase  = 'ssh.dat'

  ! Bounds checking (for output by common_hycom.f90::write_restart)
  LOGICAL :: do_physlimit=.true.
  REAL(r_size) :: max_t = 50.0d0 ! degC
  REAL(r_size) :: min_t = -4.0d0 ! degC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu
  REAL(r_size) :: max_uv = 10.0 ! m/s 
  REAL(r_size) :: min_uv =  -10.0d0 ! m/s 
  REAL(r_size) :: max_h = 1000.0d0 !m 
  REAL(r_size) :: min_h =  -1000.0d0 !m 
  REAL(r_size) :: max_ssh = 100.0d0 !m 
  REAL(r_size) :: min_ssh = -100.0d0 !m 
  
  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: ncundef = 1.267650600228229E+030 !1.0e18

  LOGICAL,SAVE :: params_model_initialized = .false.

  ! Write only the subgrid tile instead of updating the full global grid
  ! analysis file
  LOGICAL :: DO_WRITE_TILE = .true.
  INTEGER :: reclen_mult = 4

  !-----------------------------------------------------------------------------
  ! Not used for NRL NCODA-HYCOM version:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: sample_file_a = 'hycom_sample.a'
  CHARACTER(slen) :: sample_file_b = 'hycom_sample.b'
  CHARACTER(slen) :: grid_file_a   = 'regional.grid.a'
  CHARACTER(slen) :: grid_file_b   = 'regional.grid.b' !(n/a)

  !For input/output of model binary files (converted form HYCOM ab-format):
  CHARACTER(8) :: base  = '.fsd.bin' ! (HYCOM)(All variables assembled in one binary file)
  CHARACTER(2) :: base_a  = '.a'     ! (HYCOM)(All variables assembled in one file)
  CHARACTER(2) :: base_b  = '.b'     ! (HYCOM)(All variables assembled in one file)

  !model input and output file direct/sequential
  INTEGER, PARAMETER :: hycom_io_access = 1  ! 0 == direct, 1 == sequential
  

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
