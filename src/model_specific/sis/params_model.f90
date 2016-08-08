MODULE params_model

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

  ! Now definable via namelist at runtime:

  ! MOM4 ncep2012 tripolar converted to spherical
#ifdef DYNAMIC
  INTEGER :: nlon=720
  INTEGER :: nlat=410
  INTEGER :: nlev=5
#else
  INTEGER,PARAMETER :: nlon=720
  INTEGER,PARAMETER :: nlat=410
  INTEGER,PARAMETER :: nlev=5   ! These are the different ice thicknesss categories:
                                ! (0-0.1, 0.1-0.3, 0.3-0.7, 0.7-1.1, and > 1.1 m)
#endif

  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=5
  INTEGER,PARAMETER :: nv4d=0 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
  INTEGER,PARAMETER :: nv2d=1 !STEVE: update when adding ui and vi -> nv2d=3
!
  ! Initialize via subroutine below:
  INTEGER,SAVE :: nij0           ! Number of gridpoints handled by myrank processor
  INTEGER,SAVE :: nlevall        ! Total number of variables and levels (3d + 2d)
  INTEGER,SAVE :: ngpv           ! Total number of gridpoints, including nij0*nlevall


! Placeholders needed by letkf_tools.f90 to compile:
  INTEGER,PARAMETER :: iv3d_t   = -1
  INTEGER,PARAMETER :: iv2d_mld = -1

!
! 3D Variables:
  INTEGER,PARAMETER :: iv3d_hs=1
  INTEGER,PARAMETER :: iv3d_hi=2
  INTEGER,PARAMETER :: iv3d_t1=3
  INTEGER,PARAMETER :: iv3d_t2=4
  INTEGER,PARAMETER :: iv3d_ps=5    ! 3D part size, from restart file levels 2-nlev+1
! 2D Variables:
  INTEGER,PARAMETER :: iv2d_cn=1    ! 2D ice concentration as measured
  INTEGER,PARAMETER :: iv2d_ui=2    ! 2D ice drift (zonal)
  INTEGER,PARAMETER :: iv2d_vi=3    ! 2D ice drift (meridional)

  REAL(r_size) :: obs_noise_coeff = 0.01  !(SIS)  !STEVE: can set in namelist

! INTEGER,PARAMETER :: iv4d_x=1                !(OCEAN) (DRIFTERS)
! INTEGER,PARAMETER :: iv4d_y=2                !(OCEAN) (DRIFTERS)
! INTEGER,PARAMETER :: iv4d_z=3                !(OCEAN) (DRIFTERS)

  !
  ! Elements
  !
  CHARACTER(20), SAVE :: element(nv3d+nv2d+nv4d)

  !For input/output model files:
  CHARACTER(16) :: basefile = 'ice_model.res.nc'

! CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
! CHARACTER(32) :: ts_basefile = 'ocean_temp_salt.res.nc'
! CHARACTER(32) :: uv_basefile = 'ocean_velocity.res.nc'
! CHARACTER(32) :: sf_basefile = 'ocean_sbc.res.nc'
! CHARACTER(32) :: sh_basefile = 'ocean_barotropic.res.nc'
! CHARACTER(32) :: hs_basefile = 'ocean_TS.nc'

  ! For grid_spec.nc data file:
  CHARACTER(12) :: gridfile = 'grid_spec.nc'

  ! Bounds checking (for output by common_mom4.f90::write_restart)
  LOGICAL :: do_physlimit=.true.
  REAL(r_size) :: max_t = 40.0d0 ! ÂC
  REAL(r_size) :: min_t = -4.0d0 ! ÂC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu

  !STEVE: needed for letkf.f90 to compile
  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
  
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

  !
  ! Elements
  !
  element(iv3d_hs) = 'snow layer height'
  element(iv3d_hi) = ' ice layer height'
  element(iv3d_t1) = ' ice layer 1 temp'
  element(iv3d_t2) = ' ice layer 2 temp'

  params_model_initialized = .true.

END SUBROUTINE initialize_params_model

END MODULE params_model
