MODULE vars_model

USE common,       ONLY: r_size

#ifndef DYNAMIC
USE params_model, ONLY: nlon, nlat, nlev
#endif

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
#ifdef DYNAMIC
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lon !(nlon)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lat !(nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lev !(nlev)                     !(OCEAN)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lon2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lat2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dx !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dy !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: phi0 !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)

  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: fcori !(nlat)
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE      :: kmt            !(OCEAN) STEVE: the bottom topography for mom4
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: SSHclm_m       !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
  ! For AMOC computation
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE :: zb !(nlev)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE :: dz !(nlev)
!-----------------------------------------------------------------------------
#else
  REAL(r_size),DIMENSION(nlon),SAVE      :: lon !(nlon)
  REAL(r_size),DIMENSION(nlat),SAVE      :: lat !(nlat)
  REAL(r_size),DIMENSION(nlev),SAVE      :: lev !(nlev)                     !(OCEAN)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: lon2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: lat2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)

  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dx !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dy !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: phi0 !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)

  REAL(r_size),DIMENSION(nlat),SAVE      :: fcori !(nlat)
  INTEGER,DIMENSION(nlon,nlat),SAVE      :: kmt            !(OCEAN) STEVE: the bottom topography for mom4
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: SSHclm_m       !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
  ! For AMOC computation
  REAL(r_size),DIMENSION(nlev),SAVE :: zb !(nlev)
  REAL(r_size),DIMENSION(nlev),SAVE :: dz !(nlev)
!-----------------------------------------------------------------------------
#endif

!STEVE: For generalized grid
REAL(r_size),SAVE :: lon0, lonf, lat0, latf
REAL(r_size),SAVE :: wrapgap

LOGICAL,SAVE :: vars_model_initialized = .false.

CONTAINS

SUBROUTINE initialize_vars_model

USE common, ONLY: r_omega, pi
USE params_model, ONLY: nlon, nlat, nlev, params_model_initialized

IMPLICIT NONE

  if (.not. params_model_initialized) then
    WRITE(6,*) "vars_model.f90::initialize_vars_model::"
    WRITE(6,*) "ERROR: must call initialize_params_model before calling initialize_vars_model. EXITING..."
    STOP(85)
  elseif (vars_model_initialized) then
    WRITE(6,*) "initialize_vars_model:: already initialized, RETURNING..."
    RETURN
  endif

#ifdef DYNAMIC
  ALLOCATE(lon(nlon))
  ALLOCATE(lat(nlat))
  ALLOCATE(lev(nlev))
  ALLOCATE(lon2d(nlon,nlev))
  ALLOCATE(lat2d(nlon,nlev))
  ALLOCATE(dx(nlon,nlat))
  ALLOCATE(dy(nlon,nlat))
  ALLOCATE(phi0(nlon,nlat))
  ALLOCATE(kmt0(nlon,nlat))
  ALLOCATE(SSHclm_m(nlon,nlat))

  ALLOCATE(kmt(nlon,nlat))

  ALLOCATE(fcori(nlat))
  ALLOCATE(zb(nlev))
  ALLOCATE(dz(nlev))
#endif

  kmt = -1

  ! Corioris parameter
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)

  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)

  ! STEVE: for (more) generalized (longitude) grid:
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

  vars_model_initialized = .true.

END SUBROUTINE initialize_vars_model

END MODULE vars_model
