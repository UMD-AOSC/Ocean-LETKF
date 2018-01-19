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
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lon2d !(nlon,nlat)              !(2DGRID)(For irregular grids)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lat2d !(nlon,nlat)              !(2DGRID)(For irregular grids)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lev2d !(nlon,nlat)              !(2DGRID)(For irregular grids)

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dx !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dy !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: phi0 !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)

  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE     :: fcori   !(nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE   :: fcori2d !(nlon,nlat)
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE        :: kmt            !(OCEAN) STEVE: the bottom topography for mom4
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE   :: SSHclm_m       !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
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
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: lev2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)

  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dx !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dy !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: phi0 !(nlon,nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)

  REAL(r_size),DIMENSION(nlat),SAVE      :: fcori !(nlat)
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: fcori2d !(nlon,nlat)
  INTEGER,DIMENSION(nlon,nlat),SAVE      :: kmt            !(OCEAN) STEVE: the bottom topography for mom4
  REAL(r_size),DIMENSION(nlon,nlat),SAVE :: SSHclm_m       !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
  ! For AMOC computation
  REAL(r_size),DIMENSION(nlev),SAVE :: zb !(nlev)
  REAL(r_size),DIMENSION(nlev),SAVE :: dz !(nlev)
!-----------------------------------------------------------------------------
#endif

  !STEVE: For generalized grid
  REAL(r_size),SAVE :: lon0, lonf, lat0, latf, lev0, levf
  REAL(r_size),SAVE :: wrapgap

  LOGICAL,SAVE :: vars_model_initialized = .false.
  LOGICAL,SAVE :: vars_model_set         = .false.

CONTAINS
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

SUBROUTINE initialize_vars_model

  !-----------------------------------------------------------------------------
  ! Initialize the model variables
  !-----------------------------------------------------------------------------
  USE common, ONLY: r_size
  USE params_model, ONLY: nlon, nlat, nlev, params_model_initialized

  IMPLICIT NONE

  if (.not. params_model_initialized) then
    WRITE(6,*) "vars_model.f90::initialize_vars_model::"
    WRITE(6,*) "ERROR: must call initialize_params_model before calling initialize_vars_model. EXITING..."
    STOP(85)
  elseif (vars_model_initialized) then
    WRITE(6,*) "initialize_vars_model:: already initialized, RETURNING..."
    RETURN
  else
    WRITE(6,*) "initialize_vars_model..."
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

  vars_model_initialized = .true.

END SUBROUTINE initialize_vars_model


SUBROUTINE set_vars_model

  USE common, ONLY: r_size, r_omega, pi
  USE params_model, ONLY: nlon, nlat, nlev, params_model_initialized

  IMPLICIT NONE

  REAL(r_size) :: lon0_360, lonf_360

  if (.not. vars_model_initialized) then
    WRITE(6,*) "Not initialized: vars_model_initialized, EXITING..." 
    STOP(86)
  else
    WRITE(6,*) "set_vars_model..."
  endif

  ! Corioris parameter
  fcori(:)   = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
  fcori2d(:,:) = 2.0d0 * r_omega * sin(lat2d(:,:)*pi/180.0d0)

  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  lev0 = lev(1)
  levf = lev(nlev)

  if (lon0 .eq. lonf) then
    WRITE(6,*) "ERROR in vars_model.f90::initialize_vars_model(), lon0==lonf==",lon0
    STOP(25)
  endif

  ! STEVE: for (more) generalized (longitude) grid:
  ! ISSUE: may need to check additional cases
  if (lon0 .ne. 0 .and. lonf .ne. 360) then
    lon0_360 = MODULO(lon0+360.0,360.0)
    lonf_360 = MODULO(lonf+360.0,360.0)
    wrapgap = abs(lon0_360 - lonf_360)
  else
    wrapgap=0.0
  endif

  vars_model_set = .true.

END SUBROUTINE set_vars_model

END MODULE vars_model
