MODULE vars_model

USE common,       ONLY: r_size

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lon !(nlon)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lat !(nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE   :: lev !(nlev)                     !(OCEAN)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lon2d !(nlon,nlat)      !(2DGRID)(For irregular grids)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lat2d !(nlon,nlat)      !(2DGRID)(For irregular grids)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: lev2d !(nlon,nlat)      !(2DGRID)(For irregular grids)

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dx !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: dy !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: area_t !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: phi0 !(nlon,nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:),SAVE :: height !(nlon,nlat,nlev)      !(OCEAN)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: wet !(nlon,nlat)

  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE     :: fcori   !(nlat)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE   :: fcori2d !(nlon,nlat)
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE        :: kmt     !(OCEAN) STEVE: the bottom topography for mom4
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE        :: pmsk    ! JILI HYCOM mask 
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE        :: umsk
  INTEGER,ALLOCATABLE,DIMENSION(:,:),SAVE        :: vmsk
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE   :: SSHclm_m  !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
  ! For AMOC computation
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE :: zb !(nlev)
  REAL(r_size),ALLOCATABLE,DIMENSION(:),SAVE :: dz !(nlev)
!-----------------------------------------------------------------------------

  ! For reading in NCODA z-levels:
  REAL(r_size), DIMENSION(100) :: z_lvl  ! STEVE: only use the first nlev levels

  !STEVE: For generalized grid
  REAL(r_size),SAVE :: lon0, lonf, lat0, latf
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
    STOP (85)
  elseif (vars_model_initialized) then
    WRITE(6,*) "initialize_vars_model:: already initialized, RETURNING..."
    RETURN
  else
    WRITE(6,*) "initialize_vars_model..."
  endif

  ALLOCATE(lon(nlon))
  ALLOCATE(lat(nlat))
  if (.not. allocated(lev)) then
    ALLOCATE(lev(nlev))
  endif
  ALLOCATE(lon2d(nlon,nlat))
  ALLOCATE(lat2d(nlon,nlat))

  ALLOCATE(dx(nlon,nlat))
  ALLOCATE(dy(nlon,nlat))
  ALLOCATE(phi0(nlon,nlat))
  ALLOCATE(kmt0(nlon,nlat))
  ALLOCATE(SSHclm_m(nlon,nlat))

  ALLOCATE(kmt(nlon,nlat))
  ALLOCATE(pmsk(nlon,nlat))
  ALLOCATE(umsk(nlon,nlat))
  ALLOCATE(vmsk(nlon,nlat))


  ALLOCATE(fcori(nlat))
  ALLOCATE(fcori2d(nlon,nlat))
  ALLOCATE(zb(nlev))
  ALLOCATE(dz(nlev))

  kmt = -1

  vars_model_initialized = .true.

END SUBROUTINE initialize_vars_model


SUBROUTINE set_vars_model

  USE common,       ONLY: r_size, r_omega, pi
  USE params_model, ONLY: nlon, nlat, nlev

  REAL(r_size) :: lon0_360, lonf_360
  INTEGER :: i,j,k
  LOGICAL :: ex
  LOGICAL :: dodebug = .true.

  if (.not. vars_model_initialized) then
    WRITE(6,*) "Not initialized: vars_model_initialized, EXITING..." 
    STOP (86)
  else
    WRITE(6,*) "set_vars_model..."
  endif

  ! Set lon and lat, just for completeness (NEMO)
  ! These are mostly still here for backward compatibility
  lon = lon2d(:,NINT(nlat/2.0d0))
  lat = lat2d(NINT(nlon/2.0d0),:)

  !-----------------------------------------------------------------------------
  ! Assign the model grid z-levels
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Assigning NCODA oanl namelist z_lvl levels to LETKF lev() array..."
  if (dodebug) then
    WRITE(6,*) "nlev = ", nlev
    WRITE(6,*) "SHAPE(lev) = ", SHAPE(lev)
    WRITE(6,*) "SIZE(lev) = ", SIZE(lev)
    WRITE(6,*) "lev = ", lev
    WRITE(6,*) "z_lvl = ", z_lvl
  endif
  if (allocated(lev)) then
    lev(1:nlev) = z_lvl(1:nlev)
  else
    STOP ("set_vars_model:: lev array is not allocated. EXITING...")
  endif

  !
  ! Coriolis parameter
  !
  WRITE(6,*) "Computing Coriolis parameter..."
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
  do j=1,nlat
    do i=1,nlon
      fcori2d(i,j) = 2.0d0 * r_omega * sin(lat2d(i,j)*pi/180.0d0)
    enddo
  enddo

  WRITE(6,*) "Setting lon0,lonf,lat0,latf..."
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)

  if (dodebug) then
    WRITE(6,*) "common_hycom.f90::set_common_hycom:: lon0,lonf,lat0,latf = ", lon0,lonf,lat0,latf
  endif

  if (lon0 .eq. lonf) then
    WRITE(6,*) "ERROR in vars_model.f90::initialize_vars_model(), lon0==lonf==",lon0
    STOP (25)
  endif

  ! STEVE: for (more) generalized (longitude) grid:
  ! ISSUE: may need to check additional cases
  WRITE(6,*) "Setting wrapgap..."
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
