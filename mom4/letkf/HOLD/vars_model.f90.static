MODULE vars_model

USE common,       ONLY: r_size
USE params_model, ONLY: nlon, nlat, nlev

IMPLICIT NONE

PUBLIC

REAL(r_size),DIMENSION(nlon),SAVE      :: lon !(nlon)
REAL(r_size),DIMENSION(nlat),SAVE      :: lat !(nlat)
REAL(r_size),DIMENSION(nlev),SAVE      :: lev !(nlev)                     !(OCEAN)
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: lon2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: lat2d !(nlon,nlat)              !(2DGRID)(TRIPOLAR)

REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dx !(nlon,nlat)
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: dy !(nlon,nlat)
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: phi0 !(nlon,nlat)
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: kmt0 !(nlon,nlat)               !(OCEAN)

!STEVE: For generalized grid
REAL(r_size),SAVE :: lon0, lonf, lat0, latf
REAL(r_size),SAVE :: wrapgap

!-----------------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------------
REAL(r_size),DIMENSION(nlat),SAVE      :: fcori !(nlat)
INTEGER,DIMENSION(nlon,nlat),SAVE      :: kmt            !(OCEAN) STEVE: the bottom topography for mom4
REAL(r_size),DIMENSION(nlon,nlat),SAVE :: SSHclm_m       !(OCEAN)(SLA) Stores model climatology to subtract from model eta_t when assimilating SLA
! For AMOC computation
REAL(r_size),DIMENSION(nlev),SAVE :: zb !(nlev)
REAL(r_size),DIMENSION(nlev),SAVE :: dz !(nlev)

CONTAINS

SUBROUTINE initialize_vars_model

USE common, ONLY: r_omega, pi
IMPLICIT NONE

  kmt = -1

  ! Corioris parameter
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)

  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)

  ! STEVE: for (more) generalized (longitude) grid:
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

END SUBROUTINE initialize_vars_model

END MODULE vars_model
